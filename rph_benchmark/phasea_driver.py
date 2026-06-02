from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import yaml


DEFAULT_PROTOCOLS = ["ext", "full", "lite", "zero"]

STATUS_PENDING = "PENDING"
STATUS_COLLECTED = "COLLECTED"
STATUS_FAILED = "FAILED"


@dataclass(frozen=True)
class Paths:
    project_root: Path
    suite_name: str
    root_dir: Path
    plan_dir: Path
    configs_dir: Path
    runs_dir: Path
    logs_dir: Path
    cache_dir: Path
    reports_dir: Path


@dataclass(frozen=True)
class ResolvedTaskRoots:
    reaction_root: Path
    reaction_id: Optional[str]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def build_paths(project_root: Path, suite_name: str) -> Paths:
    project_root = Path(project_root)
    root_dir = project_root / suite_name
    return Paths(
        project_root=project_root,
        suite_name=suite_name,
        root_dir=root_dir,
        plan_dir=root_dir / "plan",
        configs_dir=root_dir / "configs",
        runs_dir=root_dir / "runs",
        logs_dir=root_dir / "logs",
        cache_dir=root_dir / "cache",
        reports_dir=root_dir / "reports",
    )


def ensure_plan_dirs(paths: Paths) -> None:
    for directory in (
        paths.root_dir,
        paths.plan_dir,
        paths.configs_dir,
        paths.runs_dir,
        paths.logs_dir,
        paths.cache_dir,
        paths.reports_dir,
    ):
        directory.mkdir(parents=True, exist_ok=True)


def panel_path(paths: Paths) -> Path:
    return paths.plan_dir / "panel.csv"


def task_manifest_path(paths: Paths) -> Path:
    return paths.plan_dir / "task_manifest.csv"


def rep_config_path(paths: Paths) -> Path:
    return paths.plan_dir / "rep_config.json"


def task_cache_dir(paths: Paths, task_id: str) -> Path:
    return paths.cache_dir / task_id


def write_json(path: Path, payload: Any) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
    return path


def write_csv(path: Path, rows: Sequence[Dict[str, Any]], fieldnames: Sequence[str]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})
    return path


def _read_json(path: Path, default: Any = None) -> Any:
    if not path.exists():
        return default
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def _read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def normalize_task_id(rx_id: str, protocol: str, replicate: int) -> str:
    return f"rx{str(rx_id).strip()}__{str(protocol).strip()}__rep{int(replicate):02d}"


def parse_task_id(task_id: str) -> Tuple[str, str, int]:
    parts = str(task_id).split("__")
    if len(parts) != 3 or not parts[0].startswith("rx") or not parts[2].startswith("rep"):
        raise ValueError(f"Invalid task_id: {task_id}")
    rx_id = parts[0][2:]
    protocol = parts[1]
    replicate = int(parts[2][3:])
    return rx_id, protocol, replicate


def artifact_candidates(run_dir: Path) -> Dict[str, Optional[Path]]:
    run_dir = Path(run_dir)
    reaction_roots = _candidate_reaction_roots(run_dir)
    reaction_root = reaction_roots[0] if reaction_roots else run_dir
    s1_dir = reaction_root / "S1_ConfGeneration"
    provenance_json = s1_dir / "provenance.json"
    product_candidates = [
        s1_dir / "product_min.xyz",
        s1_dir / "product" / "product_min.xyz",
    ]
    product_xyz = next((path for path in product_candidates if path.exists()), None)
    ensemble_path = s1_dir / "product" / "crest" / "ensemble.xyz"
    return {
        "reaction_root": reaction_root,
        "s1_dir": s1_dir if s1_dir.exists() else None,
        "provenance_json": provenance_json if provenance_json.exists() else None,
        "product_xyz": product_xyz,
        "ensemble_path": ensemble_path if ensemble_path.exists() else None,
    }


def validate_artifacts(run_dir: Path) -> Dict[str, bool]:
    artifacts = artifact_candidates(run_dir)
    return {
        "provenance_json": artifacts["provenance_json"] is not None,
        "product_min_xyz": artifacts["product_xyz"] is not None,
        "ensemble_found": artifacts["ensemble_path"] is not None,
    }


def artifacts_complete(run_dir: Path) -> bool:
    return all(validate_artifacts(run_dir).values())


def _is_safe_name(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text in {".", ".."}:
        return None
    if any(sep in text for sep in ("/", "\\")):
        return None
    return text


def _is_within(root: Path, candidate: Path) -> bool:
    try:
        root_resolved = root.resolve()
        candidate_resolved = candidate.resolve()
    except OSError:
        return False
    return root_resolved == candidate_resolved or root_resolved in candidate_resolved.parents


def _candidate_reaction_roots(base_root: Path) -> List[Path]:
    base_root = Path(base_root)
    if not base_root.exists() or not base_root.is_dir():
        return []

    candidates: List[Path] = []
    for child in sorted(base_root.iterdir(), key=lambda path: path.name):
        if not child.name.startswith("RXN_"):
            continue
        if not child.is_dir():
            continue
        if not _is_within(base_root, child):
            continue
        candidates.append(child.resolve())
    return candidates


def resolve_task_roots(paths: Paths, task_id: str) -> ResolvedTaskRoots:
    task_root = paths.runs_dir / task_id
    request = _read_json(task_cache_dir(paths, task_id) / "request.json", default={})
    requested_root_raw = request.get("reaction_root") if isinstance(request, dict) else None
    requested_root = Path(requested_root_raw) if requested_root_raw else task_root
    if not requested_root.is_absolute():
        requested_root = task_root / requested_root
    if not requested_root.exists():
        requested_root = task_root

    reaction_id = _is_safe_name(request.get("reaction_id") if isinstance(request, dict) else None)
    if reaction_id is not None:
        direct_candidate = requested_root / reaction_id
        if direct_candidate.exists() and _is_within(requested_root, direct_candidate):
            return ResolvedTaskRoots(reaction_root=direct_candidate.resolve(), reaction_id=reaction_id)

    nested_roots = _candidate_reaction_roots(requested_root)
    if nested_roots:
        nested_root = nested_roots[0]
        return ResolvedTaskRoots(reaction_root=nested_root, reaction_id=nested_root.name)

    return ResolvedTaskRoots(reaction_root=requested_root.resolve(), reaction_id=None)


def collect_task(paths: Paths, task_id: str) -> Dict[str, Any]:
    cache_dir = task_cache_dir(paths, task_id)
    metrics_path = cache_dir / "metrics.json"
    done_path = cache_dir / "done.json"
    failed_path = cache_dir / "failed.json"
    resolved = resolve_task_roots(paths, task_id)
    artifacts = artifact_candidates(resolved.reaction_root)

    result: Dict[str, Any] = {
        "task_id": task_id,
        "run_dir": str(paths.runs_dir / task_id),
        "reaction_root": str(resolved.reaction_root),
        "reaction_id": resolved.reaction_id,
        "provenance_json": str(artifacts["provenance_json"]) if artifacts["provenance_json"] else None,
        "product_min_xyz": str(artifacts["product_xyz"]) if artifacts["product_xyz"] else None,
        "ensemble_path": str(artifacts["ensemble_path"]) if artifacts["ensemble_path"] else None,
        "artifact_check": validate_artifacts(resolved.reaction_root),
        "metrics_path": str(metrics_path) if metrics_path.exists() else None,
        "error_message": None,
        "run_seconds": None,
    }

    if metrics_path.exists():
        metrics_payload = _read_json(metrics_path, default={})
        if isinstance(metrics_payload, dict):
            result.update(metrics_payload)

    if done_path.exists():
        done_payload = _read_json(done_path, default={})
        if isinstance(done_payload, dict):
            result["run_seconds"] = done_payload.get("run_seconds")

    if failed_path.exists():
        failed_payload = _read_json(failed_path, default={})
        if isinstance(failed_payload, dict):
            result["run_seconds"] = failed_payload.get("run_seconds", result["run_seconds"])
            result["error_message"] = failed_payload.get("error_message")

    return result


def task_status(paths: Paths, task_id: str) -> str:
    cache_dir = task_cache_dir(paths, task_id)
    if (cache_dir / "failed.json").exists():
        return STATUS_FAILED
    if (cache_dir / "done.json").exists():
        return STATUS_COLLECTED
    return STATUS_PENDING


def run_task(_paths: Paths, task_id: str, preclaimed: bool = False) -> int:
    raise RuntimeError(f"run_task is not implemented for {task_id} (preclaimed={preclaimed})")


def load_defaults(paths: Paths) -> Dict[str, Any]:
    defaults_path = paths.project_root / "config" / "defaults.yaml"
    if not defaults_path.exists():
        return {"run": {}, "step1": {}}
    with open(defaults_path, encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return data if isinstance(data, dict) else {"run": {}, "step1": {}}


def render_task_config(paths: Paths, task_id: str, rx_id: str, protocol: str) -> Path:
    defaults = dict(load_defaults(paths))
    run_cfg = dict(defaults.get("run", {}) or {})
    run_cfg["output_root"] = str(paths.runs_dir / task_id)
    single_cfg = dict(run_cfg.get("single", {}) or {})
    single_cfg["rx_id"] = rx_id
    run_cfg["single"] = single_cfg
    defaults["run"] = run_cfg

    step1_cfg = dict(defaults.get("step1", {}) or {})
    step1_cfg["protocol"] = protocol
    defaults["step1"] = step1_cfg

    config_path = paths.configs_dir / f"{task_id}.yaml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(yaml.safe_dump(defaults, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return config_path


def _extract_metrics_payload(task_metrics: Dict[str, Any]) -> Dict[str, Any]:
    excluded = {
        "task_id",
        "run_dir",
        "reaction_root",
        "reaction_id",
        "provenance_json",
        "product_min_xyz",
        "ensemble_path",
        "artifact_check",
        "metrics_path",
        "error_message",
        "run_seconds",
    }
    return {key: value for key, value in task_metrics.items() if key not in excluded}


def _comparison_to_ext(metrics: Dict[str, Any], ext_metrics: Dict[str, Any]) -> Dict[str, Optional[float]]:
    current = metrics.get("mean_native_energy")
    anchor = ext_metrics.get("mean_native_energy")
    if current is None or anchor is None:
        return {"delta_mean_native_energy": None}
    return {"delta_mean_native_energy": float(current) - float(anchor)}


def build_rx_suite_report(
    paths: Paths,
    rx_id: str,
    replicate: int,
    protocol_results: Sequence[Dict[str, Any]],
) -> Dict[str, Any]:
    report_basename = f"rx{rx_id}__rep{int(replicate):02d}_suite_report"
    json_path = paths.reports_dir / f"{report_basename}.json"
    markdown_path = paths.reports_dir / f"{report_basename}.md"
    paths.reports_dir.mkdir(parents=True, exist_ok=True)

    ext_metrics = {}
    for item in protocol_results:
        if item.get("protocol") == "ext" and isinstance(item.get("metrics"), dict):
            ext_metrics = dict(item["metrics"])
            break

    protocol_entries: List[Dict[str, Any]] = []
    for item in protocol_results:
        entry = dict(item)
        metrics = dict(entry.get("metrics") or {})
        entry["metrics"] = metrics
        entry["comparison_to_ext"] = _comparison_to_ext(metrics, ext_metrics)
        protocol_entries.append(entry)

    completed_protocols = sum(1 for item in protocol_entries if int(item.get("exit_code", 1)) == 0)
    failed_protocols = len(protocol_entries) - completed_protocols
    suite_status = "DONE" if failed_protocols == 0 else "FAILED"

    report = {
        "rx_id": str(rx_id),
        "replicate": int(replicate),
        "suite_status": suite_status,
        "protocol_order": [str(item.get("protocol")) for item in protocol_entries],
        "completed_protocols": completed_protocols,
        "failed_protocols": failed_protocols,
        "protocol_results": protocol_entries,
        "report_paths": {
            "json": str(json_path),
            "markdown": str(markdown_path),
        },
    }
    write_json(json_path, report)

    markdown_lines = [
        f"# Phase A Per-Rx Suite Report: rx{rx_id} rep{int(replicate):02d}",
        "",
        f"Suite status: **{suite_status}**",
        "",
    ]
    for item in protocol_entries:
        markdown_lines.extend(
            [
                f"### {item.get('protocol')}",
                f"- task_id: {item.get('task_id')}",
                f"- exit_code: {item.get('exit_code')}",
                f"- task_status: {item.get('task_status')}",
                f"- run_seconds: {item.get('run_seconds')}",
                f"- metrics_path: {item.get('metrics_path')}",
                f"- log_file: {item.get('log_file')}",
            ]
        )
        if item.get("error_message"):
            markdown_lines.append(f"- error: {item['error_message']}")
        for key, value in dict(item.get("metrics") or {}).items():
            markdown_lines.append(f"- {key}: {value}")
        for key, value in dict(item.get("comparison_to_ext") or {}).items():
            markdown_lines.append(f"- {key}_vs_ext: {value}")
        markdown_lines.append("")
    markdown_path.write_text("\n".join(markdown_lines).rstrip() + "\n", encoding="utf-8")
    return report


def cmd_run_rx(paths: Paths, args: argparse.Namespace) -> int:
    rx_id = str(args.rx_id)
    replicate = int(args.replicate)
    rows = _read_csv_rows(task_manifest_path(paths))
    selected_rows = [
        row for row in rows if row.get("rx_id") == rx_id and int(row.get("replicate", "0") or 0) == replicate
    ]
    if not selected_rows:
        raise RuntimeError(f"missing task_manifest rows for rx_id={rx_id} replicate={replicate}")

    order_index = {protocol: index for index, protocol in enumerate(DEFAULT_PROTOCOLS)}
    selected_rows.sort(key=lambda row: order_index.get(str(row.get("protocol")), len(DEFAULT_PROTOCOLS)))

    protocol_results: List[Dict[str, Any]] = []
    for row in selected_rows:
        protocol = str(row.get("protocol"))
        task_id = str(row.get("task_id") or normalize_task_id(rx_id, protocol, replicate))
        exit_code = int(run_task(paths, task_id, preclaimed=False))
        task_metrics = collect_task(paths, task_id)
        result = {
            "protocol": protocol,
            "task_id": task_id,
            "exit_code": exit_code,
            "task_status": task_status(paths, task_id),
            "run_seconds": task_metrics.get("run_seconds"),
            "metrics_path": task_metrics.get("metrics_path"),
            "log_file": row.get("log_relpath") or f"logs/{task_id}.log",
            "error_message": task_metrics.get("error_message"),
            "metrics": _extract_metrics_payload(task_metrics) if task_metrics.get("metrics_path") else {},
        }
        protocol_results.append(result)

    build_rx_suite_report(paths, rx_id, replicate, protocol_results)
    return 0 if all(int(item["exit_code"]) == 0 for item in protocol_results) else 1
