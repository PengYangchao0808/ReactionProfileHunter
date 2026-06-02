from __future__ import annotations
# pyright: reportDeprecated=false, reportUnusedCallResult=false, reportAny=false, reportExplicitAny=false, reportUnknownMemberType=false, reportUnknownVariableType=false, reportUnknownArgumentType=false, reportReturnType=false, reportUnnecessaryComparison=false, reportUnreachable=false, reportUnnecessaryIsInstance=false

import hashlib
import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Union

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from benchmark.dft_theory.lib.state_manager import (
    create_session_manifest,
    init_rx_manifest,
    read_rx_manifest,
    resolve_upstream_confsearch_dir,
    set_stage_status,
    write_rx_manifest,
)

KNOWN_PROTOCOLS = frozenset({"ext", "full", "lite", "zero"})
SOURCE_MANIFEST_SCHEMA_VERSION = "dft_benchmark_confsearch_import_v2"
PathLike = Union[str, Path]
ManifestUpdater = Callable[[Path, str, Dict[str, Any]], None]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return "sha256:" + digest.hexdigest()


def _resolve_input_path(value: PathLike) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return path.resolve()


def _resolve_case_path(value: PathLike, confsearch_root: Optional[PathLike]) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path.resolve()
    project_candidate = (PROJECT_ROOT / path).resolve()
    if confsearch_root is None or project_candidate.exists():
        return project_candidate
    return (_resolve_input_path(confsearch_root) / path).resolve()


def _build_case_index(cases: Optional[Sequence[Dict[str, Any]]]) -> Dict[str, Dict[str, Any]]:
    index: Dict[str, Dict[str, Any]] = {}
    for case in cases or []:
        if not isinstance(case, dict):
            continue
        rx_id = str(case.get("rx_id", "")).strip()
        if rx_id:
            index[rx_id] = dict(case)
    return index


def _select_seed_protocol(rx_id: str, case: Dict[str, Any], seed_protocol: Optional[str]) -> str:
    selected = seed_protocol
    if selected is None or not str(selected).strip():
        selected = case.get("seed_protocol")
    value = str(selected).strip() if selected is not None else ""
    if not value:
        raise ValueError("Missing seed_protocol for rx{0}".format(rx_id))
    return value


def _resolve_source_base(
    rx_id: str,
    case: Dict[str, Any],
    confsearch_root: Optional[PathLike],
) -> Dict[str, Optional[str]]:
    raw_case_source = case.get("confsearch_source")
    if raw_case_source is not None and str(raw_case_source).strip():
        resolved = _resolve_case_path(str(raw_case_source), confsearch_root)
        return {
            "case_confsearch_source": str(raw_case_source),
            "confsearch_root": str(_resolve_input_path(confsearch_root)) if confsearch_root is not None else None,
            "resolved_source_root": str(resolved),
        }
    if confsearch_root is None:
        raise ValueError("Missing confsearch source for rx{0}".format(rx_id))
    root = _resolve_input_path(confsearch_root)
    resolved_root = (root / ("rx{0}".format(rx_id))).resolve()
    return {
        "case_confsearch_source": None,
        "confsearch_root": str(root),
        "resolved_source_root": str(resolved_root),
    }


def _choose_run_dir(protocol_dir: Path) -> Path:
    direct_s1 = protocol_dir / "S1_ConfGeneration"
    if direct_s1.is_dir():
        return protocol_dir.resolve()
    candidates: List[Path] = []
    for child in sorted(protocol_dir.iterdir(), key=lambda item: item.name):
        if child.is_dir() and (child / "S1_ConfGeneration").is_dir():
            candidates.append(child.resolve())
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError("No S1_ConfGeneration subtree found under {0}".format(protocol_dir))
    raise ValueError("Ambiguous confsearch run directory under {0}".format(protocol_dir))


def _first_existing_file(candidates: Sequence[Path]) -> Optional[Path]:
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return None


def _resolve_source_layout(source_root: Path, selected_protocol: str) -> Dict[str, Optional[Path]]:
    resolved_source_root = source_root.resolve()
    protocol_dir: Optional[Path] = None
    run_dir: Optional[Path] = None
    s0_dir: Optional[Path] = None
    s1_dir: Optional[Path] = None
    rx_root: Optional[Path] = None

    if resolved_source_root.name == "S1_ConfGeneration":
        s1_dir = resolved_source_root
        run_dir = s1_dir.parent.resolve()
        protocol_dir = run_dir.parent.resolve()
        rx_root = protocol_dir.parent.resolve()
    elif (resolved_source_root / selected_protocol).is_dir():
        rx_root = resolved_source_root
        protocol_dir = (resolved_source_root / selected_protocol).resolve()
        run_dir = _choose_run_dir(protocol_dir)
        s1_dir = (run_dir / "S1_ConfGeneration").resolve()
    elif (resolved_source_root / "S1_ConfGeneration").is_dir():
        run_dir = resolved_source_root
        s1_dir = (run_dir / "S1_ConfGeneration").resolve()
        if run_dir.name in KNOWN_PROTOCOLS:
            protocol_dir = run_dir
            rx_root = protocol_dir.parent.resolve()
        elif run_dir.parent.name in KNOWN_PROTOCOLS:
            protocol_dir = run_dir.parent.resolve()
            rx_root = protocol_dir.parent.resolve()
        else:
            protocol_dir = run_dir
            rx_root = run_dir.parent.resolve()
    else:
        protocol_dir = resolved_source_root
        if protocol_dir.name in KNOWN_PROTOCOLS and protocol_dir.name != selected_protocol:
            raise ValueError(
                "Requested seed_protocol '{0}' does not match explicit source protocol '{1}'".format(
                    selected_protocol,
                    protocol_dir.name,
                )
            )
        run_dir = _choose_run_dir(protocol_dir)
        s1_dir = (run_dir / "S1_ConfGeneration").resolve()
        rx_root = protocol_dir.parent.resolve()

    if protocol_dir is None or run_dir is None or s1_dir is None or rx_root is None:
        raise RuntimeError("Failed to resolve confsearch source layout from {0}".format(source_root))
    s0_candidate = (run_dir / "S0_Mechanism").resolve()
    if s0_candidate.is_dir():
        s0_dir = s0_candidate
    if not s1_dir.is_dir():
        raise FileNotFoundError("Resolved S1 directory does not exist: {0}".format(s1_dir))

    evaluation_dir: Optional[Path] = None
    for candidate in (rx_root / "evaluation", protocol_dir / "evaluation", run_dir / "evaluation"):
        if candidate.is_dir():
            evaluation_dir = candidate.resolve()
            break

    return {
        "resolved_source_root": resolved_source_root,
        "rx_root": rx_root,
        "protocol_dir": protocol_dir,
        "run_dir": run_dir,
        "s0_dir": s0_dir,
        "s1_dir": s1_dir,
        "evaluation_dir": evaluation_dir,
    }


def _copy_file(source: Path, destination: Path, label: str, required: bool) -> Dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(source), str(destination))
    return {
        "destination": str(destination.resolve()),
        "label": label,
        "required": required,
        "sha256": _sha256(destination),
        "size_bytes": destination.stat().st_size,
        "source": str(source.resolve()),
    }


def _local_file_entry(path: Path, label: str, required: bool) -> Dict[str, Any]:
    resolved = path.resolve()
    return {
        "destination": str(resolved),
        "label": label,
        "required": required,
        "sha256": _sha256(resolved),
        "size_bytes": resolved.stat().st_size,
        "source": str(resolved),
    }


def _direct_update_rx_manifest(session_dir: Path, rx_id: str, record: Dict[str, Any]) -> None:
    init_rx_manifest(session_dir, rx_id)
    payload = read_rx_manifest(session_dir, rx_id)

    imported_files = record.get("imported_files", [])
    imported_by_label = {
        str(item.get("label")): item
        for item in imported_files
        if isinstance(item, dict) and item.get("label")
    }

    payload.setdefault("inputs", {})
    payload.setdefault("warnings", [])
    payload.setdefault("upstream", {})
    payload["upstream"]["confsearch"] = {
        "imported_at": record["imported_at"],
        "imported_files": imported_files,
        "missing_optional_files": record.get("missing_optional_files", []),
        "s0_dir": record["s0_dir"],
        "s1_dir": record["s1_dir"],
        "seed_protocol": record["selected_seed_protocol"],
        "source_manifest": record["source_manifest"],
        "source_paths": record["source_paths"],
        "status": "imported",
    }

    product_entry = imported_by_label.get("product_min_xyz")
    if isinstance(product_entry, dict):
        payload["inputs"]["product_xyz"] = {
            "hash": product_entry.get("sha256"),
            "path": product_entry.get("destination"),
            "source_step": "upstream_confsearch",
        }

    provenance_entry = imported_by_label.get("provenance_json")
    if isinstance(provenance_entry, dict):
        payload["inputs"]["confsearch_provenance"] = {
            "hash": provenance_entry.get("sha256"),
            "path": provenance_entry.get("destination"),
            "source_step": "upstream_confsearch",
        }

    source_manifest_entry = imported_by_label.get("source_manifest")
    if isinstance(source_manifest_entry, dict):
        payload["inputs"]["confsearch_source_manifest"] = {
            "hash": source_manifest_entry.get("sha256"),
            "path": source_manifest_entry.get("destination"),
            "source_step": "upstream_confsearch",
        }

    mechanism_summary_entry = imported_by_label.get("mechanism_summary_json")
    if isinstance(mechanism_summary_entry, dict):
        payload["inputs"]["mechanism_summary"] = {
            "hash": mechanism_summary_entry.get("sha256"),
            "path": mechanism_summary_entry.get("destination"),
            "source_step": "upstream_confsearch",
        }

    mechanism_graph_entry = imported_by_label.get("mechanism_graph_json")
    if isinstance(mechanism_graph_entry, dict):
        payload["inputs"]["mechanism_graph"] = {
            "hash": mechanism_graph_entry.get("sha256"),
            "path": mechanism_graph_entry.get("destination"),
            "source_step": "upstream_confsearch",
        }

    s0_status_entry = imported_by_label.get("s0_status_json")
    if isinstance(s0_status_entry, dict):
        payload["inputs"]["s0_status"] = {
            "hash": s0_status_entry.get("sha256"),
            "path": s0_status_entry.get("destination"),
            "source_step": "upstream_confsearch",
        }

    precursor_entry = imported_by_label.get("precursor_min_xyz")
    if isinstance(precursor_entry, dict):
        payload["inputs"]["precursor_xyz"] = {
            "hash": precursor_entry.get("sha256"),
            "path": precursor_entry.get("destination"),
            "source_step": "upstream_confsearch",
        }

    payload["updated_at"] = record["imported_at"]
    write_rx_manifest(session_dir, rx_id, payload)


def _update_rx_manifest(
    session_dir: Path,
    rx_id: str,
    record: Dict[str, Any],
    manifest_updater: Optional[ManifestUpdater],
) -> None:
    if manifest_updater is not None:
        manifest_updater(session_dir, rx_id, record)
        return
    _direct_update_rx_manifest(session_dir, rx_id, record)


def _build_copy_plan(layout: Dict[str, Optional[Path]], upstream_dir: Path) -> Dict[str, Any]:
    s0_dir = layout["s0_dir"]
    s1_dir = layout["s1_dir"]
    evaluation_dir = layout["evaluation_dir"]
    if s0_dir is None:
        raise FileNotFoundError("Missing resolved S0_Mechanism directory")
    if s1_dir is None:
        raise RuntimeError("Missing resolved S1 directory")

    destination_s0 = upstream_dir / "S0_Mechanism"
    destination_s1 = upstream_dir / "S1_ConfGeneration"
    product_dir = s1_dir / "product"
    precursor_dir = s1_dir / "precursor"

    required_specs = [
        (
            "mechanism_summary_json",
            [s0_dir / "mechanism_summary.json"],
            destination_s0 / "mechanism_summary.json",
        ),
        (
            "mechanism_graph_json",
            [s0_dir / "mechanism_graph.json"],
            destination_s0 / "mechanism_graph.json",
        ),
        (
            "product_min_xyz",
            [
                s1_dir / "product_min.xyz",
                product_dir / "product_min.xyz",
                product_dir / "product_global_min.xyz",
            ],
            destination_s1 / "product_min.xyz",
        ),
        (
            "provenance_json",
            [s1_dir / "provenance.json"],
            destination_s1 / "provenance.json",
        ),
        (
            "precursor_min_xyz",
            [
                s1_dir / "precursor_min.xyz",
                precursor_dir / "precursor_min.xyz",
                precursor_dir / "precursor_global_min.xyz",
            ],
            destination_s1 / "precursor" / "precursor_min.xyz",
        ),
    ]

    optional_specs = [
        (
            "s0_status_json",
            [s0_dir / "s0_status.json"],
            destination_s0 / "s0_status.json",
        ),
        (
            "conformer_energies_json",
            [
                product_dir / "finalDFT" / "conformer_energies.json",
                product_dir / "dft" / "conformer_energies.json",
                product_dir / "conformer_energies.json",
            ],
            destination_s1 / "product" / "conformer_energies.json",
        ),
        (
            "conformer_thermo_csv",
            [
                product_dir / "finalDFT" / "conformer_thermo.csv",
                product_dir / "dft" / "conformer_thermo.csv",
                product_dir / "conformer_thermo.csv",
            ],
            destination_s1 / "product" / "conformer_thermo.csv",
        ),
    ]

    if evaluation_dir is not None:
        optional_specs.extend(
            [
                (
                    "evaluation_summary_json",
                    [evaluation_dir / "summary.json"],
                    upstream_dir / "evaluation" / "summary.json",
                ),
                (
                    "evaluation_summary_md",
                    [evaluation_dir / "summary.md"],
                    upstream_dir / "evaluation" / "summary.md",
                ),
            ]
        )
    return {
        "required_specs": required_specs,
        "optional_specs": optional_specs,
    }


def _import_one_rx(
    session_dir: Path,
    rx_id: str,
    case: Dict[str, Any],
    confsearch_root: Optional[PathLike],
    seed_protocol: Optional[str],
    manifest_updater: Optional[ManifestUpdater],
) -> Dict[str, Any]:
    selected_seed_protocol = _select_seed_protocol(rx_id, case, seed_protocol)
    source_base_info = _resolve_source_base(rx_id, case, confsearch_root)
    source_root = Path(str(source_base_info["resolved_source_root"]))
    layout = _resolve_source_layout(source_root, selected_seed_protocol)

    upstream_dir = resolve_upstream_confsearch_dir(session_dir, rx_id).resolve()
    if upstream_dir.exists():
        shutil.rmtree(str(upstream_dir))
    upstream_dir.mkdir(parents=True, exist_ok=True)

    plan = _build_copy_plan(layout, upstream_dir)
    imported_files: List[Dict[str, Any]] = []
    missing_optional_files: List[str] = []

    for label, candidates, destination in plan["required_specs"]:
        source = _first_existing_file(candidates)
        if source is None:
            raise FileNotFoundError("Missing required confsearch seed artifact '{0}' for rx{1}".format(label, rx_id))
        imported_files.append(_copy_file(source, destination, label, True))

    for label, candidates, destination in plan["optional_specs"]:
        source = _first_existing_file(candidates)
        if source is None:
            missing_optional_files.append(label)
            continue
        imported_files.append(_copy_file(source, destination, label, False))

    imported_at = now_iso()
    source_paths = {
        "case_confsearch_source": source_base_info["case_confsearch_source"],
        "confsearch_root": source_base_info["confsearch_root"],
        "evaluation_dir": str(layout["evaluation_dir"]) if layout["evaluation_dir"] is not None else None,
        "protocol_dir": str(layout["protocol_dir"]),
        "resolved_source_root": str(layout["resolved_source_root"]),
        "run_dir": str(layout["run_dir"]),
        "rx_root": str(layout["rx_root"]),
        "s0_dir": str(layout["s0_dir"]) if layout["s0_dir"] is not None else None,
        "s1_dir": str(layout["s1_dir"]),
    }

    source_manifest_path = upstream_dir / "source_manifest.json"
    source_manifest_payload = {
        "imported_at": imported_at,
        "imported_files": imported_files,
        "missing_optional_files": missing_optional_files,
        "rx_id": rx_id,
        "schema_version": SOURCE_MANIFEST_SCHEMA_VERSION,
        "selected_seed_protocol": selected_seed_protocol,
        "source_paths": source_paths,
    }
    _write_json(source_manifest_path, source_manifest_payload)
    source_manifest_entry = _local_file_entry(source_manifest_path, "source_manifest", True)

    record = {
        "imported_at": imported_at,
        "imported_files": sorted(imported_files + [source_manifest_entry], key=lambda item: str(item["label"])),
        "missing_optional_files": missing_optional_files,
        "rx_id": rx_id,
        "s0_dir": str((upstream_dir / "S0_Mechanism").resolve()),
        "s1_dir": str((upstream_dir / "S1_ConfGeneration").resolve()),
        "selected_seed_protocol": selected_seed_protocol,
        "source_manifest": str(source_manifest_path.resolve()),
        "source_paths": source_paths,
        "upstream_root": str(upstream_dir),
    }
    _update_rx_manifest(session_dir, rx_id, record, manifest_updater)
    return record


def migrate_confsearch_seeds(
    session_dir: PathLike,
    rx_ids: Sequence[Union[str, int]],
    cases: Optional[Sequence[Dict[str, Any]]] = None,
    confsearch_root: Optional[PathLike] = None,
    seed_protocol: Optional[str] = None,
    manifest_updater: Optional[ManifestUpdater] = None,
) -> Dict[str, Dict[str, Any]]:
    resolved_session_dir = _resolve_input_path(session_dir)
    resolved_session_dir.mkdir(parents=True, exist_ok=True)
    cases_by_rx = _build_case_index(cases)
    create_session_manifest(
        resolved_session_dir,
        cases=[str(item.get("rx_id", "")).strip() for item in (cases or []) if isinstance(item, dict)],
        config_snapshot={"confsearch_root": str(_resolve_input_path(confsearch_root)) if confsearch_root is not None else None},
    )
    set_stage_status(resolved_session_dir, "migrate", "running")
    results: Dict[str, Dict[str, Any]] = {}
    for raw_rx_id in rx_ids:
        rx_id = str(raw_rx_id).strip()
        if not rx_id:
            raise ValueError("Encountered empty rx_id during confsearch seed migration")
        case = cases_by_rx.get(rx_id, {})
        results[rx_id] = _import_one_rx(
            session_dir=resolved_session_dir,
            rx_id=rx_id,
            case=case,
            confsearch_root=confsearch_root,
            seed_protocol=seed_protocol,
            manifest_updater=manifest_updater,
        )
    set_stage_status(resolved_session_dir, "migrate", "done")
    return results


def run_migrate_stage(
    session_dir: PathLike,
    rx_ids: Sequence[Union[str, int]],
    cases: Optional[Sequence[Dict[str, Any]]] = None,
    confsearch_root: Optional[PathLike] = None,
    seed_protocol: Optional[str] = None,
    manifest_updater: Optional[ManifestUpdater] = None,
) -> Dict[str, Dict[str, Any]]:
    return migrate_confsearch_seeds(
        session_dir=session_dir,
        rx_ids=rx_ids,
        cases=cases,
        confsearch_root=confsearch_root,
        seed_protocol=seed_protocol,
        manifest_updater=manifest_updater,
    )


def _load_cases_from_yaml(path: Path) -> List[Dict[str, Any]]:
    import yaml

    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    raw_cases = payload.get("cases", [])
    return [dict(item) for item in raw_cases if isinstance(item, dict)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Import confsearch seed artifacts into a DFT benchmark session.")
    parser.add_argument("--session-dir", required=True, help="Benchmark session directory")
    parser.add_argument("--rx-id", nargs="*", default=None, help="Optional subset of rx ids")
    parser.add_argument("--cases-config", default=str(PROJECT_ROOT / "benchmark" / "dft_theory" / "config" / "benchmark_cases.yaml"), help="Benchmark cases YAML")
    parser.add_argument("--confsearch-root", default=None, help="Optional confsearch root override")
    parser.add_argument("--seed-protocol", default=None, help="Optional protocol override for all selected rx ids")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cases = _load_cases_from_yaml(Path(args.cases_config).expanduser().resolve())
    if args.rx_id:
        selected = {str(item).strip() for item in args.rx_id if str(item).strip()}
        rx_ids = [str(item.get("rx_id", "")).strip() for item in cases if str(item.get("rx_id", "")).strip() in selected]
    else:
        rx_ids = [str(item.get("rx_id", "")).strip() for item in cases if str(item.get("rx_id", "")).strip()]
    if not rx_ids:
        raise RuntimeError("No rx ids resolved for migrate stage")
    migrate_confsearch_seeds(
        Path(args.session_dir).expanduser().resolve(),
        rx_ids,
        cases=cases,
        confsearch_root=args.confsearch_root,
        seed_protocol=args.seed_protocol,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
