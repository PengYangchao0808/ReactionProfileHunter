from __future__ import annotations
# pyright: reportDeprecated=false, reportUnusedCallResult=false, reportAny=false, reportExplicitAny=false, reportUnknownMemberType=false, reportUnknownVariableType=false, reportUnknownArgumentType=false, reportUnknownLambdaType=false, reportReturnType=false

import argparse
import json
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from benchmark.dft_theory.lib.state_manager import (
    BASELINE_SCHEMA_VERSION,
    PHASE1_RESULT_SCHEMA_VERSION,
    create_session_manifest,
    init_rx_manifest,
    read_benchmark_manifest,
    read_rx_manifest,
    record_task_status,
    resolve_method_run_dir,
    resolve_rx_phase_report_paths,
    require_canonical_baseline_manifest,
    require_canonical_stationary_points,
    require_phase1_canonical_results,
    set_phase_winner,
    set_stage_status,
    write_benchmark_manifest,
    write_rx_manifest,
)
from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.utils.optimization_config import prepare_qc_method_config
from rph_core.utils.qc_task_runner import QCTaskRunner

DEFAULTS_CONFIG = PROJECT_ROOT / "config" / "defaults.yaml"
METHODS_SP_CONFIG = PROJECT_ROOT / "benchmark" / "dft_theory" / "config" / "methods_sp.yaml"
CASES_CONFIG = PROJECT_ROOT / "benchmark" / "dft_theory" / "config" / "benchmark_cases.yaml"
PHASE1_POINT_ORDER = ("precursor", "intermediate", "ts", "product")


def _read_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid JSON payload: {path}")
    return payload


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _safe_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _load_yaml(path: Path) -> Dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid YAML payload: {path}")
    return payload


def _resolve_rx_ids(session_dir: Path, requested: Optional[Sequence[str]]) -> List[str]:
    manifest_path = session_dir / "manifests" / "benchmark_manifest.json"
    if manifest_path.exists():
        manifest = _read_json(manifest_path)
        configured = [str(item) for item in manifest.get("cases", []) if str(item).strip()]
    else:
        cases_data = _load_yaml(CASES_CONFIG)
        configured = [
            str(item.get("rx_id", "")).strip()
            for item in (cases_data.get("cases", []) or [])
            if isinstance(item, dict) and str(item.get("rx_id", "")).strip()
        ]
    if not requested:
        return configured
    selected = {str(item).strip() for item in requested if str(item).strip()}
    return [item for item in configured if item in selected]


def _resolve_methods(methods_config: Path, requested: Optional[Sequence[str]]) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Any]]:
    payload = _load_yaml(methods_config)
    methods = payload.get("methods", {}) or {}
    if not isinstance(methods, dict):
        raise ValueError(f"Invalid methods block: {methods_config}")
    if not requested:
        selected = {str(key): dict(value) for key, value in methods.items() if isinstance(value, dict)}
    else:
        requested_set = {str(item).strip() for item in requested if str(item).strip()}
        selected = {
            str(key): dict(value)
            for key, value in methods.items()
            if str(key) in requested_set and isinstance(value, dict)
        }
    solvent = payload.get("solvent", {}) or {}
    return selected, dict(solvent) if isinstance(solvent, dict) else {}


def _load_defaults(defaults_config: Path) -> Dict[str, Any]:
    payload = _load_yaml(defaults_config)
    return payload


def _baseline_paths(rx_manifest: Mapping[str, Any]) -> Tuple[Path, Dict[str, str], Dict[str, Any]]:
    baseline = rx_manifest.get("baseline")
    if not isinstance(baseline, dict):
        raise ValueError("Missing baseline block in rx manifest")
    baseline_schema_version = baseline.get("schema_version")
    if baseline_schema_version != BASELINE_SCHEMA_VERSION:
        raise ValueError(
            f"SP benchmark baseline block requires schema_version={BASELINE_SCHEMA_VERSION}, got {baseline_schema_version!r}"
        )
    baseline_root = Path(str(baseline.get("root", ""))).expanduser().resolve()
    stationary_points = baseline.get("stationary_points", {})
    if not isinstance(stationary_points, dict):
        raise ValueError("Missing stationary_points in baseline block")
    require_canonical_stationary_points(stationary_points, context="SP benchmark baseline block")
    baseline_manifest_path = Path(str(baseline.get("baseline_manifest", ""))).expanduser().resolve()
    if not baseline_manifest_path.is_file():
        raise FileNotFoundError(f"Missing baseline manifest: {baseline_manifest_path}")
    baseline_manifest = _read_json(baseline_manifest_path)
    require_canonical_baseline_manifest(baseline_manifest, context="SP benchmark")
    return baseline_root, {str(key): str(value) for key, value in stationary_points.items()}, baseline_manifest


def _path_str(path: Any) -> Optional[str]:
    if path is None:
        return None
    return str(Path(path).resolve())


def _empty_point_result(error_message: Optional[str]) -> Dict[str, Any]:
    return {
        "converged": False,
        "energy_hartree": None,
        "log_file": None,
        "fchk_file": None,
        "error_message": error_message,
    }


def _empty_phase1_point_results(error_message: Optional[str]) -> Dict[str, Dict[str, Any]]:
    return {point_name: _empty_point_result(error_message) for point_name in PHASE1_POINT_ORDER}


def _failed_phase1_result(method_id: str, error_message: str) -> Dict[str, Any]:
    return {
        "schema_version": PHASE1_RESULT_SCHEMA_VERSION,
        "method_id": method_id,
        "status": "failed",
        "runtime_seconds": None,
        "point_results": _empty_phase1_point_results(error_message),
        "error_message": error_message,
        "hard_failure": False,
        **_compute_result_metrics(_empty_phase1_point_results(error_message)),
    }


def _run_single_point(
    runner: QCTaskRunner,
    xyz_file: Path,
    output_dir: Path,
) -> Dict[str, Any]:
    result = runner.run_sp_only(xyz_file, output_dir, charge=0, spin=1)
    log_path = result.qm_output_file or result.log_file or result.output_file
    return {
        "converged": bool(result.converged),
        "energy_hartree": _safe_float(result.energy) if result.converged else None,
        "log_file": _path_str(log_path),
        "fchk_file": _path_str(result.fchk_file),
        "error_message": result.error_message,
    }


def _build_sp_config(
    defaults_config: Mapping[str, Any],
    method_spec: Mapping[str, Any],
    solvent_cfg: Mapping[str, Any],
) -> Dict[str, Any]:
    # Benchmark SP stage mandates ORCA engine for single-point calculations.
    # Gaussian lacks RI acceleration for large basis sets, making it unsuitable.
    spec = dict(method_spec)
    spec["engine"] = "orca"
    return prepare_qc_method_config(
        defaults_config,
        single_point=spec,
        solvent=solvent_cfg,
    )


def _metric_kcal(upper: Optional[float], lower: Optional[float]) -> Optional[float]:
    if upper is None or lower is None:
        return None
    return (upper - lower) * HARTREE_TO_KCAL


def _compute_result_metrics(
    point_results: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Optional[float]]:
    e_precursor = _safe_float(point_results.get("precursor", {}).get("energy_hartree"))
    e_intermediate = _safe_float(point_results.get("intermediate", {}).get("energy_hartree"))
    e_ts = _safe_float(point_results.get("ts", {}).get("energy_hartree"))
    e_prod = _safe_float(point_results.get("product", {}).get("energy_hartree"))
    return {
        "e_precursor_hartree": e_precursor,
        "e_intermediate_hartree": e_intermediate,
        "e_ts_hartree": e_ts,
        "e_product_hartree": e_prod,
        "dE_precursor_to_intermediate_kcal": _metric_kcal(e_intermediate, e_precursor),
        "dE_intermediate_to_ts_kcal": _metric_kcal(e_ts, e_intermediate),
        "dE_intermediate_to_product_kcal": _metric_kcal(e_prod, e_intermediate),
    }


def _aggregate_phase1_status(point_results: Mapping[str, Mapping[str, Any]]) -> str:
    """根据所有驻点结果计算总体状态。"""
    all_ok = True
    any_hard = False
    for pt in point_results.values():
        if not pt.get("converged") or _safe_float(pt.get("energy_hartree")) is None:
            all_ok = False
        msg = str(pt.get("error_message", "")).upper()
        if "UNRECOGNIZED" in msg or "DUPLICATED KEYWORD" in msg:
            any_hard = True
    if any_hard:
        return "failed"
    if all_ok:
        return "completed"
    return "partial_failed"


def _run_method_points_partial(
    session_dir: Path,
    rx_id: str,
    method_id: str,
    method_spec: Mapping[str, Any],
    defaults_payload: Mapping[str, Any],
    solvent_cfg: Mapping[str, Any],
    point_map: Mapping[str, Path],
    selected_points: Sequence[str],
    force_clean: bool,
) -> Dict[str, Any]:
    """仅重算指定驻点的 SP，与已有 sp_result.json 合并。"""
    run_dir = resolve_method_run_dir(session_dir, rx_id, "phase1", method_id)
    result_path = run_dir / "sp_result.json"

    if not result_path.exists():
        return _run_one_method(
            session_dir, rx_id, method_id, method_spec,
            defaults_payload, solvent_cfg, point_map, force_clean,
        )

    existing = _read_json(result_path)
    qc_dir = run_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    for point_name in selected_points:
        point_dir = qc_dir / point_name
        if force_clean and point_dir.exists():
            shutil.rmtree(point_dir)

    config = _build_sp_config(defaults_payload, method_spec, solvent_cfg)
    runner = QCTaskRunner(config)
    started = time.time()

    partial_runtimes: Dict[str, float] = {}
    for point_name in selected_points:
        point_dir = qc_dir / point_name
        point_dir.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        pt_result = _run_single_point(runner, point_map[point_name], point_dir)
        partial_runtimes[point_name] = time.time() - t0
        existing["point_results"][point_name] = pt_result

    old_runtime = _safe_float(existing.get("runtime_seconds")) or 0.0
    partial_total = sum(partial_runtimes.values())
    existing["runtime_seconds"] = old_runtime + partial_total
    existing.setdefault("partial_runtime_seconds", {})
    existing["partial_runtime_seconds"].update(partial_runtimes)

    existing.update(_compute_result_metrics(existing["point_results"]))
    existing["status"] = _aggregate_phase1_status(existing["point_results"])
    existing["partial_update"] = {
        "points": list(selected_points),
        "updated_at": datetime.now(timezone.utc).isoformat(),
    }

    _write_json(result_path, existing)
    record_task_status(
        session_dir, rx_id, "sp", method_id,
        status=existing["status"],
        run_dir=run_dir,
        exit_code=0 if existing["status"] == "completed" else 1,
        runtime_seconds=existing["runtime_seconds"],
        extra={"result_json": str(result_path.resolve())},
    )
    return existing


def _baseline_point_map(stationary_points: Mapping[str, str]) -> Dict[str, Path]:
    require_canonical_stationary_points(stationary_points, context="SP benchmark stationary_points")
    precursor = stationary_points.get("precursor_min.xyz")
    intermediate = stationary_points.get("intermediate.xyz")
    ts = stationary_points.get("ts_final.xyz")
    product = stationary_points.get("product_min.xyz")
    if not precursor or not intermediate or not ts or not product:
        raise ValueError("Baseline stationary point set is incomplete")
    return {
        "precursor": Path(precursor).expanduser().resolve(),
        "intermediate": Path(intermediate).expanduser().resolve(),
        "ts": Path(ts).expanduser().resolve(),
        "product": Path(product).expanduser().resolve(),
    }


def _run_one_method(
    session_dir: Path,
    rx_id: str,
    method_id: str,
    method_spec: Mapping[str, Any],
    defaults_payload: Mapping[str, Any],
    solvent_cfg: Mapping[str, Any],
    point_map: Mapping[str, Path],
    force_clean: bool,
) -> Dict[str, Any]:
    run_dir = resolve_method_run_dir(session_dir, rx_id, "phase1", method_id)
    qc_dir = run_dir / "qc"
    if force_clean and run_dir.exists():
        shutil.rmtree(run_dir)
    qc_dir.mkdir(parents=True, exist_ok=True)

    config = _build_sp_config(defaults_payload, method_spec, solvent_cfg)
    runner = QCTaskRunner(config)
    started = time.time()
    point_results = _empty_phase1_point_results("not_run")
    status = "completed"
    error_message: Optional[str] = None
    hard_failure = False

    normalized_spec: Optional[Dict[str, Any]] = None
    rendered_simple_keywords: Optional[str] = None
    rendered_blocks: Dict[str, str] = {}
    sp_engine = runner.sp_engine
    if hasattr(sp_engine, "method_spec") and hasattr(sp_engine.method_spec, "to_dict"):
        normalized_spec = sp_engine.method_spec.to_dict()
    if hasattr(sp_engine, "render_simple_keywords"):
        rendered_simple_keywords = sp_engine.render_simple_keywords(task_type="sp")
    if hasattr(sp_engine, "render_blocks"):
        maybe_blocks = sp_engine.render_blocks()
        if isinstance(maybe_blocks, dict):
            rendered_blocks = {str(k): str(v) for k, v in maybe_blocks.items()}

    def _is_hard_keyword_error(message: Optional[str]) -> bool:
        text = str(message or "").upper()
        return "UNRECOGNIZED" in text or "DUPLICATED KEYWORD" in text

    for point_name, xyz_file in point_map.items():
        point_output_dir = qc_dir / point_name
        point_output_dir.mkdir(parents=True, exist_ok=True)
        result = _run_single_point(runner, xyz_file, point_output_dir)
        point_results[point_name] = result
        if not result.get("converged") or _safe_float(result.get("energy_hartree")) is None:
            error_message = result.get("error_message") or f"{point_name}_single_point_failed"
            if _is_hard_keyword_error(error_message):
                status = "failed"
                hard_failure = True
                break
            if status == "completed":
                status = "partial_failed"

    runtime_seconds = time.time() - started
    metrics = _compute_result_metrics(point_results)
    payload: Dict[str, Any] = {
        "schema_version": PHASE1_RESULT_SCHEMA_VERSION,
        "method_id": method_id,
        "engine": method_spec.get("engine", "orca"),
        "method": method_spec.get("method"),
        "basis": method_spec.get("basis"),
        "aux_basis": method_spec.get("aux_basis"),
        "route_extras": method_spec.get("route_extras"),
        "normalized_spec": normalized_spec,
        "rendered_simple_keywords": rendered_simple_keywords,
        "rendered_blocks": rendered_blocks,
        "is_reference": bool(method_spec.get("is_reference", False)),
        "status": status,
        "runtime_seconds": runtime_seconds,
        "point_results": point_results,
        "error_message": error_message,
        "hard_failure": hard_failure,
    }
    payload.update(metrics)
    result_path = run_dir / "sp_result.json"
    _write_json(result_path, payload)
    record_task_status(
        session_dir,
        rx_id,
        "sp",
        method_id,
        status=status,
        run_dir=run_dir,
        exit_code=0 if status == "completed" else 1,
        runtime_seconds=runtime_seconds,
        extra={"result_json": str(result_path.resolve())},
    )
    return payload


def run_sp_benchmark_stage(
    session_dir: Path,
    rx_ids: Sequence[str],
    *,
    methods_config: Path = METHODS_SP_CONFIG,
    defaults_config: Path = DEFAULTS_CONFIG,
    selected_methods: Optional[Sequence[str]] = None,
    selected_points: Optional[Sequence[str]] = None,
    force_clean: bool = False,
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    defaults_payload = _load_defaults(defaults_config)
    methods, solvent_cfg = _resolve_methods(methods_config, selected_methods)
    if not methods:
        raise RuntimeError("No SP benchmark methods resolved")
    create_session_manifest(
        session_dir,
        methods={"sp": methods},
        config_snapshot={"defaults": str(defaults_config), "methods_sp": str(methods_config)},
    )
    set_stage_status(session_dir, "phase1_sp", "running")

    outputs: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for rx_id in rx_ids:
        init_rx_manifest(session_dir, rx_id)
        rx_manifest = read_rx_manifest(session_dir, rx_id)
        _, stationary_points, _ = _baseline_paths(rx_manifest)
        point_map = _baseline_point_map(stationary_points)
        outputs[rx_id] = {}
        for method_id, method_spec in methods.items():
            if selected_points:
                outputs[rx_id][method_id] = _run_method_points_partial(
                    session_dir, rx_id, method_id, method_spec,
                    defaults_payload, solvent_cfg, point_map,
                    selected_points, force_clean,
                )
            else:
                outputs[rx_id][method_id] = _run_one_method(
                    session_dir, rx_id, method_id, method_spec,
                    defaults_payload, solvent_cfg, point_map,
                    force_clean,
                )
            if bool(outputs[rx_id][method_id].get("hard_failure")):
                set_stage_status(session_dir, "phase1_sp", "failed")
                raise RuntimeError(
                    f"SP benchmark hard failure at rx{rx_id}/{method_id}: "
                    f"{outputs[rx_id][method_id].get('error_message', 'unknown_error')}"
                )
    set_stage_status(session_dir, "phase1_sp", "computed")
    return outputs


def _phase1_point_energy(row: Mapping[str, Any], point_name: str) -> Optional[float]:
    point_results = row.get("point_results")
    if isinstance(point_results, dict):
        point_payload = point_results.get(point_name)
        if isinstance(point_payload, dict):
            point_energy = _safe_float(point_payload.get("energy_hartree"))
            if point_energy is not None:
                return point_energy

    fallback_keys = {
        "precursor": ("e_precursor_hartree",),
        "intermediate": ("e_intermediate_hartree",),
        "ts": ("e_ts_hartree",),
        "product": ("e_product_hartree",),
    }
    for key in fallback_keys.get(point_name, ()):
        value = _safe_float(row.get(key))
        if value is not None:
            return value
    return None


def _completed_phase1_point_count(row: Mapping[str, Any]) -> int:
    return sum(1 for point_name in PHASE1_POINT_ORDER if _phase1_point_energy(row, point_name) is not None)


def _annotate_phase1_scores(
    results: Mapping[str, Mapping[str, Any]],
    reference_method: Optional[str],
) -> Dict[str, Dict[str, Any]]:
    annotated = {method_id: dict(row) for method_id, row in results.items()}
    if not reference_method:
        return annotated

    reference_row = annotated.get(reference_method)
    if not isinstance(reference_row, dict):
        return annotated

    ref_energies = {point_name: _phase1_point_energy(reference_row, point_name) for point_name in PHASE1_POINT_ORDER}
    if any(value is None for value in ref_energies.values()):
        return annotated

    if str(reference_row.get("status", "")).lower() not in {"ok", "completed"}:
        return annotated

    for row in annotated.values():
        if str(row.get("status", "")).lower() not in {"ok", "completed"}:
            continue
        point_errors: Dict[str, float] = {}
        for point_name in PHASE1_POINT_ORDER:
            point_energy = _phase1_point_energy(row, point_name)
            ref_energy = ref_energies[point_name]
            if point_energy is None or ref_energy is None:
                point_errors = {}
                break
            point_errors[point_name] = abs(point_energy - ref_energy) * HARTREE_TO_KCAL
        if len(point_errors) != len(PHASE1_POINT_ORDER):
            continue
        error_values = list(point_errors.values())
        row["absolute_sp_energy_errors_vs_ref_kcal"] = point_errors
        row["absolute_sp_energy_mae_vs_ref_kcal"] = sum(error_values) / len(error_values)
        row["absolute_sp_energy_max_error_vs_ref_kcal"] = max(error_values)
    return annotated


def _choose_phase1_winner(results: Mapping[str, Mapping[str, Any]]) -> Optional[str]:
    successful = {
        k: v
        for k, v in results.items()
        if str(v.get("status", "")).lower() in {"ok", "completed"}
    }
    if not successful:
        return None

    scored: List[Tuple[str, float, float]] = []
    for method_id, row in successful.items():
        mae = _safe_float(row.get("absolute_sp_energy_mae_vs_ref_kcal"))
        if mae is None:
            continue
        max_error = _safe_float(row.get("absolute_sp_energy_max_error_vs_ref_kcal")) or 0.0
        scored.append((method_id, mae, max_error))
    if scored:
        scored.sort(key=lambda item: (item[1], item[2], item[0]))
        return scored[0][0]

    ranked = sorted(
        successful.items(),
        key=lambda item: (-_completed_phase1_point_count(item[1]), str(item[0])),
    )
    return ranked[0][0] if ranked else None


def _render_rx_markdown(
    rx_id: str,
    winner: Optional[str],
    reference_method: Optional[str],
    results: Mapping[str, Mapping[str, Any]],
) -> str:
    lines = [
        f"# SP Benchmark Summary — rx{rx_id}",
        "",
        f"- reference_method: {reference_method or 'NA'}",
        f"- winner: {winner or 'NA'}",
        "",
        "| method | status | abs-SP MAE vs ref (kcal) | E_int-E_precursor (kcal) | E_ts-E_int (kcal) | E_prod-E_int (kcal) | runtime(s) |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for method_id, row in results.items():
        lines.append(
            "| {method} | {status} | {abs_mae} | {de_is} | {de_ti} | {de_pi} | {runtime} |".format(
                method=method_id,
                status=row.get("status", "unknown"),
                abs_mae=_format_metric(row.get("absolute_sp_energy_mae_vs_ref_kcal")),
                de_is=_format_metric(row.get("dE_precursor_to_intermediate_kcal")),
                de_ti=_format_metric(row.get("dE_intermediate_to_ts_kcal")),
                de_pi=_format_metric(row.get("dE_intermediate_to_product_kcal")),
                runtime=_format_metric(row.get("runtime_seconds")),
            )
        )
    return "\n".join(lines) + "\n"


def _format_metric(value: Any) -> str:
    parsed = _safe_float(value)
    return "NA" if parsed is None else f"{parsed:.3f}"


def evaluate_sp_benchmark_stage(session_dir: Path, rx_ids: Sequence[str]) -> Dict[str, Dict[str, Any]]:
    benchmark_manifest = read_benchmark_manifest(session_dir)
    methods = benchmark_manifest.get("methods", {}).get("sp", {})
    reference_method = None
    if isinstance(methods, dict):
        for method_id, spec in methods.items():
            if isinstance(spec, dict) and bool(spec.get("is_reference")):
                reference_method = str(method_id)
                break

    phase1_global_votes: Dict[str, int] = {}
    summaries: Dict[str, Dict[str, Any]] = {}
    for rx_id in rx_ids:
        rx_manifest = read_rx_manifest(session_dir, rx_id)
        tasks = rx_manifest.get("tasks", {}).get("phase1", {})
        if not isinstance(tasks, dict):
            tasks = {}
        results: Dict[str, Dict[str, Any]] = {}
        for method_id, task_row in tasks.items():
            result_path = Path(str((task_row or {}).get("result_json", ""))).expanduser()
            if result_path.is_file():
                result_payload = _read_json(result_path)
                require_phase1_canonical_results(result_payload, context=f"SP benchmark result rx{rx_id}/{method_id}")
                results[str(method_id)] = result_payload
            else:
                results[str(method_id)] = _failed_phase1_result(str(method_id), "missing_sp_result_json")

        scored_results = _annotate_phase1_scores(results, reference_method)
        winner = _choose_phase1_winner(scored_results)
        set_phase_winner(session_dir, rx_id, "phase1", winner)
        updated_rx_manifest = read_rx_manifest(session_dir, rx_id)
        updated_rx_manifest.setdefault("phase1", {})
        updated_rx_manifest["phase1"]["reference_method"] = reference_method
        updated_rx_manifest["phase1"]["results"] = scored_results
        report_paths = resolve_rx_phase_report_paths(session_dir, rx_id, "phase1")
        summary_json = {
            "rx_id": rx_id,
            "phase": "sp",
            "reference_method": reference_method,
            "winner": winner,
            "results": scored_results,
        }
        _write_json(report_paths["json"], summary_json)
        _write_text(report_paths["md"], _render_rx_markdown(rx_id, winner, reference_method, scored_results))
        updated_rx_manifest.setdefault("reports", {})
        updated_rx_manifest["reports"]["sp_summary_json"] = str(report_paths["json"].resolve())
        updated_rx_manifest["reports"]["sp_summary_md"] = str(report_paths["md"].resolve())
        write_rx_manifest(session_dir, rx_id, updated_rx_manifest)
        summaries[rx_id] = summary_json
        if winner:
            phase1_global_votes[winner] = phase1_global_votes.get(winner, 0) + 1

    global_winner = max(phase1_global_votes.items(), key=lambda item: item[1])[0] if phase1_global_votes else None
    benchmark_manifest["phase1_winner"] = global_winner
    benchmark_manifest.setdefault("reports", {})
    benchmark_manifest["reports"]["phase1_global_summary_json"] = str((session_dir / "reports" / "phase1_global_summary.json").resolve())
    benchmark_manifest["reports"]["phase1_global_summary_md"] = str((session_dir / "reports" / "phase1_global_summary.md").resolve())
    write_benchmark_manifest(session_dir, benchmark_manifest)

    global_summary = {
        "phase": "sp",
        "winner": global_winner,
        "votes": phase1_global_votes,
        "rx": summaries,
    }
    _write_json(session_dir / "reports" / "phase1_global_summary.json", global_summary)
    lines = ["# Phase1 Global Summary", "", f"- winner: {global_winner or 'NA'}", "", "| method | votes |", "|---|---:|"]
    for method_id, votes in sorted(phase1_global_votes.items(), key=lambda item: item[1], reverse=True):
        lines.append(f"| {method_id} | {votes} |")
    _write_text(session_dir / "reports" / "phase1_global_summary.md", "\n".join(lines) + "\n")
    set_stage_status(session_dir, "phase1_sp", "done")
    return summaries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run or evaluate SP benchmark stage.")
    parser.add_argument("--session-dir", required=True, help="Benchmark session directory")
    parser.add_argument("--mode", choices=["run", "evaluate", "all"], default="all", help="Stage mode")
    parser.add_argument("--rx-id", nargs="*", default=None, help="Optional subset of rx ids")
    parser.add_argument("--methods", nargs="*", default=None, help="Optional subset of SP method ids")
    parser.add_argument("--points", nargs="*", default=None,
                        choices=list(PHASE1_POINT_ORDER),
                        help="仅重算指定驻点 (默认: 全部四个)")
    parser.add_argument("--methods-config", default=str(METHODS_SP_CONFIG), help="SP methods YAML")
    parser.add_argument("--defaults-config", default=str(DEFAULTS_CONFIG), help="Baseline defaults YAML")
    parser.add_argument("--force-clean", action="store_true", help="Remove existing per-method directories before running")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    session_dir = Path(args.session_dir).expanduser().resolve()
    rx_ids = _resolve_rx_ids(session_dir, args.rx_id)
    if not rx_ids:
        raise RuntimeError("No rx ids resolved for SP benchmark")
    if args.mode in {"run", "all"}:
        run_sp_benchmark_stage(
            session_dir,
            rx_ids,
            methods_config=Path(args.methods_config).expanduser().resolve(),
            defaults_config=Path(args.defaults_config).expanduser().resolve(),
            selected_methods=args.methods,
            selected_points=args.points,
            force_clean=bool(args.force_clean),
        )
    if args.mode in {"evaluate", "all"}:
        evaluate_sp_benchmark_stage(session_dir, rx_ids)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
