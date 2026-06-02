from __future__ import annotations
# pyright: reportDeprecated=false, reportUnusedCallResult=false, reportAny=false, reportExplicitAny=false, reportUnknownMemberType=false, reportUnknownVariableType=false, reportUnknownArgumentType=false

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Union


PathLike = Union[str, Path]
JsonDict = Dict[str, Any]

BENCHMARK_TYPE = "dft_theory"
BENCHMARK_MANIFEST_SCHEMA_VERSION = "benchmark_v2"
RX_MANIFEST_SCHEMA_VERSION = "benchmark_rx_v2"
BASELINE_SCHEMA_VERSION = "benchmark_baseline_v2"
PHASE1_RESULT_SCHEMA_VERSION = "benchmark_phase1_sp_v2"
PHASE1_CANONICAL_POINTS = ("precursor", "intermediate", "ts", "product")
PHASE2_CANONICAL_POINTS = ("complex", "product", "ts")

DEFAULT_BENCHMARK_STAGES = {
    "migrate": "pending",
    "baseline": "pending",
    "prepare": "pending",
    "phase1_sp": "pending",
    "phase2_geo": "pending",
    "shermo_correction": "pending",
    "final_report": "pending",
}

DEFAULT_BENCHMARK_REPORTS = {
    "phase1_global_summary_json": None,
    "phase1_global_summary_md": None,
    "phase2_global_summary_json": None,
    "phase2_global_summary_md": None,
    "shermo_global_summary_json": None,
    "shermo_global_summary_md": None,
    "final_report_json": None,
    "final_report_md": None,
}

DEFAULT_RX_REPORTS = {
    "summary_json": None,
    "summary_md": None,
    "sp_summary_json": None,
    "sp_summary_md": None,
    "geo_summary_json": None,
    "geo_summary_md": None,
    "shermo_summary_json": None,
    "shermo_summary_md": None,
}

PHASE_ALIASES = {
    "sp": "phase1",
    "phase1": "phase1",
    "phase1_sp": "phase1",
    "geo": "phase2",
    "phase2": "phase2",
    "phase2_geo": "phase2",
    "shermo": "shermo",
    "shermo_correction": "shermo",
}

PHASE_MATRIX_KEYS = {
    "phase1": "phase1_methods",
    "phase2": "phase2_methods",
    "shermo": "shermo_methods",
}

PHASE_WINNER_KEYS = {
    "phase1": "current_phase1_winner",
    "phase2": "current_phase2_winner",
    "shermo": "current_shermo_winner",
}


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _session_path(session_dir: PathLike) -> Path:
    return Path(str(session_dir)).expanduser()


def _normalize_rx_id(rx_id: Any) -> str:
    value = str(rx_id).strip()
    if value.lower().startswith("rx") and len(value) > 2:
        candidate = value[2:].strip()
        if candidate:
            return candidate
    return value


def _normalize_phase(phase: str) -> str:
    key = str(phase).strip()
    if not key:
        raise ValueError("phase is required")
    return PHASE_ALIASES.get(key, key)


def resolve_benchmark_root(session_dir: PathLike) -> Path:
    session_path = _session_path(session_dir)
    if session_path.parent.name == "experiments":
        return session_path.parent.parent
    return session_path.parent


def resolve_rx_dir(session_dir: PathLike, rx_id: Any) -> Path:
    return _session_path(session_dir) / f"rx{_normalize_rx_id(rx_id)}"


def resolve_upstream_confsearch_dir(session_dir: PathLike, rx_id: Any) -> Path:
    return resolve_rx_dir(session_dir, rx_id) / "upstream" / "confsearch"


def resolve_baseline_cache_dir(session_dir: PathLike, rx_id: Any, baseline_id: str) -> Path:
    normalized_baseline_id = str(baseline_id).strip()
    if not normalized_baseline_id:
        raise ValueError("baseline_id is required")
    return resolve_benchmark_root(session_dir) / "baselines" / f"rx{_normalize_rx_id(rx_id)}" / normalized_baseline_id


def resolve_method_run_dir(session_dir: PathLike, rx_id: Any, phase: str, method_id: str) -> Path:
    phase_key = _normalize_phase(phase)
    phase_dir = "sp" if phase_key == "phase1" else "geo" if phase_key == "phase2" else phase_key
    return resolve_rx_dir(session_dir, rx_id) / phase_dir / str(method_id)


def resolve_rx_reports_dir(session_dir: PathLike, rx_id: Any) -> Path:
    return resolve_rx_dir(session_dir, rx_id) / "reports"


def resolve_rx_phase_report_paths(session_dir: PathLike, rx_id: Any, phase: str) -> Dict[str, Path]:
    phase_key = _normalize_phase(phase)
    reports_dir = resolve_rx_reports_dir(session_dir, rx_id)
    stem = "sp_summary" if phase_key == "phase1" else "geo_summary" if phase_key == "phase2" else f"{phase_key}_summary"
    return {
        "json": reports_dir / f"{stem}.json",
        "md": reports_dir / f"{stem}.md",
    }


def _phase_matrix_key(phase: str) -> str:
    return PHASE_MATRIX_KEYS.get(phase, f"{phase}_methods")


def _phase_winner_key(phase: str) -> str:
    return PHASE_WINNER_KEYS.get(phase, f"current_{phase}_winner")


def _serialize_json(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _serialize_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_serialize_json(item) for item in value]
    return value


def _ensure_mapping(container: MutableMapping[str, Any], key: str) -> JsonDict:
    value = container.get(key)
    if isinstance(value, dict):
        return value
    value = {}
    container[key] = value
    return value


def _ensure_list(container: MutableMapping[str, Any], key: str) -> List[Any]:
    value = container.get(key)
    if isinstance(value, list):
        return value
    value = []
    container[key] = value
    return value


def _read_json_file(path: Path) -> JsonDict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid manifest payload: {path}")
    return payload


def _write_json_file(path: Path, payload: Mapping[str, Any]) -> JsonDict:
    materialized = _serialize_json(payload)
    if not isinstance(materialized, dict):
        raise ValueError(f"Invalid manifest payload: {path}")
    materialized["updated_at"] = now_iso()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(materialized, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return materialized


def _require_schema_version(payload: Mapping[str, Any], expected: str, context: str) -> None:
    meta = payload.get("meta")
    actual = meta.get("schema_version") if isinstance(meta, Mapping) else None
    if actual != expected:
        raise ValueError(f"{context} requires schema_version={expected}, got {actual!r}")


def _ensure_rx_phase(payload: MutableMapping[str, Any], phase: str) -> JsonDict:
    phase_payload = _ensure_mapping(payload, phase)
    _ensure_mapping(phase_payload, "results")
    phase_payload.setdefault("winner", None)
    if phase == "phase2":
        phase_payload.setdefault("sp_reader", None)
    return phase_payload


def resolve_benchmark_manifest_path(session_dir: PathLike) -> Path:
    return _session_path(session_dir) / "manifests" / "benchmark_manifest.json"


def resolve_rx_manifest_path(session_dir: PathLike, rx_id: Any) -> Path:
    return _session_path(session_dir) / "manifests" / f"rx{_normalize_rx_id(rx_id)}.json"


def read_benchmark_manifest(session_dir: PathLike) -> JsonDict:
    payload = _read_json_file(resolve_benchmark_manifest_path(session_dir))
    _require_schema_version(payload, BENCHMARK_MANIFEST_SCHEMA_VERSION, "Benchmark session")
    return payload


def write_benchmark_manifest(session_dir: PathLike, payload: Mapping[str, Any]) -> JsonDict:
    return _write_json_file(resolve_benchmark_manifest_path(session_dir), payload)


def read_rx_manifest(session_dir: PathLike, rx_id: Any) -> JsonDict:
    payload = _read_json_file(resolve_rx_manifest_path(session_dir, rx_id))
    _require_schema_version(payload, RX_MANIFEST_SCHEMA_VERSION, f"Benchmark rx manifest rx{_normalize_rx_id(rx_id)}")
    _validate_rx_manifest_payload(payload, rx_id=rx_id)
    return payload


def write_rx_manifest(session_dir: PathLike, rx_id: Any, payload: Mapping[str, Any]) -> JsonDict:
    return _write_json_file(resolve_rx_manifest_path(session_dir, rx_id), payload)


def create_session_manifest(
    session_dir: PathLike,
    *,
    cases: Optional[Sequence[Any]] = None,
    methods: Optional[Mapping[str, Any]] = None,
    config_snapshot: Optional[Mapping[str, Any]] = None,
    reaction_type: Optional[str] = None,
    benchmark_type: str = BENCHMARK_TYPE,
    stages: Optional[Mapping[str, Any]] = None,
    reports: Optional[Mapping[str, Any]] = None,
) -> JsonDict:
    session_path = _session_path(session_dir)
    (session_path / "manifests").mkdir(parents=True, exist_ok=True)
    (session_path / "reports").mkdir(parents=True, exist_ok=True)

    manifest_path = resolve_benchmark_manifest_path(session_path)
    payload = _read_json_file(manifest_path) if manifest_path.exists() else {}

    meta = _ensure_mapping(payload, "meta")
    existing_schema_version = meta.get("schema_version")
    if existing_schema_version is not None and existing_schema_version != BENCHMARK_MANIFEST_SCHEMA_VERSION:
        raise ValueError(
            f"Benchmark session manifest schema mismatch: expected {BENCHMARK_MANIFEST_SCHEMA_VERSION}, got {existing_schema_version!r}"
        )
    meta.setdefault("schema_version", BENCHMARK_MANIFEST_SCHEMA_VERSION)
    meta.setdefault("benchmark_type", str(benchmark_type))
    meta.setdefault("created_at", now_iso())
    if reaction_type is not None:
        meta.setdefault("reaction_type", str(reaction_type))

    if isinstance(payload.get("cases"), list):
        existing_cases = payload["cases"]
    else:
        existing_cases = [_normalize_rx_id(item) for item in (cases or [])]
        payload["cases"] = existing_cases

    if cases is not None and not existing_cases:
        payload["cases"] = [_normalize_rx_id(item) for item in cases]

    methods_payload = _ensure_mapping(payload, "methods")
    method_defaults: JsonDict = {"sp": {}, "geo": {}}
    if methods:
        method_defaults.update(_serialize_json(methods))
    for key, value in method_defaults.items():
        methods_payload.setdefault(str(key), value)

    snapshot_payload = _ensure_mapping(payload, "config_snapshot")
    if config_snapshot:
        for key, value in _serialize_json(config_snapshot).items():
            snapshot_payload.setdefault(str(key), value)

    stages_payload = _ensure_mapping(payload, "stages")
    stage_defaults: JsonDict = dict(DEFAULT_BENCHMARK_STAGES)
    if stages:
        stage_defaults.update(_serialize_json(stages))
    for key, value in stage_defaults.items():
        stages_payload.setdefault(str(key), value)

    reports_payload = _ensure_mapping(payload, "reports")
    report_defaults: JsonDict = dict(DEFAULT_BENCHMARK_REPORTS)
    if reports:
        report_defaults.update(_serialize_json(reports))
    for key, value in report_defaults.items():
        reports_payload.setdefault(str(key), value)

    _ensure_list(payload, "warnings")
    return write_benchmark_manifest(session_path, payload)


def init_rx_manifest(
    session_dir: PathLike,
    rx_id: Any,
    *,
    inputs: Optional[Mapping[str, Any]] = None,
    reports: Optional[Mapping[str, Any]] = None,
) -> JsonDict:
    session_path = _session_path(session_dir)
    (session_path / "manifests").mkdir(parents=True, exist_ok=True)

    normalized_rx_id = _normalize_rx_id(rx_id)
    manifest_path = resolve_rx_manifest_path(session_path, normalized_rx_id)
    payload = _read_json_file(manifest_path) if manifest_path.exists() else {}

    meta = _ensure_mapping(payload, "meta")
    meta.setdefault("rx_id", normalized_rx_id)
    existing_schema_version = meta.get("schema_version")
    if existing_schema_version is not None and existing_schema_version != RX_MANIFEST_SCHEMA_VERSION:
        raise ValueError(
            f"Benchmark rx manifest schema mismatch for rx{normalized_rx_id}: expected {RX_MANIFEST_SCHEMA_VERSION}, got {existing_schema_version!r}"
        )
    meta.setdefault("schema_version", RX_MANIFEST_SCHEMA_VERSION)
    meta.setdefault("created_at", now_iso())

    inputs_payload = _ensure_mapping(payload, "inputs")
    if inputs:
        for key, value in _serialize_json(inputs).items():
            inputs_payload.setdefault(str(key), value)

    matrix_payload = _ensure_mapping(payload, "matrix")
    matrix_payload.setdefault("phase1_methods", [])
    matrix_payload.setdefault("phase2_methods", [])

    tasks_payload = _ensure_mapping(payload, "tasks")
    tasks_payload.setdefault("phase1", {})
    tasks_payload.setdefault("phase2", {})

    _ensure_rx_phase(payload, "phase1")
    _ensure_rx_phase(payload, "phase2")

    payload.setdefault("current_phase1_winner", None)
    payload.setdefault("current_phase2_winner", None)

    reports_payload = _ensure_mapping(payload, "reports")
    report_defaults: JsonDict = dict(DEFAULT_RX_REPORTS)
    if reports:
        report_defaults.update(_serialize_json(reports))
    for key, value in report_defaults.items():
        reports_payload.setdefault(str(key), value)

    _ensure_list(payload, "warnings")
    return write_rx_manifest(session_path, normalized_rx_id, payload)


def record_task_status(
    session_dir: PathLike,
    rx_id: Any,
    phase: str,
    method_id: str,
    *,
    status: Optional[str] = None,
    run_dir: Optional[PathLike] = None,
    exit_code: Optional[int] = None,
    runtime_seconds: Optional[float] = None,
    extra: Optional[Mapping[str, Any]] = None,
) -> JsonDict:
    payload = init_rx_manifest(session_dir, rx_id)
    phase_key = _normalize_phase(phase)

    tasks_payload = _ensure_mapping(payload, "tasks")
    phase_tasks = _ensure_mapping(tasks_payload, phase_key)
    task_payload = _ensure_mapping(phase_tasks, str(method_id))

    matrix_payload = _ensure_mapping(payload, "matrix")
    matrix_key = _phase_matrix_key(phase_key)
    methods = matrix_payload.get(matrix_key)
    if isinstance(methods, list):
        method_list = methods
    else:
        method_list = []
        matrix_payload[matrix_key] = method_list
    if method_id not in method_list:
        method_list.append(method_id)

    task_payload["status"] = status or ("done" if exit_code == 0 else "failed" if exit_code is not None else "pending")
    task_payload["run_dir"] = None if run_dir is None else str(_session_path(run_dir))
    task_payload["runtime_seconds"] = None if runtime_seconds is None else float(runtime_seconds)
    task_payload["exit_code"] = None if exit_code is None else int(exit_code)
    task_payload["updated_at"] = now_iso()

    if extra:
        task_payload.update(_serialize_json(extra))

    return write_rx_manifest(session_dir, rx_id, payload)


def set_stage_status(session_dir: PathLike, stage: str, status: str) -> JsonDict:
    payload = create_session_manifest(session_dir)
    stages_payload = _ensure_mapping(payload, "stages")
    stages_payload[str(stage)] = str(status)
    return write_benchmark_manifest(session_dir, payload)


def set_phase_winner(
    session_dir: PathLike,
    rx_id: Any,
    phase: str,
    winner: Optional[str],
    *,
    sp_reader: Optional[str] = None,
    extra: Optional[Mapping[str, Any]] = None,
) -> JsonDict:
    payload = init_rx_manifest(session_dir, rx_id)
    phase_key = _normalize_phase(phase)
    phase_payload = _ensure_rx_phase(payload, phase_key)

    phase_payload["winner"] = winner
    if sp_reader is not None:
        phase_payload["sp_reader"] = sp_reader
    if extra:
        phase_payload.update(_serialize_json(extra))

    payload[_phase_winner_key(phase_key)] = winner
    return write_rx_manifest(session_dir, rx_id, payload)


def append_warning(
    session_dir: PathLike,
    warning: Any,
    *,
    rx_id: Optional[Any] = None,
) -> JsonDict:
    if rx_id is None:
        payload = create_session_manifest(session_dir)
        warnings = _ensure_list(payload, "warnings")
        warnings.append(_serialize_json(warning))
        return write_benchmark_manifest(session_dir, payload)

    payload = init_rx_manifest(session_dir, rx_id)
    warnings = _ensure_list(payload, "warnings")
    warnings.append(_serialize_json(warning))
    return write_rx_manifest(session_dir, rx_id, payload)


def _validate_rx_manifest_payload(payload: Mapping[str, Any], *, rx_id: Any) -> None:
    context = f"Benchmark rx manifest rx{_normalize_rx_id(rx_id)}"

    baseline = payload.get("baseline")
    if isinstance(baseline, Mapping) and baseline:
        require_canonical_baseline_block(baseline, context=context)

    phase1 = payload.get("phase1")
    if not isinstance(phase1, Mapping):
        return
    results = phase1.get("results")
    if results is None:
        return
    if not isinstance(results, Mapping):
        raise ValueError(f"{context} phase1.results must be a mapping")
    for method_id, row in results.items():
        if not isinstance(row, Mapping):
            raise ValueError(f"{context} phase1.results[{method_id!r}] must be a mapping")
        require_phase1_canonical_results(row, context=f"{context} phase1.results[{method_id!r}]")


def require_phase1_canonical_results(results: Mapping[str, Any], *, context: str) -> None:
    schema_version = results.get("schema_version")
    if schema_version != PHASE1_RESULT_SCHEMA_VERSION:
        raise ValueError(f"{context} requires schema_version={PHASE1_RESULT_SCHEMA_VERSION}, got {schema_version!r}")

    point_results = results.get("point_results")
    if not isinstance(point_results, Mapping):
        raise ValueError(f"{context} missing canonical point_results payload")

    forbidden_point_names = {"substrate", "reactant"}
    present_forbidden = sorted(name for name in forbidden_point_names if name in point_results)
    if present_forbidden:
        raise ValueError(f"{context} contains legacy phase1 point names: {', '.join(present_forbidden)}")

    missing_points = [name for name in PHASE1_CANONICAL_POINTS if name not in point_results]
    if missing_points:
        raise ValueError(f"{context} missing canonical phase1 point_results entries: {', '.join(missing_points)}")

    legacy_metric_keys = (
        "e_substrate_hartree",
        "e_reactant_hartree",
        "dE_substrate_to_intermediate_kcal",
        "dE_activation_kcal",
        "dE_reaction_kcal",
    )
    present_legacy_metrics = [key for key in legacy_metric_keys if key in results]
    if present_legacy_metrics:
        raise ValueError(f"{context} contains legacy phase1 metric keys: {', '.join(present_legacy_metrics)}")

    forbidden_gibbs_keys = sorted(
        key for key in results if str(key).startswith("g_") or str(key).startswith("dG_")
    )
    if forbidden_gibbs_keys:
        raise ValueError(f"{context} contains forbidden phase1 Gibbs keys: {', '.join(forbidden_gibbs_keys)}")


def require_phase2_canonical_results(results: Mapping[str, Any], *, context: str) -> None:
    forbidden_keys = ("reactant", "g_reactant_hartree")
    present_forbidden = [key for key in forbidden_keys if key in results]
    if present_forbidden:
        raise ValueError(f"{context} contains legacy phase2 keys: {', '.join(present_forbidden)}")

    missing_points = [name for name in PHASE2_CANONICAL_POINTS if name not in results]
    if missing_points:
        raise ValueError(f"{context} missing canonical phase2 entries: {', '.join(missing_points)}")


def require_canonical_stationary_points(stationary_points: Mapping[str, Any], *, context: str) -> None:
    forbidden = (
        "reactant_complex.xyz",
        "reactant.xyz",
        "substrate.xyz",
        "reactant_min.xyz",
        "substrate_min.xyz",
        "precursor_global_min.xyz",
    )
    present_forbidden = [name for name in forbidden if str(stationary_points.get(name, "")).strip()]
    if present_forbidden:
        raise ValueError(f"{context} contains legacy stationary points: {', '.join(present_forbidden)}")

    required = (
        "precursor_min.xyz",
        "complex.xyz",
        "intermediate.xyz",
        "ts_guess.xyz",
        "ts_final.xyz",
        "product_min.xyz",
    )
    missing = [name for name in required if not str(stationary_points.get(name, "")).strip()]
    if missing:
        raise ValueError(f"{context} missing canonical stationary points: {', '.join(missing)}")


def require_canonical_baseline_block(baseline: Mapping[str, Any], *, context: str) -> None:
    schema_version = baseline.get("schema_version")
    if schema_version != BASELINE_SCHEMA_VERSION:
        raise ValueError(f"{context} requires baseline schema_version={BASELINE_SCHEMA_VERSION}, got {schema_version!r}")
    stationary_points = baseline.get("stationary_points")
    if not isinstance(stationary_points, Mapping):
        raise ValueError(f"{context} missing baseline stationary_points payload")
    require_canonical_stationary_points(stationary_points, context=f"{context} baseline block")


def require_canonical_baseline_manifest(baseline_manifest: Mapping[str, Any], *, context: str) -> None:
    schema_version = baseline_manifest.get("schema_version")
    if schema_version != BASELINE_SCHEMA_VERSION:
        raise ValueError(f"{context} requires baseline schema_version={BASELINE_SCHEMA_VERSION}, got {schema_version!r}")
    stationary_points = baseline_manifest.get("stationary_points")
    if not isinstance(stationary_points, Mapping):
        raise ValueError(f"{context} missing baseline stationary_points payload")
    require_canonical_stationary_points(stationary_points, context=f"{context} baseline_manifest")
