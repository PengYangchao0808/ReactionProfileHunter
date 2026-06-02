from __future__ import annotations
# pyright: reportDeprecated=false, reportUnusedCallResult=false, reportAny=false, reportExplicitAny=false, reportUnknownMemberType=false, reportUnknownVariableType=false, reportUnknownArgumentType=false, reportUnknownLambdaType=false, reportReturnType=false

import argparse
import csv
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from benchmark.dft_theory.lib.state_manager import (
    BASELINE_SCHEMA_VERSION,
    append_warning,
    create_session_manifest,
    init_rx_manifest,
    read_rx_manifest,
    record_task_status,
    resolve_baseline_cache_dir,
    require_canonical_baseline_manifest,
    require_canonical_stationary_points,
    set_stage_status,
    write_rx_manifest,
)
from rph_core.orchestrator import ReactionProfileHunter
from rph_core.steps.step3_opt.validator import TSValidationError, TSValidator
from rph_core.utils.gaussian_log_parser import GaussianLogParser
from rph_core.utils.layout_contract import resolve_step_dir, resolve_required_files

DEFAULTS_CONFIG = PROJECT_ROOT / "config" / "defaults.yaml"
DEFAULT_CASES_CONFIG = PROJECT_ROOT / "benchmark" / "dft_theory" / "config" / "benchmark_cases.yaml"
DEFAULT_DATASET = PROJECT_ROOT / "data" / "reaxys_cleaned.csv"
BASELINE_METHOD_ID = "BASELINE"


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _load_cases(cases_config: Path) -> List[Dict[str, Any]]:
    data = yaml.safe_load(cases_config.read_text(encoding="utf-8")) or {}
    raw_cases = data.get("cases", [])
    return [dict(item) for item in raw_cases if isinstance(item, dict)]


def _resolve_rx_ids(cases: Sequence[Mapping[str, Any]], requested: Optional[Sequence[str]]) -> List[str]:
    configured = [str(item.get("rx_id", "")).strip() for item in cases]
    configured = [item for item in configured if item]
    if not requested:
        return configured
    selected = {str(item).strip() for item in requested if str(item).strip()}
    return [item for item in configured if item in selected]


def _lookup_product_smiles(dataset_path: Path, rx_id: str, cases_by_rx: Mapping[str, Mapping[str, Any]]) -> str:
    case = cases_by_rx.get(rx_id, {})
    case_smiles = str(case.get("smiles", "")).strip()
    if case_smiles:
        return case_smiles
    with dataset_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if str(row.get("rx_id", "")).strip() != rx_id:
                continue
            product_smiles = str(row.get("product_smiles_main", "")).strip()
            if product_smiles:
                return product_smiles
            break
    raise ValueError(f"Unable to resolve product_smiles_main for rx{rx_id}")


def _copy_pre_s2_snapshot(upstream_dir: Path, baseline_root: Path) -> Path:
    required_dirs = {
        "S0_Mechanism": upstream_dir / "S0_Mechanism",
        "S1_ConfGeneration": upstream_dir / "S1_ConfGeneration",
    }
    missing_dirs = [name for name, path in required_dirs.items() if not path.is_dir()]
    if missing_dirs:
        missing_display = ", ".join(missing_dirs)
        raise FileNotFoundError(
            f"Incomplete pre-S2 snapshot in {upstream_dir}: missing required upstream directories: {missing_display}"
        )

    for step_name, source_dir in required_dirs.items():
        destination = baseline_root / step_name
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(source_dir, destination)

    return baseline_root


def _resolve_required_input_path(rx_manifest: Mapping[str, Any], input_name: str, *, rx_id: str) -> Path:
    inputs = rx_manifest.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError(f"Benchmark rx manifest rx{rx_id} missing inputs payload")
    entry = inputs.get(input_name)
    if not isinstance(entry, dict):
        raise ValueError(f"Benchmark rx manifest rx{rx_id} missing required input '{input_name}'")
    raw_path = str(entry.get("path", "")).strip()
    if not raw_path:
        raise ValueError(f"Benchmark rx manifest rx{rx_id} missing path for input '{input_name}'")
    resolved = Path(raw_path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"Benchmark rx manifest rx{rx_id} input '{input_name}' not found: {resolved}")
    return resolved


def _resolve_s3_intermediate_opt_xyz(result: Any) -> Optional[Path]:
    """优先使用 S3 DFT 优化的 intermediate 几何，而非 S2 xTB 逆向扫描的中间体。

    S2 retro-scan 产出的是 xTB 级别的 intermediate.xyz，
    该几何在 CCSD(T)/DFT 级别下不一定是真正的势能面极小点。
    S3 intermediate_opt 才是 B3LYP/def2-SVP 级别的优化几何。
    """
    if hasattr(result, 'work_dir') and result.work_dir:
        s3_int = Path(result.work_dir) / "S3_TS" / "S3_intermediate_opt" / "standard" / "intermediate_opt.xyz"
        if s3_int.is_file():
            return s3_int
    return None


def _copy_stationary_points(baseline_root: Path, precursor_xyz: Path, result: Any) -> Dict[str, str]:
    stationary_dir = baseline_root / "stationary_points"
    stationary_dir.mkdir(parents=True, exist_ok=True)
    copied: Dict[str, str] = {}

    intermediate_source = _resolve_s3_intermediate_opt_xyz(result)
    if intermediate_source is None:
        intermediate_source = result.intermediate_xyz
        import logging
        logging.getLogger(__name__).warning(
            "rx%s: S3-DFT intermediate (intermediate_opt.xyz) not found, "
            "falling back to S2-xTB retro-scan intermediate.xyz. "
            "This may cause non-physical negative barriers at higher theory levels.",
            getattr(result, 'rx_id', '?')
        )

    candidates = {
        "precursor_min.xyz": precursor_xyz,
        "complex.xyz": result.intermediate_xyz,
        "intermediate.xyz": intermediate_source,
        "ts_guess.xyz": result.ts_guess_xyz,
        "ts_final.xyz": result.ts_final_xyz,
        "product_min.xyz": result.product_xyz,
    }
    missing: List[str] = []
    for target_name, source_path in candidates.items():
        if source_path is None:
            missing.append(target_name)
            continue
        source = Path(source_path).expanduser().resolve()
        if not source.is_file():
            missing.append(target_name)
            continue
        destination = stationary_dir / target_name
        shutil.copy2(str(source), str(destination))
        copied[target_name] = str(destination.resolve())
    if missing:
        raise FileNotFoundError(
            "Baseline canonical stationary point set is incomplete: {0}".format(
                ", ".join(sorted(missing))
            )
        )
    require_canonical_stationary_points(copied, context="Baseline stage stationary_points")
    return copied


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return f"sha256:{digest.hexdigest()}"


def _serialize_sp_report(result: Any, s3_dir: Path) -> Dict[str, Any]:
    sp_report = result.sp_matrix_report
    if sp_report is not None and hasattr(sp_report, "to_dict"):
        report_payload = dict(sp_report.to_dict())
    else:
        report_payload = {}
    metadata_path = s3_dir / "sp_matrix_metadata.json"
    report_payload["sp_matrix_metadata_json"] = str(metadata_path.resolve()) if metadata_path.exists() else None
    report_payload["product_thermo"] = str(result.product_thermo.resolve()) if result.product_thermo else None
    report_payload["ts_log"] = str(result.ts_log.resolve()) if result.ts_log else None
    report_payload["intermediate_log"] = str(result.intermediate_log.resolve()) if result.intermediate_log else None
    return report_payload


def _validate_ts(result: Any) -> Dict[str, Any]:
    validation: Dict[str, Any] = {
        "status": "failed",
        "n_imag": None,
        "imag_freq": None,
        "mode_valid": None,
        "bond_lengths": [],
        "error": None,
    }
    if result.ts_log is None or result.ts_final_xyz is None:
        validation["error"] = "missing_ts_artifacts"
        return validation

    parsed = GaussianLogParser.parse_log(Path(result.ts_log))
    if parsed is None or parsed.frequencies is None:
        validation["error"] = "missing_ts_frequencies"
        return validation

    validator = TSValidator()
    frequencies = parsed.frequencies
    imag = [float(value) for value in frequencies if float(value) < 0.0]
    validation["n_imag"] = len(imag)
    validation["imag_freq"] = imag[0] if imag else None

    try:
        validator.validate_imaginary_frequencies(frequencies)
        mode_valid = True
        bond_lengths: List[float] = []
        coordinates = parsed.coordinates
        if coordinates is not None and getattr(coordinates, "size", 0) and result.forming_bonds:
            validator.validate_bond_lengths(coordinates, list(result.forming_bonds), expected_range=(2.0, 2.4))
            for i, j in result.forming_bonds:
                distance = float(((coordinates[i] - coordinates[j]) ** 2).sum() ** 0.5)
                bond_lengths.append(distance)
        validation["status"] = "pass"
        validation["mode_valid"] = mode_valid
        validation["bond_lengths"] = bond_lengths
        return validation
    except (TSValidationError, ValueError) as exc:
        validation["mode_valid"] = False
        validation["error"] = str(exc)
        return validation


def _baseline_signature_payload(
    rx_manifest: Mapping[str, Any],
    *,
    rx_id: str,
    defaults_config: Path,
    reaction_profile: str,
    product_smiles: str,
) -> Dict[str, Any]:
    inputs = rx_manifest.get("inputs", {})
    if not isinstance(inputs, Mapping):
        inputs = {}

    def _input_fingerprint(name: str) -> Dict[str, Optional[str]]:
        entry = inputs.get(name)
        if not isinstance(entry, Mapping):
            return {"hash": None}
        stored_hash = str(entry.get("hash", "")).strip() or None
        raw_path = str(entry.get("path", "")).strip()
        if raw_path:
            candidate = Path(raw_path).expanduser().resolve()
            if candidate.exists():
                if stored_hash is None and candidate.is_file():
                    stored_hash = _file_sha256(candidate)
        return {"hash": stored_hash}

    return {
        "schema_version": BASELINE_SCHEMA_VERSION,
        "rx_id": rx_id,
        "reaction_profile": reaction_profile,
        "product_smiles": product_smiles,
        "defaults_config_sha256": _file_sha256(defaults_config),
        "seed_protocol": rx_manifest.get("upstream", {}).get("confsearch", {}).get("seed_protocol") if isinstance(rx_manifest.get("upstream"), Mapping) else None,
        "inputs": {
            "confsearch_provenance": _input_fingerprint("confsearch_provenance"),
            "mechanism_summary": _input_fingerprint("mechanism_summary"),
            "mechanism_graph": _input_fingerprint("mechanism_graph"),
            "product_xyz": _input_fingerprint("product_xyz"),
            "precursor_xyz": _input_fingerprint("precursor_xyz"),
        },
    }


def _baseline_id_from_payload(payload: Mapping[str, Any]) -> str:
    normalized = json.dumps(payload, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    digest = hashlib.sha256(normalized.encode("utf-8")).hexdigest()[:16]
    return f"bl_{digest}"


def _prepare_baseline_root(session_dir: Path, rx_id: str, baseline_id: str, force_clean: bool) -> Path:
    baseline_root = resolve_baseline_cache_dir(session_dir, rx_id, baseline_id)
    if force_clean and baseline_root.exists():
        shutil.rmtree(baseline_root)
    baseline_root.parent.mkdir(parents=True, exist_ok=True)
    return baseline_root


def _load_cached_baseline(baseline_root: Path) -> Optional[Dict[str, Any]]:
    baseline_manifest_path = baseline_root / "baseline_manifest.json"
    if not baseline_manifest_path.is_file():
        return None
    payload = json.loads(baseline_manifest_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return None
    require_canonical_baseline_manifest(payload, context=f"Cached baseline {baseline_root}")
    stationary_points = payload.get("stationary_points", {})
    if not isinstance(stationary_points, Mapping):
        return None
    require_canonical_stationary_points(stationary_points, context=f"Cached baseline {baseline_root}")
    return payload


def _build_baseline_manifest(
    rx_id: str,
    baseline_id: str,
    baseline_root: Path,
    stationary_points: Mapping[str, str],
    ts_validation: Mapping[str, Any],
    thermo_reference_path: Path,
    result: Any,
    signature_payload: Mapping[str, Any],
) -> Dict[str, Any]:
    s2_dir = resolve_step_dir(baseline_root, "s2")
    s3_dir = resolve_step_dir(baseline_root, "s3")
    s4_dir = resolve_step_dir(baseline_root, "s4")
    required_s2 = resolve_required_files(baseline_root, "s2")
    required_s3 = resolve_required_files(baseline_root, "s3")
    required_s4 = resolve_required_files(baseline_root, "s4")
    return {
        "schema_version": BASELINE_SCHEMA_VERSION,
        "rx_id": rx_id,
        "baseline_id": baseline_id,
        "baseline_method_id": BASELINE_METHOD_ID,
        "baseline_root": str(baseline_root.resolve()),
        "signature": dict(signature_payload),
        "stationary_points": dict(stationary_points),
        "step_dirs": {
            "s2": str(s2_dir.resolve()),
            "s3": str(s3_dir.resolve()),
            "s4": str(s4_dir.resolve()),
        },
        "required_outputs": {
            "s2": {key: str(value.resolve()) for key, value in required_s2.items()},
            "s3": {key: str(value.resolve()) for key, value in required_s3.items()},
            "s4": {key: str(value.resolve()) for key, value in required_s4.items()},
        },
        "ts_validation": dict(ts_validation),
        "thermo_reference_json": str(thermo_reference_path.resolve()),
        "product_thermo": str(result.product_thermo.resolve()) if result.product_thermo else None,
        "forming_bonds": [list(pair) for pair in (result.forming_bonds or tuple())],
    }


def _update_rx_baseline(session_dir: Path, rx_id: str, payload: Mapping[str, Any]) -> None:
    rx_manifest = read_rx_manifest(session_dir, rx_id)
    rx_manifest["baseline"] = dict(payload)
    write_rx_manifest(session_dir, rx_id, rx_manifest)


def run_baseline_stage(
    session_dir: Path,
    rx_ids: Sequence[str],
    *,
    cases_config: Path = DEFAULT_CASES_CONFIG,
    defaults_config: Path = DEFAULTS_CONFIG,
    dataset_path: Path = DEFAULT_DATASET,
    reaction_profile: str = "[4+3]_default",
    force_clean: bool = False,
) -> Dict[str, Dict[str, Any]]:
    cases = _load_cases(cases_config)
    cases_by_rx = {str(item.get("rx_id", "")).strip(): item for item in cases}
    create_session_manifest(
        session_dir,
        cases=[str(item.get("rx_id", "")).strip() for item in cases],
        config_snapshot={
            "defaults": str(defaults_config),
            "cases": str(cases_config),
            "dataset": str(dataset_path),
        },
        reaction_type=reaction_profile,
    )
    set_stage_status(session_dir, "baseline", "running")

    results: Dict[str, Dict[str, Any]] = {}
    for rx_id in rx_ids:
        normalized_rx = str(rx_id).strip()
        init_rx_manifest(session_dir, normalized_rx)
        baseline_root: Optional[Path] = None
        try:
            rx_manifest = read_rx_manifest(session_dir, normalized_rx)
            precursor_xyz = _resolve_required_input_path(rx_manifest, "precursor_xyz", rx_id=normalized_rx)
        except (FileNotFoundError, ValueError) as exc:
            record_task_status(
                session_dir,
                normalized_rx,
                "baseline",
                BASELINE_METHOD_ID,
                status="failed",
                run_dir=baseline_root,
                exit_code=1,
                extra={"error_message": str(exc)},
            )
            raise
        product_smiles = _lookup_product_smiles(dataset_path, normalized_rx, cases_by_rx)

        signature_payload = _baseline_signature_payload(
            rx_manifest,
            rx_id=normalized_rx,
            defaults_config=defaults_config,
            reaction_profile=reaction_profile,
            product_smiles=product_smiles,
        )
        baseline_id = _baseline_id_from_payload(signature_payload)
        baseline_root = _prepare_baseline_root(session_dir, normalized_rx, baseline_id, force_clean=force_clean)
        cached_baseline = None
        try:
            cached_baseline = _load_cached_baseline(baseline_root)
        except (FileNotFoundError, ValueError, json.JSONDecodeError):
            if baseline_root.exists():
                shutil.rmtree(baseline_root)
        if cached_baseline is not None:
            baseline_manifest_path = baseline_root / "baseline_manifest.json"
            thermo_reference_path = baseline_root / "thermo_reference.json"
            stationary_points = cached_baseline.get("stationary_points", {})
            ts_validation = cached_baseline.get("ts_validation", {})
            _update_rx_baseline(
                session_dir,
                normalized_rx,
                {
                    "schema_version": BASELINE_SCHEMA_VERSION,
                    "status": "done",
                    "baseline_id": baseline_id,
                    "method_id": BASELINE_METHOD_ID,
                    "root": str(baseline_root.resolve()),
                    "baseline_manifest": str(baseline_manifest_path.resolve()),
                    "thermo_reference_json": str(thermo_reference_path.resolve()),
                    "stationary_points": dict(stationary_points) if isinstance(stationary_points, Mapping) else {},
                    "ts_validation": dict(ts_validation) if isinstance(ts_validation, Mapping) else {},
                },
            )
            record_task_status(
                session_dir,
                normalized_rx,
                "baseline",
                BASELINE_METHOD_ID,
                status="done",
                run_dir=baseline_root,
                exit_code=0,
                extra={
                    "baseline_id": baseline_id,
                    "baseline_manifest": str(baseline_manifest_path.resolve()),
                    "thermo_reference_json": str(thermo_reference_path.resolve()),
                    "cache_hit": True,
                },
            )
            results[normalized_rx] = {
                "baseline_id": baseline_id,
                "baseline_manifest": str(baseline_manifest_path.resolve()),
                "root": str(baseline_root.resolve()),
                "stationary_points": dict(stationary_points) if isinstance(stationary_points, Mapping) else {},
                "thermo_reference_json": str(thermo_reference_path.resolve()),
                "ts_validation": dict(ts_validation) if isinstance(ts_validation, Mapping) else {},
            }
            continue

        upstream_snapshot_dir = session_dir / f"rx{normalized_rx}" / "upstream" / "confsearch"
        staging_root = Path(tempfile.mkdtemp(prefix=f"{baseline_id}_", dir=str(baseline_root.parent)))
        try:
            _copy_pre_s2_snapshot(upstream_snapshot_dir, staging_root)
        except FileNotFoundError as exc:
            shutil.rmtree(staging_root, ignore_errors=True)
            record_task_status(
                session_dir,
                normalized_rx,
                "baseline",
                BASELINE_METHOD_ID,
                status="failed",
                run_dir=baseline_root,
                exit_code=1,
                extra={"error_message": str(exc)},
            )
            append_warning(session_dir, {"code": "missing_upstream_seed", "rx_id": normalized_rx}, rx_id=normalized_rx)
            raise

        hunter = ReactionProfileHunter(config_path=defaults_config)
        result = hunter.run_pipeline(product_smiles=product_smiles, work_dir=staging_root, reaction_profile=reaction_profile)
        if not result.success:
            shutil.rmtree(staging_root, ignore_errors=True)
            record_task_status(
                session_dir,
                normalized_rx,
                "baseline",
                BASELINE_METHOD_ID,
                status="failed",
                run_dir=baseline_root,
                exit_code=1,
                extra={
                    "error_step": result.error_step,
                    "error_message": result.error_message,
                },
            )
            raise RuntimeError(f"Baseline pipeline failed for rx{normalized_rx}: {result.error_step} {result.error_message}")

        staging_stationary_points = _copy_stationary_points(staging_root, precursor_xyz, result)
        ts_validation = _validate_ts(result)

        s3_dir = resolve_step_dir(staging_root, "s3")
        thermo_reference_payload = _serialize_sp_report(result, s3_dir)
        thermo_reference_path = staging_root / "thermo_reference.json"
        _write_json(thermo_reference_path, thermo_reference_payload)

        if baseline_root.exists():
            shutil.rmtree(baseline_root)
        staging_root.replace(baseline_root)
        stationary_points = {
            name: str((baseline_root / "stationary_points" / name).resolve())
            for name in staging_stationary_points.keys()
        }
        require_canonical_stationary_points(stationary_points, context=f"Published baseline {baseline_id}")
        thermo_reference_path = baseline_root / "thermo_reference.json"

        baseline_manifest_payload = _build_baseline_manifest(
            normalized_rx,
            baseline_id,
            baseline_root,
            stationary_points,
            ts_validation,
            thermo_reference_path,
            result,
            signature_payload,
        )
        baseline_manifest_path = baseline_root / "baseline_manifest.json"
        _write_json(baseline_manifest_path, baseline_manifest_payload)

        _update_rx_baseline(
            session_dir,
            normalized_rx,
            {
                "schema_version": BASELINE_SCHEMA_VERSION,
                "status": "done",
                "baseline_id": baseline_id,
                "method_id": BASELINE_METHOD_ID,
                "root": str(baseline_root.resolve()),
                "baseline_manifest": str(baseline_manifest_path.resolve()),
                "thermo_reference_json": str(thermo_reference_path.resolve()),
                "stationary_points": dict(stationary_points),
                "ts_validation": dict(ts_validation),
            },
        )
        record_task_status(
            session_dir,
            normalized_rx,
            "baseline",
            BASELINE_METHOD_ID,
            status="done",
            run_dir=baseline_root,
            exit_code=0,
            extra={
                "baseline_id": baseline_id,
                "baseline_manifest": str(baseline_manifest_path.resolve()),
                "thermo_reference_json": str(thermo_reference_path.resolve()),
                "ts_validation_status": ts_validation.get("status"),
                "cache_hit": False,
            },
        )
        results[normalized_rx] = {
            "baseline_id": baseline_id,
            "baseline_manifest": str(baseline_manifest_path.resolve()),
            "root": str(baseline_root.resolve()),
            "stationary_points": dict(stationary_points),
            "thermo_reference_json": str(thermo_reference_path.resolve()),
            "ts_validation": dict(ts_validation),
        }

    set_stage_status(session_dir, "baseline", "done")
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run baseline DFT benchmark stage.")
    parser.add_argument("--session-dir", required=True, help="Benchmark session directory")
    parser.add_argument("--rx-id", nargs="*", default=None, help="Optional subset of rx ids")
    parser.add_argument("--cases-config", default=str(DEFAULT_CASES_CONFIG), help="Benchmark cases YAML")
    parser.add_argument("--defaults-config", default=str(DEFAULTS_CONFIG), help="Baseline defaults YAML")
    parser.add_argument("--dataset", default=str(DEFAULT_DATASET), help="Dataset CSV for rx_id -> smiles")
    parser.add_argument("--reaction-profile", default="[4+3]_default", help="Reaction profile key")
    parser.add_argument("--force-clean", action="store_true", help="Remove existing baseline directory before running")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cases = _load_cases(Path(args.cases_config).expanduser().resolve())
    rx_ids = _resolve_rx_ids(cases, args.rx_id)
    if not rx_ids:
        raise RuntimeError("No rx ids resolved for baseline stage")
    run_baseline_stage(
        Path(args.session_dir).expanduser().resolve(),
        rx_ids,
        cases_config=Path(args.cases_config).expanduser().resolve(),
        defaults_config=Path(args.defaults_config).expanduser().resolve(),
        dataset_path=Path(args.dataset).expanduser().resolve(),
        reaction_profile=str(args.reaction_profile),
        force_clean=bool(args.force_clean),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
