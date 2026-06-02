from __future__ import annotations
# pyright: reportDeprecated=false, reportUnusedCallResult=false, reportAny=false, reportExplicitAny=false, reportUnknownMemberType=false, reportUnknownVariableType=false, reportUnknownArgumentType=false, reportUnknownLambdaType=false, reportReturnType=false

import argparse
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from benchmark.dft_theory.lib.state_manager import (
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
    require_phase2_canonical_results,
    set_phase_winner,
    set_stage_status,
    write_benchmark_manifest,
    write_rx_manifest,
)
import logging

logger = logging.getLogger(__name__)
from rph_core.steps.step3_opt.validator import TSValidationError, TSValidator
from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.utils.file_io import read_xyz
from rph_core.utils.optimization_config import prepare_qc_method_config
from rph_core.utils.qc_task_runner import QCTaskRunner
from rph_core.utils.shermo_runner import run_shermo

DEFAULTS_CONFIG = PROJECT_ROOT / "config" / "defaults.yaml"
METHODS_GEO_CONFIG = PROJECT_ROOT / "benchmark" / "dft_theory" / "config" / "methods_geo.yaml"
CASES_CONFIG = PROJECT_ROOT / "benchmark" / "dft_theory" / "config" / "benchmark_cases.yaml"

DEFAULT_GEO_SP_READER = "SP-1"


def _read_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid JSON payload: {path}")
    return payload


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")


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


def _sha256_file(path: Path) -> Optional[str]:
    """Compute SHA-256 hex digest of a file. Returns None if file does not exist."""
    if not path.is_file():
        return None
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


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


def _resolve_geo_methods(methods_config: Path, requested: Optional[Sequence[str]]) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Any]]:
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
    return _load_yaml(defaults_config)


def _phase1_winner_method_id(rx_manifest: Mapping[str, Any]) -> Optional[str]:
    phase1 = rx_manifest.get("phase1")
    if isinstance(phase1, dict):
        if "winner_artifact" in phase1:
            raise ValueError("Legacy phase1 winner_artifact is not supported by benchmark v2")
        for key in ("winner_method_id", "winner"):
            value = phase1.get(key)
            if isinstance(value, str) and value.strip():
                return value.strip()
            if isinstance(value, Mapping):
                for nested_key in ("method_id", "winner_method_id", "id", "method", "winner"):
                    nested_value = value.get(nested_key)
                    if isinstance(nested_value, str) and nested_value.strip():
                        return nested_value.strip()
    winner = rx_manifest.get("current_phase1_winner")
    if isinstance(winner, str) and winner.strip():
        return winner.strip()
    return None


def _resolve_sp_reader_spec(
    benchmark_manifest: Mapping[str, Any],
    requested: Optional[str],
) -> Tuple[str, Dict[str, Any]]:
    """Resolve SP reader for GEO benchmark.

    Defaults to DEFAULT_GEO_SP_READER (SP-1).  DLPNO/CCSD(T) and reference
    methods are explicitly blocked for OPT benchmark.
    """
    sp_reader = requested or DEFAULT_GEO_SP_READER
    methods = benchmark_manifest.get("methods", {}).get("sp", {})
    if not isinstance(methods, dict):
        raise ValueError("Missing SP methods in benchmark manifest")
    spec = methods.get(sp_reader)
    if not isinstance(spec, dict):
        raise ValueError(f"Missing SP reader spec for {sp_reader}")

    method_name = str(spec.get("method", "")).upper()
    family = str(spec.get("family", "")).lower()
    if "DLPNO" in method_name or "CCSD" in method_name or family == "post_hf_dlpno":
        raise ValueError(
            f"DLPNO/CCSD(T) methods are disabled for OPT benchmark. "
            f"Requested: {sp_reader} ({spec.get('method', 'N/A')})"
        )

    if spec.get("is_reference"):
        raise ValueError(
            f"SP reference method {sp_reader} is not allowed as OPT SP reader"
        )

    return sp_reader, dict(spec)


def _resolve_winner_sp_spec(benchmark_manifest: Mapping[str, Any], rx_manifest: Mapping[str, Any]) -> Tuple[str, Dict[str, Any]]:
    winner = _phase1_winner_method_id(rx_manifest)
    if not winner:
        raise ValueError("Missing phase1 winner method id in rx manifest")
    methods = benchmark_manifest.get("methods", {}).get("sp", {})
    if not isinstance(methods, dict):
        raise ValueError("Missing SP methods in benchmark manifest")
    spec = methods.get(winner)
    if not isinstance(spec, dict):
        raise ValueError(f"Missing SP winner spec for {winner}")
    return winner, dict(spec)


def _phase2_sp_reader_method_id(
    rx_manifest: Mapping[str, Any],
    results: Optional[Mapping[str, Mapping[str, Any]]] = None,
) -> str:
    phase2 = rx_manifest.get("phase2")
    if isinstance(phase2, Mapping):
        for key in ("sp_reader", "sp_reader_method_id"):
            value = phase2.get(key)
            if isinstance(value, str) and value.strip():
                return value.strip()

    if results is not None:
        for result in results.values():
            if not isinstance(result, Mapping):
                continue
            for key in ("sp_reader_method_id", "sp_reader"):
                value = result.get(key)
                if isinstance(value, str) and value.strip():
                    return value.strip()

    return DEFAULT_GEO_SP_READER


def _baseline_context(rx_manifest: Mapping[str, Any]) -> Tuple[Path, Dict[str, str], Dict[str, Any]]:
    baseline = rx_manifest.get("baseline")
    if not isinstance(baseline, dict):
        raise ValueError("Missing baseline block in rx manifest")
    baseline_root = Path(str(baseline.get("root", ""))).expanduser().resolve()
    stationary_points = baseline.get("stationary_points", {})
    if not isinstance(stationary_points, dict):
        raise ValueError("Missing baseline stationary_points")
    require_canonical_stationary_points(stationary_points, context="GEO benchmark baseline block")
    baseline_manifest_path = Path(str(baseline.get("baseline_manifest", ""))).expanduser().resolve()
    if not baseline_manifest_path.is_file():
        raise FileNotFoundError(f"Missing baseline manifest: {baseline_manifest_path}")
    baseline_manifest = _read_json(baseline_manifest_path)
    require_canonical_baseline_manifest(baseline_manifest, context="GEO benchmark")
    return baseline_root, {str(k): str(v) for k, v in stationary_points.items()}, baseline_manifest


def _build_geo_config(
    defaults_payload: Mapping[str, Any],
    geo_method_spec: Mapping[str, Any],
    sp_method_spec: Mapping[str, Any],
    solvent_cfg: Mapping[str, Any],
) -> Dict[str, Any]:
    return prepare_qc_method_config(
        defaults_payload,
        optimization=geo_method_spec,
        single_point=sp_method_spec,
        solvent=solvent_cfg,
    )


def _point_inputs(stationary_points: Mapping[str, str]) -> Dict[str, Path]:
    require_canonical_stationary_points(stationary_points, context="GEO benchmark stationary_points")
    complex_xyz = stationary_points.get("complex.xyz")
    product = stationary_points.get("product_min.xyz")
    ts_final = stationary_points.get("ts_final.xyz")
    ts_guess = stationary_points.get("ts_guess.xyz")
    if not complex_xyz or not product or not ts_final:
        raise ValueError("Incomplete baseline stationary point set for GEO benchmark")
    return {
        "complex": Path(complex_xyz).expanduser().resolve(),
        "product": Path(product).expanduser().resolve(),
        "ts_final": Path(ts_final).expanduser().resolve(),
        "ts_guess": Path(ts_guess).expanduser().resolve() if ts_guess else Path(ts_final).expanduser().resolve(),
    }


def _thermo_value(result: Any, config: Mapping[str, Any], output_dir: Path) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    l2_energy = _safe_float(getattr(result, "l2_energy", None))
    freq_log = getattr(result, "freq_log", None) or getattr(result, "log_file", None) or getattr(result, "qm_output_file", None)
    if l2_energy is None or freq_log is None:
        return None, None, "missing_l2_or_freq_log"

    executables = config.get("executables", {}) if isinstance(config.get("executables"), dict) else {}
    shermo_value = executables.get("shermo")
    if not shermo_value:
        return None, None, "missing_shermo_executable"
    if isinstance(shermo_value, dict):
        shermo_value = shermo_value.get("path")
    if not shermo_value:
        return None, None, "missing_shermo_executable_path"
    shermo_bin = Path(str(shermo_value)).expanduser()
    if not shermo_bin.exists():
        return None, None, f"missing_shermo_binary:{shermo_bin}"

    try:
        thermo = run_shermo(shermo_bin, Path(freq_log), float(l2_energy), output_dir / "thermo.sum")
    except Exception as exc:
        return None, None, str(exc)
    value = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
    return float(value), str(getattr(thermo, "output_file", output_dir / "thermo.sum")), None


def _validate_ts(ts_xyz: Path, frequencies: Any, forming_bonds: Sequence[Sequence[int]]) -> Dict[str, Any]:
    validation: Dict[str, Any] = {
        "status": "failed",
        "n_imag": None,
        "imag_freq": None,
        "mode_valid": None,
        "bond_lengths": [],
        "error": None,
    }
    if frequencies is None:
        validation["error"] = "missing_frequencies"
        return validation
    validator = TSValidator()
    imag = [float(value) for value in frequencies if float(value) < 0.0]
    validation["n_imag"] = len(imag)
    validation["imag_freq"] = imag[0] if imag else None
    try:
        validator.validate_imaginary_frequencies(frequencies)
        coords, _ = read_xyz(ts_xyz)
        bonds = [(int(pair[0]), int(pair[1])) for pair in forming_bonds if len(pair) == 2]
        bond_lengths: List[float] = []
        if bonds:
            validator.validate_bond_lengths(coords, bonds, expected_range=(2.0, 2.4))
            for i, j in bonds:
                distance = float(((coords[i] - coords[j]) ** 2).sum() ** 0.5)
                bond_lengths.append(distance)
        validation["status"] = "pass"
        validation["mode_valid"] = True
        validation["bond_lengths"] = bond_lengths
        return validation
    except (TSValidationError, ValueError) as exc:
        validation["mode_valid"] = False
        validation["error"] = str(exc)
        return validation


def _geometry_displacement(
    baseline_xyz: Path,
    optimized_xyz: Optional[str],
) -> Dict[str, Any]:
    """Compare optimized geometry against baseline input via Kabsch-aligned RMSD."""
    from rph_core.utils.geometry_tools import GeometryUtils

    result: Dict[str, Any] = {
        "rmsd_to_baseline_angstrom": None,
        "aligned_rmsd_to_baseline_angstrom": None,
        "mean_atom_displacement_angstrom": None,
        "max_atom_displacement_angstrom": None,
        "error": None,
    }
    if optimized_xyz is None or not baseline_xyz.is_file():
        result["error"] = "missing_input"
        return result

    try:
        coords_baseline, _ = read_xyz(baseline_xyz)
        opt_path = Path(optimized_xyz)
        if not opt_path.is_file():
            result["error"] = "opt_xyz_missing"
            return result
        coords_opt, _ = read_xyz(opt_path)
        if coords_baseline.shape != coords_opt.shape:
            result["error"] = (
                f"atom_count_mismatch: {coords_baseline.shape[0]}"
                f" vs {coords_opt.shape[0]}"
            )
            return result

        cb = coords_baseline.astype(np.float64)
        co = coords_opt.astype(np.float64)
        cb_centered = cb - cb.mean(axis=0)
        co_centered = co - co.mean(axis=0)

        H = co_centered.T @ cb_centered
        U, _, Vt = np.linalg.svd(H)
        rot = U @ Vt
        if np.linalg.det(rot) < 0:
            Vt[-1, :] *= -1
            rot = U @ Vt
        co_aligned = co_centered @ rot

        deltas = np.linalg.norm(cb_centered - co_aligned, axis=1)
        result["rmsd_to_baseline_angstrom"] = float(
            GeometryUtils.calculate_rmsd(cb, co)
        )
        result["aligned_rmsd_to_baseline_angstrom"] = float(
            np.sqrt(np.mean(np.sum((cb_centered - co_aligned) ** 2, axis=1)))
        )
        result["mean_atom_displacement_angstrom"] = float(deltas.mean())
        result["max_atom_displacement_angstrom"] = float(deltas.max())
    except Exception as exc:
        result["error"] = str(exc)
    return result


def _result_payload(result: Any) -> Dict[str, Any]:
    optimized_xyz_path = Path(result.optimized_xyz).resolve() if getattr(result, "optimized_xyz", None) else None
    payload: Dict[str, Any] = {
        "converged": bool(getattr(result, "converged", False)),
        "optimized_xyz": str(optimized_xyz_path) if optimized_xyz_path else None,
        "optimized_xyz_sha256": _sha256_file(optimized_xyz_path) if optimized_xyz_path else None,
        "l2_energy_hartree": _safe_float(getattr(result, "l2_energy", None)),
        "opt_energy_hartree": _safe_float(getattr(result, "opt_energy", None)),
        "imaginary_count": int(getattr(result, "imaginary_count", 0) or 0),
        "method_used": getattr(result, "method_used", None),
        "error_message": getattr(result, "error_message", None),
        "freq_log": str(Path(result.freq_log).resolve()) if getattr(result, "freq_log", None) else None,
        "log_file": str(Path(result.log_file).resolve()) if getattr(result, "log_file", None) else None,
        "fchk_file": str(Path(result.fchk_file).resolve()) if getattr(result, "fchk_file", None) else None,
        "l2_sp_log": str(Path(result.l2_sp_result.output_file).resolve()) if getattr(result, "l2_sp_result", None) and getattr(result.l2_sp_result, "output_file", None) else None,
    }
    # Canonical: record hash of ORCA sibling .xyz for cross-validation
    freq_log_val = getattr(result, "freq_log", None)
    if freq_log_val:
        freq_path = Path(freq_log_val).resolve()
        sibling = freq_path.with_suffix(".xyz")
        if sibling.exists():
            payload["freq_log_sibling_xyz_sha256"] = _sha256_file(sibling)
    l2_sp_xyz = getattr(result, "l2_sp_result", None)
    sp_input = getattr(l2_sp_xyz, "output_file", None) if l2_sp_xyz else None
    if sp_input and Path(sp_input).exists():
        sp_parent = Path(sp_input).parent
        xyz_candidates = list(sp_parent.glob("*.xyz"))
        if xyz_candidates:
            payload["sp_input_xyz_sha256"] = _sha256_file(xyz_candidates[0])
    return payload


def _metric_kcal(upper: Optional[float], lower: Optional[float]) -> Optional[float]:
    if upper is None or lower is None:
        return None
    return (upper - lower) * HARTREE_TO_KCAL


def _run_geo_method(
    session_dir: Path,
    rx_id: str,
    method_id: str,
    sp_winner: str,
    geo_method_spec: Mapping[str, Any],
    defaults_payload: Mapping[str, Any],
    sp_method_spec: Mapping[str, Any],
    solvent_cfg: Mapping[str, Any],
    inputs: Mapping[str, Path],
    baseline_manifest: Mapping[str, Any],
    force_clean: bool,
    sp_reader_policy: str = "fixed_default",
    phase1_winner: Optional[str] = None,
    shermo_winner: Optional[str] = None,
) -> Dict[str, Any]:
    run_dir = resolve_method_run_dir(session_dir, rx_id, "phase2", method_id)
    if force_clean and run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    config = _build_geo_config(defaults_payload, geo_method_spec, sp_method_spec, solvent_cfg)
    runner = QCTaskRunner(config)
    started = time.time()

    opt_cfg = config.get('theory', {}).get('optimization', {})
    sp_cfg = config.get('theory', {}).get('single_point', {})
    logger.info(
        "rx=%-2s | %-6s | %s/%s/%-6s → %s | SP: %s | %s(%s)",
        rx_id, method_id,
        opt_cfg.get('method','?'), opt_cfg.get('basis','?'), opt_cfg.get('dispersion',''),
        opt_cfg.get('engine', '?').upper(),
        sp_cfg.get('method', '?'),
        solvent_cfg.get('model', 'CPCM'), solvent_cfg.get('solvent', 'acetone'),
    )

    def _opt_label(res: Any, label: str) -> str:
        if not res.converged:
            n = int(getattr(res, 'imaginary_count', 0) or 0)
            method_used = str(getattr(res, 'method_used', '') or '')
            if "ConvergedToMinimum" in method_used:
                return f"FAIL (0 imag, converged to minimum)"
            if "WrongSaddlePoint" in method_used:
                return f"FAIL ({n} imag, wrong saddle point)"
            return f"FAIL (not converged)"
        n = int(getattr(res, 'imaginary_count', 0) or 0)
        if label == "TS":
            freqs = getattr(res, 'frequencies', None)
            if n == 1:
                imag = float(np.min(freqs)) if freqs is not None else 0.0
                return f"OK  (1 imag, ω₁={imag:.1f})"
            else:
                return f"FAIL ({n} imag)"
        else:
            return "OK" if n == 0 else f"WARN ({n} imag)"

    t0 = time.time()
    complex_result = runner.run_opt_sp_cycle(inputs["complex"], run_dir / "opt" / "complex", enable_l2_sp=True)
    t_c = time.time() - t0
    logger.info("  [1/3] %-10s opt+SP → %-28s %.1fs", "complex", _opt_label(complex_result, "complex"), t_c)

    t0 = time.time()
    product_result = runner.run_opt_sp_cycle(inputs["product"], run_dir / "opt" / "product", enable_l2_sp=True)
    t_p = time.time() - t0
    logger.info("  [2/3] %-10s opt+SP → %-28s %.1fs", "product", _opt_label(product_result, "product"), t_p)

    t0 = time.time()
    ts_result = runner.run_ts_opt_cycle(inputs["ts_guess"], run_dir / "opt" / "ts", enable_l2_sp=True)
    t_ts = time.time() - t0
    logger.info("  [3/3] %-10s OptTS  → %-28s %.1fs   input: ts_guess.xyz", "TS", _opt_label(ts_result, "TS"), t_ts)
    ts_fallback_used = False

    g_complex, complex_thermo_path, complex_thermo_error = _thermo_value(complex_result, config, run_dir / "freq" / "complex")
    g_product, product_thermo_path, product_thermo_error = _thermo_value(product_result, config, run_dir / "freq" / "product")
    g_ts, ts_thermo_path, ts_thermo_error = _thermo_value(ts_result, config, run_dir / "freq" / "ts")

    complex_displacement = _geometry_displacement(inputs["complex"], getattr(complex_result, "optimized_xyz", None))
    product_displacement = _geometry_displacement(inputs["product"], getattr(product_result, "optimized_xyz", None))
    ts_displacement = _geometry_displacement(
        inputs["ts_final"],
        getattr(ts_result, "optimized_xyz", None),
    )

    forming_bonds = baseline_manifest.get("forming_bonds", [])
    ts_validation = _validate_ts(Path(ts_result.optimized_xyz), ts_result.frequencies, forming_bonds) if getattr(ts_result, "optimized_xyz", None) else {
        "status": "failed",
        "n_imag": None,
        "imag_freq": None,
        "mode_valid": False,
        "bond_lengths": [],
        "error": "missing_optimized_ts_xyz",
    }

    def _hard_keyword_error(*messages: Optional[str]) -> bool:
        joined = "\n".join(str(m or "") for m in messages).upper()
        return "UNRECOGNIZED" in joined or "DUPLICATED KEYWORD" in joined

    hard_failure = _hard_keyword_error(
        getattr(complex_result, "error_message", None),
        getattr(product_result, "error_message", None),
        getattr(ts_result, "error_message", None),
    )

    if complex_result.converged and product_result.converged and ts_result.converged:
        status = "completed"
    elif hard_failure:
        status = "failed"
    else:
        status = "partial_failed"

    runtime_seconds = time.time() - started

    ts_method_used = str(getattr(ts_result, 'method_used', '') or '')
    ts_imag_count = int(getattr(ts_result, 'imaginary_count', 0) or 0)

    if status == "partial_failed" and "ConvergedToMinimum" in ts_method_used:
        suffix = " | TS converged to minimum (0 imaginary)"
    elif status == "partial_failed" and "WrongSaddlePoint" in ts_method_used:
        suffix = f" | TS converged to wrong saddle point ({ts_imag_count} imaginary)"
    else:
        suffix = ""

    logger.info(
        "  STATUS: %-14s | total: %.1fs%s",
        status, runtime_seconds, suffix,
    )
    complex_payload = _result_payload(complex_result)
    complex_payload["displacement"] = complex_displacement
    product_payload = _result_payload(product_result)
    product_payload["displacement"] = product_displacement
    ts_payload = _result_payload(ts_result)
    ts_payload["displacement"] = ts_displacement

    for label, station_payload in [("complex", complex_payload), ("product", product_payload), ("ts", ts_payload)]:
        opt_xyz_str = station_payload.get("optimized_xyz")
        if opt_xyz_str and status == "completed":
            opt_xyz_path = Path(opt_xyz_str)
            try:
                opt_xyz_path.relative_to(run_dir)
            except ValueError:
                logger.warning(
                    "rx%s | %s | %s optimized_xyz outside run_dir: %s (run_dir=%s)",
                    rx_id, method_id, label, opt_xyz_str, run_dir,
                )
                station_payload["optimized_xyz_under_method_run_dir"] = False
            else:
                station_payload["optimized_xyz_under_method_run_dir"] = True

    payload: Dict[str, Any] = {
        "method_id": method_id,
        "status": status,
        "runtime_seconds": runtime_seconds,
        "sp_reader": sp_winner,
        "sp_reader_method_id": sp_winner,
        "sp_reader_method": sp_method_spec.get("method"),
        "sp_reader_policy": sp_reader_policy,
        "phase1_winner": phase1_winner,
        "shermo_winner": shermo_winner,
        "geo_normalized_spec": config.get("theory", {}).get("optimization", {}).get("normalized_spec"),
        "sp_normalized_spec": config.get("theory", {}).get("single_point", {}).get("normalized_spec"),
        "complex": complex_payload,
        "product": product_payload,
        "ts": ts_payload,
        "ts_fallback_used": ts_fallback_used,
        "ts_validation": ts_validation,
        "g_complex_hartree": g_complex,
        "g_product_hartree": g_product,
        "g_ts_hartree": g_ts,
        "thermo_outputs": {
            "complex": complex_thermo_path,
            "product": product_thermo_path,
            "ts": ts_thermo_path,
        },
        "thermo_errors": {
            "complex": complex_thermo_error,
            "product": product_thermo_error,
            "ts": ts_thermo_error,
        },
        "dE_activation_kcal": _metric_kcal(_safe_float(getattr(ts_result, "l2_energy", None)), _safe_float(getattr(complex_result, "l2_energy", None))),
        "dE_reaction_kcal": _metric_kcal(_safe_float(getattr(product_result, "l2_energy", None)), _safe_float(getattr(complex_result, "l2_energy", None))),
        "dG_activation_kcal": _metric_kcal(g_ts, g_complex),
        "dG_reaction_kcal": _metric_kcal(g_product, g_complex),
        "hard_failure": hard_failure,
    }
    # P0 artifact integrity: warn when optimized_xyz diverges from freq_log sibling
    integrity_ok = True
    for label, station in [("complex", complex_payload), ("product", product_payload), ("ts", ts_payload)]:
        opt_hash = station.get("optimized_xyz_sha256")
        sib_hash = station.get("freq_log_sibling_xyz_sha256")
        if opt_hash and sib_hash and opt_hash != sib_hash:
            logger.warning(
                "rx%s | %s | %s optimized_xyz_sha256 != freq_log_sibling_xyz_sha256 "
                "(opt=%s… sib=%s…) — geometry artifact mismatch; L2-SP/dE may use wrong xyz",
                rx_id, method_id, label, opt_hash[:10], sib_hash[:10],
            )
            integrity_ok = False
    payload["artifact_integrity_ok"] = integrity_ok
    result_path = run_dir / "geo_result.json"
    _write_json(result_path, payload)
    record_task_status(
        session_dir,
        rx_id,
        "geo",
        method_id,
        status=payload["status"],
        run_dir=run_dir,
        exit_code=0 if payload["status"] == "completed" else 1,
        runtime_seconds=runtime_seconds,
        extra={"result_json": str(result_path.resolve())},
    )
    return payload


def run_geo_benchmark_stage(
    session_dir: Path,
    rx_ids: Sequence[str],
    *,
    methods_config: Path = METHODS_GEO_CONFIG,
    defaults_config: Path = DEFAULTS_CONFIG,
    selected_methods: Optional[Sequence[str]] = None,
    force_clean: bool = False,
    requested_sp_reader: Optional[str] = None,
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    benchmark_manifest = read_benchmark_manifest(session_dir)
    defaults_payload = _load_defaults(defaults_config)
    methods, solvent_cfg = _resolve_geo_methods(methods_config, selected_methods)
    if not methods:
        raise RuntimeError("No GEO benchmark methods resolved")
    create_session_manifest(
        session_dir,
        methods={"geo": methods},
        config_snapshot={"defaults": str(defaults_config), "methods_geo": str(methods_config)},
    )
    set_stage_status(session_dir, "phase2_geo", "running")

    logger.info("GEO benchmark: %d method(s) × %d reaction(s) = %d task(s)",
                len(methods), len(rx_ids), len(methods) * len(rx_ids))
    outputs: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for rx_id in rx_ids:
        init_rx_manifest(session_dir, rx_id)
        rx_manifest = read_rx_manifest(session_dir, rx_id)
        _, stationary_points, baseline_manifest = _baseline_context(rx_manifest)
        sp_winner, sp_method_spec = _resolve_sp_reader_spec(benchmark_manifest, requested_sp_reader)
        inputs = _point_inputs(stationary_points)
        sp_reader_policy = "fixed_default" if requested_sp_reader is None else "explicit"
        current_ph1_winner = _phase1_winner_method_id(rx_manifest)
        current_sh_winner = rx_manifest.get("current_shermo_winner")
        outputs[rx_id] = {}
        for method_id, method_spec in methods.items():
            outputs[rx_id][method_id] = _run_geo_method(
                session_dir,
                rx_id,
                method_id,
                sp_winner,
                method_spec,
                defaults_payload,
                sp_method_spec,
                solvent_cfg,
                inputs,
                baseline_manifest,
                force_clean,
                sp_reader_policy=sp_reader_policy,
                phase1_winner=current_ph1_winner,
                shermo_winner=current_sh_winner,
            )
            if bool(outputs[rx_id][method_id].get("hard_failure")):
                set_stage_status(session_dir, "phase2_geo", "failed")
                raise RuntimeError(
                    f"GEO benchmark hard failure at rx{rx_id}/{method_id}: {outputs[rx_id][method_id].get('ts', {}).get('error_message', 'unknown_error')}"
                )
        updated_rx_manifest = read_rx_manifest(session_dir, rx_id)
        updated_rx_manifest.setdefault("phase2", {})
        updated_rx_manifest["phase2"]["sp_reader"] = sp_winner
        write_rx_manifest(session_dir, rx_id, updated_rx_manifest)
    set_stage_status(session_dir, "phase2_geo", "computed")
    return outputs


def _choose_phase2_winner(results: Mapping[str, Mapping[str, Any]]) -> Optional[str]:
    successful = {
        k: v
        for k, v in results.items()
        if str(v.get("status", "")).lower() in {"ok", "completed"}
    }
    if not successful:
        return None
    scored: List[Tuple[str, float]] = []
    for method_id, row in successful.items():
        score = 0.0
        ts_validation = row.get("ts_validation", {})
        ts_n_imag = _safe_float(ts_validation.get("n_imag") if isinstance(ts_validation, dict) else None)
        if ts_n_imag == 1.0:
            score += 0.35
        bond_lengths = ts_validation.get("bond_lengths", []) if isinstance(ts_validation, dict) else []
        if isinstance(bond_lengths, list) and bond_lengths and all(2.0 <= float(value) <= 2.4 for value in bond_lengths):
            score += 0.25
        g_act = _safe_float(row.get("dG_activation_kcal"))
        if g_act is not None:
            score += 0.25
        runtime = _safe_float(row.get("runtime_seconds"))
        if runtime is not None:
            score += max(0.0, 0.15 - min(runtime / 100000.0, 0.15))
        scored.append((method_id, score))
    scored.sort(key=lambda item: item[1], reverse=True)
    return scored[0][0] if scored else None


def _build_comparison_rankings(
    results: Mapping[str, Mapping[str, Any]],
) -> Dict[str, List[str]]:
    completed = {
        k: v
        for k, v in results.items()
        if str(v.get("status", "")).lower() in {"ok", "completed"}
    }

    def _get_ts_disp(key: str) -> Optional[float]:
        ts = completed.get(key, {}).get("ts", {})
        if isinstance(ts, dict):
            disp = ts.get("displacement", {})
            if isinstance(disp, dict):
                return disp.get("aligned_rmsd_to_baseline_angstrom")
        return None

    by_runtime = sorted(
        completed.keys(),
        key=lambda k: _safe_float(completed[k].get("runtime_seconds")) or float("inf"),
    )
    by_rmsd = sorted(
        [k for k in completed.keys() if _get_ts_disp(k) is not None],
        key=lambda k: _get_ts_disp(k) or float("inf"),
    )
    by_dg = sorted(
        [k for k in completed.keys() if _safe_float(completed[k].get("dG_activation_kcal")) is not None],
        key=lambda k: abs(
            (_safe_float(completed[k].get("dG_activation_kcal")) or 0.0)
            - (_safe_float(completed.get("GEO-1", {}).get("dG_activation_kcal")) or 0.0)
        ),
    )

    return {
        "by_runtime": by_runtime,
        "by_aligned_rmsd": by_rmsd,
        "by_delta_g_activation_vs_geo1": by_dg,
    }


def _write_phase2_global_comparison(
    session_dir: Path,
    rx_ids: Sequence[str],
    summaries: Mapping[str, Mapping[str, Any]],
) -> None:
    benchmark_manifest = read_benchmark_manifest(session_dir)
    benchmark_manifest.setdefault("reports", {})
    benchmark_manifest["reports"]["phase2_global_summary_json"] = str(
        (session_dir / "reports" / "phase2_global_summary.json").resolve()
    )
    benchmark_manifest["reports"]["phase2_global_summary_md"] = str(
        (session_dir / "reports" / "phase2_global_summary.md").resolve()
    )
    write_benchmark_manifest(session_dir, benchmark_manifest)

    global_summary = {
        "phase": "geo",
        "benchmark_mode": "fixed_sp1_opt_comparison",
        "ranking_policy": "comparison_only_no_winner",
        "rx": dict(summaries),
    }
    _write_json(
        session_dir / "reports" / "phase2_global_summary.json",
        global_summary,
    )

    lines = [
        "# OPT Benchmark — Global Comparison (fixed SP-1)",
        "",
        "| rx_id | sp_reader | methods_completed |",
        "|---|---:|---:|",
    ]
    for rx_id in rx_ids:
        summary = summaries.get(rx_id, {})
        sp = summary.get("sp_reader", "NA")
        results = summary.get("results", {})
        n_done = sum(
            1
            for r in results.values()
            if isinstance(r, dict)
            and str(r.get("status", "")).lower() in {"ok", "completed"}
        )
        lines.append(f"| rx{rx_id} | {sp} | {n_done} / {len(results)} |")

    lines.append("")
    lines.append("Per-reaction rankings: see rx*/reports/geo_summary.md")
    _write_text(
        session_dir / "reports" / "phase2_global_summary.md",
        "\n".join(lines) + "\n",
    )


def _render_rx_markdown(rx_id: str, winner: Optional[str], sp_reader: Optional[str], results: Mapping[str, Mapping[str, Any]]) -> str:
    lines = [
        f"# GEO Benchmark Summary — rx{rx_id}",
        "",
        f"- sp_reader: {sp_reader or 'NA'}",
        f"- benchmark_mode: fixed_sp1_opt_comparison",
        f"- ranking_policy: comparison_only_no_winner",
        "",
        "| method | status | dG‡ (kcal) | dG_rxn (kcal) | TS align RMSD (Å) | TS max disp (Å) | runtime(s) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for method_id, row in results.items():
        ts = row.get("ts", {}) if isinstance(row.get("ts"), dict) else {}
        ts_disp = ts.get("displacement", {}) if isinstance(ts, dict) else {}
        lines.append(
            "| {method} | {status} | {dga} | {dgr} | {rmsd} | {maxd} | {runtime} |".format(
                method=method_id,
                status=row.get("status", "unknown"),
                dga=_format_metric(row.get("dG_activation_kcal")),
                dgr=_format_metric(row.get("dG_reaction_kcal")),
                rmsd=_format_metric(ts_disp.get("aligned_rmsd_to_baseline_angstrom") if isinstance(ts_disp, dict) else None),
                maxd=_format_metric(ts_disp.get("max_atom_displacement_angstrom") if isinstance(ts_disp, dict) else None),
                runtime=_format_metric(row.get("runtime_seconds")),
            )
        )
    return "\n".join(lines) + "\n"


def _format_metric(value: Any) -> str:
    parsed = _safe_float(value)
    return "NA" if parsed is None else f"{parsed:.3f}"


def evaluate_geo_benchmark_stage(session_dir: Path, rx_ids: Sequence[str]) -> Dict[str, Dict[str, Any]]:
    benchmark_manifest = read_benchmark_manifest(session_dir)
    summaries: Dict[str, Dict[str, Any]] = {}
    for rx_id in rx_ids:
        rx_manifest = read_rx_manifest(session_dir, rx_id)
        tasks = rx_manifest.get("tasks", {}).get("phase2", {})
        if not isinstance(tasks, dict):
            tasks = {}
        results: Dict[str, Dict[str, Any]] = {}
        for method_id, task_row in tasks.items():
            result_path = Path(str((task_row or {}).get("result_json", ""))).expanduser()
            if result_path.is_file():
                result_payload = _read_json(result_path)
                require_phase2_canonical_results(result_payload, context=f"GEO benchmark result rx{rx_id}/{method_id}")
                results[str(method_id)] = result_payload
            else:
                results[str(method_id)] = {
                    "method_id": str(method_id),
                    "status": "failed",
                    "error_message": "missing_geo_result_json",
                }
        sp_reader = _phase2_sp_reader_method_id(rx_manifest, results)
        phase1_results = rx_manifest.get("phase1", {}).get("results", {})
        if not isinstance(phase1_results, Mapping):
            raise ValueError(f"GEO benchmark requires canonical phase1 results for rx{rx_id}")
        for method_id, result_payload in phase1_results.items():
            if isinstance(result_payload, Mapping):
                require_phase1_canonical_results(result_payload, context=f"GEO benchmark phase1 cache rx{rx_id}/{method_id}")

        rankings = _build_comparison_rankings(results)

        completed_results = {
            k: v for k, v in results.items()
            if str(v.get("status", "")).lower() in {"ok", "completed"}
        }
        if len(completed_results) >= 2:
            ts_hashes: Dict[str, str] = {}
            for mid, row in completed_results.items():
                ts_block = row.get("ts", {})
                if isinstance(ts_block, dict):
                    h = ts_block.get("optimized_xyz_sha256")
                    if h:
                        ts_hashes[mid] = h
            if len(set(ts_hashes.values())) == 1 and len(ts_hashes) >= 2:
                logger.warning(
                    "rx%s: ALL %d completed GEO methods share identical TS optimized_xyz_sha256=%s. "
                    "geometry_equivalence_unverified — cannot conclude methods produce equivalent energies independently.",
                    rx_id, len(ts_hashes), next(iter(ts_hashes.values())),
                )
                summary_json_writable = True
            else:
                summary_json_writable = False
                unique_hashes = set(ts_hashes.values())
                if len(unique_hashes) < len(ts_hashes):
                    dup_groups: Dict[str, List[str]] = {}
                    for mid, h in ts_hashes.items():
                        dup_groups.setdefault(h, []).append(mid)
                    for h, group in dup_groups.items():
                        if len(group) >= 2:
                            logger.warning(
                                "rx%s: GEO methods %s share identical TS xyz hash %s",
                                rx_id, ", ".join(group), h[:16],
                            )
        else:
            summary_json_writable = False

        set_phase_winner(session_dir, rx_id, "phase2", (rankings.get("by_runtime") or [None])[0], sp_reader=sp_reader)
        updated_rx_manifest = read_rx_manifest(session_dir, rx_id)
        updated_rx_manifest.setdefault("phase2", {})
        updated_rx_manifest["phase2"]["results"] = results
        report_paths = resolve_rx_phase_report_paths(session_dir, rx_id, "phase2")
        summary_json = {
            "rx_id": rx_id,
            "phase": "geo",
            "benchmark_mode": "fixed_sp1_opt_comparison",
            "ranking_policy": "comparison_only_no_winner",
            "sp_reader": sp_reader,
            "sp_reader_method_id": sp_reader,
            "geometry_equivalence_unverified": bool(summary_json_writable),
            "results": results,
            "rankings": rankings,
        }
        _write_json(report_paths["json"], summary_json)
        _write_text(report_paths["md"], _render_rx_markdown(rx_id, None, sp_reader, results))
        updated_rx_manifest.setdefault("reports", {})
        updated_rx_manifest["reports"]["geo_summary_json"] = str(report_paths["json"].resolve())
        updated_rx_manifest["reports"]["geo_summary_md"] = str(report_paths["md"].resolve())
        write_rx_manifest(session_dir, rx_id, updated_rx_manifest)
        summaries[rx_id] = summary_json

    _write_phase2_global_comparison(session_dir, rx_ids, summaries)
    set_stage_status(session_dir, "phase2_geo", "done")
    return summaries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run or evaluate GEO benchmark stage.")
    parser.add_argument("--session-dir", required=True, help="Benchmark session directory")
    parser.add_argument("--mode", choices=["run", "evaluate", "all"], default="all", help="Stage mode")
    parser.add_argument("--rx-id", nargs="*", default=None, help="Optional subset of rx ids")
    parser.add_argument("--methods", nargs="*", default=None, help="Optional subset of GEO method ids")
    parser.add_argument("--methods-config", default=str(METHODS_GEO_CONFIG), help="GEO methods YAML")
    parser.add_argument("--defaults-config", default=str(DEFAULTS_CONFIG), help="Baseline defaults YAML")
    parser.add_argument("--force-clean", action="store_true", help="Remove existing per-method directories before running")
    parser.add_argument("--sp-reader", default=None, help="Force a specific SP reader method id (default: SP-1)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    session_dir = Path(args.session_dir).expanduser().resolve()
    rx_ids = _resolve_rx_ids(session_dir, args.rx_id)
    if not rx_ids:
        raise RuntimeError("No rx ids resolved for GEO benchmark")
    if args.mode in {"run", "all"}:
        run_geo_benchmark_stage(
            session_dir,
            rx_ids,
            methods_config=Path(args.methods_config).expanduser().resolve(),
            defaults_config=Path(args.defaults_config).expanduser().resolve(),
            selected_methods=args.methods,
            force_clean=bool(args.force_clean),
            requested_sp_reader=args.sp_reader,
        )
    if args.mode in {"evaluate", "all"}:
        evaluate_geo_benchmark_stage(session_dir, rx_ids)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
