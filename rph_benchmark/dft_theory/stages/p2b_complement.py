#!/usr/bin/env python3
"""
Phase 2b — Intermediate + Precursor Optimization Complement
============================================================

Rationale:
  GEO Phase 2 only optimizes 3 of 5 key stationary points (complex, product, TS).
  Precursor and intermediate remain at baseline geometry throughout the benchmark.
  This makes cross-method dE comparisons incomplete — we cannot know if barrier
  changes come from TS geometry variation or from the reference energy shift.

  Phase 2b completes the profile by optimizing precursor and intermediate with
  the same DFT methods used for complex/product/TS.  After Phase 2b, every
  (rx, method) combination has a full 5-point optimized profile.

Stationary points:
  Phase 2a (geo_benchmark):  complex, product, ts
  Phase 2b (this script):    precursor, intermediate

Computational cost:
  2 OPT+Freq+L2-SP per (rx, method) = 4 rx x 5 methods x 2 = 40 jobs

Output:
  p2b_result.json in each method run directory, with schema mirroring
  geo_result.json but for precursor/intermediate only.

Usage:
  # Run
  python benchmark/dft_theory/stages/p2b_complement.py \\
      --session-dir Output/.../session_opt_full_orca_sp1 \\
      --mode run --rx-id 1 3 8 15

  # Evaluate
  python benchmark/dft_theory/stages/p2b_complement.py \\
      --session-dir Output/.../session_opt_full_orca_sp1 \\
      --mode evaluate --rx-id 1 3 8 15
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import logging as _logging

_logging.basicConfig(level=_logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

from benchmark.dft_theory.lib.state_manager import (
    init_rx_manifest,
    read_benchmark_manifest,
    read_rx_manifest,
    record_task_status,
    resolve_method_run_dir,
    require_canonical_baseline_manifest,
    require_canonical_stationary_points,
    set_stage_status,
    write_rx_manifest,
)
from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.utils.file_io import read_xyz
from rph_core.utils.optimization_config import prepare_qc_method_config
from rph_core.utils.qc_task_runner import QCTaskRunner
from rph_core.utils.shermo_runner import run_shermo

DEFAULTS_CONFIG = PROJECT_ROOT / "config" / "defaults.yaml"
METHODS_GEO_CONFIG = PROJECT_ROOT / "benchmark" / "dft_theory" / "config" / "methods_geo.yaml"
P2B_POINTS = ("precursor", "intermediate")

# ────────────────────────── helpers ──────────────────────────


def _read_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid JSON payload: {path}")
    return payload


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _sha256_file(path: Path) -> Optional[str]:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


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


def _metric_kcal(upper: Optional[float], lower: Optional[float]) -> Optional[float]:
    if upper is None or lower is None:
        return None
    return (upper - lower) * HARTREE_TO_KCAL


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _result_payload(result: Any) -> Dict[str, Any]:
    optimized_xyz_path = (
        Path(result.optimized_xyz).resolve()
        if getattr(result, "optimized_xyz", None)
        else None
    )
    payload: Dict[str, Any] = {
        "converged": bool(getattr(result, "converged", False)),
        "optimized_xyz": str(optimized_xyz_path) if optimized_xyz_path else None,
        "optimized_xyz_sha256": (
            _sha256_file(optimized_xyz_path) if optimized_xyz_path else None
        ),
        "l2_energy_hartree": _safe_float(getattr(result, "l2_energy", None)),
        "opt_energy_hartree": _safe_float(getattr(result, "opt_energy", None)),
        "imaginary_count": int(getattr(result, "imaginary_count", 0) or 0),
        "method_used": getattr(result, "method_used", None),
        "error_message": getattr(result, "error_message", None),
        "freq_log": (
            str(Path(result.freq_log).resolve())
            if getattr(result, "freq_log", None)
            else None
        ),
        "log_file": (
            str(Path(result.log_file).resolve())
            if getattr(result, "log_file", None)
            else None
        ),
        "fchk_file": (
            str(Path(result.fchk_file).resolve())
            if getattr(result, "fchk_file", None)
            else None
        ),
        "l2_sp_log": (
            str(Path(result.l2_sp_result.output_file).resolve())
            if getattr(result, "l2_sp_result", None)
            and getattr(result.l2_sp_result, "output_file", None)
            else None
        ),
    }
    return payload


def _thermo_value(
    result: Any,
    config: Mapping[str, Any],
    output_dir: Path,
) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    l2_energy = _safe_float(getattr(result, "l2_energy", None))
    freq_log = (
        getattr(result, "freq_log", None)
        or getattr(result, "log_file", None)
        or getattr(result, "qm_output_file", None)
    )
    if l2_energy is None or freq_log is None:
        return None, None, "missing_l2_or_freq_log"

    executables = (
        config.get("executables", {})
        if isinstance(config.get("executables"), dict)
        else {}
    )
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
        thermo = run_shermo(
            shermo_bin, Path(freq_log), float(l2_energy), output_dir / "thermo.sum"
        )
    except Exception as exc:
        return None, None, str(exc)
    value = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
    return (
        float(value),
        str(getattr(thermo, "output_file", output_dir / "thermo.sum")),
        None,
    )


# ────────────────────────── config resolution ──────────────────────────


def _resolve_geo_methods(
    methods_config: Path,
    requested: Optional[Sequence[str]],
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Any]]:
    payload = _load_yaml(methods_config)
    methods = payload.get("methods", {}) or {}
    if not isinstance(methods, dict):
        raise ValueError(f"Invalid methods block: {methods_config}")
    if not requested:
        selected = {
            str(k): dict(v) for k, v in methods.items() if isinstance(v, dict)
        }
    else:
        requested_set = {str(s).strip() for s in requested if str(s).strip()}
        selected = {
            str(k): dict(v)
            for k, v in methods.items()
            if str(k) in requested_set and isinstance(v, dict)
        }
    solvent = payload.get("solvent", {}) or {}
    return selected, dict(solvent) if isinstance(solvent, dict) else {}


def _resolve_rx_ids(
    session_dir: Path,
    requested: Optional[Sequence[str]],
) -> List[str]:
    manifest_path = session_dir / "manifests" / "benchmark_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"No benchmark manifest: {manifest_path}")
    manifest = _read_json(manifest_path)
    configured = [str(c) for c in manifest.get("cases", []) if str(c).strip()]
    if not requested:
        return configured
    selected = {str(s).strip() for s in requested if str(s).strip()}
    return [c for c in configured if c in selected]


def _resolve_sp_reader_spec(
    benchmark_manifest: Mapping[str, Any],
    requested: Optional[str],
    geo_result_sp_reader: Optional[str],
) -> Tuple[str, Dict[str, Any]]:
    sp_reader = requested or geo_result_sp_reader or "SP-1"
    sp_methods = benchmark_manifest.get("methods", {}).get("sp", {})
    if not isinstance(sp_methods, dict):
        raise ValueError("Missing methods.sp in benchmark manifest")
    spec = sp_methods.get(sp_reader)
    if not isinstance(spec, dict):
        raise ValueError(
            f"No SP method spec for '{sp_reader}' in benchmark manifest"
        )
    return sp_reader, dict(spec)


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


def _baseline_context(
    rx_manifest: Mapping[str, Any],
) -> Tuple[Path, Dict[str, str], Dict[str, Any]]:
    baseline = rx_manifest.get("baseline")
    if not isinstance(baseline, dict):
        raise ValueError("Missing baseline block in rx manifest")
    baseline_root = Path(str(baseline.get("root", ""))).expanduser().resolve()
    stationary_points = baseline.get("stationary_points", {})
    if not isinstance(stationary_points, dict):
        raise ValueError("Missing baseline stationary_points")
    require_canonical_stationary_points(
        stationary_points, context="P2B complement baseline block"
    )
    baseline_manifest_path = (
        Path(str(baseline.get("baseline_manifest", ""))).expanduser().resolve()
    )
    if not baseline_manifest_path.is_file():
        raise FileNotFoundError(
            f"Missing baseline manifest: {baseline_manifest_path}"
        )
    baseline_manifest = _read_json(baseline_manifest_path)
    require_canonical_baseline_manifest(
        baseline_manifest, context="P2B complement"
    )
    return (
        baseline_root,
        {str(k): str(v) for k, v in stationary_points.items()},
        baseline_manifest,
    )


def _p2b_point_inputs(
    stationary_points: Mapping[str, str],
) -> Dict[str, Path]:
    require_canonical_stationary_points(
        stationary_points, context="P2B stationary_points"
    )
    precursor_xyz = stationary_points.get("precursor_min.xyz")
    intermediate_xyz = stationary_points.get("intermediate.xyz")
    if not precursor_xyz or not intermediate_xyz:
        raise ValueError(
            "Incomplete baseline stationary point set for P2B complement"
        )
    return {
        "precursor": Path(precursor_xyz).expanduser().resolve(),
        "intermediate": Path(intermediate_xyz).expanduser().resolve(),
    }


# ────────────────────────── geometry displacement ──────────────────────────


def _geometry_displacement(
    baseline_xyz: Path,
    optimized_xyz: Optional[str],
) -> Dict[str, Any]:
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
        result["aligned_rmsd_to_baseline_angstrom"] = float(
            np.sqrt(np.mean(np.sum((cb_centered - co_aligned) ** 2, axis=1)))
        )
        result["mean_atom_displacement_angstrom"] = float(deltas.mean())
        result["max_atom_displacement_angstrom"] = float(deltas.max())
    except Exception as exc:
        result["error"] = str(exc)
    return result


# ────────────────────────── run ──────────────────────────


def _run_p2b_method(
    session_dir: Path,
    rx_id: str,
    method_id: str,
    geo_method_spec: Mapping[str, Any],
    defaults_payload: Mapping[str, Any],
    sp_method_spec: Mapping[str, Any],
    solvent_cfg: Mapping[str, Any],
    p2b_inputs: Mapping[str, Path],
    force_clean: bool,
) -> Dict[str, Any]:
    run_dir = resolve_method_run_dir(session_dir, rx_id, "phase2", method_id)
    if force_clean:
        for pt_name in P2B_POINTS:
            target = run_dir / "opt" / pt_name
            if target.exists():
                shutil.rmtree(target)

    run_dir.mkdir(parents=True, exist_ok=True)
    config = _build_geo_config(
        defaults_payload, geo_method_spec, sp_method_spec, solvent_cfg
    )
    runner = QCTaskRunner(config)

    opt_cfg = config.get("theory", {}).get("optimization", {})
    sp_cfg = config.get("theory", {}).get("single_point", {})
    logger.info(
        "rx=%-2s | %-6s | %s/%s/%s → %s | SP: %s",
        rx_id,
        method_id,
        opt_cfg.get("method", "?"),
        opt_cfg.get("basis", "?"),
        opt_cfg.get("dispersion", ""),
        opt_cfg.get("engine", "?").upper(),
        sp_cfg.get("method", "?"),
    )

    started = time.time()
    overall_status = "completed"

    point_results: Dict[str, Dict[str, Any]] = {}
    for pt_name in P2B_POINTS:
        t0 = time.time()
        opt_label = "precursor" if pt_name == "precursor" else "intermediate"
        opt_dir = run_dir / "opt" / pt_name
        result = runner.run_opt_sp_cycle(
            p2b_inputs[pt_name], opt_dir, enable_l2_sp=True, charge=0, spin=1
        )
        runtime = time.time() - t0

        if not result.converged:
            overall_status = "partial_failed"
        n_imag = int(getattr(result, "imaginary_count", 0) or 0)
        status_label = "OK" if result.converged and n_imag == 0 else f"FAIL({n_imag}imag)"

        payload = _result_payload(result)
        payload["runtime_seconds"] = runtime

        disp = _geometry_displacement(
            p2b_inputs[pt_name], getattr(result, "optimized_xyz", None)
        )
        payload["displacement"] = disp

        point_results[pt_name] = payload
        logger.info(
            "  [%s] %-12s opt+SP → %-10s %.1fs  l2=%.6f",
            pt_name[0].upper(),
            opt_label,
            status_label,
            runtime,
            _safe_float(getattr(result, "l2_energy", 0.0)) or 0.0,
        )

    total_runtime = time.time() - started

    # Shermo
    g_precursor: Optional[float] = None
    g_intermediate: Optional[float] = None
    thermo_outputs: Dict[str, Any] = {}
    thermo_errors: Dict[str, Any] = {}

    if "precursor" in point_results:
        g_precursor, thermo_outputs["precursor"], thermo_errors[
            "precursor"
        ] = _thermo_value(
            point_results["precursor"], config, run_dir / "freq" / "precursor"
        )
    if "intermediate" in point_results:
        g_intermediate, thermo_outputs["intermediate"], thermo_errors[
            "intermediate"
        ] = _thermo_value(
            point_results["intermediate"],
            config,
            run_dir / "freq" / "intermediate",
        )

    e_precursor = _safe_float(
        point_results.get("precursor", {}).get("l2_energy_hartree")
    )
    e_intermediate = _safe_float(
        point_results.get("intermediate", {}).get("l2_energy_hartree")
    )

    payload: Dict[str, Any] = {
        "method_id": method_id,
        "status": overall_status,
        "runtime_seconds": total_runtime,
        "sp_reader": sp_method_spec.get("method"),
        "geo_method": geo_method_spec.get("method"),
        "geo_basis": geo_method_spec.get("basis"),
        "precursor": point_results.get("precursor", {}),
        "intermediate": point_results.get("intermediate", {}),
        "dE_precursor_to_intermediate_kcal": _metric_kcal(
            e_intermediate, e_precursor
        ),
        "g_precursor_hartree": g_precursor,
        "g_intermediate_hartree": g_intermediate,
        "dG_precursor_to_intermediate_kcal": _metric_kcal(
            g_intermediate, g_precursor
        ),
        "thermo_outputs": thermo_outputs,
        "thermo_errors": thermo_errors,
    }

    result_path = run_dir / "p2b_result.json"
    _write_json(result_path, payload)

    record_task_status(
        session_dir,
        rx_id,
        "p2b",
        method_id,
        status=payload["status"],
        run_dir=run_dir,
        exit_code=0 if payload["status"] == "completed" else 1,
        runtime_seconds=total_runtime,
        extra={"p2b_result_json": str(result_path.resolve())},
    )

    return payload


def run_p2b_stage(
    session_dir: Path,
    rx_ids: Sequence[str],
    *,
    selected_methods: Optional[Sequence[str]] = None,
    force_clean: bool = False,
    sp_reader_override: Optional[str] = None,
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    benchmark_manifest = read_benchmark_manifest(session_dir)
    defaults_payload = _load_yaml(DEFAULTS_CONFIG)
    geo_methods, solvent_cfg = _resolve_geo_methods(
        METHODS_GEO_CONFIG, selected_methods
    )

    if not geo_methods:
        raise RuntimeError("No GEO methods resolved for P2B")

    set_stage_status(session_dir, "p2b_complement", "running")

    logger.info(
        "P2B: %d method(s) x %d reaction(s) x 2 points = %d task(s)",
        len(geo_methods),
        len(rx_ids),
        len(geo_methods) * len(rx_ids) * 2,
    )

    outputs: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for rx_id in rx_ids:
        init_rx_manifest(session_dir, rx_id)
        rx_manifest = read_rx_manifest(session_dir, rx_id)
        _, stationary_points, _ = _baseline_context(rx_manifest)

        first_method = next(iter(geo_methods))
        run_dir = resolve_method_run_dir(session_dir, rx_id, "phase2", first_method)
        geo_result_path = run_dir / "geo_result.json"
        geo_sp_reader = None
        if geo_result_path.exists():
            geo_sp_reader = _read_json(geo_result_path).get("sp_reader")

        sp_reader_id, sp_method_spec = _resolve_sp_reader_spec(
            benchmark_manifest, sp_reader_override, geo_sp_reader
        )
        logger.info("P2B rx%s SP reader: %s", rx_id, sp_reader_id)

        p2b_inputs = _p2b_point_inputs(stationary_points)
        outputs[rx_id] = {}

        total = len(geo_methods)
        done = 0
        for method_id, geo_spec in geo_methods.items():
            done += 1
            logger.info("[%d/%d] rx%s %s ...", done, total, rx_id, method_id)
            outputs[rx_id][method_id] = _run_p2b_method(
                session_dir,
                rx_id,
                method_id,
                geo_spec,
                defaults_payload,
                sp_method_spec,
                solvent_cfg,
                p2b_inputs,
                force_clean,
            )

    set_stage_status(session_dir, "p2b_complement", "computed")
    return outputs


# ────────────────────────── evaluate ──────────────────────────


def _compute_full_profile(
    geo_result: Mapping[str, Any],
    p2b_result: Optional[Mapping[str, Any]],
) -> Dict[str, Any]:
    profile: Dict[str, Any] = {}

    def _point_energy(source: Optional[Mapping[str, Any]], pt_name: str) -> Optional[float]:
        if source is None:
            return None
        pt = source.get(pt_name, {})
        if not isinstance(pt, dict):
            return None
        e = pt.get("l2_energy_hartree")
        if e is not None:
            return _safe_float(e)
        rs = source.get("refresh_sp", {}).get("points", {}).get(pt_name, {})
        return rs.get("refreshed_l2_energy_hartree") if isinstance(rs, dict) else None

    e_precursor = _point_energy(p2b_result, "precursor")
    e_intermediate = _point_energy(p2b_result, "intermediate")
    e_complex = _point_energy(geo_result, "complex")
    e_ts = _point_energy(geo_result, "ts")
    e_product = _point_energy(geo_result, "product")

    profile["E_precursor"] = e_precursor
    profile["E_intermediate"] = e_intermediate
    profile["E_complex"] = e_complex
    profile["E_ts"] = e_ts
    profile["E_product"] = e_product

    profile["dE_precursor_to_intermediate"] = _metric_kcal(
        e_intermediate, e_precursor
    )
    profile["dE_intermediate_to_complex"] = _metric_kcal(e_complex, e_intermediate)
    profile["dE_complex_to_ts"] = _metric_kcal(e_ts, e_complex)
    profile["dE_intermediate_to_ts"] = _metric_kcal(e_ts, e_intermediate)
    profile["dE_precursor_to_ts"] = _metric_kcal(e_ts, e_precursor)

    profile["dE_complex_to_product"] = _metric_kcal(e_product, e_complex)
    profile["dE_intermediate_to_product"] = _metric_kcal(
        e_product, e_intermediate
    )
    profile["dE_precursor_to_product"] = _metric_kcal(e_product, e_precursor)

    return profile


def evaluate_p2b_stage(
    session_dir: Path,
    rx_ids: Sequence[str],
    selected_methods: Optional[Sequence[str]] = None,
) -> Dict[str, Dict[str, Any]]:
    geo_methods, _ = _resolve_geo_methods(METHODS_GEO_CONFIG, selected_methods)
    if not geo_methods:
        raise RuntimeError("No GEO methods resolved for evaluation")

    summaries: Dict[str, Dict[str, Any]] = {}
    for rx_id in rx_ids:
        rx_summary: Dict[str, Dict[str, Any]] = {}
        for method_id in geo_methods:
            run_dir = resolve_method_run_dir(session_dir, rx_id, "phase2", method_id)
            geo_path = run_dir / "geo_result.json"
            p2b_path = run_dir / "p2b_result.json"
            geo_result = _read_json(geo_path) if geo_path.exists() else None
            p2b_result = _read_json(p2b_path) if p2b_path.exists() else None

            if geo_result is None:
                rx_summary[method_id] = {"status": "missing_geo_result"}
                continue

            profile = _compute_full_profile(geo_result, p2b_result)
            profile["status"] = (
                "full" if p2b_result else "geo_only"
            )
            profile["method_id"] = method_id
            profile["opt_method"] = geo_methods[method_id].get("method", "?")
            rx_summary[method_id] = profile
        summaries[rx_id] = rx_summary

    report_json = session_dir / "reports" / "p2b_full_profile.json"
    report_md = session_dir / "reports" / "p2b_full_profile.md"
    _write_json(report_json, {"rx": summaries, "generated_at": _now_iso()})
    _render_p2b_markdown(summaries, report_md, geo_methods)
    logger.info("P2B profile report: %s", report_json)

    return summaries


def _render_p2b_markdown(
    summaries: Mapping[str, Mapping[str, Any]],
    path: Path,
    geo_methods: Mapping[str, Any],
) -> None:
    lines = [
        "# Phase 2b — Full 5-Point Reaction Profile",
        "",
        f"Generated: {_now_iso()}",
        "",
        "All energies in kcal/mol.  dE values computed from L2 SP (SP-1) on GEO-optimized geometries.",
        "",
    ]

    for rx_id in sorted(summaries.keys(), key=int):
        rx_data = summaries[rx_id]
        lines.extend([
            f"## rx{rx_id}",
            "",
            "| method | opt_level | E_prec | E_int | E_complex | E_ts | E_prod | dE‡(int→ts) | dE_rxn(prec→prod) | dE‡(c→ts) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ])

        for method_id in ["GEO-1", "GEO-2", "GEO-3", "GEO-4", "GEO-5"]:
            profile = rx_data.get(method_id, {})
            if not isinstance(profile, dict):
                continue
            opt_method = geo_methods.get(method_id, {}).get("method", "?")

            def f(key: str) -> str:
                v = profile.get(key)
                return f"{v:.3f}" if v is not None else "N/A"

            lines.append(
                f"| {method_id} | {opt_method} | {f('E_precursor')} | "
                f"{f('E_intermediate')} | {f('E_complex')} | {f('E_ts')} | "
                f"{f('E_product')} | {f('dE_intermediate_to_ts')} | "
                f"{f('dE_precursor_to_product')} | {f('dE_complex_to_ts')} |"
            )
        lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ────────────────────────── CLI ──────────────────────────


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Phase 2b: optimize precursor and intermediate for benchmark."
    )
    parser.add_argument(
        "--session-dir", required=True, help="Benchmark session directory"
    )
    parser.add_argument(
        "--mode",
        choices=["run", "evaluate", "all"],
        default="all",
        help="Stage mode",
    )
    parser.add_argument(
        "--rx-id", nargs="*", default=None, help="Optional subset of rx ids"
    )
    parser.add_argument(
        "--methods", nargs="*", default=None, help="Optional subset of GEO method ids"
    )
    parser.add_argument("--force-clean", action="store_true")
    parser.add_argument(
        "--sp-reader", default=None, help="Force SP reader (default: from geo_result.json or SP-1)"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    session_dir = Path(args.session_dir).expanduser().resolve()
    if not session_dir.exists():
        raise FileNotFoundError(f"Session dir not found: {session_dir}")

    rx_ids = _resolve_rx_ids(session_dir, args.rx_id)
    if not rx_ids:
        raise RuntimeError("No rx ids resolved for P2B")

    if args.mode in {"run", "all"}:
        run_p2b_stage(
            session_dir,
            rx_ids,
            selected_methods=args.methods,
            force_clean=bool(args.force_clean),
            sp_reader_override=args.sp_reader,
        )
    if args.mode in {"evaluate", "all"}:
        evaluate_p2b_stage(session_dir, rx_ids, args.methods)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
