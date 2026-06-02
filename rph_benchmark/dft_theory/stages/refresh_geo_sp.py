#!/usr/bin/env python3
"""
GEO Benchmark SP Refresh
=========================
Re-run L2 single-point calculations on the CORRECT optimized geometries
extracted from ORCA output files, then update geo_result.json with
incremental fields (never overwriting original values).

Background:
  The old _write_optimized_xyz() wrote incorrect geometries (baseline copies)
  for all GEO methods.  The actual ORCA optimization produced correct, method-
  specific final geometries stored as ``*_opt_<hash>.xyz`` in the ORCA output
  directories.  This script finds those correct files and re-runs SP on them.

Modes:
  audit   – discover + validate + produce audit report (NO QC execution)
  refresh – re-run SP + write incremental fields to geo_result.json

Design rules (per code review 2026-05-04):
  1. NEVER overwrite original geo_result.json fields — add refreshed_* fields
  2. Match correct xyz via freq_log stem (strict, no guessing)
  3. First version: only refresh dE, not dG / Shermo
  4. Per-point status tracking — partial failures are visible
  5. Backup geo_result.json before any write

Usage:
  # Step 1: audit first
  python benchmark/dft_theory/stages/refresh_geo_sp.py \
      --session-dir Output/.../session_opt_full_orca_sp1 \
      --mode audit

  # Step 2: refresh after audit passes
  python benchmark/dft_theory/stages/refresh_geo_sp.py \
      --session-dir Output/.../session_opt_full_orca_sp1 \
      --mode refresh \
      --rx-id 1

  # Step 3: re-evaluate
  python benchmark/dft_theory/stages/geo_benchmark.py \
      --session-dir ... --mode evaluate
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

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from benchmark.dft_theory.lib.state_manager import (
    read_benchmark_manifest,
    resolve_method_run_dir,
)

import logging as _logging

_logging.basicConfig(level=_logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

from rph_core.utils.constants import HARTREE_TO_KCAL

from rph_core.utils.optimization_config import prepare_qc_method_config
from rph_core.utils.qc_task_runner import QCTaskRunner

DEFAULTS_CONFIG = PROJECT_ROOT / "config" / "defaults.yaml"
METHODS_GEO_CONFIG = PROJECT_ROOT / "benchmark" / "dft_theory" / "config" / "methods_geo.yaml"

REFRESH_SCHEMA_VERSION = "refresh_v1"

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


def _atom_count(xyz_path: Path) -> int:
    """Return atom count from first line of XYZ file."""
    if not xyz_path.is_file():
        return 0
    first_line = xyz_path.read_text(encoding="utf-8").split("\n", 1)[0].strip()
    try:
        return int(first_line)
    except (ValueError, TypeError):
        return 0


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


# ────────────────────────── find correct xyz ──────────────────────────


class XYZMatchError(Exception):
    """Raised when correct ORCA final xyz cannot be unambiguously identified."""


def find_correct_orca_xyz(
    point_payload: Mapping[str, Any],
    *,
    current_xyz_hash: Optional[str] = None,
) -> Path:
    """
    Locate the correct ORCA final-geometry XYZ by matching the freq_log stem.

    Strategy:
      freq_log = ".../complex_opt_5c52d51c.out"
      → look for ".../complex_opt_5c52d51c.xyz"  (exact match, not _trj)

    Safety:
      - Rejects _trj.xyz files
      - Rejects if matched xyz hash == current (wrong) xyz hash
      - Rejects if multiple ambiguous candidates exist
      - Raises XYZMatchError instead of guessing

    Args:
        point_payload: The point block from geo_result.json (complex/product/ts).
        current_xyz_hash: SHA-256 of the current (buggy) optimized_xyz, used as
                          a guard to prevent re-selecting the same bad file.

    Returns:
        Path to the correct ORCA final-geometry XYZ file.

    Raises:
        XYZMatchError: If the file cannot be found unambiguously.
    """
    freq_log_str = point_payload.get("freq_log")
    if not freq_log_str or not isinstance(freq_log_str, str):
        raise XYZMatchError("No freq_log field in point payload")

    freq_log_path = Path(freq_log_str).expanduser().resolve()
    if not freq_log_path.is_file():
        raise XYZMatchError(f"freq_log file not found: {freq_log_path}")

    # The freq_log has a stem like "complex_opt_5c52d51c"
    freq_stem = freq_log_path.stem
    if not freq_stem:
        raise XYZMatchError(f"Empty freq_log stem: {freq_log_path}")

    # Look for exact sibling xyz
    candidate = freq_log_path.parent / f"{freq_stem}.xyz"
    if candidate.is_file():
        if "_trj" in candidate.name:
            raise XYZMatchError(f"Matched file is a trajectory (should not happen): {candidate}")
        # Guard: matched xyz must NOT be the same as the current wrong one
        matched_hash = _sha256_file(candidate)
        if matched_hash is None:
            raise XYZMatchError(f"Matched xyz file hash is None: {candidate}")
        if current_xyz_hash and matched_hash == current_xyz_hash:
            raise XYZMatchError(
                f"Matched xyz hash ({matched_hash[:16]}...) equals current "
                f"(wrong) xyz hash — would re-select the same bad file"
            )
        return candidate

    # No exact match — check for candidates with same stem prefix (non-trj)
    candidates = [
        f
        for f in freq_log_path.parent.glob(f"{freq_stem}*.xyz")
        if "_trj" not in f.name
    ]
    if len(candidates) == 1:
        matched_hash = _sha256_file(candidates[0])
        if current_xyz_hash and matched_hash == current_xyz_hash:
            raise XYZMatchError("Only candidate has same hash as current wrong xyz")
        return candidates[0]

    if len(candidates) > 1:
        names = ", ".join(c.name for c in candidates)
        raise XYZMatchError(f"Ambiguous candidates ({len(candidates)}): {names}")

    raise XYZMatchError(f"No xyz file found matching freq_log stem '{freq_stem}' in {freq_log_path.parent}")


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
        selected = {str(k): dict(v) for k, v in methods.items() if isinstance(v, dict)}
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
    """
    Resolve SP reader spec.  Priority:
      1. CLI --sp-reader override
      2. sp_reader recorded in geo_result.json
      3. DEFAULT_GEO_SP_READER = "SP-1"
    """
    sp_reader = requested or geo_result_sp_reader or "SP-1"
    sp_methods = benchmark_manifest.get("methods", {}).get("sp", {})
    if not isinstance(sp_methods, dict):
        raise ValueError("Missing methods.sp in benchmark manifest")
    spec = sp_methods.get(sp_reader)
    if not isinstance(spec, dict):
        raise ValueError(f"No SP method spec for '{sp_reader}' in benchmark manifest")
    return sp_reader, dict(spec)


# ────────────────────────── audit ──────────────────────────


def audit_one(
    session_dir: Path,
    rx_id: str,
    method_id: str,
) -> Dict[str, Any]:
    """
    Audit a single (rx, GEO method) combination.

    Returns a dict with per-point audit details:
      - matched_final_xyz path
      - matched_hash vs current_hash comparison
      - atom_count verification
      - refresh_needed flag
    """
    run_dir = resolve_method_run_dir(session_dir, rx_id, "phase2", method_id)
    geo_result_path = run_dir / "geo_result.json"
    if not geo_result_path.exists():
        return {"status": "missing_geo_result", "error": str(geo_result_path)}

    geo_result = _read_json(geo_result_path)
    if str(geo_result.get("status", "")).lower() not in {"completed", "ok"}:
        return {
            "status": "skipped_not_completed",
            "geo_status": geo_result.get("status"),
        }

    audit: Dict[str, Any] = {
        "status": "audited",
        "geo_status": geo_result.get("status"),
        "sp_reader": geo_result.get("sp_reader"),
        "points": {},
    }

    any_refresh_needed = False
    for pt_name in ["complex", "product", "ts"]:
        pt = geo_result.get(pt_name, {})
        if not isinstance(pt, dict):
            audit["points"][pt_name] = {"status": "missing_point_block"}
            continue

        # Current (potentially wrong) xyz
        current_xyz_str = pt.get("optimized_xyz")
        current_hash = None
        current_atoms = 0
        current_path_resolved: Optional[Path] = None
        if current_xyz_str:
            current_path_resolved = Path(current_xyz_str).resolve()
            current_hash = _sha256_file(current_path_resolved)
            current_atoms = _atom_count(current_path_resolved) if current_path_resolved.exists() else 0

        # L2 SP input xyz (what was actually used for SP)
        l2_sp_log_str = pt.get("l2_sp_log")
        l2_sp_input_hash: Optional[str] = None
        if l2_sp_log_str:
            l2_sp_dir = Path(l2_sp_log_str).resolve().parent
            xyz_candidates = list(l2_sp_dir.glob("*.xyz"))
            if xyz_candidates:
                l2_sp_input_hash = _sha256_file(xyz_candidates[0])

        # Find correct ORCA final xyz
        pt_audit: Dict[str, Any] = {
            "current_optimized_xyz": str(current_path_resolved) if current_path_resolved else None,
            "current_hash": current_hash,
            "current_atom_count": current_atoms,
            "l2_sp_input_hash": l2_sp_input_hash,
        }

        try:
            correct_xyz = find_correct_orca_xyz(pt, current_xyz_hash=current_hash)
            correct_hash = _sha256_file(correct_xyz)
            correct_atoms = _atom_count(correct_xyz)

            pt_audit["matched_final_xyz"] = str(correct_xyz)
            pt_audit["matched_hash"] = correct_hash
            pt_audit["matched_atom_count"] = correct_atoms
            pt_audit["match_strategy"] = "freq_log_stem_match"

            # Determine if refresh is needed
            if current_hash and correct_hash and current_hash != correct_hash:
                pt_audit["refresh_needed"] = True
                any_refresh_needed = True
            elif current_hash is None and correct_hash is not None:
                pt_audit["refresh_needed"] = True
                any_refresh_needed = True
            else:
                pt_audit["refresh_needed"] = False
                pt_audit["note"] = "current == matched (no refresh needed, or same wrong file)"

            pt_audit["status"] = "ok"
        except XYZMatchError as exc:
            pt_audit["status"] = "match_failed"
            pt_audit["error"] = str(exc)
            pt_audit["refresh_needed"] = False

        audit["points"][pt_name] = pt_audit

    audit["refresh_needed"] = any_refresh_needed
    return audit


def run_audit(
    session_dir: Path,
    rx_ids: Sequence[str],
    selected_methods: Optional[Sequence[str]],
) -> Dict[str, Dict[str, Any]]:
    """Run audit on all requested (rx, method) combinations and save report."""
    geo_methods, _ = _resolve_geo_methods(METHODS_GEO_CONFIG, selected_methods)
    if not geo_methods:
        raise RuntimeError("No GEO methods resolved for audit")

    logger.info("AUDIT: %d method(s) x %d reaction(s) = %d combination(s)",
                len(geo_methods), len(rx_ids), len(geo_methods) * len(rx_ids))

    report: Dict[str, Dict[str, Any]] = {
        "_meta": {
            "mode": "audit",
            "generated_at": _now_iso(),
            "session_dir": str(session_dir),
            "methods": list(geo_methods.keys()),
            "rx_ids": list(rx_ids),
        },
    }

    n_refreshable = 0
    n_failed = 0
    for rx_id in rx_ids:
        report[rx_id] = {}
        for method_id in geo_methods:
            result = audit_one(session_dir, rx_id, method_id)
            report[rx_id][method_id] = result
            if result.get("refresh_needed"):
                n_refreshable += 1
            if result.get("status") == "match_failed" or any(
                p.get("status") == "match_failed"
                for p in result.get("points", {}).values()
                if isinstance(p, dict)
            ):
                n_failed += 1

            # Log summary per combination
            pts = result.get("points", {})
            pt_status = {k: v.get("status", "?") for k, v in pts.items() if isinstance(v, dict)}
            logger.info("  rx%s %-6s → %s  points: %s",
                        rx_id, method_id, result.get("status", "?"), pt_status)

    report["_meta"]["refreshable_combinations"] = n_refreshable
    report["_meta"]["failed_matches"] = n_failed

    # Save report
    report_path = session_dir / "reports" / "geo_sp_refresh_audit.json"
    _write_json(report_path, report)
    logger.info("Audit report saved: %s", report_path)
    logger.info("  refreshable: %d  failed_matches: %d", n_refreshable, n_failed)

    # Also render markdown
    md_path = session_dir / "reports" / "geo_sp_refresh_audit.md"
    _render_audit_md(report, md_path)

    return report


def _render_audit_md(report: Mapping[str, Any], path: Path) -> None:
    lines = [
        "# GEO SP Refresh Audit Report",
        "",
        f"- generated: {report.get('_meta', {}).get('generated_at')}",
        f"- refreshable: {report.get('_meta', {}).get('refreshable_combinations')}",
        f"- failed_matches: {report.get('_meta', {}).get('failed_matches')}",
        "",
        "| rx | method | point | current_hash | matched_hash | refresh_needed | status |",
        "|---|---|---|---|---|---|---|",
    ]
    for rx_id, methods in report.items():
        if rx_id.startswith("_"):
            continue
        if not isinstance(methods, dict):
            continue
        for method_id, result in methods.items():
            if not isinstance(result, dict):
                continue
            pts = result.get("points", {})
            for pt_name, pt_data in pts.items():
                if not isinstance(pt_data, dict):
                    continue
                cur_h = (pt_data.get("current_hash") or "N/A")[:12]
                mat_h = (pt_data.get("matched_hash") or "N/A")[:12]
                needed = "YES" if pt_data.get("refresh_needed") else "no"
                status = pt_data.get("status", "?")
                lines.append(f"| rx{rx_id} | {method_id} | {pt_name} | {cur_h}... | {mat_h}... | {needed} | {status} |")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ────────────────────────── refresh ──────────────────────────


def refresh_one(
    session_dir: Path,
    rx_id: str,
    method_id: str,
    geo_method_spec: Mapping[str, Any],
    defaults_payload: Mapping[str, Any],
    sp_reader_spec: Mapping[str, Any],
    solvent_cfg: Mapping[str, Any],
    force: bool = False,
) -> Dict[str, Any]:
    """
    Refresh L2 SP for one (rx, GEO method) combination.

    Steps:
      1. Run audit to find correct xyz files
      2. For each point where refresh_needed=True:
         a. Re-run QCTaskRunner.run_sp_only() on the correct xyz
      3. Write incremental ``refreshed_*`` fields to geo_result.json
      4. Backup original geo_result.json first
    """
    run_dir = resolve_method_run_dir(session_dir, rx_id, "phase2", method_id)
    geo_result_path = run_dir / "geo_result.json"
    if not geo_result_path.exists():
        return {"status": "missing_geo_result"}

    geo_result = _read_json(geo_result_path)

    # Skip if already refreshed (unless force)
    existing_refresh = geo_result.get("refresh_sp")
    if existing_refresh and not force:
        if str(existing_refresh.get("schema_version")) == REFRESH_SCHEMA_VERSION:
            logger.info("  rx%s %s already refreshed at %s, skip (--force to override)",
                        rx_id, method_id, existing_refresh.get("refreshed_at", "?"))
            return {"status": "already_refreshed"}

    # Step 1: Audit to find correct xyz files
    audit_result = audit_one(session_dir, rx_id, method_id)
    if audit_result.get("status") != "audited":
        return {"status": "audit_failed", "audit": audit_result}

    # Check if any point needs refresh
    points_audit = audit_result.get("points", {})
    refreshable_pts = {
        k: v for k, v in points_audit.items()
        if isinstance(v, dict) and v.get("refresh_needed")
    }
    if not refreshable_pts:
        logger.info("  rx%s %s no points need refresh", rx_id, method_id)
        return {"status": "no_refresh_needed"}

    # Step 2: Build QCTaskRunner
    config = prepare_qc_method_config(
        defaults_payload,
        optimization=geo_method_spec,
        single_point=sp_reader_spec,
        solvent=solvent_cfg,
    )
    runner = QCTaskRunner(config)

    # Step 3: Backup
    backup_path = geo_result_path.with_suffix(".json.bak_pre_refresh")
    if not backup_path.exists():
        shutil.copy2(geo_result_path, backup_path)
        logger.info("  backup → %s", backup_path.name)

    # Step 4: Re-run SP per point
    refresh_sp: Dict[str, Any] = {
        "schema_version": REFRESH_SCHEMA_VERSION,
        "refreshed_at": _now_iso(),
        "sp_reader": audit_result.get("sp_reader"),
        "points": {},
    }

    for pt_name in ["complex", "product", "ts"]:
        pt_audit = points_audit.get(pt_name, {})
        if not isinstance(pt_audit, dict) or not pt_audit.get("refresh_needed"):
            refresh_sp["points"][pt_name] = {"status": "skipped"}
            continue

        matched_xyz_str = pt_audit.get("matched_final_xyz")
        if not matched_xyz_str:
            refresh_sp["points"][pt_name] = {
                "status": "failed",
                "error": "no_matched_xyz_in_audit",
            }
            continue

        matched_xyz = Path(matched_xyz_str).resolve()
        if not matched_xyz.exists():
            refresh_sp["points"][pt_name] = {
                "status": "failed",
                "error": f"matched_xyz_not_found: {matched_xyz}",
            }
            continue

        # Verify atom count
        expected_atoms = pt_audit.get("current_atom_count", 0)
        actual_atoms = _atom_count(matched_xyz)
        if expected_atoms > 0 and actual_atoms != expected_atoms:
            refresh_sp["points"][pt_name] = {
                "status": "failed",
                "error": f"atom_count_mismatch: expected={expected_atoms} got={actual_atoms}",
            }
            continue

        # Run SP
        sp_output_dir = run_dir / "opt" / pt_name / f"L2_SP_refreshed"
        try:
            started = time.time()
            sp_result = runner.run_sp_only(matched_xyz, sp_output_dir, charge=0, spin=1)
            sp_runtime = time.time() - started
        except Exception as exc:
            logger.error("  rx%s %s %s SP failed: %s", rx_id, method_id, pt_name, exc)
            refresh_sp["points"][pt_name] = {
                "status": "failed",
                "error": str(exc),
            }
            continue

        if not sp_result.converged:
            refresh_sp["points"][pt_name] = {
                "status": "failed",
                "error": f"sp_not_converged: {sp_result.error_message}",
            }
            continue

        new_energy = _safe_float(sp_result.energy)
        refresh_sp["points"][pt_name] = {
            "status": "completed",
            "refreshed_l2_energy_hartree": new_energy,
            "refreshed_l2_sp_log": str(sp_result.output_file.resolve()) if sp_result.output_file else None,
            "refreshed_optimized_xyz": str(matched_xyz),
            "refreshed_optimized_xyz_sha256": _sha256_file(matched_xyz),
            "source_freq_log": points_audit[pt_name].get("freq_log"),
            "refresh_sp_runtime_seconds": sp_runtime,
        }
        logger.info("  rx%s %s %s SP → %.8f Ha (%.1fs)",
                    rx_id, method_id, pt_name, new_energy or 0.0, sp_runtime)

    # Step 5: Compute refreshed dE if all 3 points completed
    refreshed_energies = {}
    for pt_name in ["complex", "product", "ts"]:
        pt_refresh = refresh_sp["points"].get(pt_name, {})
        if isinstance(pt_refresh, dict) and pt_refresh.get("status") == "completed":
            refreshed_energies[pt_name] = pt_refresh.get("refreshed_l2_energy_hartree")

    all_three = all(p in refreshed_energies for p in ["complex", "product", "ts"])
    if all_three:
        refresh_sp["refreshed_dE_activation_kcal"] = _metric_kcal(
            refreshed_energies["ts"], refreshed_energies["complex"]
        )
        refresh_sp["refreshed_dE_reaction_kcal"] = _metric_kcal(
            refreshed_energies["product"], refreshed_energies["complex"]
        )
        logger.info("  rx%s %s refreshed dE‡=%.3f  dE_rxn=%.3f kcal",
                    rx_id, method_id,
                    refresh_sp["refreshed_dE_activation_kcal"] or 0.0,
                    refresh_sp["refreshed_dE_reaction_kcal"] or 0.0)
    else:
        missing = [p for p in ["complex", "product", "ts"] if p not in refreshed_energies]
        refresh_sp["refreshed_dE_activation_kcal"] = None
        refresh_sp["refreshed_dE_reaction_kcal"] = None
        refresh_sp["incomplete_points"] = missing
        logger.warning("  rx%s %s incomplete refresh (missing: %s), dE not updated",
                       rx_id, method_id, missing)

    # Step 6: Determine overall refresh status
    completed_pts = sum(
        1 for p in refresh_sp["points"].values()
        if isinstance(p, dict) and p.get("status") == "completed"
    )
    if completed_pts == 3:
        refresh_sp["status"] = "completed"
    elif completed_pts > 0:
        refresh_sp["status"] = "partial"
    else:
        refresh_sp["status"] = "failed"

    # Step 7: Write incremental fields to geo_result.json (DO NOT overwrite originals)
    geo_result["refresh_sp"] = refresh_sp

    # P1: Compute canonical effective results for downstream report generators.
    # When refresh is complete, the effective geometry/energy should prefer the
    # refreshed values; otherwise fall back to the original (potentially stale) values.
    effective: Dict[str, Any] = {"source": "original"}
    if refresh_sp.get("status") == "completed":
        effective["source"] = "refreshed"
        for pt_name in ["complex", "product", "ts"]:
            pt_ref = refresh_sp["points"].get(pt_name, {})
            if isinstance(pt_ref, dict):
                effective[f"effective_{pt_name}_l2_energy_hartree"] = pt_ref.get("refreshed_l2_energy_hartree")
                effective[f"effective_{pt_name}_optimized_xyz"] = pt_ref.get("refreshed_optimized_xyz")
                effective[f"effective_{pt_name}_optimized_xyz_sha256"] = pt_ref.get("refreshed_optimized_xyz_sha256")
        effective["effective_dE_activation_kcal"] = refresh_sp.get("refreshed_dE_activation_kcal")
        effective["effective_dE_reaction_kcal"] = refresh_sp.get("refreshed_dE_reaction_kcal")
        effective["effective_sp_reader"] = refresh_sp.get("sp_reader")
    else:
        effective["source"] = "original_or_partial"
        orig_ts = geo_result.get("ts", {})
        if isinstance(orig_ts, dict):
            effective["effective_ts_optimized_xyz"] = orig_ts.get("optimized_xyz")
            effective["effective_ts_optimized_xyz_sha256"] = orig_ts.get("optimized_xyz_sha256")
        effective["effective_dE_activation_kcal"] = geo_result.get("dE_activation_kcal")
        effective["effective_dE_reaction_kcal"] = geo_result.get("dE_reaction_kcal")
    geo_result["effective"] = effective

    _write_json(geo_result_path, geo_result)

    return {"status": refresh_sp["status"], "refresh_sp": refresh_sp}


def run_refresh(
    session_dir: Path,
    rx_ids: Sequence[str],
    selected_methods: Optional[Sequence[str]],
    sp_reader_override: Optional[str],
    force: bool,
) -> Dict[str, Dict[str, Any]]:
    """Run refresh on all requested (rx, method) combinations."""
    benchmark_manifest = read_benchmark_manifest(session_dir)
    defaults_payload = _load_yaml(DEFAULTS_CONFIG)
    geo_methods, solvent_cfg = _resolve_geo_methods(METHODS_GEO_CONFIG, selected_methods)

    if not geo_methods:
        raise RuntimeError("No GEO methods resolved for refresh")

    logger.info("REFRESH: %d method(s) x %d reaction(s) = %d task(s)",
                len(geo_methods), len(rx_ids), len(geo_methods) * len(rx_ids))

    results: Dict[str, Dict[str, Any]] = {}
    total = len(geo_methods) * len(rx_ids)
    done = 0
    sp_reader_id = sp_reader_override or "SP-1"
    for rx_id in rx_ids:
        results[rx_id] = {}
        # Resolve SP reader from geo_result or CLI
        first_method = next(iter(geo_methods))
        run_dir = resolve_method_run_dir(session_dir, rx_id, "phase2", first_method)
        geo_result_path = run_dir / "geo_result.json"
        geo_sp_reader = None
        if geo_result_path.exists():
            geo_sp_reader = _read_json(geo_result_path).get("sp_reader")

        actual_sp_reader_id, sp_reader_spec = _resolve_sp_reader_spec(
            benchmark_manifest, sp_reader_override, geo_sp_reader
        )
        sp_reader_id = actual_sp_reader_id
        logger.info("  SP reader: %s", sp_reader_id)

        for method_id, geo_spec in geo_methods.items():
            done += 1
            logger.info("[%d/%d] rx%s %s ...", done, total, rx_id, method_id)
            result = refresh_one(
                session_dir,
                rx_id,
                method_id,
                geo_spec,
                defaults_payload,
                sp_reader_spec,
                solvent_cfg,
                force=force,
            )
            results[rx_id][method_id] = result

    # Save refresh summary
    summary_path = session_dir / "reports" / "geo_sp_refresh_summary.json"
    summary: Dict[str, Any] = {
        "generated_at": _now_iso(),
        "mode": "refresh",
        "sp_reader": sp_reader_id,
        "results": results,
    }
    _write_json(summary_path, summary)
    logger.info("Refresh summary: %s", summary_path)

    return results


# ────────────────────────── CLI ──────────────────────────


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="GEO Benchmark SP Refresh: re-run L2 SP on correct ORCA geometries."
    )
    parser.add_argument("--session-dir", required=True, help="Benchmark session directory")
    parser.add_argument(
        "--mode",
        choices=["audit", "refresh"],
        default="audit",
        help="audit: discover + validate (no QC) | refresh: re-run SP + update results",
    )
    parser.add_argument("--rx-id", nargs="*", default=None, help="Subset of rx ids")
    parser.add_argument("--methods", nargs="*", default=None, help="Subset of GEO method ids")
    parser.add_argument("--sp-reader", default=None, help="Force SP reader (default: from geo_result.json or SP-1)")
    parser.add_argument("--force", action="store_true", help="Re-refresh even if already done")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    session_dir = Path(args.session_dir).expanduser().resolve()
    if not session_dir.exists():
        raise FileNotFoundError(f"Session dir not found: {session_dir}")

    rx_ids = _resolve_rx_ids(session_dir, args.rx_id)
    if not rx_ids:
        raise RuntimeError("No rx ids resolved")

    if args.mode == "audit":
        run_audit(session_dir, rx_ids, args.methods)
    elif args.mode == "refresh":
        run_refresh(session_dir, rx_ids, args.methods, args.sp_reader, args.force)
    else:
        raise ValueError(f"Unknown mode: {args.mode}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
