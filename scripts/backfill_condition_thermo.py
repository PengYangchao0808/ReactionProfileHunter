#!/usr/bin/env python3
"""
RPH Condition Thermo Backfill — Shermo Centralization Migration
=================================================================

Walks rph_output_backup and regenerates all Shermo thermodynamic outputs
at the correct condition temperature.  Does NOT run Gaussian/ORCA/xTB —
only Shermo (~2-5 seconds per species).

Skip / idempotency strategy (v2):
  - No script-level CSV-existence skip.  ConditionThermoEngine.run()
    has its own comprehensive idempotency gate (temperature validation,
    delta-key completeness, energy-data presence, blocking-warning check).
  - ``--force`` pre-deletes thermo_calculation.json, condition_features_mlr.csv,
    merged_features.csv, merged_features.json so the engine re-runs from scratch.

Usage:
    # Phase 0 — audit only
    python scripts/backfill_condition_thermo.py \
        --rph-root external_data/rph_output_backup \
        --config config/defaults.yaml \
        --dry-run

    # Phase 1 — single reaction × single condition × single branch
    python scripts/backfill_condition_thermo.py \
        --rph-root external_data/rph_output_backup \
        --config config/defaults.yaml \
        --reaction RXN_0b71b9e9 --condition COND_8 --branch BR_MAJOR

    # Phase 2 — single reaction × all conditions
    python scripts/backfill_condition_thermo.py \
        --rph-root external_data/rph_output_backup \
        --config config/defaults.yaml \
        --reaction RXN_0b71b9e9

    # Phase 3 — full backfill with JSON summary
    python scripts/backfill_condition_thermo.py \
        --rph-root external_data/rph_output_backup \
        --config config/defaults.yaml \
        --output-summary backfill_summary.json

    # Force re-run everything (ignore existing outputs)
    python scripts/backfill_condition_thermo.py \
        --rph-root external_data/rph_output_backup \
        --config config/defaults.yaml \
        --force
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

# Ensure rph_core is importable without pip install -e
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
)
logger = logging.getLogger("backfill_thermo")


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------
@dataclass
class ConditionBranchTarget:
    reaction_id: str
    condition_id: str
    branch_id: str
    reaction_root: Path
    condition_root: Path
    branch_root: Path


def discover_targets(
    rph_root: Path,
    reaction_filter: Optional[str] = None,
    condition_filter: Optional[str] = None,
    branch_filter: Optional[str] = None,
) -> List[ConditionBranchTarget]:
    """Walk rph_root and discover all condition × branch targets.

    Directory layout expected::

        rph_root/
          RXN_*/
            conditions/
              COND_*/
                condition_manifest.json   (required)
                branches/
                  BR_*/
            branches/
              BR_*/                       (fallback if conditions/*/branches/ empty)

    """
    targets: List[ConditionBranchTarget] = []
    rph_root = Path(rph_root)

    for reaction_dir in sorted(rph_root.iterdir()):
        if not reaction_dir.is_dir() or not reaction_dir.name.startswith("RXN_"):
            continue
        reaction_id = reaction_dir.name
        if reaction_filter and reaction_id != reaction_filter:
            continue

        conditions_dir = reaction_dir / "conditions"
        if not conditions_dir.is_dir():
            logger.debug("%s: no conditions/ directory", reaction_id)
            continue

        for cond_dir in sorted(conditions_dir.iterdir()):
            if not cond_dir.is_dir() or not cond_dir.name.startswith("COND_"):
                continue
            condition_id = cond_dir.name
            if condition_filter and condition_id != condition_filter:
                continue

            manifest_path = cond_dir / "condition_manifest.json"
            if not manifest_path.exists():
                logger.warning(
                    "%s/%s: missing condition_manifest.json", reaction_id, condition_id,
                )
                continue

            # Discover branch IDs from condition-level branches/ first
            cond_branches_dir = cond_dir / "branches"
            branch_ids: List[str] = []
            if cond_branches_dir.is_dir():
                branch_ids = sorted(
                    b.name for b in cond_branches_dir.iterdir() if b.is_dir()
                )

            # Fallback: reaction-level branches/
            if not branch_ids:
                reaction_branches_dir = reaction_dir / "branches"
                if reaction_branches_dir.is_dir():
                    branch_ids = sorted(
                        b.name for b in reaction_branches_dir.iterdir() if b.is_dir()
                    )

            for branch_id in branch_ids:
                if branch_filter and branch_id != branch_filter:
                    continue
                branch_root = reaction_dir / "branches" / branch_id
                targets.append(ConditionBranchTarget(
                    reaction_id=reaction_id,
                    condition_id=condition_id,
                    branch_id=branch_id,
                    reaction_root=reaction_dir,
                    condition_root=cond_dir,
                    branch_root=branch_root,
                ))

    logger.info("Discovered %d condition×branch targets", len(targets))
    return targets


# ---------------------------------------------------------------------------
# Failure classification
# ---------------------------------------------------------------------------
class FailureCategory:
    """Classify a backfill result into A/B/C/D tiers.

    A — Complete success: dE_activation + dG_activation + dE_reaction + dG_reaction
    B — Electronic barrier only: dE_activation present, dG_activation missing
    C — S3 incomplete: missing sp_matrix_metadata.json or SP energies
    D — Anomaly: script/engine error, no output files
    """
    A = "A_complete"
    B = "B_electronic_barrier_only"
    C = "C_s3_incomplete"
    D = "D_anomaly"


def classify_result(result: Dict[str, Any]) -> str:
    """Classify a completed backfill result into a tier."""
    if result["status"] == "failed":
        return FailureCategory.D
    if result["status"] != "ok":
        return FailureCategory.D

    has_de = result.get("has_dE_activation", False)
    has_dg = result.get("has_dG_activation", False)
    has_de_rxn = result.get("has_dE_reaction", False)
    has_dg_rxn = result.get("has_dG_reaction", False)
    has_blocking = result.get("has_blocking_warning", False)

    if has_blocking:
        return FailureCategory.C
    if has_de and has_dg and has_de_rxn and has_dg_rxn:
        return FailureCategory.A
    if has_de:
        return FailureCategory.B
    # No dE at all but status was "ok" — something is off
    return FailureCategory.D


def _check_blocking_warning(warnings: List[str]) -> bool:
    """Return True if any warning is a blocking warning."""
    blocking_prefixes = ("W_MISSING_SP_MATRIX_METADATA", "W_MISSING_SP_ENERGIES")
    return any(
        any(w.startswith(prefix) for prefix in blocking_prefixes)
        for w in warnings
    )


# ---------------------------------------------------------------------------
# Backfill
# ---------------------------------------------------------------------------
def _force_clean_outputs(branch_output_dir: Path) -> None:
    """Delete existing output files to force a full re-run."""
    files_to_remove = [
        branch_output_dir / "thermo_calculation.json",
        branch_output_dir / "condition_features_mlr.csv",
        branch_output_dir / "merged_features.csv",
        branch_output_dir / "merged_features.json",
    ]
    for f in files_to_remove:
        if f.exists():
            try:
                f.unlink()
                logger.debug("Removed %s (force)", f)
            except OSError as exc:
                logger.warning("Failed to remove %s: %s", f, exc)


def backfill_target(
    target: ConditionBranchTarget,
    config: Dict[str, Any],
    dry_run: bool = False,
    force: bool = False,
) -> Dict[str, Any]:
    """Run ConditionThermoEngine + ConditionFeatureMerger for one target.

    Skip logic:
      - With ``--force``: pre-delete output files, then always run engine.
      - Without ``--force``: delegate to engine.run() internal idempotency.
      - ``--dry-run``: report what would happen without running engine.

    """
    from rph_core.steps.condition_thermo import ConditionThermoEngine
    from rph_core.steps.condition_feature_merger import ConditionFeatureMerger

    label = f"{target.reaction_id}/{target.condition_id}/{target.branch_id}"
    result: Dict[str, Any] = {
        "target": label,
        "status": "pending",
        "temperature_K": None,
        "errors": [],
        "warnings": [],
        "runtime_s": 0.0,
        "has_dE_activation": False,
        "has_dG_activation": False,
        "has_dE_reaction": False,
        "has_dG_reaction": False,
        "has_blocking_warning": False,
        "failure_category": None,
        "existing_state": None,
    }

    # --- Parse manifest for temperature ---
    manifest_path = target.condition_root / "condition_manifest.json"
    try:
        manifest_data = json.loads(manifest_path.read_text())
        result["temperature_K"] = manifest_data.get("temperature_K")
    except Exception:
        pass

    # --- Output directory ---
    branch_output_dir = target.condition_root / "branches" / target.branch_id
    thermo_json_path = branch_output_dir / "thermo_calculation.json"
    output_csv = branch_output_dir / "condition_features_mlr.csv"

    # --- Pre-flight: check existing state for dry-run reporting ---
    if thermo_json_path.exists():
        try:
            existing = json.loads(thermo_json_path.read_text())
            warnings_list = [str(w) for w in (existing.get("warnings") or [])]
            delta = existing.get("delta") or {}
            result["existing_state"] = {
                "temperature_K": existing.get("temperature_K"),
                "has_dE_activation": delta.get("dE_activation") is not None,
                "has_dG_activation": delta.get("dG_activation") is not None,
                "has_blocking_warning": _check_blocking_warning(warnings_list),
            }
        except Exception:
            result["existing_state"] = {"parse_error": True}

    # --- Dry-run: report and return ---
    if dry_run:
        # Determine what would happen
        if force:
            action = "force_rerun"
        elif result["existing_state"] is not None:
            action = "engine_idempotency_check"
        else:
            action = "fresh_run"
        result["status"] = "dry_run"
        result["proposed_action"] = action
        logger.info(
            "[DRY-RUN] %s (T=%.1f K, action=%s)",
            label, result["temperature_K"] or 298.15, action,
        )
        return result

    # --- Force: clean existing outputs ---
    if force:
        _force_clean_outputs(branch_output_dir)

    # --- Run engine (engine handles its own idempotency) ---
    start_time = time.time()
    try:
        engine = ConditionThermoEngine(
            config=config,
            reaction_root=target.reaction_root,
            condition_root=branch_output_dir,
            artifact_root=target.branch_root,
            branch_id=target.branch_id,
        )
        thermo_path = engine.run()

        # Read back thermo results for reporting
        if thermo_path and thermo_path.exists():
            try:
                thermo_data = json.loads(thermo_path.read_text())
                thermo_warnings = [str(w) for w in (thermo_data.get("warnings") or [])]
                result["warnings"].extend(thermo_warnings)
                delta = thermo_data.get("delta") or {}

                de_act = delta.get("dE_activation")
                dg_act = delta.get("dG_activation")
                de_rxn = delta.get("dE_reaction")
                dg_rxn = delta.get("dG_reaction")

                result["has_dE_activation"] = (
                    de_act is not None
                    and not (isinstance(de_act, float) and math.isnan(de_act))
                )
                result["has_dG_activation"] = (
                    dg_act is not None
                    and not (isinstance(dg_act, float) and math.isnan(dg_act))
                )
                result["has_dE_reaction"] = (
                    de_rxn is not None
                    and not (isinstance(de_rxn, float) and math.isnan(de_rxn))
                )
                result["has_dG_reaction"] = (
                    dg_rxn is not None
                    and not (isinstance(dg_rxn, float) and math.isnan(dg_rxn))
                )
                result["has_blocking_warning"] = _check_blocking_warning(thermo_warnings)
            except Exception as exc:
                result["warnings"].append(f"W_READBACK_FAILED:{exc}")

        # --- Run merger ---
        merger = ConditionFeatureMerger(
            reaction_root=target.reaction_root,
            condition_root=branch_output_dir,
            feature_root=(
                target.branch_root
                if (target.branch_root / "reaction_features").is_dir()
                else None
            ),
            branch_id=target.branch_id,
        )
        merger.run()

        result["status"] = "ok"

    except Exception as exc:
        result["status"] = "failed"
        result["errors"].append(str(exc))
        logger.error("[FAIL] %s: %s", label, exc)

    result["runtime_s"] = round(time.time() - start_time, 1)

    # --- Classify failure tier ---
    result["failure_category"] = classify_result(result)

    if result["status"] == "ok":
        logger.info(
            "[OK] %s (%.1fs, T=%.1f K, tier=%s, dE=%s, dG=%s)",
            label,
            result["runtime_s"],
            result["temperature_K"] or 298.15,
            result["failure_category"],
            "Y" if result["has_dE_activation"] else "N",
            "Y" if result["has_dG_activation"] else "N",
        )
    return result


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
def build_summary(
    results: List[Dict[str, Any]],
    rph_root: Path,
    config_path: Path,
    dry_run: bool,
    force: bool,
    elapsed_total: float,
) -> Dict[str, Any]:
    """Build a structured summary from all target results."""
    total = len(results)
    by_status: Dict[str, int] = {}
    for r in results:
        s = r["status"]
        by_status[s] = by_status.get(s, 0) + 1

    # Tier counts (only for non-dry-run)
    tier_counts: Dict[str, int] = {}
    for r in results:
        tier = r.get("failure_category")
        if tier is not None:
            tier_counts[tier] = tier_counts.get(tier, 0) + 1

    # Coverage counts
    with_de = sum(1 for r in results if r.get("has_dE_activation"))
    with_dg = sum(1 for r in results if r.get("has_dG_activation"))
    with_de_rxn = sum(1 for r in results if r.get("has_dE_reaction"))
    with_dg_rxn = sum(1 for r in results if r.get("has_dG_reaction"))
    with_blocking = sum(1 for r in results if r.get("has_blocking_warning"))

    # Warning breakdown
    warning_counts: Dict[str, int] = {}
    for r in results:
        for w in r.get("warnings", []):
            # Normalize: keep prefix before first colon/separator
            prefix = w.split(":")[0] if ":" in w else w
            warning_counts[prefix] = warning_counts.get(prefix, 0) + 1

    # Runtime stats
    runtimes = [r["runtime_s"] for r in results if r.get("runtime_s", 0) > 0]
    avg_runtime = sum(runtimes) / len(runtimes) if runtimes else 0.0

    summary: Dict[str, Any] = {
        "rph_root": str(rph_root),
        "config_path": str(config_path),
        "dry_run": dry_run,
        "force": force,
        "elapsed_total_s": round(elapsed_total, 1),
        "targets_total": total,
        "by_status": by_status,
        "tier_counts": tier_counts,
        "coverage": {
            "with_dE_activation": with_de,
            "with_dG_activation": with_dg,
            "with_dE_reaction": with_de_rxn,
            "with_dG_reaction": with_dg_rxn,
            "with_blocking_warning": with_blocking,
            "dE_rate": round(with_de / total, 3) if total > 0 else 0.0,
            "dG_rate": round(with_dg / total, 3) if total > 0 else 0.0,
        },
        "warning_counts": dict(sorted(warning_counts.items(), key=lambda x: -x[1])),
        "runtime_avg_s": round(avg_runtime, 2),
    }
    return summary


def log_summary(summary: Dict[str, Any]) -> None:
    """Log a human-readable summary to the console."""
    logger.info("\n" + "=" * 60)
    logger.info("  BACKFILL SUMMARY")
    logger.info("=" * 60)
    logger.info("RPH root:     %s", summary["rph_root"])
    logger.info("Config:       %s", summary["config_path"])
    logger.info("Dry run:      %s", summary["dry_run"])
    logger.info("Force:        %s", summary["force"])
    logger.info("Elapsed:      %.1fs", summary["elapsed_total_s"])
    logger.info("")
    logger.info("Targets:      %d", summary["targets_total"])

    by_status = summary.get("by_status", {})
    for status, count in sorted(by_status.items()):
        logger.info("  %-12s %d", status + ":", count)

    tiers = summary.get("tier_counts", {})
    if tiers:
        logger.info("")
        logger.info("Failure classification:")
        for tier, count in sorted(tiers.items()):
            label_map = {
                "A_complete": "A — Complete (dE+dG+dE_rxn+dG_rxn)",
                "B_electronic_barrier_only": "B — Electronic barrier only (dE, no dG)",
                "C_s3_incomplete": "C — S3 incomplete (blocking warning)",
                "D_anomaly": "D — Anomaly (no output / error)",
            }
            logger.info("  %s: %d", label_map.get(tier, tier), count)

    cov = summary.get("coverage", {})
    if cov:
        logger.info("")
        logger.info("Coverage:")
        logger.info("  dE_activation: %d/%d (%.1f%%)",
                     cov.get("with_dE_activation", 0), summary["targets_total"],
                     cov.get("dE_rate", 0) * 100)
        logger.info("  dG_activation: %d/%d (%.1f%%)",
                     cov.get("with_dG_activation", 0), summary["targets_total"],
                     cov.get("dG_rate", 0) * 100)
        logger.info("  dE_reaction:   %d/%d",
                     cov.get("with_dE_reaction", 0), summary["targets_total"])
        logger.info("  dG_reaction:   %d/%d",
                     cov.get("with_dG_reaction", 0), summary["targets_total"])
        logger.info("  Blocking warn: %d",
                     cov.get("with_blocking_warning", 0))

    warn_counts = summary.get("warning_counts", {})
    if warn_counts:
        logger.info("")
        logger.info("Warning breakdown:")
        for prefix, count in warn_counts.items():
            logger.info("  %-45s %d", prefix, count)

    logger.info("=" * 60)


def write_summary_json(summary: Dict[str, Any], output_path: Path) -> None:
    """Write summary to a JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    logger.info("Summary written to %s", output_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(
        description="RPH Condition Thermo Backfill",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--rph-root", type=Path, required=True,
                        help="Root directory containing RXN_* outputs")
    parser.add_argument("--config", type=Path, default=Path("config/defaults.yaml"),
                        help="Path to defaults.yaml config")
    parser.add_argument("--dry-run", action="store_true",
                        help="Audit only: report targets without running engine")
    parser.add_argument("--force", action="store_true",
                        help="Delete existing outputs before running engine")
    parser.add_argument("--reaction", type=str, default=None,
                        help="Filter to a specific reaction ID (e.g. RXN_0b71b9e9)")
    parser.add_argument("--condition", type=str, default=None,
                        help="Filter to a specific condition ID (e.g. COND_8)")
    parser.add_argument("--branch", type=str, default=None,
                        help="Filter to a specific branch ID (e.g. BR_MAJOR)")
    parser.add_argument("--output-summary", type=Path, default=None,
                        help="Write JSON summary to this path")
    args = parser.parse_args()

    rph_root = Path(args.rph_root)
    config_path = Path(args.config)

    if not rph_root.exists():
        logger.error("RPH root not found: %s", rph_root)
        return 1
    if not config_path.exists():
        logger.error("Config not found: %s", config_path)
        return 1

    try:
        from rph_core.utils.config_loader import load_config
        config = load_config(config_path)
    except Exception as exc:
        logger.error("Failed to load config: %s", exc)
        return 1

    logger.info("=== RPH Condition Thermo Backfill ===")
    logger.info("RPH root: %s", rph_root)
    logger.info("Config:   %s", config_path)
    logger.info("Dry run:  %s", args.dry_run)
    logger.info("Force:    %s", args.force)
    if args.reaction:
        logger.info("Reaction: %s", args.reaction)
    if args.condition:
        logger.info("Condition: %s", args.condition)
    if args.branch:
        logger.info("Branch:   %s", args.branch)

    targets = discover_targets(
        rph_root,
        reaction_filter=args.reaction,
        condition_filter=args.condition,
        branch_filter=args.branch,
    )
    if not targets:
        logger.warning("No targets found.")
        return 1

    logger.info("\nProcessing %d targets...\n", len(targets))

    total_start = time.time()
    results: List[Dict[str, Any]] = []

    for i, target in enumerate(targets):
        label = f"[{i + 1}/{len(targets)}]"
        logger.info(
            "%s %s/%s/%s", label,
            target.reaction_id, target.condition_id, target.branch_id,
        )
        result = backfill_target(
            target, config,
            dry_run=args.dry_run,
            force=args.force,
        )
        results.append(result)

    elapsed_total = time.time() - total_start

    # --- Build and display summary ---
    summary = build_summary(
        results, rph_root, config_path,
        dry_run=args.dry_run,
        force=args.force,
        elapsed_total=elapsed_total,
    )
    log_summary(summary)

    if args.output_summary:
        write_summary_json(summary, args.output_summary)

    # --- Report failed targets ---
    failed = [r for r in results if r["status"] == "failed"]
    if failed:
        logger.warning("\nFailed targets (%d):", len(failed))
        for r in failed:
            err = r["errors"][0] if r["errors"] else "unknown"
            logger.warning("  %s: %s", r["target"], err)

    # --- Report D-tier anomalies ---
    anomalies = [
        r for r in results
        if r.get("failure_category") == FailureCategory.D and r["status"] != "failed"
    ]
    if anomalies:
        logger.warning("\nD-tier anomalies (ok but suspicious) (%d):", len(anomalies))
        for r in anomalies:
            logger.warning(
                "  %s: dE=%s dG=%s warnings=%s",
                r["target"],
                "Y" if r.get("has_dE_activation") else "N",
                "Y" if r.get("has_dG_activation") else "N",
                len(r.get("warnings", [])),
            )

    # --- Exit code: 0 if no failures, 1 otherwise ---
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
