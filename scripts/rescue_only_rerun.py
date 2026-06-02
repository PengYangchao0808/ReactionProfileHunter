#!/usr/bin/env python3
"""
S3 TS Rescue-Only Rerun — skip primary Berny, jump directly to rescue.

Automatically infers failure kind and rescue starting geometry from the
pre-existing failed Gaussian log, then invokes QCTaskRunner.run_ts_rescue_only().

Usage:  python scripts/rescue_only_rerun.py <branch_work_dir> [<branch_work_dir> ...]
        python scripts/rescue_only_rerun.py rph_output/RXN_0b71b9e9/branches/BR_DR_001
        python scripts/rescue_only_rerun.py --clean rph_output/RXN_*/branches/BR_DR_001
"""

import sys
import argparse
import logging
import re
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from rph_core.utils.qc_task_runner import (
    QCTaskRunner,
    TS_FAILURE_OSCILLATION,
    TS_FAILURE_NO_IMAG,
    TS_FAILURE_MULTI_IMAG,
    TS_FAILURE_FATAL,
    TS_FAILURE_UNKNOWN,
)

logger = logging.getLogger("rescue_only_rerun")


def _load_config(config_path: Path):
    import yaml
    with open(config_path) as fh:
        return yaml.safe_load(fh)


def _find_failed_log(branch_dir: Path):
    rx_name = branch_dir.parent.parent.name
    br_name = branch_dir.name
    candidates = [
        REPO_ROOT / "rph_output" / "_failed_s3_backups" / f"{rx_name}_{br_name}" / "ts_guess.log",
        branch_dir / "S3_TS" / "ts_opt" / "berny" / "ts_guess.log",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _classify_from_log(log_path: Path):
    try:
        text = log_path.read_text(errors="ignore")
    except OSError:
        return TS_FAILURE_UNKNOWN

    if "Error termination" in text or "Number of steps exceeded" in text:
        nstep_match = re.search(r"NStep=\s*(\d+)", text)
        nstep = int(nstep_match.group(1)) if nstep_match else text.count("GradGradGradGrad")
        if nstep >= 60:
            return TS_FAILURE_OSCILLATION
        return TS_FAILURE_OSCILLATION

    if "Normal termination" in text:
        tail = "\n".join(text.splitlines()[-50:])
        if "imaginary" in tail.lower() or "frequencies" in tail.lower():
            return TS_FAILURE_NO_IMAG
        failed_xyz = log_path.with_suffix("").with_name(log_path.stem + "_ts_failed.xyz")
        if failed_xyz.exists():
            return TS_FAILURE_NO_IMAG
        return TS_FAILURE_NO_IMAG

    return TS_FAILURE_UNKNOWN


def _prepare_rescue_start(branch_dir: Path, failure_kind: str, old_log: Path | None):
    s2_ts_guess = branch_dir / "S2_Retro" / "ts_guess.xyz"
    s3_dir = branch_dir / "S3_TS"
    s3_dir.mkdir(parents=True, exist_ok=True)
    rescue_xyz = s3_dir / "rescue_start_from_s2_guess.xyz"
    shutil.copy2(s2_ts_guess, rescue_xyz)
    logger.info("Rescue restart from S2 initial guess (failure_kind=%s): %s", failure_kind, rescue_xyz)
    return rescue_xyz


def _rescue_one_branch(
    branch_dir: Path,
    config: dict,
    charge: int,
    spin: int,
    enable_l2_sp: bool,
    clean: bool = False,
) -> bool:
    """Rescue a single branch. Returns True on success."""
    branch_dir = branch_dir.resolve()
    if not branch_dir.exists():
        logger.error("Branch directory does not exist: %s", branch_dir)
        return False

    old_log = _find_failed_log(branch_dir)
    if old_log:
        logger.info("[%s] Found failed Gaussian log: %s", branch_dir.name, old_log)
    else:
        logger.warning("[%s] No failed Gaussian log found", branch_dir.name)

    failure_kind = _classify_from_log(old_log) if old_log else TS_FAILURE_UNKNOWN
    logger.info("[%s] Inferred failure kind: %s", branch_dir.name, failure_kind)

    if clean:
        rescue_root = branch_dir / "S3_TS" / "ts_opt" / "rescue"
        if rescue_root.exists():
            try:
                shutil.rmtree(rescue_root)
            except PermissionError:
                import subprocess
                subprocess.run(["rm", "-rf", str(rescue_root)], check=False)
            logger.info("[%s] Cleaned old rescue directory: %s", branch_dir.name, rescue_root)

    rescue_xyz = _prepare_rescue_start(branch_dir, failure_kind, old_log)
    logger.info("[%s] Rescue start geometry: %s", branch_dir.name, rescue_xyz)
    if not rescue_xyz.exists():
        logger.error("[%s] Rescue start geometry not found: %s", branch_dir.name, rescue_xyz)
        return False

    runner = QCTaskRunner(config=config)
    s3_dir = branch_dir / "S3_TS"

    result = runner.run_ts_rescue_only(
        rescue_start_xyz=rescue_xyz,
        output_dir=s3_dir,
        failure_kind=failure_kind,
        charge=charge,
        spin=spin,
        enable_l2_sp=enable_l2_sp,
    )

    if result.converged:
        logger.info("[%s] === Rescue SUCCEEDED ===", branch_dir.name)
        logger.info("[%s] Optimized TS: %s", branch_dir.name, result.optimized_xyz)
        logger.info("[%s] Method used: %s", branch_dir.name, result.method_used)
        if result.l2_energy is not None:
            logger.info("[%s] L2 SP energy: %.8f Hartree", branch_dir.name, result.l2_energy)
        return True
    else:
        logger.error("[%s] === Rescue FAILED ===", branch_dir.name)
        logger.error("[%s] Error: %s", branch_dir.name, result.error_message or "(none)")
        if result.failure_kind:
            logger.error("[%s] Failure kind: %s", branch_dir.name, result.failure_kind)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="S3 TS rescue-only rerun — skip primary Berny, jump to rescue. Supports one or more branches."
    )
    parser.add_argument(
        "branch_dirs", type=Path, nargs="+",
        help="One or more branch work directories (e.g. rph_output/RXN_*/branches/BR_DR_001)",
    )
    parser.add_argument(
        "--config", type=Path, default=REPO_ROOT / "config" / "defaults.yaml",
        help="Path to config YAML",
    )
    parser.add_argument(
        "--no-l2-sp", action="store_true",
        help="Skip L2 single-point even if rescue succeeds",
    )
    parser.add_argument("--charge", type=int, default=0, help="Molecular charge")
    parser.add_argument("--spin", type=int, default=1, help="Spin multiplicity")
    parser.add_argument(
        "--clean", action="store_true",
        help="Auto-remove old rescue directories before rerunning",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    )

    config = _load_config(args.config)
    logger.info("Loaded config from %s", args.config)

    results = {}
    for branch_dir in args.branch_dirs:
        success = _rescue_one_branch(
            branch_dir, config,
            charge=args.charge, spin=args.spin,
            enable_l2_sp=not args.no_l2_sp,
            clean=args.clean,
        )
        results[str(branch_dir)] = success

    succeeded = sum(1 for v in results.values() if v)
    failed = len(results) - succeeded
    logger.info("=== Batch rescue summary: %d succeeded, %d failed ===", succeeded, failed)
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
