#!/usr/bin/env python3
"""Backfill S1 crest/fastsp conformer_ensemble_energies.json for old raw outputs.

For each S1_ConfGeneration/product/ under a branch root:
  - If crest/crest.energies exists but conformer_ensemble_energies.json is missing,
    generate it from crest.energies.
  - If fastsp/conf_*/ directories contain SP output logs but
    conformer_ensemble_energies.json is missing, generate it.

Usage:
    python scripts/backfill_s1_ensemble_energies.py external_data/rph_output_backup
    python scripts/backfill_s1_ensemble_energies.py external_data/rph_output_backup --dry-run
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

HARTREE_TO_KCAL = 627.509608
R_KCAL_MOL_K = 1.987203e-3

SCF_DONE_RE = re.compile(r"SCF Done:\s+E\([^)]+\)\s*=\s*([-+]?\d+\.\d+)")
ORCA_SP_RE = re.compile(r"FINAL SINGLE POINT ENERGY\s+([-+]?\d+\.\d+)")


def _parse_crest_energies(energies_path: Path) -> List[float]:
    energies = []
    with open(energies_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    energies.append(float(parts[1]))
                except (ValueError, IndexError):
                    continue
            elif len(parts) == 1:
                try:
                    energies.append(float(parts[0]))
                except ValueError:
                    continue
    return energies


def _extract_last_sp_energy(log_path: Path) -> Optional[float]:
    try:
        content = log_path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    for pattern in (ORCA_SP_RE, SCF_DONE_RE):
        matches = pattern.findall(content)
        if matches:
            try:
                return float(matches[-1])
            except (ValueError, IndexError):
                continue
    return None


def _parse_fastsp_energies(fastsp_dir: Path) -> List[float]:
    energies_hartree = []
    conf_dirs = sorted(
        [d for d in fastsp_dir.iterdir() if d.is_dir()], key=lambda d: d.name
    )
    for conf_dir in conf_dirs:
        for log_path in sorted(conf_dir.glob("*.log")) + sorted(
            conf_dir.glob("*.out")
        ):
            energy = _extract_last_sp_energy(log_path)
            if energy is not None:
                energies_hartree.append(energy)
                break
    return energies_hartree


def _write_ensemble_json(
    output_path: Path,
    energies_kcal: List[float],
    source: str,
    source_file: Optional[str] = None,
    screening_method: Optional[str] = None,
    temperature_k: float = 298.15,
) -> None:
    min_e = min(energies_kcal)
    relative_energies = [e - min_e for e in energies_kcal]
    payload: Dict[str, Any] = {
        "source": source,
        "n_conformers": len(relative_energies),
        "energies": relative_energies,
        "energy_unit": "kcal/mol",
        "temperature_K": temperature_k,
    }
    if source_file:
        payload["source_file"] = source_file
    if screening_method:
        payload["screening_method"] = screening_method
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def backfill_branch(branch_root: Path, dry_run: bool = False) -> Dict[str, int]:
    stats: Dict[str, int] = {
        "crest_backfilled": 0,
        "crest_skipped_exists": 0,
        "crest_skipped_no_file": 0,
        "fastsp_backfilled": 0,
        "fastsp_skipped_exists": 0,
        "fastsp_skipped_no_data": 0,
    }
    s1_dir = branch_root / "S1_ConfGeneration" / "product"

    crest_dir = s1_dir / "crest"
    crest_energies = crest_dir / "crest.energies"
    crest_json = crest_dir / "conformer_ensemble_energies.json"
    if crest_energies.exists():
        if crest_json.exists():
            stats["crest_skipped_exists"] += 1
        else:
            energies = _parse_crest_energies(crest_energies)
            if energies:
                if not dry_run:
                    _write_ensemble_json(
                        crest_json,
                        energies,
                        source="crest_gfn2",
                        source_file=str(crest_energies),
                    )
                stats["crest_backfilled"] += 1
            else:
                stats["crest_skipped_no_file"] += 1
    else:
        stats["crest_skipped_no_file"] += 1

    fastsp_dir = s1_dir / "fastsp"
    fastsp_json = fastsp_dir / "conformer_ensemble_energies.json"
    if fastsp_dir.exists():
        if fastsp_json.exists():
            stats["fastsp_skipped_exists"] += 1
        else:
            energies_hartree = _parse_fastsp_energies(fastsp_dir)
            if energies_hartree:
                min_e = min(energies_hartree)
                energies_kcal = [(e - min_e) * HARTREE_TO_KCAL for e in energies_hartree]
                if not dry_run:
                    _write_ensemble_json(
                        fastsp_json,
                        energies_kcal,
                        source="fastsp_dft",
                        screening_method="r2SCAN-3c",
                    )
                stats["fastsp_backfilled"] += 1
            else:
                stats["fastsp_skipped_no_data"] += 1
    else:
        stats["fastsp_skipped_no_data"] += 1

    return stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Backfill S1 crest/fastsp conformer_ensemble_energies.json"
    )
    parser.add_argument(
        "root",
        type=Path,
        help="Root directory containing RXN_*/branches/BR_*/ directories",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Report what would be done without writing"
    )
    args = parser.parse_args()

    root = args.root
    if not root.is_dir():
        print(f"ERROR: root directory not found: {root}", file=sys.stderr)
        sys.exit(1)

    branch_dirs = sorted(root.glob("RXN_*/branches/BR_*"))
    if not branch_dirs:
        print(f"ERROR: no branch directories found under {root}", file=sys.stderr)
        sys.exit(1)

    total: Dict[str, int] = {}
    for br in branch_dirs:
        stats = backfill_branch(br, dry_run=args.dry_run)
        for key, value in stats.items():
            total[key] = total.get(key, 0) + value

    print(f"\nProcessed {len(branch_dirs)} branches")
    print(f"  crest backfilled:      {total.get('crest_backfilled', 0)}")
    print(f"  crest already exists:  {total.get('crest_skipped_exists', 0)}")
    print(f"  crest no crest.energies: {total.get('crest_skipped_no_file', 0)}")
    print(f"  fastsp backfilled:      {total.get('fastsp_backfilled', 0)}")
    print(f"  fastsp already exists:  {total.get('fastsp_skipped_exists', 0)}")
    print(f"  fastsp no SP data:      {total.get('fastsp_skipped_no_data', 0)}")

    if args.dry_run:
        print("(dry-run: no files written)")


if __name__ == "__main__":
    main()
