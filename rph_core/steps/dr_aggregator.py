from __future__ import annotations

import csv
import json
import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

R_KCAL = 1.987203e-3


class DRAggregator:

    def run(
        self,
        reaction_root: Path,
        temperature_k: float = 298.15,
        parent_label: str = "BR_MAJOR",
    ) -> Dict[str, Any]:
        """Aggregate branch activation barriers into a predicted DR.

        Reads ``S4_Data/features_raw.csv`` from the parent (``reaction_root``)
        and from each ``branches/{branch_id}/`` subdirectory, extracts
        ``thermo.dG_activation``, and computes the Boltzmann-weighted DR ratio.

        Returns a ``dr_prediction.json``-ready dictionary. Degraded branches
        (missing csv or barrier) are recorded in ``warnings`` without
        aborting the whole reaction.
        """

        branch_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []

        branches_root = reaction_root / "branches"
        v3_major_dir = branches_root / parent_label
        if v3_major_dir.exists():
            parent_barrier = _read_barrier(
                v3_major_dir / "S4_Data" / "features_raw.csv",
                f"S4_Data of {parent_label}",
                warnings,
            )
            branch_entries.append({
                "branch_id": parent_label,
                "source": str(v3_major_dir),
                "delta_g_act_kcal_mol": parent_barrier,
            })
        else:
            parent_barrier = _read_barrier(
                reaction_root / "S4_Data" / "features_raw.csv",
                f"S4_Data of {parent_label}",
                warnings,
            )
            branch_entries.append({
                "branch_id": parent_label,
                "source": str(reaction_root),
                "delta_g_act_kcal_mol": parent_barrier,
            })

        if branches_root.exists():
            for child_dir in sorted(branches_root.iterdir()):
                if not child_dir.is_dir():
                    continue
                child_branch_id = child_dir.name
                if child_branch_id == parent_label:
                    continue
                child_barrier = _read_barrier(
                    child_dir / "S4_Data" / "features_raw.csv",
                    f"S4_Data of {child_branch_id}",
                    warnings,
                )
                branch_entries.append({
                    "branch_id": child_branch_id,
                    "source": str(child_dir),
                    "delta_g_act_kcal_mol": child_barrier,
                })

        predicted_dr = _compute_predicted_dr(branch_entries, temperature_k, warnings)

        return {
            "version": "rph-dr-prediction-v1",
            "temperature_k": temperature_k,
            "branches": branch_entries,
            "predicted_dr": predicted_dr,
            "warnings": warnings,
        }

    def write(self, reaction_root: Path, output_path: Path, temperature_k: float = 298.15) -> Path:
        result = self.run(reaction_root, temperature_k=temperature_k)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        return output_path

    def run_for_condition(
        self,
        reaction_root: Path,
        condition_root: Path,
        temperature_k: float = 298.15,
    ) -> Dict[str, Any]:
        branch_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []

        branches_dir = condition_root / "branches"
        if branches_dir.exists():
            branch_dirs = sorted(
                (p for p in branches_dir.iterdir() if p.is_dir()),
                key=lambda p: (p.name != "BR_MAJOR", p.name),
            )
            for child_dir in branch_dirs:
                if not child_dir.is_dir():
                    continue
                child_branch_id = child_dir.name
                merged_csv = child_dir / "merged_features.csv"
                barrier = _read_barrier_from_merged(merged_csv, child_branch_id, warnings)
                if barrier is None:
                    parent_csv = reaction_root / "branches" / child_branch_id / "S4_Data" / "features_raw.csv"
                    barrier = _read_barrier(parent_csv, f"S4_Data of {child_branch_id}", warnings)
                branch_entries.append({
                    "branch_id": child_branch_id,
                    "source": str(child_dir),
                    "delta_g_act_kcal_mol": barrier,
                })

        if not any(entry.get("branch_id") == "BR_MAJOR" for entry in branch_entries):
            parent_csv = reaction_root / "branches" / "BR_MAJOR" / "S4_Data" / "features_raw.csv"
            if not parent_csv.exists():
                parent_csv = reaction_root / "S4_Data" / "features_raw.csv"
            parent_barrier = _read_barrier(parent_csv, "BR_MAJOR", warnings)
            branch_entries.insert(0, {
                "branch_id": "BR_MAJOR",
                "source": str(parent_csv.parent.parent),
                "delta_g_act_kcal_mol": parent_barrier,
            })

        predicted_dr = _compute_predicted_dr(branch_entries, temperature_k, warnings)
        result = {
            "version": "rph-dr-prediction-v1",
            "temperature_k": temperature_k,
            "branches": branch_entries,
            "predicted_dr": predicted_dr,
            "warnings": warnings,
        }

        output_path = condition_root / "dr_prediction.json"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        return result


def _read_barrier(
    csv_path: Path,
    label: str,
    warnings: List[str],
) -> Optional[float]:
    if not csv_path.exists():
        warnings.append(f"Missing features_raw.csv for {label}; cannot read barrier.")
        return None

    try:
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames or []
            rows = list(reader)
            candidates = [
                candidate
                for candidate in ("thermo.dG_activation", "thermo.dE_activation")
                if candidate in header
            ]
            if not candidates:
                warnings.append(f"No activation column found in {label} features_raw.csv.")
                return None

            for col in candidates:
                for row in rows:
                    raw = (row.get(col) or "").strip()
                    if not raw:
                        continue
                    val = float(raw)
                    if not math.isnan(val):
                        return val
    except Exception as exc:
        warnings.append(f"Failed to read barrier from {label}: {exc}")
        return None

    return None


def _read_barrier_from_merged(
    csv_path: Path,
    branch_id: str,
    warnings: List[str],
) -> Optional[float]:
    if not csv_path.exists():
        return None
    try:
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames or []
            rows = list(reader)
            candidates = [
                candidate
                for candidate in ("thermo.dG_activation", "thermo.dE_activation")
                if candidate in header
            ]
            if not candidates:
                return None
            for col in candidates:
                for row in rows:
                    raw = (row.get(col) or "").strip()
                    if not raw:
                        continue
                    val = float(raw)
                    if not math.isnan(val):
                        return val
    except Exception:
        return None
    return None


def _compute_predicted_dr(
    branch_entries: List[Dict[str, Any]],
    temperature_k: float,
    warnings: List[str],
) -> Dict[str, Any]:
    valid = [
        e for e in branch_entries
        if e["delta_g_act_kcal_mol"] is not None
        and not (
            isinstance(e["delta_g_act_kcal_mol"], float)
            and math.isnan(e["delta_g_act_kcal_mol"])
        )
    ]
    if len(valid) < 2:
        return {
            "major": None,
            "minor": None,
            "ratio": None,
            "ddg_kcal_mol": None,
            "status": "insufficient_data",
        }

    barriers = [float(e["delta_g_act_kcal_mol"]) for e in valid]
    min_barrier = min(barriers)

    rt = R_KCAL * temperature_k
    relative = [math.exp(-(b - min_barrier) / rt) for b in barriers]
    total = sum(relative)

    if total <= 0:
        return {
            "major": None,
            "minor": None,
            "ratio": None,
            "ddg_kcal_mol": None,
            "status": "degenerate",
        }

    major_id = ""
    minor_id = ""
    major_pop = 0.0
    minor_pop = 0.0
    ddg = 0.0
    status = "complete"
    ratio_val: float = float("inf")

    if len(barriers) == 2:
        idx_low = 0 if barriers[0] <= barriers[1] else 1
        idx_high = 1 - idx_low
        major_id = valid[idx_low]["branch_id"]
        minor_id = valid[idx_high]["branch_id"]
        major_pop = relative[idx_low] / total
        minor_pop = relative[idx_high] / total
        ddg = barriers[idx_high] - barriers[idx_low]

        ratio_val = float("inf")
        if minor_pop > 0:
            ratio_val = major_pop / minor_pop
    else:
        status = "partial"
        warnings.append("More than 2 branches present; returning raw relative populations.")

    return {
        "major": round(ratio_val, 2) if ratio_val != float("inf") else None,
        "minor": 1.0,
        "ratio": round(ratio_val, 2) if ratio_val is not None and ratio_val != float("inf") else None,
        "ddg_kcal_mol": round(ddg, 4),
        "major_branch_id": major_id,
        "minor_branch_id": minor_id,
        "major_relative_population": round(major_pop, 4),
        "minor_relative_population": round(minor_pop, 4),
        "status": status,
    }
