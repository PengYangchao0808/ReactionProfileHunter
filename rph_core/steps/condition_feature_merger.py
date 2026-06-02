from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

from rph_core.utils.json_io import read_json_dict, write_json

logger = logging.getLogger(__name__)


def _load_reaction_features(reaction_features_dir: Path) -> Dict[str, Any]:
    json_path = reaction_features_dir / "geo_electronic_features.json"
    if json_path.exists():
        return read_json_dict(json_path)

    csv_path = reaction_features_dir / "geo_electronic_features.csv"
    if csv_path.exists():
        import pandas as pd

        df = pd.read_csv(csv_path)
        if df.empty:
            return {}
        row = df.iloc[0].to_dict()
        return {str(k): (None if pd.isna(v) else v) for k, v in row.items()}

    return {}


@dataclass(frozen=True)
class ConditionFeatureMerger:
    reaction_root: Path
    condition_root: Path
    feature_root: Optional[Path] = None
    branch_id: Optional[str] = None

    def run(self) -> tuple[Path, Path]:
        merged_csv = self.condition_root / "merged_features.csv"
        merged_json = self.condition_root / "merged_features.json"

        if merged_csv.exists() and merged_csv.stat().st_size > 0 and merged_json.exists() and merged_json.stat().st_size > 0:
            return merged_csv, merged_json

        reaction_features_dir = (self.feature_root or self.reaction_root) / "reaction_features"
        reaction_features = _load_reaction_features(reaction_features_dir)

        s4_csv_path = (self.feature_root or self.reaction_root) / "S4_Data" / "features_raw.csv"
        if not reaction_features and s4_csv_path.exists():
            import pandas as pd
            df = pd.read_csv(s4_csv_path)
            if not df.empty:
                row = df.iloc[0].to_dict()
                reaction_features = {str(k): (None if pd.isna(v) else v) for k, v in row.items()}

        condition_manifest_path = self.condition_root / "condition_manifest.json"
        thermo_path = self.condition_root / "thermo_calculation.json"

        if not condition_manifest_path.exists():
            raise FileNotFoundError(f"Missing condition_manifest.json: {condition_manifest_path}")
        if not thermo_path.exists():
            raise FileNotFoundError(f"Missing thermo_calculation.json: {thermo_path}")

        condition_manifest = read_json_dict(condition_manifest_path)
        thermo = read_json_dict(thermo_path)

        merged: Dict[str, Any] = {}

        merged["reaction_id"] = condition_manifest.get("reaction_id")
        merged["condition_id"] = condition_manifest.get("condition_id")
        merged["row_id"] = condition_manifest.get("row_id")
        if self.branch_id is not None:
            merged["branch_id"] = self.branch_id

        for k, v in reaction_features.items():
            merged[str(k)] = v

        for k, v in condition_manifest.items():
            if k in {"reaction_id", "condition_id", "row_id"}:
                continue
            merged[f"cond.{k}"] = v

        merged["thermo.temperature_K"] = thermo.get("temperature_K")
        raw_delta = thermo.get("delta")
        delta: Dict[str, Any] = raw_delta if isinstance(raw_delta, dict) else {}
        raw_ts = thermo.get("ts")
        ts: Dict[str, Any] = raw_ts if isinstance(raw_ts, dict) else {}
        raw_intermediate = thermo.get("intermediate")
        intermediate: Dict[str, Any] = raw_intermediate if isinstance(raw_intermediate, dict) else {}
        raw_reactant = thermo.get("reactant")
        if isinstance(raw_reactant, dict) and not intermediate:
            intermediate = raw_reactant

        merged["thermo.dG_activation"] = delta.get("dG_activation")
        merged["thermo.dH_activation"] = delta.get("dH_activation")
        merged["thermo.dS_activation"] = delta.get("dS_activation")
        merged["thermo.ts_g_kcal"] = ts.get("g_kcal")
        merged["thermo.ts_h_kcal"] = ts.get("h_kcal")
        merged["thermo.ts_s_cal_mol_k"] = ts.get("s_cal_mol_k")
        merged["thermo.intermediate_g_kcal"] = intermediate.get("g_kcal")
        merged["thermo.intermediate_h_kcal"] = intermediate.get("h_kcal")
        merged["thermo.intermediate_s_cal_mol_k"] = intermediate.get("s_cal_mol_k")

        # New species and steps fields from ConditionThermoCalculator v2
        raw_species = thermo.get("species")
        species: Dict[str, Any] = raw_species if isinstance(raw_species, dict) else {}
        raw_steps = thermo.get("steps")
        steps: Dict[str, Any] = raw_steps if isinstance(raw_steps, dict) else {}

        merged["thermo.dG_oxidation"] = (steps.get("oxidation") or {}).get("dG")
        merged["thermo.dG_cycloaddition_activation"] = (steps.get("cycloaddition_activation") or {}).get("dG")
        merged["thermo.dG_activation_from_precursor"] = (steps.get("overall_activation") or {}).get("dG")
        merged["thermo.dG_total_reaction"] = (steps.get("overall_reaction") or {}).get("dG")

        reserved_species = {"ts", "intermediate", "product"}
        for sp_key, sp_data in species.items():
            if sp_key in reserved_species:
                continue
            sp = sp_data if isinstance(sp_data, dict) else {}
            merged[f"thermo.{sp_key}_g_kcal"] = sp.get("g_kcal")
            merged[f"thermo.{sp_key}_h_kcal"] = sp.get("h_kcal")
            merged[f"thermo.{sp_key}_s_cal_mol_k"] = sp.get("s_cal_mol_k")

        write_json(merged_json, merged)

        import pandas as pd

        df = pd.DataFrame([merged])
        merged_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(merged_csv, index=False)

        return merged_csv, merged_json
