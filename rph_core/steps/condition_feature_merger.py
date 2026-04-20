from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict


logger = logging.getLogger(__name__)


def _read_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def _load_reaction_features(reaction_features_dir: Path) -> Dict[str, Any]:
    json_path = reaction_features_dir / "geo_electronic_features.json"
    if json_path.exists():
        return _read_json(json_path)

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

    def run(self) -> tuple[Path, Path]:
        merged_csv = self.condition_root / "merged_features.csv"
        merged_json = self.condition_root / "merged_features.json"

        if merged_csv.exists() and merged_csv.stat().st_size > 0 and merged_json.exists() and merged_json.stat().st_size > 0:
            return merged_csv, merged_json

        reaction_features_dir = self.reaction_root / "reaction_features"
        reaction_features = _load_reaction_features(reaction_features_dir)

        condition_manifest_path = self.condition_root / "condition_manifest.json"
        thermo_path = self.condition_root / "thermo_calculation.json"

        if not condition_manifest_path.exists():
            raise FileNotFoundError(f"Missing condition_manifest.json: {condition_manifest_path}")
        if not thermo_path.exists():
            raise FileNotFoundError(f"Missing thermo_calculation.json: {thermo_path}")

        condition_manifest = _read_json(condition_manifest_path)
        thermo = _read_json(thermo_path)

        merged: Dict[str, Any] = {}

        merged["reaction_id"] = condition_manifest.get("reaction_id")
        merged["condition_id"] = condition_manifest.get("condition_id")
        merged["row_id"] = condition_manifest.get("row_id")

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
        raw_reactant = thermo.get("reactant")
        reactant: Dict[str, Any] = raw_reactant if isinstance(raw_reactant, dict) else {}

        merged["thermo.dG_activation"] = delta.get("dG_activation")
        merged["thermo.dH_activation"] = delta.get("dH_activation")
        merged["thermo.dS_activation"] = delta.get("dS_activation")
        merged["thermo.ts_g_kcal"] = ts.get("g_kcal")
        merged["thermo.ts_h_kcal"] = ts.get("h_kcal")
        merged["thermo.ts_s_cal_mol_k"] = ts.get("s_cal_mol_k")
        merged["thermo.reactant_g_kcal"] = reactant.get("g_kcal")
        merged["thermo.reactant_h_kcal"] = reactant.get("h_kcal")
        merged["thermo.reactant_s_cal_mol_k"] = reactant.get("s_cal_mol_k")

        _write_json(merged_json, merged)

        import pandas as pd

        df = pd.DataFrame([merged])
        merged_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(merged_csv, index=False)

        return merged_csv, merged_json
