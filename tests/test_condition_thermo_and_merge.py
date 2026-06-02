import json
from pathlib import Path
from typing import Any

import pandas as pd

from rph_core.steps.condition_feature_merger import ConditionFeatureMerger
from rph_core.steps.condition_thermo import ConditionThermoCalculator


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def test_condition_thermo_and_merge_smoke(tmp_path: Path, monkeypatch) -> None:
    reaction_root = tmp_path / "RXN_deadbeef"
    s3_dir = reaction_root / "S3_TS"
    s1_dir = reaction_root / "S1_ConfGeneration"
    condition_root = reaction_root / "conditions" / "COND_row_0001"
    reaction_features_dir = reaction_root / "reaction_features"
    s3_dir.mkdir(parents=True, exist_ok=True)
    s1_dir.mkdir(parents=True, exist_ok=True)
    condition_root.mkdir(parents=True, exist_ok=True)
    reaction_features_dir.mkdir(parents=True, exist_ok=True)

    ts_log = s3_dir / "ts_opt" / "standard" / "ts_freq.log"
    reactant_log = s3_dir / "reactant_opt" / "standard" / "reactant_freq.log"
    ts_log.parent.mkdir(parents=True, exist_ok=True)
    reactant_log.parent.mkdir(parents=True, exist_ok=True)
    ts_log.write_text("dummy", encoding="utf-8")
    reactant_log.write_text("dummy", encoding="utf-8")

    product_log = s1_dir / "product" / "finalDFT" / "product_freq.log"
    precursor_log = s1_dir / "precursor" / "finalDFT" / "precursor_freq.log"
    product_log.parent.mkdir(parents=True, exist_ok=True)
    precursor_log.parent.mkdir(parents=True, exist_ok=True)
    product_log.write_text("SCF Done:  E(RB3LYP) =  -99.000000", encoding="utf-8")
    precursor_log.write_text("SCF Done:  E(RB3LYP) =  -101.200000", encoding="utf-8")

    _write_json(
        s1_dir / "product" / "conformer_state.json",
        {
            "summary": {"best_conformer": "conf_000", "global_min_energy": -99.0},
            "conformers": {
                "conf_000": {
                    "record": {"log_file": str(product_log), "sp_energy": -99.0},
                }
            },
        },
    )
    _write_json(
        s1_dir / "precursor" / "conformer_state.json",
        {
            "summary": {"best_conformer": "conf_001", "global_min_energy": -101.2},
            "conformers": {
                "conf_001": {
                    "record": {"log_file": str(precursor_log), "sp_energy": -101.2},
                }
            },
        },
    )
    _write_json(
        s1_dir / "precursor" / "finalDFT" / "conformer_energies.json",
        [
            -101.2 * 627.509,
            -101.198 * 627.509,
        ],
    )

    _write_json(
        s3_dir / "sp_matrix_metadata.json",
        {
            "e_ts": -100.0,
            "e_reactant": -101.0,
            "e_product": -99.0,
            "method": "WB97X-D3BJ/def2-TZVPP",
            "solvent": "dcm",
        },
    )
    _write_json(
        s3_dir / "s3_resume.json",
        {
            "ts_opt": {"completed": True, "freq_log": str(ts_log)},
            "reactant_opt": {"completed": True, "freq_log": str(reactant_log)},
        },
    )
    _write_json(
        condition_root / "condition_manifest.json",
        {
            "version": "rph_v2.1.1",
            "condition_id": "COND_row_0001",
            "reaction_id": "RXN_deadbeef",
            "row_id": "row_0001",
            "temperature_K": 350.15,
            "solvent": "DCM",
            "yield_pct": 80.0,
        },
    )
    _write_json(reaction_features_dir / "geo_electronic_features.json", {"geom.r1": 2.1})

    from rph_core.utils import shermo_runner

    def _fake_run_shermo(*, output_file: Path, **kwargs):
        output_file.write_text("fake", encoding="utf-8")
        g_sum = {
            "ts_Shermo.sum": -100.0,
            "intermediate_Shermo.sum": -101.0,
            "product_Shermo.sum": -99.5,
            "precursor_Shermo.sum": -101.2,
        }[output_file.name]
        return shermo_runner.ThermoResult(
            g_sum=g_sum,
            h_sum=g_sum + 0.2,
            u_sum=g_sum + 0.4,
            s_total=50.0,
            g_conc=None,
            output_file=output_file,
        )

    monkeypatch.setattr(shermo_runner, "run_shermo", _fake_run_shermo)

    config = {
        "executables": {"shermo": {"path": "Shermo"}},
        "thermo": {"pressure_atm": 1.0, "scl_zpe": 0.9905, "ilowfreq": 2},
    }

    thermo_path = ConditionThermoCalculator(
        config=config,
        reaction_root=reaction_root,
        condition_root=condition_root,
    ).run()
    assert thermo_path.exists()

    thermo = json.loads(thermo_path.read_text(encoding="utf-8"))
    assert thermo["species"]["product"]["g_kcal"] is not None
    assert thermo["species"]["precursor"]["g_kcal"] is not None
    assert thermo["conformer"]["n_conformers"] == 2

    condition_mlr_csv = condition_root / "condition_features_mlr.csv"
    assert condition_mlr_csv.exists()

    merged_csv, merged_json = ConditionFeatureMerger(
        reaction_root=reaction_root,
        condition_root=condition_root,
    ).run()
    assert merged_csv.exists()
    assert merged_json.exists()

    df = pd.read_csv(merged_csv)
    assert len(df) == 1
    row = df.iloc[0].to_dict()

    assert row["reaction_id"] == "RXN_deadbeef"
    assert row["condition_id"] == "COND_row_0001"
    assert row["cond.temperature_K"] == 350.15
    assert row["geom.r1"] == 2.1

    mlr_df = pd.read_csv(condition_mlr_csv)
    assert len(mlr_df) == 1
    mlr_row = mlr_df.iloc[0].to_dict()
    assert mlr_row["sample_id"] == "RXN_deadbeef__COND_row_0001"
    assert mlr_row["thermo.temperature_K"] == 350.15
    assert pd.notna(mlr_row["thermo.product_g_kcal"])
    assert pd.notna(mlr_row["thermo.precursor_boltzmann_g_kcal"])
    assert pd.notna(mlr_row["s1_Nconf_eff"])
