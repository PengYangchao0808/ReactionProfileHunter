import json
from pathlib import Path

import pandas as pd

from rph_core.steps.condition_feature_merger import ConditionFeatureMerger
from rph_core.steps.condition_thermo import ConditionThermoCalculator


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def test_condition_thermo_and_merge_smoke(tmp_path: Path, monkeypatch) -> None:
    reaction_root = tmp_path / "RXN_deadbeef"
    s3_dir = reaction_root / "S3_TS"
    condition_root = reaction_root / "conditions" / "COND_row_0001"
    reaction_features_dir = reaction_root / "reaction_features"
    s3_dir.mkdir(parents=True, exist_ok=True)
    condition_root.mkdir(parents=True, exist_ok=True)
    reaction_features_dir.mkdir(parents=True, exist_ok=True)

    ts_log = s3_dir / "ts_opt" / "standard" / "ts_freq.log"
    reactant_log = s3_dir / "reactant_opt" / "standard" / "reactant_freq.log"
    ts_log.parent.mkdir(parents=True, exist_ok=True)
    reactant_log.parent.mkdir(parents=True, exist_ok=True)
    ts_log.write_text("dummy", encoding="utf-8")
    reactant_log.write_text("dummy", encoding="utf-8")

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
        return shermo_runner.ThermoResult(
            g_sum=-100.0,
            h_sum=-99.0,
            u_sum=-98.0,
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
