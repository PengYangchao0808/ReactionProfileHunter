from __future__ import annotations

import csv
import json
from pathlib import Path

from rph_core.steps.condition_thermo import ConditionThermoCalculator
from rph_core.steps.condition_feature_merger import ConditionFeatureMerger
from rph_core.steps.dr_aggregator import DRAggregator
from rph_core.steps.condition_thermo import _safe_path_from_json


class TestConditionThermoBranchAware:
    def test_artifact_root_overrides_reaction_root(self, tmp_path: Path):
        reaction_root = tmp_path / "reaction"
        artifact_root = tmp_path / "branch"
        condition_root = tmp_path / "condition"
        condition_root.mkdir(parents=True)
        (condition_root / "condition_manifest.json").write_text('{"temperature_K": 298.15}')

        calc = ConditionThermoCalculator(
            config={"executables": {}, "thermo": {"temperature_k": 298.15}},
            reaction_root=reaction_root,
            condition_root=condition_root,
            artifact_root=artifact_root,
        )
        assert calc.artifact_root == artifact_root

    def test_legacy_mode_no_artifact_root(self, tmp_path: Path):
        reaction_root = tmp_path / "reaction"
        condition_root = tmp_path / "condition"
        condition_root.mkdir(parents=True)
        (condition_root / "condition_manifest.json").write_text('{"temperature_K": 298.15}')

        calc = ConditionThermoCalculator(
            config={"executables": {}, "thermo": {"temperature_k": 298.15}},
            reaction_root=reaction_root,
            condition_root=condition_root,
        )
        assert calc.artifact_root is None

    def test_stale_missing_sp_metadata_placeholder_is_recomputed(self, tmp_path: Path):
        reaction_root = tmp_path / "reaction"
        artifact_root = tmp_path / "branch"
        condition_root = tmp_path / "condition"
        condition_root.mkdir(parents=True)
        (condition_root / "condition_manifest.json").write_text('{"temperature_K": 298.15}')
        (condition_root / "thermo_calculation.json").write_text(json.dumps({
            "temperature_K": 298.15,
            "species": {"ts": {"g_kcal": None}},
            "steps": {"cycloaddition_activation": {"dG": None}},
            "warnings": ["W_MISSING_SP_MATRIX_METADATA:old"],
        }))
        s3_dir = artifact_root / "S3_TS"
        s3_dir.mkdir(parents=True)
        (s3_dir / "sp_matrix_metadata.json").write_text(json.dumps({
            "method": "M062X",
            "solvent": "acetone",
            "e_ts": -100.0,
            "e_reactant": -101.0,
            "e_product": -102.0,
        }))

        thermo_path = ConditionThermoCalculator(
            config={"executables": {}, "thermo": {"temperature_k": 298.15}},
            reaction_root=reaction_root,
            condition_root=condition_root,
            artifact_root=artifact_root,
        ).run()

        data = json.loads(thermo_path.read_text())
        assert data["metadata"]["method"] == "M062X"
        assert not any(str(w).startswith("W_MISSING_SP_MATRIX_METADATA") for w in data["warnings"])


class TestConditionFeatureMergerBranchAware:
    def test_feature_root_and_branch_id(self, tmp_path: Path):
        merger = ConditionFeatureMerger(
            reaction_root=tmp_path / "reaction",
            condition_root=tmp_path / "condition",
            feature_root=tmp_path / "branch",
            branch_id="BR_DR_001",
        )
        assert merger.feature_root == tmp_path / "branch"
        assert merger.branch_id == "BR_DR_001"

    def test_legacy_mode_no_optional_fields(self, tmp_path: Path):
        merger = ConditionFeatureMerger(
            reaction_root=tmp_path / "reaction",
            condition_root=tmp_path / "condition",
        )
        assert merger.feature_root is None
        assert merger.branch_id is None


class TestDRAggregatorConditionMode:
    def test_run_for_condition_reads_merged_csv(self, tmp_path: Path):
        reaction_root = tmp_path / "reaction"
        condition_root = tmp_path / "condition"
        branches_dir = condition_root / "branches" / "BR_MAJOR"
        branches_dir.mkdir(parents=True)

        csv_content = "thermo.dG_activation,other\n15.5,test\n"
        (branches_dir / "merged_features.csv").write_text(csv_content)

        dr_branch_dir = condition_root / "branches" / "BR_DR_001"
        dr_branch_dir.mkdir(parents=True)
        (dr_branch_dir / "merged_features.csv").write_text("thermo.dG_activation,other\n18.0,minor\n")

        s4_dir = reaction_root / "S4_Data"
        s4_dir.mkdir(parents=True)
        (s4_dir / "features_raw.csv").write_text("thermo.dG_activation,other\n12.0,major\n")

        aggregator = DRAggregator()
        result = aggregator.run_for_condition(
            reaction_root=reaction_root,
            condition_root=condition_root,
            temperature_k=298.15,
        )

        assert result["temperature_k"] == 298.15
        assert len(result["branches"]) == 2
        assert result["branches"][0]["branch_id"] == "BR_MAJOR"
        assert result["branches"][0]["delta_g_act_kcal_mol"] == 15.5
        assert result["branches"][1]["branch_id"] == "BR_DR_001"
        assert result["branches"][1]["delta_g_act_kcal_mol"] == 18.0

        output = condition_root / "dr_prediction.json"
        assert output.exists()


class TestConditionThermoPathSafety:
    def test_safe_path_rejects_absolute_path_outside_s3_root(self, tmp_path: Path):
        s3_dir = tmp_path / "branch" / "S3_TS"
        s3_dir.mkdir(parents=True)
        outside = tmp_path / "outside.log"
        outside.write_text("not a freq log")

        assert _safe_path_from_json(outside, s3_dir) is None

    def test_safe_path_accepts_path_inside_s3_root(self, tmp_path: Path):
        s3_dir = tmp_path / "branch" / "S3_TS"
        log_path = s3_dir / "ts_opt" / "output.log"
        log_path.parent.mkdir(parents=True)
        log_path.write_text("freq")

        assert _safe_path_from_json("ts_opt/output.log", s3_dir) == log_path


class TestConditionThermoConditionMLR:
    def test_write_condition_features_mlr_includes_branch_sample_id(self, tmp_path: Path):
        reaction_root = tmp_path / "reaction"
        condition_root = tmp_path / "condition" / "branches" / "BR_MAJOR"
        condition_root.mkdir(parents=True)

        (condition_root / "condition_manifest.json").write_text(json.dumps({
            "reaction_id": "RXN_demo",
            "condition_id": "COND_demo",
            "temperature_K": 298.15,
        }))
        (condition_root / "thermo_calculation.json").write_text(json.dumps({
            "temperature_K": 298.15,
            "delta": {
                "dE_activation": 1.0,
                "dE_reaction": 2.0,
                "dG_activation": 3.0,
                "dG_reaction": 4.0,
                "dG_activation_gibbs": 3.0,
                "dG_reaction_gibbs": 4.0,
                "dG_act_step1": 5.0,
                "dG_act_full": 6.0,
                "Keq_act": 7.0,
            },
            "ts": {"g_kcal": 10.0, "h_kcal": 11.0, "s_cal_mol_k": 12.0},
            "intermediate": {"g_kcal": 13.0, "h_kcal": 14.0, "s_cal_mol_k": 15.0},
            "species": {
                "ts": {"g_kcal": 10.0, "h_kcal": 11.0, "s_cal_mol_k": 12.0},
                "intermediate": {"g_kcal": 13.0, "h_kcal": 14.0, "s_cal_mol_k": 15.0},
                "product": {"g_kcal": 16.0, "h_kcal": 17.0, "s_cal_mol_k": 18.0},
                "precursor": {"g_kcal": 19.0, "h_kcal": 20.0, "s_cal_mol_k": 21.0},
            },
            "conformer": {
                "g_precursor_boltzmann_avg_kcal": 18.5,
                "conformer_N_eff": 1.5,
                "conformer_S_conf_cal_mol_k": 2.5,
                "conformer_top_k_energies": [0.0, 0.3],
                "n_conformers": 2,
            },
        }))

        csv_path = ConditionThermoCalculator(
            config={"executables": {}, "thermo": {"temperature_k": 298.15}},
            reaction_root=reaction_root,
            condition_root=condition_root,
            branch_id="BR_MAJOR",
        ).write_condition_features_mlr(condition_root)

        with open(csv_path, "r", encoding="utf-8") as handle:
            row = next(csv.DictReader(handle))

        assert row["branch_id"] == "BR_MAJOR"
        assert row["sample_id"] == "RXN_demo__COND_demo__BR_MAJOR"
