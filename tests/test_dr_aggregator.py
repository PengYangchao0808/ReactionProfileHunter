import csv
import json
from pathlib import Path

from rph_core.steps.dr_aggregator import DRAggregator


def _write_features_csv(csv_path: Path, barrier: float, column: str = "thermo.dG_activation") -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["thermo.dG_activation", "thermo.dE_activation", "extra"])
        writer.writerow([str(barrier), str(barrier + 0.5), "x"])


class TestDRAggregator:
    def test_single_branch_returns_insufficient_data(self, tmp_path: Path) -> None:
        _write_features_csv(tmp_path / "S4_Data" / "features_raw.csv", 12.0)
        agg = DRAggregator()
        result = agg.run(tmp_path, temperature_k=298.15)
        assert result["predicted_dr"]["status"] == "insufficient_data"

    def test_two_branches_computes_ratio(self, tmp_path: Path) -> None:
        _write_features_csv(tmp_path / "S4_Data" / "features_raw.csv", 12.0)
        br = tmp_path / "branches" / "BR_DR_001"
        _write_features_csv(br / "S4_Data" / "features_raw.csv", 13.5)

        agg = DRAggregator()
        result = agg.run(tmp_path, temperature_k=298.15)
        pred = result["predicted_dr"]
        assert pred["status"] == "complete"
        assert pred["major_branch_id"] == "BR_MAJOR"
        assert pred["minor_branch_id"] == "BR_DR_001"
        assert pred["ddg_kcal_mol"] > 0
        assert pred["major_relative_population"] > pred.get("minor_relative_population", 0)

    def test_missing_branch_degraded_with_warning(self, tmp_path: Path) -> None:
        _write_features_csv(tmp_path / "S4_Data" / "features_raw.csv", 12.0)
        (tmp_path / "branches" / "BR_DR_001").mkdir(parents=True)

        agg = DRAggregator()
        result = agg.run(tmp_path, temperature_k=298.15)
        assert len(result["warnings"]) >= 1
        assert result["predicted_dr"]["status"] == "insufficient_data"

    def test_write_produces_valid_json(self, tmp_path: Path) -> None:
        _write_features_csv(tmp_path / "S4_Data" / "features_raw.csv", 10.0)
        br = tmp_path / "branches" / "BR_DR_001"
        _write_features_csv(br / "S4_Data" / "features_raw.csv", 11.0)

        agg = DRAggregator()
        out = tmp_path / "reaction_features" / "dr_prediction.json"
        result_path = agg.write(tmp_path, out)
        assert result_path.exists()
        with open(result_path) as f:
            data = json.load(f)
        assert data["version"] == "rph-dr-prediction-v1"
        assert data["predicted_dr"]["status"] == "complete"
