from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from rph_core.steps.mechanism_classifier.s0_record import load_s0_reaction_record
from rph_core import v4_orchestrator


REPO_ROOT = Path(__file__).resolve().parents[1]
DATASET = REPO_ROOT / "data" / "reaxys_cleaned.csv"


def test_reaxys_record_maps_forming_bonds_to_product_heavy_atom_indices():
    record = load_s0_reaction_record(DATASET, "1")

    assert record.reaction_type == "4+3"
    assert record.mapped_forming_bonds == ((5, 11), (7, 8))
    assert record.forming_bonds == ((18, 11), (15, 14))
    assert record.mapping_trusted is True


def test_csv_cli_passes_authoritative_s0_record_to_orchestrator(monkeypatch, tmp_path):
    captured = {}

    class FakeOrchestrator:
        def __init__(self, config_path, color=False):
            captured["config"] = config_path
            captured["color"] = color

        def run(self, s0_record, work_dir, stop_after):
            captured.update(
                product_smiles=s0_record.product_smiles,
                work_dir=work_dir,
                reaction_type=s0_record.reaction_type,
                stop_after=stop_after,
                s0_record=s0_record,
            )

    monkeypatch.setattr(v4_orchestrator, "V4Orchestrator", FakeOrchestrator)
    assert v4_orchestrator.main([
        "--csv", str(DATASET),
        "--rx-id", "1",
        "--output", str(tmp_path / "run"),
        "--stop-after", "s3",
    ]) == 0

    assert captured["reaction_type"] == "4+3"
    assert captured["stop_after"] == "s3"
    assert captured["s0_record"].forming_bonds == ((18, 11), (15, 14))


def test_record_driven_s0_is_written_before_s1_and_preserved_at_s1_stop(monkeypatch, tmp_path):
    record = load_s0_reaction_record(DATASET, "1")

    class FakeCensoLiteEngine:
        def __init__(self, _config, stage_dir, _molecule_name):
            self.stage_dir = stage_dir

        def run(self, _smiles):
            product_dir = self.stage_dir / "product"
            candidate_dir = product_dir / "candidates"
            candidate_dir.mkdir(parents=True)
            xyz = candidate_dir / "conf_0001.xyz"
            xyz.write_text("1\nfixture\nH 0.0 0.0 0.0\n", encoding="utf-8")
            manifest = product_dir / "manifest.json"
            manifest.write_text(
                '{"selected": "product_conf_0001", "candidates": '
                '[{"id": "product_conf_0001", "xyz": "product/candidates/conf_0001.xyz"}]}',
                encoding="utf-8",
            )
            return {"manifest": manifest, "selected_xyz": xyz, "data": __import__("json").loads(manifest.read_text())}

    monkeypatch.setattr(v4_orchestrator, "CensoLiteEngine", FakeCensoLiteEngine)
    orchestrator = object.__new__(v4_orchestrator.V4Orchestrator)
    orchestrator.config = {"step1": {"protocol": "censo_lite"}}

    result = orchestrator.run(
        record,
        tmp_path / "run",
        stop_after="s1",
    )

    s0 = __import__("json").loads((tmp_path / "run" / "S0_Mechanism" / "mechanism.json").read_text())
    assert s0["source"] == "trusted_reaction_record"
    assert s0["forming_bonds"] == [[18, 11], [15, 14]]
    assert "s0_manifest" in result
    assert result["completed_through"] == "s1"
