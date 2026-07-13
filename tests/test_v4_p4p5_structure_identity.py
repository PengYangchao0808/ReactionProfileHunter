from __future__ import annotations

import json
from collections.abc import Sequence
from importlib import import_module
from pathlib import Path
from typing import Any

v4_orchestrator = import_module("rph_core.v4_orchestrator")
S0ReactionRecord = import_module("rph_core.steps.mechanism_classifier.s0_record").S0ReactionRecord
lowlevel_module = import_module("rph_core.steps.step3_lowlevel.engine")
highlevel_module = import_module("rph_core.steps.step4_highlevel.engine")


def _record(rx_id: str = "RXN_P4P5") -> Any:
    return S0ReactionRecord(
        rx_id=rx_id,
        product_smiles="C=C",
        reaction_type="4+3",
        mapped_product_smiles="[CH2:1]=[CH2:2]",
        mapped_forming_bonds=((1, 2), (1, 2)),
        forming_bonds=((0, 1), (2, 3)),
        mapping_confidence=None,
        mapping_trusted=True,
        source_csv=Path("dataset.csv"),
        canonical_precursor_smiles="CC",
        mapped_precursor_smiles="[CH3:1][CH3:2]",
    )


def _orchestrator() -> Any:
    orchestrator = object.__new__(v4_orchestrator.V4Orchestrator)
    orchestrator.config = {
        "step1": {"protocol": "censo_lite"},
        "step2": {},
        "theory": {
            "s3_low_level": {},
            "s4_high_precision": {},
        },
    }
    return orchestrator


def _fake_s0_writer(variants: Sequence[dict[str, Any]]):
    def _write_s0(work_dir: Path, record: Any) -> Path:
        stage_dir = Path(work_dir) / "S0_Mechanism"
        stage_dir.mkdir(parents=True, exist_ok=True)
        mechanism_path = stage_dir / "mechanism.json"
        mechanism_path.write_text(
            json.dumps(
                {
                    "schema_version": "s0_mechanism_v2",
                    "stage": "S0",
                    "source": "trusted_reaction_record",
                    "rx_id": record.rx_id,
                    "forming_bonds": [list(pair) for pair in record.forming_bonds],
                    "variants": [variant["variant_id"] for variant in variants],
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        (stage_dir / "variant_registry.json").write_text(
            json.dumps(
                {
                    "schema_version": "s0_variant_registry_v1",
                    "reaction_id": record.rx_id,
                    "variants": variants,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        return mechanism_path

    return _write_s0


def _install_s1_s2_fakes(monkeypatch) -> None:
    class FakeCensoLiteEngine:
        def __init__(self, _config: object, stage_dir: Path, molecule_name: str):
            self.stage_dir = Path(stage_dir)
            self.molecule_name = molecule_name

        def run(self, smiles: str) -> dict[str, object]:
            molecule_dir = self.stage_dir / self.molecule_name
            candidate_dir = molecule_dir / "candidates"
            candidate_dir.mkdir(parents=True, exist_ok=True)
            candidate_xyz = candidate_dir / "conf_0001.xyz"
            candidate_xyz.write_text(f"1\n{smiles}\nH 0.0 0.0 0.0\n", encoding="utf-8")
            selected_xyz = molecule_dir / "selected.xyz"
            selected_xyz.write_text(candidate_xyz.read_text(encoding="utf-8"), encoding="utf-8")
            manifest = molecule_dir / "manifest.json"
            payload = {
                "selected": f"{self.molecule_name}_conf_0001",
                "selected_xyz": "selected.xyz",
                "candidates": [
                    {
                        "id": f"{self.molecule_name}_conf_0001",
                        "xyz": f"{self.molecule_name}/candidates/conf_0001.xyz",
                    }
                ],
            }
            manifest.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            return {
                "manifest": manifest,
                "selected_xyz": selected_xyz,
                "data": payload,
            }

    class FakePEBScanner:
        def __init__(self, _config: object, molecule_name: str | None = None):
            self.molecule_name = molecule_name or "unknown"

        def run(self, product_xyz: Path, output_dir: Path, forming_bonds, scan_config=None):
            del scan_config
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            ts_seed = output_dir / "ts_guess.xyz"
            intermediate_seed = output_dir / "intermediate_seed.xyz"
            profile = output_dir / "scan_profile.json"
            ts_seed.write_text(f"1\n{product_xyz}\nH 0.0 0.0 0.0\n", encoding="utf-8")
            intermediate_seed.write_text(f"1\n{product_xyz}\nH 1.0 0.0 0.0\n", encoding="utf-8")
            profile.write_text("{}", encoding="utf-8")
            return (
                ts_seed,
                intermediate_seed,
                intermediate_seed,
                tuple((int(pair[0]), int(pair[1])) for pair in forming_bonds),
                profile,
                "COMPLETE",
                "high",
                (),
            )

    monkeypatch.setattr(v4_orchestrator, "CensoLiteEngine", FakeCensoLiteEngine)
    monkeypatch.setattr(v4_orchestrator, "PEBScanner", FakePEBScanner)


def _install_fast_stage_calculators(monkeypatch) -> None:
    class FakeCalculator:
        def __init__(self, _config: dict[str, Any], _theory: dict[str, Any], event_callback=None):
            self.event_callback = event_callback

        def run_structure(self, structure: dict[str, Any], output_dir: Path) -> dict[str, Any]:
            if self.event_callback is not None:
                payload = {
                    "structure_id": structure["id"],
                    "engine": "mock",
                    "method": "mock",
                    "basis": "",
                    "solvent": None,
                    "solvent_model": None,
                }
                self.event_callback("optimization_started", payload)
                self.event_callback("optimization_finished", {**payload, "status": "complete"})
                self.event_callback("single_point_started", payload)
                self.event_callback("single_point_finished", {**payload, "status": "complete"})
            is_ts = structure.get("kind") == "ts"
            return {
                "id": structure["id"],
                "kind": structure.get("kind", "minimum"),
                "opt_xyz": str(Path(structure.get("opt_xyz") or structure["input_xyz"])),
                "opt_status": "complete",
                "sp_status": "complete",
                "frequency_status": "complete" if is_ts else "not_requested",
                "ts_frequency_valid": True if is_ts else None,
                "status": "complete",
                "usable_for_ml": True,
            }

    monkeypatch.setattr(lowlevel_module, "StageCalculator", FakeCalculator)
    monkeypatch.setattr(highlevel_module, "StageCalculator", FakeCalculator)


def test_s2_manifest_records_variant_traceability(monkeypatch, tmp_path: Path):
    variants = [
        {
            "variant_id": "product_major",
            "directory_name": "product_major",
            "branch_id": "BR_MAJOR",
            "pathway_id": "primary",
        }
    ]
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(_fake_s0_writer(variants)),
    )
    _install_s1_s2_fakes(monkeypatch)

    work_dir = tmp_path / "run"
    result = _orchestrator().run(_record(), work_dir, stop_after="s2")

    manifest = json.loads((work_dir / "S2_PEB" / "product_major" / "manifest.json").read_text(encoding="utf-8"))
    assert result["completed_through"] == "s2"
    assert manifest["variant_id"] == "product_major"
    assert manifest["parent_variant_id"] == "product_major"
    assert manifest["selected_id"] == "product_major_conf_0001"
    assert manifest["selected_xyz_ref"] == str(work_dir / "S1_ConfSearch" / "product_major" / "selected.xyz")


def test_s3_and_s4_manifests_propagate_structure_identity(monkeypatch, tmp_path: Path):
    variants = [
        {
            "variant_id": "product_major",
            "directory_name": "product_major",
            "branch_id": "BR_MAJOR",
            "pathway_id": "primary",
        }
    ]
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(_fake_s0_writer(variants)),
    )
    _install_s1_s2_fakes(monkeypatch)
    _install_fast_stage_calculators(monkeypatch)

    work_dir = tmp_path / "run"
    result = _orchestrator().run(_record(), work_dir, stop_after="s4")

    assert result["completed_through"] == "s4"

    s3_manifest = work_dir / "S3_LowLevel" / "manifest.json"
    s4_manifest = work_dir / "S4_HighLevel" / "manifest.json"
    s3_structures = {
        row["id"]: row
        for row in json.loads(s3_manifest.read_text(encoding="utf-8"))["structures"]
    }
    s4_structures = {
        row["id"]: row
        for row in json.loads(s4_manifest.read_text(encoding="utf-8"))["structures"]
    }

    assert s3_structures["precursor"]["role"] == "precursor"
    assert s3_structures["precursor"]["variant_id"] == "precursor"
    assert s3_structures["precursor"]["branch_id"] is None
    assert s3_structures["precursor"]["pathway_id"] is None
    assert s3_structures["precursor"]["parent_structure_id"] is None
    assert s3_structures["precursor"]["source_stage"] == "S1"

    assert s3_structures["product_major"]["structure_id"] == "product_major"
    assert s3_structures["product_major"]["variant_id"] == "product_major"
    assert s3_structures["product_major"]["role"] == "product"
    assert s3_structures["product_major"]["branch_id"] == "BR_MAJOR"
    assert s3_structures["product_major"]["pathway_id"] == "primary"
    assert s3_structures["product_major"]["parent_structure_id"] is None
    assert s3_structures["product_major"]["source_stage"] == "S1"

    assert s3_structures["product_major_int"]["role"] == "intermediate"
    assert s3_structures["product_major_int"]["variant_id"] == "product_major"
    assert s3_structures["product_major_int"]["branch_id"] == "BR_MAJOR"
    assert s3_structures["product_major_int"]["pathway_id"] == "primary"
    assert s3_structures["product_major_int"]["parent_structure_id"] == "product_major"
    assert s3_structures["product_major_int"]["source_stage"] == "S2"

    assert s3_structures["product_major_ts"]["role"] == "ts"
    assert s3_structures["product_major_ts"]["variant_id"] == "product_major"
    assert s3_structures["product_major_ts"]["branch_id"] == "BR_MAJOR"
    assert s3_structures["product_major_ts"]["pathway_id"] == "primary"
    assert s3_structures["product_major_ts"]["parent_structure_id"] == "product_major"
    assert s3_structures["product_major_ts"]["source_stage"] == "S2"

    assert s4_structures["product_major"]["structure_id"] == "product_major"
    assert s4_structures["product_major"]["variant_id"] == "product_major"
    assert s4_structures["product_major"]["role"] == "product"
    assert s4_structures["product_major"]["branch_id"] == "BR_MAJOR"
    assert s4_structures["product_major"]["pathway_id"] == "primary"
    assert s4_structures["product_major"]["parent_structure_id"] is None
    assert s4_structures["product_major"]["source_stage"] == "S3"
    assert s4_structures["product_major"]["source_s3"]["manifest"] == str(s3_manifest)

    assert s4_structures["product_major_ts"]["structure_id"] == "product_major_ts"
    assert s4_structures["product_major_ts"]["variant_id"] == "product_major"
    assert s4_structures["product_major_ts"]["role"] == "ts"
    assert s4_structures["product_major_ts"]["branch_id"] == "BR_MAJOR"
    assert s4_structures["product_major_ts"]["pathway_id"] == "primary"
    assert s4_structures["product_major_ts"]["parent_structure_id"] == "product_major"
    assert s4_structures["product_major_ts"]["source_stage"] == "S3"
    assert s4_structures["product_major_ts"]["source_s3"]["structure_status"] == "complete"
    assert s4_structures["product_major_ts"]["source_s3"]["opt_status"] == "complete"
    assert s4_structures["product_major_ts"]["source_s3"]["sp_status"] == "complete"
    assert s4_structures["product_major_ts"]["source_s3"]["frequency_status"] == "complete"
    assert s4_structures["product_major_ts"]["source_s3"]["ts_frequency_valid"] is True
    assert s4_structures["product_major_ts"]["source_s3"]["usable_for_ml"] is True
