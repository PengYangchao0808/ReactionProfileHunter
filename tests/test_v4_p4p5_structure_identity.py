from __future__ import annotations

import json
from collections.abc import Sequence
from importlib import import_module
from pathlib import Path
from typing import Any

from rph_core.steps.refinement import RefinementEngine
from rph_core.steps.refinement.manifest_io import write_refinement_manifest

v4_orchestrator = import_module("rph_core.v4_orchestrator")
S0ReactionRecord = import_module("rph_core.steps.mechanism_classifier.s0_record").S0ReactionRecord


def _record(rx_id: str = "RXN_P4P5") -> Any:
    return S0ReactionRecord(
        rx_id=rx_id,
        product_smiles="C1CCC1",
        reaction_type="4+3",
        mapped_product_smiles="[CH2:1]1[CH2:2][CH2:3][CH2:4]1",
        mapped_forming_bonds=((1, 2), (3, 4)),
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
        mapping_dir = stage_dir / "atom_mappings"
        mapping_dir.mkdir(parents=True, exist_ok=True)
        registry_variants = []
        for variant in variants:
            variant_id = str(variant["variant_id"])
            mapping_ref = f"atom_mappings/{variant_id}.json"
            (stage_dir / mapping_ref).write_text(json.dumps({
                "schema_version": "s0_atom_map_smiles_v2",
                "product_smiles_idx_space": "geometry_product_smiles_idx",
                "forming_bonds_map_space": [[1, 2], [3, 4]],
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
            }), encoding="utf-8")
            registry_variants.append({**variant, "atom_mapping_ref": mapping_ref})
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
                    "variants": registry_variants,
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
            candidate_xyz.write_text(
                f"4\n{smiles}\nC 0 0 0\nC 1 0 0\nC 1 1 0\nC 0 1 0\n",
                encoding="utf-8",
            )
            selected_xyz = molecule_dir / "selected.xyz"
            selected_xyz.write_text(candidate_xyz.read_text(encoding="utf-8"), encoding="utf-8")
            manifest = molecule_dir / "manifest.json"
            payload = {
                "selected": f"{self.molecule_name}_conf_0001",
                "selected_xyz": "selected.xyz",
                "atom_mapping_ref": "atom_mapping.json",
                "candidates": [
                    {
                        "id": f"{self.molecule_name}_conf_0001",
                        "xyz": f"{self.molecule_name}/candidates/conf_0001.xyz",
                    }
                ],
            }
            (molecule_dir / "atom_mapping.json").write_text(json.dumps({
                "schema_version": "s1_atom_mapping_v1",
                "mapping_source": "rdkit_generation_sidecar",
                "confidence": "high",
                "product_smiles_idx_space": "geometry_product_smiles_idx",
                "atoms": [
                    {"smiles_idx": i, "xyz_idx": i, "element": "C", "type": "organic_heavy"}
                    for i in range(4)
                ],
            }), encoding="utf-8")
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
            contents = Path(product_xyz).read_text(encoding="utf-8")
            ts_seed.write_text(contents, encoding="utf-8")
            intermediate_seed.write_text(contents, encoding="utf-8")
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


def _install_fast_refinement_engine(monkeypatch) -> None:
    class FakeRefinementEngine(RefinementEngine):
        def run(self, structures, output_dir: Path) -> Path:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            rows = []
            for structure in structures:
                is_ts = structure.get("kind") == "ts"
                rows.append(
                    {
                        "id": structure["id"],
                        "structure_id": structure.get("structure_id", structure["id"]),
                        "variant_id": structure.get("variant_id"),
                        "role": structure.get("role"),
                        "kind": structure.get("kind", "minimum"),
                        "branch_id": structure.get("branch_id"),
                        "pathway_id": structure.get("pathway_id"),
                        "parent_structure_id": structure.get("parent_structure_id"),
                        "source_stage": structure.get("source_stage"),
                        "input_xyz": structure.get("input_xyz"),
                        "opt_xyz": structure.get("opt_xyz") or structure.get("input_xyz"),
                        "forming_bonds": structure.get("forming_bonds"),
                        "opt_status": "complete",
                        "sp_status": "complete",
                        "frequency_status": "complete" if is_ts else "complete",
                        "status": "complete",
                        "usable_for_ml": True,
                    }
                )
            return write_refinement_manifest(
                output_dir / "manifest.json",
                stage=self.profile.stage,
                fidelity=self.profile.fidelity,
                profile_id=self.profile.profile_id,
                structures=rows,
                run_id=self.run_id,
                extra={"summary": {"complete": len(rows)}},
            )

    monkeypatch.setattr(v4_orchestrator, "RefinementEngine", FakeRefinementEngine)


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
    _install_fast_refinement_engine(monkeypatch)

    work_dir = tmp_path / "run"
    result = _orchestrator().run(_record(), work_dir, stop_after="s4")

    assert result["completed_through"] == "s4"

    s3_manifest = work_dir / "S3_LowLevel" / "manifest.json"
    s4_manifest = work_dir / "S4_HighLevel" / "manifest.json"
    assert json.loads(s3_manifest.read_text(encoding="utf-8"))["schema_version"] == "refinement_manifest_v1"
    assert json.loads(s4_manifest.read_text(encoding="utf-8"))["schema_version"] == "refinement_manifest_v1"
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

    assert s4_structures["product_major_ts"]["structure_id"] == "product_major_ts"
    assert s4_structures["product_major_ts"]["variant_id"] == "product_major"
    assert s4_structures["product_major_ts"]["role"] == "ts"
    assert s4_structures["product_major_ts"]["branch_id"] == "BR_MAJOR"
    assert s4_structures["product_major_ts"]["pathway_id"] == "primary"
    assert s4_structures["product_major_ts"]["parent_structure_id"] == "product_major"
    assert s4_structures["product_major_ts"]["source_stage"] == "S3"
