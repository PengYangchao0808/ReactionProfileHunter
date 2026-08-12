import json
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from rph_core.steps.conformer_search.censo_lite_runtime import CensoLiteRuntime
from rph_core.steps.mechanism_classifier.s0_artifact_writer import write_s0_artifacts
from rph_core.steps.mechanism_classifier.s0_record import load_s0_reaction_record
from rph_core.utils.atom_mapping import resolve_peb_forming_bonds


REPO_ROOT = Path(__file__).resolve().parents[1]
DATASET = REPO_ROOT / "data" / "reaxys_cleaned.csv"


def test_rx2_s0_s1_mapping_resolves_c7_c8_not_c7_c22(tmp_path: Path):
    record = load_s0_reaction_record(DATASET, "2")
    _context, variants, _manifest = write_s0_artifacts(tmp_path, record)
    major = next(variant for variant in variants if variant.variant_id == "product_major")

    runtime = CensoLiteRuntime({}, tmp_path / "S1_ConfSearch", "product_major")
    initial_xyz = runtime.embed(major.branch_product_smiles)
    result = resolve_peb_forming_bonds(
        atom_map_smiles_path=tmp_path / "S0_Mechanism" / "atom_mappings" / "product_major.json",
        smiles_to_xyz_map_path=runtime.molecule_dir / "initial_atom_mapping.json",
        product_xyz_path=initial_xyz,
    )

    assert result.forming_bonds_product_xyz_0based == ((2, 3), (7, 6))
    assert result.forming_bonds_product_xyz_1based == ((3, 4), (8, 7))
    assert (6, 21) not in result.forming_bonds_product_xyz_0based


def test_s0_variant_registry_points_to_verified_mapping_tables(tmp_path: Path):
    record = load_s0_reaction_record(DATASET, "2")
    write_s0_artifacts(tmp_path, record)
    registry = json.loads(
        (tmp_path / "S0_Mechanism" / "variant_registry.json").read_text(encoding="utf-8")
    )

    assert registry["schema_version"] == "s0_variant_registry_v2"
    for variant in registry["variants"]:
        mapping_path = tmp_path / "S0_Mechanism" / variant["atom_mapping_ref"]
        mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
        assert mapping["mapping_status"] == "verified"
        assert mapping["mapping_digest"] == variant["atom_mapping_digest"]
        assert mapping["reference_smiles"] == variant["branch_product_smiles"]
