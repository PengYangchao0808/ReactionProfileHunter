import csv
import hashlib
import json
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from rph_core.steps.mechanism_classifier.s0_artifact_writer import write_s0_artifacts
from rph_core.steps.mechanism_classifier.s0_record import load_s0_reaction_record


def _write_sample_csv(path: Path) -> Path:
    row = {
        "rx_id": "TEST_001",
        "precursor_smiles": "C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C",
        "product_smiles_main": "CC(C)(C)OC(=O)N1CCC[C@@]23C=C[C@@H](CC(=O)[C@H]12)O3",
        "rxn_smiles_mapped": "[CH2:1]=[C:2]=[CH:3][N:4]([CH2:5][CH2:6][CH2:7][c:8]1[cH:9][cH:10][cH:11][o:12]1)[C:13](=[O:14])[O:15][C:16]([CH3:17])([CH3:18])[CH3:19]>>O=[C:6]1[C@@H:5]2[N:4]([C:13](=[O:14])[O:15][C:16]([CH3:17])([CH3:18])[CH3:19])[CH2:3][CH2:2][CH2:1][C@@:11]23[CH:10]=[CH:9][C@@H:8]([CH2:7]1)[O:12]3",
        "reaction_type": "4+3",
        "mechanistic_forming_bonds": json.dumps([[5, 11], [7, 8]]),
        "mapping_trusted_for_mechanistic": "true",
        "mapping_confidence": "0.793",
        "topology": "INTRA_TYPE_I",
        "cyclo_mode": "[4+3]",
        "precursor_type": "allenamide",
        "resolver_provenance": json.dumps(
            {
                "furan_atoms": [8, 11],
                "allenamide_atoms": [1, 3],
                "raw_formed_count": 2,
                "mechanistic_formed_count": 2,
                "ring_validation": "intramolecular_bypass",
                "selection_method": "fallback",
                "component_contract_satisfied": False,
                "precursor_type": "allenamide",
                "precursor_subtype": "general",
                "intramolecular_resolution": True,
            }
        ),
        "dr_major": "",
        "dr_minor": "",
        "de": "",
        "ee": "",
        "solvent": "CH2Cl2",
        "temp_celsius": "-45",
        "has_lewis_acid": "False",
        "stereo_consistency": "NA",
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)
    return path


def test_write_s0_artifacts_emits_full_p1_context_bundle(tmp_path: Path):
    csv_path = _write_sample_csv(tmp_path / "sample.csv")
    record = load_s0_reaction_record(csv_path, "TEST_001")

    expected_row_hash = hashlib.sha256(
        json.dumps(record.raw_row, sort_keys=True).encode("utf-8")
    ).hexdigest()
    assert record.canonical_precursor_smiles == "C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C"
    assert record.mapped_precursor_smiles.startswith("[CH2:1]=[C:2]=[CH:3]")
    assert record.topology == "INTRA_TYPE_I"
    assert record.cyclo_mode == "[4+3]"
    assert record.precursor_type == "allenamide"
    assert record.source_row_hash == expected_row_hash
    assert record.raw_row["resolver_provenance"]

    reaction_context, product_variants, mechanism_path = write_s0_artifacts(tmp_path / "run", record)

    stage_dir = tmp_path / "run" / "S0_Mechanism"
    expected_files = {
        "mechanism.json",
        "reaction_context.json",
        "mechanism_graph.json",
        "variant_registry.json",
        "dr_branch_plan.json",
        "atom_map_smiles.json",
    }
    assert {path.name for path in stage_dir.iterdir() if path.is_file()} >= expected_files
    assert mechanism_path == stage_dir / "mechanism.json"

    assert reaction_context.reaction_id == "TEST_001"
    assert reaction_context.mapped_product_smiles == record.mapped_product_smiles
    assert reaction_context.forming_bonds_map_space == record.mapped_forming_bonds
    assert {variant.variant_id for variant in product_variants} == {
        "product_major",
        "product_minor_001",
    }

    mechanism = json.loads((stage_dir / "mechanism.json").read_text(encoding="utf-8"))
    reaction_context_json = json.loads((stage_dir / "reaction_context.json").read_text(encoding="utf-8"))
    mechanism_graph = json.loads((stage_dir / "mechanism_graph.json").read_text(encoding="utf-8"))
    variant_registry = json.loads((stage_dir / "variant_registry.json").read_text(encoding="utf-8"))
    atom_map_smiles = json.loads((stage_dir / "atom_map_smiles.json").read_text(encoding="utf-8"))

    assert mechanism["schema_version"] == "s0_mechanism_v2"
    assert mechanism["forming_bonds"] == [list(pair) for pair in record.forming_bonds]
    assert mechanism["mapped_forming_bonds"] == [list(pair) for pair in record.mapped_forming_bonds]
    assert mechanism["canonical_precursor_smiles"] == record.canonical_precursor_smiles
    assert mechanism["topology"] == "INTRA_TYPE_I"
    assert mechanism["cyclo_mode"] == "[4+3]"
    assert mechanism["variants"] == ["product_major", "product_minor_001"]

    assert reaction_context_json["schema_version"] == "s0_reaction_context_v1"
    assert reaction_context_json["reaction_id"] == "TEST_001"
    assert reaction_context_json["source_row_hash"] == expected_row_hash
    assert reaction_context_json["canonical_product_smiles"] == record.product_smiles
    assert reaction_context_json["canonical_precursor_smiles"] == record.canonical_precursor_smiles
    assert reaction_context_json["mapped_reaction_smiles"] == record.raw_row["rxn_smiles_mapped"]
    assert reaction_context_json["forming_bonds_product_smiles_idx"] == [list(pair) for pair in record.forming_bonds]
    assert reaction_context_json["mechanism_graph_ref"] == "mechanism_graph.json"
    assert reaction_context_json["branch_plan_ref"] == "dr_branch_plan.json"

    assert mechanism_graph["reaction_id"] == "TEST_001"
    assert mechanism_graph["cyclo_mode"] == "[4+3]"
    assert mechanism_graph["topology"] == "INTRA_TYPE_I"
    assert mechanism_graph["source_data"]["resolver_provenance"] == record.raw_row["resolver_provenance"]

    variant_ids = [variant["variant_id"] for variant in variant_registry["variants"]]
    assert variant_registry["schema_version"] == "s0_variant_registry_v1"
    assert variant_ids == ["product_major", "product_minor_001"]
    assert variant_registry["variants"][0]["branch_id"] == "BR_MAJOR"
    assert variant_registry["variants"][1]["branch_id"] == "BR_DR_001"

    assert atom_map_smiles["schema_version"] == "s0_atom_map_smiles_v1"
    assert atom_map_smiles["product_smiles_idx_space"] == "geometry_product_smiles_idx"
    assert atom_map_smiles["forming_bonds_map_space"] == [list(pair) for pair in record.mapped_forming_bonds]
    for map_pair, product_pair in zip(record.mapped_forming_bonds, record.forming_bonds):
        assert atom_map_smiles["map_to_product_smiles"][str(map_pair[0])] == product_pair[0]
        assert atom_map_smiles["map_to_product_smiles"][str(map_pair[1])] == product_pair[1]
