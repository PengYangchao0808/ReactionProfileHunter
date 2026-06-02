from pathlib import Path

from rph_core.utils.dataset_loader import load_reaction_records


def test_dataset_loader_preserves_mapped_product_smiles_without_bond_conversion(tmp_path: Path) -> None:
    dataset_file = tmp_path / "dataset.csv"
    _ = dataset_file.write_text(
        "rx_id,precursor_smiles,product_smiles_main,core_bond_changes,mapped_precursor_smiles,mapped_product_smiles,map_status,map_confidence\nrx_product_priority,CCCC,CCCC,1-4:formed,[CH3:1][CH2:2][CH2:3][CH3:4],[CH3:2][CH2:4][CH2:1][CH3:3],OK,0.98\n",
        encoding="utf-8",
    )

    records = load_reaction_records(
        dataset_cfg={
            "path": str(dataset_file),
            "delimiter": ",",
            "id_col": "rx_id",
            "precursor_smiles_col": "precursor_smiles",
            "product_smiles_col": "product_smiles_main",
        },
    )

    assert len(records) == 1
    assert records[0].raw.get("mapped_product_smiles") == "[CH3:2][CH2:4][CH2:1][CH3:3]"
    assert "formed_bond_map_pairs" not in records[0].raw
    assert "formed_bond_index_pairs" not in records[0].raw
    assert "forming_bonds" not in records[0].raw
