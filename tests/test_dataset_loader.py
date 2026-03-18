from rph_core.utils.dataset_loader import load_reaction_records


def test_dataset_loader_enriches_formed_bond_index_pairs(tmp_path) -> None:
    dataset_file = tmp_path / "dataset.csv"
    dataset_file.write_text(
        "rx_id,precursor_smiles,product_smiles_main,core_bond_changes,rxn_smiles_mapped,reaction_type,map_status\n"
        "9425282,C=C,CC,1-2:formed,[CH2:1]=[CH2:2]>>[CH3:1][CH3:2],4+3,OK\n",
        encoding="utf-8",
    )

    records = load_reaction_records(
        dataset_cfg={
            "path": str(dataset_file),
            "delimiter": ",",
            "id_col": "rx_id",
            "precursor_smiles_col": "precursor_smiles",
            "product_smiles_col": "product_smiles_main",
            "reaction_profiles": {"[4+3]_default": {"s2_strategy": "forward_scan"}},
        },
    )

    assert len(records) == 1
    assert records[0].raw.get("formed_bond_index_pairs") == "0-1"
    assert records[0].raw.get("forming_bonds") == "0-1"
    assert records[0].raw.get("forming_bonds_index_base") == "0"
    assert records[0].raw.get("index_base") == "0"
    assert records[0].raw.get("mapped_product_smiles") == "[CH3:1][CH3:2]"
    assert records[0].raw.get("reaction_profile") == "[4+3]_default"


def test_dataset_loader_supports_prefixed_filter_ids(tmp_path) -> None:
    dataset_file = tmp_path / "dataset.csv"
    dataset_file.write_text(
        "rx_id,precursor_smiles,product_smiles_main\n"
        "9422028,C=C,CC\n"
        "9425282,C=CC,CCC\n",
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
        filter_ids=["rx_9422028", "9422028"],
    )

    assert len(records) == 1
    assert records[0].rx_id == "9422028"
