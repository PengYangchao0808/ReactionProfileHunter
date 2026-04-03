from rph_core.utils.dataset_loader import _enrich_cleaner_metadata
from rph_core.utils.tsv_dataset import ReactionRecord


def test_enrich_cleaner_metadata_prefers_product_smiles_for_map_index_conversion() -> None:
    record = ReactionRecord(
        rx_id="rx_product_priority",
        precursor_smiles="CCCC",
        product_smiles_main="CCCC",
        raw={
            "core_bond_changes": "1-4:formed",
            "mapped_precursor_smiles": "[CH3:1][CH2:2][CH2:3][CH3:4]",
            "mapped_product_smiles": "[CH3:2][CH2:4][CH2:1][CH3:3]",
            "map_status": "OK",
            "map_confidence": "0.98",
        },
    )

    _enrich_cleaner_metadata(record, reaction_profiles={})

    assert record.raw.get("formed_bond_index_pairs") == "1-2"
    assert record.raw.get("formed_bond_index_pairs") != "0-3"
    assert record.raw.get("forming_bonds") == "1-2"
    assert record.raw.get("mapped_product_smiles") == "[CH3:2][CH2:4][CH2:1][CH3:3]"
