from pathlib import Path

import numpy as np

from rph_core.utils.cleaner_adapter import (
    convert_cleaner_row_to_record,
    map_pairs_to_xyz_indices,
)
from rph_core.utils.file_io import write_xyz


def test_convert_cleaner_row_preserves_loader_metadata_without_bond_derivations() -> None:
    row = {
        "rx_id": "rx_1",
        "precursor_smiles": "C=C",
        "mapped_precursor_smiles": "[CH2:1]=[CH2:2]",
        "core_bond_changes": "1-2:formed",
        "map_status": "OK",
        "map_confidence": "0.95",
        "reaction_type": "4+3",
    }

    record = convert_cleaner_row_to_record(
        row=row,
        row_index=1,
        reaction_profiles={"[4+3]_default": {"s2_strategy": "retro_scan"}},
    )

    assert record is not None
    assert record.raw.get("map_status") == "OK"
    assert record.raw.get("map_confidence") == "0.95"
    assert record.raw.get("mapped_precursor_smiles") == "[CH2:1]=[CH2:2]"
    assert "formed_bond_index_pairs" not in record.raw
    assert "forming_bonds" not in record.raw
    assert "forming_bonds_index_base" not in record.raw
    assert "index_base" not in record.raw
    assert record.raw.get("reaction_profile") == "[4+3]_default"


def test_map_pairs_to_xyz_indices_maps_atom_map_pairs(tmp_path: Path) -> None:
    xyz_file = tmp_path / "product.xyz"
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.4, 0.0, 0.0],
            [2.8, 0.0, 0.0],
            [4.2, 0.0, 0.0],
        ]
    )
    write_xyz(xyz_file, coords, ["C", "C", "C", "C"], title="product")

    pairs = map_pairs_to_xyz_indices(
        mapped_smiles="[CH3:1][CH2:2][CH2:3][CH3:4]",
        map_pairs=[(1, 2), (3, 4)],
        xyz_file=xyz_file,
    )

    assert pairs == [(0, 1), (2, 3)]
