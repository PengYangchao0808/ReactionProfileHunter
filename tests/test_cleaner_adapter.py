from pathlib import Path

import numpy as np

from rph_core.utils.cleaner_adapter import (
    convert_cleaner_row_to_record,
    map_pairs_to_xyz_indices,
)
from rph_core.utils.file_io import write_xyz


def test_convert_cleaner_row_sets_explicit_zero_based_metadata() -> None:
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
        reaction_profiles={"[4+3]_default": {"s2_strategy": "forward_scan"}},
    )

    assert record is not None
    assert record.raw.get("formed_bond_index_pairs") == "0-1"
    assert record.raw.get("forming_bonds") == "0-1"
    assert record.raw.get("forming_bonds_index_base") == "0"
    assert record.raw.get("index_base") == "0"
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
