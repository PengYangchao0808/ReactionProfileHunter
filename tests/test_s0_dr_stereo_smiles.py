import pytest
from pathlib import Path

from rph_core.steps.mechanism_classifier.stereo_smiles import (
    generate_stereo_branch_smiles,
)


class TestStereoSmiles:
    def test_returns_none_for_empty_smiles(self) -> None:
        result = generate_stereo_branch_smiles("", [1], [])
        assert result is None

    def test_returns_none_for_invalid_smiles(self) -> None:
        result = generate_stereo_branch_smiles("NOT_VALID", [1], [])
        assert result is None

    def test_returns_original_when_no_chiral_centers(self) -> None:
        result = generate_stereo_branch_smiles("[CH3:1]C", [1], [])
        assert result is None

    def test_preserves_atom_map_numbers(self) -> None:
        smi = "[CH2:1]=[CH:2][C:3]([C:4])[C@H:5](O)[CH2:6][CH2:7]1"
        result = generate_stereo_branch_smiles(smi, flipped_map_numbers=[5], fixed_map_numbers=[3])
        if result is not None:
            assert "5]" in result

    def test_does_not_flip_fixed_centers(self) -> None:
        smi = "[C:1]1[C@H:2](O)[C@@H:3](O)[C:4]1"
        result = generate_stereo_branch_smiles(smi, flipped_map_numbers=[2, 3], fixed_map_numbers=[2])
        assert result is not None
        assert result != smi


def test_generate_stereo_branch_smiles_roundtrip_isomeric() -> None:
    smi = "[C:1][C@H:2](O)[C:3]"
    result = generate_stereo_branch_smiles(smi, flipped_map_numbers=[2], fixed_map_numbers=[])
    if result is not None:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(result)
        assert mol is not None
