import pytest
from rph_core.utils.molecule_utils import (
    canonicalize_smiles,
    get_molecular_formula,
    get_molecule_key,
    is_small_molecule
)

def test_canonicalize_smiles():
    assert canonicalize_smiles("CCO") == "CCO"
    assert canonicalize_smiles("OCC") == "CCO"
    assert canonicalize_smiles("C1=CC=CC=C1") == "c1ccccc1"
    assert canonicalize_smiles("invalid") is None

def test_get_molecular_formula():
    assert get_molecular_formula("CCO") == "C2H6O"
    assert get_molecular_formula("C1=CC=CC=C1") == "C6H6"
    assert get_molecular_formula("invalid") is None

def test_get_molecule_key():
    key = get_molecule_key("CCO")
    assert key == "C2H6O_CCO"
    assert get_molecule_key("OCC") == "C2H6O_CCO"
    assert get_molecule_key("invalid") is None

def test_is_small_molecule():
    assert is_small_molecule("CCO", threshold=10) is True
    aspirin = "CC(=O)OC1=CC=CC=C1C(=O)O"
    assert is_small_molecule(aspirin, threshold=10) is False
    assert is_small_molecule(aspirin, threshold=15) is True
    assert is_small_molecule("invalid") is False
