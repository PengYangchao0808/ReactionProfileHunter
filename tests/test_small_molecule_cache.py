import pytest
from rph_core.utils.small_molecule_cache import SmallMoleculeCache

def test_cache_init(tmp_path):
    cache_root = tmp_path / "SmallMolecules"
    cache = SmallMoleculeCache(cache_root)
    assert cache.cache_root.exists()
    assert cache.cache_root.is_dir()

def test_get_path(tmp_path):
    cache_root = tmp_path / "SmallMolecules"
    cache = SmallMoleculeCache(cache_root)
    smiles = "CCO"
    path = cache.get_path(smiles)
    assert path is not None
    assert "C2H6O" in path.name
    assert "CCO" in path.name

def test_exists_empty(tmp_path):
    cache_root = tmp_path / "SmallMolecules"
    cache = SmallMoleculeCache(cache_root)
    smiles = "CCO"
    assert not cache.exists(smiles)

def test_exists_with_file(tmp_path):
    cache_root = tmp_path / "SmallMolecules"
    cache = SmallMoleculeCache(cache_root)
    smiles = "CCO"
    path = cache.get_or_create(smiles)
    (path / "molecule_min.xyz").write_text("dummy xyz")
    assert cache.exists(smiles)

def test_get_or_create(tmp_path):
    cache_root = tmp_path / "SmallMolecules"
    cache = SmallMoleculeCache(cache_root)
    smiles = "CCO"
    path = cache.get_or_create(smiles, name="Ethanol")
    assert path.exists()
    assert path.is_dir()
    
    path2 = cache.get_or_create(smiles)
    assert path == path2

def test_invalid_smiles(tmp_path):
    cache_root = tmp_path / "SmallMolecules"
    cache = SmallMoleculeCache(cache_root)
    smiles = "INVALID"
    assert cache.get_path(smiles) is None
    assert not cache.exists(smiles)
    with pytest.raises(ValueError):
        cache.get_or_create(smiles)
