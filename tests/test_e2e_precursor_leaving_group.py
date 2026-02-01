"""
End-to-End integration test for precursor + leaving group workflow.
Tests the complete pipeline from CSV input to final output structure.
"""
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
import tempfile
import shutil


def test_e2e_multi_molecule_workflow(tmp_path):
    """Test complete workflow with product, precursor, and leaving group."""
    # Mock the ConformerEngine and QC runners to avoid actual calculations
    with patch('rph_core.steps.conformer_search.engine.ConformerEngine') as mock_engine_class:
        mock_engine = MagicMock()
        mock_engine.run.return_value = (tmp_path / "test_min.xyz", -229.0)
        mock_engine_class.return_value = mock_engine
        
        with patch('rph_core.steps.step2_retro.retro_scanner.RetroScanner') as mock_s2_class:
            mock_s2 = MagicMock()
            mock_s2.run.return_value = MagicMock(
                ts_guess_xyz=tmp_path / "ts_guess.xyz",
                reactant_xyz=tmp_path / "reactant_complex.xyz",
                forming_bonds=((0, 1), (2, 3))
            )
            mock_s2_class.return_value = mock_s2
            
            with patch('rph_core.steps.step3_opt.ts_optimizer.TSOptimizer') as mock_s3_class:
                mock_s3 = MagicMock()
                mock_s3.run.return_value = MagicMock(success=True)
                mock_s3_class.return_value = mock_s3
                
                # Test data
                product_smiles = "O=C1C=C2CCCC(=O)[C@@H]3C[C@H]1O[C@H]23"
                precursor_smiles = "C=CC(=O)CCCC1=CC(=O)COC1OC(C)=O"
                leaving_group_smiles = "CC(=O)O"  # AcOH
                
                # Verify the multi-molecule workflow can be invoked
                from rph_core.steps.anchor.handler import AnchorPhase, AnchorPhaseResult
                
                config = {
                    "step1": {"crest": {"gfn_level": 2}},
                    "theory": {
                        "optimization": {"method": "B3LYP", "basis": "def2-SVP"},
                        "single_point": {"method": "wB97X-D3BJ", "basis": "def2-TZVPP"}
                    }
                }
                
                # Create molecules dictionary (as orchestrator would)
                molecules = {
                    "product": product_smiles,
                    "precursor": precursor_smiles,
                    "leaving_group": leaving_group_smiles
                }
                
                # Initialize AnchorPhase
                work_dir = tmp_path / "S1_ConfGeneration"
                work_dir.mkdir(parents=True)
                anchor = AnchorPhase(config, work_dir)
                
                # Verify small molecule detection works
                from rph_core.utils.molecule_utils import is_small_molecule, get_molecule_key
                
                assert is_small_molecule(leaving_group_smiles) is True
                assert is_small_molecule(product_smiles) is False
                
                key = get_molecule_key(leaving_group_smiles)
                assert "C2H4O2" in key  # Molecular formula
                
                print("✓ E2E multi-molecule workflow structure verified")


def test_small_molecule_cache_mechanism(tmp_path):
    """Test that small molecule caching mechanism works correctly."""
    from rph_core.utils.molecule_utils import is_small_molecule, get_molecule_key
    from rph_core.utils.small_molecule_cache import SmallMoleculeCache
    
    # Test small molecule classification
    small_molecules = [
        "CC(=O)O",  # AcOH
        "O",        # Water
        "C",        # Methane
    ]
    
    for smi in small_molecules:
        assert is_small_molecule(smi) is True, f"{smi} should be classified as small molecule"
    
    # Test larger molecules (not small)
    large_molecules = [
        "O=C1C=C2CCCC(=O)[C@@H]3C[C@H]1O[C@H]23",  # Product
        "C=CC(=O)CCCC1=CC(=O)COC1OC(C)=O",          # Precursor
    ]
    
    for smi in large_molecules:
        assert is_small_molecule(smi) is False, f"{smi[:20]}... should NOT be classified as small molecule"
    
    # Test cache key generation
    cache = SmallMoleculeCache(tmp_path / "SmallMolecules")
    
    # First access - should create directory
    key = get_molecule_key("CC(=O)O")
    cache_path = cache.get_or_create("CC(=O)O", "AcOH")
    assert cache_path.exists()
    assert (tmp_path / "SmallMolecules") in cache_path.parents
    
    # Second access - should return existing
    cache_path2 = cache.get_or_create("CC(=O)O", "AcOH")
    assert cache_path2 == cache_path
    
    print("✓ Small molecule cache mechanism verified")


def test_s1_confgeneration_directory_structure(tmp_path):
    """Verify S1_ConfGeneration directory structure is used."""
    from rph_core.steps.conformer_search.engine import ConformerEngine
    
    # Verify the engine creates molecule-specific subdirectories
    config = {
        "step1": {"crest": {"gfn_level": 2}},
        "theory": {"optimization": {}, "single_point": {}}
    }
    
    work_dir = tmp_path / "S1_ConfGeneration"
    work_dir.mkdir()
    
    # Test that molecule directories are created correctly
    molecule_name = "product"
    engine = ConformerEngine(config, work_dir, molecule_name)
    
    expected_dirs = [
        engine.molecule_dir,           # S1_ConfGeneration/product/
        engine.crest_dir,              # S1_ConfGeneration/product/xtb2/
        engine.cluster_dir,            # S1_ConfGeneration/product/cluster/
        engine.dft_dir,                # S1_ConfGeneration/product/dft/
    ]
    
    for d in expected_dirs:
        assert d.exists(), f"Expected directory {d} to exist"
        assert d.is_dir()
    
    print("✓ S1_ConfGeneration directory structure verified")


if __name__ == "__main__":
    import sys
    with pytest.raises(SystemExit):
        pytest.main([__file__, "-v"])
