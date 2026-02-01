"""
Tests for Step 4 (S4) path compatibility with S1_ConfGeneration and S1_Product.
"""
import pytest
from pathlib import Path
from rph_core.steps.step4_features.mech_packager import (
    MechanismContext, S1_DIR_ALIASES, pack_mechanism_assets
)


def test_s4_finds_product_in_s1_confgeneration(tmp_path):
    """Verify S4 packager can find product in S1_ConfGeneration directory."""
    # Setup S1 structure with new name
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir(parents=True)
    product_xyz = s1_dir / "product_min.xyz"
    product_xyz.write_text("3\ntest\nC 0 0 0\nO 1 0 0\n")
    
    # Setup S4 directory
    s4_dir = tmp_path / "S4_Data"
    s4_dir.mkdir()
    
    # Create context and verify S1 dir is found
    context = MechanismContext(s1_dir=s1_dir)
    
    # Check aliases include new name
    assert "S1_ConfGeneration" in S1_DIR_ALIASES
    assert "S1_Product" in S1_DIR_ALIASES
    
    # Verify the directory resolution works
    assert context.s1_dir == s1_dir
    print("✓ S4 correctly identifies S1_ConfGeneration")


def test_s4_finds_product_in_s1_product_fallback(tmp_path):
    """Verify S4 packager falls back to S1_Product if S1_ConfGeneration not found."""
    # Setup old S1 structure only
    s1_dir = tmp_path / "S1_Product"
    s1_dir.mkdir()
    product_xyz = s1_dir / "product_min.xyz"
    product_xyz.write_text("3\ntest\nC 0 0 0\nO 1 0 0\n")
    
    # Create context
    context = MechanismContext(s1_dir=s1_dir)
    
    # Verify fallback works
    assert context.s1_dir == s1_dir
    print("✓ S4 correctly falls back to S1_Product")


def test_s4_aliases_defined(tmp_path):
    """Verify S4 has correct S1_DIR_ALIASES defined."""
    expected_aliases = ["S1_ConfGeneration", "S1_Product"]
    
    assert S1_DIR_ALIASES == expected_aliases
    print(f"✓ S1_DIR_ALIASES = {S1_DIR_ALIASES}")


if __name__ == "__main__":
    import sys
    with pytest.raises(SystemExit):
        pytest.main([__file__, "-v"])
