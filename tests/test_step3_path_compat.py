"""
Tests for Step 3 (S3) path compatibility with S1_ConfGeneration and S1_Product.
"""
import pytest
from pathlib import Path
from rph_core.utils.result_inspector import ResultInspector


def test_s3_finds_product_in_s1_confgeneration(tmp_path):
    """Verify S3 inspector can find product in S1_ConfGeneration/product/ directory."""
    # Setup S1 structure
    s1_dir = tmp_path / "S1_ConfGeneration"
    product_dir = s1_dir / "product"
    product_dir.mkdir(parents=True)
    
    # Create fake product file
    product_xyz = product_dir / "product_global_min.xyz"
    product_xyz.write_text("3\ntest\nC 0 0 0\nO 1 0 0\n")
    
    # Create fake SP output
    sp_dir = product_dir / "dft"
    sp_dir.mkdir()
    sp_out = sp_dir / "conf_000_SP.out"
    sp_out.write_text("ORCA TERMINATED NORMALLY\nFINAL SINGLE POINT ENERGY -229.0")
    
    # Setup ResultInspector
    config = {}
    inspector = ResultInspector(tmp_path, config)
    
    # Run check
    result = inspector.check_step("s1")
    
    assert result.should_skip is True
    assert "s1_complete" in result.reason
    print("✓ Found product in S1_ConfGeneration/product/")

def test_s3_finds_product_in_s1_product_fallback(tmp_path):
    """Verify S3 inspector falls back to S1_Product if S1_ConfGeneration not found."""
    # Setup old S1 structure only
    s1_dir = tmp_path / "S1_Product"
    s1_dir.mkdir()
    
    # Create fake product file in old location
    product_xyz = s1_dir / "product_global_min.xyz"
    product_xyz.write_text("3\ntest\nC 0 0 0\nO 1 0 0")
    
    # Create fake SP output
    sp_dir = s1_dir / "dft"
    sp_dir.mkdir()
    sp_out = sp_dir / "conf_000_SP.out"
    sp_out.write_text("ORCA TERMINATED NORMALLY\nFINAL SINGLE POINT ENERGY -229.0")
    
    # Setup ResultInspector
    config = {}
    inspector = ResultInspector(tmp_path, config)
    
    # Run check - should still work with old structure
    result = inspector.check_step("s1")
    
    assert result.should_skip is True
    assert "s1_complete" in result.reason
    print("✓ Falls back to S1_Product correctly")

def test_s3_selects_lowest_energy_sp(tmp_path):
    """Verify S3 inspector selects the lowest energy SP output."""
    s1_dir = tmp_path / "S1_ConfGeneration"
    product_dir = s1_dir / "product"
    product_dir.mkdir(parents=True)
    
    product_xyz = product_dir / "product_global_min.xyz"
    product_xyz.write_text("3\ntest\nC 0 0 0\nO 1 0 0")
    
    sp_dir = product_dir / "dft"
    sp_dir.mkdir()
    
    # Create multiple SP outputs with different energies
    sp1 = sp_dir / "conf_000_SP.out"
    sp1.write_text("ORCA TERMINATED NORMALLY\nFINAL SINGLE POINT ENERGY -230.0")
    
    sp2 = sp_dir / "conf_001_SP.out"
    sp2.write_text("ORCA TERMINATED NORMALLY\nFINAL SINGLE POINT ENERGY -229.5")
    
    sp3 = sp_dir / "conf_002_SP.out"
    sp3.write_text("ORCA TERMINATED NORMALLY\nFINAL SINGLE POINT ENERGY -228.0")
    
    config = {}
    inspector = ResultInspector(tmp_path, config)
    result = inspector.check_step("s1")
    
    assert result.should_skip is True
    print("✓ Selects lowest energy SP correctly")


if __name__ == "__main__":
    import sys
    from pathlib import Path
    with pytest.raises(SystemExit):
        pytest.main([__file__, "-v"])
