import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
from rph_core.steps.step2_retro.retro_scanner import RetroScanner

@pytest.fixture
def mock_scanner():
    config = {'step2': {'xtb_settings': {'solvent': 'acetone'}}}
    scanner = RetroScanner(config)
    scanner.xtb_runner = MagicMock()
    # Mock optimize to return success
    scanner.xtb_runner.optimize.return_value = MagicMock(success=True, output_file=None)
    scanner.smarts_matcher = MagicMock()
    # Mock SMARTS match result
    match_result = MagicMock()
    match_result.matched = True
    match_result.bond_1.atom_idx_1 = 0
    match_result.bond_1.atom_idx_2 = 1
    match_result.bond_1.current_length = 1.5
    match_result.bond_2.atom_idx_1 = 0
    match_result.bond_2.atom_idx_2 = 1
    match_result.bond_2.current_length = 1.5
    match_result.pattern_name = "test"
    match_result.confidence = 1.0
    scanner.smarts_matcher.find_reactive_bonds.return_value = match_result
    
    scanner.bond_stretcher = MagicMock()
    scanner.bond_stretcher.stretch_two_bonds.return_value = [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]]
    
    return scanner

def test_s2_finds_product_min_directly(tmp_path, mock_scanner):
    """验证 S2 能在 S1_ConfGeneration 目录中直接找到 product_min.xyz (v6.1 扁平结构)"""
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir()
    product_file = s1_dir / "product_min.xyz"
    product_file.write_text("2\nTest\nC 0 0 0\nC 0 0 1.5")
    
    output_dir = tmp_path / "S2_Retro"
    
    with patch("rph_core.steps.step2_retro.retro_scanner.read_xyz") as mock_read, \
         patch("rph_core.steps.step2_retro.retro_scanner.LogParser.extract_last_converged_coords") as mock_parser:
        mock_parser.return_value = (None, None, "error") # Force fallback
        mock_read.return_value = ([[0.0, 0.0, 0.0], [0.0, 0.0, 1.5]], ["C", "C"])
        
        mock_scanner.run(s1_dir, output_dir)
        
    mock_scanner.smarts_matcher.find_reactive_bonds.assert_called_with(product_file)

def test_s2_finds_molecule_min_directly(tmp_path, mock_scanner):
    """验证 S2 能在 S1_ConfGeneration 目录中直接找到 [mol_name]_min.xyz"""
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir()
    product_file = s1_dir / "my_mol_min.xyz"
    product_file.write_text("2\nTest\nC 0 0 0\nC 0 0 1.5")
    
    output_dir = tmp_path / "S2_Retro"
    
    with patch("rph_core.steps.step2_retro.retro_scanner.read_xyz") as mock_read, \
         patch("rph_core.steps.step2_retro.retro_scanner.LogParser.extract_last_converged_coords") as mock_parser:
        mock_parser.return_value = (None, None, "error")
        mock_read.return_value = ([[0.0, 0.0, 0.0], [0.0, 0.0, 1.5]], ["C", "C"])
        
        mock_scanner.run(s1_dir, output_dir, molecule_name="my_mol")
        
    mock_scanner.smarts_matcher.find_reactive_bonds.assert_called_with(product_file)

def test_s2_finds_product_in_subdir(tmp_path, mock_scanner):
    """验证 S2 能在 S1_ConfGeneration/product/ 目录中找到 product_min.xyz (v3.0 子目录结构)"""
    s1_dir = tmp_path / "S1_ConfGeneration"
    mol_dir = s1_dir / "product"
    mol_dir.mkdir(parents=True)
    product_file = mol_dir / "product_min.xyz"
    product_file.write_text("2\nTest\nC 0 0 0\nC 0 0 1.5")
    
    output_dir = tmp_path / "S2_Retro"
    
    with patch("rph_core.steps.step2_retro.retro_scanner.read_xyz") as mock_read, \
         patch("rph_core.steps.step2_retro.retro_scanner.LogParser.extract_last_converged_coords") as mock_parser:
        mock_parser.return_value = (None, None, "error")
        mock_read.return_value = ([[0.0, 0.0, 0.0], [0.0, 0.0, 1.5]], ["C", "C"])
        
        mock_scanner.run(s1_dir, output_dir)
        
    mock_scanner.smarts_matcher.find_reactive_bonds.assert_called_with(product_file)

def test_s2_auto_detect_subdir(tmp_path, mock_scanner):
    """验证 S2 能自动检测唯一的分子子目录"""
    s1_dir = tmp_path / "S1_ConfGeneration"
    mol_dir = s1_dir / "unique_molecule"
    mol_dir.mkdir(parents=True)
    product_file = mol_dir / "global_min.xyz"
    product_file.write_text("2\nTest\nC 0 0 0\nC 0 0 1.5")
    
    output_dir = tmp_path / "S2_Retro"
    
    with patch("rph_core.steps.step2_retro.retro_scanner.read_xyz") as mock_read, \
         patch("rph_core.steps.step2_retro.retro_scanner.LogParser.extract_last_converged_coords") as mock_parser:
        mock_parser.return_value = (None, None, "error")
        mock_read.return_value = ([[0.0, 0.0, 0.0], [0.0, 0.0, 1.5]], ["C", "C"])
        
        mock_scanner.run(s1_dir, output_dir)
        
    mock_scanner.smarts_matcher.find_reactive_bonds.assert_called_with(product_file)

def test_s2_backward_compat_file_path(tmp_path, mock_scanner):
    """验证 S2 仍支持直接传入文件路径 (v2.1 结构)"""
    product_file = tmp_path / "legacy_product.xyz"
    product_file.write_text("2\nTest\nC 0 0 0\nC 0 0 1.5")
    
    output_dir = tmp_path / "S2_Retro"
    
    with patch("rph_core.steps.step2_retro.retro_scanner.read_xyz") as mock_read, \
         patch("rph_core.steps.step2_retro.retro_scanner.LogParser.extract_last_converged_coords") as mock_parser:
        mock_parser.return_value = (None, None, "error")
        mock_read.return_value = ([[0.0, 0.0, 0.0], [0.0, 0.0, 1.5]], ["C", "C"])
        
        mock_scanner.run(product_file, output_dir)
        
    mock_scanner.smarts_matcher.find_reactive_bonds.assert_called_with(product_file)

def test_s2_fail_if_no_file(tmp_path, mock_scanner):
    """验证当找不到产物文件时抛出 RuntimeError"""
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir()
    
    output_dir = tmp_path / "S2_Retro"
    
    with pytest.raises(RuntimeError, match="未找到产物文件"):
        mock_scanner.run(s1_dir, output_dir)
