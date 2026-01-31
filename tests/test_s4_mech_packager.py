"""
Test Mechanism Asset Packager (M1)
====================================

Tests for mechanism asset packaging functionality:
- Dipole source priority
- Fixed naming contract
- Missing file degradation
- Toxic path safety
- Quality flags

Author: QC Descriptors Team
Date: 2026-01-21
"""

import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

from rph_core.utils.file_io import write_xyz, read_xyz
from rph_core.steps.step4_features.mech_packager import (
    pack_mechanism_assets,
    _resolve_dipole_source,
    _copy_or_link_asset,
    _check_forming_bond_window,
    _validate_atom_count
)


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def sample_coords():
    return np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.866, 0.0]
    ])


@pytest.fixture
def sample_symbols():
    return ["C", "H", "H"]


@pytest.fixture
def step_dirs_with_all_sources(tmp_path, sample_coords, sample_symbols):
    """Fixture providing complete S1/S2/S3 directory structure."""
    s1 = tmp_path / "S1_ConfGeneration"
    s2 = tmp_path / "S2_Retro"
    s3 = tmp_path / "S3_TS"

    s1.mkdir(parents=True)
    s2.mkdir(parents=True)
    s3.mkdir(parents=True)

    # S1: product_min.xyz
    write_xyz(
        s1 / "product_min.xyz",
        sample_coords,
        sample_symbols,
        title="Product Min"
    )

    # S2: ts_guess.xyz and reactant_complex.xyz
    write_xyz(
        s2 / "ts_guess.xyz",
        sample_coords,
        sample_symbols,
        title="TS Guess"
    )
    write_xyz(
        s2 / "reactant_complex.xyz",
        sample_coords,
        sample_symbols,
        title="Reactant Complex"
    )

    # S3: ts_final.xyz and reactant_sp.xyz
    write_xyz(
        s3 / "ts_final.xyz",
        sample_coords,
        sample_symbols,
        title="TS Final"
    )
    write_xyz(
        s3 / "reactant_sp.xyz",
        sample_coords,
        sample_symbols,
        title="Reactant SP"
    )

    return {"S1": s1, "S2": s2, "S3": s3}


@pytest.fixture
def step_dirs_s3_only(tmp_path, sample_coords, sample_symbols):
    """Fixture with only S3 reactant available."""
    s1 = tmp_path / "S1_ConfGeneration"
    s2 = tmp_path / "S2_Retro"
    s3 = tmp_path / "S3_TS"

    s1.mkdir(parents=True)
    s2.mkdir(parents=True)
    s3.mkdir(parents=True)

    write_xyz(s1 / "product_min.xyz", sample_coords, sample_symbols, title="Product")
    write_xyz(s2 / "ts_guess.xyz", sample_coords, sample_symbols, title="TS Guess")

    write_xyz(s3 / "ts_final.xyz", sample_coords, sample_symbols, title="TS Final")
    write_xyz(s3 / "reactant_sp.xyz", sample_coords, sample_symbols, title="Reactant SP")

    return {"S1": s1, "S2": s2, "S3": s3}


@pytest.fixture
def step_dirs_s2_only(tmp_path, sample_coords, sample_symbols):
    """Fixture with only S2 reactant_complex available."""
    s1 = tmp_path / "S1_ConfGeneration"
    s2 = tmp_path / "S2_Retro"
    s3 = tmp_path / "S3_TS"

    s1.mkdir(parents=True)
    s2.mkdir(parents=True)
    s3.mkdir(parents=True)

    write_xyz(s1 / "product_min.xyz", sample_coords, sample_symbols, title="Product")
    write_xyz(s3 / "ts_final.xyz", sample_coords, sample_symbols, title="TS Final")

    write_xyz(s2 / "reactant_complex.xyz", sample_coords, sample_symbols, title="Reactant Complex")

    return {"S1": s1, "S2": s2, "S3": s3}


# ============================================================================
# Test Class 1: TestDipolePriority
# ============================================================================

class TestDipolePriority:
    """Test dipole source priority logic (S3 → S2 fallback)."""

    def test_dipole_priority_both_sources(
        self,
        step_dirs_with_all_sources,
        tmp_path,
        sample_coords,
        sample_symbols
    ):
        """S3 reactant should be prioritized over S2 reactant."""
        config = {
            'enabled': True,
            'copy_mode': 'copy',
            'dipole_source_priority': ['S3_reactant', 'S2_reactant_complex']
        }

        s3_reactant = step_dirs_with_all_sources['S3'] / 'reactant_sp.xyz'
        s2_reactant = step_dirs_with_all_sources['S2'] / 'reactant_complex.xyz'

        resolved_path, source_label = _resolve_dipole_source(
            s3_reactant, s2_reactant, config['dipole_source_priority']
        )

        assert source_label == 'S3_reactant'
        assert resolved_path == s3_reactant

    def test_dipole_priority_s3_only(self, step_dirs_s3_only):
        """Only S3 available - should use S3."""
        s3_reactant = step_dirs_s3_only['S3'] / 'reactant_sp.xyz'
        s2_reactant = step_dirs_s3_only['S2'] / 'reactant_complex.xyz'

        resolved_path, source_label = _resolve_dipole_source(
            s3_reactant, s2_reactant, ['S3_reactant', 'S2_reactant_complex']
        )

        assert source_label == 'S3_reactant'
        assert resolved_path == s3_reactant

    def test_dipole_priority_s2_fallback(self, step_dirs_s2_only):
        """S3 unavailable - should fallback to S2."""
        s3_reactant = step_dirs_s2_only['S3'] / 'reactant_sp.xyz'
        s2_reactant = step_dirs_s2_only['S2'] / 'reactant_complex.xyz'

        resolved_path, source_label = _resolve_dipole_source(
            s3_reactant, s2_reactant, ['S3_reactant', 'S2_reactant_complex']
        )

        assert source_label == 'S2_reactant_complex'
        assert resolved_path == s2_reactant

    def test_dipole_priority_none(self, tmp_path):
        """Both S3 and S2 unavailable - should return none."""
        s1 = tmp_path / "S1"
        s3 = tmp_path / "S3"

        s1.mkdir(parents=True)
        s3.mkdir(parents=True)

        resolved_path, source_label = _resolve_dipole_source(
            None, None, ['S3_reactant', 'S2_reactant_complex']
        )

        assert source_label == 'none'
        assert resolved_path is None

    def test_dipole_priority_end_to_end(
        self,
        step_dirs_with_all_sources,
        tmp_path
    ):
        """Verify dipole source in full mech_index.json."""
        config = {
            'enabled': True,
            'copy_mode': 'copy',
            'dipole_source_priority': ['S3_reactant', 'S2_reactant_complex']
        }

        pack_mechanism_assets(
            step_dirs=step_dirs_with_all_sources,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=((0, 1), (2, 3))
        )

        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        assert mech_index_path.exists()

        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        dipole_asset = mech_index['assets']['mech_step2_reactant_dipole']
        assert dipole_asset is not None
        assert dipole_asset['source_label'] == 'S3_reactant'


# ============================================================================
# Test Class 1b: TestDipolePriorityReversal (P0-2)
# ============================================================================

class TestDipolePriorityReversal:
    """P0-2: Test that dipole source priority order is actually respected."""

    def test_dipole_priority_reversed_order(self, step_dirs_with_all_sources):
        """When priority is [S2, S3], S2 should be chosen even if S3 exists."""
        s3_reactant = step_dirs_with_all_sources['S3'] / 'reactant_sp.xyz'
        s2_reactant = step_dirs_with_all_sources['S2'] / 'reactant_complex.xyz'

        # Priority reversed: S2 first, then S3
        resolved_path, source_label = _resolve_dipole_source(
            s3_reactant, s2_reactant, ['S2_reactant_complex', 'S3_reactant']
        )

        assert source_label == 'S2_reactant_complex'
        assert resolved_path == s2_reactant

    def test_dipole_priority_reversed_end_to_end(
        self,
        step_dirs_with_all_sources,
        tmp_path
    ):
        """End-to-end: pack_mechanism_assets should respect reversed priority."""
        config = {
            'enabled': True,
            'copy_mode': 'copy',
            # P0-2: Reversed priority - should choose S2 even though S3 exists
            'dipole_source_priority': ['S2_reactant_complex', 'S3_reactant']
        }

        pack_mechanism_assets(
            step_dirs=step_dirs_with_all_sources,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=((0, 1), (2, 3))
        )

        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        assert mech_index_path.exists()

        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        dipole_asset = mech_index['assets']['mech_step2_reactant_dipole']
        assert dipole_asset is not None
        # Should be S2 because priority is reversed
        assert dipole_asset['source_label'] == 'S2_reactant_complex'

    def test_dipole_priority_single_source_always_wins(
        self,
        step_dirs_s3_only
    ):
        """If only one source exists in priority list, it should always be chosen."""
        s3_reactant = step_dirs_s3_only['S3'] / 'reactant_sp.xyz'
        s2_reactant = step_dirs_s3_only['S2'] / 'reactant_complex.xyz'  # doesn't exist

        # Even with S2 first, S3 should be chosen because S2 doesn't exist
        resolved_path, source_label = _resolve_dipole_source(
            s3_reactant, s2_reactant, ['S2_reactant_complex', 'S3_reactant']
        )

        assert source_label == 'S3_reactant'
        assert resolved_path == s3_reactant


# ============================================================================
# Test Class 2: TestFixedNamingContract
# ============================================================================

class TestFixedNamingContract:
    """Test fixed naming contract for all assets."""

    def test_fixed_naming_contract(
        self,
        step_dirs_with_all_sources,
        tmp_path
    ):
        """All assets must use exact fixed filenames."""
        config = {'enabled': True, 'copy_mode': 'copy'}

        pack_mechanism_assets(
            step_dirs=step_dirs_with_all_sources,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=((0, 1), (2, 3))
        )

        s4_dir = tmp_path / "S4_Data"

        required_assets = [
            'mech_step2_ts2.xyz',
            'mech_step2_reactant_dipole.xyz',
            'mech_step2_product.xyz',
            'mech_step2_meta.json',
            'mech_step1_meta.json',
            'mech_index.json'
        ]

        for asset_name in required_assets:
            asset_path = s4_dir / asset_name
            assert asset_path.exists(), f"Missing required asset: {asset_name}"

        # Verify no mechanism/ subdirectory
        assert not (s4_dir / 'mechanism').exists()

    def test_mech_index_json_structure(self, step_dirs_with_all_sources, tmp_path):
        """Verify mech_index.json structure."""
        config = {'enabled': True, 'copy_mode': 'copy'}

        pack_mechanism_assets(
            step_dirs=step_dirs_with_all_sources,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=((0, 1), (2, 3))
        )

        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        # Verify top-level keys
        required_keys = ['version', 'schema_version', 'generated_at', 'mechanism_status', 'assets', 'quality_flags', 'config']
        for key in required_keys:
            assert key in mech_index, f"Missing key: {key}"

        # Verify status
        assert mech_index['mechanism_status'] in ['COMPLETE', 'INCOMPLETE']

        # Verify quality flags
        assert 'atom_count_ok' in mech_index['quality_flags']
        assert 'forming_bond_window_ok' in mech_index['quality_flags']
        assert 'suspect_optimized_to_product' in mech_index['quality_flags']
        # M4-P0: Returns None when insufficient data (not hardcoded "unknown")
        assert mech_index['quality_flags']['suspect_optimized_to_product'] is None


# ============================================================================
# Test Class 3: TestMissingFileDegradation
# ============================================================================

class TestMissingFileDegradation:
    """Test graceful degradation when inputs are missing."""

    def test_missing_s1_product(self, tmp_path, sample_coords, sample_symbols):
        """Missing S1 product should result in INCOMPLETE status."""
        s2 = tmp_path / "S2_Retro"
        s3 = tmp_path / "S3_TS"

        s2.mkdir(parents=True)
        s3.mkdir(parents=True)

        write_xyz(s3 / "ts_final.xyz", sample_coords, sample_symbols, title="TS Final")

        step_dirs = {"S1": None, "S2": s2, "S3": s3}
        config = {'enabled': True, 'copy_mode': 'copy'}

        pack_mechanism_assets(
            step_dirs=step_dirs,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=((0, 1), (2, 3))
        )

        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        assert mech_index_path.exists()

        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        assert mech_index['mechanism_status'] == 'INCOMPLETE'
        assert mech_index['assets']['mech_step2_product'] is None

    def test_missing_all_inputs(self, tmp_path):
        """All inputs missing should still create mech_index.json."""
        step_dirs = {"S1": None, "S2": None, "S3": None}
        config = {'enabled': True, 'copy_mode': 'copy'}

        pack_mechanism_assets(
            step_dirs=step_dirs,
            out_dir=tmp_path / "S4_Data",
            config=config
        )

        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        assert mech_index_path.exists()

        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        assert mech_index['mechanism_status'] == 'INCOMPLETE'
        # Most assets should be None
        asset_count = sum(1 for a in mech_index['assets'].values() if a is not None)
        assert asset_count < 5  # Only meta files created


# ============================================================================
# Test Class 4: TestToxicPathSafety
# ============================================================================

class TestToxicPathSafety:
    """Test toxic path safety using sandbox pattern."""

    def test_toxic_path_copy(self, tmp_path, sample_coords, sample_symbols):
        """Toxic path with spaces should use sandbox for copy."""
        # Create S4 directory with spaces and special chars
        s4_toxic = tmp_path / "[5+2] test space ()"
        s4_toxic.mkdir(parents=True)

        # Create S1 with a test file
        s1 = tmp_path / "S1_ConfGeneration"
        s1.mkdir(parents=True)
        test_file = s1 / "product_min.xyz"
        write_xyz(test_file, sample_coords, sample_symbols, title="Product")

        config = {'enabled': True, 'copy_mode': 'copy'}

        step_dirs = {"S1": s1, "S2": None, "S3": None}

        pack_mechanism_assets(
            step_dirs=step_dirs,
            out_dir=s4_toxic,
            config=config
        )

        # Verify assets created despite toxic path
        mech_index_path = s4_toxic / "mech_index.json"
        assert mech_index_path.exists(), "mech_index.json should be created in toxic path"

        product_asset = s4_toxic / "mech_step2_product.xyz"
        assert product_asset.exists(), "Asset should be copied via sandbox"

    def test_toxic_path_symlink_fallback(self, tmp_path, sample_coords, sample_symbols):
        """Symlink mode should fallback to copy if not supported."""
        s4_toxic = tmp_path / "[5+2] test space ()"
        s4_toxic.mkdir(parents=True)

        s1 = tmp_path / "S1_ConfGeneration"
        s1.mkdir(parents=True)
        test_file = s1 / "product_min.xyz"
        write_xyz(test_file, sample_coords, sample_symbols, title="Product")

        config = {'enabled': True, 'copy_mode': 'symlink'}

        step_dirs = {"S1": s1, "S2": None, "S3": None}

        pack_mechanism_assets(
            step_dirs=step_dirs,
            out_dir=s4_toxic,
            config=config
        )

        # Symlink might fail on some systems, so check for either symlink or copy
        product_asset = s4_toxic / "mech_step2_product.xyz"
        assert product_asset.exists(), "Asset should exist (symlink or copy)"


# ============================================================================
# Test Class 5: TestQualityFlags
# ============================================================================

class TestQualityFlags:
    """Test quality flag computation."""

    def test_atom_count_ok_match(
        self,
        step_dirs_with_all_sources,
        tmp_path
    ):
        """atom_count_ok should be True when dipole and ts2 have same atom count."""
        config = {
            'enabled': True,
            'copy_mode': 'copy',
            'write_quality_flags': True
        }

        pack_mechanism_assets(
            step_dirs=step_dirs_with_all_sources,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=((0, 1), (2, 3))
        )

        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        assert mech_index['quality_flags']['atom_count_ok'] is True

    def test_forming_bond_window_ok(self, step_dirs_with_all_sources, tmp_path):
        """forming_bond_window_ok should be computed."""
        config = {
            'enabled': True,
            'copy_mode': 'copy',
            'write_quality_flags': True
        }

        pack_mechanism_assets(
            step_dirs=step_dirs_with_all_sources,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=((0, 1), (2, 3))
        )

        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        assert 'forming_bond_window_ok' in mech_index['quality_flags']
        # Should be True or False, not None
        assert isinstance(mech_index['quality_flags']['forming_bond_window_ok'], bool)

    def test_suspect_optimized_to_product_placeholder(
        self,
        step_dirs_with_all_sources,
        tmp_path
    ):
        """M4-P0: suspect_optimized_to_product returns None when insufficient data for判断."""
        config = {
            'enabled': True,
            'copy_mode': 'copy',
            'write_quality_flags': True
        }

        pack_mechanism_assets(
            step_dirs=step_dirs_with_all_sources,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=((0, 1), (2, 3))
        )

        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        # M4-P0: Returns None when insufficient data (no product file to compare)
        assert mech_index['quality_flags']['suspect_optimized_to_product'] is None

    def test_quality_flags_disabled(
        self,
        step_dirs_with_all_sources,
        tmp_path
    ):
        """Quality flags should be default when write_quality_flags=False."""
        config = {
            'enabled': True,
            'copy_mode': 'copy',
            'write_quality_flags': False
        }

        pack_mechanism_assets(
            step_dirs=step_dirs_with_all_sources,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=((0, 1), (2, 3))
        )

        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        assert mech_index['quality_flags']['atom_count_ok'] is None
        assert mech_index['quality_flags']['forming_bond_window_ok'] is None
        # M4-P0: Returns None when insufficient data
        assert mech_index['quality_flags']['suspect_optimized_to_product'] is None


# ============================================================================
# Test Helper Functions
# ============================================================================

class TestHelperFunctions:
    """Test helper functions directly."""

    def test_validate_atom_count_match(self, tmp_path, sample_coords, sample_symbols):
        """_validate_atom_count should return True for matching counts."""
        test_file = tmp_path / "test.xyz"
        write_xyz(test_file, sample_coords, sample_symbols, title="Test")

        result = _validate_atom_count(test_file, expected_count=3)
        assert result is True

    def test_validate_atom_count_mismatch(self, tmp_path):
        """_validate_atom_count should return False for mismatch."""
        coords = np.array([[0.0, 0.0, 0.0]])
        symbols = ["C"]
        test_file = tmp_path / "test.xyz"
        write_xyz(test_file, coords, symbols, title="Test")

        result = _validate_atom_count(test_file, expected_count=5)
        assert result is False

    def test_check_forming_bond_window_valid(self, tmp_path):
        """Valid forming bonds should pass window check."""
        test_file = tmp_path / "test.xyz"

        # Create structure with 4 atoms to support bonds (0,1) and (2,3)
        valid_coords = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],  # Bond 0-1 at 2.0 Angstroms
            [0.5, 0.866, 0.0],
            [2.5, 0.866, 0.0]  # Bond 2-3 at 2.0 Angstroms
        ])
        valid_symbols = ["C", "H", "H", "H"]
        write_xyz(test_file, valid_coords, valid_symbols, title="Test")

        # Bonds (0,1) and (2,3) with distance ~2.0 (within 1.5-3.5 window)
        result = _check_forming_bond_window(test_file, forming_bonds=((0, 1), (2, 3)))
        assert result is True

    def test_check_forming_bond_window_outside(self, tmp_path):
        """Bonds outside window should fail check."""
        coords = np.array([
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],  # 4.0 Angstroms away
            [0.5, 0.866, 0.0],
            [1.5, 0.866, 0.0]
        ])
        symbols = ["C", "H", "H", "H"]
        test_file = tmp_path / "test.xyz"
        write_xyz(test_file, coords, symbols, title="Test")

        # Bond (0,1) is 4.0 Angstroms, outside 1.5-3.5 window
        result = _check_forming_bond_window(test_file, forming_bonds=((0, 1), (2, 3)))
        assert result is False

    def test_check_forming_bond_no_data(self, tmp_path, sample_coords, sample_symbols):
        """No forming bonds should return True."""
        test_file = tmp_path / "test.xyz"
        write_xyz(test_file, sample_coords, sample_symbols, title="Test")

        result = _check_forming_bond_window(test_file, forming_bonds=None)
        assert result is True  # No bonds data, assume OK


# ============================================================================
# Test Class 6: TestSchemaMigration (M3-2-4)
# ============================================================================

class TestSchemaMigration:
    """Test M3-2 schema migration from old to new format."""

    def test_migrate_old_schema_to_v1(self, tmp_path):
        """Old schema without schema_version should get migrated to v1."""
        from rph_core.steps.step4_features.mech_packager import migrate_mech_index

        old_index = {
            "version": "1.0.0",
            "timestamp": "2026-01-01T00:00:00Z",
            "status": "COMPLETE",
            "assets": {},
            "quality_flags": {},
            "config": {}
        }

        migrated = migrate_mech_index(old_index)

        assert 'schema_version' in migrated
        assert migrated['schema_version'] == 'mech_index_v1'
        assert migrated['generated_at'] == "2026-01-01T00:00:00Z"
        assert migrated['mechanism_status'] == "COMPLETE"

    def test_migrate_with_deprecated_aliases(self, tmp_path):
        """Migration should create deprecated aliases for backward compatibility when old fields exist."""
        from rph_core.steps.step4_features.mech_packager import migrate_mech_index

        old_index = {
            "version": "1.0.0",
            "timestamp": "2026-01-01T00:00:00Z",
            "status": "COMPLETE",
            "assets": {},
            "quality_flags": {},
            "config": {}
        }

        migrated = migrate_mech_index(old_index)

        assert migrated['generated_at'] == "2026-01-01T00:00:00Z"
        assert 'timestamp' in migrated

    def test_migrate_from_new_fields_only(self, tmp_path):
        """Input with only new fields should not create deprecated aliases."""
        from rph_core.steps.step4_features.mech_packager import migrate_mech_index

        new_only_index = {
            "version": "1.0.0",
            "generated_at": "2026-01-01T00:00:00Z",
            "mechanism_status": "COMPLETE",
            "assets": {},
            "quality_flags": {},
            "config": {}
        }

        migrated = migrate_mech_index(new_only_index)

        assert 'timestamp' not in migrated
        assert 'status' not in migrated
        assert migrated['generated_at'] == "2026-01-01T00:00:00Z"
        assert migrated['mechanism_status'] == "COMPLETE"

    def test_v1_schema_no_migration_needed(self, tmp_path):
        """v1 schema should pass through migration unchanged."""
        from rph_core.steps.step4_features.mech_packager import migrate_mech_index

        v1_index = {
            "version": "1.0.0",
            "schema_version": "mech_index_v1",
            "generated_at": "2026-01-01T00:00:00Z",
            "mechanism_status": "COMPLETE",
            "assets": {},
            "quality_flags": {},
            "config": {}
        }

        migrated = migrate_mech_index(v1_index)

        assert migrated == v1_index

    def test_migration_preserves_other_fields(self, tmp_path):
        """Migration should preserve all other fields unchanged."""
        from rph_core.steps.step4_features.mech_packager import migrate_mech_index

        old_index = {
            "version": "1.0.0",
            "timestamp": "2026-01-01T00:00:00Z",
            "status": "COMPLETE",
            "extra_field": "preserve_me",
            "assets": {"asset1": {"filename": "test.xyz"}},
            "quality_flags": {"atom_count_ok": True},
            "config": {"enabled": True}
        }

        migrated = migrate_mech_index(old_index)

        assert migrated['version'] == "1.0.0"
        assert migrated['extra_field'] == "preserve_me"
        assert migrated['assets'] == {"asset1": {"filename": "test.xyz"}}
        assert migrated['quality_flags'] == {"atom_count_ok": True}
        assert migrated['config'] == {"enabled": True}


# ============================================================================
# Test Class 8: M4-P0 QC Artifact Collection (Fixed Selection Strategy)
# ============================================================================

class TestM4P0QCArtifactCollection:
    """Test M4-P0 fixes for _collect_qc_artifacts selection strategy."""

    def test_collect_qc_artifacts_picks_by_mtime(self, tmp_path):
        """M4-P0: Should pick newest file by mtime, not first-found."""
        from rph_core.steps.step4_features.mech_packager import _collect_qc_artifacts

        s3_dir = tmp_path / "S3_TS"
        s3_dir.mkdir(parents=True)
        nbo_dir = s3_dir / "nbo"
        nbo_dir.mkdir(exist_ok=True)

        # Create two NBO files with different mtimes
        old_file = nbo_dir / "job_old.37"
        new_file = nbo_dir / "job_new.37"

        old_file.write_text("old data")
        new_file.write_text("new data")

        # Make new_file newer
        import time
        time.sleep(0.1)
        new_file.touch()

        result = _collect_qc_artifacts(
            s3_dir=s3_dir,
            pipeline_root=tmp_path,
            out_dir=tmp_path / "S4_Data",
            copy_mode="copy"
        )

        # Should pick the newer file
        assert "nbo_outputs" in result
        assert result["nbo_outputs"]["meta"]["reason"] == "picked_by_mtime"
        # Meta should have candidates list
        assert "candidates" in result["nbo_outputs"]["meta"]
        assert len(result["nbo_outputs"]["meta"]["candidates"]) >= 1

    def test_collect_qc_artifacts_meta_has_candidates(self, tmp_path):
        """M4-P0: Meta should record candidates with mtime, size, sha256."""
        from rph_core.steps.step4_features.mech_packager import _collect_qc_artifacts

        s3_dir = tmp_path / "S3_TS"
        s3_dir.mkdir(parents=True)
        nbo_dir = s3_dir / "nbo"
        nbo_dir.mkdir(exist_ok=True)

        # Create an NBO file
        nbo_file = s3_dir / "nbo" / "test.37"
        nbo_file.write_text("NBO data")

        result = _collect_qc_artifacts(
            s3_dir=s3_dir,
            pipeline_root=tmp_path,
            out_dir=tmp_path / "S4_Data",
            copy_mode="copy"
        )

        assert "nbo_outputs" in result
        meta = result["nbo_outputs"]["meta"]

        # Should have candidates with full metadata
        assert "candidates" in meta
        assert "picked" in meta
        assert "reason" in meta

        # Picked should have mtime and size
        assert "mtime" in meta["picked"]
        assert "size" in meta["picked"]
        assert "rel_path" in meta["picked"]

    def test_collect_qc_artifacts_empty_dir(self, tmp_path):
        """Empty S3 directory should return empty result."""
        from rph_core.steps.step4_features.mech_packager import _collect_qc_artifacts

        s3_dir = tmp_path / "S3_TS"
        s3_dir.mkdir(parents=True)

        result = _collect_qc_artifacts(
            s3_dir=s3_dir,
            pipeline_root=tmp_path,
            out_dir=tmp_path / "S4_Data",
            copy_mode="copy"
        )

        assert result == {}


# ============================================================================
# Test Class 9: M4-P0 Suspect Optimized to Product Flag
# ============================================================================

class TestM4P0SuspectOptimizedToProduct:
    """Test M4-P0 suspect_optimized_to_product flag logic."""

    def test_suspect_flag_returns_none_no_forming_bonds(self, tmp_path):
        """No forming bonds -> None (unknown)."""
        from rph_core.steps.step4_features.mech_packager import _check_suspect_optimized_to_product

        result = _check_suspect_optimized_to_product(
            ts_path=tmp_path / "ts.xyz",
            product_path=tmp_path / "product.xyz",
            forming_bonds=None
        )

        assert result is None

    def test_suspect_flag_returns_none_no_ts(self, tmp_path):
        """No TS file -> None (unknown)."""
        from rph_core.steps.step4_features.mech_packager import _check_suspect_optimized_to_product

        result = _check_suspect_optimized_to_product(
            ts_path=tmp_path / "nonexistent.xyz",
            product_path=tmp_path / "product.xyz",
            forming_bonds=((0, 1),)
        )

        assert result is None

    def test_suspect_flag_returns_ok_reasonable_distances(self, tmp_path):
        """TS in window, product shorter -> 'ok'."""
        from rph_core.steps.step4_features.mech_packager import _check_suspect_optimized_to_product

        # Create TS file with bond at 2.5 Angstroms (in window)
        ts_file = tmp_path / "ts.xyz"
        ts_coords = np.array([
            [0.0, 0.0, 0.0],
            [2.5, 0.0, 0.0],  # Bond at 2.5 Å (in 1.5-3.5 window)
            [0.5, 0.866, 0.0]
        ])
        write_xyz(ts_file, ts_coords, ["C", "H", "H"], title="TS")

        # Create product file with bond at 1.5 Angstroms (shorter)
        product_file = tmp_path / "product.xyz"
        product_coords = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],  # Bond at 1.5 Å (shorter = normal product)
            [0.5, 0.866, 0.0]
        ])
        write_xyz(product_file, product_coords, ["C", "H", "H"], title="Product")

        result = _check_suspect_optimized_to_product(
            ts_path=ts_file,
            product_path=product_file,
            forming_bonds=((0, 1),)
        )

        assert result == "ok"

    def test_suspect_flag_returns_suspect_similar_distances(self, tmp_path):
        """TS in window, product similar distance -> 'suspect'."""
        from rph_core.steps.step4_features.mech_packager import _check_suspect_optimized_to_product

        # Create TS file with bond at 2.0 Angstroms (in window)
        ts_file = tmp_path / "ts.xyz"
        ts_coords = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],  # Bond at 2.0 Å (in window)
            [0.5, 0.866, 0.0]
        ])
        write_xyz(ts_file, ts_coords, ["C", "H", "H"], title="TS")

        # Create product file with bond at 1.9 Angstroms (nearly same as TS!)
        product_file = tmp_path / "product.xyz"
        product_coords = np.array([
            [0.0, 0.0, 0.0],
            [1.9, 0.0, 0.0],  # Bond at 1.9 Å (90% of TS = suspicious)
            [0.5, 0.866, 0.0]
        ])
        write_xyz(product_file, product_coords, ["C", "H", "H"], title="Product")

        result = _check_suspect_optimized_to_product(
            ts_path=ts_file,
            product_path=product_file,
            forming_bonds=((0, 1),)
        )

        assert result == "suspect"

    def test_suspect_flag_returns_none_ts_outside_window(self, tmp_path):
        """TS outside window -> None (cannot判断)."""
        from rph_core.steps.step4_features.mech_packager import _check_suspect_optimized_to_product

        # Create TS file with bond at 1.0 Angstroms (too short, outside window)
        ts_file = tmp_path / "ts.xyz"
        ts_coords = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],  # Bond at 1.0 Å (outside 1.5-3.5 window)
            [0.5, 0.866, 0.0]
        ])
        write_xyz(ts_file, ts_coords, ["C", "H", "H"], title="TS")

        result = _check_suspect_optimized_to_product(
            ts_path=ts_file,
            product_path=tmp_path / "product.xyz",
            forming_bonds=((0, 1),)
        )

        assert result is None

    def test_suspect_flag_returns_none_no_product(self, tmp_path):
        """No product file -> None (cannot compare)."""
        from rph_core.steps.step4_features.mech_packager import _check_suspect_optimized_to_product

        # Create TS file with bond in window
        ts_file = tmp_path / "ts.xyz"
        ts_coords = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],  # Bond at 2.0 Å (in window)
            [0.5, 0.866, 0.0]
        ])
        write_xyz(ts_file, ts_coords, ["C", "H", "H"], title="TS")

        result = _check_suspect_optimized_to_product(
            ts_path=ts_file,
            product_path=None,
            forming_bonds=((0, 1),)
        )

        assert result is None
