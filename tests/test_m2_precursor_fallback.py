"""
M2-TEST-B: Precursor Fallback Chain Tests
==========================================

Tests for M2-B enhanced precursor fallback logic:
- S1 precursor priority
- S2 neutral_precursor intermediate fallback
- S2 reactant_complex final fallback

Author: QC Descriptors Team
Date: 2026-01-21
"""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from rph_core.utils.file_io import write_xyz
from rph_core.steps.step4_features.mech_packager import pack_mechanism_assets


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
def s1_only_dir(tmp_path):
    """S1 only - has precursor but no S2 neutral_precursor."""
    s1 = tmp_path / "S1_ConfGeneration"
    s1.mkdir(parents=True)

    # Create S1 precursor
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    symbols = ["C", "H", "H"]
    write_xyz(s1 / "precursor.xyz", coords, symbols, title="S1 Precursor")

    s2 = tmp_path / "S2_Retro"
    s2.mkdir(parents=True)

    # Create S2 reactant_complex (final fallback)
    write_xyz(s2 / "reactant_complex.xyz", coords, symbols, title="S2 Reactant Complex")

    return {"S1": s1, "S2": s2}


@pytest.fixture
def s1_and_s2_neutral_dir(tmp_path):
    """S1 and S2 with neutral_precursor intermediate."""
    s1 = tmp_path / "S1_ConfGeneration"
    s1.mkdir(parents=True)

    # Create S1 precursor
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    symbols = ["C", "H", "H"]
    write_xyz(s1 / "precursor.xyz", coords, symbols, title="S1 Precursor")

    s2 = tmp_path / "S2_Retro"
    s2.mkdir(parents=True)

    # Create S2 neutral_precursor (intermediate)
    write_xyz(s2 / "neutral_precursor.xyz", coords, symbols, title="S2 Neutral Precursor")

    # Create S2 reactant_complex (final fallback)
    coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    symbols = ["C", "H", "H"]
    write_xyz(s2 / "reactant_complex.xyz", coords, symbols, title="S2 Reactant Complex")

    return {"S1": s1, "S2": s2}


@pytest.fixture
def s2_only_dir(tmp_path):
    """S2 only - no S1, has reactant_complex."""
    s2 = tmp_path / "S2_Retro"
    s2.mkdir(parents=True)

    # Create S2 reactant_complex
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    symbols = ["C", "H", "H"]
    write_xyz(s2 / "reactant_complex.xyz", coords, symbols, title="S2 Reactant Complex")

    return {"S2": s2}


@pytest.fixture
def full_pipeline_dir(tmp_path):
    """Full S1+S2+S3+S4 with all assets."""
    s1 = tmp_path / "S1_ConfGeneration"
    s1.mkdir(parents=True)
    s2 = tmp_path / "S2_Retro"
    s2.mkdir(parents=True)
    s3 = tmp_path / "S3_TS"
    s3.mkdir(parents=True)

    # S1 product
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    symbols = ["C", "H", "H"]
    write_xyz(s1 / "product_min.xyz", coords, symbols, title="S1 Product")

    # S2 ts_guess and reactant_complex
    write_xyz(s2 / "ts_guess.xyz", coords, symbols, title="S2 TS Guess")
    write_xyz(s2 / "reactant_complex.xyz", coords, symbols, title="S2 Reactant Complex")

    # S3 TS final
    coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    symbols = ["C", "H", "H"]
    write_xyz(s3 / "ts_final.xyz", coords, symbols, title="S3 TS Final")

    # S3 reactant SP (for dipole)
    write_xyz(s3 / "reactant_sp.xyz", coords, symbols, title="S3 Reactant SP")

    s4 = tmp_path / "S4_Data"
    s4.mkdir(parents=True)

    # Create features_raw.csv placeholder
    with open(s4 / "features_raw.csv", "w") as f:
        f.write("reaction_id,dG_activation\n")

    return {"S1": s1, "S2": s2, "S3": s3, "S4": s4}


# ============================================================================
# Test Class: TestPrecursorPriority
# ============================================================================

class TestPrecursorPriority:
    """Test M2-B: Configurable precursor priority fallback chain."""

    def test_s1_precursor_used(self, s1_and_s2_neutral_dir, tmp_path):
        """S1 precursor should be selected when available."""
        config = {
            'enabled': True,
            'copy_mode': 'copy',
            'dipole_source_priority': ['S3_reactant', 'S2_reactant_complex'],
            'precursor_source_priority': ['S1_precursor', 'S2_neutral_precursor', 'S2_reactant_complex'],
            'write_quality_flags': True
        }

        pack_mechanism_assets(
            step_dirs=s1_and_s2_neutral_dir,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=None
        )

        # Verify mech_index.json
        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        assert mech_index_path.exists()

        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        # Check precursor asset
        precursor_asset = mech_index['assets']['mech_step1_precursor']
        assert precursor_asset is not None
        assert precursor_asset['source_label'] == 'S1_precursor'

    def test_s2_neutral_precursor_used(self, s1_and_s2_neutral_dir, tmp_path):
        """S2 neutral_precursor should be selected when S1 precursor missing."""
        config = {
            'enabled': True,
            'copy_mode': 'copy',
            'dipole_source_priority': ['S3_reactant', 'S2_reactant_complex'],
            'precursor_source_priority': ['S1_precursor', 'S2_neutral_precursor', 'S2_reactant_complex'],
            'write_quality_flags': True
        }

        # Remove S1 precursor to test fallback
        (s1_and_s2_neutral_dir['S1'] / "precursor.xyz").unlink()

        pack_mechanism_assets(
            step_dirs=s1_and_s2_neutral_dir,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=None
        )

        # Verify mech_index.json
        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        assert mech_index_path.exists()

        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        precursor_asset = mech_index['assets']['mech_step1_precursor']
        assert precursor_asset is not None
        assert precursor_asset['source_label'] == 'S2_neutral_precursor'

    def test_s2_reactant_complex_fallback(self, s2_only_dir, tmp_path):
        """S2 reactant_complex should be final fallback."""
        config = {
            'enabled': True,
            'copy_mode': 'copy',
            'dipole_source_priority': ['S3_reactant', 'S2_reactant_complex'],
            'precursor_source_priority': ['S1_precursor', 'S2_neutral_precursor', 'S2_reactant_complex'],
            'write_quality_flags': True
        }

        pack_mechanism_assets(
            step_dirs=s2_only_dir,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=None
        )

        # Verify mech_index.json
        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        assert mech_index_path.exists()

        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        precursor_asset = mech_index['assets']['mech_step1_precursor']
        assert precursor_asset is not None
        assert precursor_asset['source_label'] == 'S2_reactant_complex'

    def test_custom_priority_order(self, s1_and_s2_neutral_dir, tmp_path):
        """Custom priority order should be respected."""
        config = {
            'enabled': True,
            'copy_mode': 'copy',
            'dipole_source_priority': ['S3_reactant', 'S2_reactant_complex'],
            'precursor_source_priority': ['S2_neutral_precursor', 'S2_reactant_complex', 'S1_precursor'],  # Custom order
            'write_quality_flags': True
        }

        pack_mechanism_assets(
            step_dirs=s1_and_s2_neutral_dir,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=None
        )

        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        precursor_asset = mech_index['assets']['mech_step1_precursor']
        assert precursor_asset['source_label'] == 'S2_neutral_precursor'

    def test_all_sources_missing(self, tmp_path):
        """All sources missing should still produce mech_index.json with None."""
        empty_dirs = {
            'S1': tmp_path / "S1",
            'S2': tmp_path / "S2",
            'S3': tmp_path / "S3"
        }

        config = {
            'enabled': True,
            'copy_mode': 'copy',
            'dipole_source_priority': ['S3_reactant', 'S2_reactant_complex'],
            'precursor_source_priority': ['S1_precursor', 'S2_neutral_precursor', 'S2_reactant_complex'],
            'write_quality_flags': True
        }

        pack_mechanism_assets(
            step_dirs=empty_dirs,
            out_dir=tmp_path / "S4_Data",
            config=config,
            forming_bonds=None
        )

        mech_index_path = tmp_path / "S4_Data" / "mech_index.json"
        assert mech_index_path.exists()

        with open(mech_index_path, 'r') as f:
            mech_index = json.load(f)

        assert mech_index['mechanism_status'] == 'INCOMPLETE'
        assert mech_index['assets']['mech_step1_precursor'] is None
        assert 'mech_step1_precursor.xyz' in mech_index['missing_inputs']


# ============================================================================
# Main
# ============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
