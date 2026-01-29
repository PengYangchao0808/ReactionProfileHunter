"""
M2-TEST-D: Schema Versioning and Three-State Quality Flags Tests
==================================================

Tests for M2-D schema upgrades in mech_index.json:
- schema_version field presence
- generated_at ISO timestamp with timezone
- mechanism_status field (COMPLETE/INCOMPLETE)
- Three-state quality flags (true/false/null instead of bool)

Author: QC Descriptors Team
Date: 2026-01-21
"""

import json
import re

import pytest

from rph_core.steps.step4_features.mech_packager import pack_mechanism_assets
from rph_core.utils.file_io import write_xyz


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def sample_coords():
    """Sample coordinates for test molecules."""
    return [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.866, 0.0],
        [2.0, 0.0, 0.0]
    ]


@pytest.fixture
def sample_symbols():
    """Sample element symbols for test molecules."""
    return ["C", "H", "H"]


@pytest.fixture
def full_s4_dir_with_mech_index(tmp_path):
    """
    S4 directory with features_raw.csv and mech_index.json (M2 behavior).
    """
    s4 = tmp_path / "S4_Data"
    s4.mkdir(parents=True, exist_ok=True)

    # features_raw.csv placeholder
    (s4 / "features_raw.csv").write_text("reaction_id,dG_activation\n")

    # M2: mech_index.json with full schema
    mech_index_content = {
        "version": "1.0.0",
        "schema_version": "mech_index_v1",
        "generated_at": "2026-01-21T00:00:00.000Z",
        "mechanism_status": "COMPLETE",
        "assets": {
            "mech_step2_ts2.xyz": {
                "filename": "mech_step2_ts2.xyz",
                "source_path": "/mock/s3/ts_final.xyz",
                "source_step": "S3",
                "sha256": "mock_hash_ts2"
            },
            "mech_step2_reactant_dipole.xyz": {
                "filename": "mech_step2_reactant_dipole.xyz",
                "source_path": "/mock/s3/reactant_sp.xyz",
                "source_step": "S3",
                "sha256": "mock_hash_dipole"
            },
            "mech_step2_product.xyz": {
                "filename": "mech_step2_product.xyz",
                "source_path": "/mock/s1/product_min.xyz",
                "source_step": "S1",
                "sha256": "mock_hash_product"
            },
            "mech_step1_precursor.xyz": {
                "filename": "mech_step1_precursor.xyz",
                "source_path": "/mock/s1/precursor.xyz",
                "source_step": "S1",
                "sha256": "mock_hash_precursor"
            }
        },
        "quality_flags": {
            "atom_count_ok": True,
            "forming_bond_window_ok": True,
            "suspect_optimized_to_product": "unknown"
        },
        "missing_inputs": [],
        "degradation_reasons": [],
        "config": {
            "enabled": True,
            "schema_version": "mech_index_v1"
        }
    }

    (s4 / "mech_index.json").write_text(json.dumps(mech_index_content, indent=2))


@pytest.fixture
    def s4_with_mech_index_null_atom_count(tmp_path):
        """
        S4 directory with mech_index.json where atom_count_ok is null (M2-D three-state flag).
        """
        s4 = tmp_path / "S4_Data"
        s4.mkdir(parents=True)

        # features_raw.csv placeholder
        (s4 / "features_raw.csv").write_text("reaction_id,dG_activation\n")

        # M2: mech_index.json with null atom_count_ok (three-state flag)
        mech_index_content = {
            "version": "1.0.0",
            "schema_version": "mech_index_v1",
            "generated_at": "2026-01-21T00:00:00.000Z",
            "mechanism_status": "COMPLETE",
            "assets": {
                "mech_step2_ts2.xyz": None,
                "mech_step2_reactant_dipole.xyz": None,
                "mech_step2_product.xyz": None,
                "mech_step1_precursor.xyz": None
            },
            "quality_flags": {
                "atom_count_ok": None,
                "forming_bond_window_ok": True,
                "suspect_optimized_to_product": "unknown"
            },
            "missing_inputs": [],
            "degradation_reasons": [],
            "config": {
                "enabled": True,
                "schema_version": "mech_index_v1"
            }
        }

        (s4 / "mech_index.json").write_text(json.dumps(mech_index_content, indent=2))


@pytest.fixture
def s4_with_mech_index_false_atom_count(tmp_path):
    """
    S4 directory with mech_index.json where atom_count_ok is false (M2-D: failing quality check).
    """
    s4 = tmp_path / "S4_Data"
    s4.mkdir(parents=True, exist_ok=True)

    # features_raw.csv placeholder
    (s4 / "features_raw.csv").write_text("reaction_id,dG_activation\n")

    # M2: mech_index.json with false atom_count_ok (quality failure)
    mech_index_content = {
        "version": "1.0.0",
        "schema_version": "mech_index_v1",
        "generated_at": "2026-01-21T00:00:00.000Z",
        "mechanism_status": "INCOMPLETE",
        "assets": {
            "mech_step2_ts2.xyz": {
                "filename": "mech_step2_ts2.xyz",
                "source_path": "/mock/s3/ts_final.xyz",
                "source_step": "S3",
                "sha256": "mock_hash_ts2"
            },
            "mech_step2_reactant_dipole.xyz": {
                "filename": "mech_step2_reactant_dipole.xyz",
                "source_path": "/mock/s3/reactant_sp.xyz",
                "source_step": "S3",
                "sha256": "mock_hash_dipole"
            },
            "mech_step2_product.xyz": None,
            "mech_step1_precursor.xyz": None
        },
        "quality_flags": {
            "atom_count_ok": False,  # M2-D: failing quality check
            "forming_bond_window_ok": False,
            "suspect_optimized_to_product": "unknown"
        },
        "missing_inputs": ["mech_step2_product.xyz", "mech_step1_precursor.xyz"],
        "degradation_reasons": [
            "Atom count mismatch between dipole and TS2"
        ],
        "config": {
            "enabled": True,
            "schema_version": "mech_index_v1"
        }
    }

    (s4 / "mech_index.json").write_text(json.dumps(mech_index_content, indent=2))


@pytest.fixture
def s4_with_mech_incomplete_missing_inputs(tmp_path):
    """
    S4 directory where mech_index.json has INCOMPLETE status with missing_inputs list.
    """
    s4 = tmp_path / "S4_Data"
    s4.mkdir(parents=True, exist_ok=True)

    # features_raw.csv placeholder
    (s4 / "features_raw.csv").write_text("reaction_id,dG_activation\n")

    # M2: mech_index.json with INCOMPLETE status and missing_inputs
    mech_index_content = {
        "version": "1.0.0",
        "schema_version": "mech_index_v1",
        "generated_at": "2026-01-21T00:00:00.000Z",
        "mechanism_status": "INCOMPLETE",
        "assets": {
            "mech_step2_ts2.xyz": None,
            "mech_step2_reactant_dipole.xyz": None,
            "mech_step2_product.xyz": None,
            "mech_step1_precursor.xyz": None
        },
        "quality_flags": {
            "atom_count_ok": True,
            "forming_bond_window_ok": True,
            "suspect_optimized_to_product": "unknown"
        },
        "missing_inputs": [
            "mech_step2_ts2.xyz",
            "mech_step2_reactant_dipole.xyz",
            "mech_step2_product.xyz",
            "mech_step1_precursor.xyz"
        ],
        " degradation_reasons": [
            "All mechanism assets missing (simulation)"
        ],
        "config": {
            "enabled": True,
            "schema_version": "mech_index_v1"
        }
    }

    (s4 / "mech_index.json").write_text(json.dumps(mech_index_content, indent=2))


# ============================================================================
# Test Class: TestSchemaVersioning
# ============================================================================

class TestSchemaVersioning:
    """Test M2-D schema versioning in mech_index.json."""

    def test_schema_version_present(self, full_s4_dir_with_mech_index):
        """schema_version field must be present."""
        with open(full_s4_dir_with_mech_index["S4"] / "mech_index.json", 'r') as f:
            mech_index = json.load(f)

        assert 'schema_version' in mech_index, "schema_version field missing"
        assert mech_index['schema_version'] == 'mech_index_v1'

    def test_generated_at_iso_timestamp(self, full_s4_dir_with_mech_index):
        """generated_at should be ISO 8601 timestamp with timezone."""
        with open(full_s4_dir_with_mech_index["S4"] / "mech_index.json", 'r') as f:
            mech_index = json.load(f)

        assert 'generated_at' in mech_index, "generated_at field missing"
        
        # M2-D: Check ISO 8601 format
        iso_pattern = r'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d+Z$'
        assert re.match(iso_pattern, mech_index['generated_at']), f"Invalid ISO 8601 timestamp format: {mech_index['generated_at']}"

    def test_mechanism_status_complete_field(self, full_s4_dir_with_mech_index):
        """mechanism_status field must be present (COMPLETE or INCOMPLETE)."""
        with open(full_s4_dir_with_mech_index["S4"] / "mech_index.json", 'r') as f:
            mech_index = json.load(f)

        assert 'mechanism_status' in mech_index
        assert mech_index['mechanism_status'] in ['COMPLETE', 'INCOMPLETE'], f"Invalid mechanism_status value: {mech_index['mechanism_status']}"


# ============================================================================
# Test Class: TestThreeStateQualityFlags
# ============================================================================

class TestThreeStateQualityFlags:
    """Test M2-D three-state quality flags in mech_index.json."""

    def test_atom_count_ok_true_false_null(self, s4_with_mech_index_null_atom_count, s4_with_mech_index_false_atom_count):
        """
        atom_count_ok should be true, false, or null (three-state, not bool default).
        """
        # Test true value
        with open(s4_with_mech_index_null_atom_count["S4"] / "mech_index.json", 'r') as f:
            mech_index = json.load(f)

        assert mech_index['quality_flags']['atom_count_ok'] is True

        # Test false value
        with open(s4_with_mech_index_false_atom_count["S4"] / "mech_index.json", 'r') as f:
            mech_index = json.load(f)

        assert mech_index['quality_flags']['atom_count_ok'] is False

        # Test null value
        with open(s4_with_mech_index_null_atom_count["S4"] / "mech_index.json", 'r') as f:
            mech_index = json.load(f)

        assert mech_index['quality_flags']['atom_count_ok'] is None

    def test_forming_bond_window_ok_true_false_null(self, s4_with_mech_index_null_atom_count):
        """
        forming_bond_window_ok should be true, false, or null.
        """
        # Test true value
        with open(s4_with_mech_index_null_atom_count["S4"] / "mech_index.json", 'r') as f:
            mech_index = json.load(f)

        assert mech_index['quality_flags']['forming_bond_window_ok'] is True

        # Test false value
        with open(s4_with_mech_index_false_atom_count["S4"] / "mech_index.json", 'r') as f:
            mech_index = json.load(f)

        assert mech_index['quality_flags']['forming_bond_window_ok'] is False

        # Test null value
        with open(s4_with_mech_index_null_atom_count["S4"] / "mech_index.json", 'r') as f:
            mech_index = json.load(f)

        assert mech_index['quality_flags']['forming_bond_window_ok'] is None

    def test_suspect_optimized_to_product_unknown_placeholder(self, s4_with_mech_index_null_atom_count):
        """
        suspect_optimized_to_product should always be 'unknown' (M2 placeholder).
        """
        with open(s4_with_mech_index_null_atom_count["S4"] / "mech_index.json", 'r') as f:
            mech_index = json.load(f)

        assert mech_index['quality_flags']['suspect_optimized_to_product'] == 'unknown'


# ============================================================================
# Main
# ============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
