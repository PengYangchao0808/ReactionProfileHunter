"""
M4-E-4: Test for UpdateReason Constants
==============================================

Tests that UpdateReason class provides consistent constants for
is_mech_index_up_to_date() return values.

Author: QC Descriptors Team
Date: 2026-01-21
"""

import pytest

from rph_core.steps.step4_features.mech_packager import UpdateReason


class TestUpdateReasonConstants:
    """Test M4-E-4: UpdateReason constant values."""

    def test_update_reason_constants_defined(self):
        """Test that all UpdateReason constants are defined."""
        assert hasattr(UpdateReason, 'OK'), "OK constant should exist"
        assert hasattr(UpdateReason, 'MISSING_FILE'), "MISSING_FILE constant should exist"
        assert hasattr(UpdateReason, 'LOAD_FAILED'), "LOAD_FAILED constant should exist"
        assert hasattr(UpdateReason, 'SCHEMA_MISMATCH'), "SCHEMA_MISMATCH constant should exist"
        assert hasattr(UpdateReason, 'ASSET_MISSING'), "ASSET_MISSING constant should exist"

        print("✓ All UpdateReason constants are defined")

    def test_update_reason_values_are_strings(self):
        """Test that all constants are strings."""
        assert isinstance(UpdateReason.OK, str), "OK should be string"
        assert isinstance(UpdateReason.MISSING_FILE, str), "MISSING_FILE should be string"
        assert isinstance(UpdateReason.LOAD_FAILED, str), "LOAD_FAILED should be string"
        assert isinstance(UpdateReason.SCHEMA_MISMATCH, str), "SCHEMA_MISMATCH should be string"
        assert isinstance(UpdateReason.ASSET_MISSING, str), "ASSET_MISSING should be string"

        print("✓ All UpdateReason constants are strings")

    def test_update_reason_ok_value(self):
        """Test that OK constant has expected value."""
        assert UpdateReason.OK == "OK", f"OK should be 'OK', got '{UpdateReason.OK}'"
        print(f"✓ UpdateReason.OK = '{UpdateReason.OK}'")

    def test_update_reason_missing_file_value(self):
        """Test that MISSING_FILE constant has expected value."""
        assert "mech_index.json does not exist" in UpdateReason.MISSING_FILE
        print(f"✓ UpdateReason.MISSING_FILE = '{UpdateReason.MISSING_FILE}'")

    def test_update_reason_constants_are_unique(self):
        """Test that all constants have unique values."""
        values = [
            UpdateReason.OK,
            UpdateReason.MISSING_FILE,
            UpdateReason.LOAD_FAILED,
            UpdateReason.SCHEMA_MISMATCH,
            UpdateReason.ASSET_MISSING
        ]

        unique_values = list(set(values))
        assert len(values) == len(unique_values), "All UpdateReason values should be unique"

        print(f"✓ All {len(values)} UpdateReason constants are unique")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
