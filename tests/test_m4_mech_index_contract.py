"""
Tests for mech_index.json schema contract and quality flags.

Tests verify:
- mechanism_status field presence and valid values
- quality_flags contains ts_imag_freq_ok and asset_hash_ok
- migrate_mech_index correctly handles old schema
"""
import json
import tempfile
from pathlib import Path

import pytest

from rph_core.steps.step4_features.mech_packager import (
    _build_mech_index,
    migrate_mech_index,
    is_mech_index_up_to_date,
    ensure_mech_index_schema,  # P0-3: New function for migration + write-back
    QualityFlags,
    MechanismMetaStep1,
    MechanismMetaStep2,
    UpdateReason,
)


class TestMechIndexContract:
    """Tests for mech_index.json schema contract."""

    def test_mechanism_status_field_exists_and_valid(self):
        """mechanism_status field must be present with valid values."""
        result = _build_mech_index(
            assets={},
            quality_flags=QualityFlags(),
            step2_meta=MechanismMetaStep2(),
            step1_meta=MechanismMetaStep1(),
            config={}
        )
        assert 'mechanism_status' in result
        assert result['mechanism_status'] in ('COMPLETE', 'INCOMPLETE')

    def test_quality_flags_contains_new_fields(self):
        """quality_flags must contain ts_imag_freq_ok and asset_hash_ok."""
        result = _build_mech_index(
            assets={},
            quality_flags=QualityFlags(),
            step2_meta=MechanismMetaStep2(),
            step1_meta=MechanismMetaStep1(),
            config={}
        )
        assert 'ts_imag_freq_ok' in result['quality_flags']
        assert 'asset_hash_ok' in result['quality_flags']
        # These fields default to None (not True placeholder)
        assert result['quality_flags']['ts_imag_freq_ok'] is None
        assert result['quality_flags']['asset_hash_ok'] is None

    def test_complete_mech_index_structure(self):
        """Full mech_index should have all required top-level fields."""
        qf = QualityFlags(atom_count_ok=True, forming_bond_window_ok=True)
        result = _build_mech_index(
            assets={'mech_step2_ts2': None},
            quality_flags=qf,
            step2_meta=MechanismMetaStep2(),
            step1_meta=MechanismMetaStep1(),
            config={'enabled': True},
            missing_inputs=['mech_step2_ts2.xyz'],
            degradation_reasons=['S3 TS optimization not completed']
        )

        required_fields = [
            'version', 'schema_version', 'generated_at', 'mechanism_status',
            'assets', 'quality_flags', 'missing_inputs', 'degradation_reasons',
            'config', 'qc_artifacts'
        ]
        for field in required_fields:
            assert field in result, f"Missing required field: {field}"


class TestMigrateMechIndex:
    """Tests for mech_index migration from old schemas."""

    def test_migrate_status_to_mechanism_status(self):
        """Old status field should be migrated to mechanism_status."""
        old = {'status': 'COMPLETE', 'other_field': 'value'}
        migrated = migrate_mech_index(old)

        assert migrated['mechanism_status'] == 'COMPLETE'
        assert migrated['status'] == 'COMPLETE'  # Keep old field for backward compat

    def test_migrate_timestamp_to_generated_at(self):
        """Old timestamp field should be migrated to generated_at."""
        old = {'timestamp': '2024-01-01T00:00:00Z'}
        migrated = migrate_mech_index(old)

        assert migrated['generated_at'] == '2024-01-01T00:00:00Z'
        assert migrated['timestamp'] == '2024-01-01T00:00:00Z'

    def test_migrate_add_schema_version(self):
        """Missing schema_version should be added."""
        old = {'status': 'COMPLETE'}
        migrated = migrate_mech_index(old)

        assert migrated['schema_version'] == 'mech_index_v1'

    def test_migrate_preserves_existing_fields(self):
        """Migration should preserve fields that already exist."""
        old = {
            'schema_version': 'custom_v2',
            'status': 'INCOMPLETE',
            'custom_field': 'custom_value'
        }
        migrated = migrate_mech_index(old)

        assert migrated['schema_version'] == 'custom_v2'  # Not overwritten
        assert migrated['mechanism_status'] == 'INCOMPLETE'
        assert migrated['custom_field'] == 'custom_value'


class TestIsMechIndexUpToDate:
    """Tests for mech_index validation."""

    def test_missing_file_returns_false(self):
        """Missing mech_index.json should return False with MISSING_FILE reason."""
        with tempfile.TemporaryDirectory() as tmpdir:
            s4_dir = Path(tmpdir)
            is_up_to_date, reason = is_mech_index_up_to_date(s4_dir)

            assert is_up_to_date is False
            assert reason == UpdateReason.MISSING_FILE

    def test_valid_current_schema(self):
        """Valid mech_index_v1 should return True."""
        with tempfile.TemporaryDirectory() as tmpdir:
            s4_dir = Path(tmpdir)
            index_path = s4_dir / "mech_index.json"

            # Create a valid index
            index = {
                'schema_version': 'mech_index_v1',
                'mechanism_status': 'COMPLETE'
            }
            index_path.write_text(json.dumps(index))

            is_up_to_date, reason = is_mech_index_up_to_date(s4_dir)

            assert is_up_to_date is True
            assert reason == UpdateReason.OK

    def test_old_schema_gets_migrated(self):
        """
        P0-3 FIX: is_mech_index_up_to_date is now a PURE CHECK (no side effects).
        Migration + write-back is handled by ensure_mech_index_schema().

        This test verifies the CHECK behavior: old schema returns SCHEMA_MISMATCH.
        Use ensure_mech_index_schema() if you need migration + write-back.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            s4_dir = Path(tmpdir)
            index_path = s4_dir / "mech_index.json"

            # Create an old-style index (without schema_version)
            index = {'status': 'COMPLETE'}
            index_path.write_text(json.dumps(index))

            # is_mech_index_up_to_date should return False (pure check, no migration)
            is_up_to_date, reason = is_mech_index_up_to_date(s4_dir)

            assert is_up_to_date is False
            # Reason should indicate schema mismatch (starts with the SCHEMA_MISMATCH constant)
            assert reason.startswith(UpdateReason.SCHEMA_MISMATCH)
            assert "mech_index_v1" in reason

            # File should NOT have been modified (no side effects)
            original = json.loads(index_path.read_text())
            assert 'schema_version' not in original

    def test_ensure_mech_index_schema_migrates_and_writes_back(self):
        """
        P0-3: ensure_mech_index_schema provides migration + write-back semantics.
        This is the correct function to use when you need to migrate old files.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            s4_dir = Path(tmpdir)
            index_path = s4_dir / "mech_index.json"

            # Create an old-style index (without schema_version)
            index = {'status': 'COMPLETE'}
            index_path.write_text(json.dumps(index))

            # ensure_mech_index_schema should migrate and write back
            success, reason, migrated = ensure_mech_index_schema(s4_dir)

            assert success is True
            assert reason == UpdateReason.MIGRATED
            assert migrated is not None
            assert migrated['schema_version'] == 'mech_index_v1'

            # Verify file was actually written
            written = json.loads(index_path.read_text())
            assert written.get('schema_version') == 'mech_index_v1'

    def test_missing_mechanism_status_gets_filled(self):
        """
        P0-3 FIX: is_mech_index_up_to_date is a pure check (no side effects).
        If mechanism_status is missing from v1 schema, it returns False with SCHEMA_MISMATCH.
        Use ensure_mech_index_schema() to fill missing fields and write back.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            s4_dir = Path(tmpdir)
            index_path = s4_dir / "mech_index.json"

            # Create index with old status but no mechanism_status
            index = {'schema_version': 'mech_index_v1', 'status': 'INCOMPLETE'}
            index_path.write_text(json.dumps(index))

            # is_mech_index_up_to_date should return False (missing canonical field)
            is_up_to_date, reason = is_mech_index_up_to_date(s4_dir)

            # File should NOT have been modified
            original = json.loads(index_path.read_text())
            assert 'mechanism_status' not in original

    def test_ensure_mech_index_schema_fills_missing_fields(self):
        """
        P0-3: ensure_mech_index_schema should fill missing canonical fields.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            s4_dir = Path(tmpdir)
            index_path = s4_dir / "mech_index.json"

            # Create index with old status but no mechanism_status
            index = {'schema_version': 'mech_index_v1', 'status': 'INCOMPLETE'}
            index_path.write_text(json.dumps(index))

            # ensure_mech_index_schema should fill the missing field and write back
            success, reason, migrated = ensure_mech_index_schema(s4_dir)

            assert success is True
            # File should have been updated
            written = json.loads(index_path.read_text())
            assert written.get('mechanism_status') == 'INCOMPLETE'
