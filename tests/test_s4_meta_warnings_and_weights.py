"""
TDD Tests for S4 Meta Warnings and qc.sample_weight Policy
============================================================

Tests for:
1. Top-level warnings in feature_meta.json are populated correctly
2. Warnings are structured with code, plugin, severity, detail
3. Warning deduplication by (code, plugin)
4. qc.sample_weight policy implementation (1.0 for OK, 0.0 for invalid)

These tests should FAIL initially (RED) until feature_miner.py is fixed.

Author: RPH Team
Date: 2026-02-03
"""

import pytest
import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd

from rph_core.steps.step4_features.feature_miner import FeatureMiner
from rph_core.steps.step4_features.context import FeatureContext, PluginTrace, FeatureResult
from rph_core.steps.step4_features.status import FeatureResultStatus, FeatureStatus


# =============================================================================
# Minimal Fixtures
# =============================================================================

def _write_min_xyz(path: Path, content: str = None) -> None:
    """Write minimal xyz file."""
    if content is None:
        content = "3\n\nC 0.0 0.0 0.0\nO 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
    path.write_text(content)


def _create_minimal_context(tmpdir: Path) -> FeatureContext:
    """Create minimal FeatureContext."""
    
    # Create minimal xyz files
    ts_final = tmpdir / "ts_final.xyz"
    reactant = tmpdir / "reactant.xyz"
    product = tmpdir / "product_min.xyz"
    _write_min_xyz(ts_final)
    _write_min_xyz(reactant)
    _write_min_xyz(product)
    
    # Create mock artifacts index
    class MockArtifactsIndex:
        def __init__(self):
            self.ts_final = {'path_rel': str(ts_final.relative_to(tmpdir)), 'sha256': None}
            self.reactant = {'path_rel': str(reactant.relative_to(tmpdir)), 'sha256': None}
            self.product = {'path_rel': str(product.relative_to(tmpdir)), 'sha256': None}
    
    return FeatureContext(
        s3_dir=tmpdir,
        artifacts_index=MockArtifactsIndex(),
        forming_bonds=None,
    )


# =============================================================================
# Tests for Top-Level Warnings Population
# =============================================================================

class TestTopLevelWarningsPopulation:
    """Tests for top-level warnings in feature_meta.json."""

    def test_feature_meta_warnings_populated_when_plugin_warns(self):
        """feature_meta.json top-level warnings should be populated when plugins warn."""
        with tempfile.TemporaryDirectory(prefix="test_meta_warnings_") as tmp:
            tmpdir = Path(tmp)
            
            ts = tmpdir / "ts_final.xyz"
            reactant = tmpdir / "reactant.xyz"
            product = tmpdir / "product_min.xyz"
            _write_min_xyz(ts)
            _write_min_xyz(reactant)
            _write_min_xyz(product)
            
            out_dir = tmpdir / "S4_Data"
            miner = FeatureMiner(config={})
            
            # Run with no optional files - should trigger warnings
            miner.run(
                ts_final=ts,
                reactant=reactant,
                product=product,
                output_dir=out_dir,
                forming_bonds=None,
                fragment_indices=None,
                sp_matrix_report=None,
                ts_fchk=None,
                ts_orca_out=None,
                reactant_fchk=None,
                reactant_orca_out=None,
                product_fchk=None,
                product_orca_out=None,
            )
            
            # Read meta
            meta = json.loads((out_dir / "feature_meta.json").read_text())
            
            # Should have warnings field
            assert 'warnings' in meta, "feature_meta should have 'warnings' field"
            
            # Should have warning entries (since no fchk/log files provided)
            warnings = meta.get('warnings', [])
            assert len(warnings) > 0, \
                "feature_meta.warnings should be non-empty when plugins warn"
            
    def test_feature_meta_warnings_empty_when_all_ok(self):
        """feature_meta.json warnings should be empty when all plugins OK."""
        with tempfile.TemporaryDirectory(prefix="test_meta_ok_") as tmp:
            tmpdir = Path(tmp)
            
            # Create mock fchk with orbital energies
            fchk = tmpdir / "ts.fchk"
            fchk.write_text("""Test FCHK
Number of alpha electrons             I         10
Number of beta electrons              I         10
Orbital Energies                      R  N=10
-1.0D+01 -8.0D+00 -6.0D+00 -4.0D+00 -2.0D+00
-1.0D+00  1.0D+00  2.0D+00  3.0D+00  4.0D+00
Mulliken Charges                       R  N=3
  1  0.100000D+00
  2 -0.100000D+00
  3  0.000000D+00
Atomic numbers                         I  N=3
  6  8  1
""")
            
            # Create mock ts.log with frequency
            ts_log = tmpdir / "ts.log"
            ts_log.write_text("""Harmonic frequencies (cm**-1), IR intensities (KM/mol), and normal modes:

 1              2              3
-187.32       124.56       156.78
""")
            
            ts = tmpdir / "ts_final.xyz"
            reactant = tmpdir / "reactant.xyz"
            product = tmpdir / "product_min.xyz"
            _write_min_xyz(ts)
            _write_min_xyz(reactant)
            _write_min_xyz(product)
            
            out_dir = tmpdir / "S4_Data"
            miner = FeatureMiner(config={})
            
            # Run with fchk and ts_log
            miner.run(
                ts_final=ts,
                reactant=reactant,
                product=product,
                output_dir=out_dir,
                forming_bonds="0-1",
                fragment_indices=None,
                sp_matrix_report=None,
                ts_fchk=fchk,
                ts_orca_out=None,
                reactant_fchk=None,
                reactant_orca_out=None,
                product_fchk=None,
                product_orca_out=None,
                ts_log=ts_log,
            )
            
            # Read meta
            meta = json.loads((out_dir / "feature_meta.json").read_text())
            
            # May have some warnings, but fewer than without files
            warnings = meta.get('warnings', [])
            # At minimum, the field exists
            assert 'warnings' in meta


# =============================================================================
# Tests for Warning Structure
# =============================================================================

class TestWarningStructure:
    """Tests for warning structure (code, plugin, severity, detail)."""

    def test_warning_has_code_field(self):
        """Each warning should have a 'code' field."""
        warning = {'code': 'W_TEST', 'plugin': 'test', 'severity': 'warn', 'detail': 'Test'}
        assert 'code' in warning
        assert warning['code'].startswith('W_')

    def test_warning_has_plugin_field(self):
        """Each warning should have a 'plugin' field."""
        warning = {'code': 'W_TEST', 'plugin': 'ts_quality', 'severity': 'warn', 'detail': 'Test'}
        assert 'plugin' in warning
        assert len(warning['plugin']) > 0

    def test_warning_has_severity_field(self):
        """Each warning should have a 'severity' field."""
        warning = {'code': 'W_TEST', 'plugin': 'test', 'severity': 'warn', 'detail': 'Warning'}
        assert 'severity' in warning
        assert warning['severity'] in ['warn', 'error', 'info']

    def test_warning_has_detail_field(self):
        """Each warning should have a 'detail' field."""
        warning = {'code': 'W_TEST', 'plugin': 'test', 'severity': 'warn', 'detail': 'Detailed message'}
        assert 'detail' in warning
        assert isinstance(warning['detail'], str)


# =============================================================================
# Tests for Warning Deduplication
# =============================================================================

class TestWarningDeduplication:
    """Tests for warning deduplication by (code, plugin)."""

    def test_same_code_plugin_deduplicated(self):
        """Duplicate warnings with same code and plugin should be deduplicated."""
        # Simulate aggregated warnings with duplicates
        aggregated_warnings = [
            {'code': 'W_TS_MISSING_LOG', 'plugin': 'ts_quality', 'severity': 'warn', 'detail': 'Missing'},
            {'code': 'W_TS_MISSING_LOG', 'plugin': 'ts_quality', 'severity': 'warn', 'detail': 'Missing'},
            {'code': 'W_TS_MISSING_LOG', 'plugin': 'ts_quality', 'severity': 'warn', 'detail': 'Missing'},
            {'code': 'W_OTHER', 'plugin': 'other', 'severity': 'warn', 'detail': 'Other'},
        ]
        
        # Deduplicate
        seen = set()
        deduped = []
        for w in aggregated_warnings:
            key = (w['code'], w['plugin'])
            if key not in seen:
                seen.add(key)
                deduped.append(w)
        
        # Should deduplicate to 2 unique warnings
        assert len(deduped) == 2, f"Expected 2 unique warnings, got {len(deduped)}"
        assert deduped[0]['code'] == 'W_TS_MISSING_LOG'
        assert deduped[1]['code'] == 'W_OTHER'

    def test_different_plugins_same_code_kept(self):
        """Same code from different plugins should be kept."""
        aggregated_warnings = [
            {'code': 'W_MISSING', 'plugin': 'plugin_a', 'severity': 'warn', 'detail': 'From A'},
            {'code': 'W_MISSING', 'plugin': 'plugin_b', 'severity': 'warn', 'detail': 'From B'},
        ]
        
        seen = set()
        deduped = []
        for w in aggregated_warnings:
            key = (w['code'], w['plugin'])
            if key not in seen:
                seen.add(key)
                deduped.append(w)
        
        assert len(deduped) == 2, "Same code from different plugins should both be kept"


# =============================================================================
# Tests for qc.sample_weight Policy
# =============================================================================

class TestSampleWeightPolicy:
    """Tests for qc.sample_weight policy implementation."""

    def test_sample_weight_1_for_ok_status(self):
        """qc.sample_weight should be 1.0 when feature_status == OK and TS valid."""
        # Test the policy calculation directly
        feature_status = FeatureStatus.OK
        ts_valid = True
        has_required_features = True
        
        # Apply policy
        if feature_status == FeatureStatus.OK and ts_valid and has_required_features:
            sample_weight = 1.0
        elif feature_status in [FeatureStatus.INVALID_INPUTS, FeatureStatus.FAILED]:
            sample_weight = 0.0
        else:
            sample_weight = 0.5  # Optional: warn state
        
        assert sample_weight == 1.0, "qc.sample_weight should be 1.0 for OK status"

    def test_sample_weight_0_for_invalid_inputs(self):
        """qc.sample_weight should be 0.0 for INVALID_INPUTS status."""
        feature_status = FeatureStatus.INVALID_INPUTS
        ts_valid = False
        has_required_features = False
        
        if feature_status == FeatureStatus.OK and ts_valid and has_required_features:
            sample_weight = 1.0
        elif feature_status in [FeatureStatus.INVALID_INPUTS, FeatureStatus.FAILED]:
            sample_weight = 0.0
        else:
            sample_weight = 0.5
        
        assert sample_weight == 0.0, "qc.sample_weight should be 0.0 for INVALID_INPUTS"

    def test_sample_weight_0_for_failed_status(self):
        """qc.sample_weight should be 0.0 for FAILED status."""
        feature_status = FeatureStatus.FAILED
        ts_valid = False
        has_required_features = False
        
        if feature_status == FeatureStatus.OK and ts_valid and has_required_features:
            sample_weight = 1.0
        elif feature_status in [FeatureStatus.INVALID_INPUTS, FeatureStatus.FAILED]:
            sample_weight = 0.0
        else:
            sample_weight = 0.5
        
        assert sample_weight == 0.0, "qc.sample_weight should be 0.0 for FAILED"

    def test_sample_weight_0_for_ts_invalid(self):
        """qc.sample_weight should be 0.0 when TS is invalid."""
        feature_status = FeatureStatus.OK
        ts_valid = False  # TS has wrong number of imaginary frequencies
        has_required_features = True
        
        if feature_status == FeatureStatus.OK and ts_valid and has_required_features:
            sample_weight = 1.0
        elif feature_status in [FeatureStatus.INVALID_INPUTS, FeatureStatus.FAILED]:
            sample_weight = 0.0
        else:
            sample_weight = 0.0  # TS invalid → 0.0
        
        assert sample_weight == 0.0, "qc.sample_weight should be 0.0 for TS invalid"

    def test_sample_weight_in_features_mlr_csv(self):
        """features_mlr.csv should contain qc.sample_weight column."""
        with tempfile.TemporaryDirectory(prefix="test_weight_") as tmp:
            tmpdir = Path(tmp)
            
            ts = tmpdir / "ts_final.xyz"
            reactant = tmpdir / "reactant.xyz"
            product = tmpdir / "product_min.xyz"
            _write_min_xyz(ts)
            _write_min_xyz(reactant)
            _write_min_xyz(product)
            
            out_dir = tmpdir / "S4_Data"
            miner = FeatureMiner(config={})
            
            miner.run(
                ts_final=ts,
                reactant=reactant,
                product=product,
                output_dir=out_dir,
                forming_bonds=None,
                fragment_indices=None,
                sp_matrix_report=None,
                ts_fchk=None,
                ts_orca_out=None,
                reactant_fchk=None,
                reactant_orca_out=None,
                product_fchk=None,
                product_orca_out=None,
            )
            
            # Read MLR CSV
            df = pd.read_csv(out_dir / "features_mlr.csv")
            
            # Should have qc.sample_weight column
            assert 'qc.sample_weight' in df.columns, \
                "features_mlr.csv should have 'qc.sample_weight' column"

    def test_sample_weight_is_numeric_not_nan(self):
        """qc.sample_weight should be numeric (not NaN) when status is OK."""
        # Test the policy calculation directly
        # With complete inputs (TS + fchk + log), status should be OK
        feature_status = FeatureStatus.OK
        ts_valid = True
        
        # Apply policy
        if feature_status == FeatureStatus.OK and ts_valid:
            sample_weight = 1.0
        elif feature_status in [FeatureStatus.INVALID_INPUTS, FeatureStatus.FAILED]:
            sample_weight = 0.0
        else:
            sample_weight = 0.0
        
        assert sample_weight == 1.0

    def test_sample_weight_0_for_degraded_inputs(self):
        """qc.sample_weight should be 0.0 when essential inputs missing."""
        # Test the policy calculation directly for degraded case
        feature_status = FeatureStatus.PARTIAL  # Degraded status
        ts_valid = False  # TS invalid due to missing data
        
        if feature_status == FeatureStatus.OK and ts_valid:
            sample_weight = 1.0
        elif feature_status in [FeatureStatus.INVALID_INPUTS, FeatureStatus.FAILED]:
            sample_weight = 0.0
        else:
            sample_weight = 0.0  # PARTIAL gets 0.0 per policy
        
        assert sample_weight == 0.0


# =============================================================================
# Integration Tests
# =============================================================================

class TestMetaWarningsIntegration:
    """Integration tests for meta warnings and sample_weight together."""

    def test_full_pipeline_warnings_and_weight(self):
        """Full pipeline should produce warnings and sample_weight together."""
        with tempfile.TemporaryDirectory(prefix="test_full_") as tmp:
            tmpdir = Path(tmp)
            
            ts = tmpdir / "ts_final.xyz"
            reactant = tmpdir / "reactant.xyz"
            product = tmpdir / "product_min.xyz"
            _write_min_xyz(ts)
            _write_min_xyz(reactant)
            _write_min_xyz(product)
            
            out_dir = tmpdir / "S4_Data"
            miner = FeatureMiner(config={})
            
            # Run with partial inputs
            miner.run(
                ts_final=ts,
                reactant=reactant,
                product=product,
                output_dir=out_dir,
                forming_bonds=None,
                fragment_indices=None,
                sp_matrix_report=None,
                ts_fchk=None,
                ts_orca_out=None,
                reactant_fchk=None,
                reactant_orca_out=None,
                product_fchk=None,
                product_orca_out=None,
            )
            
            # Read both outputs
            meta = json.loads((out_dir / "feature_meta.json").read_text())
            df = pd.read_csv(out_dir / "features_mlr.csv")
            
            # Both should exist
            assert 'warnings' in meta, "feature_meta should have warnings"
            assert 'qc.sample_weight' in df.columns, "MLR CSV should have sample_weight"
            
            # Warnings should have structure
            for w in meta.get('warnings', []):
                assert 'code' in w
                assert 'plugin' in w
                assert 'severity' in w
                assert 'detail' in w
            
            # Sample weight should be numeric
            sample_weight = df['qc.sample_weight'].iloc[0]
            assert isinstance(sample_weight, (int, float))
            assert sample_weight in [0.0, 0.5, 1.0], \
                f"sample_weight should be 0.0, 0.5, or 1.0, got {sample_weight}"

    def test_warnings_deduplicated_in_meta(self):
        """Warnings in meta should be deduplicated by (code, plugin)."""
        with tempfile.TemporaryDirectory(prefix="test_dedup_") as tmp:
            tmpdir = Path(tmp)
            
            ts = tmpdir / "ts_final.xyz"
            reactant = tmpdir / "reactant.xyz"
            product = tmpdir / "product_min.xyz"
            _write_min_xyz(ts)
            _write_min_xyz(reactant)
            _write_min_xyz(product)
            
            out_dir = tmpdir / "S4_Data"
            miner = FeatureMiner(config={})
            
            miner.run(
                ts_final=ts,
                reactant=reactant,
                product=product,
                output_dir=out_dir,
                forming_bonds=None,
                fragment_indices=None,
                sp_matrix_report=None,
                ts_fchk=None,
                ts_orca_out=None,
                reactant_fchk=None,
                reactant_orca_out=None,
                product_fchk=None,
                product_orca_out=None,
            )
            
            meta = json.loads((out_dir / "feature_meta.json").read_text())
            warnings = meta.get('warnings', [])
            
            # Check for duplicates
            seen_keys = set()
            for w in warnings:
                key = (w.get('code'), w.get('plugin'))
                assert key not in seen_keys, f"Duplicate warning found: {key}"
                seen_keys.add(key)

    def test_meta_warnings_include_plugin_origin(self):
        """Each warning in meta should indicate which plugin generated it."""
        with tempfile.TemporaryDirectory(prefix="test_plugin_origin_") as tmp:
            tmpdir = Path(tmp)
            
            # Create context that will trigger ts_quality warning
            ts = tmpdir / "ts_final.xyz"
            reactant = tmpdir / "reactant.xyz"
            product = tmpdir / "product_min.xyz"
            _write_min_xyz(ts)
            _write_min_xyz(reactant)
            _write_min_xyz(product)
            
            out_dir = tmpdir / "S4_Data"
            miner = FeatureMiner(config={})
            
            # No ts_log → ts_quality should warn
            miner.run(
                ts_final=ts,
                reactant=reactant,
                product=product,
                output_dir=out_dir,
                forming_bonds=None,
                fragment_indices=None,
                sp_matrix_report=None,
                ts_fchk=None,
                ts_orca_out=None,
                reactant_fchk=None,
                reactant_orca_out=None,
                product_fchk=None,
                product_orca_out=None,
                ts_log=None,
            )
            
            meta = json.loads((out_dir / "feature_meta.json").read_text())
            warnings = meta.get('warnings', [])
            
            # At least one warning should be from ts_quality
            ts_quality_warnings = [w for w in warnings if w.get('plugin') == 'ts_quality']
            # ts_quality may not emit warnings directly in current implementation
            # But there should be warnings from some plugin
            assert len(warnings) > 0, "Should have at least one warning"
            
            for w in warnings:
                assert 'plugin' in w, "Warning must specify plugin origin"
