"""
TDD Tests for S4 Extractor Degradation Behavior
================================================

Tests to verify that extractors:
1. Emit keys even when inputs are missing (instead of being SKIPPED)
2. Record warning codes instead of message strings
3. Degrade gracefully with NaN values when files are missing

These tests should FAIL initially (RED) until extractors are fixed.

Author: RPH Team
Date: 2026-02-03
"""

import pytest
import tempfile
import json
from pathlib import Path
from unittest.mock import MagicMock, patch
from typing import Optional, Any, Dict, Union

from rph_core.steps.step4_features.context import FeatureContext, PluginTrace
from rph_core.steps.step4_features.status import FeatureResultStatus
from rph_core.steps.step4_features.feature_miner import FeatureMiner


# =============================================================================
# Minimal Fixtures
# =============================================================================

def _write_min_xyz(path: Path, content: Optional[str] = None) -> None:
    """Write minimal xyz file."""
    if content is None:
        content = "3\n\nC 0.0 0.0 0.0\nO 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
    path.write_text(content)


def _create_minimal_context(
    tmpdir: Path,
    ts_final: Optional[Union[Path, str]] = None,
    reactant: Optional[Path] = None,
    product: Optional[Path] = None,
    ts_fchk: Optional[Path] = None,
    reactant_fchk: Optional[Path] = None,
    ts_log: Optional[Path] = None,
    ts_orca_out: Optional[Path] = None,
    s1_shermo_summary: Optional[Path] = None,
    s1_hoac_thermo: Optional[Path] = None,
    forming_bonds = None,
) -> FeatureContext:
    """Create minimal FeatureContext with optional paths."""
    
    # Create minimal xyz files
    if ts_final is None:
        ts_final = tmpdir / "ts_final.xyz"
        _write_min_xyz(ts_final)
    elif isinstance(ts_final, str):
        ts_final = tmpdir / ts_final
        _write_min_xyz(ts_final)
    
    if reactant is None:
        reactant = tmpdir / "reactant.xyz"
        _write_min_xyz(reactant)
    
    if product is None:
        product = tmpdir / "product_min.xyz"
        _write_min_xyz(product)
    
    artifacts_index: Dict[str, Any] = {
        'ts_final': {'path_rel': str(Path(ts_final).relative_to(tmpdir)), 'sha256': None},
        'reactant': {'path_rel': str(Path(reactant).relative_to(tmpdir)), 'sha256': None},
        'product': {'path_rel': str(Path(product).relative_to(tmpdir)), 'sha256': None},
    }

    if ts_fchk is not None:
        artifacts_index['ts_fchk'] = {'path_rel': str(ts_fchk.relative_to(tmpdir)), 'sha256': None}
    if reactant_fchk is not None:
        artifacts_index['reactant_fchk'] = {'path_rel': str(reactant_fchk.relative_to(tmpdir)), 'sha256': None}
    if s1_shermo_summary is not None:
        artifacts_index['s1_shermo_summary_file'] = {'path_rel': str(s1_shermo_summary.relative_to(tmpdir)), 'sha256': None}
    if s1_hoac_thermo is not None:
        artifacts_index['s1_hoac_thermo_file'] = {'path_rel': str(s1_hoac_thermo.relative_to(tmpdir)), 'sha256': None}
    
    # Create context
    context = FeatureContext(
        s3_dir=tmpdir,
        artifacts_index=artifacts_index,
        forming_bonds=forming_bonds,
        ts_xyz=Path(ts_final),
        reactant_xyz=reactant,
        product_xyz=product,
        ts_fchk=ts_fchk,
        reactant_fchk=reactant_fchk,
        ts_orca_out=ts_orca_out,
        ts_log=ts_log,
        s3_ts_log=ts_log,
        s3_ts_fchk=ts_fchk,
        s3_reactant_fchk=reactant_fchk,
        s1_shermo_summary_file=s1_shermo_summary,
        s1_hoac_thermo_file=s1_hoac_thermo,
    )
    
    return context


# =============================================================================
# Tests for ts_quality Extractor
# =============================================================================

class TestTSQualityDegrade:
    """Tests for ts_quality extractor degradation behavior."""

    def test_ts_quality_emits_keys_when_ts_log_missing(self):
        """ts_quality should emit keys (NaN) when ts_log is missing, not SKIP."""
        from rph_core.steps.step4_features.extractors.ts_quality import TSQualityExtractor
        
        with tempfile.TemporaryDirectory(prefix="test_ts_quality_") as tmp:
            tmpdir = Path(tmp)
            context = _create_minimal_context(tmpdir)
            
            extractor = TSQualityExtractor()
            trace = extractor.run(context)  # Use run() not extract()
            
            # Should NOT be SKIPPED due to missing ts_log
            # Should emit keys with NaN values and record warning code
            assert trace.status != FeatureResultStatus.SKIPPED, \
                "ts_quality should not be SKIPPED when ts_log is missing"
            
            # Should emit required keys
            required_keys = ['ts.n_imag', 'ts.imag1_cm1_abs', 'ts.dipole_debye']
            for key in required_keys:
                assert key in trace._extracted_features, \
                    f"ts_quality should emit '{key}' even when ts_log missing"
            
            # ts.n_imag should be 0 (no imaginary frequencies found)
            assert trace._extracted_features.get('ts.n_imag') == 0
            # ts.imag1_cm1_abs should be NaN (no TS log to parse)
            import numpy as np
            assert np.isnan(trace._extracted_features.get('ts.imag1_cm1_abs', float('nan')))
            
            # Should have recorded a warning code
            assert len(trace.warnings) > 0, "Should record warning code when ts_log missing"

    def test_ts_quality_parses_frequency_from_ts_log(self):
        """ts_quality should parse frequency when ts_log is available."""
        from rph_core.steps.step4_features.extractors.ts_quality import TSQualityExtractor
        
        with tempfile.TemporaryDirectory(prefix="test_ts_quality_") as tmp:
            tmpdir = Path(tmp)
            
            # Create mock ts.log with frequency - use actual Gaussian format
            ts_log = tmpdir / "ts.log"
            ts_log.write_text("""Frequencies --  -187.32  -152.47   -89.23   124.56   156.78   234.12
""")
            
            context = _create_minimal_context(tmpdir, ts_log=ts_log)
            
            extractor = TSQualityExtractor()
            trace = extractor.run(context)  # Use run() not extract()
            
            # Should extract frequency
            assert trace._extracted_features.get('ts.n_imag') == 3, \
                "Should detect 3 negative frequencies"
            assert trace._extracted_features.get('ts.imag1_cm1_abs') == 187.32, \
                "Should extract first imaginary frequency value"

    def test_ts_quality_fallback_to_s3_ts_log(self):
        """ts_quality should fallback to s3_ts_log if ts_log missing."""
        from rph_core.steps.step4_features.extractors.ts_quality import TSQualityExtractor
        
        with tempfile.TemporaryDirectory(prefix="test_ts_quality_") as tmp:
            tmpdir = Path(tmp)
            
            # Create mock ts.log in S3 location - use actual Gaussian format
            s3_ts_log = tmpdir / "ts.log"
            s3_ts_log.write_text("""Frequencies --  -256.78   145.23   289.45
""")
            
            # Create context without ts_log but with s3_dir
            context = _create_minimal_context(tmpdir)
            context.s3_ts_log = s3_ts_log
            
            extractor = TSQualityExtractor()
            trace = extractor.run(context)  # Use run() not extract()
            
            # Should still extract frequency via fallback (1 negative frequency)
            assert trace._extracted_features.get('ts.n_imag') == 1, \
                "Should detect 1 imaginary frequency via s3_ts_log fallback"


# =============================================================================
# Tests for step1_activation Extractor
# =============================================================================

class TestStep1ActivationDegrade:
    """Tests for step1_activation extractor degradation behavior."""

    def test_step1_activation_not_skipped_when_shermo_missing(self):
        """step1_activation should not SKIP when shermo_summary is missing."""
        from rph_core.steps.step4_features.extractors.step1_activation import Step1ActivationExtractor
        
        with tempfile.TemporaryDirectory(prefix="test_step1_") as tmp:
            tmpdir = Path(tmp)
            context = _create_minimal_context(tmpdir)
            
            extractor = Step1ActivationExtractor()
            trace = extractor.run(context)  # Use run() not extract() to get PluginTrace with status
            
            # Should NOT be SKIPPED due to missing shermo_summary
            assert trace.status != FeatureResultStatus.SKIPPED, \
                "step1_activation should not be SKIPPED when shermo_summary missing"
            
            # Should emit s1_* keys (with NaN)
            s1_keys = [k for k in trace._extracted_features.keys() if k.startswith('s1_')]
            assert len(s1_keys) > 0, "step1_activation should emit s1_* keys"
            
            # Should record warning code for missing file
            assert len(trace.warnings) > 0, \
                "Should record warning code for missing shermo_summary"

    def test_step1_activation_not_skipped_when_hoac_thermo_missing(self):
        """step1_activation should not SKIP when hoac_thermo is missing."""
        from rph_core.steps.step4_features.extractors.step1_activation import Step1ActivationExtractor
        
        with tempfile.TemporaryDirectory(prefix="test_step1_") as tmp:
            tmpdir = Path(tmp)
            context = _create_minimal_context(tmpdir)
            
            extractor = Step1ActivationExtractor()
            trace = extractor.run(context)  # Use run() not extract() to get PluginTrace with status
            
            # Should NOT be SKIPPED due to missing hoac_thermo
            assert trace.status != FeatureResultStatus.SKIPPED, \
                "step1_activation should not be SKIPPED when hoac_thermo missing"
            
            # Should emit s1_* keys and record warning
            s1_keys = [k for k in trace._extracted_features.keys() if k.startswith('s1_')]
            assert len(s1_keys) > 0, "step1_activation should emit s1_* keys"

    def test_step1_activation_runs_with_complete_inputs(self):
        """step1_activation should calculate dG_act when all inputs present.
        
        New formula (V6.2): dG_act = G(S3_reactant) + G(leaving_group) - G(precursor)
        Where S3_reactant Gibbs is provided by Shermo (SP + freq) output.
        """
        from rph_core.steps.step4_features.extractors.step1_activation import Step1ActivationExtractor
        
        with tempfile.TemporaryDirectory(prefix="test_step1_") as tmp:
            tmpdir = Path(tmp)
            
            # Create mock S1 dir structure with precursor Shermo summary
            s1_dir = tmpdir / "S1_ConfGeneration"
            s1_dir.mkdir(parents=True)
            
            s1_shermo = s1_dir / "shermo_summary.json"
            s1_shermo.write_text(json.dumps({
                "g_precursor": -500.0,
                "temperature_K": 298.15,
                "unit": "kcal/mol"
            }))
            
            # Create mock HOAc (leaving group) thermo.json
            hoac_dir = s1_dir / "small_molecules" / "HOAc"
            hoac_dir.mkdir(parents=True)
            s1_hoac = hoac_dir / "thermo.json"
            s1_hoac.write_text(json.dumps({
                "g": -200.0,
                "temperature_K": 298.15,
                "unit": "kcal/mol"
            }))
            
            context = _create_minimal_context(
                tmpdir, 
                s1_shermo_summary=s1_shermo,
                s1_hoac_thermo=s1_hoac
            )
            context.s1_dir = s1_dir
            context.s3_reactant_g_kcal = -450.0
            
            extractor = Step1ActivationExtractor()
            trace = extractor.run(context)
            
            # Should calculate s1_dG_act
            s1_dg = trace._extracted_features.get('s1_dG_act')
            assert s1_dg is not None, "Should calculate s1_dG_act"
            # Expected: G(S3_reactant) + G(leaving_group) - G(precursor)
            # = -450.0 + (-200.0) - (-500.0) = -150.0 kcal/mol
            assert abs(s1_dg - (-150.0)) < 1.0, \
                f"s1_dG_act should be ~-150.0, got {s1_dg}"


# =============================================================================
# Tests for step2_cyclization Extractor
# =============================================================================

class TestStep2CyclizationDegrade:
    """Tests for step2_cyclization extractor degradation behavior."""

    def test_step2_cyclization_not_skipped_when_fchk_missing(self):
        """step2_cyclization should not SKIP when fchk files are missing."""
        from rph_core.steps.step4_features.extractors.step2_cyclization import Step2CyclizationExtractor
        
        with tempfile.TemporaryDirectory(prefix="test_step2_") as tmp:
            tmpdir = Path(tmp)
            context = _create_minimal_context(tmpdir, forming_bonds="0-1,2-3")
            
            extractor = Step2CyclizationExtractor()
            trace = extractor.run(context)  # Use run() not extract()
            
            # Should NOT be SKIPPED due to missing fchk
            assert trace.status != FeatureResultStatus.SKIPPED, \
                "step2_cyclization should not be SKIPPED when fchk missing"

            # De-dup contract: Step2 extractor must NOT emit TS-geometry or
            # TS-validity keys (owned by geometry/ts_quality extractors).
            assert 's2_d_forming_1' not in trace._extracted_features
            assert 's2_d_forming_2' not in trace._extracted_features
            assert 's2_asynch' not in trace._extracted_features
            assert 's2_n_imag_freq' not in trace._extracted_features
            assert 's2_imag_freq_cm1' not in trace._extracted_features
            assert 's2_ts_validity_flag' not in trace._extracted_features
            
            # Should have degraded CDFT/GEDT (NaN)
            import numpy as np
            cdft_keys = ['s2_eps_homo', 's2_eps_lumo', 's2_mu', 's2_eta', 's2_omega']
            for key in cdft_keys:
                val = trace._extracted_features.get(key)
                assert np.isnan(val) if val is not None else True, \
                    f"{key} should be NaN when fchk missing"

            gedt_val = trace._extracted_features.get('s2_gedt_value')
            assert np.isnan(gedt_val) if gedt_val is not None else True, \
                "s2_gedt_value should be NaN when charges/coords missing"

    def test_step2_cyclization_uses_ts_fchk_fallback(self):
        """step2_cyclization should use ts_fchk when reactant_fchk missing."""
        from rph_core.steps.step4_features.extractors.step2_cyclization import Step2CyclizationExtractor
        
        with tempfile.TemporaryDirectory(prefix="test_step2_") as tmp:
            tmpdir = Path(tmp)
            
            # Create mock ts.fchk
            ts_fchk = tmpdir / "ts.fchk"
            ts_fchk.write_text("""Test FCHK
Number of alpha electrons             I         10
Number of beta electrons              I         10
Orbital Energies                      R  N=10
-1.0D+01 -8.0D+00 -6.0D+00 -4.0D+00 -2.0D+00
-1.0D+00  1.0D+00  2.0D+00  3.0D+00  4.0D+00
Mulliken Charges                       R  N=4
  1  0.100000D+00
  2 -0.100000D+00
  3  0.050000D+00
  4 -0.050000D+00
Atomic numbers                         I  N=4
  6  8  1  1
""")
            
            context = _create_minimal_context(
                tmpdir, 
                ts_fchk=ts_fchk,
                forming_bonds="0-1"
            )
            extractor = Step2CyclizationExtractor()
            trace = extractor.run(context)  # Use run() not extract()
            
            # Should use ts_fchk as fallback
            # Should have parsed orbital energies
            assert trace._extracted_features.get('s2_eps_homo') is not None, \
                "Should parse eps_homo from ts_fchk fallback"


# =============================================================================
# Tests for BaseExtractor Gating Behavior
# =============================================================================

class TestBaseExtractorGating:
    """Tests for BaseExtractor.validate_inputs() behavior."""

    def test_base_extractor_respects_optional_inputs(self):
        """BaseExtractor should not skip for missing optional inputs."""
        from rph_core.steps.step4_features.extractors.base import BaseExtractor
        
        class MinimalExtractor(BaseExtractor):
            def get_plugin_name(self):
                return "minimal"
            
            def get_required_inputs(self):
                # Return empty list - all inputs optional
                return []
            
            def extract(self, context):
                return {'test_key': 42}
        
        with tempfile.TemporaryDirectory(prefix="test_base_") as tmp:
            tmpdir = Path(tmp)
            context = _create_minimal_context(tmpdir)
            
            extractor = MinimalExtractor()
            trace = extractor.run(context)  # Use run() to get PluginTrace
            
            # Should run (not SKIP)
            assert trace.status != FeatureResultStatus.SKIPPED
            assert trace._extracted_features.get('test_key') == 42

    def test_extractor_with_all_optional_files_still_runs(self):
        """Extractors should run when all files are technically optional."""
        from rph_core.steps.step4_features.extractors.base import BaseExtractor
        
        class TestExtractor(BaseExtractor):
            def get_plugin_name(self):
                return "test"
            
            def get_required_inputs(self):
                # No strictly required inputs
                return []
            
            def extract(self, context):
                return {'degraded': True}
        
        with tempfile.TemporaryDirectory(prefix="test_optional_") as tmp:
            tmpdir = Path(tmp)
            context = _create_minimal_context(tmpdir)
            
            extractor = TestExtractor()
            trace = extractor.run(context)  # Use run() to get PluginTrace
            
            # Should run and emit keys
            assert trace.status != FeatureResultStatus.SKIPPED
            assert trace._extracted_features.get('degraded') == True


# =============================================================================
# Tests for Warning Codes (Not Message Strings)
# =============================================================================

class TestWarningCodesFormat:
    """Tests to verify warnings use codes, not message strings."""

    def test_ts_quality_warning_uses_code(self):
        """ts_quality should use warning codes like 'W_TS_MISSING_LOG'."""
        from rph_core.steps.step4_features.extractors.ts_quality import TSQualityExtractor
        
        with tempfile.TemporaryDirectory(prefix="test_warning_") as tmp:
            tmpdir = Path(tmp)
            context = _create_minimal_context(tmpdir)
            
            extractor = TSQualityExtractor()
            trace = extractor.run(context)  # Use run() to get PluginTrace with warnings

            # Check warnings are codes, not full message strings
            for w in trace.warnings:
                # Should be short code, not long message
                assert isinstance(w, str), "Warning should be a string"
                assert len(w) < 50, f"Warning should be code, not message: {w}"
                # Should start with 'W_' prefix
                assert w.startswith('W_'), f"Warning should start with 'W_': {w}"

    def test_step1_activation_warning_uses_code(self):
        """step1_activation should use warning codes."""
        from rph_core.steps.step4_features.extractors.step1_activation import Step1ActivationExtractor
        
        with tempfile.TemporaryDirectory(prefix="test_warning_") as tmp:
            tmpdir = Path(tmp)
            context = _create_minimal_context(tmpdir)
            
            extractor = Step1ActivationExtractor()
            trace = extractor.run(context)  # Use run() not extract() to get PluginTrace with status
            
            # Warnings should be codes
            for w in trace.warnings:
                assert w.startswith('W_'), f"Warning should start with 'W_': {w}"

    def test_step2_cyclization_warning_uses_code(self):
        """step2_cyclization should use warning codes for degraded features."""
        from rph_core.steps.step4_features.extractors.step2_cyclization import Step2CyclizationExtractor
        
        with tempfile.TemporaryDirectory(prefix="test_warning_") as tmp:
            tmpdir = Path(tmp)
            context = _create_minimal_context(tmpdir, forming_bonds="0-1")
            
            extractor = Step2CyclizationExtractor()
            trace = extractor.run(context)  # Use run() not extract()
            
            # Warnings should be codes
            for w in trace.warnings:
                assert w.startswith('W_'), f"Warning should start with 'W_': {w}"


# =============================================================================
# Integration Tests
# =============================================================================

class TestFeatureMinerIntegration:
    """Integration tests for FeatureMiner with degraded extractors."""

    def test_featureminer_produces_output_with_missing_artifacts(self):
        """FeatureMiner should produce CSV/JSON even with missing artifacts."""
        with tempfile.TemporaryDirectory(prefix="test_miner_") as tmp:
            tmpdir = Path(tmp)
            
            ts = tmpdir / "ts_final.xyz"
            reactant = tmpdir / "reactant.xyz"
            product = tmpdir / "product_min.xyz"
            _write_min_xyz(ts)
            _write_min_xyz(reactant)
            _write_min_xyz(product)
            
            out_dir = tmpdir / "S4_Data"
            miner = FeatureMiner(config={})
            
            # Run with NO optional files
            features_raw_csv = miner.run(
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
            
            # Should produce output files
            assert features_raw_csv.exists()
            assert (out_dir / "features_mlr.csv").exists()
            assert (out_dir / "feature_meta.json").exists()

    def test_featureminer_csv_contains_ts_keys_when_log_missing(self):
        """features_mlr.csv should contain ts.* keys (NaN) when log missing."""
        with tempfile.TemporaryDirectory(prefix="test_miner_") as tmp:
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
            
            # Read CSV and check for ts keys
            import pandas as pd
            df = pd.read_csv(out_dir / "features_mlr.csv")
            
            # Should have ts columns (even if NaN)
            ts_columns = [c for c in df.columns if c.startswith('ts.')]
            assert len(ts_columns) > 0, "CSV should contain ts.* columns"

    def test_featureminer_emits_warning_codes_in_meta(self):
        """feature_meta.json should contain warning codes from extractors."""
        with tempfile.TemporaryDirectory(prefix="test_miner_") as tmp:
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
            
            # Read meta and check for warnings
            meta = json.loads((out_dir / "feature_meta.json").read_text())
            
            # Should have warnings with codes
            assert 'warnings' in meta, "feature_meta should have warnings field"
            warnings = meta.get('warnings', [])
            assert len(warnings) > 0, "Should have warning entries"
            
            # Each warning should have code
            for w in warnings:
                assert 'code' in w, "Warning should have 'code' field"
                assert w['code'].startswith('W_'), "Warning code should start with 'W_'"
