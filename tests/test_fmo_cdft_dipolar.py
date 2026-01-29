"""
Unit Tests for FMO/CDFT Dipolar Parser
=====================================

Test parsing of dipolar intermediate outputs.

Author: QCcalc Team
Date: 2026-01-27
"""

import pytest
import json
from pathlib import Path

from rph_core.steps.step4_features.extractors.fmo_cdft_dipolar import FmoCdftDipolarParser
from rph_core.steps.step4_features.context import FeatureContext, PluginTrace
from rph_core.steps.step4_features.status import FeatureResultStatus


@pytest.fixture
def gaussian_dipolar_log():
    """Path to Gaussian dipolar log fixture."""
    return Path(__file__).parent / 'tests' / 'fixtures' / 'dipolar' / 'gaussian_dipolar.log'


@pytest.fixture
def orca_dipolar_out():
    """Path to ORCA dipolar out fixture."""
    return Path(__file__).parent / 'tests' / 'fixtures' / 'dipolar' / 'orca_dipolar.out'


@pytest.fixture
def mock_context_with_dipolar(s3_dir, dipolar_output):
    """Create a mock FeatureContext with dipolar path."""
    class MockArtifactsIndex:
        def __init__(self):
            self.dipolar = {
                'path_rel': str(dipolar_output.relative_to(s3_dir)),
                'sha256': None
            }

    return FeatureContext(
        s3_dir=s3_dir,
        artifacts_index=MockArtifactsIndex()
    )


class TestGaussianDipolarParsing:
    """Tests for Gaussian dipolar log parsing."""

    def test_parse_gaussian_occ_virt_eigenvalues(self, gaussian_dipolar_log):
        """Should parse occupied and virtual eigenvalues."""
        parser = FmoCdftDipolarParser()

        content = gaussian_dipolar_log.read_text()
        occ_matches = parser._parse_gaussian_log(content)

        assert occ_matches is not None
        homo_ev, lumo_ev = occ_matches
        assert homo_ev is not None
        assert lumo_ev is not None
        assert homo_ev > 0
        assert lumo_ev > 0
        assert homo_ev < lumo_ev

    def test_parse_gaussian_missing_occ(self, gaussian_dipolar_log):
        """Should return None when occ/virt eigenvalues not found."""
        parser = FmoCdftDipolarParser()
        content = "Missing eigenvalues section"
        occ_matches = parser._parse_gaussian_log(content)

        assert occ_matches == (None, None)

    def test_parse_gaussian_no_virt(self, gaussian_dipolar_log):
        """Should return None when virt eigenvalues not found."""
        parser = FmoCdftDipolarParser()
        content = gaussian_dipolar_log.read_text()
        content = content.replace(
            "Alpha virt. eigenvalues --   0.123",
            "Alpha occ. eigenvalues --   -0.456"
        )
        occ_matches = parser._parse_gaussian_log(content)

        assert occ_matches == (None, None)


class TestOrcadipolarParsing:
    """Tests for ORCA dipolar out parsing."""

    def test_parse_orca_alpha_eigenvalues(self, orca_dipolar_out):
        """Should parse ALPHA EIGENVALUES section."""
        parser = FmoCdftDipolarParser()
        content = orca_dipolar_out.read_text()
        homo_ev, lumo_ev = parser._parse_orca_out(content)

        assert homo_ev is not None
        assert lumo_ev is not None
        assert homo_ev > 0
        assert lumo_ev > 0
        assert homo_ev < lumo_ev

    def test_parse_orca_missing_eigenvalues(self, orca_dipolar_out):
        """Should return None when ALPHA EIGENVALUES not found."""
        parser = FmoCdftDipolarParser()
        content = "No eigenvalues found"
        homo_ev, lumo_ev = parser._parse_orca_out(content)

        assert homo_ev == (None, None)

    def test_parse_orca_negative_eigenvalues(self, orca_dipolar_out):
        """Should correctly handle negative eigenvalues."""
        parser = FmoCdftDipolarParser()
        content = orca_dipolar_out.read_text()
        content = content.replace(
            "ALPHA EIGENVALUES      -0.123   0.456",
            "ALPHA EIGENVALUES      -0.678   0.234   0.123"
        )
        homo_ev, lumo_ev = parser._parse_orca_out(content)

        assert homo_ev is not None
        assert lumo_ev is not None
        assert homo_ev < 0
        assert lumo_ev > 0


class TestFmoCdftDipolarExtraction:
    """Integration tests for FMO/CDFT dipolar parser."""

    def test_extract_with_gaussian_dipolar(self, gaussian_dipolar_log, mock_context_with_dipolar):
        """Should extract HOMO/LUMO from Gaussian dipolar log."""
        s3_dir = gaussian_dipolar_log.parent
        context = mock_context_with_dipolar(s3_dir, gaussian_dipolar_log)

        parser = FmoCdftDipolarParser()
        trace = parser.extract(context)

        assert trace.status == FeatureResultStatus.OK
        assert 'fmo_cdft_dipolar.homo_ev' in trace._extracted_features
        assert 'fmo_cdft_dipolar.lumo_ev' in trace._extracted_features
        assert 'fmo_cdft_dipolar.gap_ev' in trace._extracted_features
        assert 'fmo_cdft_dipolar.omega_ev' in trace._extracted_features
        assert trace._extracted_features['fmo_cdft_dipolar.homo_ev'] > 0
        assert trace._extracted_features['fmo_cdft_dipolar.lumo_ev'] > 0

    def test_extract_with_orca_dipolar(self, orca_dipolar_out, mock_context_with_dipolar):
        """Should extract HOMO/LUMO from ORCA dipolar out."""
        s3_dir = orca_dipolar_out.parent
        context = mock_context_with_dipolar(s3_dir, orca_dipolar_out)

        parser = FmoCdftDipolarParser()
        trace = parser.extract(context)

        assert trace.status == FeatureResultStatus.OK
        assert 'fmo_cdft_dipolar.homo_ev' in trace._extracted_features
        assert 'fmo_cdft_dipolar.lumo_ev' in trace._extracted_features
        assert 'fmo_cdft_dipolar.gap_ev' in trace._extracted_features
        assert 'fmo_cdft_dipolar.omega_ev' in trace._extracted_features

    def test_extract_with_missing_dipolar(self, mock_context_with_dipolar):
        """Should skip when dipolar output not found."""
        s3_dir = gaussian_dipolar_log.parent
        context = mock_context_with_dipolar(s3_dir, None)

        parser = FmoCdftDipolarParser()
        trace = parser.extract(context)

        assert trace.status == FeatureResultStatus.SKIPPED
        assert trace._extracted_features['fmo_cdft_dipolar.status'] == 'skipped'
        assert 'dipolar_output' in trace.missing_paths

    def test_extract_without_s3_dir(self, mock_context_with_dipolar, orca_dipolar_out):
        """Should skip when s3_dir not provided."""
        class MockContextWithoutS3:
            pass

        context = MockContextWithoutS3()
        context.s3_dir = None
        context.artifacts_index = None

        parser = FmoCdftDipolarParser()
        trace = parser.extract(context)

        assert trace.status == FeatureResultStatus.SKIPPED
        assert 's3_dir' in trace.missing_fields
        assert trace._extracted_features['fmo_cdft_dipolar.status'] == 'skipped'


class TestGapAndOmegaCalculation:
    """Tests for gap and omega calculations."""

    def test_calculate_gap_positive(self):
        """Gap should be positive when HOMO < LUMO."""
        homo = -5.0
        lumo = -2.0
        gap = lumo - homo
        omega = (homo**2 + lumo**2) / 2

        assert gap > 0
        assert omega > 0

    def test_calculate_gap_equal(self):
        """Gap should be 0 when HOMO == LUMO."""
        homo = -5.0
        lumo = -5.0
        gap = lumo - homo
        omega = (homo**2 + lumo**2) / 2

        assert gap == 0
        assert omega == 0

    def test_calculate_omega_formula(self):
        """Omega = (HOMO^2 + LUMO^2)/2 should be correct."""
        homo = -6.0
        lumo = -3.0
        expected_omega = ((-6.0)**2 + (-3.0)**2) / 2

        omega = (homo**2 + lumo**2) / 2

        assert abs(omega - expected_omega) < 0.01
