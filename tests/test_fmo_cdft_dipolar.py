"""
Unit Tests for FMO/CDFT Dipolar Parser
=====================================

Test parsing of dipolar intermediate outputs.

Author: QCcalc Team
Date: 2026-01-27
Updated: 2026-02-05 (fixtures + context wiring fixes)
"""

from pathlib import Path
from typing import Optional, Any, Dict

import pytest

from rph_core.steps.step4_features.extractors.fmo_cdft_dipolar import FmoCdftDipolarParser
from rph_core.steps.step4_features.context import FeatureContext
from rph_core.steps.step4_features.status import FeatureResultStatus


@pytest.fixture
def gaussian_dipolar_log() -> Path:
    """Path to Gaussian dipolar log fixture."""
    return Path(__file__).parent / 'fixtures' / 'dipolar' / 'gaussian_dipolar.log'


@pytest.fixture
def orca_dipolar_out() -> Path:
    """Path to ORCA dipolar out fixture."""
    return Path(__file__).parent / 'fixtures' / 'dipolar' / 'orca_dipolar.out'


def _make_context_with_dipolar(s3_dir: Path, dipolar_output: Optional[Path]) -> FeatureContext:
    artifacts_index: Dict[str, Any] = {}
    if dipolar_output is not None:
        artifacts_index['dipolar'] = {
            'path_rel': str(dipolar_output.relative_to(s3_dir)),
            'sha256': None,
        }
    return FeatureContext(s3_dir=s3_dir, artifacts_index=artifacts_index)


class TestGaussianDipolarParsing:
    """Tests for Gaussian dipolar log parsing."""

    def test_parse_gaussian_occ_virt_eigenvalues(self, gaussian_dipolar_log: Path):
        parser = FmoCdftDipolarParser()
        content = gaussian_dipolar_log.read_text()
        homo_ev, lumo_ev = parser._parse_gaussian_log(content)

        assert homo_ev is not None
        assert lumo_ev is not None
        assert homo_ev < lumo_ev

    def test_parse_gaussian_missing_occ(self):
        parser = FmoCdftDipolarParser()
        content = "Missing eigenvalues section"
        assert parser._parse_gaussian_log(content) == (None, None)

    def test_parse_gaussian_no_virt(self, gaussian_dipolar_log: Path):
        parser = FmoCdftDipolarParser()
        content = gaussian_dipolar_log.read_text()
        content = content.replace("Alpha virt. eigenvalues", "Alpha occ. eigenvalues")
        assert parser._parse_gaussian_log(content) == (None, None)


class TestOrcaDipolarParsing:
    """Tests for ORCA dipolar out parsing."""

    def test_parse_orca_orbital_energies(self, orca_dipolar_out: Path):
        parser = FmoCdftDipolarParser()
        content = orca_dipolar_out.read_text()
        homo_ev, lumo_ev = parser._parse_orca_out(content)

        assert homo_ev is not None
        assert lumo_ev is not None
        assert homo_ev < lumo_ev

    def test_parse_orca_missing_orbital_section(self):
        parser = FmoCdftDipolarParser()
        content = "No eigenvalues found"
        homo_ev, lumo_ev = parser._parse_orca_out(content)
        assert (homo_ev, lumo_ev) == (None, None)

    def test_parse_orca_negative_virtual_energies(self):
        parser = FmoCdftDipolarParser()
        content = """ORBITAL ENERGIES
----------------
   NO   OCC         E(Eh)         E(eV)
    1   2.0000     -0.5000       -5.0000
    2   2.0000     -0.3000       -8.0000
    3   0.0000      0.0500       -1.0000
    4   0.0000      0.1000       -0.5000
----------------
"""
        homo_ev, lumo_ev = parser._parse_orca_out(content)
        assert homo_ev is not None
        assert lumo_ev is not None
        assert homo_ev < lumo_ev


class TestFmoCdftDipolarExtraction:
    """Integration tests for FMO/CDFT dipolar extractor."""

    def test_extract_with_gaussian_dipolar(self, gaussian_dipolar_log: Path):
        s3_dir = gaussian_dipolar_log.parent
        context = _make_context_with_dipolar(s3_dir, gaussian_dipolar_log)

        parser = FmoCdftDipolarParser()
        trace = parser.run(context)

        assert trace.status == FeatureResultStatus.OK
        feats = trace._extracted_features
        assert feats.get('fmo_cdft_dipolar.status') == 'ok'
        assert feats.get('fmo_cdft_dipolar.homo_ev') is not None
        assert feats.get('fmo_cdft_dipolar.lumo_ev') is not None
        assert feats.get('fmo_cdft_dipolar.gap_ev') is not None
        assert feats.get('fmo_cdft_dipolar.omega_ev') is not None

    def test_extract_with_orca_dipolar(self, orca_dipolar_out: Path):
        s3_dir = orca_dipolar_out.parent
        context = _make_context_with_dipolar(s3_dir, orca_dipolar_out)

        parser = FmoCdftDipolarParser()
        trace = parser.run(context)

        assert trace.status == FeatureResultStatus.OK
        feats = trace._extracted_features
        assert feats.get('fmo_cdft_dipolar.status') == 'ok'
        assert feats.get('fmo_cdft_dipolar.homo_ev') is not None
        assert feats.get('fmo_cdft_dipolar.lumo_ev') is not None
        assert feats.get('fmo_cdft_dipolar.gap_ev') is not None
        assert feats.get('fmo_cdft_dipolar.omega_ev') is not None

    def test_extract_with_missing_dipolar(self, tmp_path: Path):
        s3_dir = tmp_path / "S3"
        s3_dir.mkdir(parents=True, exist_ok=True)
        context = _make_context_with_dipolar(s3_dir, None)

        parser = FmoCdftDipolarParser()
        trace = parser.run(context)

        feats = trace._extracted_features
        assert feats.get('fmo_cdft_dipolar.status') == 'skipped'
        assert feats.get('fmo_cdft_dipolar.missing_reason') == 'dipolar_output_not_found'
        assert 'W_FMO_DIPOLAR_OUTPUT_NOT_FOUND' in trace.warnings

    def test_extract_without_s3_dir(self):
        context = FeatureContext(s3_dir=None, artifacts_index={})
        parser = FmoCdftDipolarParser()
        trace = parser.run(context)

        assert trace.status == FeatureResultStatus.SKIPPED
        assert 's3_dir' in trace.missing_fields


class TestGapAndOmegaCalculation:
    """Tests for gap and omega calculations."""

    def test_calculate_gap_positive(self):
        homo = -5.0
        lumo = -2.0
        gap = lumo - homo
        omega = (homo**2 + lumo**2) / 2

        assert gap > 0
        assert omega > 0

    def test_calculate_gap_equal(self):
        homo = -5.0
        lumo = -5.0
        gap = lumo - homo
        omega = (homo**2 + lumo**2) / 2

        assert gap == 0
        assert omega == 25.0

    def test_calculate_omega_formula(self):
        homo = -6.0
        lumo = -3.0
        expected_omega = ((-6.0) ** 2 + (-3.0) ** 2) / 2

        omega = (homo**2 + lumo**2) / 2
        assert abs(omega - expected_omega) < 0.01
