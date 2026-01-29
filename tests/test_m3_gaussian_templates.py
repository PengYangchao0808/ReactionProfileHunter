"""
M3-3-3: Pure rendering tests for Gaussian templates
=========================================================

Tests that Gaussian templates are properly rendered with:
- {{charge}} placeholder
- {{mult}} placeholder
- {{coords_block}} placeholder
- No external QC calls (pure template rendering)

Author: QC Descriptors Team
Date: 2026-01-21
"""

import pytest
from pathlib import Path


class TestGaussianTemplateRendering:
    """Test M3-3 Gaussian template rendering (pure, no QC calls)."""

    @pytest.fixture
    def sample_atoms(self):
        """Sample atomic coordinates."""
        return [
            {"symbol": "C", "x": 0.0, "y": 0.0, "z": 0.0},
            {"symbol": "H", "x": 1.0, "y": 0.0, "z": 0.0},
            {"symbol": "O", "x": 0.0, "y": 1.2, "z": 0.0},
            {"symbol": "H", "x": -1.0, "y": 0.0, "z": 0.0}
        ]

    @pytest.fixture
    def coords_block(self, sample_atoms):
        """Render coordinate block from atoms."""
        lines = []
        for atom in sample_atoms:
            lines.append(f"{atom['symbol']:<2} {atom['x']:14.8f} {atom['y']:14.8f} {atom['z']:14.8f}")
        return "\n".join(lines)

    def test_gaussian_ts_template_rendering(self, sample_atoms, coords_block):
        """TS template should include required placeholders and keywords."""
        project_root = Path(__file__).resolve().parents[1]
        template_path = project_root / "config/templates/gaussian_ts.gjf"
        template_str = template_path.read_text()

        assert "{{charge}}" in template_str
        assert "{{mult}}" in template_str
        assert "{{coords_block}}" in template_str
        assert "%chk=ts.chk" in template_str
        assert "Opt=(TS, CalcFC, NoEigenTest)" in template_str

    def test_gaussian_freq_template_rendering(self, sample_atoms, coords_block):
        """Freq template should include required placeholders and keywords."""
        project_root = Path(__file__).resolve().parents[1]
        template_path = project_root / "config/templates/gaussian_freq.gjf"
        template_str = template_path.read_text()

        assert "{{mem}}" in template_str
        assert "{{nproc}}" in template_str
        assert "{{method}}" in template_str
        assert "{{basis}}" in template_str
        assert "{{route}}" in template_str
        assert "{{charge}}" in template_str
        assert "{{mult}}" in template_str
        assert "{{coords_block}}" in template_str
        assert "%mem=" in template_str
        assert "#p" in template_str

    def test_gaussian_nmr_template_rendering(self, sample_atoms, coords_block):
        """v5.4 removes NMR template end-to-end."""
        project_root = Path(__file__).resolve().parents[1]
        template_path = project_root / "config/templates/gaussian_nmr.gjf"
        assert not template_path.exists()

    def test_gaussian_nbo_template_rendering(self, sample_atoms, coords_block):
        """NBO template should include required placeholders and keywords."""
        project_root = Path(__file__).resolve().parents[1]
        template_path = project_root / "config/templates/gaussian_nbo.gjf"
        template_str = template_path.read_text()

        assert "{{mem}}" in template_str
        assert "{{nproc}}" in template_str
        assert "{{method}}" in template_str
        assert "{{basis}}" in template_str
        assert "{{charge}}" in template_str
        assert "{{mult}}" in template_str
        assert "{{coords_block}}" in template_str
        assert "Pop=NBO" in template_str
        assert "NMR=GIAO" not in template_str


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
