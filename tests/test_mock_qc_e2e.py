"""
M4-D-5: End-to-End Mock Tests for QC Artifacts Integration
===============================================================

Tests full workflow from FakeBackend to pack_mechanism_assets() to mech_index.json
without invoking real QC binaries.

V5.4 Update: Removed Hirshfeld testing (NMR/Hirshfeld removed to reduce complexity)

Author: QC Descriptors Team
Date: 2026-01-21
"""

import pytest
import tempfile
import json
from pathlib import Path

from rph_core.utils.qc_interface import FakeBackend
from rph_core.steps.step4_features.mech_packager import pack_mechanism_assets


class TestMockQCE2E:
    """Test M4-D-5: End-to-end mock QC artifacts workflow."""

    def test_fake_backend_creates_outputs(self, tmp_path):
        """Test that FakeBackend creates fake output files."""
        backend = FakeBackend(tmp_path)

        # Mock Gaussian execution
        ts_file = tmp_path / "ts_guess.xyz"
        ts_file.write_text("6\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 2.0 0.0 0.0\n")

        result = backend.run_gaussian(ts_file, method="freq", qc_key="NBO")

        assert result["success"] is True
        assert "output_log" in result
        assert Path(result["output_log"]).exists()

        # Verify NBO output exists
        nbo_dir = tmp_path / "nbo_analysis"
        assert nbo_dir.exists()
        assert (nbo_dir / "test_job.37").exists()

        print("✓ FakeBackend creates fake QC output files (NBO only)")

    def test_packager_collects_fake_qc_artifacts(self, tmp_path):
        """Test that pack_mechanism_assets() collects fake QC artifacts."""
        # Create fake directory structure
        s4_dir = tmp_path / "S4_Data"
        s4_dir.mkdir(parents=True, exist_ok=True)

        s3_dir = tmp_path / "S3_TS"
        s3_dir.mkdir(parents=True, exist_ok=True)

        s2_dir = tmp_path / "S2_Retro"
        s2_dir.mkdir(parents=True, exist_ok=True)

        # Create fake NBO output
        nbo_dir = s3_dir / "reactant_opt/standard"
        nbo_dir.mkdir(parents=True, exist_ok=True)
        (nbo_dir / "job_reactant.37").write_text("NBO data")

        # Pack mechanism assets
        step_dirs = {
            "S1": s2_dir,
            "S2": s2_dir,
            "S3": s3_dir
        }

        config = {
            "enabled": True,
            "copy_mode": "copy",
            "use_central_resolver": False
        }

        # Create fake input files
        (s2_dir / "product_min.xyz").write_text("6\nC 0.0 0.0\nH 1.0 0.0 0.0\n")

        # Generate fake TS output
        ts_file = s3_dir / "ts_final.xyz"
        ts_file.write_text("6\nC 0.0 0.0\nH 1.0 0.0 0.0\n")

        # Generate fake reactant output
        reactant_file = s3_dir / "reactant_sp.xyz"
        reactant_file.write_text("6\nC 0.0 0.0\nH 1.0 0.0 0.0\n")

        # Forming bonds
        forming_bonds = [(0, 1), (1, 2)]

        mech_index = pack_mechanism_assets(
            step_dirs=step_dirs,
            out_dir=s4_dir,
            config=config,
            forming_bonds=forming_bonds
        )

        # Verify mech_index was created
        assert mech_index.exists()

        # Load and verify contents
        with open(mech_index, 'r') as f:
            data = json.load(f)

        assert "qc_artifacts" in data
        assert "nbo_outputs" in data["qc_artifacts"]
        assert data["qc_artifacts"]["nbo_outputs"]["filename"] == "qc_nbo.37"

        # Verify NMR/Hirshfeld are NOT present
        assert "nmr_outputs" not in data["qc_artifacts"]
        assert "hirshfeld_outputs" not in data["qc_artifacts"]

        print("✓ Packager collects NBO artifacts (NMR/Hirshfeld removed in V5.4)")


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
