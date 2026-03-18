"""
Tests for QC artifact collection structure and file operations (NBO-only, V5.4).

Tests verify:
- _collect_qc_artifacts returns correct structure with filename + meta
- Files are correctly copied/linked to S4 root
- Relative paths are calculated correctly

V5.4 Update: Removed NMR/Hirshfeld - only NBO artifacts are tested.
"""

import json
import tempfile
from pathlib import Path

import pytest

from rph_core.steps.step4_features.mech_packager import (
    _collect_qc_artifacts,
    QC_ARTIFACT_TARGETS,
    QC_ARTIFACT_PATTERNS,
)


class TestQCArtifactsCollection:
    """Tests for restricted QC artifact collection (NBO-only)."""

    def test_collect_returns_structure_with_filename_and_meta(self):
        """_collect_qc_artifacts should return relative filename + meta.source_paths."""
        with tempfile.TemporaryDirectory() as tmpdir:
            pipeline_root = Path(tmpdir)
            s4_dir = pipeline_root / "S4_Data"
            s4_dir.mkdir()

            # Create S3 directory and test NBO file in whitelisted subdir
            s3_dir = pipeline_root / "S3_TS"
            nbo_dir = s3_dir / "nbo_analysis"
            nbo_dir.mkdir(parents=True)
            nbo_file = nbo_dir / "job_nbo.37"
            nbo_file.write_text("NBO data")

            result = _collect_qc_artifacts(
                s3_dir=s3_dir,
                pipeline_root=pipeline_root,
                out_dir=s4_dir,
                copy_mode="copy"
            )

            # Verify structure
            assert 'nbo_outputs' in result
            assert result['nbo_outputs']['filename'] == 'qc_nbo.37'
            assert 'meta' in result['nbo_outputs']
            assert 'candidates' in result['nbo_outputs']['meta']
            assert 'picked' in result['nbo_outputs']['meta']
            assert 'reason' in result['nbo_outputs']['meta']
            picked = result['nbo_outputs']['meta']['picked']
            assert 'rel_path' in picked

    def test_collect_copies_file_to_s4_root(self):
        """Collected artifact should exist in S4 directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            pipeline_root = Path(tmpdir)
            s4_dir = pipeline_root / "S4_Data"
            s4_dir.mkdir()

            # Create S3 directory and test NBO file
            s3_dir = pipeline_root / "S3_TS"
            nbo_dir = s3_dir / "nbo_analysis"
            nbo_dir.mkdir(parents=True)
            nbo_file = nbo_dir / "job.37"
            nbo_file.write_text("NBO data")

            result = _collect_qc_artifacts(
                s3_dir=s3_dir,
                pipeline_root=pipeline_root,
                out_dir=s4_dir,
                copy_mode="copy"
            )

            # Verify file was copied
            target_file = s4_dir / "qc_nbo.37"
            assert target_file.exists()
            assert target_file.read_text() == "NBO data"

    def test_collect_handles_fallback_search(self):
        """Should find files in S3 root when subdir not found (NBO only)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            pipeline_root = Path(tmpdir)
            s4_dir = pipeline_root / "S4_Data"
            s4_dir.mkdir()

            # Create S3 directory with file in root (no subdir) - fallback behavior
            s3_dir = pipeline_root / "S3_TS"
            s3_dir.mkdir()
            nbo_file = s3_dir / "analysis.37"
            nbo_file.write_text("NBO data")

            result = _collect_qc_artifacts(
                s3_dir=s3_dir,
                pipeline_root=pipeline_root,
                out_dir=s4_dir,
                copy_mode="copy"
            )

            # Should still find the NBO file
            assert 'nbo_outputs' in result
            assert result['nbo_outputs']['filename'] == 'qc_nbo.37'
            assert (s4_dir / "qc_nbo.37").exists()

    def test_fixed_target_filenames_are_correct(self):
        """Target filenames should match expected fixed mapping (NBO-only)."""
        expected = {
            "nbo_outputs": "qc_nbo.37"
        }
        assert QC_ARTIFACT_TARGETS == expected

        # Verify NMR/Hirshfeld are NOT present
        assert "nmr_outputs" not in QC_ARTIFACT_TARGETS
        assert "hirshfeld_outputs" not in QC_ARTIFACT_TARGETS

    def test_qc_artifact_patterns_defined(self):
        """QC artifact patterns should be defined for NBO type."""
        # V5.4: Only nbo_outputs should have patterns defined
        for artifact_type in ["nbo_outputs"]:
            assert artifact_type in QC_ARTIFACT_PATTERNS
            patterns = QC_ARTIFACT_PATTERNS[artifact_type]
            assert len(patterns) > 0, f"No patterns defined for {artifact_type}"
            # Verify patterns look reasonable
            for pattern in patterns:
                assert '*' in pattern or pattern.endswith('.37') or pattern.endswith('.nbo'), \
                    f"Pattern {pattern} looks suspicious"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
