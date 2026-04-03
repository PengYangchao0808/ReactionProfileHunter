"""
Tests for M4-D-4: QC Artifacts Collection and mech_index.json Integration
================================================================

Tests verify:
- QC artifacts are collected from whitelisted directories
- Artifacts are copied to S4 root with fixed naming
- mech_index.json qc_artifacts field is populated correctly

V5.4 Update: Removed NMR/Hirshfeld - only NBO artifacts are collected.
"""

import json
import logging
import pytest
from pathlib import Path
import tempfile

from rph_core.steps.step4_features.mech_packager import (
    pack_mechanism_assets,
    _collect_qc_artifacts,
    _write_json_atomic,
)

logger = logging.getLogger(__name__)


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def pipeline_root(tmp_path):
    """Create a minimal pipeline directory structure for QC artifact tests."""
    root = tmp_path / "pipeline_root"
    root.mkdir(parents=True, exist_ok=True)

    # S3 directories with QC output subdirectories
    s3_ts = root / "S3_TS"
    s3_ts.mkdir(parents=True, exist_ok=True)

    nbo_dir = s3_ts / "nbo_analysis"
    nbo_dir.mkdir(parents=True, exist_ok=True)
    (nbo_dir / "job_nbo.37").write_text("NBO test data")

    # S1/S2 directories (minimal)
    s1 = root / "S1_ConfGeneration"
    s1.mkdir(parents=True, exist_ok=True)
    (s1 / "product_min.xyz").write_text("3\nProduct\nC 0 0 0\nH 1 0 0\n")

    s2 = root / "S2_Retro"
    s2.mkdir(parents=True, exist_ok=True)
    (s2 / "reactant_complex.xyz").write_text("3\nReactant\nC 0 0 0\nH 1 0 0\n")

    # S4 output directory
    s4 = root / "S4_Data"
    s4.mkdir(parents=True, exist_ok=True)

    return root


# ============================================================================
# Test Class 1: QC Artifact Collection (Pure Filesystem)
# ============================================================================

class TestQCArtifactCollection:
    """Test _collect_qc_artifacts function directly."""

    def test_collect_qc_artifacts_from_nbo_dir(self, pipeline_root):
        """Should collect NBO artifacts from whitelisted directory."""
        s3_ts = pipeline_root / "S3_TS"
        s4 = pipeline_root / "S4_Data"

        # Create NBO file in nbo_analysis directory
        nbo_dir = s3_ts / "nbo_analysis"
        nbo_dir.mkdir(exist_ok=True)
        (nbo_dir / "job_nbo.37").write_text("NBO data")

        result = _collect_qc_artifacts(
            s3_dir=s3_ts,
            pipeline_root=pipeline_root,
            out_dir=s4,
            copy_mode="copy"
        )

        assert "nbo_outputs" in result
        assert result["nbo_outputs"]["filename"] == "qc_nbo.37"
        assert (s4 / "qc_nbo.37").exists()

    def test_collect_qc_artifacts_meta_has_candidates(self, pipeline_root):
        """Meta should record candidates with mtime, size, and picked info."""
        s3_ts = pipeline_root / "S3_TS"
        s4 = pipeline_root / "S4_Data"

        # Create NBO file
        nbo_dir = s3_ts / "nbo_analysis"
        nbo_dir.mkdir(exist_ok=True)
        (nbo_dir / "job.37").write_text("NBO data")

        result = _collect_qc_artifacts(
            s3_dir=s3_ts,
            pipeline_root=pipeline_root,
            out_dir=s4,
            copy_mode="copy"
        )

        # Meta should have candidates list and picked info
        meta = result["nbo_outputs"]["meta"]
        assert "candidates" in meta
        assert "picked" in meta
        assert "reason" in meta

        # rel_path is SOURCE file relative to pipeline_root
        assert "rel_path" in meta["picked"]
        assert meta["picked"]["rel_path"] == "S3_TS/nbo_analysis/job.37"

    def test_collect_qc_artifacts_empty_dir(self, tmp_path):
        """Empty S3 directory should return empty result."""
        s3_ts = tmp_path / "S3_TS"
        s3_ts.mkdir(parents=True, exist_ok=True)
        s4 = tmp_path / "S4_Data"
        s4.mkdir(parents=True, exist_ok=True)

        result = _collect_qc_artifacts(
            s3_dir=s3_ts,
            pipeline_root=tmp_path,
            out_dir=s4,
            copy_mode="copy"
        )

        assert result == {}


# ============================================================================
# Test Class 2: QC Artifacts in mech_index.json (End-to-End)
# ============================================================================

class TestQCArtifactsInMechIndex:
    """Test qc_artifacts field in mech_index.json after pack_mechanism_assets."""

    def test_packager_creates_qc_artifacts_in_mech_index(self, pipeline_root, monkeypatch):
        """End-to-end: pack_mechanism_assets should populate qc_artifacts in mech_index.json."""
        # Prevent subprocess calls
        monkeypatch.setattr("subprocess.run", lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("No subprocess")))
        monkeypatch.setattr("subprocess.Popen", lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("No subprocess")))

        s4 = pipeline_root / "S4_Data"

        # Run packager
        mech_index = pack_mechanism_assets(
            step_dirs={
                'S1': pipeline_root / "S1_ConfGeneration",
                'S2': pipeline_root / "S2_Retro",
                'S3': pipeline_root / "S3_TS"
            },
            out_dir=s4,
            config={
                'enabled': True,
                'copy_mode': 'copy',
                'write_quality_flags': True
            },
            forming_bonds=((0, 1), (2, 3))
        )

        # Load mech_index.json and verify qc_artifacts structure
        mech_index_path = s4 / "mech_index.json"
        assert mech_index_path.exists()

        with open(mech_index_path, 'r', encoding='utf-8') as f:
            index = json.load(f)

        # qc_artifacts field should exist
        assert "qc_artifacts" in index
        qc_artifacts = index["qc_artifacts"]
        assert isinstance(qc_artifacts, dict)

        # V5.4: Only nbo_outputs should exist
        expected_keys = ["nbo_outputs"]
        for key in expected_keys:
            assert key in qc_artifacts, f"Missing expected key: {key}"

        # Verify NMR/Hirshfeld are NOT present
        assert "nmr_outputs" not in qc_artifacts
        assert "hirshfeld_outputs" not in qc_artifacts

        # Structure: {"filename": "...", "meta": {...}}
        for key, artifact in qc_artifacts.items():
            filename = artifact.get("filename", "")
            assert filename, f"Missing filename in {key}"
            assert "meta" in artifact, f"Missing meta in {key}"

        logger.info(f"✓ mech_index.json has {len(qc_artifacts)} QC artifact types: {list(qc_artifacts.keys())}")

    def test_qc_artifacts_filename_is_relative(self, pipeline_root, monkeypatch):
        """qc_artifacts filename should be relative to S4 root, not absolute."""
        monkeypatch.setattr("subprocess.run", lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("No subprocess")))

        s4 = pipeline_root / "S4_Data"

        pack_mechanism_assets(
            step_dirs={
                'S1': pipeline_root / "S1_ConfGeneration",
                'S2': pipeline_root / "S2_Retro",
                'S3': pipeline_root / "S3_TS"
            },
            out_dir=s4,
            config={'enabled': True, 'copy_mode': 'copy'},
            forming_bonds=((0, 1), (2, 3))
        )

        mech_index_path = s4 / "mech_index.json"
        with open(mech_index_path, 'r', encoding='utf-8') as f:
            index = json.load(f)

        qc_artifacts = index.get("qc_artifacts", {})
        for key, artifact in qc_artifacts.items():
            filename = artifact.get("filename", "")
            # Should NOT be an absolute path
            assert not Path(filename).is_absolute(), f"{key} filename should be relative, got: {filename}"
            # Should be just the filename (not a full path)
            assert "/" not in filename and "\\" not in filename, f"{key} filename should be simple name, got: {filename}"


# ============================================================================
# Test Class 3: Edge Cases
# ============================================================================

class TestQCArtifactsEdgeCases:
    """Test edge cases for QC artifact collection."""

    def test_missing_s3_dir_returns_empty(self, tmp_path):
        """Missing S3 directory should not error, return empty result."""
        s4 = tmp_path / "S4_Data"
        s4.mkdir(parents=True, exist_ok=True)

        result = _collect_qc_artifacts(
            s3_dir=tmp_path / "nonexistent_S3_TS",
            pipeline_root=tmp_path,
            out_dir=s4,
            copy_mode="copy"
        )

        assert result == {}

    def test_qc_files_not_in_whitelist_are_ignored(self, pipeline_root):
        """Files outside whitelisted subdirectories should be ignored."""
        s3_ts = pipeline_root / "S3_TS"
        s4 = pipeline_root / "S4_Data"

        # Create a file outside whitelist
        (s3_ts / "random_output.dat").write_text("should be ignored")

        result = _collect_qc_artifacts(
            s3_dir=s3_ts,
            pipeline_root=pipeline_root,
            out_dir=s4,
            copy_mode="copy"
        )

        # Should not contain random file
        all_filenames = []
        for artifact in result.values():
            all_filenames.append(artifact.get("filename", ""))
        assert "random_output.dat" not in all_filenames


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
