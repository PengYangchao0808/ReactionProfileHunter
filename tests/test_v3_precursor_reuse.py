from __future__ import annotations

import json
import pytest
from pathlib import Path

from rph_core.scheduling.artifact_refs import (
    link_or_copy_precursor_into_branch,
    materialize_small_molecule_refs,
    read_branch_manifest,
    write_branch_manifest,
)


class TestLinkPrecursor:
    def test_links_precursor_dir(self, tmp_path: Path):
        reaction_root = tmp_path / "RXN_001"
        precursor_dir = reaction_root / "precursor" / "S1_ConfGeneration" / "precursor"
        precursor_dir.mkdir(parents=True)
        (precursor_dir / "precursor_min.xyz").write_text("test")

        branch_root = reaction_root / "branches" / "BR_DR_001"
        link_or_copy_precursor_into_branch(reaction_root, branch_root)

        link = branch_root / "S1_ConfGeneration" / "precursor"
        assert link.exists()
        assert (link / "precursor_min.xyz").exists()

    def test_repairs_broken_precursor_symlink(self, tmp_path: Path):
        reaction_root = tmp_path / "RXN_001"
        precursor_dir = reaction_root / "precursor" / "S1_ConfGeneration" / "precursor"
        precursor_dir.mkdir(parents=True)
        (precursor_dir / "precursor_min.xyz").write_text("test")
        branch_root = reaction_root / "branches" / "BR_DR_001"
        target_dir = branch_root / "S1_ConfGeneration"
        target_dir.mkdir(parents=True)
        (target_dir / "precursor").symlink_to(Path("broken-relative-target"), target_is_directory=True)

        link_or_copy_precursor_into_branch(reaction_root, branch_root)

        link = branch_root / "S1_ConfGeneration" / "precursor"
        assert link.exists()
        assert (link / "precursor_min.xyz").exists()

    def test_no_precursor_skips(self, tmp_path: Path):
        reaction_root = tmp_path / "RXN_001"
        branch_root = reaction_root / "branches" / "BR_DR_001"
        link_or_copy_precursor_into_branch(reaction_root, branch_root)
        assert not (branch_root / "S1_ConfGeneration" / "precursor").exists()


class TestBranchManifest:
    def test_write_and_read(self, tmp_path: Path):
        branch_root = tmp_path / "BR_001"
        payload = {"branch_id": "BR_001", "status": "COMPLETE"}
        path = write_branch_manifest(branch_root, payload)
        assert path.exists()
        data = read_branch_manifest(branch_root)
        assert data["branch_id"] == "BR_001"
        assert data["status"] == "COMPLETE"

    def test_read_missing(self, tmp_path: Path):
        assert read_branch_manifest(tmp_path / "nonexistent") == {}


class TestSmallMoleculeRefsSecurity:
    def test_rejects_unsafe_small_molecule_key(self, tmp_path: Path):
        branch_root = tmp_path / "RXN" / "branches" / "BR_MAJOR"
        cache_root = tmp_path / "small_molecules"
        materialize_small_molecule_refs(branch_root, cache_root, ["../escape"])

        manifest_path = branch_root / "S1_ConfGeneration" / "small_molecules" / "small_molecule_refs.json"
        data = json.loads(manifest_path.read_text())
        assert any("Unsafe small_molecule_key" in warning for warning in data["warnings"])
        assert not (branch_root / "S1_ConfGeneration" / "escape").exists()
