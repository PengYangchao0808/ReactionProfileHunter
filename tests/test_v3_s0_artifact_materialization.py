"""Tests for S0 artifact materialization into branch workspaces.

Validates that link_or_copy_s0_into_branch correctly propagates
reaction-level S0 planning artifacts into each branch's S0_Mechanism directory.
"""
import json
import pytest
from pathlib import Path

from rph_core.scheduling.artifact_refs import link_or_copy_s0_into_branch


@pytest.fixture
def reaction_root(tmp_path: Path) -> Path:
    root = tmp_path / "RXN_test"
    root.mkdir()
    return root


@pytest.fixture
def s0_artifacts(reaction_root: Path) -> Path:
    s0_dir = reaction_root / "S0_Mechanism"
    s0_dir.mkdir()
    (s0_dir / "mechanism_summary.json").write_text(
        json.dumps({"forming_bonds": [[1, 2], [3, 4]], "reaction_type": "[4+3]"})
    )
    (s0_dir / "mechanism_graph.json").write_text(
        json.dumps({"nodes": [], "edges": []})
    )
    (s0_dir / "atom_map_smiles.json").write_text(
        json.dumps({"mapping": {}})
    )
    (s0_dir / "dr_branch_plan.json").write_text(
        json.dumps({"branches": ["BR_MAJOR", "BR_DR_001"]})
    )
    return s0_dir


def _branch_root(reaction_root: Path, branch_id: str = "BR_MAJOR") -> Path:
    br = reaction_root / "branches" / branch_id
    br.mkdir(parents=True, exist_ok=True)
    return br


class TestLinkOrCopyS0IntoBranch:
    def test_copies_all_artifacts(self, reaction_root, s0_artifacts):
        branch_root = _branch_root(reaction_root)
        link_or_copy_s0_into_branch(reaction_root, branch_root)

        dst = branch_root / "S0_Mechanism"
        assert (dst / "mechanism_summary.json").exists()
        assert (dst / "mechanism_graph.json").exists()
        assert (dst / "atom_map_smiles.json").exists()
        assert (dst / "dr_branch_plan.json").exists()

    def test_content_preserved(self, reaction_root, s0_artifacts):
        branch_root = _branch_root(reaction_root)
        link_or_copy_s0_into_branch(reaction_root, branch_root)

        dst = branch_root / "S0_Mechanism"
        data = json.loads((dst / "mechanism_summary.json").read_text())
        assert data["forming_bonds"] == [[1, 2], [3, 4]]
        assert data["reaction_type"] == "[4+3]"

    def test_no_source_dir_is_noop(self, reaction_root):
        branch_root = _branch_root(reaction_root)
        assert not (reaction_root / "S0_Mechanism").exists()
        link_or_copy_s0_into_branch(reaction_root, branch_root)
        assert not (branch_root / "S0_Mechanism").exists()

    def test_idempotent(self, reaction_root, s0_artifacts):
        branch_root = _branch_root(reaction_root)
        link_or_copy_s0_into_branch(reaction_root, branch_root)
        dst_file = branch_root / "S0_Mechanism" / "mechanism_summary.json"
        first_mtime = dst_file.stat().st_mtime

        link_or_copy_s0_into_branch(reaction_root, branch_root)
        assert dst_file.stat().st_mtime == first_mtime

    def test_partial_artifacts_only_copies_existing(self, reaction_root):
        s0_dir = reaction_root / "S0_Mechanism"
        s0_dir.mkdir()
        (s0_dir / "mechanism_summary.json").write_text('{"forming_bonds": []}')
        branch_root = _branch_root(reaction_root)

        link_or_copy_s0_into_branch(reaction_root, branch_root)

        dst = branch_root / "S0_Mechanism"
        assert (dst / "mechanism_summary.json").exists()
        assert not (dst / "mechanism_graph.json").exists()
        assert not (dst / "dr_branch_plan.json").exists()

    def test_multiple_branches_independent(self, reaction_root, s0_artifacts):
        for bid in ("BR_MAJOR", "BR_DR_001"):
            branch_root = _branch_root(reaction_root, bid)
            link_or_copy_s0_into_branch(reaction_root, branch_root)

        for bid in ("BR_MAJOR", "BR_DR_001"):
            dst = reaction_root / "branches" / bid / "S0_Mechanism"
            assert (dst / "mechanism_summary.json").exists()

    def test_does_not_overwrite_existing(self, reaction_root, s0_artifacts):
        branch_root = _branch_root(reaction_root)
        dst = branch_root / "S0_Mechanism"
        dst.mkdir(parents=True)
        existing_data = {"forming_bonds": [[99, 98]], "existing": True}
        (dst / "mechanism_summary.json").write_text(json.dumps(existing_data))

        link_or_copy_s0_into_branch(reaction_root, branch_root)

        kept = json.loads((dst / "mechanism_summary.json").read_text())
        assert kept == existing_data
