import json
from pathlib import Path

from rph_core.utils.task_builder import TaskSpec, BranchTaskSpec
from rph_core.utils.path_manager import get_branch_root, get_branches_root


def test_branch_task_spec_creates_valid_instance() -> None:
    parent = TaskSpec(
        rx_id="rx_001",
        row_id="rx_001",
        reaction_id="RXN_abcdef01",
        product_smiles="C1CCCCC1",
        meta={"precursor_smiles": "C=C"},
    )
    branch = BranchTaskSpec(
        parent=parent,
        branch_id="BR_DR_001",
        pathway_id="dr_concerted_face_flip_001",
        product_smiles="[C@H]1CCC[C@@H]1",
        generation_policy="concerted_bridgehead_pair_flip",
        flipped_map_numbers=[1, 2],
        fixed_stereocenters=[],
    )
    assert branch.branch_id == "BR_DR_001"
    assert branch.parent.product_smiles == "C1CCCCC1"
    assert branch.generation_policy == "concerted_bridgehead_pair_flip"


def test_get_branch_root_yields_correct_path() -> None:
    reaction_root = Path("/tmp/rph/RXN_abcdef01")
    root = get_branch_root(reaction_root, "BR_DR_001")
    assert root == reaction_root / "branches" / "BR_DR_001"


def test_get_branches_root_yields_correct_path() -> None:
    reaction_root = Path("/tmp/rph/RXN_abcdef01")
    root = get_branches_root(reaction_root)
    assert root == reaction_root / "branches"


def test_branch_task_serializable() -> None:
    parent = TaskSpec(
        rx_id="rx_001",
        row_id="rx_001",
        reaction_id="RXN_abcdef01",
        product_smiles="C1CCCCC1",
    )
    branch = BranchTaskSpec(
        parent=parent,
        branch_id="BR_DR_001",
        pathway_id="dr_concerted_face_flip_001",
        product_smiles="[C@H]1CCC[C@@H]1",
        generation_policy="concerted_bridgehead_pair_flip",
        flipped_map_numbers=[1, 2],
        fixed_stereocenters=[],
    )
    payload = {
        "branch_id": branch.branch_id,
        "pathway_id": branch.pathway_id,
        "product_smiles": branch.product_smiles,
        "generation_policy": branch.generation_policy,
        "flipped_map_numbers": branch.flipped_map_numbers,
        "fixed_stereocenters": branch.fixed_stereocenters,
    }
    reloaded = json.loads(json.dumps(payload))
    assert reloaded["branch_id"] == "BR_DR_001"
    assert reloaded["flipped_map_numbers"] == [1, 2]
