import json

from rph_core.steps.mechanism_classifier.clean_adapter import CleanRecord
from rph_core.steps.mechanism_classifier.graph_builder import GraphBuilder
from rph_core.steps.mechanism_classifier.dr_completion import (
    build_dr_branch_plan,
    build_disabled_dr_branch_plan,
)


def test_disabled_dr_plan_returns_single_branch() -> None:
    builder = GraphBuilder()
    record = CleanRecord(
        reaction_id="rx_disabled_test",
        reaction_type="4+3",
        precursor_smiles="C=CC=C",
        product_smiles="C1CCC1",
        topology="INTRA_TYPE_I",
        cyclo_mode="[4+3]",
        core_atom_map={1: 1, 2: 2, 3: 3, 4: 4},
        core_bond_changes={
            "forming": [(1, 4), (2, 3)],
            "breaking": [],
            "order_changed": [],
        },
        raw={
            "mapped_precursor_smiles": "[CH2:1]=[CH:2][CH:3]=[CH2:4]",
            "mapped_product_smiles": "[CH2:1]1[CH:2][CH:3][CH2:4]1",
        },
    )
    graph = builder.build(record)
    plan = build_disabled_dr_branch_plan(graph)

    assert plan.status == "disabled"
    assert len(plan.branches) == 1
    assert plan.branches[0].branch_id == "BR_MAJOR"
    assert "disabled" in plan.warnings[0].lower()


def test_disabled_plan_is_valid_json_roundtrip() -> None:
    builder = GraphBuilder()
    record = CleanRecord(
        reaction_id="rx_json_test",
        reaction_type="4+3",
        precursor_smiles="C=CC=C",
        product_smiles="C1CCC1",
        topology="INTER",
        cyclo_mode="[4+3]",
        core_atom_map={},
        core_bond_changes={"forming": [(1, 2)], "breaking": [], "order_changed": []},
        raw={},
    )
    graph = builder.build(record)
    plan = build_disabled_dr_branch_plan(graph)
    payload = plan.model_dump(mode="json")
    reloaded = json.loads(json.dumps(payload))
    assert reloaded["status"] == "disabled"
    assert len(reloaded["branches"]) == 1


def test_enabled_plan_infers_branches() -> None:
    builder = GraphBuilder()
    record = CleanRecord(
        reaction_id="rx_enabled_test",
        reaction_type="4+3",
        precursor_smiles="C=CC=C",
        product_smiles="C1CCC1",
        topology="INTRA_TYPE_I",
        cyclo_mode="[4+3]",
        core_atom_map={1: 1, 2: 2, 3: 3, 4: 4},
        core_bond_changes={
            "forming": [(1, 4), (2, 3)],
            "breaking": [],
            "order_changed": [],
        },
        raw={
            "mapped_precursor_smiles": "[CH2:1]=[CH:2][CH:3]=[CH2:4]",
            "mapped_product_smiles": "[CH2:1]1[CH:2][CH:3][CH2:4]1",
        },
    )
    graph = builder.build(record)
    plan = build_dr_branch_plan(graph)

    assert plan.status == "complete"
    assert len(plan.branches) >= 2
    assert plan.stereocenters.inferred_from_forming_bonds is True


def test_select_flip_targets_with_component_info():
    from rph_core.steps.mechanism_classifier.dr_completion import _select_flip_targets

    # [4+3] case: 3 new chiral centers, furan-side only should flip
    result = _select_flip_targets(
        [3, 7, 10],
        component_atoms={"furan": [7, 10], "allenamide": [1, 3]},
    )
    assert result == [7, 10]


def test_select_flip_targets_without_component_info_returns_all():
    from rph_core.steps.mechanism_classifier.dr_completion import _select_flip_targets

    # Without component info, all candidates returned (no cap)
    assert _select_flip_targets([1, 3, 7, 9, 12]) == [1, 3, 7, 9, 12]
    assert _select_flip_targets([5, 3, 1, 8]) == [1, 3, 5, 8]
    assert _select_flip_targets([7, 2]) == [2, 7]
    assert _select_flip_targets([4]) == [4]
    assert _select_flip_targets([]) == []
