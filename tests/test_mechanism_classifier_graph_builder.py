from rph_core.steps.mechanism_classifier.clean_adapter import CleanAdapter, CleanRecord
from rph_core.steps.mechanism_classifier.dr_completion import build_dr_branch_plan
from rph_core.steps.mechanism_classifier.graph_builder import GraphBuilder


def test_graph_builder_reuses_shared_semantic_forming_bonds_for_edges_and_annotations() -> None:
    builder = GraphBuilder()
    record = CleanRecord(
        reaction_id="rx_semantic_authority",
        reaction_type="4+3",
        precursor_smiles="C=CC=C",
        product_smiles="C1CCC1",
        topology="INTER",
        cyclo_mode="[4+3]",
        core_atom_map={1: 1, 2: 2, 3: 3, 4: 4},
        core_bond_changes={
            "forming": [(1, 4)],
            "breaking": [],
            "order_changed": [(1, 2)],
        },
        raw={
            "mapped_precursor_smiles": "[CH2:1]=[CH:2][CH:3]=[CH2:4]",
            "mapped_product_smiles": "[CH2:1]1[CH2:2][CH2:3][CH2:4]1",
        },
    )

    graph = builder.build(record)

    edge_forming_bonds = {
        tuple(pair)
        for edge in graph.edges
        for pair in edge.forming_bonds
    }
    annotated_forming_bonds = {
        tuple(annotation.map_space)
        for annotation in (graph.forming_bonds_annotated or [])
    }

    assert edge_forming_bonds == {(1, 4), (1, 2)}
    assert annotated_forming_bonds == edge_forming_bonds

    annotations_by_pair = {
        tuple(annotation.map_space): annotation
        for annotation in (graph.forming_bonds_annotated or [])
    }
    assert annotations_by_pair[(1, 4)].bond_type_precursor == "NONE"
    assert annotations_by_pair[(1, 4)].bond_type_product == "SINGLE"
    assert annotations_by_pair[(1, 2)].bond_type_precursor == "DOUBLE"
    assert annotations_by_pair[(1, 2)].bond_type_product == "SINGLE"


def test_clean_adapter_keeps_index_fallback_fields_out_of_semantic_bond_authority() -> None:
    adapter = CleanAdapter()

    record = adapter.parse_row(
        {
            "rxn_key_hash": "rx_index_fallback_only",
            "reaction_type": "4+3",
            "precursor_smiles": "C=C",
            "product_smiles_main": "CC",
            "topology": "INTER",
            "cyclo_mode": "[4+3]",
            "core_atom_map": "{}",
            "core_bond_changes": "",
            "formed_bond_index_pairs": "0-1",
            "forming_bonds": "0-1",
        }
    )

    assert record is not None
    assert record.core_bond_changes == {
        "forming": [],
        "breaking": [],
        "order_changed": [],
    }
    raw = record.raw or {}
    assert raw["formed_bond_index_pairs"] == "0-1"
    assert raw["forming_bonds"] == "0-1"


def test_dr_branch_plan_infers_unreported_cyclization_branch_from_forming_bonds() -> None:
    builder = GraphBuilder()
    record = CleanRecord(
        reaction_id="rx_dr_unreported",
        reaction_type="4+3",
        precursor_smiles="C=CC=C",
        product_smiles="C1CCC1",
        topology="INTRA_TYPE_I",
        cyclo_mode="[4+3]",
        core_atom_map={1: 1, 2: 2, 3: 3, 4: 4},
        core_bond_changes={
            "forming": [(1, 4)],
            "breaking": [],
            "order_changed": [(1, 2)],
        },
        raw={
            "mapped_precursor_smiles": "[CH2:1]=[CH:2][CH:3]=[CH2:4]",
            "mapped_product_smiles": "[CH2:1]1[CH2:2][CH2:3][CH2:4]1",
        },
    )

    graph = builder.build(record)
    plan = build_dr_branch_plan(graph)

    assert plan.status == "complete"
    assert [branch.branch_id for branch in plan.branches] == ["BR_MAJOR", "BR_DR_001"]
    assert plan.reported_dr == {}
    assert plan.stereocenters.inferred_from_forming_bonds is True
    assert set(plan.stereocenters.new_map_numbers) == {1, 2, 4}
    dr_branch = plan.branches[1]
    assert dr_branch.generation_policy == "concerted_bridgehead_pair_flip"
    assert len(dr_branch.flipped_map_numbers) == 3  # all new chiral centers (no cap)


def test_dr_branch_plan_preserves_reported_dr_as_metadata_not_branch_gate() -> None:
    builder = GraphBuilder()
    record = CleanRecord(
        reaction_id="rx_dr_reported",
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
            "dr_major": "4",
            "dr_minor": "1",
        },
    )

    plan = build_dr_branch_plan(builder.build(record))

    assert plan.reported_dr == {"dr_major": "4", "dr_minor": "1"}
    assert [branch.branch_id for branch in plan.branches] == ["BR_MAJOR", "BR_DR_001"]
