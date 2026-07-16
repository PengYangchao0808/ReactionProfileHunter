from pathlib import Path

import pytest

from rph_core.steps.mechanism_classifier.context import (
    COMPONENT_SCHEMA_VERSION,
    ReactionContext,
    SystemComponent,
    component_signature_payload,
    validate_forming_bonds_for_components,
)
from rph_core.utils.v4_checkpoint import V4Checkpoint


def _component(
    component_id: str,
    role: str,
    indices: tuple[int, ...],
    smiles: str = "CC",
) -> SystemComponent:
    return SystemComponent(
        component_id=component_id,
        role=role,
        canonical_smiles=smiles,
        mapped_smiles=smiles,
        charge=0,
        multiplicity=1,
        xyz_atom_indices=indices,
    )


def _legacy_context_payload() -> dict[str, object]:
    return {
        "schema_version": "s0_reaction_context_v1",
        "reaction_id": "RX_LEGACY",
        "source_csv": "fixture.csv",
        "source_row_hash": "row-hash",
        "mapped_reaction_smiles": "[C:1]>>[C:1]",
        "canonical_product_smiles": "CC",
        "canonical_precursor_smiles": "CC",
        "mapped_product_smiles": "[CH3:1][CH3:2]",
        "mapped_precursor_smiles": "[CH3:1][CH3:2]",
        "reaction_type": "fixture",
        "topology": "INTRA",
        "cyclo_mode": "UNKNOWN",
        "forming_bonds_map_space": [[1, 2]],
        "forming_bonds_product_smiles_idx": [[0, 1]],
        "mechanism_graph_ref": "mechanism_graph.json",
        "branch_plan_ref": "dr_branch_plan.json",
        "mapping_confidence": 1.0,
        "mapping_trusted": True,
    }


def test_component_roundtrip_single():
    context = ReactionContext.from_dict(_legacy_context_payload())

    assert len(context.components) == 1
    assert context.components[0].component_id == "substrate_0"
    assert context.components[0].role == "reactive"
    assert context.components[0].canonical_smiles == "CC"
    assert context.components[0].mapped_smiles == "[CH3:1][CH3:2]"
    assert context.to_dict()["component_schema_version"] == COMPONENT_SCHEMA_VERSION
    assert context.to_dict()["reactive_component_ids"] == ["substrate_0"]


def test_component_roundtrip_identity():
    original = _component("substrate_0", "reactive", (0, 1), smiles="C=C")

    restored = SystemComponent.from_dict(original.to_dict())

    assert restored == original
    assert restored.canonical_smiles == "C=C"
    assert restored.xyz_atom_indices == (0, 1)


def test_component_two_reactive():
    payload = component_signature_payload((
        _component("substrate_0", "reactive", (0, 1)),
        _component("substrate_1", "reactive", (2, 3), smiles="N=N"),
    ))

    assert payload["reactive_component_ids"] == ["substrate_0", "substrate_1"]
    assert payload["additive_component_ids"] == []


def test_component_la_not_reactive():
    payload = component_signature_payload((
        _component("substrate_0", "reactive", (0, 1)),
        _component("additive_0", "additive", (2, 3), smiles="[Li]Cl"),
    ))

    assert payload["reactive_component_ids"] == ["substrate_0"]
    assert payload["additive_component_ids"] == ["additive_0"]


def test_component_order_changes_signature():
    substrate_0 = _component("substrate_0", "reactive", (0, 1))
    substrate_1 = _component("substrate_1", "reactive", (2, 3), smiles="N=N")

    forward = V4Checkpoint.signature(
        component_signature_payload((substrate_0, substrate_1))
    )
    reversed_order = V4Checkpoint.signature(
        component_signature_payload((substrate_1, substrate_0))
    )

    assert forward != reversed_order


def test_component_role_changes_invalidate(tmp_path: Path):
    checkpoint = V4Checkpoint(tmp_path)
    s0_manifest = tmp_path / "S0_Mechanism" / "mechanism.json"
    s2_manifest = tmp_path / "S2_PEB" / "manifest.json"
    s0_manifest.parent.mkdir(parents=True)
    s2_manifest.parent.mkdir(parents=True)
    s0_manifest.write_text("{}", encoding="utf-8")
    s2_manifest.write_text("{}", encoding="utf-8")

    reactive = (_component("substrate_0", "reactive", (0, 1)),)
    additive = (_component("substrate_0", "additive", (0, 1)),)
    reactive_payload = component_signature_payload(reactive)
    additive_payload = component_signature_payload(additive)
    reactive_signature = checkpoint.signature(reactive_payload)
    additive_signature = checkpoint.signature(additive_payload)
    checkpoint.mark_s0(
        reactive_signature,
        s0_manifest,
        signature_payload=reactive_payload,
    )
    checkpoint.mark(
        "s2",
        checkpoint.signature({"stage": "s2"}),
        s2_manifest,
    )

    assert reactive_signature != additive_signature
    checkpoint.mark_s0(
        additive_signature,
        s0_manifest,
        signature_payload=additive_payload,
    )
    assert not checkpoint.reusable(
        "s2", checkpoint.signature({"stage": "s2"}), s2_manifest
    )


def test_forming_bonds_exclude_additive():
    components = (
        _component("substrate_0", "reactive", (0, 1)),
        _component("additive_0", "additive", (2, 3), smiles="[Li]Cl"),
    )

    with pytest.raises(ValueError, match="non-reactive component atom"):
        validate_forming_bonds_for_components(((0, 2),), components)
