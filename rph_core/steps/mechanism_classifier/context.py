from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple, cast


COMPONENT_SCHEMA_VERSION = "rph_system_components_v1"
COMPONENT_ROLES = frozenset({"reactive", "additive", "spectator"})


def _json_ready(value: Any) -> Any:
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


def _condition_signature_payload(
    solvent: Optional[str],
    temperature_celsius: Optional[float],
    has_lewis_acid: bool,
    charge: int,
    multiplicity: int,
) -> Dict[str, object]:
    return {
        "solvent": solvent,
        "temperature_celsius": temperature_celsius,
        "has_lewis_acid": has_lewis_acid,
        "charge": charge,
        "multiplicity": multiplicity,
    }


@dataclass(frozen=True)
class SystemComponent:
    """One identity-preserving molecular component in an S0 geometry system."""

    component_id: str
    role: str
    canonical_smiles: str
    mapped_smiles: str
    charge: int
    multiplicity: int
    xyz_atom_indices: Tuple[int, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if not self.component_id.strip():
            raise ValueError("System component_id must not be empty")
        if self.role not in COMPONENT_ROLES:
            raise ValueError(
                f"Unsupported system component role {self.role!r}; "
                f"expected one of {sorted(COMPONENT_ROLES)}"
            )
        if not self.canonical_smiles.strip():
            raise ValueError("System component canonical_smiles must not be empty")
        if self.multiplicity < 1:
            raise ValueError("System component multiplicity must be at least one")
        indices = tuple(int(index) for index in self.xyz_atom_indices)
        if any(index < 0 for index in indices):
            raise ValueError("System component XYZ atom indices must be 0-based and non-negative")
        if len(indices) != len(set(indices)):
            raise ValueError("System component XYZ atom indices must be unique")
        object.__setattr__(self, "xyz_atom_indices", indices)

    def to_dict(self) -> Dict[str, object]:
        return {
            "component_id": self.component_id,
            "role": self.role,
            "canonical_smiles": self.canonical_smiles,
            "mapped_smiles": self.mapped_smiles,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "xyz_atom_indices": list(self.xyz_atom_indices),
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "SystemComponent":
        return cls(
            component_id=str(payload["component_id"]),
            role=str(payload["role"]),
            canonical_smiles=str(payload["canonical_smiles"]),
            mapped_smiles=str(payload.get("mapped_smiles") or ""),
            charge=int(payload.get("charge", 0)),
            multiplicity=int(payload.get("multiplicity", 1)),
            xyz_atom_indices=tuple(
                int(index) for index in payload.get("xyz_atom_indices", ())
            ),
        )

    def signature_payload(self) -> Dict[str, object]:
        """Return an order-stable payload suitable for checkpoint hashing."""

        return self.to_dict()


def validate_system_components(
    components: Sequence[SystemComponent],
) -> Tuple[SystemComponent, ...]:
    """Validate component IDs and resolved XYZ ownership without reordering."""

    normalized = tuple(components)
    component_ids = [component.component_id for component in normalized]
    if len(component_ids) != len(set(component_ids)):
        raise ValueError("System component IDs must be unique")

    claimed_indices: Dict[int, str] = {}
    for component in normalized:
        for index in component.xyz_atom_indices:
            owner = claimed_indices.get(index)
            if owner is not None:
                raise ValueError(
                    f"XYZ atom index {index} is claimed by both {owner!r} "
                    f"and {component.component_id!r}"
                )
            claimed_indices[index] = component.component_id
    return normalized


def components_from_s0_payload(
    payload: Mapping[str, Any],
) -> Tuple[SystemComponent, ...]:
    """Load v1 components or wrap a legacy single-product S0 payload."""

    raw_components = payload.get("components")
    if isinstance(raw_components, list) and raw_components:
        return validate_system_components(
            tuple(SystemComponent.from_dict(item) for item in raw_components)
        )

    canonical_smiles = str(
        payload.get("canonical_product_smiles")
        or payload.get("product_smiles")
        or ""
    )
    if not canonical_smiles:
        raise ValueError("Legacy S0 payload has no canonical product SMILES")
    mapped_smiles = str(payload.get("mapped_product_smiles") or "")
    raw_indices = payload.get("product_xyz_atom_indices", ())
    return (
        SystemComponent(
            component_id="substrate_0",
            role="reactive",
            canonical_smiles=canonical_smiles,
            mapped_smiles=mapped_smiles,
            charge=int(payload.get("charge", 0)),
            multiplicity=int(payload.get("multiplicity", 1)),
            xyz_atom_indices=tuple(int(index) for index in raw_indices),
        ),
    )


def component_signature_payload(
    components: Sequence[SystemComponent],
) -> Dict[str, object]:
    """Preserve component order because it defines global XYZ atom identity."""

    normalized = validate_system_components(components)
    return {
        "component_schema_version": COMPONENT_SCHEMA_VERSION,
        "components": [component.signature_payload() for component in normalized],
        "reactive_component_ids": [
            component.component_id
            for component in normalized
            if component.role == "reactive"
        ],
        "additive_component_ids": [
            component.component_id
            for component in normalized
            if component.role == "additive"
        ],
    }


def validate_forming_bonds_for_components(
    forming_bonds: Sequence[Tuple[int, int]],
    components: Sequence[SystemComponent],
) -> None:
    """Reject bonds that reference additive or spectator component atoms."""

    normalized = validate_system_components(components)
    nonreactive_indices = {
        index
        for component in normalized
        if component.role != "reactive"
        for index in component.xyz_atom_indices
    }
    for raw_a, raw_b in forming_bonds:
        atom_a, atom_b = int(raw_a), int(raw_b)
        if atom_a in nonreactive_indices or atom_b in nonreactive_indices:
            raise ValueError(
                f"Forming bond ({atom_a}, {atom_b}) references a non-reactive component atom"
            )

@dataclass(frozen=True)
class ReactionContext:
    """S0-authoritative reaction identity. Written once by S0, read by S1-S4."""

    reaction_id: str
    source_csv: str
    source_row_hash: str
    mapped_reaction_smiles: str
    canonical_product_smiles: str
    canonical_precursor_smiles: str
    mapped_product_smiles: str
    mapped_precursor_smiles: str
    reaction_type: str
    topology: str
    cyclo_mode: str
    forming_bonds_map_space: Tuple[Tuple[int, int], ...]
    forming_bonds_product_smiles_idx: Tuple[Tuple[int, int], ...]
    mechanism_graph_ref: str
    branch_plan_ref: str
    mapping_confidence: Optional[float]
    mapping_trusted: bool
    components: Tuple[SystemComponent, ...] = field(default_factory=tuple)

    def normalized_components(self) -> Tuple[SystemComponent, ...]:
        if self.components:
            return validate_system_components(self.components)
        return components_from_s0_payload(self.to_legacy_dict())

    def to_legacy_dict(self) -> Dict[str, object]:
        """Return the V1 fields used to load component-free S0 contexts."""

        return cast(Dict[str, object], _json_ready({
            "reaction_id": self.reaction_id,
            "source_csv": self.source_csv,
            "source_row_hash": self.source_row_hash,
            "mapped_reaction_smiles": self.mapped_reaction_smiles,
            "canonical_product_smiles": self.canonical_product_smiles,
            "canonical_precursor_smiles": self.canonical_precursor_smiles,
            "mapped_product_smiles": self.mapped_product_smiles,
            "mapped_precursor_smiles": self.mapped_precursor_smiles,
            "reaction_type": self.reaction_type,
            "topology": self.topology,
            "cyclo_mode": self.cyclo_mode,
            "forming_bonds_map_space": self.forming_bonds_map_space,
            "forming_bonds_product_smiles_idx": self.forming_bonds_product_smiles_idx,
            "mechanism_graph_ref": self.mechanism_graph_ref,
            "branch_plan_ref": self.branch_plan_ref,
            "mapping_confidence": self.mapping_confidence,
            "mapping_trusted": self.mapping_trusted,
        }))

    def to_dict(self) -> Dict[str, object]:
        return {
            **self.to_legacy_dict(),
            **component_signature_payload(self.normalized_components()),
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ReactionContext":
        """Load current or legacy component-free reaction context payloads."""

        return cls(
            reaction_id=str(payload["reaction_id"]),
            source_csv=str(payload.get("source_csv") or ""),
            source_row_hash=str(payload.get("source_row_hash") or ""),
            mapped_reaction_smiles=str(payload.get("mapped_reaction_smiles") or ""),
            canonical_product_smiles=str(payload["canonical_product_smiles"]),
            canonical_precursor_smiles=str(
                payload.get("canonical_precursor_smiles") or ""
            ),
            mapped_product_smiles=str(payload.get("mapped_product_smiles") or ""),
            mapped_precursor_smiles=str(payload.get("mapped_precursor_smiles") or ""),
            reaction_type=str(payload.get("reaction_type") or ""),
            topology=str(payload.get("topology") or "INTER"),
            cyclo_mode=str(payload.get("cyclo_mode") or "UNKNOWN"),
            forming_bonds_map_space=tuple(
                tuple(int(index) for index in pair)
                for pair in payload.get("forming_bonds_map_space", ())
            ),
            forming_bonds_product_smiles_idx=tuple(
                tuple(int(index) for index in pair)
                for pair in payload.get("forming_bonds_product_smiles_idx", ())
            ),
            mechanism_graph_ref=str(payload.get("mechanism_graph_ref") or ""),
            branch_plan_ref=str(payload.get("branch_plan_ref") or ""),
            mapping_confidence=(
                float(payload["mapping_confidence"])
                if payload.get("mapping_confidence") is not None
                else None
            ),
            mapping_trusted=bool(payload.get("mapping_trusted", False)),
            components=components_from_s0_payload(payload),
        )

    def signature_payload(self) -> Dict[str, object]:
        payload = cast(Dict[str, object], _json_ready({
            "reaction_id": self.reaction_id,
            "source_row_hash": self.source_row_hash,
            "mapped_reaction_smiles": self.mapped_reaction_smiles,
            "canonical_product_smiles": self.canonical_product_smiles,
            "canonical_precursor_smiles": self.canonical_precursor_smiles,
            "mapped_product_smiles": self.mapped_product_smiles,
            "mapped_precursor_smiles": self.mapped_precursor_smiles,
            "reaction_type": self.reaction_type,
            "topology": self.topology,
            "cyclo_mode": self.cyclo_mode,
            "forming_bonds_map_space": self.forming_bonds_map_space,
            "forming_bonds_product_smiles_idx": self.forming_bonds_product_smiles_idx,
            "mapping_confidence": self.mapping_confidence,
            "mapping_trusted": self.mapping_trusted,
            "mechanism_graph_ref": self.mechanism_graph_ref,
            "branch_plan_ref": self.branch_plan_ref,
        }))
        payload.update(component_signature_payload(self.normalized_components()))
        return payload


@dataclass(frozen=True)
class RunContext:
    """One-per-run condition context. Written to run.manifest.json in P2; for now
    P1 only records the fields from CSV that define run conditions."""

    condition_signature: str
    solvent: Optional[str]
    temperature_celsius: Optional[float]
    has_lewis_acid: bool
    charge: int
    multiplicity: int

    def to_dict(self) -> Dict[str, object]:
        payload = _condition_signature_payload(
            self.solvent,
            self.temperature_celsius,
            self.has_lewis_acid,
            self.charge,
            self.multiplicity,
        )
        condition_signature = self.condition_signature or hashlib.sha256(
            json.dumps(payload, sort_keys=True).encode("utf-8")
        ).hexdigest()
        return {
            "condition_signature": condition_signature,
            **payload,
        }


@dataclass(frozen=True)
class ProductVariant:
    """One product-side structural variant. Created by S0 only."""

    variant_id: str
    directory_name: str
    branch_id: str
    pathway_id: str
    role: str
    branch_product_smiles: Optional[str]
    flipped_map_numbers: Tuple[int, ...] = field(default_factory=tuple)
    fixed_stereocenters: Tuple[int, ...] = field(default_factory=tuple)
    selected_structure_ref: Optional[str] = None

    def to_dict(self) -> Dict[str, object]:
        return cast(Dict[str, object], _json_ready({
            "variant_id": self.variant_id,
            "directory_name": self.directory_name,
            "branch_id": self.branch_id,
            "pathway_id": self.pathway_id,
            "role": self.role,
            "branch_product_smiles": self.branch_product_smiles,
            "flipped_map_numbers": self.flipped_map_numbers,
            "fixed_stereocenters": self.fixed_stereocenters,
            "selected_structure_ref": self.selected_structure_ref,
        }))

    def signature_payload(self) -> Dict[str, object]:
        return cast(Dict[str, object], _json_ready({
            "variant_id": self.variant_id,
            "directory_name": self.directory_name,
            "branch_id": self.branch_id,
            "pathway_id": self.pathway_id,
            "role": self.role,
            "branch_product_smiles": self.branch_product_smiles,
            "flipped_map_numbers": self.flipped_map_numbers,
            "fixed_stereocenters": self.fixed_stereocenters,
        }))
