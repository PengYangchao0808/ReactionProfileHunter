from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple, cast


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

    def to_dict(self) -> Dict[str, object]:
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

    def signature_payload(self) -> Dict[str, object]:
        return cast(Dict[str, object], _json_ready({
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
