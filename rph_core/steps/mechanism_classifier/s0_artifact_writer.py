from __future__ import annotations

import json
import hashlib
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from rph_core.steps.mechanism_classifier.context import ProductVariant, ReactionContext
from rph_core.steps.mechanism_classifier.dr_completion import (
    build_disabled_dr_branch_plan,
    build_dr_branch_plan,
)
from rph_core.steps.mechanism_classifier.models import (
    CycloMode,
    FormingBondNotation,
    GraphEdge,
    GraphNode,
    MechanismGraph,
    NodeState,
    PathwayInfo,
    TopologyType,
)
from rph_core.steps.mechanism_classifier.s0_record import (
    S0ReactionRecord,
    _build_map_to_product_smiles,
)
from rph_core.utils.json_io import write_json_atomic, write_text_atomic
from rph_core.utils.provenance import build_provenance, sha256_of_file
from rph_core.utils.run_id import RUN_ID_FIELD

logger = logging.getLogger(__name__)


def write_s0_artifacts(
    work_dir: Path,
    record: S0ReactionRecord,
    run_id: Optional[str] = None,
) -> Tuple[ReactionContext, List[ProductVariant], Path]:
    """Build and write all S0 mechanism artifacts.

    Returns (reaction_context, product_variants, mechanism_json_path).
    """

    stage_dir = Path(work_dir).resolve() / "S0_Mechanism"
    stage_dir.mkdir(parents=True, exist_ok=True)

    graph = _build_mechanism_graph(record)
    try:
        plan = build_dr_branch_plan(graph)
    except Exception as exc:
        logger.warning("Failed to build DR branch plan for rx_id=%s; using disabled fallback: %s", record.rx_id, exc)
        plan = build_disabled_dr_branch_plan(graph)

    product_variants = [_product_variant_from_branch(branch) for branch in plan.branches]
    reaction_context = ReactionContext(
        reaction_id=record.rx_id,
        source_csv=str(record.source_csv),
        source_row_hash=record.source_row_hash,
        mapped_reaction_smiles=_mapped_reaction_smiles(record),
        canonical_product_smiles=record.product_smiles,
        canonical_precursor_smiles=record.canonical_precursor_smiles,
        mapped_product_smiles=record.mapped_product_smiles,
        mapped_precursor_smiles=record.mapped_precursor_smiles,
        reaction_type=record.reaction_type,
        topology=record.topology,
        cyclo_mode=record.cyclo_mode,
        forming_bonds_map_space=record.mapped_forming_bonds,
        forming_bonds_product_smiles_idx=record.forming_bonds,
        mechanism_graph_ref="mechanism_graph.json",
        branch_plan_ref="dr_branch_plan.json",
        mapping_confidence=record.mapping_confidence,
        mapping_trusted=record.mapping_trusted,
    )

    mapping_dir = stage_dir / "atom_mappings"
    mapping_dir.mkdir(parents=True, exist_ok=True)
    registry_variants = []
    variant_mapping_payloads: Dict[str, Dict[str, object]] = {}
    for variant in product_variants:
        mapping_ref = Path("atom_mappings") / f"{variant.variant_id}.json"
        mapping_payload = _atom_map_smiles_payload(
            record,
            graph,
            reference_smiles=variant.branch_product_smiles,
            variant_id=variant.variant_id,
        )
        variant_mapping_payloads[str(variant.variant_id)] = mapping_payload
        _write_json(stage_dir / mapping_ref, mapping_payload)
        registry_variants.append(
            {
                **variant.to_dict(),
                "atom_mapping_ref": mapping_ref.as_posix(),
                "atom_mapping_digest": mapping_payload["mapping_digest"],
            }
        )

    major_mapping = next(
        (item for item in registry_variants if item["variant_id"] == "product_major"),
        registry_variants[0],
    )
    major_mapping_payload = variant_mapping_payloads[str(major_mapping["variant_id"])]
    mechanism_path = stage_dir / "mechanism.json"
    mechanism_payload = {
        "schema_version": "s0_mechanism_v3",
        "stage": "S0",
        RUN_ID_FIELD: run_id,
        "source": "trusted_reaction_record",
        "rx_id": record.rx_id,
        "source_csv": str(record.source_csv),
        "source_row_hash": record.source_row_hash,
        "reaction_type": record.reaction_type,
        "topology": record.topology,
        "cyclo_mode": record.cyclo_mode,
        "precursor_type": record.precursor_type,
        "canonical_product_smiles": record.product_smiles,
        "canonical_precursor_smiles": record.canonical_precursor_smiles,
        "mapped_product_smiles": record.mapped_product_smiles,
        "mapped_precursor_smiles": record.mapped_precursor_smiles,
        "forming_bonds": [list(item) for item in record.forming_bonds],
        "mapped_forming_bonds": [list(item) for item in record.mapped_forming_bonds],
        "index_base": 0,
        "index_space": "geometry_product_smiles_idx",
        "atom_mapping_ref": major_mapping["atom_mapping_ref"],
        "atom_mapping_digest": major_mapping["atom_mapping_digest"],
        "mapping_confidence": record.mapping_confidence,
        "mapping_trusted": record.mapping_trusted,
        "confidence": "trusted_dataset_mapping",
        "variants": [variant.variant_id for variant in product_variants],
        "mechanism_graph_ref": "mechanism_graph.json",
        "branch_plan_ref": "dr_branch_plan.json",
    }
    mechanism_payload["provenance"] = build_provenance(
        run_id=run_id,
        schema_version=str(mechanism_payload["schema_version"]),
        protocol_version=None,
        parent_manifest_paths={},
        atom_mapping_payload=major_mapping_payload,
        forming_bonds=record.forming_bonds,
        variant_manifest_path=None,
        extra={
            "source_csv_hash": sha256_of_file(record.source_csv),
            "source_row_hash": record.source_row_hash,
        },
    ).to_dict()
    _write_json(mechanism_path, mechanism_payload)
    _write_json(
        stage_dir / "reaction_context.json",
        {
            "schema_version": "s0_reaction_context_v1",
            **reaction_context.to_dict(),
        },
    )
    _write_text(stage_dir / "mechanism_graph.json", graph.model_dump_json(indent=2))
    _write_json(
        stage_dir / "variant_registry.json",
        {
            "schema_version": "s0_variant_registry_v2",
            "reaction_id": record.rx_id,
            "variants": registry_variants,
        },
    )
    _write_text(stage_dir / "dr_branch_plan.json", plan.model_dump_json(indent=2))
    _write_json(
        stage_dir / "atom_map_smiles.json",
        major_mapping_payload,
    )

    return reaction_context, product_variants, mechanism_path


def _build_mechanism_graph(record: S0ReactionRecord) -> MechanismGraph:
    annotations = [
        FormingBondNotation(
            map_space=map_pair,
            product_smiles_idx=product_pair,
        )
        for map_pair, product_pair in zip(record.mapped_forming_bonds, record.forming_bonds)
    ]
    return MechanismGraph(
        reaction_id=record.rx_id,
        reaction_type=record.reaction_type,
        cyclo_mode=_parse_cyclo_mode(record.cyclo_mode),
        topology=_parse_topology(record.topology),
        precursor_type=record.precursor_type or None,
        nodes=[
            GraphNode(
                node_id="N_reactants",
                smiles=record.mapped_precursor_smiles or None,
                state_type=NodeState.REACTANT,
                role="reactants",
            ),
            GraphNode(
                node_id="N_product",
                smiles=record.mapped_product_smiles,
                state_type=NodeState.PRODUCT,
                role="product",
            ),
        ],
        edges=[
            GraphEdge(
                source="N_reactants",
                target="N_product",
                forming_bonds=list(record.mapped_forming_bonds),
                pathway_id="primary",
            )
        ],
        pathways=[
            PathwayInfo(
                pathway_id="primary",
                description=f"Primary {record.reaction_type} pathway",
            )
        ],
        source_data=_source_data(record),
        forming_bonds_annotated=annotations,
    )


def _product_variant_from_branch(branch) -> ProductVariant:
    variant_id = _variant_id_for_branch(branch.branch_id)
    return ProductVariant(
        variant_id=variant_id,
        directory_name=variant_id,
        branch_id=branch.branch_id,
        pathway_id=branch.pathway_id,
        role=branch.role,
        branch_product_smiles=branch.product_smiles,
        flipped_map_numbers=tuple(branch.flipped_map_numbers),
        fixed_stereocenters=tuple(branch.fixed_stereocenters),
        selected_structure_ref=None,
    )


def _variant_id_for_branch(branch_id: str) -> str:
    if branch_id == "BR_MAJOR":
        return "product_major"
    if branch_id == "BR_DR_001":
        return "product_minor_001"
    try:
        return f"product_minor_{int(branch_id.split('_')[-1]):03d}"
    except (TypeError, ValueError):
        return f"product_minor_{branch_id.lower()}"


def _mapped_reaction_smiles(record: S0ReactionRecord) -> str:
    mapped = str(record.raw_row.get("rxn_smiles_mapped", "")).strip()
    if mapped:
        return mapped
    return f"{record.mapped_precursor_smiles}>>{record.mapped_product_smiles}"


def _source_data(record: S0ReactionRecord) -> Dict[str, str]:
    source_data = dict(record.raw_row)
    source_data["rxn_smiles_mapped"] = source_data.get("rxn_smiles_mapped", "") or _mapped_reaction_smiles(record)
    if "mapped_precursor_smiles" not in source_data and record.mapped_precursor_smiles:
        source_data["mapped_precursor_smiles"] = record.mapped_precursor_smiles
    if "mapped_product_smiles" not in source_data and record.mapped_product_smiles:
        source_data["mapped_product_smiles"] = record.mapped_product_smiles
    return source_data


def _parse_cyclo_mode(value: str) -> CycloMode:
    text = str(value or "").strip()
    for mode in CycloMode:
        if mode.value == text:
            return mode
    return CycloMode.UNKNOWN


def _parse_topology(value: str) -> TopologyType:
    text = str(value or "").strip()
    for topology in TopologyType:
        if topology.value == text:
            return topology
    if text.startswith("INTRA"):
        return TopologyType.INTRA_UNKNOWN
    return TopologyType.INTER


def _atom_map_smiles_payload(
    record: S0ReactionRecord,
    graph: MechanismGraph,
    *,
    reference_smiles: str | None = None,
    variant_id: str = "product_major",
) -> Dict[str, object]:
    geometry_smiles = reference_smiles or record.mapped_product_smiles
    map_to_product_smiles = _build_map_to_product_smiles(
        record.product_smiles,
        geometry_smiles,
    )
    product_smiles_to_map = {value: key for key, value in map_to_product_smiles.items()}
    forming_bonds = [
        [map_to_product_smiles[left], map_to_product_smiles[right]]
        for left, right in record.mapped_forming_bonds
    ]
    digest_payload = {
        "reference_smiles": geometry_smiles,
        "map_to_product_smiles": map_to_product_smiles,
        "forming_bonds_map_space": record.mapped_forming_bonds,
        "forming_bonds_product_smiles": forming_bonds,
    }
    mapping_digest = hashlib.sha256(
        json.dumps(digest_payload, sort_keys=True).encode("utf-8")
    ).hexdigest()
    return {
        "schema_version": "s0_atom_map_smiles_v2",
        "reaction_id": record.rx_id,
        "variant_id": variant_id,
        "reference_smiles": geometry_smiles,
        "mapped_product_smiles": record.mapped_product_smiles,
        "mapped_precursor_smiles": record.mapped_precursor_smiles,
        "canonical_product_smiles": record.product_smiles,
        "forming_bonds_map_space": [list(item) for item in record.mapped_forming_bonds],
        "forming_bonds": forming_bonds,
        "forming_bonds_index_space": "geometry_product_smiles_idx",
        "forming_bonds_annotated": [
            annotation.model_dump(mode="json") for annotation in (graph.forming_bonds_annotated or [])
        ],
        "map_to_product_smiles": {str(key): value for key, value in map_to_product_smiles.items()},
        "product_smiles_to_map": {str(key): value for key, value in product_smiles_to_map.items()},
        "smiles_atom_mapping": {
            "map_to_product_smiles": {str(key): value for key, value in map_to_product_smiles.items()},
            "product_smiles_to_map": {str(key): value for key, value in product_smiles_to_map.items()},
        },
        "product_smiles_idx_space": "geometry_product_smiles_idx",
        "mapping_digest": mapping_digest,
        "mapping_status": "verified",
        "index_spaces": {
            "map_numbers": "1-based atom map numbers from mapped SMILES",
            "geometry_product_smiles_idx": "0-based RDKit atom index in the exact mapped branch SMILES passed to S1",
        },
    }


def _write_json(path: Path, payload: Dict[str, object]) -> None:
    write_json_atomic(path, payload)


def _write_text(path: Path, content: str) -> None:
    write_text_atomic(path, content, encoding="utf-8")
