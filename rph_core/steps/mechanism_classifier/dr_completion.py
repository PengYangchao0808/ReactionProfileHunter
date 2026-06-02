from __future__ import annotations

import json
import logging
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from rph_core.steps.mechanism_classifier.dr_models import (
    ChiralCenterSet,
    ChiralCluster,
    DRBranchPlan,
    DRProductBranch,
)
from rph_core.steps.mechanism_classifier.models import MechanismGraph

logger = logging.getLogger(__name__)


def build_dr_branch_plan(graph: MechanismGraph) -> DRBranchPlan:
    """Build the S0-authoritative DR branch plan from a mechanism graph.

    The plan is intentionally conservative: it records the stereochemical branch
    intent in map-space even when RDKit cannot produce explicit `@/@@` product
    variants. Downstream task expansion should consume this file instead of
    re-interpreting DR from raw dataset columns.
    """

    raw = graph.source_data if isinstance(graph.source_data, dict) else {}
    forming_bonds = _forming_bonds_map_space(graph)
    precursor_smiles = _mapped_smiles(raw, "precursor")
    product_smiles = _mapped_smiles(raw, "product") or _product_node_smiles(graph)

    warnings: List[str] = []
    precursor_centers = _mapped_chiral_centers(precursor_smiles, warnings, "precursor")
    product_centers = _mapped_chiral_centers(product_smiles, warnings, "product")

    fixed_centers = sorted(center for center in product_centers if center in precursor_centers)
    new_centers = sorted(center for center in product_centers if center not in precursor_centers)
    inferred_from_forming_bonds = False

    if not new_centers:
        inferred = _forming_bond_center_candidates(forming_bonds)
        if inferred:
            new_centers = inferred
            inferred_from_forming_bonds = True
            warnings.append(
                "Product SMILES did not define new chiral centers; inferred DR centers from forming bonds."
            )

    stereocenters = ChiralCenterSet(
        fixed_map_numbers=fixed_centers,
        new_map_numbers=new_centers,
        inferred_from_forming_bonds=inferred_from_forming_bonds,
    )
    component_atoms = _extract_component_atoms(raw)
    clusters = _build_clusters(new_centers, forming_bonds)
    branches = _build_branches(
        product_smiles=product_smiles,
        fixed_centers=fixed_centers,
        new_centers=new_centers,
        clusters=clusters,
        component_atoms=component_atoms,
    )

    status = "complete" if len(branches) > 1 else "not_applicable"
    return DRBranchPlan(
        reaction_id=graph.reaction_id,
        status=status,
        reaction_type=graph.reaction_type,
        cyclo_mode=str(getattr(graph.cyclo_mode, "value", graph.cyclo_mode)),
        topology=str(getattr(graph.topology, "value", graph.topology)),
        reported_dr=_reported_dr(raw),
        stereocenters=stereocenters,
        chiral_clusters=clusters,
        branches=branches,
        warnings=warnings,
    )


def build_disabled_dr_branch_plan(graph: MechanismGraph) -> DRBranchPlan:
    """Return a minimal plan when s0.dr_completion.enabled is False.

    Still emits a valid plan so that layout contracts and checkpoint resume
    do not break when DR completion is toggled off.
    """

    return DRBranchPlan(
        reaction_id=graph.reaction_id,
        status="disabled",
        reaction_type=graph.reaction_type,
        cyclo_mode=str(getattr(graph.cyclo_mode, "value", graph.cyclo_mode)),
        topology=str(getattr(graph.topology, "value", graph.topology)),
        reported_dr={},
        stereocenters=ChiralCenterSet(),
        chiral_clusters=[],
        branches=[
            DRProductBranch(
                branch_id="BR_MAJOR",
                pathway_id="primary",
                role="major_or_reference",
                generation_policy="reference_product",
            )
        ],
        warnings=["DR completion disabled by config (s0.dr_completion.enabled=false)."],
    )


def attach_dr_completion(graph: MechanismGraph) -> DRBranchPlan:
    """Attach DR completion payload to the graph and return the branch plan."""

    plan = build_dr_branch_plan(graph)
    graph.dr_completion = plan.model_dump(mode="json")
    return plan


def _forming_bonds_map_space(graph: MechanismGraph) -> List[Tuple[int, int]]:
    pairs: List[Tuple[int, int]] = []
    if graph.forming_bonds_annotated:
        for annotation in graph.forming_bonds_annotated:
            pairs.append(_sorted_pair(annotation.map_space))
    else:
        for edge in graph.get_edges_for_pathway("primary"):
            for pair in edge.forming_bonds:
                pairs.append(_sorted_pair(pair))
    return _dedupe_pairs(pairs)


def _mapped_smiles(raw: Dict[str, Any], side: str) -> Optional[str]:
    rxn_mapped = raw.get("rxn_smiles_mapped")
    if isinstance(rxn_mapped, str) and ">>" in rxn_mapped:
        precursor, product = rxn_mapped.split(">>", 1)
        candidate = precursor if side == "precursor" else product
        candidate = candidate.strip()
        if candidate:
            return candidate

    keys = (
        ("mapped_precursor_smiles", "precursor_smiles_mapped", "mapped_reactant_smiles", "reactant_smiles_mapped")
        if side == "precursor"
        else ("mapped_product_smiles", "product_smiles_mapped", "mapped_product", "product_smiles_main")
    )
    for key in keys:
        value = raw.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def _product_node_smiles(graph: MechanismGraph) -> Optional[str]:
    products = graph.get_product_nodes()
    if not products:
        return None
    return products[0].smiles


def _mapped_chiral_centers(
    smiles: Optional[str],
    warnings: List[str],
    label: str,
) -> List[int]:
    if not smiles:
        return []
    try:
        from rdkit import Chem
    except ImportError:
        warnings.append("RDKit unavailable; stereocenters inferred from forming bonds only.")
        return []

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            warnings.append(f"Could not parse {label} mapped SMILES for stereocenter detection.")
            return []
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    except Exception as exc:
        warnings.append(f"Failed to inspect {label} stereocenters: {exc}")
        return []

    map_numbers: List[int] = []
    for atom_idx, _configuration in centers:
        atom = mol.GetAtomWithIdx(int(atom_idx))
        map_num = int(atom.GetAtomMapNum())
        if map_num > 0:
            map_numbers.append(map_num)
    return sorted(set(map_numbers))


def _forming_bond_center_candidates(forming_bonds: Sequence[Tuple[int, int]]) -> List[int]:
    centers: List[int] = []
    for left, right in forming_bonds:
        centers.extend([left, right])
    return sorted(set(centers))


def _extract_component_atoms(raw: Dict[str, Any]) -> Optional[Dict[str, List[int]]]:
    """Extract per-component atom map numbers from resolver_provenance.

    Parses the JSON-encoded resolver_provenance field in raw source data
    and returns a mapping from component name (e.g. "furan", "allenamide")
    to lists of atom map numbers on that component.
    """
    rp = raw.get("resolver_provenance")
    if not rp:
        return None
    try:
        data = json.loads(rp) if isinstance(rp, str) else rp
    except (json.JSONDecodeError, TypeError):
        return None

    result: Dict[str, List[int]] = {}
    for key in ("furan_atoms", "allenamide_atoms", "diene_atoms", "dienophile_atoms",
                 "dipole_atoms", "dipolarophile_atoms"):
        if key in data and isinstance(data[key], list):
            result[key.replace("_atoms", "")] = [int(a) for a in data[key]]
    return result if result else None


def _build_clusters(
    new_centers: Sequence[int],
    forming_bonds: Sequence[Tuple[int, int]],
) -> List[ChiralCluster]:
    if not new_centers:
        return []
    return [
        ChiralCluster(
            cluster_id="CL_DR_001",
            map_numbers=list(new_centers),
            forming_bonds=[list(pair) for pair in forming_bonds],
            enumeration_mode="concerted_face_flip",
        )
    ]


def _build_branches(
    *,
    product_smiles: Optional[str],
    fixed_centers: Sequence[int],
    new_centers: Sequence[int],
    clusters: Sequence[ChiralCluster],
    component_atoms: Optional[Dict[str, List[int]]] = None,
) -> List[DRProductBranch]:
    if not _requires_dr_branch(fixed_centers, new_centers):
        return [
            DRProductBranch(
                branch_id="BR_MAJOR",
                pathway_id="primary",
                role="major_or_reference",
                generation_policy="reference_product",
                product_smiles=product_smiles,
                fixed_stereocenters=list(fixed_centers),
            )
        ]

    flip_targets = _select_flip_targets(
        clusters[0].map_numbers if clusters else list(new_centers),
        component_atoms=component_atoms,
    )

    # Non-flipped new chiral centers must be treated as fixed so that
    # generate_stereo_branch_smiles does not accidentally flip them.
    non_flipped_new = sorted(set(int(c) for c in new_centers) - set(flip_targets))
    all_fixed = sorted(set(int(c) for c in fixed_centers) | set(non_flipped_new))

    dr_smiles: Optional[str] = None
    dr_notes: List[str] = []
    if product_smiles and flip_targets:
        try:
            from rph_core.steps.mechanism_classifier.stereo_smiles import (
                generate_stereo_branch_smiles,
            )

            dr_smiles = generate_stereo_branch_smiles(
                mapped_product_smiles=product_smiles,
                flipped_map_numbers=flip_targets,
                fixed_map_numbers=all_fixed,
            )
        except Exception as exc:
            logger.warning(f"DR stereo SMILES generation failed for BR_DR_001: {exc}")

    if dr_smiles is None:
        dr_notes.append(
            "Explicit stereochemical SMILES generation is deferred; "
            "branch identity is defined by flipped map numbers."
        )
    else:
        dr_notes.append(
            "Explicit stereochemical SMILES generated from mapped product SMILES."
        )

    return [
        DRProductBranch(
            branch_id="BR_MAJOR",
            pathway_id="primary",
            role="major_or_reference",
            generation_policy="reference_product",
            product_smiles=product_smiles,
            fixed_stereocenters=list(fixed_centers),
        ),
        DRProductBranch(
            branch_id="BR_DR_001",
            pathway_id="dr_concerted_face_flip_001",
            role="diastereomer_candidate",
            generation_policy="concerted_bridgehead_pair_flip",
            product_smiles=dr_smiles,
            flipped_map_numbers=flip_targets,
            fixed_stereocenters=all_fixed,
            notes=dr_notes,
        ),
    ]


def _requires_dr_branch(fixed_centers: Sequence[int], new_centers: Sequence[int]) -> bool:
    return len(new_centers) >= 2 or (len(new_centers) == 1 and len(fixed_centers) > 0)


def _select_flip_targets(
    map_numbers: Sequence[int],
    component_atoms: Optional[Dict[str, List[int]]] = None,
) -> List[int]:
    """Select chiral center map numbers to flip for DR branch generation.

    When component info from resolver_provenance is available, only
    furan-side chiral centers are selected. Allenamide-side centers are
    geometrically constrained by the ring system and must not be flipped
    independently.

    Without component info, all candidate map numbers are returned (no cap).
    """
    unique = sorted(set(int(v) for v in map_numbers))

    if component_atoms:
        furan_atoms = set(int(a) for a in component_atoms.get("furan", []))
        targets = sorted(furan_atoms & set(unique))
        if targets:
            return targets

    return unique


def _reported_dr(raw: Dict[str, Any]) -> Dict[str, Any]:
    reported: Dict[str, Any] = {}
    for key in ("dr_major", "dr_minor", "de", "ee", "stereo_consistency"):
        value = raw.get(key)
        if value not in (None, ""):
            reported[key] = value
    return reported


def _sorted_pair(pair: Iterable[Any]) -> Tuple[int, int]:
    values = [int(value) for value in pair]
    return (min(values[0], values[1]), max(values[0], values[1]))


def _dedupe_pairs(pairs: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    seen: set[Tuple[int, int]] = set()
    result: List[Tuple[int, int]] = []
    for pair in pairs:
        if pair in seen:
            continue
        seen.add(pair)
        result.append(pair)
    return result
