"""
Step2 Geometry Guard - Topology and Close-Contact Validation
===========================================================

Validates seed and intermediate geometries against product baseline
to detect unintended bond formation or topology drift.

Used by RetroScanner to guard against wrong-topology intermediates
that can arise from naive bond stretching or xTB relaxation.

Author: ReactionProfileHunter Team
Date: 2026-03-14
"""

from typing import Dict, List, Set, Tuple, Optional, Any, Sequence
from pathlib import Path
import numpy as np
import logging

from rph_core.utils.molecular_graph import build_bond_graph
from rph_core.utils.geometry_tools import GeometryUtils
from rph_core.utils.file_io import read_xyz

logger = logging.getLogger(__name__)


class TopologyGuardResult:
    """Result of topology validation."""
    
    def __init__(
        self,
        is_valid: bool,
        new_edges: List[Tuple[int, int, float]],
        lost_edges: List[Tuple[int, int]],
        forming_bonds: Set[Tuple[int, int]],
        graph_scale: float = 1.25
    ):
        self.is_valid = is_valid
        self.new_edges = new_edges  # (i, j, distance)
        self.lost_edges = lost_edges  # (i, j)
        self.forming_bonds = forming_bonds
        self.graph_scale = graph_scale


class RiskyContactResult:
    """Result of risky contact detection."""
    
    def __init__(
        self,
        risky_pairs: List[Tuple[int, int, float, float]],
        # (atom_i, atom_j, current_dist, product_dist)
        threshold_used: float
    ):
        self.risky_pairs = risky_pairs
        self.threshold_used = threshold_used


def compare_graph_topology(
    product_coords: np.ndarray,
    candidate_coords: np.ndarray,
    symbols: List[str],
    forming_bonds: Sequence[Tuple[int, int]],
    graph_scale: float = 1.25,
    min_dist: float = 0.6
) -> TopologyGuardResult:
    """
    Compare bond graph topology between product and candidate geometry.
    
    Expected changes: only the forming_bonds should be DELETED (stretched).
    Any NEW edges or loss of non-forming edges indicates topology drift.
    
    Args:
        product_coords: Product geometry (N, 3)
        candidate_coords: Candidate geometry (N, 3) - seed or xTB output
        symbols: Element symbols list
        forming_bonds: List of (i, j) tuples for bonds being broken
        graph_scale: Covalent radius multiplier for bond detection
        min_dist: Minimum distance to consider (avoid pathological overlaps)
    
    Returns:
        TopologyGuardResult with validation status and details
    """
    # Build bond graphs
    try:
        product_graph = build_bond_graph(
            product_coords, symbols, scale=graph_scale, min_dist=min_dist
        )
        candidate_graph = build_bond_graph(
            candidate_coords, symbols, scale=graph_scale, min_dist=min_dist
        )
    except ValueError as e:
        logger.warning(f"Failed to build bond graph: {e}")
        # Return invalid result on graph construction failure
        return TopologyGuardResult(
            is_valid=False,
            new_edges=[],
            lost_edges=[],
            forming_bonds=set((int(b[0]), int(b[1])) for b in forming_bonds),
            graph_scale=graph_scale
        )
    
    # Convert to edge sets
    product_edges: Set[Tuple[int, int]] = set()
    for i, neighbors in product_graph.items():
        for j in neighbors:
            if i < j:  # Avoid duplicates
                product_edges.add((i, j))
    
    candidate_edges: Set[Tuple[int, int]] = set()
    for i, neighbors in candidate_graph.items():
        for j in neighbors:
            if i < j:
                candidate_edges.add((i, j))
    
    forming_set: Set[Tuple[int, int]] = set(
        (int(b[0]), int(b[1])) for b in forming_bonds
    )
    
    # Find new edges (in candidate but not in product)
    new_edges_raw = candidate_edges - product_edges
    new_edges: List[Tuple[int, int, float]] = []
    for i, j in sorted(new_edges_raw):
        dist = GeometryUtils.calculate_distance(candidate_coords, i, j)
        new_edges.append((i, j, dist))
    
    # Find lost edges (in product but not in candidate, excluding forming bonds)
    lost_edges_raw = product_edges - candidate_edges - forming_set
    lost_edges: List[Tuple[int, int]] = sorted(lost_edges_raw)
    
    # Valid if: no new edges AND no loss of non-forming edges
    is_valid = len(new_edges) == 0 and len(lost_edges) == 0
    
    if not is_valid:
        if new_edges:
            logger.warning(
                f"Topology drift: {len(new_edges)} new edge(s) detected: "
                f"{[(i, j, f'{d:.3f}Å') for i, j, d in new_edges[:3]]}"
            )
        if lost_edges:
            logger.warning(
                f"Topology drift: {len(lost_edges)} non-forming edge(s) lost: "
                f"{lost_edges[:3]}"
            )
    
    return TopologyGuardResult(
        is_valid=is_valid,
        new_edges=new_edges,
        lost_edges=lost_edges,
        forming_bonds=forming_set,
        graph_scale=graph_scale
    )


def detect_risky_contacts(
    product_coords: np.ndarray,
    candidate_coords: np.ndarray,
    symbols: List[str],
    product_graph_scale: float = 1.25,
    near_bond_threshold_ratio: float = 0.85,
    near_bond_abs_max: float = 2.2,
    min_shrink_ratio: float = 0.75,
    max_pairs: int = 6
) -> RiskyContactResult:
    """
    Detect non-product pairs that have moved into near-bond contact.
    
    These are pairs that:
    1. Are NOT bonded in the product graph
    2. Have shrunk materially relative to product distance
    3. Are now within near-bond distance window
    
    Args:
        product_coords: Product geometry (N, 3)
        candidate_coords: Candidate geometry (N, 3)
        symbols: Element symbols
        product_graph_scale: Scale for building product bond graph
        near_bond_threshold_ratio: Ratio of covalent sum to consider "near bond"
        near_bond_abs_max: Absolute distance threshold (Å) for near-bond
        min_shrink_ratio: Pair must have shrunk to < this ratio of product distance
        max_pairs: Maximum number of risky pairs to report
    
    Returns:
        RiskyContactResult with sorted list of risky pairs
    """
    try:
        product_graph = build_bond_graph(
            product_coords, symbols, scale=product_graph_scale, min_dist=0.6
        )
    except ValueError as e:
        logger.warning(f"Failed to build product graph for risk detection: {e}")
        return RiskyContactResult(risky_pairs=[], threshold_used=near_bond_abs_max)
    
    # Get product edges
    product_edges: Set[Tuple[int, int]] = set()
    for i, neighbors in product_graph.items():
        for j in neighbors:
            if i < j:
                product_edges.add((i, j))
    
    # Covalent radii for near-bond threshold calculation
    COVALENT_RADII = {
        'H': 0.31, 'B': 0.85, 'C': 0.76, 'N': 0.71, 'O': 0.66,
        'F': 0.57, 'Si': 1.11, 'P': 1.07, 'S': 1.05,
        'Cl': 1.02, 'Br': 1.20, 'I': 1.39
    }
    
    n_atoms = len(symbols)
    risky_pairs: List[Tuple[int, int, float, float]] = []
    
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            pair = (i, j)
            
            # Skip if already bonded in product
            if pair in product_edges:
                continue
            
            # Calculate distances
            prod_dist = GeometryUtils.calculate_distance(product_coords, i, j)
            cand_dist = GeometryUtils.calculate_distance(candidate_coords, i, j)
            
            # Skip if no shrinkage
            if cand_dist >= prod_dist * min_shrink_ratio:
                continue
            
            # Calculate near-bond threshold for this pair
            r_i = COVALENT_RADII.get(symbols[i], 0.76)
            r_j = COVALENT_RADII.get(symbols[j], 0.76)
            near_bond_thresh = min(
                near_bond_abs_max,
                near_bond_threshold_ratio * (r_i + r_j)
            )
            
            # Check if in near-bond window
            if cand_dist <= near_bond_thresh:
                risky_pairs.append((i, j, cand_dist, prod_dist))
    
    # Sort by severity (closest first)
    risky_pairs.sort(key=lambda x: x[2])
    
    # Limit to max_pairs
    if len(risky_pairs) > max_pairs:
        logger.debug(
            f"Truncated risky pairs from {len(risky_pairs)} to {max_pairs}"
        )
        risky_pairs = risky_pairs[:max_pairs]
    
    return RiskyContactResult(
        risky_pairs=risky_pairs,
        threshold_used=near_bond_abs_max
    )


def compute_min_nonbonded_distance(
    coords: np.ndarray,
    bonded_pairs: Set[Tuple[int, int]],
    forming_bonds: Set[Tuple[int, int]],
    min_considered: float = 1.25
) -> Tuple[float, Optional[Tuple[int, int]]]:
    """
    Compute minimum non-bonded distance in geometry.
    
    Args:
        coords: Geometry array (N, 3)
        bonded_pairs: Set of bonded pairs to exclude
        forming_bonds: Set of forming bonds to exclude
        min_considered: Minimum distance to consider (exclude closer)
    
    Returns:
        (min_distance, pair) where pair is the closest non-bonded pair
    """
    n_atoms = len(coords)
    excluded = bonded_pairs | forming_bonds
    
    min_dist = float('inf')
    min_pair: Optional[Tuple[int, int]] = None
    
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            pair = (i, j)
            if pair in excluded:
                continue
            
            dist = GeometryUtils.calculate_distance(coords, i, j)
            if dist < min_considered:
                continue
            
            if dist < min_dist:
                min_dist = dist
                min_pair = pair
    
    return (min_dist if min_dist != float('inf') else float('nan'), min_pair)


def generate_keepaway_constraints(
    risky_pairs: List[Tuple[int, int, float, float]],
    keep_apart_floor: float = 3.0,
    force_constant: float = 0.5
) -> Dict[str, Any]:
    """
    Generate xTB constraint parameters from risky pairs.
    
    Args:
        risky_pairs: List of (i, j, cand_dist, prod_dist)
        keep_apart_floor: Minimum target distance (Å)
        force_constant: Harmonic constraint force constant
    
    Returns:
        Dict with 'distance_constraints' and 'force_constant'
    """
    distance_constraints: Dict[str, float] = {}
    
    for i, j, cand_dist, prod_dist in risky_pairs:
        # Target: max of product distance or floor
        target = max(prod_dist, keep_apart_floor)
        key = f"{i} {j}"
        distance_constraints[key] = target
    
    return {
        'distance_constraints': distance_constraints,
        'force_constant': force_constant
    }


def check_scan_trajectory(
    *,
    product_coords: np.ndarray,
    symbols: List[str],
    forming_bonds: Sequence[Tuple[int, int]],
    frame_paths: Sequence[Path],
    graph_scale: float = 1.25,
) -> Dict[str, Any]:
    total_frames = len(frame_paths)
    off_path_indices: List[int] = []
    frame_issues: List[Dict[str, Any]] = []
    prev_coords = None
    rmsd_surge_threshold = 0.5  # V5.1 Default RMSD surge threshold 

    for idx, frame_path in enumerate(frame_paths):
        try:
            frame_coords, frame_symbols = read_xyz(Path(frame_path))
            frame_coords_np = np.asarray(frame_coords, dtype=float)
        except Exception as exc:
            off_path_indices.append(idx)
            frame_issues.append({
                "frame_index": idx,
                "reason": f"frame_read_error:{exc}",
            })
            continue

        if len(frame_symbols) != len(symbols):
            off_path_indices.append(idx)
            frame_issues.append({
                "frame_index": idx,
                "reason": "symbol_count_mismatch",
            })
            continue

        # 1. Topology Drift Check 
        guard_result = compare_graph_topology(
            product_coords=product_coords,
            candidate_coords=frame_coords_np,
            symbols=frame_symbols,
            forming_bonds=forming_bonds,
            graph_scale=graph_scale,
        )
        if not guard_result.is_valid:
            if idx not in off_path_indices:
                off_path_indices.append(idx)
            frame_issues.append({
                "frame_index": idx,
                "reason": "topology_drift",
                "new_edges": len(guard_result.new_edges),
                "lost_edges": len(guard_result.lost_edges),
            })

        # 2. RMSD Surge Check (V5.1 Off-Path)
        if prev_coords is not None:
            rmsd_step = GeometryUtils.calculate_rmsd(prev_coords, frame_coords_np)
            if rmsd_step > rmsd_surge_threshold:
                if idx not in off_path_indices:
                    off_path_indices.append(idx)
                frame_issues.append({
                    "frame_index": idx,
                    "reason": f"rmsd_surge:{rmsd_step:.3f}",
                })
        
        prev_coords = frame_coords_np

    return {
        "checked": total_frames > 0,
        "total_frames": total_frames,
        "off_path_indices": off_path_indices,
        "off_path_count": len(off_path_indices),
        "frame_issues": frame_issues,
    }
