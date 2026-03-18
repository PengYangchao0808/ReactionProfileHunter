"""
Kinematic Bond Stretcher - P2 Graph-Partitioned Displacement
==============================================================

Implements P2-level algorithm: "Fragment/Subgraph-Level Rigid Body Displacement"

Key Innovation:
- Replaces "pair pulling" with topology-aware rigid body translation
- Eliminates internal stress by preserving fragment conformations
- Classifies reactions as Type-A (intermolecular) or Type-B (intramolecular)

Author: ReactionProfileHunter Team
Date: 2026-03-14
"""

import logging
from typing import Any, Dict, List, Optional, Tuple
from dataclasses import dataclass
import numpy as np

from rph_core.utils.molecular_graph import (
    build_bond_graph,
    identify_rigid_fragments,
    compute_repulsion_vector,
)

logger = logging.getLogger(__name__)


@dataclass
class KinematicParams:
    """P2 Kinematic displacement parameters."""
    target_start_distance: float = 2.2  # Target for seed generation
    fallback_to_legacy: bool = True     # Use old stretch_bonds if P2 fails


class KinematicStretcher:
    """
    P2 Graph-Partitioned Kinematic Displacement Engine

    Replaces naive "pair pulling" with topology-aware displacement:
    1. Analyze molecular graph and cut forming bonds
    2. Classify as intermolecular (2 fragments) or intramolecular (1 fragment)
    3. Apply appropriate displacement strategy:
       - Type-A (inter): Rigid body translation of entire fragment
       - Type-B (intra): Partial displacement with tether preservation
    """

    def __init__(self, params: Optional[KinematicParams] = None):
        """
        Initialize kinematic stretcher.

        Args:
            params: Kinematic parameters, uses defaults if None
        """
        self.params = params or KinematicParams()
        logger.info(
            f"P2 KinematicStretcher initialized: target={self.params.target_start_distance}Å"
        )

    def kinematic_stretch(
        self,
        coords: np.ndarray,
        symbols: List[str],
        forming_bonds: List[Tuple[int, int]],
        target_distances: Optional[List[float]] = None,
        fix_center_of_mass: bool = True
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """
        Main entry: Perform P2 kinematic displacement.

        Args:
            coords: Atomic coordinates (N, 3) in Å
            symbols: Element symbols list
            forming_bonds: List of (i, j) bonds being formed
            target_distances: Target distance for each bond (default: target_start_distance)
            fix_center_of_mass: Whether to preserve overall COM

        Returns:
            Tuple of (new_coords, metadata)
            - new_coords: Displaced coordinates (N, 3)
            - metadata: Dict with 'type', 'strategy', 'displacement_info'

        Raises:
            RuntimeError: If P2 fails and fallback is disabled
        """
        if not forming_bonds:
            logger.warning("No forming bonds provided, returning original coordinates")
            return coords.copy(), {"type": "none", "reason": "no_bonds"}

        if target_distances is None:
            target_distances = [self.params.target_start_distance] * len(forming_bonds)

        bonds_with_targets = list(zip(forming_bonds, target_distances))

        try:
            # Analyze topology
            fragment_info = identify_rigid_fragments(
                coords=coords,
                symbols=symbols,
                forming_bonds=forming_bonds
            )

            reaction_type = fragment_info['type']

            if reaction_type == 'inter':
                # Type-A: Intermolecular - pure rigid body translation
                new_coords, disp_info = self._rigid_body_translate(
                    coords=coords,
                    fragment_info=fragment_info,
                    bonds_with_targets=bonds_with_targets
                )
                metadata = {
                    'type': 'inter',
                    'strategy': 'rigid_body_translation',
                    'frag_A_size': len(fragment_info['frag_A']),
                    'frag_B_size': len(fragment_info['frag_B']),
                    'displacement': disp_info
                }
                logger.info(
                    f"P2 Type-A (inter): Rigid translation of {len(fragment_info['frag_A'])} atoms, "
                    f"fixed {len(fragment_info['frag_B'])} atoms"
                )

            else:
                # Type-B: Intramolecular - tether-aware displacement
                new_coords, disp_info = self._tether_displacement(
                    coords=coords,
                    fragment_info=fragment_info,
                    bonds_with_targets=bonds_with_targets
                )
                metadata = {
                    'type': 'intra',
                    'strategy': 'tether_aware_displacement',
                    'n_tether_paths': len(fragment_info.get('tether_paths', [])),
                    'displacement': disp_info
                }
                logger.info(
                    f"P2 Type-B (intra): Tether-aware displacement with "
                    f"{len(fragment_info.get('tether_paths', []))} paths"
                )

            # Preserve center of mass if requested
            if fix_center_of_mass:
                original_com = np.mean(coords, axis=0)
                new_com = np.mean(new_coords, axis=0)
                new_coords += (original_com - new_com)

            return new_coords, metadata

        except Exception as exc:
            if self.params.fallback_to_legacy:
                logger.warning(
                    f"P2 kinematic stretch failed: {exc}. Falling back to legacy stretch_bonds."
                )
                from .bond_stretcher import BondStretcher
                legacy = BondStretcher()
                new_coords = legacy.stretch_bonds(
                    coords,
                    [((b[0], b[1]), t) for b, t in bonds_with_targets],
                    fix_center_of_mass=fix_center_of_mass
                )
                return new_coords, {
                    'type': 'fallback',
                    'strategy': 'legacy_pair_pulling',
                    'error': str(exc)
                }
            else:
                raise RuntimeError(f"P2 kinematic stretch failed: {exc}") from exc

    def _rigid_body_translate(
        self,
        coords: np.ndarray,
        fragment_info: Dict[str, Any],
        bonds_with_targets: List[Tuple[Tuple[int, int], float]]
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """
        Type-A: Intermolecular reaction - pure rigid body translation.

        Strategy:
        1. Compute repulsion vector from interface centroids
        2. Move all atoms in fragment_A along this vector
        3. Fragment_B remains fixed
        4. This preserves internal conformations perfectly

        Args:
            coords: Atomic coordinates (N, 3)
            fragment_info: From identify_rigid_fragments()
            bonds_with_targets: Bonds with their target distances

        Returns:
            Tuple of (new_coords, displacement_info)
        """
        frag_A = fragment_info['frag_A']
        frag_B = fragment_info['frag_B']
        interface_A = fragment_info['interface_A']
        interface_B = fragment_info['interface_B']

        new_coords = coords.copy()

        # Compute repulsion vector from interface centroids
        repulsion_vec = compute_repulsion_vector(coords, interface_A, interface_B)

        # Calculate required displacement magnitude
        # Current distance between interface centroids
        center_A_current = np.mean(coords[interface_A], axis=0)
        center_B = np.mean(coords[interface_B], axis=0)
        current_interface_dist = np.linalg.norm(center_A_current - center_B)

        # Target interface distance: average of target bond distances
        target_interface_dist = np.mean([t for _, t in bonds_with_targets])

        # Displacement needed
        displacement_mag = target_interface_dist - current_interface_dist

        # Apply rigid body translation to fragment_A
        displacement = displacement_mag * repulsion_vec
        new_coords[frag_A] += displacement

        disp_info = {
            'vector': repulsion_vec.tolist(),
            'magnitude': float(displacement_mag),
            'atoms_moved': len(frag_A),
            'atoms_fixed': len(frag_B)
        }

        return new_coords, disp_info

    def _tether_displacement(
        self,
        coords: np.ndarray,
        fragment_info: Dict[str, Any],
        bonds_with_targets: List[Tuple[Tuple[int, int], float]]
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """
        Type-B: Intramolecular reaction - tether-aware partial displacement.

        Strategy:
        1. Identify rigid ends near forming bonds
        2. Compute repulsion vector for each forming bond pair
        3. Apply graded displacement: full at ends, attenuated through tether
        4. This minimizes distortion while achieving target distances

        Args:
            coords: Atomic coordinates (N, 3)
            fragment_info: From identify_rigid_fragments()
            bonds_with_targets: Bonds with their target distances

        Returns:
            Tuple of (new_coords, displacement_info)
        """
        new_coords = coords.copy()
        tether_paths = fragment_info.get('tether_paths', [])
        rigid_ends = fragment_info.get('rigid_ends', [])

        total_displacement = np.zeros_like(coords)
        displacement_weights = np.zeros(len(coords))

        # For each forming bond, compute displacement for its atoms
        for (bond, target_dist) in bonds_with_targets:
            i, j = bond

            # Current and target vectors
            vec_ij = coords[j] - coords[i]
            current_dist = np.linalg.norm(vec_ij)

            if current_dist < 1e-10:
                logger.warning(f"Atoms {i} and {j} are coincident, using arbitrary direction")
                unit_vec = np.array([1.0, 0.0, 0.0])
            else:
                unit_vec = vec_ij / current_dist

            # Displacement magnitude needed
            delta = target_dist - current_dist

            # Apply displacement to bond atoms
            # Split displacement equally between both atoms
            disp_i = -0.5 * delta * unit_vec
            disp_j = 0.5 * delta * unit_vec

            total_displacement[i] += disp_i
            total_displacement[j] += disp_j
            displacement_weights[i] += 1.0
            displacement_weights[j] += 1.0

        # Apply attenuation through tether paths
        # Atoms in tether get reduced displacement
        if tether_paths:
            for path in tether_paths:
                if len(path) > 2:
                    # Attenuate displacement for middle atoms in path
                    for idx, atom in enumerate(path[1:-1], start=1):
                        # Linear attenuation: center atoms get 30% displacement
                        attenuation = 0.3 + 0.7 * abs(idx - len(path)/2) / (len(path)/2)
                        total_displacement[atom] *= attenuation

        # Apply weighted displacements
        for atom_idx in range(len(coords)):
            if displacement_weights[atom_idx] > 0:
                new_coords[atom_idx] += (
                    total_displacement[atom_idx] / displacement_weights[atom_idx]
                )

        disp_info = {
            'n_bonds': len(bonds_with_targets),
            'n_tether_paths': len(tether_paths),
            'n_rigid_ends': len(rigid_ends),
            'max_displacement': float(np.max(np.linalg.norm(total_displacement, axis=1)))
        }

        return new_coords, disp_info


def kinematic_stretch(
    coords: np.ndarray,
    symbols: List[str],
    forming_bonds: List[Tuple[int, int]],
    target_distances: Optional[List[float]] = None,
    fix_center_of_mass: bool = True,
    fallback_to_legacy: bool = True
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Module-level function: P2 kinematic displacement.

    Args:
        coords: Atomic coordinates (N, 3)
        symbols: Element symbols list
        forming_bonds: List of (i, j) bonds being formed
        target_distances: Target distances for each bond
        fix_center_of_mass: Preserve overall COM
        fallback_to_legacy: Use old stretch_bonds if P2 fails

    Returns:
        Tuple of (new_coords, metadata)
    """
    params = KinematicParams(
        target_start_distance=target_distances[0] if target_distances else 2.2,
        fallback_to_legacy=fallback_to_legacy
    )
    stretcher = KinematicStretcher(params)
    return stretcher.kinematic_stretch(
        coords=coords,
        symbols=symbols,
        forming_bonds=forming_bonds,
        target_distances=target_distances,
        fix_center_of_mass=fix_center_of_mass
    )