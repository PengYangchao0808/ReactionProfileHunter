"""
Intramolecular Tether-Cut Fragmenter
==================================

Fragmenter for intramolecular oxidopyrylium [5+2] systems using tether-cut strategy.

Author: QCcalc Team
Date: 2026-01-27
"""

from dataclasses import dataclass
from typing import List, Tuple, Optional
import logging

import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.molecular_graph import (
    build_bond_graph,
    get_connected_components,
    find_shortest_path
)
from rph_core.utils.geometry_tools import GeometryUtils
from rph_core.utils.fragment_manipulation import h_cap_fragment

logger = logging.getLogger(__name__)


@dataclass
class FragmenterResult:
    """Result of tether-cut fragmentation."""
    fragA_indices: List[int]
    fragB_indices: List[int]
    cut_bond_indices: Tuple[int, int]
    dipole_end_indices: Tuple[int, int]
    alkene_end_indices: Tuple[int, int]
    dipole_core_indices: List[int]
    fragA_coords_R: np.ndarray
    fragB_coords_R: np.ndarray
    fragA_coords_TS: np.ndarray
    fragB_coords_TS: np.ndarray
    fragA_symbols_R: List[str]
    fragB_symbols_R: List[str]
    fragA_symbols_TS: List[str]
    fragB_symbols_TS: List[str]
    status: str
    reason: str = ""


class IntramolecularFragmenter(LoggerMixin):
    """
    Intramolecular fragmenter using tether-cut strategy.

    For oxidopyrylium [5+2] systems:
    - Fragment A: Dipole (5-atom core + substituents, no tether)
    - Fragment B: Dipolarophile + tether

    Uses forming_bonds to identify reaction atoms, then:
    1. Classify alkene vs dipole ends
    2. Find 5-atom dipole core path
    3. Locate and cut tether-dipole connecting bond
    4. H-cap both fragments
    """

    def fragment(
        self,
        reactant_coords: np.ndarray,
        reactant_symbols: List[str],
        ts_coords: np.ndarray,
        ts_symbols: List[str],
        forming_bonds: Tuple[Tuple[int, int], Tuple[int, int]],
        config: dict = None
    ) -> FragmenterResult:
        """
        Execute tether-cut fragmentation algorithm (7 steps).

        Args:
            reactant_coords: Reactant complex coordinates (N, 3) in Å
            reactant_symbols: Element symbols
            ts_coords: Transition state coordinates (N, 3) in Å
            ts_symbols: Element symbols
            forming_bonds: ((u1, v1), (u2, v2)) forming bond atom pairs
            config: Optional config dict (scale, min_dist, etc.)

        Returns:
            FragmenterResult with capped fragment geometries and metadata
        """
        if config is None:
            config = {}

        graph = self._build_graph(
            reactant_coords, reactant_symbols,
            scale=config.get('connectivity_scale', 1.25),
            min_dist=config.get('bond_min_dist_angstrom', 0.6)
        )

        alkene_pair, dipole_pair = self._classify_reaction_ends(
            graph, reactant_coords, forming_bonds
        )
        if alkene_pair is None or dipole_pair is None:
            return FragmenterResult(
                fragA_indices=[], fragB_indices=[],
                cut_bond_indices=(0, 0),
                dipole_end_indices=(0, 0),
                alkene_end_indices=(0, 0),
                dipole_core_indices=[],
                fragA_coords_R=np.array([]), fragB_coords_R=np.array([]),
                fragA_coords_TS=np.array([]), fragB_coords_TS=np.array([]),
                fragA_symbols_R=[], fragB_symbols_R=[],
                fragA_symbols_TS=[], fragB_symbols_TS=[],
                status="failed",
                reason="failed_to_classify_reaction_ends"
            )

        dipole_core = self._find_dipole_core_path(graph, dipole_pair)
        if dipole_core is None:
            return FragmenterResult(
                fragA_indices=[], fragB_indices=[],
                cut_bond_indices=(0, 0),
                dipole_end_indices=dipole_pair,
                alkene_end_indices=alkene_pair,
                dipole_core_indices=[],
                fragA_coords_R=np.array([]), fragB_coords_R=np.array([]),
                fragA_coords_TS=np.array([]), fragB_coords_TS=np.array([]),
                fragA_symbols_R=[], fragB_symbols_R=[],
                fragA_symbols_TS=[], fragB_symbols_TS=[],
                status="failed",
                reason="dipole_path_not_5_atoms"
            )

        cut_bond = self._find_cut_bond(
            graph, dipole_core, alkene_pair, reactant_symbols
        )
        if cut_bond is None:
            return FragmenterResult(
                fragA_indices=[], fragB_indices=[],
                cut_bond_indices=(0, 0),
                dipole_end_indices=dipole_pair,
                alkene_end_indices=alkene_pair,
                dipole_core_indices=dipole_core,
                fragA_coords_R=np.array([]), fragB_coords_R=np.array([]),
                fragA_coords_TS=np.array([]), fragB_coords_TS=np.array([]),
                fragA_symbols_R=[], fragB_symbols_R=[],
                fragA_symbols_TS=[], fragB_symbols_TS=[],
                status="failed",
                reason="failed_to_find_cut_bond"
            )

        fragA_indices, fragB_indices = self._cut_and_get_components(
            graph, cut_bond, dipole_pair, alkene_pair
        )

        if len(fragA_indices) == 0 or len(fragB_indices) == 0:
            return FragmenterResult(
                fragA_indices=fragA_indices, fragB_indices=fragB_indices,
                cut_bond_indices=cut_bond,
                dipole_end_indices=dipole_pair,
                alkene_end_indices=alkene_pair,
                dipole_core_indices=dipole_core,
                fragA_coords_R=np.array([]), fragB_coords_R=np.array([]),
                fragA_coords_TS=np.array([]), fragB_coords_TS=np.array([]),
                fragA_symbols_R=[], fragB_symbols_R=[],
                fragA_symbols_TS=[], fragB_symbols_TS=[],
                status="failed",
                reason="component_validation_failed"
            )

        fragA_coords_R, fragA_symbols_R = self._h_cap_geometry(
            reactant_coords, reactant_symbols, fragA_indices, cut_bond
        )
        fragB_coords_R, fragB_symbols_R = self._h_cap_geometry(
            reactant_coords, reactant_symbols, fragB_indices, cut_bond
        )

        fragA_coords_TS, fragA_symbols_TS = self._h_cap_geometry(
            ts_coords, ts_symbols, fragA_indices, cut_bond
        )
        fragB_coords_TS, fragB_symbols_TS = self._h_cap_geometry(
            ts_coords, ts_symbols, fragB_indices, cut_bond
        )

        return FragmenterResult(
            fragA_indices=fragA_indices,
            fragB_indices=fragB_indices,
            cut_bond_indices=cut_bond,
            dipole_end_indices=dipole_pair,
            alkene_end_indices=alkene_pair,
            dipole_core_indices=dipole_core,
            fragA_coords_R=fragA_coords_R,
            fragB_coords_R=fragB_coords_R,
            fragA_coords_TS=fragA_coords_TS,
            fragB_coords_TS=fragB_coords_TS,
            fragA_symbols_R=fragA_symbols_R,
            fragB_symbols_R=fragB_symbols_R,
            fragA_symbols_TS=fragA_symbols_TS,
            fragB_symbols_TS=fragB_symbols_TS,
            status="ok",
            reason=""
        )

    def _build_graph(self, coords, symbols, scale, min_dist):
        """Step A: Build topology graph using bond graph utilities."""
        return build_bond_graph(coords, symbols, scale=scale, min_dist=min_dist)

    def _classify_reaction_ends(self, graph, coords, forming_bonds):
        """
        Step B: Identify 4 reaction atoms and classify alkene vs dipole ends.

        Returns:
            (alkene_pair, dipole_pair) where each is (i, j) tuple
            Returns (None, None) if classification fails
        """
        reaction_atoms = set()
        for bond in forming_bonds:
            reaction_atoms.add(bond[0])
            reaction_atoms.add(bond[1])
        reaction_atoms = sorted(reaction_atoms)

        alkene_candidates = []
        for i in range(len(reaction_atoms)):
            for j in range(i + 1, len(reaction_atoms)):
                a, b = reaction_atoms[i], reaction_atoms[j]
                path = find_shortest_path(graph, a, b)
                if len(path) == 2:
                    dist = GeometryUtils.calculate_distance(coords, a, b)
                    if dist < 1.55:
                        alkene_candidates.append(((a, b), dist))

        if not alkene_candidates:
            return None, None

        alkene_pair = min(alkene_candidates, key=lambda x: x[1])[0]

        all_atoms = set(reaction_atoms)
        dipole_atoms = all_atoms - set(alkene_pair)
        dipole_pair = tuple(sorted(dipole_atoms))

        return alkene_pair, dipole_pair

    def _find_dipole_core_path(self, graph, dipole_pair):
        """
        Step C: Find 5-atom dipole core path.

        Returns:
            List of 5 atom indices (dipole_core_indices)
            Returns None if path length != 4 edges (5 atoms)
        """
        d1, d2 = dipole_pair
        path = find_shortest_path(graph, d1, d2)

        if len(path) != 5:
            return None

        return path

    def _find_cut_bond(self, graph, dipole_core, alkene_pair, symbols):
        """
        Step D: Find tether-dipole cut bond.

        Returns:
            (i, j) cut bond indices
            Returns None if cut would break a ring
        """
        a, b = alkene_pair

        path_to_a = self._shortest_path_to_set(graph, dipole_core, a)
        path_to_b = self._shortest_path_to_set(graph, dipole_core, b)

        if len(path_to_a) <= len(path_to_b):
            tether_path = path_to_a
        else:
            tether_path = path_to_b

        if len(tether_path) < 2:
            return None

        p0, p1 = tether_path[0], tether_path[1]

        degree_p0 = len(graph[p0])
        degree_p1 = len(graph[p1])

        if degree_p0 >= 2 and degree_p1 >= 2:
            if len(tether_path) >= 3:
                p1, p2 = tether_path[1], tether_path[2]
                if not (len(graph[p1]) >= 2 and len(graph[p2]) >= 2):
                    return (p1, p2)

        return (p0, p1)

    def _shortest_path_to_set(self, graph, source_set, target):
        """Find shortest path from any atom in source_set to target."""
        best_path = None

        for source in source_set:
            path = find_shortest_path(graph, source, target)
            if path:
                if best_path is None or len(path) < len(best_path):
                    best_path = path

        return best_path if best_path else []

    def _cut_and_get_components(self, graph, cut_bond, dipole_pair, alkene_pair):
        """
        Step E: Remove cut bond and extract components.

        Returns:
            (fragA_indices, fragB_indices)
            fragA must contain dipole_pair, fragB must contain alkene_pair
        """
        cut_graph = {k: v.copy() for k, v in graph.items()}

        i, j = cut_bond
        cut_graph[i].remove(j)
        cut_graph[j].remove(i)

        components = get_connected_components(cut_graph)

        if len(components) != 2:
            return [], []

        comp0_set = set(components[0])
        comp1_set = set(components[1])

        dipole_set = set(dipole_pair)
        alkene_set = set(alkene_pair)

        if dipole_set.issubset(comp0_set):
            fragA_indices, fragB_indices = components[0], components[1]
        else:
            fragA_indices, fragB_indices = components[1], components[0]

        if not dipole_set.issubset(set(fragA_indices)):
            return [], []
        if not alkene_set.issubset(set(fragB_indices)):
            return [], []

        return fragA_indices, fragB_indices

    def _h_cap_geometry(self, full_coords, full_symbols, frag_indices, cut_bond, cap_rules=None):
        """
        H-cap a fragment geometry with correct index transformation.

        Args:
            full_coords: Full system coordinates (N, 3)
            full_symbols: Full system element symbols
            frag_indices: Fragment atom indices (global indices)
            cut_bond: (i, j) global indices of cut bond
            cap_rules: Optional bond length by element

        Returns:
            (capped_coords, capped_symbols) - Capped fragment geometry
        """
        if cap_rules is None:
            cap_rules = {"C": 1.09, "N": 1.01, "O": 0.96, "H": 0.76}

        frag_coords = full_coords[frag_indices]
        frag_symbols = [full_symbols[i] for i in frag_indices]

        global_to_local = {g: l for l, g in enumerate(frag_indices)}

        i_global, j_global = cut_bond

        cap_positions = []
        if i_global in global_to_local:
            i_local = global_to_local[i_global]
            cap_dir = full_coords[j_global] - full_coords[i_global]
            cap_positions.append((i_local, cap_dir))
        if j_global in global_to_local:
            j_local = global_to_local[j_global]
            cap_dir = full_coords[i_global] - full_coords[j_global]
            cap_positions.append((j_local, cap_dir))

        capped_coords, capped_symbols = h_cap_fragment(
            frag_coords, frag_symbols, cap_positions, cap_bond_lengths=cap_rules
        )

        return capped_coords, capped_symbols
