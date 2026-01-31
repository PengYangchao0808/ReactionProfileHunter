"""
Intramolecular Tether-Cut Fragmenter
==================================

Fragmenter for intramolecular oxidopyrylium [5+2] systems using tether-cut strategy.

Author: QCcalc Team
Date: 2026-01-27
"""

from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Any
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
from rph_core.utils.semantic_slicer import semantic_slice, RDKIT_AVAILABLE

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
    fragA_coords_capped: Optional[np.ndarray] = None
    fragB_coords_capped: Optional[np.ndarray] = None
    fragA_symbols_capped: Optional[List[str]] = None
    fragB_symbols_capped: Optional[List[str]] = None
    status: str = "failed"
    reason: str = ""
    fragA_charge: Optional[int] = None
    fragA_mult: Optional[int] = None
    fragB_charge: Optional[int] = None
    fragB_mult: Optional[int] = None


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
        config: Optional[Dict[str, Any]] = None
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

        charge_mult = config.get("charge_multiplicity", None)
        if charge_mult is None:
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
                reason="missing_charge_multiplicity"
            )

        try:
            fragA_charge = int(charge_mult["fragA"]["charge"])
            fragA_mult = int(charge_mult["fragA"]["multiplicity"])
            fragB_charge = int(charge_mult["fragB"]["charge"])
            fragB_mult = int(charge_mult["fragB"]["multiplicity"])
        except Exception:
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
                reason="invalid_charge_multiplicity"
            )

        use_semantic_slicer = config.get('use_semantic_slicer', True)
        semantic_cut_bond = None
        
        if use_semantic_slicer and RDKIT_AVAILABLE:
            self.logger.debug("ASM 2.0: Attempting SemanticSlicer for cut bond detection")
            slice_result = semantic_slice(
                atoms=list(reactant_symbols),
                coords=reactant_coords,
                forming_bonds=forming_bonds,
                verbose=config.get('verbose', False)
            )
            if slice_result.success and slice_result.cut_bond:
                semantic_cut_bond = slice_result.cut_bond
                self.logger.info(
                    f"ASM 2.0: SemanticSlicer succeeded, cut_bond={semantic_cut_bond}"
                )
            else:
                self.logger.warning(
                    f"ASM 2.0: SemanticSlicer failed ({slice_result.reason}), "
                    "falling back to geometric algorithm"
                )
        elif use_semantic_slicer and not RDKIT_AVAILABLE:
            self.logger.warning("ASM 2.0: RDKit not available, using geometric algorithm")

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
                fragA_charge=fragA_charge,
                fragA_mult=fragA_mult,
                fragB_charge=fragB_charge,
                fragB_mult=fragB_mult,
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
                fragA_charge=fragA_charge,
                fragA_mult=fragA_mult,
                fragB_charge=fragB_charge,
                fragB_mult=fragB_mult,
                status="failed",
                reason="dipole_path_not_5_atoms"
            )

        if semantic_cut_bond is not None:
            cut_bond = semantic_cut_bond
            self.logger.debug(f"Using SemanticSlicer cut_bond: {cut_bond}")
        else:
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
                fragA_charge=fragA_charge,
                fragA_mult=fragA_mult,
                fragB_charge=fragB_charge,
                fragB_mult=fragB_mult,
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
                fragA_charge=fragA_charge,
                fragA_mult=fragA_mult,
                fragB_charge=fragB_charge,
                fragB_mult=fragB_mult,
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
            fragA_coords_capped=fragA_coords_TS,
            fragB_coords_capped=fragB_coords_TS,
            fragA_symbols_capped=fragA_symbols_TS,
            fragB_symbols_capped=fragB_symbols_TS,
            fragA_charge=fragA_charge,
            fragA_mult=fragA_mult,
            fragB_charge=fragB_charge,
            fragB_mult=fragB_mult,
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
        Step C: Find dipole core path.

        Returns:
            List of atom indices (dipole_core_indices)
            Returns None if path length not in (3, 5)
        """
        d1, d2 = dipole_pair
        path = find_shortest_path(graph, d1, d2)

        if len(path) not in (3, 5):
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
        Houk Update: Includes recursive pruning to handle ghost bonds that
        keep fragments connected after the main tether is cut.
        """
        cut_graph = {k: v.copy() for k, v in graph.items()}
        i, j = cut_bond

        if j in cut_graph[i]:
            cut_graph[i].remove(j)
        if i in cut_graph[j]:
            cut_graph[j].remove(i)

        anchor_a = dipole_pair[0]
        anchor_b = alkene_pair[0]

        max_attempts = 5
        attempt = 0

        while attempt < max_attempts:
            components = get_connected_components(cut_graph)

            if len(components) == 2:
                return self._assign_fragments(components, dipole_pair, alkene_pair)

            if len(components) == 1:
                leak_path = find_shortest_path(cut_graph, anchor_a, anchor_b)
                if not leak_path or len(leak_path) < 2:
                    break

                mid_index = (len(leak_path) - 2) // 2
                u = leak_path[mid_index]
                v = leak_path[mid_index + 1]

                if v in cut_graph[u]:
                    cut_graph[u].remove(v)
                if u in cut_graph[v]:
                    cut_graph[v].remove(u)

                attempt += 1
                continue

            return self._merge_shattered_components(components, dipole_pair, alkene_pair)

        return [], []

    def _assign_fragments(self, components, dipole_pair, alkene_pair):
        """Assign Fragment A/B based on dipole/alkene endpoints."""
        dipole_set = set(dipole_pair)
        alkene_set = set(alkene_pair)

        comp0_set = set(components[0])

        if not dipole_set.isdisjoint(comp0_set):
            fragA_indices = components[0]
            fragB_indices = components[1]
        else:
            fragA_indices = components[1]
            fragB_indices = components[0]

        if not alkene_set.issubset(set(fragB_indices)):
            return [], []

        return fragA_indices, fragB_indices

    def _merge_shattered_components(self, components, dipole_pair, alkene_pair):
        """
        Handle cases where system breaks into >2 pieces.
        Orphans are merged into the dipole-side fragment.
        """
        dipole_set = set(dipole_pair)
        alkene_set = set(alkene_pair)

        main_a = None
        main_b = None
        orphans = []

        for comp in components:
            comp_set = set(comp)
            if not dipole_set.isdisjoint(comp_set):
                main_a = comp
            elif not alkene_set.isdisjoint(comp_set):
                main_b = comp
            else:
                orphans.append(comp)

        if main_a is None or main_b is None:
            return [], []

        final_a = set(main_a)
        final_b = set(main_b)

        for orphan in orphans:
            final_a.update(orphan)

        return sorted(final_a), sorted(final_b)

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

        def _cap_direction(atom_global: int, partner_global: int) -> np.ndarray:
            """
            Houk 修正: H 放置方向沿切断键向量的**反方向**延伸
            原公式: H_pos = r_atom + k * (r_partner - r_atom) → 撞车
            正确:   H_pos = r_atom + k * (r_atom - r_partner) → 反向
            """
            bond_vec = full_coords[partner_global] - full_coords[atom_global]
            norm = np.linalg.norm(bond_vec)
            if norm < 1e-8:
                raise ValueError("Cut bond atoms are coincident")
            return -bond_vec / norm

        if i_global in global_to_local:
            i_local = global_to_local[i_global]
            cap_dir = _cap_direction(i_global, j_global)
            cap_positions.append((i_local, cap_dir))
        if j_global in global_to_local:
            j_local = global_to_local[j_global]
            cap_dir = _cap_direction(j_global, i_global)
            cap_positions.append((j_local, cap_dir))

        capped_coords, capped_symbols = h_cap_fragment(
            frag_coords, frag_symbols, cap_positions, cap_bond_lengths=cap_rules
        )

        self._validate_hcap_geometry(capped_coords, capped_symbols)

        return capped_coords, capped_symbols

    def _validate_hcap_geometry(
        self,
        coords: np.ndarray,
        symbols: List[str],
        min_distance: float = 0.7
    ) -> None:
        n_atoms = len(coords)
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < min_distance:
                    raise RuntimeError(
                        f"H-cap collision: {symbols[i]}({i})-{symbols[j]}({j}) = {dist:.2f} Å < {min_distance} Å"
                    )
