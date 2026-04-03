"""
Fragment Cutting Utility
========================

Segment molecules using forming bonds as cut edges for GEDT calculation.
Implements graph connectivity analysis to split TS into two fragments.

Author: RPH Team
Date: 2026-02-02
"""

from pathlib import Path
from typing import Optional, Tuple, List, Set, Dict, Any
import logging
import numpy as np

logger = logging.getLogger(__name__)


def _validate_forming_bonds_indices(
    n_atoms: int,
    forming_bonds: Tuple[Tuple[int, int], ...],
) -> None:
    for pair in forming_bonds:
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            raise ValueError(f"Invalid forming_bonds pair format: {pair}")
        i, j = int(pair[0]), int(pair[1])
        if i < 0 or j < 0:
            raise ValueError(f"forming_bonds contains negative index: ({i}, {j})")
        if i >= n_atoms or j >= n_atoms:
            raise ValueError(
                f"forming_bonds index out of range for atom_count={n_atoms}: ({i}, {j})"
            )


def build_adjacency_matrix(n_atoms: int, bonds: List[Tuple[int, int]]) -> np.ndarray:
    """Build adjacency matrix from bond list.

    Args:
        n_atoms: Number of atoms
        bonds: List of (i, j) bonds (0-indexed)

    Returns:
        (n_atoms, n_atoms) adjacency matrix
    """
    adj = np.zeros((n_atoms, n_atoms), dtype=int)
    for i, j in bonds:
        if 0 <= i < n_atoms and 0 <= j < n_atoms:
            adj[i, j] = 1
            adj[j, i] = 1
    return adj


def get_connected_components(
    adj: np.ndarray,
    cut_edges: List[Tuple[int, int]]
) -> Tuple[List[int], List[int]]:
    """Get two fragments after cutting forming bonds.

    Args:
        adj: Adjacency matrix
        cut_edges: List of forming bond indices to cut

    Returns:
        Tuple of (fragment_A_atoms, fragment_B_atoms) as lists of indices
    """
    n_atoms = adj.shape[0]

    # Create modified adjacency with cut edges removed
    adj_modified = adj.copy()
    for i, j in cut_edges:
        if 0 <= i < n_atoms and 0 <= j < n_atoms:
            adj_modified[i, j] = 0
            adj_modified[j, i] = 0

    # BFS to find connected components starting from atom 0
    visited = set()
    components = []

    for start in range(n_atoms):
        if start in visited:
            continue

        # Start BFS from unvisited atom
        component = []
        queue = [start]
        visited.add(start)

        while queue:
            current = queue.pop(0)
            component.append(current)

            # Check neighbors
            for neighbor in range(n_atoms):
                if adj_modified[current, neighbor] == 1 and neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

        components.append(component)

    # Return first two non-empty components
    non_empty = [c for c in components if c]
    if len(non_empty) >= 2:
        return non_empty[0], non_empty[1]
    elif len(non_empty) == 1:
        return non_empty[0], []
    else:
        return [], []


def _bfs_multi_source_dist(adj: np.ndarray, sources: List[int]) -> List[int]:
    n = adj.shape[0]
    dist = [-1] * n
    queue: List[int] = []
    for s in sources:
        if 0 <= s < n and dist[s] == -1:
            dist[s] = 0
            queue.append(s)
    head = 0
    while head < len(queue):
        cur = queue[head]
        head += 1
        for nb in range(n):
            if adj[cur, nb] == 1 and dist[nb] == -1:
                dist[nb] = dist[cur] + 1
                queue.append(nb)
    return dist


def _graph_distance(adj: np.ndarray, a: int, b: int) -> int:
    if a == b:
        return 0
    dist = _bfs_multi_source_dist(adj, [a])
    d = dist[b] if 0 <= b < len(dist) else -1
    return d if d != -1 else 10**9


def _partition_by_forming_bonds_voronoi(
    adj: np.ndarray,
    forming_bonds: Tuple[Tuple[int, int], ...],
) -> Tuple[Optional[List[int]], Optional[List[int]]]:
    """Fallback partition when forming-bond cut does not disconnect the graph.

    For intramolecular cyclizations, cutting only the forming bonds may not
    produce two connected components (a tether may keep the graph connected).
    In that case, build a deterministic 2-way partition based on graph
    distances to the forming-bond endpoints.
    """
    if forming_bonds is None or len(forming_bonds) != 2:
        return None, None

    (a1, b1), (a2, b2) = forming_bonds

    # Choose pairing that groups atoms that are closer in the covalent graph.
    d_a1_a2 = _graph_distance(adj, a1, a2)
    d_b1_b2 = _graph_distance(adj, b1, b2)
    d_a1_b2 = _graph_distance(adj, a1, b2)
    d_b1_a2 = _graph_distance(adj, b1, a2)

    if d_a1_a2 + d_b1_b2 <= d_a1_b2 + d_b1_a2:
        seeds_a = [a1, a2]
        seeds_b = [b1, b2]
    else:
        seeds_a = [a1, b2]
        seeds_b = [b1, a2]

    dist_a = _bfs_multi_source_dist(adj, seeds_a)
    dist_b = _bfs_multi_source_dist(adj, seeds_b)

    frag_a: List[int] = []
    frag_b: List[int] = []
    for i in range(adj.shape[0]):
        da = dist_a[i]
        db = dist_b[i]

        if da == -1 and db == -1:
            # Disconnected singleton: assign deterministically.
            frag_a.append(i)
        elif db == -1:
            frag_a.append(i)
        elif da == -1:
            frag_b.append(i)
        elif da < db:
            frag_a.append(i)
        elif db < da:
            frag_b.append(i)
        else:
            # Tie: assign by proximity to first seed index.
            frag_a.append(i) if i in seeds_a else frag_b.append(i)

    if not frag_a or not frag_b:
        return None, None
    return frag_a, frag_b


def generate_bonds_from_distance(
    coordinates: np.ndarray,
    threshold: float = 1.8
) -> List[Tuple[int, int]]:
    """Generate bond list from coordinates based on distance.

    Args:
        coordinates: (N, 3) coordinate array
        threshold: Distance threshold in Angstroms

    Returns:
        List of (i, j) bond indices
    """
    n_atoms = len(coordinates)
    bonds = []

    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            dist = np.linalg.norm(coordinates[i] - coordinates[j])
            if dist < threshold:
                bonds.append((i, j))

    return bonds


def generate_bonds_from_distance_with_symbols(
    coordinates: np.ndarray,
    symbols: List[str],
    scale: float = 1.25,
    max_dist: float = 2.0,
) -> List[Tuple[int, int]]:
    """Generate a bond list using covalent radii heuristics.

    A fixed distance threshold tends to over-connect TS complexes and prevents
    fragment separation after cutting forming bonds. Using covalent radii
    reduces spurious inter-fragment edges.

    Args:
        coordinates: (N, 3) coordinate array (Angstrom)
        symbols: Atom symbols length N
        scale: Multiplicative factor on (r_i + r_j)
        max_dist: Hard cap to avoid accidental long bonds

    Returns:
        List of (i, j) bond indices
    """
    # Covalent radii (Angstrom), small subset with safe fallbacks.
    rcov = {
        'H': 0.31,
        'B': 0.85,
        'C': 0.76,
        'N': 0.71,
        'O': 0.66,
        'F': 0.57,
        'P': 1.07,
        'S': 1.05,
        'Cl': 1.02,
        'Br': 1.20,
        'I': 1.39,
    }

    n_atoms = len(coordinates)
    bonds: List[Tuple[int, int]] = []

    for i in range(n_atoms):
        si = symbols[i] if i < len(symbols) else 'C'
        ri = rcov.get(si, 0.76)
        for j in range(i + 1, n_atoms):
            sj = symbols[j] if j < len(symbols) else 'C'
            rj = rcov.get(sj, 0.76)
            thresh = min(max_dist, scale * (ri + rj))
            dist = np.linalg.norm(coordinates[i] - coordinates[j])
            if dist <= thresh:
                bonds.append((i, j))

    return bonds


def identify_fragment_role(
    coordinates: np.ndarray,
    symbols: List[str],
    fragment_a: List[int],
    fragment_b: List[int],
    forming_bonds: Tuple[Tuple[int, int], ...]
) -> Tuple[str, str, str]:
    """Identify which fragment is dipole (electron-rich) vs dipolarophile (electron-poor).

    Uses deterministic heuristics:
    1. Prefer fragment containing hetero atoms (O, N, S) as dipole (oxyallyl pattern)
    2. Prefer fragment with more hetero atoms as dipole
    3. Fallback: fragment containing atom with lowest electronegativity difference

    Args:
        coordinates: (N, 3) coordinate array
        symbols: List of atom symbols
        fragment_a: List of atom indices in fragment A
        fragment_b: List of atom indices in fragment B
        forming_bonds: Tuple of (i, j) forming bond indices

    Returns:
        Tuple of (dipole_fragment_label, dipolarophile_fragment_label, label_note)
        Labels are 'A' or 'B' indicating which fragment is which
    """
    HETERO_ATOMS = {'O', 'N', 'S', 'F', 'Cl', 'Br'}

    def count_hetero_atoms(fragment: List[int]) -> int:
        return sum(1 for i in fragment if symbols[i] in HETERO_ATOMS)

    def has_carbonyl(fragment: List[int]) -> bool:
        for i in fragment:
            if symbols[i] == 'C':
                for j in fragment:
                    if i != j and symbols[j] == 'O':
                        dist = np.linalg.norm(coordinates[i] - coordinates[j])
                        if dist < 1.4:
                            return True
        return False

    hetero_a = count_hetero_atoms(fragment_a)
    hetero_b = count_hetero_atoms(fragment_b)

    # Heuristic 1: O-containing fragment is dipole
    if hetero_a != hetero_b:
        if hetero_a > hetero_b:
            return 'A', 'B', 'hetero_atom_count'
        else:
            return 'B', 'A', 'hetero_atom_count'

    # Heuristic 2: Carbonyl-containing fragment is dipolarophile
    carbonyl_a = has_carbonyl(fragment_a)
    carbonyl_b = has_carbonyl(fragment_b)
    if carbonyl_a != carbonyl_b:
        if carbonyl_a:
            return 'B', 'A', 'carbonyl_in_dipolarophile'
        else:
            return 'A', 'B', 'carbonyl_in_dipolarophile'

    # Heuristic 3: Larger fragment tends to be dipole
    if len(fragment_a) != len(fragment_b):
        if len(fragment_a) > len(fragment_b):
            return 'A', 'B', 'size_preference'
        else:
            return 'B', 'A', 'size_preference'

    # Fallback: Use first atom of forming bond
    if forming_bonds:
        first_bond = forming_bonds[0]
        atom0_fragment = 'A' if first_bond[0] in fragment_a else 'B'
        atom1_fragment = 'A' if first_bond[1] in fragment_a else 'B'

        if atom0_fragment != atom1_fragment:
            return atom0_fragment, atom1_fragment, 'forming_bond_first_atom'

    return 'A', 'B', 'fallback'


class FragmentCutter:
    """Cut molecules along forming bonds for GEDT calculation."""

    def __init__(self, config: Dict[str, Any]):
        """Initialize fragment cutter.

        Args:
            config: Configuration dict with fragment_rule
        """
        self.fragment_rule = config.get('fragment_rule', 'forming_bonds_graph_cut_v1')

    def cut_molecule(
        self,
        coordinates: np.ndarray,
        forming_bonds: Tuple[Tuple[int, int], ...],
        existing_bonds: Optional[List[Tuple[int, int]]] = None,
        symbols: Optional[List[str]] = None,
    ) -> Tuple[Optional[List[int]], Optional[List[int]]]:
        """Cut molecule along forming bonds.

        Args:
            coordinates: (N, 3) coordinate array
            forming_bonds: Tuple of (i, j) forming bond indices
            existing_bonds: Optional list of existing bonds for connectivity

        Returns:
            Tuple of (fragment_A_atoms, fragment_B_atoms)
        """
        n_atoms = len(coordinates)
        _validate_forming_bonds_indices(n_atoms=n_atoms, forming_bonds=forming_bonds)

        # Build initial bonds from geometry
        if existing_bonds:
            bonds = existing_bonds
        elif symbols:
            bonds = generate_bonds_from_distance_with_symbols(coordinates, symbols=symbols)
        else:
            bonds = generate_bonds_from_distance(coordinates, threshold=1.8)

        # Get connected components after cutting forming bonds
        fragment_a, fragment_b = get_connected_components(
            build_adjacency_matrix(n_atoms, bonds) if bonds else build_adjacency_matrix(n_atoms, []),
            list(forming_bonds)
        )

        if not fragment_a or not fragment_b:
            # Fallback: intramolecular tethered systems may not disconnect.
            adj = build_adjacency_matrix(n_atoms, bonds) if bonds else build_adjacency_matrix(n_atoms, [])
            fa, fb = _partition_by_forming_bonds_voronoi(adj, forming_bonds)
            if fa and fb:
                return fa, fb
            logger.warning("Failed to split molecule into two fragments")
            return None, None

        return fragment_a, fragment_b

    def calculate_gedt(
        self,
        charges: List[float],
        fragment_a: List[int],
        fragment_b: List[int],
        dipole_label: str = 'A',
        sign_convention: str = "gedt_positive_if_electron_flows_nuc_to_elec"
    ) -> float:
        """Calculate Global Electron Density Transfer.

        Args:
            charges: List of atomic charges (same order as atoms)
            fragment_a: List of atom indices in fragment A
            fragment_b: List of atom indices in fragment B
            dipole_label: Which fragment is dipole ('A' or 'B')
            sign_convention: Definition of GEDT sign

        Returns:
            GEDT value
        """
        if charges is None or not fragment_a or not fragment_b:
            return float('nan')

        # Identify which fragment is dipole based on label
        dipole = fragment_a if dipole_label == 'A' else fragment_b
        dipolarophile = fragment_b if dipole_label == 'A' else fragment_a

        # Calculate total charge on each fragment
        q_dipole = sum(charges[i] for i in dipole if i < len(charges))
        q_dipolarophile = sum(charges[i] for i in dipolarophile if i < len(charges))

        # Apply sign convention
        if sign_convention == "gedt_positive_if_electron_flows_nuc_to_elec":
            gedt = -q_dipole
        else:
            gedt = abs(q_dipole)

        return gedt


def cut_along_forming_bonds(
    coordinates: np.ndarray,
    forming_bonds: Tuple[Tuple[int, int], ...],
    charges: Optional[List[float]] = None,
    symbols: Optional[List[str]] = None,
    config: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """Convenience function to perform full GEDT analysis.

    Args:
        coordinates: (N, 3) coordinate array
        forming_bonds: Tuple of (i, j) forming bond indices
        charges: Optional list of atomic charges
        symbols: Optional list of atom symbols for labeling
        config: Optional configuration dict

    Returns:
        Dictionary with GEDT results and fragment information
    """
    cfg = config or {
        'fragment_rule': 'forming_bonds_graph_cut_v1',
        'sign_convention': 'gedt_positive_if_electron_flows_nuc_to_elec'
    }

    cutter = FragmentCutter(cfg)
    _validate_forming_bonds_indices(n_atoms=len(coordinates), forming_bonds=forming_bonds)

    # Cut molecule
    fragment_a, fragment_b = cutter.cut_molecule(coordinates, forming_bonds, symbols=symbols)

    result = {
        'fragment_a': fragment_a,
        'fragment_b': fragment_b,
        'gedt_value': float('nan'),
        'gedt_charge_type': 'NONE',
        'gedt_sign_convention': cfg.get('sign_convention', 'default'),
        'fragment_rule': cfg.get('fragment_rule', 'forming_bonds_graph_cut_v1'),
        'gedt_fragment_labeling': 'unknown',
        'q_fragment_dipole': float('nan'),
        'q_fragment_dipolarophile': float('nan'),
    }

    if fragment_a and fragment_b and charges is not None:
        # Get fragment charges
        q_a = sum(charges[i] for i in fragment_a if i < len(charges))
        q_b = sum(charges[i] for i in fragment_b if i < len(charges))

        result['q_fragment_dipole'] = q_a
        result['q_fragment_dipolarophile'] = q_b

        # Identify fragment roles if symbols available
        if symbols:
            dipole_label, _, label_note = identify_fragment_role(
                coordinates, symbols, fragment_a, fragment_b, forming_bonds
            )
            result['gedt_fragment_labeling'] = label_note
        else:
            dipole_label = 'A'  # Fallback
            result['gedt_fragment_labeling'] = 'no_symbols_fallback'

        # Calculate GEDT with proper labeling
        result['gedt_value'] = cutter.calculate_gedt(
            charges, fragment_a, fragment_b,
            dipole_label=dipole_label,
            sign_convention=cfg.get('sign_convention', 'gedt_positive_if_electron_flows_nuc_to_elec')
        )

    return result
