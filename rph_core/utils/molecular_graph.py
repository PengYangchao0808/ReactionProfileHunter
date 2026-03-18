"""
Molecular Graph Utilities
=========================

Build molecular topology graphs for fragmenter and analysis operations.

Author: QCcalc Team
Date: 2026-01-27
"""

from typing import Any, Optional, List, Tuple, Dict, Set
import logging

import numpy as np

from rph_core.utils.geometry_tools import GeometryUtils

logger = logging.getLogger(__name__)


# Covalent radii table (Å)
COVALENT_RADII = {
    'H': 0.31, 'B': 0.85, 'C': 0.76, 'N': 0.71, 'O': 0.66,
    'F': 0.57, 'Si': 1.11, 'P': 1.07, 'S': 1.05,
    'Cl': 1.02, 'Br': 1.20, 'I': 1.39
}


def build_bond_graph(
    coords: np.ndarray,
    symbols: List[str],
    scale: float = 1.25,
    min_dist: float = 0.6
) -> Dict[int, List[int]]:
    """
    Build adjacency list graph using covalent radius-based bond heuristic.

    Args:
        coords: (N, 3) coordinates in Å
        symbols: List of element symbols
        scale: Distance threshold multiplier
        min_dist: Minimum distance to ignore (pathological overlaps)

    Returns:
        Adjacency list: {atom_idx: [neighbor_indices]}

    Raises:
        ValueError: If unknown element encountered
    """
    n_atoms = len(coords)
    graph = {i: [] for i in range(n_atoms)}

    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            d_ij = GeometryUtils.calculate_distance(coords, i, j)

            if d_ij < min_dist:
                continue

            r_i = COVALENT_RADII.get(symbols[i])
            r_j = COVALENT_RADII.get(symbols[j])

            if r_i is None or r_j is None:
                raise ValueError(f"Unknown element radius: {symbols[i]}, {symbols[j]}")

            threshold = scale * (r_i + r_j)

            if d_ij <= threshold:
                graph[i].append(j)
                graph[j].append(i)

    return graph


def get_connected_components(graph: Dict[int, List[int]]) -> List[List[int]]:
    """
    Extract all connected components from adjacency list graph using BFS.

    Args:
        graph: Adjacency list

    Returns:
        List of components, each is a list of atom indices
    """
    visited = set()
    components = []

    for atom in graph:
        if atom not in visited:
            component = []
            queue = [atom]
            visited.add(atom)

            while queue:
                node = queue.pop(0)
                component.append(node)

                for neighbor in graph[node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)

            components.append(component)

    return components


def find_shortest_path(
    graph: Dict[int, List[int]],
    start: int,
    end: int
) -> List[int]:
    """
    Find shortest path between two atoms using BFS.

    Args:
        graph: Adjacency list
        start: Start atom index
        end: End atom index

    Returns:
        List of atom indices representing path (inclusive)
        Returns empty list if no path exists
    """
    if start == end:
        return [start]

    queue = [start]
    parent: Dict[int, Optional[int]] = {start: None}
    visited = {start}

    while queue:
        node = queue.pop(0)

        if node == end:
            path = []
            while node is not None:
                path.append(node)
                node = parent[node]
            return path[::-1]

        for neighbor in graph[node]:
            if neighbor not in visited:
                visited.add(neighbor)
                parent[neighbor] = node
                queue.append(neighbor)

    return []


def find_cycles(graph: Dict[int, List[int]]) -> List[List[int]]:
    """
    Find cycles in the graph (for ring detection).

    Args:
        graph: Adjacency list

    Returns:
        List of cycles, each is a list of atom indices
    """
    visited = set()
    cycles = []

    def dfs(node, parent, path):
        visited.add(node)
        path.append(node)

        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(neighbor, node, path)
            elif neighbor != parent and neighbor in path:
                cycle_start = path.index(neighbor)
                cycle = path[cycle_start:] + [neighbor]
                if len(cycle) >= 3:
                    min_idx = cycle.index(min(cycle))
                    normalized = cycle[min_idx:] + cycle[:min_idx]
                    if normalized not in cycles:
                        cycles.append(normalized)

        path.pop()

    for node in graph:
        if node not in visited:
            dfs(node, -1, [])

    return cycles


def get_bond_distance(
    coords: np.ndarray,
    graph: Dict[int, List[int]],
    i: int,
    j: int
) -> float:
    """
    Get bond distance between two atoms if they are bonded.

    Args:
        coords: (N, 3) coordinates in Å
        graph: Adjacency list
        i: First atom index
        j: Second atom index

    Returns:
        Distance in Å

    Raises:
        ValueError: If atoms are not bonded
    """
    if j not in graph.get(i, []):
        raise ValueError(f"Atoms {i} and {j} are not bonded")

    return GeometryUtils.calculate_distance(coords, i, j)


def graph_cut_and_components(
    graph: Dict[int, List[int]],
    bonds_to_remove: List[Tuple[int, int]]
) -> Tuple[Dict[int, List[int]], List[List[int]]]:
    """
    Cut edges from graph and return resulting components.

    This is the core algorithm for P2 kinematic displacement:
    1. Remove specified bonds (forming_bonds) from the graph
    2. Compute connected components of the resulting graph
    3. Classify as intermolecular (2 components) or intramolecular (1 component)

    Args:
        graph: Original adjacency list
        bonds_to_remove: List of (i, j) bonds to cut

    Returns:
        Tuple of (cut_graph, connected_components)
        - cut_graph: Graph with bonds removed
        - connected_components: List of component atom lists
    """
    # Create mutable copy of graph
    cut_graph = {atom: list(neighbors) for atom, neighbors in graph.items()}

    # Remove specified bonds
    for i, j in bonds_to_remove:
        if j in cut_graph.get(i, []):
            cut_graph[i].remove(j)
        if i in cut_graph.get(j, []):
            cut_graph[j].remove(i)

    # Compute connected components
    components = get_connected_components(cut_graph)

    return cut_graph, components


def identify_rigid_fragments(
    coords: np.ndarray,
    symbols: List[str],
    forming_bonds: List[Tuple[int, int]],
    graph_scale: float = 1.25
) -> Dict[str, Any]:
    """
    Identify rigid fragments and reaction type for P2 kinematic displacement.

    Algorithm:
    1. Build product bond graph
    2. Cut forming_bonds edges
    3. Analyze connected components
       - 2 components: Type-A (intermolecular) - pure rigid body translation
       - 1 component: Type-B (intramolecular) - tether relaxation needed

    Args:
        coords: (N, 3) atomic coordinates in Å
        symbols: List of element symbols
        forming_bonds: List of (i, j) bonds being formed
        graph_scale: Covalent radius multiplier for bond detection

    Returns:
        Dict with:
        - 'type': 'inter' or 'intra'
        - 'frag_A': List of atom indices in fragment A (for inter)
        - 'frag_B': List of atom indices in fragment B (for inter)
        - 'interface_A': List of atoms in fragment A at forming bonds (for inter)
        - 'interface_B': List of atoms in fragment B at forming bonds (for inter)
        - 'tether_paths': List of paths between forming bonds (for intra)
        - 'rigid_ends': List of (atom_set, center) tuples for ends (for intra)
        - 'n_components': Number of connected components

    Raises:
        ValueError: If more than 2 components detected (unexpected topology)
    """
    # Build graph and cut forming bonds
    graph = build_bond_graph(coords, symbols, scale=graph_scale)
    cut_graph, components = graph_cut_and_components(graph, forming_bonds)

    n_components = len(components)

    if n_components > 2:
        raise ValueError(
            f"Unexpected topology: {n_components} components after cutting {len(forming_bonds)} bonds. "
            f"Forming bonds: {forming_bonds}"
        )

    result: Dict[str, Any] = {
        'n_components': n_components,
        'forming_bonds': forming_bonds,
        'all_atoms': list(range(len(coords))),
    }

    if n_components == 2:
        # Type-A: Intermolecular reaction
        # Sort by size: larger = fragment A (moves), smaller = fragment B (fixed)
        components_sorted = sorted(components, key=len, reverse=True)
        frag_A = components_sorted[0]
        frag_B = components_sorted[1]

        # Identify interface atoms (those involved in forming bonds)
        interface_A = set()
        interface_B = set()
        for i, j in forming_bonds:
            if i in frag_A:
                interface_A.add(i)
                interface_B.add(j)
            else:
                interface_A.add(j)
                interface_B.add(i)

        result['type'] = 'inter'
        result['frag_A'] = frag_A
        result['frag_B'] = frag_B
        result['interface_A'] = list(interface_A)
        result['interface_B'] = list(interface_B)

    else:
        # Type-B: Intramolecular reaction (1 component)
        # Find paths between forming bond atoms to identify tether
        all_tether_paths = []
        rigid_ends = []

        for bond in forming_bonds:
            i, j = bond
            # Find shortest path in cut graph (will go through tether)
            path = find_shortest_path(cut_graph, i, j)
            if path:
                all_tether_paths.append(path)

        # Identify rigid ends: atoms near forming bonds but not in tether
        for bond in forming_bonds:
            i, j = bond
            # Define rigid end as 2-hop neighbors of forming bond atoms
            # (not crossing the tether)
            end_atoms_i = _get_local_subgraph(graph, cut_graph, i, depth=2)
            end_atoms_j = _get_local_subgraph(graph, cut_graph, j, depth=2)

            # Compute centers
            center_i = np.mean(coords[list(end_atoms_i)], axis=0) if end_atoms_i else coords[i]
            center_j = np.mean(coords[list(end_atoms_j)], axis=0) if end_atoms_j else coords[j]

            rigid_ends.append((end_atoms_i, center_i))
            rigid_ends.append((end_atoms_j, center_j))

        result['type'] = 'intra'
        result['tether_paths'] = all_tether_paths
        result['rigid_ends'] = rigid_ends

    return result


def _get_local_subgraph(
    original_graph: Dict[int, List[int]],
    cut_graph: Dict[int, List[int]],
    start_atom: int,
    depth: int = 2
) -> Set[int]:
    """
    Get local subgraph around start_atom up to specified depth.
    Uses original graph but respects cut graph connectivity.

    Args:
        original_graph: Full bond graph
        cut_graph: Graph with forming bonds removed
        start_atom: Central atom
        depth: Number of hops to include

    Returns:
        Set of atoms in local subgraph
    """
    visited = {start_atom}
    current_level = {start_atom}

    for _ in range(depth):
        next_level = set()
        for atom in current_level:
            # Use original graph for connectivity, but stay within cut component
            for neighbor in original_graph.get(atom, []):
                # Check if neighbor is reachable in cut graph (same component)
                if neighbor not in visited:
                    # Verify connectivity via BFS in cut graph
                    if _is_reachable(cut_graph, start_atom, neighbor):
                        visited.add(neighbor)
                        next_level.add(neighbor)
        current_level = next_level

    return visited


def _is_reachable(graph: Dict[int, List[int]], start: int, end: int) -> bool:
    """Check if end is reachable from start in graph using BFS."""
    if start == end:
        return True

    visited = {start}
    queue = [start]

    while queue:
        node = queue.pop(0)
        for neighbor in graph.get(node, []):
            if neighbor == end:
                return True
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)

    return False


def compute_repulsion_vector(
    coords: np.ndarray,
    interface_A: List[int],
    interface_B: List[int]
) -> np.ndarray:
    """
    Compute repulsion vector between two reaction interfaces.

    Algorithm:
    1. Compute geometric center of interface A (C_A)
    2. Compute geometric center of interface B (C_B)
    3. Return normalized vector from C_B to C_A

    This uses the reaction interface centroids rather than individual bond
    directions, avoiding conflicts when multiple bonds form simultaneously.

    Args:
        coords: (N, 3) atomic coordinates
        interface_A: Atom indices in fragment A interface
        interface_B: Atom indices in fragment B interface

    Returns:
        Normalized 3D vector pointing from B to A
    """
    if not interface_A or not interface_B:
        raise ValueError("Empty interface provided")

    center_A = np.mean(coords[interface_A], axis=0)
    center_B = np.mean(coords[interface_B], axis=0)

    vec = center_A - center_B
    norm = np.linalg.norm(vec)

    if norm < 1e-10:
        logger.warning("Interface centers are too close, using arbitrary direction")
        return np.array([1.0, 0.0, 0.0])

    return vec / norm
