"""
Molecular Graph Utilities
=========================

Build molecular topology graphs for fragmenter and analysis operations.

Author: QCcalc Team
Date: 2026-01-27
"""

from typing import Optional, List, Tuple, Dict
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
