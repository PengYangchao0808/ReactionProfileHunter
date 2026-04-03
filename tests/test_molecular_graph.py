"""
Tests for molecular_graph module
===============================

Unit tests for molecular topology graph utilities.

Author: QCcalc Team
Date: 2026-01-27
"""

import pytest
import numpy as np

from rph_core.utils.molecular_graph import (
    build_bond_graph,
    get_connected_components,
    find_shortest_path,
    find_cycles,
    get_bond_distance
)


@pytest.fixture
def methane_coords() -> np.ndarray:
    """Methane (CH4) coordinates - tetrahedral geometry."""
    return np.array([
        [0.0, 0.0, 0.0],
        [1.09, 0.0, 0.0],
        [-0.545, 0.944, 0.0],
        [-0.545, -0.472, 0.816],
        [-0.545, -0.472, -0.816]
    ])


@pytest.fixture
def methane_symbols() -> list[str]:
    """Methane (CH4) element symbols."""
    return ['C', 'H', 'H', 'H', 'H']


class TestBuildBondGraph:
    """Tests for build_bond_graph function."""

    def test_builds_simple_graph(self, methane_coords: np.ndarray, methane_symbols: list[str]) -> None:
        """Should build correct adjacency list for methane."""
        graph = build_bond_graph(methane_coords, methane_symbols)

        assert len(graph) == 5
        assert 0 in graph
        assert 1 in graph
        assert 2 in graph
        assert 3 in graph
        assert 4 in graph

        assert len(graph[0]) == 4
        assert set(graph[0]) == {1, 2, 3, 4}
        assert len(graph[1]) == 1
        assert graph[1] == [0]

    def test_ignores_overlapping_atoms(self) -> None:
        """Should ignore atoms too close together."""
        coords = np.array([
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0]
        ])
        symbols = ['H', 'C']

        graph = build_bond_graph(coords, symbols, min_dist=0.5)

        assert len(graph) == 2
        assert len(graph[0]) == 0
        assert len(graph[1]) == 0

    def test_raises_for_unknown_element(self) -> None:
        """Should raise ValueError for unknown element."""
        coords = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0]
        ])
        symbols = ['X', 'H']

        with pytest.raises((ValueError, KeyError)):
            _ = build_bond_graph(coords, symbols)

    def test_respects_scale_parameter(self, methane_coords: np.ndarray, methane_symbols: list[str]) -> None:
        """Should adjust bond detection based on scale parameter."""
        graph_default = build_bond_graph(methane_coords, methane_symbols, scale=1.25)
        graph_strict = build_bond_graph(methane_coords, methane_symbols, scale=1.0)

        assert len(graph_default[0]) == 4
        assert len(graph_strict[0]) == 0


class TestGetConnectedComponents:
    """Tests for get_connected_components function."""

    def test_single_component(self, methane_coords: np.ndarray, methane_symbols: list[str]) -> None:
        """Should find single component for methane."""
        graph = build_bond_graph(methane_coords, methane_symbols)
        components = get_connected_components(graph)

        assert len(components) == 1
        assert set(components[0]) == set(range(5))

    def test_two_separated_components(self) -> None:
        """Should find two components for separated molecules."""
        coords = np.array([
            [0.0, 0.0, 0.0],
            [1.09, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [6.09, 0.0, 0.0]
        ])
        symbols = ['C', 'H', 'C', 'H']

        graph = build_bond_graph(coords, symbols)
        components = get_connected_components(graph)

        assert len(components) == 2
        assert 0 in components[0] and 1 in components[0]
        assert 2 in components[1] and 3 in components[1]


class TestFindShortestPath:
    """Tests for find_shortest_path function."""

    def test_same_atom(self, methane_coords: np.ndarray, methane_symbols: list[str]) -> None:
        """Should return single atom for start == end."""
        graph = build_bond_graph(methane_coords, methane_symbols)

        path = find_shortest_path(graph, 0, 0)

        assert path == [0]

    def test_direct_bond(self, methane_coords: np.ndarray, methane_symbols: list[str]) -> None:
        """Should find direct bond between bonded atoms."""
        graph = build_bond_graph(methane_coords, methane_symbols)

        path = find_shortest_path(graph, 0, 1)

        assert path == [0, 1]

    def test_indirect_path(self, methane_coords: np.ndarray, methane_symbols: list[str]) -> None:
        """Should find indirect path through intermediate atoms."""
        graph = build_bond_graph(methane_coords, methane_symbols)

        path = find_shortest_path(graph, 1, 2)

        assert len(path) == 3
        assert path[0] == 1
        assert path[1] == 0
        assert path[2] == 2

    def test_no_path(self) -> None:
        """Should return empty list when no path exists."""
        graph: dict[int, list[int]] = {
            0: [],
            1: [],
            2: [],
        }

        path = find_shortest_path(graph, 0, 2)

        assert path == []


class TestFindCycles:
    """Tests for find_cycles function."""

    def test_no_cycles_in_methane(self, methane_coords: np.ndarray, methane_symbols: list[str]) -> None:
        """Should find no cycles in methane (tree structure)."""
        graph = build_bond_graph(methane_coords, methane_symbols)

        cycles = find_cycles(graph)

        assert len(cycles) == 0

    def test_finds_cyclic_structure(self) -> None:
        """Should find cycles in cyclic molecules."""
        coords = np.array([
            [0.0, 0.0, 0.0],
            [1.4, 0.0, 0.0],
            [2.4, 0.0, 0.0],
            [1.4, 0.0, 1.4]
        ])
        symbols = ['C', 'C', 'C', 'C']

        graph = build_bond_graph(coords, symbols)
        cycles = find_cycles(graph)

        assert len(cycles) >= 1


class TestGetBondDistance:
    """Tests for get_bond_distance function."""

    def test_returns_distance(self, methane_coords: np.ndarray, methane_symbols: list[str]) -> None:
        """Should return correct bond distance."""
        graph = build_bond_graph(methane_coords, methane_symbols)

        distance = get_bond_distance(methane_coords, graph, 0, 1)

        assert abs(distance - 1.09) < 0.01

    def test_raises_for_non_bonded(self, methane_coords: np.ndarray, methane_symbols: list[str]) -> None:
        """Should raise ValueError for non-bonded atoms."""
        graph = build_bond_graph(methane_coords, methane_symbols)

        with pytest.raises(ValueError, match="not bonded"):
            _ = get_bond_distance(methane_coords, graph, 1, 2)
