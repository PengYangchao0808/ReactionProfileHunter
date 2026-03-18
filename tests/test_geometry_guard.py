"""
Tests for Step2 Geometry Guard - Topology and Close-Contact Validation
"""

import numpy as np
import pytest
from pathlib import Path

from rph_core.steps.step2_retro.geometry_guard import (
    compare_graph_topology,
    detect_risky_contacts,
    generate_keepaway_constraints,
    TopologyGuardResult,
    RiskyContactResult,
)


class TestCompareGraphTopology:
    """Test graph topology comparison."""

    def test_valid_no_changes(self):
        """Test valid case: no topology changes expected."""
        coords = np.array([
            [0.0, 0.0, 0.0],  # C
            [1.5, 0.0, 0.0],  # C - bonded to 0
            [3.0, 0.0, 0.0],  # C - bonded to 1
            [4.5, 0.0, 0.0],  # C - bonded to 2
        ])
        symbols = ["C", "C", "C", "C"]
        forming_bonds = [(1, 2)]  # Bond to break

        result = compare_graph_topology(
            product_coords=coords,
            candidate_coords=coords.copy(),
            symbols=symbols,
            forming_bonds=forming_bonds,
        )

        assert result.is_valid is True
        assert len(result.new_edges) == 0
        assert len(result.lost_edges) == 0

    def test_new_edge_detected(self):
        """Test detection of unintended new bond."""
        product_coords = np.array([
            [0.0, 0.0, 0.0],  # C
            [1.5, 0.0, 0.0],  # C
            [4.0, 0.0, 0.0],  # C - far from 0 and 1
        ])
        # Candidate has atoms 0 and 2 too close (new bond)
        candidate_coords = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [1.4, 0.0, 0.0],  # Now close to atom 0
        ])
        symbols = ["C", "C", "C"]
        forming_bonds = [(0, 1)]

        result = compare_graph_topology(
            product_coords=product_coords,
            candidate_coords=candidate_coords,
            symbols=symbols,
            forming_bonds=forming_bonds,
        )

        assert result.is_valid is False
        assert len(result.new_edges) > 0
        # Check that the new edge is (0, 2)
        edge_atoms = [(e[0], e[1]) for e in result.new_edges]
        assert (0, 2) in edge_atoms

    def test_forming_bond_deletion_allowed(self):
        """Test that forming bond deletion is allowed."""
        product_coords = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
        ])
        # Candidate has bond stretched
        candidate_coords = np.array([
            [0.0, 0.0, 0.0],
            [3.5, 0.0, 0.0],  # Stretched forming bond
        ])
        symbols = ["C", "C"]
        forming_bonds = [(0, 1)]

        result = compare_graph_topology(
            product_coords=product_coords,
            candidate_coords=candidate_coords,
            symbols=symbols,
            forming_bonds=forming_bonds,
        )

        assert result.is_valid is True
        # The forming bond should be "lost" but that's allowed

    def test_lost_nonforming_edge_rejected(self):
        """Test that loss of non-forming bonds is rejected."""
        product_coords = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ])
        # Candidate has bond 0-1 broken unintentionally
        candidate_coords = np.array([
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],  # Far from 0, but forming bond is (1, 2)
            [4.5, 0.0, 0.0],
        ])
        symbols = ["C", "C", "C"]
        forming_bonds = [(1, 2)]  # Not (0, 1)

        result = compare_graph_topology(
            product_coords=product_coords,
            candidate_coords=candidate_coords,
            symbols=symbols,
            forming_bonds=forming_bonds,
        )

        assert result.is_valid is False
        assert len(result.lost_edges) > 0


class TestDetectRiskyContacts:
    """Test risky contact detection."""

    def test_no_risky_contacts(self):
        """Test case with no risky contacts."""
        coords = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ])
        symbols = ["C", "C", "C"]

        result = detect_risky_contacts(
            product_coords=coords,
            candidate_coords=coords.copy(),
            symbols=symbols,
        )

        assert len(result.risky_pairs) == 0

    def test_risky_contact_detected(self):
        """Test detection of pair that shrank into near-bond window."""
        product_coords = np.array([
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],  # Far from atom 0
            [5.5, 0.0, 0.0],  # Bonded to 1 in product
        ])
        # Candidate: atom 2 moved close to atom 0, but keep atom 1 at safe distance
        candidate_coords = np.array([
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],   # Keep far from atom 0
            [1.2, 0.0, 0.0],   # Now very close to atom 0 (within near-bond window)
        ])
        symbols = ["C", "C", "C"]

        result = detect_risky_contacts(
            product_coords=product_coords,
            candidate_coords=candidate_coords,
            symbols=symbols,
        )

        # Should detect (0, 2) as risky
        assert len(result.risky_pairs) > 0
        pair_atoms = [(p[0], p[1]) for p in result.risky_pairs]
        assert (0, 2) in pair_atoms

    def test_max_pairs_limit(self):
        """Test that max_pairs limits the result."""
        # Create geometry with many atoms
        n_atoms = 10
        product_coords = np.random.rand(n_atoms, 3) * 5.0
        # Candidate has many atoms moved closer
        candidate_coords = product_coords.copy()
        candidate_coords[1:, 0] += 0.5  # Move all atoms closer to atom 0

        symbols = ["C"] * n_atoms

        result = detect_risky_contacts(
            product_coords=product_coords,
            candidate_coords=candidate_coords,
            symbols=symbols,
            max_pairs=3,
        )

        assert len(result.risky_pairs) <= 3


class TestGenerateKeepawayConstraints:
    """Test keep-away constraint generation."""

    def test_basic_generation(self):
        """Test basic constraint generation."""
        risky_pairs = [
            (0, 2, 1.4, 3.5),  # (i, j, cand_dist, prod_dist)
            (1, 3, 1.5, 3.0),
        ]

        result = generate_keepaway_constraints(
            risky_pairs=risky_pairs,
            keep_apart_floor=3.0,
            force_constant=0.5
        )

        constraints = result["distance_constraints"]
        assert len(constraints) == 2
        # Target should be max(prod_dist, floor)
        assert constraints["0 2"] == 3.5  # max(3.5, 3.0)
        assert constraints["1 3"] == 3.0  # max(3.0, 3.0)
        assert result["force_constant"] == 0.5

    def test_floor_applied(self):
        """Test that floor is applied when product distance is smaller."""
        risky_pairs = [
            (0, 2, 1.4, 2.0),  # prod_dist is only 2.0, floor is 3.0
        ]

        result = generate_keepaway_constraints(
            risky_pairs=risky_pairs,
            keep_apart_floor=3.0,
            force_constant=0.5
        )

        constraints = result["distance_constraints"]
        assert constraints["0 2"] == 3.0  # floor applied


class TestIntegration:
    """Integration tests combining multiple functions."""

    def test_seed_to_xtb_guard_flow(self):
        """Test the full flow from seed assessment to constraint generation."""
        # Product geometry: C-C-C chain
        product_coords = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.5, 0.0, 0.0],
        ])
        symbols = ["C", "C", "C", "C"]
        forming_bonds = [(1, 2)]

        # Seed: forming bond stretched, but atoms 0 and 3 accidentally close
        seed_coords = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [4.5, 0.0, 0.0],  # Stretched
            [1.4, 0.0, 0.0],  # Close to atom 0
        ])

        # Step 1: Check seed topology
        seed_topology = compare_graph_topology(
            product_coords=product_coords,
            candidate_coords=seed_coords,
            symbols=symbols,
            forming_bonds=forming_bonds,
        )

        assert not seed_topology.is_valid  # New edge (0, 3)

        # Step 2: Detect risky contacts
        seed_risky = detect_risky_contacts(
            product_coords=product_coords,
            candidate_coords=seed_coords,
            symbols=symbols,
        )

        assert len(seed_risky.risky_pairs) > 0

        # Step 3: Generate constraints
        keepaway_cfg = generate_keepaway_constraints(
            seed_risky.risky_pairs,
            keep_apart_floor=3.0,
            force_constant=0.5
        )

        assert len(keepaway_cfg["distance_constraints"]) > 0
