"""
TDD Tests for GEDT Fragment Labeling (Task 8)
==============================================

Tests to verify that:
1. Fragment labeling is deterministic based on hetero atoms
2. O-containing fragment is identified as dipole
3. Fragment charges are output explicitly

Author: RPH Team
Date: 2026-02-03
"""

import pytest
import numpy as np
from pathlib import Path


class TestGEDTFragmentLabeling:
    """Tests for deterministic GEDT fragment labeling."""

    def test_oxyallyl_pattern_labels_o_fragment_as_dipole(self):
        """O-containing fragment should be labeled as dipole (oxyallyl pattern)."""
        from rph_core.utils.fragment_cut import identify_fragment_role

        # Create a simple [5+2] oxyallyl-like structure
        # Fragment A has O (should be dipole), Fragment B is just C
        coordinates = np.array([
            [0.0, 0.0, 0.0],   # 0: C (fragment A - dipole)
            [1.0, 0.0, 0.0],   # 1: O (fragment A - dipole)
            [2.0, 0.0, 0.0],   # 2: C (forming bond)
            [3.0, 0.0, 0.0],   # 3: C (forming bond)
            [4.0, 0.0, 0.0],   # 4: C (fragment B)
        ])
        symbols = ['C', 'O', 'C', 'C', 'C']
        fragment_a = [0, 1, 2]  # Contains O
        fragment_b = [3, 4]     # No O
        forming_bonds = ((2, 3),)

        dipole_label, _, note = identify_fragment_role(
            coordinates, symbols, fragment_a, fragment_b, forming_bonds
        )

        assert dipole_label == 'A', f"O-containing fragment should be dipole, got {dipole_label}"
        assert note == 'hetero_atom_count', f"Should use hetero_atom_count heuristic, got {note}"

    def test_labeling_with_various_heuristics(self):
        """Test that labeling uses various heuristics based on fragment composition."""
        from rph_core.utils.fragment_cut import identify_fragment_role

        # Test case: both fragments have 1 hetero atom (O), different sizes
        coordinates = np.array([
            [0.0, 0.0, 0.0],   # 0: C
            [1.0, 0.0, 0.0],   # 1: O
            [2.0, 0.0, 0.0],   # 2: C
            [3.0, 0.0, 0.0],   # 3: C
            [4.0, 0.0, 0.0],   # 4: O
            [5.0, 0.0, 0.0],   # 5: C
        ])
        symbols = ['C', 'O', 'C', 'C', 'O', 'C']
        # Both fragments have 1 O
        fragment_a = [0, 1, 2]     # 3 atoms, 1 O
        fragment_b = [3, 4, 5]     # 3 atoms, 1 O
        forming_bonds = ((2, 3),)

        dipole_label, dipolarophile_label, note = identify_fragment_role(
            coordinates, symbols, fragment_a, fragment_b, forming_bonds
        )

        # Both have 1 hetero, same size - should use size or fallback
        assert note in ['size_preference', 'forming_bond_first_atom', 'fallback']

    def test_cut_along_forming_bonds_returns_labeling_info(self):
        """cut_along_forming_bonds should return labeling info and fragment charges."""
        from rph_core.utils.fragment_cut import cut_along_forming_bonds

        coordinates = np.array([
            [0.0, 0.0, 0.0],   # 0: C (fragment A)
            [1.0, 0.0, 0.0],   # 1: O (fragment A)
            [2.0, 0.0, 0.0],   # 2: C (forming bond)
            [3.0, 0.0, 0.0],   # 3: C (forming bond)
            [4.0, 0.0, 0.0],   # 4: C (fragment B)
        ])
        symbols = ['C', 'O', 'C', 'C', 'C']
        charges = [0.1, -0.2, 0.0, 0.0, 0.1]  # O is negative, dipole loses electrons
        forming_bonds = ((2, 3),)

        result = cut_along_forming_bonds(
            coordinates=coordinates,
            forming_bonds=forming_bonds,
            charges=charges,
            symbols=symbols
        )

        # Verify result contains labeling info
        assert 'gedt_fragment_labeling' in result
        assert 'q_fragment_dipole' in result
        assert 'q_fragment_dipolarophile' in result

        # Verify fragment charges are calculated
        # Fragment A (dipole): 0.1 + (-0.2) + 0.0 = -0.1
        assert abs(result['q_fragment_dipole'] - (-0.1)) < 0.01
        # Fragment B: 0.0 + 0.1 = 0.1
        assert abs(result['q_fragment_dipolarophile'] - 0.1) < 0.01

        # Verify GEDT is calculated (negative of dipole charge)
        # GEDT = -(-0.1) = 0.1 (positive when electron flows from dipole to dipolarophile)
        assert not np.isnan(result['gedt_value'])

    def test_labeling_deterministic_across_runs(self):
        """Labeling should be deterministic (same result on repeated calls)."""
        from rph_core.utils.fragment_cut import identify_fragment_role

        coordinates = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ])
        symbols = ['C', 'O', 'C', 'C', 'C']
        fragment_a = [0, 1, 2]
        fragment_b = [3, 4]
        forming_bonds = ((2, 3),)

        # Run multiple times
        results = []
        for _ in range(10):
            result = identify_fragment_role(
                coordinates, symbols, fragment_a, fragment_b, forming_bonds
            )
            results.append(result)

        # All results should be identical
        assert all(r == results[0] for r in results), "Labeling should be deterministic"
