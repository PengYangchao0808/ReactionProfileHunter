"""
Tests for fragment manipulation module
======================================

Unit tests for fragment manipulation utilities.

Author: QCcalc Team
Date: 2026-01-27
"""

import pytest
from typing import List
import numpy as np

from rph_core.utils.fragment_manipulation import (
    h_cap_fragment,
    get_fragment_charges,
    get_fragment_multiplicities
)


@pytest.fixture
def simple_coords():
    """Simple 2-atom coordinates."""
    return np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0]
    ])


@pytest.fixture
def simple_symbols():
    """Simple element symbols."""
    return ['C', 'C']


class TestHCapFragment:
    """Tests for h_cap_fragment function."""

    def test_adds_single_cap(self, simple_coords, simple_symbols):
        """Should add one H atom."""
        cap_positions = [(0, np.array([1.0, 0.0, 0.0]))]

        capped_coords, capped_symbols = h_cap_fragment(
            simple_coords, simple_symbols, cap_positions
        )

        assert len(capped_coords) == 3
        assert len(capped_symbols) == 3
        assert capped_symbols == ['C', 'C', 'H']
        assert capped_coords[2, 0] > 1.0

    def test_adds_two_caps(self, simple_coords, simple_symbols):
        """Should add two H atoms."""
        cap_positions = [
            (0, np.array([1.0, 0.0, 0.0])),
            (1, np.array([-1.0, 0.0, 0.0]))
        ]

        capped_coords, capped_symbols = h_cap_fragment(
            simple_coords, simple_symbols, cap_positions
        )

        assert len(capped_coords) == 4
        assert capped_symbols[-2:] == ['H', 'H']

    def test_respects_custom_bond_lengths(self, simple_coords, simple_symbols):
        """Should use custom bond lengths when provided."""
        cap_positions = [(0, np.array([1.0, 0.0, 0.0]))]

        custom_lengths = {'C': 1.0}
        capped_coords, _ = h_cap_fragment(
            simple_coords, simple_symbols, cap_positions,
            cap_bond_lengths=custom_lengths
        )

        distance = np.linalg.norm(capped_coords[2] - capped_coords[0])
        assert abs(distance - 1.0) < 0.01


class TestGetFragmentCharges:
    """Tests for get_fragment_charges function."""

    def test_assigns_positive_charge_to_dipole(self):
        """Should assign +1 to fragment A (dipole)."""
        chargeA, chargeB = get_fragment_charges(
            total_charge=0, n_fragA=10, n_fragB=8, dipole_in_fragA=True
        )

        assert chargeA == 1
        assert chargeB == 0

    def test_assigns_positive_charge_to_dipole_with_total(self):
        """Should handle total charge > 0."""
        chargeA, chargeB = get_fragment_charges(
            total_charge=2, n_fragA=10, n_fragB=8, dipole_in_fragA=True
        )

        assert chargeA == 3
        assert chargeB == 0

    def test_assigns_charge_to_other_fragment(self):
        """Should assign charge to fragment B when dipole is there."""
        chargeA, chargeB = get_fragment_charges(
            total_charge=0, n_fragA=10, n_fragB=8, dipole_in_fragA=False
        )

        assert chargeA == 0
        assert chargeB == 1


class TestGetFragmentMultiplicities:
    """Tests for get_fragment_multiplicities function."""

    def test_assigns_multiplicity_to_dipole(self):
        """Should assign total multiplicity to fragment A (dipole)."""
        multA, multB = get_fragment_multiplicities(
            total_multiplicity=3, n_fragA=10, n_fragB=8, dipole_in_fragA=True
        )

        assert multA == 3
        assert multB == 1

    def test_assigns_multiplicity_to_other_fragment(self):
        """Should assign multiplicity to fragment B when dipole is there."""
        multA, multB = get_fragment_multiplicities(
            total_multiplicity=3, n_fragA=10, n_fragB=8, dipole_in_fragA=False
        )

        assert multA == 1
        assert multB == 3
