"""
Fragment Manipulation Utilities
==============================

H-cap and fragment manipulation operations for fragmenter.

Author: QCcalc Team
Date: 2026-01-27
"""

from typing import List, Tuple, Dict, Optional
import logging

import numpy as np

logger = logging.getLogger(__name__)


def h_cap_fragment(
    coords: np.ndarray,
    symbols: List[str],
    cap_positions: List[Tuple[int, np.ndarray]],
    cap_bond_lengths: Optional[Dict[str, float]] = None
) -> Tuple[np.ndarray, List[str]]:
    """
    Add hydrogen atoms to cap open valences.

    Args:
        coords: (N, 3) coordinates
        symbols: List of element symbols
        cap_positions: List of (atom_idx, cap_direction_vec)
            - atom_idx: Index of atom to cap
            - cap_direction_vec: Vector pointing from atom toward cap position
              (should be colinear with bond being cut, pointing away from cut bond)
        cap_bond_lengths: Bond length by element (Å)
            Default: {'C': 1.09, 'N': 1.01, 'O': 0.96, 'H': 0.76}

    Returns:
        (capped_coords, capped_symbols) - Extended coordinates and symbols
    """
    if cap_bond_lengths is None:
        cap_bond_lengths = {
            'C': 1.09, 'N': 1.01, 'O': 0.96, 'H': 0.76,
            'Si': 1.48, 'P': 1.42, 'S': 1.35, 'Cl': 1.27,
            'Br': 1.42, 'I': 1.54
        }

    capped_coords = coords.copy()
    capped_symbols = symbols.copy()

    for atom_idx, direction_vec in cap_positions:
        direction = direction_vec / np.linalg.norm(direction_vec)

        element = symbols[atom_idx]
        bond_length = cap_bond_lengths.get(element, 1.09)

        cap_pos = coords[atom_idx] - direction * bond_length

        capped_coords = np.vstack([capped_coords, cap_pos])
        capped_symbols.append('H')

    return capped_coords, capped_symbols


def get_fragment_charges(
    total_charge: int,
    n_fragA: int,
    n_fragB: int,
    dipole_in_fragA: bool = True
) -> Tuple[int, int]:
    """
    Distribute total charge between two fragments.

    For oxidopyrylium [5+2], assign +1 to dipole (fragment A).

    Args:
        total_charge: Total charge of the system
        n_fragA: Number of atoms in fragment A
        n_fragB: Number of atoms in fragment B
        dipole_in_fragA: If True, assign formal +1 charge to fragment A

    Returns:
        (charge_fragA, charge_fragB) tuple
    """
    if dipole_in_fragA:
        return total_charge + 1, 0
    else:
        return 0, total_charge + 1


def get_fragment_multiplicities(
    total_multiplicity: int,
    n_fragA: int,
    n_fragB: int,
    dipole_in_fragA: bool = True
) -> Tuple[int, int]:
    """
    Distribute total multiplicity between two fragments.

    Default: assign multiplicity to fragment containing dipole.

    Args:
        total_multiplicity: Total multiplicity of the system
        n_fragA: Number of atoms in fragment A
        n_fragB: Number of atoms in fragment B
        dipole_in_fragA: If True, assign multiplicity to fragment A

    Returns:
        (mult_fragA, mult_fragB) tuple
    """
    if dipole_in_fragA:
        return total_multiplicity, 1
    else:
        return 1, total_multiplicity
