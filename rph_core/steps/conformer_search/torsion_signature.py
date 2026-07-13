"""Deterministic, stereochemistry-aware torsion signatures for CENSO-LITE."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import Lipinski


BondPair = Tuple[int, int]


def circular_distance_deg(left: float, right: float) -> float:
    """Return the shortest distance between two signed dihedral angles."""
    return abs((float(left) - float(right) + 180.0) % 360.0 - 180.0)


def _dihedral(
    p0: Sequence[float],
    p1: Sequence[float],
    p2: Sequence[float],
    p3: Sequence[float],
) -> float:
    p0_array = np.asarray(p0, dtype=float)
    p1_array = np.asarray(p1, dtype=float)
    p2_array = np.asarray(p2, dtype=float)
    p3_array = np.asarray(p3, dtype=float)
    b0 = -(p1_array - p0_array)
    b1 = p2_array - p1_array
    b2 = p3_array - p2_array
    b1 /= np.linalg.norm(b1) or 1.0
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def rotatable_bonds(mol: Chem.Mol) -> Tuple[BondPair, ...]:
    """Return canonical heavy-atom rotatable bonds from the molecular graph."""
    bonds = []
    for match in mol.GetSubstructMatches(Lipinski.RotatableBondSmarts):
        i, j = sorted((int(match[0]), int(match[1])))
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond is None or bond.IsInRing():
            continue
        bonds.append((i, j))
    return tuple(sorted(set(bonds)))


def _outer_atom(mol: Chem.Mol, center: int, other: int) -> int | None:
    candidates = [atom.GetIdx() for atom in mol.GetAtomWithIdx(center).GetNeighbors() if atom.GetIdx() != other]
    heavy = [idx for idx in candidates if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1]
    return min(heavy or candidates, default=None)


def _torsion_entries(
    mol: Chem.Mol,
    coordinates: Sequence[Sequence[float]],
) -> List[Tuple[BondPair, float]]:
    coords = np.asarray(coordinates, dtype=float)
    entries: List[Tuple[BondPair, float]] = []
    for i, j in rotatable_bonds(mol):
        a = _outer_atom(mol, i, j)
        d = _outer_atom(mol, j, i)
        if a is None or d is None:
            continue
        entries.append(((i, j), _dihedral(coords[a], coords[i], coords[j], coords[d])))
    return entries


def torsion_values(mol: Chem.Mol, coordinates: Sequence[Sequence[float]]) -> Tuple[float, ...]:
    return tuple(value for _bond, value in _torsion_entries(mol, coordinates))


@dataclass(frozen=True)
class TorsionSignature:
    bonds: Tuple[BondPair, ...]
    bins: Tuple[int, ...]
    values: Tuple[float, ...]
    bond_space: str = "atom_index"

    def key(self) -> str:
        parts = [
            f"{i}-{j}:{bucket}"
            for (i, j), bucket in sorted(zip(self.bonds, self.bins), key=lambda item: item[0])
        ]
        return "|".join(parts) or "no_rotatable_bonds"


def build_signature(
    mol: Chem.Mol,
    coordinates: Sequence[Sequence[float]],
    bin_width_deg: float = 20.0,
) -> TorsionSignature:
    if bin_width_deg <= 0.0 or bin_width_deg > 180.0:
        raise ValueError("bin_width_deg must be in the interval (0, 180]")
    idx_to_mapnum = {
        atom.GetIdx(): atom.GetAtomMapNum()
        for atom in mol.GetAtoms()
        if atom.GetAtomMapNum() > 0
    }
    has_maps = bool(idx_to_mapnum)
    entries = []
    for (i, j), value in _torsion_entries(mol, coordinates):
        bond_label: BondPair = (i, j)
        if has_maps and i in idx_to_mapnum and j in idx_to_mapnum:
            bond_label = tuple(sorted((idx_to_mapnum[i], idx_to_mapnum[j])))
        entries.append((
            bond_label,
            int(math.floor((value + 180.0) / bin_width_deg)),
            value,
        ))
    entries.sort(key=lambda item: item[0])
    return TorsionSignature(
        bonds=tuple(item[0] for item in entries),
        bins=tuple(item[1] for item in entries),
        values=tuple(item[2] for item in entries),
        bond_space="atom_map" if has_maps else "atom_index",
    )


def signatures_equivalent(left: TorsionSignature, right: TorsionSignature, tolerance_deg: float = 25.0) -> bool:
    if left.bond_space != right.bond_space:
        return False
    if left.bonds != right.bonds or len(left.values) != len(right.values):
        return False
    return all(circular_distance_deg(a, b) <= tolerance_deg for a, b in zip(left.values, right.values))
