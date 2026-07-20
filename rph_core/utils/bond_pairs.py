"""Canonical helpers for unordered atom-index bond pairs."""

from __future__ import annotations

from typing import Iterable, Tuple


BondPair = Tuple[int, int]


def canonicalize_bond_pair(atom_i: int, atom_j: int) -> BondPair:
    """Return one unordered bond pair in deterministic ascending order."""

    i = int(atom_i)
    j = int(atom_j)
    if i == j:
        raise ValueError(f"A bond pair must contain two different atoms: {(i, j)!r}")
    return (i, j) if i < j else (j, i)


def canonicalize_bond_pairs(pairs: Iterable[Iterable[int]]) -> Tuple[BondPair, ...]:
    """Validate, deduplicate and sort an iterable of unordered bond pairs."""

    normalized = set()
    for raw_pair in pairs:
        pair = tuple(raw_pair)
        if len(pair) != 2:
            raise ValueError(f"A bond pair must contain exactly two indices: {pair!r}")
        normalized.add(canonicalize_bond_pair(pair[0], pair[1]))
    return tuple(sorted(normalized))
