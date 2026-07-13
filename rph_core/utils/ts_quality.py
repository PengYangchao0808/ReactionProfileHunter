"""TS quality analysis: mode displacement projection on forming bonds.

P6 of the V4 master plan — upgrades TS validation from "has one imaginary
frequency" to "chemically credible TS" by checking whether the imaginary
mode displacement aligns with forming-bond stretching.
"""

from __future__ import annotations

import logging
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

DEFAULT_PROJECTION_THRESHOLD = 0.3
_FREQUENCY_MATCH_TOLERANCE_CM1 = 5.0


def analyze_ts_quality(
    frequency_output: Optional[Path],
    optimized_xyz: Optional[Path],
    forming_bonds: Optional[Sequence[Tuple[int, int]]],
    frequencies_cm1: Sequence[float],
    imaginary_cutoff_cm1: float = -50.0,
    projection_threshold: float = DEFAULT_PROJECTION_THRESHOLD,
) -> Dict[str, Any]:
    """Analyze TS quality from frequency computation results.

    Returns dict with cleanly separated quality tags:
        frequency_count_valid — exactly one significant imaginary frequency
        mode_displacement_valid — imaginary mode aligns with forming bonds (None if unchecked)
        irc_valid — None (IRC not attempted in V4)
        significant_imaginary_count
        mode_displacement_details
        quality_summary — human-readable label
    """

    normalized_forming_bonds = _normalize_forming_bonds(forming_bonds)
    significant = [f for f in frequencies_cm1 if f <= imaginary_cutoff_cm1]
    frequency_count_valid = len(significant) == 1

    mode_displacement_valid: Optional[bool] = None
    mode_details: Dict[str, Any] = {}

    if (
        frequency_count_valid
        and normalized_forming_bonds
        and optimized_xyz
        and frequency_output
    ):
        try:
            mode_displacement_valid, mode_details = _check_mode_displacement(
                frequency_output,
                optimized_xyz,
                normalized_forming_bonds,
                significant[0],
                projection_threshold,
            )
        except Exception as exc:
            logger.warning("TS mode displacement analysis failed: %s", exc)
            mode_displacement_valid = None
            mode_details = {"error": str(exc)}

    return {
        "frequency_count_valid": frequency_count_valid,
        "mode_displacement_valid": mode_displacement_valid,
        "irc_valid": None,
        "significant_imaginary_count": len(significant),
        "mode_displacement_details": mode_details,
        "quality_summary": _quality_summary(
            frequency_count_valid, mode_displacement_valid
        ),
    }


def _check_mode_displacement(
    frequency_output: Path,
    optimized_xyz: Path,
    forming_bonds: Sequence[Tuple[int, int]],
    imaginary_frequency: float,
    threshold: float,
) -> Tuple[Optional[bool], Dict[str, Any]]:
    """Project the imaginary-mode displacement onto each forming bond axis.

    For bond (a, b) the *relative* displacement (d_b - d_a) projected onto the
    bond unit vector reveals whether the mode stretches/compresses that bond.
    A significant |projection| means the mode is "along" the reaction
    coordinate defined by the forming bonds.
    """

    geometry = _read_xyz_coordinates(optimized_xyz)
    displacements = _parse_imaginary_mode_displacements(
        frequency_output, imaginary_frequency
    )

    if displacements is None:
        return None, {"note": "Displacement parsing not available for this log format."}

    if len(displacements) != len(geometry):
        return None, {
            "error": f"Displacement count {len(displacements)} != geometry atoms {len(geometry)}",
        }

    bond_reports: List[Dict[str, Any]] = []
    all_aligned = True

    for a, b in forming_bonds:
        if a < 0 or b < 0 or a >= len(geometry) or b >= len(geometry):
            bond_reports.append({"bond": [a, b], "error": "Atom index out of range"})
            all_aligned = False
            continue

        bond_vec = [geometry[b][i] - geometry[a][i] for i in range(3)]
        norm = _vec_norm(bond_vec)
        if norm < 1e-10:
            bond_reports.append({"bond": [a, b], "error": "Degenerate bond"})
            all_aligned = False
            continue

        unit = [c / norm for c in bond_vec]
        rel = [displacements[b][i] - displacements[a][i] for i in range(3)]
        projection = sum(rel[i] * unit[i] for i in range(3))
        aligned = abs(projection) >= threshold

        bond_reports.append(
            {
                "bond": [a, b],
                "projection": round(projection, 4),
                "abs_projection": round(abs(projection), 4),
                "threshold": threshold,
                "aligned": aligned,
            }
        )
        if not aligned:
            all_aligned = False

    return all_aligned, {
        "bonds": bond_reports,
        "imaginary_frequency_cm1": imaginary_frequency,
    }


def _read_xyz_coordinates(xyz_path: Path) -> List[List[float]]:
    """Read all atom coordinates from an XYZ file."""

    lines = xyz_path.read_text(encoding="utf-8").strip().split("\n")
    if len(lines) < 2:
        raise ValueError(f"XYZ file {xyz_path} is too short")
    n = int(lines[0].strip())
    coords: List[List[float]] = []
    for line in lines[2 : 2 + n]:
        parts = line.strip().split()
        if len(parts) >= 4:
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return coords


def _parse_imaginary_mode_displacements(
    frequency_output: Path,
    imaginary_frequency: float,
) -> Optional[List[List[float]]]:
    """Best-effort parse of the imaginary mode displacement vectors.

    Gaussian and ORCA frequency logs have very different layouts.  This
    parser tries the Gaussian harmonic-displacement block first, then the
    ORCA normal-mode column block.  Returns ``None`` when the format is
    not recognised — callers treat ``None`` as "not checked" (not "failed").
    """

    text = frequency_output.read_text(encoding="utf-8", errors="replace")

    gaussian = _try_gaussian_displacements(text, imaginary_frequency)
    if gaussian is not None:
        return gaussian

    orca = _try_orca_displacements(text, imaginary_frequency)
    if orca is not None:
        return orca

    logger.warning(
        "Could not parse mode displacements from %s; skipping projection check.",
        frequency_output,
    )
    return None


def _try_gaussian_displacements(
    text: str, target_freq: float
) -> Optional[List[List[float]]]:
    """Parse Gaussian-format harmonic mode displacements.

    Gaussian prints each frequency in a column block.  After the header line
    containing the frequency value, there are rows like::

        Atom  1   X   0.12  Y   0.05  Z  -0.03

    We find the column whose header matches *target_freq* and collect the
    ``X Y Z`` triplets per atom.
    """

    lines = text.split("\n")

    index = 0
    while index < len(lines):
        freq_match = re.match(r"\s*Frequencies\s+--\s+(.*)", lines[index])
        if not freq_match:
            index += 1
            continue

        try:
            frequencies = [float(value) for value in freq_match.group(1).split()]
        except ValueError:
            index += 1
            continue
        if not frequencies:
            index += 1
            continue

        target_column = _matching_frequency_column(frequencies, target_freq)
        index += 1
        simple_displacements = _try_simple_gaussian_displacements(lines, index)
        if target_column == 0 and simple_displacements is not None:
            return simple_displacements
        while index < len(lines) and not re.match(r"\s*Atom\s+AN\b", lines[index]):
            if re.match(r"\s*Frequencies\s+--\s+", lines[index]):
                break
            index += 1
        if target_column is None or index >= len(lines) or not re.match(r"\s*Atom\s+AN\b", lines[index]):
            continue

        index += 1
        atom_rows: List[List[float]] = []
        min_expected_fields = 2 + 3 * len(frequencies)
        while index < len(lines):
            parts = lines[index].split()
            if len(parts) < min_expected_fields:
                break
            try:
                int(parts[0])
                int(parts[1])
            except ValueError:
                break
            try:
                atom_rows.append([float(value) for value in parts[2 : 2 + 3 * len(frequencies)]])
            except ValueError:
                break
            index += 1
        if not atom_rows:
            continue

        column_start = target_column * 3
        return [row[column_start : column_start + 3] for row in atom_rows]

    return None


def _try_orca_displacements(
    text: str, target_freq: float
) -> Optional[List[List[float]]]:
    """Parse ORCA normal-mode displacement block (best-effort).

    ORCA prints vibrational modes in a column matrix.  This is a simplified
    parser that returns None when the layout is not the expected one.
    """

    import re

    section_match = re.search(
        r"CARTESIAN DISPLACEMENTS\s*\n\s*-+\s*\n",
        text,
    )
    if not section_match:
        return None

    section = text[section_match.end():]
    for mode_block in re.split(r"\n\s*(?=Mode:\s*\d+\s*\n)", section):
        freq_match = re.search(r"Freq:\s*([\d\.\-+Ee]+)", mode_block)
        if not freq_match:
            continue
        try:
            frequency = float(freq_match.group(1))
        except ValueError:
            continue
        if not _frequency_matches(frequency, target_freq):
            continue

        displacements: List[List[float]] = []
        for line in mode_block.split("\n"):
            atom_match = re.match(r"\s*Atom\s+\d+\s*:\s+(.*)", line)
            if not atom_match:
                continue
            parts = atom_match.group(1).split()
            if len(parts) < 3:
                continue
            try:
                displacements.append([float(parts[0]), float(parts[1]), float(parts[2])])
            except ValueError:
                continue
        if displacements:
            return displacements

    return None


def _try_simple_gaussian_displacements(
    lines: Sequence[str],
    start_index: int,
) -> Optional[List[List[float]]]:
    displacements: List[List[float]] = []
    for index in range(start_index, min(start_index + 200, len(lines))):
        line = lines[index].strip()
        if not line:
            if displacements:
                break
            continue
        if re.match(r"\s*Frequencies\s+--\s+", line):
            break
        match = re.match(
            r"Atom\s+\d+\s+X\s+([-\d.]+)\s+Y\s+([-\d.]+)\s+Z\s+([-\d.]+)",
            line,
        )
        if match:
            displacements.append(
                [float(match.group(1)), float(match.group(2)), float(match.group(3))]
            )
    return displacements if displacements else None


def _normalize_forming_bonds(
    forming_bonds: Optional[Sequence[Tuple[int, int]]],
) -> List[Tuple[int, int]]:
    normalized: List[Tuple[int, int]] = []
    for pair in forming_bonds or ():
        if len(pair) != 2:
            continue
        normalized.append((int(pair[0]), int(pair[1])))
    return normalized


def _matching_frequency_column(
    frequencies: Sequence[float],
    target_freq: float,
) -> Optional[int]:
    candidates = [
        (index, abs(value - target_freq))
        for index, value in enumerate(frequencies)
        if _frequency_matches(value, target_freq)
    ]
    if not candidates:
        return None
    return min(candidates, key=lambda item: item[1])[0]


def _frequency_matches(value: float, target_freq: float) -> bool:
    return abs(value - target_freq) <= _FREQUENCY_MATCH_TOLERANCE_CM1


def _vec_norm(v: Sequence[float]) -> float:
    return math.sqrt(sum(c * c for c in v))


def _quality_summary(
    freq_valid: bool, mode_valid: Optional[bool]
) -> str:
    if not freq_valid:
        return "frequency_count_failed"
    if mode_valid is None:
        return "frequency_valid_mode_not_checked"
    if not mode_valid:
        return "frequency_valid_mode_misaligned"
    return "frequency_valid_mode_aligned"
