from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from rph_core.steps.refinement.engine import RefinementEngine
from rph_core.utils.orca_interface import parse_orca6_normal_mode_vectors


_ORCA6_NORMAL_MODES_HEADER = """NORMAL MODES
------------

These modes are the Cartesian displacements weighted by the diagonal matrix
M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
Thus, these vectors are normalized but *not* orthogonal

"""


def _orca6_body(mode_indices: list[int], n_coords: int, fill: float = 0.25) -> str:
    lines = []
    for start in range(0, len(mode_indices), 6):
        chunk = mode_indices[start:start + 6]
        lines.append(" " + " ".join(f"{i:>10}" for i in chunk))
        for row in range(n_coords):
            lines.append(f"{row:>7}" + "".join(f"{fill:>11.6f}" for _ in chunk))
    return "\n".join(lines) + "\n"


def test_parse_orca6_normal_modes_shape_and_normalization() -> None:
    content = _ORCA6_NORMAL_MODES_HEADER + _orca6_body([0, 1, 2, 3, 4, 5, 6], 6)
    symbols = ["H", "H"]  # 2 atoms -> 6 Cartesian coordinates

    vectors = parse_orca6_normal_mode_vectors(content, symbols)

    assert set(vectors) == {0, 1, 2, 3, 4, 5, 6}
    for mode, vector in vectors.items():
        assert vector.shape == (2, 3)
        assert np.linalg.norm(vector) == pytest.approx(1.0)


def test_parse_orca6_normal_modes_demasses_by_sqrt_mass() -> None:
    # Mass-weighted storage: heavy atoms (C) get raw values scaled by 1/sqrt(m).
    # After deweighting, per-atom Cartesian displacement is uniform.
    content = _ORCA6_NORMAL_MODES_HEADER + _orca6_body([0], 6, fill=0.25)
    symbols = ["C", "H"]

    vectors = parse_orca6_normal_mode_vectors(content, symbols)
    vector = vectors[0]

    # Raw mass-weighted rows were equal (0.25).  Deweighted C/H rows differ by
    # sqrt(m_C)/sqrt(m_H) before the shared norm rescale; the ratio is preserved.
    ratio = np.linalg.norm(vector[0]) / np.linalg.norm(vector[1])
    expected = np.sqrt(12.011 / 1.008)
    assert ratio == pytest.approx(expected, rel=1e-3)


def test_parse_orca6_normal_modes_rejects_mismatched_rows() -> None:
    content = _ORCA6_NORMAL_MODES_HEADER + _orca6_body([0], 6, fill=0.25)
    symbols = ["H", "H", "H"]  # 3 atoms -> needs 9 rows, only 6 present

    vectors = parse_orca6_normal_mode_vectors(content, symbols)

    assert vectors == {}


def test_engine_demass_weight_mode_vectors_normalizes() -> None:
    displacements = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    symbols = ["C", "H"]

    deweighted = RefinementEngine._demass_weight_mode_vectors(displacements, symbols)

    assert deweighted is not None
    assert np.linalg.norm(deweighted) == pytest.approx(1.0)
    assert deweighted.shape == (2, 3)


def test_engine_demass_weight_mode_vectors_falls_back_mass_for_unknown() -> None:
    displacements = np.asarray([[1.0, 0.0, 0.0]])
    symbols = ["Xx"]

    deweighted = RefinementEngine._demass_weight_mode_vectors(displacements, symbols)

    assert deweighted is not None
    assert deweighted.shape == (1, 3)


def test_engine_demass_weight_mode_vectors_requires_shape_match() -> None:
    displacements = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    symbols = ["C"]

    deweighted = RefinementEngine._demass_weight_mode_vectors(displacements, symbols)

    assert deweighted is None

