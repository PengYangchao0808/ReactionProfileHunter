import numpy as np
import pytest

from rph_core.steps.step4_features.context import FeatureContext
from rph_core.utils.fragment_cut import FragmentCutter, cut_along_forming_bonds


def test_feature_context_defaults_to_mol_idx_domain() -> None:
    context = FeatureContext()
    assert context.index_domain == "mol_idx"


def test_fragment_cutter_rejects_out_of_range_forming_bonds() -> None:
    cutter = FragmentCutter(config={})
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
            [2.4, 0.0, 0.0],
            [3.6, 0.0, 0.0],
        ]
    )

    with pytest.raises(ValueError, match="out of range"):
        cutter.cut_molecule(coordinates=coords, forming_bonds=((0, 8), (1, 2)))


def test_cut_along_forming_bonds_rejects_negative_indices() -> None:
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
            [2.4, 0.0, 0.0],
            [3.6, 0.0, 0.0],
        ]
    )

    with pytest.raises(ValueError, match="negative index"):
        cut_along_forming_bonds(
            coordinates=coords,
            forming_bonds=((-1, 2), (1, 3)),
            charges=[0.0, 0.0, 0.0, 0.0],
            symbols=["C", "C", "C", "C"],
        )
