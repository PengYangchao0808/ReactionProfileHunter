from __future__ import annotations

import json
import pytest
from pathlib import Path

from rph_core.utils.path_manager import (
    get_global_small_molecules_root,
    get_precursor_root,
    get_precursor_s1_dir,
    get_condition_branch_root,
    get_reaction_root,
    get_branch_root,
    get_conditions_root,
    get_branches_root,
    get_reaction_features_dir,
    get_condition_root,
)


class TestNewPathHelpers:
    def test_get_global_small_molecules_root(self):
        assert get_global_small_molecules_root(Path("/out")) == Path("/out/small_molecules")

    def test_get_precursor_root(self):
        assert get_precursor_root(Path("/out/RXN_001")) == Path("/out/RXN_001/precursor")

    def test_get_precursor_s1_dir(self):
        assert get_precursor_s1_dir(Path("/out/RXN_001")) == Path("/out/RXN_001/precursor/S1_ConfGeneration")

    def test_get_condition_branch_root(self):
        assert get_condition_branch_root(Path("/cond"), "BR_DR_001") == Path("/cond/branches/BR_DR_001")


class TestExistingPathHelpers:
    def test_get_reaction_root(self):
        assert get_reaction_root(Path("/out"), "RXN_001") == Path("/out/RXN_001")

    def test_get_branch_root(self):
        assert get_branch_root(Path("/out/RXN"), "BR_DR_001") == Path("/out/RXN/branches/BR_DR_001")

    def test_get_conditions_root(self):
        assert get_conditions_root(Path("/out/RXN")) == Path("/out/RXN/conditions")

    def test_get_condition_root(self):
        assert get_condition_root(Path("/out"), "RXN", "COND_001") == Path("/out/RXN/conditions/COND_001")
