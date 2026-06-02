from __future__ import annotations

from pathlib import Path


def get_reaction_root(output_root: Path, reaction_id: str) -> Path:
    return Path(output_root) / str(reaction_id)


def get_reaction_features_dir(reaction_root: Path) -> Path:
    return Path(reaction_root) / "reaction_features"


def get_conditions_root(reaction_root: Path) -> Path:
    return Path(reaction_root) / "conditions"


def get_branches_root(reaction_root: Path) -> Path:
    return Path(reaction_root) / "branches"


def get_branch_root(reaction_root: Path, branch_id: str) -> Path:
    return get_branches_root(reaction_root) / str(branch_id)


def get_condition_root(output_root: Path, reaction_id: str, condition_id: str) -> Path:
    return get_conditions_root(get_reaction_root(output_root, reaction_id)) / str(condition_id)


def get_global_small_molecules_root(output_root: Path) -> Path:
    return Path(output_root) / "small_molecules"


def get_precursor_root(reaction_root: Path) -> Path:
    return Path(reaction_root) / "precursor"


def get_precursor_s1_dir(reaction_root: Path) -> Path:
    return get_precursor_root(reaction_root) / "S1_ConfGeneration"


def get_condition_branch_root(condition_root: Path, branch_id: str) -> Path:
    return Path(condition_root) / "branches" / str(branch_id)
