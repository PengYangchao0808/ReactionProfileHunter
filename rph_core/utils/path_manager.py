from __future__ import annotations

from pathlib import Path


def get_reaction_root(output_root: Path, reaction_id: str) -> Path:
    return Path(output_root) / str(reaction_id)


def get_reaction_features_dir(reaction_root: Path) -> Path:
    return Path(reaction_root) / "reaction_features"


def get_conditions_root(reaction_root: Path) -> Path:
    return Path(reaction_root) / "conditions"


def get_condition_root(output_root: Path, reaction_id: str, condition_id: str) -> Path:
    return get_conditions_root(get_reaction_root(output_root, reaction_id)) / str(condition_id)
