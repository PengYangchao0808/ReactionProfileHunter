from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rph_core.utils.task_builder import TaskSpec


@dataclass
class BranchJob:
    parent_reaction_id: str
    branch_id: str
    pathway_id: str
    product_smiles: str
    branch_root: Path
    generation_policy: str = ""
    flipped_map_numbers: list[int] = field(default_factory=list)
    fixed_stereocenters: list[int] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


@dataclass
class ConditionJob:
    reaction_id: str
    condition_id: str
    condition_root: Path
    task: TaskSpec
    temperature_k: float = 298.15


@dataclass
class ReactionJob:
    reaction_id: str
    reaction_root: Path
    representative: TaskSpec
    conditions: list[ConditionJob] = field(default_factory=list)
    dr_plan_path: Path = field(default_factory=Path)
    branches: list[BranchJob] = field(default_factory=list)
