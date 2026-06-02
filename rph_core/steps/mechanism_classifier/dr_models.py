from __future__ import annotations

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field


class ChiralCenterSet(BaseModel):
    """Map-space stereocenter set used by S0 DR completion."""

    fixed_map_numbers: List[int] = Field(default_factory=list)
    new_map_numbers: List[int] = Field(default_factory=list)
    inferred_from_forming_bonds: bool = False


class ChiralCluster(BaseModel):
    """A coupled stereochemical unit that should be branched together."""

    cluster_id: str
    map_numbers: List[int] = Field(default_factory=list)
    forming_bonds: List[List[int]] = Field(default_factory=list)
    enumeration_mode: str = "concerted_face_flip"


class DRProductBranch(BaseModel):
    """One product-side stereochemical branch derived from the S0 graph."""

    branch_id: str
    pathway_id: str
    role: str
    generation_policy: str
    product_smiles: Optional[str] = None
    flipped_map_numbers: List[int] = Field(default_factory=list)
    fixed_stereocenters: List[int] = Field(default_factory=list)
    notes: List[str] = Field(default_factory=list)


class DRBranchPlan(BaseModel):
    """Serializable S0 branch plan; downstream views must derive from this."""

    version: str = "s0-dr-branch-plan-v1"
    reaction_id: str
    status: str
    reaction_type: str
    cyclo_mode: str
    topology: str
    source: str = "S0_Mechanism/mechanism_graph.json"
    reported_dr: Dict[str, Any] = Field(default_factory=dict)
    stereocenters: ChiralCenterSet = Field(default_factory=ChiralCenterSet)
    chiral_clusters: List[ChiralCluster] = Field(default_factory=list)
    branches: List[DRProductBranch] = Field(default_factory=list)
    warnings: List[str] = Field(default_factory=list)
