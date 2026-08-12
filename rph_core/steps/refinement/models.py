"""Data models for the unified refinement engine."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


@dataclass
class StructureRequest:
    """A single structure to be refined by RefinementEngine."""

    id: str
    role: str
    kind: str
    input_xyz: Path
    forming_bonds: List[Tuple[int, int]] = field(default_factory=list)
    fallback_xyz: Optional[Path] = None
    original_seed_xyz: Optional[Path] = None
    source_stage: str = "S2"
    seed_state: str = "stable_minimum_seed"
    charge: int = 0
    multiplicity: int = 1
    atom_mapping: Optional[Path] = None
    parent_structure: Optional[Dict[str, Any]] = None

    structure_id: Optional[str] = None
    variant_id: Optional[str] = None
    branch_id: Optional[str] = None
    pathway_id: Optional[str] = None
    parent_structure_id: Optional[str] = None
    mapping_audit: Optional[str] = None
    mapping_required: bool = False
    s1_manifest: Optional[str] = None
    s1_ensemble_thermodynamics: Optional[Dict[str, Any]] = None
    s1_thermochemistry_status: Optional[str] = None
    ensemble_thermochemistry_correction_hartree: Optional[float] = None

    def get(self, key: str, default: Any = None) -> Any:
        """Mapping-like accessor for stage_scheduler compatibility."""

        return getattr(self, key, default)


@dataclass
class PreflightOutcome:
    """Per-structure outcome of Pass 0 preflight."""

    structure_id: str
    status: str
    error: Optional[str] = None
    input_xyz: Optional[Path] = None
    charge: int = 0
    multiplicity: int = 1
    forming_bonds: List[Tuple[int, int]] = field(default_factory=list)


@dataclass
class Pass1Outcome:
    """Per-structure outcome of Pass 1 primary calculation."""

    structure_id: str
    role: str
    kind: str

    warmup_status: str = "not_requested"
    warmup_used: bool = False
    warmup_xyz: Optional[Path] = None
    warmup_constraints: List[Dict[str, Any]] = field(default_factory=list)
    warmup_max_cycles: Optional[int] = None
    warmup_error: Optional[str] = None

    opt_status: str = "not_run"
    opt_xyz: Optional[Path] = None
    opt_output: Optional[Path] = None
    opt_energy_hartree: Optional[float] = None
    opt_initial_hessian: str = "model"
    opt_error: Optional[str] = None
    last_geometry_path: Optional[Path] = None
    stop_geometry_path: Optional[Path] = None
    optimization_monitor: Optional[Dict[str, Any]] = None
    optimization_monitor_error: Optional[Dict[str, Any]] = None
    primary_attempt_id: Optional[str] = None

    frequency_status: str = "not_run"
    frequency_output: Optional[Path] = None
    frequency_energy_hartree: Optional[float] = None
    frequencies_cm1: List[float] = field(default_factory=list)
    imaginary_frequencies_cm1: List[float] = field(default_factory=list)
    frequency_zero_point_energy_hartree: Optional[float] = None
    frequency_thermal_energy_hartree: Optional[float] = None
    frequency_enthalpy_hartree: Optional[float] = None
    frequency_gibbs_free_energy_hartree: Optional[float] = None
    frequency_gibbs_correction_hartree: Optional[float] = None
    frequency_hessian_path: Optional[Path] = None

    ts_classification: Optional[Dict[str, Any]] = None
    int_classification: Optional[Dict[str, Any]] = None
    minimum_classification: Optional[Dict[str, Any]] = None

    attempt_history: List[Dict[str, Any]] = field(default_factory=list)
    current_attempt_id: Optional[str] = None

    pass2_rescue_attempts: List[Dict[str, Any]] = field(default_factory=list)
    rescue_winning_method: Optional[str] = None
    canonical_attempt_id: Optional[str] = None
    canonical_xyz: Optional[Path] = None
    geometry_hash: Optional[str] = None
    canonical_frequency_status: str = "not_run"
    canonical_frequency_output: Optional[Path] = None
    canonical_frequency_energy_hartree: Optional[float] = None
    canonical_frequencies_cm1: List[float] = field(default_factory=list)
    canonical_imaginary_frequencies_cm1: List[float] = field(default_factory=list)
    canonical_zero_point_energy_hartree: Optional[float] = None
    canonical_thermal_energy_hartree: Optional[float] = None
    canonical_enthalpy_hartree: Optional[float] = None
    canonical_gibbs_free_energy_hartree: Optional[float] = None
    canonical_gibbs_correction_hartree: Optional[float] = None
    canonical_hessian_path: Optional[Path] = None
    sp_status: str = "not_run"
    sp_energy_hartree: Optional[float] = None
    sp_output: Optional[Path] = None
    thermochemistry: Optional[Dict[str, Any]] = None
    ml_usability: Optional[Dict[str, bool]] = None
    irc_status: str = "not_run"
    irc_endpoints: Optional[Dict[str, Any]] = None

    resolved_kind: Optional[str] = None
    resolved_identity: Optional[str] = None
    identity_status: str = "not_checked"

    status: str = "degraded"
    error: Optional[str] = None
