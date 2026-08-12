"""Small, stage-neutral models for one quantum-chemistry job."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


@dataclass(frozen=True)
class QCJobSpec:
    engine: str
    task: str
    method: str
    basis: str = ""
    aux_basis: str = ""
    solvent: Optional[str] = None
    solvent_model: Optional[str] = None
    route: str = ""
    charge: int = 0
    multiplicity: int = 1
    nproc: Optional[int] = None
    memory: Optional[str] = None
    maxcore: Optional[int] = None
    route_extras: str = ""
    max_cycles: Optional[int] = None
    grid: Optional[str] = None
    scf: Optional[str] = None
    scf_maxiter: Optional[int] = None
    timeout: Optional[int] = None
    bond_constraints: Tuple[Tuple[int, int, Optional[float]], ...] = ()
    allow_unconverged_geometry: bool = False
    # Hessian controls (Phase 1 — V4 unified refinement)
    initial_hessian: str = "model"  # "model" | "calculate" | "read" | "hybrid"
    hessian_filename: Optional[str] = None  # for InHess Read
    ts_mode: Optional[int] = None  # for TS_Mode {M N}
    recalc_hessian: Optional[int] = None  # for Recalc_Hess N
    trust: Optional[float] = None  # for Trust X (angstrom)
    neb_end_xyz: Optional[Path] = None
    neb_ts_guess_xyz: Optional[Path] = None


@dataclass
class QCJobResult:
    status: str
    input_xyz: Path
    output_xyz: Optional[Path] = None
    output_file: Optional[Path] = None
    energy_hartree: Optional[float] = None
    error: Optional[str] = None
    frequencies_cm1: Optional[Tuple[float, ...]] = None
    zero_point_energy_hartree: Optional[float] = None
    thermal_energy_hartree: Optional[float] = None
    enthalpy_hartree: Optional[float] = None
    gibbs_free_energy_hartree: Optional[float] = None
    gibbs_correction_hartree: Optional[float] = None
    failure_type: Optional[str] = None
    failure_evidence_lines: Optional[Tuple[str, ...]] = None
    optimization_converged: Optional[bool] = None
    stop_reason: Optional[str] = None
    returncode: Optional[int] = None
    stderr_text: Optional[str] = None
    timed_out: bool = False
    orca_version: Optional[str] = None
    extra: Optional[Dict[str, Any]] = None


@dataclass(frozen=True)
class OptimizationSpec:
    """Higher-level optimization spec used by RefinementEngine (Phase 2).

    Phase 1: this is a typed bag of parameters; RefinementEngine will translate
    it to QCJobSpec at dispatch time. Phase 2 may dispatch directly.
    """

    task: str = "opt"  # "opt" | "opt_ts"
    method: str = ""
    basis: str = ""
    aux_basis: str = ""
    route: str = "Opt"
    route_extras: str = ""
    initial_hessian: str = "model"  # "model" | "calculate" | "read" | "hybrid"
    hessian_filename: Optional[str] = None
    ts_mode: Optional[int] = None
    recalc_hessian: Optional[int] = None
    trust: Optional[float] = None
    max_cycles: Optional[int] = None
    bond_constraints: Tuple[Tuple[int, int, Optional[float]], ...] = ()
    allow_unconverged_geometry: bool = False


@dataclass(frozen=True)
class IRCJobSpec:
    """ORCA IRC job spec for TS endpoint discovery."""

    direction: str = "both"  # "both" | "forward" | "backward"
    max_iter: int = 5
    init_hessian: str = "read"  # "read" | "calc_anfreq" | "calc_numfreq"
    hessian_filename: Optional[str] = None
    hessian_mode: int = 0
    init_displacement_mode: str = "energy"  # "energy" | "length"
    initial_delta_energy_mEh: float = 2.0
    scale_initial_displacement: float = 0.1
    scale_steepest_descent: float = 0.15
    adaptive_step: bool = True
    method: str = ""
    basis: str = ""


@dataclass(frozen=True)
class NEBJobSpec:
    """ORCA double-ended NEB-TS job specification."""

    method: str = ""
    basis: str = ""
    aux_basis: str = ""
    route: str = "LOOSE-NEB-TS"
    end_xyz: Optional[Path] = None
    ts_guess_xyz: Optional[Path] = None
    charge: int = 0
    multiplicity: int = 1
    max_cycles: Optional[int] = None
    timeout: Optional[int] = None


@dataclass(frozen=True)
class SurfaceScanCoordinate:
    """One ORCA relaxed-scan coordinate (all atom indices are 0-based)."""

    kind: str = "B"  # "B" | "A" | "D"
    atoms: Tuple[int, ...] = ()
    start: float = 0.0
    end: float = 0.0
    steps: int = 0


@dataclass(frozen=True)
class SurfaceScanSpec:
    """Stage-neutral specification for an ORCA relaxed surface scan.

    ``scan_ts`` is intentionally limited by callers to one local coordinate;
    multi-coordinate reactions use a simultaneous relaxed scan and hand the
    selected peak to the ordinary OptTS workflow.
    """

    method: str
    basis: str = ""
    aux_basis: str = ""
    solvent: Optional[str] = None
    solvent_model: Optional[str] = None
    route_extras: str = ""
    charge: int = 0
    multiplicity: int = 1
    nproc: Optional[int] = None
    maxcore: Optional[int] = None
    max_cycles: Optional[int] = None
    timeout: Optional[int] = None
    coordinates: Tuple[SurfaceScanCoordinate, ...] = ()
    simultaneous: bool = False
    scan_ts: bool = False
    full_scan: bool = True
