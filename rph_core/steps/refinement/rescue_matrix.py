"""Rescue method-family matrix for the unified refinement engine.

R1/R2/R3 are METHOD FAMILIES, not pipeline phases:

    R1 Restart : recompute an exact Hessian at the checkpoint and continue
                 (Calc_Hess true + Recalc_Hess 5; default step size)
    R2 Mode    : follow the imaginary-mode direction (TS_Mode / displacement)
    R3 Exact   : every-step exact Hessian (Recalc_Hess 1) final fallback

The matrix maps ``(structure kind, failure type) -> ordered method list``.
A cell is executed in order; the first method producing a valid candidate
terminates the cell.  F5 (SCF) and F6 (crash/timeout) exit directly without
rescue — each incident is investigated individually.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Sequence, Tuple


class FailureType(str, Enum):
    """Normalized task-error classification feeding the rescue matrix."""

    F1_GEOMETRY_NOT_CONVERGED = "F1"
    F2_HIGHER_ORDER_SADDLE = "F2"
    F3_TS_NO_IMAGINARY_MODE = "F3"
    F4_MINIMUM_WITH_IMAGINARY = "F4"
    F5_SCF_NOT_CONVERGED = "F5"
    F6_CRASH_TIMEOUT_OTHER = "F6"
    F7_COLLAPSED_TO_PRODUCT = "F7"

    @property
    def label(self) -> str:
        return {
            FailureType.F1_GEOMETRY_NOT_CONVERGED: "geometry not converged (maxiter)",
            FailureType.F2_HIGHER_ORDER_SADDLE: "higher-order saddle (>=2 imaginary)",
            FailureType.F3_TS_NO_IMAGINARY_MODE: "TS without imaginary mode",
            FailureType.F4_MINIMUM_WITH_IMAGINARY: "minimum with imaginary mode",
            FailureType.F5_SCF_NOT_CONVERGED: "SCF not converged",
            FailureType.F6_CRASH_TIMEOUT_OTHER: "crash / timeout / other",
            FailureType.F7_COLLAPSED_TO_PRODUCT: "INT collapsed to product minimum",
        }[self]


class StructureKind(str, Enum):
    TS = "ts"
    INT = "int"
    MINIMUM = "minimum"

    @classmethod
    def from_request(cls, kind: str, role: str) -> StructureKind:
        if str(kind).lower() == "ts":
            return cls.TS
        if str(role).lower() == "intermediate":
            return cls.INT
        return cls.MINIMUM


class MethodFamily(str, Enum):
    R1_RESTART = "R1"
    R2_MODE = "R2"
    R3_EXACT = "R3"
    R4_IRC = "R4"


class RescueMethod(str, Enum):
    FRESH_HESSIAN_RESTART = "fresh_hessian_restart"
    FRESH_HESSIAN_MODE_MONITOR = "fresh_hessian_mode_monitor"
    TS_MODE_DIRECTED = "ts_mode_directed"
    MODE_DISPLACEMENT = "mode_displacement"
    SADDLE_BREAK = "saddle_break"
    CALCALL_OPT = "calcall_opt"
    IRC_MIDPOINT_RECOVERY = "irc_midpoint_recovery"

    @property
    def family(self) -> MethodFamily:
        return {
            RescueMethod.FRESH_HESSIAN_RESTART: MethodFamily.R1_RESTART,
            RescueMethod.FRESH_HESSIAN_MODE_MONITOR: MethodFamily.R1_RESTART,
            RescueMethod.TS_MODE_DIRECTED: MethodFamily.R2_MODE,
            RescueMethod.MODE_DISPLACEMENT: MethodFamily.R2_MODE,
            RescueMethod.SADDLE_BREAK: MethodFamily.R2_MODE,
            RescueMethod.CALCALL_OPT: MethodFamily.R3_EXACT,
            RescueMethod.IRC_MIDPOINT_RECOVERY: MethodFamily.R4_IRC,
        }[self]


@dataclass(frozen=True)
class MethodParams:
    """Fixed ORCA parameters for a rescue method."""

    max_cycles: int = 30
    trust: Optional[float] = None
    recalc_hessian: Optional[int] = None
    calc_hess: bool = False
    ts_mode: bool = False
    tight_opt: bool = False
    displacement_step_angstrom: float = 0.30
    displacement_sign: Tuple[str, ...] = ("plus", "minus")
    mode_min_overlap: float = 0.35
    mode_min_overlap_margin: float = 0.08
    irc_max_iter: int = 5
    irc_direction: str = "both"
    shoulder_energy_window_kcal_mol: float = 2.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "max_cycles": self.max_cycles,
            "trust": self.trust,
            "recalc_hessian": self.recalc_hessian,
            "calc_hess": self.calc_hess,
            "ts_mode": self.ts_mode,
            "tight_opt": self.tight_opt,
            "displacement_step_angstrom": self.displacement_step_angstrom,
            "mode_min_overlap": self.mode_min_overlap,
            "mode_min_overlap_margin": self.mode_min_overlap_margin,
            "irc_max_iter": self.irc_max_iter,
            "irc_direction": self.irc_direction,
            "shoulder_energy_window_kcal_mol": self.shoulder_energy_window_kcal_mol,
        }


# Default method parameterization (config overrides via refinement.common.rescue).
DEFAULT_METHOD_PARAMS: Dict[RescueMethod, MethodParams] = {
    RescueMethod.FRESH_HESSIAN_RESTART: MethodParams(
        max_cycles=30, calc_hess=True, recalc_hessian=5
    ),
    RescueMethod.FRESH_HESSIAN_MODE_MONITOR: MethodParams(
        max_cycles=30, calc_hess=True, recalc_hessian=5, ts_mode=False
    ),
    RescueMethod.TS_MODE_DIRECTED: MethodParams(
        max_cycles=12, trust=0.15, ts_mode=True
    ),
    RescueMethod.MODE_DISPLACEMENT: MethodParams(
        max_cycles=30, calc_hess=True, tight_opt=True
    ),
    RescueMethod.SADDLE_BREAK: MethodParams(max_cycles=30, calc_hess=True),
    RescueMethod.CALCALL_OPT: MethodParams(max_cycles=30, recalc_hessian=1),
    RescueMethod.IRC_MIDPOINT_RECOVERY: MethodParams(
        max_cycles=30, irc_max_iter=5, irc_direction="both"
    ),
}

# Method-level budget: an exact-Hessian refresh every `recalc_hessian` cycles
# over `max_cycles` yields ~6 evaluations; beyond that the method is deemed
# failing and the next cell method takes over.
METHOD_BUDGET_HESSIAN_EVALUATIONS = 6


@dataclass(frozen=True)
class RescuePlan:
    """A single matrix cell: ordered methods plus metadata."""

    failure_type: FailureType
    kind: StructureKind
    methods: Tuple[RescueMethod, ...] = field(default_factory=tuple)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "failure_type": self.failure_type.value,
            "failure_label": self.failure_type.label,
            "kind": self.kind.value,
            "methods": [method.value for method in self.methods],
        }


# The rescue matrix (v3.1).  Order within a cell is execution order.
# F5/F6 intentionally have no cells: they exit directly without rescue.
_RESCUE_MATRIX: Dict[Tuple[FailureType, StructureKind], Tuple[RescueMethod, ...]] = {
    (FailureType.F1_GEOMETRY_NOT_CONVERGED, StructureKind.TS): (
        RescueMethod.FRESH_HESSIAN_RESTART,
        RescueMethod.TS_MODE_DIRECTED,
        RescueMethod.CALCALL_OPT,
    ),
    (FailureType.F1_GEOMETRY_NOT_CONVERGED, StructureKind.INT): (
        RescueMethod.FRESH_HESSIAN_RESTART,
        RescueMethod.CALCALL_OPT,
    ),
    (FailureType.F1_GEOMETRY_NOT_CONVERGED, StructureKind.MINIMUM): (
        RescueMethod.FRESH_HESSIAN_RESTART,
        RescueMethod.CALCALL_OPT,
    ),
    (FailureType.F2_HIGHER_ORDER_SADDLE, StructureKind.TS): (
        RescueMethod.SADDLE_BREAK,
        RescueMethod.TS_MODE_DIRECTED,
        RescueMethod.CALCALL_OPT,
    ),
    (FailureType.F3_TS_NO_IMAGINARY_MODE, StructureKind.TS): (
        RescueMethod.FRESH_HESSIAN_MODE_MONITOR,
        RescueMethod.TS_MODE_DIRECTED,
        RescueMethod.CALCALL_OPT,
    ),
    (FailureType.F4_MINIMUM_WITH_IMAGINARY, StructureKind.INT): (
        RescueMethod.MODE_DISPLACEMENT,
    ),
    (FailureType.F4_MINIMUM_WITH_IMAGINARY, StructureKind.MINIMUM): (
        RescueMethod.MODE_DISPLACEMENT,
    ),
    (FailureType.F7_COLLAPSED_TO_PRODUCT, StructureKind.INT): (
        RescueMethod.IRC_MIDPOINT_RECOVERY,
    ),
}


def lookup_plan(failure_type: FailureType, kind: StructureKind) -> Optional[RescuePlan]:
    """Look up the rescue plan for a (failure, kind) cell; None for F5/F6/empty."""
    methods = _RESCUE_MATRIX.get((failure_type, kind))
    if not methods:
        return None
    return RescuePlan(failure_type=failure_type, kind=kind, methods=methods)


def all_plans() -> Sequence[RescuePlan]:
    """All non-empty cells — used for exhaustive matrix testing."""
    return [
        RescuePlan(failure_type=failure, kind=kind, methods=methods)
        for (failure, kind), methods in sorted(
            _RESCUE_MATRIX.items(), key=lambda item: (item[0][0].value, item[0][1].value)
        )
    ]


def methods_params(
    config: Optional[Dict[str, Any]],
    method: RescueMethod,
) -> MethodParams:
    """Resolve method params, honoring ``refinement.common.rescue.methods``."""
    defaults = DEFAULT_METHOD_PARAMS[method]
    if not config:
        return defaults
    rescue_cfg = dict(dict(config.get("refinement", {}) or {}).get("common", {}).get("rescue", {}) or {})
    method_cfg = dict(rescue_cfg.get("methods", {}).get(method.value, {}) or {})
    if not method_cfg:
        return defaults
    return MethodParams(
        max_cycles=int(method_cfg.get("max_cycles", defaults.max_cycles)),
        trust=_optional_float(method_cfg.get("trust"), defaults.trust),
        recalc_hessian=_optional_int(method_cfg.get("recalc_hessian"), defaults.recalc_hessian),
        calc_hess=bool(method_cfg.get("calc_hess", defaults.calc_hess)),
        ts_mode=bool(method_cfg.get("ts_mode", defaults.ts_mode)),
        tight_opt=bool(method_cfg.get("tight_opt", defaults.tight_opt)),
        displacement_step_angstrom=float(
            method_cfg.get("displacement_step_angstrom", defaults.displacement_step_angstrom)
        ),
        mode_min_overlap=float(method_cfg.get("mode_min_overlap", defaults.mode_min_overlap)),
        mode_min_overlap_margin=float(
            method_cfg.get("mode_min_overlap_margin", defaults.mode_min_overlap_margin)
        ),
        irc_max_iter=int(method_cfg.get("irc_max_iter", defaults.irc_max_iter)),
        irc_direction=str(method_cfg.get("irc_direction", defaults.irc_direction)),
        shoulder_energy_window_kcal_mol=float(
            method_cfg.get(
                "shoulder_energy_window_kcal_mol",
                defaults.shoulder_energy_window_kcal_mol,
            )
        ),
    )


def _optional_float(value: Any, fallback: Optional[float]) -> Optional[float]:
    if value is None or str(value).strip() in ("", "null", "none"):
        return fallback
    try:
        return float(value)
    except (TypeError, ValueError):
        return fallback


def _optional_int(value: Any, fallback: Optional[int]) -> Optional[int]:
    if value is None or str(value).strip() in ("", "null", "none"):
        return fallback
    try:
        return int(value)
    except (TypeError, ValueError):
        return fallback


__all__ = [
    "DEFAULT_METHOD_PARAMS",
    "METHOD_BUDGET_HESSIAN_EVALUATIONS",
    "FailureType",
    "MethodFamily",
    "MethodParams",
    "RescueMethod",
    "RescuePlan",
    "StructureKind",
    "all_plans",
    "lookup_plan",
    "methods_params",
]
