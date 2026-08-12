"""S3/S4 normalized step codes — single source of truth for stage sub-steps.

Every refinement stage (S3, S4) exposes the same ordered step taxonomy:

    S*.0 preflight  → S*.1 warmup  → S*.2 primary → S*.3 freq
    → S*.4 rescue (sub-levels per method family) → S*.5 canonical
    → S*.6 properties → S*.7 finalize

All logs, manifests, UI payloads and resume records reference these codes
instead of the legacy ``Pass 0-3`` / ``wave2_pass2_pass3`` vocabulary.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional

# Legacy phase label that must not appear in V4 step metadata.
LEGACY_PHASE_LABEL = "wave2_pass2_pass3"


class StepCode(str, Enum):
    """Canonical step ordinal within a refinement stage.

    ``qualified(stage)`` renders the external id, e.g. ``S3.4``.
    """

    PREFLIGHT = "0"
    WARMUP = "1"
    PRIMARY = "2"
    FREQ = "3"
    RESCUE = "4"
    CANONICAL = "5"
    PROPERTIES = "6"
    FINALIZE = "7"

    @property
    def name(self) -> str:  # noqa: A003 - machine name for payloads/logs
        return {
            StepCode.PREFLIGHT: "preflight",
            StepCode.WARMUP: "warmup",
            StepCode.PRIMARY: "primary_opt",
            StepCode.FREQ: "freq",
            StepCode.RESCUE: "rescue",
            StepCode.CANONICAL: "canonical",
            StepCode.PROPERTIES: "properties",
            StepCode.FINALIZE: "finalize",
        }[self]

    def qualified(self, stage: str) -> str:
        """Stage-qualified code, e.g. ``S3.2`` for PRIMARY at stage S3."""
        return f"{str(stage).upper()}.{self.value}"

    @classmethod
    def ordered(cls) -> tuple[StepCode, ...]:
        return (
            cls.PREFLIGHT,
            cls.WARMUP,
            cls.PRIMARY,
            cls.FREQ,
            cls.RESCUE,
            cls.CANONICAL,
            cls.PROPERTIES,
            cls.FINALIZE,
        )


# Rescue sub-levels within S*.4 (diagnosis + method families R1/R2/R3).
RESCUE_SUB_DIAGNOSIS = "0"
RESCUE_SUB_R1_RESTART = "1"
RESCUE_SUB_R2_MODE = "2"
RESCUE_SUB_R3_EXACT = "3"


def rescue_sub_code(stage: str, sub: str) -> str:
    """Qualified rescue sub-level, e.g. ``S3.4.2`` for R2 at stage S3."""
    return f"{str(stage).upper()}.4.{sub}"


def rescue_sub_for_family(stage: str, family: str) -> str:
    """Map a method family (R1/R2/R3) to its S*.4 sub-level code."""
    if str(family).upper().startswith("R1"):
        return rescue_sub_code(stage, RESCUE_SUB_R1_RESTART)
    if str(family).upper().startswith("R2"):
        return rescue_sub_code(stage, RESCUE_SUB_R2_MODE)
    if str(family).upper().startswith("R3"):
        return rescue_sub_code(stage, RESCUE_SUB_R3_EXACT)
    return rescue_sub_code(stage, RESCUE_SUB_DIAGNOSIS)


@dataclass(frozen=True)
class StepSpec:
    """Normalized contract for a single S*.x step."""

    code: str  # qualified id, e.g. "S3.2"
    name: str  # machine name, e.g. "primary_opt"
    entry_gate: str  # condition to enter
    exit_gate: str  # branch(es) on exit
    outputs: tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "code": self.code,
            "name": self.name,
            "entry_gate": self.entry_gate,
            "exit_gate": self.exit_gate,
            "outputs": list(self.outputs),
        }


def step_specs(stage: str) -> tuple[StepSpec, ...]:
    """Ordered normalized step catalog for a refinement stage (S3/S4)."""
    s = str(stage).upper()
    return (
        StepSpec(
            code=StepCode.PREFLIGHT.qualified(s),
            name=StepCode.PREFLIGHT.name,
            entry_gate="structure request ready",
            exit_gate="validation ok -> warmup / primary",
            outputs=("provenance.json", "preflight record"),
        ),
        StepSpec(
            code=StepCode.WARMUP.qualified(s),
            name=StepCode.WARMUP.name,
            entry_gate="INT/TS role (precursor/product skip)",
            exit_gate="warmup geometry (complete or partial) -> primary",
            outputs=("warmup/opt.xyz", "bond constraints"),
        ),
        StepSpec(
            code=StepCode.PRIMARY.qualified(s),
            name=StepCode.PRIMARY.name,
            entry_gate="warmup geometry or original seed",
            exit_gate="converged -> freq; maxiter/saddle -> rescue; other -> finalize(failed)",
            outputs=("opt.xyz", "opt.out", "stop_reason"),
        ),
        StepSpec(
            code=StepCode.FREQ.qualified(s),
            name=StepCode.FREQ.name,
            entry_gate="converged primary geometry",
            exit_gate="classification done -> canonical or rescue",
            outputs=("freq.out", "hessian", "classification"),
        ),
        StepSpec(
            code=StepCode.RESCUE.qualified(s),
            name=StepCode.RESCUE.name,
            entry_gate="primary/freq failure or classification anomaly",
            exit_gate="valid candidate -> canonical; exhausted -> finalize(failed)",
            outputs=("rescue method records", "S*.4.x sub-levels"),
        ),
        StepSpec(
            code=StepCode.CANONICAL.qualified(s),
            name=StepCode.CANONICAL.name,
            entry_gate="valid candidate (primary or rescue)",
            exit_gate="canonical.xyz + canonical freq -> properties",
            outputs=("canonical.xyz", "selection record"),
        ),
        StepSpec(
            code=StepCode.PROPERTIES.qualified(s),
            name=StepCode.PROPERTIES.name,
            entry_gate="canonical geometry ready",
            exit_gate="ml_usability decided -> finalize",
            outputs=("sp.out", "thermochemistry", "ml_usability"),
        ),
        StepSpec(
            code=StepCode.FINALIZE.qualified(s),
            name=StepCode.FINALIZE.name,
            entry_gate="all prior steps settled",
            exit_gate="manifest written -> batch fail-fast check",
            outputs=("manifest.json", "summary", "health report"),
        ),
    )


@dataclass
class StepRecord:
    """Manifest record for one executed S*.x step."""

    code: str
    name: str
    status: str  # complete | failed | skipped | not_run
    duration_s: float = 0.0
    attempt_id: Optional[str] = None
    reason: Optional[str] = None
    sub_steps: List[Dict[str, Any]] = field(default_factory=list)
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        payload: Dict[str, Any] = {
            "code": self.code,
            "name": self.name,
            "status": self.status,
        }
        if self.duration_s:
            payload["duration_s"] = round(float(self.duration_s), 2)
        if self.attempt_id:
            payload["attempt_id"] = self.attempt_id
        if self.reason:
            payload["reason"] = self.reason
        if self.sub_steps:
            payload["sub_steps"] = self.sub_steps
        if self.extra:
            payload.update(self.extra)
        return payload


__all__ = [
    "LEGACY_PHASE_LABEL",
    "RESCUE_SUB_DIAGNOSIS",
    "RESCUE_SUB_R1_RESTART",
    "RESCUE_SUB_R2_MODE",
    "RESCUE_SUB_R3_EXACT",
    "StepCode",
    "StepRecord",
    "StepSpec",
    "rescue_sub_code",
    "rescue_sub_for_family",
    "step_specs",
]
