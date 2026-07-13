"""Small, stage-neutral models for one quantum-chemistry job."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple


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
    timeout: Optional[int] = None


@dataclass
class QCJobResult:
    status: str
    input_xyz: Path
    output_xyz: Optional[Path] = None
    output_file: Optional[Path] = None
    energy_hartree: Optional[float] = None
    error: Optional[str] = None
    frequencies_cm1: Optional[Tuple[float, ...]] = None
