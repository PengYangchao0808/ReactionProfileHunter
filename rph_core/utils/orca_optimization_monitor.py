"""Live, conservative detection of stalled ORCA geometry optimizations.

The monitor deliberately makes no chemical decision.  It only recognizes a
repeated optimizer pattern (energy reversals without gradient/step progress)
and asks the caller to stop the current attempt.  The caller remains
responsible for preserving the last geometry and deciding which rescue route
is appropriate.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any


_CYCLE_RE = re.compile(
    r"GEOMETRY\s+OPTIMIZATION\s+CYCLE\s+(?P<cycle>\d+)(?P<body>.*?)(?=GEOMETRY\s+OPTIMIZATION\s+CYCLE\s+\d+|\Z)",
    re.IGNORECASE | re.DOTALL,
)
_FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
_CONVERGENCE_RE = re.compile(
    r"^\s*(?P<label>RMS gradient|MAX gradient|RMS step|MAX step)\s+"
    rf"(?P<value>{_FLOAT_RE})\s+(?P<tolerance>{_FLOAT_RE})\s+"
    r"(?P<status>YES|NO)\s*$",
    re.IGNORECASE | re.MULTILINE,
)
_ENERGY_CHANGE_RE = re.compile(
    rf"(?:Actually observed|Energy)\s+change\s*(?:\.{{2,}})?\s*(?P<value>{_FLOAT_RE})",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class OptimizationMonitorConfig:
    """Configurable thresholds for a live optimization monitor."""

    enabled: bool = True
    poll_seconds: float = 20.0
    min_cycles: int = 12
    window_cycles: int = 16
    min_energy_reversals: int = 4
    energy_change_floor_hartree: float = 1.0e-5
    min_gradient_improvement_fraction: float = 0.15
    step_tolerance_ratio: float = 4.0
    min_excessive_steps: int = 6


def parse_orca_optimization_cycles(content: str) -> list[dict[str, float | int | None]]:
    """Return complete ORCA geometry-cycle convergence records.

    A partially flushed last cycle is ignored unless it includes all four
    convergence rows.  This makes it safe to call while ORCA is writing its
    output file.
    """

    records: list[dict[str, float | int | None]] = []
    for match in _CYCLE_RE.finditer(content or ""):
        values: dict[str, float] = {}
        tolerances: dict[str, float] = {}
        for row in _CONVERGENCE_RE.finditer(match.group("body")):
            label = " ".join(row.group("label").lower().split())
            values[label] = float(row.group("value"))
            tolerances[label] = float(row.group("tolerance"))
        required = {"rms gradient", "max gradient", "rms step", "max step"}
        if not required.issubset(values):
            continue
        energy_matches = list(_ENERGY_CHANGE_RE.finditer(match.group("body")))
        records.append(
            {
                "cycle": int(match.group("cycle")),
                "energy_change": (
                    float(energy_matches[-1].group("value")) if energy_matches else None
                ),
                "rms_gradient": values["rms gradient"],
                "max_gradient": values["max gradient"],
                "rms_step": values["rms step"],
                "max_step": values["max step"],
                "rms_gradient_tolerance": tolerances["rms gradient"],
                "max_gradient_tolerance": tolerances["max gradient"],
                "rms_step_tolerance": tolerances["rms step"],
                "max_step_tolerance": tolerances["max step"],
            }
        )
    return records


def detect_oscillation(
    records: list[dict[str, float | int | None]],
    config: OptimizationMonitorConfig,
) -> dict[str, Any] | None:
    """Return an auditable trigger payload only for confirmed oscillation."""

    if not config.enabled or len(records) < max(config.min_cycles, config.window_cycles):
        return None
    window = records[-config.window_cycles :]
    energy_changes = [
        float(record["energy_change"])
        for record in window
        if record.get("energy_change") is not None
        and abs(float(record["energy_change"])) >= config.energy_change_floor_hartree
    ]
    reversals = sum(
        1
        for previous, current in zip(energy_changes, energy_changes[1:])
        if previous * current < 0.0
    )
    first_gradient = float(window[0]["rms_gradient"])
    last_gradient = float(window[-1]["rms_gradient"])
    gradient_improvement = (
        (first_gradient - last_gradient) / first_gradient if first_gradient > 0.0 else 0.0
    )
    excessive_steps = sum(
        1
        for record in window
        if float(record["max_step"])
        > config.step_tolerance_ratio * float(record["max_step_tolerance"])
    )
    if (
        reversals < config.min_energy_reversals
        or gradient_improvement >= config.min_gradient_improvement_fraction
        or excessive_steps < config.min_excessive_steps
    ):
        return None
    return {
        "reason": "oscillation",
        "cycle": int(window[-1]["cycle"]),
        "window_cycles": [int(record["cycle"]) for record in window],
        "energy_reversals": reversals,
        "gradient_improvement_fraction": gradient_improvement,
        "excessive_step_count": excessive_steps,
    }


def newest_orca_output(work_dir: Path) -> Path | None:
    """Find the currently written ORCA output without assuming its basename."""

    candidates = [path for path in Path(work_dir).glob("*.out") if path.is_file()]
    return max(candidates, key=lambda path: path.stat().st_mtime) if candidates else None


__all__ = [
    "OptimizationMonitorConfig",
    "detect_oscillation",
    "newest_orca_output",
    "parse_orca_optimization_cycles",
]
