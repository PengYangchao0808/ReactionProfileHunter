# pyright: reportMissingImports=false
"""Adapters from persisted stage status payloads to UI render models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from rph_core.utils.ui_state import UiStatus, normalize_status


@dataclass
class UiTask:
    name: str
    status: UiStatus
    engine: str = ""
    method: str = ""
    basis: str = ""
    solvent: str = ""
    output: str = ""
    energy_hartree: float | None = None
    error: str = ""
    started_at: float | None = None
    finished_at: float | None = None


@dataclass
class UiBatch:
    """Normalized progress for homogeneous jobs such as S1 SP calculations."""

    id: str
    label: str
    status: UiStatus
    total: int
    done: int = 0
    failed: int = 0
    running: int = 0
    current: str = ""
    energy_hartree: float | None = None
    elapsed_seconds: float | None = None
    rate_per_minute: float | None = None
    eta_seconds: float | None = None


@dataclass
class UiAlert:
    severity: str
    message: str
    structure_id: str = ""
    action: str = ""


@dataclass
class UiStep:
    id: str
    label: str
    status: UiStatus
    index: int = 0
    total_steps: int = 0
    variant: str = ""
    purpose: str = ""
    elapsed_seconds: float = 0.0
    engine: str = ""
    method: str = ""
    output: str = ""
    error: str = ""
    detail: dict[str, Any] = field(default_factory=dict)


@dataclass
class UiActiveJob:
    id: str
    status: UiStatus
    batch: str = ""
    engine: str = ""
    method: str = ""
    attempt: int = 1
    nprocs: int = 0
    elapsed_seconds: float = 0.0
    started_at: float | None = None
    output: str = ""
    energy_hartree: float | None = None
    error: str = ""


@dataclass
class UiStructure:
    id: str
    kind: str
    source: str
    status: UiStatus
    current_task: str | None
    tasks: dict[str, UiTask] = field(default_factory=dict)
    usable_for_ml: bool = False
    error: str | None = None
    energy_hartree: float | None = None
    geometry_source: str = ""
    imaginary_frequencies: tuple[float, ...] = ()
    ts_quality_summary: str = ""
    ml_reason: str = ""
    fallback_source: str = ""


def _coerce_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _coerce_frequencies(value: Any) -> tuple[float, ...]:
    if not isinstance(value, (list, tuple)):
        return ()
    converted: list[float] = []
    for item in value:
        number = _coerce_float(item)
        if number is not None and number < 0:
            converted.append(number)
    return tuple(converted)


def adapt_batch(batch_id: str, data: dict[str, Any]) -> UiBatch:
    payload = dict(data or {})
    return UiBatch(
        id=str(batch_id),
        label=str(payload.get("label") or batch_id),
        status=normalize_status(payload.get("status")),
        total=int(payload.get("total") or 0),
        done=int(payload.get("done") or 0),
        failed=int(payload.get("failed") or 0),
        running=int(payload.get("running") or 0),
        current=str(payload.get("current") or ""),
        energy_hartree=_coerce_float(payload.get("energy_hartree")),
        elapsed_seconds=_coerce_float(payload.get("elapsed_seconds")),
        rate_per_minute=_coerce_float(payload.get("rate_per_minute")),
        eta_seconds=_coerce_float(payload.get("eta_seconds")),
    )


def adapt_step(step_id: str, data: dict[str, Any]) -> UiStep:
    payload = dict(data or {})
    known = {
        "id", "step", "step_id", "label", "status", "index", "total_steps",
        "variant", "purpose", "elapsed_seconds", "engine", "method", "output", "error",
    }
    return UiStep(
        id=str(step_id),
        label=str(payload.get("label") or step_id),
        status=normalize_status(payload.get("status")),
        index=int(payload.get("index") or 0),
        total_steps=int(payload.get("total_steps") or 0),
        variant=str(payload.get("variant") or ""),
        purpose=str(payload.get("purpose") or ""),
        elapsed_seconds=float(payload.get("elapsed_seconds") or 0.0),
        engine=str(payload.get("engine") or ""),
        method=str(payload.get("method") or ""),
        output=str(payload.get("output") or ""),
        error=str(payload.get("error") or ""),
        detail={key: value for key, value in payload.items() if key not in known},
    )


def adapt_active_job(job_id: str, data: dict[str, Any]) -> UiActiveJob:
    payload = dict(data or {})
    return UiActiveJob(
        id=str(job_id),
        status=normalize_status(payload.get("status")),
        batch=str(payload.get("batch") or ""),
        engine=str(payload.get("engine") or ""),
        method=str(payload.get("method") or ""),
        attempt=int(payload.get("attempt") or 1),
        nprocs=int(payload.get("nprocs") or 0),
        elapsed_seconds=float(payload.get("elapsed_seconds") or 0.0),
        started_at=_coerce_float(payload.get("started_at")),
        output=str(payload.get("output") or ""),
        energy_hartree=_coerce_float(payload.get("energy_hartree")),
        error=str(payload.get("error") or ""),
    )


def _task_from_structure_fields(
    name: str,
    status: Any,
    *,
    output: Any = None,
    energy_hartree: Any = None,
    error: Any = None,
) -> UiTask | None:
    if status in (None, "", "not_requested"):
        return None
    return UiTask(
        name=name,
        status=normalize_status(str(status)),
        output=str(output or ""),
        energy_hartree=_coerce_float(energy_hartree),
        error=str(error or ""),
    )


def adapt_s3_structures(raw_status: dict[str, Any]) -> list[UiStructure]:
    structures = raw_status.get("structures") or {}
    if isinstance(structures, dict):
        return [adapt_one(structure_id, data) for structure_id, data in structures.items()]
    if isinstance(structures, list):
        return [adapt_one(str(row.get("id") or row.get("structure_id") or ""), row) for row in structures if row]
    return []


def adapt_s4_structures(raw_status: dict[str, Any]) -> list[UiStructure]:
    structures = raw_status.get("structures") or []
    if isinstance(structures, list):
        return [adapt_one(str(row.get("id") or row.get("structure_id") or ""), row) for row in structures if row]
    if isinstance(structures, dict):
        return [adapt_one(structure_id, data) for structure_id, data in structures.items()]
    return []


def adapt_one(structure_id: str, data: dict[str, Any]) -> UiStructure:
    payload = dict(data or {})
    tasks: dict[str, UiTask] = {}
    for name, task in (payload.get("tasks") or {}).items():
        task_payload = dict(task or {})
        tasks[str(name)] = UiTask(
            name=str(name),
            status=normalize_status(task_payload.get("status")),
            engine=str(task_payload.get("engine") or ""),
            method=str(task_payload.get("method") or ""),
            basis=str(task_payload.get("basis") or ""),
            solvent=str(task_payload.get("solvent") or ""),
            output=str(task_payload.get("output") or ""),
            energy_hartree=_coerce_float(task_payload.get("energy_hartree")),
            error=str(task_payload.get("error") or ""),
            started_at=_coerce_float(task_payload.get("started_at")),
            finished_at=_coerce_float(task_payload.get("finished_at")),
        )
    if tasks.get("optimization") is None:
        task = _task_from_structure_fields(
            "optimization",
            payload.get("opt_status"),
            output=payload.get("opt_output"),
            error=payload.get("error"),
        )
        if task is not None:
            tasks[task.name] = task
    if tasks.get("frequency") is None:
        task = _task_from_structure_fields(
            "frequency",
            payload.get("frequency_status"),
            output=payload.get("frequency_output"),
            error=payload.get("error"),
        )
        if task is not None:
            tasks[task.name] = task
    if tasks.get("single_point") is None:
        task = _task_from_structure_fields(
            "single_point",
            payload.get("sp_status"),
            output=payload.get("sp_output"),
            energy_hartree=payload.get("sp_energy_hartree"),
            error=payload.get("error"),
        )
        if task is not None:
            tasks[task.name] = task
    energy_hartree = _coerce_float(payload.get("sp_energy_hartree"))
    if energy_hartree is None:
        energy_hartree = _coerce_float(payload.get("energy_hartree"))
    if energy_hartree is None and tasks.get("single_point") is not None:
        energy_hartree = tasks["single_point"].energy_hartree
    return UiStructure(
        id=str(structure_id),
        kind=str(payload.get("kind") or "minimum"),
        source=str(
            payload.get("input_source")
            or payload.get("source")
            or payload.get("source_stage")
            or "unknown"
        ),
        status=normalize_status(payload.get("status")),
        current_task=str(payload.get("current_task")) if payload.get("current_task") else None,
        tasks=tasks,
        usable_for_ml=bool(payload.get("usable_for_ml", False)),
        error=str(payload.get("error")) if payload.get("error") not in (None, "") else None,
        energy_hartree=energy_hartree,
        geometry_source=str(
            payload.get("geometry_source")
            or payload.get("input_source")
            or payload.get("source_stage")
            or ""
        ),
        imaginary_frequencies=_coerce_frequencies(
            payload.get("imaginary_frequencies_cm1")
            or payload.get("imaginary_frequencies")
            or payload.get("frequencies_cm1")
        ),
        ts_quality_summary=str(payload.get("ts_quality_summary") or ""),
        ml_reason=str(
            payload.get("ml_exclusion_reason")
            or payload.get("usable_for_ml_reason")
            or ""
        ),
        fallback_source=str(payload.get("fallback_source") or payload.get("fallback_xyz") or ""),
    )
