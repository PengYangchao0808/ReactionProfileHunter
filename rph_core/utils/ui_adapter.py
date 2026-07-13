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


def _coerce_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


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
    )
