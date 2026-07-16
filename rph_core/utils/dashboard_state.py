"""Thread-safe in-memory state reducer for the embedded V4 dashboard."""

from __future__ import annotations

import copy
from datetime import datetime
import threading
import time
from dataclasses import dataclass, field
from typing import Any


_TERMINAL_EVENTS = {"structure_finished", "structure_failed"}


@dataclass
class DashboardState:
    """Mutable state owned by :class:`DashboardStateReducer`."""

    stages: dict[str, dict[str, Any]] = field(default_factory=dict)
    active_stage: str | None = None
    reaction_id: str = ""
    condition_signature: str = ""
    resources: dict[str, Any] = field(default_factory=dict)
    log_path: str = ""
    started_at: float = field(default_factory=time.time)
    updated_at: float = field(default_factory=time.time)
    revision: int = 0
    final_status: str = "running"


class DashboardStateReducer:
    """Reduce durable progress events into one dashboard snapshot.

    Event callbacks may arrive from QC worker threads.  This class performs no
    terminal I/O and exposes deep-copy snapshots so rendering never observes a
    partially updated structure.
    """

    def __init__(self, *, log_path: str = "") -> None:
        self._state = DashboardState(log_path=str(log_path or ""))
        self._lock = threading.RLock()

    @property
    def revision(self) -> int:
        with self._lock:
            return self._state.revision

    def snapshot(self) -> DashboardState:
        with self._lock:
            return copy.deepcopy(self._state)

    def finish(self, status: str) -> None:
        with self._lock:
            self._state.final_status = str(status or "complete")
            self._state.updated_at = time.time()
            self._state.revision += 1

    def apply(self, event: str, record: dict[str, Any]) -> None:
        payload = dict(record or {})
        stage_name = str(payload.get("stage") or "").strip().upper()
        if not stage_name:
            return
        with self._lock:
            self._capture_run_metadata(payload)
            stage = self._state.stages.setdefault(stage_name, self._new_stage(stage_name))
            stage["updated_at"] = payload.get("timestamp") or time.time()
            stage["meta"].update(payload)

            if event.startswith("batch_job_"):
                self._apply_batch_job(stage, event, payload)
            elif event.startswith("batch_"):
                self._apply_batch(stage, event, payload)
            elif event.startswith("step_"):
                self._apply_step(stage, event, payload)
            elif event == "funnel_progress":
                variant = str(payload.get("variant") or "default")
                funnel = stage["funnel"].setdefault(variant, {})
                funnel[str(payload.get("step") or "unknown")] = payload.get("candidates")
            elif event == "science_summary":
                variant = str(payload.get("variant") or "")
                if variant:
                    stage["science"].setdefault("variants", {})[variant] = payload
                stage["science"].update(payload)
            elif event == "decision":
                stage["decisions"].append(payload)
                stage["decisions"] = stage["decisions"][-20:]
            elif event == "alert":
                stage["alerts"].append(payload)
                stage["alerts"] = stage["alerts"][-20:]
            elif self._is_stage_start(stage_name, event):
                stage["status"] = "running"
                stage["started_at"] = stage.get("started_at") or _epoch(payload.get("timestamp"))
                self._state.active_stage = stage_name
            elif self._is_stage_finish(stage_name, event):
                stage["status"] = str(payload.get("status") or ("failed" if event.endswith("failed") else "complete"))
                stage["finished_at"] = _epoch(payload.get("timestamp"))
                if self._state.active_stage == stage_name:
                    self._state.active_stage = None
            else:
                self._apply_structure(stage, event, payload)

            self._state.updated_at = time.time()
            self._state.revision += 1

    @staticmethod
    def _new_stage(name: str) -> dict[str, Any]:
        return {
            "name": name,
            "status": "pending",
            "started_at": None,
            "finished_at": None,
            "updated_at": None,
            "meta": {},
            "structures": {},
            "steps": {},
            "current_step": None,
            "batches": {},
            "funnel": {},
            "science": {},
            "alerts": [],
            "decisions": [],
        }

    def _capture_run_metadata(self, payload: dict[str, Any]) -> None:
        self._state.reaction_id = str(
            payload.get("reaction_id") or payload.get("rx_id") or self._state.reaction_id
        )
        self._state.condition_signature = str(
            payload.get("condition_signature") or self._state.condition_signature
        )
        if payload.get("resources"):
            self._state.resources = dict(payload["resources"])

    @staticmethod
    def _is_stage_start(stage: str, event: str) -> bool:
        return event == "stage_started" or event == f"{stage.lower()}_started"

    @staticmethod
    def _is_stage_finish(stage: str, event: str) -> bool:
        return event == "stage_finished" or (
            event.startswith(stage.lower()) and (event.endswith("_completed") or event.endswith("_failed"))
        )

    @staticmethod
    def _apply_step(stage: dict[str, Any], event: str, payload: dict[str, Any]) -> None:
        step_id = str(payload.get("step") or payload.get("step_id") or "step")
        variant = str(payload.get("variant") or "")
        key = f"{variant}:{step_id}" if variant else step_id
        row = dict(stage["steps"].get(key) or {})
        row.update(payload)
        row["key"] = key
        if event in {"step_started", "step_heartbeat"}:
            row["status"] = "running"
            stage["current_step"] = key
        elif event == "step_failed":
            row["status"] = "failed"
            if stage["current_step"] == key:
                stage["current_step"] = None
        elif event == "step_finished":
            row["status"] = str(payload.get("status") or "complete")
            if stage["current_step"] == key:
                stage["current_step"] = None
        stage["steps"][key] = row

    @staticmethod
    def _apply_batch(stage: dict[str, Any], event: str, payload: dict[str, Any]) -> None:
        batch_id = str(payload.get("batch") or payload.get("batch_id") or "batch")
        row = dict(stage["batches"].get(batch_id) or {"active_jobs": {}})
        active = dict(row.get("active_jobs") or {})
        row.update(payload)
        row["active_jobs"] = active
        if event == "batch_started":
            row["status"] = "running"
            row.setdefault("done", 0)
            row.setdefault("failed", 0)
        elif event == "batch_finished":
            row["status"] = str(payload.get("status") or "complete")
        else:
            row["status"] = str(payload.get("status") or row.get("status") or "running")
        stage["batches"][batch_id] = row

    @staticmethod
    def _apply_batch_job(stage: dict[str, Any], event: str, payload: dict[str, Any]) -> None:
        batch_id = str(payload.get("batch") or "batch")
        job_id = str(payload.get("job_id") or payload.get("current") or "job")
        batch = dict(stage["batches"].get(batch_id) or {"label": batch_id})
        active = dict(batch.get("active_jobs") or {})
        row = dict(active.get(job_id) or {})
        row.update(payload)
        row["id"] = job_id
        if event in {"batch_job_started", "batch_job_heartbeat", "batch_job_retry"}:
            row["status"] = "running"
            active[job_id] = row
        else:
            active.pop(job_id, None)
        batch["active_jobs"] = active
        batch["running"] = len(active)
        stage["batches"][batch_id] = batch

    @staticmethod
    def _apply_structure(stage: dict[str, Any], event: str, payload: dict[str, Any]) -> None:
        structure_id = str(payload.get("structure_id") or payload.get("id") or "").strip()
        if not structure_id:
            return
        row = dict(stage["structures"].get(structure_id) or {"tasks": {}, "status": "pending"})
        tasks = copy.deepcopy(row.get("tasks") or {})
        row.update(payload)
        row["tasks"] = tasks
        row["id"] = structure_id

        task_name, _, action = event.rpartition("_")
        if action in {"started", "finished", "skipped"} and task_name in {
            "optimization", "frequency", "single_point"
        }:
            task = dict(tasks.get(task_name) or {})
            task.update(payload)
            task["status"] = "running" if action == "started" else str(
                payload.get("status") or ("skipped" if action == "skipped" else "complete")
            )
            if action == "started":
                task["started_at"] = _epoch(payload.get("timestamp"))
                row["current_task"] = task_name
                row["status"] = "running"
            else:
                task["finished_at"] = _epoch(payload.get("timestamp"))
                if row.get("current_task") == task_name:
                    row["current_task"] = None
            tasks[task_name] = task
        elif event == "structure_started":
            row["status"] = "running"
            row.setdefault("started_at", _epoch(payload.get("timestamp")))
        elif event in _TERMINAL_EVENTS:
            row["status"] = str(payload.get("status") or ("failed" if event == "structure_failed" else "complete"))
            row["current_task"] = None
            row["finished_at"] = _epoch(payload.get("timestamp"))
        stage["structures"][structure_id] = row


def _epoch(value: Any) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str) and value:
        try:
            return datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp()
        except ValueError:
            pass
    return time.time()
