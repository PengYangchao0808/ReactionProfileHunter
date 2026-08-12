"""Unified durable progress reporting for V4 S0-S3 stages."""

from __future__ import annotations

import json
import logging
import os
import threading
import time
from pathlib import Path
from typing import Any, Callable, Dict, Optional

from rph_core.utils.run_id import RUN_ID_FIELD

logger = logging.getLogger(__name__)


def _now() -> float:
    return time.time()


def _heartbeat_interval_seconds(
    config: Optional[Dict[str, Any]],
    default_fields: Optional[Dict[str, Any]],
) -> float:
    ui_cfg: Dict[str, Any] = {}
    if default_fields is not None:
        ui_cfg.update(dict(default_fields.get("ui") or {}))
    if config is not None:
        ui_cfg.update(dict(config.get("ui", {}) or {}))
    try:
        return max(0.01, float(ui_cfg.get("heartbeat_seconds", 30)))
    except (TypeError, ValueError):
        return 30.0


def _terminal_status(value: Any) -> str:
    return str(value or "unknown").strip().lower()


class StageProgressReporter:
    """Write append-only events and a stage status snapshot for S0-S3."""

    schema_version = "rph_v4_stage_progress_v1"

    def __init__(
        self,
        stage_dir: Path,
        stage_name: str,
        default_fields: Optional[Dict[str, Any]] = None,
        event_callback: Callable[[str, dict[str, Any]], None] | None = None,
        run_id: Optional[str] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        self.stage_dir = Path(stage_dir)
        self.stage_name = stage_name
        self.stage_dir.mkdir(parents=True, exist_ok=True)
        self._events_path = self.stage_dir / "events.jsonl"
        self._status_path = self.stage_dir / "status.json"
        self._default_fields = dict(default_fields or {})
        if run_id is not None:
            self._default_fields[RUN_ID_FIELD] = run_id
        self._event_callback = event_callback
        self._structures: Dict[str, Dict[str, Any]] = {}
        self._batches: Dict[str, Dict[str, Any]] = {}
        self._steps: Dict[str, Dict[str, Any]] = {}
        self._current_step: Optional[str] = None
        self._funnel: Dict[str, Dict[str, Any]] = {}
        self._decisions: list[Dict[str, Any]] = []
        self._science_summary: Dict[str, Any] = {}
        self._alerts: list[Dict[str, Any]] = []
        self._stage_fields: Dict[str, Any] = {}
        self._status = "running"
        self._started_at = _now()
        self._updated_at = self._started_at
        self._heartbeat_at = self._started_at
        self._pid = os.getpid()
        self.heartbeat_interval_seconds = _heartbeat_interval_seconds(config, self._default_fields)
        self._finished_at: Optional[float] = None
        self._lock = threading.RLock()
        self._write_status()

    @property
    def structure_count(self) -> int:
        return len(self._structures)

    def emit(self, event: str, structure_id: Optional[str] = None, **fields: Any) -> None:
        record: Dict[str, Any] = {
            "schema_version": self.schema_version,
            "timestamp": _now(),
            "stage": self.stage_name,
            **self._default_fields,
            "event": event,
        }
        if structure_id:
            record["structure_id"] = structure_id
        record.update(fields)
        with self._lock:
            try:
                with self._events_path.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(record, ensure_ascii=False, default=str) + "\n")
                    handle.flush()
            except OSError as exc:
                logger.warning("Could not write %s progress event %s: %s", self.stage_name, event, exc)
        if self._event_callback is not None:
            try:
                self._event_callback(event, dict(record))
            except Exception as exc:  # pragma: no cover - defensive UI isolation
                logger.warning("Ignoring %s UI callback failure for %s: %s", self.stage_name, event, exc)

    def start_structure(self, structure_id: str, **meta: Any) -> None:
        entry = self._structures.setdefault(
            structure_id,
            {
                "tasks": {},
                "current_task": None,
            },
        )
        entry.update(
            {
                "status": "running",
                "started_at": entry.get("started_at") or _now(),
                "finished_at": None,
                "error": None,
                **meta,
            }
        )
        self._status = "running"
        self.emit("structure_started", structure_id, **meta)
        self._write_status()

    def register_structure(self, structure_id: str, **meta: Any) -> None:
        """Expose a queued structure before a worker starts its QC attempts."""

        entry = self._structures.setdefault(
            structure_id,
            {
                "tasks": {},
                "current_task": None,
                "status": "pending",
                "started_at": None,
                "finished_at": None,
            },
        )
        entry.update(meta)
        entry.setdefault("status", "pending")
        self.emit("structure_queued", structure_id, status=entry["status"], **meta)
        self._write_status()

    def finish_structure(self, structure_id: str, status: str, **results: Any) -> None:
        entry = self._structures.setdefault(
            structure_id,
            {
                "tasks": {},
                "current_task": None,
                "started_at": _now(),
            },
        )
        entry.update(
            {
                "status": status,
                "finished_at": _now(),
                "current_task": None,
                **results,
            }
        )
        self.emit("structure_finished", structure_id, status=status, **results)
        self._write_status()

    def emit_stage_event(self, event: str, **fields: Any) -> None:
        self._stage_fields.update(fields)
        normalized = _terminal_status(fields.get("status"))
        if event.endswith("_started") or event == "stage_started":
            self._status = "running"
            self._finished_at = None
        elif event.endswith("_failed") or normalized == "failed":
            self._status = "failed"
            self._finished_at = _now()
        elif event.endswith("_completed") or event == "stage_finished":
            self._status = str(fields.get("status") or self._stage_terminal_status())
            self._finished_at = _now()
        self.emit(event, None, **fields)
        self._write_status()

    def batch_event(self, event: str, payload: Dict[str, Any]) -> None:
        """Persist a homogeneous batch event and forward it to the UI."""

        fields = dict(payload or {})
        batch_id = str(fields.get("batch") or fields.get("batch_id") or "batch")
        state = dict(self._batches.get(batch_id) or {})
        state.update(fields)
        state["updated_at"] = _now()
        if event == "batch_progress":
            state["last_completion_at"] = state["updated_at"]
        state["label"] = str(fields.get("label") or state.get("label") or batch_id)
        if event == "batch_started":
            state.update({"status": "running", "done": 0, "failed": 0})
        elif event == "batch_finished":
            state["status"] = str(fields.get("status") or "complete")
        else:
            state["status"] = str(fields.get("status") or state.get("status") or "running")
        self._batches[batch_id] = state
        self.emit(event, None, **fields)
        self._write_status()

    def step_event(self, event: str, payload: Dict[str, Any]) -> None:
        fields = dict(payload or {})
        step_id = str(fields.get("step") or fields.get("step_id") or "step")
        variant = str(fields.get("variant") or "")
        key = f"{variant}:{step_id}" if variant else step_id
        state = dict(self._steps.get(key) or {})
        state.update(fields)
        state.update({"id": step_id, "key": key, "updated_at": _now()})
        if event == "step_started":
            state["status"] = "running"
            state["started_at"] = fields.get("started_at") or _now()
            state["finished_at"] = None
            self._current_step = key
        elif event == "step_heartbeat":
            state["status"] = "running"
        elif event == "step_finished":
            state["status"] = str(fields.get("status") or "complete")
            state["finished_at"] = fields.get("finished_at") or _now()
            if self._current_step == key:
                self._current_step = None
        elif event == "step_failed":
            state["status"] = "failed"
            state["finished_at"] = fields.get("finished_at") or _now()
            if self._current_step == key:
                self._current_step = None
        self._steps[key] = state
        fields["status"] = state["status"]
        self.emit(event, None, **fields)
        self._write_status()

    def batch_job_event(self, event: str, payload: Dict[str, Any]) -> None:
        fields = dict(payload or {})
        batch_id = str(fields.get("batch") or "batch")
        job_id = str(fields.get("job_id") or fields.get("current") or "job")
        batch: Dict[str, Any] = dict(self._batches.get(batch_id) or {"label": batch_id})
        active: Dict[str, Any] = dict(batch.get("active_jobs") or {})
        state: Dict[str, Any] = dict(active.get(job_id) or {})
        state.update(fields)
        state.update({"id": job_id, "updated_at": _now()})
        if event == "batch_job_started":
            state["status"] = "running"
            state["started_at"] = fields.get("started_at") or _now()
            active[job_id] = state
        elif event == "batch_job_heartbeat":
            state["status"] = "running"
            active[job_id] = state
        elif event == "batch_job_retry":
            state["status"] = "running"
            state["attempt"] = int(fields.get("attempt") or int(state.get("attempt") or 1) + 1)
            active[job_id] = state
            batch["retried"] = int(batch.get("retried") or 0) + 1
        else:
            state["status"] = "failed" if event == "batch_job_failed" else str(
                fields.get("status") or "complete"
            )
            state["finished_at"] = fields.get("finished_at") or _now()
            active.pop(job_id, None)
            tail: list[Dict[str, Any]] = list(batch.get("completed_job_tail") or [])
            tail.append(state)
            limit = int((self._default_fields.get("ui") or {}).get("completed_job_tail", 5))
            batch["completed_job_tail"] = tail[-max(1, limit):]
        batch["active_jobs"] = active
        batch["running"] = len(active)
        batch["updated_at"] = _now()
        self._batches[batch_id] = batch
        self.emit(event, None, **fields)
        self._write_status()

    def record_decision(self, **fields: Any) -> None:
        record = {"timestamp": _now(), **fields}
        self._decisions.append(record)
        self._decisions = self._decisions[-50:]
        self.emit("decision", None, **record)
        self._write_status()

    def external_event(self, event: str, payload: Dict[str, Any]) -> None:
        """Route stage-runtime UI events into the durable status contract."""

        fields = dict(payload or {})
        with self._lock:
            if event.startswith("batch_job_"):
                self.batch_job_event(event, fields)
            elif event.startswith("batch_"):
                self.batch_event(event, fields)
            elif event.startswith("step_"):
                self.step_event(event, fields)
            elif event == "science_summary":
                self.science_summary(**fields)
            elif event == "funnel_progress":
                variant = str(fields.get("variant") or "default")
                funnel = dict(self._funnel.get(variant) or {})
                funnel[str(fields.get("step") or "unknown")] = fields.get("candidates")
                self._funnel[variant] = funnel
                self.emit("funnel_progress", None, **fields)
                self._write_status()
            elif event == "decision":
                self.record_decision(**fields)
            elif event == "alert":
                message = str(fields.pop("message", "runtime alert"))
                self.alert(message, **fields)
            else:
                self.emit(event, None, **fields)
                self._write_status()

    def science_summary(self, **fields: Any) -> None:
        variant = str(fields.get("variant") or "")
        if variant:
            variants = dict(self._science_summary.get("variants") or {})
            variants[variant] = dict(fields)
            self._science_summary["variants"] = variants
        self._science_summary.update(fields)
        self.emit("science_summary", None, **fields)
        self._write_status()

    def alert(self, message: str, **fields: Any) -> None:
        record = {"message": message, **fields}
        self._alerts.append(record)
        limit = int((self._default_fields.get("ui") or {}).get("recent_alert_limit", 50))
        self._alerts = self._alerts[-max(1, limit):]
        self.emit("alert", None, **record)
        self._write_status()

    def calculator_event(self, event: str, payload: Dict[str, Any]) -> None:
        structure_id = str(payload.get("structure_id", "")).strip()
        event_payload = dict(payload)
        event_payload.pop("structure_id", None)
        if not structure_id:
            self.emit(event, None, **event_payload)
            return
        entry = self._structures.setdefault(
            structure_id,
            {
                "status": "running",
                "started_at": _now(),
                "tasks": {},
                "current_task": None,
            },
        )
        task, _, action = event.rpartition("_")
        if action not in {"started", "finished", "skipped"}:
            self.emit(event, structure_id, **event_payload)
            self._write_status()
            return
        task_state = dict((entry.get("tasks") or {}).get(task, {}))
        task_state.update(
            {
                "status": "running" if action == "started" else str(event_payload.get("status", "skipped")),
                "engine": event_payload.get("engine"),
                "method": event_payload.get("method"),
                "basis": event_payload.get("basis"),
                "solvent": event_payload.get("solvent"),
                "solvent_model": event_payload.get("solvent_model"),
                "output": event_payload.get("output"),
                "energy_hartree": event_payload.get("energy_hartree"),
                "error": event_payload.get("error"),
            }
        )
        if action == "started":
            task_state["started_at"] = task_state.get("started_at") or _now()
            task_state["finished_at"] = None
        else:
            task_state["finished_at"] = _now()
        entry["tasks"][task] = task_state
        entry["current_task"] = task if action == "started" else None
        if action == "started":
            entry["status"] = "running"
        self.emit(event, structure_id, **event_payload)
        self._write_status()

    def report_structure_subprocess(
        self,
        structure_id: str,
        *,
        pid: int | None,
        sandbox_path: str | None,
    ) -> None:
        """Record the active QC child process for stale-run recovery and UI."""

        entry = self._structures.setdefault(
            structure_id,
            {
                "tasks": {},
                "current_task": None,
                "status": "running",
                "started_at": _now(),
                "finished_at": None,
            },
        )
        entry.update({"pid": pid, "sandbox_path": sandbox_path})
        self.emit(
            "structure_subprocess",
            structure_id,
            pid=pid,
            sandbox_path=sandbox_path,
        )
        self._write_status()

    def touch_structure_heartbeat(self, structure_id: str) -> None:
        """Refresh the structure heartbeat while a QC child process is active."""

        entry = self._structures.get(structure_id)
        if entry is None:
            return
        entry["heartbeat_at"] = _now()
        self.emit("structure_heartbeat", structure_id, heartbeat_at=entry["heartbeat_at"])
        self._write_status()

    def reinitialize_after_archive(self) -> None:
        """Restore this reporter's live artifacts after a stage-output archive."""

        self.stage_dir.mkdir(parents=True, exist_ok=True)
        self._events_path = self.stage_dir / "events.jsonl"
        self._status_path = self.stage_dir / "status.json"
        self.emit("stage_progress_reinitialized")
        self._write_status()

    def _stage_terminal_status(self) -> str:
        if any(_terminal_status(row.get("status")) in {"failed", "degraded", "error"} for row in self._structures.values()):
            return "completed_with_failures"
        return "completed"

    def _write_status(self) -> None:
        with self._lock:
            self._write_status_unlocked()

    def _write_status_unlocked(self) -> None:
        self._heartbeat_at = _now()
        self._updated_at = self._heartbeat_at
        completed = 0
        failed = 0
        for row in self._structures.values():
            normalized = _terminal_status(row.get("status"))
            if normalized not in {"", "none", "running"}:
                completed += 1
            if normalized in {"failed", "degraded", "error"}:
                failed += 1
        payload = {
            "schema_version": self.schema_version,
            "stage": self.stage_name,
            **self._default_fields,
            "status": self._status,
            "pid": self._pid,
            "heartbeat_at": self._heartbeat_at,
            "heartbeat_interval_seconds": self.heartbeat_interval_seconds,
            "started_at": self._started_at,
            "updated_at": self._updated_at,
            "finished_at": self._finished_at,
            "total_structures": len(self._structures),
            "completed": completed,
            "failed": failed,
            "structures": dict(self._structures),
            "batches": dict(self._batches),
            "current_step": self._current_step,
            "steps": dict(self._steps),
            "funnel": dict(self._funnel),
            "decisions": list(self._decisions),
            "timing": {
                "started_at": self._started_at,
                "updated_at": self._updated_at,
                "elapsed_seconds": self._updated_at - self._started_at,
            },
            "science_summary": dict(self._science_summary),
            "alerts": list(self._alerts),
            "stage_meta": dict(self._stage_fields),
        }
        try:
            temporary = self._status_path.with_suffix(".json.tmp")
            temporary.write_text(
                json.dumps(payload, indent=2, ensure_ascii=False, default=str),
                encoding="utf-8",
            )
            temporary.replace(self._status_path)
        except OSError as exc:
            logger.warning("Could not write %s status snapshot: %s", self.stage_name, exc)


def scan_all_stages(run_dir: Path) -> Dict[str, Any]:
    run_dir = Path(run_dir)
    stages: Dict[str, Any] = {}
    for stage_name, dirname in [
        ("S0", "S0_Mechanism"),
        ("S1", "S1_ConfSearch"),
        ("S2", "S2_PEB"),
        ("S3", "S3_LowLevel"),
        ("S4", "S4_HighLevel"),
    ]:
        status_path = run_dir / dirname / "status.json"
        if not status_path.exists():
            continue
        try:
            stages[stage_name] = json.loads(status_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            logger.warning("Could not read %s progress snapshot at %s", stage_name, status_path)
    return stages
