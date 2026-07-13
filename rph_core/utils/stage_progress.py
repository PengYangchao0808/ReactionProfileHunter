"""Unified durable progress reporting for V4 S0-S3 stages."""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any, Callable, Dict, Optional

logger = logging.getLogger(__name__)


def _now() -> float:
    return time.time()


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
    ):
        self.stage_dir = Path(stage_dir)
        self.stage_name = stage_name
        self.stage_dir.mkdir(parents=True, exist_ok=True)
        self._events_path = self.stage_dir / "events.jsonl"
        self._status_path = self.stage_dir / "status.json"
        self._default_fields = dict(default_fields or {})
        self._event_callback = event_callback
        self._structures: Dict[str, Dict[str, Any]] = {}
        self._status = "running"
        self._started_at = _now()
        self._updated_at = self._started_at
        self._finished_at: Optional[float] = None
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
                "started_at": entry.get("started_at", _now()),
                "finished_at": None,
                "error": None,
                **meta,
            }
        )
        self._status = "running"
        self.emit("structure_started", structure_id, **meta)
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
        entry["tasks"][task] = task_state
        entry["current_task"] = task if action == "started" else None
        if action == "started":
            entry["status"] = "running"
        self.emit(event, structure_id, **event_payload)
        self._write_status()

    def _stage_terminal_status(self) -> str:
        if any(_terminal_status(row.get("status")) in {"failed", "degraded", "error"} for row in self._structures.values()):
            return "completed_with_failures"
        return "completed"

    def _write_status(self) -> None:
        self._updated_at = _now()
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
            "status": self._status,
            "started_at": self._started_at,
            "updated_at": self._updated_at,
            "finished_at": self._finished_at,
            "total_structures": len(self._structures),
            "completed": completed,
            "failed": failed,
            "structures": dict(self._structures),
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
