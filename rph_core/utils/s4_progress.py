"""Durable, dependency-free progress artifacts for the V4 S4 stage.

The S4 jobs can run for days.  Console output alone is not a reliable status
contract in WSL, especially when a terminal disconnects.  This module writes
an append-only event stream and an atomically replaced status snapshot next to
the S4 calculation artifacts.  Observability must never change QC behaviour:
all reporter failures are isolated from the calculation itself.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional


logger = logging.getLogger(__name__)


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class S4ProgressReporter:
    """Persist human-readable and machine-readable S4 progress artifacts."""

    schema_version = "rph_v4_s4_progress_v1"

    def __init__(
        self,
        output_dir: Path,
        structures: Iterable[Dict[str, Any]],
        event_callback: Callable[[str, dict[str, Any]], None] | None = None,
    ):
        self.output_dir = Path(output_dir)
        self.status_path = self.output_dir / "status.json"
        self.events_path = self.output_dir / "events.jsonl"
        self.log_path = self.output_dir / "s4.log"
        self._event_callback = event_callback
        structure_rows: List[Dict[str, Any]] = []
        for structure in structures:
            structure_rows.append(
                {
                    "id": str(structure["id"]),
                    "kind": str(structure.get("kind", "minimum")),
                    "input_source": "S3_OPT" if structure.get("opt_xyz") else "S3_fallback",
                    "geometry_source": "S3_OPT" if structure.get("opt_xyz") else "S3_input_fallback",
                    "fallback_source": structure.get("fallback_xyz"),
                    "source_s3": dict(structure.get("source_s3") or {}),
                    "status": "pending",
                    "current_task": None,
                    "tasks": {},
                    "usable_for_ml": False,
                    "error": None,
                }
            )
        now = _timestamp()
        self.state: Dict[str, Any] = {
            "schema_version": self.schema_version,
            "stage": "S4",
            "status": "running",
            "started_at": now,
            "updated_at": now,
            "finished_at": None,
            "structures": structure_rows,
            "summary": {},
        }
        self._refresh_summary()
        self._write_snapshot()
        self.emit("stage_started", total_structures=len(structure_rows))

    def start_structure(self, structure_id: str) -> None:
        row = self._structure(structure_id)
        if row is None:
            return
        row["status"] = "running"
        row["error"] = None
        self._touch()
        self.emit("structure_started", structure_id=structure_id, kind=row["kind"])

    def finish_structure(self, payload: Dict[str, Any]) -> None:
        structure_id = str(payload["id"])
        row = self._structure(structure_id)
        if row is None:
            return
        row["status"] = str(payload.get("status", "failed"))
        row["current_task"] = None
        row["usable_for_ml"] = bool(payload.get("usable_for_ml", False))
        row["error"] = payload.get("error")
        row["sp_energy_hartree"] = payload.get("sp_energy_hartree")
        row["ts_frequency_valid"] = payload.get("ts_frequency_valid")
        row["ts_mode_displacement_verified"] = payload.get("ts_mode_displacement_verified")
        row["frequency_count_valid"] = payload.get("frequency_count_valid")
        row["mode_displacement_valid"] = payload.get("mode_displacement_valid")
        row["irc_valid"] = payload.get("irc_valid")
        row["ts_quality_summary"] = payload.get("ts_quality_summary")
        row["imaginary_frequencies_cm1"] = payload.get("imaginary_frequencies_cm1") or payload.get(
            "imaginary_frequencies"
        )
        row["ml_exclusion_reason"] = payload.get("ml_exclusion_reason") or (
            payload.get("error") if not row["usable_for_ml"] else None
        )
        row["geometry_source"] = (
            "S4_OPT" if payload.get("opt_xyz") else row.get("geometry_source")
        )
        self._touch()
        self.emit(
            "structure_finished",
            structure_id=structure_id,
            status=row["status"],
            usable_for_ml=row["usable_for_ml"],
            error=row["error"],
            sp_energy_hartree=row["sp_energy_hartree"],
            ts_frequency_valid=row["ts_frequency_valid"],
            ts_mode_displacement_verified=row["ts_mode_displacement_verified"],
            frequency_count_valid=row["frequency_count_valid"],
            mode_displacement_valid=row["mode_displacement_valid"],
            imaginary_frequencies_cm1=row["imaginary_frequencies_cm1"],
            ts_quality_summary=row["ts_quality_summary"],
            ml_exclusion_reason=row["ml_exclusion_reason"],
            geometry_source=row["geometry_source"],
            fallback_source=row.get("fallback_source"),
        )

    def fail_structure(self, structure_id: str, error: str) -> None:
        row = self._structure(structure_id)
        if row is None:
            return
        row["status"] = "failed"
        row["current_task"] = None
        row["error"] = error
        self._touch()
        self.emit("structure_failed", structure_id=structure_id, error=error)

    def calculator_event(self, event: str, payload: Dict[str, Any]) -> None:
        """Receive lifecycle callbacks from :class:`StageCalculator`."""

        structure_id = str(payload.get("structure_id", ""))
        row = self._structure(structure_id)
        if row is None:
            return
        task, _, action = event.rpartition("_")
        if action not in {"started", "finished", "skipped"}:
            self.emit(event, **payload)
            return
        task_state = dict(row["tasks"].get(task, {}))
        task_state.update(
            {
                "status": "running" if action == "started" else str(payload.get("status", "skipped")),
                "engine": payload.get("engine"),
                "method": payload.get("method"),
                "solvent": payload.get("solvent"),
                "solvent_model": payload.get("solvent_model"),
                "output": payload.get("output"),
                "energy_hartree": payload.get("energy_hartree"),
                "error": payload.get("error"),
            }
        )
        if action == "started":
            task_state["started_at"] = task_state.get("started_at") or _timestamp()
            task_state["finished_at"] = None
        else:
            task_state["finished_at"] = _timestamp()
        row["tasks"][task] = task_state
        row["current_task"] = task if action == "started" else None
        self._touch()
        self.emit(event, **payload)

    def finish_stage(self, manifest_path: Path, results: Iterable[Dict[str, Any]]) -> None:
        materialized = list(results)
        failed = sum(1 for item in materialized if str(item.get("status")) == "failed")
        noncomplete = sum(1 for item in materialized if str(item.get("status")) != "complete")
        self.state["status"] = "completed" if noncomplete == 0 else "completed_with_failures"
        self.state["finished_at"] = _timestamp()
        self.state["manifest"] = str(Path(manifest_path))
        self._touch()
        self.emit(
            "stage_finished",
            status=self.state["status"],
            total_structures=len(materialized),
            failed_structures=failed,
            noncomplete_structures=noncomplete,
            manifest=str(Path(manifest_path)),
        )

    def emit(self, event: str, **payload: Any) -> None:
        """Append an event without allowing logging I/O to abort a QC job."""

        record = {
            "schema_version": self.schema_version,
            "timestamp": _timestamp(),
            "stage": "S4",
            "event": event,
            **payload,
        }
        try:
            with self.events_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(record, ensure_ascii=False, default=str) + "\n")
                handle.flush()
            line = self._human_line(record)
            with self.log_path.open("a", encoding="utf-8") as handle:
                handle.write(line + "\n")
                handle.flush()
            if event.endswith("_failed") or str(payload.get("status") or "").lower() == "failed":
                logger.warning("[S4] %s", line)
            elif event == "stage_finished":
                logger.info("[S4] %s", line)
            else:
                logger.debug("[S4] %s", line)
        except OSError as exc:
            logger.warning("Could not write S4 progress event %s: %s", event, exc)
        if self._event_callback is not None:
            try:
                self._event_callback(event, dict(record))
            except Exception as exc:  # pragma: no cover - defensive UI isolation
                logger.warning("Ignoring S4 UI callback failure for %s: %s", event, exc)

    def _touch(self) -> None:
        self.state["updated_at"] = _timestamp()
        self._refresh_summary()
        self._write_snapshot()

    def _refresh_summary(self) -> None:
        rows = list(self.state["structures"])
        terminal = {"complete", "ts_frequency_unverified", "degraded", "opt_failed_sp_complete", "failed"}
        self.state["summary"] = {
            "total": len(rows),
            "pending": sum(row["status"] == "pending" for row in rows),
            "running": sum(row["status"] == "running" for row in rows),
            "finished": sum(row["status"] in terminal for row in rows),
            "failed": sum(row["status"] == "failed" for row in rows),
            "degraded": sum(row["status"] in {"degraded", "opt_failed_sp_complete", "ts_frequency_unverified"} for row in rows),
            "noncomplete": sum(row["status"] not in {"pending", "running", "complete"} for row in rows),
            "usable_for_ml": sum(bool(row.get("usable_for_ml")) for row in rows),
        }

    def _write_snapshot(self) -> None:
        try:
            temporary = self.status_path.with_suffix(".json.tmp")
            temporary.write_text(
                json.dumps(self.state, indent=2, ensure_ascii=False, default=str),
                encoding="utf-8",
            )
            temporary.replace(self.status_path)
        except OSError as exc:
            logger.warning("Could not write S4 status snapshot: %s", exc)

    def _structure(self, structure_id: str) -> Optional[Dict[str, Any]]:
        for row in self.state["structures"]:
            if row["id"] == structure_id:
                return row
        logger.warning("Ignoring S4 progress event for unknown structure %s", structure_id)
        return None

    @staticmethod
    def _human_line(record: Dict[str, Any]) -> str:
        fields = [record["timestamp"], str(record["event"])]
        for key in ("structure_id", "status", "method", "engine", "error"):
            value = record.get(key)
            if value not in (None, ""):
                fields.append(f"{key}={value}")
        return " | ".join(fields)
