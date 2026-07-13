# pyright: reportMissingImports=false
"""Incremental terminal UI reporter for V4 stage and structure progress."""

from __future__ import annotations

from typing import Any

from rph_core.utils.ui_adapter import UiStructure, adapt_one
from rph_core.utils.ui_state import UiStatus, normalize_status, status_markup, status_style

try:
    from rich.box import SIMPLE as _rich_simple_box
    from rich.markup import escape as _rich_escape
    from rich.panel import Panel as _RichPanel
    from rich.table import Table as _RichTable
    from rich.text import Text as _RichText

    _has_rich = True
except ImportError:  # pragma: no cover - exercised in rich-less environments
    _rich_simple_box = None  # type: ignore[assignment]
    _RichPanel = None  # type: ignore[assignment]
    _RichTable = None  # type: ignore[assignment]
    _RichText = None  # type: ignore[assignment]
    _has_rich = False

    def _escape(value: str) -> str:
        return value

else:
    def _escape(value: str) -> str:
        return _rich_escape(value)


HAS_RICH = _has_rich


_STAGE_TITLES = {
    "S0": "Mechanism",
    "S1": "Conformer Search",
    "S2": "PEB Scan",
    "S3": "Low-level QC",
    "S4": "High-level QC",
}


class RichReporter:
    """Compact incremental UI for V4 stage and structure progress."""

    def __init__(self, console: Any, color: bool = True):
        self.console = console
        self.color = bool(color)
        self._use_rich = bool(
            color
            and HAS_RICH
            and _RichPanel is not None
            and _RichTable is not None
            and _RichText is not None
        )
        self._stage_cache: dict[str, str] = {}
        self._structure_cache: dict[tuple[str, str], str] = {}
        self._stage_meta: dict[str, dict[str, Any]] = {}
        self._stage_structures: dict[str, dict[str, dict[str, Any]]] = {}
        self._stage_headers_printed: set[str] = set()
        self._stage_tables_printed: set[str] = set()

    def event_callback(self, event: str, record: dict[str, Any]) -> None:
        self._handle_event(event, record)

    def s4_event_callback(self, event: str, record: dict[str, Any]) -> None:
        self._handle_event(event, record)

    def stage_started(self, stage: str, meta: dict[str, Any]) -> None:
        fingerprint = self._stage_fingerprint(stage, UiStatus.RUNNING, meta)
        if self._stage_cache.get(stage) == fingerprint and stage in self._stage_headers_printed:
            return
        self._stage_cache[stage] = fingerprint
        self._stage_headers_printed.add(stage)
        title = _STAGE_TITLES.get(stage, stage)
        if self._use_rich:
            assert _RichPanel is not None and _RichText is not None
            header = _RichText()
            header.append(f" {stage} ", style="step.header")
            header.append(f" {title} ", style="step.title")
            reaction_id = meta.get("reaction_id") or meta.get("rx_id")
            if reaction_id:
                header.append(f" {reaction_id} ", style="info")
            header = _RichText.assemble(header, _RichText.from_markup(f" {status_markup(UiStatus.RUNNING)} running"))
            self.console.print(_RichPanel(header, border_style="step.header", padding=(0, 1)))
            return
        reaction_id = meta.get("reaction_id") or meta.get("rx_id")
        suffix = f" {reaction_id}" if reaction_id else ""
        self.console.print(f"[{stage}] {title}{suffix} {self._status_token(UiStatus.RUNNING)} running")

    def stage_finished(self, stage: str, status: UiStatus, summary: dict[str, Any]) -> None:
        counts = self._summary_counts(stage)
        total = int(summary.get("total_structures") or summary.get("total_variants") or counts["total"])
        done = counts["done"] if total else 0
        if total and done == 0:
            done = total
        pieces = [f"{self._status_token(status)} {stage} {status.value}"]
        if total:
            pieces.append(f"{done}/{total} done")
        if counts["failed"]:
            pieces.append(f"{counts['failed']} failed")
        if counts["degraded"]:
            pieces.append(f"{counts['degraded']} degraded")
        if counts["cached"]:
            pieces.append(f"{counts['cached']} reused")
        if counts["usable_for_ml"]:
            pieces.append(f"{counts['usable_for_ml']} usable_for_ml")
        message = " | ".join(pieces)
        if self._use_rich:
            assert _RichText is not None
            self.console.print(_RichText.from_markup(message))
            return
        self.console.print(message)

    def structure_started(self, stage: str, structure: UiStructure) -> None:
        self._render_structure(stage, structure)

    def structure_updated(self, stage: str, structure: UiStructure) -> None:
        self._render_structure(stage, structure)

    def structure_finished(self, stage: str, structure: UiStructure) -> None:
        self._render_structure(stage, structure)

    def _handle_event(self, event: str, record: dict[str, Any]) -> None:
        stage = str(record.get("stage") or "").strip().upper()
        if not stage:
            return
        self._stage_meta[stage] = {**self._stage_meta.get(stage, {}), **record}
        if event == "stage_started" or (event.startswith(stage.lower()) and event.endswith("_started")):
            self.stage_started(stage, record)
            return
        if event == "stage_finished" or (
            event.startswith(stage.lower())
            and (event.endswith("_completed") or event.endswith("_failed"))
        ):
            self.stage_finished(stage, self._derive_stage_status(event, record, stage), record)
            return
        structure = self._apply_structure_event(stage, event, record)
        if structure is None:
            return
        if event == "structure_started":
            self.structure_started(stage, structure)
        elif event in {"structure_finished", "structure_failed"}:
            self.structure_finished(stage, structure)
        else:
            self.structure_updated(stage, structure)

    def _apply_structure_event(self, stage: str, event: str, record: dict[str, Any]) -> UiStructure | None:
        structure_id = str(record.get("structure_id") or record.get("id") or "").strip()
        if not structure_id:
            return None
        stage_rows = self._stage_structures.setdefault(stage, {})
        row: dict[str, Any] = dict(stage_rows.get(structure_id) or {"tasks": {}, "current_task": None})
        payload = {
            key: value
            for key, value in record.items()
            if key not in {"schema_version", "timestamp", "stage", "event", "structure_id"}
        }
        if event == "structure_started":
            row.update(payload)
            row["status"] = "running"
            row["error"] = None
        elif event in {"structure_finished", "structure_failed"}:
            row.update(payload)
            row["status"] = payload.get("status", "failed" if event == "structure_failed" else row.get("status", "complete"))
            row["current_task"] = None
        else:
            task_name, _, action = event.rpartition("_")
            if action not in {"started", "finished", "skipped"}:
                row.update(payload)
            else:
                tasks = dict(row.get("tasks") or {})
                task_state = dict(tasks.get(task_name, {}))
                task_state.update(
                    {
                        "status": "running" if action == "started" else str(payload.get("status", "skipped")),
                        "engine": payload.get("engine"),
                        "method": payload.get("method"),
                        "basis": payload.get("basis"),
                        "solvent": payload.get("solvent"),
                        "solvent_model": payload.get("solvent_model"),
                        "output": payload.get("output"),
                        "energy_hartree": payload.get("energy_hartree"),
                        "error": payload.get("error"),
                    }
                )
                tasks[task_name] = task_state
                row["tasks"] = tasks
                row["current_task"] = task_name if action == "started" else None
                if action == "started":
                    row["status"] = "running"
                if task_name == "single_point" and payload.get("energy_hartree") is not None:
                    row["sp_energy_hartree"] = payload.get("energy_hartree")
        stage_rows[structure_id] = row
        return adapt_one(structure_id, row)

    def _render_structure(self, stage: str, structure: UiStructure) -> None:
        row = dict((self._stage_structures.get(stage) or {}).get(structure.id, {}))
        fingerprint = self._fingerprint(stage, structure, row)
        key = (stage, structure.id)
        if self._structure_cache.get(key) == fingerprint:
            return
        self._structure_cache[key] = fingerprint
        if self._use_rich:
            self.console.print(self._structure_table(stage, structure, row))
        else:
            self.console.print(self._plain_structure_line(stage, structure, row))
        if structure.error and normalize_status(row.get("status")) in {UiStatus.FAILED, UiStatus.DEGRADED}:
            note = f"    error: {structure.error}"
            if self._use_rich:
                assert _RichText is not None
                self.console.print(_RichText(note, style="warning"))
                return
            self.console.print(note)

    def _structure_table(self, stage: str, structure: UiStructure, row: dict[str, Any]) -> Any:
        assert _RichTable is not None
        show_header = stage not in self._stage_tables_printed
        self._stage_tables_printed.add(stage)
        table = _RichTable(box=_rich_simple_box, show_header=show_header, expand=False)
        if stage == "S1":
            table.add_column("variant", style="step.title")
            table.add_column("smiles")
            table.add_column("status")
            table.add_column("detail")
            table.add_row(
                _escape(structure.id),
                _escape(str(row.get("smiles") or "-")),
                self._status_text(structure.status),
                _escape(self._detail_text(stage, structure, row)),
            )
            return table
        if stage == "S2":
            table.add_column("variant", style="step.title")
            table.add_column("status")
            table.add_column("detail")
            table.add_row(
                _escape(structure.id),
                self._status_text(structure.status),
                _escape(self._detail_text(stage, structure, row)),
            )
            return table
        table.add_column("id", style="step.title")
        table.add_column("kind")
        table.add_column("source")
        table.add_column("status")
        table.add_column("active")
        table.add_column("opt")
        table.add_column("freq")
        table.add_column("sp")
        table.add_column("energy", justify="right", style="energy")
        table.add_row(
            _escape(structure.id),
            _escape(structure.kind),
            _escape(structure.source),
            self._status_text(structure.status),
            _escape(structure.current_task or "-"),
            self._task_text(structure, "optimization"),
            self._task_text(structure, "frequency"),
            self._task_text(structure, "single_point"),
            _escape(self._energy_text(structure.energy_hartree)),
        )
        return table

    def _plain_structure_line(self, stage: str, structure: UiStructure, row: dict[str, Any]) -> str:
        if stage in {"S1", "S2"}:
            return " | ".join(
                [
                    stage,
                    structure.id,
                    f"{self._plain_status(structure.status)} {structure.status.value}",
                    self._detail_text(stage, structure, row),
                ]
            )
        task_summary = ", ".join(
            [
                f"opt={self._plain_task(structure, 'optimization')}",
                f"freq={self._plain_task(structure, 'frequency')}",
                f"sp={self._plain_task(structure, 'single_point')}",
            ]
        )
        return " | ".join(
            [
                stage,
                structure.id,
                structure.kind,
                structure.source,
                f"{self._plain_status(structure.status)} {structure.status.value}",
                f"task={structure.current_task or '-'}",
                task_summary,
                f"energy={self._energy_text(structure.energy_hartree)}",
            ]
        )

    def _detail_text(self, stage: str, structure: UiStructure, row: dict[str, Any]) -> str:
        pieces: list[str] = []
        if row.get("selected"):
            pieces.append(f"selected={row['selected']}")
        if row.get("selected_id"):
            pieces.append(f"selected={row['selected_id']}")
        if row.get("stage_status"):
            pieces.append(f"stage={row['stage_status']}")
        if row.get("confidence"):
            pieces.append(f"confidence={row['confidence']}")
        if row.get("reused"):
            pieces.append("reused")
        if structure.current_task:
            pieces.append(structure.current_task)
        if row.get("degraded_reasons"):
            reasons = row.get("degraded_reasons") or []
            pieces.append(f"reasons={','.join(str(item) for item in reasons)}")
        if structure.error:
            pieces.append(f"error={structure.error}")
        return " | ".join(pieces) if pieces else "-"

    def _summary_counts(self, stage: str) -> dict[str, int]:
        rows = list((self._stage_structures.get(stage) or {}).values())
        return {
            "total": len(rows),
            "done": sum(normalize_status(row.get("status")) not in {UiStatus.PENDING, UiStatus.RUNNING} for row in rows),
            "failed": sum(normalize_status(row.get("status")) == UiStatus.FAILED for row in rows),
            "degraded": sum(normalize_status(row.get("status")) == UiStatus.DEGRADED for row in rows),
            "cached": sum(bool(row.get("reused")) or normalize_status(row.get("status")) == UiStatus.CACHED for row in rows),
            "usable_for_ml": sum(bool(row.get("usable_for_ml")) for row in rows),
        }

    def _derive_stage_status(self, event: str, record: dict[str, Any], stage: str) -> UiStatus:
        if record.get("status"):
            return normalize_status(record.get("status"))
        if event.endswith("_failed"):
            return UiStatus.FAILED
        if event in {"stage_finished"} or event.endswith("_completed"):
            counts = self._summary_counts(stage)
            return UiStatus.COMPLETE if counts["total"] or stage == "S0" else UiStatus.PENDING
        return UiStatus.RUNNING

    def _fingerprint(self, stage: str, structure: UiStructure, row: dict[str, Any]) -> str:
        task_parts = "|".join(
            f"{name}:{task.status.value}:{task.energy_hartree}:{task.error}"
            for name, task in sorted(structure.tasks.items())
        )
        extras = [
            str(row.get("selected") or ""),
            str(row.get("selected_id") or ""),
            str(row.get("smiles") or ""),
            str(row.get("confidence") or ""),
            str(row.get("stage_status") or ""),
            str(row.get("reused") or ""),
            str(row.get("degraded_reasons") or ""),
        ]
        return "::".join(
            [
                stage,
                structure.id,
                structure.status.value,
                structure.current_task or "",
                task_parts,
                structure.error or "",
                str(structure.energy_hartree),
                *extras,
            ]
        )

    def _stage_fingerprint(self, stage: str, status: UiStatus, meta: dict[str, Any]) -> str:
        return "::".join(
            [
                stage,
                status.value,
                str(meta.get("reaction_id") or meta.get("rx_id") or ""),
                str(meta.get("total_structures") or meta.get("total_variants") or ""),
                str(meta.get("condition_signature") or ""),
            ]
        )

    def _status_text(self, status: UiStatus) -> Any:
        if self._use_rich:
            assert _RichText is not None
            return _RichText.from_markup(f"{status_markup(status)} {_escape(status.value)}")
        return f"{self._plain_status(status)} {status.value}"

    def _task_text(self, structure: UiStructure, task_name: str) -> Any:
        task = structure.tasks.get(task_name)
        if task is None:
            return "-"
        if self._use_rich:
            assert _RichText is not None
            return _RichText.from_markup(status_markup(task.status))
        return self._plain_status(task.status)

    def _plain_task(self, structure: UiStructure, task_name: str) -> str:
        task = structure.tasks.get(task_name)
        return self._plain_status(task.status) if task is not None else "-"

    @staticmethod
    def _energy_text(value: float | None) -> str:
        return f"{value:.6f}" if value is not None else "-"

    @staticmethod
    def _plain_status(status: UiStatus) -> str:
        return status_style(status)[0]

    def _status_token(self, status: UiStatus) -> str:
        return status_markup(status) if self._use_rich else self._plain_status(status)
