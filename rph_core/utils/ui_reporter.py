# pyright: reportMissingImports=false
"""Incremental terminal UI reporter for V4 stage and structure progress."""

from __future__ import annotations

import time
from typing import Any

from rph_core.utils.live_dashboard import LiveDashboardManager
from rph_core.utils.ui_adapter import (
    UiActiveJob,
    UiBatch,
    UiStep,
    UiStructure,
    adapt_active_job,
    adapt_batch,
    adapt_one,
    adapt_step,
)
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

    def __init__(
        self,
        console: Any,
        color: bool = True,
        mode: str = "balanced",
        *,
        dashboard_config: dict[str, Any] | None = None,
        log_path: str = "",
    ):
        self.console = console
        self.color = bool(color)
        self.mode = mode if mode in {"compact", "balanced", "verbose", "classic", "dashboard"} else "balanced"
        if self.mode == "classic":
            self.mode = "compact"
        dashboard_cfg = dict(dashboard_config or {})
        dashboard_capable = bool(
            color
            and getattr(console, "is_terminal", color)
            and int(getattr(console, "width", 120) or 120) >= int(dashboard_cfg.get("min_width", 80))
        )
        if self.mode == "dashboard" and not dashboard_capable:
            self.mode = "compact"
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
        self._batch_cache: dict[tuple[str, str], str] = {}
        self._batch_last_done: dict[tuple[str, str], int] = {}
        self._batch_last_failed: dict[tuple[str, str], int] = {}
        self._science_cache: dict[str, str] = {}
        self._step_cache: dict[tuple[str, str], str] = {}
        self._active_jobs: dict[tuple[str, str], UiActiveJob] = {}
        self._funnel: dict[tuple[str, str], dict[str, int]] = {}
        self._stage_started_at: dict[str, float] = {}
        self._dashboard_warning_printed = False
        self._dashboard = (
            LiveDashboardManager(
                console,
                log_path=log_path,
                refresh_per_second=float(dashboard_cfg.get("refresh_per_second", 4.0)),
            )
            if self.mode == "dashboard"
            else None
        )

    def start(self) -> None:
        """Start the optional embedded dashboard."""

        if self._dashboard is not None and not self._dashboard.start():
            self._activate_classic_fallback()

    def finish(self, status: str = "complete") -> None:
        if self._dashboard is not None:
            self._dashboard.finish(status)

    def close(self) -> None:
        if self._dashboard is not None:
            self._dashboard.close()

    def event_callback(self, event: str, record: dict[str, Any]) -> None:
        if self._dashboard is not None and not self._dashboard.failed:
            self._dashboard.event_callback(event, record)
            return
        if self._dashboard is not None and self._dashboard.failed:
            self._activate_classic_fallback()
        self._handle_event(event, record)

    def s4_event_callback(self, event: str, record: dict[str, Any]) -> None:
        self.event_callback(event, record)

    def _activate_classic_fallback(self) -> None:
        if not self._dashboard_warning_printed:
            reason = self._dashboard.failure_message if self._dashboard is not None else "unavailable"
            self.console.print(f"WARNING Dashboard unavailable; switching to classic UI ({reason})")
            self._dashboard_warning_printed = True
        self.mode = "compact"

    def stage_started(self, stage: str, meta: dict[str, Any]) -> None:
        self._stage_started_at.setdefault(stage, float(meta.get("timestamp") or time.time()))
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
            signature = str(meta.get("condition_signature") or "")
            if signature:
                header.append(f" cfg:{signature[:8]} ", style="dim")
            header = _RichText.assemble(header, _RichText.from_markup(f" {status_markup(UiStatus.RUNNING)} running"))
            self.console.print(_RichPanel(header, border_style="step.header", padding=(0, 1)))
            return
        reaction_id = meta.get("reaction_id") or meta.get("rx_id")
        suffix = f" {reaction_id}" if reaction_id else ""
        if len(self._stage_headers_printed) > 1:
            self.console.print("=" * 78)
        signature = str(meta.get("condition_signature") or "")
        resources = dict(meta.get("resources") or {})
        context = []
        if signature:
            context.append(f"cfg={signature[:8]}")
        if resources.get("nproc"):
            context.append(f"cores={resources['nproc']}")
        context_text = f" | {' | '.join(context)}" if context else ""
        self.console.print(
            f"{stage}  {title}{suffix}{context_text}  {self._status_token(UiStatus.RUNNING)} running"
        )
        self.console.print("-" * 78)

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
        started_at = self._stage_started_at.get(stage)
        if started_at is not None:
            pieces.append(f"elapsed {self._duration(time.time() - started_at)}")
        message = " | ".join(pieces)
        if self._use_rich:
            assert _RichText is not None
            self.console.print()
            self.console.print(_RichText.from_markup(message))
            self.console.print(_RichText("─" * 78, style="dim"))
            return
        self.console.print(message)
        self.console.print("=" * 78)

    def structure_started(self, stage: str, structure: UiStructure) -> None:
        if stage == "S1":
            self.console.print(f"\nS1 variant: {structure.id}")
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
        if event.startswith("batch_job_"):
            self._render_batch_job(stage, event, record)
            return
        if event.startswith("batch_"):
            self._render_batch(stage, event, record)
            return
        if event.startswith("step_"):
            self._render_step(stage, event, record)
            return
        if event == "decision":
            self._render_decision(stage, record)
            return
        if event in {"science_summary", "funnel_progress"}:
            self._render_science_summary(stage, record)
            return
        if event == "alert":
            self._render_alert(stage, record)
            return
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
        if (
            stage in {"S3", "S4"}
            and self.mode == "compact"
            and event.rpartition("_")[0] in {"optimization", "frequency", "single_point"}
            and normalize_status(record.get("status")) not in {UiStatus.FAILED, UiStatus.DEGRADED}
        ):
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
        table.add_column("geometry")
        table.add_column("status")
        table.add_column("validation")
        table.add_column("SP / Eh", justify="right", style="energy")
        table.add_column("ML")
        table.add_row(
            _escape(structure.id),
            _escape(structure.geometry_source or structure.source),
            self._status_text(structure.status),
            _escape(self._validation_text(structure, row)),
            _escape(self._energy_text(structure.energy_hartree)),
            "yes" if structure.usable_for_ml else "no",
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
        return " | ".join(
            [
                stage,
                structure.id,
                f"geometry={structure.geometry_source or structure.source}",
                f"{self._plain_status(structure.status)} {structure.status.value}",
                self._validation_text(structure, row),
                f"sp={self._energy_text(structure.energy_hartree)} Eh",
                f"ML={'yes' if structure.usable_for_ml else 'no'}",
            ]
        )

    def _render_batch(self, stage: str, event: str, record: dict[str, Any]) -> None:
        batch_id = str(record.get("batch") or record.get("batch_id") or "batch")
        batch: UiBatch = adapt_batch(batch_id, record)
        key = (stage, batch_id)
        fingerprint = repr((event, batch.done, batch.failed, batch.running, batch.current))
        if self._batch_cache.get(key) == fingerprint:
            return
        # Balanced mode prints start/end, failures and roughly ten snapshots.
        ui_cfg = dict(record.get("ui") or {})
        snapshots = max(1, int(ui_cfg.get("progress_snapshots", 20)))
        stride = max(1, batch.total // snapshots) if batch.total else 1
        previous = self._batch_last_done.get(key, -stride)
        previous_failed = self._batch_last_failed.get(key, 0)
        important = (
            event != "batch_progress"
            or batch.failed > previous_failed
            or batch.done - previous >= stride
        )
        if self.mode == "compact" and event == "batch_progress" and batch.failed == 0:
            important = False
        if self.mode == "verbose":
            important = True
        self._batch_cache[key] = fingerprint
        if not important:
            return
        self._batch_last_done[key] = batch.done
        self._batch_last_failed[key] = batch.failed
        pieces = [f"Progress {batch.done}/{batch.total} complete", f"{batch.running} running", f"{batch.failed} failed"]
        if batch.rate_per_minute is not None:
            pieces.append(f"{batch.rate_per_minute:.1f}/min")
        if batch.eta_seconds is not None and batch.done < batch.total:
            pieces.append(f"ETA={self._duration(batch.eta_seconds)}")
        if batch.current and self.mode == "verbose":
            pieces.append(f"latest={batch.current}")
        self.console.print(f"  {'●' if batch.done < batch.total else '✓'} {stage} {batch.label}")
        self.console.print("      " + " | ".join(pieces))
        cores = record.get("cores_per_job")
        if cores:
            self.console.print(f"      Resources {batch.running or '-'} active × {cores} cores/job")
        active = [job for (job_stage, _), job in self._active_jobs.items() if job_stage == stage and job.batch == batch_id]
        if active:
            active_limit = max(1, int(ui_cfg.get("active_job_limit", 8)))
            self.console.print("      Active " + " | ".join(
                f"{job.id} {self._duration(self._job_elapsed(job))}"
                for job in sorted(active, key=lambda item: item.id)[:active_limit]
            ))

    def _render_step(self, stage: str, event: str, record: dict[str, Any]) -> None:
        step_id = str(record.get("step") or record.get("step_id") or "step")
        step: UiStep = adapt_step(step_id, record)
        key = (stage, f"{step.variant}:{step_id}")
        fingerprint = repr((event, int(step.elapsed_seconds), step.status.value, step.error, step.detail))
        if self._step_cache.get(key) == fingerprint:
            return
        self._step_cache[key] = fingerprint
        if event == "step_heartbeat" and self.mode == "compact":
            return
        token = "●" if event in {"step_started", "step_heartbeat"} else "✗" if event == "step_failed" else "✓"
        index = f"[{step.index}/{step.total_steps}]" if step.index else ""
        self.console.print()
        self.console.print(f"  {token} {index} {step.label}")
        if event == "step_started":
            details = []
            if step.purpose:
                details.append(f"purpose={step.purpose}")
            if step.engine or step.method:
                details.append(f"method={'/'.join(filter(None, [step.engine, step.method]))}")
            for name in ("solvent", "search_mode", "energy_window_kcal", "temperature_k", "nprocs"):
                value = step.detail.get(name)
                if value not in (None, ""):
                    details.append(f"{name}={value}")
            if details:
                self.console.print("      " + " | ".join(details))
            return
        self.console.print(f"      elapsed={self._duration(step.elapsed_seconds)} | status={step.status.value}")
        if event == "step_heartbeat":
            activity = []
            if step.detail.get("last_output_age_seconds") is not None:
                activity.append(f"output updated {self._duration(float(step.detail['last_output_age_seconds']))} ago")
            if step.detail.get("output_size_bytes") is not None:
                activity.append(f"output size={self._size_text(int(step.detail['output_size_bytes']))}")
            if activity:
                self.console.print("      " + " | ".join(activity))
            active = [job for (job_stage, _), job in self._active_jobs.items() if job_stage == stage]
            if active:
                self.console.print("      active jobs: " + " | ".join(
                    f"{job.id} {self._duration(self._job_elapsed(job))}"
                    for job in sorted(active, key=lambda item: item.id)[:8]
                ))
            return
        result_names = (
            "atoms", "written", "parsed", "invalid", "input", "unique", "removed",
            "output_candidates", "valid", "failed", "ensemble_members", "representatives",
            "in_window", "selected", "selected_population", "partition_function_relative",
            "conformational_free_energy_correction_kcal", "manifest", "geometry",
        )
        results = [f"{name}={step.detail[name]}" for name in result_names if step.detail.get(name) is not None]
        if results:
            self.console.print("      " + " | ".join(results))
        if step.error:
            self.console.print(f"      error={step.error}")
        ui_cfg = dict(step.detail.get("ui") or {})
        if step.output and (self.mode == "verbose" or bool(ui_cfg.get("show_output_paths", True))):
            self.console.print(f"      output={step.output}")

    def _render_batch_job(self, stage: str, event: str, record: dict[str, Any]) -> None:
        job_id = str(record.get("job_id") or "job")
        key = (stage, job_id)
        existing = self._active_jobs.get(key)
        merged = {**(vars(existing) if existing is not None else {}), **record}
        job = adapt_active_job(job_id, merged)
        if event in {"batch_job_started", "batch_job_heartbeat", "batch_job_retry"}:
            self._active_jobs[key] = job
        else:
            self._active_jobs.pop(key, None)
        if event == "batch_job_retry":
            self.console.print(
                f"  ↻ {job_id} retry {record.get('attempt', '?')} | "
                f"GFN{record.get('gfn_level', '?')} | nprocs={record.get('nprocs', '?')} | "
                f"reason={record.get('reason', '-')}"
            )
        elif event == "batch_job_failed":
            self.console.print(
                f"  ⚠ {job_id} failed | elapsed={self._duration(job.elapsed_seconds)} | "
                f"reason={job.error or '-'} | diagnostics={job.output or '-'}"
            )
        elif event == "batch_job_finished" and self.mode == "verbose":
            energy = f" | energy={job.energy_hartree:.8f} Eh" if job.energy_hartree is not None else ""
            self.console.print(f"  ✓ {job_id} | elapsed={self._duration(job.elapsed_seconds)}{energy}")
        elif event == "batch_job_started" and self.mode == "verbose":
            self.console.print(
                f"  → {job_id} | {job.engine}/{job.method} | nprocs={job.nprocs} | output={job.output}"
            )

    def _render_decision(self, stage: str, record: dict[str, Any]) -> None:
        self.console.print()
        self.console.print(f"  ! {stage} DECISION: {record.get('message') or record.get('decision')}")
        for name in ("previous", "active", "reason"):
            if record.get(name):
                self.console.print(f"      {name}: {record[name]}")

    def _render_science_summary(self, stage: str, record: dict[str, Any]) -> None:
        fingerprint = repr(sorted((key, repr(value)) for key, value in record.items()))
        if self._science_cache.get(stage) == fingerprint:
            return
        self._science_cache[stage] = fingerprint
        if record.get("step"):
            variant = str(record.get("variant") or "default")
            funnel = self._funnel.setdefault((stage, variant), {})
            funnel[str(record["step"])] = int(record.get("candidates") or 0)
            if self.mode != "compact":
                labels = {
                    "crest_generated": "CREST generated",
                    "energy_valid": "energy valid",
                    "torsion_unique": "torsion/RMSD unique",
                    "prefiltered": "prefiltered",
                    "b97_valid": "B97-3c valid",
                    "mrrho_valid": "mRRHO valid",
                }
                self.console.print(
                    "  Funnel " + " → ".join(
                        f"{labels.get(name, name)}={count}" for name, count in funnel.items()
                    )
                )
            return
        if stage == "S1":
            pieces = [
                f"S1 science | ensemble={record.get('ensemble_members', '-')}",
                f"representatives={record.get('representatives', '-')}",
                f"selected={record.get('selected', '-')}",
                f"ranking={record.get('ranking_formula', '-')}",
            ]
            correction = record.get("conformational_free_energy_correction_kcal")
            if correction is not None:
                pieces.append(f"Gconf={float(correction):+.2f} kcal/mol")
            self.console.print(" | ".join(pieces))
            populations = list(record.get("top_population") or [])
            if populations:
                self.console.print("  Top populations")
                for row in populations:
                    population = row.get("population")
                    population_text = f"{float(population) * 100:.2f}%" if population is not None else "-"
                    delta_g = row.get("relative_free_energy_kcal")
                    delta_text = f"{float(delta_g):.2f} kcal/mol" if delta_g is not None else "-"
                    self.console.print(
                        f"      {row.get('id', '-')} | ΔG={delta_text} | population={population_text}"
                    )
            if record.get("thermochemistry_complete") is False:
                self._render_alert(stage, {"message": "thermochemistry incomplete; electronic-only ranking used"})
        elif stage == "S0":
            bonds = ", ".join(
                f"{pair[0]}-{pair[1]}" for pair in record.get("forming_bonds") or [] if len(pair) >= 2
            )
            self.console.print(
                f"S0 science | mechanism={'valid' if record.get('mechanism_valid') else 'unverified'}"
                f" | forming bonds={bonds or '-'} | variants={record.get('variants', '-')}"
            )

    def _render_alert(self, stage: str, record: dict[str, Any]) -> None:
        message = str(record.get("message") or record.get("error") or "runtime alert")
        action = str(record.get("action") or record.get("fallback") or "")
        line = f"! {stage} {message}" + (f" | action={action}" if action else "")
        if self._use_rich:
            assert _RichText is not None
            self.console.print(_RichText(line, style="warning"))
        else:
            self.console.print(line)

    @staticmethod
    def _duration(seconds: float) -> str:
        value = max(0, int(seconds))
        minutes, secs = divmod(value, 60)
        hours, minutes = divmod(minutes, 60)
        return f"{hours}h{minutes:02d}m" if hours else f"{minutes}m{secs:02d}s"

    @staticmethod
    def _size_text(size: int) -> str:
        value = float(max(0, size))
        for unit in ("B", "KB", "MB", "GB"):
            if value < 1024 or unit == "GB":
                return f"{value:.1f} {unit}"
            value /= 1024
        return f"{value:.1f} GB"

    @staticmethod
    def _job_elapsed(job: UiActiveJob) -> float:
        if job.started_at is not None:
            return max(job.elapsed_seconds, time.time() - job.started_at)
        return job.elapsed_seconds

    def _validation_text(self, structure: UiStructure, row: dict[str, Any]) -> str:
        if structure.kind != "ts":
            if structure.current_task:
                return structure.current_task
            if structure.fallback_source:
                return "fallback geometry"
            return "minimum"
        if structure.imaginary_frequencies:
            primary = min(structure.imaginary_frequencies)
            return f"{len(structure.imaginary_frequencies)} imag ({primary:.0f} cm-1)"
        valid = row.get("ts_frequency_valid")
        if valid is True:
            return "TS frequency valid"
        if valid is False:
            return "TS frequency unverified"
        return structure.ts_quality_summary or structure.current_task or "TS pending"

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
        if row.get("forming_bonds"):
            pieces.append(
                "bonds=" + ",".join(
                    f"{pair[0]}-{pair[1]}" for pair in row.get("forming_bonds") or [] if len(pair) >= 2
                )
            )
        if row.get("scan_profile"):
            pieces.append(f"scan={row['scan_profile']}")
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
