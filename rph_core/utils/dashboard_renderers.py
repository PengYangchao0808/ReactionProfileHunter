# pyright: reportMissingImports=false
"""Rich renderers shared by the embedded dashboard and terminal watcher."""

from __future__ import annotations

import time
from typing import Any

from rph_core.utils.dashboard_state import DashboardState
from rph_core.utils.ui_state import UiStatus, normalize_status, status_style

try:
    from rich.console import Group
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text

    HAS_RICH = True
except ImportError:  # pragma: no cover
    Group = Panel = Table = Text = None  # type: ignore[assignment,misc]
    HAS_RICH = False


_STAGE_TITLES = {
    "S0": "Mechanism",
    "S1": "Conformer Search",
    "S2": "PEB Scan",
    "S3": "Low-level QC",
    "S4": "High-level QC",
}


def render_dashboard(state: DashboardState, *, width: int = 120, final: bool = False) -> Any:
    """Build one bounded dashboard renderable from a state snapshot."""

    if not HAS_RICH:
        return _plain_summary(state)
    assert Group is not None and Panel is not None
    active = state.active_stage or _latest_stage(state)
    stage = state.stages.get(active or "", {})
    content = _render_stage(active or "", stage, width=width)
    alerts = _render_alerts(stage)
    parts = [_render_header(state, active, final=final), content]
    if alerts is not None:
        parts.append(alerts)
    footer = Text(style="dim")
    footer.append(f"log: {state.log_path or 'rph_v4.log'}")
    return Panel(
        Group(*parts),
        border_style="success" if final and state.final_status == "complete" else "step.header",
        padding=(0, 1),
    )


def _render_header(state: DashboardState, active: str | None, *, final: bool) -> Any:
    assert Text is not None
    elapsed = max(0.0, time.time() - state.started_at)
    header = Text()
    header.append("RPH V4", style="step.header")
    if state.reaction_id:
        header.append(f" · rx_id={state.reaction_id}", style="step.title")
    if state.condition_signature:
        header.append(f" · cfg:{state.condition_signature[:8]}", style="dim")
    header.append(f" · elapsed {_duration(elapsed)}\n", style="dim")
    for name in ("S0", "S1", "S2", "S3", "S4"):
        raw = (state.stages.get(name) or {}).get("status", "pending")
        ui_status = normalize_status(raw)
        icon, style = status_style(ui_status)
        header.append(f"{name} {icon}", style=style)
        if name == active:
            header.append(f" {_STAGE_TITLES[name]}", style="step.title")
        header.append("   ")
    progress = _stage_progress_brief(state.stages.get(active or "", {}))
    if progress:
        header.append(f"\nCurrent: {active} {progress}", style="info")
    if final:
        header.append(f"  {state.final_status}", style="success" if state.final_status == "complete" else "warning")
    return header


def _render_stage(name: str, stage: dict[str, Any], *, width: int) -> Any:
    if name == "S1":
        return _render_s1(stage, width=width)
    if name == "S2":
        return _render_s2(stage, width=width)
    if name in {"S3", "S4"}:
        return _render_qc(name, stage, width=width)
    return _render_generic(name, stage)


def _render_s1(stage: dict[str, Any], *, width: int) -> Any:
    assert Group is not None and Table is not None and Text is not None
    table = Table.grid(expand=True)
    table.add_column("Variant", ratio=2, style="step.title")
    table.add_column("Step", ratio=5)
    table.add_column("Progress", ratio=3, justify="right")
    structures = dict(stage.get("structures") or {})
    steps = dict(stage.get("steps") or {})
    variants = list(structures)
    for row in steps.values():
        variant = str(row.get("variant") or "")
        if variant and variant not in variants:
            variants.append(variant)
    for variant in variants or ["waiting"]:
        structure = structures.get(variant, {})
        variant_steps = [row for row in steps.values() if str(row.get("variant") or "") == variant]
        current = next((row for row in variant_steps if row.get("status") == "running"), None)
        latest = current or (variant_steps[-1] if variant_steps else None)
        status = normalize_status(structure.get("status"))
        icon, style = status_style(status)
        label = "waiting"
        progress = "—"
        if latest:
            label = f"[{latest.get('index', 0)}/{latest.get('total_steps', 10)}] {latest.get('label') or latest.get('step')}"
        batch = _variant_batch(stage, variant)
        if batch:
            progress = f"{int(batch.get('done') or 0)}/{int(batch.get('total') or 0)}"
        table.add_row(variant, Text(f"{icon} {label}", style=style), progress)

    current_step_key = stage.get("current_step")
    current_step = steps.get(current_step_key, {}) if current_step_key else {}
    batch = _active_batch(stage)
    current_lines = Text()
    if current_step:
        current_lines.append(
            f"[{current_step.get('index', 0)}/{current_step.get('total_steps', 10)}] "
            f"{current_step.get('label') or current_step.get('step')}\n",
            style="step.title",
        )
    if batch:
        done = int(batch.get("done") or 0)
        total = int(batch.get("total") or 0)
        current_lines.append_text(_progress_line(done, total, width=38))
        current_lines.append(
            f" · {len(batch.get('active_jobs') or {})} active · {int(batch.get('failed') or 0)} failed"
        )
        rate = batch.get("rate_per_minute")
        eta = batch.get("eta_seconds")
        if rate:
            current_lines.append(f" · {float(rate):.1f}/min")
        if eta is not None:
            current_lines.append(f" · ETA {_duration(float(eta))}")
        jobs = list((batch.get("active_jobs") or {}).values())[:6]
        if jobs:
            current_lines.append("\nactive: ", style="dim")
            current_lines.append(" · ".join(_job_label(job) for job in jobs), style="dim")
    funnel_text = _funnel_line(stage)
    sections: list[Any] = [Panel(table, title="Variants", border_style="info")]
    if current_lines.plain:
        sections.append(Panel(current_lines, title="Current work", border_style="info"))
    if funnel_text:
        sections.append(Panel(Text(funnel_text), title="Science", border_style="info"))
    return Group(*sections)


def _render_s2(stage: dict[str, Any], *, width: int) -> Any:
    assert Group is not None and Panel is not None and Table is not None and Text is not None
    table = Table.grid(expand=True)
    table.add_column("Variant", ratio=3, style="step.title")
    table.add_column("Scan", ratio=5)
    table.add_column("Status", ratio=2)
    structures = dict(stage.get("structures") or {})
    for structure_id, row in structures.items():
        batch = _variant_batch(stage, structure_id)
        done, total = int(batch.get("done") or 0), int(batch.get("total") or 0)
        status = normalize_status(row.get("status"))
        icon, style = status_style(status)
        table.add_row(structure_id, _progress_line(done, total, width=30), Text(f"{icon} {status.value}", style=style))
    if not structures:
        table.add_row("waiting", Text("○ no scan started", style="dim"), "pending")
    active = _active_batch(stage)
    parts: list[Any] = [Panel(table, title="Variants", border_style="info")]
    if active:
        done = int(active.get("done") or 0)
        total = int(active.get("total") or 0)
        current_lines = Text()
        current_lines.append(str(active.get("label") or "S2 work"), style="step.title")
        current_lines.append("\n")
        current_lines.append_text(_progress_line(done, total, width=38))
        active_jobs = list((active.get("active_jobs") or {}).values())
        current_lines.append(
            f" · {len(active_jobs)} active · {int(active.get('failed') or 0)} failed"
        )
        started = active.get("started_at")
        elapsed = (
            time.time() - float(started)
            if isinstance(started, (int, float))
            else float(active.get("elapsed_seconds") or 0.0)
        )
        current_lines.append(f" · elapsed {_duration(elapsed)}")
        rate = active.get("rate_per_minute")
        eta = active.get("eta_seconds")
        if rate:
            current_lines.append(f" · {float(rate):.1f}/min")
        if eta is not None and done < total:
            current_lines.append(f" · ETA {_duration(float(eta))}")
        if active.get("current"):
            current_lines.append(f"\nlatest: {active['current']}", style="dim")
        if active_jobs:
            current_lines.append("\nactive: ", style="dim")
            current_lines.append(
                " · ".join(_job_label(job) for job in active_jobs[:6]), style="dim"
            )
        parts.append(Panel(current_lines, title="Current work", border_style="info"))
    else:
        current_step_key = stage.get("current_step")
        current_step = (stage.get("steps") or {}).get(current_step_key, {}) if current_step_key else {}
        if current_step:
            parts.append(
                Panel(
                    Text(str(current_step.get("label") or current_step.get("step"))),
                    title="Current work",
                    border_style="info",
                )
            )
    return Group(*parts)


def _render_qc(stage_name: str, stage: dict[str, Any], *, width: int) -> Any:
    assert Group is not None and Panel is not None and Table is not None and Text is not None
    narrow = width < 105
    table = Table(expand=True, box=None, padding=(0, 1))
    table.add_column("Structure", style="step.title", no_wrap=True)
    if not narrow:
        table.add_column("Type", justify="center")
    table.add_column("OPT", no_wrap=True)
    table.add_column("FREQ", no_wrap=True)
    table.add_column("SP", no_wrap=True)
    table.add_column("E/Eh", justify="right", style="energy", no_wrap=True)
    if not narrow:
        table.add_column("Validation", no_wrap=True)
    table.add_column("ML", justify="center")
    structures = dict(stage.get("structures") or {})
    active_lines: list[str] = []
    usable = failed = degraded = terminal = 0
    for structure_id, row in structures.items():
        status = normalize_status(row.get("status"))
        if status in {UiStatus.COMPLETE, UiStatus.CACHED, UiStatus.DEGRADED, UiStatus.FAILED}:
            terminal += 1
        failed += int(status == UiStatus.FAILED)
        degraded += int(status == UiStatus.DEGRADED)
        usable += int(bool(row.get("usable_for_ml")))
        tasks = dict(row.get("tasks") or {})
        energy = row.get("sp_energy_hartree") or row.get("energy_hartree")
        if energy is None:
            energy = (tasks.get("single_point") or {}).get("energy_hartree")
        columns: list[Any] = [structure_id]
        if not narrow:
            columns.append("TS" if str(row.get("kind") or "").lower() == "ts" else "MIN")
        columns.extend([
            _task_badge(tasks.get("optimization"), default_method="OPT"),
            _frequency_badge(row, tasks.get("frequency")),
            _task_badge(tasks.get("single_point"), default_method="SP"),
            f"{float(energy):.6f}" if energy is not None else "—",
        ])
        if not narrow:
            columns.append(_validation(row, tasks))
        columns.append(Text("✓", style="success") if row.get("usable_for_ml") else Text("○", style="dim"))
        table.add_row(*columns)
        current = row.get("current_task")
        if current:
            task = tasks.get(str(current), {})
            started = task.get("started_at") or row.get("started_at")
            elapsed = time.time() - float(started) if isinstance(started, (int, float)) else 0.0
            active_lines.append(
                f"{structure_id} → {_task_method(task, str(current))} · {_duration(elapsed)}"
            )
    if not structures:
        column_count = 6 if narrow else 8
        table.add_row("waiting", *(["—"] * (column_count - 1)))
    parts: list[Any] = [Panel(table, title=f"{stage_name} Structures", border_style="info")]
    if active_lines:
        parts.append(Panel(Text("\n".join(active_lines[:6])), title="Active", border_style="info"))
    total_tasks = sum(len(row.get("tasks") or {}) for row in structures.values())
    completed_tasks = sum(
        1
        for row in structures.values()
        for task in (row.get("tasks") or {}).values()
        if normalize_status(task.get("status")) in {UiStatus.COMPLETE, UiStatus.CACHED}
    )
    summary = Text(
        f"workflow {completed_tasks}/{total_tasks or '?'} · structures {terminal}/{len(structures)}"
        f" · {len(active_lines)} active · {failed} failed · {degraded} degraded · ML usable {usable}/{len(structures)}"
    )
    if stage_name == "S4" and any(_uses_s3_fallback(row) for row in structures.values()):
        summary.append(" · S3 retained", style="warning")
    parts.append(Panel(summary, title="Summary", border_style="info"))
    return Group(*parts)


def _render_generic(name: str, stage: dict[str, Any]) -> Any:
    assert Panel is not None and Text is not None
    status = normalize_status(stage.get("status"))
    icon, style = status_style(status)
    return Panel(Text(f"{icon} {_STAGE_TITLES.get(name, name or 'Pipeline')} · {status.value}", style=style), border_style="info")


def _render_alerts(stage: dict[str, Any]) -> Any | None:
    assert Panel is not None and Text is not None
    rows = list(stage.get("alerts") or [])[-3:] + list(stage.get("decisions") or [])[-3:]
    for structure_id, structure in (stage.get("structures") or {}).items():
        error = structure.get("error")
        output = ""
        if not error:
            for task in (structure.get("tasks") or {}).values():
                if task.get("error"):
                    error = task.get("error")
                    output = task.get("output") or ""
                    break
        if error:
            rows.append({"message": f"{structure_id}: {error}", "output": output})
    if not rows:
        return None
    lines = []
    for row in rows[-3:]:
        prefix = "!" if row.get("message") else "→"
        message = str(row.get("message") or row.get("decision") or row.get("action") or "decision")
        if row.get("output"):
            message += f" · {row['output']}"
        lines.append(f"{prefix} {message}")
    return Panel(Text("\n".join(lines), style="warning"), title="Alerts / decisions", border_style="warning")


def _progress_line(done: int, total: int, *, width: int) -> Any:
    assert Text is not None
    if not total:
        return Text("pending", style="dim")
    bar_width = max(6, min(int(width), 40))
    filled = min(bar_width, max(0, round(bar_width * done / total)))
    line = Text()
    line.append("█" * filled, style="status.complete")
    line.append("░" * (bar_width - filled), style="dim")
    line.append(f" {done}/{total}")
    return line


def _active_batch(stage: dict[str, Any]) -> dict[str, Any]:
    batches = list((stage.get("batches") or {}).values())
    return next((row for row in batches if normalize_status(row.get("status")) == UiStatus.RUNNING), batches[-1] if batches else {})


def _variant_batch(stage: dict[str, Any], variant: str) -> dict[str, Any]:
    matches = [
        row
        for key, row in (stage.get("batches") or {}).items()
        if variant and variant in str(key)
    ]
    return next(
        (
            row
            for row in matches
            if normalize_status(row.get("status")) == UiStatus.RUNNING
        ),
        matches[-1] if matches else {},
    )


def _funnel_line(stage: dict[str, Any]) -> str:
    funnels = dict(stage.get("funnel") or {})
    if not funnels:
        return ""
    variant, data = next(reversed(funnels.items()))
    parts = [f"{value} {str(key).replace('_', ' ')}" for key, value in data.items()]
    return f"{variant}: " + " → ".join(parts)


def _task_method(task: dict[str, Any], fallback: str) -> str:
    return str(task.get("method") or task.get("engine") or fallback).strip()


def _task_badge(task: dict[str, Any] | None, *, default_method: str) -> Any:
    assert Text is not None
    if not task:
        return Text("○", style="dim")
    raw = str(task.get("status") or "pending").lower()
    if raw in {"skipped", "not_requested"}:
        return Text("—", style="dim")
    status = normalize_status(raw)
    icon, style = status_style(status)
    return Text(f"{icon} {_task_method(task, default_method)}", style=style)


def _frequency_badge(row: dict[str, Any], task: dict[str, Any] | None) -> Any:
    assert Text is not None
    if str(row.get("kind") or "minimum").lower() != "ts" and not task:
        return Text("—", style="dim")
    return _task_badge(task, default_method="FREQ")


def _validation(row: dict[str, Any], tasks: dict[str, Any]) -> Any:
    assert Text is not None
    kind = str(row.get("kind") or "minimum").lower()
    if kind != "ts":
        complete = normalize_status((tasks.get("optimization") or {}).get("status")) == UiStatus.COMPLETE
        return Text("✓ minimum", style="success") if complete else Text("pending", style="dim")
    valid = row.get("ts_frequency_valid")
    frequencies = row.get("significant_imaginary_frequencies_cm1") or row.get("imaginary_frequencies_cm1") or []
    if valid is True:
        return Text(f"✓ {len(frequencies) or 1} imag", style="success")
    if valid is False:
        return Text(f"✗ {len(frequencies)} imag", style="error")
    if tasks.get("frequency"):
        status = normalize_status(tasks["frequency"].get("status"))
        if status == UiStatus.RUNNING:
            return Text("● freq", style="status.running")
    return Text("unverified", style="warning")


def _uses_s3_fallback(row: dict[str, Any]) -> bool:
    source = str(row.get("geometry_source") or row.get("input_source") or row.get("fallback_source") or "")
    return "fallback" in source.lower() or bool(row.get("fallback_source"))


def _job_label(job: dict[str, Any]) -> str:
    started = job.get("started_at")
    elapsed = time.time() - float(started) if isinstance(started, (int, float)) else float(job.get("elapsed_seconds") or 0.0)
    return f"{job.get('id') or job.get('job_id') or 'job'} {_duration(elapsed)}"


def _latest_stage(state: DashboardState) -> str | None:
    present = [name for name in ("S0", "S1", "S2", "S3", "S4") if name in state.stages]
    return present[-1] if present else None


def _stage_progress_brief(stage: dict[str, Any]) -> str:
    batches = list((stage.get("batches") or {}).values())
    active = next((row for row in batches if normalize_status(row.get("status")) == UiStatus.RUNNING), None)
    if active is not None:
        done, total = int(active.get("done") or 0), int(active.get("total") or 0)
        percent = f"{100.0 * done / total:.0f}%" if total else "pending"
        return f"{done}/{total} ({percent}) · {len(active.get('active_jobs') or {})} active"
    structures = dict(stage.get("structures") or {})
    terminal = sum(
        normalize_status(row.get("status")) in {
            UiStatus.COMPLETE, UiStatus.FAILED, UiStatus.DEGRADED, UiStatus.CACHED
        }
        for row in structures.values()
    )
    if structures:
        return f"structures {terminal}/{len(structures)}"
    steps = dict(stage.get("steps") or {})
    completed = sum(
        normalize_status(row.get("status")) in {UiStatus.COMPLETE, UiStatus.CACHED}
        for row in steps.values()
    )
    return f"steps {completed}/{len(steps)}" if steps else ""


def _duration(seconds: float) -> str:
    total = max(0, int(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours}h{minutes:02d}m"
    if minutes:
        return f"{minutes}m{secs:02d}s"
    return f"{secs}s"


def _plain_summary(state: DashboardState) -> str:
    active = state.active_stage or _latest_stage(state) or "pipeline"
    stage = state.stages.get(active, {})
    return f"RPH V4 | {active} | {stage.get('status', state.final_status)}"
