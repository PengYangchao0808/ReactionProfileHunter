"""Terminal viewer for durable V4 progress artifacts.

Run from a second WSL terminal while ``bin/rph_run`` owns the QC calculations:
``python bin/rph_watch --output /tmp/rph4_rx_1 --watch``.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from rph_core.utils.shared_console import get_console
from rph_core.utils.stage_progress import scan_all_stages
from rph_core.utils.run_id import RUN_ID_FIELD
from rph_core.utils.ui_adapter import UiStructure, UiTask, adapt_s3_structures, adapt_s4_structures
from rph_core.utils.ui_state import UiStatus, normalize_status, status_markup

try:
    from rich.layout import Layout
    from rich.live import Live
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text

    _has_rich = True
except ImportError:  # pragma: no cover - exercised in rich-less environments
    _has_rich = False


HAS_RICH = _has_rich


def _read_json(path: Path) -> Optional[Dict[str, Any]]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _fallback_status(stage_dir: Path) -> Dict[str, Any]:
    manifest = _read_json(stage_dir / "manifest.json")
    if manifest is None:
        return {
            "stage": "S4",
            "status": "not_started",
            "summary": {"total": 0, "pending": 0, "running": 0, "finished": 0, "failed": 0, "usable_for_ml": 0},
            "structures": [],
        }
    rows = []
    for item in manifest.get("structures", []):
        rows.append(
            {
                "id": item.get("id", "unknown"),
                "kind": item.get("kind", "minimum"),
                "input_source": item.get("input_source", "unknown"),
                "status": item.get("status", "unknown"),
                "current_task": None,
                "tasks": {},
                "usable_for_ml": bool(item.get("usable_for_ml", False)),
                "error": item.get("error"),
                "sp_energy_hartree": item.get("sp_energy_hartree"),
            }
        )
    return {
        "stage": "S4",
        "status": "manifest_only",
        "summary": {
            "total": len(rows),
            "pending": 0,
            "running": 0,
            "finished": len(rows),
            "failed": sum(row["status"] == "failed" for row in rows),
            "usable_for_ml": sum(row["usable_for_ml"] for row in rows),
        },
        "structures": rows,
    }


def load_status(output_dir: Path) -> Dict[str, Any]:
    stage_dir = Path(output_dir) / "S4_HighLevel"
    return _read_json(stage_dir / "status.json") or _fallback_status(stage_dir)


def _scan_all_stages(run_dir: Path) -> Dict[str, Any]:
    """Load cross-stage status snapshots from S0-S4 when available."""

    return scan_all_stages(run_dir)


def _is_terminal(raw_status: Any) -> bool:
    status = normalize_status(raw_status)
    return status not in {UiStatus.PENDING, UiStatus.RUNNING}


def _all_stages_terminal(stages: Dict[str, Any]) -> bool:
    if not stages:
        return False
    return all(_is_terminal(data.get("status")) for data in stages.values())


def _active_stage(stages: Dict[str, Any]) -> str:
    """Return the first non-terminal stage, or S4 if all are terminal."""

    for stage_name in ("S0", "S1", "S2", "S3", "S4"):
        data = stages.get(stage_name)
        if data is not None and not _is_terminal(data.get("status")):
            return stage_name
    return "S4"


def _stage_dir(output_dir: Path, stage_name: str) -> Path:
    mapping = {
        "S0": "S0_Mechanism",
        "S1": "S1_ConfSearch",
        "S2": "S2_PEB",
        "S3": "S3_LowLevel",
        "S4": "S4_HighLevel",
    }
    return Path(output_dir) / mapping.get(stage_name, stage_name)


def _clip(value: Any, width: int) -> str:
    text = "-" if value in (None, "") else str(value)
    return text if len(text) <= width else text[: max(1, width - 1)] + "~"


def _payload_identity(payload: Dict[str, Any] | None) -> tuple[str, str]:
    if not isinstance(payload, dict):
        return "", ""
    reaction_id = str(payload.get("reaction_id") or payload.get("rx_id") or "")
    run_id = str(payload.get(RUN_ID_FIELD) or "")
    return reaction_id, run_id


def _run_identity(
    output_dir: Path,
    *,
    stages: Dict[str, Any] | None = None,
    status: Dict[str, Any] | None = None,
) -> tuple[str, str]:
    reaction_id, run_id = _payload_identity(status)
    for stage_name in ("S0", "S1", "S2", "S3", "S4"):
        if reaction_id and run_id:
            break
        stage_payload = (stages or {}).get(stage_name)
        stage_reaction_id, stage_run_id = _payload_identity(stage_payload)
        if not reaction_id:
            reaction_id = stage_reaction_id
        if not run_id:
            run_id = stage_run_id
    run_manifest = _read_json(Path(output_dir) / "run.manifest.json")
    manifest_reaction_id, manifest_run_id = _payload_identity(run_manifest)
    if not reaction_id:
        reaction_id = manifest_reaction_id
    if not run_id:
        run_id = manifest_run_id
    return reaction_id, run_id


def _identity_header(reaction_id: str, run_id: str) -> str:
    pieces: list[str] = []
    if reaction_id:
        pieces.append(f"RXN: {reaction_id}")
    if run_id:
        pieces.append(f"RUN: {run_id}")
    return " | ".join(pieces)


def _colour(text: str, status: str, enabled: bool) -> str:
    if not enabled:
        return text
    code = {
        "running": "36",
        "completed": "32",
        "complete": "32",
        "completed_with_failures": "33",
        "failed": "31",
        "degraded": "33",
        "ts_frequency_unverified": "33",
        "not_started": "90",
    }.get(status, "0")
    return f"\033[{code}m{text}\033[0m" if code != "0" else text


def _recent_events(stage_dir: Path, count: int) -> List[str]:
    if count <= 0:
        return []
    try:
        lines = stage_dir.joinpath("events.jsonl").read_text(encoding="utf-8").splitlines()
    except OSError:
        return []
    recent: List[str] = []
    for line in lines[-count:]:
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        fields = [event.get("timestamp", ""), event.get("event", "")]
        if event.get("structure_id"):
            fields.append(f"structure={event['structure_id']}")
        if event.get("status"):
            fields.append(f"status={event['status']}")
        if event.get("error"):
            fields.append(f"error={_clip(event['error'], 72)}")
        recent.append(" | ".join(str(field) for field in fields if field))
    return recent


def _adapt_stage_structures(stage_name: str, raw_status: Dict[str, Any]) -> List[UiStructure]:
    if stage_name in {"S0", "S1", "S2", "S3"}:
        return adapt_s3_structures(raw_status)
    return adapt_s4_structures(raw_status)


def _stage_summary(raw_status: Dict[str, Any]) -> Dict[str, Any]:
    summary = dict(raw_status.get("summary", {}) or {})
    return {
        "total": raw_status.get("total_structures", summary.get("total", 0)),
        "finished": raw_status.get("completed", summary.get("finished", 0)),
        "running": raw_status.get("running", summary.get("running", 0)),
        "failed": raw_status.get("failed", summary.get("failed", 0)),
        "noncomplete": summary.get("noncomplete", 0),
        "usable_for_ml": summary.get("usable_for_ml", 0),
    }


def _render_status_rich(status: Dict[str, Any], output_dir: Path, event_count: int) -> Any:
    assert _has_rich and Layout is not None and Panel is not None and Table is not None and Text is not None

    summary = _stage_summary(status)
    stage_status = normalize_status(status.get("status"))
    header_text = Text.from_markup(f"RPH V4 S4: {status_markup(stage_status)} {stage_status.value}")
    reaction_id, run_id = _run_identity(output_dir, status=status)
    identity = _identity_header(reaction_id, run_id)
    if identity:
        header_text.append(f" · {identity}")
    header = Panel(header_text, border_style="step.header")

    summary_text = Text(
        "structures: {finished}/{total} finished | {running} running | {failed} failed | "
        "{noncomplete} non-complete | {usable} ML-usable".format(
            finished=summary.get("finished", 0),
            total=summary.get("total", 0),
            running=summary.get("running", 0),
            failed=summary.get("failed", 0),
            noncomplete=summary.get("noncomplete", 0),
            usable=summary.get("usable_for_ml", 0),
        )
    )
    summary_panel = Panel(summary_text, title="S4 Summary", border_style="info")

    structures = _adapt_stage_structures("S4", status)
    table = Table(box=None, show_header=True, expand=False)
    table.add_column("id", style="step.title")
    table.add_column("kind")
    table.add_column("source")
    table.add_column("status")
    table.add_column("active")
    table.add_column("opt")
    table.add_column("freq")
    table.add_column("sp")
    table.add_column("energy", justify="right", style="energy")
    for structure in structures:
        table.add_row(
            _clip(structure.id, 30),
            _clip(structure.kind, 8),
            _clip(structure.source, 12),
            Text.from_markup(status_markup(structure.status)),
            _clip(structure.current_task, 17),
            Text.from_markup(status_markup(structure.tasks.get("optimization", UiTask(name="optimization", status=UiStatus.PENDING)).status))
            if structure.tasks.get("optimization")
            else "-",
            Text.from_markup(status_markup(structure.tasks.get("frequency", UiTask(name="frequency", status=UiStatus.PENDING)).status))
            if structure.tasks.get("frequency")
            else "-",
            Text.from_markup(status_markup(structure.tasks.get("single_point", UiTask(name="single_point", status=UiStatus.PENDING)).status))
            if structure.tasks.get("single_point")
            else "-",
            _clip(_energy_text(structure.energy_hartree), 14),
        )
        if structure.error:
            table.add_row(Text(f"  error: {structure.error}", style="warning"))
    structures_panel = Panel(table, title="S4 Structures", border_style="info")

    events = _recent_events(Path(output_dir) / "S4_HighLevel", event_count)
    events_text = Text("\n".join(events) if events else "no recent events")
    events_panel = Panel(events_text, title="Recent Events", border_style="info")

    layout = Layout(name="root")
    layout.split_column(Layout(header, size=3), Layout(summary_panel, size=5), Layout(structures_panel), Layout(events_panel, size=10))
    return layout




def _energy_text(value: Any) -> str:
    if value is None or value == "":
        return "-"
    try:
        return f"{float(value):.6f}"
    except (TypeError, ValueError):
        return "-"


def _render_status_plain(status: Dict[str, Any], output_dir: Path, colour: bool, event_count: int) -> str:
    summary = _stage_summary(status)
    stage_status = str(status.get("status", "unknown"))
    header = _colour(f"RPH V4 S4: {stage_status.upper()}", stage_status, colour)
    reaction_id, run_id = _run_identity(output_dir, status=status)
    identity = _identity_header(reaction_id, run_id)
    if identity:
        header += f" | {identity}"
    lines = [header]
    lines.append(
        "structures: {finished}/{total} finished | {running} running | {failed} failed | "
        "{noncomplete} non-complete | {usable} ML-usable".format(
            finished=summary.get("finished", 0),
            total=summary.get("total", 0),
            running=summary.get("running", 0),
            failed=summary.get("failed", 0),
            noncomplete=summary.get("noncomplete", 0),
            usable=summary.get("usable_for_ml", 0),
        )
    )
    lines.append(f"updated: {status.get('updated_at', '-')}")
    lines.append("")
    lines.append("ID                             KIND     SOURCE       STATUS                     ACTIVE TASK       OPT        FREQ       SP")
    lines.append("-" * 120)
    structures = _adapt_stage_structures("S4", status)
    for structure in structures:
        status_text = _colour(_clip(structure.status.value, 25), stage_status, colour)
        lines.append(
            "{id:<30} {kind:<8} {source:<12} {status:<25} {active:<17} {opt:<10} {freq:<10} {sp:<10}".format(
                id=_clip(structure.id, 30),
                kind=_clip(structure.kind, 8),
                source=_clip(structure.source, 12),
                status=status_text,
                active=_clip(structure.current_task, 17),
                opt=_clip(structure.tasks.get("optimization", UiTask(name="optimization", status=UiStatus.PENDING)).status.value, 10)
                if structure.tasks.get("optimization")
                else "-",
                freq=_clip(structure.tasks.get("frequency", UiTask(name="frequency", status=UiStatus.PENDING)).status.value, 10)
                if structure.tasks.get("frequency")
                else "-",
                sp=_clip(structure.tasks.get("single_point", UiTask(name="single_point", status=UiStatus.PENDING)).status.value, 10)
                if structure.tasks.get("single_point")
                else "-",
            )
        )
        if structure.error:
            lines.append(f"  error: {_clip(structure.error, 104)}")
    events = _recent_events(Path(output_dir) / "S4_HighLevel", event_count)
    if events:
        lines.extend(["", "Recent events:"])
        lines.extend(f"  {event}" for event in events)
    lines.extend(["", f"artifacts: {Path(output_dir) / 'S4_HighLevel'}"])
    return "\n".join(lines)


def _render_overview_plain(output_dir: Path, colour: bool) -> str:
    stages = _scan_all_stages(output_dir)
    header = _colour("RPH V4 Cross-Stage Overview", "complete", colour)
    reaction_id, run_id = _run_identity(output_dir, stages=stages)
    identity = _identity_header(reaction_id, run_id)
    if identity:
        header += f" | {identity}"
    lines = [header]
    lines.append(f"stages with status: {', '.join(sorted(stages)) or 'none'}")
    lines.append("")
    for stage_name in ("S0", "S1", "S2", "S3", "S4"):
        data = stages.get(stage_name)
        if not data:
            lines.append(f"  {stage_name}: not started")
            continue
        stage_status = str(data.get("status", "unknown"))
        summary = _stage_summary(data)
        total = summary.get("total", 0)
        completed = summary.get("finished", 0)
        failed = summary.get("failed", 0)
        stage_label = _colour(stage_status, stage_status, colour)
        if total:
            line = f"  {stage_name}: {stage_label} | {completed}/{total} done"
            if failed:
                line += f", {failed} failed"
            if summary.get("usable_for_ml"):
                line += f", {summary['usable_for_ml']} ML-usable"
            science = _science_brief(stage_name, data)
            batch = _batch_brief(data)
            lines.append(line + (f" | {science}" if science else "") + (f" | {batch}" if batch else ""))
        else:
            structures = data.get("structures", {})
            structure_count = len(structures) if isinstance(structures, (dict, list)) else 0
            science = _science_brief(stage_name, data)
            batch = _batch_brief(data)
            lines.append(
                f"  {stage_name}: {stage_label} | {structure_count} structures"
                + (f" | {science}" if science else "")
                + (f" | {batch}" if batch else "")
            )
    active_stage = _active_stage(stages)
    active_data = stages.get(active_stage) or {}
    activity = _activity_lines(active_stage, active_data)
    if activity:
        lines.extend(["", f"Current activity ({active_stage}):"])
        lines.extend(f"  {line}" for line in activity)
    run_manifest = _read_json(Path(output_dir) / "run.manifest.json")
    if run_manifest:
        lines.append("")
        lines.append(f"reaction: {run_manifest.get('reaction_id', '?')}")
        lines.append(f"condition_signature: {run_manifest.get('condition_signature', '?')}")
        lines.append(f"completed_through: {run_manifest.get('completed_through', '?')}")
        cond = run_manifest.get("conditions", {})
        if cond:
            lines.append(
                f"conditions: solvent={cond.get('solvent', '?')}, "
                f"T={cond.get('temperature_celsius', '?')}C, "
                f"LA={cond.get('has_lewis_acid', '?')}"
            )
            additive_info = cond.get("additive_info")
            if isinstance(additive_info, dict):
                catalyst = additive_info.get("catalyst", "")
                additive = additive_info.get("additive", "")
                if catalyst or additive:
                    lines.append(f"additives: catalyst={catalyst or '-'}, additive={additive or '-'}")
    return "\n".join(lines)


def _render_overview_rich(output_dir: Path, event_count: int) -> Any:
    assert _has_rich and Layout is not None and Panel is not None and Table is not None and Text is not None

    stages = _scan_all_stages(output_dir)
    active = _active_stage(stages)
    reaction_id, run_id = _run_identity(output_dir, stages=stages)
    identity = _identity_header(reaction_id, run_id)

    overview_table = Table(box=None, show_header=True, expand=False)
    overview_table.add_column("stage", style="step.title")
    overview_table.add_column("status")
    overview_table.add_column("done")
    overview_table.add_column("failed")
    overview_table.add_column("running")
    overview_table.add_column("ML")
    overview_table.add_column("scientific result")
    for stage_name in ("S0", "S1", "S2", "S3", "S4"):
        data = stages.get(stage_name)
        if data is None:
            overview_table.add_row(stage_name, Text("not started", style="dim"), "-", "-", "-", "-", "-")
            continue
        summary = _stage_summary(data)
        status = normalize_status(data.get("status"))
        overview_table.add_row(
            stage_name,
            Text.from_markup(status_markup(status)),
            f"{summary.get('finished', 0)}/{summary.get('total', 0)}",
            str(summary.get("failed", 0)),
            str(summary.get("running", 0)),
            str(summary.get("usable_for_ml", 0)),
            _clip("; ".join(filter(None, [_science_brief(stage_name, data), _batch_brief(data)])), 54),
        )
    overview_panel = Panel(
        overview_table,
        title=(f"Cross-Stage Overview · {identity}" if identity else "Cross-Stage Overview"),
        border_style="step.header",
    )

    activity_lines = _activity_lines(active, stages.get(active) or {})
    activity_panel = Panel(
        Text("\n".join(activity_lines) if activity_lines else "no active operation"),
        title=f"Current Activity ({active})",
        border_style="info",
    )

    events = _recent_events(_stage_dir(output_dir, active), event_count)
    events_text = Text("\n".join(events) if events else f"no recent events for {active}")
    events_panel = Panel(events_text, title=f"Recent Events ({active})", border_style="info")

    active_data = stages.get(active)
    structures = _adapt_stage_structures(active, active_data or {}) if active_data else []
    structure_table = Table(box=None, show_header=True, expand=False)
    structure_table.add_column("id", style="step.title")
    structure_table.add_column("geometry")
    structure_table.add_column("status")
    structure_table.add_column("validation")
    structure_table.add_column("SP / Eh", justify="right", style="energy")
    structure_table.add_column("ML")
    for structure in structures:
        structure_table.add_row(
            _clip(structure.id, 30),
            _clip(structure.geometry_source or structure.source, 16),
            Text.from_markup(status_markup(structure.status)),
            _clip(_structure_validation(structure), 28),
            _clip(_energy_text(structure.energy_hartree), 14),
            "yes" if structure.usable_for_ml else "no",
        )
        if structure.error:
            structure_table.add_row(Text(f"  error: {structure.error}", style="warning"))
    structures_panel = Panel(structure_table, title=f"{active} Structures", border_style="info")

    layout = Layout(name="root")
    layout.split_column(
        Layout(overview_panel, size=10),
        Layout(activity_panel, size=max(7, min(15, len(activity_lines) + 2))),
        Layout(events_panel, size=9),
        Layout(structures_panel),
    )
    return layout


def _science_brief(stage: str, data: Dict[str, Any]) -> str:
    science = dict(data.get("science_summary") or {})
    stage_meta = dict(data.get("stage_meta") or {})
    reused = bool(stage_meta.get("reused"))
    if stage == "S0" and science:
        bonds = [f"{pair[0]}-{pair[1]}" for pair in science.get("forming_bonds") or [] if len(pair) >= 2]
        text = "mechanism valid" + (f"; bonds {','.join(bonds)}" if bonds else "")
        return text + ("; reused" if reused else "")
    if stage == "S1" and science:
        ensemble = science.get("ensemble_members")
        representatives = science.get("representatives")
        selected = science.get("selected")
        mode = science.get("scoring_mode")
        if ensemble is not None:
            text = f"ensemble {ensemble}->{representatives}; {mode}; selected {selected}"
            rows = data.get("structures") or {}
            values = rows.values() if isinstance(rows, dict) else rows
            reused_count = sum(bool((row or {}).get("reused")) for row in values)
            return text + (f"; {reused_count} reused" if reused_count else "")
    if stage == "S2":
        rows = data.get("structures") or {}
        values = rows.values() if isinstance(rows, dict) else rows
        confidences: Dict[str, int] = {}
        s2_states: Dict[str, int] = {}
        bonds: list[str] = []
        for row in values:
            s2_state = str((row or {}).get("s2_state") or "").strip()
            if s2_state:
                s2_states[s2_state] = s2_states.get(s2_state, 0) + 1
            confidence = str((row or {}).get("confidence") or "unknown")
            confidences[confidence] = confidences.get(confidence, 0) + 1
            for pair in (row or {}).get("forming_bonds") or []:
                if len(pair) >= 2:
                    bonds.append(f"{pair[0]}-{pair[1]}")
        pieces: list[str] = []
        if s2_states:
            pieces.append("states " + ",".join(f"{name}:{count}" for name, count in sorted(s2_states.items())))
        pieces.extend(f"{name}:{count}" for name, count in sorted(confidences.items()))
        if bonds:
            pieces.append("bonds " + ",".join(dict.fromkeys(bonds)))
        return "; ".join(pieces)
    if stage in {"S3", "S4"}:
        summary = _stage_summary(data)
        return f"{summary.get('usable_for_ml', 0)} ML-usable" + ("; reused" if reused else "")
    return ""


def _structure_validation(structure: Any) -> str:
    if structure.kind != "ts":
        return structure.current_task or ("fallback geometry" if structure.fallback_source else "minimum")
    if structure.imaginary_frequencies:
        primary = min(structure.imaginary_frequencies)
        return f"{len(structure.imaginary_frequencies)} imag ({primary:.0f} cm-1)"
    return structure.ts_quality_summary or structure.current_task or "TS pending"


def _batch_brief(data: Dict[str, Any]) -> str:
    batches = dict(data.get("batches") or {})
    active = [row for row in batches.values() if str((row or {}).get("status")) == "running"]
    if not active:
        return ""
    row = active[-1]
    text = (
        f"{row.get('label', 'batch')} {row.get('done', 0)}/{row.get('total', 0)}"
        f" running={row.get('running', 0)} failed={row.get('failed', 0)}"
    )
    rate = row.get("rate_per_minute")
    if rate is not None:
        text += f" {float(rate):.1f}/min"
    threshold = float((data.get("ui") or {}).get("stalled_job_warning_seconds", 720))
    last = row.get("last_completion_at") or row.get("updated_at")
    if last is not None and time.time() - float(last) > threshold:
        text += " STALLED"
    return text


def _activity_lines(stage: str, data: Dict[str, Any]) -> list[str]:
    lines: list[str] = []
    current_key = data.get("current_step")
    steps = dict(data.get("steps") or {})
    step = dict(steps.get(current_key) or {}) if current_key else {}
    if step:
        started = step.get("started_at")
        elapsed = time.time() - float(started) if started is not None else float(step.get("elapsed_seconds") or 0)
        lines.append(
            f"step [{step.get('index', '?')}/{step.get('total_steps', '?')}] "
            f"{step.get('label', current_key)} | running {int(max(0, elapsed))}s"
        )
        detail = []
        if step.get("method") or step.get("engine"):
            detail.append("/".join(filter(None, [str(step.get("engine") or ""), str(step.get("method") or "")])))
        if step.get("purpose"):
            detail.append(str(step["purpose"]))
        if detail:
            lines.append("  " + " | ".join(detail))
        last_age = step.get("last_output_age_seconds")
        if last_age is not None:
            lines.append(f"  output updated {int(float(last_age))}s ago | {step.get('output', '-')}")
    if stage == "S2":
        selection_bits: list[str] = []
        if step.get("s2_state"):
            selection_bits.append(f"s2_state={step.get('s2_state')}")
            if step.get("seed_evidence") not in (None, ""):
                selection_bits.append(f"seed_evidence={step.get('seed_evidence')}")
            if step.get("selection_source") not in (None, ""):
                selection_bits.append(f"selection_source={step.get('selection_source')}")
        if selection_bits:
            lines.append("selection " + " | ".join(selection_bits))
        else:
            rows = data.get("structures") or {}
            values = rows.items() if isinstance(rows, dict) else [
                (str((row or {}).get("id") or (row or {}).get("structure_id") or "structure"), row)
                for row in rows
            ]
            running_states: list[str] = []
            for structure_id, row in values:
                row_data = dict(row or {})
                if str(row_data.get("status") or "") != "running":
                    continue
                s2_state = str(row_data.get("s2_state") or "").strip()
                if not s2_state:
                    continue
                piece = f"{structure_id}={s2_state}"
                seed_evidence = str(row_data.get("seed_evidence") or "").strip()
                if seed_evidence and seed_evidence != "none":
                    piece += f" ({seed_evidence})"
                running_states.append(piece)
            if running_states:
                lines.append("selection " + " | ".join(running_states[:4]))
    funnel_all = dict(data.get("funnel") or {})
    variant = str(step.get("variant") or "") if step else ""
    funnel = dict(funnel_all.get(variant) or funnel_all.get("default") or {})
    if not funnel and funnel_all:
        funnel = dict(next(reversed(funnel_all.values())))
    if funnel:
        lines.append("funnel " + " -> ".join(f"{name}={count}" for name, count in funnel.items()))
    batches = dict(data.get("batches") or {})
    for batch_id, batch in batches.items():
        if str((batch or {}).get("status")) != "running":
            continue
        lines.append(
            f"batch {(batch or {}).get('label', batch_id)} | "
            f"{(batch or {}).get('done', 0)}/{(batch or {}).get('total', 0)} complete | "
            f"{len((batch or {}).get('active_jobs') or {})} active | "
            f"{(batch or {}).get('failed', 0)} failed"
        )
        active_jobs = dict((batch or {}).get("active_jobs") or {})
        limit = max(1, int((data.get("ui") or {}).get("active_job_limit", 8)))
        for job_id, job in list(sorted(active_jobs.items()))[:limit]:
            started = (job or {}).get("started_at")
            elapsed = time.time() - float(started) if started is not None else float((job or {}).get("elapsed_seconds") or 0)
            lines.append(
                f"  {job_id} | {(job or {}).get('engine', '-')}/{(job or {}).get('method', '-')} | "
                f"attempt {(job or {}).get('attempt', 1)} | {int(max(0, elapsed))}s"
            )
    decisions = list(data.get("decisions") or [])
    if decisions:
        latest = decisions[-1]
        lines.append(f"decision {latest.get('message') or latest.get('decision')}")
    return lines


def _build_renderable(output_dir: Path, overview: bool, colour: bool, event_count: int) -> Any:
    if _has_rich and colour and sys.stdout.isatty():
        if overview:
            return _render_overview_rich(output_dir, event_count)
        return _render_status_rich(load_status(output_dir), output_dir, event_count)
    if overview:
        return _render_overview_plain(output_dir, colour)
    return _render_status_plain(load_status(output_dir), output_dir, colour, event_count)


def render_overview(output_dir: Path, colour: bool = True) -> str:
    return _render_overview_plain(output_dir, colour)


def render_status(status: Dict[str, Any], output_dir: Path, colour: bool = True, event_count: int = 5) -> str:
    return _render_status_plain(status, output_dir, colour, event_count)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Show live RPH V4 stage status in a WSL terminal")
    parser.add_argument("--output", required=True, type=Path, help="RPH run output directory")
    parser.add_argument("--watch", action="store_true", help="Refresh until all stages reach a terminal state")
    parser.add_argument("--interval", type=float, default=3.0, help="Refresh interval in seconds (default: 3)")
    parser.add_argument("--events", type=int, default=5, help="Number of recent events to display")
    parser.add_argument("--no-color", action="store_true", help="Disable colours")
    parser.add_argument("--overview", action="store_true", help="Show cross-stage overview (S0-S4) instead of S4 detail")
    args = parser.parse_args(argv)
    if args.interval <= 0:
        parser.error("--interval must be positive")
    colour = not args.no_color and sys.stdout.isatty()

    console = get_console()
    use_rich = _has_rich and colour and sys.stdout.isatty()

    if not args.watch:
        renderable = _build_renderable(args.output, args.overview, colour, args.events)
        console.print(renderable)
        return 0

    if use_rich:
        assert Live is not None
        with Live(_build_renderable(args.output, args.overview, colour, args.events), console=console, refresh_per_second=1) as live:
            while True:
                time.sleep(args.interval)
                live.update(_build_renderable(args.output, args.overview, colour, args.events))
                stages = _scan_all_stages(args.output)
                if _all_stages_terminal(stages):
                    break
        return 0

    while True:
        console.print("=" * 80)
        renderable = _build_renderable(args.output, args.overview, colour, args.events)
        console.print(renderable)
        stages = _scan_all_stages(args.output)
        if _all_stages_terminal(stages):
            break
        time.sleep(args.interval)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
