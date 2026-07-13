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
    lines = [_colour("RPH V4 Cross-Stage Overview", "complete", colour)]
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
            lines.append(f"  {stage_name}: {stage_label} | {completed}/{total} done" + (f", {failed} failed" if failed else ""))
        else:
            structures = data.get("structures", {})
            structure_count = len(structures) if isinstance(structures, (dict, list)) else 0
            lines.append(f"  {stage_name}: {stage_label} | {structure_count} structures")
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

    overview_table = Table(box=None, show_header=True, expand=False)
    overview_table.add_column("stage", style="step.title")
    overview_table.add_column("status")
    overview_table.add_column("done")
    overview_table.add_column("failed")
    overview_table.add_column("running")
    for stage_name in ("S0", "S1", "S2", "S3", "S4"):
        data = stages.get(stage_name)
        if data is None:
            overview_table.add_row(stage_name, Text("not started", style="dim"), "-", "-", "-")
            continue
        summary = _stage_summary(data)
        status = normalize_status(data.get("status"))
        overview_table.add_row(
            stage_name,
            Text.from_markup(status_markup(status)),
            f"{summary.get('finished', 0)}/{summary.get('total', 0)}",
            str(summary.get("failed", 0)),
            str(summary.get("running", 0)),
        )
    overview_panel = Panel(overview_table, title="Cross-Stage Overview", border_style="step.header")

    events = _recent_events(_stage_dir(output_dir, active), event_count)
    events_text = Text("\n".join(events) if events else f"no recent events for {active}")
    events_panel = Panel(events_text, title=f"Recent Events ({active})", border_style="info")

    active_data = stages.get(active)
    structures = _adapt_stage_structures(active, active_data or {}) if active_data else []
    structure_table = Table(box=None, show_header=True, expand=False)
    structure_table.add_column("id", style="step.title")
    structure_table.add_column("kind")
    structure_table.add_column("source")
    structure_table.add_column("status")
    structure_table.add_column("active")
    structure_table.add_column("opt")
    structure_table.add_column("freq")
    structure_table.add_column("sp")
    structure_table.add_column("energy", justify="right", style="energy")
    for structure in structures:
        structure_table.add_row(
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
            structure_table.add_row(Text(f"  error: {structure.error}", style="warning"))
    structures_panel = Panel(structure_table, title=f"{active} Structures", border_style="info")

    layout = Layout(name="root")
    layout.split_column(Layout(overview_panel, size=10), Layout(events_panel, size=10), Layout(structures_panel))
    return layout


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
