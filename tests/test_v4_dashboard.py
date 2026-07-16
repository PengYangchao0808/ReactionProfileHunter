"""Tests for the embedded, non-full-screen V4 dashboard."""

from __future__ import annotations

import threading
import time
from pathlib import Path

from rich.console import Console

from rph_core.utils.dashboard_renderers import render_dashboard
from rph_core.utils.dashboard_state import DashboardStateReducer
from rph_core.utils.live_dashboard import LiveDashboardManager
from rph_core.utils.shared_console import RPH_THEME
from rph_core.utils.ui_reporter import RichReporter


def _render_text(reducer: DashboardStateReducer, width: int = 120) -> str:
    console = Console(
        width=width,
        record=True,
        force_terminal=False,
        color_system=None,
        theme=RPH_THEME,
    )
    console.print(render_dashboard(reducer.snapshot(), width=width))
    return console.export_text()


def test_s3_dashboard_coalesces_task_events_into_one_structure_row() -> None:
    now = time.time()
    reducer = DashboardStateReducer(log_path="run/rph_v4.log")
    reducer.apply("s3_started", {"stage": "S3", "reaction_id": "1", "timestamp": now})
    reducer.apply(
        "structure_started",
        {"stage": "S3", "structure_id": "product_ts", "kind": "ts", "timestamp": now},
    )
    reducer.apply(
        "optimization_started",
        {"stage": "S3", "structure_id": "product_ts", "method": "B97-3c", "timestamp": now},
    )
    reducer.apply(
        "optimization_finished",
        {"stage": "S3", "structure_id": "product_ts", "method": "B97-3c", "status": "complete", "timestamp": now},
    )
    reducer.apply(
        "frequency_started",
        {"stage": "S3", "structure_id": "product_ts", "method": "NumFreq", "timestamp": now},
    )
    text = _render_text(reducer)
    assert text.count("product_ts") == 2  # table plus one Active line
    assert "B97-3c" in text
    assert "NumFreq" in text
    assert "freq" in text


def test_s3_validation_is_scientific_result_not_current_task() -> None:
    reducer = DashboardStateReducer()
    reducer.apply("s3_started", {"stage": "S3"})
    reducer.apply("structure_started", {"stage": "S3", "structure_id": "minimum", "kind": "minimum"})
    reducer.apply(
        "optimization_finished",
        {"stage": "S3", "structure_id": "minimum", "method": "B97-3c", "status": "complete"},
    )
    reducer.apply("structure_started", {"stage": "S3", "structure_id": "ts", "kind": "ts"})
    reducer.apply(
        "structure_finished",
        {
            "stage": "S3",
            "structure_id": "ts",
            "kind": "ts",
            "status": "ts_frequency_unverified",
            "ts_frequency_valid": False,
            "significant_imaginary_frequencies_cm1": [],
        },
    )
    text = _render_text(reducer)
    assert "minimum" in text
    assert "0 imag" in text
    assert "single_point" not in text


def test_s4_dashboard_surfaces_retained_s3_fallback() -> None:
    reducer = DashboardStateReducer()
    reducer.apply("stage_started", {"stage": "S4"})
    reducer.apply(
        "structure_finished",
        {
            "stage": "S4",
            "structure_id": "product_ts",
            "kind": "ts",
            "status": "failed",
            "geometry_source": "S3_input_fallback",
            "fallback_source": "S3/product_ts.xyz",
            "usable_for_ml": True,
        },
    )
    assert "S3 retained" in _render_text(reducer)


def test_dashboard_reducer_is_safe_for_parallel_job_events() -> None:
    reducer = DashboardStateReducer()
    reducer.apply("s1_started", {"stage": "S1"})

    def work(index: int) -> None:
        job = f"conf_{index:04d}"
        base = {"stage": "S1", "batch": "precursor:b97_sp", "job_id": job}
        reducer.apply("batch_job_started", {**base, "started_at": time.time()})
        reducer.apply("batch_job_finished", {**base, "status": "complete"})

    threads = [threading.Thread(target=work, args=(index,)) for index in range(30)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    batch = reducer.snapshot().stages["S1"]["batches"]["precursor:b97_sp"]
    assert batch["active_jobs"] == {}


def test_dashboard_mode_falls_back_for_non_tty_console() -> None:
    console = Console(force_terminal=False, theme=RPH_THEME)
    reporter = RichReporter(console, color=True, mode="dashboard")
    assert reporter.mode == "compact"
    assert reporter._dashboard is None


def test_narrow_dashboard_omits_wide_validation_column() -> None:
    reducer = DashboardStateReducer()
    reducer.apply("s3_started", {"stage": "S3"})
    reducer.apply("structure_started", {"stage": "S3", "structure_id": "product", "kind": "minimum"})
    text = _render_text(reducer, width=90)
    assert "Structure" in text
    assert "Validation" not in text


def test_live_dashboard_lifecycle_is_idempotent(monkeypatch) -> None:
    calls: list[str] = []

    class FakeLive:
        def __init__(self, _renderable, **_kwargs):
            calls.append("init")

        def start(self, refresh=True):
            calls.append("start")

        def update(self, _renderable, refresh=True):
            calls.append("update")

        def stop(self):
            calls.append("stop")

    class FakeConsole:
        width = 120
        is_terminal = True

    monkeypatch.setattr("rph_core.utils.live_dashboard.Live", FakeLive)
    manager = LiveDashboardManager(FakeConsole(), refresh_per_second=10)
    assert manager.start() is True
    manager.event_callback("s1_started", {"stage": "S1"})
    time.sleep(0.12)
    manager.finish("complete")
    manager.close()
    manager.close()
    assert calls.count("start") == 1
    assert calls.count("stop") == 1
    assert "update" in calls


def test_classic_qc_mode_suppresses_successful_task_transition_rows() -> None:
    class CaptureConsole:
        is_terminal = False
        width = 120

        def __init__(self) -> None:
            self.lines: list[str] = []

        def print(self, value="") -> None:
            self.lines.append(str(value))

    console = CaptureConsole()
    reporter = RichReporter(console, color=False, mode="classic")
    reporter.event_callback("s3_started", {"stage": "S3"})
    reporter.event_callback(
        "structure_started",
        {"stage": "S3", "structure_id": "product", "kind": "minimum"},
    )
    before = len(console.lines)
    reporter.event_callback(
        "optimization_started",
        {"stage": "S3", "structure_id": "product", "method": "B97-3c"},
    )
    reporter.event_callback(
        "optimization_finished",
        {"stage": "S3", "structure_id": "product", "method": "B97-3c", "status": "complete"},
    )
    assert len(console.lines) == before
    reporter.event_callback(
        "structure_finished",
        {"stage": "S3", "structure_id": "product", "status": "complete", "usable_for_ml": True},
    )
    assert len(console.lines) > before
