from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Callable

import pytest

import rph_core.v4_orchestrator as v4_orchestrator
from rph_core.steps.mechanism_classifier.s0_record import S0ReactionRecord
from rph_core.utils.dashboard_state import DashboardStateReducer
from rph_core.utils.live_dashboard import LiveDashboardManager
from rph_core.utils.run_id import RUN_ID_FIELD
from rph_core.utils.ui_reporter import RichReporter
from rph_core.v4_watch import render_overview


class _FakeConsole:
    is_terminal = False
    width = 120

    def print(self, value: object = "") -> None:
        del value


def _record(tmp_path: Path) -> S0ReactionRecord:
    return S0ReactionRecord(
        rx_id="RXN_UI_RUN_ID",
        product_smiles="C1CCC1",
        reaction_type="4+3",
        mapped_product_smiles="[CH2:1]1[CH2:2][CH2:3][CH2:4]1",
        mapped_forming_bonds=((1, 2), (3, 4)),
        forming_bonds=((0, 1), (2, 3)),
        mapping_confidence=0.99,
        mapping_trusted=True,
        source_csv=tmp_path / "dataset.csv",
        canonical_precursor_smiles="CC",
        mapped_precursor_smiles="[CH3:1][CH3:2]",
        topology="INTER",
        cyclo_mode="UNKNOWN",
        source_row_hash="row-hash",
    )


def _capture_reporter_events(reporter: RichReporter) -> list[tuple[str, dict[str, Any]]]:
    seen: list[tuple[str, dict[str, Any]]] = []

    def handler(event: str, record: dict[str, Any]) -> None:
        seen.append((event, dict(record)))

    callback: Callable[[str, dict[str, Any]], None] = handler
    reporter._handle_event = callback  # noqa: SLF001
    return seen


def test_rich_reporter_ignores_mismatched_run_id(caplog: pytest.LogCaptureFixture) -> None:
    reporter = RichReporter(_FakeConsole(), color=False)
    reporter.set_run_id("expected-run")
    seen = _capture_reporter_events(reporter)

    with caplog.at_level(logging.DEBUG, logger="rph_core.utils.ui_reporter"):
        reporter.event_callback("s1_started", {"stage": "S1", RUN_ID_FIELD: "stale-run"})

    assert seen == []
    assert "Ignoring event from stale run_id stale-run (expected expected-run)" in caplog.text


def test_rich_reporter_accepts_events_when_run_id_is_none() -> None:
    reporter = RichReporter(_FakeConsole(), color=False)
    seen = _capture_reporter_events(reporter)

    reporter.event_callback("s1_started", {"stage": "S1", RUN_ID_FIELD: "other-run"})

    assert seen == [("s1_started", {"stage": "S1", RUN_ID_FIELD: "other-run"})]


def test_rich_reporter_set_run_id_updates_filter() -> None:
    reporter = RichReporter(_FakeConsole(), color=False)
    seen = _capture_reporter_events(reporter)

    reporter.event_callback("s1_started", {"stage": "S1", RUN_ID_FIELD: "legacy-run"})
    reporter.set_run_id("expected-run")
    reporter.event_callback("s1_started", {"stage": "S1", RUN_ID_FIELD: "stale-run"})
    reporter.event_callback("s1_started", {"stage": "S1", RUN_ID_FIELD: "expected-run"})

    assert reporter.run_id == "expected-run"
    assert [record[RUN_ID_FIELD] for _, record in seen] == ["legacy-run", "expected-run"]


def test_live_dashboard_manager_event_callback_filters_by_run_id() -> None:
    manager = LiveDashboardManager(_FakeConsole(), run_id="expected-run")

    manager.event_callback("s1_started", {"stage": "S1", RUN_ID_FIELD: "stale-run"})
    assert manager.reducer.snapshot().stages == {}

    manager.event_callback("s1_started", {"stage": "S1", RUN_ID_FIELD: "expected-run"})
    snapshot = manager.reducer.snapshot()
    assert "S1" in snapshot.stages
    assert snapshot.run_id == "expected-run"


def test_dashboard_state_captures_first_run_id() -> None:
    reducer = DashboardStateReducer()

    reducer.apply("s1_started", {"stage": "S1", RUN_ID_FIELD: "first-run"})
    reducer.apply("s2_started", {"stage": "S2", RUN_ID_FIELD: "second-run"})

    snapshot = reducer.snapshot()
    assert snapshot.run_id == "first-run"


def test_render_overview_displays_run_id_header(tmp_path: Path) -> None:
    status_dir = tmp_path / "S1_ConfSearch"
    status_dir.mkdir(parents=True, exist_ok=True)
    (status_dir / "status.json").write_text(
        json.dumps(
            {
                "schema_version": "rph_v4_stage_progress_v1",
                "stage": "S1",
                "status": "running",
                "reaction_id": "RXN_UI_RUN_ID",
                RUN_ID_FIELD: "11111111-1111-4111-8111-111111111111",
                "structures": {},
                "batches": {},
                "steps": {},
                "funnel": {},
                "science_summary": {},
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    overview = render_overview(tmp_path, colour=False)

    assert "RXN: RXN_UI_RUN_ID" in overview
    assert "RUN: 11111111-1111-4111-8111-111111111111" in overview


def test_orchestrator_wires_run_id_to_reporter(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    run_id = "11111111-1111-4111-8111-111111111111"
    reporters: list[FakeReporter] = []

    class FakeReporter:
        def __init__(self, *_args: Any, **_kwargs: Any) -> None:
            self.order: list[str] = []
            self.run_ids: list[str | None] = []
            reporters.append(self)

        def set_run_id(self, value: str | None) -> None:
            self.order.append("set_run_id")
            self.run_ids.append(value)

        def start(self) -> None:
            self.order.append("start")

        def finish(self, status: str = "complete") -> None:
            self.order.append(f"finish:{status}")

        def close(self) -> None:
            self.order.append("close")

        def event_callback(self, event: str, record: dict[str, Any]) -> None:
            del event, record

        def s4_event_callback(self, event: str, record: dict[str, Any]) -> None:
            del event, record

    def fake_write_s0_from_record(self: v4_orchestrator.V4Orchestrator, work_dir: Path, record: S0ReactionRecord) -> Path:
        del record
        path = Path(work_dir) / "S0_Mechanism" / "mechanism.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps({"forming_bonds": [[0, 1], [2, 3]], RUN_ID_FIELD: self.run_id}, indent=2),
            encoding="utf-8",
        )
        return path

    monkeypatch.setattr(v4_orchestrator, "RichReporter", FakeReporter)
    monkeypatch.setattr(v4_orchestrator, "new_run_id", lambda: run_id)
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        fake_write_s0_from_record,
    )

    orchestrator = v4_orchestrator.V4Orchestrator(color=False)
    result = orchestrator.run(_record(tmp_path), tmp_path / "run", stop_after="s0")

    assert result[RUN_ID_FIELD] == run_id
    assert len(reporters) == 1
    assert reporters[0].run_ids == [run_id]
    assert reporters[0].order[:2] == ["set_run_id", "start"]
