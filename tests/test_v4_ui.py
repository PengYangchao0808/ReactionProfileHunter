"""Tests for V4 UI/observability layer.

These tests avoid external dependencies such as rdkit and QC binaries;
they exercise logging, status adapters, the Rich reporter, and rph_watch.
"""

from __future__ import annotations

import json
import logging
import re
import sys
import time
from pathlib import Path
from typing import Any

import pytest

from rph_core.utils.log_manager import setup_v4_logging
from rph_core.utils.ui_adapter import (
    adapt_batch,
    UiStructure,
    adapt_s3_structures,
    adapt_s4_structures,
)
from rph_core.utils.ui_reporter import RichReporter
from rph_core.utils.ui_state import UiStatus, normalize_status, status_markup
from rph_core.utils.shared_console import get_console
from rph_core.utils.stage_progress import StageProgressReporter
from rph_core.v4_watch import main as watch_main, render_overview


def _ansi_escape_count(text: str) -> int:
    return len(re.findall(r"\x1b\[[0-9;]*m", text))


def test_setup_v4_logging_plain_file_has_no_ansi(tmp_path: Path) -> None:
    log_file = tmp_path / "rph_v4.log"
    setup_v4_logging(log_file=log_file, level=logging.INFO, rich_console=False)
    logger = logging.getLogger("rph_core.test_logger")
    logger.info("plain log line")
    content = log_file.read_text(encoding="utf-8")
    assert "plain log line" in content
    assert _ansi_escape_count(content) == 0


def test_setup_v4_logging_no_duplicate_handlers(tmp_path: Path) -> None:
    log_file = tmp_path / "rph_v4.log"
    setup_v4_logging(log_file=log_file, level=logging.INFO, rich_console=False)
    setup_v4_logging(log_file=log_file, level=logging.INFO, rich_console=False)
    root = logging.getLogger()
    file_handlers = [
        h for h in root.handlers if getattr(h, "_rph_v4_log_path", None) == str(log_file.resolve())
    ]
    assert len(file_handlers) == 1


def test_setup_v4_logging_preserves_host_handlers(tmp_path: Path) -> None:
    root = logging.getLogger()
    dummy = logging.StreamHandler(sys.stdout)
    setattr(dummy, "_is_dummy_host_handler", True)
    root.addHandler(dummy)
    try:
        log_file = tmp_path / "rph_v4.log"
        setup_v4_logging(log_file=log_file, level=logging.INFO, rich_console=False)
        host_handlers = [h for h in root.handlers if getattr(h, "_is_dummy_host_handler", False)]
        assert len(host_handlers) == 1
    finally:
        root.removeHandler(dummy)


def test_v4_detail_debug_goes_to_file_not_console(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    log_file = tmp_path / "rph_v4.log"
    setup_v4_logging(log_file=log_file, level=logging.INFO, rich_console=False)
    detail_logger = logging.getLogger("rph_core.utils.resource_utils")
    detail_logger.debug("ORCA maxcore detail sentinel")
    for handler in logging.getLogger().handlers:
        handler.flush()
    assert "ORCA maxcore detail sentinel" in log_file.read_text(encoding="utf-8")
    assert "ORCA maxcore detail sentinel" not in capsys.readouterr().out


def test_ui_state_normalize_status_maps_all_existing_strings() -> None:
    cases = {
        "pending": UiStatus.PENDING,
        "running": UiStatus.RUNNING,
        "complete": UiStatus.COMPLETE,
        "completed": UiStatus.COMPLETE,
        "completed_with_failures": UiStatus.COMPLETE,
        "failed": UiStatus.FAILED,
        "degraded": UiStatus.DEGRADED,
        "opt_failed_sp_complete": UiStatus.DEGRADED,
        "ts_frequency_unverified": UiStatus.DEGRADED,
        "cached": UiStatus.CACHED,
        "reused": UiStatus.CACHED,
        "skipped": UiStatus.CACHED,
    }
    for raw, expected in cases.items():
        assert normalize_status(raw) == expected, raw
    assert normalize_status(None) == UiStatus.PENDING
    assert normalize_status("UNKNOWN") == UiStatus.PENDING


def test_ui_state_status_markup_returns_rich_markup() -> None:
    assert "[cyan]" in status_markup(UiStatus.RUNNING)
    assert "→" in status_markup(UiStatus.RUNNING)


def test_ui_adapter_s3_dict_and_s4_list_produce_same_fields() -> None:
    s3_status = {
        "structures": {
            "product_major": {
                "kind": "minimum",
                "input_source": "selected",
                "status": "complete",
                "current_task": None,
                "tasks": {
                    "optimization": {"status": "complete", "engine": "orca", "method": "B97-3c"},
                    "single_point": {
                        "status": "complete",
                        "engine": "orca",
                        "method": "r2SCAN-3c",
                        "energy_hartree": -312.456,
                    },
                },
                "usable_for_ml": True,
                "error": None,
            }
        }
    }
    s4_status = {
        "structures": [
            {
                "id": "product_major",
                "kind": "minimum",
                "input_source": "selected",
                "status": "complete",
                "current_task": None,
                "tasks": {
                    "optimization": {"status": "complete", "engine": "orca", "method": "B97-3c"},
                    "single_point": {
                        "status": "complete",
                        "engine": "orca",
                        "method": "r2SCAN-3c",
                        "energy_hartree": -312.456,
                    },
                },
                "usable_for_ml": True,
                "error": None,
            }
        ]
    }
    s3_structures = adapt_s3_structures(s3_status)
    s4_structures = adapt_s4_structures(s4_status)
    assert len(s3_structures) == len(s4_structures) == 1
    s3_first: UiStructure = s3_structures[0]
    s4_first: UiStructure = s4_structures[0]
    assert s3_first.id == s4_first.id == "product_major"
    assert s3_first.status == s4_first.status == UiStatus.COMPLETE
    assert s3_first.energy_hartree == s4_first.energy_hartree == pytest.approx(-312.456)
    assert "single_point" in s3_first.tasks
    assert "single_point" in s4_first.tasks


def test_rich_reporter_does_not_crash_on_events() -> None:
    reporter = RichReporter(get_console(), color=False)
    reporter.stage_started("S1", {"reaction_id": "RXN_UI", "total_variants": 1})
    reporter.event_callback(
        "structure_started",
        {
            "stage": "S1",
            "structure_id": "product_major",
            "kind": "minimum",
            "smiles": "CCO",
            "status": "running",
        },
    )
    reporter.event_callback(
        "structure_finished",
        {
            "stage": "S1",
            "structure_id": "product_major",
            "status": "complete",
            "selected": "conf_0001",
        },
    )
    reporter.stage_finished("S1", UiStatus.COMPLETE, {"total_structures": 1})
    assert ("S1", "product_major") in reporter._structure_cache


def test_rich_reporter_dedup_cache_ignores_duplicate_events() -> None:
    reporter = RichReporter(get_console(), color=False)
    reporter.stage_started("S0", {"reaction_id": "RXN_DEDUP"})
    reporter.event_callback(
        "structure_finished",
        {
            "stage": "S0",
            "structure_id": "product_major",
            "status": "complete",
        },
    )
    first_fingerprint = reporter._structure_cache.get(("S0", "product_major"))
    reporter.event_callback(
        "structure_finished",
        {
            "stage": "S0",
            "structure_id": "product_major",
            "status": "complete",
        },
    )
    assert reporter._structure_cache.get(("S0", "product_major")) == first_fingerprint


def test_plain_reporter_uses_stage_boundaries_and_hides_partial_science_summary() -> None:
    class CaptureConsole:
        def __init__(self) -> None:
            self.lines: list[str] = []

        def print(self, value: object = "") -> None:
            self.lines.append(str(value))

    console = CaptureConsole()
    reporter = RichReporter(console, color=False)
    reporter.event_callback(
        "s1_started",
        {
            "stage": "S1",
            "reaction_id": "1",
            "condition_signature": "1234567890abcdef",
            "resources": {"nproc": 16},
        },
    )
    reporter.event_callback(
        "funnel_progress",
        {"stage": "S1", "variant": "precursor", "step": "crest", "candidates": 191},
    )
    reporter.event_callback(
        "batch_started",
        {
            "stage": "S1",
            "batch": "precursor:b97_sp",
            "label": "B97-3c SP",
            "total": 191,
            "running": 2,
        },
    )
    output = "\n".join(console.lines)
    assert "S1  Conformer Search 1 | cfg=12345678 | cores=16" in output
    assert "-" * 78 in output
    assert "● S1 B97-3c SP" in output
    assert "Progress 0/191 complete | 2 running | 0 failed" in output
    assert "ensemble=-" not in output


def test_batch_adapter_and_stage_snapshot_preserve_scientific_progress(tmp_path: Path) -> None:
    reporter = StageProgressReporter(
        tmp_path / "S1_ConfSearch",
        "S1",
        default_fields={"reaction_id": "RXN_UI", "ui": {"stalled_job_warning_seconds": 60}},
    )
    reporter.external_event(
        "step_started",
        {
            "variant": "precursor",
            "step": "b97_sp",
            "index": 7,
            "total_steps": 10,
            "label": "ORCA B97-3c SP ranking",
            "purpose": "rank conformers",
            "started_at": 100.0,
        },
    )
    reporter.external_event(
        "batch_started",
        {"batch": "precursor:b97_sp", "label": "B97-3c SP", "total": 4, "running": 2},
    )
    reporter.external_event(
        "batch_job_started",
        {
            "batch": "precursor:b97_sp",
            "job_id": "conf_0003",
            "engine": "orca",
            "method": "B97-3c",
            "nprocs": 8,
            "started_at": 101.0,
            "output": "ranking/conf_0003",
        },
    )
    reporter.external_event(
        "batch_progress",
        {
            "batch": "precursor:b97_sp",
            "label": "B97-3c SP",
            "total": 4,
            "done": 2,
            "running": 2,
            "failed": 1,
            "rate_per_minute": 1.5,
        },
    )
    reporter.science_summary(
        variant="precursor",
        ensemble_members=4,
        representatives=3,
        scoring_mode="b97_3c_plus_mrrho",
        selected="precursor_conf_0001",
    )
    status = json.loads((tmp_path / "S1_ConfSearch" / "status.json").read_text(encoding="utf-8"))
    batch = adapt_batch("precursor:b97_sp", status["batches"]["precursor:b97_sp"])
    assert batch.done == 2
    assert batch.failed == 1
    assert batch.rate_per_minute == pytest.approx(1.5)
    assert status["science_summary"]["variants"]["precursor"]["ensemble_members"] == 4
    assert status["current_step"] == "precursor:b97_sp"
    assert status["steps"]["precursor:b97_sp"]["index"] == 7
    assert "conf_0003" in status["batches"]["precursor:b97_sp"]["active_jobs"]

    reporter.external_event(
        "batch_job_finished",
        {
            "batch": "precursor:b97_sp",
            "job_id": "conf_0003",
            "status": "complete",
            "elapsed_seconds": 22.0,
            "energy_hartree": -100.25,
        },
    )
    reporter.external_event(
        "step_finished",
        {
            "variant": "precursor",
            "step": "b97_sp",
            "index": 7,
            "total_steps": 10,
            "label": "ORCA B97-3c SP ranking",
            "elapsed_seconds": 22.0,
            "valid": 4,
        },
    )
    status = json.loads((tmp_path / "S1_ConfSearch" / "status.json").read_text(encoding="utf-8"))
    assert status["current_step"] is None
    assert not status["batches"]["precursor:b97_sp"]["active_jobs"]
    assert status["batches"]["precursor:b97_sp"]["completed_job_tail"][-1]["id"] == "conf_0003"


def test_overview_surfaces_mechanism_bonds_and_s1_science(tmp_path: Path) -> None:
    _write_terminal_stage_status(
        tmp_path,
        "S0",
        "S0_Mechanism",
        {
            "stage": "S0",
            "status": "completed",
            "structures": {},
            "science_summary": {"mechanism_valid": True, "forming_bonds": [[3, 17], [8, 22]]},
        },
    )
    _write_terminal_stage_status(
        tmp_path,
        "S1",
        "S1_ConfSearch",
        {
            "stage": "S1",
            "status": "completed",
            "structures": {},
            "science_summary": {
                "ensemble_members": 47,
                "representatives": 6,
                "scoring_mode": "b97_3c_plus_mrrho",
                "selected": "conf_0001",
            },
        },
    )
    overview = render_overview(tmp_path, colour=False)
    assert "bonds 3-17,8-22" in overview
    assert "ensemble 47->6" in overview


def test_overview_surfaces_current_step_funnel_and_active_jobs(tmp_path: Path) -> None:
    now = time.time()
    _write_terminal_stage_status(
        tmp_path,
        "S1",
        "S1_ConfSearch",
        {
            "stage": "S1",
            "status": "running",
            "current_step": "precursor:b97_sp",
            "steps": {
                "precursor:b97_sp": {
                    "variant": "precursor",
                    "label": "ORCA B97-3c SP ranking",
                    "index": 7,
                    "total_steps": 10,
                    "status": "running",
                    "started_at": now - 60,
                    "method": "B97-3c",
                    "engine": "orca",
                }
            },
            "funnel": {"precursor": {"crest_generated": 263, "torsion_unique": 191}},
            "batches": {
                "precursor:b97_sp": {
                    "label": "B97-3c SP",
                    "status": "running",
                    "total": 191,
                    "done": 37,
                    "failed": 0,
                    "active_jobs": {
                        "conf_0038": {
                            "engine": "orca",
                            "method": "B97-3c",
                            "attempt": 1,
                            "started_at": now - 20,
                        }
                    },
                }
            },
            "structures": {},
        },
    )
    overview = render_overview(tmp_path, colour=False)
    assert "step [7/10] ORCA B97-3c SP ranking" in overview
    assert "crest_generated=263 -> torsion_unique=191" in overview
    assert "conf_0038 | orca/B97-3c" in overview


def _write_terminal_stage_status(tmp_path: Path, stage_name: str, stage_dir_name: str, payload: dict[str, Any]) -> None:
    stage_dir = tmp_path / stage_dir_name
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "status.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )


def test_v4_watch_overview_exits_when_all_terminal(tmp_path: Path) -> None:
    completed_payload = {
        "schema_version": "rph_v4_stage_progress_v1",
        "stage": "S0",
        "status": "completed",
        "total_structures": 0,
        "completed": 0,
        "failed": 0,
        "structures": {},
    }
    for stage_name, dir_name in [
        ("S0", "S0_Mechanism"),
        ("S1", "S1_ConfSearch"),
        ("S2", "S2_PEB"),
        ("S3", "S3_LowLevel"),
        ("S4", "S4_HighLevel"),
    ]:
        _write_terminal_stage_status(tmp_path, stage_name, dir_name, completed_payload)

    start = time.time()
    rc = watch_main(["--output", str(tmp_path), "--overview", "--watch", "--interval", "0.01"])
    elapsed = time.time() - start
    assert rc == 0
    assert elapsed < 1.0


def test_v4_watch_s4_detail_exits_when_s4_terminal(tmp_path: Path) -> None:
    s4_status = {
        "schema_version": "rph_v4_s4_progress_v1",
        "stage": "S4",
        "status": "completed",
        "started_at": "2026-07-13T10:00:00+00:00",
        "updated_at": "2026-07-13T10:30:00+00:00",
        "finished_at": "2026-07-13T10:30:00+00:00",
        "structures": [],
        "summary": {
            "total": 0,
            "pending": 0,
            "running": 0,
            "finished": 0,
            "failed": 0,
            "noncomplete": 0,
            "usable_for_ml": 0,
        },
    }
    s4_dir = tmp_path / "S4_HighLevel"
    s4_dir.mkdir(parents=True, exist_ok=True)
    (s4_dir / "status.json").write_text(
        json.dumps(s4_status, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )

    start = time.time()
    rc = watch_main(["--output", str(tmp_path), "--watch", "--interval", "0.01", "--no-color"])
    elapsed = time.time() - start
    assert rc == 0
    assert elapsed < 1.0
