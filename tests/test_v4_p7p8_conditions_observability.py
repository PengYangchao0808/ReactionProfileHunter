from __future__ import annotations

import json
from pathlib import Path

import rph_core.v4_orchestrator as v4_orchestrator
from rph_core.steps.mechanism_classifier.s0_record import S0ReactionRecord
from rph_core.utils.stage_progress import StageProgressReporter
from rph_core.v4_watch import _scan_all_stages, render_overview


def _record(*, temp_celsius: str = "25") -> S0ReactionRecord:
    raw_row = {
        "rx_id": "RXN_P7P8",
            "solvent": "acetone",
            "temp_celsius": temp_celsius,
            "has_lewis_acid": "true",
            "catalyst": "LiCl",
            "additive": "LiCl",
    }
    return S0ReactionRecord(
        rx_id="RXN_P7P8",
        product_smiles="C=C",
        reaction_type="test",
        mapped_product_smiles="[CH2:1]=[CH2:2]",
        mapped_forming_bonds=((1, 2), (1, 2)),
        forming_bonds=((0, 1), (0, 1)),
        mapping_confidence=1.0,
        mapping_trusted=True,
        source_csv=Path("dummy.csv"),
        source_row_hash="row_hash",
        raw_row=raw_row,
    )


def _orchestrator() -> v4_orchestrator.V4Orchestrator:
    orchestrator = object.__new__(v4_orchestrator.V4Orchestrator)
    orchestrator.config = {
        "step1": {"protocol": "censo_lite"},
        "step2": {},
        "theory": {
            "s3_low_level": {},
            "s4_high_precision": {},
        },
    }
    return orchestrator


def test_run_manifest_contains_condition_metadata_and_conditions_scope_s0_checkpoint(monkeypatch, tmp_path: Path):
    writes: list[str] = []

    def fake_write_s0(work_dir: Path, record: S0ReactionRecord) -> Path:
        writes.append(record.raw_row["temp_celsius"])
        stage_dir = Path(work_dir) / "S0_Mechanism"
        stage_dir.mkdir(parents=True, exist_ok=True)
        manifest = stage_dir / "mechanism.json"
        manifest.write_text(
            json.dumps(
                {
                    "schema_version": "s0_mechanism_v2",
                    "stage": "S0",
                    "rx_id": record.rx_id,
                    "forming_bonds": [list(pair) for pair in record.forming_bonds],
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        return manifest

    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(fake_write_s0),
    )
    orchestrator = _orchestrator()
    run_dir = tmp_path / "run"

    first_record = _record(temp_celsius="25")
    orchestrator.run(first_record, run_dir, stop_after="s0")

    expected_conditions, expected_signature = orchestrator._condition_metadata(first_record)
    run_manifest = json.loads((run_dir / "run.manifest.json").read_text(encoding="utf-8"))
    assert run_manifest["condition_signature"] == expected_signature
    assert run_manifest["conditions"] == expected_conditions
    assert run_manifest["conditions"]["solvent"] == "acetone"
    assert run_manifest["conditions"]["temperature_celsius"] == "25"
    assert run_manifest["conditions"]["charge"] == 0
    assert run_manifest["conditions"]["multiplicity"] == 1
    assert run_manifest["conditions"]["additive_info"]["catalyst"] == "LiCl"

    second_record = _record(temp_celsius="80")
    orchestrator.run(second_record, run_dir, stop_after="s0")

    assert writes == ["25", "80"]


def test_stage_progress_reporter_writes_events_jsonl_and_status_json(tmp_path: Path):
    reporter = StageProgressReporter(
        tmp_path / "S1_ConfSearch",
        "S1",
        default_fields={"condition_signature": "cond_1234"},
    )

    reporter.emit_stage_event("s1_started", total_variants=1)
    reporter.start_structure("product_major", smiles="CCO")
    reporter.finish_structure("product_major", "complete", selected="conf_0001")
    reporter.emit_stage_event("s1_completed")

    events = (tmp_path / "S1_ConfSearch" / "events.jsonl").read_text(encoding="utf-8")
    status = json.loads((tmp_path / "S1_ConfSearch" / "status.json").read_text(encoding="utf-8"))
    assert '"condition_signature": "cond_1234"' in events
    assert "structure_started" in events
    assert "structure_finished" in events
    assert status["stage"] == "S1"
    assert status["status"] == "completed"
    assert status["total_structures"] == 1
    assert status["completed"] == 1


def test_stage_progress_reporter_tracks_structure_lifecycle_and_task_updates(tmp_path: Path):
    reporter = StageProgressReporter(tmp_path / "S3_LowLevel", "S3")

    reporter.start_structure("product_major_ts", kind="ts")
    reporter.calculator_event(
        "optimization_started",
        {
            "structure_id": "product_major_ts",
            "engine": "orca",
            "method": "B97-3c",
            "solvent": "acetone",
            "solvent_model": "CPCM",
        },
    )
    reporter.calculator_event(
        "optimization_finished",
        {
            "structure_id": "product_major_ts",
            "engine": "orca",
            "method": "B97-3c",
            "solvent": "acetone",
            "solvent_model": "CPCM",
            "status": "complete",
            "output": "opt/output.out",
        },
    )
    reporter.finish_structure("product_major_ts", "failed", error="opt diverged")
    reporter.emit_stage_event("s3_completed")

    status = json.loads((tmp_path / "S3_LowLevel" / "status.json").read_text(encoding="utf-8"))
    structure = status["structures"]["product_major_ts"]
    assert status["failed"] == 1
    assert structure["status"] == "failed"
    assert structure["tasks"]["optimization"]["status"] == "complete"
    assert structure["tasks"]["optimization"]["output"] == "opt/output.out"
    assert structure["error"] == "opt diverged"


def test_scan_all_stages_returns_available_stage_statuses(tmp_path: Path):
    s0_reporter = StageProgressReporter(tmp_path / "S0_Mechanism", "S0")
    s0_reporter.emit_stage_event("s0_started")
    s0_reporter.emit_stage_event("s0_completed")

    s3_reporter = StageProgressReporter(tmp_path / "S3_LowLevel", "S3")
    s3_reporter.start_structure("product_major_ts")
    s3_reporter.finish_structure("product_major_ts", "complete")
    s3_reporter.emit_stage_event("s3_completed")

    stages = _scan_all_stages(tmp_path)
    assert stages["S0"]["status"] == "completed"
    assert stages["S3"]["structures"]["product_major_ts"]["status"] == "complete"


def test_render_overview_shows_cross_stage_summary(tmp_path: Path):
    s1_reporter = StageProgressReporter(tmp_path / "S1_ConfSearch", "S1")
    s1_reporter.start_structure("product_major", smiles="CCO")
    s1_reporter.finish_structure("product_major", "complete", selected="conf_0001")
    s1_reporter.emit_stage_event("s1_completed")

    s2_reporter = StageProgressReporter(tmp_path / "S2_PEB", "S2")
    s2_reporter.start_structure("product_major")
    s2_reporter.finish_structure("product_major", "degraded", stage_status="DEGRADED")
    s2_reporter.emit_stage_event("s2_completed")

    (tmp_path / "run.manifest.json").write_text(
        json.dumps(
            {
                "reaction_id": "RXN_P7P8",
                "condition_signature": "cond_1234",
                "completed_through": "s2",
                "conditions": {
                    "solvent": "acetone",
                    "temperature_celsius": "25",
                    "has_lewis_acid": "true",
                    "charge": 0,
                    "multiplicity": 1,
                    "additive_info": {
                        "catalyst": "LiCl",
                        "additive": "LiCl",
                        "has_lewis_acid": "true",
                    },
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    rendered = render_overview(tmp_path, colour=False)
    assert "RPH V4 Cross-Stage Overview" in rendered
    assert "S1: completed | 1/1 done" in rendered
    assert "S2: completed_with_failures | 1/1 done, 1 failed" in rendered
    assert "condition_signature: cond_1234" in rendered
    assert "additives: catalyst=LiCl, additive=LiCl" in rendered
