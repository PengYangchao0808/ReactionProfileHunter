import json
import os
import time
from pathlib import Path
from unittest.mock import Mock

from rph_core.utils import stale_recovery as stale_module
from rph_core.utils.stale_recovery import (
    RecoveryAction,
    StaleRecoveryDecision,
    classify_stage_recovery,
    is_pid_alive,
    recover_stage,
)
from rph_core.utils.v4_checkpoint import V4Checkpoint
from rph_core.v4_orchestrator import V4Orchestrator


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_is_pid_alive_current_and_missing_pid() -> None:
    assert is_pid_alive(os.getpid()) is True
    assert is_pid_alive(999999) is False


def test_classify_stage_recovery_missing_status_json(tmp_path: Path) -> None:
    decision = classify_stage_recovery(tmp_path / "S3_LowLevel" / "status.json")

    assert decision.action == RecoveryAction.INSUFFICIENT_DATA
    assert "status.json not found" in decision.evidence


def test_classify_stage_recovery_terminal_status_is_complete(tmp_path: Path) -> None:
    status_path = tmp_path / "S3_LowLevel" / "status.json"
    _write_json(status_path, {"status": "complete"})

    decision = classify_stage_recovery(status_path)

    assert decision.action == RecoveryAction.MARK_COMPLETE
    assert decision.output_complete is True


def test_classify_stage_recovery_pid_alive_with_fresh_heartbeat(
    monkeypatch, tmp_path: Path
) -> None:
    status_path = tmp_path / "S3_LowLevel" / "status.json"
    _write_json(
        status_path,
        {
            "status": "running",
            "pid": 1234,
            "heartbeat_at": time.time() - 10.0,
            "heartbeat_interval_seconds": 30,
        },
    )
    monkeypatch.setattr(stale_module, "is_pid_alive", lambda pid: pid == 1234)

    decision = classify_stage_recovery(status_path, stale_threshold_seconds=300.0)

    assert decision.action == RecoveryAction.CONTINUE_MONITORING
    assert decision.pid == 1234
    assert decision.heartbeat_age_seconds is not None


def test_classify_stage_recovery_pid_alive_with_stale_heartbeat(
    monkeypatch, tmp_path: Path
) -> None:
    status_path = tmp_path / "S3_LowLevel" / "status.json"
    _write_json(
        status_path,
        {
            "status": "running",
            "pid": 1234,
            "heartbeat_at": time.time() - 1200.0,
            "heartbeat_interval_seconds": 30,
        },
    )
    monkeypatch.setattr(stale_module, "is_pid_alive", lambda pid: pid == 1234)

    decision = classify_stage_recovery(status_path, stale_threshold_seconds=300.0)

    assert decision.action == RecoveryAction.INSUFFICIENT_DATA
    assert decision.pid == 1234


def test_classify_stage_recovery_pid_dead_with_complete_manifest(
    monkeypatch, tmp_path: Path
) -> None:
    stage_dir = tmp_path / "S3_LowLevel"
    _write_json(stage_dir / "status.json", {"status": "running", "pid": 1234})
    _write_json(
        stage_dir / "manifest.json",
        {"structures": [{"status": "complete"}, {"status": "failed"}]},
    )
    monkeypatch.setattr(stale_module, "is_pid_alive", lambda _pid: False)

    decision = classify_stage_recovery(stage_dir / "status.json")

    assert decision.action == RecoveryAction.MARK_COMPLETE
    assert decision.output_complete is True


def test_classify_stage_recovery_pid_dead_with_incomplete_outputs(
    monkeypatch, tmp_path: Path
) -> None:
    stage_dir = tmp_path / "S3_LowLevel"
    _write_json(stage_dir / "status.json", {"status": "running", "pid": 1234})
    monkeypatch.setattr(stale_module, "is_pid_alive", lambda _pid: False)

    decision = classify_stage_recovery(stage_dir / "status.json")

    assert decision.action == RecoveryAction.MARK_INTERRUPTED
    assert decision.output_complete is False


def test_classify_stage_recovery_missing_pid_field(tmp_path: Path) -> None:
    stage_dir = tmp_path / "S3_LowLevel"
    _write_json(stage_dir / "status.json", {"status": "running"})

    decision = classify_stage_recovery(stage_dir / "status.json")

    assert decision.action == RecoveryAction.INSUFFICIENT_DATA
    assert decision.pid is None


def test_recover_stage_marks_complete(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        stale_module,
        "classify_stage_recovery",
        lambda *_args, **_kwargs: StaleRecoveryDecision(
            action=RecoveryAction.MARK_COMPLETE,
            structure_id=None,
            stage="s3",
            evidence=("manifest complete",),
            pid=4321,
            heartbeat_age_seconds=None,
            output_complete=True,
        ),
    )
    checkpoint = Mock()

    action = recover_stage("s3", tmp_path / "S3_LowLevel", checkpoint)

    assert action == RecoveryAction.MARK_COMPLETE
    checkpoint.mark_recovered.assert_called_once_with("s3", status="complete")
    checkpoint.invalidate_from.assert_not_called()


def test_recover_stage_invalidates_interrupted_stage(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        stale_module,
        "classify_stage_recovery",
        lambda *_args, **_kwargs: StaleRecoveryDecision(
            action=RecoveryAction.MARK_INTERRUPTED,
            structure_id=None,
            stage="s3",
            evidence=("manifest incomplete",),
            pid=4321,
            heartbeat_age_seconds=None,
            output_complete=False,
        ),
    )
    checkpoint = Mock()

    action = recover_stage("s3", tmp_path / "S3_LowLevel", checkpoint)

    assert action == RecoveryAction.MARK_INTERRUPTED
    checkpoint.mark_recovered.assert_called_once_with("s3", status="interrupted")
    checkpoint.invalidate_from.assert_called_once_with("s3")


def test_orchestrator_recovery_invalidates_on_interrupted_action(
    monkeypatch, tmp_path: Path
) -> None:
    work_dir = tmp_path / "run"
    stage_dir = work_dir / "S3_LowLevel"
    stage_dir.mkdir(parents=True)
    manifest = stage_dir / "manifest.json"
    manifest.write_text("{}", encoding="utf-8")

    checkpoint = V4Checkpoint(work_dir)
    checkpoint.mark("s3", checkpoint.signature({"stage": "s3"}), manifest)
    original_invalidate_from = checkpoint.invalidate_from
    checkpoint.invalidate_from = Mock(wraps=original_invalidate_from)

    orchestrator = object.__new__(V4Orchestrator)
    orchestrator.config = {
        "ui": {
            "stalled_job_warning_seconds": 720,
            "recovery": {
                "enabled": True,
                "stale_heartbeat_seconds": 300,
                "requeue_interrupted_structures": True,
            },
        }
    }

    monkeypatch.setattr(
        "rph_core.v4_orchestrator.recover_stage",
        lambda *_args, **_kwargs: RecoveryAction.MARK_INTERRUPTED,
    )

    orchestrator._recover_stale_resume_state(checkpoint, work_dir)

    checkpoint.invalidate_from.assert_called_with("s3")


def test_v4_checkpoint_mark_recovered_preserves_signature_and_manifest(
    tmp_path: Path,
) -> None:
    checkpoint = V4Checkpoint(tmp_path)
    manifest = tmp_path / "S4_HighLevel" / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("{}", encoding="utf-8")
    signature = checkpoint.signature({"stage": "s4"})
    checkpoint.mark("s4", signature, manifest, status="incomplete")

    checkpoint.mark_recovered("s4", status="complete")
    state = checkpoint.load()
    entry = state["stages"]["s4"]

    assert checkpoint.has_stage("s4") is True
    assert entry["signature"] == signature
    assert entry["manifest"] == str(manifest.resolve())
    assert entry["status"] == "complete"
    assert entry["recovery_status"] == "complete"
    assert entry["recovered_at"]
