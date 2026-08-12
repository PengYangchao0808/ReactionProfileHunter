"""Resume-time stale-running recovery helpers for long-lived V4 stages."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import enum
import importlib
import json
import logging
import os
from pathlib import Path
import time
from typing import TYPE_CHECKING, Any


if TYPE_CHECKING:
    from collections.abc import Mapping

    from rph_core.utils.v4_checkpoint import V4Checkpoint


logger = logging.getLogger(__name__)


_STAGE_TERMINAL_STATUSES = frozenset(
    {
        "cached",
        "complete",
        "completed",
        "completed_with_failures",
        "degraded",
        "failed",
    }
)
_STRUCTURE_TERMINAL_STATUSES = frozenset(
    {
        "complete",
        "degraded",
        "failed",
        "minimum_frequency_unverified",
        "opt_failed_sp_complete",
        "ts_frequency_unverified",
    }
)


class RecoveryAction(str, enum.Enum):
    CONTINUE_MONITORING = "continue_monitoring"
    MARK_COMPLETE = "mark_complete"
    MARK_INTERRUPTED = "mark_interrupted"
    INSUFFICIENT_DATA = "insufficient_data"


@dataclass(frozen=True)
class StaleRecoveryDecision:
    action: RecoveryAction
    structure_id: str | None
    stage: str
    evidence: tuple[str, ...]
    pid: int | None
    heartbeat_age_seconds: float | None
    output_complete: bool | None


def is_pid_alive(pid: int) -> bool:
    """Return whether the PID currently exists without blocking on the process."""

    psutil_mod = None
    try:
        psutil_mod = importlib.import_module("psutil")
    except ImportError:
        pass

    if psutil_mod is not None:
        try:
            return bool(psutil_mod.pid_exists(pid))
        except Exception as exc:  # pragma: no cover - defensive psutil fallback
            logger.warning("psutil.pid_exists(%s) failed during stale recovery: %s", pid, exc)

    try:
        os.kill(pid, 0)
        return True
    except PermissionError:
        return True
    except ProcessLookupError:
        return False
    except OSError as exc:
        logger.warning("PID probe failed for %s during stale recovery: %s", pid, exc)
        return False


def classify_stage_recovery(
    stage_status_path: Path,
    *,
    stale_threshold_seconds: float = 300.0,
) -> StaleRecoveryDecision:
    """Classify one stage-level resume decision from its durable status snapshot."""

    stage_status_path = Path(stage_status_path)
    stage_name = _stage_name_from_dir(stage_status_path.parent)
    if not stage_status_path.exists():
        return StaleRecoveryDecision(
            action=RecoveryAction.INSUFFICIENT_DATA,
            structure_id=None,
            stage=stage_name,
            evidence=("status.json not found",),
            pid=None,
            heartbeat_age_seconds=None,
            output_complete=None,
        )

    payload = _load_json_mapping(stage_status_path)
    if payload is None:
        return StaleRecoveryDecision(
            action=RecoveryAction.INSUFFICIENT_DATA,
            structure_id=None,
            stage=stage_name,
            evidence=(f"could not parse {stage_status_path.name}",),
            pid=None,
            heartbeat_age_seconds=None,
            output_complete=None,
        )

    status = _normalized_status(payload.get("status"))
    if status in _STAGE_TERMINAL_STATUSES:
        return StaleRecoveryDecision(
            action=RecoveryAction.MARK_COMPLETE,
            structure_id=None,
            stage=stage_name,
            evidence=(f"stage status already terminal: {status}",),
            pid=_coerce_pid(payload.get("pid")),
            heartbeat_age_seconds=_age_seconds(payload.get("heartbeat_at")),
            output_complete=True,
        )

    pid = _coerce_pid(payload.get("pid"))
    if pid is None:
        return StaleRecoveryDecision(
            action=RecoveryAction.INSUFFICIENT_DATA,
            structure_id=None,
            stage=stage_name,
            evidence=("PID unknown",),
            pid=None,
            heartbeat_age_seconds=_age_seconds(payload.get("heartbeat_at")),
            output_complete=None,
        )

    heartbeat_age = _age_seconds(payload.get("heartbeat_at"))
    heartbeat_interval = _coerce_float(payload.get("heartbeat_interval_seconds"))
    if is_pid_alive(pid):
        if heartbeat_age is None:
            return StaleRecoveryDecision(
                action=RecoveryAction.CONTINUE_MONITORING,
                structure_id=None,
                stage=stage_name,
                evidence=("PID alive, heartbeat missing",),
                pid=pid,
                heartbeat_age_seconds=None,
                output_complete=None,
            )
        if heartbeat_age <= stale_threshold_seconds:
            return StaleRecoveryDecision(
                action=RecoveryAction.CONTINUE_MONITORING,
                structure_id=None,
                stage=stage_name,
                evidence=(
                    _fresh_heartbeat_evidence(heartbeat_age, stale_threshold_seconds, heartbeat_interval),
                ),
                pid=pid,
                heartbeat_age_seconds=heartbeat_age,
                output_complete=None,
            )
        return StaleRecoveryDecision(
            action=RecoveryAction.INSUFFICIENT_DATA,
            structure_id=None,
            stage=stage_name,
            evidence=(
                _stale_heartbeat_evidence(heartbeat_age, heartbeat_interval),
            ),
            pid=pid,
            heartbeat_age_seconds=heartbeat_age,
            output_complete=None,
        )

    manifest_path = stage_status_path.parent / "manifest.json"
    manifest_complete, manifest_evidence = _stage_manifest_complete(manifest_path)
    if manifest_complete:
        return StaleRecoveryDecision(
            action=RecoveryAction.MARK_COMPLETE,
            structure_id=None,
            stage=stage_name,
            evidence=(f"PID {pid} dead", *manifest_evidence),
            pid=pid,
            heartbeat_age_seconds=heartbeat_age,
            output_complete=True,
        )
    return StaleRecoveryDecision(
        action=RecoveryAction.MARK_INTERRUPTED,
        structure_id=None,
        stage=stage_name,
        evidence=(f"PID {pid} dead, manifest incomplete", *manifest_evidence),
        pid=pid,
        heartbeat_age_seconds=heartbeat_age,
        output_complete=False,
    )


def classify_structure_recovery(
    structure_status: Mapping[str, Any],
    *,
    stage_output_dir: Path,
    stale_threshold_seconds: float = 300.0,
) -> StaleRecoveryDecision:
    """Classify one S4 structure-level resume decision."""

    stage_output_dir = Path(stage_output_dir)
    stage_name = _stage_name_from_dir(stage_output_dir)
    structure_id = _structure_id(structure_status)
    status = _normalized_status(structure_status.get("status"))
    if status in _STRUCTURE_TERMINAL_STATUSES:
        return StaleRecoveryDecision(
            action=RecoveryAction.MARK_COMPLETE,
            structure_id=structure_id,
            stage=stage_name,
            evidence=(f"structure status already terminal: {status}",),
            pid=_coerce_pid(structure_status.get("pid")),
            heartbeat_age_seconds=_age_seconds(structure_status.get("heartbeat_at")),
            output_complete=True,
        )

    pid = _coerce_pid(structure_status.get("pid"))
    heartbeat_age = _age_seconds(structure_status.get("heartbeat_at"))
    last_output_age = _age_seconds(structure_status.get("last_output_update"))
    if pid is None:
        evidence = ["PID unknown"]
        if last_output_age is not None:
            evidence.append(f"last_output_update age={last_output_age:.1f}s")
        return StaleRecoveryDecision(
            action=RecoveryAction.INSUFFICIENT_DATA,
            structure_id=structure_id,
            stage=stage_name,
            evidence=tuple(evidence),
            pid=None,
            heartbeat_age_seconds=heartbeat_age,
            output_complete=None,
        )

    if is_pid_alive(pid):
        if heartbeat_age is None:
            evidence = ["PID alive, heartbeat missing"]
            if last_output_age is not None:
                evidence.append(f"last_output_update age={last_output_age:.1f}s")
            return StaleRecoveryDecision(
                action=RecoveryAction.CONTINUE_MONITORING,
                structure_id=structure_id,
                stage=stage_name,
                evidence=tuple(evidence),
                pid=pid,
                heartbeat_age_seconds=None,
                output_complete=None,
            )
        if heartbeat_age <= stale_threshold_seconds:
            evidence = [
                f"PID alive, heartbeat fresh ({heartbeat_age:.1f}s <= {stale_threshold_seconds:.1f}s)"
            ]
            if last_output_age is not None:
                evidence.append(f"last_output_update age={last_output_age:.1f}s")
            return StaleRecoveryDecision(
                action=RecoveryAction.CONTINUE_MONITORING,
                structure_id=structure_id,
                stage=stage_name,
                evidence=tuple(evidence),
                pid=pid,
                heartbeat_age_seconds=heartbeat_age,
                output_complete=None,
            )
        evidence = [f"PID alive but heartbeat stale ({heartbeat_age:.1f}s old)"]
        if last_output_age is not None:
            evidence.append(f"last_output_update age={last_output_age:.1f}s")
        return StaleRecoveryDecision(
            action=RecoveryAction.INSUFFICIENT_DATA,
            structure_id=structure_id,
            stage=stage_name,
            evidence=tuple(evidence),
            pid=pid,
            heartbeat_age_seconds=heartbeat_age,
            output_complete=None,
        )

    outputs_complete, output_evidence = _structure_outputs_complete(
        structure_status,
        stage_output_dir=stage_output_dir,
    )
    if outputs_complete:
        return StaleRecoveryDecision(
            action=RecoveryAction.MARK_COMPLETE,
            structure_id=structure_id,
            stage=stage_name,
            evidence=(f"PID {pid} dead", *output_evidence),
            pid=pid,
            heartbeat_age_seconds=heartbeat_age,
            output_complete=True,
        )
    return StaleRecoveryDecision(
        action=RecoveryAction.MARK_INTERRUPTED,
        structure_id=structure_id,
        stage=stage_name,
        evidence=(f"PID {pid} dead, outputs incomplete", *output_evidence),
        pid=pid,
        heartbeat_age_seconds=heartbeat_age,
        output_complete=False,
    )


def recover_stage(
    stage_name: str,
    stage_output_dir: Path,
    checkpoint: V4Checkpoint,
    *,
    stale_threshold_seconds: float = 300.0,
    logger: logging.Logger | None = None,
) -> RecoveryAction:
    """Apply stale-running recovery before checkpoint reuse is evaluated."""

    stage_name = str(stage_name).strip().lower()
    stage_output_dir = Path(stage_output_dir)
    decision = classify_stage_recovery(
        stage_output_dir / "status.json",
        stale_threshold_seconds=stale_threshold_seconds,
    )
    if stage_name == "s4":
        decision = _merge_s4_structure_recovery(
            decision,
            stage_output_dir=stage_output_dir,
            stale_threshold_seconds=stale_threshold_seconds,
        )

    active_logger = logger if logger is not None else logging.getLogger(__name__)
    _log_decision(active_logger, decision)

    if decision.action == RecoveryAction.CONTINUE_MONITORING:
        return decision.action
    if decision.action == RecoveryAction.INSUFFICIENT_DATA:
        return decision.action
    if decision.action == RecoveryAction.MARK_COMPLETE:
        checkpoint.mark_recovered(stage_name, status="complete")
        return decision.action

    checkpoint.mark_recovered(stage_name, status="interrupted")
    checkpoint.invalidate_from(stage_name)
    return decision.action


def _merge_s4_structure_recovery(
    stage_decision: StaleRecoveryDecision,
    *,
    stage_output_dir: Path,
    stale_threshold_seconds: float,
) -> StaleRecoveryDecision:
    status_path = Path(stage_output_dir) / "status.json"
    payload = _load_json_mapping(status_path)
    if payload is None:
        return stage_decision
    raw_structures = payload.get("structures")
    if not isinstance(raw_structures, list):
        return stage_decision

    structure_decisions = [
        classify_structure_recovery(
            row,
            stage_output_dir=stage_output_dir,
            stale_threshold_seconds=stale_threshold_seconds,
        )
        for row in raw_structures
        if isinstance(row, dict)
    ]
    if not structure_decisions or stage_decision.action in {
        RecoveryAction.CONTINUE_MONITORING,
        RecoveryAction.MARK_COMPLETE,
        RecoveryAction.MARK_INTERRUPTED,
    }:
        return stage_decision

    any_continue = any(
        item.action == RecoveryAction.CONTINUE_MONITORING for item in structure_decisions
    )
    any_interrupted = any(
        item.action == RecoveryAction.MARK_INTERRUPTED for item in structure_decisions
    )
    all_complete = all(item.action == RecoveryAction.MARK_COMPLETE for item in structure_decisions)
    summary = _structure_summary_evidence(structure_decisions)
    if any_continue:
        return StaleRecoveryDecision(
            action=RecoveryAction.CONTINUE_MONITORING,
            structure_id=None,
            stage=stage_decision.stage,
            evidence=stage_decision.evidence + summary,
            pid=stage_decision.pid,
            heartbeat_age_seconds=stage_decision.heartbeat_age_seconds,
            output_complete=None,
        )
    if any_interrupted:
        return StaleRecoveryDecision(
            action=RecoveryAction.MARK_INTERRUPTED,
            structure_id=None,
            stage=stage_decision.stage,
            evidence=stage_decision.evidence + summary,
            pid=stage_decision.pid,
            heartbeat_age_seconds=stage_decision.heartbeat_age_seconds,
            output_complete=False,
        )
    if all_complete:
        return StaleRecoveryDecision(
            action=RecoveryAction.MARK_COMPLETE,
            structure_id=None,
            stage=stage_decision.stage,
            evidence=stage_decision.evidence + summary,
            pid=stage_decision.pid,
            heartbeat_age_seconds=stage_decision.heartbeat_age_seconds,
            output_complete=True,
        )
    return stage_decision


def _stage_manifest_complete(manifest_path: Path) -> tuple[bool, tuple[str, ...]]:
    if not manifest_path.exists():
        return False, ("manifest.json not found",)
    payload = _load_json_mapping(manifest_path)
    if payload is None:
        return False, ("manifest.json unreadable",)
    structures = payload.get("structures")
    if not isinstance(structures, list):
        return False, ("manifest has no structures list",)
    incomplete = []
    for row in structures:
        status = _normalized_status((row or {}).get("status") if isinstance(row, dict) else None)
        if status not in _STRUCTURE_TERMINAL_STATUSES:
            incomplete.append(status or "unknown")
    if incomplete:
        return False, (
            f"manifest has nonterminal structure statuses: {', '.join(incomplete[:5])}",
        )
    return True, (f"manifest complete with {len(structures)} terminal structures",)


def _structure_outputs_complete(
    structure_status: Mapping[str, Any],
    *,
    stage_output_dir: Path,
) -> tuple[bool, tuple[str, ...]]:
    structure_id = _structure_id(structure_status)
    if structure_id is None:
        return False, ("structure id missing",)
    structure_dir = Path(stage_output_dir) / structure_id
    if not structure_dir.exists():
        return False, (f"structure directory not found: {structure_dir}",)

    kind = _normalized_status(structure_status.get("kind")) or "minimum"
    task_definitions = {
        "optimization": ("opt", "opt_ts") if kind == "ts" else ("opt",),
        "frequency": ("freq",),
        "single_point": ("sp",),
    }
    evidence: list[str] = []
    complete = True
    for task_name, directory_names in task_definitions.items():
        task_complete, task_evidence = _task_output_complete(
            structure_status,
            task_name=task_name,
            structure_dir=structure_dir,
            directory_names=directory_names,
        )
        evidence.extend(task_evidence)
        complete = complete and task_complete
    return complete, tuple(evidence)


def _task_output_complete(
    structure_status: Mapping[str, Any],
    *,
    task_name: str,
    structure_dir: Path,
    directory_names: tuple[str, ...],
) -> tuple[bool, tuple[str, ...]]:
    explicit_paths = _task_output_paths(structure_status, task_name)
    if explicit_paths:
        valid_paths = [path for path in explicit_paths if _file_complete(path)]
        if valid_paths:
            return True, (f"{task_name} output present: {valid_paths[0]}",)
        return False, (
            f"{task_name} recorded outputs missing or empty: {', '.join(str(path) for path in explicit_paths)}",
        )

    for directory_name in directory_names:
        task_dir = structure_dir / directory_name
        if _directory_has_nonempty_file(task_dir):
            return True, (f"{task_name} output inferred from {task_dir}",)
    return False, (
        f"{task_name} output missing in {', '.join(str(structure_dir / name) for name in directory_names)}",
    )


def _task_output_paths(
    structure_status: Mapping[str, Any],
    task_name: str,
) -> tuple[Path, ...]:
    tasks = structure_status.get("tasks")
    paths: list[Path] = []
    if isinstance(tasks, dict):
        task = tasks.get(task_name)
        if isinstance(task, dict):
            output = task.get("output")
            if isinstance(output, str) and output.strip():
                paths.append(Path(output))

    top_level_key = {
        "optimization": "opt_output",
        "frequency": "frequency_output",
        "single_point": "sp_output",
    }.get(task_name)
    if top_level_key:
        output = structure_status.get(top_level_key)
        if isinstance(output, str) and output.strip():
            paths.append(Path(output))
    return tuple(paths)


def _load_json_mapping(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("Could not read stale recovery payload at %s: %s", path, exc)
        return None
    if not isinstance(payload, dict):
        logger.warning("Stale recovery payload at %s is not a JSON object", path)
        return None
    return payload


def _structure_summary_evidence(
    decisions: list[StaleRecoveryDecision],
) -> tuple[str, ...]:
    summary: dict[str, int] = {}
    for decision in decisions:
        summary[decision.action.value] = summary.get(decision.action.value, 0) + 1
    if not summary:
        return ()
    return (
        "structure recovery summary: "
        + ", ".join(f"{key}={summary[key]}" for key in sorted(summary)),
    )


def _log_decision(active_logger: logging.Logger, decision: StaleRecoveryDecision) -> None:
    log = active_logger.warning if decision.action == RecoveryAction.INSUFFICIENT_DATA else active_logger.info
    log(
        "Stale recovery decision stage=%s structure_id=%s action=%s pid=%s heartbeat_age_seconds=%s output_complete=%s evidence=%s",
        decision.stage,
        decision.structure_id,
        decision.action.value,
        decision.pid,
        _format_heartbeat_age(decision.heartbeat_age_seconds),
        decision.output_complete,
        " | ".join(decision.evidence),
    )


def _normalized_status(value: Any) -> str:
    return str(value or "").strip().lower()


def _coerce_pid(value: Any) -> int | None:
    try:
        pid = int(value)
    except (TypeError, ValueError):
        return None
    return pid if pid > 0 else None


def _coerce_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _age_seconds(value: Any) -> float | None:
    timestamp = _timestamp_to_epoch(value)
    if timestamp is None:
        return None
    return max(0.0, time.time() - timestamp)


def _timestamp_to_epoch(value: Any) -> float | None:
    if value is None:
        return None
    numeric = _coerce_float(value)
    if numeric is not None:
        return numeric
    if not isinstance(value, str):
        return None
    text = value.strip()
    if not text:
        return None
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.timestamp()


def _stage_name_from_dir(path: Path) -> str:
    name = Path(path).name
    mapping = {
        "S3_LowLevel": "s3",
        "S4_HighLevel": "s4",
    }
    return mapping.get(name, name.lower())


def _structure_id(structure_status: Mapping[str, Any]) -> str | None:
    value = structure_status.get("id")
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _file_complete(path: Path) -> bool:
    try:
        return Path(path).is_file() and Path(path).stat().st_size > 0
    except OSError:
        return False


def _directory_has_nonempty_file(path: Path) -> bool:
    try:
        return path.is_dir() and any(child.is_file() and child.stat().st_size > 0 for child in path.iterdir())
    except OSError:
        return False


def _fresh_heartbeat_evidence(
    heartbeat_age: float,
    stale_threshold_seconds: float,
    heartbeat_interval: float | None,
) -> str:
    interval_suffix = (
        f", interval={heartbeat_interval:.1f}s" if heartbeat_interval is not None else ""
    )
    return (
        f"PID alive, heartbeat fresh ({heartbeat_age:.1f}s <= {stale_threshold_seconds:.1f}s{interval_suffix})"
    )


def _stale_heartbeat_evidence(
    heartbeat_age: float,
    heartbeat_interval: float | None,
) -> str:
    interval_suffix = (
        f", interval={heartbeat_interval:.1f}s" if heartbeat_interval is not None else ""
    )
    return f"PID alive but heartbeat stale ({heartbeat_age:.1f}s old{interval_suffix})"


def _format_heartbeat_age(value: float | None) -> str | None:
    if value is None:
        return None
    return f"{value:.1f}"
