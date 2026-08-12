"""Unified S3/S4 refinement engine."""

from __future__ import annotations

import copy
import hashlib
import inspect
import json
import logging
import re
import shutil
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

from rph_core.steps.fidelity_profile import FidelityProfile
from rph_core.steps.refinement.manifest_io import REFINEMENT_MANIFEST_V1, write_refinement_manifest
from rph_core.steps.refinement.models import Pass1Outcome, PreflightOutcome, StructureRequest
from rph_core.steps.refinement.rescue_matrix import (
    FailureType,
    MethodParams,
    RescueMethod,
    RescuePlan,
    StructureKind,
    lookup_plan,
    methods_params,
)
from rph_core.steps.refinement.step_codes import (
    RESCUE_SUB_DIAGNOSIS,
    StepCode,
    StepRecord,
    rescue_sub_code,
    rescue_sub_for_family,
)
from rph_core.utils.attempt_recorder import AttemptRecord, AttemptRecorder
from rph_core.utils.file_io import read_xyz, write_xyz
from rph_core.utils.geometry_tools import GeometryUtils
from rph_core.utils import identity as identity_module
from rph_core.utils.json_io import write_json_atomic
from rph_core.utils.ml_quality import compute_legacy_usable_for_ml, compute_ml_usability
from rph_core.utils.orca_interface import (
    _ELEMENT_MASS,
    _parse_orca_displacement_vectors,
    extract_last_orca_geometry_xyz_text,
    parse_orca6_normal_mode_vectors,
)
from rph_core.utils.orca_optimization_monitor import (
    OptimizationMonitorConfig,
    detect_oscillation,
    newest_orca_output,
    parse_orca_optimization_cycles,
)
from rph_core.utils.provenance import build_provenance
from rph_core.utils.qc_jobs import run_frequency, run_irc, run_optimization, run_single_point
from rph_core.utils.qc_models import IRCJobSpec, QCJobResult, QCJobSpec
from rph_core.utils.irc_trajectory import (
    classify_irc_endpoints,
    detect_intermediate_shoulder,
    frame_to_xyz,
    parse_irc_trajectories,
)
from rph_core.utils.stage_scheduler import resolve_stage_schedule, run_structure_queue
from rph_core.utils.superseded_archive import archive_stale_stage_outputs


logger = logging.getLogger(__name__)

_MODE_DISPLACEMENT_SCALE_ANGSTROM = 0.1
# An imaginary mode on a minimum/INT is significant only below this value,
# matching rph_core.utils.identity.classify_minimum's default cutoff (-10 cm^-1).
_MINIMUM_IMAGINARY_CUTOFF_CM1 = -10.0


class _LiveOptimizationMonitor:
    """Tail an ORCA optimization and record confirmed oscillation diagnostics."""

    def __init__(self, config: OptimizationMonitorConfig, work_dir: Path):
        self.config = config
        self.work_dir = Path(work_dir)
        self.trigger: dict[str, Any] | None = None
        self.error: dict[str, Any] | None = None
        self._active_output_path: Path | None = None
        self._thread: threading.Thread | None = None
        self._stop_event = threading.Event()

    def start(self, process: Any, *, output_path: str | Path | None = None) -> None:
        if not self.config.enabled or self._thread is not None:
            return
        if output_path is not None:
            self._active_output_path = Path(output_path)

        def _loop() -> None:
            while not self._stop_event.wait(self.config.poll_seconds):
                if getattr(process, "poll", lambda: 0)() is not None:
                    return
                output_path = self._active_output_path or newest_orca_output(self.work_dir)
                if output_path is None:
                    continue
                try:
                    records = parse_orca_optimization_cycles(
                        output_path.read_text(encoding="utf-8", errors="ignore")
                    )
                    trigger = detect_oscillation(records, self.config)
                except OSError:
                    continue
                except Exception as exc:  # Defensive boundary for a background monitor.
                    self.error = {
                        "reason": "monitor_error",
                        "exception_type": type(exc).__name__,
                        "message": str(exc),
                        "output_file": str(output_path),
                    }
                    logger.exception("ORCA optimization monitor failed: %s", self.error)
                    return
                if trigger is None:
                    continue
                trigger["output_file"] = str(output_path)
                self.trigger = trigger
                logger.warning(
                    "Observed suspected ORCA optimization oscillation; retaining native ORCA checkpoint: %s",
                    trigger,
                )
                return

        self._thread = threading.Thread(
            target=_loop,
            name=f"rph-orca-monitor-{self.work_dir.name}",
            daemon=True,
        )
        self._thread.start()

    def close(self) -> None:
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=max(1.0, self.config.poll_seconds + 1.0))

class RefinementEngine:
    """Single S3/S4 implementation class."""

    stage_name: str = "Refinement"

    def __init__(
        self,
        config: Dict[str, Any],
        profile: FidelityProfile,
        run_id: Optional[str] = None,
        parent_manifest_paths: Optional[Mapping[str, Path]] = None,
    ):
        self.config = config
        self.profile = profile
        self.run_id = run_id
        self.parent_manifest_paths = {
            str(stage): (Path(path) if path is not None else None)
            for stage, path in (parent_manifest_paths or {}).items()
        }
        self.progress_reporter: Any = None
        self.event_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None
        self._recorders: Dict[str, AttemptRecorder] = {}
        self._preflight_by_id: Dict[str, PreflightOutcome] = {}
        self._requests_by_id: Dict[str, StructureRequest] = {}
        self._variant_ts_outcome: Dict[str, Pass1Outcome] = {}
        self._variant_product_outcome: Dict[str, Pass1Outcome] = {}

    def set_progress_reporter(self, reporter: Any) -> None:
        self.progress_reporter = reporter

    def set_event_callback(
        self,
        callback: Optional[Callable[[str, Dict[str, Any]], None]],
    ) -> None:
        """Backward-compat hook for the legacy S4 engine's event_callback parameter.

        Phase 3+: the unified progress reporter consumes events internally;
        this hook is preserved so callers of the old HighLevelEngine signature
        keep working.
        """
        self.event_callback = callback

    def run(
        self,
        requests: Iterable[StructureRequest | Mapping[str, Any]],
        output_dir: Path,
        event_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
        *,
        resume_incomplete: bool = False,
        structure_ids: Optional[Sequence[str]] = None,
        rescue_only: bool = False,
    ) -> Path:
        if event_callback is not None:
            self.set_event_callback(event_callback)
        output_dir = Path(output_dir)
        stage = self.profile.stage
        materialized_requests = [self._coerce_to_structure_request(request) for request in requests]
        selected_ids = {
            str(structure_id).strip()
            for structure_id in (structure_ids or ())
            if str(structure_id).strip()
        }
        partial_resume = bool(resume_incomplete or selected_ids)
        existing_payloads: Dict[str, Dict[str, Any]] = {}
        rerun_ids: set[str] = set()
        archived_to: Path | None = None

        if partial_resume:
            existing_payloads = self._load_existing_structure_payloads(output_dir)
            if not existing_payloads:
                raise RuntimeError(
                    f"{stage} partial rerun requested, but {output_dir / 'manifest.json'} has no reusable structure records"
                )
            request_ids = {request.id for request in materialized_requests}
            unknown_ids = selected_ids - request_ids
            if unknown_ids:
                raise ValueError(
                    f"{stage} partial rerun requested unknown structure ids: {sorted(unknown_ids)}"
                )
            if selected_ids:
                rerun_ids = selected_ids
            else:
                rerun_ids = {
                    request.id
                    for request in materialized_requests
                    if not self._existing_structure_is_reusable(
                        existing_payloads.get(request.id)
                    )
                }
            if rescue_only:
                collapsed_ids = self._preserved_int_collapse_ids(
                    materialized_requests, existing_payloads, rerun_ids
                )
                if collapsed_ids:
                    logger.info(
                        "%s rescue-only: %s preserved INT reclassified as collapsed -> routed to IRC rescue",
                        stage,
                        sorted(collapsed_ids),
                    )
                    rerun_ids |= collapsed_ids
            logger.info(
                "%s partial rerun: recalculating=%s preserving=%s",
                stage,
                sorted(rerun_ids),
                sorted(request_ids - rerun_ids),
            )
        else:
            archived_to = archive_stale_stage_outputs(
                output_dir,
                stage_name=stage,
                run_id=self.run_id,
            )
            if archived_to is not None:
                self._safe_reinitialize_progress_after_archive()

        output_dir.mkdir(parents=True, exist_ok=True)
        self._requests_by_id = {request.id: request for request in materialized_requests}
        for request in materialized_requests:
            self._safe_register_structure(request)

        requests_to_run = [
            request
            for request in materialized_requests
            if not partial_resume or request.id in rerun_ids
        ]
        preflight_outcomes = self._run_pass0_preflight(requests_to_run, output_dir)
        self._preflight_by_id = {
            outcome.structure_id: outcome for outcome in preflight_outcomes
        }
        active_requests = [
            request
            for request, outcome in zip(requests_to_run, preflight_outcomes)
            if outcome.status == "ok"
        ]

        schedule = resolve_stage_schedule(
            self.config,
            "step3" if stage == "S3" else "step4",
        )
        active_request_payloads = [
            {"id": request.id, "role": request.role, "kind": request.kind, "request": request}
            for request in active_requests
        ]
        if rescue_only and partial_resume:
            replay_ids, fallback_ids = self._partition_rescue_replay(
                active_requests, existing_payloads
            )
            replayed_outcomes = self._replay_failed_outcomes(existing_payloads, replay_ids)
            fallback_payloads = [
                payload
                for payload in active_request_payloads
                if payload["id"] in fallback_ids
            ]
            fallback_outcomes = (
                run_structure_queue(
                    fallback_payloads,
                    lambda payload: self._run_pass1_primary_one(payload["request"], output_dir),
                    schedule,
                    thread_name_prefix=f"rph-{stage.lower()}-pass1",
                )
                if fallback_payloads
                else []
            )
            pass1_outcomes = replayed_outcomes + fallback_outcomes
        else:
            pass1_outcomes = run_structure_queue(
                active_request_payloads,
                lambda payload: self._run_pass1_primary_one(payload["request"], output_dir),
                schedule,
                thread_name_prefix=f"rph-{stage.lower()}-pass1",
            )

        self._variant_ts_outcome.clear()
        self._variant_product_outcome.clear()
        for outcome in pass1_outcomes:
            request = self._requests_by_id.get(outcome.structure_id)
            if request is None:
                continue
            variant_key = request.variant_id or request.id
            if request.kind == "ts":
                self._variant_ts_outcome[variant_key] = outcome
            if request.role == "product":
                self._variant_product_outcome[variant_key] = outcome
        if partial_resume:
            for request in materialized_requests:
                if request.kind != "ts" or request.id in rerun_ids:
                    continue
                payload = existing_payloads.get(request.id)
                if payload is None:
                    continue
                canonical_xyz = payload.get("canonical_xyz")
                hessian = payload.get("canonical_hessian_path")
                if not canonical_xyz or not Path(str(canonical_xyz)).is_file():
                    continue
                ts_outcome = Pass1Outcome(
                    structure_id=request.id,
                    role=str(payload.get("role") or request.role),
                    kind="ts",
                )
                ts_outcome.canonical_xyz = Path(str(canonical_xyz))
                ts_outcome.canonical_hessian_path = (
                    Path(str(hessian)) if hessian and Path(str(hessian)).is_file() else None
                )
                ts_outcome.frequency_hessian_path = ts_outcome.canonical_hessian_path
                ts_outcome.canonical_frequency_output = (
                    Path(str(payload.get("canonical_frequency_output")))
                    if payload.get("canonical_frequency_output")
                    and Path(str(payload.get("canonical_frequency_output"))).is_file()
                    else None
                )
                ts_outcome.status = "complete"
                self._variant_ts_outcome.setdefault(request.variant_id or request.id, ts_outcome)

        pass2_outcomes = self._run_pass2_rescue(pass1_outcomes, output_dir)
        pass3_outcomes = self._run_pass3_canonical(pass2_outcomes, output_dir)

        pass1_by_id = {outcome.structure_id: outcome for outcome in pass1_outcomes}
        pass2_by_id = {outcome.structure_id: outcome for outcome in pass2_outcomes}
        pass3_by_id = {outcome.structure_id: outcome for outcome in pass3_outcomes}

        structures_payload = []
        preflight_by_id = {outcome.structure_id: outcome for outcome in preflight_outcomes}
        for request in materialized_requests:
            if partial_resume and request.id not in rerun_ids:
                preserved_payload = copy.deepcopy(existing_payloads[request.id])
                structures_payload.append(preserved_payload)
                self._safe_finish_structure(request.id, preserved_payload)
                continue

            preflight = preflight_by_id[request.id]
            if preflight.status != "ok":
                structure_payload = self._build_failed_preflight_payload(request, preflight)
                structures_payload.append(structure_payload)
                self._safe_finish_structure(request.id, structure_payload)
                continue

            pass1 = pass1_by_id[request.id]
            pass2 = pass2_by_id.get(request.id, pass1)
            pass3 = pass3_by_id.get(request.id, pass2)
            structure_payload = self._build_structure_payload(
                request,
                preflight,
                pass1,
                pass2,
                pass3,
            )
            structures_payload.append(structure_payload)
            self._safe_finish_structure(request.id, structure_payload)

        manifest_path = self._write_manifest(
            output_dir,
            structures_payload,
            archived_to=archived_to,
            schedule_record=schedule.manifest_record(),
            source_requests=materialized_requests,
            partial_rerun_ids=sorted(rerun_ids) if partial_resume else None,
            preserved_structure_ids=(
                sorted({request.id for request in materialized_requests} - rerun_ids)
                if partial_resume
                else None
            ),
        )
        logger.info(
            "RefinementEngine (%s, profile=%s) wrote manifest: %s",
            stage,
            self.profile.profile_id,
            manifest_path,
        )
        return manifest_path

    def _heartbeat_interval_seconds(self) -> float:
        interval = getattr(self.progress_reporter, "heartbeat_interval_seconds", None)
        if interval is not None:
            try:
                return max(0.01, float(interval))
            except (TypeError, ValueError):
                pass
        try:
            ui_cfg = dict(self.config.get("ui", {}) or {})
            return max(0.01, float(ui_cfg.get("heartbeat_seconds", 30)))
        except (TypeError, ValueError):
            return 30.0

    def _supports_subprocess_tracking(self) -> bool:
        reporter = self.progress_reporter
        return bool(
            reporter is not None
            and hasattr(reporter, "report_structure_subprocess")
            and hasattr(reporter, "touch_structure_heartbeat")
        )

    def _safe_report_structure_subprocess(
        self,
        structure_id: str,
        *,
        pid: int | None,
        sandbox_path: str | None,
    ) -> None:
        reporter = self.progress_reporter
        if reporter is None or not hasattr(reporter, "report_structure_subprocess"):
            return
        try:
            reporter.report_structure_subprocess(
                structure_id,
                pid=pid,
                sandbox_path=sandbox_path,
            )
        except Exception as exc:  # pragma: no cover - defensive observability isolation
            logger.warning(
                "Ignoring subprocess progress callback failure for %s: %s",
                structure_id,
                exc,
            )

    def _safe_touch_structure_heartbeat(self, structure_id: str) -> None:
        reporter = self.progress_reporter
        if reporter is None or not hasattr(reporter, "touch_structure_heartbeat"):
            return
        try:
            reporter.touch_structure_heartbeat(structure_id)
        except Exception as exc:  # pragma: no cover - defensive observability isolation
            logger.warning(
                "Ignoring structure heartbeat failure for %s: %s",
                structure_id,
                exc,
            )

    def _safe_reinitialize_progress_after_archive(self) -> None:
        reporter = self.progress_reporter
        if reporter is None or not hasattr(reporter, "reinitialize_after_archive"):
            return
        try:
            reporter.reinitialize_after_archive()
        except Exception as exc:  # pragma: no cover - defensive observability isolation
            logger.warning("Ignoring progress reinitialization failure: %s", exc)

    def _safe_register_structure(self, request: StructureRequest) -> None:
        reporter = self.progress_reporter
        if reporter is None or not hasattr(reporter, "register_structure"):
            return
        try:
            reporter.register_structure(
                request.id,
                role=request.role,
                kind=request.kind,
                source_stage=request.source_stage,
            )
        except Exception as exc:  # pragma: no cover - defensive observability isolation
            logger.warning("Ignoring structure registration failure for %s: %s", request.id, exc)

    def _safe_start_structure(self, request: StructureRequest) -> None:
        reporter = self.progress_reporter
        if reporter is None or not hasattr(reporter, "start_structure"):
            return
        try:
            reporter.start_structure(
                request.id,
                role=request.role,
                kind=request.kind,
                source_stage=request.source_stage,
            )
        except Exception as exc:  # pragma: no cover - defensive observability isolation
            logger.warning("Ignoring structure start failure for %s: %s", request.id, exc)

    def _safe_finish_structure(self, structure_id: str, payload: Mapping[str, Any]) -> None:
        reporter = self.progress_reporter
        if reporter is None or not hasattr(reporter, "finish_structure"):
            return
        results = dict(payload)
        results.pop("id", None)
        results.pop("structure_id", None)
        status = str(results.pop("status", "complete"))
        try:
            reporter.finish_structure(structure_id, status, **results)
        except TypeError:
            # S4's legacy reporter accepts one payload mapping instead.  Keep
            # reporter incompatibilities isolated from the QC calculation.
            return
        except Exception as exc:  # pragma: no cover - defensive observability isolation
            logger.warning("Ignoring structure completion failure for %s: %s", structure_id, exc)

    def _safe_calculator_event(
        self,
        event: str,
        structure_id: str,
        payload: Mapping[str, Any],
    ) -> None:
        reporter = self.progress_reporter
        if reporter is None or not hasattr(reporter, "calculator_event"):
            return
        event_payload = {"structure_id": structure_id, **dict(payload)}
        try:
            reporter.calculator_event(event, event_payload)
        except Exception as exc:  # pragma: no cover - defensive observability isolation
            logger.warning("Ignoring QC task progress failure for %s: %s", structure_id, exc)

    def _start_heartbeat_thread(
        self,
        structure_id: str,
        work_dir: Path,
        process: Any,
        stop_event: threading.Event,
    ) -> threading.Thread:
        interval_seconds = self._heartbeat_interval_seconds()

        def _heartbeat_loop() -> None:
            while not stop_event.is_set():
                self._safe_touch_structure_heartbeat(structure_id)
                if getattr(process, "poll", lambda: 0)() is not None:
                    break
                if stop_event.wait(interval_seconds):
                    break
            self._safe_touch_structure_heartbeat(structure_id)

        heartbeat_thread = threading.Thread(
            target=_heartbeat_loop,
            name=f"rph-structure-heartbeat-{structure_id}-{Path(work_dir).name}",
            daemon=True,
        )
        heartbeat_thread.start()
        return heartbeat_thread

    def _run_job_with_subprocess_tracking(
        self,
        structure_id: str,
        work_dir: Path,
        runner: Callable[..., QCJobResult],
        *runner_args: Any,
        optimization_monitor: _LiveOptimizationMonitor | None = None,
    ) -> QCJobResult:
        accepts_callback = self._runner_accepts_subprocess_callback(runner)
        if not self._supports_subprocess_tracking() and optimization_monitor is None:
            return runner(*runner_args)
        if not accepts_callback:
            # Test doubles and third-party runners are allowed to omit the
            # callback.  They simply cannot participate in live monitoring.
            return runner(*runner_args)

        stop_event = threading.Event()
        heartbeat_threads: List[threading.Thread] = []

        def subprocess_callback(process: Any) -> None:
            self._safe_report_structure_subprocess(
                structure_id,
                pid=getattr(process, "pid", None),
                sandbox_path=str(Path(work_dir)),
            )
            self._safe_touch_structure_heartbeat(structure_id)
            if not heartbeat_threads:
                heartbeat_threads.append(
                    self._start_heartbeat_thread(
                        structure_id,
                        Path(work_dir),
                        process,
                        stop_event,
                    )
                )
            if optimization_monitor is not None:
                optimization_monitor.start(
                    process,
                    output_path=getattr(process, "_rph_orca_output_path", None),
                )

        try:
            result = runner(*runner_args, subprocess_callback=subprocess_callback)
        finally:
            stop_event.set()
            if heartbeat_threads:
                heartbeat_threads[0].join(timeout=1.0)
            if optimization_monitor is not None:
                optimization_monitor.close()
        if optimization_monitor is not None and optimization_monitor.trigger is not None:
            result.extra = {
                **dict(result.extra or {}),
                "optimization_monitor": dict(optimization_monitor.trigger),
            }
        if optimization_monitor is not None and optimization_monitor.error is not None:
            result.extra = {
                **dict(result.extra or {}),
                "optimization_monitor_error": dict(optimization_monitor.error),
            }
        return result

    @staticmethod
    def _runner_accepts_subprocess_callback(runner: Callable[..., QCJobResult]) -> bool:
        try:
            signature = inspect.signature(runner)
        except (TypeError, ValueError):
            return False
        return "subprocess_callback" in signature.parameters or any(
            parameter.kind == inspect.Parameter.VAR_KEYWORD
            for parameter in signature.parameters.values()
        )

    def _optimization_monitor_for(
        self,
        *,
        spec: QCJobSpec,
        work_dir: Path,
    ) -> _LiveOptimizationMonitor | None:
        if self.profile.stage != "S3" or spec.task.lower() not in {"opt_ts", "ts", "opt-ts"}:
            return None
        return _LiveOptimizationMonitor(
            OptimizationMonitorConfig(
                enabled=self.profile.live_monitor_enabled,
                poll_seconds=self.profile.live_monitor_poll_seconds,
                min_cycles=self.profile.live_monitor_min_cycles,
                window_cycles=self.profile.live_monitor_window_cycles,
                min_energy_reversals=self.profile.live_monitor_min_energy_reversals,
                energy_change_floor_hartree=self.profile.live_monitor_energy_change_floor_hartree,
                min_gradient_improvement_fraction=(
                    self.profile.live_monitor_min_gradient_improvement_fraction
                ),
                step_tolerance_ratio=self.profile.live_monitor_step_tolerance_ratio,
                min_excessive_steps=self.profile.live_monitor_min_excessive_steps,
            ),
            work_dir,
        )

    def _coerce_to_structure_request(
        self,
        request: Any,
    ) -> StructureRequest:
        if isinstance(request, StructureRequest):
            return request
        if not isinstance(request, Mapping):
            raise TypeError(f"Unsupported refinement request type: {type(request)!r}")

        role = str(request.get("role") or request.get("kind") or "product").strip().lower()
        kind = str(request.get("kind") or ("ts" if role == "ts" else "minimum")).strip().lower()
        structure_id = str(request.get("id") or request.get("structure_id") or "").strip()
        if not structure_id:
            raise ValueError("Structure request is missing required field: id")
        input_xyz_value = request.get("input_xyz")
        if input_xyz_value is None:
            raise ValueError(f"Structure request {structure_id!r} is missing input_xyz")

        return StructureRequest(
            id=structure_id,
            role=role,
            kind=kind,
            input_xyz=Path(input_xyz_value),
            forming_bonds=self._coerce_forming_bonds(request.get("forming_bonds")),
            fallback_xyz=self._coerce_optional_path(request.get("fallback_xyz")),
            original_seed_xyz=self._coerce_optional_path(request.get("original_seed_xyz")),
            source_stage=str(request.get("source_stage") or "S2").strip().upper(),
            seed_state=str(request.get("seed_state") or "stable_minimum_seed"),
            charge=int(request.get("charge", 0)),
            multiplicity=int(request.get("multiplicity", 1)),
            atom_mapping=self._coerce_optional_path(request.get("atom_mapping")),
            parent_structure=(
                dict(request.get("parent_structure", {}) or {})
                if request.get("parent_structure") is not None
                else None
            ),
            structure_id=self._coerce_optional_str(request.get("structure_id")),
            variant_id=self._coerce_optional_str(request.get("variant_id")),
            branch_id=self._coerce_optional_str(request.get("branch_id")),
            pathway_id=self._coerce_optional_str(request.get("pathway_id")),
            parent_structure_id=self._coerce_optional_str(request.get("parent_structure_id")),
            mapping_audit=self._coerce_optional_str(request.get("mapping_audit")),
            mapping_required=bool(request.get("mapping_required", False)),
            s1_manifest=self._coerce_optional_str(request.get("s1_manifest")),
            s1_ensemble_thermodynamics=(
                dict(request.get("s1_ensemble_thermodynamics", {}) or {})
                if request.get("s1_ensemble_thermodynamics") is not None
                else None
            ),
            s1_thermochemistry_status=self._coerce_optional_str(
                request.get("s1_thermochemistry_status")
            ),
            ensemble_thermochemistry_correction_hartree=self._coerce_optional_float(
                request.get("ensemble_thermochemistry_correction_hartree")
            ),
        )

    def _run_pass0_preflight(
        self,
        requests: Sequence[StructureRequest],
        output_dir: Path,
    ) -> List[PreflightOutcome]:
        outcomes: List[PreflightOutcome] = []
        self._recorders = {}

        for request in requests:
            try:
                if request.mapping_required and request.atom_mapping is None:
                    raise ValueError(f"atom_mapping required for {request.id}")
                if request.atom_mapping is not None and not Path(request.atom_mapping).exists():
                    raise ValueError(f"atom_mapping not found: {request.atom_mapping}")

                resolved_input = self._resolve_input_xyz(request)
                atom_count = self._read_atom_count(resolved_input)
                for atom_i, atom_j in request.forming_bonds:
                    if atom_i < 0 or atom_j < 0 or atom_i >= atom_count or atom_j >= atom_count:
                        raise ValueError(
                            f"forming_bond {(atom_i, atom_j)} out of range for atom_count={atom_count}"
                        )

                structure_dir = output_dir / request.id
                structure_dir.mkdir(parents=True, exist_ok=True)
                (structure_dir / "attempts").mkdir(parents=True, exist_ok=True)
                (structure_dir / "diagnostics").mkdir(parents=True, exist_ok=True)
                recorder = AttemptRecorder(structure_dir, structure_id=request.id)
                self._recorders[request.id] = recorder
                self._write_provenance(request, structure_dir, resolved_input)

                outcomes.append(
                    PreflightOutcome(
                        structure_id=request.id,
                        status="ok",
                        input_xyz=resolved_input,
                        charge=int(request.charge),
                        multiplicity=int(request.multiplicity),
                        forming_bonds=list(request.forming_bonds),
                    )
                )
            except Exception as exc:
                logger.warning("Preflight failed for %s: %s", request.id, exc)
                outcomes.append(
                    PreflightOutcome(
                        structure_id=request.id,
                        status="failed_preflight",
                        error=str(exc),
                    )
                )
        return outcomes

    def _run_pass1_primary_one(
        self,
        request: StructureRequest,
        output_dir: Path,
    ) -> Pass1Outcome:
        preflight = self._preflight_by_id[request.id]
        if preflight.input_xyz is None:
            raise ValueError(f"Preflight input missing for {request.id}")

        self._safe_start_structure(request)

        structure_dir = Path(output_dir) / request.id
        structure_dir.mkdir(parents=True, exist_ok=True)
        recorder = self._recorders.get(request.id)
        if recorder is None:
            recorder = AttemptRecorder(structure_dir, structure_id=request.id)
            self._recorders[request.id] = recorder

        outcome = Pass1Outcome(
            structure_id=request.id,
            role=request.role,
            kind=request.kind,
            status="degraded",
        )

        warmup_input = Path(preflight.input_xyz)
        warmup_xyz, warmup_payload = self._maybe_warmup(
            request,
            warmup_input,
            structure_dir,
            recorder,
        )
        self._apply_warmup_payload(outcome, warmup_payload)

        opt_outcome, _opt_spec, initial_hessian = self._run_primary_opt(
            request,
            warmup_xyz,
            recorder,
        )
        outcome.opt_initial_hessian = initial_hessian
        outcome.primary_attempt_id = opt_outcome["record"].attempt_id
        outcome.current_attempt_id = opt_outcome["record"].attempt_id
        if warmup_payload.get("history_entry") is not None:
            outcome.attempt_history.append(warmup_payload["history_entry"])
        outcome.attempt_history.append(opt_outcome["history_entry"])
        opt_result = opt_outcome["result"]
        outcome.opt_status = "complete" if opt_result.status == "complete" else "failed"
        outcome.opt_energy_hartree = opt_result.energy_hartree
        outcome.opt_error = opt_result.error
        outcome.last_geometry_path = (
            Path(opt_outcome["last_geometry_path"])
            if opt_outcome["last_geometry_path"] is not None
            else None
        )
        outcome.stop_geometry_path = (
            Path(opt_outcome["stop_geometry_path"])
            if opt_outcome.get("stop_geometry_path") is not None
            else None
        )
        outcome.optimization_monitor = opt_outcome.get("optimization_monitor")
        outcome.optimization_monitor_error = opt_outcome.get("optimization_monitor_error")

        legacy_opt_dir = structure_dir / "opt"
        self._populate_legacy_attempt_dir(legacy_opt_dir, opt_outcome, xyz_target_name="opt.xyz")
        outcome.opt_xyz = self._latest_existing_path(
            legacy_opt_dir / "opt.xyz",
            Path(opt_result.output_xyz) if opt_result.output_xyz else None,
            outcome.last_geometry_path,
        )
        outcome.opt_output = self._latest_existing_path(
            legacy_opt_dir / "output.out",
            opt_outcome["record"].output_out_path,
        )

        frequency_outcome = None
        if outcome.opt_status == "complete":
            frequency_outcome, _frequency_spec = self._run_frequency(
                request,
                outcome.opt_xyz,
                recorder,
            )
        if frequency_outcome is None:
            outcome.frequency_status = "skipped" if outcome.opt_status != "complete" else "not_run"
        else:
            frequency_result = frequency_outcome["result"]
            outcome.current_attempt_id = frequency_outcome["record"].attempt_id
            outcome.attempt_history.append(frequency_outcome["history_entry"])
            outcome.frequency_status = "complete" if frequency_result.status == "complete" else "failed"
            legacy_freq_dir = structure_dir / "freq"
            self._populate_legacy_attempt_dir(legacy_freq_dir, frequency_outcome)
            outcome.frequency_output = self._latest_existing_path(
                legacy_freq_dir / "output.out",
                frequency_outcome["record"].output_out_path,
            )
            outcome.frequency_energy_hartree = frequency_result.energy_hartree
            outcome.frequencies_cm1 = [
                float(value) for value in (frequency_result.frequencies_cm1 or ())
            ]
            outcome.imaginary_frequencies_cm1 = [
                float(value) for value in outcome.frequencies_cm1 if float(value) < 0.0
            ]
            outcome.frequency_zero_point_energy_hartree = (
                frequency_result.zero_point_energy_hartree
            )
            outcome.frequency_thermal_energy_hartree = frequency_result.thermal_energy_hartree
            outcome.frequency_enthalpy_hartree = frequency_result.enthalpy_hartree
            outcome.frequency_gibbs_free_energy_hartree = (
                frequency_result.gibbs_free_energy_hartree
            )
            outcome.frequency_gibbs_correction_hartree = (
                frequency_result.gibbs_correction_hartree
            )
            outcome.frequency_hessian_path = self._resolve_hessian_path(frequency_outcome)
            if outcome.frequency_status != "complete" and frequency_result.error and outcome.opt_error is None:
                outcome.error = frequency_result.error

        (
            outcome.ts_classification,
            outcome.int_classification,
            outcome.minimum_classification,
            _alignment_score,
        ) = self._classify_structure(
            request,
            opt_xyz=outcome.opt_xyz if outcome.opt_status == "complete" else None,
            frequency_output=outcome.frequency_output,
            frequencies_cm1=outcome.frequencies_cm1,
        )

        if outcome.opt_status != "complete":
            outcome.error = outcome.opt_error or f"R0 {'OptTS' if request.kind == 'ts' else 'Opt'} failed"
        return outcome

    def _maybe_warmup(
        self,
        request: StructureRequest,
        input_xyz: Path,
        structure_dir: Path,
        recorder: AttemptRecorder,
    ) -> Tuple[Path, Dict[str, Any]]:
        if not self._warmup_applies(request):
            return input_xyz, {
                "warmup_status": "not_requested",
                "warmup_used": False,
                "warmup_xyz": None,
                "warmup_constraints": [],
                "warmup_max_cycles": None,
                "warmup_error": None,
                "history_entry": None,
            }

        constraint_source_xyz = (
            input_xyz
            if str(request.source_stage).upper() == "S3"
            else Path(request.original_seed_xyz or input_xyz)
        )
        constraints, details = self._build_bond_constraints(
            constraint_source_xyz,
            request.forming_bonds,
        )
        max_cycles = (
            self.profile.warmup_max_cycles_int
            if request.role == "intermediate"
            else self.profile.warmup_max_cycles_ts
        )
        warmup_spec = QCJobSpec(
            engine="orca",
            task="opt",
            method=self.profile.geometry_method,
            basis=self.profile.geometry_basis,
            aux_basis=self.profile.geometry_aux_basis,
            solvent=self.profile.solvent,
            solvent_model=self.profile.solvent_model,
            route=self.profile.route_minimum,
            route_extras="LooseOpt",
            max_cycles=max_cycles,
            charge=int(request.charge),
            multiplicity=int(request.multiplicity),
            bond_constraints=constraints,
            allow_unconverged_geometry=True,
            grid=self.profile.geometry_grid,
            scf=self.profile.geometry_scf,
            timeout=self.profile.timeout_seconds,
        )
        warmup_outcome = self._execute_attempt(
            request.id,
            recorder,
            kind="warmup_opt",
            retry_level=0,
            runner=run_optimization,
            spec=warmup_spec,
            input_xyz=input_xyz,
        )
        legacy_warmup_dir = structure_dir / "warmup"
        self._populate_legacy_attempt_dir(legacy_warmup_dir, warmup_outcome, xyz_target_name="opt.xyz")
        warmup_result = warmup_outcome["result"]
        warmup_status = warmup_result.status if warmup_result.status in {"complete", "partial"} else "failed"
        warmup_xyz = self._latest_existing_path(
            legacy_warmup_dir / "opt.xyz",
            Path(warmup_result.output_xyz) if warmup_result.output_xyz else None,
            Path(warmup_outcome["last_geometry_path"]) if warmup_outcome["last_geometry_path"] else None,
        )
        next_input = input_xyz
        warmup_used = False
        if warmup_status in {"complete", "partial"} and warmup_xyz is not None:
            next_input = Path(warmup_xyz)
            warmup_used = True

        return next_input, {
            "warmup_status": warmup_status,
            "warmup_used": warmup_used,
            "warmup_xyz": warmup_xyz if warmup_used else None,
            "warmup_constraints": details,
            "warmup_max_cycles": max_cycles,
            "warmup_error": warmup_result.error,
            "history_entry": warmup_outcome["history_entry"],
        }

    def _run_primary_opt(
        self,
        request: StructureRequest,
        warmup_xyz: Path,
        recorder: AttemptRecorder,
    ) -> Tuple[Dict[str, Any], QCJobSpec, str]:
        initial_hessian = self.profile.initial_hessian_for_role(request.role)
        is_ts = request.kind == "ts"
        opt_spec = QCJobSpec(
            engine="orca",
            task="opt_ts" if is_ts else "opt",
            method=self.profile.geometry_method,
            basis=self.profile.geometry_basis,
            aux_basis=self.profile.geometry_aux_basis,
            solvent=self.profile.solvent,
            solvent_model=self.profile.solvent_model,
            route=self.profile.route_ts if is_ts else self.profile.route_minimum,
            route_extras="",
            max_cycles=(
                self.profile.max_cycles_ts
                if is_ts
                else (
                    self.profile.max_cycles_intermediate
                    if request.role == "intermediate"
                    else self.profile.max_cycles_minimum
                )
            ),
            charge=int(request.charge),
            multiplicity=int(request.multiplicity),
            grid=self.profile.geometry_grid,
            scf=self.profile.geometry_scf,
            timeout=self.profile.timeout_seconds,
            initial_hessian=initial_hessian,
            trust=self.profile.ts_trust_radius if is_ts else None,
        )
        opt_outcome = self._execute_attempt(
            request.id,
            recorder,
            kind="opt_ts" if is_ts else "opt",
            retry_level=0,
            runner=run_optimization,
            spec=opt_spec,
            input_xyz=warmup_xyz,
        )
        return opt_outcome, opt_spec, initial_hessian

    def _run_frequency(
        self,
        request: StructureRequest,
        opt_xyz: Optional[Path],
        recorder: AttemptRecorder,
    ) -> Tuple[Optional[Dict[str, Any]], Optional[QCJobSpec]]:
        if opt_xyz is None or not Path(opt_xyz).exists():
            return None, None

        freq_spec = QCJobSpec(
            engine="orca",
            task="freq",
            method=self.profile.frequency_method,
            basis=self.profile.frequency_basis,
            aux_basis=self.profile.geometry_aux_basis,
            solvent=self.profile.solvent,
            solvent_model=self.profile.solvent_model,
            charge=int(request.charge),
            multiplicity=int(request.multiplicity),
            grid=self.profile.geometry_grid,
            scf=self.profile.geometry_scf,
            timeout=self.profile.timeout_seconds,
        )
        freq_outcome = self._execute_attempt(
            request.id,
            recorder,
            kind="freq",
            retry_level=0,
            runner=run_frequency,
            spec=freq_spec,
            input_xyz=Path(opt_xyz),
        )
        return freq_outcome, freq_spec

    def _build_optimization_spec(
        self,
        request: StructureRequest,
        *,
        task: Optional[str] = None,
        route: Optional[str] = None,
        route_extras: str = "",
        initial_hessian: str = "model",
        hessian_filename: str | None = None,
        ts_mode: int | None = None,
        recalc_hessian: int | None = None,
        trust: float | None = None,
        max_cycles: int | None = None,
        grid: str | None = None,
        scf: str | None = None,
        allow_unconverged_geometry: bool = False,
        bond_constraints: Tuple[Tuple[int, int, Optional[float]], ...] = (),
    ) -> QCJobSpec:
        is_ts = (task or ("opt_ts" if request.kind == "ts" else "opt")).lower() in {
            "opt_ts",
            "ts",
            "opt-ts",
        }
        return QCJobSpec(
            engine="orca",
            task="opt_ts" if is_ts else "opt",
            method=self.profile.geometry_method,
            basis=self.profile.geometry_basis,
            aux_basis=self.profile.geometry_aux_basis,
            solvent=self.profile.solvent,
            solvent_model=self.profile.solvent_model,
            route=route or (self.profile.route_ts if is_ts else self.profile.route_minimum),
            route_extras=route_extras,
            max_cycles=max_cycles or (
                self.profile.max_cycles_ts
                if is_ts
                else (
                    self.profile.max_cycles_intermediate
                    if request.role == "intermediate"
                    else self.profile.max_cycles_minimum
                )
            ),
            charge=int(request.charge),
            multiplicity=int(request.multiplicity),
            grid=grid if grid is not None else self.profile.geometry_grid,
            scf=scf if scf is not None else self.profile.geometry_scf,
            timeout=self.profile.timeout_seconds,
            allow_unconverged_geometry=allow_unconverged_geometry,
            bond_constraints=bond_constraints,
            initial_hessian=initial_hessian,
            hessian_filename=hessian_filename,
            ts_mode=ts_mode,
            recalc_hessian=recalc_hessian,
            trust=trust,
        )

    def _build_frequency_spec(
        self,
        request: StructureRequest,
        *,
        grid: str | None = None,
        scf: str | None = None,
    ) -> QCJobSpec:
        return QCJobSpec(
            engine="orca",
            task="freq",
            method=self.profile.frequency_method,
            basis=self.profile.frequency_basis,
            aux_basis=self.profile.geometry_aux_basis,
            solvent=self.profile.solvent,
            solvent_model=self.profile.solvent_model,
            charge=int(request.charge),
            multiplicity=int(request.multiplicity),
            grid=grid if grid is not None else self.profile.geometry_grid,
            scf=scf if scf is not None else self.profile.geometry_scf,
            timeout=self.profile.timeout_seconds,
        )

    def _build_single_point_spec(self, request: StructureRequest) -> QCJobSpec:
        return QCJobSpec(
            engine="orca",
            task="sp",
            method=self.profile.sp_method,
            basis=self.profile.sp_basis,
            aux_basis=self.profile.sp_aux_basis,
            solvent=self.profile.solvent,
            solvent_model=self.profile.solvent_model,
            charge=int(request.charge),
            multiplicity=int(request.multiplicity),
            timeout=self.profile.timeout_seconds,
        )

    def _classify_structure(
        self,
        request: StructureRequest,
        *,
        opt_xyz: Optional[Path],
        frequency_output: Optional[Path],
        frequencies_cm1: Sequence[float],
        hessian_path: Optional[Path] = None,
    ) -> Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]], Dict[str, Any], Optional[float]]:
        ts_classification: Optional[Dict[str, Any]] = None
        int_classification: Optional[Dict[str, Any]] = None
        alignment_score: Optional[float] = None

        if request.kind == "ts":
            normal_modes = self._build_ts_normal_mode_payload(
                opt_xyz,
                frequency_output,
                hessian_path=hessian_path,
                frequencies_cm1=frequencies_cm1,
            )
            if normal_modes is not None:
                alignment_score = self._mode_alignment_score(
                    np.asarray(normal_modes["coordinates"], dtype=float),
                    np.asarray(normal_modes["displacements"], dtype=float),
                    request.forming_bonds,
                )
            ts_classification = identity_module.classify_ts(
                frequencies_cm1=frequencies_cm1,
                forming_bonds=request.forming_bonds,
                normal_modes=normal_modes,
            )
            if ts_classification is not None and alignment_score is not None:
                ts_classification = dict(ts_classification)
                ts_classification["alignment_score"] = alignment_score

        if request.role == "intermediate":
            int_classification = identity_module.classify_int(
                opt_xyz=opt_xyz,
                forming_bonds=request.forming_bonds,
                precursor_ref=self._resolve_reference_path(request, "precursor"),
                product_ref=self._resolve_reference_path(request, "product"),
                atom_mapping=self._load_atom_mapping(request.atom_mapping),
                frequencies_cm1=frequencies_cm1,
            )

        minimum_classification = identity_module.classify_minimum(
            frequencies_cm1=frequencies_cm1,
            expected_role=request.role,
            rmsd_to_expected=None,
        )
        return ts_classification, int_classification, minimum_classification, alignment_score

    @staticmethod
    def _read_orca_final_energy(output_path: Optional[Path]) -> Optional[float]:
        if output_path is None or not Path(output_path).is_file():
            return None
        try:
            lines = Path(output_path).read_text(
                encoding="utf-8", errors="replace"
            ).splitlines()
        except OSError:
            return None
        for line in reversed(lines):
            if "FINAL SINGLE POINT ENERGY" in line:
                try:
                    return float(line.split()[-1])
                except (ValueError, IndexError):
                    return None
        return None

    def _opt_energy_for(self, outcome: Optional[Pass1Outcome]) -> Optional[float]:
        if outcome is None:
            return None
        energy = getattr(outcome, "opt_energy_hartree", None)
        if energy is not None:
            return float(energy)
        geometry = getattr(outcome, "canonical_xyz", None) or getattr(outcome, "opt_xyz", None)
        if geometry is None:
            return None
        geometry_path = Path(geometry)
        if not geometry_path.is_file():
            return None
        energy = self._read_orca_final_energy(geometry_path.with_name("output.out"))
        if energy is not None:
            return energy
        structure_dir = geometry_path.parents[1] if geometry_path.parent.name == "opt" else geometry_path.parent
        attempt_dir = structure_dir / "diagnostics"
        if attempt_dir.is_dir():
            for attempt in sorted(attempt_dir.glob("attempt_*_opt"), reverse=True):
                energy = self._read_orca_final_energy(attempt / "output.out")
                if energy is not None:
                    return energy
        return None

    def _int_identity_thresholds(self) -> Dict[str, Any]:
        return dict(
            (self.config or {})
            .get("refinement", {})
            .get("common", {})
            .get("identity", {})
            .get("int_identity_v2", {})
            or {}
        )

    def _classify_int_v2(
        self,
        request: StructureRequest,
        *,
        opt_xyz: Path,
        frequencies_cm1: Sequence[float],
        energy_int_hartree: Optional[float],
        converged: Optional[bool] = None,
    ) -> Dict[str, Any]:
        variant_key = request.variant_id or request.id
        ts_outcome = self._variant_ts_outcome.get(variant_key)
        product_outcome = self._variant_product_outcome.get(variant_key)
        ts_xyz = None
        if ts_outcome is not None:
            ts_xyz = ts_outcome.canonical_xyz or ts_outcome.opt_xyz
        product_xyz = None
        if product_outcome is not None:
            product_xyz = product_outcome.canonical_xyz or product_outcome.opt_xyz
        if product_xyz is None:
            product_xyz = self._resolve_reference_path(request, "product")
        if ts_xyz is None:
            ts_xyz = self._resolve_reference_path(request, "ts")
        result = identity_module.classify_int_v2(
            opt_xyz=Path(opt_xyz),
            product_ref=Path(product_xyz) if product_xyz is not None else None,
            ts_ref=Path(ts_xyz) if ts_xyz is not None else None,
            forming_bonds=request.forming_bonds,
            frequencies_cm1=frequencies_cm1,
            energy_int_hartree=energy_int_hartree,
            energy_ts_hartree=self._opt_energy_for(ts_outcome),
            energy_product_hartree=self._opt_energy_for(product_outcome),
            atom_mapping=self._load_atom_mapping(request.atom_mapping),
            converged=converged,
            thresholds=self._int_identity_thresholds(),
        )
        return result.as_dict()

    def _build_ts_normal_mode_payload(
        self,
        opt_xyz: Optional[Path],
        frequency_output: Optional[Path],
        *,
        hessian_path: Optional[Path] = None,
        frequencies_cm1: Sequence[float] = (),
    ) -> Optional[Dict[str, Any]]:
        if opt_xyz is None:
            return None
        try:
            coordinates, symbols = read_xyz(Path(opt_xyz))
            if hessian_path is not None:
                hessian_frequencies, mode_matrix = self._parse_orca_hessian_vibrations(
                    Path(hessian_path)
                )
                frequencies = list(frequencies_cm1) or hessian_frequencies
                imaginary_indices = [
                    index for index, value in enumerate(frequencies) if float(value) < 0.0
                ]
                if mode_matrix is None or len(imaginary_indices) != 1:
                    return None
                mode_index = imaginary_indices[0]
                if mode_index >= mode_matrix.shape[1]:
                    return None
                displacement = mode_matrix[:, mode_index].reshape((-1, 3))
            else:
                if frequency_output is None:
                    return None
                content = Path(frequency_output).read_text(encoding="utf-8", errors="ignore")
                symbol_list: List[str] = [str(symbol) for symbol in symbols]
                displacement_vectors = _parse_orca_displacement_vectors(content)
                if not displacement_vectors:
                    displacement_vectors = parse_orca6_normal_mode_vectors(
                        content, symbol_list
                    )
                if not displacement_vectors:
                    return None
                displacement = displacement_vectors[sorted(displacement_vectors)[0]]
        except (OSError, ValueError):
            return None
        if len(coordinates) != len(displacement):
            return None
        # ORCA Hessian ``$normal_modes`` are mass-weighted; deweight and
        # renormalise so alignment scores stay Cartesian-comparable.
        if hessian_path is not None:
            deweighted = self._demass_weight_mode_vectors(
                np.asarray(displacement, dtype=float), symbols
            )
            if deweighted is None:
                return None
            displacement = deweighted
        return {
            "coordinates": np.asarray(coordinates, dtype=float).tolist(),
            "displacements": np.asarray(displacement, dtype=float).tolist(),
        }

    @staticmethod
    def _demass_weight_mode_vectors(
        displacements: np.ndarray,
        symbols: Sequence[str],
    ) -> Optional[np.ndarray]:
        if displacements.ndim != 2 or displacements.shape[0] != len(symbols):
            return None
        masses = np.asarray(
            [_ELEMENT_MASS.get(str(symbol).strip().capitalize(), 12.011) for symbol in symbols],
            dtype=float,
        )
        deweighted = displacements * np.sqrt(masses)[:, None]
        norm = float(np.linalg.norm(deweighted))
        if norm <= 1e-12:
            return None
        return deweighted / norm

    @staticmethod
    def _recover_normal_modes_from_hessian_section(
        content: str,
    ) -> Optional[np.ndarray]:
        """Recover mass-weighted normal modes from the force-constant matrix.

        Some ORCA 6.1.1 builds write an all-zero ``$normal_modes`` block in the
        ``.hess`` file while ``$hessian`` still holds the complete force
        constant matrix.  Mass-weighted diagonalisation recovers the mode
        vectors at zero recalculation cost.  Columns are mass-weighted (ORCA
        ``$normal_modes`` semantics, i.e. each column is a normalised vector of
        M^-1/2-weighted displacements) and sorted ascending by eigenvalue, so
        column ``i`` corresponds to frequency index ``i`` in
        ``$vibrational_frequencies``.
        """
        atoms_section = re.search(r"(?ms)^\$atoms\s*$\n(.*?)(?=^\$|\Z)", content)
        hessian_section = re.search(r"(?ms)^\$hessian\s*$\n(.*?)(?=^\$|\Z)", content)
        if atoms_section is None or hessian_section is None:
            return None
        try:
            atom_lines = atoms_section.group(1).splitlines()
            n_atoms = int(atom_lines[0].strip())
            masses = np.asarray(
                [float(line.split()[1]) for line in atom_lines[1 : 1 + n_atoms]],
                dtype=float,
            )
        except (ValueError, IndexError):
            return None
        hessian_lines = hessian_section.group(1).splitlines()
        try:
            size = int(hessian_lines[0].split()[0])
        except (ValueError, IndexError):
            return None
        if size != 3 * n_atoms or size <= 0:
            return None
        matrix = np.zeros((size, size), dtype=float)
        active_columns: List[int] = []
        for line in hessian_lines[1:]:
            tokens = line.split()
            if not tokens:
                continue
            if all(re.fullmatch(r"[+-]?\d+", token) for token in tokens):
                active_columns = [int(token) for token in tokens]
                continue
            if not active_columns or len(tokens) != len(active_columns) + 1:
                continue
            try:
                row = int(tokens[0])
                values = [float(token.replace("D", "E")) for token in tokens[1:]]
            except ValueError:
                continue
            if 0 <= row < size:
                for column, value in zip(active_columns, values):
                    if 0 <= column < size:
                        matrix[row, column] = value
        if not np.any(matrix):
            return None
        matrix = (matrix + matrix.T) / 2.0
        inv_sqrt_m = np.repeat(1.0 / np.sqrt(masses), 3)
        mass_weighted = (
            matrix * inv_sqrt_m[:, None] * inv_sqrt_m[None, :]
        )
        try:
            _eigenvalues, eigenvectors = np.linalg.eigh(mass_weighted)
        except np.linalg.LinAlgError:
            return None
        if eigenvectors.shape[1] != size:
            return None
        return eigenvectors

    @staticmethod
    def _parse_orca_hessian_vibrations(
        hessian_path: Path,
    ) -> Tuple[List[float], Optional[np.ndarray]]:
        """Read frequencies and normal modes written by ORCA ``CalcAll``.

        ORCA's Hessian file contains the final optimization Hessian, including
        ``$vibrational_frequencies`` and ``$normal_modes``.  This is the same
        data a separate numerical frequency job would consume, so reuse it for
        stationary-point validation instead of submitting a redundant FREQ.
        """
        content = Path(hessian_path).read_text(encoding="utf-8", errors="ignore")

        def section(name: str) -> List[str]:
            match = re.search(
                rf"(?ms)^\${re.escape(name)}\s*$\n(.*?)(?=^\$|\Z)", content
            )
            return match.group(1).splitlines() if match else []

        frequency_lines = section("vibrational_frequencies")
        frequencies_by_index: Dict[int, float] = {}
        for line in frequency_lines[1:]:
            tokens = line.split()
            if len(tokens) < 2:
                continue
            try:
                frequencies_by_index[int(tokens[0])] = float(tokens[1].replace("D", "E"))
            except ValueError:
                continue
        frequencies = [frequencies_by_index[index] for index in sorted(frequencies_by_index)]

        mode_lines = section("normal_modes")
        if not mode_lines:
            return frequencies, None
        try:
            dimensions = [int(token) for token in mode_lines[0].split()[:2]]
            rows, columns = dimensions
        except (TypeError, ValueError):
            return frequencies, None
        if rows <= 0 or columns <= 0:
            return frequencies, None
        matrix = np.zeros((rows, columns), dtype=float)
        active_columns: List[int] = []
        for line in mode_lines[1:]:
            tokens = line.split()
            if not tokens:
                continue
            if all(re.fullmatch(r"[+-]?\d+", token) for token in tokens):
                active_columns = [int(token) for token in tokens]
                continue
            if not active_columns or len(tokens) != len(active_columns) + 1:
                continue
            try:
                row = int(tokens[0])
                values = [float(token.replace("D", "E")) for token in tokens[1:]]
            except ValueError:
                continue
            if 0 <= row < rows:
                for column, value in zip(active_columns, values):
                    if 0 <= column < columns:
                        matrix[row, column] = value
        if not np.any(matrix):
            recovered = RefinementEngine._recover_normal_modes_from_hessian_section(
                content
            )
            if recovered is not None:
                matrix = recovered
                rows, columns = matrix.shape
            else:
                return frequencies, None
        return frequencies, matrix

    @staticmethod
    def _mode_alignment_score(
        coordinates: np.ndarray,
        displacements: np.ndarray,
        forming_bonds: Sequence[Tuple[int, int]],
    ) -> Optional[float]:
        if coordinates.shape != displacements.shape or coordinates.size == 0:
            return None
        projections: List[float] = []
        for atom_i, atom_j in forming_bonds:
            if atom_i < 0 or atom_j < 0 or atom_i >= len(coordinates) or atom_j >= len(coordinates):
                return None
            bond_vector = coordinates[int(atom_j)] - coordinates[int(atom_i)]
            norm = float(np.linalg.norm(bond_vector))
            if norm <= 1e-12:
                return None
            relative_displacement = displacements[int(atom_j)] - displacements[int(atom_i)]
            projections.append(
                abs(float(np.dot(relative_displacement, bond_vector / norm)))
            )
        if not projections:
            return None
        return float(min(projections))

    @staticmethod
    def _resolve_hessian_path(attempt_outcome: Dict[str, Any] | None) -> Optional[Path]:
        if attempt_outcome is None:
            return None
        result = attempt_outcome.get("result")
        record = attempt_outcome.get("record")
        candidate_output = Path(result.output_file) if result and result.output_file else None
        candidates: List[Path] = []
        if candidate_output is not None:
            candidates.append(candidate_output.with_suffix(".hess"))
        if record is not None:
            candidates.extend(sorted(Path(record.directory).glob("*.hess")))
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return None

    @staticmethod
    def _target_reaction_displacement(
        coordinates: np.ndarray,
        forming_bonds: Sequence[Tuple[int, int]],
    ) -> Optional[np.ndarray]:
        """Build a normalized displacement that shortens every forming bond."""
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            return None
        target = np.zeros_like(coordinates, dtype=float)
        for atom_i, atom_j in forming_bonds:
            if atom_i < 0 or atom_j < 0 or atom_i >= len(coordinates) or atom_j >= len(coordinates):
                return None
            axis = coordinates[int(atom_j)] - coordinates[int(atom_i)]
            length = float(np.linalg.norm(axis))
            if length <= 1.0e-12:
                return None
            axis /= length
            target[int(atom_i)] += axis
            target[int(atom_j)] -= axis
        norm = float(np.linalg.norm(target))
        return target / norm if norm > 1.0e-12 else None

    def _select_target_ts_mode(
        self,
        request: StructureRequest,
        coordinates: np.ndarray,
        displacement_vectors: Mapping[int, np.ndarray],
    ) -> Dict[str, Any]:
        """Select the imaginary mode with the strongest forming-bond overlap.

        ``TS_Mode`` accepts a mode ordinal, not a chemistry label.  Selecting
        it here prevents the non-convergence ladder from treating ORCA's first
        available imaginary mode as the target reaction coordinate by default.
        """
        target = self._target_reaction_displacement(coordinates, request.forming_bonds)
        if target is None:
            return {"selection_status": "unavailable", "reason": "forming_bonds_unavailable", "mode_scores": []}

        scores: List[Dict[str, float | int]] = []
        for mode_index, displacement in displacement_vectors.items():
            vector = np.asarray(displacement, dtype=float)
            if vector.shape != target.shape:
                continue
            norm = float(np.linalg.norm(vector))
            if norm <= 1.0e-12:
                continue
            overlap = abs(float(np.vdot(target.reshape(-1), (vector / norm).reshape(-1))))
            scores.append({"mode_index": int(mode_index), "overlap": overlap})

        scores.sort(key=lambda item: float(item["overlap"]), reverse=True)
        if not scores:
            return {"selection_status": "unavailable", "reason": "no_parseable_imaginary_modes", "mode_scores": []}

        best = scores[0]
        best_overlap = float(best["overlap"])
        margin = best_overlap - float(scores[1]["overlap"]) if len(scores) > 1 else best_overlap
        if best_overlap < self.profile.ts_nonconvergence_mode_min_overlap:
            reason = "target_overlap_below_threshold"
        elif len(scores) > 1 and margin < self.profile.ts_nonconvergence_mode_min_overlap_margin:
            reason = "target_mode_ambiguous"
        else:
            return {
                "selection_status": "selected",
                "selected_mode": int(best["mode_index"]),
                "selected_overlap": best_overlap,
                "overlap_margin": margin,
                "mode_scores": scores,
            }
        return {
            "selection_status": "ambiguous",
            "reason": reason,
            "selected_overlap": best_overlap,
            "overlap_margin": margin,
            "mode_scores": scores,
        }

    def _run_local_hessian_analysis(
        self,
        request: StructureRequest,
        input_xyz: Path,
        recorder: AttemptRecorder,
        *,
        level: int,
        kind: str,
    ) -> Dict[str, Any]:
        """Calculate a checkpoint Hessian; TS requests also select a target mode."""
        frequency_spec = self._build_frequency_spec(request)
        frequency_outcome = self._execute_attempt(
            request.id,
            recorder,
            kind=f"{kind}_mode_analysis",
            retry_level=level,
            runner=run_frequency,
            spec=frequency_spec,
            input_xyz=Path(input_xyz),
        )
        result: QCJobResult = frequency_outcome["result"]
        output_path = self._latest_existing_path(
            Path(result.output_file) if result.output_file else None,
            Path(frequency_outcome["record"].output_out_path),
        )
        hessian_path = self._resolve_hessian_path(frequency_outcome)
        frequencies = [float(value) for value in (result.frequencies_cm1 or ())]
        analysis: Dict[str, Any] = {
            "attempt_id": frequency_outcome["record"].attempt_id,
            "strategy": f"{kind}_mode_analysis",
            "level": level,
            "status": result.status,
            "input_geometry": str(input_xyz),
            "frequency_output": str(output_path) if output_path else None,
            "hessian_filename": str(hessian_path) if hessian_path else None,
            "frequencies_cm1": frequencies,
            "imaginary_frequencies_cm1": [value for value in frequencies if value < 0.0],
        }
        if result.status != "complete" or output_path is None or hessian_path is None:
            analysis.update({"selection_status": "unavailable", "reason": "local_hessian_unavailable", "mode_scores": []})
            return analysis

        if request.kind != "ts":
            analysis.update(
                {
                    "selection_status": "not_applicable",
                    "reason": "minimum_optimization_does_not_use_ts_mode",
                    "mode_scores": [],
                }
            )
            return analysis

        try:
            coordinates, _symbols = read_xyz(Path(input_xyz))
            displacement_vectors = _parse_orca_displacement_vectors(
                Path(output_path).read_text(encoding="utf-8", errors="ignore")
            )
        except (OSError, ValueError):
            analysis.update({"selection_status": "unavailable", "reason": "local_mode_parse_failed", "mode_scores": []})
            return analysis
        analysis.update(self._select_target_ts_mode(request, coordinates, displacement_vectors))
        return analysis

    def _rescue_matrix_enabled(self) -> bool:
        """Master switch for the S3.4 rescue matrix (refinement.common.rescue.matrix.enabled)."""
        rescue_cfg = dict(dict(self.config.get("refinement", {}) or {}).get("common", {}).get("rescue", {}) or {})
        matrix_cfg = dict(rescue_cfg.get("matrix", {}) or {})
        return bool(matrix_cfg.get("enabled", True))

    def _run_pass2_rescue(
        self,
        pass1_outcomes: Sequence[Pass1Outcome],
        output_dir: Path,
    ) -> List[Pass1Outcome]:
        """Run the S3.4 rescue matrix for non-viable outcomes.

        Every outcome is classified (F1-F6), its matrix cell looked up and
        the cell methods executed in order until a valid candidate emerges.
        F5 (SCF) and F6 (crash/timeout) exit directly without rescue — each
        incident is investigated individually.
        """
        outcomes = [copy.deepcopy(outcome) for outcome in pass1_outcomes]
        matrix_enabled = self._rescue_matrix_enabled()
        for outcome in outcomes:
            request = self._requests_by_id[outcome.structure_id]
            recorder = self._recorders[outcome.structure_id]
            failure = self._classify_rescue_failure(request, outcome)
            self._record_rescue_diagnosis(outcome, failure)
            if failure is None:
                continue
            if failure in (FailureType.F5_SCF_NOT_CONVERGED, FailureType.F6_CRASH_TIMEOUT_OTHER):
                logger.info(
                    "[%s.4] %s diagnosis=%s -> direct failure exit (no rescue)",
                    self.profile.stage,
                    outcome.structure_id,
                    failure.value,
                )
                continue
            if not matrix_enabled:
                logger.info(
                    "[%s.4] %s diagnosis=%s -> rescue disabled (primary-only mode)",
                    self.profile.stage,
                    outcome.structure_id,
                    failure.value,
                )
                continue
            plan = lookup_plan(
                failure, StructureKind.from_request(request.kind, request.role)
            )
            if plan is None:
                continue
            logger.info(
                "[%s.4] %s diagnosis=%s plan=%s",
                self.profile.stage,
                outcome.structure_id,
                failure.value,
                [method.value for method in plan.methods],
            )
            self._run_rescue_cell(request, outcome, recorder, Path(output_dir), plan)
        return outcomes

    def _classify_rescue_failure(
        self,
        request: StructureRequest,
        outcome: Pass1Outcome,
    ) -> Optional[FailureType]:
        """Classify a non-viable outcome into the F1-F6 failure taxonomy."""
        if outcome.opt_status == "complete":
            if (
                request.kind != "ts"
                and outcome.imaginary_frequencies_cm1
                and min(outcome.imaginary_frequencies_cm1) < _MINIMUM_IMAGINARY_CUTOFF_CM1
            ):
                return FailureType.F4_MINIMUM_WITH_IMAGINARY
            if (
                request.role == "intermediate"
                and str(
                    (outcome.int_classification or {}).get("identity") or ""
                )
                == "collapsed_to_product"
            ):
                return FailureType.F7_COLLAPSED_TO_PRODUCT
            return None
        history = outcome.attempt_history[-1] if outcome.attempt_history else {}
        failure_type = str(history.get("failure_type") or "").lower()
        stop_reason = str(history.get("stop_reason") or "").lower()
        error_text = str(outcome.opt_error or "").lower()
        if "scf" in failure_type or "scf" in error_text:
            return FailureType.F5_SCF_NOT_CONVERGED
        if (
            bool(history.get("timed_out"))
            or "timed_out" in stop_reason
            or any(marker in error_text for marker in ("timed out", "crash", "killed"))
        ):
            return FailureType.F6_CRASH_TIMEOUT_OTHER
        if not self._is_geometry_nonconvergence(outcome):
            return FailureType.F6_CRASH_TIMEOUT_OTHER
        if request.kind == "ts":
            saddle_counts = self._parse_saddle_eigenvalue_counts(outcome)
            if saddle_counts:
                last = saddle_counts[-1]
                if last >= 2:
                    return FailureType.F2_HIGHER_ORDER_SADDLE
                if last == 0:
                    return FailureType.F3_TS_NO_IMAGINARY_MODE
        return FailureType.F1_GEOMETRY_NOT_CONVERGED

    @staticmethod
    def _parse_saddle_eigenvalue_counts(outcome: Pass1Outcome) -> List[int]:
        """Count ``Hessian has N negative eigenvalues`` at each recalc point."""
        output = outcome.opt_output
        if output is None or not Path(output).exists():
            return []
        try:
            text = Path(output).read_text(encoding="utf-8", errors="ignore")
        except OSError:
            return []
        counts: List[int] = []
        for match in re.finditer(r"Hessian has\s+(\d+)\s+negative eigenvalues", text):
            counts.append(int(match.group(1)))
        return counts

    def _record_rescue_diagnosis(
        self,
        outcome: Pass1Outcome,
        failure: Optional[FailureType],
    ) -> None:
        if failure is None:
            return
        record: Dict[str, Any] = {
            "strategy": "S3.4.0_diagnosis",
            "failure_type": failure.value,
            "failure_label": failure.label,
            "status": "diagnosed",
        }
        if failure is FailureType.F4_MINIMUM_WITH_IMAGINARY and outcome.imaginary_frequencies_cm1:
            record["imaginary_frequencies_cm1"] = list(outcome.imaginary_frequencies_cm1)
        if failure in (
            FailureType.F2_HIGHER_ORDER_SADDLE,
            FailureType.F3_TS_NO_IMAGINARY_MODE,
        ):
            saddle_counts = self._parse_saddle_eigenvalue_counts(outcome)
            if saddle_counts:
                record["saddle_eigenvalue_history"] = saddle_counts
        outcome.pass2_rescue_attempts.append(record)
        if failure in (FailureType.F5_SCF_NOT_CONVERGED, FailureType.F6_CRASH_TIMEOUT_OTHER):
            if not outcome.error:
                outcome.error = f"Direct failure exit: {failure.label} (no rescue)"

    def _rescue_trigger(self, outcome: Pass1Outcome) -> Dict[str, Any]:
        history = outcome.attempt_history[-1] if outcome.attempt_history else {}
        return {
            "failure_type": str(history.get("failure_type") or "unknown"),
            "stop_reason": str(history.get("stop_reason") or "unknown"),
        }

    @staticmethod
    def _rescue_input_geometry(outcome: Pass1Outcome) -> Optional[Path]:
        for candidate in (outcome.stop_geometry_path, outcome.last_geometry_path):
            if candidate is not None and Path(candidate).exists():
                return Path(candidate)
        return None

    def _run_rescue_cell(
        self,
        request: StructureRequest,
        outcome: Pass1Outcome,
        recorder: AttemptRecorder,
        structure_dir: Path,
        plan: RescuePlan,
    ) -> None:
        """Execute a matrix cell: ordered methods until a valid candidate."""
        for method in plan.methods:
            params = methods_params(self.config, method)
            logger.info(
                "[%s.4] %s executing method=%s (family=%s)",
                self.profile.stage,
                outcome.structure_id,
                method.value,
                method.family.value,
            )
            attempt = self._run_rescue_method(
                request, outcome, recorder, structure_dir, method, params
            )
            if attempt is not None and self._candidate_is_valid_for_request(request, attempt):
                outcome.error = None
                logger.info(
                    "[%s.4] %s method=%s produced valid candidate",
                    self.profile.stage,
                    outcome.structure_id,
                    method.value,
                )
                return

    def _run_rescue_method(
        self,
        request: StructureRequest,
        outcome: Pass1Outcome,
        recorder: AttemptRecorder,
        structure_dir: Path,
        method: RescueMethod,
        params: MethodParams,
    ) -> Optional[Dict[str, Any]]:
        try:
            if method is RescueMethod.FRESH_HESSIAN_RESTART:
                return self._rescue_fresh_hessian(request, outcome, recorder, params, monitor=False)
            if method is RescueMethod.FRESH_HESSIAN_MODE_MONITOR:
                return self._rescue_fresh_hessian(request, outcome, recorder, params, monitor=True)
            if method is RescueMethod.TS_MODE_DIRECTED:
                return self._rescue_ts_mode_directed(request, outcome, recorder, params)
            if method is RescueMethod.MODE_DISPLACEMENT:
                return self._rescue_mode_displacement(
                    request, outcome, recorder, structure_dir, params
                )
            if method is RescueMethod.SADDLE_BREAK:
                return self._rescue_saddle_break(
                    request, outcome, recorder, structure_dir, params
                )
            if method is RescueMethod.CALCALL_OPT:
                return self._rescue_calcall(request, outcome, recorder, params)
            if method is RescueMethod.IRC_MIDPOINT_RECOVERY:
                return self._rescue_int_via_irc(request, outcome, recorder, structure_dir, params)
        except Exception as exc:
            logger.warning(
                "%s rescue method %s failed for %s: %s",
                self.profile.stage,
                method.value,
                request.id,
                exc,
            )
        return None

    def _append_rescue_skipped(
        self,
        outcome: Pass1Outcome,
        family: str,
        method: str,
        reason: str,
        extra: Optional[Dict[str, Any]] = None,
    ) -> None:
        record: Dict[str, Any] = {
            "strategy": method,
            "method_family": family,
            "method": method,
            "status": "skipped",
            "error": reason,
            "trigger": self._rescue_trigger(outcome),
        }
        if extra:
            record.update(extra)
        outcome.pass2_rescue_attempts.append(record)

    def _rescue_fresh_hessian(
        self,
        request: StructureRequest,
        outcome: Pass1Outcome,
        recorder: AttemptRecorder,
        params: MethodParams,
        *,
        monitor: bool,
    ) -> Optional[Dict[str, Any]]:
        input_xyz = self._rescue_input_geometry(outcome)
        if input_xyz is None:
            self._append_rescue_skipped(
                outcome, "R1", "fresh_hessian_restart", "no_usable_stop_geometry"
            )
            return None
        is_ts = request.kind == "ts"
        method_name = "fresh_hessian_mode_monitor" if monitor else "fresh_hessian_restart"
        attempt = self._run_rescue_attempt(
            request,
            recorder,
            strategy=method_name,
            level=1,
            kind=f"r1_{method_name}",
            input_xyz=input_xyz,
            opt_spec=self._build_optimization_spec(
                request,
                task="opt_ts" if is_ts else "opt",
                route=self.profile.route_ts if is_ts else self.profile.route_minimum,
                initial_hessian="calculate",
                recalc_hessian=params.recalc_hessian,
                max_cycles=params.max_cycles,
            ),
        )
        attempt.update(
            {
                "method_family": "R1",
                "method": method_name,
                "trigger": self._rescue_trigger(outcome),
                "use_default_trust": True,
            }
        )
        outcome.pass2_rescue_attempts.append(attempt)
        return attempt

    def _rescue_ts_mode_directed(
        self,
        request: StructureRequest,
        outcome: Pass1Outcome,
        recorder: AttemptRecorder,
        params: MethodParams,
    ) -> Optional[Dict[str, Any]]:
        input_xyz = self._rescue_input_geometry(outcome)
        if input_xyz is None:
            self._append_rescue_skipped(
                outcome, "R2", "ts_mode_directed", "no_usable_stop_geometry"
            )
            return None
        analysis = self._run_local_hessian_analysis(
            request, input_xyz, recorder, level=1, kind="ts_mode_directed"
        )
        mode = analysis.get("selected_mode")
        hessian = analysis.get("hessian_filename")
        if analysis.get("selection_status") != "selected" or mode is None or not hessian:
            outcome.pass2_rescue_attempts.append(
                {
                    "strategy": "ts_mode_directed",
                    "method_family": "R2",
                    "method": "ts_mode_directed",
                    "status": "skipped",
                    "error": "No unambiguous target imaginary mode for mode-directed OptTS",
                    "trigger": self._rescue_trigger(outcome),
                    "local_hessian_analysis": analysis,
                    "mode_analysis": analysis,
                }
            )
            return None
        attempt = self._run_rescue_attempt(
            request,
            recorder,
            strategy="ts_mode_directed",
            level=1,
            kind="r2_ts_mode_directed",
            input_xyz=input_xyz,
            opt_spec=self._build_optimization_spec(
                request,
                task="opt_ts",
                route=self.profile.route_ts,
                initial_hessian="read",
                hessian_filename=str(hessian),
                ts_mode=int(mode),
                trust=params.trust,
                max_cycles=params.max_cycles,
            ),
        )
        attempt.update(
            {
                "method_family": "R2",
                "method": "ts_mode_directed",
                "trigger": self._rescue_trigger(outcome),
                "mode_analysis": analysis,
                "local_hessian_analysis": analysis,
            }
        )
        outcome.pass2_rescue_attempts.append(attempt)
        return attempt

    def _rescue_mode_displacement(
        self,
        request: StructureRequest,
        outcome: Pass1Outcome,
        recorder: AttemptRecorder,
        structure_dir: Path,
        params: MethodParams,
    ) -> Optional[Dict[str, Any]]:
        source = outcome.opt_xyz or self._rescue_input_geometry(outcome)
        if source is None or not Path(source).exists():
            self._append_rescue_skipped(
                outcome, "R2", "mode_displacement", "no_usable_geometry"
            )
            return None
        freq_output = outcome.frequency_output
        if freq_output is None or not Path(freq_output).exists():
            analysis = self._run_local_hessian_analysis(
                request, Path(source), recorder, level=1, kind="mode_displacement"
            )
            freq_output = analysis.get("frequency_output")
        if freq_output is None or not Path(freq_output).exists():
            self._append_rescue_skipped(
                outcome, "R2", "mode_displacement", "no_imaginary_mode_output"
            )
            return None
        best: Optional[Dict[str, Any]] = None
        for sign, label in ((1.0, "plus"), (-1.0, "minus")):
            try:
                displaced = self._build_mode_displaced_xyz(
                    Path(source),
                    Path(freq_output),
                    structure_dir / f"r2_mode_displacement_{label}.xyz",
                    sign=sign,
                )
            except ValueError as exc:
                logger.debug("mode displacement %s skipped for %s: %s", label, request.id, exc)
                continue
            attempt = self._run_rescue_attempt(
                request,
                recorder,
                strategy=f"mode_displacement_{label}",
                level=1,
                kind=f"r2_mode_displacement_{label}",
                input_xyz=displaced,
                opt_spec=self._build_optimization_spec(
                    request,
                    task="opt",
                    route=self.profile.route_minimum,
                    route_extras="TightOpt",
                    initial_hessian="calculate",
                    grid=self._stricter_grid(),
                    scf=self._stricter_scf(),
                    max_cycles=params.max_cycles,
                ),
            )
            attempt.update(
                {
                    "method_family": "R2",
                    "method": "mode_displacement",
                    "displacement_sign": label,
                    "trigger": self._rescue_trigger(outcome),
                }
            )
            outcome.pass2_rescue_attempts.append(attempt)
            if best is None:
                best = attempt
        return best

    def _rescue_saddle_break(
        self,
        request: StructureRequest,
        outcome: Pass1Outcome,
        recorder: AttemptRecorder,
        structure_dir: Path,
        params: MethodParams,
    ) -> Optional[Dict[str, Any]]:
        input_xyz = self._rescue_input_geometry(outcome)
        if input_xyz is None:
            self._append_rescue_skipped(
                outcome, "R2", "saddle_break", "no_usable_stop_geometry"
            )
            return None
        analysis = self._run_local_hessian_analysis(
            request, input_xyz, recorder, level=1, kind="saddle_break"
        )
        freq_output = analysis.get("frequency_output")
        if freq_output is None or not Path(freq_output).exists():
            self._append_rescue_skipped(
                outcome, "R2", "saddle_break", "local_hessian_unavailable"
            )
            return None
        best: Optional[Dict[str, Any]] = None
        for sign, label in ((1.0, "plus"), (-1.0, "minus")):
            try:
                displaced = self._build_mode_displaced_xyz(
                    input_xyz,
                    Path(freq_output),
                    structure_dir / f"r2_saddle_break_{label}.xyz",
                    sign=sign,
                    mode_index=1,
                )
            except ValueError as exc:
                logger.debug("saddle break %s skipped for %s: %s", label, request.id, exc)
                continue
            attempt = self._run_rescue_attempt(
                request,
                recorder,
                strategy=f"saddle_break_{label}",
                level=1,
                kind=f"r2_saddle_break_{label}",
                input_xyz=displaced,
                opt_spec=self._build_optimization_spec(
                    request,
                    task="opt_ts",
                    route=self.profile.route_ts,
                    initial_hessian="calculate",
                    max_cycles=params.max_cycles,
                ),
            )
            attempt.update(
                {
                    "method_family": "R2",
                    "method": "saddle_break",
                    "displacement_sign": label,
                    "mode_index": 1,
                    "trigger": self._rescue_trigger(outcome),
                }
            )
            outcome.pass2_rescue_attempts.append(attempt)
            if best is None:
                best = attempt
        return best

    def _rescue_calcall(
        self,
        request: StructureRequest,
        outcome: Pass1Outcome,
        recorder: AttemptRecorder,
        params: MethodParams,
    ) -> Optional[Dict[str, Any]]:
        input_xyz = self._rescue_input_geometry(outcome)
        if input_xyz is None:
            self._append_rescue_skipped(
                outcome, "R3", "calcall_opt", "no_usable_stop_geometry"
            )
            return None
        is_ts = request.kind == "ts"
        ts_mode: Optional[int] = None
        hessian: Optional[str] = None
        mode_analysis: Optional[Dict[str, Any]] = None
        if is_ts:
            analysis = self._run_local_hessian_analysis(
                request, input_xyz, recorder, level=2, kind="calcall"
            )
            mode_analysis = analysis
            if (
                analysis.get("selection_status") == "selected"
                and analysis.get("selected_mode") is not None
                and analysis.get("hessian_filename")
            ):
                ts_mode = int(analysis["selected_mode"])
                hessian = str(analysis["hessian_filename"])
        attempt = self._run_rescue_attempt(
            request,
            recorder,
            strategy="calcall_opt",
            level=2,
            kind="r3_calcall",
            input_xyz=input_xyz,
            opt_spec=self._build_optimization_spec(
                request,
                task="opt_ts" if is_ts else "opt",
                route=self.profile.route_ts if is_ts else self.profile.route_minimum,
                initial_hessian="read" if hessian else "calculate",
                hessian_filename=hessian,
                ts_mode=ts_mode,
                recalc_hessian=params.recalc_hessian,
                max_cycles=params.max_cycles,
            ),
            reuse_final_hessian=True,
        )
        attempt.update(
            {
                "method_family": "R3",
                "method": "calcall_opt",
                "calc_all": True,
                "use_default_trust": True,
                "trigger": self._rescue_trigger(outcome),
            }
        )
        if mode_analysis is not None:
            attempt["mode_analysis"] = mode_analysis
        outcome.pass2_rescue_attempts.append(attempt)
        return attempt

    @staticmethod
    def _is_geometry_nonconvergence(outcome: Pass1Outcome) -> bool:
        """Return true only for ORCA's normal geometry-iteration exhaustion."""

        history = outcome.attempt_history[-1] if outcome.attempt_history else {}
        failure_type = str(history.get("failure_type") or "").strip().lower()
        stop_reason = str(history.get("stop_reason") or "").strip().lower()
        error_text = str(outcome.opt_error or "").lower()
        return (
            failure_type == "geometry_optimization_not_converged"
            or stop_reason == "native_maxiter_checkpoint"
            or "geometry_optimization_not_converged" in error_text
            or "geometry optimization not converged" in error_text
            or "geometry optimization did not converge" in error_text
        )

    def _run_rescue_attempt(
        self,
        request: StructureRequest,
        recorder: AttemptRecorder,
        *,
        strategy: str,
        level: int,
        kind: str,
        input_xyz: Path,
        opt_spec: QCJobSpec,
        frequency_grid: str | None = None,
        frequency_scf: str | None = None,
        reuse_final_hessian: bool = False,
    ) -> Dict[str, Any]:
        opt_outcome = self._execute_attempt(
            request.id,
            recorder,
            kind=f"{kind}_opt",
            retry_level=level,
            runner=run_optimization,
            spec=opt_spec,
            input_xyz=Path(input_xyz),
        )
        opt_result = opt_outcome["result"]
        opt_xyz = self._latest_existing_path(
            Path(opt_result.output_xyz) if opt_result.output_xyz else None,
            Path(opt_outcome["last_geometry_path"]) if opt_outcome["last_geometry_path"] else None,
        )
        freq_outcome: Optional[Dict[str, Any]] = None
        freq_result: Optional[QCJobResult] = None
        freq_spec: Optional[QCJobSpec] = None
        hessian_path: Optional[Path] = None
        frequency_source = "separate_frequency"
        frequencies_cm1: List[float] = []
        frequency_status = "not_run"
        if opt_result.status == "complete" and opt_xyz is not None:
            if reuse_final_hessian:
                hessian_path = self._resolve_hessian_path(opt_outcome)
                if hessian_path is not None:
                    try:
                        frequencies_cm1, _mode_matrix = self._parse_orca_hessian_vibrations(
                            hessian_path
                        )
                    except (OSError, ValueError):
                        frequencies_cm1 = []
                if frequencies_cm1:
                    frequency_status = "complete"
                    frequency_source = "calcall_final_hessian"
                else:
                    frequency_status = "failed"
                    frequency_source = "calcall_hessian_unavailable"
            else:
                freq_spec = self._build_frequency_spec(
                    request,
                    grid=frequency_grid,
                    scf=frequency_scf,
                )
                freq_outcome = self._execute_attempt(
                    request.id,
                    recorder,
                    kind=f"{kind}_freq",
                    retry_level=level,
                    runner=run_frequency,
                    spec=freq_spec,
                    input_xyz=Path(opt_xyz),
                )
                freq_result = freq_outcome["result"]
                frequency_status = freq_result.status
                frequencies_cm1 = [float(value) for value in (freq_result.frequencies_cm1 or ())]
                hessian_path = self._resolve_hessian_path(freq_outcome)

        frequency_output = self._latest_existing_path(
            Path(freq_result.output_file) if freq_result and freq_result.output_file else None,
            Path(freq_outcome["record"].output_out_path) if freq_outcome else None,
        )
        imaginary_frequencies_cm1 = [
            float(value) for value in frequencies_cm1 if float(value) < 0.0
        ]
        (
            ts_classification,
            int_classification,
            minimum_classification,
            alignment_score,
        ) = self._classify_structure(
            request,
            opt_xyz=opt_xyz,
            frequency_output=frequency_output,
            frequencies_cm1=frequencies_cm1,
            hessian_path=hessian_path if frequency_source == "calcall_final_hessian" else None,
        )

        error_messages = [
            message
            for message in (
                opt_result.error,
                freq_result.error if freq_result is not None else None,
            )
            if message
        ]
        return {
            "attempt_id": opt_outcome["record"].attempt_id,
            "strategy": strategy,
            "level": level,
            "status": (
                "complete"
                if opt_result.status == "complete" and frequency_status == "complete"
                else (frequency_status if opt_result.status == "complete" else opt_result.status)
            ),
            "error": "; ".join(error_messages) or None,
            "opt_status": opt_result.status,
            "opt_xyz": str(opt_xyz) if opt_xyz else None,
            "stop_geometry_path": (
                str(opt_outcome["stop_geometry_path"])
                if opt_outcome.get("stop_geometry_path")
                else None
            ),
            "optimization_monitor": opt_outcome.get("optimization_monitor"),
            "opt_output": str(opt_outcome["record"].output_out_path),
            "opt_energy_hartree": opt_result.energy_hartree,
            "frequency_status": frequency_status,
            "frequency_source": frequency_source,
            "frequency_output": str(frequency_output) if frequency_output else None,
            "frequency_energy_hartree": freq_result.energy_hartree if freq_result is not None else None,
            "frequencies_cm1": frequencies_cm1,
            "imaginary_frequencies_cm1": imaginary_frequencies_cm1,
            "zero_point_energy_hartree": (
                freq_result.zero_point_energy_hartree if freq_result is not None else None
            ),
            "thermal_energy_hartree": (
                freq_result.thermal_energy_hartree if freq_result is not None else None
            ),
            "enthalpy_hartree": freq_result.enthalpy_hartree if freq_result is not None else None,
            "gibbs_free_energy_hartree": (
                freq_result.gibbs_free_energy_hartree if freq_result is not None else None
            ),
            "gibbs_correction_hartree": (
                freq_result.gibbs_correction_hartree if freq_result is not None else None
            ),
            "hessian_filename": str(hessian_path) if hessian_path is not None else None,
            "ts_classification": ts_classification,
            "int_classification": int_classification,
            "minimum_classification": minimum_classification,
            "alignment_score": alignment_score,
            "gradient_norm": (
                float(opt_result.extra["gradient_norm"])
                if opt_result.extra
                and isinstance(opt_result.extra, Mapping)
                and opt_result.extra.get("gradient_norm") is not None
                else None
            ),
            "theory_signature": {
                "method": opt_spec.method,
                "basis": opt_spec.basis,
                "aux_basis": opt_spec.aux_basis,
                "route": opt_spec.route,
                "route_extras": opt_spec.route_extras,
                "initial_hessian": opt_spec.initial_hessian,
                "hessian_filename": opt_spec.hessian_filename,
                "ts_mode": opt_spec.ts_mode,
                "recalc_hessian": opt_spec.recalc_hessian,
                "trust": opt_spec.trust,
                "grid": opt_spec.grid,
                "scf": opt_spec.scf,
            },
        }

    def _run_irc_attempt(
        self,
        request: StructureRequest,
        recorder: AttemptRecorder,
        *,
        irc_spec: IRCJobSpec,
        input_xyz: Path,
        output_dir: Path,
        level: int = 1,
    ) -> Dict[str, Any]:
        irc_dir = Path(output_dir) / "diagnostics" / recorder.next_attempt(
            "irc", retry_level=level
        ).attempt_id
        irc_dir.mkdir(parents=True, exist_ok=True)
        result = run_irc(
            irc_spec,
            input_xyz,
            irc_dir,
            self.config,
            charge=int(request.charge),
            spin=int(request.multiplicity),
        )
        return {
            "status": result.status,
            "error": result.error,
            "irc_dir": str(irc_dir),
            "endpoint_a_xyz": (
                str(Path(result.extra["endpoint_a_xyz"]))
                if result.extra and result.extra.get("endpoint_a_xyz")
                else None
            ),
            "endpoint_b_xyz": (
                str(Path(result.extra["endpoint_b_xyz"]))
                if result.extra and result.extra.get("endpoint_b_xyz")
                else None
            ),
        }

    def _rescue_int_via_irc(
        self,
        request: StructureRequest,
        outcome: Pass1Outcome,
        recorder: AttemptRecorder,
        structure_dir: Path,
        params: MethodParams,
    ) -> Optional[Dict[str, Any]]:
        """Recover the dipolar intermediate by IRC-stretched-end re-optimization.

        A collapsed INT means the S2 midpoint seed slid to the product; the
        true intermediate lies on the stretched (forming-bond-open) side of the
        TS.  We run a short IRC (``irc_max_iter`` steps per direction, reusing
        the TS Hessian), take the terminal frame of the stretched trajectory as
        a fresh seed, and re-optimize it without constraints.
        """
        ts_outcome = self._variant_ts_outcome.get(request.variant_id or request.id)
        if ts_outcome is None:
            self._append_rescue_skipped(
                outcome, "R4", "irc_midpoint_recovery", "no_complete_sibling_ts"
            )
            return None
        ts_xyz = ts_outcome.canonical_xyz
        if ts_xyz is None or not Path(ts_xyz).exists():
            self._append_rescue_skipped(
                outcome, "R4", "irc_midpoint_recovery", "ts_geometry_missing"
            )
            return None
        ts_hessian = self._resolve_ts_hessian(ts_outcome)
        if ts_hessian is None:
            self._append_rescue_skipped(
                outcome, "R4", "irc_midpoint_recovery", "ts_hessian_missing"
            )
            return None

        irc_spec = IRCJobSpec(
            direction=str(params.irc_direction or "both"),
            max_iter=int(params.irc_max_iter or 5),
            init_hessian="read",
            hessian_filename=str(ts_hessian),
            method=self.profile.geometry_method,
            basis=self.profile.geometry_basis,
        )
        irc_attempt = self._run_irc_attempt(
            request,
            recorder,
            irc_spec=irc_spec,
            input_xyz=Path(ts_xyz),
            output_dir=Path(structure_dir),
        )
        if irc_attempt["status"] != "complete":
            record = {
                "strategy": "irc_midpoint_recovery",
                "method_family": "R4",
                "method": "irc_midpoint_recovery",
                "status": "failed",
                "error": irc_attempt.get("error"),
                "irc_dir": irc_attempt.get("irc_dir"),
                "trigger": self._rescue_trigger(outcome),
            }
            outcome.pass2_rescue_attempts.append(record)
            return record

        trajectories = parse_irc_trajectories(Path(irc_attempt["irc_dir"]))
        endpoints = classify_irc_endpoints(trajectories, request.forming_bonds)
        if endpoints is None:
            self._append_rescue_skipped(
                outcome, "R4", "irc_midpoint_recovery", "no_usable_irc_endpoints"
            )
            return None

        stretch_frame = endpoints["stretch_frame"]
        _coords, symbols = read_xyz(Path(ts_xyz))
        seed_xyz = frame_to_xyz(
            endpoints["stretch_trajectory"],
            stretch_frame.frame_index,
            Path(structure_dir) / "irc_stretch_seed.xyz",
            symbols,
        )
        if seed_xyz is None:
            self._append_rescue_skipped(
                outcome, "R4", "irc_midpoint_recovery", "stretch_frame_unusable"
            )
            return None

        remedy_outcome, _remedy_spec, _initial_hessian = self._run_primary_opt(
            request,
            Path(seed_xyz),
            recorder,
        )
        if remedy_outcome["result"].status != "complete":
            self._append_rescue_skipped(
                outcome, "R4", "irc_midpoint_recovery", "stretch_opt_failed"
            )
            return None

        opt_xyz = self._latest_existing_path(
            Path(remedy_outcome["result"].output_xyz)
            if remedy_outcome["result"].output_xyz
            else None,
            Path(remedy_outcome["last_geometry_path"])
            if remedy_outcome.get("last_geometry_path")
            else None,
        )
        if opt_xyz is None or not Path(opt_xyz).exists():
            self._append_rescue_skipped(
                outcome, "R4", "irc_midpoint_recovery", "stretch_opt_geometry_missing"
            )
            return None

        freq_outcome = self._run_frequency(request, Path(opt_xyz), recorder)
        freq_result = freq_outcome[0] if freq_outcome else None
        frequencies_cm1 = (
            [float(value) for value in (freq_result["result"].frequencies_cm1 or ())]
            if freq_result and freq_result.get("result")
            else []
        )
        int_classification = self._classify_int_v2(
            request,
            opt_xyz=Path(opt_xyz),
            frequencies_cm1=frequencies_cm1,
            energy_int_hartree=(
                remedy_outcome["result"].energy_hartree
                if remedy_outcome["result"].energy_hartree is not None
                else self._opt_energy_for(outcome)
            ),
            converged=(
                True
                if remedy_outcome["result"].status == "complete"
                else None
            ),
        )
        identity = str(int_classification.get("identity") or "")
        significant_imag = [
            value for value in frequencies_cm1 if value < _MINIMUM_IMAGINARY_CUTOFF_CM1
        ]
        if identity in {"distinct_intermediate", "dipolar_intermediate"} and not significant_imag:
            outcome.opt_xyz = Path(opt_xyz)
            outcome.opt_status = "complete"
            outcome.irc_status = "complete"
            outcome.int_classification = int_classification
            outcome.irc_endpoints = {
                "mechanism": "intermediate_recovered",
                "seed_source": "irc_stretched_endpoint",
                "stretch_mean_fb_distance": endpoints["stretch_mean_fb_distance"],
                "product_mean_fb_distance": endpoints["product_mean_fb_distance"],
                "endpoint_a_xyz": irc_attempt.get("endpoint_a_xyz"),
                "endpoint_b_xyz": irc_attempt.get("endpoint_b_xyz"),
            }
            record = {
                "strategy": "irc_midpoint_recovery",
                "method_family": "R4",
                "method": "irc_midpoint_recovery",
                "status": "complete",
                "mechanism": "intermediate_recovered",
                "opt_xyz": str(opt_xyz),
                "stretch_mean_fb_distance": endpoints["stretch_mean_fb_distance"],
                "irc_dir": irc_attempt.get("irc_dir"),
                "trigger": self._rescue_trigger(outcome),
            }
            outcome.pass2_rescue_attempts.append(record)
            return record

        outcome.int_classification = int_classification
        record = {
            "strategy": "irc_midpoint_recovery",
            "method_family": "R4",
            "method": "irc_midpoint_recovery",
            "status": "complete",
            "mechanism": "stretch_endpoint_not_stationary",
            "opt_xyz": str(opt_xyz),
            "identity": identity,
            "imaginary_frequencies_cm1": significant_imag,
            "stretch_mean_fb_distance": endpoints["stretch_mean_fb_distance"],
            "irc_dir": irc_attempt.get("irc_dir"),
            "trigger": self._rescue_trigger(outcome),
        }
        outcome.irc_status = "complete"
        outcome.irc_endpoints = {
            "mechanism": "stretch_endpoint_not_stationary",
            "identity": identity,
            "stretch_mean_fb_distance": endpoints["stretch_mean_fb_distance"],
            "endpoint_a_xyz": irc_attempt.get("endpoint_a_xyz"),
            "endpoint_b_xyz": irc_attempt.get("endpoint_b_xyz"),
        }
        outcome.pass2_rescue_attempts.append(record)
        return record

    def _resolve_ts_hessian(self, ts_outcome: Pass1Outcome) -> Optional[Path]:
        for candidate in (
            ts_outcome.canonical_hessian_path,
            ts_outcome.frequency_hessian_path,
            ts_outcome.canonical_frequency_output,
        ):
            if candidate is not None and Path(candidate).exists():
                path = Path(candidate)
                if path.suffix == ".hess":
                    return path
                hess = path.with_suffix(".hess")
                if hess.exists():
                    return hess
        return None

    @staticmethod
    def _freq_output_fallback(pass3: Pass1Outcome) -> Optional[str]:
        if pass3.canonical_hessian_path is not None and Path(pass3.canonical_hessian_path).exists():
            hess = Path(pass3.canonical_hessian_path)
            out = hess.with_suffix(".out")
            if out.exists():
                return str(out)
            return str(hess)
        return None

    def _run_pass3_canonical(
        self,
        pass2_outcomes: Sequence[Pass1Outcome],
        output_dir: Path,
    ) -> List[Pass1Outcome]:
        stage = self.profile.stage
        schedule = resolve_stage_schedule(
            self.config,
            "step3" if stage == "S3" else "step4",
        )
        payloads = [
            {"outcome": copy.deepcopy(outcome), "output_dir": Path(output_dir)}
            for outcome in pass2_outcomes
        ]
        return run_structure_queue(
            payloads,
            lambda payload: self._run_pass3_canonical_one(payload["outcome"], payload["output_dir"]),
            schedule,
            thread_name_prefix=f"rph-{stage.lower()}-pass3",
        )

    def _run_pass3_canonical_one(
        self,
        outcome: Pass1Outcome,
        output_dir: Path,
    ) -> Pass1Outcome:
        request = self._requests_by_id[outcome.structure_id]
        recorder = self._recorders[outcome.structure_id]
        structure_dir = Path(output_dir) / outcome.structure_id

        if not self._pass3_enabled_for_request(request):
            return outcome

        canonical_attempt_id = self._select_canonical_attempt(outcome)
        outcome.canonical_attempt_id = canonical_attempt_id
        outcome.rescue_winning_method = self._rescue_winning_method(outcome, canonical_attempt_id)
        resolved_candidate = self._lookup_candidate(outcome, canonical_attempt_id)
        if resolved_candidate is None:
            outcome.resolved_kind = None
            outcome.resolved_identity = None
            outcome.identity_status = "not_checked"
            outcome.ml_usability = compute_ml_usability(
                requested_role=request.role,
                requested_kind=request.kind,
                opt_converged=False,
                sp_complete=False,
                freq_available=False,
                identity_status=outcome.identity_status,
                ts_frequency_valid=None,
                geometry_consistent=None,
            )
            outcome.status = "failed" if outcome.opt_status != "complete" else outcome.status
            return outcome

        source_xyz = resolved_candidate.get("opt_xyz")
        if source_xyz is None:
            outcome.ml_usability = compute_ml_usability(
                requested_role=request.role,
                requested_kind=request.kind,
                opt_converged=False,
                sp_complete=False,
                freq_available=False,
                identity_status="not_checked",
                ts_frequency_valid=None,
                geometry_consistent=None,
            )
            outcome.status = "failed"
            return outcome

        canonical_xyz = structure_dir / "canonical.xyz"
        copied_canonical = self._copy_file(Path(source_xyz), canonical_xyz)
        if copied_canonical is None:
            outcome.status = "failed"
            outcome.error = outcome.error or f"Unable to materialize canonical geometry for {outcome.structure_id}"
            return outcome
        outcome.canonical_xyz = canonical_xyz
        outcome.geometry_hash = self._file_sha256(canonical_xyz)

        primary_geometry_hash = self._file_sha256(outcome.opt_xyz) if outcome.opt_xyz else None
        canonical_geometry_hash = self._file_sha256(canonical_xyz)
        reuse_primary_frequency = bool(
            canonical_attempt_id == (outcome.primary_attempt_id or outcome.current_attempt_id)
            and outcome.frequency_status == "complete"
            and primary_geometry_hash is not None
            and primary_geometry_hash == canonical_geometry_hash
        )
        # R1 CalcAll already wrote the final Hessian at this exact geometry.
        # Re-submitting canonical FREQ here would defeat the rescue's purpose
        # and duplicate its expensive Hessian evaluation.
        reuse_calcall_hessian = bool(
            resolved_candidate.get("frequency_status") == "complete"
            and resolved_candidate.get("frequency_source") == "calcall_final_hessian"
            and resolved_candidate.get("hessian_filename")
        )

        candidate_bundle = copy.deepcopy(resolved_candidate)
        freq_geometry_hash: Optional[str] = None
        if reuse_primary_frequency or reuse_calcall_hessian:
            frequency_candidate: Mapping[str, Any]
            if reuse_primary_frequency:
                frequency_candidate = {
                    "frequency_status": outcome.frequency_status,
                    "frequency_output": outcome.frequency_output,
                    "frequency_energy_hartree": outcome.frequency_energy_hartree,
                    "frequencies_cm1": outcome.frequencies_cm1,
                    "imaginary_frequencies_cm1": outcome.imaginary_frequencies_cm1,
                    "zero_point_energy_hartree": outcome.frequency_zero_point_energy_hartree,
                    "thermal_energy_hartree": outcome.frequency_thermal_energy_hartree,
                    "enthalpy_hartree": outcome.frequency_enthalpy_hartree,
                    "gibbs_free_energy_hartree": outcome.frequency_gibbs_free_energy_hartree,
                    "gibbs_correction_hartree": outcome.frequency_gibbs_correction_hartree,
                    "hessian_filename": outcome.frequency_hessian_path,
                }
            else:
                frequency_candidate = resolved_candidate
            outcome.canonical_frequency_status = str(frequency_candidate.get("frequency_status") or "not_run")
            outcome.canonical_frequency_output = self._latest_existing_path(
                Path(frequency_candidate["frequency_output"])
                if frequency_candidate.get("frequency_output")
                else None,
            )
            outcome.canonical_frequency_energy_hartree = frequency_candidate.get("frequency_energy_hartree")
            outcome.canonical_frequencies_cm1 = [
                float(value) for value in (frequency_candidate.get("frequencies_cm1") or ())
            ]
            outcome.canonical_imaginary_frequencies_cm1 = [
                float(value)
                for value in (frequency_candidate.get("imaginary_frequencies_cm1") or ())
            ]
            outcome.canonical_zero_point_energy_hartree = frequency_candidate.get("zero_point_energy_hartree")
            outcome.canonical_thermal_energy_hartree = frequency_candidate.get("thermal_energy_hartree")
            outcome.canonical_enthalpy_hartree = frequency_candidate.get("enthalpy_hartree")
            outcome.canonical_gibbs_free_energy_hartree = frequency_candidate.get("gibbs_free_energy_hartree")
            outcome.canonical_gibbs_correction_hartree = frequency_candidate.get("gibbs_correction_hartree")
            raw_hessian = frequency_candidate.get("hessian_filename")
            outcome.canonical_hessian_path = Path(raw_hessian) if raw_hessian else None
            freq_geometry_hash = (
                primary_geometry_hash if reuse_primary_frequency else canonical_geometry_hash
            )
        else:
            freq_spec = self._build_frequency_spec(request)
            freq_outcome = self._execute_attempt(
                request.id,
                recorder,
                kind="canonical_freq",
                retry_level=0,
                runner=run_frequency,
                spec=freq_spec,
                input_xyz=canonical_xyz,
            )
            self._populate_legacy_attempt_dir(structure_dir / "canonical_freq", freq_outcome)
            freq_result = freq_outcome["result"]
            outcome.canonical_frequency_status = freq_result.status
            outcome.canonical_frequency_output = self._latest_existing_path(
                Path(freq_result.output_file) if freq_result.output_file else None,
                freq_outcome["record"].output_out_path,
            )
            outcome.canonical_frequency_energy_hartree = freq_result.energy_hartree
            outcome.canonical_frequencies_cm1 = [
                float(value) for value in (freq_result.frequencies_cm1 or ())
            ]
            outcome.canonical_imaginary_frequencies_cm1 = [
                float(value)
                for value in outcome.canonical_frequencies_cm1
                if float(value) < 0.0
            ]
            outcome.canonical_zero_point_energy_hartree = freq_result.zero_point_energy_hartree
            outcome.canonical_thermal_energy_hartree = freq_result.thermal_energy_hartree
            outcome.canonical_enthalpy_hartree = freq_result.enthalpy_hartree
            outcome.canonical_gibbs_free_energy_hartree = freq_result.gibbs_free_energy_hartree
            outcome.canonical_gibbs_correction_hartree = freq_result.gibbs_correction_hartree
            outcome.canonical_hessian_path = self._resolve_hessian_path(freq_outcome)
            freq_geometry_hash = canonical_geometry_hash
            (
                candidate_bundle["ts_classification"],
                candidate_bundle["int_classification"],
                candidate_bundle["minimum_classification"],
                candidate_bundle["alignment_score"],
            ) = self._classify_structure(
                request,
                opt_xyz=canonical_xyz,
                frequency_output=outcome.canonical_frequency_output,
                frequencies_cm1=outcome.canonical_frequencies_cm1,
            )
            if request.role == "intermediate":
                candidate_bundle["int_classification"] = self._classify_int_v2(
                    request,
                    opt_xyz=canonical_xyz,
                    frequencies_cm1=outcome.canonical_frequencies_cm1,
                    energy_int_hartree=self._opt_energy_for(outcome),
                    converged=True if outcome.opt_status == "complete" else None,
                )
            candidate_bundle["frequency_status"] = outcome.canonical_frequency_status
            candidate_bundle["frequency_output"] = (
                str(outcome.canonical_frequency_output) if outcome.canonical_frequency_output else None
            )
            candidate_bundle["frequency_energy_hartree"] = outcome.canonical_frequency_energy_hartree
            candidate_bundle["frequencies_cm1"] = list(outcome.canonical_frequencies_cm1)
            candidate_bundle["imaginary_frequencies_cm1"] = list(
                outcome.canonical_imaginary_frequencies_cm1
            )
            candidate_bundle["enthalpy_hartree"] = outcome.canonical_enthalpy_hartree
            candidate_bundle["gibbs_free_energy_hartree"] = (
                outcome.canonical_gibbs_free_energy_hartree
            )
            candidate_bundle["gibbs_correction_hartree"] = outcome.canonical_gibbs_correction_hartree
            candidate_bundle["hessian_filename"] = (
                str(outcome.canonical_hessian_path) if outcome.canonical_hessian_path else None
            )

        outcome.ts_classification = candidate_bundle.get("ts_classification")
        outcome.int_classification = candidate_bundle.get("int_classification")
        outcome.minimum_classification = candidate_bundle.get("minimum_classification")
        outcome.resolved_kind, outcome.resolved_identity = self._resolved_kind_identity(
            request,
            candidate_bundle,
        )
        outcome.identity_status = self._infer_identity_status(request, candidate_bundle)

        sp_spec = self._build_single_point_spec(request)
        sp_outcome = self._execute_attempt(
            request.id,
            recorder,
            kind="sp",
            retry_level=0,
            runner=run_single_point,
            spec=sp_spec,
            input_xyz=canonical_xyz,
        )
        self._populate_legacy_attempt_dir(structure_dir / "sp", sp_outcome)
        sp_result = sp_outcome["result"]
        outcome.sp_status = "complete" if sp_result.status == "complete" else "failed"
        outcome.sp_energy_hartree = sp_result.energy_hartree
        outcome.sp_output = self._latest_existing_path(
            Path(sp_result.output_file) if sp_result.output_file else None,
            sp_outcome["record"].output_out_path,
        )
        sp_geometry_hash = canonical_geometry_hash
        outcome.thermochemistry = self._compute_composite_thermochemistry(
            request,
            outcome,
            freq_geometry_hash=freq_geometry_hash,
            sp_geometry_hash=sp_geometry_hash,
        )

        geometry_consistent = (
            outcome.thermochemistry.get("geometry_consistent")
            if outcome.thermochemistry is not None
            else None
        )
        outcome.ml_usability = compute_ml_usability(
            requested_role=request.role,
            requested_kind=request.kind,
            opt_converged=bool(canonical_xyz.exists()),
            sp_complete=outcome.sp_status == "complete",
            freq_available=outcome.canonical_frequency_status == "complete",
            identity_status=outcome.identity_status,
            ts_frequency_valid=self._final_ts_frequency_valid(candidate_bundle),
            geometry_consistent=geometry_consistent,
        )

        error_messages = [
            message
            for message in (
                outcome.error,
                sp_result.error,
            )
            if message
        ]
        outcome.error = "; ".join(dict.fromkeys(error_messages)) or None
        if outcome.sp_status == "complete" and outcome.canonical_frequency_status == "complete":
            outcome.status = "complete"
        elif canonical_xyz.exists():
            outcome.status = "degraded"
        else:
            outcome.status = "failed"
        return outcome

    def _select_canonical_attempt(self, outcome: Pass1Outcome) -> Optional[str]:
        candidates = []
        primary_candidate = self._candidate_from_primary(outcome)
        if primary_candidate.get("opt_status") == "complete":
            candidates.append(primary_candidate)
        candidates.extend(
            attempt
            for attempt in outcome.pass2_rescue_attempts
            if attempt.get("attempt_id") is not None and attempt.get("opt_status") == "complete"
        )
        if not candidates:
            return None
        candidates.sort(key=self._candidate_sort_key)
        return str(candidates[0].get("attempt_id") or "") or None

    @staticmethod
    def _rescue_winning_method(
        outcome: Pass1Outcome,
        canonical_attempt_id: Optional[str],
    ) -> Optional[str]:
        if not canonical_attempt_id:
            return None
        if canonical_attempt_id == (outcome.primary_attempt_id or outcome.current_attempt_id):
            return None
        for attempt in outcome.pass2_rescue_attempts:
            if attempt.get("attempt_id") == canonical_attempt_id:
                return str(attempt.get("strategy") or attempt.get("method") or "") or None
        return None

    @staticmethod
    def _pass3_enabled_for_request(request: StructureRequest) -> bool:
        if request.role in {"precursor", "product"}:
            return bool(
                request.s1_ensemble_thermodynamics is not None
                or request.ensemble_thermochemistry_correction_hartree is not None
                or request.s1_thermochemistry_status is not None
            )
        return True

    def _candidate_sort_key(self, candidate: Mapping[str, Any]) -> Tuple[Any, ...]:
        ts_classification = dict(candidate.get("ts_classification") or {})
        int_classification = dict(candidate.get("int_classification") or {})
        minimum_classification = dict(candidate.get("minimum_classification") or {})
        stationary_point_class = str(ts_classification.get("stationary_point_class") or "")
        mode_identity = str(ts_classification.get("mode_identity") or "unavailable")
        hessian_index = ts_classification.get("hessian_index")
        if hessian_index is None:
            hessian_rank = 99
        else:
            hessian_rank = 0 if int(hessian_index) == 1 else abs(int(hessian_index) - 1) + 1

        stationary_rank_map = {
            "valid_target_ts": 0,
            "soft_target_ts": 1,
            "minimum_after_optts": 2,
            "first_order_wrong_mode": 3,
            "higher_order_saddle": 4,
        }
        if stationary_point_class:
            stationary_rank = stationary_rank_map.get(stationary_point_class, 5)
        elif int_classification:
            stationary_rank = {
                "distinct_intermediate": 0,
                "collapsed_to_precursor": 1,
                "collapsed_to_product": 1,
                "topology_ambiguous": 2,
                "imaginary_frequency": 3,
                "opt_failed": 4,
            }.get(str(int_classification.get("identity") or ""), 5)
        else:
            stationary_rank = {
                "valid_minimum": 0,
                "not_checked": 0,
                "identity_drift": 2,
                "imaginary_frequency": 3,
            }.get(str(minimum_classification.get("identity") or "not_checked"), 4)

        mode_rank = {"target": 0, "unavailable": 1, "unrelated": 2}.get(mode_identity, 3)
        converged_rank = 0 if candidate.get("frequency_status") == "complete" else 1
        alignment_score = candidate.get("alignment_score")
        if alignment_score is None:
            alignment_score = ts_classification.get("alignment_score")
        if alignment_score is None:
            alignment_score = {"target": 1.0, "unavailable": 0.0, "unrelated": -1.0}.get(mode_identity, 0.0)
        gradient_norm = candidate.get("gradient_norm")
        if gradient_norm is None:
            gradient_norm = float("inf")
        return (
            stationary_rank,
            mode_rank,
            hessian_rank,
            converged_rank,
            -float(alignment_score),
            float(gradient_norm),
        )

    def _lookup_candidate(
        self,
        outcome: Pass1Outcome,
        attempt_id: Optional[str],
    ) -> Optional[Dict[str, Any]]:
        if attempt_id is None:
            return None
        primary = self._candidate_from_primary(outcome)
        if primary.get("attempt_id") == attempt_id:
            return primary
        for candidate in outcome.pass2_rescue_attempts:
            if candidate.get("attempt_id") == attempt_id:
                return copy.deepcopy(candidate)
        return None

    def _candidate_from_primary(self, outcome: Pass1Outcome) -> Dict[str, Any]:
        """Build a candidate record from the primary (pass 1) outcome."""
        return {
            "attempt_id": outcome.primary_attempt_id,
            "strategy": "primary",
            "level": 0,
            "status": "complete" if outcome.opt_status == "complete" else outcome.opt_status,
            "opt_status": outcome.opt_status,
            "opt_xyz": str(outcome.opt_xyz) if outcome.opt_xyz else None,
            "opt_output": str(outcome.opt_output) if outcome.opt_output else None,
            "opt_energy_hartree": outcome.opt_energy_hartree,
            "frequency_status": outcome.frequency_status,
            "frequency_output": (
                str(outcome.frequency_output) if outcome.frequency_output else None
            ),
            "frequency_energy_hartree": outcome.frequency_energy_hartree,
            "frequencies_cm1": list(outcome.frequencies_cm1),
            "imaginary_frequencies_cm1": list(outcome.imaginary_frequencies_cm1),
            "frequency_zero_point_energy_hartree": outcome.frequency_zero_point_energy_hartree,
            "frequency_thermal_energy_hartree": outcome.frequency_thermal_energy_hartree,
            "frequency_enthalpy_hartree": outcome.frequency_enthalpy_hartree,
            "frequency_gibbs_free_energy_hartree": outcome.frequency_gibbs_free_energy_hartree,
            "frequency_gibbs_correction_hartree": outcome.frequency_gibbs_correction_hartree,
            "frequency_hessian_path": (
                str(outcome.frequency_hessian_path)
                if outcome.frequency_hessian_path
                else None
            ),
            "ts_classification": copy.deepcopy(outcome.ts_classification),
            "int_classification": copy.deepcopy(outcome.int_classification),
            "minimum_classification": copy.deepcopy(outcome.minimum_classification),
        }

    def _candidate_is_valid_for_request(
        self,
        request: StructureRequest,
        candidate: Mapping[str, Any],
    ) -> bool:
        # A partial geometry must never terminate a rescue ladder.  In
        # particular, a failed optimization has no frequency result; treating
        # its default imaginary-count payload as zero previously skipped INT
        # R2 entirely.
        if str(candidate.get("opt_status") or "") != "complete":
            return False
        if str(candidate.get("frequency_status") or "") != "complete":
            return False
        if request.kind == "ts":
            return str((candidate.get("ts_classification") or {}).get("stationary_point_class") or "") == "valid_target_ts"
        if request.role == "intermediate":
            return (
                self._candidate_has_no_significant_imag(candidate)
                and str((candidate.get("int_classification") or {}).get("identity") or "")
                == "distinct_intermediate"
            )
        return self._candidate_has_no_significant_imag(candidate)

    @staticmethod
    def _candidate_has_no_significant_imag(candidate: Mapping[str, Any]) -> bool:
        minimum_classification = dict(candidate.get("minimum_classification") or {})
        imag_count = int(minimum_classification.get("imaginary_count") or 0)
        if imag_count != 0:
            return False
        int_classification = dict(candidate.get("int_classification") or {})
        return str(int_classification.get("identity") or "") != "imaginary_frequency"

    def _build_mode_displaced_xyz(
        self,
        source_xyz: Path,
        frequency_output: Path,
        destination: Path,
        *,
        sign: float,
        mode_index: int = 0,
    ) -> Path:
        if not source_xyz.exists() or not frequency_output.exists():
            raise ValueError("Mode-displacement rescue requires both geometry and frequency output")
        content = frequency_output.read_text(encoding="utf-8", errors="ignore")
        displacement_vectors = _parse_orca_displacement_vectors(content)
        if not displacement_vectors:
            raise ValueError(f"No imaginary-mode displacement vectors found in {frequency_output}")
        coordinates, symbols = read_xyz(source_xyz)
        ordered = sorted(displacement_vectors)
        if mode_index >= len(ordered):
            raise ValueError(
                f"Imaginary-mode index {mode_index} unavailable (only {len(ordered)} modes parsed)"
            )
        mode_vector = np.asarray(displacement_vectors[ordered[mode_index]], dtype=float)
        if coordinates.shape != mode_vector.shape:
            raise ValueError("Imaginary-mode displacement vector shape mismatch")
        max_norm = float(np.max(np.linalg.norm(mode_vector, axis=1)))
        if max_norm <= 1e-12:
            raise ValueError("Imaginary-mode displacement vector is zero-length")
        displaced = np.asarray(coordinates, dtype=float) + (
            float(sign) * (mode_vector / max_norm) * _MODE_DISPLACEMENT_SCALE_ANGSTROM
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        write_xyz(destination, displaced, symbols, title=f"mode_displaced_{sign:+.1f}")
        return destination

    def _stricter_grid(self) -> str:
        if str(self.profile.geometry_grid or "").strip():
            return str(self.profile.geometry_grid)
        return "DefGrid3"

    def _stricter_scf(self) -> str:
        if str(self.profile.geometry_scf or "").strip():
            return str(self.profile.geometry_scf)
        return "TightSCF"

    def _minimum_rescue_enabled(self) -> bool:
        refinement_cfg = dict(self.config.get("refinement", {}) or {})
        common_cfg = dict(refinement_cfg.get("common", {}) or {})
        rescue_cfg = dict(common_cfg.get("rescue", {}) or {})
        minimum_cfg = dict(rescue_cfg.get("minimum_rescue", {}) or {})
        return bool(minimum_cfg.get("enabled", True))

    def _compute_composite_thermochemistry(
        self,
        request: StructureRequest,
        outcome: Pass1Outcome,
        *,
        freq_geometry_hash: Optional[str],
        sp_geometry_hash: Optional[str],
    ) -> Optional[Dict[str, Any]]:
        if outcome.sp_status != "complete" or outcome.sp_energy_hartree is None:
            return None
        if outcome.canonical_frequency_status != "complete":
            return None
        if outcome.canonical_frequency_energy_hartree is None:
            return None
        if outcome.canonical_enthalpy_hartree is None or outcome.canonical_gibbs_free_energy_hartree is None:
            return None

        ensemble_gibbs, ensemble_enthalpy, ensemble_correction_id = self._extract_ensemble_corrections(request)
        gibbs_correction = (
            float(outcome.canonical_gibbs_free_energy_hartree)
            - float(outcome.canonical_frequency_energy_hartree)
        )
        enthalpy_correction = (
            float(outcome.canonical_enthalpy_hartree)
            - float(outcome.canonical_frequency_energy_hartree)
        )
        geometry_consistent = (
            freq_geometry_hash == sp_geometry_hash
            if freq_geometry_hash is not None and sp_geometry_hash is not None
            else None
        )
        return {
            "gibbs_free_energy_hartree": (
                float(outcome.sp_energy_hartree) + gibbs_correction + float(ensemble_gibbs)
            ),
            "enthalpy_hartree": (
                float(outcome.sp_energy_hartree) + enthalpy_correction + float(ensemble_enthalpy)
            ),
            "single_point_energy_hartree": float(outcome.sp_energy_hartree),
            "frequency_energy_hartree": float(outcome.canonical_frequency_energy_hartree),
            "frequency_gibbs_free_energy_hartree": float(outcome.canonical_gibbs_free_energy_hartree),
            "frequency_enthalpy_hartree": float(outcome.canonical_enthalpy_hartree),
            "ensemble_gibbs_correction_hartree": float(ensemble_gibbs),
            "ensemble_enthalpy_correction_hartree": float(ensemble_enthalpy),
            "geometry_level": self._format_theory_level(
                self.profile.geometry_method,
                self.profile.geometry_basis,
                self.profile.geometry_aux_basis,
            ),
            "frequency_level": self._format_theory_level(
                self.profile.frequency_method,
                self.profile.frequency_basis,
                self.profile.geometry_aux_basis,
            ),
            "single_point_level": self._format_theory_level(
                self.profile.sp_method,
                self.profile.sp_basis,
                self.profile.sp_aux_basis,
            ),
            "temperature_K": float(self.profile.temperature_k),
            "standard_state": self.profile.standard_state,
            "qrrho": bool(self.profile.qrrho),
            "ensemble_correction_id": ensemble_correction_id,
            "geometry_consistent": geometry_consistent,
        }

    def _extract_ensemble_corrections(
        self,
        request: StructureRequest,
    ) -> Tuple[float, float, Optional[str]]:
        payload = dict(request.s1_ensemble_thermodynamics or {})
        gibbs_candidates = [
            payload.get("ensemble_gibbs_correction_hartree"),
            payload.get("ensemble_thermochemistry_correction_hartree"),
            request.ensemble_thermochemistry_correction_hartree,
        ]
        enthalpy_candidates = [
            payload.get("ensemble_enthalpy_correction_hartree"),
            payload.get("ensemble_thermochemistry_enthalpy_correction_hartree"),
            payload.get("enthalpy_correction_hartree"),
            request.ensemble_thermochemistry_correction_hartree,
        ]
        ensemble_gibbs = next(
            (float(value) for value in gibbs_candidates if value is not None),
            0.0,
        )
        ensemble_enthalpy = next(
            (float(value) for value in enthalpy_candidates if value is not None),
            ensemble_gibbs,
        )
        ensemble_correction_id = (
            str(payload.get("ensemble_correction_id") or payload.get("reference_conformer_id") or self.profile.ensemble_correction_source)
            if payload or self.profile.ensemble_correction_source
            else None
        )
        return ensemble_gibbs, ensemble_enthalpy, ensemble_correction_id

    def _resolved_kind_identity(
        self,
        request: StructureRequest,
        candidate: Mapping[str, Any],
    ) -> Tuple[Optional[str], Optional[str]]:
        if request.kind == "ts":
            ts_classification = dict(candidate.get("ts_classification") or {})
            stationary_point_class = str(ts_classification.get("stationary_point_class") or "")
            if stationary_point_class == "minimum_after_optts":
                return "minimum", stationary_point_class
            if stationary_point_class:
                return "ts", stationary_point_class
            return "ts", ts_classification.get("mode_identity")
        if request.role == "intermediate":
            int_classification = dict(candidate.get("int_classification") or {})
            return "minimum", int_classification.get("identity") or request.role
        minimum_classification = dict(candidate.get("minimum_classification") or {})
        return "minimum", minimum_classification.get("identity") or request.role

    def _infer_identity_status(
        self,
        request: StructureRequest,
        candidate: Mapping[str, Any],
    ) -> str:
        if request.kind == "ts":
            stationary_point_class = str((candidate.get("ts_classification") or {}).get("stationary_point_class") or "")
            if stationary_point_class == "valid_target_ts":
                return "role_matched"
            if stationary_point_class == "soft_target_ts":
                return "ambiguous"
            if stationary_point_class in {"minimum_after_optts", "first_order_wrong_mode", "higher_order_saddle"}:
                return "mismatched"
            return "not_checked"
        if request.role == "intermediate":
            identity = str((candidate.get("int_classification") or {}).get("identity") or "")
            if identity in {"distinct_intermediate", "dipolar_intermediate"}:
                return "role_matched"
            if identity in {
                "collapsed_to_precursor",
                "collapsed_to_product",
                "merged_with_ts",
                "ts_energy_degenerate",
                "above_ts_not_minimum",
                "product_like_geometry_energy_inconsistent",
            }:
                return "mismatched"
            if identity in {
                "imaginary_frequency",
                "topology_ambiguous",
                "candidate_intermediate_unverified_energy",
            }:
                return "ambiguous"
            return "not_checked"
        minimum_classification = dict(candidate.get("minimum_classification") or {})
        if int(minimum_classification.get("imaginary_count") or 0) > 0:
            return "ambiguous"
        if str(minimum_classification.get("identity") or "") == "identity_drift":
            return "mismatched"
        return "role_matched"

    @staticmethod
    def _final_ts_frequency_valid(candidate: Mapping[str, Any]) -> Optional[bool]:
        stationary_point_class = str((candidate.get("ts_classification") or {}).get("stationary_point_class") or "")
        if not stationary_point_class:
            return None
        return stationary_point_class == "valid_target_ts"

    @staticmethod
    def _format_theory_level(method: str, basis: str, aux_basis: str = "") -> str:
        level = str(method).strip()
        basis_text = str(basis).strip()
        aux_text = str(aux_basis).strip()
        if basis_text:
            level = f"{level}/{basis_text}"
        if aux_text:
            level = f"{level} ({aux_text})"
        return level

    def _build_failed_preflight_payload(
        self,
        request: StructureRequest,
        preflight: PreflightOutcome,
    ) -> Dict[str, Any]:
        return {
            "id": request.id,
            "role": request.role,
            "kind": request.kind,
            "preflight_status": "failed_preflight",
            "status": "failed",
            "error": preflight.error,
            "pass2_rescue_attempts": [],
            "canonical_attempt_id": None,
            "canonical_xyz": None,
            "sp_status": "not_run",
            "thermochemistry": None,
            "ml_usability": None,
            "irc_status": "not_run",
            "irc_endpoints": None,
            "resolved_kind": None,
            "resolved_identity": None,
            "identity_status": "not_checked",
        }

    def _build_steps_payload(
        self,
        request: StructureRequest,
        preflight: PreflightOutcome,
        pass1: Pass1Outcome,
        pass2: Pass1Outcome,
        pass3: Pass1Outcome,
    ) -> List[Dict[str, Any]]:
        """Reconstruct the normalized S3.0-S3.7 steps[] record for a structure."""
        stage = self.profile.stage

        def status_map(value: str, skip_values: frozenset[str]) -> str:
            if str(value or "").strip().lower() in skip_values:
                return "skipped"
            return str(value or "not_run")

        rescue_attempts = pass2.pass2_rescue_attempts
        rescue_sub_steps: List[Dict[str, Any]] = []
        rescue_status = "not_run"
        if rescue_attempts:
            valid = [
                attempt
                for attempt in rescue_attempts
                if self._candidate_is_valid_for_request(request, attempt)
            ]
            rescue_status = (
                "complete"
                if valid
                else (
                    "partial"
                    if any(str(attempt.get("status")) == "complete" for attempt in rescue_attempts)
                    else "failed"
                )
            )
            for attempt in rescue_attempts:
                strategy = str(attempt.get("strategy") or "")
                if strategy.startswith("S3.4.0") or (
                    attempt.get("failure_type") and not attempt.get("method_family")
                ):
                    sub_code = rescue_sub_code(stage, RESCUE_SUB_DIAGNOSIS)
                else:
                    sub_code = rescue_sub_for_family(
                        stage, str(attempt.get("method_family") or "R1")
                    )
                sub: Dict[str, Any] = {
                    "code": sub_code,
                    "method_family": attempt.get("method_family"),
                    "method": attempt.get("method") or attempt.get("strategy"),
                    "status": attempt.get("status"),
                }
                if attempt.get("attempt_id"):
                    sub["attempt_id"] = attempt["attempt_id"]
                if attempt.get("failure_type"):
                    sub["failure_type"] = attempt["failure_type"]
                if attempt.get("trigger"):
                    sub["trigger"] = attempt["trigger"]
                if attempt.get("error"):
                    sub["error"] = attempt["error"]
                rescue_sub_steps.append(sub)

        properties_extra: Dict[str, Any] = {}
        if pass3.sp_status == "complete" and pass3.sp_energy_hartree is not None:
            properties_extra["sp_energy_hartree"] = pass3.sp_energy_hartree

        steps: List[Dict[str, Any]] = [
            StepRecord(
                code=StepCode.PREFLIGHT.qualified(stage),
                name=StepCode.PREFLIGHT.name,
                status="complete" if preflight.status == "ok" else "failed",
                reason=preflight.error,
            ).to_dict(),
            StepRecord(
                code=StepCode.WARMUP.qualified(stage),
                name=StepCode.WARMUP.name,
                status=status_map(pass1.warmup_status, frozenset({"not_requested"})),
                attempt_id=(
                    pass1.primary_attempt_id if pass1.warmup_used else None
                ),
                reason=pass1.warmup_error,
            ).to_dict(),
            StepRecord(
                code=StepCode.PRIMARY.qualified(stage),
                name=StepCode.PRIMARY.name,
                status=status_map(pass1.opt_status, frozenset({"not_run"})),
                attempt_id=pass1.primary_attempt_id,
                reason=pass1.opt_error,
            ).to_dict(),
            StepRecord(
                code=StepCode.FREQ.qualified(stage),
                name=StepCode.FREQ.name,
                status=status_map(pass1.frequency_status, frozenset({"not_run"})),
            ).to_dict(),
            StepRecord(
                code=StepCode.RESCUE.qualified(stage),
                name=StepCode.RESCUE.name,
                status=rescue_status,
                sub_steps=rescue_sub_steps,
            ).to_dict(),
            StepRecord(
                code=StepCode.CANONICAL.qualified(stage),
                name=StepCode.CANONICAL.name,
                status="complete" if pass3.canonical_xyz is not None else "skipped",
                attempt_id=pass3.canonical_attempt_id,
            ).to_dict(),
            StepRecord(
                code=StepCode.PROPERTIES.qualified(stage),
                name=StepCode.PROPERTIES.name,
                status=status_map(pass3.sp_status, frozenset({"not_run"})),
                extra=properties_extra,
            ).to_dict(),
            StepRecord(
                code=StepCode.FINALIZE.qualified(stage),
                name=StepCode.FINALIZE.name,
                status=pass3.status,
            ).to_dict(),
        ]
        return steps

    def _build_structure_payload(
        self,
        request: StructureRequest,
        preflight: PreflightOutcome,
        pass1: Pass1Outcome,
        pass2: Pass1Outcome,
        pass3: Pass1Outcome,
    ) -> Dict[str, Any]:
        return {
            "id": request.id,
            "role": request.role,
            "kind": request.kind,
            "preflight_status": preflight.status,
            "warmup_status": pass1.warmup_status,
            "warmup_used": pass1.warmup_used,
            "warmup_xyz": str(pass1.warmup_xyz) if pass1.warmup_xyz else None,
            "warmup_constraints": pass1.warmup_constraints,
            "warmup_max_cycles": pass1.warmup_max_cycles,
            "warmup_error": pass1.warmup_error,
            "opt_status": pass1.opt_status,
            "opt_xyz": str(pass1.opt_xyz) if pass1.opt_xyz else None,
            "opt_output": str(pass1.opt_output) if pass1.opt_output else None,
            "opt_energy_hartree": pass1.opt_energy_hartree,
            "opt_initial_hessian": pass1.opt_initial_hessian,
            "opt_error": pass1.opt_error,
            "last_geometry_path": (
                str(pass1.last_geometry_path) if pass1.last_geometry_path else None
            ),
            "stop_geometry_path": (
                str(pass1.stop_geometry_path) if pass1.stop_geometry_path else None
            ),
            "optimization_monitor": self._json_safe(pass1.optimization_monitor),
            "optimization_monitor_error": self._json_safe(pass1.optimization_monitor_error),
            "frequency_status": pass1.frequency_status,
            "frequency_output": (
                str(pass1.frequency_output) if pass1.frequency_output else None
            ),
            "frequency_energy_hartree": pass1.frequency_energy_hartree,
            "frequencies_cm1": pass1.frequencies_cm1,
            "imaginary_frequencies_cm1": pass1.imaginary_frequencies_cm1,
            "frequency_zero_point_energy_hartree": pass1.frequency_zero_point_energy_hartree,
            "frequency_thermal_energy_hartree": pass1.frequency_thermal_energy_hartree,
            "frequency_enthalpy_hartree": pass1.frequency_enthalpy_hartree,
            "frequency_gibbs_free_energy_hartree": pass1.frequency_gibbs_free_energy_hartree,
            "frequency_gibbs_correction_hartree": pass1.frequency_gibbs_correction_hartree,
            "frequency_hessian_path": (
                str(pass1.frequency_hessian_path) if pass1.frequency_hessian_path else None
            ),
            "ts_classification": self._json_safe(pass3.ts_classification),
            "int_classification": self._json_safe(pass3.int_classification),
            "minimum_classification": self._json_safe(pass3.minimum_classification),
            "pass2_rescue_attempts": self._json_safe(pass2.pass2_rescue_attempts),
            "rescue_winning_method": pass3.rescue_winning_method,
            "steps": self._build_steps_payload(request, preflight, pass1, pass2, pass3),
            "canonical_attempt_id": pass3.canonical_attempt_id,
            "canonical_xyz": str(pass3.canonical_xyz) if pass3.canonical_xyz else None,
            "geometry_hash": pass3.geometry_hash,
            "canonical_frequency_status": pass3.canonical_frequency_status,
            "canonical_frequency_output": (
                str(pass3.canonical_frequency_output) if pass3.canonical_frequency_output else None
            )
            if pass3.canonical_frequency_output is not None
            else self._freq_output_fallback(pass3),
            "canonical_frequency_energy_hartree": pass3.canonical_frequency_energy_hartree,
            "canonical_frequencies_cm1": pass3.canonical_frequencies_cm1,
            "canonical_imaginary_frequencies_cm1": pass3.canonical_imaginary_frequencies_cm1,
            "canonical_zero_point_energy_hartree": pass3.canonical_zero_point_energy_hartree,
            "canonical_thermal_energy_hartree": pass3.canonical_thermal_energy_hartree,
            "canonical_enthalpy_hartree": pass3.canonical_enthalpy_hartree,
            "canonical_gibbs_free_energy_hartree": pass3.canonical_gibbs_free_energy_hartree,
            "canonical_gibbs_correction_hartree": pass3.canonical_gibbs_correction_hartree,
            "canonical_hessian_path": (
                str(pass3.canonical_hessian_path) if pass3.canonical_hessian_path else None
            ),
            "sp_status": pass3.sp_status,
            "sp_energy_hartree": pass3.sp_energy_hartree,
            "sp_output": str(pass3.sp_output) if pass3.sp_output else None,
            "thermochemistry": self._json_safe(pass3.thermochemistry),
            "ml_usability": self._json_safe(pass3.ml_usability),
            "usable_for_ml": (
                compute_legacy_usable_for_ml(pass3.ml_usability)
                if pass3.ml_usability is not None
                else None
            ),
            "irc_status": pass3.irc_status,
            "irc_endpoints": self._json_safe(pass3.irc_endpoints),
            "resolved_kind": pass3.resolved_kind,
            "resolved_identity": pass3.resolved_identity,
            "identity_status": pass3.identity_status,
            "status": pass3.status,
            "error": pass3.error,
            "attempt_history": pass1.attempt_history,
            "current_attempt_id": pass1.current_attempt_id,
            "forming_bonds": [list(pair) for pair in preflight.forming_bonds],
            "structure_id": request.structure_id,
            "variant_id": request.variant_id,
            "branch_id": request.branch_id,
            "pathway_id": request.pathway_id,
            "parent_structure_id": request.parent_structure_id,
            "parent_structure": self._json_safe(request.parent_structure),
            "source_stage": request.source_stage,
            "seed_state": request.seed_state,
            "input_xyz": str(preflight.input_xyz) if preflight.input_xyz else None,
            "charge": preflight.charge,
            "multiplicity": preflight.multiplicity,
            "atom_mapping": str(request.atom_mapping) if request.atom_mapping else None,
            "mapping_audit": self._json_safe(request.mapping_audit),
            "mapping_required": request.mapping_required,
            "s1_manifest": self._json_safe(request.s1_manifest),
            "s1_ensemble_thermodynamics": self._json_safe(
                request.s1_ensemble_thermodynamics
            ),
            "s1_thermochemistry_status": request.s1_thermochemistry_status,
            "ensemble_thermochemistry_correction_hartree": (
                request.ensemble_thermochemistry_correction_hartree
            ),
        }

    def _write_manifest(
        self,
        output_dir: Path,
        structures_payload: List[Dict[str, Any]],
        *,
        archived_to: Path | None,
        schedule_record: Dict[str, Any],
        source_requests: Sequence[StructureRequest],
        partial_rerun_ids: Optional[Sequence[str]] = None,
        preserved_structure_ids: Optional[Sequence[str]] = None,
    ) -> Path:
        manifest_path = output_dir / "manifest.json"
        summary: Dict[str, int] = {}
        for structure in structures_payload:
            status = str(structure.get("status", "unknown"))
            summary[status] = summary.get(status, 0) + 1

        provenance = build_provenance(
            run_id=self.run_id,
            schema_version=REFINEMENT_MANIFEST_V1,
            protocol_version=None,
            parent_manifest_paths=self._relevant_parent_manifest_paths(),
            atom_mapping_payload=None,
            forming_bonds=self._aggregate_forming_bonds(source_requests),
            variant_manifest_path=self._variant_parent_manifest_path(),
            extra={"phase": StepCode.FINALIZE.qualified(self.profile.stage)},
        ).to_dict()

        return write_refinement_manifest(
            manifest_path,
            stage=self.profile.stage,
            fidelity=self.profile.fidelity,
            profile_id=self.profile.profile_id,
            structures=structures_payload,
            run_id=self.run_id,
            extra={
                "stale_outputs_archived_to": str(archived_to) if archived_to else None,
                "scheduling": schedule_record,
                "summary": summary,
                "provenance": provenance,
                "partial_rerun": {
                    "enabled": partial_rerun_ids is not None,
                    "rerun_structure_ids": list(partial_rerun_ids or ()),
                    "preserved_structure_ids": list(preserved_structure_ids or ()),
                },
            },
        )

    @staticmethod
    def _load_existing_structure_payloads(output_dir: Path) -> Dict[str, Dict[str, Any]]:
        """Return prior structure records keyed by id without trusting their status."""

        manifest_path = Path(output_dir) / "manifest.json"
        if not manifest_path.exists():
            return {}
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise RuntimeError(f"Unable to read existing refinement manifest {manifest_path}: {exc}") from exc
        rows = payload.get("structures", [])
        if not isinstance(rows, list):
            raise RuntimeError(f"Existing refinement manifest {manifest_path} has invalid structures payload")
        return {
            str(row.get("id") or row.get("structure_id")): copy.deepcopy(dict(row))
            for row in rows
            if isinstance(row, Mapping) and str(row.get("id") or row.get("structure_id") or "").strip()
        }

    @staticmethod
    def _existing_structure_is_reusable(payload: Optional[Mapping[str, Any]]) -> bool:
        """Accept only complete records whose essential produced files remain present."""

        if not payload or str(payload.get("status") or "").lower() != "complete":
            return False
        required = [payload.get("canonical_xyz"), payload.get("opt_output")]
        if str(payload.get("canonical_frequency_status") or "") == "complete":
            freq_output = payload.get("canonical_frequency_output")
            if not freq_output:
                freq_output = RefinementEngine._freq_output_fallback_from_payload(payload)
            required.append(freq_output)
        if str(payload.get("sp_status") or "") == "complete":
            required.append(payload.get("sp_output"))
        return all(value and Path(str(value)).is_file() for value in required)

    @staticmethod
    def _freq_output_fallback_from_payload(payload: Mapping[str, Any]) -> Optional[str]:
        hessian = payload.get("canonical_hessian_path")
        if not hessian or not Path(str(hessian)).is_file():
            return None
        hess = Path(str(hessian))
        out = hess.with_suffix(".out")
        return str(out) if out.exists() else str(hess)

    def _preserved_int_collapse_ids(
        self,
        requests: Sequence[StructureRequest],
        existing_payloads: Mapping[str, Dict[str, Any]],
        rerun_ids: Set[str],
    ) -> Set[str]:
        """Reclassify preserved INTs and return those collapsed to product.

        INT structures preserved from a previous run still carry the old
        ``distinct_intermediate`` identity label (computed before the
        ``parent_structure`` references were injected).  In rescue-only mode
        we re-run the identity classifier on the preserved geometry so that a
        collapsed INT is routed into the rescue matrix (F7 -> IRC mid-point
        recovery) instead of being silently carried forward.
        """
        collapsed: Set[str] = set()
        for request in requests:
            if request.role != "intermediate" or request.id in rerun_ids:
                continue
            payload = existing_payloads.get(request.id)
            if payload is None:
                continue
            canonical_xyz = payload.get("canonical_xyz")
            if not canonical_xyz or not Path(str(canonical_xyz)).is_file():
                continue
            int_classification = self._classify_int_v2(
                request,
                opt_xyz=Path(str(canonical_xyz)),
                frequencies_cm1=payload.get("imaginary_frequencies_cm1") or (),
                energy_int_hartree=(
                    float(payload["opt_energy_hartree"])
                    if payload.get("opt_energy_hartree")
                    else self._read_orca_final_energy(
                        Path(str(canonical_xyz)).with_name("output.out")
                    )
                ),
                converged=True if str(payload.get("opt_status")) == "complete" else None,
            )
            if str(int_classification.get("identity") or "") in {
                "collapsed_to_product",
                "product_like_geometry_energy_inconsistent",
            }:
                payload["int_classification"] = int_classification
                collapsed.add(request.id)
        return collapsed

    def _partition_rescue_replay(
        self,
        requests: Sequence[StructureRequest],
        existing_payloads: Mapping[str, Dict[str, Any]],
    ) -> Tuple[Set[str], Set[str]]:
        replay_ids: Set[str] = set()
        fallback_ids: Set[str] = set()
        for request in requests:
            payload = existing_payloads.get(request.id)
            if payload is None:
                fallback_ids.add(request.id)
                continue
            if str(payload.get("opt_status") or "").lower() != "failed":
                if (
                    request.role == "intermediate"
                    and str(
                        (payload.get("int_classification") or {}).get("identity") or ""
                    )
                    == "collapsed_to_product"
                ):
                    replay_ids.add(request.id)
                    continue
                fallback_ids.add(request.id)
                continue
            history = payload.get("attempt_history") or []
            if not history:
                fallback_ids.add(request.id)
                continue
            geometry = self._replay_geometry(payload)
            if geometry is None:
                fallback_ids.add(request.id)
                continue
            replay_ids.add(request.id)
        return replay_ids, fallback_ids

    @staticmethod
    def _replay_geometry(payload: Mapping[str, Any]) -> Optional[Path]:
        for key in ("stop_geometry_path", "last_geometry_path"):
            value = payload.get(key)
            if value and Path(str(value)).is_file():
                return Path(str(value))
        return None

    def _replay_failed_outcomes(
        self,
        existing_payloads: Mapping[str, Dict[str, Any]],
        replay_ids: Set[str],
    ) -> List[Pass1Outcome]:
        """Reconstruct Pass1Outcome from prior payloads, carrying only the fields the rescue matrix consumes."""
        outcomes: List[Pass1Outcome] = []
        for request_id in sorted(replay_ids):
            request = self._requests_by_id[request_id]
            payload = existing_payloads.get(request_id) or {}
            history = payload.get("attempt_history") or []
            outcome = Pass1Outcome(
                structure_id=request_id,
                role=str(payload.get("role") or request.role),
                kind=str(payload.get("kind") or request.kind),
            )
            outcome.opt_status = str(payload.get("opt_status") or "failed")
            outcome.opt_error = payload.get("opt_error")
            outcome.opt_output = (
                Path(str(payload.get("opt_output"))) if payload.get("opt_output") else None
            )
            outcome.imaginary_frequencies_cm1 = list(
                payload.get("imaginary_frequencies_cm1") or ()
            )
            outcome.attempt_history = list(history)
            outcome.current_attempt_id = payload.get("current_attempt_id")
            outcome.primary_attempt_id = payload.get("primary_attempt_id")
            stop = payload.get("stop_geometry_path")
            last = payload.get("last_geometry_path")
            outcome.stop_geometry_path = (
                Path(str(stop)) if stop and Path(str(stop)).is_file() else None
            )
            outcome.last_geometry_path = (
                Path(str(last)) if last and Path(str(last)).is_file() else None
            )
            outcome.ts_classification = copy.deepcopy(payload.get("ts_classification"))
            outcome.int_classification = copy.deepcopy(payload.get("int_classification"))
            outcome.minimum_classification = copy.deepcopy(
                payload.get("minimum_classification")
            )
            canonical_xyz = payload.get("canonical_xyz")
            if (
                request.role == "intermediate"
                and canonical_xyz
                and Path(str(canonical_xyz)).is_file()
            ):
                outcome.opt_xyz = Path(str(canonical_xyz))
                outcome.int_classification = self._classify_int_v2(
                    request,
                    opt_xyz=Path(str(canonical_xyz)),
                    frequencies_cm1=outcome.imaginary_frequencies_cm1,
                    energy_int_hartree=(
                        float(payload.get("opt_energy_hartree"))
                        if payload.get("opt_energy_hartree")
                        else self._read_orca_final_energy(
                            Path(str(canonical_xyz)).with_name("output.out")
                        )
                    ),
                    converged=True if str(payload.get("opt_status")) == "complete" else None,
                )
            outcome.status = "failed" if outcome.opt_status != "complete" else "complete"
            logger.info(
                "[%s] rescue-only replay: %s opt_status=%s int_identity=%s",
                self.profile.stage,
                request_id,
                outcome.opt_status,
                (outcome.int_classification or {}).get("identity"),
            )
            outcomes.append(outcome)
        return outcomes

    def _relevant_parent_manifest_paths(self) -> Dict[str, Optional[Path]]:
        stage = self.profile.stage
        if stage == "S3":
            return {"s2": self.parent_manifest_paths.get("s2")}
        if stage == "S4":
            return {"s3": self.parent_manifest_paths.get("s3")}
        return dict(self.parent_manifest_paths)

    def _variant_parent_manifest_path(self) -> Optional[Path]:
        stage = self.profile.stage
        if stage == "S3":
            return self.parent_manifest_paths.get("s2")
        if stage == "S4":
            return self.parent_manifest_paths.get("s3")
        return None

    @staticmethod
    def _aggregate_forming_bonds(
        structures: Iterable[StructureRequest | Mapping[str, Any]],
    ) -> Optional[List[List[int]]]:
        aggregated = [
            [int(pair[0]), int(pair[1])]
            for structure in structures
            for pair in (
                structure.forming_bonds
                if isinstance(structure, StructureRequest)
                else (structure.get("forming_bonds") or [])
            )
        ]
        return aggregated or None

    @staticmethod
    def _coerce_optional_path(value: Any) -> Optional[Path]:
        if value in {None, ""}:
            return None
        return Path(value)

    @staticmethod
    def _coerce_optional_str(value: Any) -> Optional[str]:
        if value is None:
            return None
        text = str(value).strip()
        return text or None

    @staticmethod
    def _coerce_optional_float(value: Any) -> Optional[float]:
        if value is None:
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _coerce_forming_bonds(value: Any) -> List[Tuple[int, int]]:
        bonds: List[Tuple[int, int]] = []
        for raw_pair in value or []:
            if len(raw_pair) < 2:
                raise ValueError(f"Invalid forming-bond record: {raw_pair!r}")
            bonds.append((int(raw_pair[0]), int(raw_pair[1])))
        return bonds

    @staticmethod
    def _resolve_input_xyz(request: StructureRequest) -> Path:
        input_xyz = Path(request.input_xyz)
        if input_xyz.exists():
            return input_xyz
        if request.fallback_xyz is not None and Path(request.fallback_xyz).exists():
            return Path(request.fallback_xyz)
        raise ValueError(f"input_xyz not found: {request.input_xyz}")

    @staticmethod
    def _read_atom_count(input_xyz: Path) -> int:
        coordinates, _symbols = read_xyz(Path(input_xyz))
        return int(len(coordinates))

    def _write_provenance(
        self,
        request: StructureRequest,
        structure_dir: Path,
        input_xyz: Path,
    ) -> None:
        write_json_atomic(
            structure_dir / "provenance.json",
            {
                "structure_id": request.id,
                "kind": request.kind,
                "run_id": self.run_id
                or self.config.get("run_id")
                or dict(self.config.get("run", {}) or {}).get("id"),
                "profile_id": self.profile.profile_id,
                "stage": self.profile.stage,
                "input_hashes": {
                    "input_xyz_sha256": self._file_sha256(input_xyz),
                    "atom_mapping_sha256": self._file_sha256(request.atom_mapping),
                },
                "forming_bonds_hash": self._payload_sha256(request.forming_bonds),
                "input_xyz": str(input_xyz),
                "fallback_xyz": str(request.fallback_xyz) if request.fallback_xyz else None,
                "original_seed_xyz": (
                    str(request.original_seed_xyz) if request.original_seed_xyz else None
                ),
                "atom_mapping": str(request.atom_mapping) if request.atom_mapping else None,
                "source_stage": request.source_stage,
                "seed_state": request.seed_state,
            },
        )

    def _warmup_applies(self, request: StructureRequest) -> bool:
        warmup_cfg = (
            dict(self.config.get("refinement", {}) or {})
            .get("common", {})
            .get("workflow", {})
            .get("warmup", {})
        )
        enabled_roles = {
            str(role).strip().lower()
            for role in warmup_cfg.get("enabled_roles", ("intermediate", "ts"))
        }
        return request.role in enabled_roles and bool(request.forming_bonds)

    @staticmethod
    def _build_bond_constraints(
        source_xyz: Path,
        forming_bonds: Sequence[Tuple[int, int]],
    ) -> Tuple[Tuple[Tuple[int, int, Optional[float]], ...], List[Dict[str, Any]]]:
        coordinates, _symbols = read_xyz(Path(source_xyz))
        atom_count = len(coordinates)
        constraints: List[Tuple[int, int, Optional[float]]] = []
        details: List[Dict[str, Any]] = []
        for atom_i, atom_j in forming_bonds:
            if atom_i < 0 or atom_j < 0 or atom_i >= atom_count or atom_j >= atom_count:
                raise ValueError(
                    f"forming_bond {(atom_i, atom_j)} out of range for atom_count={atom_count}"
                )
            distance = float(GeometryUtils.calculate_distance(coordinates, atom_i, atom_j))
            constraints.append((int(atom_i), int(atom_j), distance))
            details.append(
                {
                    "atoms": [int(atom_i), int(atom_j)],
                    "target_distance_angstrom": distance,
                }
            )
        return tuple(constraints), details

    @staticmethod
    def _copy_file(source: Path | None, destination: Path) -> Path | None:
        if source is None:
            return None
        source_path = Path(source)
        if not source_path.exists():
            return None
        destination.parent.mkdir(parents=True, exist_ok=True)
        try:
            if source_path.resolve() == destination.resolve():
                return destination
        except OSError:
            if source_path == destination:
                return destination
        shutil.copy2(source_path, destination)
        return destination

    @staticmethod
    def _file_sha256(path: Path | None) -> str | None:
        if path is None or not Path(path).exists():
            return None
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    @staticmethod
    def _payload_sha256(payload: Any) -> str | None:
        if payload is None:
            return None
        serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(serialized.encode("utf-8")).hexdigest()

    @staticmethod
    def _attempt_runtime_payload(
        record: AttemptRecord,
        spec: QCJobSpec,
        result: QCJobResult,
        started_at: datetime,
        finished_at: datetime,
    ) -> Dict[str, Any]:
        runtime_seconds = max(0.0, (finished_at - started_at).total_seconds())
        return {
            "started_at": started_at.astimezone(timezone.utc).isoformat(),
            "finished_at": finished_at.astimezone(timezone.utc).isoformat(),
            "runtime_seconds": runtime_seconds,
            "returncode": result.returncode,
            "orca_version": result.orca_version,
            "nproc": spec.nproc,
            "mem": spec.memory,
            "theory_signature": {
                "engine": spec.engine,
                "task": spec.task,
                "method": spec.method,
                "basis": spec.basis,
                "aux_basis": spec.aux_basis,
                "route": spec.route,
                "route_extras": spec.route_extras,
                "grid": spec.grid,
                "scf": spec.scf,
                "scf_maxiter": spec.scf_maxiter,
                "max_cycles": spec.max_cycles,
                "charge": spec.charge,
                "multiplicity": spec.multiplicity,
                "solvent": spec.solvent,
                "solvent_model": spec.solvent_model,
                "initial_hessian": spec.initial_hessian,
                "hessian_filename": spec.hessian_filename,
                "ts_mode": spec.ts_mode,
                "recalc_hessian": spec.recalc_hessian,
                "trust": spec.trust,
            },
            "attempt_kind": record.attempt_kind,
            "retry_level": record.retry_level,
            "status": result.status,
            "timed_out": bool(result.timed_out),
            "optimization_converged": result.optimization_converged,
            "stop_reason": result.stop_reason,
        }

    @staticmethod
    def _attempt_error_summary(result: QCJobResult) -> Dict[str, Any] | None:
        if result.status == "complete":
            return None
        return {
            "failure_type": result.failure_type,
            "evidence_lines": list(result.failure_evidence_lines or ()),
            "summary": result.error or result.failure_type or "QC job failed",
            "optimization_converged": result.optimization_converged,
            "stop_reason": result.stop_reason,
        }

    @staticmethod
    def _resolve_generated_input_file(
        result: QCJobResult,
        attempt_dir: Path,
    ) -> Path | None:
        output_file = Path(result.output_file) if result.output_file else None
        if output_file is not None:
            candidate = output_file.with_suffix(".inp")
            if candidate.exists():
                return candidate
        candidates = sorted(
            path
            for path in attempt_dir.iterdir()
            if path.is_file() and path.suffix.lower() == ".inp" and path.name != "input.inp"
        )
        return candidates[-1] if candidates else None

    @staticmethod
    def _resolve_generated_output_file(
        result: QCJobResult,
        attempt_dir: Path,
    ) -> Path | None:
        output_file = Path(result.output_file) if result.output_file else None
        if output_file is not None and output_file.exists():
            return output_file
        candidates = sorted(
            path
            for path in attempt_dir.iterdir()
            if path.is_file() and path.suffix.lower() == ".out" and path.name != "output.out"
        )
        return candidates[-1] if candidates else None

    @staticmethod
    def _resolve_last_geometry_text(result: QCJobResult, input_xyz: Path) -> str | None:
        if result.output_xyz is not None:
            output_xyz = Path(result.output_xyz)
            if output_xyz.exists():
                try:
                    text = output_xyz.read_text(encoding="utf-8")
                    return text if text.endswith("\n") else f"{text}\n"
                except OSError:
                    pass
        return extract_last_orca_geometry_xyz_text(result.output_file, input_xyz)

    def _execute_attempt(
        self,
        structure_id: str,
        recorder: AttemptRecorder,
        *,
        kind: str,
        retry_level: int,
        runner: Callable[..., QCJobResult],
        spec: QCJobSpec,
        input_xyz: Path,
    ) -> Dict[str, Any]:
        record = recorder.next_attempt(kind, retry_level)
        self._copy_file(Path(input_xyz), record.input_xyz_path)
        started_at = datetime.now(timezone.utc)
        task_name = self._progress_task_name(kind)
        task_payload = {
            "engine": spec.engine,
            "method": spec.method,
            "basis": spec.basis,
            "solvent": spec.solvent,
            "solvent_model": spec.solvent_model,
            "attempt_id": record.attempt_id,
            "retry_level": retry_level,
            "output": str(record.directory),
        }
        self._safe_calculator_event(f"{task_name}_started", structure_id, task_payload)
        optimization_monitor = self._optimization_monitor_for(spec=spec, work_dir=record.directory)
        try:
            result = self._run_job_with_subprocess_tracking(
                structure_id,
                record.directory,
                runner,
                spec,
                Path(input_xyz),
                record.directory,
                self.config,
                optimization_monitor=optimization_monitor,
            )
        except Exception as exc:
            self._safe_calculator_event(
                f"{task_name}_finished",
                structure_id,
                {**task_payload, "status": "failed", "error": str(exc)},
            )
            raise
        finished_at = datetime.now(timezone.utc)

        generated_input = self._resolve_generated_input_file(result, record.directory)
        generated_output = self._resolve_generated_output_file(result, record.directory)
        self._copy_file(generated_input, record.input_inp_path)
        copied_output = self._copy_file(generated_output, record.output_out_path)
        if copied_output is None and not record.output_out_path.exists():
            record.output_out_path.touch(exist_ok=True)
        record.stderr_log_path.write_text(result.stderr_text or "", encoding="utf-8")

        last_geometry_text = self._resolve_last_geometry_text(result, Path(input_xyz))
        monitor_payload = dict((result.extra or {}).get("optimization_monitor") or {})
        monitor_error_payload = dict(
            (result.extra or {}).get("optimization_monitor_error") or {}
        )
        stop_geometry_path: Path | None = None
        if monitor_payload:
            monitor_path = record.directory / "monitor.json"
            if last_geometry_text is not None:
                stop_geometry_path = record.directory / "stop.xyz"
                stop_geometry_path.write_text(last_geometry_text, encoding="utf-8")
                monitor_payload["stop_geometry"] = str(stop_geometry_path)
            write_json_atomic(monitor_path, monitor_payload)
        if monitor_error_payload:
            write_json_atomic(record.directory / "monitor_error.json", monitor_error_payload)
        runtime = self._attempt_runtime_payload(record, spec, result, started_at, finished_at)
        error_summary = self._attempt_error_summary(result)
        recorder.finalize_attempt(
            record,
            runtime=runtime,
            error_summary=error_summary,
            last_geometry_xyz=last_geometry_text,
        )
        self._safe_calculator_event(
            f"{task_name}_finished",
            structure_id,
            {
                **task_payload,
                "status": result.status,
                "energy_hartree": result.energy_hartree,
                "output": str(record.output_out_path),
                "error": result.error,
            },
        )
        return {
            "record": record,
            "result": result,
            "runtime": runtime,
            "error_summary": error_summary,
            "last_geometry_path": (
                record.last_geometry_xyz_path if last_geometry_text is not None else None
            ),
            "stop_geometry_path": stop_geometry_path,
            "optimization_monitor": monitor_payload or None,
            "optimization_monitor_error": monitor_error_payload or None,
            "history_entry": {
                "attempt_id": record.attempt_id,
                "retry_level": record.retry_level,
                "status": result.status,
                "failure_type": result.failure_type,
                "optimization_converged": result.optimization_converged,
                "stop_reason": result.stop_reason,
                "runtime_seconds": runtime["runtime_seconds"],
                "directory": str(record.directory),
                "optimization_monitor": monitor_payload or None,
                "optimization_monitor_error": monitor_error_payload or None,
            },
        }

    @staticmethod
    def _progress_task_name(kind: str) -> str:
        normalized = str(kind).strip().lower()
        if "freq" in normalized:
            return "frequency"
        if normalized == "sp":
            return "single_point"
        if normalized == "irc":
            return "irc"
        if normalized.startswith("warmup"):
            return "warmup"
        return "optimization"

    @staticmethod
    def _latest_existing_path(*paths: Path | None) -> Path | None:
        for candidate in paths:
            if candidate is not None and Path(candidate).exists():
                return Path(candidate)
        return None

    def _populate_legacy_attempt_dir(
        self,
        legacy_dir: Path,
        outcome: Dict[str, Any],
        *,
        xyz_target_name: str | None = None,
    ) -> None:
        if legacy_dir.exists():
            shutil.rmtree(legacy_dir)
        legacy_dir.mkdir(parents=True, exist_ok=True)
        record = outcome["record"]
        for source_path in (
            record.input_inp_path,
            record.input_xyz_path,
            record.output_out_path,
            record.stderr_log_path,
            record.runtime_json_path,
            record.error_summary_path,
            record.last_geometry_xyz_path,
        ):
            if source_path is None or not Path(source_path).exists():
                continue
            self._copy_file(Path(source_path), legacy_dir / Path(source_path).name)
        if xyz_target_name:
            preferred_xyz = self._latest_existing_path(
                Path(outcome["result"].output_xyz) if outcome["result"].output_xyz else None,
                record.last_geometry_xyz_path,
            )
            if preferred_xyz is not None:
                self._copy_file(preferred_xyz, legacy_dir / xyz_target_name)

    @staticmethod
    def _load_atom_mapping(atom_mapping_path: Path | None) -> Optional[Dict[int, int]]:
        if atom_mapping_path is None or not Path(atom_mapping_path).exists():
            return None
        try:
            payload = json.loads(Path(atom_mapping_path).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        source = payload
        if isinstance(payload, Mapping):
            for candidate_key in ("map_to_xyz", "smiles_to_xyz"):
                candidate = payload.get(candidate_key)
                if isinstance(candidate, Mapping):
                    source = candidate
                    break
        if not isinstance(source, Mapping):
            return None
        mapping: Dict[int, int] = {}
        for raw_key, raw_value in source.items():
            try:
                mapping[int(raw_key)] = int(raw_value)
            except (TypeError, ValueError):
                continue
        return mapping or None

    @staticmethod
    def _resolve_reference_path(
        request: StructureRequest,
        role: str,
    ) -> Optional[Path]:
        parent_structure = request.parent_structure or {}
        direct_value = parent_structure.get(f"{role}_ref")
        if direct_value:
            return Path(direct_value)
        for key in (
            f"{role}_xyz",
            f"{role}_opt_xyz",
            f"{role}_input_xyz",
        ):
            value = parent_structure.get(key)
            if value:
                return Path(value)
        reference_payload = parent_structure.get(role)
        if isinstance(reference_payload, Mapping):
            for key in ("opt_xyz", "input_xyz", "xyz"):
                value = reference_payload.get(key)
                if value:
                    return Path(value)
        return None

    @staticmethod
    def _apply_warmup_payload(outcome: Pass1Outcome, payload: Dict[str, Any]) -> None:
        outcome.warmup_status = str(payload.get("warmup_status", outcome.warmup_status))
        outcome.warmup_used = bool(payload.get("warmup_used", False))
        warmup_xyz = payload.get("warmup_xyz")
        outcome.warmup_xyz = Path(warmup_xyz) if warmup_xyz is not None else None
        outcome.warmup_constraints = list(payload.get("warmup_constraints", []))
        outcome.warmup_max_cycles = payload.get("warmup_max_cycles")
        outcome.warmup_error = payload.get("warmup_error")

    @staticmethod
    def _json_safe(value: Any) -> Any:
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, Mapping):
            return {str(key): RefinementEngine._json_safe(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [RefinementEngine._json_safe(item) for item in value]
        return value


__all__ = ["RefinementEngine", "StructureRequest", "PreflightOutcome", "Pass1Outcome"]
