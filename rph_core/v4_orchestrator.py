"""V4 S0–S4 orchestrator with a fixed CENSO-LITE S1."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from rph_core.steps.conformer_search.censo_lite import (
    CensoLiteEngine,
    S1_MANIFEST_SCHEMA_VERSION,
)
from rph_core.steps.mechanism_classifier.s0_artifact_writer import write_s0_artifacts
from rph_core.steps.mechanism_classifier.s0_record import S0ReactionRecord, load_s0_reaction_record
from rph_core.steps.step2_retro.peb_scanner import PEBScanner
from rph_core.steps.step3_lowlevel import LowLevelEngine
from rph_core.steps.step4_highlevel import HighLevelEngine
from rph_core.utils.config_loader import load_config
from rph_core.utils.log_manager import setup_v4_logging
from rph_core.utils.shared_console import get_console
from rph_core.utils.stage_progress import StageProgressReporter
from rph_core.utils.ui_reporter import RichReporter
from rph_core.utils.v4_checkpoint import V4Checkpoint

logger = logging.getLogger(__name__)
SCHEMA_VERSION = "rph_v4_s0_s4_v1"


class ResumeCheckpointError(RuntimeError):
    """Raised before QC starts when requested resume state is unsafe."""


class V4Orchestrator:
    """Run S0–S4 using manifest-based handoffs and per-structure failures."""

    def __init__(
        self,
        config_path: Optional[Path] = None,
        color: bool = True,
        ui_mode: Optional[str] = None,
        resume_policy: Optional[str] = None,
        start_from: Optional[str] = None,
        recompute_from: Optional[str] = None,
    ):
        self.config = load_config(config_path) if config_path else load_config()
        self.color = bool(color)
        self.ui_mode = ui_mode
        run_config = self.config.get("run", {}) or {}
        self.resume_policy = str(
            resume_policy or run_config.get("resume_policy", "strict")
        ).strip().lower()
        self.start_from = str(start_from).strip().lower() if start_from else None
        self.recompute_from = (
            str(recompute_from).strip().lower() if recompute_from else None
        )
        if self.resume_policy not in {
            "strict",
            "recompute",
            "use-existing-upstream",
        }:
            raise ValueError(
                "resume_policy must be strict, recompute, or use-existing-upstream"
            )
        if self.start_from and self.start_from not in V4Checkpoint.stage_order:
            raise ValueError("start_from must be one of s0, s1, s2, s3, s4")
        if self.recompute_from and self.recompute_from not in V4Checkpoint.stage_order:
            raise ValueError("recompute_from must be one of s0, s1, s2, s3, s4")
        if self.start_from and self.resume_policy != "use-existing-upstream":
            raise ValueError(
                "start_from requires resume_policy=use-existing-upstream"
            )
        if self.resume_policy == "use-existing-upstream" and not self.start_from:
            raise ValueError(
                "resume_policy=use-existing-upstream requires start_from"
            )
        if self.start_from and self.recompute_from:
            raise ValueError("start_from and recompute_from cannot be combined")
        protocol = str((self.config.get("step1", {}) or {}).get("protocol", "censo_lite"))
        if protocol != "censo_lite":
            raise ValueError("RPH V4 only supports step1.protocol=censo_lite")

    def run(
        self,
        s0_record: S0ReactionRecord,
        work_dir: Path,
        stop_after: str = "s4",
    ) -> Dict[str, Any]:
        """Run the pipeline and always restore the terminal UI."""

        self._ui_reporter: RichReporter | None = None
        try:
            result = self._run_pipeline(s0_record, work_dir, stop_after=stop_after)
        except KeyboardInterrupt:
            if self._ui_reporter is not None:
                self._ui_reporter.finish("interrupted")
            raise
        except Exception:
            if self._ui_reporter is not None:
                self._ui_reporter.finish("failed")
            raise
        else:
            if self._ui_reporter is not None:
                self._ui_reporter.finish("complete")
            return result
        finally:
            if self._ui_reporter is not None:
                self._ui_reporter.close()

    def _run_pipeline(
        self,
        s0_record: S0ReactionRecord,
        work_dir: Path,
        stop_after: str = "s4",
    ) -> Dict[str, Any]:
        stop_after = str(stop_after or "s4").strip().lower()
        if stop_after not in {"s0", "s1", "s2", "s3", "s4"}:
            raise ValueError("stop_after must be one of: s0, s1, s2, s3, s4")
        work_dir = Path(work_dir).resolve()
        work_dir.mkdir(parents=True, exist_ok=True)
        self._reused_signature_overrides: Dict[Tuple[str, Optional[str]], str] = {}
        self._write_schema(work_dir)
        checkpoint = V4Checkpoint(work_dir)
        checkpoint.write_config_snapshot(self.config)
        checkpoint.set_reaction_id(s0_record.rx_id)
        recompute_from = getattr(self, "recompute_from", None)
        if recompute_from:
            checkpoint.invalidate_from(recompute_from)
            logger.warning(
                "[V4] Explicitly invalidated checkpoint from %s", recompute_from
            )
        product_smiles = s0_record.product_smiles
        condition_fields, condition_signature = self._condition_metadata(s0_record)
        progress_fields = {
            "reaction_id": s0_record.rx_id,
            "condition_signature": condition_signature,
            "conditions": condition_fields,
            "resources": dict(self.config.get("resources", {}) or {}),
            "ui": dict(self.config.get("ui", {}) or {}),
        }
        requested_ui_mode = str(
            getattr(self, "ui_mode", None)
            or (self.config.get("ui", {}) or {}).get("mode", "auto")
        ).strip().lower()
        if requested_ui_mode == "auto":
            requested_ui_mode = "dashboard" if bool(getattr(self, "color", True)) else "compact"
        ui_config = dict(self.config.get("ui", {}) or {})
        ui_reporter = RichReporter(
            get_console(),
            color=bool(getattr(self, "color", True)),
            mode=requested_ui_mode,
            dashboard_config=dict(ui_config.get("dashboard") or {}),
            log_path=str(work_dir / "rph_v4.log"),
        )
        self._ui_reporter = ui_reporter
        ui_reporter.start()
        s0_reporter = StageProgressReporter(
            work_dir / "S0_Mechanism",
            "S0",
            default_fields=progress_fields,
            event_callback=ui_reporter.event_callback,
        )
        s0_reporter.emit_stage_event("s0_started", rx_id=s0_record.rx_id)

        # A trusted reaction record is the original S0 contract: mapped
        # mechanistic bonds are resolved before geometry generation and remain
        # authoritative through S2.  Do not turn map numbers into XYZ indices.
        s0_manifest = work_dir / "S0_Mechanism" / "mechanism.json"
        s0_signature_payload = {
            "signature_schema": "s0_signature_v2",
            "stage": "s0",
            "source": "trusted_reaction_record",
            "condition": condition_signature,
            "record": s0_record.signature_payload(),
        }
        s0_signature = V4Checkpoint.signature(s0_signature_payload)
        s0_reused = self._checkpoint_reusable(
            checkpoint,
            "s0",
            s0_signature,
            s0_manifest,
            signature_payload=s0_signature_payload,
        )
        if s0_reused:
            s0_data = self._read_json(s0_manifest)
            forming_bonds = tuple(tuple(pair) for pair in s0_data.get("forming_bonds", []))
            logger.info("[V4] Reusing record-driven S0 manifest")
        else:
            forming_bonds = s0_record.forming_bonds
            s0_manifest = self._write_s0_from_record(work_dir, s0_record)
            checkpoint.mark_s0(
                s0_signature,
                s0_manifest,
                signature_payload=s0_signature_payload,
            )
            logger.info("[V4] S0 restored from trusted reaction record rx_id=%s", s0_record.rx_id)

        variants = self._load_variants(work_dir, product_smiles)
        s0_reporter.science_summary(
            mechanism_valid=True,
            forming_bonds=[list(pair) for pair in forming_bonds],
            variants=len(variants),
            source="trusted_reaction_record",
        )
        s0_reporter.emit_stage_event(
            "s0_completed",
            status="completed",
            reused=s0_reused,
            variants=[str(variant["variant_id"]) for variant in variants],
        )
        precursor_smiles = getattr(s0_record, "canonical_precursor_smiles", "") or ""
        precursor_s1_result: Optional[Dict[str, Any]] = None
        all_s1_results: Dict[str, Dict[str, Any]] = {}
        all_s2_results: Dict[str, Dict[str, Any]] = {}
        if stop_after == "s0":
            self._write_run_manifest(
                work_dir,
                s0_record,
                variants,
                all_s1_results,
                all_s2_results,
                "s0",
                s0_manifest,
            )
            return self._write_pipeline_result(
                work_dir,
                stop_after,
                {"s0": s0_manifest},
            )

        s1_reporter = StageProgressReporter(
            work_dir / "S1_ConfSearch",
            "S1",
            default_fields=progress_fields,
            event_callback=ui_reporter.event_callback,
        )
        s1_reporter.emit_stage_event(
            "s1_started",
            total_variants=len(variants),
            has_precursor=bool(precursor_smiles),
        )
        if precursor_smiles:
            precursor_s1_manifest = work_dir / "S1_ConfSearch" / "precursor" / "manifest.json"
            precursor_s1_payload = self._s1_signature_payload(
                precursor_smiles, "precursor", "precursor"
            )
            precursor_s1_signature = V4Checkpoint.signature(precursor_s1_payload)
            if self._checkpoint_reusable(
                checkpoint,
                "s1",
                precursor_s1_signature,
                precursor_s1_manifest,
                scope="precursor_s1",
                signature_payload=precursor_s1_payload,
            ):
                precursor_s1_data = self._read_json(precursor_s1_manifest)
                precursor_s1_result = {
                    "status": "complete",
                    "manifest": precursor_s1_manifest,
                    "selected_xyz": self._selected_s1_xyz(
                        work_dir,
                        precursor_s1_data,
                        "precursor",
                        precursor_s1_manifest,
                    ),
                    "data": precursor_s1_data,
                    "variant_id": "precursor",
                    "directory_name": "precursor",
                    "checkpoint_signature": self._effective_checkpoint_signature(
                        "s1", "precursor_s1", precursor_s1_signature
                    ),
                }
                s1_reporter.finish_structure(
                    "precursor",
                    "complete",
                    role="precursor",
                    smiles=precursor_smiles,
                    manifest=str(precursor_s1_manifest),
                    reused=True,
                    **self._s1_progress_fields(precursor_s1_data),
                )
                s1_reporter.science_summary(
                    variant="precursor", reused=True, **self._s1_progress_fields(precursor_s1_data)
                )
                logger.info("[V4] Reusing precursor S1 manifest")
            else:
                s1_reporter.start_structure("precursor", role="precursor", smiles=precursor_smiles)
                try:
                    precursor_engine = CensoLiteEngine(
                        self.config, work_dir / "S1_ConfSearch", "precursor"
                    )
                    if hasattr(precursor_engine, "set_event_callback"):
                        precursor_engine.set_event_callback(s1_reporter.external_event)
                    raw_precursor_s1 = precursor_engine.run(precursor_smiles)
                    manifest_path = Path(raw_precursor_s1["manifest"])
                    precursor_s1_data = dict(raw_precursor_s1.get("data") or self._read_json(manifest_path))
                    selected_xyz = Path(
                        raw_precursor_s1.get("selected_xyz")
                        or self._selected_s1_xyz(work_dir, precursor_s1_data, "precursor", manifest_path)
                    )
                    precursor_s1_result = {
                        "status": "complete",
                        "manifest": manifest_path,
                        "selected_xyz": selected_xyz,
                        "data": precursor_s1_data,
                        "variant_id": "precursor",
                        "directory_name": "precursor",
                        "checkpoint_signature": precursor_s1_signature,
                    }
                    checkpoint.mark_scope(
                        "precursor_s1",
                        "s1",
                        precursor_s1_signature,
                        manifest_path,
                        signature_payload=precursor_s1_payload,
                    )
                    s1_reporter.finish_structure(
                        "precursor",
                        "complete",
                        role="precursor",
                        smiles=precursor_smiles,
                        manifest=str(manifest_path),
                        **self._s1_progress_fields(precursor_s1_data),
                    )
                except Exception as exc:
                    checkpoint.mark_failed(
                        "precursor_s1",
                        "s1",
                        precursor_s1_signature,
                        str(exc),
                        signature_payload=precursor_s1_payload,
                    )
                    logger.error("[V4] Precursor S1 failed: %s", exc)
                    precursor_s1_result = {
                        "status": "failed",
                        "error": str(exc),
                        "variant_id": "precursor",
                        "directory_name": "precursor",
                        "checkpoint_signature": precursor_s1_signature,
                    }
                    s1_reporter.finish_structure(
                        "precursor",
                        "failed",
                        role="precursor",
                        smiles=precursor_smiles,
                        error=str(exc),
                    )

        shared_legacy_s1: Optional[Dict[str, Any]] = None
        for variant in variants:
            variant_id = str(variant["variant_id"])
            directory_name = str(variant["directory_name"])
            variant_smiles = str(variant.get("branch_product_smiles") or product_smiles)
            s1_manifest = work_dir / "S1_ConfSearch" / directory_name / "manifest.json"
            s1_payload = self._s1_signature_payload(
                variant_smiles, variant_id, directory_name
            )
            s1_signature = V4Checkpoint.signature(s1_payload)
            if shared_legacy_s1 is not None:
                shared_result = dict(shared_legacy_s1)
                shared_data = dict(shared_result.get("data") or {})
                shared_result.update(
                    status="complete",
                    variant_id=variant_id,
                    directory_name=directory_name,
                    checkpoint_signature=s1_signature,
                )
                checkpoint.mark_scope(
                    variant_id,
                    "s1",
                    s1_signature,
                    Path(shared_result["manifest"]),
                    signature_payload=s1_payload,
                )
                all_s1_results[variant_id] = shared_result
                s1_reporter.finish_structure(
                    variant_id,
                    "complete",
                    smiles=variant_smiles,
                    manifest=str(shared_result["manifest"]),
                    reused=True,
                    **self._s1_progress_fields(shared_data),
                )
                s1_reporter.science_summary(
                    variant=variant_id,
                    reused=True,
                    **self._s1_progress_fields(shared_data),
                )
                continue
            if self._checkpoint_reusable(
                checkpoint,
                "s1",
                s1_signature,
                s1_manifest,
                scope=variant_id,
                signature_payload=s1_payload,
            ):
                s1_data = self._read_json(s1_manifest)
                all_s1_results[variant_id] = {
                    "status": "complete",
                    "manifest": s1_manifest,
                    "selected_xyz": self._selected_s1_xyz(work_dir, s1_data, directory_name, s1_manifest),
                    "data": s1_data,
                    "variant_id": variant_id,
                    "directory_name": directory_name,
                    "checkpoint_signature": self._effective_checkpoint_signature(
                        "s1", variant_id, s1_signature
                    ),
                }
                s1_reporter.finish_structure(
                    variant_id,
                    "complete",
                    smiles=variant_smiles,
                    manifest=str(s1_manifest),
                    reused=True,
                    **self._s1_progress_fields(s1_data),
                )
                s1_reporter.science_summary(
                    variant=variant_id, reused=True, **self._s1_progress_fields(s1_data)
                )
                logger.info("[V4] Reusing S1 manifest for variant %s", variant_id)
                continue
            s1_reporter.start_structure(variant_id, smiles=variant_smiles)
            try:
                s1_engine = CensoLiteEngine(
                    self.config, work_dir / "S1_ConfSearch", directory_name
                )
                if hasattr(s1_engine, "set_event_callback"):
                    s1_engine.set_event_callback(s1_reporter.external_event)
                raw_s1 = s1_engine.run(variant_smiles)
                manifest_path = Path(raw_s1["manifest"])
                s1_data = dict(raw_s1.get("data") or self._read_json(manifest_path))
                selected_xyz = Path(
                    raw_s1.get("selected_xyz")
                    or self._selected_s1_xyz(work_dir, s1_data, directory_name, manifest_path)
                )
                s1_result = {
                    "status": "complete",
                    "manifest": manifest_path,
                    "selected_xyz": selected_xyz,
                    "data": s1_data,
                    "variant_id": variant_id,
                    "directory_name": directory_name,
                    "checkpoint_signature": s1_signature,
                }
                checkpoint.mark_scope(
                    variant_id,
                    "s1",
                    s1_signature,
                    manifest_path,
                    signature_payload=s1_payload,
                )
                all_s1_results[variant_id] = s1_result
                s1_reporter.finish_structure(
                    variant_id,
                    "complete",
                    smiles=variant_smiles,
                    manifest=str(manifest_path),
                    **self._s1_progress_fields(s1_data),
                )
                if manifest_path.parent.resolve() != s1_manifest.parent.resolve():
                    shared_legacy_s1 = dict(s1_result)
                    logger.warning(
                        "[V4] S1 for variant %s returned legacy path %s; remaining variants will reuse the shared output",
                        variant_id,
                        manifest_path,
                    )
            except Exception as exc:
                checkpoint.mark_failed(
                    variant_id,
                    "s1",
                    s1_signature,
                    str(exc),
                    signature_payload=s1_payload,
                )
                logger.error("[V4] Variant %s S1 failed: %s", variant_id, exc)
                all_s1_results[variant_id] = {
                    "status": "failed",
                    "error": str(exc),
                    "variant_id": variant_id,
                    "directory_name": directory_name,
                    "checkpoint_signature": s1_signature,
                }
                s1_reporter.finish_structure(
                    variant_id,
                    "failed",
                    smiles=variant_smiles,
                    error=str(exc),
                )

        s1_stage_manifest = self._representative_stage_manifest(variants, all_s1_results)
        s1_stage_signature = self._variant_stage_signature("s1", variants, all_s1_results)
        s1_reporter.emit_stage_event(
            "s1_completed",
            total_structures=s1_reporter.structure_count,
            manifest=str(s1_stage_manifest) if s1_stage_manifest is not None else None,
        )
        if (
            s1_stage_manifest is not None
            and not checkpoint.reusable(
                "s1", s1_stage_signature, s1_stage_manifest
            )
        ):
            checkpoint.mark("s1", s1_stage_signature, s1_stage_manifest)
        if stop_after == "s1":
            manifests = {"s0": s0_manifest}
            if s1_stage_manifest is not None:
                manifests["s1"] = s1_stage_manifest
            self._write_run_manifest(
                work_dir,
                s0_record,
                variants,
                all_s1_results,
                all_s2_results,
                "s1",
                s0_manifest,
            )
            return self._write_pipeline_result(
                work_dir,
                stop_after,
                manifests,
            )

        s2_reporter = StageProgressReporter(
            work_dir / "S2_PEB",
            "S2",
            default_fields=progress_fields,
            event_callback=ui_reporter.event_callback,
        )
        s2_reporter.emit_stage_event("s2_started", total_variants=len(variants))
        for variant in variants:
            variant_id = str(variant["variant_id"])
            s1_result = all_s1_results.get(variant_id, {})
            if s1_result.get("status") != "complete":
                all_s2_results[variant_id] = {
                    "status": "skipped",
                    "reason": "s1_failed",
                }
                s2_reporter.finish_structure(variant_id, "skipped", reason="s1_failed")
                continue
            product_xyz = Path(s1_result["selected_xyz"])
            s2_dir = work_dir / "S2_PEB" / variant_id
            s2_manifest = s2_dir / "manifest.json"
            s2_payload = {
                "signature_schema": "s2_signature_v3",
                "stage": "s2",
                "variant_id": variant_id,
                "s0": s0_signature,
                "selected_id": (s1_result.get("data") or {}).get("selected"),
                "selected_xyz_ref": str(work_dir / "S1_ConfSearch" / variant_id / "selected.xyz"),
                "product_xyz": V4Checkpoint.file_signature(product_xyz),
                "forming_bonds": forming_bonds,
                "config": self.config.get("step2", {}),
            }
            s2_signature = V4Checkpoint.signature(s2_payload)
            scan_batch = f"{variant_id}:peb_scan"
            try:
                if self._checkpoint_reusable(
                    checkpoint,
                    "s2",
                    s2_signature,
                    s2_manifest,
                    scope=variant_id,
                    signature_payload=s2_payload,
                ):
                    s2_data = self._read_json(s2_manifest)
                    normalized_stage_status = self._normalize_progress_status(s2_data.get("status", "complete"))
                    all_s2_results[variant_id] = {
                        "status": "complete",
                        "manifest": s2_manifest,
                        "ts_seed": Path(s2_data["peb_peak"]),
                        "intermediate_seed": Path(s2_data["intermediate_seed"]),
                        "variant_id": s2_data.get("variant_id"),
                        "parent_variant_id": s2_data.get("parent_variant_id"),
                        "selected_id": s2_data.get("selected_id"),
                        "selected_xyz_ref": s2_data.get("selected_xyz_ref"),
                        "forming_bonds": tuple(tuple(pair) for pair in s2_data.get("forming_bonds", [])),
                        "scan_profile": Path(s2_data["scan_profile"]),
                        "stage_status": s2_data.get("status", "reused"),
                        "confidence": s2_data.get("confidence", "unknown"),
                        "degraded_reasons": tuple(s2_data.get("degraded_reasons", [])),
                        "checkpoint_signature": s2_signature,
                    }
                    s2_reporter.finish_structure(
                        variant_id,
                        normalized_stage_status,
                        manifest=str(s2_manifest),
                        stage_status=s2_data.get("status", "reused"),
                        confidence=s2_data.get("confidence", "unknown"),
                        degraded_reasons=list(s2_data.get("degraded_reasons", [])),
                        reused=True,
                    )
                    s2_reporter.science_summary(
                        variant=variant_id,
                        forming_bonds=s2_data.get("forming_bonds", []),
                        confidence=s2_data.get("confidence", "unknown"),
                        degraded_reasons=list(s2_data.get("degraded_reasons", [])),
                        reused=True,
                    )
                    logger.info("[V4] Reusing S2 manifest for variant %s", variant_id)
                    continue
                s2_reporter.start_structure(
                    variant_id,
                    selected_id=(s1_result.get("data") or {}).get("selected"),
                    product_xyz=str(product_xyz),
                )
                scanner = PEBScanner(self.config, molecule_name=variant_id)
                scan_steps = int(
                    ((self.config.get("step2", {}) or {}).get("scan", {}) or {}).get(
                        "scan_steps", 24
                    )
                )
                s2_reporter.external_event(
                    "batch_started",
                    {
                        "batch": scan_batch,
                        "label": "xTB PEB scan",
                        "total": scan_steps,
                        "running": 1,
                        "forming_bonds": [list(pair) for pair in forming_bonds],
                    },
                )
                result = scanner.run(product_xyz, s2_dir, forming_bonds)
                ts_seed, _legacy, intermediate_seed, returned_bonds, profile, stage_status, confidence, degraded = result[:8]
                s2_manifest = self._write_s2(
                    s2_dir,
                    ts_seed,
                    intermediate_seed,
                    returned_bonds,
                    profile,
                    stage_status,
                    confidence,
                    degraded,
                    variant_id=variant_id,
                    selected_id=(s1_result.get("data") or {}).get("selected"),
                    selected_xyz_ref=work_dir / "S1_ConfSearch" / variant_id / "selected.xyz",
                )
                checkpoint.mark_scope(
                    variant_id,
                    "s2",
                    s2_signature,
                    s2_manifest,
                    signature_payload=s2_payload,
                )
                s2_reporter.external_event(
                    "batch_finished",
                    {
                        "batch": scan_batch,
                        "label": "xTB PEB scan",
                        "total": scan_steps,
                        "done": scan_steps,
                        "running": 0,
                        "failed": 0,
                        "status": "complete",
                    },
                )
                s2_reporter.science_summary(
                    variant=variant_id,
                    forming_bonds=[list(pair) for pair in returned_bonds],
                    confidence=confidence,
                    degraded_reasons=list(degraded),
                    ts_seed=str(ts_seed),
                    scan_profile=str(profile),
                )
                all_s2_results[variant_id] = {
                    "status": "complete",
                    "manifest": s2_manifest,
                    "ts_seed": Path(ts_seed),
                    "intermediate_seed": Path(intermediate_seed),
                    "variant_id": variant_id,
                    "parent_variant_id": variant_id,
                    "selected_id": (s1_result.get("data") or {}).get("selected"),
                    "selected_xyz_ref": str(work_dir / "S1_ConfSearch" / variant_id / "selected.xyz"),
                    "forming_bonds": tuple(tuple(pair) for pair in returned_bonds),
                    "scan_profile": Path(profile),
                    "stage_status": stage_status,
                    "confidence": confidence,
                    "degraded_reasons": tuple(degraded),
                    "checkpoint_signature": s2_signature,
                }
                s2_reporter.finish_structure(
                    variant_id,
                    self._normalize_progress_status(stage_status),
                    manifest=str(s2_manifest),
                    stage_status=stage_status,
                    confidence=confidence,
                    degraded_reasons=list(degraded),
                    ts_seed=str(ts_seed),
                    intermediate_seed=str(intermediate_seed),
                    forming_bonds=[list(pair) for pair in returned_bonds],
                    scan_profile=str(profile),
                )
            except ResumeCheckpointError:
                raise
            except Exception as exc:
                s2_reporter.external_event(
                    "batch_finished",
                    {
                        "batch": scan_batch,
                        "label": "xTB PEB scan",
                        "running": 0,
                        "failed": 1,
                        "status": "failed",
                    },
                )
                checkpoint.mark_failed(
                    variant_id,
                    "s2",
                    s2_signature,
                    str(exc),
                    signature_payload=s2_payload,
                )
                logger.error("[V4] Variant %s S2 failed: %s", variant_id, exc)
                all_s2_results[variant_id] = {
                    "status": "failed",
                    "error": str(exc),
                    "checkpoint_signature": s2_signature,
                }
                s2_reporter.finish_structure(variant_id, "failed", error=str(exc))

        s2_stage_manifest = self._representative_stage_manifest(variants, all_s2_results)
        s2_stage_signature = self._variant_stage_signature("s2", variants, all_s2_results)
        s2_reporter.emit_stage_event(
            "s2_completed",
            total_structures=s2_reporter.structure_count,
            manifest=str(s2_stage_manifest) if s2_stage_manifest is not None else None,
        )
        if (
            s2_stage_manifest is not None
            and not checkpoint.reusable(
                "s2", s2_stage_signature, s2_stage_manifest
            )
        ):
            checkpoint.mark("s2", s2_stage_signature, s2_stage_manifest)
        if stop_after == "s2":
            manifests = {"s0": s0_manifest}
            if s1_stage_manifest is not None:
                manifests["s1"] = s1_stage_manifest
            if s2_stage_manifest is not None:
                manifests["s2"] = s2_stage_manifest
            self._write_run_manifest(
                work_dir,
                s0_record,
                variants,
                all_s1_results,
                all_s2_results,
                "s2",
                s0_manifest,
            )
            return self._write_pipeline_result(
                work_dir,
                stop_after,
                manifests,
            )

        structures: List[Dict[str, Any]] = []
        if precursor_s1_result and precursor_s1_result.get("status") == "complete":
            structures.append({
                "id": "precursor",
                "kind": "minimum",
                "input_xyz": str(Path(precursor_s1_result["selected_xyz"])),
                "structure_id": "precursor",
                "variant_id": "precursor",
                "role": "precursor",
                "branch_id": None,
                "pathway_id": None,
                "parent_structure_id": None,
                "source_stage": "S1",
                **self._s1_ensemble_fields(precursor_s1_result),
            })
        for variant in variants:
            variant_id = str(variant["variant_id"])
            branch_id = variant.get("branch_id")
            pathway_id = variant.get("pathway_id")
            s1_result = all_s1_results.get(variant_id, {})
            s2_result = all_s2_results.get(variant_id, {})
            if s1_result.get("status") != "complete" or s2_result.get("status") != "complete":
                continue
            structures.extend([
                {
                    "id": variant_id,
                    "kind": "minimum",
                    "input_xyz": str(Path(s1_result["selected_xyz"])),
                    "structure_id": variant_id,
                    "variant_id": variant_id,
                    "role": "product",
                    "branch_id": branch_id,
                    "pathway_id": pathway_id,
                    "parent_structure_id": None,
                    "source_stage": "S1",
                    **self._s1_ensemble_fields(s1_result),
                },
                {
                    "id": f"{variant_id}_int",
                    "kind": "minimum",
                    "input_xyz": str(Path(s2_result["intermediate_seed"])),
                    "structure_id": f"{variant_id}_int",
                    "variant_id": variant_id,
                    "role": "intermediate",
                    "branch_id": branch_id,
                    "pathway_id": pathway_id,
                    "parent_structure_id": variant_id,
                    "source_stage": "S2",
                    "forming_bonds": [list(b) for b in forming_bonds],
                },
                {
                    "id": f"{variant_id}_ts",
                    "kind": "ts",
                    "input_xyz": str(Path(s2_result["ts_seed"])),
                    "structure_id": f"{variant_id}_ts",
                    "variant_id": variant_id,
                    "role": "ts",
                    "branch_id": branch_id,
                    "pathway_id": pathway_id,
                    "parent_structure_id": variant_id,
                    "source_stage": "S2",
                    "forming_bonds": [list(b) for b in forming_bonds],
                },
            ])
        s3_dir = work_dir / "S3_LowLevel"
        s3_manifest = s3_dir / "manifest.json"
        s3_reporter = StageProgressReporter(
            s3_dir,
            "S3",
            default_fields=progress_fields,
            event_callback=ui_reporter.event_callback,
        )
        s3_reporter.emit_stage_event("s3_started", total_structures=len(structures))
        s3_payload = {
            "signature_schema": "s3_signature_v2",
            "stage": "s3",
            "structures": self._structure_signatures(structures),
            "config": self.config.get("theory", {}).get("s3_low_level", {}),
        }
        s3_signature = V4Checkpoint.signature(s3_payload)
        s3_reused = self._checkpoint_reusable(
            checkpoint,
            "s3",
            s3_signature,
            s3_manifest,
            signature_payload=s3_payload,
        )
        if s3_reused:
            logger.info("[V4] Reusing S3 manifest")
        else:
            s3_engine = LowLevelEngine(self.config)
            if hasattr(s3_engine, "set_progress_reporter"):
                s3_engine.set_progress_reporter(s3_reporter)
            s3_manifest = s3_engine.run(structures, s3_dir)
            checkpoint.mark(
                "s3",
                s3_signature,
                s3_manifest,
                signature_payload=s3_payload,
            )
        if s3_reporter.structure_count == 0:
            self._sync_stage_reporter_from_stage_manifest(s3_reporter, s3_manifest)
        s3_reporter.emit_stage_event(
            "s3_completed",
            total_structures=s3_reporter.structure_count,
            manifest=str(s3_manifest),
            reused=s3_reused,
        )
        if stop_after == "s3":
            manifests = {"s0": s0_manifest, "s3": s3_manifest}
            if s1_stage_manifest is not None:
                manifests["s1"] = s1_stage_manifest
            if s2_stage_manifest is not None:
                manifests["s2"] = s2_stage_manifest
            self._write_run_manifest(
                work_dir,
                s0_record,
                variants,
                all_s1_results,
                all_s2_results,
                "s3",
                s0_manifest,
                s3_manifest=s3_manifest,
            )
            return self._write_pipeline_result(
                work_dir,
                stop_after,
                manifests,
            )
        s3_data = json.loads(s3_manifest.read_text(encoding="utf-8"))
        s4_structures = [
            {
                "id": item["id"],
                "kind": item.get("kind", "minimum"),
                "input_xyz": item["input_xyz"],
                "opt_xyz": item.get("opt_xyz"),
                "fallback_xyz": item.get("input_xyz"),
                "charge": item.get("charge"),
                "multiplicity": item.get("multiplicity"),
                "structure_id": item.get("structure_id", item["id"]),
                "variant_id": item.get("variant_id"),
                "role": item.get("role"),
                "branch_id": item.get("branch_id"),
                "pathway_id": item.get("pathway_id"),
                "parent_structure_id": item.get("parent_structure_id"),
                "source_stage": "S3",
                "s1_manifest": item.get("s1_manifest"),
                "s1_ensemble_thermodynamics": item.get("s1_ensemble_thermodynamics"),
                "s1_thermochemistry_status": item.get("s1_thermochemistry_status"),
                "ensemble_thermochemistry_correction_hartree": item.get(
                    "ensemble_thermochemistry_correction_hartree"
                ),
                "source_s3": {
                    "structure_status": item.get("status"),
                    "opt_status": item.get("opt_status"),
                    "sp_status": item.get("sp_status"),
                    "frequency_status": item.get("frequency_status"),
                    "ts_frequency_valid": item.get("ts_frequency_valid"),
                    "ts_mode_displacement_verified": item.get("ts_mode_displacement_verified"),
                    "frequency_count_valid": item.get("frequency_count_valid"),
                    "mode_displacement_valid": item.get("mode_displacement_valid"),
                    "irc_valid": item.get("irc_valid"),
                    "ts_quality_summary": item.get("ts_quality_summary"),
                    "usable_for_ml": item.get("usable_for_ml"),
                    "manifest": str(s3_manifest),
                },
                **(
                    {"forming_bonds": item.get("forming_bonds")}
                    if item.get("kind", "minimum") == "ts" and item.get("forming_bonds") is not None
                    else {}
                ),
            }
            for item in s3_data.get("structures", [])
        ]
        s4_dir = work_dir / "S4_HighLevel"
        s4_manifest = s4_dir / "manifest.json"
        s4_payload = {
            "signature_schema": "s4_signature_v3",
            "stage": "s4",
            "structures": self._structure_signatures(s4_structures),
            "config": self.config.get("theory", {}).get("s4_high_precision", {}),
        }
        s4_signature = V4Checkpoint.signature(s4_payload)
        if self._checkpoint_reusable(
            checkpoint,
            "s4",
            s4_signature,
            s4_manifest,
            signature_payload=s4_payload,
        ):
            logger.info("[V4] Reusing S4 manifest")
        else:
            s4_manifest = HighLevelEngine(self.config).run(
                s4_structures,
                s4_dir,
                event_callback=ui_reporter.s4_event_callback,
            )
            s4_data = self._read_json(s4_manifest)
            s4_checkpoint_status = (
                "complete" if s4_data.get("status") == "complete" else "incomplete"
            )
            checkpoint.mark(
                "s4",
                s4_signature,
                s4_manifest,
                status=s4_checkpoint_status,
                signature_payload=s4_payload,
            )
        manifests = {"s0": s0_manifest, "s3": s3_manifest, "s4": s4_manifest}
        if s1_stage_manifest is not None:
            manifests["s1"] = s1_stage_manifest
        if s2_stage_manifest is not None:
            manifests["s2"] = s2_stage_manifest
        self._write_run_manifest(
            work_dir,
            s0_record,
            variants,
            all_s1_results,
            all_s2_results,
            "s4",
            s0_manifest,
            s3_manifest=s3_manifest,
            s4_manifest=s4_manifest,
        )
        return self._write_pipeline_result(
            work_dir,
            stop_after,
            manifests,
        )

    @staticmethod
    def _s1_ensemble_fields(s1_result: Dict[str, Any]) -> Dict[str, Any]:
        """Carry the S1 conformational correction into downstream minima."""

        data = dict(s1_result.get("data") or {})
        thermo = dict(data.get("ensemble_thermodynamics") or {})
        if not thermo:
            return {}
        return {
            "s1_manifest": str(s1_result.get("manifest")) if s1_result.get("manifest") else None,
            "s1_ensemble_thermodynamics": thermo,
            "s1_thermochemistry_status": thermo.get(
                "thermochemistry_status", thermo.get("status")
            ),
            "ensemble_thermochemistry_correction_hartree": thermo.get(
                "ensemble_thermochemistry_correction_hartree"
            ),
        }

    @staticmethod
    def _s1_progress_fields(data: Dict[str, Any]) -> Dict[str, Any]:
        candidates = list(data.get("candidates") or [])
        representatives = list(
            data.get("reactivity_screening_candidates")
            or data.get("representative_candidates")
            or []
        )
        thermo = dict(data.get("ensemble_thermodynamics") or {})
        pipeline_log = list(data.get("pipeline_log") or [])
        funnel = {
            str(row.get("step")): row.get("candidates")
            for row in pipeline_log
            if row.get("step") and row.get("candidates") is not None
        }
        top_population = [
            {
                "id": row.get("id"),
                "population": row.get("boltzmann_population"),
                "relative_free_energy_kcal": row.get("relative_free_energy_kcal"),
            }
            for row in sorted(
                candidates,
                key=lambda item: float(item.get("boltzmann_population") or 0.0),
                reverse=True,
            )[:5]
        ]
        return {
            "funnel": funnel,
            "ensemble_members": len(candidates),
            "representatives": len(representatives),
            "ranking_formula": data.get("ranking_formula"),
            "scoring_mode": data.get("scoring_mode"),
            "selected": data.get("selected"),
            "thermodynamic_rank1": data.get("thermodynamic_rank1"),
            "top_population": top_population,
            "conformational_free_energy_correction_kcal": thermo.get(
                "conformational_free_energy_correction_kcal"
            ),
            "thermochemistry_complete": thermo.get("thermochemistry_complete"),
            "thermochemistry_status": thermo.get(
                "thermochemistry_status", thermo.get("status")
            ),
        }

    @staticmethod
    def _write_pipeline_result(
        work_dir: Path,
        completed_through: str,
        manifests: Dict[str, Path],
    ) -> Dict[str, Any]:
        result: Dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            "success": True,
            "completed_through": completed_through,
            "work_dir": str(work_dir),
        }
        for stage in ("s0", "s1", "s2", "s3", "s4"):
            manifest = manifests.get(stage)
            if manifest is not None:
                result[f"{stage}_manifest"] = str(Path(manifest))
        (Path(work_dir) / "pipeline.result.json").write_text(
            json.dumps(result, indent=2),
            encoding="utf-8",
        )
        return result

    @staticmethod
    def _read_json(path: Path) -> Dict[str, Any]:
        return json.loads(Path(path).read_text(encoding="utf-8"))

    @staticmethod
    def _selected_s1_xyz(
        work_dir: Path,
        manifest: Dict[str, Any],
        variant_id: str = "product",
        manifest_path: Optional[Path] = None,
    ) -> Path:
        selected_xyz = manifest.get("selected_xyz")
        if isinstance(selected_xyz, str) and selected_xyz.strip():
            if manifest_path is not None:
                return Path(manifest_path).parent / selected_xyz
            return work_dir / "S1_ConfSearch" / variant_id / selected_xyz
        selected = manifest.get("selected")
        for candidate in manifest.get("candidates", []):
            if candidate.get("id") == selected:
                return work_dir / "S1_ConfSearch" / str(candidate["xyz"])
        raise RuntimeError("S1 manifest has no selected candidate")

    @staticmethod
    def _load_variants(work_dir: Path, product_smiles: str) -> List[Dict[str, Any]]:
        registry_path = work_dir / "S0_Mechanism" / "variant_registry.json"
        if not registry_path.exists():
            return [{
                "variant_id": "product_major",
                "directory_name": "product_major",
                "branch_product_smiles": product_smiles,
            }]
        payload = json.loads(registry_path.read_text(encoding="utf-8"))
        variants = payload.get("variants", []) or []
        normalized: List[Dict[str, Any]] = []
        for index, variant in enumerate(variants, start=1):
            variant_id = str(variant.get("variant_id") or f"product_{index:03d}")
            normalized.append({
                **variant,
                "variant_id": variant_id,
                "directory_name": str(variant.get("directory_name") or variant_id),
                "branch_product_smiles": variant.get("branch_product_smiles") or product_smiles,
            })
        if normalized:
            return normalized
        return [{
                "variant_id": "product_major",
                "directory_name": "product_major",
                "branch_product_smiles": product_smiles,
            }]

    def _s1_signature_payload(
        self,
        smiles: str,
        variant_id: str,
        directory_name: str,
    ) -> Dict[str, Any]:
        """Build a stable, explainable signature for scientific S1 identity."""

        executables = self.config.get("executables", {}) or {}

        def executable_identity(name: str) -> Optional[Dict[str, Any]]:
            value = executables.get(name, {}) or {}
            path = value.get("path") if isinstance(value, Mapping) else None
            if not path:
                return None
            configured = Path(str(path))
            identity: Dict[str, Any] = {"path": str(path)}
            try:
                stat = configured.stat()
            except OSError:
                return identity
            identity.update(size=stat.st_size, mtime_ns=stat.st_mtime_ns)
            return identity

        return {
            "signature_schema": "s1_signature_v3",
            "stage": "s1",
            "protocol_version": "censo_light_ranking_v3",
            "manifest_schema": S1_MANIFEST_SCHEMA_VERSION,
            "config": self._scientific_s1_config(
                self.config.get("step1", {}) or {}
            ),
            "executables": {
                "crest": executable_identity("crest"),
                "xtb": executable_identity("xtb"),
                "orca": executable_identity("orca"),
            },
            "smiles": smiles,
            "variant_id": variant_id,
            "directory_name": directory_name,
        }

    @staticmethod
    def _scientific_s1_config(config: Mapping[str, Any]) -> Dict[str, Any]:
        """Remove scheduling-only settings from the S1 scientific signature."""

        operational_keys = {
            "parallel_jobs",
            "cores_per_job",
            "max_parallel_jobs",
            "nproc",
            "retry_nproc",
            "fallback_nproc",
            "timeout",
            "omp_stacksize",
            "omp_max_active_levels",
            "adaptive_tail",
            "queue_drain",
            "scheduling",
            "orca_math_threads_per_rank",
        }

        def normalize(value: Any) -> Any:
            if isinstance(value, Mapping):
                return {
                    str(key): normalize(item)
                    for key, item in sorted(value.items(), key=lambda row: str(row[0]))
                    if str(key) not in operational_keys
                }
            if isinstance(value, (list, tuple)):
                return [normalize(item) for item in value]
            return value

        return dict(normalize(config))

    def _checkpoint_reusable(
        self,
        checkpoint: V4Checkpoint,
        stage: str,
        signature: str,
        manifest: Path,
        *,
        scope: Optional[str] = None,
        signature_payload: Optional[Mapping[str, Any]] = None,
    ) -> bool:
        """Apply the configured resume policy to one checkpoint decision."""

        decision = checkpoint.check_reuse(
            stage,
            signature,
            manifest,
            scope=scope,
            signature_payload=signature_payload,
        )
        if decision.reusable:
            return True

        start_from = getattr(self, "start_from", None)
        is_required_upstream = bool(
            start_from
            and V4Checkpoint.stage_order.index(stage)
            < V4Checkpoint.stage_order.index(start_from)
        )
        policy = str(getattr(self, "resume_policy", "strict") or "strict")

        if (
            policy == "use-existing-upstream"
            and is_required_upstream
            and stage in ("s0", "s1", "s2")
            and decision.reason == "signature_mismatch"
        ):
            accepted = checkpoint.check_reuse(
                stage,
                signature,
                manifest,
                scope=scope,
                signature_payload=signature_payload,
                accept_signature_mismatch=True,
            )
            if accepted.reusable:
                if accepted.recorded_signature:
                    overrides = getattr(
                        self, "_reused_signature_overrides", None
                    )
                    if overrides is None:
                        overrides = {}
                        self._reused_signature_overrides = overrides
                    overrides[(stage, scope)] = accepted.recorded_signature
                logger.warning(
                    "[V4] Explicitly reusing completed upstream %s:%s despite "
                    "signature mismatch; recorded=%s current=%s",
                    stage,
                    scope or "aggregate",
                    accepted.recorded_signature,
                    accepted.requested_signature,
                )
                return True

        if is_required_upstream:
            raise ResumeCheckpointError(
                f"Cannot start from {start_from}: required upstream "
                f"{decision.describe()}"
            )

        if decision.reason in {
            "signature_mismatch",
            "manifest_mismatch",
            "manifest_missing",
            "manifest_unrecorded",
        }:
            # S0 is a cheap record-driven reconstruction. A changed trusted
            # reaction record or condition set must regenerate S0 and
            # invalidate downstream state rather than blocking an ordinary
            # rerun under the strict policy.
            if stage == "s0":
                logger.warning("[V4] Recomputing S0 after %s", decision.describe())
                return False
            if policy == "strict":
                raise ResumeCheckpointError(
                    f"{decision.describe()}. Use --recompute-from {stage} to "
                    "recalculate, or explicitly start from a downstream stage "
                    "with --resume-policy use-existing-upstream."
                )
            logger.warning("[V4] Recomputing after %s", decision.describe())
        return False

    def _effective_checkpoint_signature(
        self,
        stage: str,
        scope: Optional[str],
        requested_signature: str,
    ) -> str:
        overrides = getattr(self, "_reused_signature_overrides", {}) or {}
        return str(overrides.get((stage, scope), requested_signature))

    @staticmethod
    def _representative_stage_manifest(
        variants: Sequence[Dict[str, Any]],
        results: Dict[str, Dict[str, Any]],
    ) -> Optional[Path]:
        for variant in variants:
            variant_id = str(variant["variant_id"])
            result = results.get(variant_id, {})
            if result.get("status") == "complete" and result.get("manifest") is not None:
                return Path(result["manifest"])
        return None

    @staticmethod
    def _variant_stage_signature(
        stage: str,
        variants: Sequence[Dict[str, Any]],
        results: Dict[str, Dict[str, Any]],
    ) -> str:
        payload = {
            "stage": stage,
            "variants": [],
        }
        for variant in variants:
            variant_id = str(variant["variant_id"])
            result = results.get(variant_id, {})
            payload["variants"].append({
                "variant_id": variant_id,
                "status": result.get("status", "pending"),
                "checkpoint_signature": result.get("checkpoint_signature"),
                "manifest": str(result["manifest"]) if result.get("manifest") is not None else None,
                "error": result.get("error"),
                "reason": result.get("reason"),
            })
        return V4Checkpoint.signature(payload)

    @staticmethod
    def _condition_metadata(record: S0ReactionRecord) -> Tuple[Dict[str, Any], str]:
        raw_row = getattr(record, "raw_row", {}) or {}
        condition_fields: Dict[str, Any] = {
            "solvent": raw_row.get("solvent", ""),
            "temperature_celsius": raw_row.get("temp_celsius", ""),
            "has_lewis_acid": raw_row.get("has_lewis_acid", "false"),
            "charge": 0,
            "multiplicity": 1,
        }
        additive_info = {
            "has_lewis_acid": raw_row.get("has_lewis_acid", "false"),
            "catalyst": raw_row.get("catalyst", ""),
            "additive": raw_row.get("additive", ""),
        }
        if any(str(value).strip() for key, value in additive_info.items() if key != "has_lewis_acid"):
            condition_fields["additive_info"] = additive_info
        condition_signature = hashlib.sha256(
            json.dumps(condition_fields, sort_keys=True, default=str).encode("utf-8")
        ).hexdigest()[:16]
        return condition_fields, condition_signature

    @staticmethod
    def _normalize_progress_status(status: Any) -> str:
        text = str(status or "unknown").strip().lower()
        if text in {"complete", "completed", "success", "ok"}:
            return "complete"
        if text in {"degraded", "opt_failed_sp_complete", "ts_frequency_unverified"}:
            return text
        if text in {"failed", "error"}:
            return "failed"
        if text == "skipped":
            return "skipped"
        return text or "unknown"

    @staticmethod
    def _sync_stage_reporter_from_stage_manifest(
        reporter: StageProgressReporter,
        manifest_path: Path,
    ) -> None:
        payload = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
        for row in payload.get("structures", []):
            structure_id = str(row.get("id") or row.get("structure_id") or "").strip()
            if not structure_id:
                continue
            structure_payload = dict(row)
            structure_payload.pop("id", None)
            structure_payload.pop("structure_id", None)
            structure_payload.pop("status", None)
            reporter.finish_structure(
                structure_id,
                V4Orchestrator._normalize_progress_status(row.get("status", "complete")),
                **structure_payload,
            )

    def _write_run_manifest(
        self,
        work_dir: Path,
        record: S0ReactionRecord,
        variants: Sequence[Dict[str, Any]],
        all_s1_results: Dict[str, Dict[str, Any]],
        all_s2_results: Dict[str, Dict[str, Any]],
        completed_through: str,
        s0_manifest: Path,
        s3_manifest: Optional[Path] = None,
        s4_manifest: Optional[Path] = None,
    ) -> Path:
        condition_fields, condition_signature = self._condition_metadata(record)
        manifest: Dict[str, Any] = {
            "schema_version": "rph_run_manifest_v1",
            "reaction_id": record.rx_id,
            "completed_through": completed_through,
            "condition_signature": condition_signature,
            "conditions": condition_fields,
            "variants": [str(variant["variant_id"]) for variant in variants],
            "variant_status": {},
            "summary": {
                "total_variants": len(variants),
                "s1_complete": 0,
                "s2_complete": 0,
                "failed": 0,
            },
            "stages": {
                "s0": {
                    "status": "complete",
                    "manifest": str(Path(s0_manifest)),
                },
                "s1": {},
                "s2": {},
            },
        }
        for variant in variants:
            variant_id = str(variant["variant_id"])
            s1_result = all_s1_results.get(variant_id, {})
            s2_result = all_s2_results.get(variant_id, {})
            if s1_result:
                row = {"status": s1_result.get("status", "pending")}
                if s1_result.get("manifest") is not None:
                    row["manifest"] = str(Path(s1_result["manifest"]))
                if s1_result.get("error"):
                    row["error"] = s1_result["error"]
                manifest["stages"]["s1"][variant_id] = row
            if s2_result:
                row = {"status": s2_result.get("status", "pending")}
                if s2_result.get("manifest") is not None:
                    row["manifest"] = str(Path(s2_result["manifest"]))
                if s2_result.get("stage_status") is not None:
                    row["stage_status"] = s2_result["stage_status"]
                if s2_result.get("reason"):
                    row["reason"] = s2_result["reason"]
                if s2_result.get("error"):
                    row["error"] = s2_result["error"]
                manifest["stages"]["s2"][variant_id] = row

            if s1_result.get("status") == "complete":
                manifest["summary"]["s1_complete"] += 1
            if s2_result.get("status") == "complete":
                manifest["summary"]["s2_complete"] += 1

            if s1_result.get("status") == "failed" or s2_result.get("status") == "failed":
                manifest["variant_status"][variant_id] = "failed"
                manifest["summary"]["failed"] += 1
            elif s2_result.get("status") == "complete":
                manifest["variant_status"][variant_id] = "s2_complete"
            elif s1_result.get("status") == "complete":
                manifest["variant_status"][variant_id] = "s1_complete"
            else:
                manifest["variant_status"][variant_id] = "pending"

        if s3_manifest is not None:
            manifest["stages"]["s3"] = {
                "status": "complete",
                "manifest": str(Path(s3_manifest)),
            }
        if s4_manifest is not None:
            manifest["stages"]["s4"] = {
                "status": "complete",
                "manifest": str(Path(s4_manifest)),
            }
        manifest["overall_status"] = "partial_failure" if manifest["summary"]["failed"] else "complete"
        path = Path(work_dir) / "run.manifest.json"
        path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        return path

    @staticmethod
    def _structure_signatures(structures: Sequence[Dict[str, Any]]) -> Sequence[Dict[str, Any]]:
        payload = []
        for structure in structures:
            input_value = structure.get("input_xyz") or structure.get("fallback_xyz")
            if not input_value:
                raise ValueError(f"Structure {structure.get('id')} has no input XYZ")
            path = Path(input_value)
            record = {
                "id": structure.get("id"),
                "kind": structure.get("kind", "minimum"),
                "xyz": str(path.resolve()),
                "xyz_hash": V4Checkpoint.file_signature(path),
                "structure_id": structure.get("structure_id"),
                "variant_id": structure.get("variant_id"),
                "role": structure.get("role"),
                "branch_id": structure.get("branch_id"),
                "pathway_id": structure.get("pathway_id"),
                "parent_structure_id": structure.get("parent_structure_id"),
                "source_stage": structure.get("source_stage"),
                "ensemble_thermochemistry_correction_hartree": structure.get(
                    "ensemble_thermochemistry_correction_hartree"
                ),
                "s1_ensemble_thermodynamics": structure.get(
                    "s1_ensemble_thermodynamics"
                ),
                "s1_thermochemistry_status": structure.get(
                    "s1_thermochemistry_status"
                ),
            }
            s1_manifest_value = structure.get("s1_manifest")
            if s1_manifest_value:
                s1_manifest_path = Path(s1_manifest_value)
                record["s1_manifest"] = str(s1_manifest_path.resolve())
                record["s1_manifest_hash"] = (
                    V4Checkpoint.file_signature(s1_manifest_path)
                    if s1_manifest_path.exists()
                    else None
                )
            opt_value = structure.get("opt_xyz")
            if opt_value:
                opt_path = Path(opt_value)
                record["opt_xyz"] = str(opt_path.resolve())
                record["opt_xyz_hash"] = V4Checkpoint.file_signature(opt_path) if opt_path.exists() else None
            if structure.get("forming_bonds") is not None:
                record["forming_bonds"] = [
                    [int(pair[0]), int(pair[1])]
                    for pair in structure.get("forming_bonds", [])
                ]
            payload.append(record)
        return payload

    @staticmethod
    def _write_schema(work_dir: Path) -> None:
        (work_dir / ".rph_schema.json").write_text(json.dumps({"schema_version": SCHEMA_VERSION}, indent=2), encoding="utf-8")

    @staticmethod
    def _write_s0_from_record(work_dir: Path, record: S0ReactionRecord) -> Path:
        _reaction_context, _product_variants, mechanism_path = write_s0_artifacts(work_dir, record)
        return mechanism_path

    @staticmethod
    def _write_s2(
        directory: Path,
        ts_seed: Path,
        intermediate_seed: Path,
        forming_bonds: Sequence[Tuple[int, int]],
        profile: Path,
        status: str,
        confidence: str,
        degraded: Sequence[str],
        variant_id: Optional[str] = None,
        selected_id: Optional[str] = None,
        selected_xyz_ref: Optional[Path] = None,
    ) -> Path:
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / "manifest.json"
        profile_data: Dict[str, Any] = {}
        try:
            profile_data = json.loads(Path(profile).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            profile_data = {}
        selections = profile_data.get("selections", {}) or {}
        ts_selection = selections.get("ts_guess", {}) or {}
        intermediate_selection = selections.get("intermediate", {}) or {}

        payload: Dict[str, Any] = {
            "schema_version": "s2_peb_manifest_v2",
            "stage": "S2",
            "peb_peak": str(ts_seed),
            "intermediate_seed": str(intermediate_seed),
            "intermediate_guess": str(intermediate_seed),
            "forming_bonds": [list(item) for item in forming_bonds],
            "scan_profile": str(profile),
            "scan_plot": profile_data.get("scan_plot"),
            "ts_guess_index": ts_selection.get("index"),
            "intermediate_index": intermediate_selection.get("index"),
            "intermediate_selection_method": intermediate_selection.get("rule"),
            "intermediate_confidence": (profile_data.get("scan_quality", {}) or {}).get(
                "intermediate_confidence"
            ),
            "status": status,
            "confidence": confidence,
            "degraded_reasons": list(degraded),
        }
        if variant_id:
            payload["variant_id"] = variant_id
            payload["parent_variant_id"] = variant_id
        if selected_id:
            payload["selected_id"] = selected_id
        if selected_xyz_ref:
            payload["selected_xyz_ref"] = str(selected_xyz_ref)
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return path


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="RPH V4 fixed CENSO-LITE S0–S4 pipeline")
    parser.add_argument("--csv", type=Path, required=True, help="Trusted Reaxys-cleaned CSV")
    parser.add_argument("--rx-id", required=True, help="Dataset rx_id")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
        help="Console and file log level (default: INFO)",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        help="Pipeline log file; default: <output>/rph_v4.log",
    )
    parser.add_argument(
        "--no-color",
        action="store_true",
        help="Disable colored console output",
    )
    ui_group = parser.add_mutually_exclusive_group()
    ui_group.add_argument(
        "--ui",
        choices=("auto", "dashboard", "classic", "verbose"),
        help="Terminal UI mode; default comes from config (auto)",
    )
    ui_group.add_argument(
        "--compact",
        action="store_true",
        help="Print stage summaries and failures only",
    )
    ui_group.add_argument(
        "--verbose-ui",
        action="store_true",
        help="Print every structured progress update",
    )
    parser.add_argument(
        "--stop-after",
        choices=("s1", "s2", "s3", "s4"),
        default="s4",
        help="Stop after the selected stage; default: s4",
    )
    parser.add_argument(
        "--resume-policy",
        choices=("strict", "recompute", "use-existing-upstream"),
        help="Resume mismatch policy; default comes from config (strict)",
    )
    parser.add_argument(
        "--start-from",
        choices=("s1", "s2", "s3", "s4"),
        help="Use completed upstream artifacts and begin new work at this stage",
    )
    parser.add_argument(
        "--recompute-from",
        choices=("s0", "s1", "s2", "s3", "s4"),
        help="Explicitly invalidate this stage and all downstream stages",
    )
    args = parser.parse_args(argv)
    if args.start_from and args.resume_policy != "use-existing-upstream":
        parser.error("--start-from requires --resume-policy use-existing-upstream")
    if args.resume_policy == "use-existing-upstream" and not args.start_from:
        parser.error("--resume-policy use-existing-upstream requires --start-from")
    if args.start_from and args.recompute_from:
        parser.error("--start-from and --recompute-from cannot be combined")
    if args.start_from and V4Checkpoint.stage_order.index(args.start_from) > V4Checkpoint.stage_order.index(args.stop_after):
        parser.error("--start-from cannot be later than --stop-after")
    if args.recompute_from and V4Checkpoint.stage_order.index(args.recompute_from) > V4Checkpoint.stage_order.index(args.stop_after):
        parser.error("--recompute-from cannot be later than --stop-after")
    log_file = args.log_file or args.output / "rph_v4.log"
    color = not args.no_color and sys.stdout.isatty()
    setup_v4_logging(
        log_file=log_file,
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        rich_console=color,
    )
    logger.info("[V4] Pipeline log: %s", log_file)
    s0_record = load_s0_reaction_record(args.csv, args.rx_id)
    ui_mode = args.ui or ("classic" if args.compact else "verbose" if args.verbose_ui else None)
    orchestrator_kwargs: Dict[str, Any] = {"color": color}
    if ui_mode is not None:
        orchestrator_kwargs["ui_mode"] = ui_mode
    if args.resume_policy is not None:
        orchestrator_kwargs["resume_policy"] = args.resume_policy
    if args.start_from is not None:
        orchestrator_kwargs["start_from"] = args.start_from
    if args.recompute_from is not None:
        orchestrator_kwargs["recompute_from"] = args.recompute_from
    V4Orchestrator(args.config, **orchestrator_kwargs).run(
        s0_record,
        args.output,
        stop_after=args.stop_after,
    )
    return 0
