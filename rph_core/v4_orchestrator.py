"""V4 S0–S4 orchestrator with a fixed CENSO-LITE S1."""

from __future__ import annotations

import argparse
import hashlib
import inspect
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
from rph_core.steps.fidelity_profile import FidelityProfile
from rph_core.steps.refinement.engine import RefinementEngine
from rph_core.steps.refinement.manifest_io import read_refinement_manifest
from rph_core.steps.step3_lowlevel import LowLevelEngine  # backward-compat alias
from rph_core.steps.step4_highlevel import HighLevelEngine  # backward-compat alias
from rph_core.utils.config_loader import load_config
from rph_core.utils.bond_pairs import canonicalize_bond_pairs
from rph_core.utils.atom_mapping import (
    PebResolutionError,
    bind_atom_mapping_to_xyz,
    resolve_peb_forming_bonds,
    write_peb_mapping,
)
from rph_core.utils.artifact_reconciler import ArtifactReconciler
from rph_core.utils.json_io import write_json_atomic
from rph_core.utils.log_manager import setup_v4_logging
from rph_core.utils.provenance import (
    build_provenance,
    sha256_of_file,
    verify_provenance_chain,
)
from rph_core.utils.run_id import RUN_ID_FIELD, new_run_id
from rph_core.utils.shared_console import get_console
from rph_core.utils.stale_recovery import RecoveryAction, recover_stage
from rph_core.utils.stage_progress import StageProgressReporter
from rph_core.utils.ui_reporter import RichReporter
from rph_core.utils.v4_checkpoint import V4Checkpoint

logger = logging.getLogger(__name__)
SCHEMA_VERSION = "rph_v4_s0_s4_v1"
_BACKWARD_COMPAT_STAGE_ENGINES = (LowLevelEngine, HighLevelEngine)


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
        refresh_stages: Optional[Sequence[str]] = None,
        rerun_failed_s3: bool = False,
        s3_structure_ids: Optional[Sequence[str]] = None,
        rescue_only: bool = False,
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
        self.refresh_stages = {
            str(stage).strip().lower()
            for stage in (refresh_stages or ())
            if str(stage).strip()
        }
        self.rerun_failed_s3 = bool(rerun_failed_s3)
        self.rescue_only = bool(rescue_only)
        self.s3_structure_ids = tuple(
            structure_id
            for structure_id in (
                str(value).strip() for value in (s3_structure_ids or ())
            )
            if structure_id
        )
        unknown_refresh = self.refresh_stages - set(V4Checkpoint.stage_order)
        if unknown_refresh:
            raise ValueError(f"refresh_stages contains unknown stages: {sorted(unknown_refresh)}")
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
        if self.refresh_stages and (self.start_from or self.recompute_from):
            raise ValueError("refresh_stages cannot be combined with start_from or recompute_from")
        if (self.rerun_failed_s3 or self.s3_structure_ids) and self.refresh_stages:
            raise ValueError("S3 partial rerun cannot be combined with refresh_stages")
        if self.rerun_failed_s3 and self.s3_structure_ids:
            raise ValueError("rerun_failed_s3 and s3_structure_ids are mutually exclusive")
        protocol = str((self.config.get("step1", {}) or {}).get("protocol", "censo_lite"))
        if protocol != "censo_lite":
            raise ValueError("RPH V4 only supports step1.protocol=censo_lite")
        self.run_id: Optional[str] = None

    def run(
        self,
        s0_record: S0ReactionRecord,
        work_dir: Path,
        stop_after: str = "s4",
    ) -> Dict[str, Any]:
        """Run the pipeline and always restore the terminal UI."""

        self.run_id = new_run_id()
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
        if (self.rerun_failed_s3 or self.s3_structure_ids) and stop_after != "s3":
            raise ValueError("S3 partial rerun requires stop_after=s3")
        work_dir = Path(work_dir).resolve()
        work_dir.mkdir(parents=True, exist_ok=True)
        self._reused_signature_overrides: Dict[Tuple[str, Optional[str]], str] = {}
        self._recomputed_stages: set[str] = set()
        self._write_schema(work_dir, run_id=self.run_id)
        checkpoint = V4Checkpoint(work_dir)
        self._artifact_reconciler = ArtifactReconciler(
            work_dir,
            S1_MANIFEST_SCHEMA_VERSION,
        )
        checkpoint.write_config_snapshot(self.config)
        checkpoint.set_run_id(self._current_run_id())
        checkpoint.set_reaction_id(s0_record.rx_id)
        self._recover_stale_resume_state(checkpoint, work_dir)
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
        ui_reporter.set_run_id(self._current_run_id())
        ui_reporter.start()
        s0_reporter = StageProgressReporter(
            work_dir / "S0_Mechanism",
            "S0",
            default_fields=progress_fields,
            event_callback=ui_reporter.event_callback,
            run_id=self._current_run_id(),
            config=self.config,
        )
        s0_reporter.emit_stage_event("s0_started", rx_id=s0_record.rx_id)

        # A trusted reaction record is the original S0 contract: mapped
        # mechanistic bonds are resolved before geometry generation and remain
        # authoritative through S2.  Do not turn map numbers into XYZ indices.
        s0_manifest = work_dir / "S0_Mechanism" / "mechanism.json"
        s0_signature_payload = {
            "signature_schema": "s0_signature_v3",
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
            self._recomputed_stages.add("s0")
            logger.info("[V4] S0 restored from trusted reaction record rx_id=%s", s0_record.rx_id)
        stage_manifest_paths: Dict[str, Path] = {"s0": s0_manifest}

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
        # Precursor S1 has no PEB bond dependency.  Preserve its established
        # canonical geometry order so existing expensive S1 results remain
        # deterministically migratable; its own S1 atom table is authoritative.
        precursor_smiles = (
            getattr(s0_record, "canonical_precursor_smiles", "")
            or getattr(s0_record, "mapped_precursor_smiles", "")
            or ""
        )
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
                stage_manifest_paths=stage_manifest_paths,
            )
            return self._write_pipeline_result(
                work_dir,
                stop_after,
                dict(stage_manifest_paths),
                run_id=self._current_run_id(),
            )

        s1_reporter = StageProgressReporter(
            work_dir / "S1_ConfSearch",
            "S1",
            default_fields=progress_fields,
            event_callback=ui_reporter.event_callback,
            run_id=self._current_run_id(),
            config=self.config,
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
                    precursor_engine = self._instantiate_with_optional_run_id(
                        CensoLiteEngine,
                        self.config,
                        work_dir / "S1_ConfSearch",
                        "precursor",
                        optional_kwargs={
                            "s0_mechanism_json_path": s0_manifest,
                            "forming_bonds": forming_bonds,
                        },
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
                    self._recomputed_stages.add("s1")
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
                self._recomputed_stages.add("s1")
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
                s1_engine = self._instantiate_with_optional_run_id(
                    CensoLiteEngine,
                    self.config,
                    work_dir / "S1_ConfSearch",
                    directory_name,
                    optional_kwargs={
                        "s0_mechanism_json_path": s0_manifest,
                        "forming_bonds": forming_bonds,
                    },
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
                self._recomputed_stages.add("s1")
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
        if s1_stage_manifest is not None:
            stage_manifest_paths["s1"] = s1_stage_manifest
        if stop_after == "s1":
            self._write_run_manifest(
                work_dir,
                s0_record,
                variants,
                all_s1_results,
                all_s2_results,
                "s1",
                s0_manifest,
                stage_manifest_paths=stage_manifest_paths,
            )
            return self._write_pipeline_result(
                work_dir,
                stop_after,
                dict(stage_manifest_paths),
                run_id=self._current_run_id(),
            )

        s2_reporter = StageProgressReporter(
            work_dir / "S2_PEB",
            "S2",
            default_fields=progress_fields,
            event_callback=ui_reporter.event_callback,
            run_id=self._current_run_id(),
            config=self.config,
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
            s0_mapping_ref = str(
                variant.get("atom_mapping_ref")
                or (Path("atom_mappings") / f"{variant_id}.json").as_posix()
            )
            s0_mapping_path = work_dir / "S0_Mechanism" / s0_mapping_ref
            s1_manifest_path = Path(s1_result["manifest"])
            s1_mapping_ref = str(
                (s1_result.get("data") or {}).get("atom_mapping_ref")
                or "atom_mapping.json"
            )
            s1_mapping_path = s1_manifest_path.parent / s1_mapping_ref
            resolved_mapping = resolve_peb_forming_bonds(
                atom_map_smiles_path=s0_mapping_path,
                smiles_to_xyz_map_path=s1_mapping_path,
                product_xyz_path=product_xyz,
            )
            peb_mapping_path = write_peb_mapping(
                resolved_mapping,
                s2_dir / "atom_mapping.json",
                product_xyz=product_xyz,
            )
            peb_mapping_payload = self._read_json(peb_mapping_path)
            s2_product_mapping_table = bind_atom_mapping_to_xyz(
                s1_mapping_path,
                product_xyz,
                s2_dir / "product_atom_mapping.json",
            )
            variant_forming_bonds = resolved_mapping.forming_bonds_product_xyz_0based
            s2_payload = {
                "signature_schema": "s2_signature_v9",
                "schema_version": "s2_peb_manifest_v11",
                "stage": "s2",
                "variant_id": variant_id,
                "s0": s0_signature,
                "selected_id": (s1_result.get("data") or {}).get("selected"),
                "selected_xyz_ref": str(work_dir / "S1_ConfSearch" / variant_id / "selected.xyz"),
                "product_xyz": V4Checkpoint.file_signature(product_xyz),
                "forming_bonds": variant_forming_bonds,
                "s0_atom_mapping": V4Checkpoint.file_signature(s0_mapping_path),
                "s1_atom_mapping": V4Checkpoint.file_signature(s1_mapping_path),
                "peb_atom_mapping": V4Checkpoint.file_signature(peb_mapping_path),
                "config": self.config.get("step2", {}),
            }
            s2_signature = V4Checkpoint.signature(s2_payload)
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
                    intermediate_value = (
                        s2_data.get("intermediate_xyz")
                        or s2_data.get("intermediate_seed")
                    )
                    normalized_stage_status = self._normalize_progress_status(s2_data.get("status", "complete"))
                    all_s2_results[variant_id] = {
                        "status": "complete",
                        "manifest": s2_manifest,
                        "ts_seed": (
                            Path(s2_data["peb_peak"])
                            if s2_data.get("peb_peak")
                            else None
                        ),
                        "intermediate_xyz": Path(intermediate_value) if intermediate_value else None,
                        "intermediate_selection_status": s2_data.get(
                            "intermediate_selection_status"
                        ),
                        "s3_dispatch": dict(s2_data.get("s3_dispatch", {}) or {}),
                        "selection_source": s2_data.get("selection_source"),
                        "s2_state": s2_data.get("s2_state")
                        or (s2_data.get("s3_dispatch", {}) or {}).get("resolution")
                        or "unresolved",
                        "seed_evidence": s2_data.get("seed_evidence"),
                        "ts_search_seed": None
                        if s2_data.get("ts_search_seed") is None
                        else dict(s2_data.get("ts_search_seed") or {}),
                        "int_search_seed": None
                        if s2_data.get("int_search_seed") is None
                        else dict(s2_data.get("int_search_seed") or {}),
                        "has_independent_int": bool(s2_data.get("has_independent_int", False)),
                        "rejection_reason": s2_data.get("rejection_reason"),
                        "variant_id": s2_data.get("variant_id"),
                        "parent_variant_id": s2_data.get("parent_variant_id"),
                        "selected_id": s2_data.get("selected_id"),
                        "selected_xyz_ref": s2_data.get("selected_xyz_ref"),
                        "forming_bonds": tuple(tuple(pair) for pair in s2_data.get("forming_bonds", [])),
                        "atom_mapping": Path(s2_data.get("atom_mapping", peb_mapping_path)),
                        "mapping_status": s2_data.get("mapping_status", "verified"),
                        "atom_mapping_tables": {
                            key: Path(value)
                            for key, value in (s2_data.get("atom_mapping_tables") or {}).items()
                        },
                        "scan_profile": Path(s2_data["scan_profile"]) if s2_data.get("scan_profile") else None,
                        "stage_status": s2_data.get("status", "reused"),
                        "confidence": s2_data.get("confidence", "unknown"),
                        "degraded_reasons": tuple(s2_data.get("degraded_reasons", [])),
                        "checkpoint_signature": s2_signature,
                        "reused": True,
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
                s2_extra_metadata: Optional[Dict[str, Any]] = None
                scanner = PEBScanner(self.config, molecule_name=variant_id)
                scanner.event_callback = s2_reporter.external_event
                result = scanner.run(product_xyz, s2_dir, variant_forming_bonds)
                ts_seed, _legacy, intermediate_xyz, returned_bonds, profile, stage_status, confidence, degraded = result[:8]
                s2_method_actual = "legacy_peb"
                expected_bonds = canonicalize_bond_pairs(variant_forming_bonds)
                returned_bonds = canonicalize_bond_pairs(returned_bonds)
                if returned_bonds != expected_bonds:
                    raise RuntimeError(
                        f"S2 scanner changed verified forming bonds for {variant_id}: "
                        f"expected={variant_forming_bonds} returned={returned_bonds}"
                    )
                s2_mapping_tables = {
                    "product": s2_product_mapping_table,
                }
                if ts_seed is not None:
                    s2_mapping_tables["ts"] = bind_atom_mapping_to_xyz(
                        s1_mapping_path,
                        Path(ts_seed),
                        s2_dir / "ts_atom_mapping.json",
                    )
                if intermediate_xyz is not None:
                    s2_mapping_tables["intermediate"] = bind_atom_mapping_to_xyz(
                        s1_mapping_path,
                        Path(intermediate_xyz),
                        s2_dir / "intermediate_atom_mapping.json",
                    )
                s2_manifest = self._write_s2(
                    s2_dir,
                    ts_seed,
                    intermediate_xyz,
                    returned_bonds,
                    profile,
                    stage_status,
                    confidence,
                    degraded,
                    variant_id=variant_id,
                    selected_id=(s1_result.get("data") or {}).get("selected"),
                    selected_xyz_ref=work_dir / "S1_ConfSearch" / variant_id / "selected.xyz",
                    atom_mapping=peb_mapping_path,
                    atom_mapping_payload=peb_mapping_payload,
                    atom_mapping_tables=s2_mapping_tables,
                    s0_manifest_path=s0_manifest,
                    profile_data=getattr(scanner, "last_profile_payload", None),
                    s1_manifest_path=s1_manifest_path,
                    run_id=self._current_run_id(),
                    s2_method=s2_method_actual,
                    extra_metadata=s2_extra_metadata,
                )
                s2_manifest_data = self._read_json(s2_manifest)
                checkpoint.mark_scope(
                    variant_id,
                    "s2",
                    s2_signature,
                    s2_manifest,
                    signature_payload=s2_payload,
                )
                self._recomputed_stages.add("s2")
                s2_reporter.science_summary(
                    variant=variant_id,
                    forming_bonds=[list(pair) for pair in returned_bonds],
                    confidence=confidence,
                    degraded_reasons=list(degraded),
                    ts_seed=str(ts_seed) if ts_seed else None,
                    scan_profile=str(profile) if profile else None,
                )
                all_s2_results[variant_id] = {
                    "status": "complete",
                    "manifest": s2_manifest,
                    "ts_seed": Path(ts_seed) if ts_seed else None,
                    "intermediate_xyz": Path(intermediate_xyz) if intermediate_xyz else None,
                    "intermediate_selection_status": s2_manifest_data.get(
                        "intermediate_selection_status"
                    ),
                    "s3_dispatch": dict(s2_manifest_data.get("s3_dispatch", {}) or {}),
                    "selection_source": s2_manifest_data.get("selection_source"),
                    "s2_state": s2_manifest_data.get("s2_state")
                    or (s2_manifest_data.get("s3_dispatch", {}) or {}).get("resolution")
                    or "unresolved",
                    "seed_evidence": s2_manifest_data.get("seed_evidence"),
                    "ts_search_seed": None
                    if s2_manifest_data.get("ts_search_seed") is None
                    else dict(s2_manifest_data.get("ts_search_seed") or {}),
                    "int_search_seed": None
                    if s2_manifest_data.get("int_search_seed") is None
                    else dict(s2_manifest_data.get("int_search_seed") or {}),
                    "has_independent_int": bool(
                        s2_manifest_data.get("has_independent_int", False)
                    ),
                    "rejection_reason": s2_manifest_data.get("rejection_reason"),
                    "variant_id": variant_id,
                    "parent_variant_id": variant_id,
                    "selected_id": (s1_result.get("data") or {}).get("selected"),
                    "selected_xyz_ref": str(work_dir / "S1_ConfSearch" / variant_id / "selected.xyz"),
                    "forming_bonds": tuple(tuple(pair) for pair in returned_bonds),
                    "atom_mapping": peb_mapping_path,
                    "mapping_status": "verified",
                    "atom_mapping_tables": s2_mapping_tables,
                    "scan_profile": Path(profile) if profile else None,
                    "stage_status": stage_status,
                    "confidence": confidence,
                    "degraded_reasons": tuple(degraded),
                    "checkpoint_signature": s2_signature,
                    "reused": False,
                }
                s2_reporter.finish_structure(
                    variant_id,
                    self._normalize_progress_status(stage_status),
                    manifest=str(s2_manifest),
                    stage_status=stage_status,
                    confidence=confidence,
                    degraded_reasons=list(degraded),
                    ts_seed=str(ts_seed) if ts_seed else None,
                    intermediate_xyz=str(intermediate_xyz) if intermediate_xyz else None,
                    forming_bonds=[list(pair) for pair in returned_bonds],
                    scan_profile=str(profile) if profile else None,
                )
            except ResumeCheckpointError:
                raise
            except Exception as exc:
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
        if s2_stage_manifest is not None:
            stage_manifest_paths["s2"] = s2_stage_manifest
        if stop_after == "s2":
            self._write_run_manifest(
                work_dir,
                s0_record,
                variants,
                all_s1_results,
                all_s2_results,
                "s2",
                s0_manifest,
                stage_manifest_paths=stage_manifest_paths,
            )
            return self._write_pipeline_result(
                work_dir,
                stop_after,
                dict(stage_manifest_paths),
                run_id=self._current_run_id(),
            )

        self._assert_stage_results_run_id(all_s2_results, "S2", "S3")
        for variant in variants:
            variant_id = str(variant["variant_id"])
            s2_result = all_s2_results.get(variant_id, {})
            s1_result = all_s1_results.get(variant_id, {})
            if s2_result.get("status") != "complete" or s2_result.get("manifest") is None:
                continue
            expected_parent_paths = {"s0": s0_manifest}
            if s1_result.get("manifest") is not None:
                expected_parent_paths["s1"] = Path(s1_result["manifest"])
            self._warn_provenance_mismatch(
                Path(s2_result["manifest"]),
                boundary_label=f"S2→S3 ({variant_id})",
                expected_parent_paths=expected_parent_paths,
            )
        structures: List[Dict[str, Any]] = []
        if precursor_s1_result and precursor_s1_result.get("status") == "complete":
            precursor_mapping = self._s1_atom_mapping_path(precursor_s1_result)
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
                "atom_mapping": str(precursor_mapping),
                "mapping_required": True,
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
            atom_mapping_table = self._s1_atom_mapping_path(s1_result)
            verified_bonds = [list(b) for b in s2_result.get("forming_bonds", ())]
            s3_dispatch = dict(s2_result.get("s3_dispatch", {}) or {})
            # Generated S2 manifests always carry an explicit dispatch
            # resolution.  Only an explicit unresolved dispatch suppresses
            # S3; accepting an absent field keeps older/test doubles
            # backwards-compatible while they migrate to the manifest
            # contract.
            if s3_dispatch.get("resolution") == "unresolved":
                logger.info(
                    "[V4] Skipping S3 for variant %s: S2 selector/rescue is unresolved",
                    variant_id,
                )
                continue
            submit_ts = bool(s3_dispatch.get("submit_ts", True)) and bool(
                s2_result.get("ts_seed")
            )
            submit_intermediate = bool(s3_dispatch.get("submit_intermediate", False))
            variant_structures: List[Dict[str, Any]] = [
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
                    "atom_mapping": str(
                        (s2_result.get("atom_mapping_tables") or {}).get("product")
                        or atom_mapping_table
                    ),
                    "mapping_required": True,
                    **self._s1_ensemble_fields(s1_result),
                },
            ]
            if submit_ts:
                variant_structures.append({
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
                    "forming_bonds": verified_bonds,
                    "atom_mapping": str(
                        (s2_result.get("atom_mapping_tables") or {}).get("ts")
                        or atom_mapping_table
                    ),
                    "mapping_audit": str(s2_result["atom_mapping"]),
                    "mapping_required": True,
                })
            intermediate_path = s2_result.get("intermediate_xyz")
            if submit_intermediate and intermediate_path:
                variant_structures.insert(
                    1,
                    {
                        "id": f"{variant_id}_int",
                        "kind": "minimum",
                        "input_xyz": str(Path(intermediate_path)),
                        "structure_id": f"{variant_id}_int",
                        "variant_id": variant_id,
                        "role": "intermediate",
                        "branch_id": branch_id,
                        "pathway_id": pathway_id,
                        "parent_structure_id": variant_id,
                        "source_stage": "S2",
                        "seed_state": "stable_minimum_seed",
                        "forming_bonds": verified_bonds,
                        "atom_mapping": str(
                            (s2_result.get("atom_mapping_tables") or {}).get("intermediate")
                            or atom_mapping_table
                        ),
                        "mapping_audit": str(s2_result["atom_mapping"]),
                        "mapping_required": True,
                        "parent_structure": self._reference_geometry_bundle(
                            variant_id,
                            s1_result,
                            precursor_s1_result,
                        ),
                    },
                )
            structures.extend(variant_structures)
        s3_dir = work_dir / "S3_LowLevel"
        s3_manifest = s3_dir / "manifest.json"
        s3_reporter = StageProgressReporter(
            s3_dir,
            "S3",
            default_fields=progress_fields,
            event_callback=ui_reporter.event_callback,
            run_id=self._current_run_id(),
            config=self.config,
        )
        s3_reporter.emit_stage_event("s3_started", total_structures=len(structures))
        s3_payload = {
            "signature_schema": "s3_signature_v3",
            "stage": "s3",
            "structures": self._structure_signatures(structures),
            "config": self.config.get("refinement", {}).get("s3", {}),
        }
        s3_signature = V4Checkpoint.signature(s3_payload)
        s3_partial_rerun = bool(self.rerun_failed_s3 or self.s3_structure_ids)
        s3_reused = False if s3_partial_rerun else self._checkpoint_reusable(
            checkpoint,
            "s3",
            s3_signature,
            s3_manifest,
            signature_payload=s3_payload,
        )
        if s3_reused:
            logger.info("[V4] Reusing S3 manifest")
        else:
            s3_profile = FidelityProfile.from_config(self.config, "S3")
            s3_engine = RefinementEngine(
                config=self.config,
                profile=s3_profile,
                run_id=self._current_run_id(),
                parent_manifest_paths=dict(stage_manifest_paths),
            )
            s3_engine.set_progress_reporter(s3_reporter)
            s3_manifest =             s3_engine.run(
                structures,
                s3_dir,
                resume_incomplete=self.rerun_failed_s3,
                structure_ids=self.s3_structure_ids,
                rescue_only=self.rescue_only,
            )
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
        stage_manifest_paths["s3"] = s3_manifest
        if stop_after == "s3":
            self._write_run_manifest(
                work_dir,
                s0_record,
                variants,
                all_s1_results,
                all_s2_results,
                "s3",
                s0_manifest,
                stage_manifest_paths=stage_manifest_paths,
                s3_manifest=s3_manifest,
            )
            return self._write_pipeline_result(
                work_dir,
                stop_after,
                dict(stage_manifest_paths),
                run_id=self._current_run_id(),
            )
        self._assert_manifest_run_id(
            s3_manifest,
            "S3",
            "S4",
            allow_reused=s3_reused,
        )
        self._warn_provenance_mismatch(
            s3_manifest,
            boundary_label="S3→S4",
            expected_parent_paths={"s2": stage_manifest_paths.get("s2")},
        )
        s3_manifest_data = read_refinement_manifest(s3_manifest)
        s4_structures = [
            {
                "id": item["id"],
                "kind": item.get("kind", "minimum"),
                "input_xyz": (
                    item.get("opt_xyz")
                    or item.get("fallback_xyz")
                    or item["input_xyz"]
                ),
                "fallback_xyz": item.get("fallback_xyz") or item.get("input_xyz"),
                "original_seed_xyz": (
                    item.get("original_seed_xyz")
                    or item.get("fallback_xyz")
                    or item.get("input_xyz")
                ),
                "charge": item.get("charge"),
                "multiplicity": item.get("multiplicity"),
                "structure_id": item.get("structure_id", item["id"]),
                "variant_id": item.get("variant_id"),
                "role": item.get("role"),
                "branch_id": item.get("branch_id"),
                "pathway_id": item.get("pathway_id"),
                "parent_structure_id": item.get("parent_structure_id"),
                "source_stage": "S3",
                "source_s3_manifest": str(s3_manifest),
                "atom_mapping": item.get("atom_mapping"),
                "mapping_audit": item.get("mapping_audit"),
                "mapping_status": item.get("mapping_status"),
                "mapping_required": item.get("mapping_required", True),
                "parent_structure": item.get("parent_structure"),
                "s1_manifest": item.get("s1_manifest"),
                "s1_ensemble_thermodynamics": item.get("s1_ensemble_thermodynamics"),
                "s1_thermochemistry_status": item.get("s1_thermochemistry_status"),
                "ensemble_thermochemistry_correction_hartree": item.get(
                    "ensemble_thermochemistry_correction_hartree"
                ),
                **(
                    {"forming_bonds": item.get("forming_bonds")}
                    if item.get("forming_bonds") is not None
                    else {}
                ),
            }
            for item in s3_manifest_data.get("structures", [])
        ]
        s4_dir = work_dir / "S4_HighLevel"
        s4_manifest = s4_dir / "manifest.json"
        s4_payload = {
            "signature_schema": "s4_signature_v4",
            "stage": "s4",
            "structures": self._structure_signatures(s4_structures),
            "config": self.config.get("refinement", {}).get("s4", {}),
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
            s4_profile = FidelityProfile.from_config(self.config, "S4")
            s4_engine = RefinementEngine(
                config=self.config,
                profile=s4_profile,
                run_id=self._current_run_id(),
                parent_manifest_paths=dict(stage_manifest_paths),
            )
            s4_engine.set_event_callback(ui_reporter.s4_event_callback)
            s4_manifest = s4_engine.run(s4_structures, s4_dir)
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
        stage_manifest_paths["s4"] = s4_manifest
        self._write_run_manifest(
            work_dir,
            s0_record,
            variants,
            all_s1_results,
            all_s2_results,
            "s4",
            s0_manifest,
            stage_manifest_paths=stage_manifest_paths,
            s3_manifest=s3_manifest,
            s4_manifest=s4_manifest,
        )
        return self._write_pipeline_result(
            work_dir,
            stop_after,
            dict(stage_manifest_paths),
            run_id=self._current_run_id(),
        )

    def _current_run_id(self) -> str:
        run_id = self.run_id
        if not isinstance(run_id, str) or not run_id.strip():
            raise RuntimeError("V4Orchestrator.run_id is not initialized")
        return run_id

    def _recover_stale_resume_state(
        self,
        checkpoint: V4Checkpoint,
        work_dir: Path,
    ) -> None:
        ui_config = dict(self.config.get("ui", {}) or {})
        recovery_config = dict(ui_config.get("recovery") or {})
        if not bool(recovery_config.get("enabled", True)):
            return
        stale_threshold_seconds = float(
            recovery_config.get(
                "stale_heartbeat_seconds",
                ui_config.get("stalled_job_warning_seconds", 720),
            )
        )
        requeue_interrupted = bool(
            recovery_config.get("requeue_interrupted_structures", True)
        )
        for stage_name in ("s3", "s4"):
            stage_dir = Path(work_dir) / self._stage_output_dir_name(stage_name)
            if not stage_dir.exists():
                continue
            had_stage_entry = checkpoint.has_stage(stage_name)
            action = recover_stage(
                stage_name,
                stage_dir,
                checkpoint,
                stale_threshold_seconds=stale_threshold_seconds,
            )
            logger.info(
                "[V4] Resume recovery stage=%s action=%s stage_dir=%s",
                stage_name,
                action.value,
                stage_dir,
            )
            if (
                action == RecoveryAction.MARK_INTERRUPTED
                and requeue_interrupted
                and had_stage_entry
                and checkpoint.has_stage(stage_name)
            ):
                checkpoint.invalidate_from(stage_name)
                logger.info(
                    "[V4] Resume recovery invalidated stale checkpoint stage=%s",
                    stage_name,
                )

    @staticmethod
    def _stage_output_dir_name(stage_name: str) -> str:
        mapping = {
            "s3": "S3_LowLevel",
            "s4": "S4_HighLevel",
        }
        try:
            return mapping[str(stage_name).strip().lower()]
        except KeyError as exc:
            raise ValueError(f"Unsupported recovery stage: {stage_name}") from exc

    def _instantiate_with_optional_run_id(
        self,
        factory: Any,
        *args: Any,
        optional_kwargs: Mapping[str, Any] | None = None,
    ) -> Any:
        kwargs: Dict[str, Any] = {}
        try:
            signature = inspect.signature(factory)
        except (TypeError, ValueError):
            signature = None
        accepts_var_kwargs = bool(
            signature is not None
            and any(
                parameter.kind == inspect.Parameter.VAR_KEYWORD
                for parameter in signature.parameters.values()
            )
        )
        if signature is not None and (
            RUN_ID_FIELD in signature.parameters
            or accepts_var_kwargs
        ):
            kwargs[RUN_ID_FIELD] = self.run_id
        for key, value in (optional_kwargs or {}).items():
            if signature is None or key in signature.parameters or accepts_var_kwargs:
                kwargs[key] = value
        return factory(*args, **kwargs)

    def _warn_provenance_mismatch(
        self,
        manifest_path: Path,
        *,
        boundary_label: str,
        expected_parent_paths: Mapping[str, Path | None],
    ) -> None:
        comparable_paths = {
            str(stage): Path(path)
            for stage, path in expected_parent_paths.items()
            if path is not None
        }
        if not comparable_paths:
            return
        ok, mismatches = verify_provenance_chain(
            manifest_path,
            expected_parent_paths=comparable_paths,
        )
        if not ok:
            logger.warning(
                "[V4] Soft provenance check failed at %s for %s: %s",
                boundary_label,
                manifest_path,
                "; ".join(mismatches),
            )

    def _assert_stage_results_run_id(
        self,
        results: Mapping[str, Mapping[str, Any]],
        stage_name: str,
        next_stage_name: str,
    ) -> None:
        for result in results.values():
            if result.get("status") != "complete" or result.get("manifest") is None:
                continue
            self._assert_manifest_run_id(
                Path(result["manifest"]),
                stage_name,
                next_stage_name,
                allow_reused=bool(result.get("reused", False)),
            )

    def _assert_manifest_run_id(
        self,
        manifest_path: Path,
        stage_name: str,
        next_stage_name: str,
        *,
        allow_reused: bool = False,
    ) -> None:
        payload = self._read_json(manifest_path)
        manifest_run_id = payload.get(RUN_ID_FIELD)
        current_run_id = self._current_run_id()
        if manifest_run_id != current_run_id:
            if allow_reused:
                logger.info(
                    "[V4] Allowing validated reused %s manifest from run_id %s while current run_id is %s",
                    stage_name,
                    manifest_run_id,
                    current_run_id,
                )
                return
            raise RuntimeError(
                f"{stage_name} manifest run_id {manifest_run_id} does not match current run_id {current_run_id}. "
                f"Stale {stage_name} output detected. Refusing to start {next_stage_name}."
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
    def _reference_geometry_bundle(
        variant_id: str,
        s1_result: Dict[str, Any],
        precursor_s1_result: Optional[Dict[str, Any]],
    ) -> Dict[str, str]:
        """Reference geometries consumed by ``classify_int`` collapse detection.

        ``classify_int`` compares an intermediate against the precursor and
        product minima (mapped RMSD + per-bond progress).  The product
        reference is the S1 selected conformer of the same variant (shared
        atom ordering with the PEB intermediate seed); the precursor
        reference is the S1 selected conformer of the precursor.
        """
        bundle: Dict[str, str] = {}
        product_xyz = s1_result.get("selected_xyz")
        if product_xyz:
            bundle["product_ref"] = str(Path(product_xyz))
        if precursor_s1_result is not None:
            precursor_xyz = precursor_s1_result.get("selected_xyz")
            if precursor_xyz:
                bundle["precursor_ref"] = str(Path(precursor_xyz))
        return bundle

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
        run_id: Optional[str] = None,
    ) -> Dict[str, Any]:
        result: Dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            RUN_ID_FIELD: run_id,
            "success": True,
            "completed_through": completed_through,
            "work_dir": str(work_dir),
        }
        for stage in ("s0", "s1", "s2", "s3", "s4"):
            manifest = manifests.get(stage)
            if manifest is not None:
                result[f"{stage}_manifest"] = str(Path(manifest))
        result["provenance"] = build_provenance(
            run_id=run_id,
            schema_version=SCHEMA_VERSION,
            protocol_version=None,
            parent_manifest_paths={
                stage: manifests[stage]
                for stage in ("s0", "s1", "s2", "s3", "s4")
                if stage in manifests
            },
            atom_mapping_payload=None,
            forming_bonds=None,
            variant_manifest_path=None,
            extra={"completed_through": completed_through},
        ).to_dict()
        write_json_atomic(Path(work_dir) / "pipeline.result.json", result)
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
    def _s1_atom_mapping_path(result: Mapping[str, Any]) -> Path:
        manifest_path = Path(result["manifest"])
        data = dict(result.get("data") or {})
        mapping_ref = str(data.get("atom_mapping_ref") or "atom_mapping.json")
        mapping_path = manifest_path.parent / mapping_ref
        if not mapping_path.is_file():
            raise ResumeCheckpointError(f"S1 atom mapping sidecar is missing: {mapping_path}")
        return mapping_path

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
            "protocol_version": "censo_light_ranking_v4",
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

        if stage in set(getattr(self, "refresh_stages", set()) or set()):
            logger.warning(
                "[V4] Refresh requested for %s%s; other validated artifacts remain reusable",
                stage,
                f":{scope}" if scope else "",
            )
            return False

        decision = checkpoint.check_reuse(
            stage,
            signature,
            manifest,
            scope=scope,
            signature_payload=signature_payload,
        )
        explicit_recompute = getattr(self, "recompute_from", None)
        reconciliation_allowed = not (
            explicit_recompute
            and V4Checkpoint.stage_order.index(stage)
            >= V4Checkpoint.stage_order.index(explicit_recompute)
        )
        if stage == "s1" and scope and signature_payload is not None and reconciliation_allowed:
            state = checkpoint.load()
            stage_record = ((state.get("stages") or {}).get("s1") or {})
            entry = ((stage_record.get("scopes") or {}).get(scope) or {})
            recorded_payload = dict(entry.get("signature_payload") or {})
            recorded_smiles = recorded_payload.get("smiles")
            reconciler = getattr(self, "_artifact_reconciler", None)
            if reconciler is None:
                reconciler = ArtifactReconciler(checkpoint.work_dir, S1_MANIFEST_SCHEMA_VERSION)
            assessment = reconciler.assess_s1(
                manifest,
                expected_smiles=str(signature_payload.get("smiles") or ""),
                recorded_smiles=str(recorded_smiles) if recorded_smiles else None,
                changed_paths=[str(row.get("path", "")) for row in decision.differences],
            )
            if assessment.status == "reusable_after_migration":
                try:
                    assessment = reconciler.migrate_s1(assessment)
                except (OSError, ValueError, PebResolutionError) as exc:
                    logger.warning("[V4] S1 artifact migration rejected for %s: %s", scope, exc)
            if assessment.reusable:
                if not decision.reusable or assessment.status != "reusable_exact":
                    checkpoint.mark_scope(
                        scope,
                        "s1",
                        signature,
                        manifest,
                        signature_payload=signature_payload,
                    )
                    logger.info(
                        "[V4] Rebuilt S1 checkpoint from validated artifacts for %s (%s)",
                        scope,
                        assessment.status,
                    )
                return True
        if decision.reusable:
            return True

        recomputed_stages = set(getattr(self, "_recomputed_stages", set()) or set())
        current_index = V4Checkpoint.stage_order.index(stage)
        recomputed_upstream = [
            upstream
            for upstream in recomputed_stages
            if V4Checkpoint.stage_order.index(upstream) < current_index
        ]
        if recomputed_upstream and decision.reason in {
            "signature_mismatch",
            "manifest_mismatch",
            "manifest_missing",
            "manifest_unrecorded",
        }:
            logger.warning(
                "[V4] Recomputing %s%s after upstream %s changed: %s",
                stage,
                f":{scope}" if scope else "",
                ", ".join(sorted(recomputed_upstream)),
                decision.describe(),
            )
            return False

        migration_fields = {
            "signature_schema",
            "manifest_schema",
            "protocol_version",
            "schema_version",
        }
        changed_paths = {
            str(row.get("path", "")).rsplit(".", 1)[-1]
            for row in decision.differences
        }
        if decision.reason == "signature_mismatch" and changed_paths & migration_fields:
            logger.warning(
                "[V4] Recomputing %s%s for checkpoint contract migration: %s",
                stage,
                f":{scope}" if scope else "",
                decision.describe(),
            )
            return False
        if decision.reason == "signature_mismatch" and stage == "s1" and "smiles" in changed_paths:
            logger.warning(
                "[V4] Recomputing %s%s because its S0-derived SMILES input changed: %s",
                stage,
                f":{scope}" if scope else "",
                decision.describe(),
            )
            return False

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
        stage_manifest_paths: Mapping[str, Path] | None = None,
        s3_manifest: Optional[Path] = None,
        s4_manifest: Optional[Path] = None,
    ) -> Path:
        condition_fields, condition_signature = self._condition_metadata(record)
        manifest: Dict[str, Any] = {
            "schema_version": "rph_run_manifest_v1",
            RUN_ID_FIELD: self._current_run_id(),
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
        available_stage_manifests: Dict[str, Path] = {"s0": Path(s0_manifest)}
        for stage, manifest_path in (stage_manifest_paths or {}).items():
            available_stage_manifests[str(stage)] = Path(manifest_path)
        if s3_manifest is not None:
            available_stage_manifests["s3"] = Path(s3_manifest)
        if s4_manifest is not None:
            available_stage_manifests["s4"] = Path(s4_manifest)
        manifest["provenance"] = build_provenance(
            run_id=self._current_run_id(),
            schema_version=str(manifest["schema_version"]),
            protocol_version=None,
            parent_manifest_paths={
                stage: available_stage_manifests[stage]
                for stage in ("s0", "s1", "s2", "s3", "s4")
                if stage in available_stage_manifests
            },
            atom_mapping_payload=None,
            forming_bonds=None,
            variant_manifest_path=None,
            extra={"condition_signature": condition_signature},
        ).to_dict()
        path = Path(work_dir) / "run.manifest.json"
        write_json_atomic(path, manifest)
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
    def _write_schema(work_dir: Path, run_id: Optional[str] = None) -> None:
        write_json_atomic(
            work_dir / ".rph_schema.json",
            {"schema_version": SCHEMA_VERSION, RUN_ID_FIELD: run_id},
        )

    def _write_s0_from_record(self, work_dir: Path, record: S0ReactionRecord) -> Path:
        _reaction_context, _product_variants, mechanism_path = write_s0_artifacts(
            work_dir,
            record,
            run_id=self.run_id,
        )
        return mechanism_path

    @staticmethod
    def _write_s2(
        directory: Path,
        ts_seed: Optional[Path],
        intermediate_xyz: Optional[Path],
        forming_bonds: Sequence[Tuple[int, int]],
        profile: Optional[Path],
        status: str,
        confidence: str,
        degraded: Sequence[str],
        variant_id: Optional[str] = None,
        selected_id: Optional[str] = None,
        selected_xyz_ref: Optional[Path] = None,
        atom_mapping: Optional[Path] = None,
        atom_mapping_payload: Optional[Mapping[str, Any]] = None,
        atom_mapping_tables: Optional[Mapping[str, Path]] = None,
        s0_manifest_path: Optional[Path] = None,
        s1_manifest_path: Optional[Path] = None,
        run_id: Optional[str] = None,
        s2_method: str = "legacy_peb",
        extra_metadata: Optional[Mapping[str, Any]] = None,
        profile_data: Optional[Mapping[str, Any]] = None,
    ) -> Path:
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / "manifest.json"
        resolved_profile_data: Dict[str, Any]
        if profile_data is not None:
            resolved_profile_data = dict(profile_data)
        elif profile is not None:
            try:
                resolved_profile_data = json.loads(Path(profile).read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                resolved_profile_data = {}
        else:
            resolved_profile_data = {}
        profile_data = resolved_profile_data
        selections = profile_data.get("selections", {}) or {}
        ts_selection = selections.get("ts_guess", {}) or {}
        intermediate_selection = selections.get("intermediate", {}) or {}
        path_selection = profile_data.get("selection_decision", {}) or {}
        s3_dispatch = profile_data.get("s3_dispatch", {}) or {}
        resolution = s3_dispatch.get("resolution") or "unresolved"
        s2_state = profile_data.get("s2_state")
        if s2_state is None:
            s2_state = resolution
        shared_ts_int = intermediate_selection.get("selection_mode") == "shared_ts_fallback"

        payload: Dict[str, Any] = {
            "schema_version": "s2_peb_manifest_v11",
            "stage": "S2",
            RUN_ID_FIELD: run_id,
            "method": s2_method,
            "peb_peak": str(ts_seed) if ts_seed else None,
            "intermediate_xyz": str(intermediate_xyz) if intermediate_xyz else None,
            "forming_bonds": [list(item) for item in forming_bonds],
            "mapping_status": "verified" if atom_mapping else "missing",
            "scan_profile": str(profile) if profile is not None else None,
            "scan_plot": profile_data.get("scan_plot"),
            "ts_guess_index": ts_selection.get("index"),
            "ts_energy_peak_index": ts_selection.get("energy_peak_index", ts_selection.get("index")),
            "ts_seed_index": ts_selection.get("seed_index", ts_selection.get("index")),
            "ts_seed_reactant_backoff_A": ts_selection.get("seed_backoff_applied_A", 0.0),
            "ts_selection_method": ts_selection.get("rule"),
            "ts_selection_configured_method": ts_selection.get("configured_method"),
            "ts_selection_actual_method": ts_selection.get("actual_method"),
            "selected_path_branch": path_selection.get("selected_branch"),
            "path_selection_rule": path_selection.get("rule"),
            "intermediate_index": intermediate_selection.get("index"),
            "intermediate_selection_method": intermediate_selection.get("rule"),
            "intermediate_selection_status": intermediate_selection.get(
                "selection_status",
                (profile_data.get("scan_quality", {}) or {}).get(
                    "intermediate_selection_status"
                ),
            ),
            "selection_energy_source": (profile_data.get("selection_policy", {}) or {}).get(
                "actual_source",
                (profile_data.get("scan_quality", {}) or {}).get("selection_energy_source"),
            ),
            "scan_profile_schema": profile_data.get("profile_schema_version"),
            "scan_attempt_count": len(profile_data.get("attempts", []) or []),
            "composite_point_count": len(
                ((profile_data.get("composite_profile", {}) or {}).get("points", []) or [])
            ),
            "intermediate_confidence": (profile_data.get("scan_quality", {}) or {}).get(
                "intermediate_confidence"
            ),
            "selection_source": profile_data.get("selection_source"),
            "resolution": resolution,
            "s2_state": s2_state,
            "seed_evidence": profile_data.get("seed_evidence"),
            "ts_search_seed": profile_data.get("ts_search_seed"),
            "int_search_seed": profile_data.get("int_search_seed"),
            "has_independent_int": bool(profile_data.get("has_independent_int", False)),
            "rejection_reason": profile_data.get("rejection_reason"),
            "s3_dispatch": s3_dispatch,
            "rescue": profile_data.get("rescue", {}),
            "status": status,
            "confidence": confidence,
            "degraded_reasons": list(degraded),
            "nodes": [
                {
                    "role": "ts",
                    "xyz": str(ts_seed) if ts_seed else None,
                    "frame_index": ts_selection.get("index"),
                    "point_id": ts_selection.get("point_id"),
                    "source_attempt": ts_selection.get("source_attempt"),
                    "selection_rule": ts_selection.get("rule"),
                    "selection_actual_method": ts_selection.get("actual_method"),
                    "energy_source": ts_selection.get("energy_source"),
                    "confidence": confidence,
                    "usable_for_s3": bool(s3_dispatch.get("submit_ts", False)),
                },
                {
                    "role": "intermediate",
                    "xyz": str(intermediate_xyz) if intermediate_xyz else None,
                    "frame_index": intermediate_selection.get("index"),
                    "point_id": intermediate_selection.get("point_id"),
                    "source_attempt": intermediate_selection.get("source_attempt"),
                    "selection_rule": intermediate_selection.get("rule"),
                    "energy_source": intermediate_selection.get("energy_source"),
                    "status": (
                        "selected"
                        if intermediate_xyz
                        else "shared_with_ts" if shared_ts_int else "not_found"
                    ),
                    "selection_status": intermediate_selection.get("selection_status"),
                    "reason": intermediate_selection.get("reason") or intermediate_selection.get("error"),
                    "shared_with_ts": shared_ts_int,
                    "confidence": (profile_data.get("scan_quality", {}) or {}).get(
                        "intermediate_confidence"
                    ),
                    "usable_for_s3": bool(
                        intermediate_xyz and s3_dispatch.get("submit_intermediate", False)
                    ),
                },
            ],
        }
        if variant_id:
            payload["variant_id"] = variant_id
            payload["parent_variant_id"] = variant_id
        if selected_id:
            payload["selected_id"] = selected_id
        if selected_xyz_ref:
            payload["selected_xyz_ref"] = str(selected_xyz_ref)
        if atom_mapping:
            payload["atom_mapping"] = str(atom_mapping)
        if atom_mapping_tables:
            payload["atom_mapping_tables"] = {
                str(key): str(value) for key, value in atom_mapping_tables.items()
            }
        if extra_metadata:
            payload.update(
                {key: value for key, value in extra_metadata.items() if value is not None}
            )
        payload["provenance"] = build_provenance(
            run_id=run_id,
            schema_version=str(payload["schema_version"]),
            protocol_version=None,
            parent_manifest_paths={"s0": s0_manifest_path, "s1": s1_manifest_path},
            atom_mapping_payload=atom_mapping_payload,
            forming_bonds=forming_bonds,
            variant_manifest_path=s1_manifest_path,
            extra={
                "ts_guess_xyz_hash": sha256_of_file(ts_seed) if ts_seed else "",
                "intermediate_xyz_hash": sha256_of_file(intermediate_xyz) if intermediate_xyz else "",
            },
        ).to_dict()
        write_json_atomic(path, payload)
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
    parser.add_argument(
        "--refresh-stages",
        help="Comma-separated stages to rerun while retaining other validated artifacts, e.g. s0,s2,s3",
    )
    s3_partial_group = parser.add_mutually_exclusive_group()
    s3_partial_group.add_argument(
        "--rerun-failed-s3",
        action="store_true",
        help="Preserve complete S3 structures and rerun only failed or artifact-missing S3 structures",
    )
    s3_partial_group.add_argument(
        "--s3-structures",
        help="Comma-separated S3 structure ids to rerun while preserving all other S3 structures",
    )
    parser.add_argument(
        "--rescue-only",
        action="store_true",
        help="S3 partial rerun: replay failed structures directly into the rescue matrix, skipping a full pass1 re-optimization when a usable stop geometry exists",
    )
    args = parser.parse_args(argv)
    if args.start_from and args.resume_policy != "use-existing-upstream":
        parser.error("--start-from requires --resume-policy use-existing-upstream")
    if args.resume_policy == "use-existing-upstream" and not args.start_from:
        parser.error("--resume-policy use-existing-upstream requires --start-from")
    if args.start_from and args.recompute_from:
        parser.error("--start-from and --recompute-from cannot be combined")
    if args.refresh_stages and (args.start_from or args.recompute_from):
        parser.error("--refresh-stages cannot be combined with --start-from or --recompute-from")
    if (args.rerun_failed_s3 or args.s3_structures) and args.refresh_stages:
        parser.error("S3 partial rerun cannot be combined with --refresh-stages")
    if (args.rerun_failed_s3 or args.s3_structures) and args.stop_after != "s3":
        parser.error("S3 partial rerun requires --stop-after s3")
    if args.rescue_only and not (args.rerun_failed_s3 or args.s3_structures):
        parser.error("--rescue-only requires --rerun-failed-s3 or --s3-structures")
    if args.rescue_only and args.stop_after != "s3":
        parser.error("--rescue-only requires --stop-after s3")
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
    if args.rerun_failed_s3:
        orchestrator_kwargs["rerun_failed_s3"] = True
    if args.rescue_only:
        orchestrator_kwargs["rescue_only"] = True
    if args.s3_structures:
        orchestrator_kwargs["s3_structure_ids"] = tuple(
            item.strip() for item in args.s3_structures.split(",") if item.strip()
        )
    if args.start_from is not None:
        orchestrator_kwargs["start_from"] = args.start_from
    if args.recompute_from is not None:
        orchestrator_kwargs["recompute_from"] = args.recompute_from
    if args.refresh_stages:
        orchestrator_kwargs["refresh_stages"] = tuple(
            stage.strip() for stage in args.refresh_stages.split(",") if stage.strip()
        )
    V4Orchestrator(args.config, **orchestrator_kwargs).run(
        s0_record,
        args.output,
        stop_after=args.stop_after,
    )
    return 0
