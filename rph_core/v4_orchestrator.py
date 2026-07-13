"""V4 S0–S4 orchestrator with a fixed CENSO-LITE S1."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from rph_core.steps.conformer_search.censo_lite import CensoLiteEngine
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
from rph_core.utils.variant_checkpoint import VariantCheckpoint
from rph_core.utils.v4_checkpoint import V4Checkpoint

logger = logging.getLogger(__name__)
SCHEMA_VERSION = "rph_v4_s0_s4_v1"


class V4Orchestrator:
    """Run S0–S4 using manifest-based handoffs and per-structure failures."""

    def __init__(self, config_path: Optional[Path] = None, color: bool = True):
        self.config = load_config(config_path) if config_path else load_config()
        self.color = bool(color)
        protocol = str((self.config.get("step1", {}) or {}).get("protocol", "censo_lite"))
        if protocol != "censo_lite":
            raise ValueError("RPH V4 only supports step1.protocol=censo_lite")

    def run(
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
        self._write_schema(work_dir)
        checkpoint = VariantCheckpoint(work_dir)
        legacy_checkpoint = V4Checkpoint(work_dir)
        checkpoint.set_reaction_id(s0_record.rx_id)
        product_smiles = s0_record.product_smiles
        _, condition_signature = self._condition_metadata(s0_record)
        progress_fields = {
            "reaction_id": s0_record.rx_id,
            "condition_signature": condition_signature,
        }
        ui_reporter = RichReporter(get_console(), color=bool(getattr(self, "color", True)))
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
        s0_signature = V4Checkpoint.signature({
            "stage": "s0",
            "source": "trusted_reaction_record",
            "condition": condition_signature,
            "record": s0_record.signature_payload(),
        })
        s0_reused = checkpoint.reusable("s0", "s0", s0_signature, s0_manifest)
        if s0_reused:
            s0_data = self._read_json(s0_manifest)
            forming_bonds = tuple(tuple(pair) for pair in s0_data.get("forming_bonds", []))
            logger.info("[V4] Reusing record-driven S0 manifest")
        else:
            forming_bonds = s0_record.forming_bonds
            s0_manifest = self._write_s0_from_record(work_dir, s0_record)
            checkpoint.mark("s0", "s0", s0_signature, s0_manifest)
            logger.info("[V4] S0 restored from trusted reaction record rx_id=%s", s0_record.rx_id)
        if not legacy_checkpoint.reusable("s0", s0_signature, s0_manifest):
            legacy_checkpoint.mark_s0(s0_signature, s0_manifest)

        variants = self._load_variants(work_dir, product_smiles)
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
            precursor_s1_signature = V4Checkpoint.signature({
                "stage": "s1",
                "config": self.config.get("step1", {}),
                "smiles": precursor_smiles,
                "variant_id": "precursor",
                "directory_name": "precursor",
            })
            if checkpoint.reusable("precursor_s1", "s1", precursor_s1_signature, precursor_s1_manifest):
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
                    "checkpoint_signature": precursor_s1_signature,
                }
                s1_reporter.finish_structure(
                    "precursor",
                    "complete",
                    role="precursor",
                    smiles=precursor_smiles,
                    selected=precursor_s1_data.get("selected"),
                    manifest=str(precursor_s1_manifest),
                    reused=True,
                )
                logger.info("[V4] Reusing precursor S1 manifest")
            else:
                s1_reporter.start_structure("precursor", role="precursor", smiles=precursor_smiles)
                try:
                    raw_precursor_s1 = CensoLiteEngine(self.config, work_dir / "S1_ConfSearch", "precursor").run(precursor_smiles)
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
                    checkpoint.mark("precursor_s1", "s1", precursor_s1_signature, manifest_path)
                    s1_reporter.finish_structure(
                        "precursor",
                        "complete",
                        role="precursor",
                        smiles=precursor_smiles,
                        selected=precursor_s1_data.get("selected"),
                        manifest=str(manifest_path),
                    )
                except Exception as exc:
                    checkpoint.mark_failed("precursor_s1", "s1", precursor_s1_signature, str(exc))
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
            s1_signature = V4Checkpoint.signature({
                "stage": "s1",
                "config": self.config.get("step1", {}),
                "smiles": variant_smiles,
                "variant_id": variant_id,
                "directory_name": directory_name,
            })
            if shared_legacy_s1 is not None:
                shared_result = dict(shared_legacy_s1)
                shared_result.update(
                    status="complete",
                    variant_id=variant_id,
                    directory_name=directory_name,
                    checkpoint_signature=s1_signature,
                )
                checkpoint.mark(variant_id, "s1", s1_signature, Path(shared_result["manifest"]))
                all_s1_results[variant_id] = shared_result
                s1_reporter.finish_structure(
                    variant_id,
                    "complete",
                    smiles=variant_smiles,
                    selected=(shared_result.get("data") or {}).get("selected"),
                    manifest=str(shared_result["manifest"]),
                    reused=True,
                )
                continue
            if checkpoint.reusable(variant_id, "s1", s1_signature, s1_manifest):
                s1_data = self._read_json(s1_manifest)
                all_s1_results[variant_id] = {
                    "status": "complete",
                    "manifest": s1_manifest,
                    "selected_xyz": self._selected_s1_xyz(work_dir, s1_data, directory_name, s1_manifest),
                    "data": s1_data,
                    "variant_id": variant_id,
                    "directory_name": directory_name,
                    "checkpoint_signature": s1_signature,
                }
                s1_reporter.finish_structure(
                    variant_id,
                    "complete",
                    smiles=variant_smiles,
                    selected=s1_data.get("selected"),
                    manifest=str(s1_manifest),
                    reused=True,
                )
                logger.info("[V4] Reusing S1 manifest for variant %s", variant_id)
                continue
            s1_reporter.start_structure(variant_id, smiles=variant_smiles)
            try:
                raw_s1 = CensoLiteEngine(self.config, work_dir / "S1_ConfSearch", directory_name).run(variant_smiles)
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
                checkpoint.mark(variant_id, "s1", s1_signature, manifest_path)
                all_s1_results[variant_id] = s1_result
                s1_reporter.finish_structure(
                    variant_id,
                    "complete",
                    smiles=variant_smiles,
                    selected=s1_data.get("selected"),
                    manifest=str(manifest_path),
                )
                if manifest_path.parent.resolve() != s1_manifest.parent.resolve():
                    shared_legacy_s1 = dict(s1_result)
                    logger.warning(
                        "[V4] S1 for variant %s returned legacy path %s; remaining variants will reuse the shared output",
                        variant_id,
                        manifest_path,
                    )
            except Exception as exc:
                checkpoint.mark_failed(variant_id, "s1", s1_signature, str(exc))
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
        if s1_stage_manifest is not None and not legacy_checkpoint.reusable("s1", s1_stage_signature, s1_stage_manifest):
            legacy_checkpoint.mark("s1", s1_stage_signature, s1_stage_manifest)
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
            s2_signature = V4Checkpoint.signature({
                "stage": "s2",
                "variant_id": variant_id,
                "s0": s0_signature,
                "selected_id": (s1_result.get("data") or {}).get("selected"),
                "selected_xyz_ref": str(work_dir / "S1_ConfSearch" / variant_id / "selected.xyz"),
                "product_xyz": V4Checkpoint.file_signature(product_xyz),
                "forming_bonds": forming_bonds,
                "config": self.config.get("step2", {}),
            })
            try:
                if checkpoint.reusable(variant_id, "s2", s2_signature, s2_manifest):
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
                    logger.info("[V4] Reusing S2 manifest for variant %s", variant_id)
                    continue
                s2_reporter.start_structure(
                    variant_id,
                    selected_id=(s1_result.get("data") or {}).get("selected"),
                    product_xyz=str(product_xyz),
                )
                scanner = PEBScanner(self.config, molecule_name=variant_id)
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
                checkpoint.mark(variant_id, "s2", s2_signature, s2_manifest)
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
                )
            except Exception as exc:
                checkpoint.mark_failed(variant_id, "s2", s2_signature, str(exc))
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
        if s2_stage_manifest is not None and not legacy_checkpoint.reusable("s2", s2_stage_signature, s2_stage_manifest):
            legacy_checkpoint.mark("s2", s2_stage_signature, s2_stage_manifest)
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
        s3_signature = V4Checkpoint.signature({
            "stage": "s3",
            "structures": self._structure_signatures(structures),
            "config": self.config.get("theory", {}).get("s3_low_level", {}),
        })
        s3_reused = checkpoint.reusable("s3", "s3", s3_signature, s3_manifest)
        if s3_reused:
            logger.info("[V4] Reusing S3 manifest")
        else:
            s3_engine = LowLevelEngine(self.config)
            if hasattr(s3_engine, "set_progress_reporter"):
                s3_engine.set_progress_reporter(s3_reporter)
            s3_manifest = s3_engine.run(structures, s3_dir)
            checkpoint.mark("s3", "s3", s3_signature, s3_manifest)
        if s3_reporter.structure_count == 0:
            self._sync_stage_reporter_from_stage_manifest(s3_reporter, s3_manifest)
        s3_reporter.emit_stage_event(
            "s3_completed",
            total_structures=s3_reporter.structure_count,
            manifest=str(s3_manifest),
            reused=s3_reused,
        )
        if not legacy_checkpoint.reusable("s3", s3_signature, s3_manifest):
            legacy_checkpoint.mark("s3", s3_signature, s3_manifest)
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
        s4_signature = V4Checkpoint.signature({
            "stage": "s4",
            "structures": self._structure_signatures(s4_structures),
            "config": self.config.get("theory", {}).get("s4_high_precision", {}),
        })
        if checkpoint.reusable("s4", "s4", s4_signature, s4_manifest):
            logger.info("[V4] Reusing S4 manifest")
        else:
            s4_manifest = HighLevelEngine(self.config).run(
                s4_structures,
                s4_dir,
                event_callback=ui_reporter.s4_event_callback,
            )
            checkpoint.mark("s4", "s4", s4_signature, s4_manifest)
        if not legacy_checkpoint.reusable("s4", s4_signature, s4_manifest):
            legacy_checkpoint.mark("s4", s4_signature, s4_manifest)
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
            }
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
        payload: Dict[str, Any] = {
            "schema_version": "s2_peb_manifest_v1",
            "stage": "S2",
            "peb_peak": str(ts_seed),
            "intermediate_seed": str(intermediate_seed),
            "forming_bonds": [list(item) for item in forming_bonds],
            "scan_profile": str(profile),
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
    parser.add_argument(
        "--stop-after",
        choices=("s1", "s2", "s3", "s4"),
        default="s4",
        help="Stop after the selected stage; default: s4",
    )
    args = parser.parse_args(argv)
    log_file = args.log_file or args.output / "rph_v4.log"
    color = not args.no_color and sys.stdout.isatty()
    setup_v4_logging(
        log_file=log_file,
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        rich_console=color,
    )
    logger.info("[V4] Pipeline log: %s", log_file)
    s0_record = load_s0_reaction_record(args.csv, args.rx_id)
    V4Orchestrator(args.config, color=color).run(
        s0_record,
        args.output,
        stop_after=args.stop_after,
    )
    return 0
