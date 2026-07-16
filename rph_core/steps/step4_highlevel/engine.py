# pyright: reportMissingImports=false
"""S4 high-precision OPT/FREQ/SP stage with per-structure failure isolation."""

from __future__ import annotations

import json
import threading
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List

from rph_core.steps.stage_calculator import StageCalculator
from rph_core.utils.s4_progress import S4ProgressReporter
from rph_core.utils.stage_scheduler import (
    resolve_stage_schedule,
    run_structure_queue,
    worker_config,
    worker_theory,
)


class HighLevelEngine:
    stage_name = "S4"

    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.stage_config = dict((config.get("theory", {}) or {}).get("s4_high_precision", {}) or {})
        for section in ("optimization", "single_point"):
            engine = str(
                (self.stage_config.get(section, {}) or {}).get("engine", "orca")
            ).strip().lower()
            if engine != "orca":
                raise ValueError(
                    f"V4 S4 {section} requires engine=orca, got {engine!r}"
                )

    def _stage_config(self) -> Dict[str, Any]:
        return dict(self.config)

    def run(
        self,
        structures: Iterable[Dict[str, Any]],
        output_dir: Path,
        event_callback: Callable[[str, dict[str, Any]], None] | None = None,
    ) -> Path:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        materialized = list(structures)
        reporter = S4ProgressReporter(output_dir, materialized, event_callback=event_callback)
        schedule = resolve_stage_schedule(self.config, "step4")
        calculator_config = worker_config(self.config, schedule)
        calculator_theory = worker_theory(self.stage_config, schedule)
        reporter_lock = threading.RLock()

        def calculator_event(event: str, payload: Dict[str, Any]) -> None:
            with reporter_lock:
                reporter.calculator_event(event, payload)

        def run_one(structure: Dict[str, Any]) -> Dict[str, Any]:
            calculator = StageCalculator(
                calculator_config,
                calculator_theory,
                event_callback=calculator_event,
            )
            return self._run_one(
                structure,
                output_dir,
                calculator,
                reporter,
                reporter_lock,
            )

        results = run_structure_queue(
            materialized,
            run_one,
            schedule,
            thread_name_prefix="rph-s4",
        )
        summary = self._summary(results)
        stage_status = (
            "complete"
            if summary["complete"] == summary["total_structures"]
            else "incomplete"
        )
        path = output_dir / "manifest.json"
        path.write_text(
            json.dumps(
                {
                    "schema_version": "s4_high_level_v4",
                    "stage": "S4",
                    "status": stage_status,
                    "theory": self.stage_config,
                    "scheduling": schedule.manifest_record(),
                    "summary": summary,
                    "structures": results,
                },
                indent=2,
                default=str,
            ),
            encoding="utf-8",
        )
        reporter.finish_stage(path, results)
        return path

    def _run_one(
        self,
        structure: Dict[str, Any],
        output_dir: Path,
        calculator: StageCalculator,
        reporter: S4ProgressReporter,
        reporter_lock: Any = None,
    ) -> Dict[str, Any]:
        structure_id = str(structure["id"])
        source = Path(structure.get("opt_xyz") or structure.get("fallback_xyz") or structure["input_xyz"])
        target = output_dir / structure_id
        target.mkdir(parents=True, exist_ok=True)
        payload: Dict[str, Any] = {
            "id": structure_id,
            "kind": structure.get("kind", "minimum"),
            "input_xyz": str(source),
            "input_source": "S3_OPT" if structure.get("opt_xyz") else "S3_fallback",
            "source_s3": dict(structure.get("source_s3", {}) or {}),
            "status": "failed",
            "usable_for_ml": False,
        }
        source_s3 = dict(structure.get("source_s3", {}) or {})
        for field in (
            "structure_id",
            "variant_id",
            "role",
            "branch_id",
            "pathway_id",
            "parent_structure_id",
            "source_stage",
            "s1_manifest",
            "s1_ensemble_thermodynamics",
            "s1_thermochemistry_status",
            "ensemble_thermochemistry_correction_hartree",
        ):
            if field in structure:
                payload[field] = structure[field]
            elif field in source_s3:
                payload[field] = source_s3[field]
        if payload["kind"] == "ts" and structure.get("forming_bonds") is not None:
            payload["forming_bonds"] = [
                [int(pair[0]), int(pair[1])]
                for pair in structure.get("forming_bonds", [])
            ]
        lock = reporter_lock or threading.RLock()
        with lock:
            reporter.start_structure(structure_id)
        try:
            payload.update(calculator.run_structure(structure, target))
            correction = payload.get("ensemble_thermochemistry_correction_hartree")
            if payload.get("sp_energy_hartree") is not None and correction is not None:
                payload["ensemble_corrected_sp_free_energy_hartree"] = (
                    float(payload["sp_energy_hartree"]) + float(correction)
                )
                payload["ensemble_free_energy_formula"] = (
                    "E_S4_SP + G_RRHO(reference)_S1 + G_conf_rel_S1"
                )
                payload["ensemble_free_energy_status"] = "complete"
            elif payload.get("sp_energy_hartree") is not None:
                payload["ensemble_corrected_sp_free_energy_hartree"] = None
                payload["ensemble_free_energy_status"] = (
                    "thermochemistry_incomplete"
                    if payload.get("s1_thermochemistry_status") == "incomplete"
                    else "not_available"
                )
            payload["usable_for_ml"] = bool(payload.get("usable_for_ml", False))
        except (OSError, RuntimeError, ValueError) as exc:
            payload["error"] = str(exc)
        payload["s4_status"] = payload["status"]
        payload["s4_usable_for_ml"] = bool(payload["usable_for_ml"])
        with lock:
            reporter.finish_structure(payload)
        return payload

    @staticmethod
    def _summary(results: Iterable[Dict[str, Any]]) -> Dict[str, int]:
        rows = list(results)
        return {
            "total_structures": len(rows),
            "complete": sum(row.get("status") == "complete" for row in rows),
            "ts_frequency_unverified": sum(row.get("status") == "ts_frequency_unverified" for row in rows),
            "minimum_frequency_unverified": sum(row.get("status") == "minimum_frequency_unverified" for row in rows),
            "degraded": sum(row.get("status") in {"degraded", "opt_failed_sp_complete", "minimum_frequency_unverified"} for row in rows),
            "failed": sum(row.get("status") == "failed" for row in rows),
            "usable_for_ml": sum(bool(row.get("usable_for_ml")) for row in rows),
        }
