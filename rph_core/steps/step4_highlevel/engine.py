# pyright: reportMissingImports=false
"""S4 high-precision OPT/SP stage with per-structure failure isolation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List

from rph_core.steps.stage_calculator import StageCalculator
from rph_core.utils.s4_progress import S4ProgressReporter


class HighLevelEngine:
    stage_name = "S4"

    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.stage_config = dict((config.get("theory", {}) or {}).get("s4_high_precision", {}) or {})

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
        calculator = StageCalculator(
            self.config,
            self.stage_config,
            event_callback=reporter.calculator_event,
        )
        results = [self._run_one(structure, output_dir, calculator, reporter) for structure in materialized]
        path = output_dir / "manifest.json"
        path.write_text(
            json.dumps(
                {
                    "schema_version": "s4_high_level_v2",
                    "stage": "S4",
                    "theory": self.stage_config,
                    "summary": self._summary(results),
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
        reporter.start_structure(structure_id)
        try:
            payload.update(calculator.run_structure(structure, target))
            payload["usable_for_ml"] = bool(payload.get("sp_status") == "complete")
        except (OSError, RuntimeError, ValueError) as exc:
            payload["error"] = str(exc)
        payload["s4_status"] = payload["status"]
        payload["s4_usable_for_ml"] = bool(payload["usable_for_ml"])
        reporter.finish_structure(payload)
        return payload

    @staticmethod
    def _summary(results: Iterable[Dict[str, Any]]) -> Dict[str, int]:
        rows = list(results)
        return {
            "total_structures": len(rows),
            "complete": sum(row.get("status") == "complete" for row in rows),
            "ts_frequency_unverified": sum(row.get("status") == "ts_frequency_unverified" for row in rows),
            "degraded": sum(row.get("status") in {"degraded", "opt_failed_sp_complete"} for row in rows),
            "failed": sum(row.get("status") == "failed" for row in rows),
            "usable_for_ml": sum(bool(row.get("usable_for_ml")) for row in rows),
        }
