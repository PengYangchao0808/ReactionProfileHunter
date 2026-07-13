"""S3 low-level OPT/SP stage for V4."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from rph_core.steps.stage_calculator import StageCalculator
from rph_core.utils.stage_progress import StageProgressReporter


class LowLevelEngine:
    stage_name = "S3"

    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.stage_config = dict((config.get("theory", {}) or {}).get("s3_low_level", {}) or {})
        self.progress_reporter: Optional[StageProgressReporter] = None

    def _stage_config(self) -> Dict[str, Any]:
        return dict(self.config)

    def set_progress_reporter(self, reporter: StageProgressReporter) -> None:
        self.progress_reporter = reporter

    def run(self, structures: Iterable[Dict[str, Any]], output_dir: Path) -> Path:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        materialized = list(structures)
        reporter = self.progress_reporter
        manifest: List[Dict[str, Any]] = []
        calculator = StageCalculator(
            self.config,
            self.stage_config,
            event_callback=reporter.calculator_event if reporter is not None else None,
        )
        for structure in materialized:
            manifest.append(self._run_one(structure, output_dir, calculator, reporter))
        path = output_dir / "manifest.json"
        path.write_text(json.dumps({"schema_version": "s3_low_level_v1", "stage": "S3", "structures": manifest}, indent=2, default=str), encoding="utf-8")
        return path

    def _run_one(
        self,
        structure: Dict[str, Any],
        output_dir: Path,
        calculator: StageCalculator,
        reporter: Optional[StageProgressReporter],
    ) -> Dict[str, Any]:
        structure_id = str(structure["id"])
        source = Path(structure["input_xyz"])
        target = output_dir / structure_id
        target.mkdir(parents=True, exist_ok=True)
        payload: Dict[str, Any] = {
            "id": structure_id,
            "kind": structure.get("kind", "minimum"),
            "input_xyz": str(source),
            "status": "degraded",
            "usable_for_ml": False,
        }
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
        if payload["kind"] == "ts" and structure.get("forming_bonds") is not None:
            payload["forming_bonds"] = [
                [int(pair[0]), int(pair[1])]
                for pair in structure.get("forming_bonds", [])
            ]
        if reporter is not None:
            reporter.start_structure(
                structure_id,
                kind=payload["kind"],
                input_xyz=str(source),
                variant_id=payload.get("variant_id"),
                role=payload.get("role"),
                source_stage=payload.get("source_stage"),
            )
        try:
            payload.update(calculator.run_structure(structure, target))
        except (OSError, RuntimeError, ValueError) as exc:
            payload["error"] = str(exc)
        if reporter is not None:
            reporter.finish_structure(
                structure_id,
                str(payload.get("status", "failed")),
                usable_for_ml=bool(payload.get("usable_for_ml", False)),
                error=payload.get("error"),
                opt_status=payload.get("opt_status"),
                sp_status=payload.get("sp_status"),
                frequency_status=payload.get("frequency_status"),
                ts_frequency_valid=payload.get("ts_frequency_valid"),
                ts_mode_displacement_verified=payload.get("ts_mode_displacement_verified"),
                frequency_count_valid=payload.get("frequency_count_valid"),
                mode_displacement_valid=payload.get("mode_displacement_valid"),
                irc_valid=payload.get("irc_valid"),
                ts_quality_summary=payload.get("ts_quality_summary"),
                sp_energy_hartree=payload.get("sp_energy_hartree"),
            )
        return payload
