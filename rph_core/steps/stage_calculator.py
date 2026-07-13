"""Compose exactly one optimization and one single-point job per structure."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Callable, Dict, Optional

from rph_core.utils.qc_jobs import run_frequency, run_optimization, run_single_point
from rph_core.utils.qc_models import QCJobResult, QCJobSpec
from rph_core.utils.ts_quality import analyze_ts_quality


logger = logging.getLogger(__name__)


class StageCalculator:
    def __init__(
        self,
        config: Dict[str, Any],
        theory: Dict[str, Any],
        event_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
    ):
        self.config = config
        self.theory = theory
        self.event_callback = event_callback

    def run_structure(self, structure: Dict[str, Any], output_dir: Path) -> Dict[str, Any]:
        structure_id = str(structure["id"])
        kind = str(structure.get("kind", "minimum"))
        input_value = structure.get("input_xyz") or structure.get("fallback_xyz")
        if not input_value:
            raise ValueError(f"Structure {structure_id} has no XYZ input")
        input_xyz = Path(input_value)
        opt_cfg = dict(self.theory.get("optimization", {}) or {})
        sp_cfg = dict(self.theory.get("single_point", {}) or {})
        opt_route = opt_cfg.get("route_ts" if kind == "ts" else "route_minimum", "OptTS" if kind == "ts" else "Opt")
        charge_value = structure.get("charge")
        multiplicity_value = structure.get("multiplicity")
        charge = int(opt_cfg.get("charge", 0) if charge_value is None else charge_value)
        multiplicity = int(opt_cfg.get("multiplicity", 1) if multiplicity_value is None else multiplicity_value)
        frequency_cfg = dict(opt_cfg.get("frequency", {}) or {})
        frequency_required = kind == "ts" and bool(frequency_cfg.get("enabled_for_ts", False))
        frequency_engine = str(opt_cfg.get("engine", "orca")).lower()
        if frequency_required and frequency_engine not in {"orca", "gaussian"}:
            raise ValueError("V4 TS frequency validation requires an ORCA or Gaussian optimization engine")
        opt_spec = QCJobSpec(
            engine=str(opt_cfg.get("engine", "orca")),
            task="opt_ts" if kind == "ts" else "opt",
            method=str(opt_cfg.get("method", "B97-3c")),
            basis=str(opt_cfg.get("basis", "")),
            aux_basis=str(opt_cfg.get("aux_basis", "")),
            solvent=opt_cfg.get("solvent"),
            solvent_model=opt_cfg.get("solvent_model"),
            route=str(opt_route),
            charge=charge,
            multiplicity=multiplicity,
            nproc=opt_cfg.get("nproc"),
            memory=opt_cfg.get("mem"),
            maxcore=opt_cfg.get("maxcore"),
            route_extras=str(opt_cfg.get("route_extras", "")),
            max_cycles=opt_cfg.get("max_cycles"),
            grid=opt_cfg.get("grid"),
            scf=opt_cfg.get("scf"),
            timeout=opt_cfg.get("timeout"),
        )
        self._emit("optimization_started", structure_id, opt_spec)
        opt = run_optimization(opt_spec, input_xyz, output_dir / "opt", self.config)
        self._emit("optimization_finished", structure_id, opt_spec, opt)
        frequency = None
        frequency_spec = None
        if frequency_required:
            frequency_task = str(frequency_cfg.get("task", "numfreq")).strip().lower()
            frequency_spec = QCJobSpec(
                engine=str(opt_cfg.get("engine", "orca")),
                task=frequency_task,
                method=str(opt_cfg.get("method", "B97-3c")),
                basis=str(opt_cfg.get("basis", "")),
                aux_basis=str(opt_cfg.get("aux_basis", "")),
                solvent=opt_cfg.get("solvent"),
                solvent_model=opt_cfg.get("solvent_model"),
                charge=charge,
                multiplicity=multiplicity,
                nproc=opt_cfg.get("nproc"),
                memory=opt_cfg.get("mem"),
                maxcore=opt_cfg.get("maxcore"),
                route_extras=str(frequency_cfg.get("route_extras", "")),
                grid=opt_cfg.get("grid"),
                scf=opt_cfg.get("scf"),
                timeout=frequency_cfg.get("timeout", opt_cfg.get("timeout")),
            )
        if frequency_spec is not None and opt.status == "complete":
            self._emit("frequency_started", structure_id, frequency_spec)
            frequency = run_frequency(
                frequency_spec,
                opt.output_xyz or input_xyz,
                output_dir / "freq",
                self.config,
            )
            self._emit("frequency_finished", structure_id, frequency_spec, frequency)
        elif frequency_spec is not None:
            self._emit(
                "frequency_skipped",
                structure_id,
                frequency_spec,
                QCJobResult(
                    "skipped",
                    opt.output_xyz or input_xyz,
                    error="Skipped because the prerequisite optimization did not converge",
                ),
            )
        frequencies_cm1 = tuple(frequency.frequencies_cm1 or ()) if frequency is not None else tuple()
        imaginary_frequencies_cm1 = tuple(value for value in frequencies_cm1 if value < 0.0)
        imaginary_cutoff_cm1 = float(frequency_cfg.get("imaginary_cutoff_cm1", -50.0))
        significant_imaginary_cm1 = tuple(
            value for value in frequencies_cm1 if value <= imaginary_cutoff_cm1
        )
        require_exactly_one = bool(frequency_cfg.get("require_exactly_one", True))
        if not frequency_required:
            frequency_status = "not_requested"
            ts_frequency_valid = None
        elif frequency is None or frequency.status != "complete" or not frequencies_cm1:
            frequency_status = "failed"
            ts_frequency_valid = False
        else:
            frequency_status = "complete"
            ts_frequency_valid = (
                len(significant_imaginary_cm1) == 1
                if require_exactly_one
                else bool(significant_imaginary_cm1)
            )
        sp_input = opt.output_xyz or input_xyz
        ts_quality_result = None
        if frequency_required and kind == "ts":
            ts_quality_result = analyze_ts_quality(
                frequency_output=Path(frequency.output_file) if frequency and frequency.output_file else None,
                optimized_xyz=opt.output_xyz,
                forming_bonds=structure.get("forming_bonds"),
                frequencies_cm1=frequencies_cm1,
                imaginary_cutoff_cm1=imaginary_cutoff_cm1,
            )
        sp_spec = QCJobSpec(
            engine=str(sp_cfg.get("engine", "orca")),
            task="sp",
            method=str(sp_cfg.get("method", "r2SCAN-3c")),
            basis=str(sp_cfg.get("basis", "")),
            aux_basis=str(sp_cfg.get("aux_basis", "")),
            solvent=sp_cfg.get("solvent"),
            solvent_model=sp_cfg.get("solvent_model"),
            charge=charge,
            multiplicity=multiplicity,
            nproc=sp_cfg.get("nproc"),
            memory=sp_cfg.get("mem"),
            maxcore=sp_cfg.get("maxcore"),
            route_extras=str(sp_cfg.get("route_extras", "")),
            timeout=sp_cfg.get("timeout"),
        )
        self._emit("single_point_started", structure_id, sp_spec)
        sp = run_single_point(sp_spec, sp_input, output_dir / "sp", self.config)
        self._emit("single_point_finished", structure_id, sp_spec, sp)
        if opt.status == "complete" and sp.status == "complete":
            status = "complete" if ts_frequency_valid is not False else "ts_frequency_unverified"
        else:
            status = "opt_failed_sp_complete" if sp.status == "complete" else "degraded"
        return {
            "id": structure_id,
            "kind": kind,
            "input_xyz": str(input_xyz),
            "opt_xyz": str(opt.output_xyz) if opt.output_xyz else None,
            "sp_input_xyz": str(sp_input),
            "opt_output": str(opt.output_file) if opt.output_file else None,
            "sp_output": str(sp.output_file) if sp.output_file else None,
            "opt_status": opt.status,
            "sp_status": sp.status,
            "frequency_status": frequency_status,
            "frequency_output": str(frequency.output_file) if frequency and frequency.output_file else None,
            "frequencies_cm1": list(frequencies_cm1),
            "imaginary_frequencies_cm1": list(imaginary_frequencies_cm1),
            "significant_imaginary_frequencies_cm1": list(significant_imaginary_cm1),
            "imaginary_cutoff_cm1": imaginary_cutoff_cm1 if frequency_required else None,
            "ts_frequency_valid": ts_frequency_valid,
            "ts_mode_displacement_verified": ts_quality_result["mode_displacement_valid"] if ts_quality_result else (False if frequency_required else None),
            "frequency_count_valid": ts_quality_result["frequency_count_valid"] if ts_quality_result else ts_frequency_valid,
            "mode_displacement_valid": ts_quality_result["mode_displacement_valid"] if ts_quality_result else None,
            "irc_valid": ts_quality_result["irc_valid"] if ts_quality_result else None,
            "ts_quality_summary": ts_quality_result["quality_summary"] if ts_quality_result else None,
            "ts_mode_displacement_details": ts_quality_result["mode_displacement_details"] if ts_quality_result else {},
            "opt_energy_hartree": opt.energy_hartree,
            "sp_energy_hartree": sp.energy_hartree,
            "status": status,
            "usable_for_ml": (
                sp.status == "complete"
                and (kind != "ts" or ts_frequency_valid is True)
            ),
            "error": "; ".join(
                item for item in (opt.error, frequency.error if frequency else None, sp.error) if item
            ) or None,
        }

    def _emit(
        self,
        event: str,
        structure_id: str,
        spec: QCJobSpec,
        result: Optional[QCJobResult] = None,
    ) -> None:
        """Forward a best-effort QC lifecycle event to the owning stage."""

        if self.event_callback is None:
            return
        payload: Dict[str, Any] = {
            "structure_id": structure_id,
            "engine": spec.engine,
            "method": spec.method,
            "basis": spec.basis,
            "solvent": spec.solvent,
            "solvent_model": spec.solvent_model,
        }
        if result is not None:
            payload.update(
                {
                    "status": result.status,
                    "output": str(result.output_file) if result.output_file else None,
                    "energy_hartree": result.energy_hartree,
                    "error": result.error,
                }
            )
        try:
            self.event_callback(event, payload)
        except (OSError, RuntimeError, ValueError, TypeError) as exc:
            logger.warning("Ignoring progress callback failure for %s: %s", event, exc)
