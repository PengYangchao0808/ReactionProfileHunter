"""Compose one optimization, frequency and single-point sequence per structure."""

from __future__ import annotations

import logging
from dataclasses import replace
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Tuple

from rph_core.utils.file_io import read_xyz
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

    @staticmethod
    def _bond_constraints_from_seed(
        input_xyz: Path,
        forming_bonds: Any,
    ) -> Tuple[
        Tuple[Tuple[int, int, Optional[float]], ...],
        List[Dict[str, Any]],
    ]:
        """Freeze forming bonds at their S2 seed distances during DFT warm-up."""

        coordinates, _ = read_xyz(input_xyz)
        atom_count = len(coordinates)
        constraints: List[Tuple[int, int, Optional[float]]] = []
        details: List[Dict[str, Any]] = []
        for raw_pair in forming_bonds or ():
            if len(raw_pair) < 2:
                raise ValueError(f"Invalid forming-bond record: {raw_pair!r}")
            atom_i, atom_j = int(raw_pair[0]), int(raw_pair[1])
            if atom_i < 0 or atom_j < 0 or atom_i >= atom_count or atom_j >= atom_count:
                raise ValueError(
                    f"Forming bond {(atom_i, atom_j)} is outside XYZ atom range 0..{atom_count - 1}"
                )
            if atom_i == atom_j:
                raise ValueError(f"Forming bond cannot reference one atom twice: {(atom_i, atom_j)}")
            delta = coordinates[atom_i] - coordinates[atom_j]
            distance = float(sum(float(value) ** 2 for value in delta) ** 0.5)
            constraints.append((atom_i, atom_j, distance))
            details.append(
                {
                    "atoms": [atom_i, atom_j],
                    "target_distance_angstrom": distance,
                }
            )
        return tuple(constraints), details

    @staticmethod
    def _warmup_applies(structure: Dict[str, Any], warmup_cfg: Dict[str, Any]) -> bool:
        if not bool(warmup_cfg.get("enabled", False)):
            return False
        if str(structure.get("source_stage", "")).upper() != "S2":
            return False
        role = str(structure.get("role") or structure.get("kind", "minimum")).lower()
        roles = {str(value).lower() for value in warmup_cfg.get("roles", ("intermediate", "ts"))}
        return role in roles and bool(structure.get("forming_bonds"))

    @staticmethod
    def _warmup_max_cycles(structure: Dict[str, Any], warmup_cfg: Dict[str, Any]) -> int:
        """Resolve a role-specific warm-up budget while accepting the legacy scalar form."""

        configured = warmup_cfg.get("max_cycles", 8)
        if isinstance(configured, Mapping):
            role = str(structure.get("role") or structure.get("kind", "minimum")).lower()
            configured = configured.get(role, configured.get("default", 8))
        cycles = int(configured)
        if cycles <= 0:
            raise ValueError("S3 warm-up max_cycles must be a positive integer")
        return cycles

    @staticmethod
    def _warmup_route_extras(opt_spec: QCJobSpec, warmup_cfg: Dict[str, Any]) -> Tuple[str, str]:
        """Add an optimization convergence keyword only to the constrained warm-up."""

        convergence = str(warmup_cfg.get("convergence", "normal")).strip().lower()
        keywords = {
            "normal": None,
            "default": None,
            "loose": "LooseOpt",
        }
        if convergence not in keywords:
            raise ValueError(
                "S3 warm-up convergence must be one of: normal, default, loose"
            )
        tokens = str(opt_spec.route_extras or "").split()
        keyword = keywords[convergence]
        if keyword:
            optimization_convergence_keywords = {
                "looseopt",
                "normalopt",
                "tightopt",
                "verytightopt",
            }
            tokens = [
                token
                for token in tokens
                if token.lower() not in optimization_convergence_keywords
            ]
            tokens.append(keyword)
        return " ".join(tokens), convergence

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
        is_ts = kind == "ts"
        frequency_required = bool(
            frequency_cfg.get("enabled_for_ts", False)
            if is_ts
            else frequency_cfg.get("enabled_for_minima", False)
        )
        frequency_engine = str(opt_cfg.get("engine", "orca")).lower()
        if frequency_required and frequency_engine not in {"orca", "gaussian"}:
            raise ValueError("V4 frequency validation requires an ORCA or Gaussian optimization engine")
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
        warmup_cfg = dict(opt_cfg.get("warmup", {}) or {})
        warmup_status = "not_requested"
        warmup_used = False
        warmup_error = None
        warmup_output = None
        warmup_xyz = None
        warmup_energy_hartree = None
        warmup_constraints: List[Dict[str, Any]] = []
        warmup_convergence = None
        warmup_max_cycles = None
        optimization_input_xyz = input_xyz
        if self._warmup_applies(structure, warmup_cfg):
            if opt_spec.engine.lower() != "orca":
                warmup_status = "skipped"
                warmup_error = "S3 constrained warm-up currently requires ORCA"
            else:
                bond_constraints, warmup_constraints = self._bond_constraints_from_seed(
                    input_xyz,
                    structure.get("forming_bonds"),
                )
                warmup_max_cycles = self._warmup_max_cycles(structure, warmup_cfg)
                warmup_route_extras, warmup_convergence = self._warmup_route_extras(
                    opt_spec,
                    warmup_cfg,
                )
                warmup_spec = replace(
                    opt_spec,
                    task="opt",
                    route=str(opt_cfg.get("route_minimum", "Opt")),
                    route_extras=warmup_route_extras,
                    max_cycles=warmup_max_cycles,
                    bond_constraints=bond_constraints,
                    allow_unconverged_geometry=bool(
                        warmup_cfg.get("accept_partial_geometry", True)
                    ),
                )
                self._emit("warmup_started", structure_id, warmup_spec)
                warmup = run_optimization(
                    warmup_spec,
                    input_xyz,
                    output_dir / "warmup",
                    self.config,
                )
                self._emit("warmup_finished", structure_id, warmup_spec, warmup)
                warmup_status = warmup.status
                warmup_error = warmup.error
                warmup_output = str(warmup.output_file) if warmup.output_file else None
                warmup_xyz = str(warmup.output_xyz) if warmup.output_xyz else None
                warmup_energy_hartree = warmup.energy_hartree
                if warmup.output_xyz is not None and warmup.status in {"complete", "partial"}:
                    optimization_input_xyz = warmup.output_xyz
                    warmup_used = True
        elif bool(warmup_cfg.get("enabled", False)):
            warmup_status = "not_applicable"
        self._emit("optimization_started", structure_id, opt_spec)
        opt = run_optimization(
            opt_spec,
            optimization_input_xyz,
            output_dir / "opt",
            self.config,
        )
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
                opt.output_xyz or optimization_input_xyz,
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
                    opt.output_xyz or optimization_input_xyz,
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
        minimum_frequency_valid = None
        if not frequency_required:
            frequency_status = "not_requested"
            ts_frequency_valid = None
        elif frequency is None or frequency.status != "complete" or not frequencies_cm1:
            frequency_status = "failed"
            ts_frequency_valid = False if is_ts else None
            minimum_frequency_valid = False if not is_ts else None
        else:
            frequency_status = "complete"
            if is_ts:
                ts_frequency_valid = (
                    len(significant_imaginary_cm1) == 1
                    if require_exactly_one
                    else bool(significant_imaginary_cm1)
                )
            else:
                ts_frequency_valid = None
                minimum_frequency_valid = len(significant_imaginary_cm1) == 0
        sp_input = opt.output_xyz or optimization_input_xyz
        ts_quality_result = None
        if frequency_required and is_ts:
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
        frequency_energy_hartree = frequency.energy_hartree if frequency else None
        enthalpy_correction_hartree = (
            float(frequency.enthalpy_hartree) - float(frequency_energy_hartree)
            if frequency
            and frequency.enthalpy_hartree is not None
            and frequency_energy_hartree is not None
            else None
        )
        thermal_energy_correction_hartree = (
            float(frequency.thermal_energy_hartree) - float(frequency_energy_hartree)
            if frequency
            and frequency.thermal_energy_hartree is not None
            and frequency_energy_hartree is not None
            else None
        )
        composite_enthalpy_hartree = (
            float(sp.energy_hartree) + enthalpy_correction_hartree
            if sp.energy_hartree is not None and enthalpy_correction_hartree is not None
            else None
        )
        composite_gibbs_free_energy_hartree = (
            float(sp.energy_hartree) + float(frequency.gibbs_correction_hartree)
            if frequency
            and sp.energy_hartree is not None
            and frequency.gibbs_correction_hartree is not None
            else None
        )
        frequency_valid_for_role = (
            ts_frequency_valid if is_ts else minimum_frequency_valid
        )
        if opt.status == "complete" and sp.status == "complete":
            if not frequency_required or frequency_valid_for_role is True:
                status = "complete"
            else:
                status = "ts_frequency_unverified" if is_ts else "minimum_frequency_unverified"
        else:
            status = "opt_failed_sp_complete" if sp.status == "complete" else "degraded"
        return {
            "id": structure_id,
            "kind": kind,
            "input_xyz": str(input_xyz),
            "optimization_input_xyz": str(optimization_input_xyz),
            "warmup_status": warmup_status,
            "warmup_used": warmup_used,
            "warmup_xyz": warmup_xyz,
            "warmup_output": warmup_output,
            "warmup_energy_hartree": warmup_energy_hartree,
            "warmup_constraints": warmup_constraints,
            "warmup_convergence": warmup_convergence,
            "warmup_max_cycles": warmup_max_cycles,
            "warmup_error": warmup_error,
            "opt_xyz": str(opt.output_xyz) if opt.output_xyz else None,
            "sp_input_xyz": str(sp_input),
            "opt_output": str(opt.output_file) if opt.output_file else None,
            "sp_output": str(sp.output_file) if sp.output_file else None,
            "opt_status": opt.status,
            "sp_status": sp.status,
            "frequency_status": frequency_status,
            "frequency_output": str(frequency.output_file) if frequency and frequency.output_file else None,
            "frequency_energy_hartree": frequency_energy_hartree,
            "zero_point_energy_hartree": frequency.zero_point_energy_hartree if frequency else None,
            "thermal_energy_hartree": frequency.thermal_energy_hartree if frequency else None,
            "enthalpy_hartree": frequency.enthalpy_hartree if frequency else None,
            "gibbs_free_energy_hartree": frequency.gibbs_free_energy_hartree if frequency else None,
            "gibbs_correction_hartree": frequency.gibbs_correction_hartree if frequency else None,
            "thermal_energy_correction_hartree": thermal_energy_correction_hartree,
            "enthalpy_correction_hartree": enthalpy_correction_hartree,
            "composite_enthalpy_hartree": composite_enthalpy_hartree,
            "composite_gibbs_free_energy_hartree": composite_gibbs_free_energy_hartree,
            "frequencies_cm1": list(frequencies_cm1),
            "imaginary_frequencies_cm1": list(imaginary_frequencies_cm1),
            "significant_imaginary_frequencies_cm1": list(significant_imaginary_cm1),
            "imaginary_cutoff_cm1": imaginary_cutoff_cm1 if frequency_required else None,
            "ts_frequency_valid": ts_frequency_valid,
            "minimum_frequency_valid": minimum_frequency_valid,
            "ts_mode_displacement_verified": ts_quality_result["mode_displacement_valid"] if ts_quality_result else (False if frequency_required and is_ts else None),
            "frequency_count_valid": ts_quality_result["frequency_count_valid"] if ts_quality_result else frequency_valid_for_role,
            "mode_displacement_valid": ts_quality_result["mode_displacement_valid"] if ts_quality_result else None,
            "irc_valid": ts_quality_result["irc_valid"] if ts_quality_result else None,
            "ts_quality_summary": ts_quality_result["quality_summary"] if ts_quality_result else None,
            "ts_mode_displacement_details": ts_quality_result["mode_displacement_details"] if ts_quality_result else {},
            "opt_energy_hartree": opt.energy_hartree,
            "sp_energy_hartree": sp.energy_hartree,
            "status": status,
            "usable_for_ml": (
                opt.status == "complete"
                and sp.status == "complete"
                and (not frequency_required or frequency_valid_for_role is True)
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
