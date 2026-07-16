"""One-shot OPT/OPTTS/FREQ/SP jobs.

This module deliberately contains no rescue, NBO, checkpoint or stage policy
logic.  All external execution remains behind existing QC interfaces.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from rph_core.utils.file_io import read_xyz, write_xyz
from rph_core.utils.keyword_translator import KeywordTranslator
from rph_core.utils.qc_interface import GaussianInterface
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.utils.qc_models import QCJobResult, QCJobSpec


def _resources(config: Dict[str, Any], spec: QCJobSpec) -> tuple[int, str]:
    resources = dict(config.get("resources", {}) or {})
    return int(spec.nproc or resources.get("nproc", 1)), str(spec.memory or resources.get("mem", "1GB"))


def _orca_route_extras(spec: QCJobSpec) -> str:
    """Return normalized ORCA controls, excluding the task keyword."""

    tokens = spec.route_extras.split()
    tokens.extend(str(spec.route or "").replace("=", " ").split())
    if spec.grid:
        tokens.append(str(spec.grid).strip())
    if spec.scf:
        tokens.append(str(spec.scf).strip())
    ignored = {"opt", "optts", "ts", "tightscf"}
    normalized = []
    seen = set()
    for token in tokens:
        lowered = token.lower()
        if not token or lowered in ignored or lowered in seen:
            continue
        normalized.append(token)
        seen.add(lowered)
    return " ".join(normalized)


def _frequency_route_extras(spec: QCJobSpec) -> str:
    """Add the configured ORCA Hessian keyword to a frequency job."""

    task_name = spec.task.strip().lower()
    keyword_by_task = {"freq": "Freq", "numfreq": "NumFreq"}
    try:
        hessian_keyword = keyword_by_task[task_name]
    except KeyError as exc:
        raise ValueError(f"Unsupported ORCA frequency task: {spec.task}") from exc
    tokens = _orca_route_extras(spec).split()
    if hessian_keyword.lower() not in {token.lower() for token in tokens}:
        tokens.append(hessian_keyword)
    return " ".join(tokens)


def _orca_geom_block(spec: QCJobSpec) -> str:
    if spec.max_cycles is None and not spec.bond_constraints:
        return ""
    lines = ["%geom"]
    if spec.max_cycles is not None:
        lines.append(f"  MaxIter {int(spec.max_cycles)}")
    if spec.bond_constraints:
        lines.append("  Constraints")
        for atom_i, atom_j, target in spec.bond_constraints:
            if int(atom_i) < 0 or int(atom_j) < 0 or int(atom_i) == int(atom_j):
                raise ValueError(f"Invalid ORCA bond constraint: {(atom_i, atom_j, target)!r}")
            if target is None:
                lines.append(f"    {{ B {int(atom_i)} {int(atom_j)} C }}")
            else:
                target_value = float(target)
                if target_value <= 0.0:
                    raise ValueError(f"Invalid ORCA bond-constraint distance: {target_value}")
                lines.append(
                    f"    {{ B {int(atom_i)} {int(atom_j)} {target_value:.8f} C }}"
                )
        lines.append("  end")
    lines.append("end")
    return "\n".join(lines)


def _gaussian_route(spec: QCJobSpec, require_frequency: bool = False) -> str:
    basis = KeywordTranslator.to_gaussian_basis(spec.basis) if spec.basis else ""
    route = f"{spec.method}/{basis} {spec.route}".strip() if basis else f"{spec.method} {spec.route}".strip()
    if require_frequency and "freq" not in route.lower():
        route += " Freq"
    if spec.grid and "int=" not in route.lower():
        route += f" Int={spec.grid}"
    if spec.scf and "scf=" not in route.lower():
        route += f" SCF={spec.scf}"
    if spec.max_cycles is not None and "maxcycle" not in route.lower():
        route += f" MaxCycle={int(spec.max_cycles)}"
    if spec.route_extras:
        route += f" {spec.route_extras.strip()}"
    if spec.solvent and "scrf=" not in route.lower():
        solvent_model = str(spec.solvent_model or "").strip()
        if solvent_model:
            route += f" SCRF=({solvent_model},Solvent={spec.solvent})"
        else:
            route += f" SCRF=(Solvent={spec.solvent})"
    return route


def _write_coordinates(input_xyz: Path, coordinates: Any, output_xyz: Path) -> Optional[Path]:
    if coordinates is None:
        return None
    try:
        _, symbols = read_xyz(input_xyz)
        output_xyz.parent.mkdir(parents=True, exist_ok=True)
        write_xyz(output_xyz, coordinates, symbols, title="optimized")
        return output_xyz
    except (OSError, TypeError, ValueError, IndexError):
        return None


def _frequencies_cm1(result: Any) -> Optional[Tuple[float, ...]]:
    """Normalize interface-native frequency arrays for manifest serialization."""

    values = getattr(result, "frequencies", None)
    if values is None:
        return None
    try:
        return tuple(float(value) for value in values)
    except (TypeError, ValueError):
        return None


def run_optimization(spec: QCJobSpec, input_xyz: Path, output_dir: Path, config: Dict[str, Any]) -> QCJobResult:
    input_xyz = Path(input_xyz).resolve()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    nproc, memory = _resources(config, spec)
    try:
        if spec.engine.lower() == "orca":
            interface = ORCAInterface(
                method=spec.method,
                basis=spec.basis,
                aux_basis=spec.aux_basis,
                nprocs=nproc,
                maxcore=spec.maxcore,
                solvent=spec.solvent or "",
                solvent_model=spec.solvent_model,
                route_extras=_orca_route_extras(spec),
                config=config,
            )
            task_name = spec.task.lower()
            task_type = "ts_freq" if task_name in {"ts_freq", "opt_ts_freq", "ts-freq"} else (
                "ts" if task_name in {"opt_ts", "ts", "opt-ts"} else "opt"
            )
            result = interface.optimize_with_task(
                input_xyz,
                output_dir,
                task_type=task_type,
                charge=spec.charge,
                spin=spec.multiplicity,
                geom_block=_orca_geom_block(spec),
                timeout=spec.timeout,
            )
        elif spec.engine.lower() == "gaussian":
            route = _gaussian_route(spec)
            interface = GaussianInterface(charge=spec.charge, multiplicity=spec.multiplicity, nprocshared=nproc, mem=memory, config=config)
            result = interface.optimize(input_xyz, output_dir, route=route, timeout=spec.timeout, charge=spec.charge, spin=spec.multiplicity)
        else:
            raise ValueError(f"Unsupported optimization engine: {spec.engine}")
        output_file = getattr(result, "output_file", None)
        converged = bool(getattr(result, "converged", False))
        may_export_geometry = converged or bool(spec.allow_unconverged_geometry)
        output_xyz = (
            _write_coordinates(
                input_xyz,
                getattr(result, "coordinates", None),
                output_dir / "opt.xyz",
            )
            if may_export_geometry
            else None
        )
        status = "complete" if converged else ("partial" if output_xyz is not None else "failed")
        return QCJobResult(
            status,
            input_xyz,
            output_xyz,
            output_file,
            getattr(result, "energy", None),
            getattr(result, "error_message", None),
            frequencies_cm1=_frequencies_cm1(result),
        )
    except (OSError, RuntimeError, ValueError) as exc:
        return QCJobResult("failed", input_xyz, error=str(exc))


def run_single_point(spec: QCJobSpec, input_xyz: Path, output_dir: Path, config: Dict[str, Any]) -> QCJobResult:
    input_xyz = Path(input_xyz).resolve()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    nproc, _ = _resources(config, spec)
    try:
        if spec.engine.lower() != "orca":
            raise ValueError("V4 single-point jobs currently require ORCA")
        interface = ORCAInterface(
            method=spec.method,
            basis=spec.basis,
            aux_basis=spec.aux_basis,
            nprocs=nproc,
            maxcore=spec.maxcore,
            solvent=spec.solvent or "",
            solvent_model=spec.solvent_model,
            route_extras=spec.route_extras,
            config=config,
        )
        result = interface.single_point(input_xyz, output_dir, timeout=spec.timeout, charge=spec.charge, spin=spec.multiplicity)
        return QCJobResult("complete" if getattr(result, "converged", False) else "failed", input_xyz, output_file=getattr(result, "output_file", None), energy_hartree=getattr(result, "energy", None), error=getattr(result, "error_message", None))
    except (OSError, RuntimeError, ValueError) as exc:
        return QCJobResult("failed", input_xyz, error=str(exc))


def run_frequency(spec: QCJobSpec, input_xyz: Path, output_dir: Path, config: Dict[str, Any]) -> QCJobResult:
    """Run an ORCA or Gaussian frequency calculation on one geometry."""

    input_xyz = Path(input_xyz).resolve()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    nproc, _ = _resources(config, spec)
    try:
        engine = spec.engine.lower()
        if engine == "orca":
            interface = ORCAInterface(
                method=spec.method,
                basis=spec.basis,
                aux_basis=spec.aux_basis,
                nprocs=nproc,
                maxcore=spec.maxcore,
                solvent=spec.solvent or "",
                solvent_model=spec.solvent_model,
                route_extras=_frequency_route_extras(spec),
                config=config,
            )
            result = interface.single_point(
                input_xyz,
                output_dir,
                timeout=spec.timeout,
                charge=spec.charge,
                spin=spec.multiplicity,
            )
            if getattr(result, "frequencies", None) is None:
                output_file = getattr(result, "output_file", None)
                parser = getattr(interface, "extract_frequencies_from_output", None)
                if output_file is not None and callable(parser):
                    result.frequencies = parser(Path(output_file))
        elif engine == "gaussian":
            _, memory = _resources(config, spec)
            interface = GaussianInterface(
                charge=spec.charge,
                multiplicity=spec.multiplicity,
                nprocshared=nproc,
                mem=memory,
                config=config,
            )
            result = interface.optimize(
                input_xyz,
                output_dir,
                route=_gaussian_route(spec, require_frequency=True),
                timeout=spec.timeout,
                charge=spec.charge,
                spin=spec.multiplicity,
            )
        else:
            raise ValueError(f"Unsupported V4 frequency engine: {spec.engine}")
        frequencies_cm1 = _frequencies_cm1(result)
        if getattr(result, "converged", False) and frequencies_cm1:
            return QCJobResult(
                "complete",
                input_xyz,
                output_file=getattr(result, "output_file", None),
                energy_hartree=getattr(result, "energy", None),
                frequencies_cm1=frequencies_cm1,
            )
        return QCJobResult(
            "failed",
            input_xyz,
            output_file=getattr(result, "output_file", None),
            energy_hartree=getattr(result, "energy", None),
            frequencies_cm1=frequencies_cm1,
            error=getattr(result, "error_message", None) or "ORCA frequency output contains no parsed frequencies",
        )
    except (OSError, RuntimeError, ValueError) as exc:
        return QCJobResult("failed", input_xyz, error=str(exc))
