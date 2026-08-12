"""One-shot OPT/OPTTS/FREQ/SP jobs.

This module deliberately contains no rescue, NBO, checkpoint or stage policy
logic.  All external execution remains behind existing QC interfaces.
"""

from __future__ import annotations

from dataclasses import replace
import re
import subprocess
from pathlib import Path
from typing import Any, Callable, Dict, Mapping, Optional, Tuple

from rph_core.utils.file_io import read_xyz, write_xyz
from rph_core.utils.keyword_translator import KeywordTranslator
from rph_core.utils.qc_interface import GaussianInterface
from rph_core.utils.orca_interface import ORCAInterface, extract_orca_version_from_output
from rph_core.utils.qc_models import (
    IRCJobSpec,
    NEBJobSpec,
    QCJobResult,
    QCJobSpec,
    SurfaceScanSpec,
)


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
    has_hessian_controls = (
        spec.initial_hessian in {"calculate", "read"}
        or spec.ts_mode is not None
        or spec.recalc_hessian is not None
        or spec.trust is not None
    )
    if spec.max_cycles is None and not has_hessian_controls and not spec.bond_constraints:
        return ""
    lines = ["%geom"]
    if spec.max_cycles is not None:
        lines.append(f"  MaxIter {int(spec.max_cycles)}")
    if spec.initial_hessian == "calculate":
        lines.append("  Calc_Hess true")
    elif spec.initial_hessian == "read":
        lines.append("  InHess Read")
        if spec.hessian_filename:
            lines.append(f'  InHessName "{spec.hessian_filename}"')
    if spec.ts_mode is not None:
        lines.append(f"  TS_Mode {{M {int(spec.ts_mode)}}}")
    if spec.recalc_hessian is not None:
        lines.append(f"  Recalc_Hess {int(spec.recalc_hessian)}")
    if spec.trust is not None:
        lines.append(f"  Trust {float(spec.trust)}")
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


def _merge_orca_route_extras(*chunks: Any) -> str:
    tokens = []
    seen = set()
    for chunk in chunks:
        for token in str(chunk or "").replace("=", " ").split():
            lowered = token.lower()
            if not token or lowered in seen:
                continue
            tokens.append(token)
            seen.add(lowered)
    return " ".join(tokens)


def _s3_optimization_config(config: Mapping[str, Any]) -> Mapping[str, Any]:
    theory = config.get("theory", {}) if isinstance(config, Mapping) else {}
    if not isinstance(theory, Mapping):
        return {}
    s3_low_level = theory.get("s3_low_level", {})
    if not isinstance(s3_low_level, Mapping):
        return {}
    optimization = s3_low_level.get("optimization", {})
    return optimization if isinstance(optimization, Mapping) else {}


def _resolve_irc_method(config: Mapping[str, Any]) -> str:
    method = str(_s3_optimization_config(config).get("method", "") or "").strip()
    return method or "B97-3c"


def _resolve_orca_interface_for_irc(
    config: Mapping[str, Any],
    *,
    method: str,
    basis: str,
) -> ORCAInterface:
    optimization = _s3_optimization_config(config)
    resources = dict(config.get("resources", {}) or {}) if isinstance(config, Mapping) else {}
    route_extras = _merge_orca_route_extras(
        optimization.get("route_extras", ""),
        optimization.get("grid", ""),
        optimization.get("scf", ""),
    )
    return ORCAInterface(
        method=method,
        basis=basis,
        aux_basis=str(optimization.get("aux_basis", "") or ""),
        nprocs=int(resources.get("nproc", 1)),
        solvent=str(optimization.get("solvent", "") or ""),
        solvent_model=(
            str(optimization.get("solvent_model", "") or "").strip() or None
        ),
        route_extras=route_extras,
        scf_maxiter=(
            int(optimization["scf_maxiter"])
            if str(optimization.get("scf_maxiter", "")).strip()
            else None
        ),
        config=dict(config),
    )


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


def _orca_thermochemistry(output_file: Optional[Path]) -> Dict[str, Optional[float]]:
    """Extract the final ORCA thermochemistry summary from a FREQ output."""

    fields = {
        "zero_point_energy_hartree": r"Zero point energy",
        "thermal_energy_hartree": r"Total thermal energy",
        "enthalpy_hartree": r"Total Enthalpy",
        "gibbs_free_energy_hartree": r"Final Gibbs free energy",
        "gibbs_correction_hartree": r"G-E\(el\)",
    }
    values: Dict[str, Optional[float]] = {key: None for key in fields}
    if output_file is None:
        return values
    try:
        content = Path(output_file).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return values
    for key, label in fields.items():
        matches = re.findall(
            rf"{label}\s*(?:\.\.\.)?\s*([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s+Eh",
            content,
            flags=re.IGNORECASE,
        )
        if matches:
            values[key] = float(matches[-1])
    return values


def run_optimization(
    spec: QCJobSpec,
    input_xyz: Path,
    output_dir: Path,
    config: Dict[str, Any],
    subprocess_callback: Callable[[subprocess.Popen[Any]], None] | None = None,
) -> QCJobResult:
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
                scf_maxiter=spec.scf_maxiter,
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
                max_cycles_opt=spec.max_cycles,
                timeout=spec.timeout,
                subprocess_callback=subprocess_callback,
            )
        elif spec.engine.lower() == "gaussian":
            route = _gaussian_route(spec)
            interface = GaussianInterface(charge=spec.charge, multiplicity=spec.multiplicity, nprocshared=nproc, mem=memory, config=config)
            result = interface.optimize(input_xyz, output_dir, route=route, timeout=spec.timeout, charge=spec.charge, spin=spec.multiplicity)
        else:
            raise ValueError(f"Unsupported optimization engine: {spec.engine}")
        output_file = getattr(result, "output_file", None)
        run_metadata = getattr(interface, "_last_orca_run_metadata", None) if spec.engine.lower() == "orca" else None
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
            failure_type=getattr(result, "failure_type", None),
            failure_evidence_lines=getattr(result, "failure_evidence_lines", None),
            optimization_converged=getattr(result, "optimization_converged", None),
            stop_reason=getattr(result, "stop_reason", None),
            returncode=getattr(run_metadata, "returncode", None),
            stderr_text=getattr(run_metadata, "stderr_text", None),
            timed_out=bool(getattr(run_metadata, "timed_out", False)),
            orca_version=extract_orca_version_from_output(output_file) if output_file else None,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        return QCJobResult("failed", input_xyz, error=str(exc))


def run_single_point(
    spec: QCJobSpec,
    input_xyz: Path,
    output_dir: Path,
    config: Dict[str, Any],
    subprocess_callback: Callable[[subprocess.Popen[Any]], None] | None = None,
) -> QCJobResult:
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
            scf_maxiter=spec.scf_maxiter,
            config=config,
        )
        result = interface.single_point(
            input_xyz,
            output_dir,
            timeout=spec.timeout,
            charge=spec.charge,
            spin=spec.multiplicity,
            subprocess_callback=subprocess_callback,
        )
        run_metadata = getattr(interface, "_last_orca_run_metadata", None)
        return QCJobResult(
            "complete" if getattr(result, "converged", False) else "failed",
            input_xyz,
            output_file=getattr(result, "output_file", None),
            energy_hartree=getattr(result, "energy", None),
            error=getattr(result, "error_message", None),
            failure_type=getattr(result, "failure_type", None),
            failure_evidence_lines=getattr(result, "failure_evidence_lines", None),
            returncode=getattr(run_metadata, "returncode", None),
            stderr_text=getattr(run_metadata, "stderr_text", None),
            timed_out=bool(getattr(run_metadata, "timed_out", False)),
            orca_version=extract_orca_version_from_output(getattr(result, "output_file", None)),
        )
    except (OSError, RuntimeError, ValueError) as exc:
        return QCJobResult("failed", input_xyz, error=str(exc))


def run_frequency(
    spec: QCJobSpec,
    input_xyz: Path,
    output_dir: Path,
    config: Dict[str, Any],
    subprocess_callback: Callable[[subprocess.Popen[Any]], None] | None = None,
) -> QCJobResult:
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
                scf_maxiter=spec.scf_maxiter,
                config=config,
            )
            result = interface.single_point(
                input_xyz,
                output_dir,
                timeout=spec.timeout,
                charge=spec.charge,
                spin=spec.multiplicity,
                subprocess_callback=subprocess_callback,
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
        output_file = getattr(result, "output_file", None)
        run_metadata = getattr(interface, "_last_orca_run_metadata", None) if engine == "orca" else None
        thermochemistry = (
            _orca_thermochemistry(Path(output_file) if output_file else None)
            if engine == "orca"
            else {}
        )
        if getattr(result, "converged", False) and frequencies_cm1:
            return QCJobResult(
                status="complete",
                input_xyz=input_xyz,
                output_file=output_file,
                energy_hartree=getattr(result, "energy", None),
                frequencies_cm1=frequencies_cm1,
                zero_point_energy_hartree=thermochemistry.get("zero_point_energy_hartree"),
                thermal_energy_hartree=thermochemistry.get("thermal_energy_hartree"),
                enthalpy_hartree=thermochemistry.get("enthalpy_hartree"),
                gibbs_free_energy_hartree=thermochemistry.get("gibbs_free_energy_hartree"),
                gibbs_correction_hartree=thermochemistry.get("gibbs_correction_hartree"),
                returncode=getattr(run_metadata, "returncode", None),
                stderr_text=getattr(run_metadata, "stderr_text", None),
                timed_out=bool(getattr(run_metadata, "timed_out", False)),
                orca_version=extract_orca_version_from_output(output_file) if output_file else None,
            )
        return QCJobResult(
            status="failed",
            input_xyz=input_xyz,
            output_file=output_file,
            energy_hartree=getattr(result, "energy", None),
            frequencies_cm1=frequencies_cm1,
            error=getattr(result, "error_message", None) or "Frequency output contains no parsed frequencies",
            failure_type=getattr(result, "failure_type", None),
            failure_evidence_lines=getattr(result, "failure_evidence_lines", None),
            zero_point_energy_hartree=thermochemistry.get("zero_point_energy_hartree"),
            thermal_energy_hartree=thermochemistry.get("thermal_energy_hartree"),
            enthalpy_hartree=thermochemistry.get("enthalpy_hartree"),
            gibbs_free_energy_hartree=thermochemistry.get("gibbs_free_energy_hartree"),
            gibbs_correction_hartree=thermochemistry.get("gibbs_correction_hartree"),
            returncode=getattr(run_metadata, "returncode", None),
            stderr_text=getattr(run_metadata, "stderr_text", None),
            timed_out=bool(getattr(run_metadata, "timed_out", False)),
            orca_version=extract_orca_version_from_output(output_file) if output_file else None,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        return QCJobResult("failed", input_xyz, error=str(exc))


def run_irc(
    spec: IRCJobSpec,
    input_xyz: Path,
    output_dir: Path,
    config: Mapping[str, Any],
    subprocess_callback: Optional[Callable[[Any], None]] = None,
    charge: int = 0,
    spin: int = 1,
) -> QCJobResult:
    """Execute an ORCA IRC calculation for TS endpoint discovery."""

    input_xyz = Path(input_xyz).resolve()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    optimization = _s3_optimization_config(config)
    engine = str(optimization.get("engine", getattr(spec, "engine", "orca")) or "orca").strip().lower()
    if engine != "orca":
        return QCJobResult(
            status="failed",
            input_xyz=input_xyz,
            error=f"Unsupported IRC engine: {engine}",
        )

    resolved_method = str(spec.method or "").strip() or _resolve_irc_method(config)
    resolved_basis = str(spec.basis or "").strip() or str(optimization.get("basis", "") or "").strip()
    resolved_spec = replace(spec, method=resolved_method, basis=resolved_basis)
    interface = _resolve_orca_interface_for_irc(
        config,
        method=resolved_method,
        basis=resolved_basis,
    )
    return interface.run_irc(
        spec=resolved_spec,
        input_xyz=input_xyz,
        output_dir=output_dir,
        charge=int(charge),
        spin=int(spin),
        timeout=getattr(spec, "timeout", None),
        subprocess_callback=subprocess_callback,
    )


def run_neb_ts(
    spec: QCJobSpec,
    input_xyz: Path,
    output_dir: Path,
    config: Mapping[str, Any],
    subprocess_callback: Optional[Callable[[Any], None]] = None,
) -> QCJobResult:
    """Execute an ORCA double-ended NEB-TS calculation."""
    input_xyz = Path(input_xyz).resolve()
    optimization = _s3_optimization_config(config)
    if str(optimization.get("engine", "orca") or "orca").strip().lower() != "orca":
        return QCJobResult(status="failed", input_xyz=input_xyz, error="NEB-TS requires ORCA")
    method = str(spec.method or _resolve_irc_method(config)).strip()
    basis = str(spec.basis or optimization.get("basis", "") or "").strip()
    interface = _resolve_orca_interface_for_irc(config, method=method, basis=basis)
    return interface.run_neb_ts(spec, input_xyz, output_dir, subprocess_callback=subprocess_callback)


def run_surface_scan(
    spec: SurfaceScanSpec,
    input_xyz: Path,
    output_dir: Path,
    config: Mapping[str, Any],
    subprocess_callback: Optional[Callable[[Any], None]] = None,
) -> QCJobResult:
    """Execute an ORCA relaxed surface scan through the shared QC boundary."""

    spec = resolve_surface_scan_spec(spec, config)
    input_xyz = Path(input_xyz).resolve()
    output_dir = Path(output_dir)
    if not spec.coordinates:
        return QCJobResult(status="failed", input_xyz=input_xyz, error="Surface scan requires coordinates")
    if len(spec.coordinates) > 3:
        return QCJobResult(status="failed", input_xyz=input_xyz, error="ORCA supports at most three scan coordinates")
    if spec.scan_ts and len(spec.coordinates) != 1:
        return QCJobResult(
            status="failed",
            input_xyz=input_xyz,
            error="ScanTS is restricted to one local scan coordinate",
        )
    optimization = _s3_optimization_config(config)
    if str(optimization.get("engine", "orca") or "orca").strip().lower() != "orca":
        return QCJobResult(status="failed", input_xyz=input_xyz, error="Relaxed surface scans require ORCA")
    nproc, _memory = _resources(dict(config), QCJobSpec(
        engine="orca", task="scan", method=spec.method, nproc=spec.nproc,
    ))
    interface = ORCAInterface(
        method=spec.method,
        basis=spec.basis,
        aux_basis=spec.aux_basis,
        nprocs=nproc,
        maxcore=spec.maxcore,
        solvent=spec.solvent or "",
        solvent_model=spec.solvent_model,
        route_extras=spec.route_extras,
        config=dict(config),
    )
    return interface.run_surface_scan(
        spec,
        input_xyz,
        output_dir,
        subprocess_callback=subprocess_callback,
    )


def resolve_surface_scan_spec(
    spec: SurfaceScanSpec,
    config: Mapping[str, Any],
) -> SurfaceScanSpec:
    """Fill an ORCA scan spec from the shared low-level ORCA profile.

    Surface scans and NEB are QC capabilities, not S2-specific engines.  A
    caller may provide a deliberate per-job override, but leaving a field
    empty must consistently inherit the project-wide B97-3c/solvation,
    resource and timeout settings.  This keeps S2 limited to reaction
    coordinates and selection policy.
    """

    optimization = _s3_optimization_config(config)
    resources = dict(config.get("resources", {}) or {}) if isinstance(config, Mapping) else {}
    method = str(spec.method or optimization.get("method", "") or "").strip()
    basis = str(spec.basis or optimization.get("basis", "") or "").strip()
    aux_basis = str(spec.aux_basis or optimization.get("aux_basis", "") or "").strip()
    solvent = str(spec.solvent or optimization.get("solvent", "") or "").strip() or None
    solvent_model = (
        str(spec.solvent_model or optimization.get("solvent_model", "") or "").strip()
        or None
    )
    route_extras = _merge_orca_route_extras(
        optimization.get("route_extras", ""),
        spec.route_extras,
    )
    nproc = spec.nproc if spec.nproc is not None else resources.get("nproc", 1)
    max_cycles = (
        spec.max_cycles
        if spec.max_cycles is not None
        else optimization.get("max_cycles")
    )
    timeout = (
        spec.timeout
        if spec.timeout is not None
        else optimization.get("timeout")
    )
    return replace(
        spec,
        method=method,
        basis=basis,
        aux_basis=aux_basis,
        solvent=solvent,
        solvent_model=solvent_model,
        route_extras=route_extras,
        nproc=max(1, int(nproc or 1)),
        max_cycles=(None if max_cycles is None else max(1, int(max_cycles))),
        timeout=(None if timeout is None else max(1, int(timeout))),
    )
