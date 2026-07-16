"""
xTB Single-Point Hessian + MRRHO Thermo Interface
==================================================

Implements the CENSO-style ``xtb --bhess --enso`` workflow for
thermochemical property calculation via xTB SPH + mRRHO.

Used by the conformer search funnel to compute the thermostatistical
free-energy correction (G_RRHO) for each conformer candidate.

**Important**: The ``G(T)`` field in ``xtb_enso.json`` is the
*thermostatistical free-energy correction* (G_RRHO contribution),
**not** the total free energy.  The total free energy would be::

    G_total = E_electronic + G(T)

This follows the xTB official ENSO documentation which states:
    "free energy is obtained as the sum of the low-level DFT energy E,
     G_RRHO (thermostatistical contribution to free energy) and dG_solv"
"""

import os
import subprocess
import json
import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

from rph_core.utils.constants import BOHR_TO_ANGSTROM

logger = logging.getLogger(__name__)

__all__ = ["run_xtb_enso", "XTBThermoResult"]


def _tail(text: str, limit: int = 2000) -> str:
    """Return a compact tail of a subprocess stream for diagnostics."""
    if not text:
        return ""
    return text[-limit:]


def _failure_context(result: subprocess.CompletedProcess, output_dir: Path) -> str:
    """Collect xTB failure context from streams and marker files."""
    chunks: List[str] = []
    if result.stderr:
        chunks.append(f"stderr: {_tail(result.stderr)}")
    if result.stdout:
        chunks.append(f"stdout: {_tail(result.stdout)}")

    markers = [
        ".sccnotconverged",
        "NOT_CONVERGED",
        "xtbopt.log",
        "xtb_enso.json",
    ]
    present = [name for name in markers if (output_dir / name).exists()]
    if present:
        chunks.append(f"markers/files: {', '.join(present)}")

    xtbopt_log = output_dir / "xtbopt.log"
    if xtbopt_log.exists():
        try:
            chunks.append(f"xtbopt.log tail: {_tail(xtbopt_log.read_text(encoding='utf-8', errors='replace'))}")
        except OSError:
            pass

    return " | ".join(chunks) if chunks else "(no stdout/stderr/marker context)"


@dataclass(frozen=True)
class XTBThermoResult:
    """Result from xTB SPH + mRRHO thermochemistry calculation.

    The ``g_rrho_correction_hartree`` field is the *thermostatistical
    free-energy correction* (G_RRHO) parsed from ``xtb_enso.json["G(T)"]``.
    It is **not** the total free energy.  To obtain a ranked free energy::

        G_ranked = E_electronic(DFT) + g_rrho_correction_hartree

    Attributes:
        g_rrho_correction_hartree: Thermostatistical free-energy
            correction G(T) from xtb_enso.json (Hartree).  This is the
            G(RRHO) contribution, not the total free energy.
        zpve: Zero-point vibrational energy (Hartree), optional.
        h_total: Enthalpy correction H(T) (Hartree), optional.
        success: Whether the calculation completed successfully.
        error: Error message if the calculation failed, else None.
    """
    # Failed jobs carry ``None`` rather than numerical zero so missing
    # thermochemistry cannot silently become a valid correction.
    g_rrho_correction_hartree: Optional[float]
    xtb_electronic_energy_hartree: Optional[float] = None
    xtb_total_free_energy_hartree: Optional[float] = None
    energy_ledger_residual_hartree: Optional[float] = None
    xtb_solvation_energy_hartree: Optional[float] = None
    temperature_k: Optional[float] = None
    gfn_level: Optional[int] = None
    solvent_model: Optional[str] = None
    explicit_solvation_correction_added: bool = False
    zpve: Optional[float] = None
    h_total: Optional[float] = None
    success: bool = True
    error: Optional[str] = None

    @property
    def g_total(self) -> Optional[float]:
        """Backward-compatible alias for *g_rrho_correction_hartree*.

        .. deprecated::
            Use ``g_rrho_correction_hartree`` directly to avoid confusion
            with the total free energy.
        """
        return self.g_rrho_correction_hartree


def _xyz_to_coord(xyz_path: Path, coord_path: Path) -> None:
    """Convert an XYZ file to xTB .coord format.

    Plain Turbomole ``$coord`` data use Bohr by default. XYZ coordinates are
    therefore converted from Angstrom before writing lowercase atom symbols.

    Args:
        xyz_path: Path to the input XYZ file.
        coord_path: Path to write the .coord file.
    """
    xyz_path = Path(xyz_path)
    coord_path = Path(coord_path)

    with open(xyz_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    if len(lines) < 3:
        raise ValueError(f"XYZ file too short: {xyz_path}")

    natoms = int(lines[0].strip())
    atom_lines = lines[2:2 + natoms]

    with open(coord_path, 'w', encoding='utf-8') as f:
        f.write("$coord\n")
        for line in atom_lines:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            symbol = parts[0].lower()
            # Plain Turbomole $coord values are Bohr; XYZ values are Angstrom.
            x = float(parts[1]) / BOHR_TO_ANGSTROM
            y = float(parts[2]) / BOHR_TO_ANGSTROM
            z = float(parts[3]) / BOHR_TO_ANGSTROM
            f.write(f"{x:16.10f}  {y:16.10f}  {z:16.10f}  {symbol}\n")
        f.write("$end\n")


def _write_xcontrol(
    xcontrol_path: Path,
    temperature_k: float,
    sthr: float,
    imagthr: float,
) -> None:
    """Write the .xcontrol file with $thermo and $symmetry blocks.

    Args:
        xcontrol_path: Path to the .xcontrol file to write.
        temperature_k: Temperature in Kelvin.
        sthr: Rotational/vibrational entropy threshold (cm⁻¹ / K).
        imagthr: Imaginary frequency threshold for thermo (cm⁻¹).
    """
    content = (
        f"$thermo\n"
        f"    temp={temperature_k}\n"
        f"    sthr={sthr}\n"
        f"    imagthr={imagthr}\n"
        f"$symmetry\n"
        f"     maxat=1000\n"
        f"$end\n"
    )
    xcontrol_path.write_text(content, encoding='utf-8')


def run_xtb_enso(
    xtb_bin: Path,
    coord_file: Path,
    output_dir: Path,
    *,
    nproc: int = 1,
    gfn_level: int = 2,
    temperature_k: float = 298.15,
    sthr: float = 50.0,
    imagthr: float = -100.0,
    charge: int = 0,
    unpaired: int = 0,
    solvent: Optional[str] = None,
    timeout: Optional[int] = None,
    bhess_level: Optional[str] = "normal",
    omp_stacksize: Optional[str] = None,
    omp_max_active_levels: Optional[int] = 1,
    max_scc_iterations: Optional[int] = None,
    ledger_tolerance_hartree: float = 1.0e-7,
) -> XTBThermoResult:
    """Run xTB single-point Hessian + mRRHO (--bhess --enso) calculation.

    Writes a ``.xcontrol`` file into *output_dir*, invokes xTB, and
    parses the resulting ``xtb_enso.json`` file to extract the G(T)
    thermostatistical free-energy correction, ZPVE, and H(T).

    **Semantics note**: The ``G(T)`` field in ``xtb_enso.json`` is the
    *thermostatistical free-energy correction* (G_RRHO), not the total
    free energy.  The returned ``g_rrho_correction_hartree`` should be
    added to an electronic energy to obtain a ranked free energy.

    Args:
        xtb_bin: Path to the xTB binary.
        coord_file: Path to the .coord file (input geometry).
        output_dir: Directory for all output files (created if needed).
        nproc: Number of threads/cores. Sets ``--parallel`` flag and
            ``OMP_NUM_THREADS`` / ``MKL_NUM_THREADS`` /
            ``OPENBLAS_NUM_THREADS`` environment variables. Default 1.
        gfn_level: GFN-xTB level (0, 1, or 2). Default 2.
        temperature_k: Temperature in Kelvin. Default 298.15.
        sthr: Rotational/vibrational entropy threshold. Default 50.0.
        imagthr: Imaginary frequency threshold. Default -100.0.
        charge: Molecular charge. Default 0.
        unpaired: Number of unpaired electrons. Default 0.
        solvent: ALPB solvent name (e.g. ``"toluene"``). Default None.
        timeout: Subprocess timeout in seconds. Default 600.
        bhess_level: Numerical Hessian level passed after ``--bhess``.
        omp_stacksize: Optional per-thread OpenMP stack size. ``None`` removes
            any inherited oversized setting instead of forcing a multi-GB
            stack for every xTB thread.
        omp_max_active_levels: OpenMP nesting limit, or ``None`` to inherit.
        max_scc_iterations: Optional maximum SCC iterations for convergence
            retries. ``None`` retains the xTB default.
        ledger_tolerance_hartree: Maximum allowed residual for
            ``free energy - electronic energy - G(T)``.

    Returns:
        XTBThermoResult with parsed thermochemical data on success,
        or error details on failure.
    """
    output_dir = Path(output_dir)
    coord_file = Path(coord_file)
    xtb_bin = Path(xtb_bin)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write .xcontrol
    xcontrol_name = "xcontrol"
    xcontrol_path = output_dir / xcontrol_name
    _write_xcontrol(xcontrol_path, temperature_k, sthr, imagthr)

    # Build command
    cmd: List[str] = [
        str(xtb_bin),
        str(coord_file),
        "--gfn", str(gfn_level),
        "--bhess",
        "--enso",
        "--chrg", str(charge),
        "--uhf", str(unpaired),
        "--parallel", str(nproc),
        "-I", xcontrol_name,
    ]
    if bhess_level:
        cmd.insert(5, str(bhess_level))

    if solvent:
        cmd.extend(["--alpb", solvent])
    if max_scc_iterations is not None:
        cmd.extend(["--iterations", str(int(max_scc_iterations))])

    # Thread environment
    xtb_env = os.environ.copy()
    xtb_env["OMP_NUM_THREADS"] = str(nproc)
    xtb_env["MKL_NUM_THREADS"] = str(nproc)
    xtb_env["OPENBLAS_NUM_THREADS"] = str(nproc)
    if omp_stacksize:
        xtb_env["OMP_STACKSIZE"] = str(omp_stacksize)
    else:
        xtb_env.pop("OMP_STACKSIZE", None)
    if omp_max_active_levels is not None:
        xtb_env["OMP_MAX_ACTIVE_LEVELS"] = str(omp_max_active_levels)

    logger.debug(
        "Running xTB SPH+MRRHO: %s in %s",
        " ".join(cmd), output_dir,
    )

    # Run subprocess
    try:
        result = subprocess.run(
            cmd,
            cwd=output_dir,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=xtb_env,
        )

        (output_dir / "xtb_stdout.log").write_text(result.stdout or "", encoding='utf-8')
        (output_dir / "xtb_stderr.log").write_text(result.stderr or "", encoding='utf-8')
    except subprocess.TimeoutExpired:
        error_msg = f"xTB SPH+MRRHO timed out after {timeout}s"
        logger.error(error_msg)
        return XTBThermoResult(
            g_rrho_correction_hartree=None,
            success=False,
            error=error_msg,
        )
    except FileNotFoundError:
        error_msg = f"xTB binary not found: {xtb_bin}"
        logger.error(error_msg)
        return XTBThermoResult(
            g_rrho_correction_hartree=None,
            success=False,
            error=error_msg,
        )

    # Check return code
    if result.returncode != 0:
        context = _failure_context(result, output_dir)
        error_msg = f"xTB SPH+MRRHO failed (rc={result.returncode}): {context}"
        logger.error(error_msg)
        return XTBThermoResult(
            g_rrho_correction_hartree=None,
            success=False,
            error=error_msg,
        )

    # Parse xtb_enso.json
    json_path = output_dir / "xtb_enso.json"
    if not json_path.exists():
        error_msg = f"xTB SPH+MRRHO completed but xtb_enso.json not found in {output_dir}"
        logger.error(error_msg)
        return XTBThermoResult(
            g_rrho_correction_hartree=None,
            success=False,
            error=error_msg,
        )

    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        error_msg = f"Failed to parse xtb_enso.json: {e}"
        logger.error(error_msg)
        return XTBThermoResult(
            g_rrho_correction_hartree=None,
            success=False,
            error=error_msg,
        )

    # Extract G(T) — this is the G(RRHO) thermostatistical correction,
    # NOT the total free energy.
    try:
        g_rrho_correction_hartree = float(data["G(T)"])
        xtb_electronic_energy_hartree = float(data["energy"])
        xtb_total_free_energy_hartree = float(data["free energy"])
    except (KeyError, TypeError, ValueError) as e:
        error_msg = f"Incomplete or invalid xTB energy ledger in xtb_enso.json: {e}"
        logger.error(error_msg)
        return XTBThermoResult(
            g_rrho_correction_hartree=None,
            success=False,
            error=error_msg,
        )

    ledger_residual = (
        xtb_total_free_energy_hartree
        - xtb_electronic_energy_hartree
        - g_rrho_correction_hartree
    )
    tolerance = float(ledger_tolerance_hartree)
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("ledger_tolerance_hartree must be greater than zero")
    if not all(
        math.isfinite(value)
        for value in (
            g_rrho_correction_hartree,
            xtb_electronic_energy_hartree,
            xtb_total_free_energy_hartree,
            ledger_residual,
        )
    ):
        error_msg = "Non-finite value in xTB energy ledger"
        logger.error(error_msg)
        return XTBThermoResult(
            g_rrho_correction_hartree=None,
            xtb_electronic_energy_hartree=xtb_electronic_energy_hartree,
            xtb_total_free_energy_hartree=xtb_total_free_energy_hartree,
            energy_ledger_residual_hartree=ledger_residual,
            gfn_level=int(gfn_level),
            solvent_model=f"ALPB({solvent})" if solvent else None,
            success=False,
            error=error_msg,
        )
    if abs(ledger_residual) > tolerance:
        error_msg = (
            "Inconsistent xTB energy ledger: free energy - energy - G(T) = "
            f"{ledger_residual:.12g} Eh exceeds {tolerance:.3g} Eh"
        )
        logger.error(error_msg)
        return XTBThermoResult(
            g_rrho_correction_hartree=None,
            xtb_electronic_energy_hartree=xtb_electronic_energy_hartree,
            xtb_total_free_energy_hartree=xtb_total_free_energy_hartree,
            energy_ledger_residual_hartree=ledger_residual,
            temperature_k=float(data.get("temperature", temperature_k)),
            gfn_level=int(gfn_level),
            solvent_model=f"ALPB({solvent})" if solvent else None,
            success=False,
            error=error_msg,
        )

    # Extract optional fields
    zpve = None
    if "ZPVE" in data:
        try:
            zpve = float(data["ZPVE"])
        except (TypeError, ValueError):
            logger.warning("Invalid 'ZPVE' value in xtb_enso.json, ignoring")

    h_total = None
    if "H(T)" in data:
        try:
            h_total = float(data["H(T)"])
        except (TypeError, ValueError):
            logger.warning("Invalid 'H(T)' value in xtb_enso.json, ignoring")

    return XTBThermoResult(
        g_rrho_correction_hartree=g_rrho_correction_hartree,
        xtb_electronic_energy_hartree=xtb_electronic_energy_hartree,
        xtb_total_free_energy_hartree=xtb_total_free_energy_hartree,
        energy_ledger_residual_hartree=ledger_residual,
        # xtb_enso.json does not expose an independent ALPB solvation term.
        xtb_solvation_energy_hartree=None,
        temperature_k=float(data.get("temperature", temperature_k)),
        gfn_level=int(gfn_level),
        solvent_model=f"ALPB({solvent})" if solvent else None,
        explicit_solvation_correction_added=False,
        zpve=zpve,
        h_total=h_total,
        success=True,
        error=None,
    )
