import logging
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Tuple

from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.utils.path_compat import is_toxic_path
from rph_core.utils.thermo.schema import ThermoRecord
from rph_core.utils.thermo.options import ShermoOptions

logger = logging.getLogger(__name__)

# Shermo occasionally pauses for one or more bare Enter presses even in batch
# mode. Feed a small newline buffer explicitly so the subprocess never inherits
# interactive terminal stdin from the parent pipeline.
_SHERMO_AUTO_CONTINUE_INPUT = "\n\n"


def _extract_last_float(line: str) -> Optional[float]:
    tokens = [t for t in line.strip().split() if t]
    for token in reversed(tokens):
        try:
            return float(token)
        except ValueError:
            continue
    return None


def _extract_entropy_cal_mol_k(line: str) -> Optional[float]:
    """Extract entropy in cal/mol/K from a 'Total S' line.

    The Total S line format is:
        Total S:      300.011 J/mol/K      71.704 cal/mol/K    -TS:   -21.379 kcal/mol
    The cal/mol/K value is the second numeric token; the last token is -TS.
    """
    numeric_tokens = []
    for part in line.replace(":", " ").split():
        try:
            numeric_tokens.append(float(part))
        except ValueError:
            continue
    # Pattern: J/mol/K  cal/mol/K  -TS(kcal/mol)
    if len(numeric_tokens) >= 2:
        return numeric_tokens[1]
    return None


def parse_shermo_sum(sum_file: Path) -> ThermoRecord:
    """Parse a Shermo .sum file into a canonical ThermoRecord.

    Supports Shermo 2.0+ stdout format (primary) and legacy G(total) format (fallback).
    """
    content = sum_file.read_text(errors="ignore").splitlines()

    g_sum = None
    h_sum = None
    u_sum = None
    g_conc = None
    s_total = None
    parsed_temperature = 298.15

    # Primary: Shermo 2.0+ stdout "Sum of electronic energy and thermal correction to X"
    for line in content:
        # Parse the header temperature when present; fall back to 298.15 K if absent.
        temperature_match = re.search(r"Temperature:\s+(\d+\.?\d*)\s+K", line)
        if temperature_match:
            parsed_temperature = float(temperature_match.group(1))
        if "Sum of electronic energy and thermal correction to U" in line:
            u_sum = _extract_last_float(line)
        if "Sum of electronic energy and thermal correction to H" in line:
            h_sum = _extract_last_float(line)
        if "Sum of electronic energy and thermal correction to G" in line:
            g_sum = _extract_last_float(line)
        if "Gibbs free energy at specified concentration" in line:
            g_conc = _extract_last_float(line)
        if "Total S" in line:
            s_total = _extract_entropy_cal_mol_k(line)

    # Fallback: legacy "G(total)" / "H(total)" / "S(total)" format
    if g_sum is None or h_sum is None or u_sum is None:
        g_sum_fb, h_sum_fb, u_sum_fb, s_total_fb = _parse_legacy_format(content)
        if g_sum is None:
            g_sum = g_sum_fb
        if h_sum is None:
            h_sum = h_sum_fb
        if u_sum is None:
            u_sum = u_sum_fb
        if s_total is None:
            s_total = s_total_fb

    if g_sum is None or h_sum is None or u_sum is None:
        raise RuntimeError(f"Shermo output missing thermodynamic fields: {sum_file}")

    g_value = g_conc if g_conc is not None else g_sum

    return ThermoRecord(
        g_kcal=g_value * HARTREE_TO_KCAL,
        h_kcal=h_sum * HARTREE_TO_KCAL,
        u_kcal=u_sum * HARTREE_TO_KCAL,
        s_cal_mol_k=s_total,
        g_sum_hartree=g_sum,
        h_sum_hartree=h_sum,
        u_sum_hartree=u_sum,
        g_conc_hartree=g_conc,
        temperature_K=parsed_temperature,
        pressure_atm=None,
        concentration=None,
        source="shermo_sum",
        source_file=sum_file,
    )


def _parse_legacy_format(
    content: list[str],
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
    """Parse legacy G(total)/H(total)/S(total) format as fallback."""
    g_sum = None
    h_sum = None
    u_sum = None
    s_total = None

    for line in content:
        parts = line.split()
        if len(parts) < 2:
            continue
        label = parts[0].lower()
        try:
            val = float(parts[-1])
        except (ValueError, IndexError):
            continue
        if "G(total)" in line or label.startswith("g("):
            g_sum = val
        elif "H(total)" in line or label.startswith("h("):
            h_sum = val
        elif "U(total)" in line or label.startswith("u("):
            u_sum = val
        elif "S(total)" in line or label.startswith("s("):
            s_total = val

    return g_sum, h_sum, u_sum, s_total


def run_shermo_record(
    shermo_bin: Path,
    freq_output: Path,
    sp_energy: float,
    output_file: Path,
    options: ShermoOptions,
) -> ThermoRecord:
    """Invoke Shermo subprocess and return a ThermoRecord."""
    args = [str(shermo_bin), str(freq_output), "-E", f"{sp_energy:.12f}"]

    if options.temperature_k is not None:
        args.extend(["-T", f"{options.temperature_k}"])
    if options.pressure_atm is not None:
        args.extend(["-P", f"{options.pressure_atm}"])
    if options.scl_zpe is not None:
        args.extend(["-sclZPE", f"{options.scl_zpe}"])
    if options.ilowfreq is not None:
        args.extend(["-ilowfreq", f"{options.ilowfreq}"])
    if options.imagreal is not None:
        args.extend(["-imagreal", f"{options.imagreal}"])
    if options.conc:
        _cleaned = re.sub(r"\s+", "", str(options.conc))
        if _cleaned and _cleaned != "0":
            args.extend(["-conc", _cleaned])

    output_file.parent.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["Shermopath"] = str(shermo_bin.parent)

    def _run_shermo_subprocess(*, cwd: Optional[str] = None):
        return subprocess.run(
            args,
            input=_SHERMO_AUTO_CONTINUE_INPUT,
            capture_output=True,
            text=True,
            cwd=cwd,
            env=env,
        )

    if is_toxic_path(freq_output):
        temp_dir = Path(tempfile.mkdtemp(prefix="RPH_Shermo_", dir="/tmp"))
        try:
            temp_freq = temp_dir / "freq.log"
            shutil.copy2(freq_output, temp_freq)
            args[1] = str(temp_freq)
            result = _run_shermo_subprocess(cwd=str(temp_dir))
            if result.returncode != 0:
                raise RuntimeError(f"Shermo failed: {result.stderr}")
            output_file.write_text(result.stdout)
            if "Error:" in result.stdout:
                raise RuntimeError(f"Shermo error: {result.stdout[:500]}")
            record = parse_shermo(output_file)
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)
    else:
        result = _run_shermo_subprocess()
        if result.returncode != 0:
            raise RuntimeError(f"Shermo failed: {result.stderr}")
        output_file.write_text(result.stdout)
        if "Error:" in result.stdout:
            raise RuntimeError(f"Shermo error: {result.stdout[:500]}")
        record = parse_shermo(output_file)

    # Override record fields with actual invocation parameters
    return ThermoRecord(
        g_kcal=record.g_kcal,
        h_kcal=record.h_kcal,
        u_kcal=record.u_kcal,
        s_cal_mol_k=record.s_cal_mol_k,
        g_sum_hartree=record.g_sum_hartree,
        h_sum_hartree=record.h_sum_hartree,
        u_sum_hartree=record.u_sum_hartree,
        g_conc_hartree=record.g_conc_hartree,
        temperature_K=options.temperature_k,
        pressure_atm=options.pressure_atm,
        concentration=options.conc,
        source="shermo_run",
        source_file=output_file,
    )


def parse_shermo(output_file: Path) -> ThermoRecord:
    """Alias for parse_shermo_sum used internally by run_shermo_record."""
    return parse_shermo_sum(output_file)
