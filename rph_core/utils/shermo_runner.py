import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class ThermoResult:
    g_sum: float
    h_sum: float
    u_sum: float
    s_total: Optional[float]
    g_conc: Optional[float]
    output_file: Path


def _extract_last_float(line: str) -> Optional[float]:
    tokens = [token for token in line.strip().split() if token]
    for token in reversed(tokens):
        try:
            return float(token)
        except ValueError:
            continue
    return None


def _parse_sum_file(sum_file: Path) -> ThermoResult:
    content = sum_file.read_text(errors="ignore").splitlines()
    g_sum = None
    h_sum = None
    u_sum = None
    g_conc = None
    s_total = None

    for line in content:
        if "Sum of electronic energy and thermal correction to U" in line:
            u_sum = _extract_last_float(line)
        if "Sum of electronic energy and thermal correction to H" in line:
            h_sum = _extract_last_float(line)
        if "Sum of electronic energy and thermal correction to G" in line:
            g_sum = _extract_last_float(line)
        if "Gibbs free energy at specified concentration" in line:
            g_conc = _extract_last_float(line)
        if "Total S" in line or "Vibrational entropy" in line:
            s_total = _extract_last_float(line)

    if g_sum is None or h_sum is None or u_sum is None:
        raise RuntimeError(f"Shermo 输出缺少热力学字段: {sum_file}")

    return ThermoResult(
        g_sum=g_sum,
        h_sum=h_sum,
        u_sum=u_sum,
        s_total=s_total,
        g_conc=g_conc,
        output_file=sum_file
    )


def _clean_conc(conc: Optional[str]) -> Optional[str]:
    if not conc:
        return None
    cleaned = re.sub(r"\s+", "", str(conc))
    if cleaned == "" or cleaned == "0":
        return None
    return cleaned


def run_shermo(
    shermo_bin: Path,
    freq_output: Path,
    sp_energy: float,
    output_file: Path,
    temperature_k: Optional[float] = None,
    pressure_atm: Optional[float] = None,
    scl_zpe: Optional[float] = None,
    ilowfreq: Optional[int] = None,
    imagreal: Optional[float] = None,
    conc: Optional[str] = None
) -> ThermoResult:
    args = [str(shermo_bin), str(freq_output), "-E", f"{sp_energy:.12f}"]

    if temperature_k is not None:
        args.extend(["-T", f"{temperature_k}"])
    if pressure_atm is not None:
        args.extend(["-P", f"{pressure_atm}"])
    if scl_zpe is not None:
        args.extend(["-sclZPE", f"{scl_zpe}"])
    if ilowfreq is not None:
        args.extend(["-ilowfreq", f"{ilowfreq}"])
    if imagreal is not None:
        args.extend(["-imagreal", f"{imagreal}"])

    conc_value = _clean_conc(conc)
    if conc_value:
        args.extend(["-conc", conc_value])

    output_file.parent.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        args,
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"Shermo 运行失败: {result.stderr}")

    output_file.write_text(result.stdout)
    return _parse_sum_file(output_file)
