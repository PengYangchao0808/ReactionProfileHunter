import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Any

from rph_core.utils.constants import HARTREE_TO_KCAL


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


def derive_shermo_summary_from_sum(
    sum_file: Path,
    output_json: Path,
    molecule_type: str = "precursor",
    temperature_K: float = 298.15
) -> Dict[str, Any]:
    """Derive shermo_summary.json from an existing .sum file.

    This enables Step1 activation to work even when only .sum files exist
    (extract-only approach, no QC execution).

    Args:
        sum_file: Path to Shermo .sum file
        output_json: Path to write derived JSON summary
        molecule_type: Type of molecule ('precursor', 'ylide', 'hoac')
        temperature_K: Temperature for thermodynamics (default 298.15 K)

    Returns:
        Dictionary containing the derived summary data
    """
    thermo = _parse_sum_file(sum_file)

    # Shermo reports energies in a.u. (Hartree). Convert to kcal/mol to match
    # Step4/Step1 activation expectations.
    g_value_au = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
    g_conc_au = thermo.g_conc

    summary = {
        "unit": "kcal/mol",
        "temperature_K": temperature_K,
        f"g_{molecule_type}": g_value_au * HARTREE_TO_KCAL,
        "g_sum": thermo.g_sum * HARTREE_TO_KCAL,
        "g_conc": (g_conc_au * HARTREE_TO_KCAL) if g_conc_au is not None else None,
        "h_sum": thermo.h_sum * HARTREE_TO_KCAL,
        "u_sum": thermo.u_sum * HARTREE_TO_KCAL,
        "s_total": thermo.s_total,
        "derived_from_sum": str(sum_file),
        "derived_artifacts": True,
    }

    # Write JSON for S4 extractors
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, 'w') as f:
        json.dump(summary, f, indent=2)

    return summary


def derive_hoac_thermo_from_sum(
    sum_file: Path,
    output_json: Path,
    temperature_K: float = 298.15
) -> Dict[str, Any]:
    """Derive HOAc thermo.json from Shermo .sum file.

    Creates a thermo.json compatible with step1_activation extractor,
    which expects 'g' or 'G' key (not 'g_hoac').

    Args:
        sum_file: Path to HOAc Shermo .sum file
        output_json: Path to write thermo.json
        temperature_K: Temperature for thermodynamics

    Returns:
        Dictionary with 'g', 'temperature_K', 'unit' keys
    """
    thermo = _parse_sum_file(sum_file)
    g_value = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum

    # Convert from Hartree to kcal/mol if needed
    # Shermo outputs in Hartree (a.u.) by default
    g_kcal = g_value * HARTREE_TO_KCAL

    thermo_data = {
        "g": g_kcal,
        "G": g_kcal,
        "temperature_K": temperature_K,
        "unit": "kcal/mol",
        "derived_from_sum": str(sum_file),
        "derived_artifacts": True
    }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, 'w') as f:
        json.dump(thermo_data, f, indent=2)

    return thermo_data


def find_shermo_sum_files(s1_dir: Path) -> Dict[str, Optional[Path]]:
    """Find Shermo .sum files in S1 directory.

    Searches for *_Shermo.sum files in common locations:
    - product/dft/ directory
    - Any subdirectory

    Args:
        s1_dir: S1 working directory

    Returns:
        Dictionary mapping molecule types to .sum file paths
    """
    import glob

    results: Dict[str, Optional[Path]] = {
        "precursor": None,
        "ylide": None,
        "hoac": None,
        "leaving_group": None,
    }

    # Search patterns for .sum files
    patterns = [
        str(s1_dir / "**" / "*_Shermo.sum"),
        str(s1_dir / "**" / "*Shermo*.sum"),
    ]

    def _is_hoac_like(token: str) -> bool:
        t = (token or "").lower()
        return any(key in t for key in ("hoac", "acoh", "acetic", "cc(=o)o", "c2h4o2"))

    # Track best (lowest) G value per type to prefer global-min conformer.
    best: Dict[str, tuple[float, Path]] = {}

    for pattern in patterns:
        for match in glob.glob(pattern, recursive=True):
            path = Path(match)
            parts_lower = [p.lower() for p in path.parts]
            stem_lower = path.stem.lower()

            mol_type: Optional[str] = None

            # Prefer directory-based classification (stable across naming schemes)
            # CRITICAL: leaving_group/ is the leaving group (e.g., HOAc), NOT the ylide.
            # The "ylide" (reactive intermediate) comes from S3_reactant, not S1.
            if "precursor" in parts_lower:
                mol_type = "precursor"
            elif "leaving_group" in parts_lower:
                # leaving_group/ contains the leaving group (e.g., acetic acid)
                # This is NOT the ylide - classify as "leaving_group" separately
                mol_type = "leaving_group"
            elif "ylide" in parts_lower:
                # Only explicit "ylide" directory maps to ylide
                mol_type = "ylide"
            elif "hoac" in parts_lower:
                mol_type = "hoac"
            elif "smallmolecules" in parts_lower or "small_molecules" in parts_lower:
                # Small-molecule cache keys vary; use heuristics
                for part in parts_lower:
                    if _is_hoac_like(part):
                        mol_type = "hoac"
                        break

            # Fallback to filename token classification
            if mol_type is None:
                if _is_hoac_like(stem_lower):
                    mol_type = "hoac"
                elif "ylide" in stem_lower or "leaving" in stem_lower:
                    mol_type = "ylide"
                elif "precursor" in stem_lower or "product" in stem_lower:
                    mol_type = "precursor"

            if mol_type is None:
                continue

            g_val = None
            try:
                thermo = _parse_sum_file(path)
                g_val = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
            except Exception:
                g_val = None

            # If parse fails, keep the first match.
            if g_val is None:
                if mol_type == "precursor" and results["precursor"] is None:
                    results["precursor"] = path
                elif mol_type == "ylide" and results["ylide"] is None:
                    results["ylide"] = path
                elif mol_type == "hoac" and results["hoac"] is None:
                    results["hoac"] = path
                elif mol_type == "leaving_group" and results["leaving_group"] is None:
                    results["leaving_group"] = path
                continue

            prev = best.get(mol_type)
            if prev is None or float(g_val) < float(prev[0]):
                best[mol_type] = (float(g_val), path)
                if mol_type == "precursor":
                    results["precursor"] = path
                elif mol_type == "ylide":
                    results["ylide"] = path
                elif mol_type == "hoac":
                    results["hoac"] = path
                elif mol_type == "leaving_group":
                    results["leaving_group"] = path

    return results
