import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from rph_core.utils.qc_interface import QCResult

logger = logging.getLogger(__name__)


class GaussianLogParser:
    @staticmethod
    def parse_log(log_file: Path) -> Optional[QCResult]:
        content = ""
        try:
            with open(log_file, "r") as f:
                content = f.read()
        except OSError as exc:
            logger.error(f"无法读取日志文件 {log_file}: {exc}")
            return None

        converged = "Normal termination" in content
        energy_match = re.findall(r"SCF Done:\s+E\(\w+\)\s+=\s+(-?\d+\.\d+)", content)
        energy = float(energy_match[-1]) if energy_match else 0.0

        gibbs_match = re.findall(r"Sum of electronic and thermal Free Energies=\s+(-?\d+\.\d+)", content)
        if gibbs_match:
            energy = float(gibbs_match[-1])

        coords = GaussianLogParser._extract_coordinates(content)
        homo, lumo = GaussianLogParser._extract_orbitals(content)
        frequencies = GaussianLogParser._extract_frequencies(content)

        return QCResult(
            energy=energy,
            converged=converged,
            coordinates=coords,
            homo=homo,
            lumo=lumo,
            gap=(lumo - homo) if (lumo and homo) else None,
            frequencies=np.array(frequencies) if frequencies else None,
            output_file=log_file,
        )

    @staticmethod
    def _extract_coordinates(content: str) -> np.ndarray:
        sections = re.findall(
            r"Standard orientation:.*?---------------------------------------------------------------------.*?---------------------------------------------------------------------.*?---------------------------------------------------------------------",
            content,
            re.DOTALL,
        )
        if not sections:
            return np.array([])

        last_section = sections[-1]
        lines = last_section.split("\n")[5:-1]
        coords = []
        for line in lines:
            parts = line.split()
            if len(parts) == 6:
                coords.append([float(parts[3]), float(parts[4]), float(parts[5])])
        return np.array(coords)

    @staticmethod
    def _extract_orbitals(content: str) -> Tuple[Optional[float], Optional[float]]:
        h_to_ev = 27.2114
        occ_matches = re.findall(r"Alpha  occ\. eigenvalues -- (.*)", content)
        vir_matches = re.findall(r"Alpha virt\. eigenvalues -- (.*)", content)
        if not occ_matches or not vir_matches:
            return None, None

        last_occ_line = occ_matches[-1].split()
        homo = float(last_occ_line[-1]) * h_to_ev
        first_virt_line = vir_matches[0].split()
        lumo = float(first_virt_line[0]) * h_to_ev
        return homo, lumo

    @staticmethod
    def _extract_frequencies(content: str) -> List[float]:
        freq_matches = re.findall(r"Frequencies --\s+(.*)", content)
        freqs = []
        for line in freq_matches:
            freqs.extend([float(f) for f in line.split()])
        return freqs

    @staticmethod
    def _extract_displacement_vectors(content: str) -> Dict[int, np.ndarray]:
        """
        Parse displacement vectors for imaginary-frequency modes from a Gaussian log.

        Gaussian prints normal-mode displacement vectors in blocks of up to 3
        frequencies per table.  Each block layout (*n_modes* ∈ {1,2,3})::

                             1                      2                      3
                           ?A                     ?A                     ?A
            Frequencies --  -123.4567              45.6789               100.1234
            Red. masses --     5.6789               6.7890                 8.9012
            Frc consts  --     ...
            IR Inten    --     ...
            Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z
               1   6     0.12   0.34  -0.56     0.01   0.02   0.03      ...
               2   6     0.00   0.00   0.01     ...

        Every atom line: index (skip), atomic number (skip), then *n_modes*×3
        displacement components.  Only **imaginary** (negative) frequencies are
        collected; positive frequencies are skipped.

        Returns
        -------
        dict[int, np.ndarray]
            *mode_index* (0 for the first imaginary mode) → (*N_atoms*, 3) array.
            Empty dict if no imaginary modes or no frequency section found.
        """
        result: Dict[int, np.ndarray] = {}
        lines = content.split("\n")

        global_freq_idx = 0
        imaginary_count = 0

        i = 0
        while i < len(lines):
            line = lines[i]
            freq_match = re.match(r"\s*Frequencies\s+--\s+(.*)", line)
            if not freq_match:
                i += 1
                continue

            freq_str = freq_match.group(1)
            try:
                freqs = [float(f) for f in freq_str.split()]
            except ValueError:
                logger.warning("Failed to parse Frequencies line: %r", line.strip())
                i += 1
                continue
            n_modes = len(freqs)

            imaginary_cols = [col for col, freq in enumerate(freqs) if freq < 0]

            i += 1
            header_found = False
            while i < len(lines):
                if re.match(r"\s*Atom\s+AN\b", lines[i]):
                    header_found = True
                    break
                i += 1

            if not header_found:
                global_freq_idx += n_modes
                continue

            i += 1

            min_expected_fields = 2 + 3 * n_modes
            atom_rows: List[List[float]] = []

            while i < len(lines):
                atom_line = lines[i].strip()
                if not atom_line:
                    break
                parts = atom_line.split()
                if len(parts) < min_expected_fields:
                    break
                try:
                    _ = int(parts[0]), int(parts[1])
                except ValueError:
                    break
                try:
                    disp = [float(p) for p in parts[2:2 + 3 * n_modes]]
                except ValueError:
                    break
                atom_rows.append(disp)
                i += 1

            if not atom_rows:
                global_freq_idx += n_modes
                continue

            n_atoms = len(atom_rows)

            for col in imaginary_cols:
                array = np.zeros((n_atoms, 3))
                col_start = col * 3
                for atom_j, disp in enumerate(atom_rows):
                    array[atom_j, :] = disp[col_start:col_start + 3]
                result[imaginary_count] = array
                imaginary_count += 1

            global_freq_idx += n_modes

        if result:
            logger.debug(
                "Extracted %d imaginary-mode displacement vector(s) from Gaussian log",
                len(result),
            )
        return result


def extract_displacement_vectors(
    content: str, engine: str = "gaussian"
) -> Dict[int, np.ndarray]:
    """Parse imaginary-mode displacement vectors from a QC log.

    Parameters
    ----------
    content : str
        Raw log / output file text.
    engine : str
        ``"gaussian"`` or ``"orca"``.  Unknown engines log a warning and
        return an empty dict.

    Returns
    -------
    dict[int, np.ndarray]
        *mode_index* (0-based for imaginary modes) → (*N_atoms*, 3) array.
    """
    engine_lower = engine.lower()

    if engine_lower == "gaussian":
        return GaussianLogParser._extract_displacement_vectors(content)

    if engine_lower == "orca":
        from rph_core.utils.orca_interface import _parse_orca_displacement_vectors

        return _parse_orca_displacement_vectors(content)

    logger.warning(
        "extract_displacement_vectors: unknown engine '%s'; returning empty dict",
        engine,
    )
    return {}

