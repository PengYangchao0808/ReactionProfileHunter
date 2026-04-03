import logging
import re
from pathlib import Path
from typing import List, Optional, Tuple

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
            r"Standard orientation:.*?---------------------------------------------------------------------.*?---------------------------------------------------------------------",
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
