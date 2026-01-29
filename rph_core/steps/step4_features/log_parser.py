"""
Gaussian Log Parser
===================

从 Gaussian 输出文件 (.log/.out) 中提取能量、频率和描述符

Author: QCcalc Team
Date: 2026-01-10
"""

import re
import logging
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import numpy as np

from rph_core.utils.qc_interface import QCResult

logger = logging.getLogger(__name__)


class GaussianLogParser:
    """Gaussian 日志解析器"""

    @staticmethod
    def parse_log(log_file: Path) -> QCResult:
        """
        解析 Gaussian 日志文件

        Args:
            log_file: 日志文件路径

        Returns:
            QCResult 对象
        """
        content = ""
        try:
            with open(log_file, 'r') as f:
                content = f.read()
        except Exception as e:
            logger.error(f"无法读取日志文件 {log_file}: {e}")
            return None

        # 1. 检查收敛
        converged = "Normal termination" in content

        # 2. 提取能量 (SCF Done)
        energy_match = re.findall(r"SCF Done:\s+E\(\w+\)\s+=\s+(-?\d+\.\d+)", content)
        energy = float(energy_match[-1]) if energy_match else 0.0

        # 3. 提取 Gibbs 自由能 (Sum of electronic and thermal Free Energies)
        gibbs_match = re.findall(r"Sum of electronic and thermal Free Energies=\s+(-?\d+\.\d+)", content)
        if gibbs_match:
            energy = float(gibbs_match[-1])

        # 4. 提取坐标 (Last occurrence of Standard orientation)
        coords = GaussianLogParser._extract_coordinates(content)

        # 5. 提取 HOMO/LUMO
        homo, lumo = GaussianLogParser._extract_orbitals(content)

        # 6. 提取频率
        frequencies = GaussianLogParser._extract_frequencies(content)

        return QCResult(
            energy=energy,
            converged=converged,
            coordinates=coords,
            homo=homo,
            lumo=lumo,
            gap=(lumo - homo) if (lumo and homo) else None,
            frequencies=np.array(frequencies) if frequencies else None,
            output_file=log_file
        )

    @staticmethod
    def _extract_coordinates(content: str) -> np.ndarray:
        """提取标准取向坐标"""
        sections = re.findall(r"Standard orientation:.*?---------------------------------------------------------------------.*?---------------------------------------------------------------------", content, re.DOTALL)
        if not sections:
            return np.array([])
        
        last_section = sections[-1]
        lines = last_section.split('\n')[5:-1]
        
        coords = []
        for line in lines:
            parts = line.split()
            if len(parts) == 6:
                coords.append([float(parts[3]), float(parts[4]), float(parts[5])])
        
        return np.array(coords)

    @staticmethod
    def _extract_orbitals(content: str) -> Tuple[Optional[float], Optional[float]]:
        """提取 HOMO 和 LUMO (eV)"""
        # Hartree to eV
        H_TO_EV = 27.2114
        
        occ_matches = re.findall(r"Alpha  occ\. eigenvalues -- (.*)", content)
        vir_matches = re.findall(r"Alpha virt\. eigenvalues -- (.*)", content)
        
        if not occ_matches or not vir_matches:
            return None, None
            
        # 最后一个 occ 是 HOMO
        last_occ_line = occ_matches[-1].split()
        homo = float(last_occ_line[-1]) * H_TO_EV
        
        # 第一个 virt 是 LUMO
        first_virt_line = vir_matches[0].split()
        lumo = float(first_virt_line[0]) * H_TO_EV
        
        return homo, lumo

    @staticmethod
    def _extract_frequencies(content: str) -> List[float]:
        """提取振动频率"""
        freq_matches = re.findall(r"Frequencies --\s+(.*)", content)
        freqs = []
        for line in freq_matches:
            freqs.extend([float(f) for f in line.split()])
        return freqs
