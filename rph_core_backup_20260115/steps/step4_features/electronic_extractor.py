"""
Electronic Descriptor Extractor
==============================

从量子化学计算结果中提取电子结构描述符 (HOMO, LUMO, Charges, etc.)

Author: QCcalc Team
Date: 2026-01-10
"""

import numpy as np
from typing import Dict, Optional
from pathlib import Path
import logging
import re

from rph_core.utils.qc_interface import QCResult

logger = logging.getLogger(__name__)


class ElectronicExtractor:
    """电子结构描述符提取器"""

    @staticmethod
    def extract(qc_result: QCResult, nbo_file: Optional[Path] = None) -> Dict[str, float]:
        """
        提取电子结构描述符

        Args:
            qc_result: 量子化学计算结果
            nbo_file: 可选的 NBO 输出文件路径 [Phase 2.2 新增]

        Returns:
            描述符字典
        """
        descriptors = {}

        # 1. 轨道能量 (Energy levels)
        if qc_result.homo is not None:
            descriptors['homo'] = qc_result.homo
        if qc_result.lumo is not None:
            descriptors['lumo'] = qc_result.lumo
        if qc_result.gap is not None:
            descriptors['gap'] = qc_result.gap

        # 2. 偶极矩 (Dipole)
        if qc_result.dipole_moment is not None:
            dipole_magnitude = np.linalg.norm(qc_result.dipole_moment)
            descriptors['dipole_magnitude'] = float(dipole_magnitude)
            descriptors['dipole_x'] = float(qc_result.dipole_moment[0])
            descriptors['dipole_y'] = float(qc_result.dipole_moment[1])
            descriptors['dipole_z'] = float(qc_result.dipole_moment[2])

        # 3. 电荷统计 (Charge statistics)
        if qc_result.charges is not None:
            charges = np.array(list(qc_result.charges.values()))
            descriptors['max_charge'] = float(np.max(charges))
            descriptors['min_charge'] = float(np.min(charges))
            descriptors['charge_variance'] = float(np.var(charges))
            descriptors['charge_std'] = float(np.std(charges))

        # 4. 总能量 (Energy)
        descriptors['total_energy'] = qc_result.energy

        # 5. [Phase 2.2 新增] NBO 分析
        if nbo_file and nbo_file.exists():
            nbo_descriptors = ElectronicExtractor._extract_nbo(nbo_file)
            descriptors.update(nbo_descriptors)
            logger.info(f"✓ NBO 分析完成: {len(nbo_descriptors)} 个描述符")

        return descriptors

    @staticmethod
    def calculate_reactivity_indices(homo: float, lumo: float) -> Dict[str, float]:
        """
        计算基于 Conceptual DFT 的反应活性指数

        Args:
            homo: HOMO 能级 (eV)
            lumo: LUMO 能级 (eV)

        Returns:
            活性指数字典
        """
        indices = {}
        
        # Koopmans 定理估算 IP 和 EA
        ip = -homo
        ea = -lumo

        # 电负性 (Electronegativity)
        chi = (ip + ea) / 2.0
        indices['electronegativity'] = chi

        # 化学硬度 (Hardness)
        eta = (ip - ea) / 2.0
        indices['chemical_hardness'] = eta

        # 化学软度 (Softness)
        if eta > 0:
            indices['chemical_softness'] = 1.0 / eta

        # 亲电指数 (Electrophilicity index)
        indices['electrophilicity'] = (chi ** 2) / (2.0 * eta) if eta > 0 else 0.0

        # 亲核指数 (Nucleophilicity index - 基于 HOMO)
        indices['nucleophilicity'] = -homo

        return indices

    @staticmethod
    def _extract_nbo(nbo_file: Path) -> Dict[str, float]:
        """
        [Phase 2.2 新增] 从 Gaussian NBO 输出文件提取描述符

        提取的 NBO 描述符:
        - Natural charges: 原子自然电荷
        - Wiberg bond orders: Wiberg 键序
        - Second-order perturbation: E(2) 能量分析

        Args:
            nbo_file: NBO 输出文件路径（通常为 .nbo 或 .nbo7）

        Returns:
            NBO 描述符字典
        """
        descriptors = {}

        try:
            content = nbo_file.read_text()

            # 1. 提取 Natural Charges (第 4 列通常是自然电荷)
            # 格式示例: "C   1  2.12345  1  -0.12345  ..."
            charge_pattern = re.compile(r'^\s*[A-Z][a-z]?\s+\d+\s+[\d\.\-]+\s+\d+\s+([\d\.\-]+)')
            charges = charge_pattern.findall(content)

            if charges:
                charge_values = [float(c[0]) for c in charges]
                descriptors['nbo_max_charge'] = float(np.max(charge_values))
                descriptors['nbo_min_charge'] = float(np.min(charge_values))
                descriptors['nbo_mean_charge'] = float(np.mean(charge_values))
                descriptors['nbo_charge_std'] = float(np.std(charge_values))
                logger.debug(f"提取到 {len(charges)} 个 NBO 原子电荷")

            # 2. 提取 Wiberg Bond Orders (键序)
            # 格式示例: "   1   2  1.0567  ..."
            bond_order_pattern = re.compile(r'^\s*\d+\s+\d+\s+([\d\.\-]+)')
            bond_orders = bond_order_pattern.findall(content)

            if bond_orders:
                descriptors['nbo_max_bond_order'] = float(np.max(bond_orders))
                descriptors['nbo_mean_bond_order'] = float(np.mean(bond_orders))
                descriptors['nbo_min_bond_order'] = float(np.min(bond_orders))
                logger.debug(f"提取到 {len(bond_orders)} 个 Wiberg 键序")

            # 3. 提取第二阶微扰能量 (E2 stabilization energies)
            # 格式示例: "   1    2   15.23 kcal/mol  ..."
            e2_pattern = re.compile(r'^\s*\d+\s+\d+\s+([\d\.\-]+)\s+kcal/mol')
            e2_energies = e2_pattern.findall(content)

            if e2_energies:
                descriptors['nbo_max_e2'] = float(np.max(e2_energies))
                descriptors['nbo_mean_e2'] = float(np.mean(e2_energies))
                descriptors['nbo_total_e2'] = float(np.sum(e2_energies))
                logger.debug(f"提取到 {len(e2_energies)} 个 E2 能量")

        except Exception as e:
            logger.warning(f"NBO 解析失败: {e}")

        return descriptors
