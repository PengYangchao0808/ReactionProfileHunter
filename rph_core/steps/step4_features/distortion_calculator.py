"""
Distortion/Interaction Calculator
==================================

评估反应过渡态的畸变能 (Distortion Energy) 和相互作用能 (Interaction Energy)

Author: QCcalc Team
Date: 2026-01-10
"""

import logging
from typing import Dict, List, Tuple
import numpy as np

logger = logging.getLogger(__name__)


class DistortionCalculator:
    """畸变/相互作用分析引擎"""

    HARTREE_TO_KCAL = 627.509  # Hartree to kcal/mol conversion factor

    @staticmethod
    def calculate_activation_energy(e_ts: float, e_reactant: float) -> float:
        """计算活化能 (kcal/mol)"""
        return (e_ts - e_reactant) * DistortionCalculator.HARTREE_TO_KCAL

    @staticmethod
    def calculate_reaction_energy(e_product: float, e_reactant: float) -> float:
        """计算反应能 (kcal/mol)"""
        return (e_product - e_reactant) * DistortionCalculator.HARTREE_TO_KCAL

    @staticmethod
    def calculate_distortion_interaction(
        e_ts: float,
        e_fragment_a_at_ts: float,
        e_fragment_b_at_ts: float,
        e_fragment_a_relaxed: float,
        e_fragment_b_relaxed: float
    ) -> Dict[str, float]:
        """
        计算 Distortion/Interaction 特征

        公式:
        - E_dist_A = E(A_at_ts) - E(A_relaxed)
        - E_dist_B = E(B_at_ts) - E(B_relaxed)
        - E_dist_total = E_dist_A + E_dist_B
        - E_int = (E_ts - e_reactant) - E_dist_total

        Returns:
            Dict containing distortion and interaction energies
        """
        e_dist_a = (e_fragment_a_at_ts - e_fragment_a_relaxed) * DistortionCalculator.HARTREE_TO_KCAL
        e_dist_b = (e_fragment_b_at_ts - e_fragment_b_relaxed) * DistortionCalculator.HARTREE_TO_KCAL
        e_dist_total = e_dist_a + e_dist_b
        
        # 活化能
        e_activation = (e_ts - (e_fragment_a_relaxed + e_fragment_b_relaxed)) * DistortionCalculator.HARTREE_TO_KCAL
        e_int = e_activation - e_dist_total

        return {
            "e_distortion_a": e_dist_a,
            "e_distortion_b": e_dist_b,
            "e_distortion_total": e_dist_total,
            "e_interaction": e_int,
            "e_activation": e_activation
        }

    @staticmethod
    def calculate_asynchronicity(bond_1_len: float, bond_2_len: float) -> float:
        """
        计算非同步性指数 (Asynchronicity)

        公式: |r1 - r2| / (r1 + r2)
        """
        if bond_1_len + bond_2_len == 0:
            return 0.0
        return abs(bond_1_len - bond_2_len) / (bond_1_len + bond_2_len)

    @staticmethod
    def calculate_bond_lengths(coords: np.ndarray, bond_1_indices: Tuple[int, int], 
                              bond_2_indices: Tuple[int, int]) -> Tuple[float, float]:
        """从坐标计算两根键的长度"""
        r1 = np.linalg.norm(coords[bond_1_indices[0]] - coords[bond_1_indices[1]])
        r2 = np.linalg.norm(coords[bond_2_indices[0]] - coords[bond_2_indices[1]])
        return float(r1), float(r2)
       