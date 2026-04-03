"""
Bond Stretcher - Geometry Manipulation
=======================================

实现分子几何操作，特别是用于 TS 初猜生成的键长拉伸

Author: QCcalc Team (based on original by RetroTS Team)
Date: 2026-01-10
"""

import logging
from typing import List, Tuple, Optional
from dataclasses import dataclass
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class StretchingParams:
    """键拉伸参数"""
    target_length_A: float = 2.2  # 目标键长 (Å) - ✅ PROMOTE.md 标准
    min_length: float = 2.00       # 下限
    max_length: float = 2.30       # 上限
    step_size: float = 0.05        # 扫描步长


class BondStretcher:
    """
    键长拉伸引擎

    核心算法:
    1. 计算键向量 v = pos_B - pos_A
    2. 计算拉伸量 delta = target - current
    3. 沿键向量方向移动原子 (对称分配)
    """

    def __init__(self, params: Optional[StretchingParams] = None):
        """
        初始化键拉伸器

        Args:
            params: 拉伸参数，默认使用标准 TS 参数
        """
        self.params = params or StretchingParams()
        logger.info(f"键拉伸器初始化: 目标长度 = {self.params.target_length_A} Å")

    def stretch_bond(
        self,
        coords: np.ndarray,
        atom_i: int,
        atom_j: int,
        target_length: Optional[float] = None,
        fix_center_of_mass: bool = True
    ) -> np.ndarray:
        """
        拉伸指定键到目标长度

        Args:
            coords: 坐标数组 (N, 3)
            atom_i: 第一个原子索引
            atom_j: 第二个原子索引
            target_length: 目标键长
            fix_center_of_mass: 是否固定质心

        Returns:
            修改后的坐标数组 (N, 3)
        """
        if target_length is None:
            target_length = self.params.target_length_A

        # 复制坐标
        new_coords = coords.copy()

        # 获取两个原子的位置
        pos_i = new_coords[atom_i]
        pos_j = new_coords[atom_j]

        # 计算当前键长
        bond_vector = pos_j - pos_i
        current_length = np.linalg.norm(bond_vector)

        if current_length == 0:
            raise ValueError(f"原子 {atom_i} 和 {atom_j} 重合，无法计算键向量")

        # 计算需要拉伸的距离
        delta = target_length - current_length

        # 单位向量
        unit_vector = bond_vector / current_length

        # 对称移动: 每个原子移动 delta/2
        displacement = (delta / 2.0) * unit_vector

        new_coords[atom_i] -= displacement
        new_coords[atom_j] += displacement

        # 可选: 固定质心
        if fix_center_of_mass:
            original_com = np.mean(coords, axis=0)
            new_com = np.mean(new_coords, axis=0)
            new_coords += (original_com - new_com)

        return new_coords

    def stretch_bonds(
        self,
        coords: np.ndarray,
        bonds: List[Tuple[Tuple[int, int], float]],
        fix_center_of_mass: bool = True
    ) -> np.ndarray:
        """
        拉伸任意数量的键到各自目标长度

        Args:
            coords: 坐标数组 (N, 3)
            bonds: 键列表，每项为 ((atom_i, atom_j), target_length)
            fix_center_of_mass: 是否固定质心

        Returns:
            修改后的坐标数组 (N, 3)

        Example:
            >>> coords = np.array([[0,0,0], [1,0,0], [2,0,0], [3,0,0]])
            >>> bonds = [((0, 1), 1.5), ((2, 3), 1.8)]
            >>> new_coords = stretcher.stretch_bonds(coords, bonds)
        """
        if not bonds:
            return coords

        new_coords = coords.copy()
        bond_descriptions = []

        for bond, target_length in bonds:
            atom_i, atom_j = bond
            new_coords = self.stretch_bond(
                new_coords, atom_i, atom_j, target_length, fix_center_of_mass=False
            )
            bond_descriptions.append(f"{bond} → {target_length} Å")

        # 只在最后固定一次质心
        if fix_center_of_mass:
            original_com = np.mean(coords, axis=0)
            new_com = np.mean(new_coords, axis=0)
            new_coords += (original_com - new_com)

        logger.info(f"✓ 拉伸了 {len(bonds)} 根键: {bond_descriptions}")

        return new_coords

    def stretch_two_bonds(
        self,
        coords: np.ndarray,
        bond1: Tuple[int, int],
        bond2: Tuple[int, int],
        target_length: Optional[float] = None
    ) -> np.ndarray:
        """
        同时拉伸两根键 (用于 [5+2] Retro-TS)

        Args:
            coords: 坐标数组
            bond1: 第一根键 (atom_i, atom_j)
            bond2: 第二根键 (atom_k, atom_l)
            target_length: 目标键长

        Returns:
            修改后的坐标数组
        """
        if target_length is None:
            target_length = self.params.target_length_A

        # 使用 stretch_bonds 统一处理
        bonds = [
            (bond1, target_length),
            (bond2, target_length)
        ]
        return self.stretch_bonds(coords, bonds)


def stretch_bonds(
    coords: np.ndarray,
    bonds: List[Tuple[Tuple[int, int], float]],
    fix_center_of_mass: bool = True
) -> np.ndarray:
    """
    模块级函数：拉伸任意数量的键

    Args:
        coords: 坐标数组 (N, 3)
        bonds: 键列表，每项为 ((atom_i, atom_j), target_length)
        fix_center_of_mass: 是否固定质心

    Returns:
        修改后的坐标数组 (N, 3)
    """
    stretcher = BondStretcher()
    return stretcher.stretch_bonds(coords, bonds, fix_center_of_mass)
