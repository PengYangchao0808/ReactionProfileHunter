"""
TS Validator - 虚频和IRC验证
==============================

过渡态验证模块

Author: QCcalc Team
Date: 2026-01-09
"""

import re
import logging
from pathlib import Path
from typing import List, Tuple, Optional
from dataclasses import dataclass

import numpy as np

from rph_core.utils.geometry_tools import GeometryUtils


logger = logging.getLogger(__name__)


class TSValidationError(Exception):
    """TS验证错误"""
    pass


class TSValidator:
    """
    TS 验证器 (Step 3)

    功能:
    - 虚频验证（必须恰好1个）
    - IRC 验证（可选）
    - 键长验证
    """

    def __init__(self):
        """初始化验证器"""
        self.logger = logger

    def validate_imaginary_frequencies(self, frequencies: np.ndarray):
        """
        验证虚频

        要求: 恰好 1 个虚频

        Args:
            frequencies: 频率数组

        Raises:
            TSValidationError: 虚频数量不符合要求
        """
        imaginary_freqs = [f for f in frequencies if f < 0]

        if len(imaginary_freqs) == 0:
            raise TSValidationError(
                "未找到虚频，这可能不是过渡态。"
            )

        if len(imaginary_freqs) > 1:
            raise TSValidationError(
                f"找到 {len(imaginary_freqs)} 个虚频，期望恰好 1 个。"
                f"虚频值: {imaginary_freqs}"
            )

        self.logger.info(f"✓ 虚频验证通过: {imaginary_freqs[0]:.1f} cm⁻¹")

    def validate_imaginary_freq(self, frequencies: List[float]) -> bool:
        """
        验证虚频（旧接口，保持兼容）

        Args:
            frequencies: 频率列表

        Returns:
            是否有效
        """
        try:
            self.validate_imaginary_frequencies(np.array(frequencies))
            return True
        except TSValidationError:
            return False

    def validate_bond_lengths(
        self,
        coordinates: np.ndarray,
        bond_indices: List[Tuple[int, int]],
        expected_range: Tuple[float, float] = (2.0, 2.4)
    ):
        """
        验证键长是否在合理范围内

        Args:
            coordinates: 坐标数组 (N, 3)
            bond_indices: 键索引列表 [(i, j), ...]
            expected_range: 期望的键长范围 (min, max)

        Raises:
            TSValidationError: 键长超出范围
        """
        for i, j in bond_indices:
            dist = np.linalg.norm(coordinates[i] - coordinates[j])

            if not (expected_range[0] <= dist <= expected_range[1]):
                raise TSValidationError(
                    f"键长 {i}-{j} = {dist:.2f} Å "
                    f"超出范围 {expected_range}"
                )

        self.logger.info(f"✓ 键长验证通过")

    def validate_irc_endpoints(
        self,
        irc_result: 'IRCResult',
        reactant: Path,
        product: Path,
        tolerance: float = 2.0  # Å
    ):
        """
        验证 IRC 端点是否正确连接

        Args:
            irc_result: IRC 计算结果
            reactant: 底物结构文件
            product: 产物结构文件
            tolerance: RMSD 容忍度

        Raises:
            TSValidationError: IRC 端点不匹配
        """
        # 提取 IRC 终点坐标
        irc_forward_end = irc_result.forward_endpoint
        irc_reverse_end = irc_result.reverse_endpoint

        # 读取 reactant 和 product 坐标
        from rph_core.utils.file_io import read_xyz

        coords_r, _ = read_xyz(reactant)
        coords_p, _ = read_xyz(product)

        # 计算 RMSD
        rmsd_forward = self._calculate_rmsd(irc_forward_end, coords_r)
        rmsd_reverse = self._calculate_rmsd(irc_reverse_end, coords_p)

        self.logger.info(f"IRC RMSD: forward={rmsd_forward:.2f} Å, reverse={rmsd_reverse:.2f} Å")

        if rmsd_forward > tolerance or rmsd_reverse > tolerance:
            raise TSValidationError(
                f"IRC 端点不匹配: RMSD (forward={rmsd_forward:.2f}, "
                f"reverse={rmsd_reverse:.2f}) > {tolerance} Å"
            )

        self.logger.info(f"✓ IRC 端点验证通过")

    def _calculate_rmsd(self, coords1: np.ndarray, coords2: np.ndarray) -> float:
        """
        计算 RMSD（均方根偏差）

        Args:
            coords1: 第一个坐标集 (N, 3)
            coords2: 第二个坐标集 (N, 3)

        Returns:
            RMSD (Å)
        """
        if coords1.shape != coords2.shape:
            raise ValueError("坐标数组形状不匹配")

        diff = coords1 - coords2
        mse = np.mean(np.sum(diff**2, axis=1))
        rmsd = np.sqrt(mse)

        return rmsd

    def validate_no_atom_clash(
        self,
        coordinates: np.ndarray,
        min_distance: float = 0.8  # Å
    ):
        """
        验证没有原子重叠

        Args:
            coordinates: 坐标数组 (N, 3)
            min_distance: 最小允许距离

        Raises:
            TSValidationError: 存在原子重叠
        """
        n_atoms = len(coordinates)

        distance_matrix = GeometryUtils.compute_distance_matrix(coordinates)

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                dist = distance_matrix[i, j]

                if dist < min_distance:
                    raise TSValidationError(
                        f"原子 {i} 和 {j} 距离过近: {dist:.2f} Å < {min_distance} Å"
                    )

        self.logger.debug("✓ 无原子重叠")


@dataclass
class IRCResult:
    """IRC 计算结果"""
    converged: bool
    forward_endpoint: np.ndarray
    reverse_endpoint: np.ndarray
    forward_path: Optional[np.ndarray] = None
    reverse_path: Optional[np.ndarray] = None
    output_file: Optional[Path] = None
