"""
振荡检测模块 (Oscillation Detector)
====================================

检测几何优化过程中的能量振荡模式，并提供救援策略建议

Author: QC Descriptors Team
Date: 2026-01-13
Purpose: ReactionProfileHunter v2.1 - 优化控制与振荡救援
"""

from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Tuple
import numpy as np
import logging

logger = logging.getLogger(__name__)


class RescueLevel(Enum):
    """
    救援策略等级

    从低到高依次升级，每一步比前一步更激进但也更昂贵
    """
    NONE = "none"  # 无需救援
    RECALC_HESS = "recalc_hess"  # 增加Hessian重算频率
    REDUCE_STEP = "reduce_step"  # 减小步长和信任半径
    QST2 = "qst2"  # 使用 QST2 救援（端点优化）
    CALCALL = "calcall"  # 每步都计算Hessian（最昂贵）


@dataclass
class OscillationResult:
    """
    振荡检测结果

    Attributes:
        is_oscillating: 是否检测到振荡
        oscillation_type: 振荡类型 (energy, gradient, both)
        severity: 严重程度 (mild, moderate, severe)
        oscillation_count: 检测到的振荡次数
        energy_variance: 能量方差
        recommended_rescue: 推荐的救援策略
        rescue_params: 建议的救援参数调整
    """
    is_oscillating: bool
    oscillation_type: str  # 'energy', 'gradient', 'both'
    severity: str  # 'mild', 'moderate', 'severe'
    oscillation_count: int
    energy_variance: float
    recommended_rescue: RescueLevel
    rescue_params: dict


class OscillationDetector:
    """
    几何优化振荡检测器

    功能：
    1. 监测能量和梯度的历史变化
    2. 检测振荡模式（符号交替、能量聚集）
    3. 评估振荡严重程度
    4. 推荐适当的救援策略和参数调整
    """

    def __init__(
        self,
        window_size: int = 10,
        energy_tolerance: float = 1e-4,  # Hartree
        gradient_tolerance: float = 1e-3,  # Hartree/Bohr
        max_oscillation_count: int = 3
    ):
        """
        初始化振荡检测器

        Args:
            window_size: 检测窗口大小（步数）
            energy_tolerance: 能量变化容忍阈值（Hartree）
            gradient_tolerance: 梯度变化容忍阈值
            max_oscillation_count: 最大允许振荡次数
        """
        self.window_size = window_size
        self.energy_tolerance = energy_tolerance
        self.gradient_tolerance = gradient_tolerance
        self.max_oscillation_count = max_oscillation_count

        # 历史数据
        self.energy_history: List[float] = []
        self.gradient_history: List[float] = []  # 梯度的最大分量
        self.step_history: List[int] = []

        logger.debug(
            f"OscillationDetector 初始化: "
            f"window={window_size}, E_tol={energy_tolerance}, "
            f"G_tol={gradient_tolerance}"
        )

    def record_step(self, step: int, energy: float, gradient_max: float):
        """
        记录优化步骤

        Args:
            step: 步骤编号
            energy: 当前能量
            gradient_max: 最大梯度分量
        """
        self.step_history.append(step)
        self.energy_history.append(energy)
        self.gradient_history.append(gradient_max)

        # 限制历史长度
        if len(self.energy_history) > self.window_size * 2:
            self.energy_history.pop(0)
            self.gradient_history.pop(0)
            self.step_history.pop(0)

    def detect(self) -> OscillationResult:
        """
        检测振荡

        Returns:
            OscillationResult 对象
        """
        if len(self.energy_history) < 3:
            # 数据不足，无法检测
            return OscillationResult(
                is_oscillating=False,
                oscillation_type='none',
                severity='none',
                oscillation_count=0,
                energy_variance=0.0,
                recommended_rescue=RescueLevel.NONE,
                rescue_params={}
            )

        # 1. 检测能量振荡
        energy_oscillating, energy_count = self._detect_energy_oscillation()
        energy_var = self._calculate_energy_variance()

        # 2. 检测梯度振荡
        gradient_oscillating, gradient_count = self._detect_gradient_oscillation()

        # 3. 判断振荡类型和严重程度
        if not energy_oscillating and not gradient_oscillating:
            # 无振荡
            return OscillationResult(
                is_oscillating=False,
                oscillation_type='none',
                severity='none',
                oscillation_count=0,
                energy_variance=energy_var,
                recommended_rescue=RescueLevel.NONE,
                rescue_params={}
            )

        # 有振荡，确定类型和严重程度
        oscillation_type = 'both' if energy_oscillating and gradient_oscillating else \
                          'energy' if energy_oscillating else 'gradient'

        total_count = max(energy_count, gradient_count)

        # 根据振荡次数和方差判断严重程度
        if total_count >= self.max_oscillation_count:
            severity = 'severe'
        elif total_count >= self.max_oscillation_count // 2:
            severity = 'moderate'
        else:
            severity = 'mild'

        # 4. 推荐救援策略
        rescue_level, rescue_params = self._recommend_rescue(
            severity, oscillation_type, total_count, energy_var
        )

        logger.warning(
            f"检测到振荡: type={oscillation_type}, severity={severity}, "
            f"count={total_count}, rescue={rescue_level.value}"
        )

        return OscillationResult(
            is_oscillating=True,
            oscillation_type=oscillation_type,
            severity=severity,
            oscillation_count=total_count,
            energy_variance=energy_var,
            recommended_rescue=rescue_level,
            rescue_params=rescue_params
        )

    def _detect_energy_oscillation(self) -> Tuple[bool, int]:
        """
        检测能量振荡模式

        Returns:
            (是否振荡, 振荡次数)
        """
        oscillations = 0
        energies = self.energy_history[-self.window_size:] if len(self.energy_history) >= self.window_size else self.energy_history

        # 检测符号交替
        for i in range(2, len(energies)):
            delta1 = energies[i-1] - energies[i-2]
            delta2 = energies[i] - energies[i-1]

            # 如果两次能量变化方向相反且都超过容忍阈值
            if abs(delta1) > self.energy_tolerance and abs(delta2) > self.energy_tolerance:
                if (delta1 > 0 and delta2 < 0) or (delta1 < 0 and delta2 > 0):
                    oscillations += 1

        return oscillations >= 2, oscillations

    def _detect_gradient_oscillation(self) -> Tuple[bool, int]:
        """
        检测梯度振荡模式

        Returns:
            (是否振荡, 振荡次数)
        """
        oscillations = 0
        gradients = self.gradient_history[-self.window_size:] if len(self.gradient_history) >= self.window_size else self.gradient_history

        # 检测符号交替
        for i in range(2, len(gradients)):
            delta1 = gradients[i-1] - gradients[i-2]
            delta2 = gradients[i] - gradients[i-1]

            if abs(delta1) > self.gradient_tolerance and abs(delta2) > self.gradient_tolerance:
                if (delta1 > 0 and delta2 < 0) or (delta1 < 0 and delta2 > 0):
                    oscillations += 1

        return oscillations >= 2, oscillations

    def _calculate_energy_variance(self) -> float:
        """
        计算能量方差

        Returns:
            能量方差
        """
        energies = self.energy_history[-self.window_size:] if len(self.energy_history) >= self.window_size else self.energy_history
        if len(energies) < 2:
            return 0.0
        return float(np.var(energies))

    def _recommend_rescue(
        self,
        severity: str,
        oscillation_type: str,
        oscillation_count: int,
        energy_variance: float
    ) -> Tuple[RescueLevel, dict]:
        """
        推荐救援策略和参数调整

        Args:
            severity: 严重程度
            oscillation_type: 振荡类型
            oscillation_count: 振荡次数
            energy_variance: 能量方差

        Returns:
            (救援等级, 参数调整字典)
        """
        params = {}

        # 根据严重程度推荐不同策略
        if severity == 'mild':
            # 轻微振荡：增加Hessian重算频率
            params = {
                'recalc_hess_every': 5,  # 从10减到5
                'trust_radius': 0.2  # 降低信任半径
            }
            return RescueLevel.RECALC_HESS, params

        elif severity == 'moderate':
            # 中度振荡：减小步长，更频繁重算Hessian
            params = {
                'recalc_hess_every': 3,
                'max_step': 20,  # Gaussian: 0.20 Bohr (从30降低)
                'trust_radius': 0.15
            }
            return RescueLevel.REDUCE_STEP, params

        else:  # severe
            # 严重振荡：需要激进策略
            if oscillation_type == 'energy' and energy_variance > 1e-3:
                # 能量方差大，使用 QST2
                params = {}
                return RescueLevel.QST2, params
            else:
                # 梯度振荡或无法用QST2，使用 CalcAll
                params = {
                    'initial_hessian': 'calcall',
                    'recalc_hess_every': 1
                }
                return RescueLevel.CALCALL, params

    def reset(self):
        """重置历史数据"""
        self.energy_history.clear()
        self.gradient_history.clear()
        self.step_history.clear()
        logger.debug("OscillationDetector 历史数据已重置")

    def get_status_summary(self) -> str:
        """
        获取状态摘要

        Returns:
            状态摘要字符串
        """
        if len(self.energy_history) == 0:
            return "无数据"

        last_energy = self.energy_history[-1]
        last_gradient = self.gradient_history[-1]
        energy_range = max(self.energy_history) - min(self.energy_history)

        return (
            f"步骤数: {len(self.energy_history)}, "
            f"最新能量: {last_energy:.6f} Hartree, "
            f"能量范围: {energy_range:.6f} Hartree, "
            f"最大梯度: {last_gradient:.6f}"
        )


# 便捷函数
def detect_oscillation_from_output(
    output_lines: List[str],
    window_size: int = 10,
    energy_tolerance: float = 1e-4
) -> OscillationResult:
    """
    从输出文件中检测振荡

    Args:
        output_lines: 输出文件行列表
        window_size: 检测窗口
        energy_tolerance: 能量容忍阈值

    Returns:
        OscillationResult 对象
    """
    detector = OscillationDetector(
        window_size=window_size,
        energy_tolerance=energy_tolerance
    )

    # 解析输出文件提取能量和梯度
    import re

    for i, line in enumerate(output_lines):
        # Gaussian: SCF Done:  E(RB3LYP) =  -XXX.XXXXX
        match = re.search(r'SCF Done:\s+E\([^)]+\)\s+=\s+([-\d.]+)', line)
        if match:
            energy = float(match.group(1))
            # 梯度通常在下一行或附近的收敛信息中
            detector.record_step(i, energy, 0.0)  # 梯度暂设为0

    return detector.detect()
