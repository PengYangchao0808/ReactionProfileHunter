"""
Gaussian TS Driver - Berny TS Optimization
============================================

Berny过渡态优化驱动器

Author: QCcalc Team
Date: 2026-01-09
"""

import subprocess
import re
import logging
from pathlib import Path
from typing import Optional, TYPE_CHECKING
from dataclasses import dataclass
import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.qc_interface import GaussianInterface

if TYPE_CHECKING:
    from rph_core.utils.optimization_config import OptimizationConfig

logger = logging.getLogger(__name__)


@dataclass
class TSOptResult:
    """TS优化结果"""
    converged: bool
    coordinates: 'np.ndarray'
    energy: float
    imaginary_freq: float
    frequencies: 'np.ndarray'
    output_file: Path
    method_used: str = "Berny"


class BernyTSDriver(LoggerMixin):
    """
    Berny TS 优化驱动器 (Step 3 - 主策略)

    功能:
    - 生成 Gaussian TS 输入文件
    - 提交 Gaussian 计算
    - 解析输出文件
    - 验证收敛和虚频
    """

    def __init__(self, method: str = "B3LYP", basis: str = "def2-SVP",
                 dispersion: str = "GD3BJ", nprocshared: int = 16,
                 mem: str = "32GB", opt_config: Optional['OptimizationConfig'] = None,
                 config: dict = {}):
        from rph_core.utils.optimization_config import OptimizationConfig

        self.method = method
        self.basis = basis
        self.dispersion = dispersion
        self.nprocshared = nprocshared
        self.mem = mem
        self.config = config

        if opt_config:
            self.route = opt_config.to_gaussian_route(method, basis, dispersion, is_ts=True)
        else:
            self.route = (
                f"# {method}/{basis} EmpiricalDispersion={dispersion} "
                f"Opt=(TS, CalcFC, NoEigenTest) Freq"
            )

        self.gaussian = GaussianInterface(
            nprocshared=nprocshared,
            mem=mem,
            config=self.config
        )

        self.logger.info(f"BernyTSDriver 初始化: {self.route}")

    def optimize(self, ts_guess: Path, output_dir: Path,
                 old_checkpoint: Optional[Path] = None) -> TSOptResult:
        """
        执行 Berny TS 优化

        Args:
            ts_guess: TS 初猜 XYZ 文件
            output_dir: 输出目录
            old_checkpoint: 可选的checkpoint文件 (P1: 复用轨道)

        Returns:
            TSOptResult 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"开始 Berny TS 优化: {ts_guess.name}")

        if old_checkpoint and old_checkpoint.exists():
            self.logger.info(f"  P1: 复用checkpoint {old_checkpoint.name}")

        # 直接调用 Gaussian 接口的优化方法
        res = self.gaussian.optimize(
            xyz_file=ts_guess,
            output_dir=output_dir,
            route=self.route,
            old_checkpoint=old_checkpoint  # P1: 传递checkpoint
        )

        if not res.converged:
            # 即使没完全收敛，我们也尝试获取当前结果
            self.logger.warning("Gaussian 优化未收敛")
            
        # 提取第一个虚频
        imag_freq = 0.0
        if res.frequencies is not None:
            imag_freqs = [f for f in res.frequencies if f < 0]
            if imag_freqs:
                imag_freq = imag_freqs[0]

        return TSOptResult(
            converged=res.converged,
            coordinates=res.coordinates,
            energy=res.energy,
            imaginary_freq=imag_freq,
            frequencies=res.frequencies,
            output_file=res.output_file,
            method_used="Berny"
        )
