"""
Quantum Chemistry Task Runner (QCTaskRunner)
=============================================
统一计算中枢 - 负责所有几何优化和单点能计算

职责:
1. Normal 模式：基态优化 → 频率验证 → (失败则) 救援 → L2 高精度 SP
2. TS 模式：过渡态优化 → 虚频验证 → (失败则) 救援 → L2 高精度 SP

优势: 一处修改，全案通用（底物、产物、过渡态全部调用此逻辑）

Author: QCcalc Team
Date: 2026-01-13
Version: v2.2-Rigor-Efficient
"""

import logging
from pathlib import Path
from typing import Optional, Tuple, Union
from dataclasses import dataclass
import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.optimization_config import OptimizationConfig
from rph_core.utils.qc_interface import (
    XTBInterface,
    GaussianInterface,
    QCResult,
    QCInterfaceFactory
)
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.utils.file_io import write_xyz, read_xyz

logger = logging.getLogger(__name__)


@dataclass
class QCSPResult:
    """单点能计算结果"""
    energy: float  # Hartree
    converged: bool
    output_file: Optional[Path] = None
    error_message: Optional[str] = None


@dataclass
class QCOptimizationResult:
    """完整优化结果（含频率验证和 L2 SP）"""
    optimized_xyz: Path  # 优化后的结构 XYZ
    l2_energy: Optional[float] = None  # L2 高精度能量 (Hartree)
    opt_energy: Optional[float] = None  # 几何优化能量 (Hartree)
    converged: bool = False
    frequencies: Optional[np.ndarray] = None
    imaginary_count: int = 0
    method_used: str = ""  # 使用的方法 (Normal/Berny/Berny_Rescue/QST2)
    error_message: Optional[str] = None
    checkpoint_file: Optional[Path] = None  # .chk 文件 (可选)


class QCTaskRunner(LoggerMixin):
    """
    统一计算中枢

    职责:
    1. Normal 模式：基态优化 + 频率验证 + 救援 + L2 SP
    2. TS 模式：过渡态优化 + 虚频验证 + 救援 + L2 SP

    设计原则:
    - 单一职责：只负责几何优化和单点能计算
    - 统一接口：所有优化流程通过此模块
    - 智能救援：标准策略失败时自动触发救援
    """

    def __init__(self, config: dict):
        """
        初始化 QCTaskRunner

        Args:
            config: 配置字典，包含:
                - theory.optimization: 几何优化配置
                - theory.single_point: L2 SP 配置
                - optimization_control: 优化控制参数
        """
        self.config = config
        self.theory_opt = config.get('theory', {}).get('optimization', {})
        self.theory_sp = config.get('theory', {}).get('single_point', {})

        # 资源配置
        self.nprocshared = self.theory_opt.get('nproc', 16)
        self.mem = self.theory_opt.get('mem', '32GB')

        # 优化配置
        self.method = self.theory_opt.get('method', 'B3LYP')
        self.basis = self.theory_opt.get('basis', 'def2-SVP')
        self.dispersion = self.theory_opt.get('dispersion', 'GD3BJ')
        self.engine_type = self.theory_opt.get('engine', 'gaussian').lower()

        # 默认优化配置
        self.opt_config = OptimizationConfig.from_config(config)

        # 初始化几何优化接口
        if self.engine_type == 'gaussian':
            self.qc_engine = GaussianInterface(
                charge=0,
                multiplicity=1,
                nprocshared=self.nprocshared,
                mem=self.mem,
                config=self.config
            )
        elif self.engine_type == 'orca':
            self.qc_engine = ORCAInterface(
                method=self.method,
                basis=self.basis,
                nprocs=self.nprocshared,
                config=self.config
            )
        else:
            raise ValueError(f"不支持的引擎: {self.engine_type}")

        # 初始化 ORCA L2 SP 接口
        self.orca = ORCAInterface(
            method=self.theory_sp.get('method', 'WB97M-V'),
            basis=self.theory_sp.get('basis', 'def2-TZVPP'),
            aux_basis=self.theory_sp.get('aux_basis', 'def2/J'),
            nprocs=self.theory_sp.get('nproc', 16),
            maxcore=self.theory_sp.get('maxcore', 4000),
            solvent=self.theory_sp.get('solvent', 'acetone'),
            config=self.config
        )

        self.logger.info(f"QCTaskRunner 初始化: {self.engine_type} {self.method}/{self.basis}")
        self.logger.info(f"  L2 SP: {self.orca.method}/{self.orca.basis}")

    def run_opt_sp_cycle(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int = 0,
        spin: int = 1,
        enable_l2_sp: bool = True,
        old_checkpoint: Optional[Path] = None
    ) -> QCOptimizationResult:
        """
        执行基态优化 + 频率验证 + L2 SP 完整流程

        策略:
        1. 标准几何优化
        2. 频率分析验证无虚频
        3. 失败则救援: CalcFC + MaxStep=10
        4. L2 高精度单点能

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录
            charge: 分子电荷
            spin: 自旋多重度
            enable_l2_sp: 是否启用 L2 SP
            old_checkpoint: 复用之前的 checkpoint

        Returns:
            QCOptimizationResult 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"=== Normal 模式优化: {xyz_file.name} ===")

        # 尝试标准优化
        opt_result = self._try_normal_optimization(
            xyz_file, output_dir, charge, spin, old_checkpoint
        )

        # 如果标准优化失败，尝试救援
        if not opt_result.converged:
            self.logger.warning("标准优化失败，启动救援策略...")
            opt_result = self._try_normal_rescue(
                xyz_file, output_dir, charge, spin, old_checkpoint
            )

        # L2 高精度 SP
        l2_energy = None
        if enable_l2_sp and opt_result.converged:
            self.logger.info("执行 L2 高精度单点能...")
            l2_energy = self._run_l2_sp(opt_result.optimized_xyz, output_dir / "L2_SP")
            opt_result.l2_energy = l2_energy

        return opt_result

    def run_ts_opt_cycle(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int = 0,
        spin: int = 1,
        enable_l2_sp: bool = True,
        old_checkpoint: Optional[Path] = None
    ) -> QCOptimizationResult:
        """
        执行 TS 优化 + 虚频验证 + L2 SP 完整流程

        策略:
        1. 标准 Berny TS 优化 (Opt=TS, CalcFC, NoEigenTest)
        2. 验证恰好 1 个虚频
        3. 失败则救援: Recalc=5 + NoEigenTest + MaxStep=10
        4. L2 高精度单点能

        Args:
            xyz_file: TS 猜想 XYZ 文件
            output_dir: 输出目录
            charge: 分子电荷
            spin: 自旋多重度
            enable_l2_sp: 是否启用 L2 SP
            old_checkpoint: 复用之前的 checkpoint

        Returns:
            QCOptimizationResult 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"=== TS 模式优化: {xyz_file.name} ===")

        # 尝试标准 TS 优化
        opt_result = self._try_ts_optimization(
            xyz_file, output_dir, charge, spin, old_checkpoint
        )

        # 如果 TS 优化失败，尝试救援
        if not opt_result.converged or opt_result.imaginary_count != 1:
            self.logger.warning("TS 优化失败或虚频不符，启动救援策略...")
            opt_result = self._try_ts_rescue(
                xyz_file, output_dir, charge, spin, old_checkpoint
            )

        # L2 高精度 SP
        l2_energy = None
        if enable_l2_sp and opt_result.converged:
            self.logger.info("执行 L2 高精度单点能...")
            l2_energy = self._run_l2_sp(opt_result.optimized_xyz, output_dir / "L2_SP")
            opt_result.l2_energy = l2_energy

        return opt_result

    def _count_imaginary(self, frequencies: Optional[np.ndarray]) -> int:
        """计算虚频数量"""
        if frequencies is None:
            return 0
        return int(np.sum(frequencies < 0))

    def _try_normal_optimization(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int,
        spin: int,
        old_checkpoint: Optional[Path]
    ) -> QCOptimizationResult:
        """尝试标准基态优化"""
        self.logger.info("尝试标准基态优化...")

        # 生成 route card
        route = self.opt_config.to_gaussian_route(
            self.method, self.basis, self.dispersion, is_ts=False
        )
        route += " Freq"

        # 执行优化
        result = self.qc_engine.optimize(
            xyz_file,
            output_dir / "standard",
            route=route,
            old_checkpoint=old_checkpoint
        )

        # 验证结果
        imaginary_count = self._count_imaginary(result.frequencies)

        if result.converged and imaginary_count == 0:
            self.logger.info("✓ 标准优化成功，无虚频")
            return QCOptimizationResult(
                optimized_xyz=result.output_file if result.output_file else xyz_file,
                opt_energy=result.energy,
                converged=True,
                frequencies=result.frequencies,
                imaginary_count=0,
                method_used="Normal"
            )
        else:
            self.logger.warning(
                f"标准优化失败 (converged={result.converged}, imaginary={imaginary_count})"
            )
            return QCOptimizationResult(
                optimized_xyz=result.output_file if result.output_file else xyz_file,
                opt_energy=result.energy,
                converged=False,
                frequencies=result.frequencies,
                imaginary_count=imaginary_count,
                method_used="Normal_Failed"
            )

    def _try_normal_rescue(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int,
        spin: int,
        old_checkpoint: Optional[Path]
    ) -> QCOptimizationResult:
        """尝试基态救援 (CalcFC + MaxStep=10)"""
        self.logger.info("救援: CalcFC + MaxStep=10")

        # 创建救援配置
        rescue_config = OptimizationConfig()
        rescue_config.update_for_rescue({
            'initial_hessian': 'calcfc',
            'max_step': 10,
            'recalc_hess_every': 5
        })

        route = rescue_config.to_gaussian_route(
            self.method, self.basis, self.dispersion, is_ts=False
        )
        route += " Freq"

        result = self.qc_engine.optimize(
            xyz_file,
            output_dir / "rescue",
            route=route,
            old_checkpoint=old_checkpoint
        )

        imaginary_count = self._count_imaginary(result.frequencies)

        if result.converged and imaginary_count == 0:
            self.logger.info("✓ 救援成功")
            return QCOptimizationResult(
                optimized_xyz=result.output_file if result.output_file else xyz_file,
                opt_energy=result.energy,
                converged=True,
                frequencies=result.frequencies,
                imaginary_count=0,
                method_used="Normal_Rescue"
            )
        else:
            self.logger.error(f"救援失败 (converged={result.converged}, imaginary={imaginary_count})")
            return QCOptimizationResult(
                optimized_xyz=result.output_file if result.output_file else xyz_file,
                opt_energy=result.energy,
                converged=False,
                frequencies=result.frequencies,
                imaginary_count=imaginary_count,
                method_used="Normal_Rescue_Failed"
            )

    def _try_ts_optimization(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int,
        spin: int,
        old_checkpoint: Optional[Path]
    ) -> QCOptimizationResult:
        """尝试标准 TS 优化"""
        self.logger.info(f"尝试标准 TS 优化 (引擎: {self.engine_type})...")

        if self.engine_type == 'gaussian':
            return self._try_ts_optimization_gaussian(
                xyz_file, output_dir, charge, spin, old_checkpoint
            )
        elif self.engine_type == 'orca':
            return self._try_ts_optimization_orca(
                xyz_file, output_dir, charge, spin, old_checkpoint
            )
        else:
            raise ValueError(f"不支持的 TS 优化引擎: {self.engine_type}")

    def _try_ts_optimization_gaussian(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int,
        spin: int,
        old_checkpoint: Optional[Path]
    ) -> QCOptimizationResult:
        """Gaussian TS 优化"""
        route = self.opt_config.to_gaussian_route(
            self.method, self.basis, self.dispersion, is_ts=True
        )
        route += " Freq"

        result = self.qc_engine.optimize(
            xyz_file,
            output_dir / "berny",
            route=route,
            old_checkpoint=old_checkpoint
        )

        imaginary_count = self._count_imaginary(result.frequencies)

        if result.converged and imaginary_count == 1:
            imaginary_freq = float(np.min(result.frequencies)) if result.frequencies is not None else 0.0
            self.logger.info(f"✓ Berny TS 优化成功，虚频 = {imaginary_freq:.1f} cm⁻¹")
            return QCOptimizationResult(
                optimized_xyz=result.output_file if result.output_file else xyz_file,
                opt_energy=result.energy,
                converged=True,
                frequencies=result.frequencies,
                imaginary_count=1,
                method_used="Berny"
            )
        else:
            self.logger.warning(
                f"Berny TS 优化失败或虚频不符 "
                f"(converged={result.converged}, imaginary={imaginary_count})"
            )
            return QCOptimizationResult(
                optimized_xyz=result.output_file if result.output_file else xyz_file,
                opt_energy=result.energy,
                converged=False,
                frequencies=result.frequencies,
                imaginary_count=imaginary_count,
                method_used="Berny_Failed"
            )

    def _try_ts_optimization_orca(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int,
        spin: int,
        old_checkpoint: Optional[Path]
    ) -> QCOptimizationResult:
        """ORCA TS 优化 (OptTS + Calc_Hess)"""
        # ORCA 使用内置的 ts_optimization 方法
        if not hasattr(self.qc_engine, 'ts_optimization'):
            # 简化实现：使用通用优化方法（无频率）
            self.logger.warning("ORCA 接口不支持 ts_optimization，使用简化模式")
            return QCOptimizationResult(
                optimized_xyz=xyz_file,
                opt_energy=0.0,
                converged=False,
                frequencies=None,
                imaginary_count=0,
                method_used="ORCA_Simplified"
            )

        opt_config = OptimizationConfig()
        result = self.qc_engine.ts_optimization(
            xyz_file,
            output_dir / "orcasts",
            opt_config=opt_config,
            timeout=None,
            charge=charge,
            spin=spin
        )

        imaginary_count = self._count_imaginary(result.frequencies)

        if result.converged and imaginary_count == 1:
            imaginary_freq = float(np.min(result.frequencies)) if result.frequencies is not None else 0.0
            self.logger.info(f"✓ ORCA OptTS 成功，虚频 = {imaginary_freq:.1f} cm⁻¹")
            return QCOptimizationResult(
                optimized_xyz=result.output_file if result.output_file else xyz_file,
                opt_energy=result.energy,
                converged=True,
                frequencies=result.frequencies,
                imaginary_count=1,
                method_used="ORCA_OptTS"
            )
        else:
            self.logger.warning(
                f"ORCA OptTS 失败或虚频不符 "
                f"(converged={result.converged}, imaginary={imaginary_count})"
            )
            return QCOptimizationResult(
                optimized_xyz=result.output_file if result.output_file else xyz_file,
                opt_energy=result.energy,
                converged=False,
                frequencies=result.frequencies,
                imaginary_count=imaginary_count,
                method_used="ORCA_OptTS_Failed"
            )

    def _try_ts_rescue(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int,
        spin: int,
        old_checkpoint: Optional[Path]
    ) -> QCOptimizationResult:
        """尝试 TS 救援 (Recalc=5 + NoEigenTest + MaxStep=10)"""
        self.logger.info("救援: Recalc=5 + NoEigenTest + MaxStep=10")

        # 创建救援配置
        rescue_config = OptimizationConfig()
        rescue_config.update_for_rescue({
            'recalc_hess_every': 5,
            'max_step': 10,
            'trust_radius': 0.15
        })
        rescue_config.ts_eigentest = False

        route = rescue_config.to_gaussian_route(
            self.method, self.basis, self.dispersion, is_ts=True
        )
        route += " Freq"

        result = self.qc_engine.optimize(
            xyz_file,
            output_dir / "rescue",
            route=route,
            old_checkpoint=old_checkpoint
        )

        imaginary_count = self._count_imaginary(result.frequencies)

        if result.converged and imaginary_count == 1:
            imaginary_freq = float(np.min(result.frequencies)) if result.frequencies is not None else 0.0
            self.logger.info(f"✓ TS 救援成功，虚频 = {imaginary_freq:.1f} cm⁻¹")
            return QCOptimizationResult(
                optimized_xyz=result.output_file if result.output_file else xyz_file,
                opt_energy=result.energy,
                converged=True,
                frequencies=result.frequencies,
                imaginary_count=1,
                method_used="TS_Rescue"
            )
        else:
            self.logger.error(f"TS 救援失败 (converged={result.converged}, imaginary={imaginary_count})")
            return QCOptimizationResult(
                optimized_xyz=result.output_file if result.output_file else xyz_file,
                opt_energy=result.energy,
                converged=False,
                frequencies=result.frequencies,
                imaginary_count=imaginary_count,
                method_used="TS_Rescue_Failed"
            )

    def _run_l2_sp(self, xyz_file: Path, output_dir: Path) -> Optional[float]:
        """
        运行 L2 高精度单点能

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录

        Returns:
            L2 能量 (Hartree)，失败返回 None
        """
        try:
            result = self.orca.single_point(xyz_file, output_dir)

            if not result.converged:
                self.logger.error(f"L2 SP 未收敛: {result.error_message}")
                return None

            self.logger.info(f"✓ L2 SP 成功: {result.energy:.8f} Hartree")
            return result.energy

        except Exception as e:
            self.logger.error(f"L2 SP 计算失败: {e}")
            return None

    def run_sp_only(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int = 0,
        spin: int = 1
    ) -> QCSPResult:
        """
        仅执行单点能计算（不进行优化）

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录
            charge: 分子电荷
            spin: 自旋多重度

        Returns:
            QCSPResult 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"执行单点能计算: {xyz_file.name}")

        try:
            result = self.orca.single_point(xyz_file, output_dir)

            return QCSPResult(
                energy=result.energy,
                converged=result.converged,
                output_file=result.output_file,
                error_message=result.error_message
            )

        except Exception as e:
            self.logger.error(f"SP 计算失败: {e}")
            return QCSPResult(
                energy=0.0,
                converged=False,
                error_message=str(e)
            )
