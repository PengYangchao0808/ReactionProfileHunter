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

# pyright: ignore

import logging
from pathlib import Path
from typing import Any, Optional, Tuple, Union
from dataclasses import dataclass
import re
import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.optimization_config import OptimizationConfig, build_gaussian_route_from_config
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
    log_file: Optional[Path] = None  # Gaussian .log or ORCA .out
    chk_file: Optional[Path] = None
    fchk_file: Optional[Path] = None
    qm_output_file: Optional[Path] = None  # Generic QM output (Gaussian=.log, ORCA=.out)

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
    freq_log: Optional[Path] = None
    log_file: Optional[Path] = None  # Gaussian .log or ORCA .out
    chk_file: Optional[Path] = None
    fchk_file: Optional[Path] = None
    qm_output_file: Optional[Path] = None  # Generic QM output (Gaussian=.log, ORCA=.out)


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

    def _get_normal_route(self, rescue: bool = False) -> str:
        route = self.theory_opt.get('rescue_route') if rescue else self.theory_opt.get('route')
        if route:
            return route
        return build_gaussian_route_from_config(self.config, rescue=rescue)

    def _write_optimized_xyz(
        self,
        coords: Optional[np.ndarray],
        source_xyz: Path,
        output_dir: Path,
        name: str
    ) -> Path:
        if coords is None or coords.size == 0:
            return source_xyz
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        optimized_xyz = output_dir / f"{name}.xyz"
        symbols = None
        if source_xyz.suffix.lower() == ".xyz":
            _, symbols = read_xyz(source_xyz)
        if symbols is None:
            n_atoms = coords.shape[0]
            symbols = ["X"] * n_atoms
        write_xyz(optimized_xyz, coords, symbols, title=f"{name} optimized")
        return optimized_xyz

    def _get_ts_route(self, rescue: bool = False) -> str:
        step3_keywords = self.config.get('step3', {}).get('gaussian_keywords', {})
        if rescue:
            route = step3_keywords.get('ts_rescue')
            if not route:
                route = step3_keywords.get('berny')
        else:
            route = step3_keywords.get('berny')
        if not route:
            route = "# B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt=(TS, CalcFC, NoEigenTest) Freq"
        return self._ensure_ts_force_constants(route)

    def _ensure_ts_force_constants(self, route: str) -> str:
        normalized = route
        upper_route = normalized.upper()
        if "OPT" not in upper_route or "TS" not in upper_route:
            return normalized
        required_tokens = ("CALCFC", "CALCALL", "READFC", "MODREDUNDANT")
        if any(token in upper_route for token in required_tokens):
            return normalized
        opt_start = upper_route.find("OPT=")
        if opt_start == -1:
            return normalized
        paren_start = normalized.find("(", opt_start)
        if paren_start == -1 and "OPT=TS" in upper_route:
            normalized = normalized.replace("Opt=TS", "Opt=(TS,CalcFC)")
            normalized = normalized.replace("OPT=TS", "Opt=(TS,CalcFC)")
            self.logger.warning("TS 路由缺少力常数关键词，已自动补充 CalcFC")
            return normalized
        paren_end = normalized.find(")", paren_start + 1) if paren_start != -1 else -1
        if paren_start != -1 and paren_end != -1:
            opt_body = normalized[paren_start + 1:paren_end]
            opt_body = f"{opt_body},CalcFC"
            normalized = f"{normalized[:paren_start + 1]}{opt_body}{normalized[paren_end:]}"
        else:
            normalized = normalized.replace("Opt=TS", "Opt=(TS,CalcFC)")
            normalized = normalized.replace("OPT=TS", "Opt=(TS,CalcFC)")
        self.logger.warning("TS 路由缺少力常数关键词，已自动补充 CalcFC")
        return normalized

    def run_opt_sp_cycle(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int = 0,
        spin: int = 1,
        enable_l2_sp: bool = True,
        enable_nbo: bool = False,
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
            enable_nbo: 是否在 OPT+Freq 中集成 NBO (默认 Pop=NBO)
            old_checkpoint: 复用之前的 checkpoint

        Returns:
            QCOptimizationResult 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"=== Normal 模式优化: {xyz_file.name} ===")

        self.normal_route = self._get_normal_route()
        self.normal_rescue_route = self._get_normal_route(rescue=True)

        # 尝试标准优化
        opt_result = self._try_normal_optimization(
            xyz_file, output_dir, charge, spin, old_checkpoint, enable_nbo
        )

        # 如果标准优化失败，尝试救援
        if not opt_result.converged:
            if opt_result.error_message and opt_result.error_message.startswith("FATAL:"):
                self.logger.error(f"Fatal error detected; skipping rescue: {opt_result.error_message}")
                return opt_result
            self.logger.warning("标准优化失败，启动救援策略...")
            opt_result = self._try_normal_rescue(
                xyz_file, output_dir, charge, spin, old_checkpoint, enable_nbo
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

        self.ts_route = self._get_ts_route()
        self.ts_rescue_route = self._get_ts_route(rescue=True)

        # 尝试标准 TS 优化
        opt_result = self._try_ts_optimization(
            xyz_file, output_dir, charge, spin, old_checkpoint
        )

        # 如果 TS 优化失败，尝试救援
        if not opt_result.converged or opt_result.imaginary_count != 1:
            if opt_result.error_message and opt_result.error_message.startswith("FATAL:"):
                self.logger.error(f"Fatal error detected; skipping TS rescue: {opt_result.error_message}")
                return opt_result
            self.logger.warning("TS 优化失败或虚频不符，启动救援策略...")
            opt_result = self._try_ts_rescue(
                xyz_file, output_dir, charge, spin, old_checkpoint
            )

        l2_energy = None
        if enable_l2_sp and opt_result.converged:
            self.logger.info("执行 L2 高精度单点能...")
            l2_energy = self._run_l2_sp(opt_result.optimized_xyz, output_dir / "L2_SP")
            opt_result.l2_energy = l2_energy

        return opt_result

    def _count_imaginary(self, frequencies: Optional[np.ndarray]) -> int:
        if frequencies is None:
            return 0
        return int(np.sum(frequencies < 0))

    def _coalesce_qc_error(self, result: Any) -> Optional[str]:
        """Prefer explicit error_message; otherwise point to output/log."""

        msg = getattr(result, "error_message", None)
        if msg:
            return str(msg)
        out = getattr(result, "output_file", None) or getattr(result, "log_file", None)
        if out:
            return f"See log: {out}"
        return None

    def _analyze_log_for_fatal_errors(self, log_path: Optional[Path]) -> Tuple[bool, Optional[str]]:
        """Detect fatal Gaussian errors that rescue cannot fix.

        This prevents repeated OPT/rescue reruns for route/keyword/system errors
        (e.g., QPErr due to unsupported Pop= options), which can otherwise
        accumulate large scratch/output and stall WSL.
        """

        if log_path is None:
            return True, "FATAL: log file not available"

        try:
            if not log_path.exists():
                return True, f"FATAL: log file not found ({log_path})"

            tail = "\n".join(log_path.read_text(errors="ignore").splitlines()[-200:])
            fatal_patterns = [
                "QPErr",
                "Atomic number out of range",
                "Error termination via Lnk1e",
            ]
            for pat in fatal_patterns:
                if pat in tail:
                    return True, f"FATAL: {pat} ({log_path})"
            return False, None

        except Exception as exc:
            return True, f"FATAL: log read failed: {exc} ({log_path})"

    def _apply_nbo_route(self, route: str) -> str:
        """Enable NBO in a Gaussian route (G16 Rev A compatible).

        Default: Pop=NBO  (built-in NBO 3.1)
        If a keylist is provided via config, use Pop=(NBORead) and rely on
        GaussianInterface to inject the $NBO ... $END block.

        Note: Do not use Pop=NBO7 / Pop=NBO7Read on G16 Rev A.*; these options
        are introduced in Rev C and will trigger QPErr at Pop1.

        Also normalizes legacy invalid forms like Pop=NBO7Read,6D by removing
        the NBO7 tokens and splitting 6D out of Pop=.
        """

        normalized = re.sub(r"\s+", " ", (route or "").strip())
        lower = normalized.lower()

        # Normalize known-bad legacy forms from older revisions.
        normalized = re.sub(r"Pop\s*=\s*NBO7Read\s*,\s*6D", "6D Pop=(NBORead)", normalized, flags=re.IGNORECASE)
        normalized = re.sub(r"Pop\s*=\s*NBO7Read\b", "Pop=(NBORead)", normalized, flags=re.IGNORECASE)
        normalized = re.sub(r"Pop\s*=\s*NBO7\b", "Pop=NBO", normalized, flags=re.IGNORECASE)
        lower = normalized.lower()

        # If route already has a Pop= directive, don't try to outsmart it.
        if "pop=" in lower:
            return normalized

        nbo_keylist = (self.theory_opt.get("nbo_keylist") or self.config.get("nbo", {}).get("keylist"))
        if nbo_keylist:
            return f"{normalized} Pop=(NBORead)"
        return f"{normalized} Pop=NBO"

    def _try_normal_optimization(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int,
        spin: int,
        old_checkpoint: Optional[Path],
        enable_nbo: bool = False
    ) -> QCOptimizationResult:
        """尝试标准基态优化"""
        self.logger.info("尝试标准基态优化...")

        route = self.normal_route if hasattr(self, "normal_route") and self.normal_route else self._get_normal_route()
        if "freq" not in route.lower():
            route = f"{route} Freq"
        if enable_nbo:
            route = self._apply_nbo_route(route)

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
            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "standard",
                f"{xyz_file.stem}_opt"
            )
            

            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=True,
                frequencies=result.frequencies,
                imaginary_count=0,
                method_used="Normal",
                freq_log=result.output_file,
                log_file=result.log_file,
                chk_file=result.chk_file,
                fchk_file=result.fchk_file,
                qm_output_file=result.qm_output_file
            )
        else:
            self.logger.warning(
                f"标准优化失败 (converged={result.converged}, imaginary={imaginary_count})"
            )

            log_path = result.log_file or result.output_file
            is_fatal, fatal_msg = self._analyze_log_for_fatal_errors(log_path)
            error_message = fatal_msg if is_fatal else self._coalesce_qc_error(result)

            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "standard",
                f"{xyz_file.stem}_opt_failed"
            )
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=False,
                frequencies=result.frequencies,
                imaginary_count=imaginary_count,
                method_used="Normal_Failed" if not is_fatal else "Normal_Fatal",
                error_message=error_message,
                freq_log=result.output_file,
                log_file=result.log_file,
                chk_file=result.chk_file,
                fchk_file=result.fchk_file,
                qm_output_file=result.qm_output_file
            )

    def _try_normal_rescue(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int,
        spin: int,
        old_checkpoint: Optional[Path],
        enable_nbo: bool = False
    ) -> QCOptimizationResult:
        """尝试基态救援 (CalcFC + MaxStep=10)"""
        self.logger.info("救援: CalcFC + MaxStep=10")

        route = self.normal_rescue_route if hasattr(self, "normal_rescue_route") and self.normal_rescue_route else self._get_normal_route(rescue=True)
        if "freq" not in route.lower():
            route = f"{route} Freq"
        if enable_nbo:
            route = self._apply_nbo_route(route)

        result = self.qc_engine.optimize(
            xyz_file,
            output_dir / "rescue",
            route=route,
            old_checkpoint=old_checkpoint
        )

        imaginary_count = self._count_imaginary(result.frequencies)

        if result.converged and imaginary_count == 0:
            self.logger.info("✓ 救援成功")
            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "rescue",
                f"{xyz_file.stem}_opt_rescue"
            )
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=True,
                frequencies=result.frequencies,
                imaginary_count=0,
                method_used="Normal_Rescue",
                freq_log=result.output_file,
                log_file=result.log_file,
                chk_file=result.chk_file,
                fchk_file=result.fchk_file,
                qm_output_file=result.qm_output_file
            )
        else:
            self.logger.error(f"救援失败 (converged={result.converged}, imaginary={imaginary_count})")

            log_path = result.log_file or result.output_file
            is_fatal, fatal_msg = self._analyze_log_for_fatal_errors(log_path)
            error_message = fatal_msg if is_fatal else self._coalesce_qc_error(result)

            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "rescue",
                f"{xyz_file.stem}_opt_rescue_failed"
            )
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=False,
                frequencies=result.frequencies,
                imaginary_count=imaginary_count,
                method_used="Normal_Rescue_Failed" if not is_fatal else "Normal_Rescue_Fatal",
                error_message=error_message,
                freq_log=result.output_file,
                log_file=result.log_file,
                chk_file=result.chk_file,
                fchk_file=result.fchk_file,
                qm_output_file=result.qm_output_file
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
        route = self.ts_route if hasattr(self, "ts_route") and self.ts_route else self._get_ts_route()
        if "freq" not in route.lower():
            route = f"{route} Freq"
        route = self._ensure_ts_force_constants(route)

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
            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "berny",
                f"{xyz_file.stem}_ts"
            )
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=True,
                frequencies=result.frequencies,
                imaginary_count=1,
                method_used="Berny",
                freq_log=result.output_file,
                log_file=result.log_file,
                chk_file=result.chk_file,
                fchk_file=result.fchk_file,
                qm_output_file=result.qm_output_file
            )
        else:
            self.logger.warning(
                f"Berny TS 优化失败或虚频不符 "
                f"(converged={result.converged}, imaginary={imaginary_count})"
            )

            log_path = result.log_file or result.output_file
            is_fatal, fatal_msg = self._analyze_log_for_fatal_errors(log_path)
            error_message = fatal_msg if is_fatal else self._coalesce_qc_error(result)

            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "berny",
                f"{xyz_file.stem}_ts_failed"
            )
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=False,
                frequencies=result.frequencies,
                imaginary_count=imaginary_count,
                method_used="Berny_Failed" if not is_fatal else "Berny_Fatal",
                error_message=error_message,
                freq_log=result.output_file,
                log_file=result.log_file,
                chk_file=result.chk_file,
                fchk_file=result.fchk_file,
                qm_output_file=result.qm_output_file
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
        orca_ts_opt = getattr(self.qc_engine, "ts_optimization", None)
        if orca_ts_opt is None:
            # 简化实现：使用通用优化方法（无频率）
            self.logger.warning("ORCA 接口不支持 ts_optimization，使用简化模式")
            return QCOptimizationResult(
                optimized_xyz=xyz_file,
                opt_energy=0.0,
                converged=False,
                frequencies=None,
                imaginary_count=0,
                method_used="ORCA_Simplified",
                freq_log=None
            )

        opt_config = OptimizationConfig()
        result = orca_ts_opt(
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
            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "orcasts",
                f"{xyz_file.stem}_ts"
            )
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=True,
                frequencies=result.frequencies,
                imaginary_count=1,
                method_used="ORCA_OptTS",
                freq_log=result.output_file,
                log_file=getattr(result, 'log_file', None),
                chk_file=getattr(result, 'chk_file', None),
                fchk_file=getattr(result, 'fchk_file', None),
                qm_output_file=getattr(result, 'qm_output_file', None)
            )
        else:
            self.logger.warning(
                f"ORCA OptTS 失败或虚频不符 "
                f"(converged={result.converged}, imaginary={imaginary_count})"
            )
            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "orcasts",
                f"{xyz_file.stem}_ts_failed"
            )
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=False,
                frequencies=result.frequencies,
                imaginary_count=imaginary_count,
                method_used="ORCA_OptTS_Failed",
                error_message=self._coalesce_qc_error(result),
                freq_log=result.output_file,
                log_file=getattr(result, 'log_file', None),
                chk_file=getattr(result, 'chk_file', None),
                fchk_file=getattr(result, 'fchk_file', None),
                qm_output_file=getattr(result, 'qm_output_file', None)
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

        route = self.ts_rescue_route if hasattr(self, "ts_rescue_route") and self.ts_rescue_route else self._get_ts_route(rescue=True)
        if "freq" not in route.lower():
            route = f"{route} Freq"
        route = self._ensure_ts_force_constants(route)

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
            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "rescue",
                f"{xyz_file.stem}_ts_rescue"
            )
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=True,
                frequencies=result.frequencies,
                imaginary_count=1,
                method_used="TS_Rescue",
                freq_log=result.output_file,
                log_file=getattr(result, 'log_file', None),
                chk_file=getattr(result, 'chk_file', None),
                fchk_file=getattr(result, 'fchk_file', None),
                qm_output_file=getattr(result, 'qm_output_file', None)
            )
        else:
            self.logger.error(f"TS 救援失败 (converged={result.converged}, imaginary={imaginary_count})")

            log_path = getattr(result, 'log_file', None) or result.output_file
            is_fatal, fatal_msg = self._analyze_log_for_fatal_errors(log_path)
            error_message = fatal_msg if is_fatal else self._coalesce_qc_error(result)

            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "rescue",
                f"{xyz_file.stem}_ts_rescue_failed"
            )
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=False,
                frequencies=result.frequencies,
                imaginary_count=imaginary_count,
                method_used="TS_Rescue_Failed" if not is_fatal else "TS_Rescue_Fatal",
                error_message=error_message,
                freq_log=result.output_file,
                log_file=getattr(result, 'log_file', None),
                chk_file=getattr(result, 'chk_file', None),
                fchk_file=getattr(result, 'fchk_file', None),
                qm_output_file=getattr(result, 'qm_output_file', None)
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

            energy = result.energy if result.energy is not None else 0.0
            return QCSPResult(
                energy=energy,
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
