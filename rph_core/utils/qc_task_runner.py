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
from typing import Any, Dict, Optional, Tuple, cast
from dataclasses import dataclass
import re
import shutil
import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.optimization_config import OptimizationConfig, build_gaussian_route_from_config
from rph_core.utils.keyword_translator import KeywordTranslator
from rph_core.utils.data_types import QCResult
from rph_core.utils.qc_interface import (
    GaussianInterface,
    QCInterfaceFactory
)
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
    l2_sp_result: Optional[QCSPResult] = None  # L2 SP provenance/artifacts
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
    failure_kind: Optional[str] = None
    rescue_start: Optional[Path] = None
    rescue_route: Optional[str] = None


TS_FAILURE_OSCILLATION = "oscillation_or_step_exceeded"
TS_FAILURE_NO_IMAG = "converged_to_minimum_no_imag"
TS_FAILURE_MULTI_IMAG = "wrong_saddle_multi_imag"
TS_FAILURE_FATAL = "gaussian_fatal_or_no_geometry"
TS_FAILURE_UNKNOWN = "unknown_ts_failure"


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

    def __init__(
        self,
        config: dict[str, Any],
        resource_override: Optional[Dict[str, Any]] = None,
    ):
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

        # 资源配置 - 统一从 resources 读取
        resources = config.get('resources', {})
        self.nprocshared = resources.get('nproc', 16)
        self.mem = resources.get('mem', '32GB')

        # Apply resource overrides for intra-reaction parallel scheduling
        if resource_override:
            if 'nproc' in resource_override:
                self.nprocshared = resource_override['nproc']
            if 'mem' in resource_override:
                self.mem = resource_override['mem']

        # 优化配置
        self.method = self.theory_opt.get('method', 'B3LYP')
        self.basis = self.theory_opt.get('basis', 'def2-SVP')
        self.dispersion = self.theory_opt.get('dispersion', 'GD3BJ')
        self.engine_type = self.theory_opt.get('engine', 'gaussian').lower()
        self.sp_engine_type = self.theory_sp.get('engine', 'orca').lower()
        self.opt_method_spec = self.theory_opt.get('normalized_spec') if isinstance(self.theory_opt.get('normalized_spec'), dict) else None
        self.sp_method_spec = self.theory_sp.get('normalized_spec') if isinstance(self.theory_sp.get('normalized_spec'), dict) else None

        # 默认优化配置
        self.opt_config = OptimizationConfig.from_config(config)
        self.qc_engine: Any

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
            self.qc_engine = QCInterfaceFactory.create_interface(
                'orca',
                charge=0,
                multiplicity=1,
                nprocshared=self.nprocshared,
                mem=self.mem,
                method=self.theory_opt.get('method', self.method),
                basis=self.theory_opt.get('basis', self.basis),
                aux_basis=self.theory_opt.get('aux_basis', 'def2/J'),
                maxcore=self.theory_opt.get('maxcore', 4000),
                solvent=self.theory_opt.get('solvent', 'acetone'),
                route_extras=self.theory_opt.get('route_extras', ''),
                normalized_spec=self.opt_method_spec,
                config=self.config,
            )
        else:
            raise ValueError(f"不支持的引擎: {self.engine_type}")

        self.sp_engine: Any = QCInterfaceFactory.create_interface(
            self.sp_engine_type,
            charge=0,
            multiplicity=1,
            nprocshared=self.nprocshared,
            mem=self.mem,
            method=self.theory_sp.get('method', 'WB97M-V'),
            basis=self.theory_sp.get('basis', 'def2-TZVPP'),
            aux_basis=self.theory_sp.get('aux_basis', 'def2/J'),
            maxcore=self.theory_sp.get('maxcore', 4000),
            solvent=self.theory_sp.get('solvent', 'acetone'),
            route_extras=self.theory_sp.get('route_extras', ''),
            normalized_spec=self.sp_method_spec,
            config=self.config,
        )

        self.normal_route = None
        self.normal_rescue_route = None
        self.ts_route = None
        self.ts_rescue_route = None

        self.logger.info(f"QCTaskRunner 初始化: {self.engine_type} {self.method}/{self.basis}")
        self.logger.info(
            f"  L2 SP: {self.sp_engine_type} "
            f"{self.theory_sp.get('method', 'WB97M-V')}/{self.theory_sp.get('basis', 'def2-TZVPP')}"
        )

    @classmethod
    def create_with_override(
        cls,
        config: Dict[str, Any],
        resource_override: Optional[Dict[str, Any]] = None,
    ) -> "QCTaskRunner":
        """Create QCTaskRunner with optional resource override for parallel scheduling."""
        return cls(config=config, resource_override=resource_override)

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
    ) -> Tuple[Path, bool]:
        if coords is None or coords.size == 0:
            self.logger.warning(
                "Optimized coordinates missing or empty for %s; geometry file will point to input %s",
                name,
                source_xyz,
            )
            return source_xyz, False
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        optimized_xyz = output_dir / f"{name}.xyz"
        symbols = None
        if source_xyz.suffix.lower() == ".xyz":
            _, symbols = read_xyz(source_xyz)
        if symbols is not None and coords.shape[0] != len(symbols):
            self.logger.error(
                "原子数不一致: optimized coords %d vs source xyz %d symbols. "
                "回退使用原始结构以避免下游损坏。source=%s",
                coords.shape[0], len(symbols), source_xyz,
            )
            return source_xyz, False
        if symbols is None:
            n_atoms = coords.shape[0]
            symbols = ["X"] * n_atoms
        write_xyz(optimized_xyz, coords, symbols, title=f"{name} optimized")
        return optimized_xyz, True

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
        if not self._route_has_method_basis(route):
            route = route.strip()
            if route.startswith("#"):
                route = route.lstrip("#").strip()
                if route.lower().startswith("p "):
                    route = route[1:].strip()
            route = f"{self._build_ts_route_prefix()} {route}".strip()
        return self._ensure_ts_force_constants(route)

    def _route_has_method_basis(self, route: str) -> bool:
        pattern = r"\b[A-Za-z0-9][A-Za-z0-9+_.-]*\s*/\s*[A-Za-z0-9][A-Za-z0-9+_.()\-]*"
        return re.search(pattern, route) is not None

    def _build_ts_route_prefix(self) -> str:
        method = self.theory_opt.get('method', 'B3LYP')
        basis = KeywordTranslator.to_gaussian_basis(self.theory_opt.get('basis', 'def2SVP'))
        dispersion = KeywordTranslator.to_gaussian_dispersion(self.theory_opt.get('dispersion'))
        solvent = KeywordTranslator.to_gaussian_solvent(self.theory_opt.get('solvent'))
        parts = [f"#p {method}/{basis}"]
        if dispersion:
            parts.append(dispersion)
        if solvent:
            parts.append(solvent.strip())
        return " ".join(parts)

    def _ensure_ts_force_constants(self, route: str) -> str:
        normalized = route
        upper_route = normalized.upper()
        if "OPT" not in upper_route or "TS" not in upper_route:
            return normalized
        # RecalcFC alone is NOT sufficient for Opt=TS — Gaussian requires an
        # initial Hessian via CalcFC/CalcAll/ReadFC or a ModRedundant section.
        # Use word-boundary match to avoid "CalcFC" matching inside "RecalcFC".
        sufficient_tokens = ("CALCFC", "CALCALL", "READFC", "MODREDUNDANT")
        if any(re.search(r'\b' + token + r'\b', upper_route) for token in sufficient_tokens):
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
...
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"[QCTaskRunner] === Normal 模式优化: {xyz_file.name} ===")

        # 【NEW】Step 0: 几何预处理 (检测重叠 → xTB 预优化)

        from rph_core.utils.geometry_preprocessor import GeometryPreprocessor
        
        preprocessor = GeometryPreprocessor(self.config)
        preopt_result = preprocessor.preprocess(
            xyz_file=xyz_file,
            output_dir=output_dir / "preoptimization",
            charge=charge,
            uhf=spin - 1
        )
        
        # 如果预优化成功，使用预优化后的结构；否则继续使用原始结构
        if preopt_result.success and preopt_result.coordinates:
            working_xyz = preopt_result.coordinates
            # 通过 error_message 判断预优化状态
            status = preopt_result.error_message or ""
            if status == "preopt_success":
                self.logger.info(f"✓ xTB 预优化已完成，将使用预优化结构: {working_xyz}")
            elif status in ("no_overlap", "preopt_disabled"):
                self.logger.info(f"无需预优化（{status}），继续使用原始结构: {working_xyz}")
            elif status.startswith("xtb_preopt_failed"):
                self.logger.warning(f"xTB 预优化失败但继续：{status}")
            else:
                self.logger.info(f"预优化状态: {status}, 使用结构: {working_xyz}")
        else:
            working_xyz = xyz_file
            self.logger.warning(
                f"预优化失败或被跳过，将直接使用原始结构进行 DFT 优化: {working_xyz}"
            )

        self.normal_route = self._get_normal_route()
        self.normal_rescue_route = self._get_normal_route(rescue=True)

        # 尝试标准优化（使用预处理后的结构）
        opt_result = self._try_normal_optimization(
            working_xyz, output_dir, charge, spin, old_checkpoint, enable_nbo
        )

        # 如果标准优化失败，尝试救援
        if not opt_result.converged:
            if opt_result.error_message and opt_result.error_message.startswith("FATAL:"):
                self.logger.error(f"Fatal error detected; skipping rescue: {opt_result.error_message}")
                return opt_result
            self.logger.warning("标准优化失败，启动救援策略...")
            opt_result = self._try_normal_rescue(
                working_xyz, output_dir, charge, spin, old_checkpoint, enable_nbo
            )

        # L2 高精度 SP
        l2_energy = None
        if enable_l2_sp and opt_result.converged:
            self.logger.info("[QCTaskRunner] 执行 L2 高精度单点能...")

            l2_result = self._run_l2_sp(
                opt_result.optimized_xyz,
                output_dir / "L2_SP",
                charge=charge,
                spin=spin
            )
            opt_result.l2_sp_result = l2_result
            opt_result.l2_energy = l2_result.energy if l2_result.converged else None

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
        """Execute TS optimization + frequency validation + L2 SP workflow.

        This method orchestrates the TS optimization flow:
        1. Attempt standard TS optimization (Berny for Gaussian, OptTS for ORCA)
        2. If fails, attempt TS rescue (QST2 or modified strategy)
        3. If converged and L2 SP requested, run high-level single-point

        Args:
            xyz_file: Input XYZ file with TS guess geometry
            output_dir: Directory for all output files
            charge: Molecular charge (default 0)
            spin: Spin multiplicity (default 1 for singlet)
            enable_l2_sp: Whether to run L2 single-point calculation
            old_checkpoint: Optional checkpoint file from previous calculation

        Returns:
            QCOptimizationResult with converged TS geometry and optional L2 SP result
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"[QCTaskRunner] === TS optimization cycle: {xyz_file.name} ===")

        # Step 1: Attempt TS optimization
        ts_output_dir = output_dir / "ts_opt"
        ts_output_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info("[QCTaskRunner] Attempting TS optimization (primary)...")
        ts_result = self._try_ts_optimization(xyz_file, ts_output_dir, charge, spin, old_checkpoint)

        # Step 2: If primary TS attempt failed, try rescue from S2 initial guess
        if not ts_result.converged:
            policy = self.config.get("step3", {}).get("ts_rescue_policy", {})
            if policy.get("enabled", True):
                failure_kind = self._classify_ts_failure(ts_result, xyz_file)
                ts_result.failure_kind = failure_kind
                self.logger.warning(f"[QCTaskRunner] Primary TS failed, classified as: {failure_kind}")
                if failure_kind == TS_FAILURE_FATAL:
                    self.logger.error("[QCTaskRunner] Fatal TS failure; skipping rescue.")
                    return ts_result
                rescue_xyz = self._prepare_ts_rescue_start(
                    result=ts_result,
                    original_xyz=xyz_file,
                    output_dir=ts_output_dir,
                    failure_kind=failure_kind,
                )
                ts_result.rescue_start = rescue_xyz
                ts_result = self._try_ts_rescue(
                    rescue_xyz,
                    ts_output_dir,
                    charge,
                    spin,
                    old_checkpoint,
                    failure_kind=failure_kind,
                )
            else:
                self.logger.warning("[QCTaskRunner] Primary TS failed and rescue policy is disabled. Returning failed result.")
                return ts_result

        # Step 3: If we have converged TS and L2 SP requested, run L2 SP and attach
        if enable_l2_sp and ts_result.converged:
            self.logger.info("[QCTaskRunner] Running L2 single-point on TS geometry...")
            try:
                l2_result = self._run_l2_sp(ts_result.optimized_xyz, ts_output_dir / "L2_SP", charge=charge, spin=spin)
                ts_result.l2_sp_result = l2_result
                ts_result.l2_energy = l2_result.energy if l2_result.converged else None
            except Exception as e:
                self.logger.error(f"L2 SP execution raised exception: {e}")
                ts_result.l2_sp_result = None
                ts_result.l2_energy = None

        return ts_result

    def run_ts_rescue_only(
        self,
        rescue_start_xyz: Path,
        output_dir: Path,
        failure_kind: str,
        charge: int = 0,
        spin: int = 1,
        old_checkpoint: Optional[Path] = None,
        enable_l2_sp: bool = True,
    ) -> QCOptimizationResult:
        """跳过主 Berny TS 优化，直接进入救援阶段。

        适用场景：主 Berny 已跑过且明确失败（如 NStep 超限 / 收敛至极小值），
        无需浪费 CPU 重复主优化。救援路由来自 config 中的 ts_rescue 配置。

        Args:
            rescue_start_xyz: 救援起始结构 XYZ（由旧日志推断的起点）
            output_dir: S3 输出目录（通常是 S3_TS/）
            failure_kind: 失败分类常量（oscillation_or_step_exceeded 等）
            charge: 分子电荷
            spin: 自旋多重度
            old_checkpoint: 可选的旧 checkpoint 文件
            enable_l2_sp: 救援成功后是否执行 L2 高精度单点能

        Returns:
            QCOptimizationResult，若救援成功则 converged=True
        """
        ts_output_dir = output_dir / "ts_opt"
        ts_output_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(
            "[QCTaskRunner] Rescue-only mode: skipping primary TS opt, "
            "directly invoking Berny rescue for failure kind: %s",
            failure_kind,
        )

        rescue_result = self._try_ts_rescue(
            rescue_start_xyz,
            ts_output_dir,
            charge,
            spin,
            old_checkpoint,
            failure_kind=failure_kind,
        )

        if enable_l2_sp and rescue_result.converged:
            self.logger.info("[QCTaskRunner] Rescue succeeded; running L2 single-point on TS geometry...")
            try:
                l2_result = self._run_l2_sp(
                    rescue_result.optimized_xyz,
                    ts_output_dir / "L2_SP",
                    charge=charge,
                    spin=spin,
                )
                rescue_result.l2_sp_result = l2_result
                rescue_result.l2_energy = l2_result.energy if l2_result.converged else None
            except Exception as exc:
                self.logger.error("L2 SP after rescue raised exception: %s", exc)
                rescue_result.l2_sp_result = None
                rescue_result.l2_energy = None

        return rescue_result

    def _count_imaginary(self, frequencies: Optional[np.ndarray]) -> int:
        if frequencies is None:
            return 0
        return int(np.sum(frequencies < 0))

    def _parse_gaussian_opt_step_count(self, log_path: Optional[Path]) -> int:
        if not log_path or not log_path.exists():
            return 0
        try:
            text = log_path.read_text(errors="ignore")
            match = re.search(r"NStep=\s*(\d+)", text)
            if match:
                return int(match.group(1))
            return text.count("GradGradGradGrad")
        except OSError as exc:
            self.logger.debug("Could not read Gaussian optimization step count from %s: %s", log_path, exc)
            return 0

    def _same_path(self, first: Path, second: Path) -> bool:
        try:
            return first.resolve() == second.resolve()
        except OSError:
            return first == second

    def _classify_ts_failure(
        self,
        result: "QCOptimizationResult",
        original_xyz: Path,
    ) -> str:
        _ = original_xyz
        imag = result.imaginary_count
        log_path = result.log_file or result.freq_log or result.qm_output_file
        step_count = self._parse_gaussian_opt_step_count(log_path)
        policy = self.config.get("step3", {}).get("ts_rescue_policy", {})
        trigger_steps = policy.get("trigger_after_steps", 60)

        if log_path and log_path.exists():
            try:
                tail = "\n".join(log_path.read_text(errors="ignore").splitlines()[-200:])
                if "Number of steps exceeded" in tail and step_count >= trigger_steps:
                    return TS_FAILURE_OSCILLATION
            except OSError as exc:
                self.logger.debug("Could not inspect Gaussian TS failure log %s: %s", log_path, exc)

        error_text = result.error_message or ""
        if imag == 0 and (
            result.converged
            or result.method_used in {
                "Berny_ConvergedToMinimum",
                "ORCA_OptTS_ConvergedToMinimum",
            }
            or "0 imaginary" in error_text
            or "converged to minimum" in error_text.lower()
            or "虚频数为 0" in error_text
        ):
            return TS_FAILURE_NO_IMAG
        if imag > 1:
            return TS_FAILURE_MULTI_IMAG
        if result.error_message and result.error_message.startswith("FATAL"):
            return TS_FAILURE_FATAL
        return TS_FAILURE_UNKNOWN

    def _prepare_ts_rescue_start(
        self,
        result: "QCOptimizationResult",
        original_xyz: Path,
        output_dir: Path,
        failure_kind: str,
    ) -> Path:
        """Always restart rescue from S2 initial guess, not from failed-run geometry.

        The failed optimization may have fallen into a product well (narrow barrier);
        restarting from the last geometry would perpetuate the failure.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        rescue_xyz = output_dir / "rescue_start_from_s2_guess.xyz"
        shutil.copy2(original_xyz, rescue_xyz)
        self.logger.info(
            "Rescue restart from S2 initial guess (failure_kind=%s): %s",
            failure_kind, rescue_xyz,
        )
        return rescue_xyz

    def _coalesce_qc_error(self, result: Any) -> Optional[str]:
        """Prefer explicit error_message; otherwise point to output/log."""

        msg = getattr(result, "error_message", None)
        if msg:
            return str(msg)
        out = getattr(result, "output_file", None) or getattr(result, "log_file", None)
        if out:
            return f"See log: {out}"
        return None

    def _to_sp_result(self, result: QCResult) -> QCSPResult:
        qm_output_file = (
            getattr(result, 'qm_output_file', None)
            or getattr(result, 'log_file', None)
            or getattr(result, 'output_file', None)
        )
        output_file = getattr(result, 'output_file', None) or qm_output_file
        log_file = getattr(result, 'log_file', None) or qm_output_file
        energy = float(result.energy) if result.energy is not None else 0.0
        converged = bool(getattr(result, 'converged', False) and result.energy is not None)
        return QCSPResult(
            energy=energy,
            converged=converged,
            output_file=output_file,
            error_message=getattr(result, 'error_message', None),
            log_file=log_file,
            chk_file=getattr(result, 'chk_file', None),
            fchk_file=getattr(result, 'fchk_file', None),
            qm_output_file=qm_output_file,
        )

    def _execute_single_point(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int = 0,
        spin: int = 1,
    ) -> QCSPResult:
        try:
            sp_engine = cast(Any, self.sp_engine)
            if not hasattr(sp_engine, 'single_point'):
                raise AttributeError(f"{type(sp_engine).__name__} does not support single_point()")
            result = sp_engine.single_point(
                xyz_file,
                output_dir,
                charge=charge,
                spin=spin,
            )
            return self._to_sp_result(result)
        except Exception as e:
            return QCSPResult(
                energy=0.0,
                converged=False,
                error_message=str(e),
            )

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
            ]
            for pat in fatal_patterns:
                if pat in tail:
                    return True, f"FATAL: {pat} ({log_path})"
            if re.search(r'Error termination via Lnk1e\b(?!.*\bl9999\b)', tail):
                return True, f"FATAL: Error termination via Lnk1e ({log_path})"
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
        self.logger.info("[QCTaskRunner] 尝试标准基态优化...")

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
            old_checkpoint=old_checkpoint,
            charge=charge,
            spin=spin,
        )

        # 验证结果
        imaginary_count = self._count_imaginary(result.frequencies)

        if result.converged and imaginary_count == 0:
            self.logger.info("✓ 标准优化成功，无虚频")
            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz, coords_ok = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "standard",
                f"{xyz_file.stem}_opt"
            )
            error_message = None if coords_ok else f"Optimized coordinates not available for {xyz_file.stem}_opt"
            

            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=coords_ok,
                frequencies=result.frequencies,
                imaginary_count=0,
                method_used="Normal",
                error_message=error_message,
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
            optimized_xyz, _coords_ok = self._write_optimized_xyz(
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
        self.logger.info("[QCTaskRunner] 救援: CalcFC + MaxStep=10")

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
            optimized_xyz, coords_ok = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "rescue",
                f"{xyz_file.stem}_opt_rescue"
            )
            error_message = None if coords_ok else f"Optimized coordinates not available for {xyz_file.stem}_opt_rescue"
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=coords_ok,
                frequencies=result.frequencies,
                imaginary_count=0,
                method_used="Normal_Rescue",
                error_message=error_message,
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
            optimized_xyz, _coords_ok = self._write_optimized_xyz(
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
        self.logger.info(f"[S3] 尝试标准 TS 优化 (引擎: {self.engine_type})...")

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
            old_checkpoint=old_checkpoint,
            charge=charge,
            spin=spin,
        )

        imaginary_count = self._count_imaginary(result.frequencies)

        if result.converged and imaginary_count == 1:
            imaginary_freq = float(np.min(result.frequencies)) if result.frequencies is not None else 0.0
            self.logger.info(f"✓ Berny TS 优化成功，虚频 = {imaginary_freq:.1f} cm⁻¹")
            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz, coords_ok = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "berny",
                f"{xyz_file.stem}_ts"
            )
            error_message = None if coords_ok else f"Optimized coordinates not available for {xyz_file.stem}_ts"
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=coords_ok,
                frequencies=result.frequencies,
                imaginary_count=1,
                method_used="Berny",
                error_message=error_message,
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
            if is_fatal:
                method_used = "Berny_Fatal"
                error_message = fatal_msg
            elif result.converged and imaginary_count == 0:
                method_used = "Berny_ConvergedToMinimum"
                error_message = "TS optimization converged to minimum: 0 imaginary frequencies"
            elif result.converged and imaginary_count > 1:
                method_used = "Berny_WrongSaddle"
                error_message = f"TS optimization found high-order saddle: {imaginary_count} imaginary frequencies"
            else:
                if log_path and log_path.exists():
                    try:
                        tail_text = log_path.read_text(errors="ignore")
                        if "Number of steps exceeded" in tail_text:
                            method_used = "Berny_StepExceeded"
                        else:
                            method_used = "Berny_Failed"
                    except Exception:
                        method_used = "Berny_Failed"
                else:
                    method_used = "Berny_Failed"
                error_message = self._coalesce_qc_error(result)

            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz, _coords_ok = self._write_optimized_xyz(
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
                method_used=method_used,
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

        opt_config = self.opt_config
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
            # Prefer the ORCA-resolved final .xyz path; fallback to parsed coords
            final_xyz_str = getattr(result, 'final_xyz_path', None)
            if final_xyz_str and Path(final_xyz_str).exists():
                optimized_xyz = Path(final_xyz_str).resolve()
                coords_ok = True
            else:
                coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
                optimized_xyz, coords_ok = self._write_optimized_xyz(
                    coords, xyz_file, output_dir / "orcasts", f"{xyz_file.stem}_ts"
                )
            error_message = None if coords_ok else f"Optimized coordinates not available for {xyz_file.stem}_ts"
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=coords_ok,
                frequencies=result.frequencies,
                imaginary_count=1,
                method_used="ORCA_OptTS",
                error_message=error_message,
                freq_log=result.output_file,
                log_file=getattr(result, 'log_file', None),
                chk_file=getattr(result, 'chk_file', None),
                fchk_file=getattr(result, 'fchk_file', None),
                qm_output_file=getattr(result, 'qm_output_file', None)
            )
        else:
            self.logger.warning(
                f"ORCA OptTS 失败 (%s): converged=%s, imaginary=%s",
                xyz_file.name, result.converged, imaginary_count
            )
            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz, _coords_ok = self._write_optimized_xyz(
                coords,
                xyz_file,
                output_dir / "orcasts",
                f"{xyz_file.stem}_ts_failed"
            )

            # 区分失效模式：ORCA 几何收敛状态与虚频数独立判定
            if result.converged and imaginary_count == 0:
                method_used = "ORCA_OptTS_ConvergedToMinimum"
                error_message = "几何已收敛但虚频数为 0（找到极小值点，非过渡态）"
            elif result.converged and imaginary_count > 1:
                method_used = "ORCA_OptTS_WrongSaddlePoint"
                error_message = f"几何已收敛但虚频数={imaginary_count}（期望 1 个，找到非目标鞍点）"
            else:
                method_used = "ORCA_OptTS_Failed"
                error_message = self._coalesce_qc_error(result)

            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=False,
                frequencies=result.frequencies,
                imaginary_count=imaginary_count,
                method_used=method_used,
                error_message=error_message,
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
        old_checkpoint: Optional[Path],
        failure_kind: str = TS_FAILURE_UNKNOWN,
    ) -> QCOptimizationResult:
        """尝试 TS 救援 (RecalcFC=5 + NoEigenTest + MaxStep=15)"""
        self.logger.info(f"[S3] Berny rescue for failure kind: {failure_kind}")

        route = self.ts_rescue_route if hasattr(self, "ts_rescue_route") and self.ts_rescue_route else self._get_ts_route(rescue=True)
        if "freq" not in route.lower():
            route = f"{route} Freq"
        route = self._ensure_ts_force_constants(route)
        rescue_dir = output_dir / "rescue" / failure_kind
        rescue_dir.mkdir(parents=True, exist_ok=True)

        result = self.qc_engine.optimize(
            xyz_file,
            rescue_dir,
            route=route,
            old_checkpoint=old_checkpoint,
            charge=charge,
            spin=spin,
        )

        import json
        ts_rescue_metadata = {
            "failure_kind": failure_kind,
            "rescue_start": str(xyz_file),
            "rescue_route": route,
        }
        try:
            with open(rescue_dir / "ts_rescue_metadata.json", "w") as f:
                json.dump(ts_rescue_metadata, f, indent=2)
        except OSError as exc:
            self.logger.warning("Could not write TS rescue metadata in %s: %s", rescue_dir, exc)

        imaginary_count = self._count_imaginary(result.frequencies)

        if result.converged and imaginary_count == 1:
            imaginary_freq = float(np.min(result.frequencies)) if result.frequencies is not None else 0.0
            self.logger.info(f"✓ TS 救援成功，虚频 = {imaginary_freq:.1f} cm⁻¹")
            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz, coords_ok = self._write_optimized_xyz(
                coords,
                xyz_file,
                rescue_dir,
                f"{xyz_file.stem}_ts_rescue"
            )
            error_message = None if coords_ok else f"Optimized coordinates not available for {xyz_file.stem}_ts_rescue"
            return QCOptimizationResult(
                optimized_xyz=optimized_xyz,
                opt_energy=result.energy,
                converged=coords_ok,
                frequencies=result.frequencies,
                imaginary_count=1,
                method_used="TS_Rescue",
                error_message=error_message,
                freq_log=result.output_file,
                log_file=getattr(result, 'log_file', None),
                chk_file=getattr(result, 'chk_file', None),
                fchk_file=getattr(result, 'fchk_file', None),
                qm_output_file=getattr(result, 'qm_output_file', None),
                failure_kind=failure_kind,
                rescue_start=xyz_file,
                rescue_route=route,
            )
        else:
            self.logger.error(f"TS 救援失败 (converged={result.converged}, imaginary={imaginary_count})")

            log_path = getattr(result, 'log_file', None) or result.output_file
            is_fatal, fatal_msg = self._analyze_log_for_fatal_errors(log_path)
            error_message = fatal_msg if is_fatal else self._coalesce_qc_error(result)

            coords = result.coordinates if isinstance(result.coordinates, np.ndarray) else None
            optimized_xyz, _coords_ok = self._write_optimized_xyz(
                coords,
                xyz_file,
                rescue_dir,
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
                qm_output_file=getattr(result, 'qm_output_file', None),
                failure_kind=failure_kind,
                rescue_start=xyz_file,
                rescue_route=route,
            )

    def _run_l2_sp(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int = 0,
        spin: int = 1
    ) -> QCSPResult:
        """
        运行 L2 高精度单点能

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录

        Returns:
            QCSPResult (includes energy and provenance)
        """
        result = self._execute_single_point(
            xyz_file,
            output_dir,
            charge=charge,
            spin=spin,
        )

        if not result.converged:
            self.logger.error(f"L2 SP 未收敛: {result.error_message}")
            return result

        self.logger.info(f"[QCTaskRunner] ✓ L2 SP 成功: {result.energy:.8f} Hartree")
        return result

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

        self.logger.info(f"[QCTaskRunner] 执行单点能计算: {xyz_file.name}")
        result = self._execute_single_point(
            xyz_file,
            output_dir,
            charge=charge,
            spin=spin,
        )
        if not result.converged:
            self.logger.error(f"SP 计算失败: {result.error_message}")
        return result
