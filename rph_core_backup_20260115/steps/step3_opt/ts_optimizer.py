"""
Step 3: Transition Analyzer
==========================

反应中心全分析模块 - TS优化验证 + Reactant SP + 片段处理

Author: QCcalc Team
Date: 2026-01-09
Updated: 2026-01-13 (合并 S3.5，集成 SP 和片段处理)
"""

import logging
from pathlib import Path
from typing import Optional, Tuple, List
from dataclasses import dataclass
import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import write_xyz, read_xyz
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.utils.qc_interface import XTBInterface
from rph_core.utils.geometry_tools import LogParser
# v3.0: SPMatrixReport 已移除（批处理 SP 功能被 ConformerEngine 替代）
# 创建简化的占位符类以保持接口兼容性
from dataclasses import dataclass
from typing import Optional
@dataclass
class SPMatrixReport:
    """简化的 SP 报告（v3.0 占位符）"""
    e_ts: float = 0.0
    e_reactant: float = 0.0
    e_product: float = 0.0
    e_reactant_l2: float = 0.0
    e_product_l2: float = 0.0
    e_ts_final: float = 0.0  # TS 最终能量（从 l2_energy 赋值）
    # 片段能量（v3.0 可能不使用，占位符）
    e_frag_a_ts: float = 0.0
    e_frag_b_ts: float = 0.0
    e_frag_a_relaxed: float = 0.0
    e_frag_b_relaxed: float = 0.0
    # 方法和溶剂（占位符）
    method: str = "Berny"
    solvent: str = "acetone"

    def get_activation_energy(self) -> float:
        """返回活化能"""
        return (self.e_ts - self.e_reactant) * 627.509  # Hartree → kcal/mol
    def get_reaction_energy(self) -> float:
        """返回反应能"""
        return (self.e_product - self.e_reactant) * 627.509  # Hartree → kcal/mol

from rph_core.steps.step4_features.fragment_extractor import FragmentExtractor
from .berny_driver import BernyTSDriver, TSOptResult
from .qst2_rescue import QST2RescueDriver, QST2Result
from .validator import TSValidator, TSValidationError
from .irc_driver import IRCDriver, IRCResult


logger = logging.getLogger(__name__)


@dataclass
class TransitionAnalysisResult:
    """
    Transition Analyzer 完整结果 (S3 输出)

    Attributes:
        ts_final_xyz: 优化后的 TS 结构 XYZ 文件
        ts_checkpoint: TS 的 checkpoint 文件 (可选)
        sp_report: SP 矩阵报告 (SPMatrixReport)
        method_used: 使用的方法 (Berny 或 QST2)
    """
    ts_final_xyz: Path
    ts_checkpoint: Optional[Path]
    sp_report: SPMatrixReport
    method_used: str


class TSOptimizer(LoggerMixin):
    """
    TS 优化引擎 (Step 3) - v3.0

    策略:
    1. Berny TS 优化 (Opt=TS, CalcFC, NoEigenTest)
    2. 虚频检验 (必须恰好 1 个虚频)
    3. 失败 → QST2 救援
    4. 可选 IRC 验证

    输入:
    - ts_guess: TS 初猜 (来自 Step 2，可能是文件或目录)
    - reactant: 底物复合物 (来自 Step 2，用于 QST2 救援)
    - product: 产物 (来自 Step 1，用于 QST2 救援)

    输出:
    - ts_final_xyz: 优化后的 TS 结构
    """

    def __init__(self, config: dict, molecule_name: Optional[str] = None):
        """
        初始化反应分析引擎 - v3.0

        Args:
            config: 配置字典，包含:
                - theory.optimization: 几何优化配置
                - theory.single_point: L2 SP 配置
                - step3: 步骤特定配置 (IRC, 片段处理等)
            molecule_name: 分子名称（用于定位 v3.0 目录结构）
        """
        self.molecule_name = molecule_name  # v3.0
        self.config = config
        self.step3_config = config.get('step3', {})
        self.theory_opt = config.get('theory', {}).get('optimization', {})
        self.theory_sp = config.get('theory', {}).get('single_point', {})

        # 几何优化配置
        self.method = self.theory_opt.get('method', 'B3LYP')
        self.basis = self.theory_opt.get('basis', 'def2-SVP')
        self.dispersion = self.theory_opt.get('dispersion', 'GD3BJ')

        # IRC 配置
        self.verify_irc = self.step3_config.get('verify_irc', False)
        self.irc_max_points = self.step3_config.get('irc_max_points', 50)
        self.irc_step_size = self.step3_config.get('irc_step_size', 10)

        # 资源配置
        self.nprocshared = self.theory_opt.get('nproc', 16)
        self.mem = self.theory_opt.get('mem', '32GB')

        # 片段处理配置
        frag_config = self.step3_config.get('fragments', {})
        self.calculate_relaxed_fragments = frag_config.get('calculate_relaxed', True)
        self.fragment_opt_level = frag_config.get('opt_level', 'xtb')

        # 初始化几何优化驱动器
        self.berny_driver = BernyTSDriver(
            method=self.method,
            basis=self.basis,
            dispersion=self.dispersion,
            nprocshared=self.nprocshared,
            mem=self.mem,
            config=self.config
        )

        self.qst2_driver = QST2RescueDriver(
            method=self.method,
            basis=self.basis,
            dispersion=self.dispersion,
            nprocshared=self.nprocshared,
            mem=self.mem,
            config=self.config
        )

        self.validator = TSValidator()

        # IRC 驱动器（仅在需要时初始化）
        self.irc_driver: Optional[IRCDriver] = None

        # 初始化 ORCA SP 接口
        self.orca = ORCAInterface(
            method=self.theory_sp.get('method', 'WB97M-V'),
            basis=self.theory_sp.get('basis', 'def2-TZVPP'),
            aux_basis=self.theory_sp.get('aux_basis', 'def2/J'),
            nprocs=self.theory_sp.get('nproc', 16),
            maxcore=self.theory_sp.get('maxcore', 4000),
            solvent=self.theory_sp.get('solvent', 'acetone'),
            config=self.config
        )

        # 初始化 XTB 接口
        self.xtb = XTBInterface(
            gfn_level=2,
            solvent=self.theory_sp.get('solvent', 'acetone'),
            nproc=self.theory_sp.get('nproc', 8),
            config=self.config
        )

        # 初始化 QCTaskRunner（统一计算中枢）
        from rph_core.utils.qc_task_runner import QCTaskRunner
        self.qc_runner = QCTaskRunner(config=config)

        # 初始化片段提取器
        self.frag_extractor = FragmentExtractor(
            config={
                'method': self.method,
                'basis': self.basis,
                'dispersion': self.dispersion,
                'nprocshared': self.nprocshared,
                'mem': self.mem,
                'solvent': self.theory_sp.get('solvent', 'acetone')
            }
        )

        self.logger.info(f"TransitionAnalyzer 初始化: {self.method}/{self.basis} D3={self.dispersion}")
        self.logger.info(f"  L2 SP: {self.orca.method}/{self.orca.basis}")
        self.logger.info(f"  片段处理: {'启用' if self.calculate_relaxed_fragments else '禁用'} ({self.fragment_opt_level})")
        self.logger.info("  QCTaskRunner 已初始化（统一计算中枢）")
    
    def run_with_qctaskrunner(
        self,
        ts_guess: Path,
        reactant: Path,
        product: Path,
        output_dir: Path,
        e_product_l2: float,
        forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]] = None,
        old_checkpoint: Optional[Path] = None
    ) -> TransitionAnalysisResult:
        """
        使用 QCTaskRunner 执行 TS 优化（统一计算中枢）- v3.0
        """
        def safe_extract_coords(file_path: Path, name: str):
            """安全提取坐标的辅助函数"""
            coords, symbols, error = LogParser.extract_last_converged_coords(
                file_path,
                engine_type='auto'
            )
            if coords is None:
                self.logger.warning(f"无法从 {name} ({file_path}) 提取坐标: {error}")
                # 回退
                from rph_core.utils.file_io import read_xyz
                return read_xyz(file_path)
            self.logger.info(f"从 {name} 成功提取 {len(coords)} 个原子坐标")
            return coords, symbols

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("=" * 60)
        self.logger.info("Step 3: TS 优化（使用 QCTaskRunner 统一计算中枢）- v3.0")
        self.logger.info("=" * 60)

        # 提取输入坐标（使用 LogParser）
        self.logger.info("验证输入文件...")
        ts_coords, ts_symbols = safe_extract_coords(ts_guess, "TS 初猜")
        reactant_coords, reactant_symbols = safe_extract_coords(reactant, "底物复合物")
        product_coords, product_symbols = safe_extract_coords(product, "产物")

        # 1. TS 优化 + L2 SP
        self.logger.info("步骤 1/3: TS 优化 + 虚频验证 + L2 SP")
        ts_result = self.qc_runner.run_ts_opt_cycle(
            xyz_file=ts_guess,
            output_dir=output_dir / "ts_opt",
            charge=0,
            spin=1,
            enable_l2_sp=True,
            old_checkpoint=old_checkpoint
        )

        if not ts_result.converged:
            raise RuntimeError(f"TS 优化失败: {ts_result.error_message}")

        # 2. Reactant L2 SP 计算
        self.logger.info("步骤 2/3: Reactant (Complex) L2 SP 计算...")
        e_reactant_l2 = self._run_sp(reactant, output_dir / "reactant_sp")
        
        # 3. 片段处理 (基于 TS 优化后的几何)
        self.logger.info("步骤 3/3: 构建 SP 矩阵报告 (含片段)...")
        sp_report = self._build_sp_matrix(
            ts_final_xyz=ts_result.optimized_xyz,
            reactant=reactant,
            e_product_l2=e_product_l2,
            output_dir=output_dir / "SP_Matrix",
            forming_bonds=forming_bonds
        )
        # 更新报告中的 Reactant 能量（默认是复合物能量）
        sp_report.e_reactant = e_reactant_l2
        sp_report.e_ts_final = ts_result.l2_energy

        # 返回结果
        result = TransitionAnalysisResult(
            ts_final_xyz=ts_result.optimized_xyz,
            ts_checkpoint=ts_result.optimized_xyz.with_suffix('.chk') if ts_result.optimized_xyz.with_suffix('.chk').exists() else None,
            sp_report=sp_report,
            method_used=ts_result.method_used
        )
        
        self.logger.info("=" * 60)
        self.logger.info(f"ΔG‡ (Complex-based): {sp_report.get_activation_energy():.2f} kcal/mol")
        self.logger.info("Step 3 完成")
        self.logger.info("=" * 60)
        
        return result

    def _build_sp_matrix(
        self,
        ts_final_xyz: Path,
        reactant: Path,
        e_product_l2: float,
        output_dir: Path,
        forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]]
    ) -> SPMatrixReport:
        """
        构建完整 SP 矩阵 (TS SP + Reactant SP + Fragments SP)
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 1. TS L2 SP
        self.logger.info("  [TS_Final] 执行 L2 SP...")
        e_ts_l2 = self._run_sp(ts_final_xyz, output_dir / "ts_sp")

        # 2. Reactant L2 SP
        self.logger.info("  [Reactant] 执行 L2 SP...")
        e_reactant_l2 = self._run_sp(reactant, output_dir / "reactant_sp")

        # 3. 提取片段并计算能量
        e_frag_a_l2 = 0.0
        e_frag_b_l2 = 0.0
        e_frag_a_relaxed = None
        e_frag_b_relaxed = None

        if forming_bonds:
            self.logger.info("  提取片段并计算能量...")
            # 简化片段分配：假设 forming_bonds[0] 将体系分为两部分
            (bond1, bond2) = forming_bonds
            split_idx = bond1[1] # 简单切分点
            
            coords, symbols = read_xyz(ts_final_xyz)
            frag_a_indices = list(range(0, split_idx))
            frag_b_indices = list(range(split_idx, len(symbols)))
            
            frag_results = self.frag_extractor.extract_and_calculate(
                ts_xyz=ts_final_xyz,
                fragment_indices=(frag_a_indices, frag_b_indices),
                output_dir=output_dir / "fragments"
            )
            
            e_frag_a_l2 = frag_results['e_fragment_a_ts']
            e_frag_b_l2 = frag_results['e_fragment_b_ts']
            e_frag_a_relaxed = frag_results['e_fragment_a_relaxed']
            e_frag_b_relaxed = frag_results['e_fragment_b_relaxed']

        # 构建报告
        report = SPMatrixReport(
            e_product=e_product_l2,
            e_reactant=e_reactant_l2,
            e_ts_final=e_ts_l2,
            e_frag_a_ts=e_frag_a_l2,
            e_frag_b_ts=e_frag_b_l2,
            e_frag_a_relaxed=e_frag_a_relaxed,
            e_frag_b_relaxed=e_frag_b_relaxed,
            method=f"{self.orca.method}/{self.orca.basis}",
            solvent=self.orca.solvent
        )

        return report

    def _run_sp(self, xyz_file: Path, output_dir: Path) -> float:
        """运行 ORCA SP 计算"""
        result = self.orca.single_point(xyz_file, output_dir)
        if not result.converged:
             raise RuntimeError(f"SP 计算失败: {result.error_message}")
        return result.energy

    def run(self, *args, **kwargs):
        # 保持旧接口兼容性（内部调用新逻辑）
        return self.run_with_qctaskrunner(*args, **kwargs)
