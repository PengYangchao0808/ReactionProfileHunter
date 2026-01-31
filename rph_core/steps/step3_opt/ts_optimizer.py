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
from typing import Optional, Tuple, List, cast
from dataclasses import dataclass
import numpy as np
import hashlib
import json
import datetime

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import write_xyz, read_xyz
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.utils.qc_interface import XTBInterface
from rph_core.utils.geometry_tools import LogParser
# ASM 功能已禁用 (2026-01-31) - 以下导入不再使用
# from rph_core.steps.step3_opt.post_qc_enrichment import run_post_qc_enrichment
# v3.0: SPMatrixReport 已移除（批处理 SP 功能被 ConformerEngine 替代）
# 创建简化的占位符类以保持接口兼容性
from dataclasses import dataclass
from typing import Optional
@dataclass
class SPMatrixReport:
    e_ts: float = 0.0
    e_reactant: Optional[float] = None
    e_product: Optional[float] = None
    e_reactant_l2: Optional[float] = None
    e_product_l2: Optional[float] = None
    e_ts_final: Optional[float] = None
    g_ts: Optional[float] = None
    g_reactant: Optional[float] = None
    g_product: Optional[float] = None
    # 片段能量（v3.0 可能不使用，占位符）
    e_frag_a_ts: float = 0.0
    e_frag_b_ts: float = 0.0
    e_frag_a_relaxed: float = 0.0
    e_frag_b_relaxed: float = 0.0
    fragment_split_source: Optional[str] = None
    fragment_split_reason: Optional[str] = None
    fragment_indices: Optional[Tuple[List[int], List[int]]] = None
    # 方法和溶剂（占位符）
    method: str = "Berny"
    solvent: str = "acetone"

    def get_activation_energy(self) -> Optional[float]:
        g_ts = self.g_ts
        g_reactant = self.g_reactant
        if g_ts is not None and g_reactant is not None:
            return float(g_ts) - float(g_reactant)
        e_ts_final = self.e_ts_final
        e_reactant = self.e_reactant
        if e_ts_final is not None and e_reactant is not None:
            return (e_ts_final - e_reactant) * 627.509
        e_ts = self.e_ts
        if e_reactant is not None:
            return (e_ts - e_reactant) * 627.509
        return None

    def get_reaction_energy(self) -> Optional[float]:
        g_product = self.g_product
        g_reactant = self.g_reactant
        if g_product is not None and g_reactant is not None:
            return float(g_product) - float(g_reactant)
        e_product = self.e_product
        e_reactant = self.e_reactant
        if e_product is not None and e_reactant is not None:
            return (e_product - e_reactant) * 627.509
        return None
    
    @staticmethod
    def _write_artifacts_index(s3_dir: Path, config: dict) -> None:
        """
        Write S3 artifacts_index.json with dipolar artifact info.

        Args:
            s3_dir: Step3 output directory
            config: Configuration dict
        """
        import json
        import logging
        from pathlib import Path
        
        enrichment_write_dir = config.get('step3', {}).get('enrichment', {}).get('write_dirname', 'S3_PostQCEnrichment')
        index_path = s3_dir / 'artifacts_index.json'
        
        dipolar_info = {}
        dipolar_dir = s3_dir / 'dipolar'
        
        if dipolar_dir.exists():
            dipolar_patterns = ['*dipolar*.log', '*dipolar*.out']
            for pattern in dipolar_patterns:
                for dipolar_output in sorted(dipolar_dir.glob(pattern)):
                    dipolar_info = {
                        'path_rel': f'dipolar/{dipolar_output.name}',
                        'sha256': hashlib.sha256(dipolar_output.read_bytes()).hexdigest()
                    }
                    break
                if dipolar_info:
                    break
        
        index_data = {
            'schema_version': 'artifacts_index_v1',
            'created_at': __import__('datetime').datetime.utcnow().isoformat() + 'Z',
            'dipolar': dipolar_info if dipolar_info else None
        }
        
        with open(index_path, 'w', encoding='utf-8') as f:
            json.dump(index_data, f, indent=2)
        
        logging.getLogger(__name__).info(f"Wrote artifacts index to {index_path}")

# ASM 功能已禁用 (2026-01-31) - 以下导入不再使用
# from rph_core.steps.step4_features.fragment_extractor import FragmentExtractor
from typing import Any
import json
import os
import shutil
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
    ts_fchk: Optional[Path] = None
    ts_log: Optional[Path] = None
    ts_qm_output: Optional[Path] = None
    reactant_fchk: Optional[Path] = None
    reactant_log: Optional[Path] = None
    reactant_qm_output: Optional[Path] = None


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

        # ASM 功能已禁用 (2026-01-31) - FragmentExtractor 不再使用
        # self.frag_extractor = FragmentExtractor(
        #     config={'solvent': self.theory_sp.get('solvent', 'acetone')},
        #     sp_engine=self.orca
        # )

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
        product_thermo: Optional[Path] = None,
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

        # 2. Reactant OPT + Freq + L2 SP (完整工作流)
        self.logger.info("步骤 2/3: Reactant 优化 + 频率 + L2 SP...")
        
        # Get charge/multiplicity from config or defaults
        reactant_opt_config = self.step3_config.get('reactant_opt', {})
        reactant_charge = reactant_opt_config.get('charge', 0)
        reactant_mult = reactant_opt_config.get('multiplicity', 1)
        enable_nbo = reactant_opt_config.get('enable_nbo', False)
        
        # Standard optimization attempt
        self.logger.info(f"尝试标准 Reactant 优化 (charge={reactant_charge}, mult={reactant_mult})...")
        reactant_opt_result = self.qc_runner.run_opt_sp_cycle(
            xyz_file=reactant,
            output_dir=output_dir / "reactant_opt" / "standard",
            charge=reactant_charge,
            spin=reactant_mult,
            enable_nbo=enable_nbo,
            old_checkpoint=None
        )
        
        # Rescue if standard fails
        if not reactant_opt_result.converged:
            self.logger.warning(f"Reactant 标准优化失败: {reactant_opt_result.error_message}")
            self.logger.info("尝试 Reactant 救援优化（由 QCTaskRunner 内部处理）...")
            reactant_opt_result = self.qc_runner.run_opt_sp_cycle(
                xyz_file=reactant,
                output_dir=output_dir / "reactant_opt" / "rescue",
                charge=reactant_charge,
                spin=reactant_mult,
                enable_nbo=enable_nbo,
                old_checkpoint=None
            )
            if not reactant_opt_result.converged:
                raise RuntimeError(f"Reactant 优化失败（含rescue）: {reactant_opt_result.error_message}")
        
        assert reactant_opt_result.l2_energy is not None, "Reactant L2 能量缺失"
        e_reactant_l2 = reactant_opt_result.l2_energy
        reactant_optimized_xyz = reactant_opt_result.optimized_xyz
        self.logger.info(f"✓ Reactant 优化完成，L2 能量: {e_reactant_l2:.6f} Hartree")
        
        if e_product_l2 is None:
            self.logger.warning("Product L2 能量缺失，执行补算 SP。")
            product_sp_result = self._run_sp(product, output_dir / "product_sp")
            assert product_sp_result.energy is not None
            e_product_l2 = cast(float, product_sp_result.energy)

        # 3. 片段处理 (基于 TS 优化后的几何)
        self.logger.info("步骤 3/3: 构建 SP 矩阵报告 (含片段)...")
        ts_l2_result = ts_result.l2_sp_result
        if ts_l2_result is None or not ts_l2_result.converged:
            raise RuntimeError("TS L2 SP 结果缺失或未收敛，无法构建 SP 矩阵报告。")

        sp_report = self._build_sp_matrix(
            ts_final_xyz=ts_result.optimized_xyz,
            reactant=reactant,
            e_product_l2=e_product_l2,
            e_reactant_l2=e_reactant_l2,
            output_dir=output_dir / "ASM_SP_Mat",
            forming_bonds=forming_bonds,
            ts_l2_result=ts_l2_result,
            ts_l2_dir=output_dir / "ts_opt" / "L2_SP"
        )
        # 更新报告中的 Reactant 能量（默认是复合物能量）
        sp_report.e_reactant = e_reactant_l2
        
        # ASM 功能已禁用 (2026-01-31) - Post-QC enrichment 不再运行
        # if forming_bonds is not None:
        #     self.logger.info("Running post-QC enrichment...")
        #     enrichment_config = self.config.get('step3', {}).get('enrichment', {})
        #     run_post_qc_enrichment(
        #         s3_dir=output_dir,
        #         config=enrichment_config,
        #         sp_report=sp_report,
        #         reactant_complex_xyz=reactant,
        #         forming_bonds=forming_bonds
        #     )

        dipolar_output = self._ensure_dipolar_output(output_dir, output_dir / "ts_opt" / "L2_SP")
        if dipolar_output is None:
            self.logger.warning("Dipolar output not found; artifacts_index will omit dipolar entry")

        # V6.3: Write artifacts_index.json
        SPMatrixReport._write_artifacts_index(output_dir, self.config)
        
        # 返回结果
        result = TransitionAnalysisResult(
            ts_final_xyz=ts_result.optimized_xyz,
            ts_checkpoint=ts_result.optimized_xyz.with_suffix('.chk') if ts_result.optimized_xyz.with_suffix('.chk').exists() else None,
            sp_report=sp_report,
            method_used=ts_result.method_used,
            ts_fchk=ts_result.fchk_file,
            ts_log=ts_result.log_file,
            ts_qm_output=ts_result.qm_output_file,
            reactant_fchk=reactant_opt_result.fchk_file,
            reactant_log=reactant_opt_result.freq_log or reactant_opt_result.log_file,
            reactant_qm_output=reactant_opt_result.qm_output_file
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
        e_reactant_l2: float,
        output_dir: Path,
        forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]],
        ts_l2_result: Any,
        ts_l2_dir: Path
    ) -> SPMatrixReport:
        """
        构建完整 SP 矩阵 (TS SP + Reactant SP + Fragments SP)
        """
        output_dir = Path(output_dir)
        # Don't create output_dir yet - only create when we actually write files

        # 1. TS L2 SP (authoritative): reuse ts_opt/L2_SP artifacts
        self.logger.info("  [TS_Final] 复用 ts_opt/L2_SP 作为权威 L2 SP...")
        e_ts_l2 = ts_l2_result.energy
        self._alias_ts_sp_dir(output_dir / "ts_sp", ts_l2_dir)

        # ASM 功能已禁用 (2026-01-31)
        # 原因: 片段分割算法未能正常工作
        e_frag_a_l2 = 0.0
        e_frag_b_l2 = 0.0
        e_frag_a_relaxed = 0.0
        e_frag_b_relaxed = 0.0
        fragment_split_source = None
        fragment_split_reason = None
        fragment_indices = None

        # if forming_bonds:
        #     self.logger.info("  提取片段并计算能量...")
        #     split_config = self.step3_config.get("fragment_split", None)
        #     if split_config is None:
        #         split_config = self.step3_config.get("fragments", {})
        #
        #     frag_results = self.frag_extractor.extract_and_calculate(
        #         ts_xyz=ts_final_xyz,
        #         fragment_indices=None,
        #         output_dir=output_dir,
        #         forming_bonds=forming_bonds,
        #         reactant_xyz=reactant,
        #         split_config=split_config,
        #         system_charge=0,
        #         system_mult=1
        #     )
        #
        #     e_frag_a_l2 = frag_results.get("e_fragment_a_ts", 0.0)
        #     e_frag_b_l2 = frag_results.get("e_fragment_b_ts", 0.0)
        #     e_frag_a_relaxed = frag_results.get("e_fragment_a_relaxed", 0.0)
        #     e_frag_b_relaxed = frag_results.get("e_fragment_b_relaxed", 0.0)
        #     fragment_split_source = frag_results.get("fragment_split_source")
        #     fragment_split_reason = frag_results.get("fragment_split_reason")
        #     fragment_indices = frag_results.get("fragment_indices")

        # 构建报告
        report = SPMatrixReport(
            e_product=e_product_l2,
            e_reactant=e_reactant_l2,
            e_ts_final=e_ts_l2,
            e_frag_a_ts=e_frag_a_l2,
            e_frag_b_ts=e_frag_b_l2,
            e_frag_a_relaxed=e_frag_a_relaxed,
            e_frag_b_relaxed=e_frag_b_relaxed,
            fragment_split_source=fragment_split_source,
            fragment_split_reason=fragment_split_reason,
            fragment_indices=fragment_indices,
            method=f"{self.orca.method}/{self.orca.basis}",
            solvent=self.orca.solvent
        )

        return report

    def _alias_ts_sp_dir(self, ts_sp_dir: Path, l2_sp_dir: Path) -> None:
        """Make ASM_SP_Mat/ts_sp a view of ts_opt/L2_SP without recomputation."""
        ts_sp_dir = Path(ts_sp_dir)
        l2_sp_dir = Path(l2_sp_dir)
        ts_sp_dir.parent.mkdir(parents=True, exist_ok=True)

        if ts_sp_dir.exists():
            # If already present, do not delete; assume resume.
            return

        try:
            os.symlink(l2_sp_dir, ts_sp_dir, target_is_directory=True)
            (ts_sp_dir / "provenance.json").write_text(
                '{"mode": "symlink", "source_dir": ' + json.dumps(str(l2_sp_dir)) + '}'
            )
            return
        except Exception:
            pass

        try:
            shutil.copytree(l2_sp_dir, ts_sp_dir)
            (ts_sp_dir / "provenance.json").write_text(
                '{"mode": "copy", "copied_from": ' + json.dumps(str(l2_sp_dir)) + '}'
            )
        except Exception as e:
            self.logger.warning(f"Failed to alias ts_sp dir: {e}")

    def _run_sp(self, xyz_file: Path, output_dir: Path):
        """运行 ORCA SP 计算"""
        result = self.orca.single_point(xyz_file, output_dir)
        if not result.converged:
             raise RuntimeError(f"SP 计算失败: {result.error_message}")
        if result.energy is None:
            raise RuntimeError("SP 计算失败: 未返回能量值")
        return result

    def _ensure_dipolar_output(self, s3_dir: Path, ts_l2_dir: Path) -> Optional[Path]:
        """Ensure a dipolar output file exists under S3/dipolar."""
        s3_dir = Path(s3_dir)
        dipolar_dir = s3_dir / "dipolar"
        
        # Check if dipolar directory already exists with files
        if dipolar_dir.exists():
            existing = list(dipolar_dir.glob("*dipolar*.out")) + list(dipolar_dir.glob("*dipolar*.log"))
            if existing:
                return existing[0]

        candidate_dirs = [
            ts_l2_dir,
            s3_dir / "ASM_SP_Mat" / "ts_sp",
            s3_dir / "ts_opt" / "L2_SP"
        ]

        for candidate_dir in candidate_dirs:
            candidate_dir = Path(candidate_dir)
            if not candidate_dir.exists():
                continue
            candidates = list(candidate_dir.glob("*.out")) + list(candidate_dir.glob("*.log"))
            if not candidates:
                continue
            source = candidates[0]
            # Only create directory when we actually have a file to copy
            dipolar_dir.mkdir(parents=True, exist_ok=True)
            target = dipolar_dir / f"ts_dipolar{source.suffix}"
            try:
                shutil.copy2(source, target)
                return target
            except Exception as e:
                self.logger.warning(f"Failed to copy dipolar output from {source}: {e}")
                return None

        return None

    def run(self, *args, **kwargs):
        # 保持旧接口兼容性（内部调用新逻辑）
        return self.run_with_qctaskrunner(*args, **kwargs)
