"""
Step 3: Transition Analyzer
==========================

Reaction center analysis module - TS optimization validation + Intermediate SP + fragment processing

Author: QCcalc Team
Date: 2026-01-09
Updated: 2026-01-13 (merged S3.5, integrated SP and fragment processing)
"""

import logging
from pathlib import Path
from typing import Optional, Tuple, List, cast, Dict, Any
from dataclasses import dataclass, asdict
import numpy as np
import hashlib
import json
import datetime

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import write_xyz, read_xyz
from rph_core.utils.naming_compat import (
    INTERMEDIATE_XYZ,
    REACTANT_COMPLEX_XYZ,
    resolve_intermediate_path,
    DIR_INTERMEDIATE_OPT,
)
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.utils.qc_interface import XTBInterface
from rph_core.utils.geometry_tools import LogParser
from rph_core.utils.ui import get_progress_manager
from rph_core.utils.constants import HARTREE_TO_KCAL


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
    # Shermo 来源与错误追踪
    g_ts_source: Optional[str] = None
    g_ts_error: Optional[str] = None
    g_reactant_source: Optional[str] = None
    g_reactant_error: Optional[str] = None
    # 片段能量（v3.0 可能不使用，占位符）
    e_frag_a_ts: float = 0.0
    e_frag_b_ts: float = 0.0
    e_frag_a_relaxed: Optional[float] = None
    e_frag_b_relaxed: Optional[float] = None
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
            return (e_ts_final - e_reactant) * HARTREE_TO_KCAL
        e_ts = self.e_ts
        if e_reactant is not None:
            return (e_ts - e_reactant) * HARTREE_TO_KCAL
        return None

    def get_reaction_energy(self) -> Optional[float]:
        g_product = self.g_product
        g_reactant = self.g_reactant
        if g_product is not None and g_reactant is not None:
            return float(g_product) - float(g_reactant)
        e_product = self.e_product
        e_reactant = self.e_reactant
        if e_product is not None and e_reactant is not None:
            return (e_product - e_reactant) * HARTREE_TO_KCAL
        return None

    def to_dict(self) -> dict[str, Any]:
        """Serialize to dictionary using dataclass field names."""
        import dataclasses
        return dataclasses.asdict(self)

    def to_json(self, indent: int = 2) -> str:
        """Serialize to JSON string."""
        import json
        data = self.to_dict()
        return json.dumps(data, indent=indent, ensure_ascii=False)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "SPMatrixReport":
        """Create instance from dictionary, filtering unknown keys."""
        import dataclasses
        valid_fields = {f.name for f in dataclasses.fields(cls)}
        filtered = {k: v for k, v in data.items() if k in valid_fields}
        return cls(**filtered)

    def validate(self) -> bool:
        """Validate data integrity. Check numeric/string fields."""
        import dataclasses
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if value is None:
                continue
            type_str = str(field.type)
            if 'float' in type_str:
                if not isinstance(value, (int, float)):
                    return False
            elif 'str' in type_str:
                if not isinstance(value, str):
                    return False
        return True

    def __str__(self) -> str:
        """Formatted string representation with key energy info."""
        lines = [f"SPMatrixReport(method={self.method}, solvent={self.solvent})"]
        lines.append(f"  e_reactant  = {self.e_reactant}")
        lines.append(f"  e_product   = {self.e_product}")
        lines.append(f"  e_ts_final  = {self.e_ts_final}")

        activation = self.get_activation_energy()
        if activation is not None:
            lines.append(f"  ΔG‡ (activation) = {activation:.4f} kcal/mol")

        reaction = self.get_reaction_energy()
        if reaction is not None:
            lines.append(f"  ΔG_rxn (reaction) = {reaction:.4f} kcal/mol")

        return "\n".join(lines)
    
    @staticmethod
    def _write_artifacts_index(s3_dir: Path, config: dict[str, Any]) -> None:
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
    intermediate_fchk: Optional[Path] = None
    intermediate_log: Optional[Path] = None
    intermediate_qm_output: Optional[Path] = None


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
    - intermediate: 反应中间体 (来自 Step 2，用于 QST2 救援)
    - product: 产物 (来自 Step 1，用于 QST2 救援)

    输出:
    - ts_final_xyz: 优化后的 TS 结构
    """

    def __init__(self, config: dict[str, Any], molecule_name: Optional[str] = None):
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

        # 中间体 DFT 优化驱动
        int_cfg = self.step3_config.get('intermediate_opt', {})
        if int_cfg.get('enabled', True):
            from rph_core.steps.step3_opt.intermediate_driver import IntermediateDriver
            self.intermediate_driver = IntermediateDriver(config)
            self.logger.info(f"  中间体 DFT: 启用")
        else:
            self.intermediate_driver = None
            self.logger.info(f"  中间体 DFT: 禁用")

        self.logger.info(f"TransitionAnalyzer 初始化: {self.method}/{self.basis} D3={self.dispersion}")
        self.logger.info(f"  L2 SP: {self.orca.method}/{self.orca.basis}")
        self.logger.info(f"  片段处理: {'启用' if self.calculate_relaxed_fragments else '禁用'} ({self.fragment_opt_level})")
        self.logger.info("  QCTaskRunner 已初始化（统一计算中枢）")

    def _load_s3_resume_state(self, output_dir: Path) -> Dict[str, Any]:
        resume_file = output_dir / "s3_resume.json"
        if not resume_file.exists():
            return {}
        try:
            with open(resume_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            self.logger.warning(f"Failed to load s3_resume.json: {e}")
            return {}

    def _save_s3_resume_state(self, output_dir: Path, state: Dict[str, Any]) -> None:
        resume_file = output_dir / "s3_resume.json"
        try:
            with open(resume_file, 'w') as f:
                json.dump(state, f, indent=2)
        except Exception as e:
            self.logger.warning(f"Failed to save s3_resume.json: {e}")

    def _verify_ts_result(self, ts_result, output_dir: Path) -> bool:
        if not ts_result or not ts_result.converged:
            return False
        if not ts_result.optimized_xyz or not ts_result.optimized_xyz.exists():
            return False
        l2_sp_candidates = [
            output_dir / "L2_SP",
            output_dir / "ts_opt" / "L2_SP",
        ]
        l2_sp_dir = None
        for candidate in l2_sp_candidates:
            if candidate.exists():
                out_files = list(candidate.glob("*.out")) + list(candidate.glob("*.log"))
                if out_files:
                    l2_sp_dir = candidate
                    break
        if l2_sp_dir is None:
            return False
        return True

    def _verify_reactant_result(self, reactant_result, output_dir: Path) -> bool:
        if not reactant_result or not reactant_result.converged:
            return False
        if not reactant_result.optimized_xyz or not reactant_result.optimized_xyz.exists():
            return False
        if reactant_result.imaginary_count != 0:
            return False
        l2_sp_candidates = [
            output_dir / "L2_SP",
            output_dir / "standard" / "L2_SP",
        ]
        l2_sp_dir = None
        for candidate in l2_sp_candidates:
            if candidate.exists():
                out_files = list(candidate.glob("*.out")) + list(candidate.glob("*.log"))
                if out_files:
                    l2_sp_dir = candidate
                    break
        if l2_sp_dir is None:
            return False
        return True

    def run_with_qctaskrunner(
        self,
        ts_guess: Path,
        intermediate: Path,
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
        pm = get_progress_manager()
        ui_step_id = "s3"

        def ui_log(message: str) -> None:
            pm.log_event("S3", message)

        pm.update_step("s3", description="S3: Starting...")
        pm.enter_phase(ui_step_id, "Initialization")
        pm.set_subtask(ui_step_id, "S3 Workflow", 0, 3)
        ui_log("Transition analysis started")

        def safe_extract_coords(file_path: Path, name: str):
            """安全提取坐标的辅助函数"""
            coords, symbols, error = LogParser.extract_last_converged_coords(
                file_path,
                engine_type='auto'
            )
            if coords is None:
                warning_msg = f"无法从 {name} ({file_path}) 提取坐标: {error}"
                self.logger.warning(warning_msg)
                ui_log(warning_msg)
                # 回退
                from rph_core.utils.file_io import read_xyz
                return read_xyz(file_path)
            info_msg = f"从 {name} 成功提取 {len(coords)} 个原子坐标"
            self.logger.info(info_msg)
            ui_log(info_msg)
            return coords, symbols

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        resume_state = self._load_s3_resume_state(output_dir)

        with open(output_dir / ".rph_step_status.json", "w") as f:
            json.dump({
                "step": "s3",
                "status": "in_progress",
                "start_time": datetime.datetime.now().isoformat(),
                "description": "S3: Transition Analysis Started"
            }, f, indent=2)

        self.logger.info("=" * 60)
        self.logger.info("Step 3: TS 优化（使用 QCTaskRunner 统一计算中枢）- v3.0")
        self.logger.info("=" * 60)
        ui_log("Step 3 execution initialized")

        self.logger.info("验证输入文件...")
        pm.enter_phase(ui_step_id, "Input Validation")
        ui_log("Validating input structures")
        ts_coords, ts_symbols = safe_extract_coords(ts_guess, "TS 初猜")
        intermediate_coords, intermediate_symbols = safe_extract_coords(intermediate, "中间体")
        product_coords, product_symbols = safe_extract_coords(product, "产物")

        # ========================================================================
        # Phase 0: 中间体 DFT 优化 (可选)
        # ========================================================================
        intermediate_result = None
        intermediate_optimized_xyz = None
        driver = getattr(self, 'intermediate_driver', None)
        if driver is not None:
            self.logger.info("步骤 0/3: 中间体 DFT 优化")
            pm.enter_phase(ui_step_id, "Intermediate Optimization")
            pm.set_subtask(ui_step_id, "Intermediate Opt", 0, 3)
            pm.update_step("s3", description="S3: Intermediate Opt")
            ui_log("Starting intermediate DFT optimization")

            intermediate_result = driver.run(
                intermediate_xyz=intermediate,
                output_dir=output_dir / DIR_INTERMEDIATE_OPT,
            )

            if intermediate_result.get('converged'):
                intermediate_optimized_xyz = intermediate_result.get('optimized_xyz')
                self.logger.info(f"中间体优化完成: {intermediate_optimized_xyz}")
                self.logger.info(f"  -> 将作为 Reactant 优化的输入，避免重复优化")
                ui_log("Intermediate optimization converged")
            elif intermediate_result.get('skipped'):
                ui_log("Intermediate optimization skipped (disabled)")
            else:
                self.logger.warning(f"中间体优化失败: {intermediate_result.get('error')}")
                ui_log(f"Intermediate optimization failed: {intermediate_result.get('error')}")

        ts_result = None
        ts_output_dir = output_dir / "ts_opt"

        if resume_state.get("ts_opt", {}).get("completed"):
            self.logger.info("检测到已完成的 TS 优化结果，尝试复用...")
            ui_log("Found resumable TS optimization result")
            ts_result_data = resume_state["ts_opt"]
            from rph_core.utils.qc_task_runner import QCOptimizationResult
            ts_result = QCOptimizationResult(
                optimized_xyz=Path(ts_result_data["optimized_xyz"]),
                l2_energy=ts_result_data.get("l2_energy"),
                opt_energy=ts_result_data.get("opt_energy"),
                converged=ts_result_data.get("converged", False),
                imaginary_count=ts_result_data.get("imaginary_count", 0),
                method_used=ts_result_data.get("method_used", "Berny"),
                freq_log=Path(ts_result_data["freq_log"]) if ts_result_data.get("freq_log") else None,
                log_file=Path(ts_result_data["log_file"]) if ts_result_data.get("log_file") else None,
                chk_file=Path(ts_result_data["chk_file"]) if ts_result_data.get("chk_file") else None,
                fchk_file=Path(ts_result_data["fchk_file"]) if ts_result_data.get("fchk_file") else None,
                qm_output_file=Path(ts_result_data["qm_output_file"]) if ts_result_data.get("qm_output_file") else None,
            )
            if not self._verify_ts_result(ts_result, ts_output_dir):
                self.logger.warning("TS 优化结果验证失败，将重新计算")
                ui_log("Resumed TS result validation failed, recomputing")
                ts_result = None
            else:
                self.logger.info("TS 优化结果验证通过，成功复用")
                ui_log("Reused TS optimization result")

        if ts_result is None:
            self.logger.info("步骤 1/3: TS 优化 + 虚频验证 + L2 SP")
            pm.enter_phase(ui_step_id, "TS Optimization")
            pm.set_subtask(ui_step_id, "TS Opt", 1, 3)
            pm.update_step("s3", description="S3: TS Opt")
            ui_log("Starting TS optimization + frequency validation + L2 SP")
            ts_result = self.qc_runner.run_ts_opt_cycle(
                xyz_file=ts_guess,
                output_dir=output_dir,
                charge=0,
                spin=1,
                enable_l2_sp=True,
                old_checkpoint=None
            )

            if ts_result.converged and self._verify_ts_result(ts_result, ts_output_dir):
                ui_log("TS optimization converged and validated")
                resume_state["ts_opt"] = {
                    "completed": True,
                    "optimized_xyz": str(ts_result.optimized_xyz),
                    "l2_energy": ts_result.l2_energy,
                    "opt_energy": ts_result.opt_energy,
                    "imaginary_count": ts_result.imaginary_count,
                    "method_used": ts_result.method_used,
                    "freq_log": str(ts_result.freq_log) if ts_result.freq_log else None,
                    "log_file": str(ts_result.log_file) if ts_result.log_file else None,
                    "chk_file": str(ts_result.chk_file) if ts_result.chk_file else None,
                    "fchk_file": str(ts_result.fchk_file) if ts_result.fchk_file else None,
                    "qm_output_file": str(ts_result.qm_output_file) if ts_result.qm_output_file else None,
                }
                self._save_s3_resume_state(output_dir, resume_state)

        if not ts_result.converged:
            ui_log(f"TS optimization failed: {ts_result.error_message}")
            with open(output_dir / ".rph_step_status.json", "w") as f:
                json.dump({
                    "step": "s3",
                    "status": "failed",
                    "error": ts_result.error_message,
                    "end_time": datetime.datetime.now().isoformat()
                }, f, indent=2)
            raise RuntimeError(f"TS 优化失败: {ts_result.error_message}")

        reactant_opt_result = None
        reactant_output_dir = output_dir / "reactant_opt"

        if resume_state.get("reactant_opt", {}).get("completed"):
            self.logger.info("检测到已完成的 Reactant 优化结果，尝试复用...")
            ui_log("Found resumable Reactant optimization result")
            reactant_result_data = resume_state["reactant_opt"]
            from rph_core.utils.qc_task_runner import QCOptimizationResult
            reactant_opt_result = QCOptimizationResult(
                optimized_xyz=Path(reactant_result_data["optimized_xyz"]),
                l2_energy=reactant_result_data.get("l2_energy"),
                opt_energy=reactant_result_data.get("opt_energy"),
                converged=reactant_result_data.get("converged", False),
                imaginary_count=reactant_result_data.get("imaginary_count", 0),
                method_used=reactant_result_data.get("method_used", "Normal"),
                freq_log=Path(reactant_result_data["freq_log"]) if reactant_result_data.get("freq_log") else None,
                log_file=Path(reactant_result_data["log_file"]) if reactant_result_data.get("log_file") else None,
                chk_file=Path(reactant_result_data["chk_file"]) if reactant_result_data.get("chk_file") else None,
                fchk_file=Path(reactant_result_data["fchk_file"]) if reactant_result_data.get("fchk_file") else None,
                qm_output_file=Path(reactant_result_data["qm_output_file"]) if reactant_result_data.get("qm_output_file") else None,
            )
            if not self._verify_reactant_result(reactant_opt_result, reactant_output_dir / "standard"):
                self.logger.warning("Reactant 优化结果验证失败，将重新计算")
                ui_log("Resumed Reactant result validation failed, recomputing")
                reactant_opt_result = None
            else:
                self.logger.info("Reactant 优化结果验证通过，成功复用")
                ui_log("Reused Reactant optimization result")

        if reactant_opt_result is None:
            # 复用中间体结果作为 Reactant，避免重复 DFT 优化
            if intermediate_result is not None and intermediate_result.get('converged') and intermediate_optimized_xyz:
                self.logger.info("=" * 60)
                self.logger.info("复用中间体优化结果作为 Reactant (跳过重复 DFT 优化)")
                self.logger.info(f"  来源: intermediate_optimized_xyz = {intermediate_optimized_xyz}")
                self.logger.info(f"  L2 Energy: {intermediate_result.get('l2_energy'):.6f} Hartree")
                self.logger.info("=" * 60)
                ui_log("Reusing intermediate optimization result as Reactant")
                
                from rph_core.utils.qc_task_runner import QCOptimizationResult
                reactant_opt_result = QCOptimizationResult(
                    optimized_xyz=intermediate_optimized_xyz,
                    l2_energy=intermediate_result.get('l2_energy'),
                    opt_energy=intermediate_result.get('l2_energy'),
                    converged=True,
                    imaginary_count=0,
                    method_used="Intermediate_Reuse",
                    freq_log=intermediate_result.get('freq_output'),
                    log_file=intermediate_result.get('opt_output'),
                    chk_file=None,
                    fchk_file=None,
                    qm_output_file=intermediate_result.get('sp_output'),
                )
                
                resume_state["reactant_opt"] = {
                    "completed": True,
                    "source": "intermediate_reuse",
                    "optimized_xyz": str(intermediate_optimized_xyz),
                    "l2_energy": intermediate_result.get('l2_energy'),
                    "opt_energy": intermediate_result.get('l2_energy'),
                    "imaginary_count": 0,
                    "method_used": "Intermediate_Reuse",
                    "freq_log": str(intermediate_result.get('freq_output')) if intermediate_result.get('freq_output') else None,
                    "log_file": str(intermediate_result.get('opt_output')) if intermediate_result.get('opt_output') else None,
                    "chk_file": None,
                    "fchk_file": None,
                    "qm_output_file": str(intermediate_result.get('sp_output')) if intermediate_result.get('sp_output') else None,
                }
                self._save_s3_resume_state(output_dir, resume_state)
            else:
                self.logger.info("步骤 2/3: Reactant 优化 + 频率 + L2 SP...")
                pm.enter_phase(ui_step_id, "Reactant Optimization")
                pm.set_subtask(ui_step_id, "Reactant Opt", 2, 3)
                pm.update_step("s3", description="S3: Reactant Opt")
                ui_log("Starting Reactant optimization + frequency + L2 SP")
                
                reactant_opt_config = self.step3_config.get('reactant_opt', {})
                reactant_charge = reactant_opt_config.get('charge', 0)
                reactant_mult = reactant_opt_config.get('multiplicity', 1)
                enable_nbo = reactant_opt_config.get('enable_nbo', False)
                
                self.logger.info(f"尝试标准 Reactant 优化 (charge={reactant_charge}, mult={reactant_mult})...")
                ui_log(f"Reactant optimization parameters: charge={reactant_charge}, mult={reactant_mult}")
                
                reactant_input = intermediate_optimized_xyz if intermediate_optimized_xyz else intermediate
                if intermediate_optimized_xyz:
                    self.logger.info(f"  -> 使用中间体优化结果: {intermediate_optimized_xyz}")
                else:
                    self.logger.info(f"  -> 使用原始输入 (中间体优化未成功): {intermediate}")
                
                reactant_opt_result = self.qc_runner.run_opt_sp_cycle(
                    xyz_file=reactant_input,
                    output_dir=reactant_output_dir,
                    charge=reactant_charge,
                    spin=reactant_mult,
                    enable_nbo=enable_nbo,
                    old_checkpoint=None
                )

                if reactant_opt_result.converged and self._verify_reactant_result(reactant_opt_result, reactant_output_dir / "standard"):
                    ui_log("Reactant optimization converged and validated")
                    resume_state["reactant_opt"] = {
                        "completed": True,
                        "optimized_xyz": str(reactant_opt_result.optimized_xyz),
                        "l2_energy": reactant_opt_result.l2_energy,
                        "opt_energy": reactant_opt_result.opt_energy,
                        "imaginary_count": reactant_opt_result.imaginary_count,
                        "method_used": reactant_opt_result.method_used,
                        "freq_log": str(reactant_opt_result.freq_log) if reactant_opt_result.freq_log else None,
                        "log_file": str(reactant_opt_result.log_file) if reactant_opt_result.log_file else None,
                        "chk_file": str(reactant_opt_result.chk_file) if reactant_opt_result.chk_file else None,
                        "fchk_file": str(reactant_opt_result.fchk_file) if reactant_opt_result.fchk_file else None,
                        "qm_output_file": str(reactant_opt_result.qm_output_file) if reactant_opt_result.qm_output_file else None,
                    }
                    self._save_s3_resume_state(output_dir, resume_state)

        if not reactant_opt_result.converged:
            ui_log(f"Reactant optimization failed: {reactant_opt_result.error_message}")
            with open(output_dir / ".rph_step_status.json", "w") as f:
                json.dump({
                    "step": "s3",
                    "status": "failed",
                    "error": reactant_opt_result.error_message,
                    "end_time": datetime.datetime.now().isoformat()
                }, f, indent=2)
                raise RuntimeError(f"Reactant 优化失败（含rescue）: {reactant_opt_result.error_message}")
        
        assert reactant_opt_result.l2_energy is not None, "Reactant L2 能量缺失"
        e_reactant_l2 = reactant_opt_result.l2_energy
        reactant_optimized_xyz = reactant_opt_result.optimized_xyz
        self.logger.info(f"✓ Reactant 优化完成，L2 能量: {e_reactant_l2:.6f} Hartree")
        ui_log(f"Reactant optimization completed, L2={e_reactant_l2:.6f} Hartree")
        
        if e_product_l2 is None:
            self.logger.warning("Product L2 能量缺失，执行补算 SP。")
            product_sp_result = self._run_sp(product, output_dir / "product_sp")
            assert product_sp_result.energy is not None
            e_product_l2 = cast(float, product_sp_result.energy)

        # 3. 片段处理 (基于 TS 优化后的几何)
        self.logger.info("步骤 3/3: 构建 SP 矩阵报告 (含片段)...")
        pm.enter_phase(ui_step_id, "SP Matrix Build")
        pm.set_subtask(ui_step_id, "SP Matrix", 3, 3)
        pm.update_step("s3", description="S3: Building SP Matrix")
        ui_log("Building SP matrix report")
        ts_l2_result = ts_result.l2_sp_result
        if ts_l2_result is None or not ts_l2_result.converged:
            raise RuntimeError("TS L2 SP 结果缺失或未收敛，无法构建 SP 矩阵报告。")

        sp_report = self._build_sp_matrix(
            ts_final_xyz=ts_result.optimized_xyz,
            intermediate=intermediate,
            e_product_l2=e_product_l2,
            e_intermediate_l2=e_reactant_l2,
            output_dir=output_dir / "ASM_SP_Mat",
            forming_bonds=forming_bonds,
            ts_l2_result=ts_l2_result,
            ts_l2_dir=output_dir / "ts_opt" / "L2_SP"
        )
        # 更新报告中的 Reactant 能量（默认是复合物能量）
        sp_report.e_reactant = e_reactant_l2

        # 计算 Shermo Gibbs 自由能 (TS 和 Reactant)
        self.logger.info("计算 Shermo Gibbs 自由能 (SP energy + freq log)...")
        pm.enter_phase(ui_step_id, "Thermochemistry")
        ui_log("Computing Shermo Gibbs energies")
        self._compute_shermo_gibbs(
            sp_report=sp_report,
            ts_result=ts_result,
            reactant_result=reactant_opt_result,
            ts_l2_energy=ts_l2_result.energy,
            reactant_l2_energy=e_reactant_l2,
            output_dir=output_dir
        )

        with open(output_dir / "sp_matrix_metadata.json", "w") as f:
            json.dump({
                "e_ts": sp_report.e_ts_final,
                "e_reactant": sp_report.e_reactant,
                "e_product": sp_report.e_product,
                "g_ts": sp_report.g_ts,
                "g_reactant": sp_report.g_reactant,
                "g_ts_source": sp_report.g_ts_source,
                "g_ts_error": sp_report.g_ts_error,
                "g_reactant_source": sp_report.g_reactant_source,
                "g_reactant_error": sp_report.g_reactant_error,
                "activation_energy_kcal": sp_report.get_activation_energy(),
                "reaction_energy_kcal": sp_report.get_reaction_energy(),
                "method": sp_report.method,
                "solvent": sp_report.solvent,
                "timestamp": datetime.datetime.now().isoformat()
            }, f, indent=2)

        dipolar_output = self._ensure_dipolar_output(output_dir, output_dir / "ts_opt" / "L2_SP")
        if dipolar_output is None:
            self.logger.warning("Dipolar output not found; artifacts_index will omit dipolar entry")
            ui_log("Dipolar output missing; artifacts index omits dipolar entry")

        # V6.3: Write artifacts_index.json
        SPMatrixReport._write_artifacts_index(output_dir, self.config)
        ui_log("Wrote artifacts_index.json")

        ts_final_canonical: Optional[Path] = None
        try:
            import shutil
            ts_final_canonical = output_dir / "ts_final.xyz"
            shutil.copy2(ts_result.optimized_xyz, ts_final_canonical)
            if isinstance(resume_state.get("ts_opt"), dict):
                resume_state["ts_opt"]["optimized_xyz"] = str(ts_final_canonical)
                self._save_s3_resume_state(output_dir, resume_state)

            from rph_core.utils.checkpoint_manager import CheckpointManager
            cm = CheckpointManager(output_dir.parent)
            sig = cm._compute_step3_signature(self.config)
            hashes = {
                "ts_guess": cm.compute_file_hash(ts_guess) or "",
                "intermediate": cm.compute_file_hash(intermediate) or "",
                "product": cm.compute_file_hash(product) or ""
            }
            with open(output_dir / "step3_signature.json", "w") as f:
                json.dump({"step3_signature": sig, "input_hashes": hashes}, f, indent=2)
        except Exception as e:
            self.logger.warning(f"Failed to write rehydration artifacts: {e}")
            ui_log(f"Failed to write rehydration artifacts: {e}")

        with open(output_dir / ".rph_step_status.json", "w") as f:
            json.dump({
                "step": "s3",
                "status": "completed",
                "end_time": datetime.datetime.now().isoformat(),
                "description": "S3: Transition Analysis Completed Successfully"
            }, f, indent=2)

        # 返回结果
        result = TransitionAnalysisResult(
            ts_final_xyz=ts_final_canonical if ts_final_canonical is not None else ts_result.optimized_xyz,
            ts_checkpoint=ts_result.optimized_xyz.with_suffix('.chk') if ts_result.optimized_xyz.with_suffix('.chk').exists() else None,
            sp_report=sp_report,
            method_used=ts_result.method_used,
            ts_fchk=ts_result.fchk_file,
            ts_log=ts_result.log_file,
            ts_qm_output=ts_result.qm_output_file,
            intermediate_fchk=reactant_opt_result.fchk_file,
            intermediate_log=reactant_opt_result.freq_log or reactant_opt_result.log_file,
            intermediate_qm_output=reactant_opt_result.qm_output_file
        )
        
        self.logger.info("=" * 60)
        self.logger.info(f"ΔG‡ (Complex-based): {sp_report.get_activation_energy():.2f} kcal/mol")
        self.logger.info("Step 3 完成")
        self.logger.info("=" * 60)
        pm.enter_phase(ui_step_id, "Completed")
        pm.update_step("s3", completed=100, description="S3: Completed")
        ui_log("Step 3 completed successfully")
        
        return result

    def _build_sp_matrix(
        self,
        ts_final_xyz: Path,
        intermediate: Path,
        e_product_l2: float,
        e_intermediate_l2: float,
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

        e_frag_a_l2 = 0.0
        e_frag_b_l2 = 0.0
        e_frag_a_relaxed = None
        e_frag_b_relaxed = None
        fragment_split_source = None
        fragment_split_reason = None
        fragment_indices = None

        # 构建报告
        report = SPMatrixReport(
            e_product=e_product_l2,
            e_reactant=e_intermediate_l2,
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
        except Exception as exc:
            self.logger.debug(f"Symlink alias failed, fallback to copytree: {exc}")

        try:
            shutil.copytree(l2_sp_dir, ts_sp_dir)
            (ts_sp_dir / "provenance.json").write_text(
                '{"mode": "copy", "copied_from": ' + json.dumps(str(l2_sp_dir)) + '}'
            )
        except Exception as e:
            self.logger.warning(f"Failed to alias ts_sp dir (from {l2_sp_dir} to {ts_sp_dir}): {e}")

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

    def _compute_shermo_gibbs(
        self,
        sp_report: SPMatrixReport,
        ts_result: Any,
        reactant_result: Any,
        ts_l2_energy: float,
        reactant_l2_energy: float,
        output_dir: Path
    ) -> None:
        from rph_core.utils.shermo_runner import run_shermo
        from pathlib import Path

        shermo_dir = output_dir / "shermo"
        shermo_dir.mkdir(parents=True, exist_ok=True)

        shermo_bin_path = self.config.get("executables", {}).get("shermo", {}).get("path", "Shermo")
        shermo_bin = Path(shermo_bin_path) if shermo_bin_path else Path("Shermo")

        thermo_config = self.config.get("thermo", {})
        temperature_k = thermo_config.get("temperature_k", 298.15)
        pressure_atm = thermo_config.get("pressure_atm", 1.0)
        scl_zpe = thermo_config.get("scl_zpe", 0.9905)
        ilowfreq = thermo_config.get("ilowfreq", 2)
        imagreal = thermo_config.get("imagreal", None)
        conc = thermo_config.get("conc", None)

        sp_report.g_ts = None
        sp_report.g_ts_source = None
        sp_report.g_ts_error = None

        freq_log_ts = ts_result.freq_log if hasattr(ts_result, "freq_log") and ts_result.freq_log else ts_result.log_file
        if freq_log_ts and freq_log_ts.exists():
            try:
                ts_shermo_out = shermo_dir / "ts_Shermo.sum"
                thermo_ts = run_shermo(
                    shermo_bin=shermo_bin,
                    freq_output=freq_log_ts,
                    sp_energy=ts_l2_energy,
                    output_file=ts_shermo_out,
                    temperature_k=temperature_k,
                    pressure_atm=pressure_atm,
                    scl_zpe=scl_zpe,
                    ilowfreq=ilowfreq,
                    imagreal=imagreal,
                    conc=conc
                )
                g_ts_h = thermo_ts.g_conc if thermo_ts.g_conc is not None else thermo_ts.g_sum
                sp_report.g_ts = g_ts_h * HARTREE_TO_KCAL
                sp_report.g_ts_source = "shermo_sp_freq"
                self.logger.info(
                    f"✓ TS Gibbs 计算成功: {sp_report.g_ts:.6f} kcal/mol (来源: Shermo SP+Freq)"
                )
            except Exception as e:
                sp_report.g_ts_error = str(e)
                sp_report.g_ts_source = "shermo_failed"
                self.logger.warning(f"Shermo TS Gibbs 计算失败: {e}")
        else:
            sp_report.g_ts_error = "freq_log not found"
            sp_report.g_ts_source = "missing_freq_log"
            self.logger.warning(f"TS freq log 未找到，跳过 Gibbs 计算")

        sp_report.g_reactant = None
        sp_report.g_reactant_source = None
        sp_report.g_reactant_error = None

        freq_log_reactant = (
            reactant_result.freq_log
            if hasattr(reactant_result, "freq_log") and reactant_result.freq_log
            else reactant_result.log_file
        )
        if freq_log_reactant and freq_log_reactant.exists():
            try:
                reactant_shermo_out = shermo_dir / "reactant_Shermo.sum"
                thermo_reactant = run_shermo(
                    shermo_bin=shermo_bin,
                    freq_output=freq_log_reactant,
                    sp_energy=reactant_l2_energy,
                    output_file=reactant_shermo_out,
                    temperature_k=temperature_k,
                    pressure_atm=pressure_atm,
                    scl_zpe=scl_zpe,
                    ilowfreq=ilowfreq,
                    imagreal=imagreal,
                    conc=conc
                )
                g_reactant_h = (
                    thermo_reactant.g_conc
                    if thermo_reactant.g_conc is not None
                    else thermo_reactant.g_sum
                )
                sp_report.g_reactant = g_reactant_h * HARTREE_TO_KCAL
                sp_report.g_reactant_source = "shermo_sp_freq"
                self.logger.info(
                    f"✓ Reactant Gibbs 计算成功: {sp_report.g_reactant:.6f} kcal/mol (来源: Shermo SP+Freq)"
                )
            except Exception as e:
                sp_report.g_reactant_error = str(e)
                sp_report.g_reactant_source = "shermo_failed"
                self.logger.warning(f"Shermo Reactant Gibbs 计算失败: {e}")
        else:
            sp_report.g_reactant_error = "freq_log not found"
            sp_report.g_reactant_source = "missing_freq_log"
            self.logger.warning(f"Reactant freq log 未找到，跳过 Gibbs 计算")

    def run(self, *args, **kwargs):
        # 保持旧接口兼容性（内部调用新逻辑）
        return self.run_with_qctaskrunner(*args, **kwargs)
