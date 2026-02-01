"""
Reaction Profile Hunter Orchestrator
======================================

串行四步走架构的总指挥 + Step 3.5 SP矩阵集成

Author: QCcalc Team
Date: 2026-01-09
Session: #13 - 集成 S3.5 SP矩阵
"""

import logging
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple
from datetime import datetime

from rph_core.utils.log_manager import setup_logger
from rph_core.utils.path_compat import normalize_path, is_toxic_path
from rph_core.utils.task_builder import build_tasks_from_run_config, sanitize_rx_id
from rph_core.utils.config_loader import load_config
from rph_core.utils.small_molecule_catalog import SmallMoleculeCatalog
from rph_core.utils.checkpoint_manager import CheckpointManager
from rph_core.utils.checkpoint_manager import PipelineState
from rph_core.utils.optimization_config import normalize_qc_config
from rph_core.version import __version__
from rph_core.utils import ui, notify


@dataclass
class PipelineResult:
    """流水线结果"""
    success: bool
    product_smiles: Optional[str] = None
    work_dir: Optional[Path] = None

    # Step outputs
    product_xyz: Optional[Path] = None
    e_product_l2: Optional[float] = None
    product_checkpoint: Optional[Path] = None
    product_thermo: Optional[Path] = None
    product_fchk: Optional[Path] = None
    product_log: Optional[Path] = None
    product_qm_output: Optional[Path] = None
    ts_guess_xyz: Optional[Path] = None
    reactant_xyz: Optional[Path] = None
    ts_final_xyz: Optional[Path] = None
    features_csv: Optional[Path] = None

    # Metadata from steps
    forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]] = None
    sp_matrix_report: Optional[object] = None  # 来自 S3 的完整 SP 报告
    ts_fchk: Optional[Path] = None
    ts_log: Optional[Path] = None
    ts_qm_output: Optional[Path] = None
    reactant_fchk: Optional[Path] = None
    reactant_log: Optional[Path] = None
    reactant_qm_output: Optional[Path] = None

    # Error tracking
    error_step: Optional[str] = None
    error_message: Optional[str] = None

    def __str__(self):
        if self.success:
            l2_info = f", L2: {self.e_product_l2:.6f} Ha" if self.e_product_l2 else ""
            return f"✅ Pipeline 成功: {self.product_smiles}\n" \
                   f"   Product: {self.product_xyz}{l2_info}\n" \
                   f"   TS Final: {self.ts_final_xyz}\n" \
                   f"   Features: {self.features_csv}"
        else:
            return f"❌ Pipeline 失败: {self.error_step}\n" \
                   f"   错误: {self.error_message}"


class ReactionProfileHunter:
    """
    Reaction Profile Hunter v3.0 (分子自治架构)

    设计模式:
    - Orchestrator 只是"调度员"，不是"工人"
    - 每个 Step 是独立的"工人"，有明确的输入输出
    - v3.0 核心改进：分子自治目录结构 + OPT-SP 耦合循环

    串行流程:
    S1: AnchorPhase (分子锚定 + CREST + DFT OPT-SP 耦合)
      → S2: RetroScanner (逆向扫描，从 S1_ConfGeneration/[Molecule]/dft 读取)
      → S3: TSOptimizer (TS优化，使用 LogParser 提取坐标)
      → S4: FeatureMiner (特征提取)

    v3.0 目录结构:
    S1_ConfGeneration/[Molecule_Name]/
        ├── crest/          # CREST 搜索结果
        ├── dft/            # DFT OPT + SP (无子目录，扁平结构)
        └── [Molecule_Name]_global_min.xyz
    """

    def __init__(self, config_path: Optional[Path] = None, log_level: Optional[str] = None):
        """
        初始化 Reaction Profile Hunter

        Args:
            config_path: 配置文件路径（可选，默认使用defaults.yaml）
        """
        # 加载配置
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "defaults.yaml"

        self.config = load_config(config_path)
        if log_level:
            self.config.setdefault("global", {})["log_level"] = log_level
        self.config, qc_fixes = normalize_qc_config(self.config, auto_fix=True)
        self.logger = setup_logger(
            "ReactionProfileHunter",
            level=self.config.get('global', {}).get('log_level', 'INFO')
        )
        if qc_fixes:
            self.logger.warning(f"QC config normalized: {len(qc_fixes)} change(s)")
            for fix in qc_fixes:
                self.logger.debug(
                    f"QC config fix [{fix['field']}]: {fix['original']} -> {fix['updated']}"
                )

        ui.print_pipeline_header(__version__)
        self.logger.info(f"Reaction Profile Hunter v{__version__} 初始化 (含 S3.5)")
        
        # Initialize small molecule catalog
        self.small_mol_catalog = SmallMoleculeCatalog(self.config)

        # 延迟初始化各步骤引擎（懒加载）
        self._s1_engine = None
        self._s1_engine_type = "ProductAnchor"  # 默认类型
        self._s2_engine = None
        self._s3_engine = None
        self._s4_engine = None

    @property
    def s1_engine(self):
        """Step 1 引擎（懒加载）- v3.0: 分子自治架构"""
        if self._s1_engine is None:
            # v3.0 强制使用新的 AnchorPhase
            try:
                from rph_core.steps.anchor.handler import AnchorPhase
                self._s1_engine = AnchorPhase(
                    config=self.config,
                    base_work_dir=Path.cwd()  # 默认值，将在运行时更新
                )
                self._s1_engine_type = "AnchorPhase_v3"
                self.logger.debug("✓ Step 1 使用 AnchorPhase v3.0（分子自治架构）")
                return self._s1_engine
            except ImportError as e:
                self.logger.error(f"无法导入 AnchorPhase v3.0: {e}")
                self.logger.error("v3.0 要求必须使用新的 AnchorPhase")
                raise RuntimeError(
                    "ReactionProfileHunter v3.0 要求必须使用 AnchorPhase v3.0。"
                    "请确保 rph_core/steps/anchor/handler.py 存在且可导入。"
                )
        return self._s1_engine

    @property
    def s2_engine(self):
        """Step 2 引擎 (懒加载)"""
        if self._s2_engine is None:
            from rph_core.steps.step2_retro import RetroScanner
            self._s2_engine = RetroScanner(self.config)
            self.logger.debug("✓ Step 2 (RetroScanner) 已初始化")
        return self._s2_engine

    @property
    def s3_engine(self):
        """Step 3 引擎 (懒加载)"""
        if self._s3_engine is None:
            from rph_core.steps.step3_opt import TSOptimizer
            self._s3_engine = TSOptimizer(self.config)
            self.logger.debug("✓ Step 3 (TSOptimizer) 已初始化")
        return self._s3_engine



    @property
    def s4_engine(self):
        """Step 4 引擎 (懒加载)"""
        if self._s4_engine is None:
            from rph_core.steps.step4_features import FeatureMiner
            self._s4_engine = FeatureMiner(self.config)
            self.logger.debug("✓ Step 4 (FeatureMiner) 已初始化")
        return self._s4_engine

    def run_pipeline(
        self,
        product_smiles: str,
        work_dir: Path,
        skip_steps: Optional[list] = None,
        precursor_smiles: Optional[str] = None,
        leaving_group_key: Optional[str] = None
    ) -> PipelineResult:
        """
        执行 v2.1 串行四步走架构

        数据流:
        S1_Output (Product_Min) ──┬──> S2_Input
                                  │
        S2_Output (TS_Guess, R) ──┼──> S3_Input
                                  │
        S3_Output (TS_Final) ─────┼──> S4_Input
                                  │
        S1 + S2 + S3 ─────────────┴──> S4_Input

        Args:
            product_smiles: 产物 SMILES
            work_dir: 工作目录
            skip_steps: 要跳过的步骤列表（用于调试）

        Returns:
            PipelineResult 对象
        """
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"🚀 任务启动: {product_smiles}")
        self.logger.info(f"📁 工作目录: {work_dir}")

        result = PipelineResult(
            success=False,
            product_smiles=product_smiles,
            work_dir=work_dir
        )

        def _notify(success: bool, error_step: Optional[str] = None, error_message: Optional[str] = None) -> None:
            title = "RPH 任务完成" if success else "RPH 任务失败"
            if success:
                message = f"{product_smiles} 完成: {work_dir}"
            else:
                step = error_step or "Unknown"
                message = f"{product_smiles} 失败于 {step}: {error_message or 'unknown error'}"
            notify.notify_completion(title, message, self.config)

        skip_steps = skip_steps or []

        resume_enabled = bool((self.config.get("run", {}) or {}).get("resume", True))
        checkpoint_mgr = CheckpointManager(work_dir)
        if resume_enabled:
            state = checkpoint_mgr.load_state()
            if state is None:
                state = PipelineState(
                    product_smiles=product_smiles,
                    work_dir=str(work_dir),
                    start_time=datetime.now().isoformat(),
                    last_update=datetime.now().isoformat(),
                    steps={},
                    config_snapshot={},
                )
                checkpoint_mgr.save_state(state)

        # Resume: reuse Step1/Step2 outputs to avoid repeated heavy QC on S3 failures.
        if resume_enabled and 's1' not in skip_steps and checkpoint_mgr.is_step_completed('s1'):
            product_xyz = checkpoint_mgr.get_step_output('s1', 'product_xyz')
            if product_xyz and Path(product_xyz).exists():
                result.product_xyz = Path(product_xyz)
                product_sp = checkpoint_mgr.get_step_metadata('s1', 'e_product_sp')
                if product_sp is not None:
                    result.e_product_l2 = float(product_sp)
                fchk = checkpoint_mgr.get_step_output('s1', 'product_fchk')
                if fchk and Path(fchk).exists():
                    result.product_fchk = Path(fchk)
                logp = checkpoint_mgr.get_step_output('s1', 'product_log')
                if logp and Path(logp).exists():
                    result.product_log = Path(logp)
                qmp = checkpoint_mgr.get_step_output('s1', 'product_qm_output')
                if qmp and Path(qmp).exists():
                    result.product_qm_output = Path(qmp)
                chk = checkpoint_mgr.get_step_output('s1', 'product_checkpoint')
                if chk and Path(chk).exists():
                    result.product_checkpoint = Path(chk)
                thermo = checkpoint_mgr.get_step_output('s1', 'product_thermo')
                if thermo and Path(thermo).exists():
                    result.product_thermo = Path(thermo)

                self.logger.info(f"✅ Resume: Step1 already complete, reuse product: {result.product_xyz}")
                skip_steps.append('s1')

        if resume_enabled and 's2' not in skip_steps and checkpoint_mgr.is_step_completed('s2'):
            ts_guess = checkpoint_mgr.get_step_output('s2', 'ts_guess_xyz')
            reactant_xyz = checkpoint_mgr.get_step_output('s2', 'reactant_xyz')
            if ts_guess and reactant_xyz and Path(ts_guess).exists() and Path(reactant_xyz).exists():
                result.ts_guess_xyz = Path(ts_guess)
                result.reactant_xyz = Path(reactant_xyz)
                self.logger.info(
                    f"✅ Resume: Step2 already complete, reuse ts_guess/reactant: {result.ts_guess_xyz}, {result.reactant_xyz}"
                )
                skip_steps.append('s2')

        # === Step 1: 产物锚定 (Product Anchor - v3.0 分子自治架构) ===
        if 's1' not in skip_steps:
            try:
                ui.print_step_header("Step 1", "Product Anchor", "Global Minimum Search (v3.0)")
                self.logger.info(">>> Step 1: 寻找产物全局最低构象 (v3.0 OPT-SP 耦合)...")

                # v3.0: 使用 AnchorPhase 处理产物
                s1_work_dir = work_dir / "S1_ConfGeneration"

                # 设置 AnchorPhase 的工作目录
                self.s1_engine.base_work_dir = s1_work_dir

                # Resolve leaving group if provided
                resolved_lg_smiles = None
                if leaving_group_key:
                    mol_obj = self.small_mol_catalog.get(leaving_group_key)
                    if mol_obj:
                        resolved_lg_smiles = mol_obj.smiles
                    else:
                        self.logger.warning(f"Leaving group key '{leaving_group_key}' not found in catalog. Skipping.")

                # Build molecules dictionary for S1
                molecules = {"product": product_smiles}
                if precursor_smiles:
                    molecules["precursor"] = precursor_smiles
                if resolved_lg_smiles:
                    molecules["leaving_group"] = resolved_lg_smiles

                with ui.status("Running AnchorPhase (CREST + DFT)..."):
                    anchor_result = self.s1_engine.run(
                        molecules=molecules
                    )

                # 检查执行结果
                if not anchor_result.success:
                    raise RuntimeError(f"AnchorPhase v3.0 失败: {anchor_result.error_message}")

                # 从 AnchorPhaseResult 中提取产物数据（v3.0 结构）
                product_data = anchor_result.anchored_molecules.get("product", {})

                # v3.0 新结构：xyz 是 SP 输出文件路径，e_sp 是 SP 能量
                product_sp_out = product_data.get("xyz")
                e_product_sp = product_data.get("e_sp")

                if product_sp_out is None or e_product_sp is None:
                    raise RuntimeError(
                        f"AnchorPhase 未返回完整的产物数据。"
                        f"product_data = {product_data}"
                    )

                # 使用 LogParser 从 SP 输出中提取最终坐标
                from rph_core.utils.geometry_tools import LogParser
                coords, symbols, error = LogParser.extract_last_converged_coords(
                    product_sp_out,
                    engine_type='auto'
                )

                if coords is None:
                    self.logger.warning(f"无法从 {product_sp_out} 提取坐标: {error}")
                    # 回退：直接使用 SP 输出路径
                    product_min_xyz = product_sp_out
                else:
                    # 创建最终的产物 XYZ 文件
                    product_min_xyz = s1_work_dir / "product_min.xyz"
                    from rph_core.utils.file_io import write_xyz
                    # 确保 symbols 不为 None
                    if symbols is None:
                        self.logger.warning("未提取到符号，从 SP 输出文件读取")
                        from rph_core.utils.file_io import read_xyz
                        _, fallback_symbols = read_xyz(product_sp_out)
                        symbols = fallback_symbols
                    write_xyz(product_min_xyz, coords, symbols, title=f"Product SP E={e_product_sp:.6f}")

                # 保存结果（v3.0 使用 e_sp 而不是 e_l2）
                result.product_xyz = product_min_xyz
                result.e_product_l2 = e_product_sp  # 保持向后兼容，使用 SP 能量
                result.product_fchk = product_data.get("fchk")
                result.product_log = product_data.get("log")
                result.product_qm_output = product_data.get("qm_output")

                product_thermo_file = s1_work_dir / "product" / "dft" / "conformer_thermo.csv"
                if product_thermo_file.exists():
                    result.product_thermo = product_thermo_file

                # 保存 checkpoint 路径（如果存在）
                product_checkpoint = product_data.get("chk")
                if product_checkpoint and product_checkpoint.exists():
                    result.product_checkpoint = product_checkpoint
                    self.logger.info(f"    ✓ S1 checkpoint 可用: {product_checkpoint}")

                self.logger.info(f"    ✓ 产物锚定完成: {product_min_xyz}")
                self.logger.info(f"    ✓ SP 能量: {e_product_sp:.8f} Hartree")
                if product_sp_out.suffix == ".xyz":
                    self.logger.info(f"    ✓ SP 输出文件(已是global_min): {product_sp_out}")
                else:
                    self.logger.info(f"    ✓ SP 输出文件: {product_sp_out}")

                if resume_enabled:
                    checkpoint_mgr.mark_step_completed(
                        "s1",
                        output_files={
                            "product_xyz": str(result.product_xyz),
                            "product_fchk": str(result.product_fchk) if result.product_fchk else "",
                            "product_log": str(result.product_log) if result.product_log else "",
                            "product_qm_output": str(result.product_qm_output) if result.product_qm_output else "",
                            "product_checkpoint": str(result.product_checkpoint) if result.product_checkpoint else "",
                            "product_thermo": str(result.product_thermo) if result.product_thermo else "",
                        },
                        metadata={"e_product_sp": float(e_product_sp)},
                    )

            except Exception as e:
                result.error_step = "Step1_ProductAnchor_v3"
                result.error_message = str(e)
                self.logger.error(f"Step 1 v3.0 失败: {e}", exc_info=True)
                _notify(False, result.error_step, result.error_message)
                return result
        else:
            self.logger.warning("⚠️  跳过 Step 1")
            s1_candidates = [
                work_dir / "S1_ConfGeneration",
                work_dir.parent / "S1_test" / "S1_ConfGeneration",
                work_dir.parent / "S1_ConfGeneration",
                work_dir / "S1_Product",
                work_dir.parent / "S1_test" / "S1_Product",
                work_dir.parent / "S1_Product"
            ]
            for candidate in s1_candidates:
                if candidate.exists():
                    result.product_xyz = candidate
                    self.logger.info(f"    ✓ 复用 S1 输出目录: {candidate}")
                    break

        # === Step 2: 逆向扫描 (Retro Scanner) ===
        if 's2' not in skip_steps and result.product_xyz:
            try:
                ui.print_step_header("Step 2", "Retro Scanner", "Generating TS Guess & Reactant Complex")
                self.logger.info(">>> Step 2: 逆向生成 TS 初猜与底物...")
                
                with ui.status("Scanning bond coordinates (xTB)..."):
                    ts_guess_xyz, reactant_xyz, bonds = self.s2_engine.run(
                        product_xyz=result.product_xyz,
                        output_dir=work_dir / "S2_Retro"
                    )
                result.ts_guess_xyz = ts_guess_xyz
                result.reactant_xyz = reactant_xyz
                result.forming_bonds = bonds
                self.logger.info(f"    ✓ TS 初猜: {ts_guess_xyz}")
                self.logger.info(f"    ✓ 底物复合物: {reactant_xyz}")

                if resume_enabled:
                    checkpoint_mgr.mark_step_completed(
                        "s2",
                        output_files={
                            "ts_guess_xyz": str(result.ts_guess_xyz),
                            "reactant_xyz": str(result.reactant_xyz),
                        },
                        metadata={"forming_bonds": str(bonds)},
                    )
            except Exception as e:
                result.error_step = "Step2_RetroScanner"
                result.error_message = str(e)
                self.logger.error(f"Step 2 失败: {e}", exc_info=True)
                _notify(False, result.error_step, result.error_message)
                return result
        else:
            self.logger.warning("⚠️  跳过 Step 2")

        # === Step 3: 反应分析 (Transition Analyzer) ===
        if 's3' not in skip_steps and result.ts_guess_xyz:
            try:
                ui.print_step_header("Step 3", "Transition Analyzer", "TS Optimization & Verification")
                self.logger.info(">>> Step 3: 反应中心全分析 (TS优化 + Reactant/Fragments SP)...")

                # 传递S1的checkpoint以复用轨道
                old_checkpoint = result.product_checkpoint
                if old_checkpoint:
                    self.logger.info(f"  复用S1 checkpoint: {old_checkpoint.name}")

                if result.product_xyz and result.product_xyz.is_dir():
                    product_dir = result.product_xyz
                    candidates = [
                        product_dir / "product_min.xyz",
                        product_dir / "product" / "product_global_min.xyz",
                        product_dir / "product_global_min.xyz",
                        product_dir / "product" / "global_min.xyz",
                        product_dir / "global_min.xyz"
                    ]
                    product_file = next((p for p in candidates if p.exists()), None)
                    if product_file is None:
                        raise RuntimeError(
                            f"无法在 {product_dir} 中找到产物结构文件用于 S3。"
                            f"已尝试: {[str(p) for p in candidates]}"
                        )
                    result.product_xyz = product_file
                    self.logger.info(f"  ✓ 使用产物文件: {result.product_xyz}")

                if result.reactant_xyz is None or result.product_xyz is None:
                    raise RuntimeError("Step3 输入缺失: reactant 或 product 为 None")

                with ui.status("Optimizing Transition State (Berny/QST2)..."):
                    s3_result = self.s3_engine.run(
                        ts_guess=result.ts_guess_xyz,
                        reactant=result.reactant_xyz,
                        product=result.product_xyz,
                        output_dir=work_dir / "S3_TransitionAnalysis",
                        e_product_l2=result.e_product_l2,
                        product_thermo=result.product_thermo,
                        forming_bonds=result.forming_bonds,
                        old_checkpoint=old_checkpoint
                    )

                result.ts_final_xyz = s3_result.ts_final_xyz
                result.sp_matrix_report = s3_result.sp_report
                
                result.ts_fchk = s3_result.ts_fchk
                result.ts_log = s3_result.ts_log
                result.ts_qm_output = s3_result.ts_qm_output
                result.reactant_fchk = s3_result.reactant_fchk
                result.reactant_log = s3_result.reactant_log
                result.reactant_qm_output = s3_result.reactant_qm_output
                
                self.logger.info(f"    ✓ S3 完成")
                dg_act = s3_result.sp_report.get_activation_energy()
                dg_rxn = s3_result.sp_report.get_reaction_energy()
                if dg_act is not None:
                    self.logger.info(f"      ΔG‡ = {dg_act:.3f} kcal/mol")
                else:
                    self.logger.info("      ΔG‡ = N/A")
                if dg_rxn is not None:
                    self.logger.info(f"      ΔG_rxn = {dg_rxn:.3f} kcal/mol")
                else:
                    self.logger.info("      ΔG_rxn = N/A")
            except Exception as e:
                result.error_step = "Step3_TransitionAnalyzer"
                result.error_message = str(e)
                self.logger.error(f"Step 3 失败: {e}", exc_info=True)
                _notify(False, result.error_step, result.error_message)
                return result
        else:
            self.logger.warning("⚠️  跳过 Step 3")

        # === Step 4: 特征挖掘 (Feature Miner) ===
        if 's4' not in skip_steps and result.ts_final_xyz:
            try:
                ui.print_step_header("Step 4", "Feature Miner", "Extracting Features (Extract-Only)")
                self.logger.info(">>> Step 4: 提取物理有机特征...")

                # Check S3 artifacts for S4 and warn if missing
                if result.ts_fchk is None:
                    self.logger.warning("TS fchk missing: formchk failed or not produced. S4 will degrade.")
                if result.reactant_fchk is None:
                    self.logger.warning("Reactant fchk missing: formchk failed or not produced. S4 will degrade.")
                if result.product_fchk is None:
                    self.logger.warning("Product fchk missing: formchk failed or not produced. S4 will degrade.")
                if result.ts_log is None and result.ts_qm_output is None:
                    self.logger.warning("TS log/out missing: Gaussian .log or ORCA .out not available. S4 will degrade.")
                if result.reactant_log is None and result.reactant_qm_output is None:
                    self.logger.warning("Reactant log/out missing: Gaussian .log or ORCA .out not available. S4 will degrade.")
                if result.product_log is None and result.product_qm_output is None:
                    self.logger.warning("Product log/out missing: Gaussian .log or ORCA .out not available. S4 will degrade.")

                if result.reactant_xyz is None or result.product_xyz is None:
                    raise RuntimeError("Step4 输入缺失: reactant 或 product 为 None")

                with ui.status("Extracting features (thermo, geom, qc)..."):
                    features_csv = self.s4_engine.run(
                        ts_final=result.ts_final_xyz,
                        reactant=result.reactant_xyz,
                        product=result.product_xyz,
                        output_dir=work_dir / "S4_Data",
                        forming_bonds=result.forming_bonds,
                        sp_matrix_report=result.sp_matrix_report,
                        ts_fchk=result.ts_fchk,
                        reactant_fchk=result.reactant_fchk,
                        product_fchk=result.product_fchk,
                        ts_orca_out=result.ts_qm_output,
                        reactant_orca_out=result.reactant_qm_output,
                        product_orca_out=result.product_qm_output
                    )
                result.features_csv = features_csv
                self.logger.info(f"    ✓ 特征提取完成: {features_csv}")
            except Exception as e:
                result.error_step = "Step4_FeatureMiner"
                result.error_message = str(e)
                self.logger.error(f"Step 4 失败: {e}", exc_info=True)
                _notify(False, result.error_step, result.error_message)
                return result
        else:
            self.logger.warning("⚠️  跳过 Step 4")

        # 成功完成
        result.success = True
        ui.print_result_summary(result)
        self.logger.info(f"✅ 任务完成! 数据已保存至: {work_dir}")
        _notify(True)

        return result

    def run_batch(
        self,
        smiles_list: list,
        work_dir: Path,
        max_workers: int = 1,
        skip_steps: Optional[list] = None
    ):
        if max_workers <= 0:
            max_workers = int((self.config.get("run", {}) or {}).get("max_workers", 1) or 1)
        skip_steps = list(skip_steps) if skip_steps is not None else []
        self.logger.info(f"开始批量处理: {len(smiles_list)} 个分子")
        self.logger.info(f"并发数: {max_workers}")

        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        from concurrent.futures import ProcessPoolExecutor, as_completed
        from tqdm import tqdm

        results = []
        tasks = []

        for idx, smiles in enumerate(smiles_list):
            task_work_dir = work_dir / f"task_{idx}_{smiles[:10]}"
            tasks.append((idx, smiles, task_work_dir, skip_steps))

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(self.run_pipeline, smiles, task_dir, skip_steps): (idx, smiles, task_dir)
                for idx, smiles, task_dir in tasks
            }

            # 使用 tqdm 显示进度
            with tqdm(total=len(tasks), desc="处理进度") as pbar:
                for future in as_completed(futures):
                    idx, smiles, _ = futures[future]
                    try:
                        result = future.result()
                        results.append(result)
                        status = "✅" if result.success else "❌"
                        pbar.set_postfix_str(f"{status} {idx}: {smiles[:15]}")
                    except Exception as e:
                        self.logger.error(f"Task {idx} 失败: {e}")
                        pbar.set_postfix_str(f"❌ {idx}: Error")

        self.logger.info(f"批量处理完成: {len(results)}/{len(smiles_list)} 成功")
        return results


# =============================================================================
# 命令行接口
# =============================================================================

def _resolve_run_config(config: dict, args) -> dict:
    run_cfg = dict(config.get("run", {}) or {})
    global_cfg = config.get("global", {}) or {}
    run_cfg.setdefault("source", "single")
    run_cfg.setdefault("output_root", global_cfg.get("work_dir_base", "./rph_output"))
    run_cfg.setdefault("workdir_naming", "rx_{rx_id}")
    run_cfg.setdefault("resume", True)
    run_cfg.setdefault("dry_run", False)
    run_cfg.setdefault("max_tasks", 0)
    run_cfg.setdefault("filter_ids", [])

    if args.output:
        run_cfg["output_root"] = args.output

    if args.smiles:
        run_cfg["source"] = "single"
        run_cfg["single"] = {
            "rx_id": "manual",
            "product_smiles": args.smiles,
        }

    return run_cfg


def _run_tasks(hunter: ReactionProfileHunter, run_cfg: dict) -> list[PipelineResult]:
    tasks = build_tasks_from_run_config(run_cfg)
    output_root = normalize_path(str(run_cfg.get("output_root", "./rph_output")))

    if is_toxic_path(output_root):
        hunter.logger.warning(
            f"Output path contains toxic characters: {output_root}. "
            "Consider using a safe directory without spaces/brackets."
        )

    output_root.mkdir(parents=True, exist_ok=True)

    results = []
    for task in tasks:
        rx_id = sanitize_rx_id(task.rx_id)
        work_dir = output_root / run_cfg["workdir_naming"].format(rx_id=rx_id)

        if run_cfg.get("dry_run", False):
            hunter.logger.info(f"[dry-run] {task.rx_id} -> {work_dir}")
            continue

        if run_cfg.get("resume", True):
            s4_dir = work_dir / "S4_Data"
            if s4_dir.exists():
                checkpoint_mgr = CheckpointManager(work_dir)
                if checkpoint_mgr.is_step4_complete(s4_dir, hunter.config):
                    hunter.logger.info(f"Skip {task.rx_id}: Step4 already complete")
                    continue

        result = hunter.run_pipeline(
            product_smiles=task.product_smiles,
            work_dir=work_dir,
            skip_steps=[],
            precursor_smiles=task.meta.get("precursor_smiles"),
            leaving_group_key=task.meta.get("leaving_small_molecule_key")
        )
        results.append(result)

    return results


def main():
    """命令行主入口"""
    import argparse

    parser = argparse.ArgumentParser(
        description="ReactionProfileHunter v2.1 - 过渡态搜索与特征提取"
    )
    parser.add_argument(
        '--smiles',
        type=str,
        default=None,
        help='产物 SMILES 字符串（可选；默认使用 config.run）'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='输出目录（覆盖 config.run.output_root）'
    )
    parser.add_argument(
        '--config',
        type=str,
        default=None,
        help='配置文件路径'
    )
    parser.add_argument(
        '--log-level',
        type=str,
        default=None,
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='日志级别（覆盖 config.global.log_level）'
    )

    args = parser.parse_args()

    try:
        hunter = ReactionProfileHunter(
            config_path=Path(args.config) if args.config else None,
            log_level=args.log_level
        )

        run_cfg = _resolve_run_config(hunter.config, args)
        results = _run_tasks(hunter, run_cfg)

        if not results and run_cfg.get("dry_run", False):
            return 0

        success_count = sum(1 for r in results if r.success)
        hunter.logger.info(f"批量处理完成: {success_count}/{len(results)} 成功")

        return 0 if success_count == len(results) else 1

    except Exception as e:
        logging.error(f"程序异常: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
