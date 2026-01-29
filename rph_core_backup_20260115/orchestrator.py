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
from dataclasses import dataclass, field
from typing import Optional, Tuple

from rph_core.utils.log_manager import setup_logger
from rph_core.utils.config_loader import load_config
from rph_core.utils.checkpoint_manager import CheckpointManager, load_checkpoint_state


@dataclass
class PipelineResult:
    """流水线结果"""
    success: bool
    product_smiles: str = None
    work_dir: Path = None

    # Step outputs
    product_xyz: Optional[Path] = None
    e_product_l2: Optional[float] = None
    product_checkpoint: Optional[Path] = None
    ts_guess_xyz: Optional[Path] = None
    reactant_xyz: Optional[Path] = None
    ts_final_xyz: Optional[Path] = None
    features_csv: Optional[Path] = None

    # Metadata from steps
    forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]] = None
    sp_matrix_report: Optional[object] = None  # 来自 S3 的完整 SP 报告

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
      → S2: RetroScanner (逆向扫描，从 S1_Product/[Molecule]/dft 读取)
      → S3: TSOptimizer (TS优化，使用 LogParser 提取坐标)
      → S4: FeatureMiner (特征提取)

    v3.0 目录结构:
    S1_Product/[Molecule_Name]/
        ├── crest/          # CREST 搜索结果
        ├── dft/            # DFT OPT + SP (无子目录，扁平结构)
        └── [Molecule_Name]_global_min.xyz
    """

    def __init__(self, config_path: Optional[Path] = None):
        """
        初始化 Reaction Profile Hunter

        Args:
            config_path: 配置文件路径（可选，默认使用defaults.yaml）
        """
        # 加载配置
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "defaults.yaml"

        self.config = load_config(config_path)
        self.logger = setup_logger(
            "ReactionProfileHunter",
            level=self.config.get('global', {}).get('log_level', 'INFO')
        )

        self.logger.info("=" * 60)
        self.logger.info("Reaction Profile Hunter v2.1 初始化 (含 S3.5)")
        self.logger.info("=" * 60)

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
        skip_steps: list = None
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

        skip_steps = skip_steps or []

        # === Step 1: 产物锚定 (Product Anchor - v3.0 分子自治架构) ===
        if 's1' not in skip_steps:
            try:
                self.logger.info(">>> Step 1: 寻找产物全局最低构象 (v3.0 OPT-SP 耦合)...")

                # v3.0: 使用 AnchorPhase 处理产物
                s1_work_dir = work_dir / "S1_Product"

                # 设置 AnchorPhase 的工作目录
                self.s1_engine.base_work_dir = s1_work_dir

                anchor_result = self.s1_engine.run(
                    molecules={"product": product_smiles}
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

                # 保存 checkpoint 路径（如果存在）
                product_checkpoint = product_min_xyz.with_suffix('.chk')
                if product_checkpoint.exists():
                    result.product_checkpoint = product_checkpoint
                    self.logger.info(f"    ✓ S1 checkpoint 可用: {product_checkpoint}")

                self.logger.info(f"    ✓ 产物锚定完成: {product_min_xyz}")
                self.logger.info(f"    ✓ SP 能量: {e_product_sp:.8f} Hartree")
                self.logger.info(f"    ✓ SP 输出文件: {product_sp_out}")

            except Exception as e:
                result.error_step = "Step1_ProductAnchor_v3"
                result.error_message = str(e)
                self.logger.error(f"Step 1 v3.0 失败: {e}", exc_info=True)
                return result
        else:
            self.logger.warning("⚠️  跳过 Step 1")

        # === Step 2: 逆向扫描 (Retro Scanner) ===
        if 's2' not in skip_steps and result.product_xyz:
            try:
                self.logger.info(">>> Step 2: 逆向生成 TS 初猜与底物...")
                ts_guess_xyz, reactant_xyz, bonds = self.s2_engine.run(
                    product_xyz=result.product_xyz,
                    output_dir=work_dir / "S2_Retro"
                )
                result.ts_guess_xyz = ts_guess_xyz
                result.reactant_xyz = reactant_xyz
                result.forming_bonds = bonds
                self.logger.info(f"    ✓ TS 初猜: {ts_guess_xyz}")
                self.logger.info(f"    ✓ 底物复合物: {reactant_xyz}")
            except Exception as e:
                result.error_step = "Step2_RetroScanner"
                result.error_message = str(e)
                self.logger.error(f"Step 2 失败: {e}", exc_info=True)
                return result
        else:
            self.logger.warning("⚠️  跳过 Step 2")

        # === Step 3: 反应分析 (Transition Analyzer) ===
        if 's3' not in skip_steps and result.ts_guess_xyz:
            try:
                self.logger.info(">>> Step 3: 反应中心全分析 (TS优化 + Reactant/Fragments SP)...")

                # 传递S1的checkpoint以复用轨道
                old_checkpoint = result.product_checkpoint
                if old_checkpoint:
                    self.logger.info(f"  复用S1 checkpoint: {old_checkpoint.name}")

                s3_result = self.s3_engine.run(
                    ts_guess=result.ts_guess_xyz,
                    reactant=result.reactant_xyz,
                    product=result.product_xyz,
                    output_dir=work_dir / "S3_TransitionAnalysis",
                    e_product_l2=result.e_product_l2,
                    forming_bonds=result.forming_bonds,
                    old_checkpoint=old_checkpoint
                )

                result.ts_final_xyz = s3_result.ts_final_xyz
                result.sp_matrix_report = s3_result.sp_report
                self.logger.info(f"    ✓ S3 完成")
                self.logger.info(f"      ΔG‡ = {s3_result.sp_report.get_activation_energy():.3f} kcal/mol")
                self.logger.info(f"      ΔG_rxn = {s3_result.sp_report.get_reaction_energy():.3f} kcal/mol")
            except Exception as e:
                result.error_step = "Step3_TransitionAnalyzer"
                result.error_message = str(e)
                self.logger.error(f"Step 3 失败: {e}", exc_info=True)
                return result
        else:
            self.logger.warning("⚠️  跳过 Step 3")

        # === Step 4: 特征挖掘 (Feature Miner) ===
        if 's4' not in skip_steps and result.ts_final_xyz:
            try:
                self.logger.info(">>> Step 4: 提取物理有机特征...")
                features_csv = self.s4_engine.run(
                    ts_final=result.ts_final_xyz,
                    reactant=result.reactant_xyz,
                    product=result.product_xyz,
                    output_dir=work_dir / "S4_Data",
                    forming_bonds=result.forming_bonds,
                    sp_matrix_report=result.sp_matrix_report
                )
                result.features_csv = features_csv
                self.logger.info(f"    ✓ 特征提取完成: {features_csv}")
            except Exception as e:
                result.error_step = "Step4_FeatureMiner"
                result.error_message = str(e)
                self.logger.error(f"Step 4 失败: {e}", exc_info=True)
                return result
        else:
            self.logger.warning("⚠️  跳过 Step 4")

        # 成功完成
        result.success = True
        self.logger.info("=" * 60)
        self.logger.info(f"✅ 任务完成! 数据已保存至: {work_dir}")
        self.logger.info("=" * 60)

        return result

    def run_batch(
        self,
        smiles_list: list,
        work_dir: Path,
        max_workers: int = 4,
        skip_steps: list = None
    ):
        """
        [Phase 2.4 新增] 批量处理多个 SMILES

        Args:
            smiles_list: SMILES 字符串列表
            work_dir: 工作目录（每个任务将有子目录）
            max_workers: 最大并发数（默认 4）
            skip_steps: 要跳过的步骤列表

        Returns:
            PipelineResult 对象列表
        """
        self.logger.info(f"开始批量处理: {len(smiles_list)} 个分子")
        self.logger.info(f"并发数: {max_workers}")

        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        from concurrent.futures import ProcessPoolExecutor, as_completed
        from tqdm import tqdm

        results = []
        tasks = []

        # 为每个 SMILES 创建任务
        for idx, smiles in enumerate(smiles_list):
            task_work_dir = work_dir / f"task_{idx}_{smiles[:10]}"
            tasks.append((idx, smiles, task_work_dir, skip_steps))

        # 使用 ProcessPoolExecutor 并行执行
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
                        pbar.set_postfix(f"{status} Task {idx}: {smiles[:15]}")
                    except Exception as e:
                        self.logger.error(f"Task {idx} 失败: {e}")
                        pbar.set_postfix(f"❌ Task {idx}: Error")

        self.logger.info(f"批量处理完成: {len(results)}/{len(smiles_list)} 成功")
        return results


# =============================================================================
# 命令行接口
# =============================================================================

def main():
    """命令行主入口"""
    import argparse

    parser = argparse.ArgumentParser(
        description="ReactionProfileHunter v2.1 - 过渡态搜索与特征提取"
    )
    parser.add_argument(
        '--smiles',
        type=str,
        required=True,
        help='产物 SMILES 字符串'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='./rph_output',
        help='输出目录'
    )
    parser.add_argument(
        '--config',
        type=str,
        default=None,
        help='配置文件路径'
    )
    parser.add_argument(
        '--skip-steps',
        type=str,
        default=None,
        help='要跳过的步骤（逗号分隔，如: s1,s2）'
    )
    parser.add_argument(
        '--log-level',
        type=str,
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='日志级别'
    )

    args = parser.parse_args()

    # 解析跳过的步骤
    skip_steps = []
    if args.skip_steps:
        skip_steps = [s.strip().lower() for s in args.skip_steps.split(',')]

    # 初始化并运行
    try:
        hunter = ReactionProfileHunter(
            config_path=Path(args.config) if args.config else None
        )

        result = hunter.run_pipeline(
            product_smiles=args.smiles,
            work_dir=Path(args.output),
            skip_steps=skip_steps
        )

        print("\n" + str(result))

        return 0 if result.success else 1

    except Exception as e:
        logging.error(f"程序异常: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
