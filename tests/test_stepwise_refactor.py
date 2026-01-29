"""
Step-wise 重构测试
==================

为重构后的每一步创建独立测试，使用 reaxys_cleaned.csv 的第一个例子

Author: QCcalc Team
Date: 2026-01-13
"""

import pytest
import logging
from pathlib import Path
import sys

# 添加项目路径
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from rph_core.orchestrator import ReactionProfileHunter

logging.basicConfig(level=logging.INFO)

# 第一个产物的 SMILES (来自 reaxys_cleaned.csv)
FIRST_PRODUCT_SMILES = "O=C1C=C2CCCC(=O)[C@@H]3C[C@H]1O[C@H]23"


class TestStepWiseRefactor:
    """分步重构测试类"""

    @pytest.fixture
    def hunter(self):
        """初始化 RPH"""
        config_path = PROJECT_ROOT / "config" / "defaults.yaml"
        return ReactionProfileHunter(config_path=config_path)

    @pytest.fixture
    def test_smiles(self):
        """测试用的产物 SMILES"""
        return FIRST_PRODUCT_SMILES

    @pytest.fixture
    def tmp_path(self, tmp_path_factory):
        """临时目录 fixture"""
        return tmp_path_factory.mktemp("rph_test_")

    def test_initialization(self, hunter):
        """测试初始化 - 验证新配置结构"""
        assert hunter is not None
        assert hunter.config is not None

        # 验证 theory 配置存在
        assert 'theory' in hunter.config
        assert 'optimization' in hunter.config['theory']
        assert 'single_point' in hunter.config['theory']

        print("✓ 初始化测试通过 - 新配置结构加载正确")

    def test_step1_product_anchor_only(self, hunter, test_smiles, tmp_path):
        """
        测试 Step 1: Product Anchor
        ================================
        目的: 验证 S1 能够正常工作，包括 L2 SP 计算
        """
        result = hunter.run_pipeline(
            product_smiles=test_smiles,
            work_dir=tmp_path / "test_s1_only",
            skip_steps=['s2', 's3', 's4']
        )

        # 验证 S1 输出
        assert result.success, f"S1 失败: {result.error_message}"
        assert result.product_xyz is not None
        assert result.product_xyz.exists(), "产物 XYZ 文件不存在"

        # 验证 L2 能量被计算
        assert result.e_product_l2 is not None, "L2 能量未计算"
        assert result.e_product_l2 < 0, "L2 能量应为负值"

        print("✓ Step 1 (Product Anchor) 测试通过")
        print(f"  Product XYZ: {result.product_xyz}")
        print(f"  L2 Energy: {result.e_product_l2:.8f} Hartree")

    def test_step1_to_step2(self, hunter, test_smiles, tmp_path):
        """
        测试 Step 1 + Step 2: Retro Scanner
        =====================================
        目的: 验证 S1 → S2 的数据传递
        """
        result = hunter.run_pipeline(
            product_smiles=test_smiles,
            work_dir=tmp_path / "test_s1_s2",
            skip_steps=['s3', 's4']
        )

        # 验证 S1 输出
        assert result.success, f"S1+S2 失败: {result.error_message}"
        assert result.product_xyz.exists()
        assert result.e_product_l2 is not None

        # 验证 S2 输出
        assert result.ts_guess_xyz is not None, "TS 初猜未生成"
        assert result.ts_guess_xyz.exists(), "TS 初猜文件不存在"
        assert result.reactant_xyz is not None, "底物复合物未生成"
        assert result.reactant_xyz.exists(), "底物复合物文件不存在"
        assert result.forming_bonds is not None, "形成键信息未识别"

        print("✓ Step 1 + 2 (Retro Scanner) 测试通过")
        print(f"  Product: {result.product_xyz}")
        print(f"  TS Guess: {result.ts_guess_xyz}")
        print(f"  Reactant: {result.reactant_xyz}")
        print(f"  Forming bonds: {result.forming_bonds}")

    def test_step3_transition_analyzer_only(self, hunter, test_smiles, tmp_path):
        """
        测试 Step 3: Transition Analyzer (需要手动准备输入)
        =======================================================
        目的: 验证 S3 的核心功能 - TS 优化 + SP 矩阵构建

        注意: 此测试需要手动准备输入文件，模拟 S1+S2 的输出
        """
        # 手动创建模拟的 S1+S2 输出 (用于独立测试 S3)
        from rph_core.utils.file_io import write_xyz
        import numpy as np

        # 创建模拟产物
        mock_product_dir = tmp_path / "mock_s1"
        mock_product_dir.mkdir(parents=True, exist_ok=True)
        mock_product_xyz = mock_product_dir / "mock_product.xyz"
        write_xyz(mock_product_xyz,
                 np.array([[0.0, 0.0, 0.0]]),
                 ['C'], title="Mock Product")

        # 创建模拟底物
        mock_reactant_xyz = mock_product_dir / "mock_reactant.xyz"
        write_xyz(mock_reactant_xyz,
                 np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
                 ['C', 'C'], title="Mock Reactant")

        # 创建模拟 TS 初猜
        mock_ts_guess_xyz = mock_product_dir / "mock_ts_guess.xyz"
        write_xyz(mock_ts_guess_xyz,
                 np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]),
                 ['C', 'C'], title="Mock TS Guess")

        # 准备模拟 forming_bonds
        mock_forming_bonds = ((0, 1),)

        # 创建模拟的 PipelineResult 用于传递给 S3
        from rph_core.orchestrator import PipelineResult
        from rph_core.steps.step3_5_sp.sp_report import SPMatrixReport

        mock_pipeline_result = PipelineResult(
            success=True,
            product_smiles=test_smiles,
            work_dir=mock_product_dir,
            product_xyz=mock_product_xyz,
            e_product_l2=-150.0,  # 模拟 L2 能量
            product_checkpoint=None,
            ts_guess_xyz=mock_ts_guess_xyz,
            reactant_xyz=mock_reactant_xyz,
            ts_final_xyz=None,
            features_csv=None,
            forming_bonds=mock_forming_bonds,
            sp_matrix_report=None,
            error_step=None,
            error_message=None
        )

        # 直接测试 S3 引擎
        from rph_core.steps.step3_opt.ts_optimizer import TSOptimizer

        s3_engine = TSOptimizer(hunter.config)

        output_dir = tmp_path / "test_s3_direct"
        output_dir.mkdir(parents=True, exist_ok=True)

        print("\n>>> 直接测试 Step 3 引擎...")
        print(f"  输入 TS Guess: {mock_ts_guess_xyz}")
        print(f"  输入 Reactant: {mock_reactant_xyz}")
        print(f"  输入 Product: {mock_product_xyz}")
        print(f"  输入 Product L2 Energy: {-150.0:.8f} Hartree")
        print(f"  Forming bonds: {mock_forming_bonds}")

        try:
            s3_result = s3_engine.run(
                ts_guess=mock_ts_guess_xyz,
                reactant=mock_reactant_xyz,
                product=mock_product_xyz,
                output_dir=output_dir,
                e_product_l2=mock_pipeline_result.e_product_l2,
                forming_bonds=mock_forming_bonds,
                old_checkpoint=None
            )

            # 验证 S3 输出
            assert s3_result.ts_final_xyz.exists(), "TS 最终结构不存在"
            assert s3_result.sp_report is not None, "SP 报告未生成"

            print("✓ Step 3 (Transition Analyzer) 引擎测试通过")
            print(f"  TS Final: {s3_result.ts_final_xyz}")
            print(f"  Method: {s3_result.method_used}")

            # 验证 SP 报告
            sp_report = s3_result.sp_report
            print(f"  SP Report 生成:")
            print(f"    E_TS = {sp_report.e_ts_final:.8f} Hartree")
            print(f"    E_reactant = {sp_report.e_reactant:.8f} Hartree")
            print(f"    E_product = {sp_report.e_product:.8f} Hartree")
            print(f"    ΔG‡ = {sp_report.get_activation_energy():.3f} kcal/mol")
            print(f"    ΔG_rxn = {sp_report.get_reaction_energy():.3f} kcal/mol")

        except Exception as e:
            pytest.fail(f"S3 引擎测试失败: {e}")

    def test_full_pipeline_refactored(self, hunter, test_smiles, tmp_path):
        """
        测试完整流程: S1 → S2 → S3 → S4
        ==========================================
        目的: 验证重构后的完整串行流程

        注意: 此测试需要真实的 QC 软件 (ORCA/Gaussian/XBT)
        """
        result = hunter.run_pipeline(
            product_smiles=test_smiles,
            work_dir=tmp_path / "test_full_pipeline",
            skip_steps=[]  # 运行所有步骤
        )

        # 验证完整输出
        assert result.success, f"完整流程失败: {result.error_message}"

        # S1 输出
        assert result.product_xyz.exists(), "产物 XYZ 不存在"
        assert result.e_product_l2 is not None, "L2 能量未计算"

        # S2 输出
        assert result.ts_guess_xyz.exists(), "TS 初猜不存在"
        assert result.reactant_xyz.exists(), "底物复合物不存在"
        assert result.forming_bonds is not None, "形成键信息缺失"

        # S3 输出 (重构后应包含 SP 报告)
        assert result.ts_final_xyz.exists(), "TS 最终结构不存在"
        assert result.sp_matrix_report is not None, "SP 矩阵报告未生成 (S3 应生成)"

        # S4 输出
        assert result.features_csv.exists(), "特征 CSV 不存在"

        print("✓ 完整流程测试通过")
        print(f"  Work Dir: {result.work_dir}")
        print(f"  Product: {result.product_xyz}")
        print(f"  TS Final: {result.ts_final_xyz}")
        print(f"  Features: {result.features_csv}")

        # 打印 SP 报告摘要
        if result.sp_matrix_report:
            print(f"\n  SP Report Summary:")
            print(f"    E_product = {result.sp_matrix_report.e_product:.8f} Ha")
            print(f"    E_reactant = {result.sp_matrix_report.e_reactant:.8f} Ha")
            print(f"    E_TS = {result.sp_matrix_report.e_ts_final:.8f} Ha")
            print(f"    ΔG‡ = {result.sp_matrix_report.get_activation_energy():.3f} kcal/mol")
            print(f"    ΔG_rxn = {result.sp_matrix_report.get_reaction_energy():.3f} kcal/mol")


# 运行标记
# ============
# 运行特定测试:
#
#   pytest tests/test_stepwise_refactor.py -v -k "test_step1_product_anchor_only"
#   pytest tests/test_stepwise_refactor.py -v -k "test_step1_to_step2"
#   pytest tests/test_stepwise_refactor.py -v -k "test_step3_transition_analyzer_only"
#   pytest tests/test_stepwise_refactor.py -v -k "test_full_pipeline_refactored"
