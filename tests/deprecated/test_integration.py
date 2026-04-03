"""
示例测试：ReactionProfileHunter v2.1 集成测试
==============================================

注意: 这是一个示例测试，展示如何测试完整的串行工作流
实际使用时需要：
1. 安装 XTB 和 Gaussian
2. 准备测试分子
3. 根据计算资源调整超时时间

Author: QCcalc Team
Date: 2026-01-09
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


class TestReactionProfileHunter:
    """ReactionProfileHunter 集成测试"""

    @pytest.fixture
    def hunter(self):
        """初始化 RPH"""
        config_path = PROJECT_ROOT / "config" / "defaults.yaml"
        return ReactionProfileHunter(config_path=config_path)

    @pytest.fixture
    def test_smiles(self):
        """测试用的 [5+2] 产物 SMILES"""
        # 氧杂桥环产物
        return "O=C1CC2CCCCC2O1"

    def test_initialization(self, hunter):
        """测试初始化"""
        assert hunter is not None
        assert hunter.config is not None
        print("✓ 初始化测试通过")

    def test_step1_only(self, hunter, test_smiles, tmp_path):
        """测试仅运行 Step 1"""
        result = hunter.run_pipeline(
            product_smiles=test_smiles,
            work_dir=tmp_path / "test_s1",
            skip_steps=['s2', 's3', 's4']
        )

        assert result.success
        assert result.product_xyz.exists()
        print(f"✓ Step 1 测试通过: {result.product_xyz}")

    def test_step1_step2_only(self, hunter, test_smiles, tmp_path):
        """测试运行 Step 1 + Step 2"""
        result = hunter.run_pipeline(
            product_smiles=test_smiles,
            work_dir=tmp_path / "test_s1_s2",
            skip_steps=['s3', 's4']
        )

        assert result.success
        assert result.product_xyz.exists()
        assert result.ts_guess_xyz.exists()
        assert result.reactant_xyz.exists()
        print(f"✓ Step 1+2 测试通过")
        print(f"  Product: {result.product_xyz}")
        print(f"  TS Guess: {result.ts_guess_xyz}")
        print(f"  Reactant: {result.reactant_xyz}")

    def test_full_pipeline_mock(self, hunter, test_smiles, tmp_path):
        """测试完整流程（使用模拟 Step 3/4）"""
        result = hunter.run_pipeline(
            product_smiles=test_smiles,
            work_dir=tmp_path / "test_full"
        )

        # 注意: Step 3 和 4 使用了简化实现，不会实际运行 Gaussian
        assert result.success
        assert result.product_xyz.exists()
        assert result.ts_final_xyz.exists()
        print(f"✓ 完整流程测试通过")
        print(f"  所有文件已生成")


if __name__ == "__main__":
    # 快速测试
    print("运行 ReactionProfileHunter 快速测试...")
    print("=" * 60)

    test = TestReactionProfileHunter()

    # 初始化
    print("\n[1/4] 测试初始化...")
    hunter = test.hunter()
    test.test_initialization(hunter)

    # Step 1
    print("\n[2/4] 测试 Step 1...")
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as tmpdir:
        test.test_step1_only(hunter, test.test_smiles(), Path(tmpdir))

    # Step 1+2
    print("\n[3/4] 测试 Step 1+2...")
    with TemporaryDirectory() as tmpdir:
        test.test_step1_step2_only(hunter, test.test_smiles(), Path(tmpdir))

    # 完整流程
    print("\n[4/4] 测试完整流程...")
    with TemporaryDirectory() as tmpdir:
        test.test_full_pipeline_mock(hunter, test.test_smiles(), Path(tmpdir))

    print("\n" + "=" * 60)
    print("✅ 所有测试通过！")
    print("\n注意:")
    print("- Step 3 和 4 使用简化实现，未实际调用 Gaussian")
    print("- 完整实现需要安装 XTB 和 Gaussian")
    print("- 真实计算需要较长时间（小时级）")
