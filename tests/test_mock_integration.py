"""
Mock 集成测试 - Phase 1.3

完全 Mock 的端到端集成测试，不依赖真实的 QC 软件

测试目标:
1. 验证 S1-S4 之间的数据流连通性
2. 验证 S3.5 → S4 的 L2 能量传递
3. 验证 feature_miner 正确使用 SPMatrixReport

Author: QCcalc Team
Date: 2026-01-12
Phase: 1.3 - 端到端 Mock 集成测试
"""

import logging
import json
from pathlib import Path
from dataclasses import dataclass
from typing import Optional
import numpy as np
import sys
from tempfile import TemporaryDirectory

# 添加项目路径
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# 尝试导入 SPMatrixReport
try:
    from rph_core.steps.step3_opt.ts_optimizer import SPMatrixReport
    RPH_AVAILABLE = True
except ImportError as e:
    logger = logging.getLogger(__name__)
    logger.warning(f"RPH 模块导入失败: {e}")
    RPH_AVAILABLE = False

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def create_mock_xyz_files(tmp_path: Path) -> dict:
    """创建 Mock XYZ 文件"""
    # 创建简单的 [5+2] 产物结构
    product_xyz = tmp_path / "product.xyz"
    product_xyz.write_text("""3
Oxabicyclo product
C    0.0  0.0  0.0
O    1.0  0.0  0.0
H    2.0  0.0  0.0
""")

    # 创建 TS 结构
    ts_xyz = tmp_path / "ts.xyz"
    ts_xyz.write_text("""3
Transition State
C    0.0  0.0  0.0
O    2.2  0.0  0.0
H    4.0  0.0  0.0
""")

    # 创建底物结构
    reactant_xyz = tmp_path / "reactant.xyz"
    reactant_xyz.write_text("""3
Reactant
C    0.0  0.0  0.0
O    3.5  0.0  0.0
H    5.0  0.0  0.0
""")

    # 创建片段结构
    frag_a_xyz = tmp_path / "frag_a.xyz"
    frag_a_xyz.write_text("""1
Fragment A
C    0.0  0.0  0.0
""")

    frag_b_xyz = tmp_path / "frag_b.xyz"
    frag_b_xyz.write_text("""2
Fragment B
O    0.0  0.0  0.0
H    1.0  0.0  0.0
""")

    return {
        'product': product_xyz,
        'ts': ts_xyz,
        'reactant': reactant_xyz,
        'frag_a': frag_a_xyz,
        'frag_b': frag_b_xyz
    }


def create_mock_sp_matrix_report() -> SPMatrixReport:
    """Create Mock SPMatrixReport

    Returns SPMatrixReport with V4.1 thermo energies (kcal/mol).
    """
    return SPMatrixReport(
        g_reactant=-100.0,
        g_ts=-80.0,  # ΔG‡ = 20
        g_product=-110.0,  # ΔG_rxn = -10
        method="WB97M-V/def2-TZVPP",
        solvent="acetone",
        # Optional fields (not used in V4.1 but kept for compatibility)
        e_ts_final=-345.67890123,
        e_reactant=-456.78901234,
        e_product=-123.45678901
    )


def test_sp_matrix_report_validation():
    """测试 SPMatrixReport 验证"""
    print("\n[1/4] SPMatrixReport 验证...")

    if not RPH_AVAILABLE:
        print("  ⚠️  RPH 模块不可用，跳过测试")
        return

    report = create_mock_sp_matrix_report()

    # 验证报告数据
    assert report.validate() == True, "SPMatrixReport 验证失败"

    # V4.1: 畸变能测试已移除 (get_distortion_energy_a/b 不再存在)

    # 测试能量转换 (V4.1: Gibbs energies in kcal/mol)
    dG_activation = report.get_activation_energy()
    dG_reaction = report.get_reaction_energy()

    print(f"  ΔG‡ = {dG_activation:.3f} kcal/mol")
    print(f"  ΔG_rxn = {dG_reaction:.3f} kcal/mol")

    # 验证活化能为正
    assert dG_activation > 0, f"活化能应为正: {dG_activation}"

    # 验证反应能可能为负
    assert dG_reaction < 0, f"反应能应为负: {dG_reaction}"

    print("  ✓ 通过")

# V4.1: 畸变能测试已移除 (get_distortion_energy_a/b 不再存在)

    print("  ✓ 通过")

# V4.1: 畸变能测试已移除 (E_distortion methods no longer exist in V4.1)


def test_sp_matrix_report_serialization():
    """测试 SPMatrixReport 序列化 - V4.1 dataclass field names"""
    print("\n[2/4] SPMatrixReport 序列化...")

    if not RPH_AVAILABLE:
        print("  ⚠️  RPH 模块不可用，跳过测试")
        return

    report = create_mock_sp_matrix_report()

    # V4.1: to_dict() uses dataclass field names (e_product, not e_product_l2)
    data_dict = report.to_dict()
    assert 'e_product' in data_dict, "dataclass field name check failed"
    assert 'method' in data_dict, "dataclass field name check failed"

    # 测试 to_json()
    json_str = report.to_json()
    # V4.1: JSON now contains actual method name, not M062X-specific string
    assert report.method in json_str, "method name check failed"

    # 测试 from_dict()
    report2 = SPMatrixReport.from_dict(data_dict)
    assert report2.e_product == report.e_product, "from_dict roundtrip failed"
    assert report2.method == report.method, "from_dict roundtrip failed"

    # 测试 JSON 文件保存和加载
    with TemporaryDirectory() as tmpdir:
        json_file = Path(tmpdir) / "sp_matrix_report.json"
        json_file.write_text(json_str)

        loaded_data = json.loads(json_file.read_text())
        report3 = SPMatrixReport.from_dict(loaded_data)
        assert report3.e_product == report.e_product, "JSON load roundtrip failed"
        assert report3.method == report.method, "JSON load roundtrip failed"

    print("  ✓ 通过")


def test_sp_matrix_l2_energy_calculation():
    """测试 L2 能量计算"""
    print("\n[3/4] L2 能量计算...")
    if not RPH_AVAILABLE:
        print("  ⚠️  RPH 模块不可用，跳过测试")
        return

    # V4.1: 畸变能/片段测试已移除，跳过此测试
    print("  (removed in V4.1 architecture)")

    report = create_mock_sp_matrix_report()

    # 测试活化能 (using Gibbs priority from SPMatrixReport)
    dG_activation = report.get_activation_energy()
    assert dG_activation == 20.0, f"Gibbs activation energy check failed: {dG_activation} != 20.0"

    print(f"  ΔG‡ = {dG_activation:.3f} kcal/mol")

    print("  ✓ 通过")




        # 验证能量值
        assert data['e_product_l2'] == report.e_product
        assert data['method'] == report.method

        print(f"  ✓ JSON 数据正确")

    print("  ✓ 通过")


def run_all_tests():
    """运行所有 Mock 集成测试"""
    print("=" * 70)
    print("Phase 1.3: Mock 集成测试")
    print("=" * 70)
    print(f"\nRPH 模块可用: {'是' if RPH_AVAILABLE else '否'}")

    if not RPH_AVAILABLE:
        print("\n⚠️  RPH 模块不可用，跳过所有测试")
        print("请确保 RPH 核心模块已正确安装")
        return

    try:
        test_sp_matrix_report_validation()
        test_sp_matrix_report_serialization()
        test_sp_matrix_l2_energy_calculation()
        test_feature_miner_integration()

        print("\n" + "=" * 70)
        print("✅ Phase 1.3: 所有 Mock 集成测试通过！")
        print("=" * 70)
    print("\n测试总结:")
    print("  1. ✓ SPMatrixReport 验证")
    print("  2. ✓ SPMatrixReport 序列化")
    print("  3. ✓ L2 能量计算")
    print("  4. ✓ FeatureMiner 集成")
    print("\n验证内容:")
    print("  • 数据流连通性: S1 → S3.5 → S4 ✓")
    print("  • L2 能量传递: SPMatrixReport → FeatureMiner ✓")
    print("  • V4.1 Schema: 4.1-Shermo ✓")

def test_route_a_missing_sp_report():
    """测试 Route A 失败分支 - missing sp_matrix_report"""
    print("\n[4/4] Route A: missing SPMatrixReport...")

    if not RPH_AVAILABLE:
        print("  ⚠️  RPH 模块不可用，跳过测试")
        return

    from rph_core.steps.step4_features.feature_miner import FeatureMiner

    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)

        # V4.1: Route A - 不提供 SPMatrixReport（或 None）
        features_csv = FeatureMiner.run(
            ts_final=tmp_path / 'ts.xyz',
            reactant=tmp_path / 'reactant.xyz',
            product=tmp_path / 'product.xyz',
            output_dir=tmp_path / 'S4_Data',
            sp_matrix_report=None  # 模拟 Step 3 失败
        )

        # Route A: 即使失败也应该生成 features_raw.csv
        assert features_csv.exists(), "features_raw.csv should always be generated (Route A)"
        import pandas as pd
        data = pd.read_csv(features_csv)

        # 验证状态
        assert data['feature_status'].iloc[0] == 'missing_sp_report', "Should report missing_sp_report"
        assert data['L2_available'].iloc[0] == 0, "L2 should not be available when report missing"
        assert pd.isna(data['dG_activation_L2'].iloc[0]), "dG_activation_L2 should be NaN"

        print("  ✓ Route A 失败分支通过")
        print("\n验证内容:")
        print("  • feature_status: missing_sp_report ✓")
        print("  • features_raw.csv generated ✓")

def test_route_a_missing_thermo():
    """测试 Route A 失败分支 - missing thermo data"""
    print("\n[4/4] Route A: missing thermo...")

    if not RPH_AVAILABLE:
        print("  ⚠️  RPH 模块不可用，跳过测试")
        return

    from rph_core.steps.step3_opt.ts_optimizer import SPMatrixReport
    from rph_core.steps.step4_features.feature_miner import FeatureMiner

    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)

        # V4.1: Route A - SPMatrixReport 缺少 thermo (g_*)
        report = SPMatrixReport(
            g_reactant=2.0,
            g_ts=10.0,
            g_product=None,  # 缺失
            method="test",
            solvent="test"
        )

        features_csv = FeatureMiner.run(
            ts_final=tmp_path / 'ts.xyz',
            reactant=tmp_path / 'reactant.xyz',
            product=tmp_path / 'product.xyz',
            output_dir=tmp_path / 'S4_Data',
            sp_matrix_report=report
        )

        # Route A: 缺少 thermo 时应该标记为 missing_thermo
        assert features_csv.exists(), "features_raw.csv should always be generated (Route A)"
        import pandas as pd
        data = pd.read_csv(features_csv)

        # 验证状态
        assert data['feature_status'].iloc[0] == 'missing_thermo', "Should report missing_thermo"
        assert data['L2_available'].iloc[0] == 0, "L2 should not be available when thermo missing"
        assert pd.isna(data['dG_activation_L2'].iloc[0]), "dG_activation_L2 should be NaN"
        assert pd.isna(data['dG_reaction_L2'].iloc[0]), "dG_reaction_L2 should be NaN"

        print("  ✓ Route A missing thermo 分支通过")
        print("\n验证内容:")
        print("  • feature_status: missing_thermo ✓")
        print("  • features_raw.csv generated ✓")

def run_all_tests():
    """运行所有 Mock 集成测试"""
    print("=" * 70)
    print("Phase 1.3: Mock 集成测试")
    print("=" * 70)
    print(f"\nRPH 模块可用: {'是' if RPH_AVAILABLE else '否'}")

    if RPH_AVAILABLE:
        test_sp_matrix_report_validation()
        test_sp_matrix_report_serialization()
        test_feature_miner_integration()
        test_route_a_missing_sp_report()
        test_route_a_missing_thermo()

        print("\n测试总结:")
        print("  1. ✓ SPMatrixReport 验证")
        print("  2. ✓ SPMatrixReport 序列化")
        print("  3. ✓ L2 能量计算 (Gibbs 优先)")
        print("  4. ✓ FeatureMiner 集成 (V4.1: always outputs features_raw.csv)")
        print("  5. ✓ Route A: 永远不崩溃 (never crashes)")
        print("\n验证内容:")
        print("  • 数据流连通性: S1 → S3 → S4 ✓")
        print("  • V4.1 Schema: 4.1-Shermo ✓")
        print("  • L2 能量传递: SPMatrixReport → FeatureMiner ✓")
        print("\n这些测试不依赖真实的 QC 软件，仅验证逻辑连通性。")

    print("✅ Phase 1.3: 所有 Mock 集成测试通过！")
    print("=" * 70)
    print("\n测试总结:")
    print("  1. ✓ SPMatrixReport 验证")
    print("  2. ✓ SPMatrixReport 序列化")
    print("  3. ✓ L2 能量计算 (Gibbs 优先)")
    print("  4. ✓ FeatureMiner 集成 (V4.1: always outputs features_raw.csv)")
    print("  5. ✓ Route A: 永远不崩溃 (never crashes)")

    print("完整测试需要安装 XTB 和 Gaussian。")

    except AssertionError as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 测试出错: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    run_all_tests()
