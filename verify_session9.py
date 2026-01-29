"""
Session #9 验证脚本 - 测试 SPMatrixBuilder 核心逻辑
"""
import sys
import tempfile
from pathlib import Path
from unittest.mock import Mock, MagicMock, patch

# 添加项目路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

# Mock 所有 rdkit 相关导入
sys.modules['rdkit'] = MagicMock()
sys.modules['rdkit.Chem'] = MagicMock()

# 直接导入，避免触发 steps/__init__.py 的依赖
import importlib.util

# 导入 SPMatrixReport
spec_report = importlib.util.spec_from_file_location(
    "sp_report",
    project_root / "rph_core" / "steps" / "step3_5_sp" / "sp_report.py"
)
sp_report = importlib.util.module_from_spec(spec_report)
spec_report.loader.exec_module(sp_report)
SPMatrixReport = sp_report.SPMatrixReport

# 导入 SPMatrixBuilder - 需要提前 mock SPMatrixReport
spec_matrix = importlib.util.spec_from_file_location(
    "sp_matrix",
    project_root / "rph_core" / "steps" / "step3_5_sp" / "sp_matrix.py"
)
sp_matrix = importlib.util.module_from_spec(spec_matrix)

# 设置已导入的 SPMatrixReport
sp_matrix.SPMatrixReport = SPMatrixReport

# Mock ORCAInterface 导入
mock_qc_result = MagicMock()
mock_qc_result.QCResult = MagicMock
sys.modules['rph_core.utils.orca_interface'] = MagicMock()
sys.modules['rph_core.utils.orca_interface'].ORCAInterface = Mock
sys.modules['rph_core.utils.orca_interface'].QCResult = mock_qc_result.QCResult

# Mock GaussianInterface 导入 (可能被 fallback 使用)
sys.modules['rph_core.utils.qc_interface'] = MagicMock()

spec_matrix.loader.exec_module(sp_matrix)
SPMatrixBuilder = sp_matrix.SPMatrixBuilder


def create_mock_xyz(directory: Path, name: str) -> Path:
    """创建模拟 XYZ 文件"""
    xyz_file = directory / f"{name}.xyz"
    xyz_file.write_text("""3
comment
C    0.0  0.0  0.0
H    1.0  0.0  0.0
H   -1.0  0.0  0.0
""")
    return xyz_file


def test_sp_matrix_builder_init():
    """测试 SPMatrixBuilder 初始化"""
    print("=" * 60)
    print("测试: SPMatrixBuilder.__init__()")
    print("=" * 60)

    config = {
        'method': 'M062X',
        'basis': 'def2-TZVPP',
        'aux_basis': 'def2/J',
        'nprocs': 16,
        'maxcore': 4000,
        'solvent': 'acetone'
    }

    # Mock ORCAInterface
    with patch('rph_core.utils.orca_interface.ORCAInterface') as mock_orca_class:
        mock_orca = Mock()
        mock_orca.method = 'M062X'
        mock_orca.basis = 'def2-TZVPP'
        mock_orca.solvent = 'acetone'
        mock_orca_class.return_value = mock_orca

        builder = SPMatrixBuilder(config)

        tests = [
            (builder.config == config, "配置正确保存"),
            (builder.orca.method == 'M062X', "ORCA method 设置正确"),
            (builder.orca.basis == 'def2-TZVPP', "ORCA basis 设置正确"),
            (builder.orca.solvent == 'acetone', "ORCA solvent 设置正确"),
        ]

        all_passed = True
        for passed, description in tests:
            status = "✓ PASS" if passed else "✗ FAIL"
            print(f"  {status}: {description}")
            all_passed = all_passed and passed

    print("\n" + "=" * 60)
    if all_passed:
        print("✓ 测试通过!")
    else:
        print("✗ 测试失败")
    print("=" * 60)

    return all_passed


def test_sp_matrix_builder_run_with_energy_reuse():
    """测试 SPMatrixBuilder.run() 复用 S1 能量"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixBuilder.run() - 复用 S1 能量")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # 创建模拟 XYZ 文件
        product_xyz = create_mock_xyz(tmp_path, "product")
        reactant_xyz = create_mock_xyz(tmp_path, "reactant")
        ts_xyz = create_mock_xyz(tmp_path, "ts")
        frag_a_xyz = create_mock_xyz(tmp_path, "frag_a")
        frag_b_xyz = create_mock_xyz(tmp_path, "frag_b")
        output_dir = tmp_path / "S35_SP"

        # Mock ORCAInterface
        with patch('rph_core.utils.orca_interface.ORCAInterface') as mock_orca_class:
            # Mock ORCA instance
            mock_orca = Mock()
            mock_orca.method = 'M062X'
            mock_orca.basis = 'def2-TZVPP'
            mock_orca.solvent = 'acetone'

            # Mock single_point 返回结果
            mock_result = Mock()
            mock_result.converged = True
            mock_result.energy = -456.78901234
            mock_orca.single_point.return_value = mock_result

            mock_orca_class.return_value = mock_orca

            # 创建 SPMatrixBuilder
            config = {'method': 'M062X', 'basis': 'def2-TZVPP'}
            builder = SPMatrixBuilder(config)

            # 运行 (提供 S1 已计算的产物能量)
            e_product_l2_s1 = -123.45678901
            report = builder.run(
                product_xyz=product_xyz,
                reactant_xyz=reactant_xyz,
                ts_final_xyz=ts_xyz,
                frag_a_xyz=frag_a_xyz,
                frag_b_xyz=frag_b_xyz,
                output_dir=output_dir,
                e_product_l2=e_product_l2_s1
            )

            # 验证结果
            tests = [
                (report.e_product == e_product_l2_s1, f"产物能量复用 S1: {report.e_product}"),
                (report.e_reactant == -456.78901234, f"底物能量正确: {report.e_reactant}"),
                (report.e_ts_final == -456.78901234, f"TS 能量正确: {report.e_ts_final}"),
                (report.e_frag_a_ts == -456.78901234, f"片段A能量正确: {report.e_frag_a_ts}"),
                (report.e_frag_b_ts == -456.78901234, f"片段B能量正确: {report.e_frag_b_ts}"),
                (mock_orca.single_point.call_count == 4, "调用 4 次 single_point (跳过产物)"),
            ]

            all_passed = True
            for passed, description in tests:
                status = "✓ PASS" if passed else "✗ FAIL"
                print(f"  {status}: {description}")
                all_passed = all_passed and passed

    print("\n" + "=" * 60)
    if all_passed:
        print("✓ 测试通过!")
    else:
        print("✗ 测试失败")
    print("=" * 60)

    return all_passed


def test_sp_matrix_builder_run_full_calculation():
    """测试 SPMatrixBuilder.run() 完整计算"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixBuilder.run() - 完整计算")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # 创建模拟 XYZ 文件
        product_xyz = create_mock_xyz(tmp_path, "product")
        reactant_xyz = create_mock_xyz(tmp_path, "reactant")
        ts_xyz = create_mock_xyz(tmp_path, "ts")
        frag_a_xyz = create_mock_xyz(tmp_path, "frag_a")
        frag_b_xyz = create_mock_xyz(tmp_path, "frag_b")
        output_dir = tmp_path / "S35_SP"

        # Mock ORCAInterface
        with patch('rph_core.utils.orca_interface.ORCAInterface') as mock_orca_class:
            # Mock ORCA instance
            mock_orca = Mock()
            mock_orca.method = 'M062X'
            mock_orca.basis = 'def2-TZVPP'
            mock_orca.solvent = 'acetone'

            # Mock single_point 返回不同能量
            def mock_single_point(xyz, out_dir):
                result = Mock()
                result.converged = True
                # 根据文件名返回不同能量
                if 'product' in str(xyz):
                    result.energy = -123.45678901
                elif 'reactant' in str(xyz):
                    result.energy = -456.78901234
                elif 'ts' in str(xyz):
                    result.energy = -345.67890123
                elif 'frag_a' in str(xyz):
                    result.energy = -234.56789012
                elif 'frag_b' in str(xyz):
                    result.energy = -123.45678901
                else:
                    result.energy = -100.0
                return result

            mock_orca.single_point = Mock(side_effect=mock_single_point)
            mock_orca_class.return_value = mock_orca

            # 创建 SPMatrixBuilder
            config = {'method': 'M062X', 'basis': 'def2-TZVPP'}
            builder = SPMatrixBuilder(config)

            # 运行 (不提供 S1 能量)
            report = builder.run(
                product_xyz=product_xyz,
                reactant_xyz=reactant_xyz,
                ts_final_xyz=ts_xyz,
                frag_a_xyz=frag_a_xyz,
                frag_b_xyz=frag_b_xyz,
                output_dir=output_dir,
                e_product_l2=None  # 不复用
            )

            # 验证结果
            tests = [
                (report.e_product == -123.45678901, f"产物能量正确: {report.e_product}"),
                (report.e_reactant == -456.78901234, f"底物能量正确: {report.e_reactant}"),
                (report.e_ts_final == -345.67890123, f"TS 能量正确: {report.e_ts_final}"),
                (report.e_frag_a_ts == -234.56789012, f"片段A能量正确: {report.e_frag_a_ts}"),
                (report.e_frag_b_ts == -123.45678901, f"片段B能量正确: {report.e_frag_b_ts}"),
                (mock_orca.single_point.call_count == 5, "调用 5 次 single_point"),
            ]

            all_passed = True
            for passed, description in tests:
                status = "✓ PASS" if passed else "✗ FAIL"
                print(f"  {status}: {description}")
                all_passed = all_passed and passed

    print("\n" + "=" * 60)
    if all_passed:
        print("✓ 测试通过!")
    else:
        print("✗ 测试失败")
    print("=" * 60)

    return all_passed


def test_sp_matrix_builder_report_saved():
    """测试报告文件保存"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixBuilder 报告保存")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # 创建模拟 XYZ 文件
        product_xyz = create_mock_xyz(tmp_path, "product")
        reactant_xyz = create_mock_xyz(tmp_path, "reactant")
        ts_xyz = create_mock_xyz(tmp_path, "ts")
        frag_a_xyz = create_mock_xyz(tmp_path, "frag_a")
        frag_b_xyz = create_mock_xyz(tmp_path, "frag_b")
        output_dir = tmp_path / "S35_SP"

        # Mock ORCAInterface
        with patch('rph_core.utils.orca_interface.ORCAInterface') as mock_orca_class:
            mock_orca = Mock()
            mock_orca.method = 'PWPB95'
            mock_orca.basis = 'def2-TZVPP'
            mock_orca.solvent = 'water'

            mock_result = Mock()
            mock_result.converged = True
            mock_result.energy = -200.0
            mock_orca.single_point.return_value = mock_result

            mock_orca_class.return_value = mock_orca

            # 创建 SPMatrixBuilder
            config = {
                'method': 'PWPB95',
                'basis': 'def2-TZVPP',
                'solvent': 'water'
            }
            builder = SPMatrixBuilder(config)

            # 运行
            report = builder.run(
                product_xyz=product_xyz,
                reactant_xyz=reactant_xyz,
                ts_final_xyz=ts_xyz,
                frag_a_xyz=frag_a_xyz,
                frag_b_xyz=frag_b_xyz,
                output_dir=output_dir,
                e_product_l2=-100.0  # 复用产物能量
            )

            # 验证报告文件
            report_file = output_dir / "sp_matrix_report.json"

            tests = [
                (report_file.exists(), "报告文件已创建"),
                (report.method == "PWPB95/def2-TZVPP", f"方法正确: {report.method}"),
                (report.solvent == "water", f"溶剂正确: {report.solvent}"),
                (report.e_product == -100.0, f"产物能量正确: {report.e_product}"),
            ]

            all_passed = True
            for passed, description in tests:
                status = "✓ PASS" if passed else "✗ FAIL"
                print(f"  {status}: {description}")
                all_passed = all_passed and passed

            # 验证文件内容
            if report_file.exists():
                import json
                content = json.loads(report_file.read_text())
                tests.append(
                    (content["e_product_l2"] == -100.0, "JSON 文件内容正确")
                )

    print("\n" + "=" * 60)
    if all_passed:
        print("✓ 测试通过!")
    else:
        print("✗ 测试失败")
    print("=" * 60)

    return all_passed


if __name__ == "__main__":
    print("\n")
    print("*" * 60)
    print("*" + " " * 58 + "*")
    print("*" + " " * 8 + "Session #9: SPMatrixBuilder 测试套件" + " " * 18 + "*")
    print("*" + " " * 58 + "*")
    print("*" * 60)

    results = []
    results.append(("初始化", test_sp_matrix_builder_init()))
    results.append(("复用 S1 能量", test_sp_matrix_builder_run_with_energy_reuse()))
    results.append(("完整计算", test_sp_matrix_builder_run_full_calculation()))
    results.append(("报告保存", test_sp_matrix_builder_report_saved()))

    print("\n\n")
    print("*" * 60)
    print("*" + " " * 20 + "测试总结" + " " * 28 + "*")
    print("*" * 60)

    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {name}")

    all_passed = all(r[1] for r in results)

    print("*" * 60)
    if all_passed:
        print("✓✓✓ 所有测试通过! Session #9 完成 ✓✓✓")
    else:
        print("✗✗✗ 部分测试失败 ✗✗✗")
    print("*" * 60)
    print()

    sys.exit(0 if all_passed else 1)
