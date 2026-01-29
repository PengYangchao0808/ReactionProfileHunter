"""
Session #9 简化验证脚本 - 测试 SPMatrixBuilder 结构
"""
import sys
from pathlib import Path

# 添加项目路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

def test_sp_matrix_file_exists():
    """测试 SPMatrixBuilder 文件存在"""
    print("=" * 60)
    print("测试: SPMatrixBuilder 文件结构")
    print("=" * 60)

    sp_matrix_file = project_root / "rph_core" / "steps" / "step3_5_sp" / "sp_matrix.py"

    tests = [
        (sp_matrix_file.exists(), f"sp_matrix.py 存在: {sp_matrix_file}"),
    ]

    # 读取文件内容
    if sp_matrix_file.exists():
        content = sp_matrix_file.read_text()

        tests.append((
            "class SPMatrixBuilder" in content,
            "包含 SPMatrixBuilder 类定义"
        ))
        tests.append((
            "def __init__(self, config: dict)" in content,
            "包含 __init__ 方法"
        ))
        tests.append((
            "def run(" in content,
            "包含 run 方法"
        ))
        tests.append((
            "def _run_sp(" in content,
            "包含 _run_sp 方法"
        ))
        tests.append((
            "def _save_report(" in content,
            "包含 _save_report 方法"
        ))
        tests.append((
            "ORCAInterface" in content,
            "导入 ORCAInterface"
        ))
        tests.append((
            "SPMatrixReport" in content,
            "导入 SPMatrixReport"
        ))

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


def test_sp_matrix_init_structure():
    """测试 __init__ 方法结构"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixBuilder.__init__() 结构")
    print("=" * 60)

    sp_matrix_file = project_root / "rph_core" / "steps" / "step3_5_sp" / "sp_matrix.py"
    content = sp_matrix_file.read_text()

    # 检查 __init__ 方法的实现
    tests = [
        ("self.config = config" in content, "保存 config"),
        ("self.orca = ORCAInterface(" in content, "初始化 ORCAInterface"),
        ("self.logger" in content, "初始化 logger"),
        ("config.get('method'" in content, "支持 method 配置"),
        ("config.get('basis'" in content, "支持 basis 配置"),
        ("config.get('solvent'" in content, "支持 solvent 配置"),
        ("config.get('nprocs'" in content, "支持 nprocs 配置"),
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


def test_sp_matrix_run_structure():
    """测试 run 方法结构"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixBuilder.run() 结构")
    print("=" * 60)

    sp_matrix_file = project_root / "rph_core" / "steps" / "step3_5_sp" / "sp_matrix.py"
    content = sp_matrix_file.read_text()

    # 检查 run 方法的实现
    tests = [
        ("def run(" in content and "product_xyz" in content, "接受 product_xyz 参数"),
        ("def run(" in content and "reactant_xyz" in content, "接受 reactant_xyz 参数"),
        ("def run(" in content and "ts_final_xyz" in content, "接受 ts_final_xyz 参数"),
        ("def run(" in content and "frag_a_xyz" in content, "接受 frag_a_xyz 参数"),
        ("def run(" in content and "frag_b_xyz" in content, "接受 frag_b_xyz 参数"),
        ("e_product_l2: Optional[float]" in content, "支持 e_product_l2 参数"),
        ("-> SPMatrixReport" in content, "返回 SPMatrixReport"),
        ("if e_product_l2 is not None:" in content, "支持复用 S1 能量"),
        ("self._run_sp(" in content, "调用 _run_sp 方法"),
        ("SPMatrixReport(" in content, "创建 SPMatrixReport"),
        ("self._save_report(" in content, "保存报告"),
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


def test_sp_matrix_run_sp_structure():
    """测试 _run_sp 方法结构"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixBuilder._run_sp() 结构")
    print("=" * 60)

    sp_matrix_file = project_root / "rph_core" / "steps" / "step3_5_sp" / "sp_matrix.py"
    content = sp_matrix_file.read_text()

    # 检查 _run_sp 方法的实现
    tests = [
        ("def _run_sp(self, xyz_file: Path" in content, "接受 xyz_file 参数"),
        ("def _run_sp(" in content and "output_dir: Path" in content, "接受 output_dir 参数"),
        ("self.orca.single_point(" in content, "调用 ORCAInterface.single_point"),
        ("result.converged" in content, "检查收敛状态"),
        ("result.energy" in content, "返回能量"),
        ("fallback_to_gaussian" in content, "支持 Gaussian fallback"),
        ("RuntimeError" in content, "抛出 RuntimeError"),
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


def test_sp_matrix_save_report_structure():
    """测试 _save_report 方法结构"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixBuilder._save_report() 结构")
    print("=" * 60)

    sp_matrix_file = project_root / "rph_core" / "steps" / "step3_5_sp" / "sp_matrix.py"
    content = sp_matrix_file.read_text()

    # 检查 _save_report 方法的实现
    tests = [
        ("def _save_report(self, report: SPMatrixReport" in content, "接受 report 参数"),
        ("output_path.write_text(" in content, "写入文件"),
        ("report.to_json(" in content, "调用 to_json 方法"),
        ("self.logger.info" in content, "记录日志"),
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


def test_sp_matrix_docstrings():
    """测试文档字符串"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixBuilder 文档字符串")
    print("=" * 60)

    sp_matrix_file = project_root / "rph_core" / "steps" / "step3_5_sp" / "sp_matrix.py"
    content = sp_matrix_file.read_text()

    # 检查文档字符串
    tests = [
        ('"""' in content, "包含文档字符串"),
        ("SP 矩阵构建器" in content or "SP Matrix Builder" in content, "类/模块描述"),
        ("Args:" in content, "包含 Args 说明"),
        ("Returns:" in content, "包含 Returns 说明"),
        ("Raises:" in content, "包含 Raises 说明"),
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


def test_sp_matrix_code_quality():
    """测试代码质量"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixBuilder 代码质量")
    print("=" * 60)

    sp_matrix_file = project_root / "rph_core" / "steps" / "step3_5_sp" / "sp_matrix.py"
    content = sp_matrix_file.read_text()

    # 统计代码行数
    lines = content.split('\n')
    code_lines = [l for l in lines if l.strip() and not l.strip().startswith('#')]
    comment_lines = [l for l in lines if l.strip().startswith('"""') or l.strip().startswith("#")]

    tests = [
        (len(lines) > 200, f"代码行数 > 200 (实际: {len(lines)})"),
        (len(code_lines) > 100, f"有效代码行数 > 100 (实际: {len(code_lines)})"),
        ("from typing import" in content, "使用 typing 注解"),
        (": Path" in content or "-> Path" in content, "使用 Path 类型标注"),
        (": float" in content or "-> float" in content, "使用 float 类型标注"),
        (": Optional[" in content, "使用 Optional 类型标注"),
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


if __name__ == "__main__":
    print("\n")
    print("*" * 60)
    print("*" + " " * 58 + "*")
    print("*" + " " * 8 + "Session #9: SPMatrixBuilder 结构测试" + " " * 18 + "*")
    print("*" + " " * 58 + "*")
    print("*" * 60)

    results = []
    results.append(("文件结构", test_sp_matrix_file_exists()))
    results.append(("初始化方法", test_sp_matrix_init_structure()))
    results.append(("run 方法", test_sp_matrix_run_structure()))
    results.append(("_run_sp 方法", test_sp_matrix_run_sp_structure()))
    results.append(("_save_report 方法", test_sp_matrix_save_report_structure()))
    results.append(("文档字符串", test_sp_matrix_docstrings()))
    results.append(("代码质量", test_sp_matrix_code_quality()))

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
