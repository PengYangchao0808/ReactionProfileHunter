"""
Session #2 验证脚本 - 测试 ORCAInterface._parse_output()
"""
import sys
import tempfile
from pathlib import Path

# 添加项目路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from rph_core.utils.orca_interface import ORCAInterface, QCResult


def test_parse_output_energy():
    """测试成功解析输出能量"""
    print("=" * 60)
    print("测试: ORCAInterface._parse_output() - 成功案例")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # 创建成功的 ORCA 输出文件
        out_file = tmp_path / "orca_success.out"
        out_file.write_text("""
-------------------
ORCA CALCULATION DONE
-------------------

Total Energy        :     -123.45678901
FINAL SINGLE POINT ENERGY      -123.45678901

ORCA TERMINATED NORMALLY
""")

        # 解析输出
        orca = ORCAInterface()
        result = orca._parse_output(out_file)

        # 验证结果
        print("\n验证结果:")
        tests = [
            (result.converged == True, "converged = True"),
            (result.energy == -123.45678901, f"energy = {result.energy}"),
            (result.error_message is None, "error_message = None"),
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


def test_parse_output_failed():
    """测试解析失败输出"""
    print("\n\n")
    print("=" * 60)
    print("测试: ORCAInterface._parse_output() - 失败案例")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # 创建失败的 ORCA 输出文件
        out_file = tmp_path / "orca_failed.out"
        out_file.write_text("""
ORCA calculation encountered an error

ERROR: SCF failed to converge

ORCA TERMINATED WITH ERRORS
""")

        # 解析输出
        orca = ORCAInterface()
        result = orca._parse_output(out_file)

        # 验证结果
        print("\n验证结果:")
        tests = [
            (result.converged == False, "converged = False"),
            (result.energy == 0.0, f"energy = {result.energy}"),
            (result.error_message is not None, "error_message 存在"),
            ("未正常终止" in result.error_message, "错误消息包含'未正常终止'"),
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


def test_parse_output_no_energy():
    """测试缺少能量信息"""
    print("\n\n")
    print("=" * 60)
    print("测试: ORCAInterface._parse_output() - 缺少能量")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # 创建缺少能量信息的 ORCA 输出文件
        out_file = tmp_path / "orca_no_energy.out"
        out_file.write_text("""
ORCA calculation completed

ORCA TERMINATED NORMALLY
""")

        # 解析输出
        orca = ORCAInterface()
        result = orca._parse_output(out_file)

        # 验证结果
        print("\n验证结果:")
        tests = [
            (result.converged == False, "converged = False"),
            (result.energy == 0.0, f"energy = {result.energy}"),
            ("无法找到能量信息" in result.error_message, "错误消息包含'无法找到能量信息'"),
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
    print("*" + " " * 10 + "Session #2: _parse_output() 测试套件" + " " * 18 + "*")
    print("*" + " " * 58 + "*")
    print("*" * 60)

    results = []
    results.append(("成功解析能量", test_parse_output_energy()))
    results.append(("失败输出处理", test_parse_output_failed()))
    results.append(("缺少能量信息", test_parse_output_no_energy()))

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
        print("✓✓✓ 所有测试通过! Session #2 完成 ✓✓✓")
    else:
        print("✗✗✗ 部分测试失败 ✗✗✗")
    print("*" * 60)
    print()

    sys.exit(0 if all_passed else 1)
