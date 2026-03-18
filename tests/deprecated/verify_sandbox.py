"""
简单的沙盒机制验证脚本
======================

直接运行验证，不依赖 pytest
"""

import sys
import tempfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT))

from rph_core.utils.qc_interface import is_path_toxic, run_in_sandbox


def test_toxic_path_detection():
    """测试有毒路径检测"""
    print("=" * 60)
    print("测试有毒路径检测")
    print("=" * 60)

    test_cases = [
        ("/home/user/my calculations/test", True, "包含空格"),
        ("/home/user/[test]/results", True, "包含方括号"),
        ("/home/user/test(dir)", True, "包含圆括号"),
        ("/home/user/calc/test", False, "干净路径"),
        ("/home/user/test{dir}", True, "包含花括号"),
    ]

    all_passed = True
    for path_str, expected, desc in test_cases:
        path = Path(path_str)
        result = is_path_toxic(path)
        passed = result == expected
        status = "✓ PASS" if passed else "✗ FAIL"
        result_str = str(result)
        expected_str = str(expected)
        print(f"{status}: {desc:30s} -> {result_str:5s} (期望: {expected_str})")
        if not passed:
            all_passed = False

    print()
    if all_passed:
        print("✓ 有毒路径检测测试全部通过!")
    else:
        print("✗ 有毒路径检测测试有失败!")
    return all_passed


def test_sandbox_basic():
    """测试基本沙盒执行"""
    print("=" * 60)
    print("测试基本沙盒执行")
    print("=" * 60)

    with tempfile.TemporaryDirectory(prefix='test_sandbox_') as temp_dir:
        base = Path(temp_dir)
        input_file = base / "input.txt"
        input_file.write_text("test content")

        output_dir = base / "output"
        output_dir.mkdir()

        # 创建一个简单的命令来测试沙盒
        # 使用 python 来创建输出文件
        import os
        python_cmd = [sys.executable, "-c", "with open('output.txt', 'w') as f: f.write('success')"]

        try:
            result = run_in_sandbox(
                cmd=python_cmd,
                input_files={"input.txt": input_file},
                output_files=["output.txt"],
                output_dir=output_dir,
                timeout=10
            )

            if result.returncode == 0 and (output_dir / "output.txt").exists():
                print("✓ 基本沙盒执行测试通过!")
                return True
            else:
                print(f"✗ 沙盒执行失败: returncode={result.returncode}")
                return False
        except Exception as e:
            print(f"✗ 沙盒执行异常: {e}")
            return False


def test_sandbox_with_toxic_path():
    """测试在有毒路径下使用沙盒"""
    print("=" * 60)
    print("测试有毒路径下的沙盒执行")
    print("=" * 60)

    with tempfile.TemporaryDirectory(prefix='test_toxic_path_') as base_dir:
        base = Path(base_dir)

        # 创建包含特殊字符的路径
        toxic_base = base / "[5+2] Test Space"
        toxic_base.mkdir()

        input_file = base / "input.txt"
        input_file.write_text("test content")

        output_dir = toxic_base / "output"
        output_dir.mkdir()

        import os
        python_cmd = [sys.executable, "-c", "with open('output.txt', 'w') as f: f.write('success from sandbox')"]

        try:
            print(f"有毒路径: {toxic_base}")
            print(f"输出目录: {output_dir}")

            result = run_in_sandbox(
                cmd=python_cmd,
                input_files={"input.txt": input_file},
                output_files=["output.txt"],
                output_dir=output_dir,
                timeout=10
            )

            if result.returncode == 0:
                output_file = output_dir / "output.txt"
                if output_file.exists():
                    content = output_file.read_text()
                    print(f"✓ 有毒路径沙盒执行测试通过!")
                    print(f"  输出文件内容: {content}")
                    return True
                else:
                    print(f"✗ 输出文件未创建: {output_file}")
                    return False
            else:
                print(f"✗ 沙盒执行失败: returncode={result.returncode}")
                return False
        except Exception as e:
            print(f"✗ 沙盒执行异常: {e}")
            import traceback
            traceback.print_exc()
            return False


def main():
    """运行所有测试"""
    print("\n" + "=" * 60)
    print("沙盒机制验证测试套件")
    print("=" * 60 + "\n")

    results = []

    # 运行测试
    results.append(("有毒路径检测", test_toxic_path_detection()))
    print()
    results.append(("基本沙盒执行", test_sandbox_basic()))
    print()
    results.append(("有毒路径沙盒执行", test_sandbox_with_toxic_path()))
    print()

    # 总结
    print("=" * 60)
    print("测试总结")
    print("=" * 60)

    passed = sum(1 for _, r in results if r)
    total = len(results)

    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {name}")

    print()
    print(f"通过: {passed}/{total}")

    if passed == total:
        print("\n✓ 所有测试通过!")
        return 0
    else:
        print(f"\n✗ {total - passed} 个测试失败!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
