"""
临时验证脚本 - 测试 ORCAInterface._generate_input()
"""
import sys
import tempfile
from pathlib import Path

# 添加项目路径
sys.path.insert(0, r'E:\Calculations\[5+2] Mechain learning\Scripts\ReactionProfileHunter\ReactionProfileHunter')

from rph_core.utils.orca_interface import ORCAInterface

def test_generate_input_m062x():
    """测试 M062X 输入生成"""
    print("=" * 60)
    print("测试: ORCAInterface._generate_input() - M062X")
    print("=" * 60)

    # 创建临时目录
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # 创建测试 XYZ 文件
        xyz_file = tmp_path / "test.xyz"
        xyz_file.write_text("""3
comment
C    0.0  0.0  0.0
H    1.0  0.0  0.0
H   -1.0  0.0  0.0
""")

        # 创建 ORCAInterface 实例
        orca = ORCAInterface(method="M062X", basis="def2-TZVPP")

        # 生成输入文件
        inp_file = orca._generate_input(xyz_file, tmp_path)

        # 读取并显示内容
        content = inp_file.read_text()

        print("\n✓ 输入文件已生成:")
        print(f"  路径: {inp_file}")
        print(f"  内容:\n")
        print("-" * 60)
        print(content)
        print("-" * 60)

        # 验证内容
        print("\n验证结果:")
        tests = [
            ("! M062X def2-TZVPP def2/J RIJCOSX" in content, "包含 M062X 路由行"),
            ("%maxcore 4000" in content, "包含 maxcore 设置"),
            ("%pal nprocs 16" in content, "包含并行设置"),
        ]

        all_passed = True
        for passed, description in tests:
            status = "✓ PASS" if passed else "✗ FAIL"
            print(f"  {status}: {description}")
            all_passed = all_passed and passed

        print("\n" + "=" * 60)
        if all_passed:
            print("✓ 所有测试通过!")
        else:
            print("✗ 部分测试失败")
        print("=" * 60)

        return all_passed

def test_generate_input_with_solvent():
    """测试溶剂模型输入生成"""
    print("\n\n")
    print("=" * 60)
    print("测试: ORCAInterface._generate_input() - 溶剂模型")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        xyz_file = tmp_path / "test.xyz"
        xyz_file.write_text("""3
comment
C    0.0  0.0  0.0
H    1.0  0.0  0.0
H   -1.0  0.0  0.0
""")

        orca = ORCAInterface(solvent="acetone")
        inp_file = orca._generate_input(xyz_file, tmp_path)

        content = inp_file.read_text()

        print("\n✓ 输入文件已生成:")
        print(f"  路径: {inp_file}")
        print(f"  内容:\n")
        print("-" * 60)
        print(content)
        print("-" * 60)

        print("\n验证结果:")
        tests = [
            ("%cpcm" in content, "包含 CPCM 块"),
            ('SMDsolvent "acetone"' in content, "包含 SMD 溶剂设置"),
            ("smd true" in content, "启用 SMD 模型"),
        ]

        all_passed = True
        for passed, description in tests:
            status = "✓ PASS" if passed else "✗ FAIL"
            print(f"  {status}: {description}")
            all_passed = all_passed and passed

        print("\n" + "=" * 60)
        if all_passed:
            print("✓ 所有测试通过!")
        else:
            print("✗ 部分测试失败")
        print("=" * 60)

        return all_passed

def test_double_hybrid():
    """测试双杂化泛函"""
    print("\n\n")
    print("=" * 60)
    print("测试: ORCAInterface._generate_input() - 双杂化泛函")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        xyz_file = tmp_path / "test.xyz"
        xyz_file.write_text("""3
comment
C    0.0  0.0  0.0
H    1.0  0.0  0.0
H   -1.0  0.0  0.0
""")

        orca = ORCAInterface(method="PWPB95", basis="def2-TZVPP")
        inp_file = orca._generate_input(xyz_file, tmp_path)

        content = inp_file.read_text()

        print("\n✓ 输入文件已生成:")
        print(f"  路径: {inp_file}")
        print(f"  内容:\n")
        print("-" * 60)
        print(content)
        print("-" * 60)

        print("\n验证结果:")
        tests = [
            ("def2-TZVPP/C" in content, "自动添加 /C 辅助基组"),
            ("PWPB95" in content, "包含双杂化泛函名称"),
        ]

        all_passed = True
        for passed, description in tests:
            status = "✓ PASS" if passed else "✗ FAIL"
            print(f"  {status}: {description}")
            all_passed = all_passed and passed

        print("\n" + "=" * 60)
        if all_passed:
            print("✓ 所有测试通过!")
        else:
            print("✗ 部分测试失败")
        print("=" * 60)

        return all_passed

if __name__ == "__main__":
    print("\n")
    print("*" * 60)
    print("*" + " " * 58 + "*")
    print("*" + " " * 15 + "ORCAInterface 测试套件" + " " * 23 + "*")
    print("*" + " " * 58 + "*")
    print("*" * 60)

    results = []
    results.append(("M062X 输入生成", test_generate_input_m062x()))
    results.append(("溶剂模型", test_generate_input_with_solvent()))
    results.append(("双杂化泛函", test_double_hybrid()))

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
        print("✓✓✓ 所有测试通过! Session #1 完成 ✓✓✓")
    else:
        print("✗✗✗ 部分测试失败 ✗✗✗")
    print("*" * 60)
    print()

    sys.exit(0 if all_passed else 1)
