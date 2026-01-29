"""
Session #8 验证脚本 - 测试 SPMatrixReport 数据结构
"""
import sys
import json
from pathlib import Path

# 添加项目路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

# 直接导入，避免触发 steps/__init__.py 的依赖
import importlib.util
spec = importlib.util.spec_from_file_location(
    "sp_report",
    project_root / "rph_core" / "steps" / "step3_5_sp" / "sp_report.py"
)
sp_report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sp_report)

SPMatrixReport = sp_report.SPMatrixReport


def test_sp_report_creation():
    """测试 SPMatrixReport 创建"""
    print("=" * 60)
    print("测试: SPMatrixReport 创建")
    print("=" * 60)

    report = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901
    )

    # 验证基本属性
    tests = [
        (report.e_product == -123.45678901, "e_product = -123.45678901"),
        (report.e_reactant == -456.78901234, "e_reactant = -456.78901234"),
        (report.e_ts_final == -345.67890123, "e_ts_final = -345.67890123"),
        (report.e_frag_a_ts == -234.56789012, "e_frag_a_ts = -234.56789012"),
        (report.e_frag_b_ts == -123.45678901, "e_frag_b_ts = -123.45678901"),
        (report.e_frag_a_relaxed is None, "e_frag_a_relaxed = None"),
        (report.e_frag_b_relaxed is None, "e_frag_b_relaxed = None"),
        (report.method == "M062X/def2-TZVPP", "method = M062X/def2-TZVPP"),
        (report.solvent == "acetone", "solvent = acetone"),
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


def test_sp_report_to_dict():
    """测试转换为字典"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixReport.to_dict()")
    print("=" * 60)

    report = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901,
        e_frag_a_relaxed=-235.00000000
    )

    data = report.to_dict()

    # 验证字典内容
    tests = [
        (data["e_product_l2"] == -123.45678901, "e_product_l2 = -123.45678901"),
        (data["e_reactant_l2"] == -456.78901234, "e_reactant_l2 = -456.78901234"),
        (data["e_ts_final_l2"] == -345.67890123, "e_ts_final_l2 = -345.67890123"),
        (data["e_frag_a_ts_l2"] == -234.56789012, "e_frag_a_ts_l2 = -234.56789012"),
        (data["e_frag_b_ts_l2"] == -123.45678901, "e_frag_b_ts_l2 = -123.45678901"),
        (data["e_frag_a_relaxed_l2"] == -235.00000000, "e_frag_a_relaxed_l2 = -235.00000000"),
        (data["e_frag_b_relaxed_l2"] is None, "e_frag_b_relaxed_l2 = None"),
        (len(data) == 9, f"字典键数量 = 9 (实际: {len(data)})"),
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


def test_sp_report_to_json():
    """测试转换为 JSON"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixReport.to_json()")
    print("=" * 60)

    report = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901
    )

    json_str = report.to_json()

    # 验证 JSON 格式
    data = json.loads(json_str)

    tests = [
        (data["e_product_l2"] == -123.45678901, "JSON 包含正确的 e_product_l2"),
        (data["e_reactant_l2"] == -456.78901234, "JSON 包含正确的 e_reactant_l2"),
        ("\n" in json_str, "JSON 包含缩进格式"),
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


def test_sp_report_from_dict():
    """测试从字典创建"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixReport.from_dict()")
    print("=" * 60)

    data = {
        "e_product_l2": -123.45678901,
        "e_reactant_l2": -456.78901234,
        "e_ts_final_l2": -345.67890123,
        "e_frag_a_ts_l2": -234.56789012,
        "e_frag_b_ts_l2": -123.45678901,
        "e_frag_a_relaxed_l2": -235.00000000,
        "e_frag_b_relaxed_l2": None,
        "method": "PWPB95/def2-TZVPP",
        "solvent": "toluene"
    }

    report = SPMatrixReport.from_dict(data)

    tests = [
        (report.e_product == -123.45678901, "e_product = -123.45678901"),
        (report.e_reactant == -456.78901234, "e_reactant = -456.78901234"),
        (report.e_ts_final == -345.67890123, "e_ts_final = -345.67890123"),
        (report.e_frag_a_ts == -234.56789012, "e_frag_a_ts = -234.56789012"),
        (report.e_frag_b_ts == -123.45678901, "e_frag_b_ts = -123.45678901"),
        (report.e_frag_a_relaxed == -235.00000000, "e_frag_a_relaxed = -235.00000000"),
        (report.e_frag_b_relaxed is None, "e_frag_b_relaxed = None"),
        (report.method == "PWPB95/def2-TZVPP", "method = PWPB95/def2-TZVPP"),
        (report.solvent == "toluene", "solvent = toluene"),
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


def test_sp_report_validate():
    """测试数据验证"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixReport.validate()")
    print("=" * 60)

    # 有效数据
    valid_report = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901
    )

    # 无效数据 (包含字符串)
    invalid_report = SPMatrixReport(
        e_product="invalid",  # type: ignore
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901
    )

    tests = [
        (valid_report.validate() is True, "有效报告通过验证"),
        (invalid_report.validate() is False, "无效报告验证失败"),
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


def test_sp_report_calculation_methods():
    """测试计算方法"""
    print("\n\n")
    print("=" * 60)
    print("测试: SPMatrixReport 计算方法")
    print("=" * 60)

    report = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901,
        e_frag_a_relaxed=-235.00000000
    )

    # 计算活化能、反应能、畸变能
    activation = report.get_activation_energy()
    reaction = report.get_reaction_energy()
    distortion_a = report.get_distortion_energy_a()

    tests = [
        (isinstance(activation, float), "活化能为浮点数"),
        (isinstance(reaction, float), "反应能为浮点数"),
        (isinstance(distortion_a, float), "畸变能为浮点数"),
        (abs(activation - 69723.4) < 1.0, f"活化能 ≈ 69723.4 kcal/mol (实际: {activation:.1f})"),
        (abs(reaction - 209168.7) < 1.0, f"反应能 ≈ 209168.7 kcal/mol (实际: {reaction:.1f})"),
        (abs(distortion_a - 271.1) < 1.0, f"畸变能A ≈ 271.1 kcal/mol (实际: {distortion_a:.1f})"),
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


def test_sp_report_serialization_roundtrip():
    """测试序列化/反序列化往返"""
    print("\n\n")
    print("=" * 60)
    print("测试: 序列化/反序列化往返")
    print("=" * 60)

    original = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901,
        e_frag_a_relaxed=-235.00000000,
        method="PWPB95/def2-TZVPP",
        solvent="water"
    )

    # 序列化
    json_str = original.to_json()
    data = json.loads(json_str)

    # 反序列化
    restored = SPMatrixReport.from_dict(data)

    # 验证所有字段一致
    tests = [
        (restored.e_product == original.e_product, "e_product 一致"),
        (restored.e_reactant == original.e_reactant, "e_reactant 一致"),
        (restored.e_ts_final == original.e_ts_final, "e_ts_final 一致"),
        (restored.e_frag_a_ts == original.e_frag_a_ts, "e_frag_a_ts 一致"),
        (restored.e_frag_b_ts == original.e_frag_b_ts, "e_frag_b_ts 一致"),
        (restored.e_frag_a_relaxed == original.e_frag_a_relaxed, "e_frag_a_relaxed 一致"),
        (restored.method == original.method, "method 一致"),
        (restored.solvent == original.solvent, "solvent 一致"),
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
    print("*" + " " * 10 + "Session #8: SPMatrixReport 测试套件" + " " * 17 + "*")
    print("*" + " " * 58 + "*")
    print("*" * 60)

    results = []
    results.append(("SPMatrixReport 创建", test_sp_report_creation()))
    results.append(("to_dict() 方法", test_sp_report_to_dict()))
    results.append(("to_json() 方法", test_sp_report_to_json()))
    results.append(("from_dict() 方法", test_sp_report_from_dict()))
    results.append(("validate() 方法", test_sp_report_validate()))
    results.append(("计算方法", test_sp_report_calculation_methods()))
    results.append(("序列化/反序列化", test_sp_report_serialization_roundtrip()))

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
        print("✓✓✓ 所有测试通过! Session #8 完成 ✓✓✓")
    else:
        print("✗✗✗ 部分测试失败 ✗✗✗")
    print("*" * 60)
    print()

    sys.exit(0 if all_passed else 1)
