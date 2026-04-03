"""
测试 SPMatrixReport 数据结构
==============================

Session #8: 测试 SPMatrixReport 类

Author: QC Descriptors Team
Date: 2026-01-11
"""

import pytest
import json
from rph_core.steps.step3_opt.ts_optimizer import SPMatrixReport


def test_sp_report_creation():
    """测试 SPMatrixReport 创建"""
    report = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901,
        method="M062X/def2-TZVPP",
        solvent="acetone"
    )

    # 验证基本属性
    assert report.e_product == -123.45678901
    assert report.e_reactant == -456.78901234
    assert report.e_ts_final == -345.67890123
    assert report.e_frag_a_ts == -234.56789012
    assert report.e_frag_b_ts == -123.45678901

    # 验证默认值
    assert report.e_frag_a_relaxed is None
    assert report.e_frag_b_relaxed is None
    assert report.method == "M062X/def2-TZVPP"
    assert report.solvent == "acetone"


def test_sp_report_with_optional_fields():
    """测试包含可选字段的 SPMatrixReport"""
    report = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901,
        e_frag_a_relaxed=-235.00000000,
        e_frag_b_relaxed=-124.00000000,
        method="PWPB95/def2-TZVPP",
        solvent="water"
    )

    assert report.e_frag_a_relaxed == -235.00000000
    assert report.e_frag_b_relaxed == -124.00000000
    assert report.method == "PWPB95/def2-TZVPP"
    assert report.solvent == "water"


def test_sp_report_to_dict():
    """测试转换为字典 - V4.1 dataclass field names"""
    report = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901,
        e_frag_a_relaxed=-235.00000000
    )

    data = report.to_dict()

    # V4.1: to_dict() now uses dataclass field names (e_product, not e_product_l2)
    assert "e_product" in data
    assert "e_reactant" in data
    assert "e_ts_final" in data
    assert "e_frag_a_ts" in data
    assert "e_frag_b_ts" in data
    assert "e_frag_a_relaxed" in data
    assert "method" in data
    assert "solvent" in data
    # g_* fields should exist (even if None)
    assert "g_ts" in data
    assert "g_reactant" in data
    assert "g_product" in data


def test_sp_report_to_json():
    """测试转换为 JSON"""
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
    assert data["e_product"] == -123.45678901
    assert data["e_reactant"] == -456.78901234

    # 验证缩进
    assert "\n" in json_str  # 有换行符表示有缩进


def test_sp_report_from_dict():
    """测试从字典创建 - V4.1 dataclass field names"""
    # V4.1: from_dict() filters unknown keys, so input should use valid field names
    data = {
        "e_product": -123.45678901,
        "e_reactant": -456.78901234,
        "e_ts_final": -345.67890123,
        "e_frag_a_ts": -234.56789012,
        "e_frag_b_ts": -123.45678901,
        "e_frag_a_relaxed": -235.00000000,
        "e_frag_b_relaxed": None,
        "method": "PWPB95/def2-TZVPP",
        "solvent": "toluene"
    }

    report = SPMatrixReport.from_dict(data)

    assert report.e_product == -123.45678901
    assert report.e_reactant == -456.78901234
    assert report.e_ts_final == -345.67890123
    assert report.e_frag_a_ts == -234.56789012
    assert report.e_frag_b_ts == -123.45678901
    assert report.e_frag_a_relaxed == -235.00000000
    assert report.e_frag_b_relaxed is None
    assert report.method == "PWPB95/def2-TZVPP"
    assert report.solvent == "toluene"


def test_sp_report_validate():
    """测试数据验证"""
    # 有效数据
    valid_report = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901
    )
    assert valid_report.validate() is True

    # 无效数据 (包含字符串)
    invalid_report = SPMatrixReport(
        e_product="invalid",  # pyright: ignore[reportArgumentType]
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901
    )
    assert invalid_report.validate() is False


def test_sp_report_get_activation_energy():
    """测试计算活化能 - V4.1 Gibbs priority"""
    # V4.1: Test with Gibbs energies (kcal/mol)
    report = SPMatrixReport(
        g_ts=10.0,
        g_reactant=2.0,
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123
    )

    # Gibbs priority: g_ts - g_reactant (kcal/mol)
    delta_g = report.get_activation_energy()
    assert delta_g == 8.0  # 10.0 - 2.0

def test_sp_report_get_activation_energy_fallback():
    """测试计算活化能 - Fallback to electronic energy"""
    # V4.1: Test fallback when Gibbs energies are missing
    report = SPMatrixReport(
        g_ts=None,
        g_reactant=None,
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123
    )

    # Fallback: (e_ts_final - e_reactant) * 627.509
    expected = (-345.67890123 - (-456.78901234)) * 627.509
    expected = 111.11011111 * 627.509  # ≈ 69723.4
    delta_g = report.get_activation_energy()

    assert isinstance(delta_g, float)
    assert abs(delta_g - expected) < 0.1  # 允许小误差


def test_sp_report_get_reaction_energy():
    """测试计算反应能"""
    report = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901
    )

    # ΔG_rxn = (E_product - E_reactant) * 627.509
    # = (-123.45678901 - (-456.78901234)) * 627.509
    # = 333.33222333 * 627.509
    # ≈ 209168.7 kcal/mol
    delta_g = report.get_reaction_energy()

    assert isinstance(delta_g, float)
    assert abs(delta_g - 209168.7) < 1.0  # Allow small tolerance





def test_sp_report_str():
    """测试字符串表示"""
    report = SPMatrixReport(
        e_product=-123.45678901,
        e_reactant=-456.78901234,
        e_ts_final=-345.67890123,
        e_frag_a_ts=-234.56789012,
        e_frag_b_ts=-123.45678901
    )

    str_repr = str(report)

    # 验证字符串包含关键信息
    assert "SPMatrixReport" in str_repr
    # V4.1: Default method is "Berny", not M062X/def2-TZVPP
    assert "-123.45678901" in str_repr
    assert "-456.78901234" in str_repr
    assert "ΔG‡" in str_repr or "activation" in str_repr.lower()
    assert "ΔG_rxn" in str_repr or "reaction" in str_repr.lower()


def test_sp_report_serialization_roundtrip():
    """测试序列化/反序列化往返 - V4.1 dataclass field names"""
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

    # V4.1: to_dict/from_dict use dataclass field names (e_product, not e_product_l2)
    assert restored.e_product == original.e_product
    assert restored.e_reactant == original.e_reactant
    assert restored.e_ts_final == original.e_ts_final
    assert restored.e_frag_a_ts == original.e_frag_a_ts
    assert restored.e_frag_b_ts == original.e_frag_b_ts
    assert restored.e_frag_a_relaxed == original.e_frag_a_relaxed
    assert restored.method == original.method
    assert restored.solvent == original.solvent
