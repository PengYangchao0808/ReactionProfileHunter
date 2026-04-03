from pathlib import Path

from rph_core.utils.result_inspector import ResultInspector


def test_s3_finds_product_in_canonical_s1(tmp_path: Path) -> None:
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir(parents=True)

    product_xyz = s1_dir / "product_min.xyz"
    product_xyz.write_text("3\ntest\nC 0 0 0\nO 1 0 0\n", encoding="utf-8")

    sp_dir = s1_dir / "product" / "dft"
    sp_dir.mkdir(parents=True)
    sp_out = sp_dir / "conf_000_SP.out"
    sp_out.write_text("ORCA TERMINATED NORMALLY\nFINAL SINGLE POINT ENERGY -229.0", encoding="utf-8")

    inspector = ResultInspector(tmp_path, {})
    result = inspector.check_step("s1")

    assert result.should_skip is True
    assert "s1_complete" in result.reason


def test_s3_rejects_legacy_s1_product_layout(tmp_path: Path) -> None:
    legacy_dir = tmp_path / "S1_Product"
    legacy_dir.mkdir(parents=True)
    (legacy_dir / "product_global_min.xyz").write_text("3\ntest\nC 0 0 0\nO 1 0 0\n", encoding="utf-8")

    inspector = ResultInspector(tmp_path, {})
    result = inspector.check_step("s1")

    assert result.should_skip is False
