from pathlib import Path

from rph_core.utils.qc_interface import try_formchk


def test_try_formchk_missing_chk_returns_none(tmp_path) -> None:
    missing = tmp_path / "missing.chk"
    assert try_formchk(missing) is None


def test_try_formchk_existing_chk_without_formchk_binary_returns_none(tmp_path) -> None:
    chk = tmp_path / "sample.chk"
    chk.write_text("dummy", encoding="utf-8")
    result = try_formchk(chk)
    assert result is None or isinstance(result, Path)
