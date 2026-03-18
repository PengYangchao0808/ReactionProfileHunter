from pathlib import Path

from rph_core.utils.qc_interface import is_path_toxic


def test_path_with_space_is_toxic() -> None:
    assert is_path_toxic(Path("/tmp/my calculations/test")) is True


def test_path_with_brackets_is_toxic() -> None:
    assert is_path_toxic(Path("/tmp/[4+3]/test")) is True


def test_clean_path_is_not_toxic() -> None:
    assert is_path_toxic(Path("/tmp/rph/test")) is False
