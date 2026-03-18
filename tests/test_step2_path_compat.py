from pathlib import Path

import pytest

from rph_core.orchestrator import ReactionProfileHunter


@pytest.fixture
def hunter() -> ReactionProfileHunter:
    return ReactionProfileHunter(config_path=Path("config/defaults.yaml"))


def test_resolve_product_xyz_direct_file(tmp_path: Path, hunter: ReactionProfileHunter) -> None:
    product = tmp_path / "product.xyz"
    product.write_text("2\nproduct\nC 0 0 0\nC 0 0 1\n", encoding="utf-8")

    resolved = hunter._resolve_product_xyz_for_s2(product)
    assert resolved == product


def test_resolve_product_xyz_v61_flat(tmp_path: Path, hunter: ReactionProfileHunter) -> None:
    s1 = tmp_path / "S1_ConfGeneration"
    s1.mkdir()
    product = s1 / "product_min.xyz"
    product.write_text("2\nproduct\nC 0 0 0\nC 0 0 1\n", encoding="utf-8")

    resolved = hunter._resolve_product_xyz_for_s2(s1)
    assert resolved == product


def test_resolve_product_xyz_v30_subdir(tmp_path: Path, hunter: ReactionProfileHunter) -> None:
    s1 = tmp_path / "S1_ConfGeneration"
    product_dir = s1 / "product"
    product_dir.mkdir(parents=True)
    product = product_dir / "global_min.xyz"
    product.write_text("2\nproduct\nC 0 0 0\nC 0 0 1\n", encoding="utf-8")

    resolved = hunter._resolve_product_xyz_for_s2(s1)
    assert resolved == product


def test_resolve_product_xyz_missing_raises(tmp_path: Path, hunter: ReactionProfileHunter) -> None:
    s1 = tmp_path / "S1_ConfGeneration"
    s1.mkdir()

    with pytest.raises(RuntimeError, match="无法在"):
        hunter._resolve_product_xyz_for_s2(s1)


def test_resolve_product_xyz_nonexistent_raises(tmp_path: Path, hunter: ReactionProfileHunter) -> None:
    with pytest.raises(FileNotFoundError, match="S2 输入产物路径不存在"):
        hunter._resolve_product_xyz_for_s2(tmp_path / "does_not_exist")
