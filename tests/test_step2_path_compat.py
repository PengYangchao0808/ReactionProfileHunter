import json
from pathlib import Path

import pytest

from rph_core.orchestrator import ReactionProfileHunter


@pytest.fixture
def hunter() -> ReactionProfileHunter:
    return ReactionProfileHunter(config_path=Path("config/defaults.yaml"))


def _write_product_xyz(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "4\nproduct\nC 0.0 0.0 0.0\nC 0.0 0.0 3.0\nC 2.0 0.0 0.0\nC 2.0 0.0 3.0\n",
        encoding="utf-8",
    )


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


def test_resolve_product_xyz_missing_raises(tmp_path: Path, hunter: ReactionProfileHunter) -> None:
    s1 = tmp_path / "S1_ConfGeneration"
    s1.mkdir()

    with pytest.raises(RuntimeError, match="无法在"):
        hunter._resolve_product_xyz_for_s2(s1)


def test_resolve_product_xyz_nonexistent_raises(tmp_path: Path, hunter: ReactionProfileHunter) -> None:
    with pytest.raises(FileNotFoundError, match="S2 输入产物路径不存在"):
        hunter._resolve_product_xyz_for_s2(tmp_path / "does_not_exist")


def test_resolve_forming_bonds_for_s2_prefers_xyz_mapping_artifact(
    tmp_path: Path,
    hunter: ReactionProfileHunter,
) -> None:
    product = tmp_path / "S1_ConfGeneration" / "product_min.xyz"
    _write_product_xyz(product)

    s0_dir = tmp_path / "S0_Mechanism"
    s0_dir.mkdir(parents=True, exist_ok=True)
    (s0_dir / "mechanism_graph.json").write_text(
        json.dumps({"edges": [{"pathway_id": "primary", "forming_bonds": [[9, 10], [11, 12]]}]}, indent=2),
        encoding="utf-8",
    )

    (tmp_path / "S1_ConfGeneration" / "atom_map_xyz.json").write_text(
        json.dumps(
            {
                "map_to_product_xyz_1based": {"5": 1, "7": 2, "9": 3, "11": 4},
                "forming_bonds_full_notation": [
                    {"map_space": [5, 9], "product_xyz_1based": [1, 3]},
                    {"map_space": [7, 11], "product_xyz_1based": [2, 4]},
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    resolved = hunter._resolve_forming_bonds_for_s2(
        cleaner_data={"formed_bond_xyz_pairs": [[0, 1], [2, 3]]},
        product_xyz_file=product,
        work_dir=tmp_path,
    )

    assert resolved == ((0, 2), (1, 3))


def test_resolve_forming_bonds_for_s2_uses_graph_edges_with_xyz_mapping(
    tmp_path: Path,
    hunter: ReactionProfileHunter,
) -> None:
    product = tmp_path / "S1_ConfGeneration" / "product_min.xyz"
    _write_product_xyz(product)

    s0_dir = tmp_path / "S0_Mechanism"
    s0_dir.mkdir(parents=True, exist_ok=True)
    (s0_dir / "mechanism_graph.json").write_text(
        json.dumps(
            {
                "edges": [{"pathway_id": "primary", "forming_bonds": [[0, 2], [1, 3]]}],
                "smiles_atom_mapping": {
                    "product_smiles_to_map": {"0": 5, "1": 7, "2": 9, "3": 11},
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (tmp_path / "S1_ConfGeneration" / "atom_map_xyz.json").write_text(
        json.dumps(
            {
                "map_to_product_xyz_1based": {"5": 1, "7": 2, "9": 3, "11": 4},
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    resolved = hunter._resolve_forming_bonds_for_s2(
        cleaner_data={"forming_bonds": [[0, 1], [2, 3]]},
        product_xyz_file=product,
        work_dir=tmp_path,
    )

    assert resolved == ((0, 2), (1, 3))
