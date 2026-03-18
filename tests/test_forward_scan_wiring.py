from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock
from typing import Any, cast
import logging

import pytest

from rph_core.orchestrator import ReactionProfileHunter, _resolve_run_config
from rph_core.steps.runners import run_step2


def _write_product_xyz(path: Path, natoms: int = 6) -> None:
    lines = [str(natoms), "product"]
    for idx in range(natoms):
        lines.append(f"C {float(idx):.3f} 0.000 0.000")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_resolve_run_config_includes_reaction_profile() -> None:
    config = {"global": {}, "run": {}}
    args = SimpleNamespace(
        output=None,
        smiles=None,
        reaction_type="[4+3]_default",
    )

    run_cfg = _resolve_run_config(config, args)

    assert run_cfg["reaction_profile"] == "[4+3]_default"
    assert run_cfg["reaction_type"] == "[4+3]_default"


def test_profile_drives_scan_parameters() -> None:
    hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))

    scan_cfg = hunter._resolve_forward_scan_config(reaction_profile="[4+3]_default")

    assert scan_cfg["scan_start_distance"] == 3.5
    assert scan_cfg["scan_end_distance"] == 1.8
    assert scan_cfg["scan_steps"] == 20
    assert scan_cfg["scan_mode"] == "concerted"
    assert scan_cfg["scan_force_constant"] == 0.5
    assert scan_cfg["reject_boundary_maximum"] is True
    assert scan_cfg["boundary_retry_once"] is True
    assert scan_cfg["allow_boundary_degradation"] is True


def test_forming_bonds_prefers_cleaner_data() -> None:
    hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))

    forming_bonds = cast(Any, hunter)._resolve_forming_bonds_for_s2(
        cleaner_data={"formed_bond_index_pairs": "0-1;2-3"},
    )

    assert forming_bonds == ((0, 1), (2, 3))


def test_forming_bonds_parses_comma_delimited_string() -> None:
    hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))

    forming_bonds = cast(Any, hunter)._resolve_forming_bonds_for_s2(
        cleaner_data={"formed_bond_index_pairs": "0-1,2-3"},
    )

    assert forming_bonds == ((0, 1), (2, 3))


def test_forming_bonds_resolved_from_map_pairs_and_product_mapping(tmp_path: Path) -> None:
    hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))
    product_xyz = tmp_path / "product.xyz"
    _write_product_xyz(product_xyz, natoms=4)

    forming_bonds = cast(Any, hunter)._resolve_forming_bonds_for_s2(
        cleaner_data={
            "formed_bond_map_pairs": "1-2;3-4",
            "mapped_product_smiles": "[CH3:1][CH2:2][CH2:3][CH3:4]",
        },
        product_xyz_file=product_xyz,
    )

    assert forming_bonds == ((0, 1), (2, 3))


def test_forming_bonds_one_based_requires_explicit_base_and_converts(tmp_path: Path) -> None:
    hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))
    product_xyz = tmp_path / "product.xyz"
    _write_product_xyz(product_xyz, natoms=8)

    forming_bonds = cast(Any, hunter)._resolve_forming_bonds_for_s2(
        cleaner_data={
            "formed_bond_index_pairs": "1-2;3-4",
            "index_base": 1,
        },
        product_xyz_file=product_xyz,
    )

    assert forming_bonds == ((0, 1), (2, 3))


def test_forming_bonds_level3_fallback_warns_for_long_product_distance(tmp_path: Path, caplog) -> None:
    hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))
    product_xyz = tmp_path / "product.xyz"
    _write_product_xyz(product_xyz, natoms=4)

    lines = [
        "4",
        "product",
        "C 0.000 0.000 0.000",
        "C 1.000 0.000 0.000",
        "C 0.000 1.000 0.000",
        "C 4.000 0.000 0.000",
    ]
    product_xyz.write_text("\n".join(lines) + "\n", encoding="utf-8")

    caplog.set_level(logging.WARNING)
    forming_bonds = cast(Any, hunter)._resolve_forming_bonds_for_s2(
        cleaner_data={"formed_bond_index_pairs": "0-3;1-2", "index_base": 0},
        product_xyz_file=product_xyz,
    )

    assert forming_bonds == ((0, 3), (1, 2))
    assert any(
        "possible index mapping error" in rec.message and "(0, 3)" in rec.message
        for rec in caplog.records
    )


def test_forming_bonds_ambiguous_without_index_base_raises(tmp_path: Path) -> None:
    hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))
    product_xyz = tmp_path / "product.xyz"
    _write_product_xyz(product_xyz, natoms=8)

    with pytest.raises(RuntimeError, match="ambiguous"):
        cast(Any, hunter)._resolve_forming_bonds_for_s2(
            cleaner_data={"formed_bond_index_pairs": "1-2;3-4"},
            product_xyz_file=product_xyz,
        )


def test_forming_bonds_missing_raises() -> None:
    hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))
    hunter.config.setdefault("step2", {}).pop("forming_bonds", None)

    with pytest.raises(RuntimeError, match="requires forming_bonds"):
        cast(Any, hunter)._resolve_forming_bonds_for_s2(cleaner_data={})


def test_resolve_run_config_applies_rx_id_filter() -> None:
    config = {"global": {}, "run": {"source": "dataset", "filter_ids": []}}
    args = SimpleNamespace(
        output=None,
        smiles=None,
        reaction_type=None,
        rx_id="rx_9422028",
    )

    run_cfg = _resolve_run_config(config, args)

    assert run_cfg["filter_rx_id"] == "9422028"
    assert "9422028" in run_cfg["filter_ids"]
    assert "rx_9422028" in run_cfg["filter_ids"]


def test_run_step2_uses_single_workflow(tmp_path: Path) -> None:
    hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))
    hunter._resolve_product_xyz_for_s2 = MagicMock(return_value=tmp_path / "product.xyz")
    cast(Any, hunter)._resolve_forming_bonds_for_s2 = MagicMock(return_value=((0, 1), (2, 3)))
    hunter._resolve_forward_scan_config = MagicMock(return_value={"scan_start_distance": 3.2, "scan_end_distance": 1.8})
    hunter._s2_engine = MagicMock()
    runner_return = (
        tmp_path / "ts_guess.xyz",
        tmp_path / "reactant_complex.xyz",
        tmp_path / "dipolar_intermediate.xyz",
        ((0, 1), (2, 3)),
        tmp_path / "scan_profile.json",
        "COMPLETE",
        "high",
        tuple(),
    )
    hunter._s2_engine.run.return_value = runner_return
    hunter._s2_engine.run_forward_scan.return_value = runner_return

    result = run_step2(
        hunter=hunter,
        product_xyz=tmp_path / "product.xyz",
        work_dir=tmp_path,
        reaction_profile="[4+3]_default",
        cleaner_data={"formed_bond_index_pairs": "0-1;2-3"},
    )

    assert getattr(result, "generation_method") == "forward_scan"
    assert getattr(result, "dipolar_intermediate_xyz") == tmp_path / "dipolar_intermediate.xyz"
    assert result.forming_bonds == ((0, 1), (2, 3))
    assert result.status == "COMPLETE"
    assert result.ts_guess_confidence == "high"
    hunter._s2_engine.run_forward_scan.assert_called_once()


def test_run_step2_uses_profile_strategy_for_generation_method(tmp_path: Path) -> None:
    hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))
    hunter.config.setdefault("reaction_profiles", {})["[4+3]_default"] = {
        "s2_strategy": "retro_scan",
        "scan": {
            "scan_start_distance": 3.5,
            "scan_end_distance": 1.8,
            "scan_steps": 20,
            "scan_mode": "concerted",
            "scan_force_constant": 0.5,
        },
    }
    hunter._resolve_product_xyz_for_s2 = MagicMock(return_value=tmp_path / "product.xyz")
    cast(Any, hunter)._resolve_forming_bonds_for_s2 = MagicMock(return_value=((0, 1), (2, 3)))
    hunter._resolve_forward_scan_config = MagicMock(return_value={"scan_start_distance": 3.5, "scan_end_distance": 1.8})
    hunter._s2_engine = MagicMock()
    runner_return = (
        tmp_path / "ts_guess.xyz",
        tmp_path / "reactant_complex.xyz",
        tmp_path / "dipolar_intermediate.xyz",
        ((0, 1), (2, 3)),
        tmp_path / "scan_profile.json",
        "COMPLETE",
        "high",
        tuple(),
    )
    hunter._s2_engine.run.return_value = runner_return
    hunter._s2_engine.run_retro_scan.return_value = runner_return

    result = run_step2(
        hunter=hunter,
        product_xyz=tmp_path / "product.xyz",
        work_dir=tmp_path,
        reaction_profile="[4+3]_default",
        cleaner_data={"formed_bond_index_pairs": "0-1;2-3"},
    )

    assert result.generation_method == "retro_scan"
    hunter._s2_engine.run_retro_scan.assert_called_once()
