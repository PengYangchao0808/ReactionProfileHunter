import json
import logging
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from rph_core.steps.runners import run_step2
from rph_core.steps.runners import _adapt_product_xyz_for_s2_if_needed


def _write_xyz(path: Path) -> None:
    path.write_text("3\nproduct\nC 0 0 0\nH 0 0 1\nH 0 1 0\n", encoding="utf-8")


def test_pes_adapter_not_triggered_for_dft_geometry(tmp_path: Path) -> None:
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir(parents=True, exist_ok=True)
    (s1_dir / "provenance.json").write_text(
        json.dumps({"protocol": "lite", "schema_version": "s1_provenance_v1", "has_geometry_optimization": True}, indent=2),
        encoding="utf-8",
    )

    product = tmp_path / "product_min.xyz"
    _write_xyz(product)

    hunter = SimpleNamespace(
        config={"step2": {"pes_adapter": {"enabled": True, "mode": "fallback"}}},
        logger=logging.getLogger("test_step2_pes_adapter"),
    )

    adapted = _adapt_product_xyz_for_s2_if_needed(
        hunter=hunter,
        work_dir=tmp_path,
        product_xyz_file=product,
    )

    assert adapted == product
    assert not (tmp_path / "S2_Retro" / "pes_adapter_fallback").exists()


def test_pes_adapter_triggered_for_non_dft_geometry(tmp_path: Path) -> None:
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir(parents=True, exist_ok=True)
    (s1_dir / "provenance.json").write_text(
        json.dumps({"protocol": "zero", "schema_version": "s1_provenance_v1", "has_geometry_optimization": False}, indent=2),
        encoding="utf-8",
    )

    product = tmp_path / "product_min.xyz"
    _write_xyz(product)

    hunter = SimpleNamespace(
        config={"step2": {"pes_adapter": {"enabled": True, "mode": "fallback"}}},
        logger=logging.getLogger("test_step2_pes_adapter"),
    )

    adapted = _adapt_product_xyz_for_s2_if_needed(
        hunter=hunter,
        work_dir=tmp_path,
        product_xyz_file=product,
    )

    assert adapted != product
    assert adapted.exists()
    assert adapted.name == "product_relaxed.xyz"

    meta_file = tmp_path / "S2_Retro" / "pes_adapter_fallback" / "adapter_meta.json"
    assert meta_file.exists()
    meta = json.loads(meta_file.read_text(encoding="utf-8"))
    assert meta["mode"] == "fallback"
    assert meta["triggered"] is True
    assert meta["reason"] == "s1_geometry_not_dft_optimized"
    assert meta["protocol"] == "zero"
    assert meta["provenance_schema_version"] == "s1_provenance_v1"
    assert meta["trigger_source"] == "s1_provenance.has_geometry_optimization=false"
    assert meta["trigger_class"] == "exception_path"


def test_run_step2_uses_fallback_adapted_product_path(tmp_path: Path) -> None:
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir(parents=True, exist_ok=True)
    (s1_dir / "provenance.json").write_text(
        json.dumps({"protocol": "zero", "schema_version": "s1_provenance_v1", "has_geometry_optimization": False}, indent=2),
        encoding="utf-8",
    )

    product = tmp_path / "product_min.xyz"
    _write_xyz(product)

    ts_guess = tmp_path / "ts_guess.xyz"
    substrate = tmp_path / "reactant_complex.xyz"
    intermediate = tmp_path / "intermediate.xyz"
    _write_xyz(ts_guess)
    _write_xyz(substrate)
    _write_xyz(intermediate)

    hunter = SimpleNamespace()
    hunter.config = {
        "step2": {
            "pes_adapter": {"enabled": True, "mode": "fallback"},
            "path_search": {"enabled": False},
        }
    }
    hunter.logger = logging.getLogger("test_step2_runner")
    hunter._resolve_product_xyz_for_s2 = MagicMock(return_value=product)
    hunter._resolve_profile_key = MagicMock(return_value="[4+3]_default")
    hunter._resolve_forming_bonds_for_s2 = MagicMock(return_value=((0, 1), (2, 3)))
    hunter._build_step2_signature = MagicMock(return_value={"ok": True})
    hunter.s2_engine = SimpleNamespace()
    hunter.s2_engine.run_retro_scan = MagicMock(
        return_value=(
            ts_guess,
            substrate,
            intermediate,
            ((0, 1), (2, 3)),
            tmp_path / "scan_profile.json",
            "COMPLETE",
            "high",
            tuple(),
            None,
        )
    )

    result = run_step2(
        hunter=hunter,
        product_xyz=product,
        work_dir=tmp_path,
        reaction_profile="[4+3]_default",
        cleaner_data=None,
    )

    assert result.ts_guess_xyz == ts_guess
    assert "pes_adapter_fallback_applied" in result.degraded_reasons
    assert "pes_adapter_fallback_exception_path" in result.degraded_reasons
    called_product = hunter.s2_engine.run_retro_scan.call_args.kwargs["product_xyz"]
    assert Path(called_product).name == "product_relaxed.xyz"


def test_run_step2_passes_resolved_xyz_bonds_without_map_override(tmp_path: Path) -> None:
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir(parents=True, exist_ok=True)
    (s1_dir / "provenance.json").write_text(
        json.dumps({"protocol": "lite", "schema_version": "s1_provenance_v1", "has_geometry_optimization": True}, indent=2),
        encoding="utf-8",
    )

    product = tmp_path / "product_min.xyz"
    _write_xyz(product)

    ts_guess = tmp_path / "ts_guess.xyz"
    substrate = tmp_path / "reactant_complex.xyz"
    intermediate = tmp_path / "intermediate.xyz"
    _write_xyz(ts_guess)
    _write_xyz(substrate)
    _write_xyz(intermediate)

    hunter = SimpleNamespace()
    hunter.config = {
        "step2": {
            "pes_adapter": {"enabled": True, "mode": "fallback"},
            "path_search": {"enabled": False},
        }
    }
    hunter.logger = logging.getLogger("test_step2_runner_contract")
    hunter._resolve_product_xyz_for_s2 = MagicMock(return_value=product)
    hunter._resolve_profile_key = MagicMock(return_value="[4+3]_default")
    hunter._resolve_forming_bonds_for_s2 = MagicMock(return_value=((0, 2), (1, 3)))
    hunter._build_step2_signature = MagicMock(return_value={"ok": True})
    hunter.s2_engine = SimpleNamespace()
    hunter.s2_engine.run_retro_scan = MagicMock(
        return_value=(
            ts_guess,
            substrate,
            intermediate,
            ((0, 2), (1, 3)),
            tmp_path / "scan_profile.json",
            "COMPLETE",
            "high",
            tuple(),
            None,
        )
    )

    run_step2(
        hunter=hunter,
        product_xyz=product,
        work_dir=tmp_path,
        reaction_profile="[4+3]_default",
        cleaner_data={"formed_bond_map_pairs": "9-10;11-12"},
    )

    retro_kwargs = hunter.s2_engine.run_retro_scan.call_args.kwargs
    assert retro_kwargs["forming_bonds"] == ((0, 2), (1, 3))
    assert "atom_map" not in retro_kwargs


def test_run_step2_path_search_reuses_resolved_xyz_bonds(tmp_path: Path) -> None:
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir(parents=True, exist_ok=True)
    (s1_dir / "provenance.json").write_text(
        json.dumps({"protocol": "lite", "schema_version": "s1_provenance_v1", "has_geometry_optimization": True}, indent=2),
        encoding="utf-8",
    )

    product = tmp_path / "product_min.xyz"
    _write_xyz(product)

    ts_guess = tmp_path / "ts_guess.xyz"
    substrate = tmp_path / "reactant_complex.xyz"
    intermediate = tmp_path / "intermediate.xyz"
    path_ts_guess = tmp_path / "path_ts_guess.xyz"
    _write_xyz(ts_guess)
    _write_xyz(substrate)
    _write_xyz(intermediate)
    _write_xyz(path_ts_guess)

    hunter = SimpleNamespace()
    hunter.config = {
        "step2": {
            "pes_adapter": {"enabled": True, "mode": "fallback"},
            "path_search": {"enabled": True},
        }
    }
    hunter.logger = logging.getLogger("test_step2_runner_path_contract")
    hunter._resolve_product_xyz_for_s2 = MagicMock(return_value=product)
    hunter._resolve_profile_key = MagicMock(return_value="[4+3]_default")
    hunter._resolve_forming_bonds_for_s2 = MagicMock(return_value=((0, 2), (1, 3)))
    hunter._build_step2_signature = MagicMock(return_value={"ok": True})
    hunter.s2_engine = SimpleNamespace()
    hunter.s2_engine.run_retro_scan = MagicMock(
        return_value=(
            ts_guess,
            substrate,
            intermediate,
            ((0, 2), (1, 3)),
            tmp_path / "scan_profile.json",
            "COMPLETE",
            "high",
            tuple(),
            None,
        )
    )
    hunter.s2_engine.run_path_search = MagicMock(
        return_value=(
            path_ts_guess,
            substrate,
            intermediate,
            ((0, 2), (1, 3)),
            tmp_path / "path_scan_profile.json",
            "COMPLETE",
            "high",
            tuple(),
            None,
        )
    )

    run_step2(
        hunter=hunter,
        product_xyz=product,
        work_dir=tmp_path,
        reaction_profile="[4+3]_default",
        cleaner_data={"formed_bond_map_pairs": "9-10;11-12"},
    )

    path_kwargs = hunter.s2_engine.run_path_search.call_args.kwargs
    assert path_kwargs["forming_bonds"] == ((0, 2), (1, 3))
