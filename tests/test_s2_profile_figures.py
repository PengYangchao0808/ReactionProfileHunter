import json
from pathlib import Path

import pytest

from rph_core.utils.scan_profile_plotter import HARTREE_TO_KCAL, plot_scan_profile


pytest.importorskip("matplotlib")


def _energies(relative, baseline):
    return [baseline + value / HARTREE_TO_KCAL for value in relative]


def _profile(*, fallback_int: bool = True) -> dict:
    relative = [0.0, 1.0, 3.0, -1.0]
    selections = {"ts_guess": {"index": 2, "energy_source": "B97-3c"}}
    if fallback_int:
        selections["intermediate"] = {"index": 1, "selection_mode": "midpoint_fallback"}
    return {
        "profile_schema_version": "s2_scan_profile_v9",
        "xtb_path": {"path_arclength": [0.0, 0.1, 0.2, 0.3]},
        "energy_curves": {
            "xtb": {"energies_hartree": _energies(relative, -10.0), "relative_energies_kcal_mol": relative},
            "b973c": {"energies_hartree": _energies(relative, -100.0), "relative_energies_kcal_mol": relative},
        },
        "selection_policy": {"actual_source": "B97-3c"},
        "selections": selections,
        "trajectory_quality": {"off_path_indices": []},
    }


def test_profile_figure_is_the_single_root_png(tmp_path: Path):
    profile_path = tmp_path / "scan_profile.json"
    profile_path.write_text(json.dumps(_profile()), encoding="utf-8")

    output = plot_scan_profile(profile_path)

    assert output == tmp_path / "scan_profile.png"
    assert output.is_file()
    assert output.stat().st_size > 0
    assert not (tmp_path / "outputs").exists()


def test_profile_figure_allows_ts_only(tmp_path: Path):
    profile_path = tmp_path / "scan_profile.json"
    profile_path.write_text(json.dumps(_profile(fallback_int=False)), encoding="utf-8")

    output = plot_scan_profile(profile_path)

    assert output is not None
    assert output.is_file()
