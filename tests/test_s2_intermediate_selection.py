import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from rph_core.steps.step2_retro.peb_engine import PEBScanEngine
from rph_core.utils.scan_profile_plotter import HARTREE_TO_KCAL, plot_scan_profile
from rph_core.v4_orchestrator import V4Orchestrator


def _hartree_profile(relative_kcal):
    return [float(value) / HARTREE_TO_KCAL for value in relative_kcal]


def _selection_config(**overrides):
    config = {
        "method": "stable_midpoint",
        "endpoint_exclusion_points": 1,
        "stable_quantile": 0.5,
        "min_candidate_points": 2,
        "require_topology_valid": True,
        "allow_degraded_fallback": False,
        "endpoint_jump_warning_ratio": 5.0,
        "endpoint_jump_warning_min_kcal": 5.0,
    }
    config.update(overrides)
    return config


def test_stable_midpoint_selects_major_frame_21():
    energies = _hartree_profile(
        [-20.0] * 19 + [0.0, -2.051, -2.289, -2.459, -61.100]
    )
    frames = [Path(f"frame_{index:03d}.xyz") for index in range(len(energies))]

    selection = PEBScanEngine._select_intermediate_frame(
        energies,
        frames,
        peak_index=19,
        off_path_indices=[],
        selection_config=_selection_config(),
    )

    assert selection["index"] == 21
    assert selection["endpoint_index"] == 23
    assert 23 not in selection["candidate_indices"]


def test_stable_midpoint_selects_minor_frame_18():
    energies = _hartree_profile(
        [-25.0] * 14
        + [0.0, -32.434, -35.798, -37.443, -37.995, -37.782, -36.847, -35.248, -32.860, -29.740]
    )
    frames = [Path(f"frame_{index:03d}.xyz") for index in range(len(energies))]

    selection = PEBScanEngine._select_intermediate_frame(
        energies,
        frames,
        peak_index=14,
        off_path_indices=[],
        selection_config=_selection_config(),
    )

    assert selection["index"] == 18
    assert selection["midpoint_index"] == pytest.approx(18.5)


def test_stable_midpoint_rejects_off_path_candidates():
    energies = _hartree_profile([-2.0, -1.0, 3.0, 1.0, 0.9, -5.0])
    frames = [Path(f"frame_{index:03d}.xyz") for index in range(len(energies))]

    selection = PEBScanEngine._select_intermediate_frame(
        energies,
        frames,
        peak_index=2,
        off_path_indices=[3],
        selection_config=_selection_config(min_candidate_points=1),
    )

    assert selection["index"] == 4
    assert 3 not in selection["stable_indices"]
    assert selection["index"] != len(energies) - 1


def test_stable_midpoint_fails_without_valid_interior_frame():
    energies = _hartree_profile([-2.0, -1.0, 3.0, 1.0, 0.9, -5.0])
    frames = [Path(f"frame_{index:03d}.xyz") for index in range(len(energies))]

    with pytest.raises(RuntimeError, match="no topology-valid interior frames"):
        PEBScanEngine._select_intermediate_frame(
            energies,
            frames,
            peak_index=2,
            off_path_indices=[3, 4],
            selection_config=_selection_config(),
        )


def test_peb_run_copies_selected_scan_frame_and_survives_plot_failure(tmp_path, monkeypatch):
    product = tmp_path / "product.xyz"
    product.write_text("2\nproduct\nH 0 0 0\nH 1 0 0\n", encoding="utf-8")
    frames = []
    for index in range(6):
        frame = tmp_path / f"frame_{index:03d}.xyz"
        frame.write_text(
            f"2\nframe {index}\nH 0 0 0\nH {1.0 + 0.2 * index:.2f} 0 0\n",
            encoding="utf-8",
        )
        frames.append(frame)

    energies = _hartree_profile([-2.0, -1.0, 5.0, 2.0, 1.9, -10.0])
    engine = PEBScanEngine(
        {
            "step2": {
                "scan": {
                    "topology_guard_enabled": False,
                    "scan_start_distance": 4.0,
                    "scan_end_distance": 1.5,
                    "scan_steps": 6,
                    "intermediate_selection": _selection_config(),
                }
            }
        }
    )

    def fake_execute_scan(**_kwargs):
        result = SimpleNamespace(geometries=frames, ts_guess_xyz=frames[2])
        return result, energies, 2, False, True

    monkeypatch.setattr(engine, "_execute_scan", fake_execute_scan)
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.plot_scan_profile",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("plot failed")),
    )

    result = engine.run(product, tmp_path / "s2", [(0, 1)])
    ts_guess, _, intermediate, _, profile_path = result[:5]
    profile = json.loads(profile_path.read_text(encoding="utf-8"))

    assert ts_guess.read_text(encoding="utf-8") == frames[2].read_text(encoding="utf-8")
    assert intermediate.read_text(encoding="utf-8") == frames[3].read_text(encoding="utf-8")
    assert profile["selections"]["intermediate"]["index"] == 3
    assert profile["trajectory_quality"]["endpoint_excluded"] is True
    assert profile["scan_plot"] is None


def test_plot_scan_profile_marks_ts_intermediate_and_endpoint(tmp_path):
    pytest.importorskip("matplotlib")
    profile_path = tmp_path / "scan_profile.json"
    profile_path.write_text(
        json.dumps(
            {
                "energies_hartree": [-10.0, -9.99, -9.98, -9.985, -9.987, -10.01],
                "reaction_coordinate_angstrom": [1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
                "forming_bonds": [[0, 1]],
                "scan_parameters": {
                    "scan_start_distance": 4.0,
                    "scan_end_distance": 1.5,
                    "scan_steps": 6,
                },
                "selections": {
                    "ts_guess": {"index": 2},
                    "intermediate": {"index": 4},
                },
                "trajectory_quality": {
                    "endpoint_index": 5,
                    "endpoint_excluded": True,
                },
            }
        ),
        encoding="utf-8",
    )

    plot_path = plot_scan_profile(profile_path)

    assert plot_path == profile_path.with_suffix(".png")
    assert plot_path.exists()
    assert plot_path.stat().st_size > 0


def test_s2_manifest_archives_plot_and_selection_metadata(tmp_path):
    ts_guess = tmp_path / "ts_guess.xyz"
    intermediate = tmp_path / "intermediate.xyz"
    plot = tmp_path / "scan_profile.png"
    for path in (ts_guess, intermediate, plot):
        path.write_text("artifact", encoding="utf-8")
    profile = tmp_path / "scan_profile.json"
    profile.write_text(
        json.dumps(
            {
                "scan_plot": str(plot),
                "scan_quality": {"intermediate_confidence": "high"},
                "selections": {
                    "ts_guess": {"index": 19},
                    "intermediate": {"index": 21, "rule": "stable_midpoint"},
                },
            }
        ),
        encoding="utf-8",
    )

    manifest_path = V4Orchestrator._write_s2(
        tmp_path,
        ts_guess,
        intermediate,
        [(11, 18), (14, 15)],
        profile,
        "COMPLETE",
        "high",
        [],
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert manifest["schema_version"] == "s2_peb_manifest_v2"
    assert manifest["scan_plot"] == str(plot)
    assert manifest["ts_guess_index"] == 19
    assert manifest["intermediate_index"] == 21
    assert manifest["intermediate_selection_method"] == "stable_midpoint"
