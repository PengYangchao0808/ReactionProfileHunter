import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from rph_core.steps.step2_retro.energy_refinement import ScanEnergyRefiner
from rph_core.steps.step2_retro.peb_engine import PEBScanEngine
from rph_core.steps.step2_retro.relaxed_scan_rescue import B97CRelaxedScanRescuer
from rph_core.steps.step2_retro.scan_trajectory import ScanAttempt
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.utils.bond_pairs import canonicalize_bond_pairs
from rph_core.utils.qc_models import QCJobResult
from rph_core.utils.xtb_runner import XTBRunner
from rph_core.utils.scan_profile_plotter import compute_scan_distances


def _write_xyz(path: Path, distance: float = 1.0) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"2\nframe\nH 0 0 0\nH {distance:.6f} 0 0\n",
        encoding="utf-8",
    )
    return path


def _run_profile_with_b973c(tmp_path, monkeypatch, *, fail_int_segment=False):
    product = _write_xyz(tmp_path / "product.xyz", 1.5)
    engine = PEBScanEngine(
        {
            "step2": {
                "scan": {
                    "topology_guard_enabled": False,
                    "scan_start_distance": 3.4,
                    "scan_end_distance": 1.5,
                    "scan_steps": "auto",
                    "coarse_step_A": 0.20,
                    "candidate_refinement": {
                        "enabled": True,
                        "ts": {"half_window_A": 0.25, "step_A": 0.05},
                        "intermediate": {"half_window_A": 0.25, "step_A": 0.05},
                        "min_overlap_points": 2,
                    },
                    "endpoint_extension": {"enabled": False},
                    "selection": {
                        "preferred_energy_source": "b973c",
                    },
                },
                "xtb_path": {"enabled": False},
                "energy_refinement": {
                    "enabled": True,
                    "scope": "accepted_candidate_refinement_frames",
                },
            }
        }
    )
    def relative_energy(distance):
        return 44.0 / (1.0 + math.exp(-8.0 * (float(distance) - 2.15)))

    def fake_scan(*, output_dir, params, direction="outward", **_kwargs):
        coordinates = compute_scan_distances(
            params["scan_start_distance"],
            params["scan_end_distance"],
            params["scan_steps"],
            direction=direction,
        )
        frames = [
            _write_xyz(Path(output_dir) / "scan_frames" / f"frame_{index:03d}.xyz", q)
            for index, q in enumerate(coordinates)
        ]
        energies = [-10.0 + relative_energy(q) / 627.509 for q in coordinates]
        return SimpleNamespace(geometries=frames), energies

    monkeypatch.setattr(engine, "_execute_scan", fake_scan)

    def fake_refine(_self, requested_frames, point_ids=None):
        energies = []
        for frame in requested_frames:
            distance = float(Path(frame).read_text(encoding="utf-8").splitlines()[3].split()[1])
            value = -100.0 + relative_energy(distance) / 627.509
            energies.append(None if fail_int_segment and distance > 2.75 else value)
        return {
            "enabled": True,
            "status": "complete" if all(value is not None for value in energies) else "partial",
            "energies_hartree": energies,
            "records": [],
        }
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.ScanEnergyRefiner.refine",
        fake_refine,
    )
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.plot_scan_profile",
        lambda *_args, **_kwargs: None,
    )
    result = engine.run(product, tmp_path / "s2", [(0, 1)])
    return json.loads(Path(result[4]).read_text(encoding="utf-8"))


def test_bond_pairs_are_unordered_and_deduplicated():
    assert canonicalize_bond_pairs([(7, 6), (2, 3), (6, 7)]) == ((2, 3), (6, 7))


def test_only_a_fully_distorted_path_requires_relaxed_scan():
    assert PEBScanEngine._path_requires_relaxed_scan(
        {
            "checked": True,
            "persistent_off_path_start": 0,
            "usable_end_index": -1,
        }
    )
    assert PEBScanEngine._path_requires_relaxed_scan(
        {
            "checked": True,
            "topology_drift_index": 0,
            "last_valid_before_drift_index": None,
        }
    )
    assert not PEBScanEngine._path_requires_relaxed_scan(
        {
            "checked": True,
            "persistent_off_path_start": 3,
            "usable_end_index": 2,
        }
    )
    assert not PEBScanEngine._path_requires_relaxed_scan(
        {
            "checked": True,
            "persistent_off_path_start": None,
            "usable_end_index": 23,
        }
    )


def test_clean_path_dispatches_shared_ts_int_seed_to_s3(tmp_path, monkeypatch):
    profile = _run_profile_with_b973c(tmp_path, monkeypatch)

    assert profile["s3_dispatch"]["resolution"] == "path_seeded"
    assert profile["s3_dispatch"]["submit_ts"] is True
    assert profile["s3_dispatch"]["submit_intermediate"] is True
    assert profile["s3_dispatch"]["neb_eligible"] is False


def test_b973c_relaxed_scan_rescue_accepts_positive_barrier_direct_ts(
    tmp_path, monkeypatch
):
    product = _write_xyz(tmp_path / "product.xyz", 1.5)
    frames = [
        _write_xyz(
            tmp_path / "rescue_frames" / f"frame_{index:03d}.xyz",
            1.5 + index * 0.1,
        )
        for index in range(17)
    ]
    events = []
    captured = {}

    def fake_surface_scan(spec, *_args, **_kwargs):
        captured["spec"] = spec
        return QCJobResult(
            "complete",
            product,
            output_file=tmp_path / "rescue.out",
            extra={
                "frames": [str(frame) for frame in frames],
                # Product -> stretched-end ordering with a strict internal peak.
                # The reactant-side endpoint is the final frame, so the barrier
                # is evaluated against the tail of this ledger, not frame 0.
                "energies_hartree": [
                    -100.0 + value / 627.509
                    for value in [
                        1.0,
                        1.2,
                        1.4,
                        1.8,
                        2.0,
                        2.2,
                        2.5,
                        3.0,
                        5.5,
                        2.8,
                        1.5,
                        0.8,
                        0.4,
                        0.2,
                        0.1,
                        0.05,
                        0.0,
                    ]
                ],
                "energy_source": "orca.relaxscanact.dat",
                "scan_ts_candidate_xyz": str(frames[8]),
            },
        )

    monkeypatch.setattr(
        "rph_core.steps.step2_retro.relaxed_scan_rescue.run_surface_scan",
        fake_surface_scan,
    )
    rescuer = B97CRelaxedScanRescuer(
        {
            "resources": {"nproc": 7},
            "theory": {"s3_low_level": {"optimization": {
                "method": "B97-3c",
                "solvent": "acetone",
                "solvent_model": "CPCM",
                "max_cycles": 88,
                "timeout": 12345,
            }}},
            "step2": {"rescue": {"enabled": True, "relaxed_scan": {
                "points": 17,
                "stretch_end_A": 2.5,
                "ts_min_prominence_kcal_mol": 1.0,
                "ts_min_reactant_barrier_kcal_mol": 3.0,
            }}},
        },
        event_callback=lambda event, record: events.append((event, record)),
        variant="product_major",
    )

    result = rescuer.run(product, tmp_path / "rescue", [(0, 1)], trigger_reasons=["weak_peak"])

    assert result["status"] == "complete"
    assert result["resolution"] == "rescue_seeded"
    assert result["s2_state"] == "rescue_seeded"
    assert result["seed_evidence"] == "knee_shifted"
    ts_index = int(result["ts_search_seed"]["frame_index"])
    int_index = int(result["int_search_seed"]["frame_index"])
    knee_index = int(result["ts_search_seed"]["knee_frame_index"])
    assert ts_index > knee_index
    assert int_index > ts_index
    assert result["ts_search_seed"]["stationary_point_claimed"] is False
    assert result["int_search_seed"]["selection_mode"] == (
        "ts_to_effective_endpoint_midpoint"
    )
    rescue_root = tmp_path / "rescue"
    assert result["ts_xyz"] == str(rescue_root / frames[ts_index].name)
    assert result["intermediate_xyz"] == str(rescue_root / frames[int_index].name)
    assert "direct_ts_only" not in json.dumps(result, sort_keys=True)
    rescue_profile = Path(result["scan_profile"])
    assert rescue_profile.is_file()
    payload = json.loads(rescue_profile.read_text(encoding="utf-8"))
    assert payload["profile_kind"] == "relaxed_scan_rescue"
    assert payload["selection_policy"]["actual_source"] == "B97-3c"
    assert payload["rescue"]["resolution"] == "rescue_seeded"
    assert payload["rescue"]["s2_state"] == "rescue_seeded"
    assert payload["rescue"]["seed_evidence"] == "knee_shifted"
    assert result["scan_plot"] is None
    assert not (rescue_root / "scan_profile.png").exists()
    assert not (rescue_root / "scan_frames").exists()
    assert captured["spec"].nproc == 7
    assert captured["spec"].max_cycles == 88
    assert captured["spec"].timeout == 12345
    assert captured["spec"].solvent == "acetone"
    assert [event for event, _record in events] == [
        "step_started",
        "batch_started",
        "batch_job_started",
        "batch_job_finished",
        "batch_finished",
        "step_finished",
    ]
    started = events[0][1]
    assert started["label"] == "ORCA B97-3c relaxed-scan rescue"
    assert started["scan_mode"] == "ScanTS"
    assert started["point_count"] == 17
    assert started["nprocs"] == 7
    assert all(record["variant"] == "product_major" for _event, record in events)


def test_orca_relaxed_scan_energy_ledger_uses_optimized_point_energies(tmp_path):
    ledger = tmp_path / "scan.relaxscanact.dat"
    ledger.write_text(
        "1.50 1.55 -100.00000000\n"
        "1.80 1.85 -99.95000000\n"
        "2.10 2.15 -100.02000000\n",
        encoding="utf-8",
    )

    assert ORCAInterface._parse_relaxed_scan_energy_ledger(
        ledger,
        coordinate_count=2,
    ) == [-100.0, -99.95, -100.02]


def test_b973c_monotonic_profile_uses_refined_curve_fallback(tmp_path, monkeypatch):
    profile = _run_profile_with_b973c(
        tmp_path,
        monkeypatch,
    )

    assert profile["selection_policy"]["actual_source"] == "B97-3c"
    assert profile["selection_policy"]["xtb_fallback_used"] is False
    ts_index = int(profile["selections"]["ts_guess"]["index"])
    assert 0 <= ts_index < len(profile["reaction_coordinate_angstrom"])
    assert profile["selections"]["ts_guess"]["rule"] in {
        "refined_curve_local_maximum",
        "refined_curve_weak_local_maximum",
        "refined_curve_internal_maximum_fallback",
        "anchor_fallback_no_curve_extremum",
        "reactant_side_backoff_from_refined_curve_local_maximum",
        "reactant_side_backoff_from_refined_curve_weak_local_maximum",
        "reactant_side_backoff_from_refined_curve_internal_maximum_fallback",
        "reactant_side_backoff_from_anchor_fallback_no_curve_extremum",
        "unified_selector_monotonic_shoulder",
        "unified_selector_knee_shifted_ts",
    }
    assert profile["energy_refinement"]["scope"] == "path_full_coverage"


def test_refined_curve_selects_pre_ts_basin_before_midpoint_fallback(tmp_path):
    engine = PEBScanEngine({"step2": {}})
    frames = [_write_xyz(tmp_path / f"frame_{index:03d}.xyz") for index in range(6)]
    energies = [-100.0 + value / 627.509 for value in [0.0, -1.0, 1.0, 3.0, 0.5, -1.0]]

    ts, intermediate, ts_index, int_index, status = engine._select_refined_path_nodes(
        anchors={},
        frame_paths=frames,
        reaction_coordinate=[3.0, 2.8, 2.6, 2.4, 2.2, 2.0],
        method_energies=energies,
        off_path_indices=[],
        path_arclength=np.arange(6, dtype=float),
        selection_config={
            "ts_min_prominence_kcal_mol": 0.15,
            "int_min_basin_prominence_kcal_mol": 0.50,
        },
    )

    assert ts_index == 3
    assert ts["rule"] == "refined_curve_local_maximum"
    assert int_index == 1
    assert intermediate["selection_mode"] == "stable_basin_candidate"
    assert intermediate["stationary_point_claimed"] is False
    assert status == "selected"


def test_ts_seed_reactant_side_backoff_preserves_energy_peak():
    seed_index, metadata = PEBScanEngine._select_reactant_side_ts_seed(
        peak_index=3,
        forming_bond_distances=[
            (3.0, 3.0),
            (2.8, 2.8),
            (2.6, 2.6),
            (2.4, 2.4),
            (2.2, 2.2),
        ],
        invalid_indices=[],
        backoff_A=0.20,
    )

    assert seed_index == 2
    assert metadata["energy_peak_index"] == 3
    assert metadata["seed_index"] == 2
    assert metadata["seed_backoff_applied_A"] == pytest.approx(0.20)
    assert metadata["seed_rule"] == "reactant_side_mean_distance_backoff"


def test_ts_seed_boundary_fallback_never_uses_product_side_frames():
    seed_index, metadata = PEBScanEngine._select_reactant_side_ts_seed(
        peak_index=1,
        forming_bond_distances=[
            (1.47, 1.47),
            (1.48, 1.48),
            (1.50, 1.50),
            (1.54, 1.54),
        ],
        invalid_indices=[],
        backoff_A=0.20,
        path_arclength=[0.0, 0.1, 0.2, 0.3],
    )

    assert seed_index == 0
    assert seed_index < 1
    assert metadata["seed_rule"] == "path_order_boundary_fallback"
    assert metadata["seed_backoff_status"] == "no_bond_consistent_pre_peak_frame"


def test_product_connected_valid_indices_excludes_disconnected_path_prefix():
    assert PEBScanEngine._product_connected_valid_indices(8, [0, 1, 5]) == [6, 7]
    assert PEBScanEngine._product_connected_valid_indices(5, []) == [0, 1, 2, 3, 4]
    assert PEBScanEngine._product_connected_valid_indices(4, [3]) == []


def test_full_endpoint_path_is_retained_but_selection_excludes_its_distorted_prefix(
    tmp_path, monkeypatch
):
    product = _write_xyz(tmp_path / "product.xyz", 1.5)
    engine = PEBScanEngine(
        {
            "step2": {
                "scan": {
                    "topology_guard_enabled": True,
                    "scan_start_distance": 3.4,
                    "scan_end_distance": 1.5,
                    "scan_steps": 6,
                    "endpoint_extension": {"enabled": False},
                    "selection": {"preferred_energy_source": "b973c"},
                },
                "xtb_path": {
                    "enabled": True,
                    "endpoint_strategy": "full_endpoint_with_valid_corridor",
                    "retain_valid_corridor_path": True,
                },
            }
        }
    )
    coarse_frames = [
        _write_xyz(tmp_path / "coarse" / f"frame_{index:03d}.xyz", 1.5 + index * 0.2)
        for index in range(6)
    ]
    monkeypatch.setattr(
        engine,
        "_execute_scan",
        lambda **_kwargs: (SimpleNamespace(geometries=coarse_frames), [-10.0, -9.9, -9.8, -9.7, -9.6, -9.5]),
    )
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.check_scan_trajectory",
        lambda **_kwargs: {"checked": True, "total_frames": 6, "off_path_indices": [4, 5], "off_path_count": 2},
    )

    def fake_path(*, path_name, output_dir, **_kwargs):
        distances = (
            [3.4, 3.0, 2.7, 2.4, 2.1, 1.8, 1.5]
            if path_name == "full_endpoint"
            else [2.7, 2.4, 2.1, 1.8, 1.5]
        )
        frames = tuple(
            _write_xyz(Path(output_dir) / f"{path_name}_{index:03d}.xyz", distance)
            for index, distance in enumerate(distances)
        )
        energies = tuple(-100.0 + value / 627.509 for value in range(len(distances)))
        off_path = (0, 1) if path_name == "full_endpoint" else ()
        return (
            ScanAttempt(
                attempt_id=f"NN_xtb_path_{path_name}", kind="xtb_path", directory=Path(output_dir),
                frame_paths=frames, target_coordinates_A=tuple(distances),
                xtb_energies_hartree=energies, off_path_indices=off_path,
                trajectory_quality={"off_path_indices": list(off_path)}, scan_policy="xtb_path",
            ),
            {"path_name": path_name, "path_arclength": list(np.arange(len(frames), dtype=float))},
        )

    monkeypatch.setattr(engine, "_execute_xtb_path", fake_path)

    def fake_refine(_self, frames, point_ids=None):
        branch = "valid_corridor" if "valid_corridor" in Path(frames[0]).name else "full_endpoint"
        relative = (
            [0.0, 0.0, 0.0, 1.0, 3.0, 1.0, 0.0]
            if branch == "full_endpoint"
            else [0.0, 1.0, 4.0, 1.0, 0.0]
        )
        return {
            "status": "complete",
            "energies_hartree": [-100.0 + value / 627.509 for value in relative],
            "records": [],
        }

    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.ScanEnergyRefiner.refine",
        fake_refine,
    )
    monkeypatch.setattr("rph_core.steps.step2_retro.peb_engine.plot_scan_profile", lambda *_args, **_kwargs: None)

    result = engine.run(product, tmp_path / "s2", [(0, 1)])
    profile = json.loads(Path(result[4]).read_text(encoding="utf-8"))

    assert set(profile["path_branches"]) == {"full_endpoint", "valid_corridor"}
    assert profile["path_branches"]["valid_corridor"]["selected_for_s3"] is True
    assert profile["path_branches"]["full_endpoint"]["trajectory_quality"]["selection_excluded_indices"] == [0, 1]
    assert profile["selection_decision"]["rule"] == "unified_selector_unresolved"
    # The unified selector applies the endpoint and neighbor gates to the
    # valid corridor as well.  This short corridor has no admissible TS seed;
    # it must remain unresolved rather than revive the legacy boundary seed.
    assert profile["s2_state"] == "unresolved"
    assert profile["selection_source"] is None
    assert profile["s3_dispatch"]["resolution"] == "unresolved"
    assert profile["s3_dispatch"]["submit_ts"] is False
    assert profile["selections"]["ts_guess"]["index"] is None
    assert profile["selections"]["intermediate"]["index"] is None
    assert result[0] is None
    assert result[1] is None
    assert result[2] is None


def test_refined_curve_marks_unresolved_int_as_shared_ts_seed(tmp_path):
    engine = PEBScanEngine({"step2": {}})
    frames = [_write_xyz(tmp_path / f"frame_{index:03d}.xyz") for index in range(6)]
    energies = [-100.0 + value / 627.509 for value in [0.0, 0.5, 1.0, 1.5, 3.0, 0.0]]

    _, intermediate, ts_index, int_index, status = engine._select_refined_path_nodes(
        anchors={},
        frame_paths=frames,
        reaction_coordinate=[3.0, 2.8, 2.6, 2.4, 2.2, 2.0],
        method_energies=energies,
        off_path_indices=[],
        path_arclength=np.arange(6, dtype=float),
        selection_config={
            "ts_min_prominence_kcal_mol": 0.15,
            "int_min_basin_prominence_kcal_mol": 0.50,
        },
    )

    assert ts_index == 4
    assert int_index is None
    assert intermediate["index"] == ts_index
    assert intermediate["selection_mode"] == "shared_ts_fallback"
    assert intermediate["rule"] == "no_resolved_pre_ts_basin_shared_ts_seed"
    assert intermediate["reason"] == "no_resolved_pre_ts_basin"
    assert intermediate["s3_job_required"] is False
    assert status == "shared_with_ts"


def test_late_pre_ts_platform_seed_is_preferred_to_shared_ts(tmp_path):
    engine = PEBScanEngine({"step2": {}})
    frames = [_write_xyz(tmp_path / f"frame_{index:03d}.xyz") for index in range(7)]
    energies = [
        -100.0 + value / 627.509
        for value in [0.0, 0.1, 0.2, 0.3, 0.5, 4.0, 0.0]
    ]

    index, selection = engine._select_late_pre_ts_platform_seed(
        frame_paths=frames,
        reaction_coordinate=[0.1 * value for value in range(7)],
        method_energies=energies,
        off_path_indices=[],
        path_arclength=np.asarray([0.1 * value for value in range(7)]),
        ts_peak_index=5,
        ts_seed_index=5,
        selection_config={
            "int_plateau_fallback_enabled": True,
            "int_plateau_min_consecutive_frames": 3,
            "int_plateau_min_ts_separation_A": 0.10,
            "int_plateau_energy_window_kcal_mol": 2.0,
            "int_plateau_barrier_fraction": 0.25,
            "int_plateau_max_slope_kcal_mol_A": 15.0,
        },
    )

    assert index == 4
    assert selection["selection_mode"] == "late_pre_ts_platform_fallback"
    assert selection["s3_job_required"] is True
    assert selection["stationary_point_claimed"] is False


def test_b973c_failure_in_m_d_domain_uses_method_consistent_xtb_fallback(tmp_path, monkeypatch):
    profile = _run_profile_with_b973c(
        tmp_path,
        monkeypatch,
        fail_int_segment=True,
    )

    # With a unified refinement corridor, a B97-3c SP failure inside the
    # corridor cascades to both TS and INT selections because they share the
    # same refinement frames.
    assert profile["energy_refinement"]["status"] == "partial"
    assert profile["selection_policy"]["actual_source"] == "xtb_fallback"
    assert profile["selection_policy"]["intermediate_source"] == "xtb_fallback"
    assert profile["selection_policy"]["xtb_fallback_used"] is True


def test_geometry_guard_can_reclassify_only_an_isolated_low_jump_suspect(
    tmp_path, monkeypatch
):
    product = _write_xyz(tmp_path / "product.xyz", 1.5)
    frames = [
        _write_xyz(tmp_path / f"frame_{index:03d}.xyz", 1.5 + index * 0.2)
        for index in range(6)
    ]
    relative = [-2.0, -1.0, 5.0, 4.0, 3.0, 2.0]
    energies = [-10.0 + value / 627.509 for value in relative]
    engine = PEBScanEngine(
        {
            "step2": {
                "scan": {
                    "topology_guard_enabled": True,
                    "reclassify_isolated_topology_suspects": True,
                    "isolated_topology_max_energy_jump_kcal": 5.0,
                    "isolated_topology_max_reaction_core_rmsd_A": 0.25,
                    "scan_start_distance": 3.4,
                    "scan_end_distance": 1.5,
                    "scan_steps": 6,
                    "candidate_refinement": {"enabled": False},
                    "endpoint_extension": {"enabled": False},
                    "selection": {
                        "preferred_energy_source": "xtb",
                    },
                },
                "xtb_path": {"enabled": False},
                "energy_refinement": {"enabled": False},
            }
        }
    )
    monkeypatch.setattr(
        engine,
        "_execute_scan",
        lambda **_kwargs: (SimpleNamespace(geometries=frames), energies),
    )
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.check_scan_trajectory",
        lambda **_kwargs: {
            "checked": True,
            "total_frames": 6,
            "off_path_indices": [3],
            "off_path_count": 1,
            "frame_issues": [{"frame_index": 3}],
        },
    )
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.plot_scan_profile",
        lambda *_args, **_kwargs: None,
    )

    result = engine.run(product, tmp_path / "s2", [(0, 1)])
    profile = json.loads(Path(result[4]).read_text(encoding="utf-8"))
    quality = profile["attempts"][0]["trajectory_quality"]
    assert quality["raw_off_path_indices"] == [3]
    assert quality["isolated_suspect_indices"] == [3]
    assert quality["reclassified_isolated_indices"] == [3]
    assert quality["off_path_indices"] == []


def test_xtb_scan_keeps_fixed_constraints_out_of_scan_block(tmp_path):
    runner = object.__new__(XTBRunner)
    runner.config = {}
    runner.work_dir = tmp_path

    scan_input = runner._write_scan_input(
        constraints={"2 3": 1.5, "6 7": 1.5},
        fixed_constraints={"0 9": 3.0},
        scan_range=(1.5, 3.4),
        scan_steps=25,
        scan_force_constant=0.05,
    )
    text = scan_input.read_text(encoding="utf-8")

    assert "distance: 1, 10, 3.000" in text
    scan_block = text.split("$scan", 1)[1].split("$opt", 1)[0]
    assert "1: 1.500, 3.400, 25" in scan_block
    assert "2: 1.500, 3.400, 25" in scan_block
    assert "3:" not in scan_block


def test_frame_energy_refinement_is_cached_per_geometry(tmp_path, monkeypatch):
    frames = [_write_xyz(tmp_path / f"frame_{index:03d}.xyz", 1.0 + index) for index in range(3)]
    calls = []

    def fake_single_point(spec, input_xyz, output_dir, config):
        calls.append(Path(input_xyz).name)
        output = Path(output_dir) / "job.out"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("complete", encoding="utf-8")
        index = int(Path(input_xyz).stem.split("_")[-1])
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_file=output,
            energy_hartree=-100.0 + index,
        )

    monkeypatch.setattr(
        "rph_core.steps.step2_retro.energy_refinement.run_single_point",
        fake_single_point,
    )
    config = {
        "resources": {"nproc": 2, "mem": "2GB"},
        "step2": {
            "energy_refinement": {
                "enabled": True,
                "parallel_jobs": 2,
                "cores_per_job": 1,
                "method": "B97-3c",
                "solvent": "acetone",
            }
        },
    }
    refiner = ScanEnergyRefiner(config, tmp_path / "sp")

    first = refiner.refine(frames)
    second = refiner.refine(frames)
    reordered = refiner.refine(list(reversed(frames)))

    assert first["status"] == "complete"
    assert first["energies_hartree"] == [-100.0, -99.0, -98.0]
    assert len(calls) == 3
    assert all(record["reused"] for record in second["records"])
    assert reordered["energies_hartree"] == [-98.0, -99.0, -100.0]
    assert all(record["reused"] for record in reordered["records"])


def test_frame_energy_refinement_reports_dynamic_batch_and_orca_jobs(tmp_path, monkeypatch):
    frames = [_write_xyz(tmp_path / f"frame_{index:03d}.xyz", 1.0 + index) for index in range(3)]
    events = []

    def fake_single_point(spec, input_xyz, output_dir, config):
        output = Path(output_dir) / "job.out"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("complete", encoding="utf-8")
        index = int(Path(input_xyz).stem.split("_")[-1])
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_file=output,
            energy_hartree=-100.0 + index,
        )

    monkeypatch.setattr(
        "rph_core.steps.step2_retro.energy_refinement.run_single_point",
        fake_single_point,
    )
    config = {
        "resources": {"nproc": 2, "mem": "4GB"},
        "step2": {
            "energy_refinement": {
                "enabled": True,
                "parallel_jobs": 2,
                "cores_per_job": 1,
                "method": "B97-3c",
            }
        },
    }
    result = ScanEnergyRefiner(
        config,
        tmp_path / "sp",
        variant="product_major",
        event_callback=lambda event, payload: events.append((event, payload)),
    ).refine(frames, point_ids=["p_000", "p_001", "p_002"])

    started = next(payload for event, payload in events if event == "batch_started")
    finished = next(payload for event, payload in events if event == "batch_finished")
    job_starts = [payload for event, payload in events if event == "batch_job_started"]
    progresses = [payload for event, payload in events if event == "batch_progress"]
    assert started["batch"] == "product_major:b973c_sp"
    assert started["label"] == "ORCA B97-3c SP refinement"
    assert started["total"] == 3
    assert started["parallel_jobs"] == 2
    assert len(job_starts) == 3
    assert all(payload["engine"] == "orca" for payload in job_starts)
    assert progresses[-1]["done"] == 3
    assert progresses[-1]["total"] == 3
    assert finished["done"] == 3
    assert finished["status"] == "complete"
    assert result["status"] == "complete"


def test_geometry_cache_migrates_matching_legacy_frame_cache(tmp_path, monkeypatch):
    frame = _write_xyz(tmp_path / "source.xyz", 1.75)
    config = {
        "resources": {"nproc": 1, "mem": "2GB"},
        "step2": {
            "energy_refinement": {
                "enabled": True,
                "parallel_jobs": 1,
                "cores_per_job": 1,
                "method": "B97-3c",
                "solvent": "acetone",
            }
        },
    }
    output_dir = tmp_path / "sp"
    refiner = ScanEnergyRefiner(config, output_dir)
    spec_signature = refiner._signature(refiner._spec())
    legacy_path = output_dir / "frame_007" / "result.json"
    legacy_path.parent.mkdir(parents=True)
    legacy_path.write_text(
        json.dumps(
            {
                "schema_version": "s2_frame_sp_cache_v1",
                "frame_index": 7,
                "input_xyz": str(frame),
                "xyz_sha256": refiner._file_hash(frame),
                "spec_sha256": spec_signature,
                "status": "complete",
                "energy_hartree": -55.25,
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.energy_refinement.run_single_point",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("QC reran")),
    )

    result = refiner.refine([frame], point_ids=["p_0000"])

    assert result["energies_hartree"] == [-55.25]
    assert result["records"][0]["reused"] is True
    assert result["records"][0]["migrated_from"] == str(legacy_path)
    assert list((output_dir / "cache").glob("*/result.json"))


def test_failed_intermediate_selection_archives_stale_published_nodes(tmp_path, monkeypatch):
    product = _write_xyz(tmp_path / "product.xyz")
    output_dir = tmp_path / "s2"
    output_dir.mkdir()
    for name in ("manifest.json", "intermediate.xyz", "reactant_complex.xyz"):
        (output_dir / name).write_text("stale", encoding="utf-8")
    frames = [_write_xyz(tmp_path / f"scan_{index:03d}.xyz", 1.0 + index * 0.1) for index in range(6)]
    energies = [-10.0, -9.9, -9.0, -9.2, -9.4, -9.6]
    engine = PEBScanEngine(
        {
            "step2": {
                "scan": {
                    "topology_guard_enabled": True,
                    "scan_start_distance": 3.4,
                    "scan_end_distance": 1.5,
                    "scan_steps": 6,
                    "selection": {
                        "preferred_energy_source": "xtb",
                    },
                },
                "xtb_path": {"enabled": False},
                "energy_refinement": {"enabled": False},
            }
        }
    )

    monkeypatch.setattr(
        engine,
        "_execute_scan",
        lambda **_kwargs: (SimpleNamespace(geometries=frames), energies),
    )
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.check_scan_trajectory",
        lambda **_kwargs: {
            "checked": True,
            "total_frames": 6,
            "off_path_indices": [3, 4, 5],
            "off_path_count": 3,
            "frame_issues": [
                {"frame_index": 3},
                {"frame_index": 4},
                {"frame_index": 5},
            ],
        },
    )
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.plot_scan_profile",
        lambda *_args, **_kwargs: None,
    )

    result = engine.run(product, output_dir, [(0, 1)])
    profile = json.loads(Path(result[4]).read_text(encoding="utf-8"))

    # The unified selector does not fabricate a seed from a short prefix that
    # cannot satisfy its neighborhood and endpoint gates.
    assert result[0] is None
    assert result[1] is None
    assert result[2] is None
    assert result[5] == "DEGRADED"
    assert not (output_dir / "ts_guess.xyz").exists()
    assert not (output_dir / "intermediate.xyz").exists()
    assert not (output_dir / "reactant_complex.xyz").exists()
    # PEBEngine owns scan_profile.json; the orchestrator publishes manifest.json.
    assert not (output_dir / "manifest.json").exists()
    assert profile["artifact_lifecycle"]["legacy_reactant_complex_written"] is False
    assert profile["profile_schema_version"] == "s2_scan_profile_v10"
    assert profile["attempts"][0]["kind"] == "coarse"
    assert profile["composite_profile"]["coverage"]["coordinate_min_A"] == 1.5
    assert profile["composite_profile"]["coverage"]["coordinate_max_A"] == 3.4
    assert profile["energy_curves"]["xtb"]["status"] == "complete"
    assert profile["s2_state"] == "unresolved"
    assert profile["s3_dispatch"]["resolution"] == "unresolved"
    assert profile["s3_dispatch"]["submit_ts"] is False
    assert profile["s3_dispatch"]["submit_intermediate"] is False
    assert list((output_dir / "superseded").glob("*/manifest.json"))
