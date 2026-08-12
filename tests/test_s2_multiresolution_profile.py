import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from rph_core.steps.step2_retro.peb_engine import PEBScanEngine
from rph_core.steps.step2_retro.scan_trajectory import (
    CompositeProfileBuilder,
    ScanAttempt,
)
from rph_core.utils.qc_models import QCJobResult
from rph_core.utils.scan_profile_plotter import plot_scan_profile
from rph_core.utils.scan_profile_plotter import compute_scan_distances


def _xyz(path: Path, distance: float) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"2\nframe\nH 0 0 0\nH {distance:.6f} 0 0\n",
        encoding="utf-8",
    )
    return path


def _attempt(
    tmp_path: Path,
    attempt_id: str,
    kind: str,
    coordinates,
    energies,
    *,
    off_path=(),
) -> ScanAttempt:
    frames = tuple(
        _xyz(tmp_path / attempt_id / f"frame_{index:03d}.xyz", coordinate)
        for index, coordinate in enumerate(coordinates)
    )
    return ScanAttempt(
        attempt_id=attempt_id,
        kind=kind,
        directory=tmp_path / attempt_id,
        frame_paths=frames,
        target_coordinates_A=tuple(coordinates),
        xtb_energies_hartree=tuple(energies),
        off_path_indices=tuple(off_path),
    )


def _hartree_profile(relative_kcal, reference_hartree=-100.0):
    return [reference_hartree + float(value) / 627.509 for value in relative_kcal]


def test_composite_profile_preserves_full_coarse_coverage_and_inserts_refinement(tmp_path):
    coarse = _attempt(
        tmp_path,
        "00_coarse",
        "coarse",
        [1.5, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5],
        [-10.0, -9.9, -9.6, -9.2, -8.95, -8.8, -9.1, -9.4],
    )
    refined = _attempt(
        tmp_path,
        "01_ts_refinement",
        "ts_refinement",
        [2.5, 2.6, 2.75, 2.9, 3.0],
        [-9.21, -9.08, -8.96, -8.88, -8.81],
    )

    profile = CompositeProfileBuilder(
        coordinate_tolerance_A=0.01,
        min_overlap_points=3,
    ).build([coarse, refined])

    coordinates = [point["target_coordinate_A"] for point in profile["points"]]
    assert min(coordinates) == 1.5
    assert max(coordinates) == 3.5
    assert "00_coarse" in profile["accepted_attempt_ids"]
    assert "01_ts_refinement" in profile["accepted_attempt_ids"]
    assert any(point["source_attempt"] == "01_ts_refinement" for point in profile["points"])
    assert profile["coverage"]["complete_xtb_curve"] is True


def test_composite_uses_reaction_core_rmsd_and_keeps_global_rmsd_as_diagnostic(
    tmp_path, monkeypatch
):
    coarse = _attempt(
        tmp_path,
        "00_coarse",
        "coarse",
        [1.5, 2.0, 2.5, 3.0],
        [-10.0, -9.8, -9.5, -9.7],
    )
    refined = _attempt(
        tmp_path,
        "01_ts_refinement",
        "ts_refinement",
        [2.0, 2.25, 2.5, 2.75, 3.0],
        [-9.8, -9.65, -9.5, -9.6, -9.7],
    )
    monkeypatch.setattr(
        CompositeProfileBuilder,
        "_aligned_rmsd",
        staticmethod(
            lambda _left, _right, atom_indices=None: 0.10 if atom_indices is not None else 1.20
        ),
    )

    profile = CompositeProfileBuilder(
        min_overlap_points=3,
        max_overlap_rmsd_A=0.75,
        max_reaction_core_rmsd_A=0.50,
        forming_bonds=[(0, 1)],
    ).build([coarse, refined])

    check = profile["continuity_checks"][0]
    assert check["accepted"] is True
    assert check["inserted_into_composite"] is True
    assert check["global_rmsd_warning"] is True
    assert check["p95_reaction_core_rmsd_A"] == 0.10


def test_composite_rejects_alternate_conformer_energy_offset(tmp_path, monkeypatch):
    coarse = _attempt(
        tmp_path,
        "00_coarse",
        "coarse",
        [1.5, 2.0, 2.5, 3.0],
        [-10.0, -9.8, -9.5, -9.7],
    )
    refined = _attempt(
        tmp_path,
        "01_ts_refinement",
        "ts_refinement",
        [2.0, 2.25, 2.5, 2.75, 3.0],
        [-9.70, -9.55, -9.40, -9.50, -9.60],
    )
    monkeypatch.setattr(
        CompositeProfileBuilder,
        "_aligned_rmsd",
        staticmethod(lambda *_args, **_kwargs: 0.10),
    )

    profile = CompositeProfileBuilder(
        min_overlap_points=3,
        max_overlap_energy_gap_kcal=25.0,
        forming_bonds=[(0, 1)],
    ).build([coarse, refined])

    check = profile["continuity_checks"][0]
    assert check["accepted"] is False
    assert "overlap_energy_offset_exceeded" in check["rejection_reasons"]
    assert check["alternate_path_candidate"] is True


def test_composite_profile_retains_topology_drift_for_diagnostics(tmp_path):
    coarse = _attempt(
        tmp_path,
        "00_coarse",
        "coarse",
        [1.5, 2.0, 2.5, 3.0],
        [-10.0, -9.0, -8.0, -9.0],
        off_path=(3,),
    )

    profile = CompositeProfileBuilder().build([coarse])

    assert profile["coverage"]["point_count"] == 4
    assert profile["coverage"]["topology_valid_point_count"] == 3
    assert profile["points"][3]["topology_valid"] is False


def test_dual_curve_plot_uses_full_xtb_profile(tmp_path):
    profile = tmp_path / "scan_profile.json"
    profile.write_text(
        json.dumps(
            {
                "profile_schema_version": "s2_scan_profile_v4",
                "reaction_coordinate_angstrom": [1.5, 2.0, 2.5, 3.0, 3.5],
                "energy_curves": {
                    "xtb": {
                        "energies_hartree": [-10.0, -9.9, -9.0, -9.4, -9.8],
                        "relative_energies_kcal_mol": [0.0, 62.75, 627.51, 376.51, 125.50],
                    },
                    "b973c": {
                        "energies_hartree": [None, -100.0, -99.0, -99.5, None],
                        "relative_energies_kcal_mol": [None, 0.0, 627.51, 313.75, None],
                    },
                },
                "selection_policy": {"actual_source": "B97-3c"},
                "selections": {
                    "ts_guess": {"index": 2},
                    "intermediate": {"index": 3},
                },
                "trajectory_quality": {"off_path_indices": [4]},
                "forming_bonds": [[0, 1]],
            }
        ),
        encoding="utf-8",
    )

    output = plot_scan_profile(profile)
    expected_png = profile.parent / "scan_profile.png"

    assert output == expected_png
    assert output is not None
    assert expected_png.is_file()
    assert expected_png.stat().st_size > 0


def test_engine_endpoint_extension_and_refinement_keep_full_profile(tmp_path, monkeypatch):
    product = _xyz(tmp_path / "product.xyz", 1.5)
    progress_events = []
    engine = PEBScanEngine(
        {
            "step2": {
                "scan": {
                    "topology_guard_enabled": False,
                    "scan_start_distance": 3.4,
                    "scan_end_distance": 1.5,
                    "scan_steps": "auto",
                    "coarse_step_A": 0.20,
                    "endpoint_extension": {
                        "enabled": True,
                        "increment_A": 0.25,
                        "maximum_coordinate_A": 4.0,
                        "min_dissociation_side_span_A": 0.80,
                        "min_valid_post_ts_points": 3,
                        "max_extensions": 3,
                    },
                    "candidate_refinement": {
                        "enabled": True,
                        "ts": {"enabled": True, "half_window_A": 0.25, "step_A": 0.05},
                        "intermediate": {"enabled": True, "half_window_A": 0.25, "step_A": 0.05},
                        "min_overlap_points": 2,
                        "coordinate_tolerance_A": 0.002,
                    },
                    "selection": {
                        "preferred_energy_source": "xtb",
                        "ts": {"method": "maximum_positive_gradient"},
                        "intermediate": {"minimum_candidate_points": 2},
                    },
                },
                "energy_refinement": {"enabled": False},
            }
        },
        molecule_name="product_major",
        event_callback=lambda event, payload: progress_events.append((event, payload)),
    )

    scan_starts = []

    def fake_scan(*, start_xyz, output_dir, params, direction="outward", **_kwargs):
        scan_starts.append((Path(output_dir).name, Path(start_xyz)))
        coordinates = compute_scan_distances(
            params["scan_start_distance"],
            params["scan_end_distance"],
            params["scan_steps"],
            direction=direction,
        )
        frames = [
            _xyz(Path(output_dir) / "scan_frames" / f"frame_{index:03d}.xyz", coordinate)
            for index, coordinate in enumerate(coordinates)
        ]
        energies = [-10.0 - (float(coordinate) - 3.15) ** 2 for coordinate in coordinates]
        return SimpleNamespace(geometries=frames), energies

    monkeypatch.setattr(engine, "_execute_scan", fake_scan)
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.plot_scan_profile",
        lambda *_args, **_kwargs: None,
    )

    result = engine.run(product, tmp_path / "s2", [(0, 1)])
    profile = json.loads(Path(result[4]).read_text(encoding="utf-8"))

    coverage = profile["composite_profile"]["coverage"]
    assert coverage["coordinate_min_A"] <= 1.5
    assert coverage["coordinate_max_A"] >= 3.65
    assert profile["coarse_scan"]["reaction_coordinate_angstrom"][0] == 1.5
    assert profile["coarse_scan"]["reaction_coordinate_angstrom"][-1] == 3.4
    assert {attempt["kind"] for attempt in profile["attempts"]} >= {
        "coarse",
        "endpoint_extension",
    }
    assert profile["energy_curves"]["xtb"]["status"] == "complete"
    completed_attempts = [
        payload
        for event, payload in progress_events
        if event == "batch_finished" and payload.get("phase") == "xtb_scan"
    ]
    assert {payload["attempt_kind"] for payload in completed_attempts} >= {
        "coarse",
        "endpoint_extension",
    }
    assert all(payload["done"] == payload["total"] > 0 for payload in completed_attempts)
    assert any(
        event == "step_finished"
        and payload.get("step") == "node_selection_render"
        for event, payload in progress_events
    )
    assert scan_starts[0][1] == product
    assert all("attempts" in str(start) for _, start in scan_starts[1:])
    warm_started = [attempt for attempt in profile["attempts"] if attempt["parent_attempt_id"]]
    assert warm_started
    assert all(attempt["seed_source_attempt"] for attempt in warm_started)
    assert all(attempt["seed_source_frame_index"] is not None for attempt in warm_started)


def _run_forced_rescue_engine(
    tmp_path: Path,
    monkeypatch,
    *,
    rescue_energies_hartree,
    expected_scan_ts_candidate_index=8,
):
    product = _xyz(tmp_path / "product.xyz", 1.5)
    engine = PEBScanEngine(
        {
            "step2": {
                "scan": {
                    "topology_guard_enabled": False,
                    "scan_start_distance": 3.4,
                    "scan_end_distance": 1.5,
                    "scan_steps": 8,
                    "selection": {
                        "preferred_energy_source": "b973c",
                        "allow_monotonic_shoulder": True,
                        "shoulder_max_abs_slope_kcal_mol_per_A": 5.0,
                    },
                },
                "xtb_path": {
                    "enabled": True,
                    "retain_valid_corridor_path": True,
                },
                "energy_refinement": {"enabled": True},
                "rescue": {
                    "enabled": True,
                    "relaxed_scan": {
                        "points": 17,
                        "stretch_end_A": 3.1,
                        "single_coordinate_use_scants": True,
                    },
                },
            },
            "theory": {
                "s3_low_level": {
                    "optimization": {
                        "method": "B97-3c",
                        "solvent": "acetone",
                        "solvent_model": "CPCM",
                    }
                }
            },
        }
    )

    def fake_scan(*, output_dir, params, direction="outward", **_kwargs):
        coordinates = compute_scan_distances(
            params["scan_start_distance"],
            params["scan_end_distance"],
            params["scan_steps"],
            direction=direction,
        )
        frames = [
            _xyz(Path(output_dir) / "scan_frames" / f"frame_{index:03d}.xyz", coordinate)
            for index, coordinate in enumerate(coordinates)
        ]
        energies = [-10.0 - (float(coordinate) - 2.6) ** 2 for coordinate in coordinates]
        return SimpleNamespace(geometries=frames), energies

    monkeypatch.setattr(engine, "_execute_scan", fake_scan)

    def fake_path(*, output_dir, path_name, **_kwargs):
        coordinates = tuple(3.4 - 0.2 * index for index in range(10))
        frames = tuple(
            _xyz(Path(output_dir) / path_name / f"frame_{index:03d}.xyz", coordinate)
            for index, coordinate in enumerate(coordinates)
        )
        attempt = ScanAttempt(
            attempt_id=f"NN_xtb_path_{path_name}",
            kind="xtb_path",
            directory=Path(output_dir) / path_name,
            frame_paths=frames,
            target_coordinates_A=coordinates,
            xtb_energies_hartree=tuple(_hartree_profile([0.0, 0.3, 0.8, 1.2, 2.0, 4.0, 2.5, 1.0, 0.2, 0.0])),
            off_path_indices=(0, 1, 2),
            trajectory_quality={
                "checked": True,
                "persistent_off_path_start": 0,
                "usable_end_index": -1,
                "topology_drift_index": 0,
                "last_valid_before_drift_index": None,
            },
        )
        return attempt, {
            "path_name": path_name,
            "path_arclength": [float(index) for index in range(len(frames))],
        }

    monkeypatch.setattr(engine, "_execute_xtb_path", fake_path)
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.ScanEnergyRefiner.refine",
        lambda _self, frames, point_ids=None: {
            "status": "complete",
            "energies_hartree": _hartree_profile(
                [0.0, 0.3, 0.8, 1.2, 2.0, 4.0, 2.5, 1.0, 0.2, 0.0]
            ),
            "records": [],
        },
    )
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.plot_scan_profile",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.relaxed_scan_rescue.plot_scan_profile",
        lambda *_args, **_kwargs: None,
    )

    rescue_frames = [
        _xyz(
            tmp_path / "rescue_frames" / f"frame_{index:03d}.xyz",
            1.5 + 0.1 * index,
        )
        for index in range(17)
    ]

    def fake_surface_scan(spec, *_args, **_kwargs):
        assert spec.method == "B97-3c"
        return QCJobResult(
            "complete",
            product,
            output_file=tmp_path / "rescue.out",
            extra={
                "frames": [str(frame) for frame in rescue_frames],
                "energies_hartree": list(rescue_energies_hartree),
                "energy_source": "orca.relaxscanact.dat",
                "scan_ts_candidate_xyz": str(
                    rescue_frames[expected_scan_ts_candidate_index]
                ),
            },
        )

    monkeypatch.setattr(
        "rph_core.steps.step2_retro.relaxed_scan_rescue.run_surface_scan",
        fake_surface_scan,
    )

    result = engine.run(product, tmp_path / "s2", [(0, 1)])
    profile = json.loads(Path(result[4]).read_text(encoding="utf-8"))
    return result, profile


@pytest.mark.parametrize(
    ("relative_kcal", "expected_evidence"),
    [
        (
            [
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
            ],
                "knee_shifted",
        ),
        (
            [
                5.1,
                5.0,
                4.9,
                4.8,
                4.7,
                4.6,
                4.25,
                4.25,
                4.25,
                4.25,
                4.25,
                4.2,
                3.8,
                1.5,
                0.0,
                0.1,
                0.0,
            ],
                "knee_shifted",
        ),
    ],
)
def test_engine_rescue_dispatches_only_normal_ts_int_jobs(
    tmp_path,
    monkeypatch,
    relative_kcal,
    expected_evidence,
):
    _result, profile = _run_forced_rescue_engine(
        tmp_path,
        monkeypatch,
        rescue_energies_hartree=_hartree_profile(relative_kcal),
    )

    assert profile["selection_source"] == "orca_relaxed_scan"
    assert profile["s2_state"] == "rescue_seeded"
    assert profile["seed_evidence"] == expected_evidence
    assert profile["ts_search_seed"] is not None
    assert profile["s3_dispatch"] == {
        "resolution": "rescue_seeded",
        "submit_ts": True,
        "submit_intermediate": True,
        "neb_eligible": False,
        "source": "b973c_relaxed_scan",
        "path_fully_distorted": True,
    }
    assert profile["selections"]["ts_guess"]["source"] == "orca_relaxed_scan"
    assert profile["selections"]["intermediate"]["source"] == "orca_relaxed_scan"


@pytest.mark.parametrize(
    ("relative_kcal", "expected_reason"),
    [
        (
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            {"no_valid_knee_point"},
        ),
        (
            [1.0, 1.2, 1.4, 1.8, 2.0, 2.2],
            None,
        ),
    ],
)
def test_engine_unresolved_rescue_never_dispatches_fabricated_s3_jobs(
    tmp_path,
    monkeypatch,
    relative_kcal,
    expected_reason,
):
    energies_hartree = _hartree_profile(relative_kcal)
    _result, profile = _run_forced_rescue_engine(
        tmp_path,
        monkeypatch,
        rescue_energies_hartree=energies_hartree,
    )

    # An unresolved rescue may have been attempted with ORCA, but it must not
    # advertise a stale or fabricated seed as canonical output.
    assert profile["selection_source"] is None
    assert profile["s2_state"] == "unresolved"
    assert profile["s3_dispatch"]["resolution"] == "unresolved"
    assert profile["s3_dispatch"]["submit_ts"] is False
    assert profile["s3_dispatch"]["submit_intermediate"] is False
    assert profile["s3_dispatch"]["neb_eligible"] is False
    if expected_reason is None:
        assert profile["rescue"]["error"] == "Relaxed scan lacks complete frame-energy coverage"
        assert profile["rejection_reason"] is None
    else:
        assert profile["rejection_reason"] in expected_reason
