import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from rph_core.steps.step2_retro.energy_refinement import ScanEnergyRefiner
from rph_core.steps.step2_retro.peb_engine import PEBScanEngine
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


def test_refined_curve_uses_pre_ts_arclength_midpoint_when_no_basin(tmp_path):
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
    assert int_index == 2
    assert intermediate["selection_mode"] == "midpoint_fallback"
    assert intermediate["rule"] == "pre_ts_arclength_midpoint"
    assert intermediate["reason"] == "no_resolved_pre_ts_basin"
    assert status == "selected"


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

    assert result[2] is None
    assert result[5] == "DEGRADED"
    assert not (output_dir / "intermediate.xyz").exists()
    assert not (output_dir / "reactant_complex.xyz").exists()
    assert not (output_dir / "manifest.json").exists()
    assert profile["artifact_lifecycle"]["legacy_reactant_complex_written"] is False
    assert profile["profile_schema_version"] == "s2_scan_profile_v9"
    assert profile["attempts"][0]["kind"] == "coarse"
    assert profile["composite_profile"]["coverage"]["coordinate_min_A"] == 1.5
    assert profile["composite_profile"]["coverage"]["coordinate_max_A"] == 3.4
    assert profile["energy_curves"]["xtb"]["status"] == "complete"
    assert list((output_dir / "superseded").glob("*/manifest.json"))
