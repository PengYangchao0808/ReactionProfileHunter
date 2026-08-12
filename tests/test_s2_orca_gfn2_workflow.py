from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from rph_core.steps.step2_retro.peb_engine import PEBScanEngine
from rph_core.utils.orca_interface import ORCAInterface


def _xyz(path: Path, distance: float) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"3\nframe\nH 0.0 0.0 0.0\nH {distance:.6f} 0.0 0.0\nH 0.0 2.0 0.0\n",
        encoding="utf-8",
    )
    return path


def _config() -> dict:
    return {
        "step2": {
            "method": "orca_gfn2_relaxed_scan",
            "scan": {
                "scan_start_distance": 3.4,
                "scan_end_distance": 1.2,
                "scan_steps": "auto",
                "selection": {
                    "min_reaction_progress": 0.0,
                    "endpoint_exclusion_frames": 1,
                    "minimum_clean_frames_after_knee": 2,
                    "endpoint": {"guard_frames": 1, "min_valid_frames": 3},
                    "knee": {"enabled": True, "smoothing_window": 5},
                    "ts_seed": {"right_shift": {"base_A": 0.15, "min_A": 0.05, "max_A": 0.40}},
                },
            },
            "orca_gfn2_scan": {
                "enabled": True,
                "relaxed_scan": {
                    "stretch_end_A": 3.4,
                    "points": 17,
                    "single_coordinate_use_scants": True,
                    "solvent": "acetone",
                    "solvent_model": "ALPB",
                },
            },
            "rescue": {"enabled": False, "relaxed_scan": {"stretch_end_A": 3.4, "points": 17}},
            "energy_refinement": {
                "enabled": True,
                "engine": "orca",
                "method": "B97-3c",
                "parallel_jobs": 1,
                "cores_per_job": 1,
                "memory_per_job": "2GB",
            },
        },
        "theory": {
            "s3_low_level": {
                "optimization": {
                    "engine": "orca",
                    "method": "B97-3c",
                    "solvent": "acetone",
                    "solvent_model": "CPCM",
                }
            }
        },
        "resources": {"nproc": 1},
    }


def test_orca_renderer_emits_alpb_for_gfn2_without_changing_cpcm():
    gfn2 = ORCAInterface(method="GFN2-xTB", solvent="acetone", solvent_model="ALPB")
    b973c = ORCAInterface(method="B97-3c", solvent="acetone", solvent_model="CPCM")

    gfn2_route = gfn2._render_route(task_type="opt")
    b973c_route = b973c._render_route(task_type="opt")

    assert "ALPB(Acetone)" in gfn2_route
    assert "CPCM(" not in gfn2_route
    assert "CPCM(Acetone)" in b973c_route


def test_orca_renderer_translates_cpcm_or_smd_to_alpb_for_gfn_xtb():
    gfn2_cpcm = ORCAInterface(method="GFN2-xTB", solvent="acetone", solvent_model="CPCM")
    gfn2_default = ORCAInterface(method="GFN2-xTB", solvent="acetone")
    b973c = ORCAInterface(method="B97-3c", solvent="acetone", solvent_model="SMD")

    assert gfn2_cpcm.solvent_model == "ALPB"
    assert gfn2_default.solvent_model == "ALPB"
    assert "ALPB(Acetone)" in gfn2_cpcm._render_route(task_type="opt")
    assert "CPCM(" not in gfn2_cpcm._render_route(task_type="opt")
    assert gfn2_cpcm._render_cpcm_block() == ""
    assert b973c.solvent_model == "SMD"


def test_surface_scan_preflight_fails_early_when_xtb_missing(tmp_path, monkeypatch):
    from rph_core.utils.qc_models import SurfaceScanCoordinate, SurfaceScanSpec

    product = _xyz(tmp_path / "product.xyz", 1.5)
    orca_dir = tmp_path / "orca_bin"
    orca_dir.mkdir()
    monkeypatch.setattr(
        "rph_core.utils.orca_interface.ORCAInterface._find_orca_binary",
        lambda *args, **kwargs: orca_dir / "orca",
    )
    interface = ORCAInterface(method="GFN2-xTB", solvent="acetone", solvent_model="ALPB")
    spec = SurfaceScanSpec(
        method="GFN2-xTB",
        solvent="acetone",
        solvent_model="ALPB",
        coordinates=(SurfaceScanCoordinate(kind="B", atoms=(0, 1), start=1.5, end=3.4, steps=17),),
        simultaneous=False,
        scan_ts=True,
        full_scan=True,
    )
    result = interface.run_surface_scan(spec, product, tmp_path / "scan")
    assert result.status == "failed"
    assert "xtb binary in the ORCA directory" in (result.error or "")
    assert str(orca_dir / "xtb") in (result.error or "")


def test_v4_s2_uses_orca_gfn2_then_b973c_sp_without_rescue(tmp_path, monkeypatch):
    product = _xyz(tmp_path / "product.xyz", 1.5)
    methods: list[str] = []
    specs = []

    def fake_scan(spec, input_xyz, output_dir, config, **_kwargs):
        methods.append(spec.method)
        specs.append(spec)
        frames = [_xyz(Path(output_dir) / "scan_frames" / f"frame_{i:04d}.xyz", 1.5 + 0.11 * i) for i in range(17)]
        return SimpleNamespace(
            status="complete",
            output_file=Path(output_dir) / "gfn2_scan.out",
            output_xyz=None,
            error=None,
            extra={
                "frames": [str(frame) for frame in frames],
                "energies_hartree": [-100.0 + value / 627.509 for value in [0, .1, .5, 1.5, 3, 6, 10, 8, 5, 3, 2, 1.5, 1, .5, .2, .1, 0]],
                "energy_source": "gfn2_scan_ledger",
            },
        )

    monkeypatch.setattr("rph_core.steps.step2_retro.relaxed_scan_rescue.run_surface_scan", fake_scan)
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.ScanEnergyRefiner.refine",
        lambda _self, frames, point_ids=None: {
            "status": "complete",
            "method": "B97-3c",
            "energies_hartree": [-100.0 + value / 627.509 for value in [0, .1, .5, 1.5, 3, 6, 10, 8, 5, 3, 2, 1.5, 1, .5, .2, .1, 0]],
        },
    )
    monkeypatch.setattr("rph_core.steps.step2_retro.peb_engine.plot_scan_profile", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.path_profile.assess_path_topology",
        lambda *_args, **_kwargs: {"off_path_indices": [], "topology_reason_by_frame": []},
    )

    result = PEBScanEngine(_config(), molecule_name="product").run(product, tmp_path / "s2", [(0, 1)])
    profile = json.loads(Path(result[4]).read_text(encoding="utf-8"))

    assert methods == ["GFN2-xTB"]
    assert specs[0].solvent == "acetone"
    assert specs[0].solvent_model == "ALPB"
    assert profile["scan_engine"] == "orca"
    assert profile["scan_method"] == "GFN2-xTB"
    assert profile["energy_refinement_method"] == "B97-3c"
    assert profile["selection_source"] == "orca_gfn2_b973c_sp"
    assert profile["s2_state"] == "gfn2_seeded"
    assert profile["s3_dispatch"]["submit_ts"] is True
    assert not (tmp_path / "s2" / "rescue").exists()


def test_v4_s2_enters_b973c_rescue_only_after_gfn2_distortion(tmp_path, monkeypatch):
    product = _xyz(tmp_path / "product.xyz", 1.5)
    methods: list[str] = []
    specs = []
    relative = [0, .1, .5, 1.5, 3, 6, 10, 8, 5, 3, 2, 1.5, 1, .5, .2, .1, 0]

    def fake_scan(spec, input_xyz, output_dir, config, **_kwargs):
        methods.append(spec.method)
        specs.append(spec)
        frames = [_xyz(Path(output_dir) / "scan_frames" / f"frame_{i:04d}.xyz", 1.5 + 0.11 * i) for i in range(17)]
        return SimpleNamespace(
            status="complete",
            output_file=Path(output_dir) / "scan.out",
            output_xyz=None,
            error=None,
            extra={
                "frames": [str(frame) for frame in frames],
                "energies_hartree": [-100.0 + value / 627.509 for value in relative],
                "energy_source": "scan_ledger",
            },
        )

    monkeypatch.setattr("rph_core.steps.step2_retro.relaxed_scan_rescue.run_surface_scan", fake_scan)
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.ScanEnergyRefiner.refine",
        lambda _self, frames, point_ids=None: {
            "status": "complete",
            "method": "B97-3c",
            "energies_hartree": [-100.0 + value / 627.509 for value in relative],
        },
    )
    monkeypatch.setattr("rph_core.steps.step2_retro.peb_engine.plot_scan_profile", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.path_profile.assess_path_topology",
        lambda _product, _symbols, _bonds, frames: {
            "off_path_indices": [] if any(Path(frame).parent.name == "rescue" for frame in frames) else list(range(len(frames))),
            "topology_reason_by_frame": [],
        },
    )
    config = _config()
    config["step2"]["rescue"]["enabled"] = True

    result = PEBScanEngine(config, molecule_name="product").run(product, tmp_path / "s2", [(0, 1)])
    profile = json.loads(Path(result[4]).read_text(encoding="utf-8"))

    assert methods == ["GFN2-xTB", "B97-3c"]
    assert [spec.solvent_model for spec in specs] == ["ALPB", "CPCM"]
    assert profile["selection_source"] == "orca_b973c_relaxed_scan"
    assert profile["s2_state"] == "rescue_seeded"
    assert profile["s3_dispatch"]["submit_ts"] is True
    assert (tmp_path / "s2" / "rescue").is_dir()


def test_v4_s2_retains_late_tail_distortion_without_rescue(tmp_path, monkeypatch):
    product = _xyz(tmp_path / "product.xyz", 1.5)
    methods: list[str] = []
    relative = [0, .1, .5, 1.5, 3, 6, 10, 8, 5, 3, 2, 1.5, 1, .5, .2, .1, 0]

    def fake_scan(spec, input_xyz, output_dir, config, **_kwargs):
        methods.append(spec.method)
        frames = [_xyz(Path(output_dir) / "scan_frames" / f"frame_{i:04d}.xyz", 1.5 + 0.11 * i) for i in range(17)]
        return SimpleNamespace(
            status="complete",
            output_file=Path(output_dir) / "scan.out",
            output_xyz=None,
            error=None,
            extra={
                "frames": [str(frame) for frame in frames],
                "energies_hartree": [-100.0 + value / 627.509 for value in relative],
                "energy_source": "scan_ledger",
            },
        )

    monkeypatch.setattr("rph_core.steps.step2_retro.relaxed_scan_rescue.run_surface_scan", fake_scan)
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.ScanEnergyRefiner.refine",
        lambda _self, frames, point_ids=None: {
            "status": "complete",
            "method": "B97-3c",
            "energies_hartree": [-100.0 + value / 627.509 for value in relative],
        },
    )
    monkeypatch.setattr("rph_core.steps.step2_retro.peb_engine.plot_scan_profile", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.path_profile.assess_path_topology",
        lambda _product, _symbols, _bonds, frames: {
            "off_path_indices": [15, 16]
            if any(Path(frame).parent.name == "scan" for frame in frames)
            else [],
            "topology_reason_by_frame": [],
        },
    )
    config = _config()
    config["step2"]["rescue"]["enabled"] = True

    result = PEBScanEngine(config, molecule_name="product").run(product, tmp_path / "s2", [(0, 1)])
    profile = json.loads(Path(result[4]).read_text(encoding="utf-8"))

    assert methods == ["GFN2-xTB"]
    assert profile["selection_source"] == "orca_gfn2_b973c_sp"
    assert profile["s2_state"] == "gfn2_seeded"
    assert profile["trajectory_quality"]["topology_state"] == "tail_distorted_after_knee"
    decision = profile["trajectory_quality"]["topology_rescue_decision"]
    assert decision["distortion_frame_index"] == 15
    assert decision["knee_frame_index"] is not None
    assert decision["knee_frame_index"] < decision["distortion_frame_index"]
    assert decision["clean_frames_after_knee"] >= 2
    assert decision["minimum_clean_frames_after_knee"] == 2
    assert decision["rescue_required"] is False
    assert profile["s3_dispatch"]["path_fully_distorted"] is False
    assert not (tmp_path / "s2" / "rescue").exists()


def test_topology_rescue_requires_knee_and_two_clean_post_knee_frames():
    profile = SimpleNamespace(
        excluded_frames=(9,),
        frames=tuple(
            SimpleNamespace(frame_index=index, topology_valid=index != 9)
            for index in range(12)
        ),
    )
    selection = SimpleNamespace(diagnostics={"knee_frame_index": 6})

    decision = PEBScanEngine._topology_rescue_decision(profile, selection, 2)

    assert decision["post_knee_frame_indices"] == [7, 8]
    assert decision["clean_frames_after_knee"] == 2
    assert decision["rescue_required"] is False
    assert decision["decision"] == "retain_primary_knee_supported"


def test_topology_rescue_triggers_when_distortion_leaves_one_post_knee_frame():
    profile = SimpleNamespace(
        excluded_frames=(8,),
        frames=tuple(
            SimpleNamespace(frame_index=index, topology_valid=index != 8)
            for index in range(12)
        ),
    )
    selection = SimpleNamespace(diagnostics={"knee_frame_index": 6})

    decision = PEBScanEngine._topology_rescue_decision(profile, selection, 2)

    assert decision["post_knee_frame_indices"] == [7]
    assert decision["clean_frames_after_knee"] == 1
    assert decision["rescue_required"] is True
    assert decision["decision"] == "rescue_insufficient_post_knee_support"


def test_v4_s2_tolerates_rescue_warranted_drift_when_rescue_disabled(tmp_path, monkeypatch):
    """With step2.rescue.enabled=false, a drift that would warrant the B97-3c
    rescue is tolerated: the primary GFN2 selection is retained and no rescue
    scan is launched (GFN2-only robustness mode)."""
    product = _xyz(tmp_path / "product.xyz", 1.5)
    methods: list[str] = []
    relative = [0, .1, .5, 1.5, 3, 6, 10, 8, 5, 3, 2, 1.5, 1, .5, .2, .1, 0]

    def fake_scan(spec, input_xyz, output_dir, config, **_kwargs):
        methods.append(spec.method)
        frames = [_xyz(Path(output_dir) / "scan_frames" / f"frame_{i:04d}.xyz", 1.5 + 0.11 * i) for i in range(17)]
        return SimpleNamespace(
            status="complete",
            output_file=Path(output_dir) / "scan.out",
            output_xyz=None,
            error=None,
            extra={
                "frames": [str(frame) for frame in frames],
                "energies_hartree": [-100.0 + value / 627.509 for value in relative],
                "energy_source": "scan_ledger",
            },
        )

    monkeypatch.setattr("rph_core.steps.step2_retro.relaxed_scan_rescue.run_surface_scan", fake_scan)
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.peb_engine.ScanEnergyRefiner.refine",
        lambda _self, frames, point_ids=None: {
            "status": "complete",
            "method": "B97-3c",
            "energies_hartree": [-100.0 + value / 627.509 for value in relative],
        },
    )
    monkeypatch.setattr("rph_core.steps.step2_retro.peb_engine.plot_scan_profile", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "rph_core.steps.step2_retro.path_profile.assess_path_topology",
        lambda _product, _symbols, _bonds, frames: {
            "off_path_indices": [0, 1, 2]
            if any(Path(frame).parent.name == "scan" for frame in frames)
            else [],
            "topology_reason_by_frame": [],
        },
    )
    config = _config()
    config["step2"]["rescue"]["enabled"] = False

    result = PEBScanEngine(config, molecule_name="product").run(product, tmp_path / "s2", [(0, 1)])
    profile = json.loads(Path(result[4]).read_text(encoding="utf-8"))

    assert methods == ["GFN2-xTB"]
    decision = profile["trajectory_quality"]["topology_rescue_decision"]
    assert decision["rescue_required"] is True
    assert decision["decision"] == "rescue_knee_not_before_distortion"
    assert profile["selection_source"] == "orca_gfn2_b973c_sp"
    assert profile["s2_state"] == "gfn2_seeded"
    assert not (tmp_path / "s2" / "rescue").exists()
