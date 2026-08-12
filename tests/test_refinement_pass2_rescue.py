from __future__ import annotations

from collections import deque
from dataclasses import replace
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest

import rph_core.steps.refinement.engine as engine_module
from rph_core.steps.fidelity_profile import FidelityProfile
from rph_core.steps.refinement import RefinementEngine
from rph_core.steps.refinement.models import Pass1Outcome, StructureRequest
from rph_core.utils.attempt_recorder import AttemptRecorder
from rph_core.utils.config_loader import load_config
from rph_core.utils.qc_models import QCJobResult


def _engine(stage: str = "S3", **profile_overrides: object) -> RefinementEngine:
    config = load_config()
    profile = FidelityProfile.from_config(config, stage)
    if profile_overrides:
        profile = replace(profile, **profile_overrides)
    return RefinementEngine(config, profile)


def _write_xyz(path: Path, distance: float = 1.8) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"2\nmock\nH 0.0 0.0 0.0\nH {distance:.6f} 0.0 0.0\n",
        encoding="utf-8",
    )
    return path


def _request(
    tmp_path: Path,
    *,
    structure_id: str,
    role: str,
    kind: str,
    seed_distance: float | None = None,
) -> StructureRequest:
    inputs = tmp_path / "inputs"
    input_xyz = _write_xyz(inputs / f"{structure_id}.xyz")
    original_seed_xyz = _write_xyz(inputs / f"{structure_id}_seed.xyz", seed_distance) if seed_distance is not None else None
    return StructureRequest(
        id=structure_id,
        role=role,
        kind=kind,
        input_xyz=input_xyz,
        original_seed_xyz=original_seed_xyz,
        source_stage="S2",
        forming_bonds=[(0, 1)] if role in {"intermediate", "ts"} else [],
        s1_ensemble_thermodynamics={
            "ensemble_thermochemistry_correction_hartree": 0.0,
            "ensemble_enthalpy_correction_hartree": 0.0,
            "ensemble_correction_id": "test",
        },
        ensemble_thermochemistry_correction_hartree=0.0,
    )


def _manifest_for(engine: RefinementEngine, requests: list[StructureRequest], output_dir: Path) -> dict:
    import json

    return json.loads(engine.run(requests, output_dir).read_text(encoding="utf-8"))


def _install_qc_mocks(
    monkeypatch: pytest.MonkeyPatch,
    *,
    frequency_sequence: list[tuple[float, ...]],
    fail_optimization: Callable[[dict[str, object]], bool] | None = None,
):
    opt_calls: list[dict[str, Any]] = []
    freq_calls: list[dict[str, Any]] = []
    sp_calls: list[dict[str, Any]] = []
    frequencies = deque(frequency_sequence)
    fallback_frequencies = frequency_sequence[-1] if frequency_sequence else (50.0, 150.0, 250.0)

    def fake_run_optimization(spec, input_xyz, output_dir, config):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        call = {
            "spec": spec,
            "input_xyz": Path(input_xyz),
            "output_dir": output_dir,
            "config": config,
        }
        opt_calls.append(call)
        if fail_optimization is not None and fail_optimization(call):
            return QCJobResult(
                "failed",
                Path(input_xyz),
                output_file=output_dir / "opt.out",
                error="mocked optimization failure",
            )
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
        out_file = output_dir / "opt.out"
        out_file.write_text("mock opt\n", encoding="utf-8")
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_xyz=out_xyz,
            output_file=out_file,
            energy_hartree=-1.0,
        )

    def fake_run_frequency(spec, input_xyz, output_dir, config):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / "freq.out"
        hess_file = out_file.with_suffix(".hess")
        hess_file.write_text("mock hess\n", encoding="utf-8")
        freqs = frequencies.popleft() if frequencies else fallback_frequencies
        mode_blocks = []
        imaginary_mode = 0
        for frequency in freqs:
            if float(frequency) >= 0.0:
                continue
            mode_blocks.append(
                "\n".join(
                    (
                        f"Mode:   {imaginary_mode}",
                        f"Freq:   {float(frequency):.2f} cm**-1 (imaginary mode)",
                        "        dx          dy          dz",
                        "Atom 0:  1.000000   0.000000   0.000000",
                        "Atom 1: -1.000000   0.000000   0.000000",
                    )
                )
            )
            imaginary_mode += 1
        out_file.write_text(
            "CARTESIAN DISPLACEMENTS\n-----------------------\n" + "\n".join(mode_blocks) + "\n",
            encoding="utf-8",
        )
        call = {
            "spec": spec,
            "input_xyz": Path(input_xyz),
            "output_dir": output_dir,
            "config": config,
            "frequencies": freqs,
        }
        freq_calls.append(call)
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_file=out_file,
            energy_hartree=-0.9,
            frequencies_cm1=freqs,
            enthalpy_hartree=-0.85,
            gibbs_free_energy_hartree=-0.88,
            gibbs_correction_hartree=0.02,
        )

    def fake_run_single_point(spec, input_xyz, output_dir, config):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / "sp.out"
        out_file.write_text("mock sp\n", encoding="utf-8")
        sp_calls.append({"spec": spec, "input_xyz": Path(input_xyz), "output_dir": output_dir})
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_file=out_file,
            energy_hartree=-1.2,
        )

    monkeypatch.setattr(engine_module, "run_optimization", fake_run_optimization)
    monkeypatch.setattr(engine_module, "run_frequency", fake_run_frequency)
    monkeypatch.setattr(engine_module, "run_single_point", fake_run_single_point)
    return opt_calls, freq_calls, sp_calls


def _patch_ts_classifier(monkeypatch: pytest.MonkeyPatch, mapping: dict[tuple[float, ...], dict]):
    def fake_classify_ts(**kwargs):
        key = tuple(float(value) for value in kwargs["frequencies_cm1"])
        return dict(
            mapping.get(
                key,
                {
                    "hessian_index": 0,
                    "mode_identity": "unavailable",
                    "stationary_point_class": "unclassifiable",
                },
            )
        )

    monkeypatch.setattr(engine_module.identity_module, "classify_ts", fake_classify_ts)


def _patch_int_classifier(monkeypatch: pytest.MonkeyPatch, mapping: dict[tuple[float, ...], dict]):
    def fake_classify_int(**kwargs):
        key = tuple(float(value) for value in kwargs["frequencies_cm1"])
        return dict(mapping.get(key, {"identity": "opt_failed"}))

    monkeypatch.setattr(engine_module.identity_module, "classify_int", fake_classify_int)


def _patch_minimum_classifier(monkeypatch: pytest.MonkeyPatch, mapping: dict[tuple[float, ...], dict]):
    def fake_classify_minimum(**kwargs):
        key = tuple(float(value) for value in kwargs["frequencies_cm1"])
        return dict(mapping.get(key, {"identity": "not_checked", "imaginary_count": 0}))

    monkeypatch.setattr(engine_module.identity_module, "classify_minimum", fake_classify_minimum)


@pytest.mark.xfail(reason="R0/R1 supersedes the removed mode-directed rescue ladder")
def test_ts_rescue_l1_read_hessian_spec(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls, _sp_calls = _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[(-20.0, 120.0, 250.0), (-150.0, 120.0, 250.0)],
    )
    _patch_ts_classifier(
        monkeypatch,
        {
            (-20.0, 120.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "soft_target_ts",
            },
            (-150.0, 120.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "valid_target_ts",
            },
        },
    )
    engine = _engine("S3", irc_enabled=False)
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts", seed_distance=2.5)

    _manifest_for(engine, [request], tmp_path / "stage")

    assert any(
        call["spec"].initial_hessian == "read"
        and call["spec"].hessian_filename
        and call["spec"].ts_mode == 0
        for call in opt_calls
    )


@pytest.mark.xfail(reason="R0/R1 supersedes the removed mode-directed rescue ladder")
def test_ts_rescue_l2_recalc_hess_spec(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls, _sp_calls = _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[(-20.0, 120.0, 250.0), (-20.0, 120.0, 250.0), (-150.0, 120.0, 250.0)],
    )
    _patch_ts_classifier(
        monkeypatch,
        {
            (-20.0, 120.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "soft_target_ts",
            },
            (-150.0, 120.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "valid_target_ts",
            },
        },
    )
    engine = _engine("S3", irc_enabled=False)
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts", seed_distance=2.5)

    _manifest_for(engine, [request], tmp_path / "stage")

    assert any(
        call["spec"].initial_hessian == "calculate"
        and call["spec"].recalc_hessian == 5
        and call["spec"].trust == pytest.approx(0.3)
        for call in opt_calls
    )


def test_target_mode_selection_uses_forming_bond_overlap(tmp_path: Path) -> None:
    engine = _engine("S3", irc_enabled=False)
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")
    coordinates = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    displacements = {
        0: np.asarray([[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]]),
        1: np.asarray([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]),
    }

    selection = engine._select_target_ts_mode(request, coordinates, displacements)

    assert selection["selection_status"] == "selected"
    assert selection["selected_mode"] == 1
    assert selection["selected_overlap"] == pytest.approx(1.0)


def test_partial_int_candidate_cannot_short_circuit_rescue(tmp_path: Path):
    engine = _engine("S3", irc_enabled=False)
    request = _request(tmp_path, structure_id="int", role="intermediate", kind="minimum")

    assert not engine._candidate_is_valid_for_request(
        request,
        {
            "opt_status": "failed",
            "frequency_status": "not_run",
            "minimum_classification": {"imaginary_count": 0},
            "int_classification": {"identity": "distinct_intermediate"},
        },
    )


def test_nonconvergence_trigger_accepts_orca_normal_termination_wording(tmp_path: Path):
    engine = _engine("S3", irc_enabled=False)
    outcome = Pass1Outcome(
        structure_id="ts",
        role="ts",
        kind="ts",
        opt_status="failed",
        opt_error="ORCA terminated normally, but the geometry optimization did not converge",
        attempt_history=[{"failure_type": None, "stop_reason": None}],
    )

    assert engine._is_geometry_nonconvergence(outcome)


def test_ts_rescue_not_run_for_valid_ts(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch, frequency_sequence=[(-150.0, 120.0, 250.0)])
    _patch_ts_classifier(
        monkeypatch,
        {
            (-150.0, 120.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "valid_target_ts",
            },
        },
    )
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [_request(tmp_path, structure_id="ts", role="ts", kind="ts", seed_distance=2.5)],
        tmp_path / "stage",
    )

    assert manifest["structures"][0]["pass2_rescue_attempts"] == []


def test_ts_rescue_l2_not_run_if_l1_succeeds(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls, _sp_calls = _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[(-20.0, 120.0, 250.0), (-150.0, 120.0, 250.0)],
    )
    _patch_ts_classifier(
        monkeypatch,
        {
            (-20.0, 120.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "soft_target_ts",
            },
            (-150.0, 120.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "valid_target_ts",
            },
        },
    )
    engine = _engine("S3", irc_enabled=False)
    _manifest_for(
        engine,
        [_request(tmp_path, structure_id="ts", role="ts", kind="ts", seed_distance=2.5)],
        tmp_path / "stage",
    )

    assert not any(call["spec"].recalc_hessian == 5 for call in opt_calls)


@pytest.mark.xfail(reason="INT now has one CalcAll R1, without TightOpt")
def test_int_rescue_l1_tightopt_spec(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls, _sp_calls = _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[(-80.0, 120.0, 250.0), (50.0, 120.0, 250.0)],
    )
    _patch_int_classifier(
        monkeypatch,
        {
            (-80.0, 120.0, 250.0): {"identity": "imaginary_frequency"},
            (50.0, 120.0, 250.0): {"identity": "distinct_intermediate"},
        },
    )
    _patch_minimum_classifier(
        monkeypatch,
        {
            (-80.0, 120.0, 250.0): {"identity": "imaginary_frequency", "imaginary_count": 1},
            (50.0, 120.0, 250.0): {"identity": "valid_minimum", "imaginary_count": 0},
        },
    )
    engine = _engine("S3", irc_enabled=False)
    _manifest_for(
        engine,
        [_request(tmp_path, structure_id="int", role="intermediate", kind="minimum", seed_distance=2.4)],
        tmp_path / "stage",
    )

    assert any(call["spec"].route_extras == "TightOpt" for call in opt_calls)


def test_int_rescue_l2_mode_displacement(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls, _sp_calls = _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[(-80.0, 120.0, 250.0), (-80.0, 120.0, 250.0), (50.0, 120.0, 250.0)],
    )
    _patch_int_classifier(
        monkeypatch,
        {
            (-80.0, 120.0, 250.0): {"identity": "imaginary_frequency"},
            (50.0, 120.0, 250.0): {"identity": "distinct_intermediate"},
        },
    )
    _patch_minimum_classifier(
        monkeypatch,
        {
            (-80.0, 120.0, 250.0): {"identity": "imaginary_frequency", "imaginary_count": 1},
            (50.0, 120.0, 250.0): {"identity": "valid_minimum", "imaginary_count": 0},
        },
    )
    monkeypatch.setattr(
        engine_module,
        "_parse_orca_displacement_vectors",
        lambda _content: {0: np.array([[0.2, 0.0, 0.0], [-0.2, 0.0, 0.0]])},
    )
    engine = _engine("S3", irc_enabled=False)
    _manifest_for(
        engine,
        [_request(tmp_path, structure_id="int", role="intermediate", kind="minimum", seed_distance=2.4)],
        tmp_path / "stage",
    )

    displaced_calls = [
        call for call in opt_calls if "int_rescue_l2" in str(call["output_dir"])
    ]
    assert len(displaced_calls) >= 1
    assert any(call["spec"].initial_hessian == "calculate" for call in displaced_calls)


@pytest.mark.xfail(reason="INT now has one CalcAll R1 rather than the old L3")
def test_int_rescue_l3_recalc_hess_5(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls, _sp_calls = _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[
            (-80.0, 120.0, 250.0),
            (-80.0, 120.0, 250.0),
            (-80.0, 120.0, 250.0),
            (-80.0, 120.0, 250.0),
            (50.0, 120.0, 250.0),
        ],
    )
    _patch_int_classifier(
        monkeypatch,
        {
            (-80.0, 120.0, 250.0): {"identity": "imaginary_frequency"},
            (50.0, 120.0, 250.0): {"identity": "distinct_intermediate"},
        },
    )
    _patch_minimum_classifier(
        monkeypatch,
        {
            (-80.0, 120.0, 250.0): {"identity": "imaginary_frequency", "imaginary_count": 1},
            (50.0, 120.0, 250.0): {"identity": "valid_minimum", "imaginary_count": 0},
        },
    )
    monkeypatch.setattr(
        engine_module,
        "_parse_orca_displacement_vectors",
        lambda _content: {0: np.array([[0.2, 0.0, 0.0], [-0.2, 0.0, 0.0]])},
    )
    engine = _engine("S3", irc_enabled=False)
    _manifest_for(
        engine,
        [_request(tmp_path, structure_id="int", role="intermediate", kind="minimum", seed_distance=2.4)],
        tmp_path / "stage",
    )

    assert any(call["spec"].recalc_hessian == 5 for call in opt_calls)


@pytest.mark.xfail(reason="minimum rescue now only follows a failed R0")
def test_minimum_rescue_for_precursor_with_imaginary(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls, _sp_calls = _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[(-60.0, 120.0, 250.0), (50.0, 120.0, 250.0)],
    )
    _patch_minimum_classifier(
        monkeypatch,
        {
            (-60.0, 120.0, 250.0): {"identity": "imaginary_frequency", "imaginary_count": 1},
            (50.0, 120.0, 250.0): {"identity": "valid_minimum", "imaginary_count": 0},
        },
    )
    engine = _engine("S3", irc_enabled=False)
    _manifest_for(
        engine,
        [_request(tmp_path, structure_id="precursor", role="precursor", kind="minimum")],
        tmp_path / "stage",
    )

    assert any(call["spec"].route_extras == "TightOpt" for call in opt_calls)


@pytest.mark.xfail(reason="a successful R0 intentionally has no R1 attempt")
def test_rescue_attempts_appended_to_pass2_rescue_attempts(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[(-20.0, 120.0, 250.0), (-150.0, 120.0, 250.0)],
    )
    _patch_ts_classifier(
        monkeypatch,
        {
            (-20.0, 120.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "soft_target_ts",
            },
            (-150.0, 120.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "valid_target_ts",
            },
        },
    )
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [_request(tmp_path, structure_id="ts", role="ts", kind="ts", seed_distance=2.5)],
        tmp_path / "stage",
    )

    attempts = manifest["structures"][0]["pass2_rescue_attempts"]
    assert attempts
    assert {"attempt_id", "strategy", "level", "status", "error"} <= set(attempts[0])


@pytest.mark.xfail(reason="covers removed ts_rescue_l1 scheduling")
def test_rescue_does_not_block_other_structures(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    def fail_optimization(call: dict[str, Any]) -> bool:
        return "ts_a" in str(call["output_dir"]) and "ts_rescue_l1" in str(call["output_dir"])

    _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[
            (-20.0, 120.0, 250.0),
            (-20.0, 120.0, 250.0),
            (-150.0, 120.0, 250.0),
            (-20.0, 120.0, 250.0),
            (-150.0, 120.0, 250.0),
        ],
        fail_optimization=fail_optimization,
    )
    _patch_ts_classifier(
        monkeypatch,
        {
            (-20.0, 120.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "soft_target_ts",
            },
            (-150.0, 120.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "valid_target_ts",
            },
        },
    )
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [
            _request(tmp_path, structure_id="ts_a", role="ts", kind="ts", seed_distance=2.5),
            _request(tmp_path, structure_id="ts_b", role="ts", kind="ts", seed_distance=2.5),
        ],
        tmp_path / "stage",
    )

    structures = {entry["id"]: entry for entry in manifest["structures"]}
    assert structures["ts_a"]["pass2_rescue_attempts"]
    assert structures["ts_b"]["pass2_rescue_attempts"]
