"""Tests for the IRC mid-point INT rescue (FailureType F7, R4 family)."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

import rph_core.steps.refinement.engine as engine_module
from rph_core.steps.fidelity_profile import FidelityProfile
from rph_core.steps.refinement import RefinementEngine
from rph_core.steps.refinement.models import Pass1Outcome, StructureRequest
from rph_core.steps.refinement.rescue_matrix import (
    FailureType,
    RescueMethod,
    StructureKind,
    all_plans,
    lookup_plan,
    methods_params,
)
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


def _int_request(tmp_path: Path) -> StructureRequest:
    inputs = tmp_path / "inputs"
    input_xyz = _write_xyz(inputs / "product_minor_001_int.xyz")
    return StructureRequest(
        id="product_minor_001_int",
        role="intermediate",
        kind="minimum",
        input_xyz=input_xyz,
        variant_id="product_minor_001",
        source_stage="S2",
        forming_bonds=[(0, 1)],
        charge=0,
        multiplicity=1,
    )


def _collapsed_int_outcome(tmp_path: Path) -> Pass1Outcome:
    return Pass1Outcome(
        structure_id="product_minor_001_int",
        role="intermediate",
        kind="minimum",
        opt_status="complete",
        status="complete",
        int_classification={"identity": "collapsed_to_product"},
    )


def _dual_end_irc_xyz(tmp_path: Path, n_frames: int = 5) -> tuple[Path, Path]:
    """Build both IRC trajectories: Fwd stretches to ~3.5 A, Bwd contracts to ~1.5 A."""
    fwd = tmp_path / "job_IRC_Fwd_Traj.xyz"
    bwd = tmp_path / "job_IRC_Bwd_Traj.xyz"
    fwd_lines = []
    bwd_lines = []
    base_e = -1052.70
    for frame in range(n_frames):
        fwd_lines.append("2")
        fwd_lines.append(f"# energy {base_e - 0.001 * frame:.6f} E_h")
        fwd_lines.append("H 0.00000000 0.00000000 0.00000000")
        fwd_lines.append(f"H {2.8 + 0.15 * frame:.8f} 0.00000000 0.00000000")
        bwd_lines.append("2")
        bwd_lines.append(f"# energy {base_e - 0.02 * frame:.6f} E_h")
        bwd_lines.append("H 0.00000000 0.00000000 0.00000000")
        bwd_lines.append(f"H {2.0 - 0.1 * frame:.8f} 0.00000000 0.00000000")
    fwd.write_text("\n".join(fwd_lines) + "\n", encoding="utf-8")
    bwd.write_text("\n".join(bwd_lines) + "\n", encoding="utf-8")
    return fwd, bwd


def _install_irc_mock(monkeypatch: pytest.MonkeyPatch, fwd: Path, bwd: Path):
    calls: list[dict[str, Any]] = []

    def fake_run_irc(spec, input_xyz, output_dir, config, subprocess_callback=None, charge=0, spin=1):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        calls.append(
            {
                "spec": spec,
                "input_xyz": Path(input_xyz),
                "charge": charge,
                "spin": spin,
            }
        )
        (output_dir / "job_IRC_Fwd_Traj.xyz").write_text(
            fwd.read_text(encoding="utf-8"), encoding="utf-8"
        )
        (output_dir / "job_IRC_Bwd_Traj.xyz").write_text(
            bwd.read_text(encoding="utf-8"), encoding="utf-8"
        )
        endpoint_a = output_dir / "irc_endpoint_a.xyz"
        endpoint_a.write_text("2\nendpoint\nH 0 0 0\nH 3.5 0 0\n", encoding="utf-8")
        endpoint_b = output_dir / "irc_endpoint_b.xyz"
        endpoint_b.write_text("2\nendpoint\nH 0 0 0\nH 1.5 0 0\n", encoding="utf-8")
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_file=output_dir / "irc.out",
            extra={
                "endpoint_a_xyz": endpoint_a,
                "endpoint_b_xyz": endpoint_b,
            },
        )

    monkeypatch.setattr(engine_module, "run_irc", fake_run_irc)
    return calls


def _install_opt_mock(monkeypatch: pytest.MonkeyPatch, distance: float = 2.9):
    def fake_run_optimization(spec, input_xyz, output_dir, config):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(
            f"2\nmock\nH 0.0 0.0 0.0\nH {distance:.6f} 0.0 0.0\n",
            encoding="utf-8",
        )
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_xyz=out_xyz,
            output_file=output_dir / "opt.out",
            energy_hartree=-1.0,
        )

    monkeypatch.setattr(engine_module, "run_optimization", fake_run_optimization)


def _install_freq_mock(monkeypatch: pytest.MonkeyPatch, frequencies: tuple[float, ...] = (50.0, 150.0)):
    def fake_run_frequency(spec, input_xyz, output_dir, config):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / "freq.out"
        out_file.write_text("mock freq\n", encoding="utf-8")
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_file=out_file,
            energy_hartree=-0.9,
            frequencies_cm1=list(frequencies),
        )

    monkeypatch.setattr(engine_module, "run_frequency", fake_run_frequency)


def test_f7_classification_collapsed_int() -> None:
    engine = _engine("S3")
    request = _int_request(Path("/tmp"))
    outcome = _collapsed_int_outcome(Path("/tmp"))

    failure = engine._classify_rescue_failure(request, outcome)

    assert failure is FailureType.F7_COLLAPSED_TO_PRODUCT


def test_f7_not_classified_for_intact_int() -> None:
    engine = _engine("S3")
    request = _int_request(Path("/tmp"))
    outcome = _collapsed_int_outcome(Path("/tmp"))
    outcome.int_classification = {"identity": "distinct_intermediate"}

    assert engine._classify_rescue_failure(request, outcome) is None


def test_matrix_has_f7_irc_cell() -> None:
    plan = lookup_plan(FailureType.F7_COLLAPSED_TO_PRODUCT, StructureKind.INT)

    assert plan is not None
    assert plan.methods == (RescueMethod.IRC_MIDPOINT_RECOVERY,)
    assert any(
        cell.failure_type is FailureType.F7_COLLAPSED_TO_PRODUCT
        and RescueMethod.IRC_MIDPOINT_RECOVERY in cell.methods
        for cell in all_plans()
    )


def test_irc_method_params_defaults() -> None:
    engine = _engine("S3")

    params = methods_params(engine.config, RescueMethod.IRC_MIDPOINT_RECOVERY)

    assert params.irc_max_iter == 5
    assert params.irc_direction == "both"
    assert params.shoulder_energy_window_kcal_mol == pytest.approx(2.0)


def test_rescue_irc_recovers_stretched_intermediate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fwd, bwd = _dual_end_irc_xyz(tmp_path)
    calls = _install_irc_mock(monkeypatch, fwd, bwd)
    _install_opt_mock(monkeypatch, distance=2.9)
    _install_freq_mock(monkeypatch, frequencies=(50.0, 150.0))
    engine = _engine("S3")
    request = _int_request(tmp_path)
    request.parent_structure = {
        "product_ref": str(_write_xyz(tmp_path / "product_ref.xyz", distance=1.55)),
        "precursor_ref": str(_write_xyz(tmp_path / "precursor_ref.xyz", distance=3.5)),
    }
    outcome = _collapsed_int_outcome(tmp_path)
    ts_outcome = Pass1Outcome(
        structure_id="product_minor_001_ts",
        role="ts",
        kind="ts",
        canonical_xyz=_write_xyz(tmp_path / "ts_canonical.xyz"),
        canonical_hessian_path=tmp_path / "ts.hess",
        opt_energy_hartree=-0.99,
    )
    (tmp_path / "ts.hess").write_text("$vibrational_frequencies\n0  -54.0\n", encoding="utf-8")
    engine._variant_ts_outcome = {"product_minor_001": ts_outcome}
    engine._variant_product_outcome = {
        "product_minor_001": Pass1Outcome(
            structure_id="product_minor_001",
            role="product",
            kind="minimum",
            canonical_xyz=_write_xyz(tmp_path / "product_canonical.xyz", distance=1.55),
            opt_energy_hartree=-1.5,
        )
    }
    structure_dir = tmp_path / "structure_dir"
    structure_dir.mkdir(parents=True, exist_ok=True)
    recorder = engine_module.AttemptRecorder(structure_dir, structure_id=request.id)

    params = methods_params(engine.config, RescueMethod.IRC_MIDPOINT_RECOVERY)
    record = engine._rescue_int_via_irc(
        request, outcome, recorder, structure_dir, params
    )

    assert record is not None
    assert record["status"] == "complete"
    assert record["mechanism"] == "intermediate_recovered"
    assert outcome.opt_status == "complete"
    assert outcome.irc_status == "complete"
    assert outcome.irc_endpoints["mechanism"] == "intermediate_recovered"
    assert outcome.irc_endpoints["seed_source"] == "irc_stretched_endpoint"
    assert (structure_dir / "irc_stretch_seed.xyz").exists()
    assert calls and calls[0]["spec"].max_iter == 5
    assert calls[0]["spec"].direction == "both"
    assert calls[0]["charge"] == 0
    assert calls[0]["spin"] == 1


def test_rescue_irc_stretch_endpoint_not_stationary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fwd, bwd = _dual_end_irc_xyz(tmp_path)
    _install_irc_mock(monkeypatch, fwd, bwd)
    _install_opt_mock(monkeypatch, distance=1.55)
    _install_freq_mock(monkeypatch, frequencies=(50.0, 150.0))
    engine = _engine("S3")
    request = _int_request(tmp_path)
    request.parent_structure = {
        "product_ref": str(_write_xyz(tmp_path / "product_ref.xyz", distance=1.55)),
        "precursor_ref": str(_write_xyz(tmp_path / "precursor_ref.xyz", distance=3.5)),
    }
    outcome = _collapsed_int_outcome(tmp_path)
    ts_outcome = Pass1Outcome(
        structure_id="product_minor_001_ts",
        role="ts",
        kind="ts",
        canonical_xyz=_write_xyz(tmp_path / "ts_canonical.xyz"),
        canonical_hessian_path=tmp_path / "ts.hess",
    )
    (tmp_path / "ts.hess").write_text("$vibrational_frequencies\n0  -54.0\n", encoding="utf-8")
    engine._variant_ts_outcome = {"product_minor_001": ts_outcome}
    structure_dir = tmp_path / "structure_dir"
    structure_dir.mkdir(parents=True, exist_ok=True)
    recorder = engine_module.AttemptRecorder(structure_dir, structure_id=request.id)

    params = methods_params(engine.config, RescueMethod.IRC_MIDPOINT_RECOVERY)
    record = engine._rescue_int_via_irc(
        request, outcome, recorder, structure_dir, params
    )

    assert record is not None
    assert record["mechanism"] == "stretch_endpoint_not_stationary"
    assert outcome.irc_status == "complete"
    assert outcome.irc_endpoints["mechanism"] == "stretch_endpoint_not_stationary"


def test_classify_irc_endpoints_identifies_stretch_side(tmp_path: Path) -> None:
    from rph_core.utils.irc_trajectory import (
        classify_irc_endpoints,
        parse_irc_trajectory_file,
    )

    fwd, bwd = _dual_end_irc_xyz(tmp_path)
    trajectories = [
        parse_irc_trajectory_file(fwd),
        parse_irc_trajectory_file(bwd),
    ]

    result = classify_irc_endpoints(trajectories, [(0, 1)])

    assert result is not None
    assert result["stretch_mean_fb_distance"] > 3.0
    assert result["product_mean_fb_distance"] < 1.8
    assert result["stretch_mean_fb_distance"] > result["product_mean_fb_distance"]


def _collapsed_int_payload(tmp_path: Path, structure_id: str = "product_minor_001_int") -> dict:
    canonical = _write_xyz(tmp_path / f"{structure_id}_canonical.xyz", distance=1.55)
    return {
        "id": structure_id,
        "role": "intermediate",
        "kind": "minimum",
        "status": "complete",
        "opt_status": "complete",
        "canonical_xyz": str(canonical),
        "opt_output": str(canonical),
        "imaginary_frequencies_cm1": [],
        "int_classification": {"identity": "distinct_intermediate"},
    }


def test_preserved_collapsed_int_rerouted_to_irc(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fwd, bwd = _dual_end_irc_xyz(tmp_path)
    _install_irc_mock(monkeypatch, fwd, bwd)
    engine = _engine("S3")
    request = _int_request(tmp_path)
    request.parent_structure = {
        "product_ref": str(_write_xyz(tmp_path / "product_ref.xyz", distance=1.55)),
        "precursor_ref": str(_write_xyz(tmp_path / "precursor_ref.xyz", distance=3.5)),
    }
    payloads = {"product_minor_001_int": _collapsed_int_payload(tmp_path)}
    ts_payload = {
        "id": "product_minor_001_ts",
        "role": "ts",
        "kind": "ts",
        "status": "complete",
        "canonical_xyz": str(_write_xyz(tmp_path / "ts_canonical.xyz")),
        "canonical_hessian_path": str(tmp_path / "ts.hess"),
    }
    (tmp_path / "ts.hess").write_text("$vibrational_frequencies\n0  -54.0\n", encoding="utf-8")
    payloads["product_minor_001_ts"] = ts_payload
    engine._requests_by_id = {
        "product_minor_001_int": request,
        "product_minor_001_ts": StructureRequest(
            id="product_minor_001_ts", role="ts", kind="ts",
            input_xyz=_write_xyz(tmp_path / "ts_input.xyz"),
            variant_id="product_minor_001",
        ),
    }

    collapsed = engine._preserved_int_collapse_ids(
        list(engine._requests_by_id.values()), payloads, set()
    )

    assert collapsed == {"product_minor_001_int"}


def test_preserved_intact_int_not_rerouted(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _int_request(tmp_path)
    request.parent_structure = {
        "product_ref": str(_write_xyz(tmp_path / "product_ref.xyz", distance=1.55)),
        "precursor_ref": str(_write_xyz(tmp_path / "precursor_ref.xyz", distance=3.5)),
    }
    canonical = _write_xyz(tmp_path / "intact_canonical.xyz", distance=2.9)
    payloads = {
        "product_minor_001_int": {
            "id": "product_minor_001_int",
            "role": "intermediate",
            "kind": "minimum",
            "status": "complete",
            "opt_status": "complete",
            "canonical_xyz": str(canonical),
            "opt_output": str(canonical),
            "imaginary_frequencies_cm1": [],
            "int_classification": {"identity": "distinct_intermediate"},
        }
    }
    engine._requests_by_id = {"product_minor_001_int": request}

    collapsed = engine._preserved_int_collapse_ids(
        list(engine._requests_by_id.values()), payloads, set()
    )

    assert collapsed == set()


def test_reusable_ts_with_missing_freq_output_falls_back_to_hessian(tmp_path: Path) -> None:
    engine = _engine("S3")
    canonical = _write_xyz(tmp_path / "ts_canonical.xyz")
    hess = tmp_path / "ts.hess"
    hess.write_text("$vibrational_frequencies\n0  -54.0\n", encoding="utf-8")
    payload = {
        "id": "product_major_ts",
        "role": "ts",
        "kind": "ts",
        "status": "complete",
        "opt_status": "failed",
        "canonical_xyz": str(canonical),
        "opt_output": str(canonical),
        "canonical_frequency_status": "complete",
        "canonical_frequency_output": None,
        "canonical_hessian_path": str(hess),
        "sp_status": "complete",
        "sp_output": str(canonical),
    }

    assert engine._existing_structure_is_reusable(payload) is True
