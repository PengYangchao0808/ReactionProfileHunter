"""Tests for the S3 rescue method-family matrix (v3.1).

Covers matrix exhaustiveness, cell ordering, method parameter resolution,
F1-F6 failure classification, and the matrix-driven method executors.
"""

from __future__ import annotations

from collections import deque
from dataclasses import replace
from pathlib import Path
from typing import Any, Callable

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
) -> StructureRequest:
    inputs = tmp_path / "inputs"
    input_xyz = _write_xyz(inputs / f"{structure_id}.xyz")
    return StructureRequest(
        id=structure_id,
        role=role,
        kind=kind,
        input_xyz=input_xyz,
        source_stage="S2",
        forming_bonds=[(0, 1)] if role in {"intermediate", "ts"} else [],
    )


def _failed_outcome(
    *,
    structure_id: str,
    role: str,
    kind: str,
    opt_output: Path | None = None,
    stop_geometry: Path | None = None,
    failure_type: str | None = "geometry_optimization_not_converged",
    error: str | None = None,
) -> Pass1Outcome:
    return Pass1Outcome(
        structure_id=structure_id,
        role=role,
        kind=kind,
        opt_status="failed",
        opt_output=opt_output,
        opt_error=error or (
            "ORCA terminated normally, but the geometry optimization did not converge"
            if failure_type == "geometry_optimization_not_converged"
            else error
        ),
        stop_geometry_path=stop_geometry,
        attempt_history=[{"failure_type": failure_type, "stop_reason": None}],
    )


def _install_qc_mocks(
    monkeypatch: pytest.MonkeyPatch,
    *,
    frequency_sequence: list[tuple[float, ...]],
):
    opt_calls: list[dict[str, Any]] = []
    freq_calls: list[dict[str, Any]] = []
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
        freq_calls.append({"spec": spec, "input_xyz": Path(input_xyz), "frequencies": freqs})
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_file=out_file,
            energy_hartree=-0.9,
            frequencies_cm1=freqs,
        )

    monkeypatch.setattr(engine_module, "run_optimization", fake_run_optimization)
    monkeypatch.setattr(engine_module, "run_frequency", fake_run_frequency)
    return opt_calls, freq_calls


def _patch_ts_classifier(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_classify_ts(**kwargs):
        return {
            "hessian_index": 1,
            "mode_identity": "target",
            "stationary_point_class": "valid_target_ts",
        }

    monkeypatch.setattr(engine_module.identity_module, "classify_ts", fake_classify_ts)


# ---------------------------------------------------------------------------
# Matrix structure
# ---------------------------------------------------------------------------


def test_matrix_cells_are_exhaustive() -> None:
    plans = all_plans()
    assert len(plans) == 8
    cells = {(plan.failure_type, plan.kind) for plan in plans}
    assert (FailureType.F1_GEOMETRY_NOT_CONVERGED, StructureKind.TS) in cells
    assert (FailureType.F1_GEOMETRY_NOT_CONVERGED, StructureKind.INT) in cells
    assert (FailureType.F1_GEOMETRY_NOT_CONVERGED, StructureKind.MINIMUM) in cells
    assert (FailureType.F2_HIGHER_ORDER_SADDLE, StructureKind.TS) in cells
    assert (FailureType.F3_TS_NO_IMAGINARY_MODE, StructureKind.TS) in cells
    assert (FailureType.F4_MINIMUM_WITH_IMAGINARY, StructureKind.INT) in cells
    assert (FailureType.F4_MINIMUM_WITH_IMAGINARY, StructureKind.MINIMUM) in cells
    assert (FailureType.F7_COLLAPSED_TO_PRODUCT, StructureKind.INT) in cells


def test_matrix_f5_f6_have_no_cells() -> None:
    assert lookup_plan(FailureType.F5_SCF_NOT_CONVERGED, StructureKind.TS) is None
    assert lookup_plan(FailureType.F6_CRASH_TIMEOUT_OTHER, StructureKind.INT) is None
    assert lookup_plan(FailureType.F6_CRASH_TIMEOUT_OTHER, StructureKind.MINIMUM) is None


def test_matrix_cell_orders() -> None:
    plan = lookup_plan(FailureType.F2_HIGHER_ORDER_SADDLE, StructureKind.TS)
    assert plan is not None
    assert plan.methods[0] is RescueMethod.SADDLE_BREAK
    plan = lookup_plan(FailureType.F1_GEOMETRY_NOT_CONVERGED, StructureKind.TS)
    assert plan is not None
    assert [m.value for m in plan.methods] == [
        "fresh_hessian_restart",
        "ts_mode_directed",
        "calcall_opt",
    ]
    plan = lookup_plan(FailureType.F1_GEOMETRY_NOT_CONVERGED, StructureKind.INT)
    assert plan is not None
    assert [m.value for m in plan.methods] == ["fresh_hessian_restart", "calcall_opt"]
    plan = lookup_plan(FailureType.F4_MINIMUM_WITH_IMAGINARY, StructureKind.MINIMUM)
    assert plan is not None
    assert [m.value for m in plan.methods] == ["mode_displacement"]


def test_structure_kind_from_request() -> None:
    assert StructureKind.from_request("ts", "ts") is StructureKind.TS
    assert StructureKind.from_request("minimum", "intermediate") is StructureKind.INT
    assert StructureKind.from_request("minimum", "precursor") is StructureKind.MINIMUM


# ---------------------------------------------------------------------------
# Method parameters (v3.1 final values)
# ---------------------------------------------------------------------------


def test_method_params_final_values() -> None:
    fresh = methods_params(None, RescueMethod.FRESH_HESSIAN_RESTART)
    assert fresh.calc_hess is True
    assert fresh.recalc_hessian == 5
    assert fresh.trust is None
    directed = methods_params(None, RescueMethod.TS_MODE_DIRECTED)
    assert directed.trust == pytest.approx(0.15)
    assert directed.max_cycles == 12
    calcall = methods_params(None, RescueMethod.CALCALL_OPT)
    assert calcall.recalc_hessian == 1
    assert calcall.trust is None


def test_method_params_config_override() -> None:
    config = {
        "refinement": {
            "common": {
                "rescue": {
                    "methods": {
                        "ts_mode_directed": {"trust": 0.2, "max_cycles": 8},
                        "calcall_opt": {"recalc_hessian": 2},
                    }
                }
            }
        }
    }
    directed = methods_params(config, RescueMethod.TS_MODE_DIRECTED)
    assert directed.trust == pytest.approx(0.2)
    assert directed.max_cycles == 8
    assert methods_params(config, RescueMethod.CALCALL_OPT).recalc_hessian == 2


# ---------------------------------------------------------------------------
# Failure classification
# ---------------------------------------------------------------------------


def test_classify_f1_geometry_not_converged(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="int", role="intermediate", kind="minimum")
    outcome = _failed_outcome(structure_id="int", role="intermediate", kind="minimum")
    assert (
        engine._classify_rescue_failure(request, outcome)
        is FailureType.F1_GEOMETRY_NOT_CONVERGED
    )


def test_classify_f5_scf_not_converged(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="int", role="intermediate", kind="minimum")
    outcome = _failed_outcome(
        structure_id="int",
        role="intermediate",
        kind="minimum",
        failure_type="scf_not_converged",
    )
    assert engine._classify_rescue_failure(request, outcome) is FailureType.F5_SCF_NOT_CONVERGED


def test_classify_f6_timed_out(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")
    outcome = Pass1Outcome(
        structure_id="ts",
        role="ts",
        kind="ts",
        opt_status="failed",
        attempt_history=[{"failure_type": None, "stop_reason": "timed_out"}],
    )
    assert engine._classify_rescue_failure(request, outcome) is FailureType.F6_CRASH_TIMEOUT_OTHER


def test_classify_f2_higher_order_saddle_from_output(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")
    opt_output = tmp_path / "opt.out"
    opt_output.write_text(
        "Hessian has   2 negative eigenvalues\n...\nHessian has   2 negative eigenvalues\n",
        encoding="utf-8",
    )
    outcome = _failed_outcome(
        structure_id="ts",
        role="ts",
        kind="ts",
        opt_output=opt_output,
    )
    assert (
        engine._classify_rescue_failure(request, outcome)
        is FailureType.F2_HIGHER_ORDER_SADDLE
    )


def test_classify_f3_ts_without_imaginary_mode(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")
    opt_output = tmp_path / "opt.out"
    opt_output.write_text("Hessian has   0 negative eigenvalues\n", encoding="utf-8")
    outcome = _failed_outcome(
        structure_id="ts",
        role="ts",
        kind="ts",
        opt_output=opt_output,
    )
    assert (
        engine._classify_rescue_failure(request, outcome)
        is FailureType.F3_TS_NO_IMAGINARY_MODE
    )


def test_classify_ts_uses_last_hessian_state_not_max(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")
    opt_output = tmp_path / "opt.out"
    # trajectory wandered through a 2nd-order region but ENDED as 1st-order:
    # must classify F1 (slow convergence), not F2.
    opt_output.write_text(
        "Hessian has   1 negative eigenvalues\n"
        "Hessian has   2 negative eigenvalues\n"
        "Hessian has   2 negative eigenvalues\n"
        "Hessian has   1 negative eigenvalues\n",
        encoding="utf-8",
    )
    outcome = _failed_outcome(
        structure_id="ts",
        role="ts",
        kind="ts",
        opt_output=opt_output,
    )
    assert (
        engine._classify_rescue_failure(request, outcome)
        is FailureType.F1_GEOMETRY_NOT_CONVERGED
    )


def test_classify_f4_converged_minimum_with_imaginary(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="precursor", role="precursor", kind="minimum")
    outcome = Pass1Outcome(
        structure_id="precursor",
        role="precursor",
        kind="minimum",
        opt_status="complete",
        imaginary_frequencies_cm1=[-40.0],
    )
    assert (
        engine._classify_rescue_failure(request, outcome)
        is FailureType.F4_MINIMUM_WITH_IMAGINARY
    )


def test_classify_f4_ignores_noise_imaginary_within_cutoff(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="precursor", role="precursor", kind="minimum")
    outcome = Pass1Outcome(
        structure_id="precursor",
        role="precursor",
        kind="minimum",
        opt_status="complete",
        imaginary_frequencies_cm1=[-4.7],
    )
    assert engine._classify_rescue_failure(request, outcome) is None


def test_saddle_eigenvalue_parsing(tmp_path: Path) -> None:
    engine = _engine("S3")
    opt_output = tmp_path / "opt.out"
    opt_output.write_text(
        "Hessian has   1 negative eigenvalues\n"
        "Hessian has   2 negative eigenvalues\n"
        "Hessian has   2 negative eigenvalues\n",
        encoding="utf-8",
    )
    outcome = _failed_outcome(
        structure_id="ts", role="ts", kind="ts", opt_output=opt_output
    )
    assert engine._parse_saddle_eigenvalue_counts(outcome) == [1, 2, 2]


# ---------------------------------------------------------------------------
# Matrix-driven execution
# ---------------------------------------------------------------------------


def test_f2_ts_cell_runs_saddle_break_first(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    opt_calls, _freq_calls = _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[
            (-80.0, -30.0, 120.0),  # saddle_break local Hessian (2 imag modes)
            (-80.0, -30.0, 120.0),  # saddle_break candidate validation
        ],
    )
    _patch_ts_classifier(monkeypatch)
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")
    stop_xyz = _write_xyz(tmp_path / "stop.xyz", distance=2.2)
    opt_output = tmp_path / "opt.out"
    opt_output.write_text("Hessian has   2 negative eigenvalues\n", encoding="utf-8")
    outcome = _failed_outcome(
        structure_id="ts",
        role="ts",
        kind="ts",
        opt_output=opt_output,
        stop_geometry=stop_xyz,
    )
    failure = engine._classify_rescue_failure(request, outcome)
    assert failure is FailureType.F2_HIGHER_ORDER_SADDLE
    plan = lookup_plan(failure, StructureKind.from_request(request.kind, request.role))
    assert plan is not None

    engine._run_rescue_cell(
        request,
        outcome,
        AttemptRecorder(tmp_path / "stage" / "ts", structure_id="ts"),
        tmp_path / "rescue",
        plan,
    )

    assert len(opt_calls) >= 1
    first_rescue = opt_calls[0]
    assert "r2_saddle_break" in str(first_rescue["input_xyz"])
    assert first_rescue["spec"].task == "opt_ts"
    assert first_rescue["spec"].initial_hessian == "calculate"
    strategies = [
        str(record.get("strategy")) for record in outcome.pass2_rescue_attempts
    ]
    assert any("saddle_break" in strategy for strategy in strategies)


def test_f1_int_cell_runs_fresh_hessian_restart(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    opt_calls, _freq_calls = _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[(50.0, 120.0, 250.0)],
    )
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="int", role="intermediate", kind="minimum")
    stop_xyz = _write_xyz(tmp_path / "stop_int.xyz", distance=2.2)
    outcome = _failed_outcome(
        structure_id="int",
        role="intermediate",
        kind="minimum",
        stop_geometry=stop_xyz,
    )
    failure = engine._classify_rescue_failure(request, outcome)
    assert failure is FailureType.F1_GEOMETRY_NOT_CONVERGED
    plan = lookup_plan(failure, StructureKind.from_request(request.kind, request.role))
    assert plan is not None

    engine._run_rescue_cell(
        request,
        outcome,
        AttemptRecorder(tmp_path / "stage" / "int", structure_id="int"),
        tmp_path / "rescue",
        plan,
    )

    assert opt_calls, "fresh_hessian_restart must run an optimization"
    first = opt_calls[0]
    assert first["input_xyz"] == stop_xyz
    assert first["spec"].initial_hessian == "calculate"
    assert first["spec"].recalc_hessian == 5
    assert first["spec"].trust is None


def test_f5_scf_failure_exits_directly_without_rescue(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="int", role="intermediate", kind="minimum")
    outcome = _failed_outcome(
        structure_id="int",
        role="intermediate",
        kind="minimum",
        failure_type="scf_not_converged",
    )
    failure = engine._classify_rescue_failure(request, outcome)
    assert failure is FailureType.F5_SCF_NOT_CONVERGED
    assert lookup_plan(failure, StructureKind.from_request(request.kind, request.role)) is None


# ---------------------------------------------------------------------------
# Mode displacement with explicit mode index (saddle_break uses mode 1)
# ---------------------------------------------------------------------------


def test_build_mode_displaced_xyz_selects_second_mode(tmp_path: Path) -> None:
    engine = _engine("S3")
    source = _write_xyz(tmp_path / "source.xyz", distance=2.0)
    freq_out = tmp_path / "freq.out"
    freq_out.write_text(
        "CARTESIAN DISPLACEMENTS\n-----------------------\n"
        "Mode:   0\nFreq:   -80.00 cm**-1 (imaginary mode)\n"
        "        dx          dy          dz\n"
        "Atom 0:  1.000000   0.000000   0.000000\n"
        "Atom 1: -1.000000   0.000000   0.000000\n"
        "Mode:   1\nFreq:   -30.00 cm**-1 (imaginary mode)\n"
        "        dx          dy          dz\n"
        "Atom 0:  0.000000   1.000000   0.000000\n"
        "Atom 1:  0.000000  -1.000000   0.000000\n",
        encoding="utf-8",
    )
    displaced = engine._build_mode_displaced_xyz(
        source, freq_out, tmp_path / "displaced.xyz", sign=1.0, mode_index=1
    )
    content = displaced.read_text(encoding="utf-8").splitlines()
    atom0 = content[2].split()[1:4]
    assert abs(float(atom0[0])) < 1e-6  # mode 1 displaces along y, not x
    assert abs(float(atom0[1])) > 1e-3


# ---------------------------------------------------------------------------
# steps[] payload contract
# ---------------------------------------------------------------------------


def test_steps_payload_contains_s3_0_to_s3_7(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="precursor", role="precursor", kind="minimum")
    from rph_core.steps.refinement.models import PreflightOutcome

    preflight = PreflightOutcome(
        structure_id="precursor", status="ok", input_xyz=request.input_xyz, charge=0, multiplicity=1
    )
    pass1 = Pass1Outcome(
        structure_id="precursor",
        role="precursor",
        kind="minimum",
        warmup_status="not_requested",
        opt_status="complete",
        frequency_status="complete",
        primary_attempt_id="attempt_001_opt",
        sp_status="complete",
        sp_energy_hartree=-1.2,
    )
    pass2 = Pass1Outcome(structure_id="precursor", role="precursor", kind="minimum")
    pass3 = Pass1Outcome(
        structure_id="precursor",
        role="precursor",
        kind="minimum",
        canonical_xyz=tmp_path / "canonical.xyz",
        canonical_attempt_id="attempt_003_canonical",
        sp_status="complete",
        sp_energy_hartree=-1.2,
        status="complete",
    )

    steps = engine._build_steps_payload(request, preflight, pass1, pass2, pass3)

    assert [step["code"] for step in steps] == [
        "S3.0",
        "S3.1",
        "S3.2",
        "S3.3",
        "S3.4",
        "S3.5",
        "S3.6",
        "S3.7",
    ]
    assert steps[0]["status"] == "complete"
    assert steps[1]["status"] == "skipped"  # precursor has no warmup
    assert steps[2]["status"] == "complete"
    assert steps[3]["status"] == "complete"
    assert steps[4]["status"] == "not_run"  # no rescue needed
    assert steps[5]["status"] == "complete"
    assert steps[6]["status"] == "complete"
    assert steps[6]["sp_energy_hartree"] == pytest.approx(-1.2)
    assert steps[7]["status"] == "complete"


def test_steps_payload_rescue_sub_levels(tmp_path: Path) -> None:
    engine = _engine("S3")
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")
    from rph_core.steps.refinement.models import PreflightOutcome

    preflight = PreflightOutcome(structure_id="ts", status="ok", input_xyz=request.input_xyz)
    pass1 = Pass1Outcome(
        structure_id="ts",
        role="ts",
        kind="ts",
        opt_status="failed",
        frequency_status="not_run",
    )
    pass2 = Pass1Outcome(
        structure_id="ts",
        role="ts",
        kind="ts",
        pass2_rescue_attempts=[
            {
                "strategy": "S3.4.0_diagnosis",
                "failure_type": "F2",
                "status": "diagnosed",
            },
            {
                "strategy": "saddle_break_plus",
                "method_family": "R2",
                "method": "saddle_break",
                "status": "complete",
                "attempt_id": "attempt_004_r2_saddle_break",
            },
        ],
    )
    pass3 = Pass1Outcome(
        structure_id="ts",
        role="ts",
        kind="ts",
        canonical_xyz=tmp_path / "canonical.xyz",
        status="complete",
    )

    steps = engine._build_steps_payload(request, preflight, pass1, pass2, pass3)
    rescue_step = steps[4]
    assert rescue_step["code"] == "S3.4"
    assert rescue_step["status"] == "partial"  # method complete but no valid candidate in fixture
    sub_codes = [sub["code"] for sub in rescue_step["sub_steps"]]
    assert "S3.4.0" in sub_codes
    assert "S3.4.2" in sub_codes  # R2 family
