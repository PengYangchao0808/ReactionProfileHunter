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


def _write_failed_manifest(
    output_dir: Path,
    *,
    structure_id: str,
    role: str,
    kind: str,
    last_geometry: Path | None,
    stop_geometry: Path | None = None,
) -> None:
    payload = {
        "schema_version": "refinement_manifest_v1",
        "stage": "S3",
        "fidelity": "low",
        "profile_id": "b97_3c_r2scan_3c_v1",
        "run_id": "test-rescue-only",
        "structures": [
            {
                "id": structure_id,
                "role": role,
                "kind": kind,
                "opt_status": "failed",
                "opt_error": "ORCA terminated normally, but the geometry optimization did not converge",
                "opt_output": str(last_geometry.with_suffix(".out")) if last_geometry else None,
                "last_geometry_path": str(last_geometry) if last_geometry else None,
                "stop_geometry_path": str(stop_geometry) if stop_geometry else None,
                "imaginary_frequencies_cm1": [],
                "attempt_history": [
                    {
                        "attempt_id": "attempt_002_opt_ts",
                        "retry_level": 0,
                        "status": "failed",
                        "failure_type": "geometry_optimization_not_converged",
                        "stop_reason": "native_maxiter_checkpoint",
                        "optimization_converged": False,
                        "runtime_seconds": 100.0,
                        "directory": str(output_dir / structure_id / "diagnostics" / "attempt_002_opt_ts"),
                    }
                ],
                "current_attempt_id": "attempt_002_opt_ts",
                "primary_attempt_id": "attempt_002_opt_ts",
                "ts_classification": {
                    "hessian_index": 1,
                    "mode_identity": "unrelated",
                    "stationary_point_class": "first_order_wrong_mode",
                },
                "status": "failed",
            }
        ],
        "summary": {"complete": 0, "failed": 1},
        "provenance": {},
    }
    manifest = output_dir / "manifest.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(
        __import__("json").dumps(payload, indent=2),
        encoding="utf-8",
    )


def test_partition_rescue_replay_last_geometry_fallback(tmp_path: Path) -> None:
    engine = _engine("S3", irc_enabled=False)
    output_dir = tmp_path / "stage"
    last_xyz = _write_xyz(output_dir / "structure" / "diagnostics" / "attempt_002_opt_ts" / "last_geometry.xyz")
    _write_failed_manifest(
        output_dir,
        structure_id="ts",
        role="ts",
        kind="ts",
        last_geometry=last_xyz,
    )
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")

    replay_ids, fallback_ids = engine._partition_rescue_replay(
        [request],
        engine._load_existing_structure_payloads(output_dir),
    )

    assert replay_ids == {"ts"}
    assert fallback_ids == set()


def test_partition_rescue_replay_missing_geometry_falls_back(tmp_path: Path) -> None:
    engine = _engine("S3", irc_enabled=False)
    output_dir = tmp_path / "stage"
    missing = output_dir / "structure" / "diagnostics" / "attempt_002_opt_ts" / "last_geometry.xyz"
    _write_failed_manifest(
        output_dir,
        structure_id="ts",
        role="ts",
        kind="ts",
        last_geometry=missing,
    )
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")

    replay_ids, fallback_ids = engine._partition_rescue_replay(
        [request],
        engine._load_existing_structure_payloads(output_dir),
    )

    assert replay_ids == set()
    assert fallback_ids == {"ts"}


def test_replay_failed_outcomes_carries_rescue_fields(tmp_path: Path) -> None:
    engine = _engine("S3", irc_enabled=False)
    output_dir = tmp_path / "stage"
    last_xyz = _write_xyz(output_dir / "structure" / "diagnostics" / "attempt_002_opt_ts" / "last_geometry.xyz")
    opt_out = output_dir / "structure" / "opt" / "output.out"
    opt_out.parent.mkdir(parents=True, exist_ok=True)
    opt_out.write_text("Hessian has 1 negative eigenvalues\n", encoding="utf-8")
    _write_failed_manifest(
        output_dir,
        structure_id="ts",
        role="ts",
        kind="ts",
        last_geometry=last_xyz,
    )
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")
    engine._requests_by_id = {"ts": request}
    payloads = engine._load_existing_structure_payloads(output_dir)

    outcomes = engine._replay_failed_outcomes(payloads, {"ts"})
    assert len(outcomes) == 1
    outcome = outcomes[0]
    assert outcome.structure_id == "ts"
    assert outcome.opt_status == "failed"
    assert outcome.attempt_history[-1]["failure_type"] == "geometry_optimization_not_converged"
    assert outcome.last_geometry_path == last_xyz
    assert outcome.stop_geometry_path is None
    assert outcome.ts_classification["stationary_point_class"] == "first_order_wrong_mode"


def test_rescue_only_skips_pass1_and_rescues(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    opt_calls, _freq_calls, _sp_calls = _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[(-150.0, 120.0, 250.0)],
    )
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
    output_dir = tmp_path / "stage"
    last_xyz = _write_xyz(output_dir / "structure" / "diagnostics" / "attempt_002_opt_ts" / "last_geometry.xyz")
    _write_failed_manifest(
        output_dir,
        structure_id="ts",
        role="ts",
        kind="ts",
        last_geometry=last_xyz,
    )
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")

    manifest = engine.run(
        [request],
        output_dir,
        resume_incomplete=True,
        rescue_only=True,
    )
    record = next(
        row for row in __import__("json").loads(manifest.read_text(encoding="utf-8"))["structures"]
        if row["id"] == "ts"
    )

    assert record["status"] == "complete"
    assert record["usable_for_ml"] is True
    assert any(a["strategy"] == "fresh_hessian_restart" for a in record["pass2_rescue_attempts"])
    assert record["rescue_winning_method"] == "fresh_hessian_restart"
    assert record["canonical_attempt_id"] != "attempt_002_opt_ts"


def test_rescue_winning_method_none_when_primary_wins(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _install_qc_mocks(
        monkeypatch,
        frequency_sequence=[(50.0, 150.0, 250.0)],
    )
    _patch_ts_classifier(
        monkeypatch,
        {
            (50.0, 150.0, 250.0): {
                "hessian_index": 1,
                "mode_identity": "target",
                "stationary_point_class": "valid_target_ts",
            },
        },
    )
    engine = _engine("S3", irc_enabled=False)
    request = _request(tmp_path, structure_id="ts", role="ts", kind="ts")

    manifest = engine.run([request], tmp_path / "stage")
    record = next(
        row for row in __import__("json").loads(manifest.read_text(encoding="utf-8"))["structures"]
        if row["id"] == "ts"
    )

    assert record["status"] == "complete"
    assert record["rescue_winning_method"] is None
    assert record["pass2_rescue_attempts"] == []
