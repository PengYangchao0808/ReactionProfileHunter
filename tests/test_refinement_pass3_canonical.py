from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

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
    ensemble_g: float = -0.02,
    ensemble_h: float = -0.01,
) -> StructureRequest:
    input_xyz = _write_xyz(tmp_path / "inputs" / f"{structure_id}.xyz")
    return StructureRequest(
        id=structure_id,
        role=role,
        kind=kind,
        input_xyz=input_xyz,
        forming_bonds=[(0, 1)] if role in {"intermediate", "ts"} else [],
        s1_ensemble_thermodynamics={
            "ensemble_thermochemistry_correction_hartree": ensemble_g,
            "ensemble_enthalpy_correction_hartree": ensemble_h,
            "ensemble_correction_id": "s1-test-correction",
        },
        ensemble_thermochemistry_correction_hartree=ensemble_g,
    )


def _manifest_for(engine: RefinementEngine, requests: list[StructureRequest], output_dir: Path) -> dict:
    import json

    return json.loads(engine.run(requests, output_dir).read_text(encoding="utf-8"))


def _install_qc_mocks(
    monkeypatch: pytest.MonkeyPatch,
    *,
    frequencies: tuple[float, ...] = (50.0, 120.0, 250.0),
    freq_energy: float = -100.0,
    enthalpy: float = -99.6,
    gibbs: float = -99.5,
    sp_energy: float = -200.0,
    fail_opt: bool = False,
):
    opt_calls: list[dict[str, Any]] = []
    freq_calls: list[dict[str, Any]] = []
    sp_calls: list[dict[str, Any]] = []

    def fake_run_optimization(spec, input_xyz, output_dir, config):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        opt_calls.append({"spec": spec, "input_xyz": Path(input_xyz), "output_dir": output_dir})
        if fail_opt:
            return QCJobResult("failed", Path(input_xyz), output_file=output_dir / "opt.out", error="opt failed")
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
        out_file = output_dir / "opt.out"
        out_file.write_text("mock opt\n", encoding="utf-8")
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_xyz=out_xyz,
            output_file=out_file,
            energy_hartree=-110.0,
        )

    def fake_run_frequency(spec, input_xyz, output_dir, config):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / "freq.out"
        out_file.write_text("mock freq\n", encoding="utf-8")
        out_file.with_suffix(".hess").write_text("mock hess\n", encoding="utf-8")
        freq_calls.append({"spec": spec, "input_xyz": Path(input_xyz), "output_dir": output_dir})
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_file=out_file,
            energy_hartree=freq_energy,
            frequencies_cm1=frequencies,
            enthalpy_hartree=enthalpy,
            gibbs_free_energy_hartree=gibbs,
            gibbs_correction_hartree=gibbs - freq_energy,
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
            energy_hartree=sp_energy,
        )

    monkeypatch.setattr(engine_module, "run_optimization", fake_run_optimization)
    monkeypatch.setattr(engine_module, "run_frequency", fake_run_frequency)
    monkeypatch.setattr(engine_module, "run_single_point", fake_run_single_point)
    return opt_calls, freq_calls, sp_calls


def test_final_sp_uses_profile_sp_method(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _opt_calls, _freq_calls, sp_calls = _install_qc_mocks(monkeypatch)
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [_request(tmp_path, structure_id="product", role="product", kind="minimum")],
        tmp_path / "stage",
    )

    assert sp_calls[0]["spec"].method == engine.profile.sp_method
    assert manifest["structures"][0]["sp_status"] == "complete"


@pytest.mark.parametrize(
    ("stage", "basis", "aux_basis"),
    [("S3", "", ""), ("S4", "def2-TZVPP", "def2/J")],
)
def test_final_sp_uses_profile_sp_basis(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    stage: str,
    basis: str,
    aux_basis: str,
):
    _opt_calls, _freq_calls, sp_calls = _install_qc_mocks(monkeypatch)
    engine = _engine(stage, irc_enabled=False)
    _manifest_for(
        engine,
        [_request(tmp_path, structure_id=f"product_{stage}", role="product", kind="minimum")],
        tmp_path / f"stage_{stage}",
    )

    assert sp_calls[0]["spec"].basis == basis
    assert sp_calls[0]["spec"].aux_basis == aux_basis


def test_composite_thermochemistry_structure(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch)
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [_request(tmp_path, structure_id="product", role="product", kind="minimum")],
        tmp_path / "stage",
    )

    thermo = manifest["structures"][0]["thermochemistry"]
    expected_keys = {
        "gibbs_free_energy_hartree",
        "enthalpy_hartree",
        "geometry_level",
        "frequency_level",
        "single_point_level",
        "temperature_K",
        "standard_state",
        "qrrho",
        "ensemble_correction_id",
        "geometry_consistent",
    }
    assert expected_keys <= set(thermo)


def test_composite_gibbs_formula(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(
        monkeypatch,
        freq_energy=-100.0,
        gibbs=-99.5,
        enthalpy=-99.6,
        sp_energy=-200.0,
    )
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [_request(tmp_path, structure_id="product", role="product", kind="minimum", ensemble_g=-0.02)],
        tmp_path / "stage",
    )

    thermo = manifest["structures"][0]["thermochemistry"]
    assert thermo["gibbs_free_energy_hartree"] == pytest.approx(-200.0 + ((-99.5) - (-100.0)) - 0.02)


def test_composite_enthalpy_formula(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(
        monkeypatch,
        freq_energy=-100.0,
        gibbs=-99.5,
        enthalpy=-99.6,
        sp_energy=-200.0,
    )
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [_request(tmp_path, structure_id="product", role="product", kind="minimum", ensemble_h=-0.01)],
        tmp_path / "stage",
    )

    thermo = manifest["structures"][0]["thermochemistry"]
    assert thermo["enthalpy_hartree"] == pytest.approx(-200.0 + ((-99.6) - (-100.0)) - 0.01)


def test_ml_usability_dict_has_9_keys(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch)
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [_request(tmp_path, structure_id="product", role="product", kind="minimum")],
        tmp_path / "stage",
    )

    assert set(manifest["structures"][0]["ml_usability"]) == {
        "geometry",
        "electronic_energy",
        "enthalpy",
        "gibbs_free_energy",
        "frequency_descriptors",
        "ts_descriptors",
        "intermediate_descriptors",
        "mechanism_label",
        "multifidelity_pair",
    }


def test_ml_usability_ts_descriptors_true_for_valid_ts(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch, frequencies=(-150.0, 120.0, 250.0))
    monkeypatch.setattr(
        engine_module.identity_module,
        "classify_ts",
        lambda **kwargs: {
            "hessian_index": 1,
            "mode_identity": "target",
            "stationary_point_class": "valid_target_ts",
        },
    )
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [_request(tmp_path, structure_id="ts", role="ts", kind="ts")],
        tmp_path / "stage",
    )

    assert manifest["structures"][0]["ml_usability"]["ts_descriptors"] is True


def test_ml_usability_intermediate_descriptors_role_matched(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch)
    monkeypatch.setattr(
        engine_module.identity_module,
        "classify_int",
        lambda **kwargs: {"identity": "distinct_intermediate"},
    )
    monkeypatch.setattr(
        engine_module.identity_module,
        "classify_minimum",
        lambda **kwargs: {"identity": "valid_minimum", "imaginary_count": 0},
    )
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [_request(tmp_path, structure_id="int", role="intermediate", kind="minimum")],
        tmp_path / "stage",
    )

    assert manifest["structures"][0]["ml_usability"]["intermediate_descriptors"] is True


def test_ml_usability_mechanism_label_always_true(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch, fail_opt=True)
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [_request(tmp_path, structure_id="broken", role="product", kind="minimum")],
        tmp_path / "stage",
    )

    entry = manifest["structures"][0]
    assert entry["canonical_attempt_id"] is None
    assert entry["ml_usability"]["mechanism_label"] is True


def test_canonical_selection_picks_best_attempt():
    engine = _engine("S3", irc_enabled=False)
    outcome = Pass1Outcome(
        structure_id="ts",
        role="ts",
        kind="ts",
        opt_status="complete",
        frequency_status="complete",
        primary_attempt_id="primary",
        ts_classification={
            "stationary_point_class": "soft_target_ts",
            "mode_identity": "target",
            "hessian_index": 1,
        },
    )
    outcome.pass2_rescue_attempts = [
        {
            "attempt_id": "rescue_wrong",
            "opt_status": "complete",
            "frequency_status": "complete",
            "ts_classification": {
                "stationary_point_class": "first_order_wrong_mode",
                "mode_identity": "unrelated",
                "hessian_index": 1,
            },
            "alignment_score": 0.1,
            "gradient_norm": 0.2,
        },
        {
            "attempt_id": "rescue_valid",
            "opt_status": "complete",
            "frequency_status": "complete",
            "ts_classification": {
                "stationary_point_class": "valid_target_ts",
                "mode_identity": "target",
                "hessian_index": 1,
            },
            "alignment_score": 0.9,
            "gradient_norm": 0.1,
        },
    ]

    assert engine._select_canonical_attempt(outcome) == "rescue_valid"


def test_canonical_none_when_no_candidates():
    engine = _engine("S3", irc_enabled=False)
    outcome = Pass1Outcome(structure_id="x", role="product", kind="minimum", opt_status="failed")

    assert engine._select_canonical_attempt(outcome) is None


def test_manifest_has_all_pass3_fields(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch)
    engine = _engine("S3", irc_enabled=False)
    manifest = _manifest_for(
        engine,
        [_request(tmp_path, structure_id="product", role="product", kind="minimum")],
        tmp_path / "stage",
    )

    entry = manifest["structures"][0]
    required = {
        "pass2_rescue_attempts",
        "canonical_attempt_id",
        "canonical_xyz",
        "geometry_hash",
        "canonical_frequency_status",
        "sp_status",
        "thermochemistry",
        "ml_usability",
        "irc_status",
        "resolved_kind",
        "resolved_identity",
        "identity_status",
        "usable_for_ml",
    }
    assert required <= set(entry)
