import json
from dataclasses import replace
from pathlib import Path
from typing import cast

import pytest

import rph_core.steps.refinement.engine as engine_module
from rph_core.steps.fidelity_profile import FidelityProfile
from rph_core.steps.refinement import RefinementEngine
from rph_core.steps.step3_lowlevel import LowLevelEngine
from rph_core.steps.step4_highlevel import HighLevelEngine
from rph_core.utils.config_loader import load_config
from rph_core.utils.qc_models import QCJobResult, QCJobSpec


def _config() -> dict:
    return load_config()


def _engine(stage: str, **profile_overrides: object) -> RefinementEngine:
    config = _config()
    profile = FidelityProfile.from_config(config, stage)
    if profile_overrides:
        profile = replace(profile, **profile_overrides)
    return RefinementEngine(config, profile)


def _write_xyz(path: Path, distance: float = 2.4) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"2\nmock\nH 0.0 0.0 0.0\nH {distance:.6f} 0.0 0.0\n",
        encoding="utf-8",
    )
    return path


def _structure(
    tmp_path: Path,
    *,
    structure_id: str,
    role: str,
    kind: str,
    input_distance: float = 2.4,
    seed_distance: float | None = None,
    source_stage: str = "S2",
    **extra: object,
) -> dict[str, object]:
    input_xyz = _write_xyz(tmp_path / "inputs" / f"{structure_id}.xyz", input_distance)
    payload: dict[str, object] = {
        "id": structure_id,
        "role": role,
        "kind": kind,
        "input_xyz": str(input_xyz),
        "source_stage": source_stage,
    }
    if role in {"intermediate", "ts"}:
        payload["forming_bonds"] = [[0, 1]]
    if seed_distance is not None:
        payload["original_seed_xyz"] = str(
            _write_xyz(tmp_path / "inputs" / f"{structure_id}_seed.xyz", seed_distance)
        )
    payload.update(extra)
    return payload


def _read_manifest(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _install_default_qc_mocks(monkeypatch: pytest.MonkeyPatch):
    calls: list[tuple[str, object, Path, Path]] = []

    def fake_opt(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del config
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
        if subprocess_callback is not None:
            subprocess_callback(type("Process", (), {"pid": 101, "poll": lambda self: 0})())
        calls.append(("opt", spec, Path(input_xyz), output_dir))
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_xyz=out_xyz,
            output_file=output_dir / "opt.out",
            energy_hartree=-1.0,
        )

    def fake_frequency(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / "freq.out"
        out_file.write_text("mock freq\n", encoding="utf-8")
        calls.append(("freq", spec, Path(input_xyz), output_dir))
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=out_file,
            energy_hartree=-1.0,
            frequencies_cm1=(-321.4, 120.0, 311.2) if "ts" in str(output_dir) else (25.0, 125.0, 325.0),
            zero_point_energy_hartree=0.12,
            thermal_energy_hartree=-0.97,
            enthalpy_hartree=-0.98,
            gibbs_free_energy_hartree=-1.04,
            gibbs_correction_hartree=-0.04,
        )

    def fake_sp(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del config
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        if subprocess_callback is not None:
            subprocess_callback(type("Process", (), {"pid": 202, "poll": lambda self: 0})())
        calls.append(("sp", spec, Path(input_xyz), output_dir))
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=output_dir / "sp.out",
            energy_hartree=-1.1,
        )

    monkeypatch.setattr(engine_module, "run_optimization", fake_opt)
    monkeypatch.setattr(engine_module, "run_frequency", fake_frequency)
    monkeypatch.setattr(engine_module, "run_single_point", fake_sp)
    return calls


def test_stage_calculator_maps_s4_orca_controls(monkeypatch, tmp_path: Path):
    calls = _install_default_qc_mocks(monkeypatch)
    manifest = _read_manifest(
        _engine("S4", irc_enabled=False).run(
            [_structure(tmp_path, structure_id="ts", role="ts", kind="ts")],
            tmp_path / "stage",
        )
    )

    opt = cast(
        QCJobSpec,
        next(spec for kind, spec, *_rest in calls if kind == "opt" and getattr(spec, "task") == "opt_ts"),
    )
    sp = cast(QCJobSpec, next(spec for kind, spec, *_rest in calls if kind == "sp"))
    assert opt.engine == "orca"
    assert opt.task == "opt_ts"
    assert opt.route == "OptTS"
    assert opt.aux_basis == "def2/J"
    assert opt.grid == "DefGrid3"
    assert opt.scf == "TightSCF"
    assert opt.max_cycles == 200
    assert opt.timeout == 864000
    assert sp.engine == "orca"
    assert sp.method == "wB97M-V"
    assert sp.aux_basis == "def2/J"
    assert sp.timeout == 864000
    assert manifest["structures"][0]["status"] == "complete"


def test_stage_calculator_runs_and_validates_ts_frequency(monkeypatch, tmp_path: Path):
    _install_default_qc_mocks(monkeypatch)
    monkeypatch.setattr(
        engine_module.identity_module,
        "classify_ts",
        lambda **kwargs: {
            "hessian_index": 1,
            "curvature_class": "strict",
            "mode_identity": "target",
            "stationary_point_class": "valid_target_ts",
        },
    )
    manifest = _read_manifest(
        _engine("S3", irc_enabled=False).run(
            [_structure(tmp_path, structure_id="ts", role="ts", kind="ts")],
            tmp_path / "stage",
        )
    )

    result = manifest["structures"][0]
    assert result["frequency_status"] == "complete"
    assert result["imaginary_frequencies_cm1"] == [-321.4]
    assert result["ts_classification"]["stationary_point_class"] == "valid_target_ts"
    assert result["status"] == "complete"


def test_stage_calculator_marks_ts_unusable_when_frequency_fails(monkeypatch, tmp_path: Path):
    def fake_opt(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_xyz=out_xyz,
            output_file=output_dir / "opt.out",
            energy_hartree=-1.0,
        )

    monkeypatch.setattr(engine_module, "run_optimization", fake_opt)
    monkeypatch.setattr(
        engine_module,
        "run_frequency",
        lambda spec, input_xyz, output_dir, config, subprocess_callback=None: QCJobResult(
            status="failed",
            input_xyz=Path(input_xyz),
            output_file=Path(output_dir) / "freq.out",
            error="parse failure",
        ),
    )
    monkeypatch.setattr(
        engine_module,
        "run_single_point",
        lambda spec, input_xyz, output_dir, config, subprocess_callback=None: QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=Path(output_dir) / "sp.out",
            energy_hartree=-1.1,
        ),
    )
    manifest = _read_manifest(
        _engine("S3", irc_enabled=False).run(
            [_structure(tmp_path, structure_id="ts", role="ts", kind="ts")],
            tmp_path / "stage",
        )
    )

    result = manifest["structures"][0]
    assert result["frequency_status"] == "failed"
    assert result["ml_usability"]["frequency_descriptors"] is False
    assert result["ml_usability"]["ts_descriptors"] is False


def test_stage_calculator_runs_minimum_frequency_and_archives_thermochemistry(
    monkeypatch,
    tmp_path: Path,
):
    _install_default_qc_mocks(monkeypatch)
    manifest = _read_manifest(
        _engine("S3", irc_enabled=False).run(
            [
                _structure(
                    tmp_path,
                    structure_id="minimum",
                    role="product",
                    kind="minimum",
                    s1_thermochemistry_status="complete",
                )
            ],
            tmp_path / "stage",
        )
    )

    result = manifest["structures"][0]
    assert result["frequency_status"] == "complete"
    assert result["minimum_classification"]["imaginary_count"] == 0
    assert result["frequency_zero_point_energy_hartree"] == 0.12
    assert result["frequency_enthalpy_hartree"] == -0.98
    assert result["frequency_gibbs_free_energy_hartree"] == -1.04
    assert result["frequency_gibbs_correction_hartree"] == -0.04
    assert result["thermochemistry"]["enthalpy_hartree"] == pytest.approx(-1.08)
    assert result["thermochemistry"]["gibbs_free_energy_hartree"] == pytest.approx(-1.14)
    assert result["status"] == "complete"
    assert result["usable_for_ml"] is True


def test_stage_calculator_rejects_minimum_with_significant_imaginary_mode(
    monkeypatch,
    tmp_path: Path,
):
    calls = []

    def fake_opt(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
        calls.append("opt")
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_xyz=out_xyz,
            output_file=output_dir / "opt.out",
            energy_hartree=-1.0,
        )

    def fake_frequency(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        calls.append("freq")
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=output_dir / "freq.out",
            frequencies_cm1=(-120.0, 100.0),
            energy_hartree=-1.0,
            enthalpy_hartree=-0.98,
            gibbs_free_energy_hartree=-1.04,
        )

    monkeypatch.setattr(engine_module, "run_optimization", fake_opt)
    monkeypatch.setattr(engine_module, "run_frequency", fake_frequency)
    monkeypatch.setattr(
        engine_module,
        "run_single_point",
        lambda spec, input_xyz, output_dir, config, subprocess_callback=None: QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=Path(output_dir) / "sp.out",
            energy_hartree=-1.1,
        ),
    )
    manifest = _read_manifest(
        _engine("S3", irc_enabled=False).run(
            [
                _structure(
                    tmp_path,
                    structure_id="minimum",
                    role="product",
                    kind="minimum",
                    s1_thermochemistry_status="complete",
                )
            ],
            tmp_path / "stage_bad",
        )
    )

    result = manifest["structures"][0]
    assert result["minimum_classification"]["identity"] == "imaginary_frequency"
    assert result["sp_status"] == "complete"
    assert calls[:2] == ["opt", "freq"]
    assert len(calls) >= 2


def test_s3_and_s4_apply_s1_ensemble_thermochemistry_correction(monkeypatch, tmp_path: Path):
    def fake_opt(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
        return QCJobResult("complete", Path(input_xyz), out_xyz, output_dir / "opt.out", -100.0)

    def fake_frequency(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=output_dir / "freq.out",
            energy_hartree=-100.0,
            frequencies_cm1=(25.0, 125.0, 325.0),
            enthalpy_hartree=-100.0,
            gibbs_free_energy_hartree=-100.0,
            gibbs_correction_hartree=0.0,
        )

    def fake_sp(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        return QCJobResult("complete", Path(input_xyz), output_file=output_dir / "sp.out", energy_hartree=-100.0)

    monkeypatch.setattr(engine_module, "run_optimization", fake_opt)
    monkeypatch.setattr(engine_module, "run_frequency", fake_frequency)
    monkeypatch.setattr(engine_module, "run_single_point", fake_sp)
    structure = _structure(
        tmp_path,
        structure_id="product",
        role="product",
        kind="minimum",
        ensemble_thermochemistry_correction_hartree=-0.0125,
        s1_ensemble_thermodynamics={"partition_function_relative": 2.0},
    )
    s3_manifest = _read_manifest(_engine("S3", irc_enabled=False).run([structure], tmp_path / "s3"))
    s4_manifest = _read_manifest(_engine("S4", irc_enabled=False).run([structure], tmp_path / "s4"))

    assert s3_manifest["structures"][0]["thermochemistry"]["gibbs_free_energy_hartree"] == pytest.approx(-100.0125)
    assert s3_manifest["structures"][0]["s1_ensemble_thermodynamics"]["partition_function_relative"] == 2.0
    assert s4_manifest["structures"][0]["thermochemistry"]["gibbs_free_energy_hartree"] == pytest.approx(-100.0125)


def test_s3_and_s4_do_not_treat_incomplete_s1_thermochemistry_as_zero(tmp_path: Path, monkeypatch):
    _install_default_qc_mocks(monkeypatch)
    structure = _structure(
        tmp_path,
        structure_id="product",
        role="product",
        kind="minimum",
        s1_thermochemistry_status="incomplete",
        ensemble_thermochemistry_correction_hartree=None,
        s1_ensemble_thermodynamics={
            "status": "incomplete",
            "partition_function_relative": None,
        },
    )
    s3_manifest = _read_manifest(_engine("S3", irc_enabled=False).run([structure], tmp_path / "s3_incomplete"))
    s4_manifest = _read_manifest(_engine("S4", irc_enabled=False).run([structure], tmp_path / "s4_incomplete"))

    assert s3_manifest["structures"][0]["s1_thermochemistry_status"] == "incomplete"
    assert s3_manifest["structures"][0]["s1_ensemble_thermodynamics"]["status"] == "incomplete"
    assert s4_manifest["structures"][0]["s1_thermochemistry_status"] == "incomplete"
    assert s4_manifest["structures"][0]["s1_ensemble_thermodynamics"]["status"] == "incomplete"


def test_s3_uses_sequential_constrained_warmup_for_s2_seed(monkeypatch, tmp_path: Path):
    warmup_xyz = tmp_path / "stage" / "intermediate" / "warmup" / "opt.xyz"
    calls = []

    def fake_opt(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del config, subprocess_callback
        calls.append((spec, Path(input_xyz), Path(output_dir)))
        if spec.bond_constraints:
            warmup_xyz.parent.mkdir(parents=True, exist_ok=True)
            warmup_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
            return QCJobResult(
                status="partial",
                input_xyz=Path(input_xyz),
                output_xyz=warmup_xyz,
                output_file=Path(output_dir) / "warmup.out",
                energy_hartree=-10.0,
                error="maximum warm-up cycles reached",
            )
        out_xyz = Path(output_dir) / "opt.xyz"
        out_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_xyz=out_xyz,
            output_file=Path(output_dir) / "opt.out",
            energy_hartree=-10.1,
        )

    monkeypatch.setattr(engine_module, "run_optimization", fake_opt)
    monkeypatch.setattr(
        engine_module,
        "run_frequency",
        lambda spec, input_xyz, output_dir, config, subprocess_callback=None: QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=Path(output_dir) / "freq.out",
            energy_hartree=-10.1,
            frequencies_cm1=(25.0, 125.0, 325.0),
            enthalpy_hartree=-10.0,
            gibbs_free_energy_hartree=-10.0,
        ),
    )
    monkeypatch.setattr(
        engine_module,
        "run_single_point",
        lambda spec, input_xyz, output_dir, config, subprocess_callback=None: QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=Path(output_dir) / "sp.out",
            energy_hartree=-10.2,
        ),
    )
    manifest = _read_manifest(
        _engine("S3", irc_enabled=False).run(
            [
                _structure(
                    tmp_path,
                    structure_id="intermediate",
                    role="intermediate",
                    kind="minimum",
                    input_distance=2.4,
                    seed_distance=2.4,
                    s1_thermochemistry_status="complete",
                )
            ],
            tmp_path / "stage",
        )
    )

    result = manifest["structures"][0]
    assert len(calls) >= 2
    assert calls[0][0].task == "opt"
    assert calls[0][0].max_cycles == 40
    assert calls[0][0].route_extras == "LooseOpt"
    assert calls[0][0].allow_unconverged_geometry is True
    assert calls[0][0].bond_constraints == ((0, 1, 2.4),)
    assert calls[1][1] == warmup_xyz
    assert not calls[1][0].bond_constraints
    assert "LooseOpt" not in calls[1][0].route_extras
    assert result["warmup_status"] == "partial"
    assert result["warmup_used"] is True
    assert result["warmup_constraints"][0]["target_distance_angstrom"] == 2.4
    assert result["warmup_max_cycles"] == 40


def test_s3_warmup_resolves_role_specific_cycle_budgets():
    profile = FidelityProfile.from_config(_config(), "S3")

    assert profile.warmup_max_cycles_int == 40
    assert profile.warmup_max_cycles_ts == 50


def test_s3_engine_archives_prior_outputs_and_records_manifest_path(monkeypatch, tmp_path: Path):
    _install_default_qc_mocks(monkeypatch)
    structures = [
        _structure(
            tmp_path,
            structure_id="product",
            role="product",
            kind="minimum",
            s1_thermochemistry_status="complete",
        )
    ]
    engine = LowLevelEngine(_config())
    stage_dir = tmp_path / "S3_LowLevel"

    first_manifest = _read_manifest(engine.run(structures, stage_dir))
    second_manifest = _read_manifest(engine.run(structures, stage_dir))

    assert first_manifest["stale_outputs_archived_to"] is None
    archived_to = Path(second_manifest["stale_outputs_archived_to"])
    assert archived_to.is_dir()
    assert (archived_to / "manifest.json").is_file()
    assert (archived_to / "product").is_dir()
    assert stage_dir.joinpath("product").is_dir()


def test_s4_engine_manifest_records_null_archived_path_on_first_run(monkeypatch, tmp_path: Path):
    _install_default_qc_mocks(monkeypatch)
    manifest = _read_manifest(
        HighLevelEngine(_config()).run(
            [
                _structure(
                    tmp_path,
                    structure_id="product",
                    role="product",
                    kind="minimum",
                    s1_thermochemistry_status="complete",
                )
            ],
            tmp_path / "S4_HighLevel",
        )
    )

    assert manifest["schema_version"] == "refinement_manifest_v1"
    assert manifest["stale_outputs_archived_to"] is None
