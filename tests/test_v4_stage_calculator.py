from pathlib import Path

from rph_core.steps import stage_calculator as stage_module
from rph_core.steps.stage_calculator import StageCalculator
from rph_core.utils.qc_models import QCJobResult


def test_stage_calculator_maps_orca_and_gaussian_controls(monkeypatch, tmp_path: Path):
    input_xyz = tmp_path / "seed.xyz"
    input_xyz.write_text("1\nH\nH 0.0 0.0 0.0\n", encoding="utf-8")
    calls = []

    def fake_opt(spec, input_path, output_dir, config):
        calls.append(("opt", spec))
        return QCJobResult("complete", input_path, input_path, output_dir / "opt.out", -1.0)

    def fake_sp(spec, input_path, output_dir, config):
        calls.append(("sp", spec))
        return QCJobResult("complete", input_path, None, output_dir / "sp.out", -1.1)

    monkeypatch.setattr(stage_module, "run_optimization", fake_opt)
    monkeypatch.setattr(stage_module, "run_single_point", fake_sp)

    config = {"resources": {"nproc": 8, "mem": "16GB"}}
    theory = {
        "optimization": {
            "engine": "gaussian",
            "method": "M062X",
            "basis": "def2-SVP",
            "route_ts": "Opt=(TS,CalcFC,NoEigenTest)",
            "grid": "UltraFine",
            "scf": "XQC",
            "max_cycles": 200,
            "timeout": 123,
        },
        "single_point": {
            "engine": "orca",
            "method": "wB97M-V",
            "basis": "def2-TZVPP",
            "aux_basis": "def2/J",
            "timeout": 456,
        },
    }

    result = StageCalculator(config, theory).run_structure(
        {"id": "ts", "kind": "ts", "input_xyz": str(input_xyz)},
        tmp_path / "stage",
    )

    opt = calls[0][1]
    sp = calls[1][1]
    assert opt.engine == "gaussian"
    assert opt.task == "opt_ts"
    assert opt.route == "Opt=(TS,CalcFC,NoEigenTest)"
    assert opt.grid == "UltraFine"
    assert opt.scf == "XQC"
    assert opt.max_cycles == 200
    assert opt.timeout == 123
    assert sp.engine == "orca"
    assert sp.method == "wB97M-V"
    assert sp.aux_basis == "def2/J"
    assert sp.timeout == 456
    assert result["status"] == "complete"


def test_stage_calculator_runs_and_validates_ts_frequency(monkeypatch, tmp_path: Path):
    input_xyz = tmp_path / "ts.xyz"
    input_xyz.write_text("1\nH\nH 0.0 0.0 0.0\n", encoding="utf-8")
    calls = []

    def fake_opt(spec, input_path, output_dir, config):
        calls.append(("opt", spec))
        return QCJobResult(
            "complete",
            input_path,
            input_path,
            output_dir / "opt.out",
            -1.0,
        )

    def fake_frequency(spec, input_path, output_dir, config):
        calls.append(("freq", spec))
        return QCJobResult(
            "complete",
            input_path,
            output_file=output_dir / "numfreq.out",
            frequencies_cm1=(-321.4, 120.0, 311.2),
        )

    def fake_sp(spec, input_path, output_dir, config):
        return QCJobResult("complete", input_path, None, output_dir / "sp.out", -1.1)

    monkeypatch.setattr(stage_module, "run_optimization", fake_opt)
    monkeypatch.setattr(stage_module, "run_frequency", fake_frequency)
    monkeypatch.setattr(stage_module, "run_single_point", fake_sp)
    theory = {
        "optimization": {
            "engine": "orca",
            "method": "B97-3c",
            "frequency": {
                "enabled_for_ts": True,
                "task": "numfreq",
                "imaginary_cutoff_cm1": -50.0,
                "require_exactly_one": True,
            },
        },
        "single_point": {"engine": "orca", "method": "r2SCAN-3c"},
    }

    result = StageCalculator({}, theory).run_structure(
        {"id": "ts", "kind": "ts", "input_xyz": str(input_xyz)},
        tmp_path / "stage",
    )

    assert calls[0][1].task == "opt_ts"
    assert calls[1][1].task == "numfreq"
    assert result["frequency_status"] == "complete"
    assert result["imaginary_frequencies_cm1"] == [-321.4]
    assert result["significant_imaginary_frequencies_cm1"] == [-321.4]
    assert result["ts_frequency_valid"] is True
    assert result["status"] == "complete"


def test_stage_calculator_marks_ts_unusable_when_frequency_fails(monkeypatch, tmp_path: Path):
    input_xyz = tmp_path / "ts.xyz"
    input_xyz.write_text("1\nH\nH 0.0 0.0 0.0\n", encoding="utf-8")

    monkeypatch.setattr(
        stage_module,
        "run_optimization",
        lambda spec, input_path, output_dir, config: QCJobResult(
            "complete", input_path, input_path, output_dir / "opt.out", -1.0
        ),
    )
    monkeypatch.setattr(
        stage_module,
        "run_frequency",
        lambda spec, input_path, output_dir, config: QCJobResult(
            "failed", input_path, output_file=output_dir / "freq.out", error="parse failure"
        ),
    )
    monkeypatch.setattr(
        stage_module,
        "run_single_point",
        lambda spec, input_path, output_dir, config: QCJobResult(
            "complete", input_path, None, output_dir / "sp.out", -1.1
        ),
    )

    result = StageCalculator(
        {},
        {
            "optimization": {
                "engine": "orca",
                "method": "B97-3c",
                "frequency": {
                    "enabled_for_ts": True,
                    "task": "freq",
                    "require_exactly_one": True,
                },
            },
            "single_point": {"engine": "orca", "method": "r2SCAN-3c"},
        },
    ).run_structure({"id": "ts", "kind": "ts", "input_xyz": str(input_xyz)}, tmp_path / "stage")

    assert result["frequency_status"] == "failed"
    assert result["ts_frequency_valid"] is False
    assert result["usable_for_ml"] is False
