from pathlib import Path

from rph_core.steps import stage_calculator as stage_module
from rph_core.steps.stage_calculator import StageCalculator
from rph_core.steps.step3_lowlevel.engine import LowLevelEngine
from rph_core.steps.step4_highlevel.engine import HighLevelEngine
from rph_core.utils.qc_models import QCJobResult


def test_stage_calculator_maps_s4_orca_controls(monkeypatch, tmp_path: Path):
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
            "engine": "orca",
            "method": "M062X",
            "basis": "def2-SVP",
            "aux_basis": "def2/J",
            "route_ts": "OptTS",
            "grid": "DefGrid3",
            "scf": "TightSCF",
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
    assert opt.engine == "orca"
    assert opt.task == "opt_ts"
    assert opt.route == "OptTS"
    assert opt.aux_basis == "def2/J"
    assert opt.grid == "DefGrid3"
    assert opt.scf == "TightSCF"
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


def test_s3_and_s4_apply_s1_ensemble_thermochemistry_correction(tmp_path: Path):
    input_xyz = tmp_path / "minimum.xyz"
    input_xyz.write_text("1\nH\nH 0.0 0.0 0.0\n", encoding="utf-8")
    structure = {
        "id": "product",
        "kind": "minimum",
        "input_xyz": str(input_xyz),
        "ensemble_thermochemistry_correction_hartree": -0.0125,
        "s1_ensemble_thermodynamics": {"partition_function_relative": 2.0},
    }

    class FakeCalculator:
        def run_structure(self, _structure, _target):
            return {
                "sp_energy_hartree": -100.0,
                "sp_status": "complete",
                "status": "complete",
                "usable_for_ml": True,
            }

    s3 = LowLevelEngine({"theory": {"s3_low_level": {}}})._run_one(
        structure, tmp_path / "s3", FakeCalculator(), None
    )
    assert s3["ensemble_corrected_sp_free_energy_hartree"] == -100.0125
    assert s3["s1_ensemble_thermodynamics"]["partition_function_relative"] == 2.0

    class FakeReporter:
        def start_structure(self, _structure_id):
            pass

        def finish_structure(self, _payload):
            pass

    s4 = HighLevelEngine({"theory": {"s4_high_precision": {}}})._run_one(
        structure, tmp_path / "s4", FakeCalculator(), FakeReporter()
    )
    assert s4["ensemble_corrected_sp_free_energy_hartree"] == -100.0125


def test_s3_and_s4_do_not_treat_incomplete_s1_thermochemistry_as_zero(tmp_path: Path):
    input_xyz = tmp_path / "minimum_incomplete.xyz"
    input_xyz.write_text("1\nH\nH 0.0 0.0 0.0\n", encoding="utf-8")
    structure = {
        "id": "product",
        "kind": "minimum",
        "input_xyz": str(input_xyz),
        "s1_thermochemistry_status": "incomplete",
        "ensemble_thermochemistry_correction_hartree": None,
        "s1_ensemble_thermodynamics": {
            "status": "incomplete",
            "partition_function_relative": None,
        },
    }

    class FakeCalculator:
        def run_structure(self, _structure, _target):
            return {
                "sp_energy_hartree": -100.0,
                "sp_status": "complete",
                "status": "complete",
                "usable_for_ml": True,
            }

    s3 = LowLevelEngine({"theory": {"s3_low_level": {}}})._run_one(
        structure, tmp_path / "s3_incomplete", FakeCalculator(), None
    )
    assert s3["ensemble_corrected_sp_free_energy_hartree"] is None
    assert s3["ensemble_free_energy_status"] == "thermochemistry_incomplete"

    class FakeReporter:
        def start_structure(self, _structure_id):
            pass

        def finish_structure(self, _payload):
            pass

    s4 = HighLevelEngine({"theory": {"s4_high_precision": {}}})._run_one(
        structure, tmp_path / "s4_incomplete", FakeCalculator(), FakeReporter()
    )
    assert s4["ensemble_corrected_sp_free_energy_hartree"] is None
    assert s4["ensemble_free_energy_status"] == "thermochemistry_incomplete"


def test_s3_uses_sequential_constrained_warmup_for_s2_seed(monkeypatch, tmp_path: Path):
    input_xyz = tmp_path / "intermediate.xyz"
    input_xyz.write_text(
        "2\nS2 intermediate\nC 0.0 0.0 0.0\nC 2.4 0.0 0.0\n",
        encoding="utf-8",
    )
    warmup_xyz = tmp_path / "stage" / "warmup" / "opt.xyz"
    calls = []

    def fake_opt(spec, input_path, output_dir, config):
        calls.append((spec, Path(input_path), Path(output_dir)))
        if spec.bond_constraints:
            warmup_xyz.parent.mkdir(parents=True, exist_ok=True)
            warmup_xyz.write_text(input_xyz.read_text(encoding="utf-8"), encoding="utf-8")
            return QCJobResult(
                "partial",
                Path(input_path),
                warmup_xyz,
                Path(output_dir) / "warmup.out",
                -10.0,
                "maximum warm-up cycles reached",
            )
        return QCJobResult(
            "complete",
            Path(input_path),
            Path(input_path),
            Path(output_dir) / "opt.out",
            -10.1,
        )

    monkeypatch.setattr(stage_module, "run_optimization", fake_opt)
    monkeypatch.setattr(
        stage_module,
        "run_single_point",
        lambda spec, input_path, output_dir, config: QCJobResult(
            "complete", Path(input_path), output_file=Path(output_dir) / "sp.out", energy_hartree=-10.2
        ),
    )

    result = StageCalculator(
        {},
        {
            "optimization": {
                "engine": "orca",
                "method": "B97-3c",
                "warmup": {
                    "enabled": True,
                    "roles": ["intermediate", "ts"],
                    "convergence": "loose",
                    "max_cycles": 8,
                    "accept_partial_geometry": True,
                },
            },
            "single_point": {"engine": "orca", "method": "r2SCAN-3c"},
        },
    ).run_structure(
        {
            "id": "product_major_int",
            "kind": "minimum",
            "role": "intermediate",
            "source_stage": "S2",
            "forming_bonds": [[0, 1]],
            "input_xyz": str(input_xyz),
        },
        tmp_path / "stage",
    )

    assert len(calls) == 2
    assert calls[0][0].task == "opt"
    assert calls[0][0].max_cycles == 8
    assert calls[0][0].route_extras == "LooseOpt"
    assert calls[0][0].allow_unconverged_geometry is True
    assert calls[0][0].bond_constraints == ((0, 1, 2.4),)
    assert calls[1][1] == warmup_xyz
    assert not calls[1][0].bond_constraints
    assert "LooseOpt" not in calls[1][0].route_extras
    assert result["warmup_status"] == "partial"
    assert result["warmup_used"] is True
    assert result["optimization_input_xyz"] == str(warmup_xyz)
    assert result["warmup_constraints"][0]["target_distance_angstrom"] == 2.4
    assert result["warmup_convergence"] == "loose"
    assert result["warmup_max_cycles"] == 8


def test_s3_warmup_resolves_role_specific_cycle_budgets():
    warmup_cfg = {
        "max_cycles": {
            "intermediate": 8,
            "ts": 12,
        }
    }

    assert StageCalculator._warmup_max_cycles(
        {"role": "intermediate", "kind": "minimum"}, warmup_cfg
    ) == 8
    assert StageCalculator._warmup_max_cycles(
        {"role": "ts", "kind": "ts"}, warmup_cfg
    ) == 12
    assert StageCalculator._warmup_max_cycles(
        {"role": "intermediate"}, {"max_cycles": 6}
    ) == 6
