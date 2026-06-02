"""
Tests for S3 TS rescue policy: failure classification, rescue start selection, and route configuration.

Tests cover the four failure modes defined in ts_rescue_policy:
- oscillation_or_step_exceeded
- converged_to_minimum_no_imag
- wrong_saddle_multi_imag
- gaussian_fatal_or_no_geometry
"""
from pathlib import Path
from unittest import mock
import builtins

import numpy as np
import pytest

from rph_core.utils.data_types import QCResult
from rph_core.utils.qc_task_runner import (
    QCOptimizationResult,
    TS_FAILURE_FATAL,
    TS_FAILURE_MULTI_IMAG,
    TS_FAILURE_NO_IMAG,
    TS_FAILURE_OSCILLATION,
    TS_FAILURE_UNKNOWN,
)


def _make_runner(config):
    from rph_core.utils.qc_task_runner import QCTaskRunner
    with mock.patch.object(QCTaskRunner, "__init__", lambda self, *a, **kw: None):
        runner = QCTaskRunner.__new__(QCTaskRunner)
        runner.config = config
    return runner


@pytest.fixture
def tmp_dir(tmp_path):
    return tmp_path


@pytest.fixture
def base_config():
    return {
        "theory": {
            "optimization": {
                "method": "M062X",
                "basis": "def2-SVP",
                "dispersion": None,
                "engine": "gaussian",
                "solvent": "acetone",
            },
            "single_point": {
                "method": "wB97M-V",
                "basis": "def2-TZVPP",
                "engine": "orca",
                "solvent": "acetone",
                "maxcore": 2600,
            },
        },
        "resources": {"mem": "48GB", "nproc": 16},
        "step3": {
            "gaussian_keywords": {
                "berny": "Opt=(TS, CalcFC, NoEigenTest) Freq",
                "ts_rescue": "Opt=(TS, NoEigenTest, MaxStep=15, RecalcFC=5, MaxCycles=120) Freq",
                "qst2": "Opt=(QST2, CalcFC) Freq",
                "irc": "IRC=(CalcFC, MaxPoints=50, StepSize=10)",
            },
            "disable_qst2_rescue": True,
            "ts_rescue_policy": {
                "enabled": True,
                "trigger_after_steps": 60,
                "max_step_rescue": 15,
                "recalc_fc_every": 5,
                "max_cycles_rescue": 120,
                "oscillation_uses_last_geometry": True,
                "step_exceeded_uses_last_geometry": True,
                "no_imag_uses_original_guess": True,
                "multi_imag_uses_failed_geometry": True,
            },
        },
    }


def _make_ts_guess_xyz(path: Path) -> Path:
    xyz_file = path / "ts_guess.xyz"
    xyz_file.write_text(
        "3\nTS guess\n"
        "C     0.000000     0.000000     0.000000\n"
        "H     0.000000     0.000000     1.089000\n"
        "H     0.000000     1.089000     0.000000\n"
    )
    return xyz_file


def _make_gaussian_log_with_step_exceeded(path: Path, nstep: int = 162) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    log_path = path / "ts_guess.log"
    lines = []
    for i in range(min(nstep, 10)):
        lines.append(f" SCF Done:  E(RB3LYP) =  -{380 + i * 0.001:.6f}     A.U. after    2 cycles")
    lines.append(" Opt Step 162")
    lines.append(f" NStep=  {nstep}")
    lines.append(" Optimization stopped.")
    lines.append(f" -- Number of steps exceeded,  NStep=  {nstep}")
    lines.append(" Error termination via Lnk1e in /g16/l502.exe at Thu Jan  1 00:00:00 2026")
    log_path.write_text("\n".join(lines))
    return log_path


def _make_gaussian_log_with_orientation(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    log_path = path / "ts_guess.log"
    log_path.write_text(
        " Standard orientation:\n"
        " ---------------------------------------------------------------------\n"
        " Center     Atomic      Atomic             Coordinates (Angstroms)\n"
        " Number     Number       Type             X           Y           Z\n"
        " ---------------------------------------------------------------------\n"
        "      1          6           0        2.000000    0.000000    0.000000\n"
        "      2          1           0        2.000000    0.000000    1.089000\n"
        "      3          1           0        2.000000    1.089000    0.000000\n"
        " ---------------------------------------------------------------------\n"
        " Optimization stopped.\n"
        " -- Number of steps exceeded,  NStep=  162\n"
    )
    return log_path


class TestFailureClassification:
    def test_oscillation_step_exceeded(self, base_config, tmp_dir):
        runner = _make_runner(base_config)
        log_path = _make_gaussian_log_with_step_exceeded(tmp_dir / "ts_opt" / "berny", nstep=162)

        result = QCOptimizationResult(
            optimized_xyz=_make_ts_guess_xyz(tmp_dir),
            converged=False,
            imaginary_count=0,
            method_used="Berny_StepExceeded",
            error_message="Number of steps exceeded",
            log_file=log_path,
        )

        assert runner._classify_ts_failure(result, tmp_dir / "ts_guess.xyz") == TS_FAILURE_OSCILLATION

    def test_converged_no_imaginary(self, base_config, tmp_dir):
        runner = _make_runner(base_config)

        result = QCOptimizationResult(
            optimized_xyz=_make_ts_guess_xyz(tmp_dir),
            converged=False,
            imaginary_count=0,
            method_used="Berny_ConvergedToMinimum",
            error_message="TS optimization converged to minimum: 0 imaginary frequencies",
            log_file=None,
        )

        assert runner._classify_ts_failure(result, tmp_dir / "ts_guess.xyz") == TS_FAILURE_NO_IMAG

    def test_multi_imaginary(self, base_config, tmp_dir):
        runner = _make_runner(base_config)

        result = QCOptimizationResult(
            optimized_xyz=_make_ts_guess_xyz(tmp_dir),
            converged=True,
            imaginary_count=2,
            method_used="Berny_WrongSaddle",
            error_message="TS optimization found high-order saddle: 2 imaginary frequencies",
            log_file=None,
        )

        assert runner._classify_ts_failure(result, tmp_dir / "ts_guess.xyz") == TS_FAILURE_MULTI_IMAG

    def test_fatal_error(self, base_config, tmp_dir):
        runner = _make_runner(base_config)

        result = QCOptimizationResult(
            optimized_xyz=_make_ts_guess_xyz(tmp_dir),
            converged=False,
            imaginary_count=0,
            method_used="Berny_Fatal",
            error_message="FATAL: QPErr (gaussian log)",
            log_file=None,
        )

        assert runner._classify_ts_failure(result, tmp_dir / "ts_guess.xyz") == TS_FAILURE_FATAL

    def test_unknown_failure(self, base_config, tmp_dir):
        runner = _make_runner(base_config)

        result = QCOptimizationResult(
            optimized_xyz=_make_ts_guess_xyz(tmp_dir),
            converged=False,
            imaginary_count=0,
            method_used="Berny_Failed",
            error_message="Gaussian did not terminate normally",
            log_file=None,
        )

        assert runner._classify_ts_failure(result, tmp_dir / "ts_guess.xyz") == TS_FAILURE_UNKNOWN


class TestRescueStartSelection:
    def test_oscillation_uses_last_geometry(self, base_config, tmp_dir):
        runner = _make_runner(base_config)
        original = _make_ts_guess_xyz(tmp_dir)
        failed_xyz = tmp_dir / "ts_opt" / "berny" / "ts_guess_ts_failed.xyz"
        failed_xyz.parent.mkdir(parents=True, exist_ok=True)
        failed_xyz.write_text(
            "3\nfailed\n"
            "C     1.000000     0.000000     0.000000\n"
            "H     1.000000     0.000000     1.089000\n"
            "H     1.000000     1.089000     0.000000\n"
        )

        result = QCOptimizationResult(
            optimized_xyz=failed_xyz,
            converged=False,
            imaginary_count=0,
            method_used="Berny_StepExceeded",
            log_file=None,
        )

        rescue_xyz = runner._prepare_ts_rescue_start(
            result=result,
            original_xyz=original,
            output_dir=tmp_dir / "ts_opt",
            failure_kind=TS_FAILURE_OSCILLATION,
        )
        assert rescue_xyz.name == "rescue_start_from_oscillation.xyz"
        assert rescue_xyz.exists()

    def test_oscillation_original_geometry_falls_back_to_log(self, base_config, tmp_dir):
        runner = _make_runner(base_config)
        original = _make_ts_guess_xyz(tmp_dir)
        log_path = _make_gaussian_log_with_orientation(tmp_dir / "ts_opt" / "berny")

        result = QCOptimizationResult(
            optimized_xyz=original,
            converged=False,
            imaginary_count=0,
            method_used="Berny_StepExceeded",
            log_file=log_path,
        )

        rescue_xyz = runner._prepare_ts_rescue_start(
            result=result,
            original_xyz=original,
            output_dir=tmp_dir / "ts_opt",
            failure_kind=TS_FAILURE_OSCILLATION,
        )
        assert rescue_xyz.name == "rescue_start_from_oscillation.xyz"
        assert rescue_xyz.exists()
        assert "2.000000" in rescue_xyz.read_text()

    def test_no_imag_uses_original_guess(self, base_config, tmp_dir):
        runner = _make_runner(base_config)
        original = _make_ts_guess_xyz(tmp_dir)

        result = QCOptimizationResult(
            optimized_xyz=tmp_dir / "failed.xyz",
            converged=False,
            imaginary_count=0,
            method_used="Berny_ConvergedToMinimum",
        )

        rescue_xyz = runner._prepare_ts_rescue_start(
            result=result,
            original_xyz=original,
            output_dir=tmp_dir / "ts_opt",
            failure_kind=TS_FAILURE_NO_IMAG,
        )
        assert "no_imag_original_guess" in rescue_xyz.name

    def test_multi_imag_uses_failed_geometry(self, base_config, tmp_dir):
        runner = _make_runner(base_config)
        original = _make_ts_guess_xyz(tmp_dir)
        failed_xyz = tmp_dir / "wrong_saddle.xyz"
        failed_xyz.write_text("3\nsaddle\nC 1 0 0\nH 1 0 1\nH 1 1 0\n")

        result = QCOptimizationResult(
            optimized_xyz=failed_xyz,
            converged=True,
            imaginary_count=2,
            method_used="Berny_WrongSaddle",
        )

        rescue_xyz = runner._prepare_ts_rescue_start(
            result=result,
            original_xyz=original,
            output_dir=tmp_dir / "ts_opt",
            failure_kind=TS_FAILURE_MULTI_IMAG,
        )
        assert "wrong_saddle" in rescue_xyz.name
        assert rescue_xyz.exists()


class TestRescueRoute:
    def test_rescue_route_contains_expected_keywords(self, base_config, tmp_dir):
        runner = _make_runner(base_config)
        runner.theory_opt = base_config["theory"]["optimization"]

        route = runner._get_ts_route(rescue=True)
        assert "MaxStep=15" in route
        assert "RecalcFC=5" in route
        assert "MaxCycles=120" in route
        assert ",CalcFC" not in route
        assert "TS" in route.upper()


class TestRunTsRescueControlFlow:
    def test_fatal_failure_does_not_trigger_rescue(self, base_config, tmp_dir):
        runner = _make_runner(base_config)
        original = _make_ts_guess_xyz(tmp_dir)
        fatal_result = QCOptimizationResult(
            optimized_xyz=original,
            converged=False,
            imaginary_count=0,
            method_used="Berny_Fatal",
            error_message="FATAL: QPErr (gaussian log)",
        )
        runner._try_ts_optimization = mock.Mock(return_value=fatal_result)
        runner._try_ts_rescue = mock.Mock()

        result = runner.run_ts_opt_cycle(original, tmp_dir / "S3_TS", enable_l2_sp=False)

        assert result is fatal_result
        assert result.failure_kind == TS_FAILURE_FATAL
        runner._try_ts_rescue.assert_not_called()

    def test_unknown_failure_without_new_geometry_does_not_repeat_original_guess(self, base_config, tmp_dir):
        runner = _make_runner(base_config)
        original = _make_ts_guess_xyz(tmp_dir)
        unknown_result = QCOptimizationResult(
            optimized_xyz=original,
            converged=False,
            imaginary_count=0,
            method_used="Berny_Failed",
            error_message="Gaussian did not terminate normally",
        )
        runner._try_ts_optimization = mock.Mock(return_value=unknown_result)
        runner._try_ts_rescue = mock.Mock()

        result = runner.run_ts_opt_cycle(original, tmp_dir / "S3_TS", enable_l2_sp=False)

        assert result is unknown_result
        assert result.failure_kind == TS_FAILURE_UNKNOWN
        assert result.rescue_start == original
        runner._try_ts_rescue.assert_not_called()

    def test_ts_rescue_metadata_write_failure_does_not_drop_qc_result(self, base_config, tmp_dir, monkeypatch):
        runner = _make_runner(base_config)
        runner.ts_rescue_route = base_config["step3"]["gaussian_keywords"]["ts_rescue"]
        runner.qc_engine = mock.Mock()
        runner.qc_engine.optimize.return_value = QCResult(
            success=True,
            energy=-100.0,
            coordinates=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
            converged=True,
            frequencies=np.array([-100.0, 50.0, 100.0]),
            output_file=tmp_dir / "ts.log",
            log_file=tmp_dir / "ts.log",
        )
        original = _make_ts_guess_xyz(tmp_dir)
        real_open = builtins.open

        def _failing_metadata_open(file, *args, **kwargs):
            if Path(file).name == "ts_rescue_metadata.json":
                raise OSError("metadata disk error")
            return real_open(file, *args, **kwargs)

        monkeypatch.setattr(builtins, "open", _failing_metadata_open)

        result = runner._try_ts_rescue(
            original,
            tmp_dir / "S3_TS" / "ts_opt",
            charge=0,
            spin=1,
            old_checkpoint=None,
            failure_kind=TS_FAILURE_NO_IMAG,
        )

        assert result.converged is True
        assert result.failure_kind == TS_FAILURE_NO_IMAG
        assert result.rescue_route is not None
        assert "RecalcFC=5" in result.rescue_route


class TestRescueOnlyEntryPoint:
    def test_runs_rescue_and_l2_sp_on_success(self, base_config, tmp_dir, monkeypatch):
        runner = _make_runner(base_config)
        runner.ts_rescue_route = base_config["step3"]["gaussian_keywords"]["ts_rescue"]
        runner.qc_engine = mock.Mock()
        runner.qc_engine.optimize.return_value = QCResult(
            success=True, energy=-100.0,
            coordinates=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
            converged=True,
            frequencies=np.array([-100.0, 50.0, 100.0]),
            output_file=tmp_dir / "rescue.log",
            log_file=tmp_dir / "rescue.log",
        )

        class _FakeSP:
            def single_point(self, xyz_file, output_dir, charge=None, spin=None):
                return QCResult(success=True, energy=-200.0, converged=True,
                               output_file=output_dir / "sp.out")

        monkeypatch.setattr(
            "rph_core.utils.qc_task_runner.QCInterfaceFactory.create_interface",
            lambda engine_type, **kw: _FakeSP(),
        )
        runner.sp_engine = _FakeSP()

        rescue_xyz = _make_ts_guess_xyz(tmp_dir)
        result = runner.run_ts_rescue_only(
            rescue_start_xyz=rescue_xyz,
            output_dir=tmp_dir / "S3_TS",
            failure_kind=TS_FAILURE_NO_IMAG,
            enable_l2_sp=True,
        )

        assert result.converged is True
        assert result.method_used == "TS_Rescue"
        assert result.l2_energy == -200.0

    def test_skips_l2_sp_when_disabled(self, base_config, tmp_dir):
        runner = _make_runner(base_config)
        runner.ts_rescue_route = base_config["step3"]["gaussian_keywords"]["ts_rescue"]
        runner.qc_engine = mock.Mock()
        runner.qc_engine.optimize.return_value = QCResult(
            success=True, energy=-100.0,
            coordinates=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
            converged=True,
            frequencies=np.array([-100.0, 50.0, 100.0]),
            output_file=tmp_dir / "rescue.log",
            log_file=tmp_dir / "rescue.log",
        )

        rescue_xyz = _make_ts_guess_xyz(tmp_dir)
        result = runner.run_ts_rescue_only(
            rescue_start_xyz=rescue_xyz,
            output_dir=tmp_dir / "S3_TS",
            failure_kind=TS_FAILURE_NO_IMAG,
            enable_l2_sp=False,
        )

        assert result.converged is True
        assert result.l2_energy is None
        assert result.l2_sp_result is None


class TestStepCountParser:
    def test_nstep_from_log(self, base_config, tmp_dir):
        runner = _make_runner(base_config)
        log_path = tmp_dir / "test.log"
        log_path.write_text("Some log content\nNStep=  162\nMore content")

        assert runner._parse_gaussian_opt_step_count(log_path) == 162

    def test_fallback_grad_count(self, base_config, tmp_dir):
        runner = _make_runner(base_config)
        log_path = tmp_dir / "test.log"
        log_path.write_text("GradGradGradGrad\nGradGradGradGrad\nGradGradGradGrad")

        assert runner._parse_gaussian_opt_step_count(log_path) == 3

    def test_missing_log_returns_zero(self, base_config, tmp_dir):
        runner = _make_runner(base_config)

        assert runner._parse_gaussian_opt_step_count(None) == 0
        assert runner._parse_gaussian_opt_step_count(tmp_dir / "nonexistent.log") == 0


class TestOscillationDetector:
    def test_posthoc_oscillation_detection_insufficient_steps(self, tmp_dir):
        from rph_core.utils.oscillation_detector import analyze_gaussian_ts_log_for_oscillation

        log_path = tmp_dir / "short.log"
        lines = []
        for i in range(10):
            lines.append(f" SCF Done:  E(RB3LYP) =  -{380 + i * 0.001:.6f}     A.U.")
        log_path.write_text("\n".join(lines))

        result = analyze_gaussian_ts_log_for_oscillation(log_path, min_steps=60)
        assert not result.is_oscillating

    def test_posthoc_missing_log(self, tmp_dir):
        from rph_core.utils.oscillation_detector import analyze_gaussian_ts_log_for_oscillation

        result = analyze_gaussian_ts_log_for_oscillation(tmp_dir / "nonexistent.log")
        assert not result.is_oscillating

    def test_posthoc_oscillation_detected(self, tmp_dir):
        from rph_core.utils.oscillation_detector import analyze_gaussian_ts_log_for_oscillation

        log_path = tmp_dir / "oscillating.log"
        lines = []
        energy = -380.0
        for i in range(100):
            if i > 50:
                energy += 0.001 * (1 if i % 2 == 0 else -1)
            else:
                energy -= 0.0001
            lines.append(f" SCF Done:  E(RB3LYP) =  {energy:.6f}     A.U.")
        log_path.write_text("\n".join(lines))

        result = analyze_gaussian_ts_log_for_oscillation(log_path, min_steps=60)
        assert result.is_oscillating


class TestQCOptimizationResultFields:
    def test_new_fields_default_none(self):
        result = QCOptimizationResult(
            optimized_xyz=Path("/tmp/test.xyz"),
        )
        assert result.failure_kind is None
        assert result.rescue_start is None
        assert result.rescue_route is None

    def test_new_fields_populated(self):
        result = QCOptimizationResult(
            optimized_xyz=Path("/tmp/test.xyz"),
            failure_kind=TS_FAILURE_OSCILLATION,
            rescue_start=Path("/tmp/rescue.xyz"),
            rescue_route="Opt=(TS, NoEigenTest, MaxStep=15, RecalcFC=5) Freq",
        )
        assert result.failure_kind == "oscillation_or_step_exceeded"
        assert result.rescue_start == Path("/tmp/rescue.xyz")
        assert result.rescue_route is not None
        assert "MaxStep=15" in result.rescue_route
