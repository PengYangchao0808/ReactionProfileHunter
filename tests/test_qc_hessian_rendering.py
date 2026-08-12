from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

from rph_core.utils import qc_jobs
from rph_core.utils.qc_jobs import _orca_geom_block, run_irc
from rph_core.utils.qc_models import IRCJobSpec, OptimizationSpec, QCJobResult, QCJobSpec


def _spec(**overrides):
    return QCJobSpec(engine="orca", task="opt", method="B97-3c", **overrides)


def test_orca_geom_block_default_unchanged():
    assert _orca_geom_block(_spec(max_cycles=8)) == "%geom\n  MaxIter 8\nend"


def test_orca_geom_block_calc_hess_true():
    block = _orca_geom_block(_spec(initial_hessian="calculate"))

    assert "Calc_Hess true" in block


def test_orca_geom_block_inhess_read():
    block = _orca_geom_block(
        _spec(initial_hessian="read", hessian_filename="prelim.hess")
    )

    assert "InHess Read" in block
    assert 'InHessName "prelim.hess"' in block


def test_orca_geom_block_ts_mode():
    assert "TS_Mode {M 0}" in _orca_geom_block(_spec(ts_mode=0))


def test_orca_geom_block_recalc_hess():
    assert "Recalc_Hess 5" in _orca_geom_block(_spec(recalc_hessian=5))


def test_orca_geom_block_trust():
    assert "Trust 0.3" in _orca_geom_block(_spec(trust=0.3))


def test_orca_geom_block_combined():
    block = _orca_geom_block(
        _spec(
            max_cycles=8,
            initial_hessian="calculate",
            ts_mode=0,
            recalc_hessian=5,
            trust=0.3,
        )
    )

    assert block.splitlines() == [
        "%geom",
        "  MaxIter 8",
        "  Calc_Hess true",
        "  TS_Mode {M 0}",
        "  Recalc_Hess 5",
        "  Trust 0.3",
        "end",
    ]


def test_orca_geom_block_bond_constraints_preserved():
    assert _orca_geom_block(_spec(bond_constraints=((0, 1, 1.5),))) == (
        "%geom\n"
        "  Constraints\n"
        "    { B 0 1 1.50000000 C }\n"
        "  end\n"
        "end"
    )


def test_irc_job_spec_defaults():
    spec = IRCJobSpec()

    assert spec.direction == "both"
    assert spec.max_iter == 50
    assert spec.init_hessian == "read"
    assert spec.hessian_filename is None
    assert spec.hessian_mode == 0
    assert spec.init_displacement_mode == "energy"
    assert spec.initial_delta_energy_mEh == 2.0
    assert spec.scale_initial_displacement == 0.1
    assert spec.scale_steepest_descent == 0.15
    assert spec.adaptive_step is True
    assert spec.method == ""
    assert spec.basis == ""


def test_optimization_spec_defaults():
    spec = OptimizationSpec(method="B97-3c")

    assert spec.task == "opt"
    assert spec.method == "B97-3c"
    assert spec.basis == ""
    assert spec.aux_basis == ""
    assert spec.route == "Opt"
    assert spec.route_extras == ""
    assert spec.initial_hessian == "model"
    assert spec.hessian_filename is None
    assert spec.ts_mode is None
    assert spec.recalc_hessian is None
    assert spec.trust is None
    assert spec.max_cycles is None
    assert spec.bond_constraints == ()
    assert spec.allow_unconverged_geometry is False


def test_qc_job_result_extra_field():
    result = QCJobResult(
        status="complete",
        input_xyz=Path("/tmp/x.xyz"),
        extra={"foo": "bar"},
    )

    payload = asdict(result)

    assert payload["input_xyz"] == Path("/tmp/x.xyz")
    assert payload["extra"] == {"foo": "bar"}


def test_run_irc_is_callable(tmp_path: Path):
    expected = QCJobResult(status="complete", input_xyz=Path("/tmp/x.xyz"), output_xyz=Path("/tmp/x.xyz"))

    with patch.object(qc_jobs.ORCAInterface, "run_irc", return_value=expected) as mocked_run_irc:
        result = run_irc(IRCJobSpec(method="B97-3c"), Path("/tmp/x.xyz"), tmp_path, {})

    assert result is expected
    mocked_run_irc.assert_called_once()
