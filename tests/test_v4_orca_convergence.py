import logging
from pathlib import Path

from rph_core.utils.qc_jobs import (
    _frequency_route_extras,
    _orca_geom_block,
    _orca_route_extras,
)
from rph_core.utils.method_registry import MethodRegistry
from rph_core.utils.orca_input_renderer import OrcaInputRenderer
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.utils.orca_failure_classifier import OrcaFailureType, classify_orca_failure
from rph_core.utils.qc_models import QCJobSpec


def _parser() -> ORCAInterface:
    interface = object.__new__(ORCAInterface)
    interface.logger = logging.getLogger("test.orca.parser")
    return interface


def test_orca_optimization_requires_explicit_convergence_marker(tmp_path: Path):
    output = tmp_path / "not_converged.out"
    output.write_text(
        "FINAL SINGLE POINT ENERGY      -100.123456\n"
        "ORCA TERMINATED NORMALLY\n",
        encoding="utf-8",
    )

    result = _parser()._parse_output(
        output,
        require_optimization_convergence=True,
    )

    assert result.converged is False
    assert result.energy == -100.123456
    assert "did not converge" in str(result.error_message)


def test_orca_maxiter_checkpoint_is_classified_and_persisted(tmp_path: Path):
    output_text = (
        "The optimization did not converge but reached the maximum number of\n"
        "optimization cycles.\n"
        "FINAL SINGLE POINT ENERGY      -100.123456\n"
        "ORCA TERMINATED NORMALLY\n"
    )
    classification = classify_orca_failure(
        output_text=output_text,
        returncode=0,
        max_cycles_opt=30,
    )
    assert classification is not None
    assert classification.failure_type is OrcaFailureType.GEOMETRY_OPTIMIZATION_NOT_CONVERGED

    output = tmp_path / "maxiter.out"
    output.write_text(output_text, encoding="utf-8")
    result = _parser()._parse_output(
        output,
        require_optimization_convergence=True,
        returncode=0,
        max_cycles_opt=30,
    )

    assert result.converged is False
    assert result.failure_type == "geometry_optimization_not_converged"
    assert result.optimization_converged is False
    assert result.stop_reason == "native_maxiter_checkpoint"


def test_orca_single_point_does_not_require_optimization_marker(tmp_path: Path):
    output = tmp_path / "sp.out"
    output.write_text(
        "FINAL SINGLE POINT ENERGY      -100.123456\n"
        "ORCA TERMINATED NORMALLY\n",
        encoding="utf-8",
    )

    result = _parser()._parse_output(output)

    assert result.converged is True
    assert result.energy == -100.123456


def test_orca_optimization_accepts_explicit_convergence_marker(tmp_path: Path):
    output = tmp_path / "converged.out"
    output.write_text(
        "THE OPTIMIZATION HAS CONVERGED\n"
        "FINAL SINGLE POINT ENERGY      -100.123456\n"
        "ORCA TERMINATED NORMALLY\n",
        encoding="utf-8",
    )

    result = _parser()._parse_output(
        output,
        require_optimization_convergence=True,
    )

    assert result.converged is True
    assert result.error_message is None


def test_orca_geom_block_renders_target_bond_constraints():
    block = _orca_geom_block(
        QCJobSpec(
            engine="orca",
            task="opt",
            method="B97-3c",
            max_cycles=8,
            bond_constraints=((11, 18, 2.3456789), (14, 15, 2.1)),
        )
    )

    assert "MaxIter 8" in block
    assert "{ B 11 18 2.34567890 C }" in block
    assert "{ B 14 15 2.10000000 C }" in block


def test_orca_route_maps_grid_and_deduplicates_renderer_owned_controls():
    spec = QCJobSpec(
        engine="orca",
        task="opt_ts",
        method="M062X",
        basis="def2-SVP",
        route="OptTS",
        route_extras="DefGrid3 TightSCF DefGrid3",
        grid="DefGrid3",
        scf="TightSCF",
    )

    assert _orca_route_extras(spec) == "DefGrid3"


def test_orca_frequency_inherits_grid_and_adds_analytic_frequency_keyword():
    spec = QCJobSpec(
        engine="orca",
        task="freq",
        method="M062X",
        basis="def2-SVP",
        grid="DefGrid3",
        scf="TightSCF",
    )

    assert _frequency_route_extras(spec) == "DefGrid3 Freq"


def test_s4_orca_ts_route_renders_expected_method_and_controls():
    job = QCJobSpec(
        engine="orca",
        task="opt_ts",
        method="M062X",
        basis="def2-SVP",
        aux_basis="def2/J",
        route="OptTS",
        grid="DefGrid3",
        scf="TightSCF",
    )
    normalized = MethodRegistry.normalize_spec(
        {
            "engine": job.engine,
            "method": job.method,
            "basis": job.basis,
            "aux_basis": job.aux_basis,
            "route_extras": _orca_route_extras(job),
        },
        default_engine="orca",
        default_basis=job.basis,
    )

    route = OrcaInputRenderer(normalized).render_simple_keywords(task_type="ts")

    assert "M062X" in route
    assert "def2-SVP" in route
    assert "def2/J" in route
    assert "RIJCOSX" in route
    assert "OptTS" in route
    assert "tightscf" in route.lower()
    assert "DefGrid3" in route


def test_normal_termination_with_benign_gradient_warning_is_success():
    """A converged SP whose SCF was slow to converge must NOT be flagged.

    ORCA prints "WARNING: the maximum gradient error descreased on average
    only by a factor 0.9" on every slow SCF; the old loose
    (?:grad|gradient).{0,40}(?:abort|fail|error) regex misclassified it as
    scf_gradient_aborted even though the run terminated normally with rc=0.
    """
    output_text = (
        "SCF CONVERGED AFTER  46 CYCLES\n"
        "WARNING: the maximum gradient error descreased on average only by a "
        "factor   0.9\n"
        "FINAL SINGLE POINT ENERGY      -938.996043165520\n"
        "****ORCA TERMINATED NORMALLY****\n"
    )
    assert classify_orca_failure(output_text=output_text, returncode=0) is None


def test_normal_termination_after_soscf_recovery_is_success():
    """ORCA's non-fatal '* ABORTING THE RUN *' SOSCF banner must NOT flag.

    ORCA prints the banner when SOSCF struggles, then recovers via the
    TRAH-SCF procedure and terminates normally with a final energy.
    """
    output_text = (
        "SERIOUS PROBLEM IN SOSCF\n"
        "LARGE STEP WAS ABOUT TO BE TAKEN\n"
        "ABORTING THE RUN\n"
        "Rediagonalizing the Fockian in SOSCF/NRSCF\n"
        "SCF CONVERGED AFTER  62 CYCLES\n"
        "FINAL SINGLE POINT ENERGY      -938.994112345678\n"
        "****ORCA TERMINATED NORMALLY****\n"
    )
    assert classify_orca_failure(output_text=output_text, returncode=0) is None


def test_genuine_scf_subprocess_abort_is_classified():
    output_text = (
        "SCF ITERATIONS\n"
        "ORCA finished by error termination in SCF\n"
        "[file orca_tools/qcmsg.cpp, line 165]: .... aborting the run\n"
    )
    classification = classify_orca_failure(output_text=output_text, returncode=1)
    assert classification is not None
    assert classification.failure_type is OrcaFailureType.SCF_SUBPROCESS_ABORTED


def test_genuine_scf_gradient_abort_is_classified_as_gradient():
    """'error termination in SCF gradient' must map to SCF_GRADIENT_ABORTED,
    not SCF_SUBPROCESS_ABORTED."""
    output_text = (
        "SCF ITERATIONS\n"
        "ORCA finished by error termination in SCF gradient\n"
        "[file orca_tools/qcmsg.cpp, line 465]: .... aborting the run\n"
    )
    classification = classify_orca_failure(output_text=output_text, returncode=1)
    assert classification is not None
    assert classification.failure_type is OrcaFailureType.SCF_GRADIENT_ABORTED


def test_scf_not_converged_survives_normal_termination_guard():
    """An SP whose SCF did not converge may still print TERMINATED NORMALLY
    and exit 0; the guard must NOT swallow this genuine failure."""
    output_text = (
        "SCF NOT CONVERGED AFTER 125 CYCLES\n"
        "FINAL SINGLE POINT ENERGY      -100.0 (Wavefunction not fully "
        "converged!)\n"
        "****ORCA TERMINATED NORMALLY****\n"
    )
    classification = classify_orca_failure(output_text=output_text, returncode=0)
    assert classification is not None
    assert classification.failure_type is OrcaFailureType.SCF_NOT_CONVERGED


def test_gradient_pattern_rejects_degradation_false_positive():
    """'degrad'/'gradual' words must not trigger the gradient-abort regex."""
    output_text = (
        "SCF CONVERGED AFTER  30 CYCLES\n"
        "WARNING: degradation of convergence rate detected\n"
        "FINAL SINGLE POINT ENERGY      -100.0\n"
        "****ORCA TERMINATED NORMALLY****\n"
    )
    assert classify_orca_failure(output_text=output_text, returncode=0) is None


def test_xtb_binary_missing_next_to_orca_is_classified():
    """ORCA GFN-xTB aborts when it cannot find xtb in the ORCA directory."""
    output_text = (
        "Your calculation utilizes the semiempirical GFN2-xTB method\n"
        "WARNING: otool_xtb not available.\n"
        "         Trying xtb instead   ... failed!\n"
        "  ===> : Please provide the xtb executables in the same path as where "
        "the orca binaries are located.\n"
        "         Aborting the run.\n"
        ": Error (ORCA_MAIN): ... aborting the run\n"
    )
    classification = classify_orca_failure(output_text=output_text, returncode=52)
    assert classification is not None
    assert classification.failure_type is OrcaFailureType.XTB_BINARY_UNAVAILABLE
    assert classification.retry_hint == "provide_xtb_next_to_orca"
