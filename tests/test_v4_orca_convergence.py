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
