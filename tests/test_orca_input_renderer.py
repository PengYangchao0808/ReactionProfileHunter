from rph_core.utils.method_registry import MethodRegistry
from rph_core.utils.capability_validator import CapabilityValidator
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.utils.orca_input_renderer import OrcaInputRenderer
from rph_core.utils.optimization_config import OptimizationConfig, prepare_qc_method_config


def test_renderer_double_hybrid_includes_aux_and_dispersion():
    spec = MethodRegistry.normalize_spec(
        {
            "engine": "orca",
            "family": "double_hybrid_dft",
            "method": "PWPB95",
            "basis": "def2-TZVPP",
            "scf": {"accel": "RIJCOSX", "aux_j": "def2/J"},
            "correlation": {"model": "RI-MP2", "aux_c": "def2-TZVPP/C"},
            "dispersion": {"mode": "external", "keyword": "D4"},
        }
    )
    renderer = OrcaInputRenderer(spec)
    line = renderer.render_simple_keywords(task_type="sp")
    blocks = renderer.render_blocks()

    assert "PWPB95" in line
    assert "RIJCOSX" in line
    assert "def2/J" in line
    assert "def2-TZVPP/C" in line
    assert "D4" in line
    assert blocks == {}


def test_renderer_3c_omits_basis_and_aux_blocks():
    spec = MethodRegistry.normalize_spec(
        {
            "engine": "orca",
            "family": "composite_3c",
            "method": "r2SCAN-3c",
            "basis": "mTZVPP",
            "dispersion": {"mode": "forbidden"},
        }
    )
    renderer = OrcaInputRenderer(spec)
    line = renderer.render_simple_keywords(task_type="sp")
    blocks = renderer.render_blocks()

    assert "r2SCAN-3c" in line
    assert "mTZVPP" not in line
    assert "RIJCOSX" not in line
    assert blocks == {}


def test_capability_validator_rejects_vv10_with_external_dispersion():
    spec = MethodRegistry.normalize_spec(
        {
            "engine": "orca",
            "family": "vv10_family",
            "method": "wB97M-V",
            "basis": "def2-TZVPP",
            "dispersion": {"mode": "external", "keyword": "D4"},
        }
    )
    errors = CapabilityValidator.validate(spec, task_type="sp")
    assert errors
    assert any("VV10" in item for item in errors)


def test_composite_3c_normalize_clears_aux_basis_when_not_in_input():
    """
    Regression test: composite_3c (r2SCAN-3c) should have aux_basis=None
    when not specified in input, regardless of profile defaults.
    
    Bug: defaults.yaml aux_basis=def2/J was leaking into SP-6 (r2SCAN-3c)
    because normalize_spec didn't clear aux_basis for self_contained methods.
    """
    spec = MethodRegistry.normalize_spec(
        {
            "engine": "orca",
            "family": "composite_3c",
            "method": "r2SCAN-3c",
            "basis": "mTZVPP",
            "dispersion": {"mode": "forbidden"},
        }
    )
    # aux_basis must be cleared for self-contained composite_3c
    assert spec.aux_basis is None, f"Expected aux_basis=None for composite_3c, got {spec.aux_basis}"
    # Verify the method is properly marked as self-contained
    assert spec.scf.aux_j is None, f"Expected scf.aux_j=None for composite_3c, got {spec.scf.aux_j}"
    assert spec.correlation.aux_c is None, f"Expected correlation.aux_c=None for composite_3c, got {spec.correlation.aux_c}"


def test_composite_3c_prepare_qc_method_config_clears_inherited_aux_basis():
    """
    Regression test: when prepare_qc_method_config merges defaults.yaml
    containing aux_basis=def2/J into composite_3c spec, aux_basis must
    still be cleared after normalization.
    
    Bug: prepare_qc_method_config -> normalize_spec path was passing
    through defaults.yaml aux_basis=def2/J into r2SCAN-3c specs.
    """
    base_config = {
        "theory": {
            "single_point": {
                "method": "wB97X-D4",
                "basis": "def2-TZVPP",
                "aux_basis": "def2/J",
                "engine": "orca",
            }
        }
    }
    
    result = prepare_qc_method_config(
        base_config,
        single_point={
            "method": "r2SCAN-3c",
            "family": "composite_3c",
            "basis": "mTZVPP",
            "dispersion": {"mode": "forbidden"},
        }
    )
    
    normalized_sp = result.get("theory", {}).get("single_point", {}).get("normalized_spec")
    assert normalized_sp is not None, "normalized_spec missing from result"
    
    # aux_basis must be cleared in normalized spec (dict form)
    aux_basis_val = normalized_sp.get("aux_basis")
    assert aux_basis_val is None, (
        f"Expected aux_basis=None for composite_3c after prepare_qc, got {aux_basis_val}"
    )
    scf_dict = normalized_sp.get("scf", {})
    aux_j_val = scf_dict.get("aux_j")
    assert aux_j_val is None, (
        f"Expected scf.aux_j=None for composite_3c, got {aux_j_val}"
    )


def test_prepare_qc_method_config_preserves_benchmark_solvent_model():
    base_config = {"theory": {"optimization": {}, "single_point": {}}}

    result = prepare_qc_method_config(
        base_config,
        optimization={
            "engine": "orca",
            "family": "hybrid_gga",
            "method": "B3LYP",
            "basis": "def2-SVP",
        },
        single_point={
            "engine": "orca",
            "family": "vv10_family",
            "method": "wB97M-V",
            "basis": "def2-TZVPP",
        },
        solvent={"model": "CPCM", "solvent": "acetone"},
    )

    theory = result["theory"]
    assert theory["solvent"]["model"] == "CPCM"
    assert theory["optimization"]["solvent_model"] == "CPCM"
    assert theory["single_point"]["solvent_model"] == "CPCM"
    assert theory["optimization"]["normalized_spec"]["raw"]["solvent_model"] == "CPCM"


def test_non_composite_3c_preserves_aux_basis_from_defaults():
    """
    Verify that non-self-contained methods (double_hybrid) correctly
    preserve aux_basis from defaults.yaml when not explicitly overridden.
    """
    spec = MethodRegistry.normalize_spec(
        {
            "engine": "orca",
            "family": "double_hybrid_dft",
            "method": "PWPB95",
            "basis": "def2-TZVPP",
            # No aux_basis in input - should get profile default
        }
    )
    # Double hybrid should have aux_basis from profile default
    assert spec.scf.aux_j == "def2/J", f"Expected scf.aux_j=def2/J, got {spec.scf.aux_j}"


def test_orca_xyzfile_inputs_do_not_emit_closing_star(tmp_path, monkeypatch):
    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nh2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    interface = ORCAInterface(
        method="B3LYP",
        basis="def2-SVP",
        nprocs=1,
        solvent="none",
        method_spec={
            "engine": "orca",
            "family": "hybrid_gga",
            "method": "B3LYP",
            "basis": "def2-SVP",
            "scf": {"accel": "RIJCOSX", "aux_j": "def2/J"},
        },
    )

    captured = {}

    def fake_run_orca(inp_file, output_dir, timeout=None):
        captured["normal"] = inp_file.read_text(encoding="utf-8")
        raise RuntimeError("stop after input render")

    monkeypatch.setattr(interface, "_run_orca", fake_run_orca)
    _ = interface._run_normal_optimization(xyz, tmp_path / "normal")

    ts_dir = tmp_path / "ts"
    ts_dir.mkdir(parents=True)
    ts_inp = interface._generate_ts_input(
        xyz,
        ts_dir,
        OptimizationConfig(initial_hessian="calcfc"),
        0,
        1,
    )
    captured["ts"] = ts_inp.read_text(encoding="utf-8")

    for rendered in captured.values():
        assert "* xyzfile" in rendered
        assert "AuxJ" not in rendered
        assert "AuxC" not in rendered
        assert all(line.strip() != "*" for line in rendered.splitlines())
