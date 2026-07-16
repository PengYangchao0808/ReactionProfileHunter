from rph_core.steps.conformer_search.protocols import resolve_protocol_spec
from rph_core.utils.config_loader import load_config


def test_v4_accepts_only_censo_lite():
    spec = resolve_protocol_spec({"step1": {"protocol": "censo_lite"}}, "censo_lite")
    assert spec.name == "censo_lite"
    assert spec.funnel_policy.clustering_mode == "torsion_signature"


def test_v4_rejects_legacy_protocols():
    try:
        resolve_protocol_spec({"step1": {"protocol": "lite"}}, "lite")
    except RuntimeError as exc:
        assert "censo_lite" in str(exc)
    else:
        raise AssertionError("legacy protocol was accepted")


def test_default_config_is_v4_only():
    config = load_config()
    assert config["schema_version"] == "rph_v4_config_v1"
    assert config["step1"]["allowed_protocols"] == ["censo_lite"]
    assert "protocol_stack" not in config["step1"]
    assert "optimization" not in config["theory"]
    assert config["theory"]["s3_low_level"]["optimization"]["engine"] == "orca"
    assert config["theory"]["s4_high_precision"]["optimization"]["engine"] == "orca"
    optimization = config["theory"]["s3_low_level"]["optimization"]
    assert optimization["route_minimum"] == "Opt"
    assert optimization["route_ts"] == "OptTS"
    assert optimization["route_extras"] == ""
    assert optimization["solvent_model"] == "CPCM"
    assert optimization["frequency"]["enabled_for_ts"] is True
    assert optimization["frequency"]["task"] == "freq"
    assert optimization["frequency"]["require_exactly_one"] is True
    assert config["theory"]["s3_low_level"]["single_point"]["solvent_model"] == "CPCM"
    s4_optimization = config["theory"]["s4_high_precision"]["optimization"]
    assert s4_optimization["method"] == "M062X"
    assert s4_optimization["basis"] == "def2-SVP"
    assert s4_optimization["aux_basis"] == "def2/J"
    assert s4_optimization["route_minimum"] == "Opt"
    assert s4_optimization["route_ts"] == "OptTS"
    assert s4_optimization["grid"] == "DefGrid3"
    assert s4_optimization["scf"] == "TightSCF"
    assert s4_optimization["solvent_model"] == "CPCM"
    assert s4_optimization["frequency"]["enabled_for_ts"] is True
    assert s4_optimization["frequency"]["task"] == "freq"
    assert config["theory"]["s4_high_precision"]["single_point"]["engine"] == "orca"
    assert config["theory"]["s4_high_precision"]["single_point"]["solvent_model"] == "CPCM"
    scan = config["step2"]["scan"]
    assert scan["scan_start_distance"] > scan["scan_end_distance"]
