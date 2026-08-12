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
    assert optimization["frequency"]["enabled_for_minima"] is True
    assert optimization["frequency"]["task"] == "freq"
    assert optimization["frequency"]["require_exactly_one"] is True
    assert optimization["frequency"]["soft_mode_window_cm1"] == [-50.0, -10.0]
    assert optimization["frequency"]["soft_mode_review"]["enabled"] is True
    assert optimization["frequency"]["soft_mode_review"]["perturbation_scale_bohr"] == 0.1
    assert optimization["frequency"]["soft_mode_review"]["projection_threshold"] == 0.3
    assert optimization["frequency"]["soft_mode_review"]["max_perturbation_attempts"] == 2
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
    assert s4_optimization["frequency"]["enabled_for_minima"] is True
    assert s4_optimization["frequency"]["task"] == "freq"
    assert s4_optimization["frequency"]["soft_mode_window_cm1"] == [-50.0, -10.0]
    assert s4_optimization["frequency"]["soft_mode_review"]["enabled"] is True
    assert config["theory"]["s4_high_precision"]["single_point"]["engine"] == "orca"
    assert config["theory"]["s4_high_precision"]["single_point"]["solvent_model"] == "CPCM"
    assert config["step3"]["scheduling"]["max_workers"] == 1
    assert config["step3"]["scheduling"]["nproc_per_job"] == 16
    assert config["step3"]["scheduling"]["memory_per_job"] == "32GB"
    assert config["step4"]["scheduling"]["max_workers"] == 1
    assert config["step4"]["scheduling"]["nproc_per_job"] == 16
    assert config["step4"]["scheduling"]["memory_per_job"] == "32GB"
    scan = config["step2"]["scan"]
    assert scan["scan_start_distance"] > scan["scan_end_distance"]
    assert config["step2"]["method"] == "orca_gfn2_relaxed_scan"
    assert config["step2"]["orca_gfn2_scan"]["relaxed_scan"]["points"] == 17
    assert config["step2"]["orca_gfn2_scan"]["relaxed_scan"]["solvent_model"] == "ALPB"
    assert config["theory"]["s3_low_level"]["optimization"]["solvent_model"] == "CPCM"
    assert "xtb_path" not in config["step2"]
    assert scan["selection"]["preferred_energy_source"] == "b973c"
    assert config["ui"]["recovery"]["enabled"] is True
    assert config["ui"]["recovery"]["stale_heartbeat_seconds"] == 300
    assert config["ui"]["recovery"]["requeue_interrupted_structures"] is True
    assert config["step2"]["energy_refinement"]["method"] == "B97-3c"
    assert config["refinement"]["common"]["workflow"]["warmup"]["enabled_roles"] == [
        "intermediate",
        "ts",
    ]
    assert config["refinement"]["s3"]["profile_id"] == "b97_3c_r2scan_3c_v1"
    assert config["refinement"]["s4"]["profile_id"] == "m062x_wb97mv_v1"
    assert config["refinement"]["s3"]["geometry"]["method"] == "B97-3c"
    assert config["refinement"]["s4"]["geometry"]["method"] == "M062X"
