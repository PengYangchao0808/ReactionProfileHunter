from __future__ import annotations

from dataclasses import FrozenInstanceError
from pathlib import Path
from typing import Any

import pytest
import yaml

from rph_core.steps.fidelity_profile import FidelityProfile
from rph_core.utils.config_loader import load_config


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULTS_PATH = REPOSITORY_ROOT / "config" / "defaults.yaml"


def _load_defaults_yaml() -> dict[str, Any]:
    validated = load_config(DEFAULTS_PATH)
    raw = yaml.safe_load(DEFAULTS_PATH.read_text(encoding="utf-8"))
    assert raw["schema_version"] == validated["schema_version"]
    return raw


def _make_profile(**overrides: Any) -> FidelityProfile:
    values: dict[str, Any] = {
        "stage": "S3",
        "fidelity": "low",
        "profile_id": "b97_3c_r2scan_3c_v1",
        "geometry_method": "B97-3c",
        "geometry_basis": "",
        "geometry_aux_basis": "",
        "geometry_grid": None,
        "geometry_scf": None,
        "route_minimum": "Opt",
        "route_ts": "OptTS",
        "max_cycles_minimum": 60,
        "max_cycles_intermediate": 60,
        "max_cycles_ts": 60,
        "ts_trust_radius": 0.15,
        "initial_hessian_precursor": "model",
        "initial_hessian_product": "model",
        "initial_hessian_intermediate": "calculate",
        "initial_hessian_ts": "calculate",
        "warmup_max_cycles_int": 40,
        "warmup_max_cycles_ts": 50,
        "frequency_method": "B97-3c",
        "frequency_basis": "",
        "sp_method": "r2SCAN-3c",
        "sp_basis": "",
        "sp_aux_basis": "",
        "solvent": "acetone",
        "solvent_model": "CPCM",
        "ts_rescue_enabled": True,
        "ts_rescue_recalc_hessian_interval": 5,
        "ts_rescue_trust_radius": 0.3,
        "live_monitor_enabled": True,
        "live_monitor_poll_seconds": 20.0,
        "live_monitor_min_cycles": 12,
        "live_monitor_window_cycles": 16,
        "live_monitor_min_energy_reversals": 4,
        "live_monitor_energy_change_floor_hartree": 1.0e-5,
        "live_monitor_min_gradient_improvement_fraction": 0.15,
        "live_monitor_step_tolerance_ratio": 4.0,
        "live_monitor_min_excessive_steps": 6,
        "ts_nonconvergence_rescue_enabled": True,
        "ts_nonconvergence_l1_trust_radius": 0.10,
        "ts_nonconvergence_l1_max_cycles": 30,
        "ts_nonconvergence_mode_min_overlap": 0.35,
        "ts_nonconvergence_mode_min_overlap_margin": 0.08,
        "ts_nonconvergence_l2_trust_radius": 0.05,
        "ts_nonconvergence_l2_use_default_trust": True,
        "ts_nonconvergence_l2_recalc_hessian_interval": 5,
        "ts_nonconvergence_l2_max_cycles": 30,
        "ts_nonconvergence_reseed_from_warmup_on_missing_target_mode": True,
        "int_rescue_enabled": True,
        "int_rescue_mode_displacement": True,
        "int_nonconvergence_rescue_enabled": True,
        "int_nonconvergence_l1_trust_radius": 0.10,
        "int_nonconvergence_l1_max_cycles": 30,
        "int_nonconvergence_l2_trust_radius": 0.05,
        "int_nonconvergence_l2_use_default_trust": True,
        "int_nonconvergence_l2_recalc_hessian_interval": 1,
        "int_nonconvergence_l2_max_cycles": 30,
        "irc_enabled": True,
        "irc_max_iter": 50,
        "irc_direction": "both",
        "cores_per_worker": 16,
        "memory_gb_per_worker": 32,
        "max_workers": 1,
        "timeout_seconds": 864000,
        "temperature_k": 247.55,
        "standard_state": "1M",
        "qrrho": True,
        "ensemble_correction_source": "s1",
    }
    values.update(overrides)
    return FidelityProfile(**values)


def test_fidelity_profile_construction() -> None:
    profile = _make_profile()

    assert profile.stage == "S3"
    assert profile.fidelity == "low"
    assert profile.profile_id == "b97_3c_r2scan_3c_v1"
    assert isinstance(hash(profile), int)


def test_fidelity_profile_s3_from_legacy_config() -> None:
    config = _load_defaults_yaml()

    profile = FidelityProfile.from_config(config, "S3")

    assert profile.stage == "S3"
    assert profile.fidelity == "low"
    assert profile.profile_id == "b97_3c_r2scan_3c_v1"
    assert profile.geometry_method == "B97-3c"
    assert profile.geometry_basis == ""
    assert profile.geometry_aux_basis == ""
    assert profile.geometry_grid is None
    assert profile.geometry_scf is None
    assert profile.route_minimum == "Opt"
    assert profile.route_ts == "OptTS"
    assert profile.max_cycles_minimum == 60
    assert profile.max_cycles_intermediate == 60
    assert profile.max_cycles_ts == 60
    assert profile.ts_trust_radius == pytest.approx(0.15)
    assert profile.warmup_max_cycles_int == 40
    assert profile.warmup_max_cycles_ts == 50
    assert profile.frequency_method == "B97-3c"
    assert profile.frequency_basis == ""
    assert profile.sp_method == "r2SCAN-3c"
    assert profile.sp_basis == ""
    assert profile.sp_aux_basis == ""
    assert profile.solvent == "acetone"
    assert profile.solvent_model == "CPCM"
    assert profile.cores_per_worker == 16
    assert profile.memory_gb_per_worker == 32
    assert profile.max_workers == 1
    assert profile.timeout_seconds == 864000


def test_fidelity_profile_s4_from_legacy_config() -> None:
    config = _load_defaults_yaml()

    profile = FidelityProfile.from_config(config, "S4")

    assert profile.stage == "S4"
    assert profile.fidelity == "high"
    assert profile.profile_id == "m062x_wb97mv_v1"
    assert profile.geometry_method == "M062X"
    assert profile.geometry_basis == "def2-SVP"
    assert profile.geometry_aux_basis == "def2/J"
    assert profile.geometry_grid == "DefGrid3"
    assert profile.geometry_scf == "TightSCF"
    assert profile.route_minimum == "Opt"
    assert profile.route_ts == "OptTS"
    assert profile.max_cycles_minimum == 200
    assert profile.max_cycles_ts == 200
    assert profile.warmup_max_cycles_int == 4
    assert profile.warmup_max_cycles_ts == 6
    assert profile.frequency_method == "M062X"
    assert profile.frequency_basis == "def2-SVP"
    assert profile.sp_method == "wB97M-V"
    assert profile.sp_basis == "def2-TZVPP"
    assert profile.sp_aux_basis == "def2/J"
    assert profile.solvent == "acetone"
    assert profile.solvent_model == "CPCM"
    assert profile.cores_per_worker == 16
    assert profile.memory_gb_per_worker == 32
    assert profile.max_workers == 1
    assert profile.timeout_seconds == 864000


def test_fidelity_profile_s3_from_refinement_config() -> None:
    config = {
        "resources": {"mem": "64GB", "nproc": 32},
        "step3": {
            "scheduling": {
                "max_workers": 5,
                "nproc_per_job": 20,
                "memory_per_job": "40GB",
            }
        },
        "refinement": {
            "common": {
                "workflow": {
                    "initial_hessian": {
                        "precursor": "model",
                        "product": "model",
                        "intermediate": "calculate",
                        "ts": "calculate",
                    }
                },
                "rescue": {
                    "ts_rescue": {
                        "enabled": False,
                        "recalc_hessian_interval": 7,
                        "trust_radius": 0.2,
                    },
                    "int_rescue": {
                        "enabled": False,
                        "mode_displacement_first": False,
                    },
                },
                "irc": {"enabled": False, "max_iter": 25, "direction": "backward"},
                "thermochemistry": {
                    "temperature_K": 247.55,
                    "standard_state": "1M",
                    "qrrho": True,
                    "ensemble_correction_source": "s1",
                },
            },
            "s3": {
                "fidelity": "low",
                "profile_id": "custom_s3_profile",
                "warmup": {"max_cycles": {"intermediate": 10, "ts": 14}},
                "geometry": {
                    "method": "PBEh-3c",
                    "basis": "",
                    "aux_basis": "",
                    "route_minimum": "Opt",
                    "route_ts": "OptTS",
                    "max_cycles_minimum": 111,
                    "max_cycles_ts": 166,
                    "ts_trust_radius": 0.12,
                    "timeout": 654321,
                },
                "frequency": {"method": "PBEh-3c", "basis": ""},
                "single_point": {
                    "method": "B2PLYP",
                    "basis": "def2-TZVP",
                    "aux_basis": "def2/J",
                },
                "solvent": {"solvent": "acetone", "solvent_model": "CPCM"},
                "resources": {
                    "cores_per_worker": 24,
                    "memory_gb_per_worker": 48,
                    "max_workers": 3,
                },
            },
        },
    }

    profile = FidelityProfile.from_config(config, "S3")

    assert profile.profile_id == "custom_s3_profile"
    assert profile.geometry_method == "PBEh-3c"
    assert profile.max_cycles_minimum == 111
    assert profile.max_cycles_intermediate == 111
    assert profile.max_cycles_ts == 166
    assert profile.ts_trust_radius == pytest.approx(0.12)
    assert profile.warmup_max_cycles_int == 10
    assert profile.warmup_max_cycles_ts == 14
    assert profile.frequency_method == "PBEh-3c"
    assert profile.sp_method == "B2PLYP"
    assert profile.sp_basis == "def2-TZVP"
    assert profile.sp_aux_basis == "def2/J"
    assert profile.ts_rescue_enabled is False
    assert profile.ts_rescue_recalc_hessian_interval == 7
    assert profile.ts_rescue_trust_radius == pytest.approx(0.2)
    assert profile.int_rescue_enabled is False
    assert profile.int_rescue_mode_displacement is False
    assert profile.irc_enabled is False
    assert profile.irc_max_iter == 25
    assert profile.irc_direction == "backward"
    assert profile.cores_per_worker == 24
    assert profile.memory_gb_per_worker == 48
    assert profile.max_workers == 3
    assert profile.timeout_seconds == 654321


def test_fidelity_profile_case_insensitive_stage() -> None:
    config = _load_defaults_yaml()

    assert FidelityProfile.from_config(config, "s3") == FidelityProfile.from_config(config, "S3")


def test_fidelity_profile_is_frozen() -> None:
    profile = _make_profile()

    with pytest.raises(FrozenInstanceError):
        setattr(profile, "stage", "S4")


def test_initial_hessian_for_role_all_roles() -> None:
    profile = _make_profile()

    assert profile.initial_hessian_for_role("precursor") == "model"
    assert profile.initial_hessian_for_role("product") == "model"
    assert profile.initial_hessian_for_role("intermediate") == "calculate"
    assert profile.initial_hessian_for_role("ts") == "calculate"
    assert profile.initial_hessian_for_role("TS") == "calculate"


def test_initial_hessian_for_role_unknown_raises() -> None:
    profile = _make_profile()

    with pytest.raises(ValueError, match="Unknown refinement role"):
        profile.initial_hessian_for_role("unknown")


def test_fidelity_profile_defaults_match_plan() -> None:
    config = _load_defaults_yaml()

    s3_profile = FidelityProfile.from_config(config, "S3")
    s4_profile = FidelityProfile.from_config(config, "S4")

    assert s3_profile.warmup_max_cycles_int == 40
    assert s4_profile.warmup_max_cycles_int == 4
    assert s3_profile.warmup_max_cycles_ts == 50
    assert s4_profile.warmup_max_cycles_ts == 6
    assert s3_profile.temperature_k == pytest.approx(247.55)
    assert s4_profile.temperature_k == pytest.approx(247.55)
    assert s3_profile.standard_state == "1M"
    assert s4_profile.standard_state == "1M"
    assert s3_profile.qrrho is True
    assert s4_profile.qrrho is True
    assert s3_profile.ensemble_correction_source == "s1"
    assert s4_profile.ensemble_correction_source == "s1"
    assert s3_profile.initial_hessian_for_role("precursor") == "model"
    assert s3_profile.initial_hessian_for_role("product") == "model"
    assert s3_profile.initial_hessian_for_role("intermediate") == "calculate"
    assert s3_profile.initial_hessian_for_role("ts") == "calculate"
    assert s3_profile.ts_trust_radius == pytest.approx(0.15)
    assert s4_profile.ts_trust_radius is None
    assert s4_profile.ts_rescue_recalc_hessian_interval == 5
    assert s4_profile.ts_rescue_trust_radius == pytest.approx(0.3)
    assert s4_profile.irc_max_iter == 50
    assert s4_profile.irc_direction == "both"
