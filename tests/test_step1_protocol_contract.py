from pathlib import Path
from dataclasses import asdict

import pytest

from rph_core.steps.anchor.engine_factory import create_s1_engine
from rph_core.steps.conformer_search.protocols import resolve_protocol_spec


def _base_config() -> dict[str, object]:
    return {
        "step1": {
            "protocol": "ext",
            "conformer_search": {
                "two_stage_enabled": True,
                "ngeom_default": 3,
                "ngeom_max": 6,
            },
            "protocol_stack": {
                "lite": {
                    "two_stage_enabled": False,
                    "ngeom_default": 1,
                    "ngeom_max": 1,
                    "final_opt_sp": {
                        "enabled": True,
                        "freq": True,
                        "selection_mode": "single_lowest_energy",
                    },
                },
                "zero": {
                    "two_stage_enabled": False,
                    "ngeom_default": 1,
                    "ngeom_max": 1,
                    "final_opt_sp": {
                        "enabled": True,
                        "freq": True,
                        "selection_mode": "single_lowest_energy",
                    },
                },
            },
        },
        "theory": {
            "optimization": {"engine": "gaussian"},
            "single_point": {"engine": "orca"},
        },
        "solvent": {"name": "acetone"},
    }


@pytest.mark.parametrize("protocol", ["ext", "default", "full", "lite", "zero"])
def test_protocol_spec_resolution_contract(protocol: str) -> None:
    spec = resolve_protocol_spec(_base_config(), protocol)

    assert spec.name == protocol
    assert isinstance(spec.ngeom_default, int)
    assert isinstance(spec.ngeom_max, int)
    assert isinstance(spec.selection_mode, str)


def test_protocol_spec_rejects_unknown_protocol() -> None:
    with pytest.raises(RuntimeError, match="Unknown step1.protocol"):
        resolve_protocol_spec(_base_config(), "unknown")


def test_protocol_spec_full_is_supported() -> None:
    spec = resolve_protocol_spec(_base_config(), "full")
    assert spec.name == "full"
    assert spec.handoff_policy.mode == "optimize_all_survivors_within_window"


def test_spec_to_engine_config_mapping_for_lite_zero(tmp_path: Path) -> None:
    config = _base_config()

    lite_engine = create_s1_engine(
        protocol="lite",
        config=config,
        work_dir=tmp_path,
        molecule_name="product",
    )
    zero_engine = create_s1_engine(
        protocol="zero",
        config=config,
        work_dir=tmp_path,
        molecule_name="product2",
    )

    assert lite_engine.protocol_spec.name == "lite"
    assert lite_engine.two_stage_enabled is False
    assert lite_engine.ngeom_default == 1
    assert lite_engine.max_conformers == 1

    assert zero_engine.protocol_spec.name == "zero"
    assert zero_engine.two_stage_enabled is False
    assert zero_engine.ngeom_default == 1
    assert zero_engine.max_conformers == 1


def test_lite_zero_spec_contract_captures_freq_and_selection_mode() -> None:
    config = _base_config()

    lite_spec = resolve_protocol_spec(config, "lite")
    zero_spec = resolve_protocol_spec(config, "zero")

    assert lite_spec.final_opt_sp_enabled is True
    assert zero_spec.final_opt_sp_enabled is True
    assert lite_spec.freq_enabled is True
    assert zero_spec.freq_enabled is True
    assert lite_spec.handoff_policy.mode == "optimize_rank1"
    assert zero_spec.handoff_policy.mode == "optimize_rank1"
    assert lite_spec.handoff_policy.fallback_mode == "optimize_top2_if_gap_small"
    assert zero_spec.handoff_policy.fallback_mode == "optimize_all_within_0p5_kcal"


@pytest.mark.parametrize("protocol", ["ext", "default", "full", "lite", "zero"])
def test_protocol_spec_is_serializable_and_stable(protocol: str) -> None:
    spec = resolve_protocol_spec(_base_config(), protocol)
    payload = asdict(spec)

    assert payload["name"] == protocol
    assert "two_stage_enabled" in payload
    assert "ngeom_default" in payload
    assert "ngeom_max" in payload
    assert "funnel_policy" in payload
    assert "handoff_policy" in payload
    assert "final_opt_sp_enabled" in payload
    assert "freq_enabled" in payload
    assert "final_sp_enabled" in payload
    assert "selection_mode" in payload
    assert "provenance_flags" in payload
    assert isinstance(payload["provenance_flags"], dict)
