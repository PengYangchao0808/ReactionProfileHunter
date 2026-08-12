"""Immutable S3/S4 refinement-profile parameters.

This module stores configuration-only values. It intentionally carries no
filesystem paths; callers should keep using :class:`pathlib.Path` at the edges
where path handling is required.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import ceil
from typing import Any, Optional

_RESOURCE_DEFAULTS: dict[str, int] = {
    "cores_per_worker": 16,
    "memory_gb_per_worker": 32,
    "max_workers": 1,
    "timeout_seconds": 864000,
}

_INITIAL_HESSIAN_DEFAULTS: dict[str, str] = {
    "precursor": "model",
    "product": "model",
    "intermediate": "calculate",
    "ts": "calculate",
}

_TS_RESCUE_DEFAULTS: dict[str, Any] = {
    "enabled": True,
    "recalc_hessian_interval": 5,
    "trust_radius": 0.3,
}

_TS_NONCONVERGENCE_RESCUE_DEFAULTS: dict[str, Any] = {
    "enabled": True,
    # R1 is a deliberately short, mode-directed OptTS block. It starts from
    # a fresh local Hessian and never burns a long generic optimization tail.
    "level_1_trust_radius": 0.03,
    "level_1_max_cycles": 12,
    "mode_min_overlap": 0.35,
    "mode_min_overlap_margin": 0.08,
    "level_2_trust_radius": 0.05,
    "level_2_use_default_trust": True,
    # ORCA has no standalone CalcAll geometry keyword. Recalc_Hess 1 is the
    # corresponding every-step exact-Hessian fallback.
    "level_2_recalc_hessian_interval": 1,
    "level_2_max_cycles": 30,
    "reseed_from_warmup_on_missing_target_mode": True,
}

_LIVE_MONITOR_DEFAULTS: dict[str, Any] = {
    "enabled": True,
    "poll_seconds": 20.0,
    "min_cycles": 12,
    "window_cycles": 16,
    "min_energy_reversals": 4,
    "energy_change_floor_hartree": 1.0e-5,
    "min_gradient_improvement_fraction": 0.15,
    "step_tolerance_ratio": 4.0,
    "min_excessive_steps": 6,
}

_INT_RESCUE_DEFAULTS: dict[str, Any] = {
    "enabled": True,
    "mode_displacement_first": True,
}

_INT_NONCONVERGENCE_RESCUE_DEFAULTS: dict[str, Any] = {
    "enabled": True,
    "level_1_trust_radius": 0.10,
    "level_1_max_cycles": 30,
    "level_2_trust_radius": 0.05,
    "level_2_use_default_trust": True,
    # The shared stationary-point R2 is the final every-step-Hessian fallback.
    "level_2_recalc_hessian_interval": 1,
    "level_2_max_cycles": 30,
}

_IRC_DEFAULTS: dict[str, Any] = {
    "enabled": True,
    "max_iter": 50,
    "direction": "both",
}

_THERMOCHEMISTRY_DEFAULTS: dict[str, Any] = {
    "temperature_k": 247.55,
    "standard_state": "1M",
    "qrrho": True,
    "ensemble_correction_source": "s1",
}

_STAGE_DEFAULTS: dict[str, dict[str, Any]] = {
    "S3": {
        "fidelity": "low",
        "profile_id": "b97_3c_r2scan_3c_v1",
        "geometry_method": "B97-3c",
        "geometry_basis": "",
        "geometry_aux_basis": "",
        "geometry_grid": None,
        "geometry_scf": None,
        "route_minimum": "Opt",
        "route_ts": "OptTS",
        "max_cycles_minimum": 100,
        "max_cycles_intermediate": 30,
        "max_cycles_ts": 30,
        "ts_trust_radius": None,
        "warmup_max_cycles_int": 8,
        "warmup_max_cycles_ts": 12,
        "frequency_method": "B97-3c",
        "frequency_basis": "",
        "sp_method": "r2SCAN-3c",
        "sp_basis": "",
        "sp_aux_basis": "",
        "solvent": "acetone",
        "solvent_model": "CPCM",
        "timeout_seconds": _RESOURCE_DEFAULTS["timeout_seconds"],
    },
    "S4": {
        "fidelity": "high",
        "profile_id": "m062x_wb97mv_v1",
        "geometry_method": "M062X",
        "geometry_basis": "def2-SVP",
        "geometry_aux_basis": "def2/J",
        "geometry_grid": "DefGrid3",
        "geometry_scf": "TightSCF",
        "route_minimum": "Opt",
        "route_ts": "OptTS",
        "max_cycles_minimum": 200,
        "max_cycles_ts": 200,
        "ts_trust_radius": None,
        "warmup_max_cycles_int": 4,
        "warmup_max_cycles_ts": 6,
        "frequency_method": "M062X",
        "frequency_basis": "def2-SVP",
        "sp_method": "wB97M-V",
        "sp_basis": "def2-TZVPP",
        "sp_aux_basis": "def2/J",
        "solvent": "acetone",
        "solvent_model": "CPCM",
        "timeout_seconds": _RESOURCE_DEFAULTS["timeout_seconds"],
    },
}


def _as_mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _mapping_path(root: Mapping[str, Any], *keys: str) -> Mapping[str, Any]:
    current: Any = root
    for key in keys:
        if not isinstance(current, Mapping):
            return {}
        current = current.get(key)
    return current if isinstance(current, Mapping) else {}


def _coerce_str(value: Any, default: str) -> str:
    if value is None:
        return default
    text = str(value).strip()
    if not text and default:
        return default
    return text


def _coerce_optional_str(value: Any, default: Optional[str] = None) -> Optional[str]:
    if value is None:
        return default
    text = str(value).strip()
    return text or default


def _coerce_int(value: Any, default: int) -> int:
    if value is None or isinstance(value, bool):
        return default
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        try:
            return int(value.strip())
        except ValueError:
            return default
    return default


def _coerce_bool(value: Any, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "yes", "on"}:
            return True
        if normalized in {"0", "false", "no", "off"}:
            return False
        return default
    if isinstance(value, (int, float)):
        return bool(value)
    return default


def _coerce_float(value: Any, default: float) -> float:
    if value is None or isinstance(value, bool):
        return default
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return default
    return default


def _coerce_optional_float(value: Any, default: Optional[float]) -> Optional[float]:
    if value is None or isinstance(value, bool):
        return default
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return default
    return default


def _coerce_memory_gb(value: Any, default: int) -> int:
    if value is None or isinstance(value, bool):
        return default
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if not isinstance(value, str):
        return default

    text = value.strip().upper().replace(" ", "")
    try:
        if text.endswith("GB"):
            return int(float(text[:-2]))
        if text.endswith("G"):
            return int(float(text[:-1]))
        if text.endswith("MB"):
            return max(1, ceil(float(text[:-2]) / 1024.0))
        if text.endswith("M"):
            return max(1, ceil(float(text[:-1]) / 1024.0))
        return int(float(text))
    except ValueError:
        return default


def _normalize_stage(stage: str) -> str:
    normalized = stage.strip().upper()
    if normalized not in _STAGE_DEFAULTS:
        raise ValueError(f"Unsupported refinement stage: {stage!r}")
    return normalized


@dataclass(frozen=True)
class FidelityProfile:
    """Immutable snapshot of all S3/S4 refinement-stage differences."""

    stage: str
    fidelity: str
    profile_id: str

    geometry_method: str
    geometry_basis: str
    geometry_aux_basis: str
    geometry_grid: Optional[str]
    geometry_scf: Optional[str]
    route_minimum: str
    route_ts: str
    max_cycles_minimum: int
    max_cycles_intermediate: int
    max_cycles_ts: int
    ts_trust_radius: Optional[float]

    initial_hessian_precursor: str
    initial_hessian_product: str
    initial_hessian_intermediate: str
    initial_hessian_ts: str

    warmup_max_cycles_int: int
    warmup_max_cycles_ts: int

    frequency_method: str
    frequency_basis: str

    sp_method: str
    sp_basis: str
    sp_aux_basis: str

    solvent: str
    solvent_model: str

    ts_rescue_enabled: bool
    ts_rescue_recalc_hessian_interval: int
    ts_rescue_trust_radius: Optional[float]

    live_monitor_enabled: bool
    live_monitor_poll_seconds: float
    live_monitor_min_cycles: int
    live_monitor_window_cycles: int
    live_monitor_min_energy_reversals: int
    live_monitor_energy_change_floor_hartree: float
    live_monitor_min_gradient_improvement_fraction: float
    live_monitor_step_tolerance_ratio: float
    live_monitor_min_excessive_steps: int

    ts_nonconvergence_rescue_enabled: bool
    ts_nonconvergence_l1_trust_radius: float
    ts_nonconvergence_l1_max_cycles: int
    ts_nonconvergence_mode_min_overlap: float
    ts_nonconvergence_mode_min_overlap_margin: float
    ts_nonconvergence_l2_trust_radius: float
    ts_nonconvergence_l2_use_default_trust: bool
    ts_nonconvergence_l2_recalc_hessian_interval: int
    ts_nonconvergence_l2_max_cycles: int
    ts_nonconvergence_reseed_from_warmup_on_missing_target_mode: bool

    int_rescue_enabled: bool
    int_rescue_mode_displacement: bool

    int_nonconvergence_rescue_enabled: bool
    int_nonconvergence_l1_trust_radius: float
    int_nonconvergence_l1_max_cycles: int
    int_nonconvergence_l2_trust_radius: float
    int_nonconvergence_l2_use_default_trust: bool
    int_nonconvergence_l2_recalc_hessian_interval: int
    int_nonconvergence_l2_max_cycles: int

    irc_enabled: bool
    irc_max_iter: int
    irc_direction: str

    cores_per_worker: int
    memory_gb_per_worker: int
    max_workers: int
    timeout_seconds: int

    temperature_k: float
    standard_state: str
    qrrho: bool
    ensemble_correction_source: str

    @classmethod
    def from_config(cls, config: Mapping[str, Any], stage: str) -> FidelityProfile:
        """Build a profile from new ``refinement.*`` config or legacy theory keys."""

        normalized_stage = _normalize_stage(stage)
        stage_defaults = _STAGE_DEFAULTS[normalized_stage]
        legacy_key = "s3_low_level" if normalized_stage == "S3" else "s4_high_precision"
        legacy_stage_cfg = _mapping_path(config, "theory", legacy_key)

        refinement_cfg = _mapping_path(config, "refinement")
        common_cfg = _mapping_path(refinement_cfg, "common")
        stage_name = normalized_stage.lower()
        # Phase 1 dual-write period: defaults.yaml contains both legacy theory.*
        # and refinement.* blocks. Until orchestrator cutover, callers loading the
        # canonical mixed config should continue to receive the legacy profile.
        has_refinement_stage = isinstance(refinement_cfg.get(stage_name), Mapping) and not legacy_stage_cfg
        refinement_stage_cfg = _as_mapping(refinement_cfg.get(stage_name))

        if has_refinement_stage:
            stage_cfg = refinement_stage_cfg
            geometry_cfg = _mapping_path(stage_cfg, "geometry")
            frequency_cfg = _mapping_path(stage_cfg, "frequency")
            single_point_cfg = _mapping_path(stage_cfg, "single_point")
            solvent_cfg = _mapping_path(stage_cfg, "solvent")
            warmup_cycles_cfg = _mapping_path(stage_cfg, "warmup", "max_cycles")
            stage_resources_cfg = _mapping_path(stage_cfg, "resources")
            timeout_source = geometry_cfg.get("timeout")
            if timeout_source is None:
                timeout_source = single_point_cfg.get("timeout")
            max_cycles_minimum = _coerce_int(
                geometry_cfg.get("max_cycles_minimum"),
                _coerce_int(geometry_cfg.get("max_cycles"), stage_defaults["max_cycles_minimum"]),
            )
            max_cycles_ts = _coerce_int(
                geometry_cfg.get("max_cycles_ts"),
                _coerce_int(geometry_cfg.get("max_cycles"), stage_defaults["max_cycles_ts"]),
            )
            max_cycles_intermediate = _coerce_int(
                geometry_cfg.get("max_cycles_intermediate"),
                max_cycles_minimum,
            )
        else:
            stage_cfg = legacy_stage_cfg
            geometry_cfg = _mapping_path(stage_cfg, "optimization")
            frequency_cfg = {}
            single_point_cfg = _mapping_path(stage_cfg, "single_point")
            solvent_cfg = {}
            warmup_cycles_cfg = _mapping_path(geometry_cfg, "warmup", "max_cycles")
            stage_resources_cfg = {}
            timeout_source = geometry_cfg.get("timeout")
            if timeout_source is None:
                timeout_source = single_point_cfg.get("timeout")
            max_cycles_minimum = _coerce_int(
                geometry_cfg.get("max_cycles_minimum"),
                _coerce_int(geometry_cfg.get("max_cycles"), stage_defaults["max_cycles_minimum"]),
            )
            max_cycles_ts = _coerce_int(
                geometry_cfg.get("max_cycles_ts"),
                stage_defaults["max_cycles_ts"],
            )
            max_cycles_intermediate = _coerce_int(
                geometry_cfg.get("max_cycles_intermediate"),
                max_cycles_minimum,
            )

        workflow_cfg = _mapping_path(common_cfg, "workflow")
        initial_hessian_cfg = _mapping_path(workflow_cfg, "initial_hessian")
        rescue_cfg = _mapping_path(common_cfg, "rescue")
        ts_rescue_cfg = _mapping_path(rescue_cfg, "ts_rescue")
        ts_nonconvergence_cfg = _mapping_path(rescue_cfg, "ts_nonconvergence_rescue")
        live_monitor_cfg = _mapping_path(rescue_cfg, "live_monitor")
        int_rescue_cfg = _mapping_path(rescue_cfg, "int_rescue")
        int_nonconvergence_cfg = _mapping_path(rescue_cfg, "int_nonconvergence_rescue")
        irc_cfg = _mapping_path(common_cfg, "irc")
        thermo_cfg = _mapping_path(common_cfg, "thermochemistry")

        trust_radius_value = ts_rescue_cfg.get("trust_radius")
        if trust_radius_value is None:
            trust_radius_value = ts_rescue_cfg.get("trust")

        step_key = "step3" if normalized_stage == "S3" else "step4"
        scheduling_cfg = _mapping_path(config, step_key, "scheduling")
        global_resources_cfg = _mapping_path(config, "resources")

        cores_per_worker = _coerce_int(
            stage_resources_cfg.get("cores_per_worker"),
            _coerce_int(
                scheduling_cfg.get("nproc_per_job"),
                _coerce_int(
                    scheduling_cfg.get("nproc"),
                    _coerce_int(
                        global_resources_cfg.get("nproc"),
                        _RESOURCE_DEFAULTS["cores_per_worker"],
                    ),
                ),
            ),
        )
        memory_gb_per_worker = _coerce_memory_gb(
            stage_resources_cfg.get("memory_gb_per_worker"),
            _coerce_memory_gb(
                scheduling_cfg.get("memory_per_job"),
                _coerce_memory_gb(
                    global_resources_cfg.get("mem"),
                    _RESOURCE_DEFAULTS["memory_gb_per_worker"],
                ),
            ),
        )
        max_workers = _coerce_int(
            stage_resources_cfg.get("max_workers"),
            _coerce_int(
                scheduling_cfg.get("max_workers"),
                _RESOURCE_DEFAULTS["max_workers"],
            ),
        )
        timeout_seconds = _coerce_int(
            timeout_source,
            stage_defaults["timeout_seconds"],
        )

        solvent = _coerce_str(
            solvent_cfg.get("solvent", geometry_cfg.get("solvent", single_point_cfg.get("solvent"))),
            stage_defaults["solvent"],
        )
        solvent_model = _coerce_str(
            solvent_cfg.get(
                "solvent_model",
                geometry_cfg.get("solvent_model", single_point_cfg.get("solvent_model")),
            ),
            stage_defaults["solvent_model"],
        )

        temperature_value = thermo_cfg.get("temperature_K")
        if temperature_value is None:
            temperature_value = thermo_cfg.get("temperature_k")

        return cls(
            stage=normalized_stage,
            fidelity=_coerce_str(stage_cfg.get("fidelity"), stage_defaults["fidelity"]),
            profile_id=_coerce_str(stage_cfg.get("profile_id"), stage_defaults["profile_id"]),
            geometry_method=_coerce_str(geometry_cfg.get("method"), stage_defaults["geometry_method"]),
            geometry_basis=_coerce_str(geometry_cfg.get("basis"), stage_defaults["geometry_basis"]),
            geometry_aux_basis=_coerce_str(
                geometry_cfg.get("aux_basis"),
                stage_defaults["geometry_aux_basis"],
            ),
            geometry_grid=_coerce_optional_str(
                geometry_cfg.get("grid"),
                stage_defaults["geometry_grid"],
            ),
            geometry_scf=_coerce_optional_str(
                geometry_cfg.get("scf"),
                stage_defaults["geometry_scf"],
            ),
            route_minimum=_coerce_str(
                geometry_cfg.get("route_minimum"),
                stage_defaults["route_minimum"],
            ),
            route_ts=_coerce_str(geometry_cfg.get("route_ts"), stage_defaults["route_ts"]),
            max_cycles_minimum=max_cycles_minimum,
            max_cycles_intermediate=max_cycles_intermediate,
            max_cycles_ts=max_cycles_ts,
            ts_trust_radius=_coerce_optional_float(
                geometry_cfg.get("ts_trust_radius"),
                stage_defaults["ts_trust_radius"],
            ),
            initial_hessian_precursor=_coerce_str(
                initial_hessian_cfg.get("precursor"),
                _INITIAL_HESSIAN_DEFAULTS["precursor"],
            ),
            initial_hessian_product=_coerce_str(
                initial_hessian_cfg.get("product"),
                _INITIAL_HESSIAN_DEFAULTS["product"],
            ),
            initial_hessian_intermediate=_coerce_str(
                initial_hessian_cfg.get("intermediate"),
                _INITIAL_HESSIAN_DEFAULTS["intermediate"],
            ),
            initial_hessian_ts=_coerce_str(
                initial_hessian_cfg.get("ts"),
                _INITIAL_HESSIAN_DEFAULTS["ts"],
            ),
            warmup_max_cycles_int=_coerce_int(
                warmup_cycles_cfg.get("intermediate"),
                stage_defaults["warmup_max_cycles_int"],
            ),
            warmup_max_cycles_ts=_coerce_int(
                warmup_cycles_cfg.get("ts"),
                stage_defaults["warmup_max_cycles_ts"],
            ),
            frequency_method=_coerce_str(
                frequency_cfg.get("method"),
                stage_defaults["frequency_method"],
            ),
            frequency_basis=_coerce_str(
                frequency_cfg.get("basis"),
                stage_defaults["frequency_basis"],
            ),
            sp_method=_coerce_str(single_point_cfg.get("method"), stage_defaults["sp_method"]),
            sp_basis=_coerce_str(single_point_cfg.get("basis"), stage_defaults["sp_basis"]),
            sp_aux_basis=_coerce_str(
                single_point_cfg.get("aux_basis"),
                stage_defaults["sp_aux_basis"],
            ),
            solvent=solvent,
            solvent_model=solvent_model,
            ts_rescue_enabled=_coerce_bool(
                ts_rescue_cfg.get("enabled"),
                _TS_RESCUE_DEFAULTS["enabled"],
            ),
            ts_rescue_recalc_hessian_interval=_coerce_int(
                ts_rescue_cfg.get("recalc_hessian_interval"),
                _TS_RESCUE_DEFAULTS["recalc_hessian_interval"],
            ),
            ts_rescue_trust_radius=_coerce_optional_float(
                trust_radius_value,
                _TS_RESCUE_DEFAULTS["trust_radius"],
            ),
            live_monitor_enabled=_coerce_bool(
                live_monitor_cfg.get("enabled"), _LIVE_MONITOR_DEFAULTS["enabled"]
            ),
            live_monitor_poll_seconds=_coerce_float(
                live_monitor_cfg.get("poll_seconds"), _LIVE_MONITOR_DEFAULTS["poll_seconds"]
            ),
            live_monitor_min_cycles=_coerce_int(
                live_monitor_cfg.get("min_cycles"), _LIVE_MONITOR_DEFAULTS["min_cycles"]
            ),
            live_monitor_window_cycles=_coerce_int(
                live_monitor_cfg.get("window_cycles"), _LIVE_MONITOR_DEFAULTS["window_cycles"]
            ),
            live_monitor_min_energy_reversals=_coerce_int(
                live_monitor_cfg.get("min_energy_reversals"),
                _LIVE_MONITOR_DEFAULTS["min_energy_reversals"],
            ),
            live_monitor_energy_change_floor_hartree=_coerce_float(
                live_monitor_cfg.get("energy_change_floor_hartree"),
                _LIVE_MONITOR_DEFAULTS["energy_change_floor_hartree"],
            ),
            live_monitor_min_gradient_improvement_fraction=_coerce_float(
                live_monitor_cfg.get("min_gradient_improvement_fraction"),
                _LIVE_MONITOR_DEFAULTS["min_gradient_improvement_fraction"],
            ),
            live_monitor_step_tolerance_ratio=_coerce_float(
                live_monitor_cfg.get("step_tolerance_ratio"),
                _LIVE_MONITOR_DEFAULTS["step_tolerance_ratio"],
            ),
            live_monitor_min_excessive_steps=_coerce_int(
                live_monitor_cfg.get("min_excessive_steps"),
                _LIVE_MONITOR_DEFAULTS["min_excessive_steps"],
            ),
            ts_nonconvergence_rescue_enabled=_coerce_bool(
                ts_nonconvergence_cfg.get("enabled"),
                _TS_NONCONVERGENCE_RESCUE_DEFAULTS["enabled"],
            ),
            ts_nonconvergence_l1_trust_radius=_coerce_float(
                ts_nonconvergence_cfg.get("level_1_trust_radius"),
                _TS_NONCONVERGENCE_RESCUE_DEFAULTS["level_1_trust_radius"],
            ),
            ts_nonconvergence_l1_max_cycles=_coerce_int(
                ts_nonconvergence_cfg.get("level_1_max_cycles"),
                _TS_NONCONVERGENCE_RESCUE_DEFAULTS["level_1_max_cycles"],
            ),
            ts_nonconvergence_mode_min_overlap=_coerce_float(
                ts_nonconvergence_cfg.get("mode_min_overlap"),
                _TS_NONCONVERGENCE_RESCUE_DEFAULTS["mode_min_overlap"],
            ),
            ts_nonconvergence_mode_min_overlap_margin=_coerce_float(
                ts_nonconvergence_cfg.get("mode_min_overlap_margin"),
                _TS_NONCONVERGENCE_RESCUE_DEFAULTS["mode_min_overlap_margin"],
            ),
            ts_nonconvergence_l2_trust_radius=_coerce_float(
                ts_nonconvergence_cfg.get("level_2_trust_radius"),
                _TS_NONCONVERGENCE_RESCUE_DEFAULTS["level_2_trust_radius"],
            ),
            ts_nonconvergence_l2_use_default_trust=_coerce_bool(
                ts_nonconvergence_cfg.get("level_2_use_default_trust"),
                _TS_NONCONVERGENCE_RESCUE_DEFAULTS["level_2_use_default_trust"],
            ),
            ts_nonconvergence_l2_recalc_hessian_interval=_coerce_int(
                ts_nonconvergence_cfg.get("level_2_recalc_hessian_interval"),
                _TS_NONCONVERGENCE_RESCUE_DEFAULTS["level_2_recalc_hessian_interval"],
            ),
            ts_nonconvergence_l2_max_cycles=_coerce_int(
                ts_nonconvergence_cfg.get("level_2_max_cycles"),
                _TS_NONCONVERGENCE_RESCUE_DEFAULTS["level_2_max_cycles"],
            ),
            ts_nonconvergence_reseed_from_warmup_on_missing_target_mode=_coerce_bool(
                ts_nonconvergence_cfg.get("reseed_from_warmup_on_missing_target_mode"),
                _TS_NONCONVERGENCE_RESCUE_DEFAULTS[
                    "reseed_from_warmup_on_missing_target_mode"
                ],
            ),
            int_rescue_enabled=_coerce_bool(
                int_rescue_cfg.get("enabled"),
                _INT_RESCUE_DEFAULTS["enabled"],
            ),
            int_rescue_mode_displacement=_coerce_bool(
                int_rescue_cfg.get("mode_displacement_first"),
                _INT_RESCUE_DEFAULTS["mode_displacement_first"],
            ),
            int_nonconvergence_rescue_enabled=_coerce_bool(
                int_nonconvergence_cfg.get("enabled"),
                _INT_NONCONVERGENCE_RESCUE_DEFAULTS["enabled"],
            ),
            int_nonconvergence_l1_trust_radius=_coerce_float(
                int_nonconvergence_cfg.get("level_1_trust_radius"),
                _INT_NONCONVERGENCE_RESCUE_DEFAULTS["level_1_trust_radius"],
            ),
            int_nonconvergence_l1_max_cycles=_coerce_int(
                int_nonconvergence_cfg.get("level_1_max_cycles"),
                _INT_NONCONVERGENCE_RESCUE_DEFAULTS["level_1_max_cycles"],
            ),
            int_nonconvergence_l2_trust_radius=_coerce_float(
                int_nonconvergence_cfg.get("level_2_trust_radius"),
                _INT_NONCONVERGENCE_RESCUE_DEFAULTS["level_2_trust_radius"],
            ),
            int_nonconvergence_l2_use_default_trust=_coerce_bool(
                int_nonconvergence_cfg.get("level_2_use_default_trust"),
                _INT_NONCONVERGENCE_RESCUE_DEFAULTS["level_2_use_default_trust"],
            ),
            int_nonconvergence_l2_recalc_hessian_interval=_coerce_int(
                int_nonconvergence_cfg.get("level_2_recalc_hessian_interval"),
                _INT_NONCONVERGENCE_RESCUE_DEFAULTS[
                    "level_2_recalc_hessian_interval"
                ],
            ),
            int_nonconvergence_l2_max_cycles=_coerce_int(
                int_nonconvergence_cfg.get("level_2_max_cycles"),
                _INT_NONCONVERGENCE_RESCUE_DEFAULTS["level_2_max_cycles"],
            ),
            irc_enabled=_coerce_bool(irc_cfg.get("enabled"), _IRC_DEFAULTS["enabled"]),
            irc_max_iter=_coerce_int(irc_cfg.get("max_iter"), _IRC_DEFAULTS["max_iter"]),
            irc_direction=_coerce_str(irc_cfg.get("direction"), _IRC_DEFAULTS["direction"]),
            cores_per_worker=cores_per_worker,
            memory_gb_per_worker=memory_gb_per_worker,
            max_workers=max_workers,
            timeout_seconds=timeout_seconds,
            temperature_k=_coerce_float(
                temperature_value,
                _THERMOCHEMISTRY_DEFAULTS["temperature_k"],
            ),
            standard_state=_coerce_str(
                thermo_cfg.get("standard_state"),
                _THERMOCHEMISTRY_DEFAULTS["standard_state"],
            ),
            qrrho=_coerce_bool(thermo_cfg.get("qrrho"), _THERMOCHEMISTRY_DEFAULTS["qrrho"]),
            ensemble_correction_source=_coerce_str(
                thermo_cfg.get("ensemble_correction_source"),
                _THERMOCHEMISTRY_DEFAULTS["ensemble_correction_source"],
            ),
        )

    def initial_hessian_for_role(self, role: str) -> str:
        """Return the initial-Hessian policy for a chemical role."""

        normalized_role = role.strip().lower()
        if normalized_role == "precursor":
            return self.initial_hessian_precursor
        if normalized_role == "product":
            return self.initial_hessian_product
        if normalized_role == "intermediate":
            return self.initial_hessian_intermediate
        if normalized_role == "ts":
            return self.initial_hessian_ts
        raise ValueError(f"Unknown refinement role: {role!r}")


__all__ = ["FidelityProfile"]
