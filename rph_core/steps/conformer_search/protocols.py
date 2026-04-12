from dataclasses import asdict, dataclass
from typing import Any, Dict, Optional


SUPPORTED_PROTOCOLS = {"ext", "default", "full", "lite", "zero"}


def _optional_float(value: Any, default: Optional[float] = None) -> Optional[float]:
    if value is None:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _optional_int(value: Any, default: Optional[int] = None) -> Optional[int]:
    if value is None:
        return default
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


@dataclass(frozen=True)
class FunnelPolicy:
    search_mode: str
    clustering_mode: str
    prescreen_mode: str
    rerank_mode: str
    use_mrrho_like_correction: bool
    survivor_window_kcal: Optional[float]
    narrow_window_kcal: Optional[float]
    prescreen_window_kcal: Optional[float]
    screening_window_kcal: Optional[float]
    optimize_limit: Optional[int]
    top2_fallback_enabled: bool
    boltzmann_cutoff: Optional[float]


@dataclass(frozen=True)
class HandoffPolicy:
    enabled: bool
    mode: str
    fallback_mode: Optional[str]
    small_gap_kcal: Optional[float]
    ranking_after_handoff: str


@dataclass(frozen=True)
class ProtocolSpec:
    name: str
    two_stage_enabled: bool
    ngeom_default: int
    ngeom_max: int
    funnel_policy: FunnelPolicy
    handoff_policy: HandoffPolicy
    final_opt_sp_enabled: bool
    freq_enabled: bool
    final_sp_enabled: bool
    selection_mode: str
    provenance_flags: Dict[str, bool]
    raw_config: Dict[str, Any]


def resolve_protocol_spec(config: Dict[str, Any], protocol: str) -> ProtocolSpec:
    normalized = (protocol or "ext").strip().lower()

    if normalized not in SUPPORTED_PROTOCOLS:
        raise RuntimeError(
            f"Unknown step1.protocol '{protocol}'. "
            "Allowed values currently: ext, default, full, lite, zero."
        )

    step1 = config.get("step1", {})
    conformer_search = step1.get("conformer_search", {}) if isinstance(step1, dict) else {}
    protocol_stack = step1.get("protocol_stack", {}) if isinstance(step1, dict) else {}
    raw_cfg = protocol_stack.get(normalized, {}) if isinstance(protocol_stack, dict) else {}
    proto_cfg = raw_cfg if isinstance(raw_cfg, dict) else {}
    final_opt_sp = proto_cfg.get("final_opt_sp", {}) if isinstance(proto_cfg, dict) else {}
    final_opt_sp_cfg = final_opt_sp if isinstance(final_opt_sp, dict) else {}
    funnel_cfg = proto_cfg.get("funnel", {}) if isinstance(proto_cfg, dict) else {}
    funnel_cfg = funnel_cfg if isinstance(funnel_cfg, dict) else {}
    handoff_cfg = proto_cfg.get("handoff", {}) if isinstance(proto_cfg, dict) else {}
    handoff_cfg = handoff_cfg if isinstance(handoff_cfg, dict) else {}
    shared_handoff = step1.get("shared_handoff", {}) if isinstance(step1, dict) else {}
    shared_handoff = shared_handoff if isinstance(shared_handoff, dict) else {}

    is_light = normalized in {"lite", "zero"}
    is_full = normalized == "full"
    default_two_stage = bool(conformer_search.get("two_stage_enabled", True))
    default_ngeom_default = int(conformer_search.get("ngeom_default", 3))
    default_ngeom_max = int(conformer_search.get("ngeom_max", 6))

    final_opt_sp_enabled = bool(
        handoff_cfg.get(
            "enabled",
            shared_handoff.get("enabled", final_opt_sp_cfg.get("enabled", True)),
        )
    )
    default_freq_enabled = bool(shared_handoff.get("dft_optfreq", final_opt_sp_cfg.get("freq", True)))
    freq_enabled = bool(final_opt_sp_cfg.get("freq", default_freq_enabled))
    final_sp_enabled = bool(
        shared_handoff.get(
            "final_sp",
            final_opt_sp_cfg.get("final_sp_enabled", final_opt_sp_enabled),
        )
    )
    selection_mode = str(
        handoff_cfg.get("mode", final_opt_sp_cfg.get("selection_mode", proto_cfg.get("selection_mode", "default")))
    )

    default_search_mode = {
        "ext": "crest_two_stage_gfn0_to_gfn2",
        "default": "crest_two_stage_gfn0_to_gfn2",
        "full": "crest_gfn2",
        "lite": "crest_gfn2",
        "zero": "crest_gfn2_or_skip",
    }[normalized]
    default_handoff_mode = {
        "ext": "optimize_all_candidates",
        "default": "optimize_all_candidates",
        "full": "optimize_all_survivors_within_window",
        "lite": "optimize_rank1",
        "zero": "optimize_rank1",
    }[normalized]
    default_fallback = {
        "ext": None,
        "default": None,
        "full": None,
        "lite": "optimize_top2_if_gap_small",
        "zero": "optimize_all_within_0p5_kcal",
    }[normalized]

    funnel_policy = FunnelPolicy(
        search_mode=str(funnel_cfg.get("search_mode", default_search_mode)),
        clustering_mode=str(funnel_cfg.get("clustering_mode", "isostat")),
        prescreen_mode=str(funnel_cfg.get("prescreen_mode", "none")),
        rerank_mode=str(funnel_cfg.get("rerank_mode", "none")),
        use_mrrho_like_correction=bool(funnel_cfg.get("use_mrrho_like_correction", is_full or normalized == "lite")),
        survivor_window_kcal=_optional_float(funnel_cfg.get("survivor_window_kcal"), 3.0 if is_full else None),
        narrow_window_kcal=_optional_float(funnel_cfg.get("narrow_window_kcal"), 0.5 if normalized == "zero" else None),
        prescreen_window_kcal=_optional_float(funnel_cfg.get("prescreen_window_kcal"), 4.0 if is_full else None),
        screening_window_kcal=_optional_float(funnel_cfg.get("screening_window_kcal"), 3.5 if is_full else None),
        optimize_limit=_optional_int(funnel_cfg.get("optimize_limit"), 1 if normalized in {"lite", "zero"} else None),
        top2_fallback_enabled=bool(funnel_cfg.get("top2_fallback_enabled", normalized == "lite")),
        boltzmann_cutoff=_optional_float(funnel_cfg.get("boltzmann_cutoff"), 0.90 if normalized == "lite" else None),
    )

    handoff_policy = HandoffPolicy(
        enabled=final_opt_sp_enabled,
        mode=str(handoff_cfg.get("mode", default_handoff_mode)),
        fallback_mode=(
            str(handoff_cfg.get("fallback_mode"))
            if handoff_cfg.get("fallback_mode") is not None
            else default_fallback
        ),
        small_gap_kcal=_optional_float(handoff_cfg.get("small_gap_kcal"), 1.0 if normalized == "lite" else None),
        ranking_after_handoff=str(
            handoff_cfg.get(
                "ranking_after_handoff",
                shared_handoff.get("ranking_after_handoff", "final_sp_minimum"),
            )
        ),
    )

    provenance_flags = {
        "write_stage_meta": bool(final_opt_sp_cfg.get("write_stage_meta", True)),
        "write_provenance_flags": bool(final_opt_sp_cfg.get("write_provenance_flags", True)),
        "protocol_is_lightweight": is_light,
        "shared_handoff": True,
    }

    return ProtocolSpec(
        name=normalized,
        two_stage_enabled=bool(
            proto_cfg.get(
                "two_stage_enabled",
                default_two_stage if normalized in {"ext", "default"} else False,
            )
        ),
        ngeom_default=int(
            proto_cfg.get(
                "ngeom_default",
                12 if normalized == "full" else (default_ngeom_default if normalized in {"ext", "default"} else 1),
            )
        ),
        ngeom_max=int(
            proto_cfg.get(
                "ngeom_max",
                24 if normalized == "full" else (default_ngeom_max if normalized in {"ext", "default"} else 1),
            )
        ),
        funnel_policy=funnel_policy,
        handoff_policy=handoff_policy,
        final_opt_sp_enabled=final_opt_sp_enabled,
        freq_enabled=freq_enabled,
        final_sp_enabled=final_sp_enabled,
        selection_mode=selection_mode,
        provenance_flags=provenance_flags,
        raw_config=proto_cfg,
    )


def as_provenance_payload(spec: ProtocolSpec) -> Dict[str, Any]:
    payload = asdict(spec)
    return payload


def is_ext_protocol(spec: ProtocolSpec) -> bool:
    return spec.name in {"ext", "default"}


def is_full_protocol(spec: ProtocolSpec) -> bool:
    return spec.name == "full"


def is_lite_protocol(spec: ProtocolSpec) -> bool:
    return spec.name == "lite"


def is_zero_protocol(spec: ProtocolSpec) -> bool:
    return spec.name == "zero"
