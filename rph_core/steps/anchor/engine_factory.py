from copy import deepcopy
from pathlib import Path
from typing import Any, Dict

from rph_core.steps.conformer_search.engine import ConformerEngine
from rph_core.steps.conformer_search.protocols import ProtocolSpec, SUPPORTED_PROTOCOLS, resolve_protocol_spec


def create_s1_engine(
    *,
    protocol: str,
    config: Dict[str, Any],
    work_dir: Path,
    molecule_name: str,
) -> ConformerEngine:
    normalized = (protocol or "ext").strip().lower()
    protocol_spec = resolve_protocol_spec(config, normalized)

    if normalized in SUPPORTED_PROTOCOLS:
        effective_config = _build_protocol_config(config, normalized, protocol_spec)
        return ConformerEngine(
            config=effective_config,
            work_dir=work_dir,
            molecule_name=molecule_name,
            protocol_spec=protocol_spec,
        )

    raise RuntimeError(f"Unexpected unsupported protocol after resolution: {protocol}")


def _build_protocol_config(config: Dict[str, Any], protocol: str, protocol_spec: ProtocolSpec) -> Dict[str, Any]:
    merged = deepcopy(config)
    step1 = merged.setdefault("step1", {})
    step1["protocol"] = protocol
    proto_cfg = step1.get("protocol_stack", {}).get(protocol, {}) or {}
    conf = step1.setdefault("conformer_search", {})

    conf["two_stage_enabled"] = bool(proto_cfg.get("two_stage_enabled", protocol_spec.two_stage_enabled))
    conf["ngeom_default"] = int(proto_cfg.get("ngeom_default", protocol_spec.ngeom_default))
    conf["ngeom_max"] = int(proto_cfg.get("ngeom_max", protocol_spec.ngeom_max))
    conf_energy = float(conf.get("energy_window_kcal", 3.0))
    proto_energy = float(proto_cfg.get("conformer_energy_window_kcal", conf_energy))
    conf["energy_window_kcal"] = proto_energy

    crest = step1.setdefault("crest", {})
    crest_energy = float(crest.get("energy_window", 6.0))
    proto_crest_energy = float(proto_cfg.get("crest_energy_window_kcal", crest_energy))
    crest["energy_window"] = proto_crest_energy

    step1["enable_dft_opt"] = True
    return merged
