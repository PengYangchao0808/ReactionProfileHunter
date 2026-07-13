"""V4-only S1 protocol contract.

The former protocol stack selected several conformer workflows.  V4 has one
fixed CENSO-LITE implementation; this small public contract remains for
validation and callers that need to assert that invariant.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict


SUPPORTED_PROTOCOLS = {"censo_lite"}


@dataclass(frozen=True)
class FunnelPolicy:
    """The fixed conformer-selection policy used by CENSO-LITE."""

    clustering_mode: str = "torsion_signature"


@dataclass(frozen=True)
class ProtocolSpec:
    """Minimal immutable description of the sole supported S1 protocol."""

    name: str
    funnel_policy: FunnelPolicy


def resolve_protocol_spec(config: Dict[str, Any], protocol: str) -> ProtocolSpec:
    """Validate and describe the only S1 protocol available in V4."""

    normalized = str(protocol or "").strip().lower()
    if normalized not in SUPPORTED_PROTOCOLS:
        raise RuntimeError(
            f"Unknown step1.protocol '{protocol}'. Allowed value: censo_lite."
        )
    return ProtocolSpec(name="censo_lite", funnel_policy=FunnelPolicy())
