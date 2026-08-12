"""Pure helpers for property-level ML usability classification."""

from __future__ import annotations

from typing import Dict, Optional


def compute_ml_usability(
    *,
    requested_role: str,
    requested_kind: str,
    opt_converged: bool,
    sp_complete: bool,
    freq_available: bool,
    identity_status: str,
    ts_frequency_valid: Optional[bool],
    geometry_consistent: Optional[bool],
) -> Dict[str, bool]:
    """Compute property-level ML usability."""

    _ = geometry_consistent

    return {
        "geometry": bool(opt_converged),
        "electronic_energy": bool(sp_complete),
        "enthalpy": bool(sp_complete and freq_available),
        "gibbs_free_energy": bool(sp_complete and freq_available),
        "frequency_descriptors": bool(freq_available),
        "ts_descriptors": bool(requested_kind == "ts" and ts_frequency_valid is True),
        "intermediate_descriptors": bool(
            requested_role == "intermediate" and identity_status == "role_matched"
        ),
        "mechanism_label": True,
        "multifidelity_pair": False,
    }


def compute_legacy_usable_for_ml(ml_usability: Dict[str, bool]) -> bool:
    """Map property-level usability back to the legacy single boolean."""

    return bool(ml_usability.get("geometry") and ml_usability.get("electronic_energy"))


__all__ = ["compute_ml_usability", "compute_legacy_usable_for_ml"]
