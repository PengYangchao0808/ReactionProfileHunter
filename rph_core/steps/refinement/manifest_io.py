"""Read/write helpers for refinement manifests.

Supports the unified `refinement_manifest_v1` schema AND legacy
`s3_low_level_v3` / `s4_high_level_v4` schemas via read adapters.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict

from rph_core.utils.json_io import write_json_atomic


logger = logging.getLogger(__name__)


REFINEMENT_MANIFEST_V1 = "refinement_manifest_v1"
LEGACY_S3_SCHEMA = "s3_low_level_v3"
LEGACY_S4_SCHEMA = "s4_high_level_v4"


def write_refinement_manifest(
    path: Path,
    *,
    stage: str,
    fidelity: str,
    profile_id: str,
    structures: list[Dict[str, Any]],
    run_id: str | None = None,
    extra: Dict[str, Any] | None = None,
) -> Path:
    """Write a refinement_manifest_v1 manifest."""

    path = Path(path)
    payload: Dict[str, Any] = {
        "schema_version": REFINEMENT_MANIFEST_V1,
        "stage": stage,
        "fidelity": fidelity,
        "profile_id": profile_id,
        "run_id": run_id,
        "structures": structures,
    }
    if extra:
        payload.update(extra)
    write_json_atomic(path, payload)
    return path


def read_refinement_manifest(path: Path) -> Dict[str, Any]:
    """Read a refinement manifest, normalizing legacy schemas to v1 shape.

    Accepts:
    - `refinement_manifest_v1` (pass-through)
    - `s3_low_level_v3` (legacy S3) → adapt: stage="S3", fidelity="low", profile_id inferred
    - `s4_high_level_v4` (legacy S4) → adapt: stage="S4", fidelity="high", profile_id inferred

    The returned dict always contains:
    - `schema_version`: the ORIGINAL schema_version (caller can check legacy vs v1)
    - `stage`, `fidelity`, `profile_id`: normalized
    - `structures`: list (legacy S3 stored as dict-keyed — convert to list)
    """

    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    schema = str(data.get("schema_version", ""))

    if schema == REFINEMENT_MANIFEST_V1:
        return data

    if schema == LEGACY_S3_SCHEMA:
        return _adapt_legacy_s3(data)

    if schema == LEGACY_S4_SCHEMA:
        return _adapt_legacy_s4(data)

    logger.warning(
        "Unknown refinement manifest schema_version=%r in %s; returning raw payload",
        schema,
        path,
    )
    return data


def _adapt_legacy_s3(data: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize s3_low_level_v3 → refinement_manifest_v1 shape."""

    normalized = dict(data)
    normalized.setdefault("stage", "S3")
    normalized.setdefault("fidelity", "low")
    normalized.setdefault("profile_id", "b97_3c_r2scan_3c_v1_legacy")
    structures = data.get("structures")
    if isinstance(structures, dict):
        normalized["structures"] = list(structures.values())
    elif structures is None:
        normalized["structures"] = []
    return normalized


def _adapt_legacy_s4(data: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize s4_high_level_v4 → refinement_manifest_v1 shape."""

    normalized = dict(data)
    normalized.setdefault("stage", "S4")
    normalized.setdefault("fidelity", "high")
    normalized.setdefault("profile_id", "m062x_wb97mv_v1_legacy")
    normalized.setdefault("structures", data.get("structures") or [])
    return normalized


__all__ = [
    "LEGACY_S3_SCHEMA",
    "LEGACY_S4_SCHEMA",
    "REFINEMENT_MANIFEST_V1",
    "read_refinement_manifest",
    "write_refinement_manifest",
]
