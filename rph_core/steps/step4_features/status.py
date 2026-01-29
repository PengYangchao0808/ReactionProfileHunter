"""
Step 4: Feature Status Enums (V4.2 Phase A)
==========================================

FeatureStatus (global) and FeatureResultStatus (per-plugin) enums
with helper functions.

Author: QC Descriptors Team
Date: 2026-01-18
"""

from enum import Enum
from typing import Optional


class NICSTriggerPolicy(Enum):
    """NICS (Nucleus Independent Chemical Shift) trigger policy.

    Phase B/C: determines whether to generate NICS descriptors.
    Phase A: unused but locked in interface.
    """
    AUTO_RUN = "auto"
    GENERATE_ONLY = "gen"
    SKIP = "skip"


class FeatureStatus(Enum):
    """Global feature extraction status for entire Step4 pipeline.

    This status is written to `features_raw.csv` as `feature_status` column.
    """
    OK = "ok"
    PARTIAL = "partial"
    FAILED = "failed"
    MISSING_INPUTS = "missing_inputs"
    INVALID_INPUTS = "invalid_inputs"
    SKIPPED = "skipped"


class FeatureResultStatus(Enum):
    """Per-plugin feature extraction result status.

    Each plugin reports its own status, which is aggregated into
    `feature_meta.json` under `trace.plugins.<plugin_name>.status`.
    """
    OK = "ok"
    SKIPPED = "skipped"
    FAILED = "failed"
    PARTIAL = "partial"


def is_global_success(status: FeatureStatus) -> bool:
    """Check if global status indicates successful extraction.

    Args:
        status: FeatureStatus enum value

    Returns:
        True if status is OK or PARTIAL, False otherwise
    """
    return status in (FeatureStatus.OK, FeatureStatus.PARTIAL)


def is_plugin_success(status: FeatureResultStatus) -> bool:
    """Check if plugin status indicates successful extraction.

    Args:
        status: FeatureResultStatus enum value

    Returns:
        True if status is OK or PARTIAL, False otherwise
    """
    return status in (FeatureResultStatus.OK, FeatureResultStatus.PARTIAL)


def aggregate_plugin_status(plugin_statuses: dict) -> FeatureStatus:
    """Aggregate per-plugin statuses into global FeatureStatus.

    Aggregation rules:
    - If all plugins are OK → OK
    - If any plugins are PARTIAL and none FAILED → PARTIAL
    - If any plugins are FAILED → FAILED
    - If all plugins are SKIPPED → SKIPPED

    Args:
        plugin_statuses: Dictionary mapping plugin_name -> FeatureResultStatus

    Returns:
        Aggregated global FeatureStatus
    """
    if not plugin_statuses:
        return FeatureStatus.SKIPPED

    statuses = list(plugin_statuses.values())

    # Check for FAILED
    if FeatureResultStatus.FAILED in statuses:
        return FeatureStatus.FAILED

    # Check for SKIPPED only
    if all(s == FeatureResultStatus.SKIPPED for s in statuses):
        return FeatureStatus.SKIPPED

    # Check for PARTIAL
    if FeatureResultStatus.PARTIAL in statuses:
        return FeatureStatus.PARTIAL

    # All OK
    return FeatureStatus.OK


def status_from_error(error: Exception) -> FeatureResultStatus:
    """Convert an exception to a FeatureResultStatus.

    Args:
        error: Exception that occurred during plugin execution

    Returns:
        FeatureResultStatus.FAILED for any exception
    """
    return FeatureResultStatus.FAILED


def status_from_missing_inputs(missing_keys: list) -> FeatureStatus:
    """Create FeatureStatus based on missing input keys.

    Args:
        missing_keys: List of missing context keys

    Returns:
        FeatureStatus.MISSING_INPUTS if any keys are missing
    """
    if missing_keys:
        return FeatureStatus.MISSING_INPUTS
    return FeatureStatus.OK


__all__ = [
    "FeatureStatus",
    "FeatureResultStatus",
    "NICSTriggerPolicy",
    "is_global_success",
    "is_plugin_success",
    "aggregate_plugin_status",
    "status_from_error",
    "status_from_missing_inputs",
]
