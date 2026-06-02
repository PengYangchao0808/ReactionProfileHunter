"""
Unified JSON I/O helpers.

Single source of truth for reading/writing JSON files across RPH.
Replaces duplicated _read_json / _write_json in orchestrator, condition_thermo,
condition_feature_merger, and other modules.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


def read_json(path: Path, default: Any = None) -> Any:
    """Read JSON from path, returning `default` on any failure.

    Args:
        path: JSON file path.
        default: Value to return if file missing or parse fails.

    Returns:
        Parsed JSON object, or `default`.
    """
    if not path.exists():
        return default
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as exc:
        logger.warning(f"Failed to read JSON {path}: {exc}")
        return default


def read_json_dict(path: Path) -> Dict[str, Any]:
    """Read JSON, guarantee a dict return (empty dict on failure)."""
    result = read_json(path, default={})
    return result if isinstance(result, dict) else {}


def write_json(path: Path, payload: Any) -> None:
    """Write JSON to path atomically (write-then-rename where safe).

    Creates parent directories automatically.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        tmp = path.with_suffix(path.suffix + ".tmp")
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False, default=str)
        tmp.replace(path)
    except Exception as exc:
        logger.warning(f"Failed to write JSON {path}: {exc}")
