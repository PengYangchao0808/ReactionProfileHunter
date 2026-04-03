"""
Cache Key Utilities
===================

Shared cache key generation for S3 enrichment and S4 caching.

Author: QCcalc Team
Date: 2026-01-27
"""

import hashlib
import json
from pathlib import Path
from typing import Dict, Any, Optional, List


def generate_step4_plugin_cache_key(
    plugin_name: str,
    input_files: Dict[str, Path],
    params: Optional[Dict[str, Any]] = None
) -> str:
    """
    Generate cache key for Step 4 plugin results (Contract 5).

    Phase B/C: used for caching expensive QM calculations.
    Phase A: interface is locked but caching not yet active.

    Args:
        plugin_name: Name of plugin generating cache
        input_files: Dictionary mapping canonical key -> Path for input files
        params: Additional parameters affecting computation

    Returns:
        Cache key string (sha1[:16] of normalized input)
    """
    key_data = {
        "plugin": plugin_name,
        "files": {},
        "params": params or {}
    }

    sorted_file_keys = sorted(input_files.keys())
    for key in sorted_file_keys:
        file_path = input_files[key]
        if file_path is None or not file_path.exists():
            continue

        try:
            stat = file_path.stat()
            size = stat.st_size
            mtime = stat.st_mtime

            with open(file_path, "rb") as f:
                content = f.read()

            content_hash = hashlib.sha256(content).hexdigest()[:8]

            key_data["files"][key] = {
                "size": size,
                "mtime": mtime,
                "sha256_prefix": content_hash
            }
        except Exception:
            continue

    key_str = json.dumps(key_data, sort_keys=True, default=str)
    full_hash = hashlib.sha1(key_str.encode('utf-8')).hexdigest()
    return full_hash[:16]


def generate_enrichment_cache_key(
    orca_template_hash: str,
    geometry_hashes: Dict[str, str],
    fragment_hash: str
) -> str:
    """
    Generate enrichment cache key from S3 enrichment inputs.

    Content-hash based (no mtime):
    - orca_template_hash: SHA256 hash of rendered ORCA template
    - geometry_hashes: JSON with sha256 hashes of input XYZ files
    - fragment_hash: SHA256 hash of JSON with fragmenter metadata (indices, cut_bond, cap_rules)

    Args:
        orca_template_hash: SHA256 hash of ORCA .inp file content
        geometry_hashes: Dictionary with sha256 hashes of input XYZ files
            - reactant_complex_xyz_sha256: hash of reactant complex XYZ bytes
            - ts_final_xyz_sha256: hash of TS final XYZ bytes
        fragment_hash: SHA256 hash of JSON with fragmenter metadata

    Returns:
        Full sha256 hex string
    """
    payload = {
        "orca_template_hash": orca_template_hash,
        "geometry_hashes": geometry_hashes,
        "fragment_hash": fragment_hash
    }

    payload_str = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(payload_str.encode('utf-8')).hexdigest()
