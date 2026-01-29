"""
Step 4: Path Accessor (V4.2 Phase A)
======================================

Path validation and access utilities for input files.

Author: QC Descriptors Team
Date: 2026-01-18
"""

import hashlib
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any
import numpy as np

from rph_core.utils.file_io import read_xyz

logger = logging.getLogger(__name__)


def compute_file_fingerprint(file_path: Path) -> Optional[Dict[str, Any]]:
    """Compute file fingerprint including hash, size, and mtime.

    Args:
        file_path: Path to file

    Returns:
        Dictionary with keys:
            - hash: SHA1 hex digest
            - size: File size in bytes
            - mtime: Modification time (timestamp)
            or None if file cannot be read
    """
    file_path = Path(file_path)
    if not file_path.exists():
        return None

    try:
        stat = file_path.stat()
        size = stat.st_size
        mtime = stat.st_mtime

        with open(file_path, "rb") as f:
            content = f.read()

        file_hash = hashlib.sha1(content).hexdigest()

        return {
            "hash": file_hash,
            "size": size,
            "mtime": mtime,
        }
    except Exception as e:
        logger.warning(f"Failed to compute fingerprint for {file_path}: {e}")
        return None


def compute_sp_report_fingerprint(sp_report: Any) -> Optional[str]:
    """Compute fingerprint of SPMatrixReport key fields.

    Args:
        sp_report: SPMatrixReport instance

    Returns:
        SHA1 hex digest of key fields JSON or None
    """
    if sp_report is None:
        return None

    try:
        import json

        # Extract key fields for fingerprinting
        key_fields = {
            "method": getattr(sp_report, "method", ""),
            "solvent": getattr(sp_report, "solvent", ""),
            "e_ts": getattr(sp_report, "e_ts", 0.0),
            "e_reactant": getattr(sp_report, "e_reactant", 0.0),
            "e_product": getattr(sp_report, "e_product", 0.0),
            "g_ts": getattr(sp_report, "g_ts", None),
            "g_reactant": getattr(sp_report, "g_reactant", None),
            "g_product": getattr(sp_report, "g_product", None),
        }

        # Serialize to JSON with sorted keys
        json_str = json.dumps(key_fields, sort_keys=True, default=str)

        # Compute SHA1 hash
        return hashlib.sha1(json_str.encode("utf-8")).hexdigest()
    except Exception as e:
        logger.warning(f"Failed to compute sp_report fingerprint: {e}")
        return None


class PathAccessor:
    """Path validation and access utility for input files.

    Validates that required files exist and optionally computes fingerprints
    for reproducibility tracking.
    """

    def __init__(self, required_paths: Dict[str, Path], compute_fingerprints: bool = True):
        """Initialize PathAccessor.

        Args:
            required_paths: Dictionary mapping canonical key -> Path object
            compute_fingerprints: Whether to compute file fingerprints
        """
        self.required_paths = required_paths
        self.compute_fingerprints = compute_fingerprints
        self.fingerprints: Dict[str, Optional[Dict[str, Any]]] = {}
        self.missing_paths: List[str] = []

        if compute_fingerprints:
            self._compute_all_fingerprints()

    def _compute_all_fingerprints(self) -> None:
        for key, path_obj in self.required_paths.items():
            if path_obj is not None and path_obj.exists():
                fp = compute_file_fingerprint(path_obj)
                self.fingerprints[key] = fp
            else:
                self.fingerprints[key] = None
                if path_obj is not None:
                    self.missing_paths.append(key)

    def validate(self) -> bool:
        """Validate that all required paths exist.

        Returns:
            True if all paths exist, False otherwise
        """
        self.missing_paths = []

        for key, path_obj in self.required_paths.items():
            if path_obj is None or not path_obj.exists():
                self.missing_paths.append(key)

        return len(self.missing_paths) == 0

    def get_path(self, key: str) -> Optional[Path]:
        """Get a path by canonical key.

        Args:
            key: Canonical key name (e.g., "ts_xyz")

        Returns:
            Path object or None if not found/missing
        """
        return self.required_paths.get(key, None)

    def get_fingerprint(self, key: str) -> Optional[Dict[str, Any]]:
        """Get fingerprint for a path.

        Args:
            key: Canonical key name

        Returns:
            Dictionary with keys {hash, size, mtime} or None if not computed
        """
        return self.fingerprints.get(key, None)

    def get_fingerprint_hash(self, key: str) -> Optional[str]:
        """Get hash-only fingerprint for backward compatibility.

        Args:
            key: Canonical key name

        Returns:
            SHA1 hex digest or None if not computed
        """
        fp = self.fingerprints.get(key)
        return fp["hash"] if fp else None

    def read_xyz_file(self, key: str) -> Optional[tuple]:
        """Read an XYZ file using file_io.read_xyz.

        Args:
            key: Canonical key name for the path

        Returns:
            Tuple (symbols, coordinates) where:
                - symbols: List[str] - element symbols
                - coordinates: np.ndarray - (N, 3) coordinate array
            Returns None if read fails
        """
        path_obj = self.get_path(key)
        if path_obj is None or not path_obj.exists():
            logger.warning(f"Cannot read XYZ for missing key: {key}")
            return None

        try:
            # read_xyz returns (coordinates, symbols) - 2-tuple
            coords, symbols = read_xyz(path_obj)
            # Return (symbols, coords) so symbols[0] gives list for natoms calculation
            return (symbols, np.array(coords))
        except Exception as e:
            logger.error(f"Failed to read XYZ file {path_obj}: {e}")
            return None

    def get_missing_fields(self, required_fields: List[str]) -> List[str]:
        """Get list of missing path keys.

        Args:
            required_fields: List of canonical field names

        Returns:
            List of missing field names
        """
        missing = []

        for field in required_fields:
            if field not in self.required_paths:
                missing.append(field)
            elif self.required_paths[field] is None or not self.required_paths[field].exists():
                missing.append(field)

        return missing

    def __repr__(self) -> str:
        """String representation of PathAccessor."""
        return f"PathAccessor(paths={list(self.required_paths.keys())}, missing={self.missing_paths})"
