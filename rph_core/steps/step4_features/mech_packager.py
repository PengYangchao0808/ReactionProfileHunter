"""
Mechanism Asset Packager (M1)
==============================

Packages mechanism assets from S1/S2/S3 into S4 root directory with fixed naming.
Creates mech_index.json as single source of truth for downstream consumption.

Author: QC Descriptors Team
Date: 2026-01-21
"""

import json
import logging
import os
import shutil
import tempfile
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np

from rph_core.utils.file_io import read_xyz, write_xyz
from rph_core.utils.qc_interface import is_path_toxic
from rph_core.utils.naming_compat import (
    normalize_source_label,
    get_intermediate_source_priority,
    SOURCE_S2_INTERMEDIATE,
    SOURCE_S3_INTERMEDIATE,
)
from rph_core.utils.forming_bonds_resolver import resolve_forming_bonds

logger = logging.getLogger(__name__)


# ============================================================================
# Constants
# ============================================================================

S3_DIR_ALIASES = ["S3_TS"]
S2_DIR_ALIASES = ["S2_Retro"]
S1_DIR_ALIASES = ["S1_ConfGeneration"]

# ============================================================================
# Dataclasses
# ============================================================================

@dataclass
class MechanismContext:
    """
    M4-C: Unified mechanism context with resolved asset paths.

    Contains all resolved paths for S1/S2/S3 assets with their source labels.
    Used by pack_mechanism_assets() for centralized asset resolution.

    All Optional[Path] fields default to None for safe construction.
    """
    s1_dir: Optional[Path] = None
    s2_dir: Optional[Path] = None
    s3_dir: Optional[Path] = None

    s1_product: Optional[Path] = None
    s1_precursor: Optional[Path] = None
    s2_ts_guess: Optional[Path] = None
    s2_intermediate: Optional[Path] = None

    s3_ts_final: Optional[Path] = None
    s3_reactant_sp: Optional[Path] = None

    s1_precursor_source: str = "none"


@dataclass
class AssetInfo:
    """Information about a mechanism asset."""
    filename: str
    source_path: Optional[Path]
    source_step: str
    source_label: Optional[str] = None
    sha256: Optional[str] = None


@dataclass
class QualityFlags:
    """
    M2-D: Three-state quality flags for mechanism assets.

    Uses null for "not applicable/unknown" instead of false for "missing data".
    Avoids "missing input → false" misinterpretation in ML filtering.
    """
    atom_count_ok: Optional[bool] = None
    forming_bond_window_ok: Optional[bool] = None
    suspect_optimized_to_product: Optional[str] = None
    # M4-D-3: QC artifact quality flags
    ts_imag_freq_ok: Optional[bool] = None
    asset_hash_ok: Optional[bool] = None


@dataclass
class MechanismMetaStep2:
    """Metadata for forming bonds derived from optimized geometries."""
    forming_bonds: Optional[Tuple[Tuple[int, int], ...]] = None
    source_summary: str = ""


@dataclass
class MechanismMetaStep1:
    """Metadata for Step 1 mechanism assets."""
    ts_computed: bool = False
    degradation_reason: str = ""


@dataclass
class MechanismIndex:
    """
    M2-D: Complete mechanism index (mech_index.json structure).

    Contains schema versioning, timestamp, status tracking, quality flags,
    asset inventory, and metadata for downstream ML training and review.
    """
    version: str = "1.0.0"
    schema_version: str = "mech_index_v1"
    generated_at: str = ""
    status: str = "COMPLETE"
    assets: Dict[str, Optional[AssetInfo]] = field(default_factory=dict)
    quality_flags: QualityFlags = field(default_factory=QualityFlags)
    config: Dict[str, Any] = field(default_factory=dict)
    # M4-D-2: Optional QC artifacts field (NMR, Hirshfeld, NBO outputs)
    qc_artifacts: Dict[str, Any] = field(default_factory=dict)


# ============================================================================
# Helper Functions
# ============================================================================

# M4-E-1: REASON constants for is_mech_index_up_to_date()
class UpdateReason:
    """Constant reasons for mech_index.json update decisions."""
    OK = "OK"
    MISSING_FILE = "mech_index.json does not exist"
    LOAD_FAILED = "Failed to load mech_index.json"
    SCHEMA_MISMATCH = "Schema version mismatch"
    ASSET_MISSING = "Asset file missing"
    # P0-3: Migration-related reasons
    MIGRATION_FAILED = "Schema migration failed"
    WRITE_FAILED = "Failed to write migrated mech_index.json"
    MIGRATED = "Migrated to new schema version"


def resolve_mechanism_context(
    s4_dir: Path,
    config: Dict[str, Any],
    max_recursion_depth: int = 3
) -> MechanismContext:
    """
    M4-C-1: Unified S1/S2/S3 asset path resolution with alias support.

    Resolves all S1/S2/S3 asset paths using directory alias tables.
    Provides single source of truth for mechanism packaging.

    P1 FIX: Updated docstring - currently uses alias table lookup only
    (max_recursion_depth parameter is reserved for future recursive search).

    Args:
        s4_dir: S4_Data directory (parent of S1, S2, S3 subdirs)
        config: Configuration dictionary
        max_recursion_depth: Reserved for future recursive search (currently unused)

    Returns:
        MechanismContext with all resolved paths and source labels

    Raises:
        RuntimeError: If work_dir cannot be determined from s4_dir
    """
    work_dir = _find_work_dir(s4_dir, max_recursion_depth)
    if work_dir is None:
        return MechanismContext(
            s1_dir=None,
            s2_dir=None,
            s3_dir=None,
            s1_product=None,
            s1_precursor=None,
            s1_precursor_source="none",
            s2_ts_guess=None,
            s2_intermediate=None,
            s3_ts_final=None,
            s3_reactant_sp=None,
        )

    s1_dir = _find_step_dir(work_dir, S1_DIR_ALIASES)
    s2_dir = _find_step_dir(work_dir, S2_DIR_ALIASES)

    s1_product = _resolve_s1_product(s1_dir)
    s2_ts_guess, s2_reactant = _resolve_s2_assets(s2_dir)

    s3_dir = _find_s3_dir(work_dir, max_recursion_depth)
    s3_ts_final, s3_intermediate = _resolve_s3_assets(s3_dir)

    precursor_priority = config.get('precursor_source_priority', None)
    if precursor_priority is not None:
        s1_precursor, s1_precursor_source = _resolve_s1_precursor(
            s1_dir,
            s2_dir,
            precursor_priority
        )
    else:
        s1_precursor = None
        s1_precursor_source = "none"

    return MechanismContext(
        s1_dir=s1_dir,
        s2_dir=s2_dir,
        s3_dir=s3_dir,

        s1_product=s1_product,
        s1_precursor=s1_precursor,
        s1_precursor_source=s1_precursor_source,

        s2_ts_guess=s2_ts_guess,
        s2_intermediate=s2_reactant,

        s3_ts_final=s3_ts_final,
        s3_reactant_sp=s3_intermediate,
    )


def _find_s3_dir(work_dir: Path, max_recursion_depth: int) -> Optional[Path]:
    """
    Find S3 directory using alias table.

    P1 FIX: Updated docstring - uses alias table lookup only
    (max_recursion_depth parameter is reserved for future recursive search).

    Args:
        work_dir: Base directory containing S3 subdirectories
        max_recursion_depth: Reserved for future recursive search (currently unused)

    Returns:
        Resolved S3 directory path or None if not found
    """
    s3_path = work_dir / "S3_TS"
    if s3_path.exists():
        return s3_path
    logger.warning("S3 directory not found: S3_TS")
    return None


def _find_step_dir(work_dir: Path, aliases: List[str]) -> Optional[Path]:
    if not aliases:
        return None
    candidate = work_dir / aliases[0]
    if candidate.exists():
        return candidate
    return None


def _has_step_dirs(work_dir: Path) -> bool:
    for alias in S1_DIR_ALIASES + S2_DIR_ALIASES + S3_DIR_ALIASES:
        if (work_dir / alias).exists():
            return True
    return False


def _find_work_dir(s4_dir: Path, max_recursion_depth: int) -> Optional[Path]:
    base_dir = s4_dir.parent
    if base_dir.exists() and _has_step_dirs(base_dir):
        return base_dir
    if max_recursion_depth <= 1:
        return None
    max_depth = max_recursion_depth - 1
    queue = [(base_dir, 0)]
    while queue:
        current, depth = queue.pop(0)
        if depth >= max_depth:
            continue
        children = [p for p in current.iterdir() if p.is_dir()]
        for child in sorted(children):
            if _has_step_dirs(child):
                return child
            queue.append((child, depth + 1))
    return None


def _first_existing(candidates: List[Path]) -> Optional[Path]:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None
def _resolve_s1_product(s1_dir: Optional[Path]) -> Optional[Path]:
    """Resolve S1 product file.

    Args:
        s1_dir: S1_ConfGeneration directory path

    Returns:
        Path to product file or None
    """
    if s1_dir is None or not s1_dir.exists():
        return None

    # Search in order of preference
    candidates = [
        s1_dir / "product_min.xyz",
    ]

    return _first_existing(candidates)


def _copy_or_link_asset(src: Path, target: Path, mode: str = "copy") -> bool:
    """
    Copy or create a hard link/symlink for an asset in a safe, idempotent way.

    M1-P0: Accept both "symlink" and "link" as valid values for mode.
    Falls back to copy if linking fails (e.g., cross-filesystem, permission).
    M1-P0: Force copy mode if source or target path is toxic (contains spaces/brackets).

    Args:
        src: Source file path
        target: Target file path
        mode: "copy", "symlink", or "link" (symlink and link are treated identically)

    Returns:
        True on success, False otherwise.
    """
    try:
        target.parent.mkdir(parents=True, exist_ok=True)

        # M1-P0: Check for toxic paths - force copy mode if detected
        src_toxic = is_path_toxic(src)
        dst_toxic = is_path_toxic(target)
        if src_toxic or dst_toxic:
            if src_toxic:
                logger.debug(f"Toxic source path detected, forcing copy: {src}")
            if dst_toxic:
                logger.debug(f"Toxic target path detected, forcing copy: {target}")
            mode = "copy"

        # M1-P0: Accept both "symlink" and "link" as linking modes
        if mode in ("link", "symlink"):
            # prefer hard link when possible, fall back to copy
            try:
                if target.exists():
                    target.unlink()
                os.link(src, target)
            except Exception:
                logger.debug(f"Hard link failed, falling back to copy: {src} -> {target}")
                shutil.copy2(src, target)
        else:
            shutil.copy2(src, target)
        return True
    except Exception as e:
        logger.warning(f"Failed to copy/link asset {src} -> {target}: {e}")
        return False


def _compute_sha256(path: Path) -> Optional[str]:
    try:
        h = sha256()
        with open(path, 'rb') as f:
            for chunk in iter(lambda: f.read(8192), b''):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None


def _write_json_atomic(data: Dict[str, Any], target_path: Path) -> None:
    """
    Write JSON data atomically using tmp file + rename.

    Args:
        data: Dictionary to serialize as JSON
        target_path: Final target path (will be atomically replaced)
    """
    try:
        target_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = target_path.with_suffix(target_path.suffix + ".tmp")
        tmp_path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
        # Atomic rename on POSIX; on Windows this may not be truly atomic
        tmp_path.replace(target_path)
        logger.debug(f"Atomically wrote: {target_path}")
    except Exception as e:
        logger.error(f"Failed to write {target_path}: {e}")
        raise


def _resolve_dipole_source(s3_intermediate: Optional[Path], s2_reactant: Optional[Path], priority: List[str]) -> Tuple[Optional[Path], str]:
    """
    Resolve intermediate source based on configurable priority order.

    P0-2 FIX: Now actually respects the priority list instead of hardcoding S3 > S2.

    Args:
        s3_reactant: Path to S3 reactant_sp.xyz (or None)
        s2_reactant: Path to S2 intermediate.xyz (or None)
        priority: List of source labels in preference order (e.g., ['S3_intermediate', 'S2_intermediate'])

    Returns:
        Tuple of (resolved_path, source_label). Returns (None, "none") if no valid source found.
    """
    # Normalize priority labels to new standard
    normalized_priority = [normalize_source_label(l) for l in priority]
    
    # Build candidate map for O(1) lookup
    candidates = {
        SOURCE_S3_INTERMEDIATE: s3_intermediate,
        SOURCE_S2_INTERMEDIATE: s2_reactant,
    }
    
    # Add legacy aliases for backward compatibility
    candidates['S3_reactant'] = s3_intermediate
    candidates['S2_reactant_complex'] = s2_reactant

    # Iterate through priority list and return first valid .xyz path
    for label in normalized_priority:
        path = candidates.get(label)
        if path and path.exists() and path.suffix == ".xyz":
            return path, label

    # No valid source found
    return None, "none"


def _resolve_s2_assets(s2_dir: Optional[Path]) -> Tuple[Optional[Path], Optional[Path]]:
    """Resolve S2 assets (ts_guess and intermediate)."""
    if s2_dir is None or not s2_dir.exists():
        return None, None

    ts_guess = _first_existing([s2_dir / "ts_guess.xyz"])
    intermediate = _first_existing([s2_dir / "intermediate.xyz"])

    return ts_guess, intermediate


def _resolve_s3_assets(s3_dir: Optional[Path]) -> Tuple[Optional[Path], Optional[Path]]:
    """Resolve S3 assets (ts_final and reactant SP).

    Args:
        s3_dir: S3_TransitionAnalysis directory path

    Returns:
        Tuple of (ts_final_path, reactant_sp_path)
    """
    if s3_dir is None or not s3_dir.exists():
        return None, None

    ts_final = _first_existing([s3_dir / "ts_final.xyz"])

    reactant_candidates = [
        s3_dir / "reactant_sp.xyz",
        s3_dir / "reactant_sp.out",
        s3_dir / "reactant_sp.log"
    ]

    reactant_opt_dir = s3_dir / "reactant_opt"
    if reactant_opt_dir.exists():
        reactant_candidates.extend([
            reactant_opt_dir / "final_output.log",
            reactant_opt_dir / "final_output.xyz",
            reactant_opt_dir / "reactant_sp.xyz"
        ])

    xyz_candidates = [c for c in reactant_candidates if c.suffix == ".xyz"]
    reactant_sp = _first_existing(xyz_candidates)

    return ts_final, reactant_sp


def _resolve_s1_precursor(
    s1_dir: Optional[Path],
    s2_dir: Optional[Path],
    priority: Optional[List[str]] = None
) -> Tuple[Optional[Path], str]:
    """
    M2-B2: Resolve S1 precursor state with enhanced fallback chain.

    Priority: S1 precursor → S2 neutral_precursor → S2 intermediate

    Priority order controlled by config.precursor_source_priority:
    - Default: ["S1_precursor", "S2_neutral_precursor", "S2_intermediate"]
    - Respects Step2 contract: intermediate.xyz is canonical fallback

    Args:
        s1_dir: S1_ConfGeneration directory path (MUST be a directory, not a file)
        s2_dir: S2_Retro directory path (MUST be a directory, not a file)
        priority: Precursor source priority list from config

    Returns:
        Tuple of (precursor_path, source_label)

    Raises:
        TypeError: If s1_dir or s2_dir is not None and not a directory
    """
    # M1-P0: Validate that we receive directories, not file paths
    for param_name, param_value in [('s1_dir', s1_dir), ('s2_dir', s2_dir)]:
        if param_value is not None:
            if not isinstance(param_value, Path):
                raise TypeError(f"{param_name} must be a Path, got {type(param_value).__name__}")
            if param_value.exists() and not param_value.is_dir():
                raise TypeError(f"{param_name} must be a directory, got file path: {param_value}")

    # M2-B2: Use priority list from config
    if priority is None:
        priority = ["S1_precursor", "S2_neutral_precursor", "S2_intermediate"]

    # M2-B2: Iterate through priority list and return first valid source
    for label in priority:
        if label == "S1_precursor":
            if s1_dir is not None and s1_dir.exists():
                s1_candidates = [
                    s1_dir / "precursor.xyz",
                    s1_dir / "neutral_precursor.xyz"
                ]
                s1_candidate = _first_existing(s1_candidates)
                if s1_candidate is not None:
                    return s1_candidate, "S1_precursor"
        
        elif label == "S2_neutral_precursor":
            if s2_dir is not None and s2_dir.exists():
                s2_precursor = _first_existing([s2_dir / "neutral_precursor.xyz"])
                if s2_precursor is not None:
                    return s2_precursor, "S2_neutral_precursor"
        
        elif label == "S2_intermediate":
            if s2_dir is not None and s2_dir.exists():
                s2_reactant = _first_existing([s2_dir / "intermediate.xyz"])
                if s2_reactant is not None:
                    return s2_reactant, "S2_intermediate"

    return None, "none"


def _check_forming_bond_window(
    file_path: Path,
    forming_bonds: Optional[Tuple[Tuple[int, int], ...]]
) -> bool:
    """Check if forming bonds are within reasonable distance window.

    Args:
        file_path: Path to structure file (TS or dipole)
        forming_bonds: Forming bond atom indices

    Returns:
        True if bonds are within reasonable window, False otherwise
    """
    if forming_bonds is None:
        # No forming bonds data, assume OK
        return True

    if not file_path.exists():
        logger.warning(f"File not found for bond window check: {file_path}")
        return False

    try:
        coords, symbols = read_xyz(file_path)

        # Check each forming bond
        for bond in forming_bonds:
            idx1, idx2 = bond

            # Validate indices
            if idx1 >= len(coords) or idx2 >= len(coords):
                logger.warning(
                    f"Forming bond indices out of range: {bond} "
                    f"(natoms={len(coords)})"
                )
                return False

            # Calculate distance
            coord1 = coords[idx1]
            coord2 = coords[idx2]
            distance = np.linalg.norm(coord1 - coord2)

            # Reasonable window: 1.5Å to 3.5Å
            # (TS bonds are typically stretched)
            if not (1.5 <= distance <= 3.5):
                logger.warning(
                    f"Forming bond distance outside window: {bond} "
                    f"distance={distance:.2f}Å"
                )
                return False

        return True
    except Exception as e:
        logger.error(f"Error checking forming bond window: {e}")
        return False


# M4-P0: Forming bond distance window constants
_FORMING_BOND_TS_WINDOW = (1.5, 3.5)  # TS bonds are typically stretched
_FORMING_BOND_PRODUCT_THRESHOLD = 2.5  # Typical product bond length


def _compute_bond_distances(
    file_path: Path,
    forming_bonds: Optional[Tuple[Tuple[int, int], ...]]
) -> Optional[Dict[Tuple[int, int], float]]:
    """
    Compute forming bond distances from an XYZ file.

    M4-P0: Helper for suspect_optimized_to_product check.

    Args:
        file_path: Path to structure file
        forming_bonds: Forming bond atom indices

    Returns:
        Dict mapping bond tuples to distances, or None on error
    """
    if forming_bonds is None or not file_path.exists():
        return None

    try:
        coords, _ = read_xyz(file_path)
        distances: Dict[Tuple[int, int], float] = {}

        for bond in forming_bonds:
            idx1, idx2 = bond

            if idx1 >= len(coords) or idx2 >= len(coords):
                continue

            coord1 = coords[idx1]
            coord2 = coords[idx2]
            distances[bond] = float(np.linalg.norm(coord1 - coord2))

        return distances if distances else None
    except Exception as e:
        logger.debug(f"Could not compute bond distances from {file_path}: {e}")
        return None


def _check_suspect_optimized_to_product(
    ts_path: Optional[Path],
    product_path: Optional[Path],
    forming_bonds: Optional[Tuple[Tuple[int, int], ...]]
) -> Optional[str]:
    """
    M4-P0: Check if optimization may have converged to product instead of TS.

    Logic:
    1. TS bond distance should be in [1.5, 3.5] window (stretched TS geometry)
    2. Product bond distance should be < 2.5 (typical product bond)
    3. If TS is in window but product distance is significantly smaller,
       suspect_optimized_to_product = "suspect"
    4. If any required data is missing, return None (unknown)

    Args:
        ts_path: Path to TS structure (optimized or guess)
        product_path: Path to product structure (optional)
        forming_bonds: Forming bond atom indices

    Returns:
        "suspect" if optimization likely converged to product,
        "ok" if distances look reasonable,
        None if insufficient data to判断
    """
    if not forming_bonds:
        return None

    if not ts_path or not ts_path.exists():
        return None

    ts_distances = _compute_bond_distances(ts_path, forming_bonds)
    if not ts_distances:
        return None

    # Check if TS distances are in expected window
    ts_in_window = all(
        _FORMING_BOND_TS_WINDOW[0] <= d <= _FORMING_BOND_TS_WINDOW[1]
        for d in ts_distances.values()
    )

    if not ts_in_window:
        # TS not in expected window - cannot判断
        logger.debug(f"TS distances not in window, cannot判断 suspect_optimized_to_product")
        return None

    if not product_path:
        # No product to compare against
        logger.debug(f"No product file, cannot判断 suspect_optimized_to_product")
        return None

    product_distances = _compute_bond_distances(product_path, forming_bonds)
    if not product_distances:
        return None

    # Compare: if product bonds are significantly shorter than TS,
    # this suggests optimization went in the right direction
    for bond in forming_bonds:
        if bond in ts_distances and bond in product_distances:
            ts_d = ts_distances[bond]
            prod_d = product_distances[bond]

            # If product distance is much smaller than TS, that's normal
            # If product distance is close to or larger than TS, that's suspicious
            if prod_d >= ts_d * 0.9:  # Product is >= 90% of TS distance
                logger.warning(
                    f"  ⚠ suspect_optimized_to_product: bond {bond} "
                    f"TS={ts_d:.2f}Å, product={prod_d:.2f}Å "
                    f"(product not significantly shorter)"
                )
                return "suspect"

    logger.debug(f"  ✓ suspect_optimized_to_product: distances look reasonable")
    return "ok"


def _validate_atom_count(xyz_path: Path, expected_count: int) -> bool:
    """
    Validate that an XYZ file contains the expected number of atoms.

    Args:
        xyz_path: Path to XYZ file
        expected_count: Expected number of atoms

    Returns:
        True if count matches, False otherwise
    """
    if not xyz_path.exists():
        logger.warning(f"XYZ file not found: {xyz_path}")
        return False

    try:
        coords, symbols = read_xyz(xyz_path)
        actual_count = len(coords)
        if actual_count == expected_count:
            return True
        else:
            logger.warning(
                f"Atom count mismatch: expected {expected_count}, got {actual_count}"
            )
            return False
    except Exception as e:
        logger.error(f"Error reading XYZ file {xyz_path}: {e}")
        return False


def _build_mech_index(
    assets: Dict[str, Optional[AssetInfo]],
    quality_flags: QualityFlags,
    step2_meta: MechanismMetaStep2,
    step1_meta: MechanismMetaStep1,
    config: Dict[str, Any],
    missing_inputs: Optional[List[str]] = None,
    degradation_reasons: Optional[List[str]] = None,
    qc_artifacts: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    M2-D: Build mech_index.json structure with schema versioning.

    Args:
        assets: Dictionary of asset info (key -> AssetInfo)
        quality_flags: Quality flags
        step2_meta: Step 2 metadata
        step1_meta: Step 1 metadata
        config: Packager configuration
        missing_inputs: List of missing inputs (for INCOMPLETE status)
        degradation_reasons: List of degradation reasons (for INCOMPLETE status)
        qc_artifacts: QC artifact paths (NMR, Hirshfeld, NBO) - M4-D-4

    Returns:
        Dictionary ready for JSON serialization
    """
    # M2-D: Build assets dict, excluding None values
    assets_dict = {}
    for key, asset_info in assets.items():
        if asset_info is None:
            assets_dict[key] = None
        else:
            assets_dict[key] = {
                "filename": asset_info.filename,
                "source_path": str(asset_info.source_path) if asset_info.source_path else None,
                "source_step": asset_info.source_step,
                "sha256": asset_info.sha256
            }
            if asset_info.source_label:
                assets_dict[key]["source_label"] = asset_info.source_label

    # Build step2 metadata
    step2_meta_dict: Dict[str, Any] = {
        "filename": "mech_step2_meta.json",
        "forming_bonds": None,
        "source_summary": step2_meta.source_summary
    }
    if step2_meta.forming_bonds:
        step2_meta_dict["forming_bonds"] = [list(b) for b in step2_meta.forming_bonds]

    # Build step1 metadata
    step1_meta_dict = {
        "filename": "mech_step1_meta.json",
        "ts_computed": step1_meta.ts_computed,
        "degradation_reason": step1_meta.degradation_reason
    }

    # M2-D2: Track missing inputs
    missing_inputs_list = []
    if missing_inputs:
        missing_inputs_list = missing_inputs
    else:
        # Auto-detect missing inputs
        if not assets.get('mech_step2_ts2'):
            missing_inputs_list.append("mech_step2_ts2.xyz")
        if not assets.get('mech_step2_reactant_dipole'):
            missing_inputs_list.append("mech_step2_reactant_dipole.xyz")
        if not assets.get('mech_step2_product'):
            missing_inputs_list.append("mech_step2_product.xyz")
        if not assets.get('mech_step1_precursor'):
            missing_inputs_list.append("mech_step1_precursor.xyz")

    # M2-D2: Track degradation reasons
    degradation_reasons_list = degradation_reasons if degradation_reasons else []
    if not assets.get('mech_step2_ts2'):
        degradation_reasons_list.append("S3 TS optimization not completed")
    if quality_flags.atom_count_ok is False:
        degradation_reasons_list.append("Atom count mismatch between dipole and TS2")
    if quality_flags.forming_bond_window_ok is False:
        degradation_reasons_list.append("Forming bond distances outside expected window")

    # Determine overall mechanism_status
    any_missing = any(asset is None for asset in assets.values())
    mechanism_status = "INCOMPLETE" if any_missing else "COMPLETE"

    # M2-D3: Get schema_version from config
    schema_version = config.get('schema_version', 'mech_index_v1')

    # Build complete index
    mech_index = {
        "version": "1.0.0",
        "schema_version": schema_version,  # M2-D3
        "generated_at": datetime.now(timezone.utc).isoformat(),  # M2-D3
        "mechanism_status": mechanism_status,  # M2-D3
        "assets": {
            **assets_dict,
            "mech_step2_meta": step2_meta_dict,
            "mech_step1_meta": step1_meta_dict
        },
        "quality_flags": {
            "atom_count_ok": quality_flags.atom_count_ok,
            "forming_bond_window_ok": quality_flags.forming_bond_window_ok,
            "suspect_optimized_to_product": quality_flags.suspect_optimized_to_product,
            # M4-D-3: QC artifact quality flags (default None, calculated in pack_mechanism_assets)
            "ts_imag_freq_ok": getattr(quality_flags, 'ts_imag_freq_ok', None),
            "asset_hash_ok": getattr(quality_flags, 'asset_hash_ok', None)
        },
        "missing_inputs": missing_inputs_list,  # M2-D2
        "degradation_reasons": degradation_reasons_list,  # M2-D2
        "config": config,
        "qc_artifacts": qc_artifacts or {}  # M4-D-4
    }

    return mech_index


def is_mech_index_up_to_date(
    s4_dir: Path,
    expected_schema: str = "mech_index_v1",
    validate_assets: bool = False
) -> Tuple[bool, str]:
    """
    M3-1: Check if mech_index.json is up-to-date with expected schema.

    P0-3 FIX: This function is now explicitly a PURE CHECK with NO side effects.
    It does NOT write back migrated results to disk.

    For migration + write-back, use ensure_mech_index_schema() instead.

    Args:
        s4_dir: S4_Data directory path
        expected_schema: Expected schema_version (default: mech_index_v1)
        validate_assets: If True, also check asset existence (optional, default False)

    Returns:
        Tuple of (is_up_to_date, reason)
        is_up_to_date: True if all checks pass, False otherwise
        reason: String explaining why not up-to-date, or "OK" if up-to-date
    """
    mech_index_path = s4_dir / "mech_index.json"
    if not mech_index_path.exists():
        return False, UpdateReason.MISSING_FILE

    try:
        with open(mech_index_path, 'r', encoding='utf-8') as f:
            mech_index = json.load(f)
    except Exception as e:
        return False, f"{UpdateReason.LOAD_FAILED}: {e}"

    # Check schema version (NO migration, just check)
    actual_schema = mech_index.get('schema_version')
    if actual_schema != expected_schema:
        return False, f"{UpdateReason.SCHEMA_MISMATCH}: expected {expected_schema}, got {actual_schema}"

    if validate_assets:
        assets = mech_index.get('assets', {})
        for asset_key, asset_info in assets.items():
            if asset_info is None:
                continue

            filename = asset_info.get('filename')
            if filename:
                asset_path = s4_dir / filename
                if not asset_path.exists():
                    return False, f"{UpdateReason.ASSET_MISSING}: {filename}"

    return True, UpdateReason.OK


def ensure_mech_index_schema(
    s4_dir: Path,
    expected_schema: str = "mech_index_v1"
) -> Tuple[bool, str, Optional[Dict[str, Any]]]:
    """
    P0-3: Ensure mech_index.json has expected schema, with optional write-back.

    This function provides migration + write-back semantics that is_mech_index_up_to_date() lacks.
    It also ensures canonical fields (like mechanism_status) are present.

    Args:
        s4_dir: S4_Data directory path
        expected_schema: Expected schema_version (default: mech_index_v1)

    Returns:
        Tuple of (success, reason, migrated_index)
        - success: True if already up-to-date or migration/fix succeeded
        - reason: UpdateReason constant or error message
        - migrated_index: The migrated/fixed dict if changes occurred, None otherwise
    """
    mech_index_path = s4_dir / "mech_index.json"

    if not mech_index_path.exists():
        return False, UpdateReason.MISSING_FILE, None

    try:
        with open(mech_index_path, 'r', encoding='utf-8') as f:
            mech_index = json.load(f)
    except Exception as e:
        return False, f"{UpdateReason.LOAD_FAILED}: {e}", None

    actual_schema = mech_index.get('schema_version')
    needs_write = False

    # Check schema version mismatch
    if actual_schema != expected_schema:
        logger.info(f"Migrating mech_index.json from '{actual_schema}' to '{expected_schema}'")
        try:
            mech_index = migrate_mech_index(mech_index)
        except Exception as e:
            return False, f"{UpdateReason.MIGRATION_FAILED}: {e}", None
        needs_write = True

    # Ensure canonical fields are present (even if schema version matches)
    # P0-3: Fill mechanism_status from status if missing
    if 'mechanism_status' not in mech_index and 'status' in mech_index:
        mech_index['mechanism_status'] = mech_index['status']
        needs_write = True

    # Write back if changes were made
    if needs_write:
        try:
            _write_json_atomic(mech_index, mech_index_path)
            logger.info(f"  ✓ mech_index.json updated and written to {mech_index_path}")
        except Exception as e:
            return False, f"{UpdateReason.WRITE_FAILED}: {e}", mech_index
        return True, UpdateReason.MIGRATED if actual_schema != expected_schema else UpdateReason.OK, mech_index

    return True, UpdateReason.OK, None


def migrate_mech_index(old_index: Dict[str, Any]) -> Dict[str, Any]:
    """
    M3-2-1: Migrate old mech_index.json to new schema (v1).

    Handles:
    - timestamp/status → generated_at/mechanism_status
    - Missing schema_version → add mech_index_v1
    - Optional deprecated alias fields for backward compatibility

    Args:
        old_index: Old mech_index.json dictionary

    Returns:
        Migrated mech_index.json dictionary with v1 schema
    """
    migrated = old_index.copy()

    if 'schema_version' not in migrated:
        migrated['schema_version'] = 'mech_index_v1'

    if 'timestamp' in old_index and 'generated_at' not in migrated:
        migrated['generated_at'] = old_index['timestamp']

    if 'status' in old_index and 'mechanism_status' not in migrated:
        migrated['mechanism_status'] = old_index['status']

    return migrated


# ============================================================================
# QC Artifact Collection Constants
# ============================================================================

QC_ARTIFACT_SUBDIRS = {
    "nbo": ["nbo_analysis/", "nbo/", "reactant_opt/standard/", "reactant_opt/rescue/"]
}

QC_ARTIFACT_TARGETS = {
    "nbo_outputs": "qc_nbo.37"
}

QC_ARTIFACT_PATTERNS = {
    "nbo_outputs": ["*.37", "*.nbo", "*.nbo7"]
}


def _collect_qc_artifacts(
    s3_dir: Path,
    pipeline_root: Path,
    out_dir: Path,
    copy_mode: str = "copy"
) -> Dict[str, Dict[str, Any]]:
    """
    M4-D-4: Collect QC artifacts with restricted search + whitelist + relative paths.

    Searches only in whitelisted subdirectories of s3_dir, copies artifacts to S4 root
    with fixed naming, and returns structure with relative filename + meta.source_paths.

    M4-P0 FIX: Now collects ALL candidates and picks by mtime (not first-found).
    Meta records: candidates list with metadata, picked file, and accurate reason.

    Args:
        s3_dir: S3 directory path (S3_TS or S3_TransitionAnalysis)
        pipeline_root: Pipeline root directory (for relative path calculation)
        out_dir: S4_Data output directory
        copy_mode: Copy mode ("copy" or "link")

    Returns:
        Dictionary mapping artifact type to structure:
        {
            "filename": "qc_nmrdat",
            "meta": {
                "candidates": [
                    {"rel_path": "S3_TS/nmr/output.nmrdat", "mtime": 1234567890.0, "size": 1024, "sha256": "..."}
                ],
                "picked": {"rel_path": "...", "mtime": 1234567890.0, "size": 1024},
                "reason": "picked_by_mtime"
            }
        }
        Empty result for types not found.
    """
    artifacts: Dict[str, Dict[str, Any]] = {}
    glob_cache: Dict[Tuple[str, str, str], List[Path]] = {}
    import hashlib

    def _compute_file_hash(file_path: Path, limit_bytes: int = 4096) -> str:
        """Compute SHA256 hash of file (first N bytes for performance)."""
        try:
            h = hashlib.sha256()
            with open(file_path, 'rb') as f:
                # Read only first chunk for performance (identification, not security)
                data = f.read(limit_bytes)
                h.update(data)
            return h.hexdigest()[:16]  # First 16 chars sufficient for identification
        except Exception:
            return "unknown"

    def _glob_cached(search_root: Path, pattern: str, recursive: bool = False) -> List[Path]:
        key = (str(search_root), pattern, "r" if recursive else "g")
        cached = glob_cache.get(key)
        if cached is not None:
            return cached
        results = list(search_root.rglob(pattern)) if recursive else list(search_root.glob(pattern))
        glob_cache[key] = results
        return results

    if not s3_dir or not s3_dir.exists():
        logger.debug(f"S3 directory not found, skipping QC artifact collection: {s3_dir}")
        return artifacts

    for artifact_type, target_filename in QC_ARTIFACT_TARGETS.items():
        patterns = QC_ARTIFACT_PATTERNS.get(artifact_type, [])
        if not patterns:
            continue

        # Collect ALL candidates first (M4-P0: fix contract fraud)
        candidates: List[Dict[str, Any]] = []
        source_subdirs = QC_ARTIFACT_SUBDIRS.get(artifact_type.split("_")[0], [])

        # Search in whitelisted subdirectories
        for subdir in source_subdirs:
            search_root = s3_dir / subdir
            if not search_root.exists():
                continue

            for pattern in patterns:
                for candidate in _glob_cached(search_root, pattern, recursive=False):
                    if candidate.is_file():
                        try:
                            stat = candidate.stat()
                            candidates.append({
                                "rel_path": str(candidate.relative_to(pipeline_root)),
                                "abs_path": str(candidate),
                                "mtime": stat.st_mtime,
                                "size": stat.st_size,
                                "sha256": _compute_file_hash(candidate)
                            })
                        except (ValueError, OSError) as e:
                            logger.debug(f"Could not stat candidate {candidate}: {e}")

        # Fallback: search in S3 root with very restricted patterns (max depth 2)
        if not candidates:
            for pattern in patterns:
                for candidate in _glob_cached(s3_dir, pattern, recursive=True):
                    if candidate.is_file():
                        try:
                            rel_depth = len(candidate.relative_to(s3_dir).parts)
                            if rel_depth <= 2:
                                stat = candidate.stat()
                                candidates.append({
                                    "rel_path": str(candidate.relative_to(pipeline_root)),
                                    "abs_path": str(candidate),
                                    "mtime": stat.st_mtime,
                                    "size": stat.st_size,
                                    "sha256": _compute_file_hash(candidate)
                                })
                        except (ValueError, OSError) as e:
                            logger.debug(f"Could not stat candidate {candidate}: {e}")

        if not candidates:
            logger.debug(f"No {artifact_type} found in {s3_dir}")
            continue

        # M4-P0: Select by mtime (not first-found)
        # Sort by mtime descending, pick the newest
        candidates.sort(key=lambda x: x["mtime"], reverse=True)
        picked = candidates[0]

        # Minimal validation: size > 0 and file exists
        if picked["size"] == 0:
            logger.warning(f"  ⚠ Picked file has zero size: {picked['rel_path']}")
            # Check if there are other non-empty candidates
            non_empty = [c for c in candidates if c["size"] > 0]
            if non_empty:
                non_empty.sort(key=lambda x: x["mtime"], reverse=True)
                picked = non_empty[0]
                logger.info(f"  ↳ Falling back to non-zero file: {picked['rel_path']}")
            else:
                logger.warning(f"  ✗ All candidates for {artifact_type} are empty, skipping")
                continue

        # Copy to S4 root with fixed naming
        source_path = Path(picked["abs_path"])
        target = out_dir / target_filename

        try:
            if _copy_or_link_asset(source_path, target, copy_mode):
                # Record only top 5 candidates to avoid bloating meta
                meta_candidates = [
                    {"rel_path": c["rel_path"], "mtime": c["mtime"], "size": c["size"], "sha256": c["sha256"]}
                    for c in candidates[:5]
                ]

                artifacts[artifact_type] = {
                    "filename": target_filename,
                    "meta": {
                        "candidates": meta_candidates,
                        "picked": {
                            "rel_path": picked["rel_path"],
                            "mtime": picked["mtime"],
                            "size": picked["size"]
                        },
                        "reason": "picked_by_mtime"
                    }
                }
                logger.info(f"  ✓ Collected {artifact_type}: {target} (from {picked['rel_path']})")
            else:
                logger.warning(f"  ✗ Failed to copy {artifact_type}: {source_path} -> {target}")
        except Exception as e:
            logger.warning(f"  ✗ Error collecting {artifact_type}: {e}")

    collected_count = len(artifacts)
    logger.debug(f"QC artifacts collected: {collected_count}/{len(QC_ARTIFACT_TARGETS)}")

    return artifacts


# ============================================================================
# Main Function
# ============================================================================

def pack_mechanism_assets(
    step_dirs: Dict[str, Optional[Path]],
    out_dir: Path,
    config: Dict[str, Any],
    forming_bonds: Optional[Tuple[Tuple[int, int], ...]] = None,
    use_central_resolver: bool = True
) -> Dict[str, Any]:
    """
    Package mechanism assets from S1/S2/S3 into S4 root directory with fixed naming.

    Creates fixed-named asset files and mech_index.json in S4 root.

    M4-C: Updated to use centralized resolve_mechanism_context() for asset resolution.

    Args:
        step_dirs: Dictionary with keys 'S1', 'S2', 'S3' pointing to step directories
        out_dir: S4 output directory (S4_Data/)
        config: Configuration dictionary for mechanism packaging
        forming_bonds: Optional override for forming bonds
        use_central_resolver: Use centralized asset resolution (default: True)

    Returns:
        Dictionary containing mech_index.json content
    """
    logger.info("=" * 60)
    logger.info("Mechanism Asset Packager (M1) with Central Resolver")
    logger.info("=" * 60)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    context = resolve_mechanism_context(out_dir, config)

    resolved_forming_bonds = forming_bonds
    if resolved_forming_bonds is None:
        forming_cfg = config.get('forming_bonds', {}) or {}
        resolved = resolve_forming_bonds(
            product_xyz=context.s1_product,
            ts_xyz=context.s3_ts_final,
            s3_dir=context.s3_dir,
            s4_dir=out_dir,
            config=forming_cfg,
            write_meta=forming_cfg.get('write_meta', True)
        )
        resolved_forming_bonds = resolved.forming_bonds
        if resolved.warnings:
            for w in resolved.warnings:
                logger.warning(f"Forming bonds resolver: {w}")

    # Derive work_dir from step_dirs when available (prefer explicit over guessing)
    work_dir = None
    for key in ('S1', 'S2', 'S3'):
        p = step_dirs.get(key)
        if p:
            work_dir = Path(p)
            break
    if work_dir is None:
        work_dir = out_dir.parent

    assets: Dict[str, Optional[AssetInfo]] = {}

    # Resolve intermediate source
    intermediate_source, intermediate_label = _resolve_dipole_source(
        context.s3_reactant_sp,
        context.s2_intermediate,
        get_intermediate_source_priority(config)
    )
    logger.debug(f"Intermediate source: {intermediate_label} -> {intermediate_source}")

    # Resolve precursor
    # M1-P0: Pass directories (context.s1_dir, context.s2_dir), not file paths
    precursor, precursor_source = _resolve_s1_precursor(
        context.s1_dir,
        context.s2_dir,
        config.get('precursor_source_priority', None)
    )
    logger.debug(f"Precursor: {precursor_source} -> {precursor}")

    # M3-2-2: Migrate old mech_index.json if needed
    mech_index_path = out_dir / "mech_index.json"
    if mech_index_path.exists():
        try:
            with open(mech_index_path, 'r', encoding='utf-8') as f:
                existing_index = json.load(f)
            if existing_index.get('schema_version') != 'mech_index_v1':
                logger.info(f"Migrating old mech_index.json to v1 schema...")
                migrated_index = migrate_mech_index(existing_index)
                with open(mech_index_path, 'w', encoding='utf-8') as f:
                    json.dump(migrated_index, f, indent=2, default=str)
                logger.info(f"  ✓ Migrated mech_index.json to v1 schema")
        except Exception as e:
            logger.warning(f"Failed to migrate mech_index.json: {e}, will regenerate")

    # Extract configuration
    enabled = config.get('enabled', False)
    copy_mode = config.get('copy_mode', 'copy')

    if not enabled:
        logger.info("Mechanism packaging disabled (enabled=False)")
        return {}

    logger.info(f"Copy mode: {copy_mode}")
    logger.info(f"Intermediate source priority: {get_intermediate_source_priority(config)}")

    # Copy/link assets to S4 root with fixed naming
    assets: Dict[str, Optional[AssetInfo]] = {}

    # Resolve commonly used context variables (avoid relying on re-used local names)
    s3_ts_final = context.s3_ts_final
    s1_product = context.s1_product
    precursor = context.s1_precursor
    precursor_label = context.s1_precursor_source
    write_quality_flags = config.get('write_quality_flags', True)

    # Asset 1: mech_step2_ts2.xyz (from S3 ts_final)
    if s3_ts_final:
        target = out_dir / "mech_step2_ts2.xyz"
        if _copy_or_link_asset(s3_ts_final, target, copy_mode):
            assets['mech_step2_ts2'] = AssetInfo(
                filename="mech_step2_ts2.xyz",
                source_path=s3_ts_final,
                source_step="S3",
                sha256=_compute_sha256(target)
            )
            logger.info(f"  ✓ Created: {target}")
        else:
            logger.warning(f"  ✗ Failed to create: {target}")
            assets['mech_step2_ts2'] = None
    else:
        logger.warning("  ✗ S3 ts_final not available")
        assets['mech_step2_ts2'] = None

    # Asset 2: mech_step2_reactant_intermediate.xyz (intermediate source)
    if intermediate_source:
        target = out_dir / "mech_step2_reactant_intermediate.xyz"
        if _copy_or_link_asset(intermediate_source, target, copy_mode):
            assets['mech_step2_reactant_intermediate'] = AssetInfo(
                filename="mech_step2_reactant_intermediate.xyz",
                source_path=intermediate_source,
                source_step="S3" if "S3" in intermediate_label else "S2",
                source_label=intermediate_label,
                sha256=_compute_sha256(target)
            )
            logger.info(f"  ✓ Created: {target}")
        else:
            logger.warning(f"  ✗ Failed to create: {target}")
            assets['mech_step2_reactant_intermediate'] = None
    else:
        if intermediate_label == "none":
            logger.warning("  ✗ Intermediate source not available: no .xyz file found in S3_intermediate or S2_intermediate")
        else:
            logger.warning(f"  ✗ Intermediate source not available: {intermediate_label}")
        assets['mech_step2_reactant_intermediate'] = None

    # Asset 3: mech_step2_product.xyz (from S1 product)
    if s1_product:
        target = out_dir / "mech_step2_product.xyz"
        if _copy_or_link_asset(s1_product, target, copy_mode):
            assets['mech_step2_product'] = AssetInfo(
                filename="mech_step2_product.xyz",
                source_path=s1_product,
                source_step="S1",
                sha256=_compute_sha256(target)
            )
            logger.info(f"  ✓ Created: {target}")
        else:
            logger.warning(f"  ✗ Failed to create: {target}")
            assets['mech_step2_product'] = None
    else:
        logger.warning("  ✗ S1 product not available")
        assets['mech_step2_product'] = None

    # Asset 5: mech_step1_precursor.xyz
    if precursor:
        target = out_dir / "mech_step1_precursor.xyz"
        if _copy_or_link_asset(precursor, target, copy_mode):
            assets['mech_step1_precursor'] = AssetInfo(
                filename="mech_step1_precursor.xyz",
                source_path=precursor,
                source_step="S1" if "S1" in precursor_label else "S2",
                source_label=precursor_label,
                sha256=_compute_sha256(target)
            )
            logger.info(f"  ✓ Created: {target}")
        else:
            logger.warning(f"  ✗ Failed to create: {target}")
            assets['mech_step1_precursor'] = None
    else:
        logger.warning("  ✗ Precursor not available")
        assets['mech_step1_precursor'] = None

    # Perform quality checks
    quality_flags = QualityFlags()

    # M4-P0: Extract ts2_asset for both quality checks and suspect_optimized_to_product
    ts2_asset = assets.get('mech_step2_ts2')
    dipole_asset = assets.get('mech_step2_reactant_dipole')

    if write_quality_flags:
        logger.info("Performing lightweight quality checks...")

        # atom_count_ok: Compare dipole and ts2
        if ts2_asset and ts2_asset.source_path and dipole_asset and dipole_asset.source_path:
            ts2_natoms = len(read_xyz(ts2_asset.source_path)[1]) if ts2_asset.source_path.exists() else 0
            dipole_natoms = len(read_xyz(dipole_asset.source_path)[1]) if dipole_asset.source_path.exists() else 0

            if ts2_natoms > 0 and ts2_natoms == dipole_natoms:
                quality_flags.atom_count_ok = True
                logger.info(f"  ✓ atom_count_ok: {ts2_natoms} atoms")
            else:
                quality_flags.atom_count_ok = False
                logger.warning(f"  ✗ atom_count_ok: ts2={ts2_natoms}, dipole={dipole_natoms}")
        else:
            # M1-P0: Missing assets -> None (unknown/not applicable), not False
            quality_flags.atom_count_ok = None
            missing_assets = []
            if not (ts2_asset and ts2_asset.source_path):
                missing_assets.append("ts2")
            if not (dipole_asset and dipole_asset.source_path):
                missing_assets.append("dipole")
            logger.warning(f"  ⊘ atom_count_ok: Missing assets ({', '.join(missing_assets)})")

        # forming_bond_window_ok: Check TS or dipole
        if ts2_asset and ts2_asset.source_path:
            if resolved_forming_bonds:
                quality_flags.forming_bond_window_ok = _check_forming_bond_window(
                    ts2_asset.source_path, resolved_forming_bonds
                )
                if quality_flags.forming_bond_window_ok:
                    logger.info("  ✓ forming_bond_window_ok: Within window")
                else:
                    logger.warning("  ✗ forming_bond_window_ok: Outside window")
            else:
                quality_flags.forming_bond_window_ok = True
                logger.info("  ✓ forming_bond_window_ok: No forming bonds data")
        elif dipole_asset and dipole_asset.source_path:
            if resolved_forming_bonds:
                quality_flags.forming_bond_window_ok = _check_forming_bond_window(
                    dipole_asset.source_path, resolved_forming_bonds
                )
                if quality_flags.forming_bond_window_ok:
                    logger.info("  ✓ forming_bond_window_ok: Within window")
                else:
                    logger.warning("  ✗ forming_bond_window_ok: Outside window")
            else:
                quality_flags.forming_bond_window_ok = True
                logger.info("  ✓ forming_bond_window_ok: No forming bonds data")
        else:
            # M1-P0: Missing assets -> None (unknown/not applicable), not False
            quality_flags.forming_bond_window_ok = None
            logger.warning("  ⊘ forming_bond_window_ok: Missing TS and dipole assets")

    # M4-P0: Check suspect_optimized_to_product based on bond distances
    # This is NOT dependent on write_quality_flags - always run if data available
    ts_path = None
    product_path = None

    if ts2_asset and ts2_asset.source_path:
        ts_path = ts2_asset.source_path
    product_asset = assets.get('mech_step2_product')
    if product_asset and product_asset.source_path:
        product_path = product_asset.source_path

    suspect_result = _check_suspect_optimized_to_product(
        ts_path, product_path, resolved_forming_bonds
    )

    if suspect_result is None:
        quality_flags.suspect_optimized_to_product = None
        logger.debug("  suspect_optimized_to_product: None (insufficient data)")
    else:
        quality_flags.suspect_optimized_to_product = suspect_result
        if suspect_result == "suspect":
            logger.warning("  ⚠ suspect_optimized_to_product: suspect")
        else:
            logger.info("  ✓ suspect_optimized_to_product: ok")

    # M2-D2: Track missing inputs for degradation reporting
    missing_inputs_list = []
    degradation_reasons_list = []

    # Determine missing inputs
    if not assets.get('mech_step2_ts2'):
        missing_inputs_list.append("mech_step2_ts2.xyz")
    if not assets.get('mech_step2_reactant_dipole'):
        missing_inputs_list.append("mech_step2_reactant_dipole.xyz")
    if not assets.get('mech_step2_product'):
        missing_inputs_list.append("mech_step2_product.xyz")
    if not assets.get('mech_step1_precursor'):
        missing_inputs_list.append("mech_step1_precursor.xyz")

    # M2-D2: Track degradation reasons
    if not assets.get('mech_step2_ts2'):
        degradation_reasons_list.append("S3 TS optimization not completed")
    # M1-P0: Only add quality failure reason if check actually ran and FAILED (False), not if unknown (None)
    if quality_flags.atom_count_ok is False:
        degradation_reasons_list.append("Atom count mismatch between dipole and TS2")
    if quality_flags.forming_bond_window_ok is False:
        degradation_reasons_list.append("Forming bond distances outside expected window")

    # Determine overall mechanism_status
    any_missing = any(asset is None for asset in assets.values())
    mechanism_status = "INCOMPLETE" if any_missing else "COMPLETE"

    # M2-E: Enhanced visibility for INCOMPLETE status (one-time warning with missing items)
    if mechanism_status == "INCOMPLETE":
        missing_summary = ", ".join(missing_inputs_list)
        degradation_summary = ", ".join(degradation_reasons_list)
        logger.warning(
            f"⚠️  Mechanism packaging INCOMPLETE:"
            f"\n  Missing inputs: {missing_summary}"
            f"\n  Degradation reasons: {degradation_summary}"
        )

    step2_meta = MechanismMetaStep2(
        forming_bonds=resolved_forming_bonds,
        source_summary="derived_from_s1_s3"
    )

    step1_meta = MechanismMetaStep1(
        ts_computed=False,
        degradation_reason=(
            f"Using {precursor_label} as precursor"
            if precursor and precursor_label != "none"
            else "Precursor not available"
        )
    )

    # M4-D-4: Collect and copy QC artifacts
    # P0-1 FIX: Use out_dir.parent as pipeline_root (single source of truth for job root)
    # This ensures candidate.relative_to(pipeline_root) works correctly
    pipeline_root = out_dir.parent

    # Derive s3_dir from context, or fallback to pipeline_root/S3_TS
    s3_dir = context.s3_dir if context.s3_dir else (pipeline_root / "S3_TS")
    qc_artifacts = _collect_qc_artifacts(
        s3_dir=s3_dir,
        pipeline_root=pipeline_root,
        out_dir=out_dir,
        copy_mode=copy_mode
    )

    # Build mech_index.json
    mech_index = _build_mech_index(
        assets=assets,
        quality_flags=quality_flags,
        step2_meta=step2_meta,
        step1_meta=step1_meta,
        config=config,
        missing_inputs=missing_inputs_list,
        degradation_reasons=degradation_reasons_list,
        qc_artifacts=qc_artifacts
    )

    # M1-P0: Write all output files to disk with atomic rename
    # Write mech_index.json (atomic rename)
    mech_index_path = out_dir / "mech_index.json"
    _write_json_atomic(mech_index, mech_index_path)

    # Write mech_step2_meta.json
    mech_step2_meta = {
        "filename": "mech_step2_meta.json",
        "forming_bonds": [list(b) for b in step2_meta.forming_bonds] if step2_meta.forming_bonds else None,
        "source_summary": step2_meta.source_summary
    }
    _write_json_atomic(mech_step2_meta, out_dir / "mech_step2_meta.json")

    # Write mech_step1_meta.json
    mech_step1_meta = {
        "filename": "mech_step1_meta.json",
        "ts_computed": step1_meta.ts_computed,
        "degradation_reason": step1_meta.degradation_reason
    }
    _write_json_atomic(mech_step1_meta, out_dir / "mech_step1_meta.json")

    logger.info("=" * 60)
    logger.info(f"Mechanism assets packaged to: {out_dir}")
    logger.info(f"  Status: {mech_index.get('mechanism_status', 'UNKNOWN')}")
    logger.info(f"  Assets: {sum(1 for a in assets.values() if a is not None)}/{len(assets)}")
    logger.info(f"  Written: mech_index.json, mech_step1_meta.json, mech_step2_meta.json")
    logger.info("=" * 60)

    return mech_index
