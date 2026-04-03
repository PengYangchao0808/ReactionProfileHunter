"""
S3 Artifact Resolver
===================

Resolves S3 inputs from S2 artifacts with support for both new role-based
naming and legacy naming conventions.

New mode:
    - start_structure (reactant_complex) + product -> ts_guess

Legacy mode:
- ts_guess.xyz + reactant_complex.xyz + product.xyz
"""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


class S3ArtifactSource:
    """Source tracking for S3 inputs"""
    TS_GUESS = "ts_guess"
    REACTANT = "reactant"
    PRODUCT = "product"


@dataclass
class S3InputPaths:
    ts_guess: Path
    reactant: Path
    product: Path
    source_ts_guess: str
    source_reactant: str
    source_product: str
    metadata: Dict[str, Any]


def resolve_s3_inputs(
    s2_dir: Path,
    product_xyz: Path,
) -> S3InputPaths:
    """
    Resolve S3 input paths from S2 artifacts.
    
    Priority:
    1. Check S2 metadata for new role-based naming
    2. Fall back to legacy files (ts_guess.xyz, reactant_complex.xyz)
    
    Args:
        s2_dir: S2 output directory (S2_Retro)
        product_xyz: Product XYZ path (from S1)
        
    Returns:
        S3InputPaths with resolved paths and source tracking
    """
    s2_dir = Path(s2_dir)
    
    metadata = _load_s2_metadata(s2_dir)
    
    if metadata and metadata.get("generation_method") == "xtb_path_search":
        return _resolve_new_mode(s2_dir, product_xyz, metadata)
    else:
        return _resolve_legacy_mode(s2_dir, product_xyz)


def _load_s2_metadata(s2_dir: Path) -> Optional[Dict[str, Any]]:
    profile_path = s2_dir / "scan_profile.json"
    if not profile_path.exists():
        return None
    
    try:
        with open(profile_path) as f:
            return json.load(f)
    except Exception as e:
        logger.warning(f"Failed to load S2 metadata: {e}")
        return None


def _resolve_new_mode(
    s2_dir: Path,
    product_xyz: Path,
    metadata: Dict[str, Any],
) -> S3InputPaths:
    logger.info("Resolving S3 inputs: new mode (xtb_path_search)")
    
    start_role = metadata.get("start_structure_role", "intermediate")
    
    ts_guess = s2_dir / "ts_guess.xyz"
    if not ts_guess.exists():
        raise FileNotFoundError(f"ts_guess.xyz not found in {s2_dir}")
    
    reactant_path = s2_dir / "intermediate.xyz"
    is_alias = False
    
    if is_alias:
        logger.info("Using intermediate as reactant")
    
    product_path = Path(product_xyz)
    if not product_path.exists():
        raise FileNotFoundError(f"Product xyz not found: {product_xyz}")
    
    return S3InputPaths(
        ts_guess=ts_guess,
        reactant=reactant_path,
        product=product_path,
        source_ts_guess=S3ArtifactSource.TS_GUESS,
        source_reactant="intermediate",
        source_product=S3ArtifactSource.PRODUCT,
        metadata={
            "generation_method": metadata.get("generation_method"),
            "start_structure_role": start_role,
            "is_reactant_alias": is_alias,
        },
    )


def _resolve_legacy_mode(
    s2_dir: Path,
    product_xyz: Path,
) -> S3InputPaths:
    """Resolve inputs from legacy scan mode"""
    logger.info("Resolving S3 inputs: legacy mode")
    
    ts_guess = s2_dir / "ts_guess.xyz"
    if not ts_guess.exists():
        raise FileNotFoundError(f"ts_guess.xyz not found in {s2_dir}")
    
    reactant_candidates = [
        s2_dir / "intermediate.xyz",
    ]
    
    reactant = None
    for candidate in reactant_candidates:
        if candidate.exists():
            reactant = candidate
            break
    
    if reactant is None:
        raise FileNotFoundError(
            f"No reactant found in {s2_dir}. "
            f"Checked: {[c.name for c in reactant_candidates]}"
        )
    
    product_path = Path(product_xyz)
    if not product_path.exists():
        raise FileNotFoundError(f"Product xyz not found: {product_xyz}")
    
    return S3InputPaths(
        ts_guess=ts_guess,
        reactant=reactant,
        product=product_path,
        source_ts_guess=S3ArtifactSource.TS_GUESS,
        source_reactant=S3ArtifactSource.REACTANT,
        source_product=S3ArtifactSource.PRODUCT,
        metadata={"generation_method": "legacy_scan"},
    )


def check_s2_artifacts(s2_dir: Path) -> Dict[str, bool]:
    """Check which S2 artifacts exist"""
    s2_dir = Path(s2_dir)
    
    return {
        "ts_guess_exists": (s2_dir / "ts_guess.xyz").exists(),
        "intermediate_exists": (s2_dir / "intermediate.xyz").exists(),
        "scan_profile_exists": (s2_dir / "scan_profile.json").exists(),
    }
