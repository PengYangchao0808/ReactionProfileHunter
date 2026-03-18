"""
Unified Naming Conventions
=========================

Standardized naming for ReactionProfileHunter across all reaction types.
This module defines the single source of truth for all naming.

Naming Schema:
- S1: substrate (starting material before reaction)
- S2: intermediate (reaction intermediate/zwitterion)
- S3: transition state + optimized intermediate
- S4: features extracted from above

Author: QCcalc Team
Date: 2026-03-16
"""

from pathlib import Path
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)


# =============================================================================
# File Names
# =============================================================================

INTERMEDIATE_XYZ = "intermediate.xyz"
TS_GUESS_XYZ = "ts_guess.xyz"
PRODUCT_XYZ = "product_min.xyz"
PRECURSOR_XYZ = "precursor_min.xyz"
NEUTRAL_PRECURSOR_XYZ = "neutral_precursor.xyz"
SUBSTRATE_XYZ = "substrate.xyz"
REACTANT_COMPLEX_XYZ = "reactant_complex.xyz"


# =============================================================================
# Directory Names
# =============================================================================

DIR_S1_CONFORMATION = "S1_ConfGeneration"
DIR_S2_RETRO = "S2_Retro"
DIR_S3_TRANSITION = "S3_TransitionAnalysis"
DIR_S4_FEATURES = "S4_Data"


# =============================================================================
# S3 Subdirectories
# =============================================================================

DIR_INTERMEDIATE_OPT = "S3_intermediate_opt"
DIR_TS_OPT = "ts_opt"
DIR_SP_MATRIX = "ASM_SP_Mat"


# =============================================================================
# Source Labels for S4
# =============================================================================

SOURCE_S1_SUBSTRATE = "S1_substrate"
SOURCE_S2_INTERMEDIATE = "S2_intermediate"
SOURCE_S3_INTERMEDIATE = "S3_intermediate"
SOURCE_S3_TS = "S3_ts"


# =============================================================================
# Feature Name Prefixes
# =============================================================================

PREFIX_INTERMEDIATE = "intermediate_"
PREFIX_SUBSTRATE = "substrate_"


# =============================================================================
# Resolution Functions
# =============================================================================

def resolve_intermediate_path(directory: Path) -> Tuple[Path, str]:
    """Resolve intermediate file path."""
    directory = Path(directory)
    path = directory / INTERMEDIATE_XYZ
    if path.exists():
        return path, "intermediate"
    return path, "intermediate"


def resolve_s3_intermediate_opt_dir(s3_dir: Path) -> Path:
    """Resolve S3 intermediate optimization directory."""
    return s3_dir / DIR_INTERMEDIATE_OPT


def normalize_source_label(label: str) -> str:
    """Normalize source label to standard."""
    mapping = {
        "s2_reactant_complex": SOURCE_S2_INTERMEDIATE,
        "reactant_complex": SOURCE_S2_INTERMEDIATE,
        "s3_reactant": SOURCE_S3_INTERMEDIATE,
        "s3_reactant_sp": SOURCE_S3_INTERMEDIATE,
        "dipole": SOURCE_S2_INTERMEDIATE,
        "dipolar_intermediate": SOURCE_S2_INTERMEDIATE,
    }
    return mapping.get(label.lower(), label)


def get_intermediate_source_priority(config: dict) -> list:
    """Get intermediate source priority from config."""
    default = [SOURCE_S3_INTERMEDIATE, SOURCE_S2_INTERMEDIATE]
    
    priority = config.get("intermediate_source_priority")
    if priority and isinstance(priority, list):
        return [normalize_source_label(l) for l in priority]
    
    return default


def create_intermediate_alias(source_xyz: Path, output_dir: Path) -> Path:
    import shutil
    output_dir = Path(output_dir)
    alias_path = output_dir / REACTANT_COMPLEX_XYZ
    shutil.copy2(source_xyz, alias_path)
    logger.debug(f"Created backward compatibility alias: {alias_path}")
    return alias_path
