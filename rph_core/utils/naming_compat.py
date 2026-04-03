from pathlib import Path
from typing import Any, Dict, List, Tuple

INTERMEDIATE_XYZ = "intermediate.xyz"
TS_GUESS_XYZ = "ts_guess.xyz"
PRODUCT_XYZ = "product_min.xyz"
PRECURSOR_XYZ = "precursor_min.xyz"
NEUTRAL_PRECURSOR_XYZ = "neutral_precursor.xyz"
SUBSTRATE_XYZ = "substrate.xyz"
REACTANT_COMPLEX_XYZ = INTERMEDIATE_XYZ

DIR_S0_MECHANISM = "S0_Mechanism"
DIR_S1_CONFORMATION = "S1_ConfGeneration"
DIR_S2_RETRO = "S2_Retro"
DIR_S3_TRANSITION = "S3_TS"
DIR_S4_FEATURES = "S4_Data"

DIR_INTERMEDIATE_OPT = "S3_intermediate_opt"
DIR_TS_OPT = "ts_opt"
DIR_SP_MATRIX = "ASM_SP_Mat"

SOURCE_S1_SUBSTRATE = "S1_substrate"
SOURCE_S2_INTERMEDIATE = "S2_intermediate"
SOURCE_S3_INTERMEDIATE = "S3_intermediate"
SOURCE_S3_TS = "S3_ts"

PREFIX_INTERMEDIATE = "intermediate_"
PREFIX_SUBSTRATE = "substrate_"


def resolve_intermediate_path(directory: Path) -> Tuple[Path, str]:
    path = Path(directory) / INTERMEDIATE_XYZ
    return path, "intermediate"


def resolve_s3_intermediate_opt_dir(s3_dir: Path) -> Path:
    return Path(s3_dir) / DIR_INTERMEDIATE_OPT


def normalize_source_label(label: str) -> str:
    normalized = str(label).strip().lower()
    if normalized in {"s2_intermediate", "intermediate", "dipole", "dipolar_intermediate"}:
        return SOURCE_S2_INTERMEDIATE
    if normalized in {"s3_intermediate", "s3_reactant_sp", "s3_reactant"}:
        return SOURCE_S3_INTERMEDIATE
    return label


def get_intermediate_source_priority(config: Dict[str, Any]) -> List[str]:
    default = [SOURCE_S3_INTERMEDIATE, SOURCE_S2_INTERMEDIATE]
    priority = config.get("intermediate_source_priority")
    if isinstance(priority, list):
        return [normalize_source_label(item) for item in priority]
    return default
