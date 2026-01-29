"""
Step 4: Feature Schema (V4.2 Phase A)
====================================

Fixed columns + dynamic columns support + schema signature + CSV writer.

Author: QC Descriptors Team
Date: 2026-01-18
"""

import hashlib
import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import pandas as pd
import numpy as np

SCHEMA_VERSION = "6.1"


def validate_forming_bonds(
    forming_bonds: Optional[Tuple[Tuple[int, int], ...]],
    natoms: Optional[int] = None
) -> Tuple[bool, List[str]]:
    """Validate forming_bonds tuple structure and values.

    Performs strict validation:
    # TS quality features (ts.*)
    "ts.n_imag",             # number of imaginary frequencies
    "ts.imag1_cm1_abs",      # absolute value of most negative frequency (cm-1)
    "ts.dipole_debye",       # dipole magnitude (Debye)
    
    # Strict-QC features (qc.*)
    "qc.sample_weight",       # sample weight (1.0 default, 0.0 if strict-QC downgrade)

    - Must be tuple/list of (int, int) pairs
    - Indices must be 0-based (>= 0)
    - Indices must be within natoms if provided
    - No i == j (no self-bonds)
    - No duplicate pairs (order-insensitive)
    - Minimum/maximum length validation (2 pairs max for typical reactions)

    Args:
        forming_bonds: Tuple of (atom_idx, atom_idx) pairs, or None
        natoms: Total number of atoms for bounds checking (optional)

    Returns:
        (is_valid, error_messages) tuple
    """
    errors = []

    if forming_bonds is None:
        return True, []

    if not isinstance(forming_bonds, (tuple, list)):
        errors.append(f"forming_bonds must be tuple or list, got {type(forming_bonds).__name__}")
        return False, errors

    if len(forming_bonds) == 0:
        errors.append("forming_bonds must not be empty")
        return False, errors

    if len(forming_bonds) > 2:
        errors.append(f"forming_bonds has {len(forming_bonds)} pairs; expected at most 2 for typical reactions")

    seen_pairs = set()

    for idx, pair in enumerate(forming_bonds):
        pair_name = f"pair {idx}"

        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            errors.append(f"{pair_name}: must be (int, int), got {pair}")
            continue

        i, j = pair

        if not (isinstance(i, int) and isinstance(j, int)):
            errors.append(f"{pair_name}: both indices must be ints, got ({type(i).__name__}, {type(j).__name__})")
            continue

        if i < 0 or j < 0:
            errors.append(f"{pair_name}: indices must be 0-based (>= 0), got ({i}, {j})")

        if i == j:
            errors.append(f"{pair_name}: self-bond detected (i==j={i})")

        if natoms is not None:
            if i >= natoms or j >= natoms:
                errors.append(f"{pair_name}: indices out of bounds (natoms={natoms}), got ({i}, {j})")

        canonical = tuple(sorted((i, j)))
        if canonical in seen_pairs:
            errors.append(f"{pair_name}: duplicate pair {canonical}")
        seen_pairs.add(canonical)

    is_valid = len(errors) == 0
    return is_valid, errors

# Fixed columns (in order) - Phase A
FIXED_COLUMNS = [
    # Schema metadata
    "schema_version",
    "schema_signature",
    "feature_status",

    # Thermo features (thermo.*)
    "thermo.dG_activation",  # kcal/mol, Gibbs or electronic activation energy (best available)
    "thermo.dG_reaction",   # kcal/mol, Gibbs or electronic reaction energy (best available)
    "thermo.dG_activation_gibbs",  # kcal/mol, Gibbs activation energy (V6.1)
    "thermo.dG_reaction_gibbs",   # kcal/mol, Gibbs reaction energy (V6.1)
    "thermo.dE_activation",  # kcal/mol, electronic activation energy
    "thermo.dE_reaction",   # kcal/mol, electronic reaction energy
    "thermo.energy_source_activation",  # "gibbs" or "electronic" (describes dG_*)
    "thermo.energy_source_reaction",    # "gibbs" or "electronic" (describes dG_*)
    "thermo.method",          # method string from sp_report
    "thermo.solvent",         # solvent string from sp_report

    # Geometry features (geom.*)
    "geom.natoms_ts",         # number of atoms in TS
    "geom.r1",                # forming bond 1 distance (Å)
    "geom.r2",                # forming bond 2 distance (Å)
    "geom.asynch",            # asynchronicity = |r1 - r2| (Å)
    "geom.asynch_index",      # asynchronicity index: asynch / (r1 + r2)
    "geom.rg_ts",             # radius of gyration of TS (Å)
    "geom.min_nonbonded",     # minimum non-bonded distance (Å)
    "geom.close_contacts",     # number of close contacts (< default 2.2 Å)
    # Geometry V6.1 additions
    "geom.r_avg",  # (r1+r2)/2, average forming bond distance
    "geom.dr",              # (r1-r2), signed asynchronicity
    "geom.close_contacts_density",  # close_contacts / natoms_ts
    

    # QC validation features (qc.*)
    "qc.has_gibbs",               # 1 if g_* present, else 0
    "qc.used_fallback_electronic", # 1 if fallback to electronic energy, else 0
    "qc.sp_report_validated",     # 1 if sp_report.validate() succeeded, else 0
    "qc.forming_bonds_valid",      # 1 if forming_bonds provided and valid, else 0
    "qc.warnings_count",           # total warnings across all plugins
]

# Default close contacts cutoff (Å)
DEFAULT_CLOSE_CONTACTS_CUTOFF = 2.2


def get_schema_signature(*, features: Dict[str, Any], enabled_plugins: List[str]) -> str:
    """Compute schema signature hash from schema metadata (not feature values).

    The signature represents the schema itself: schema_version, enabled_plugins,
    and column names. This ensures the same schema produces the same signature
    regardless of actual feature values.

    Args:
        features: Feature dictionary with all keys (including dynamic)
        enabled_plugins: List of enabled plugin names for this run

    Returns:
        sha1 hex digest representing the schema signature

    Notes:
        - Hash includes: schema_version, sorted enabled_plugins, sorted column names
        - Does NOT hash feature values (those change per run)
        - Plugin ordering is insensitive (sorted before hashing)
    """
    columns = sorted(features.keys())
    sorted_plugins = sorted(enabled_plugins)

    schema_data = {
        "schema_version": SCHEMA_VERSION,
        "enabled_plugins": sorted_plugins,
        "columns": columns,
    }

    # Serialize to JSON with sorted keys
    signature_data = json.dumps(schema_data, sort_keys=True, default=str)

    # Compute SHA1 hash
    return hashlib.sha1(signature_data.encode('utf-8')).hexdigest()


def write_features_csv(output_path: Path, row_data: Dict[str, Any]) -> None:
    """Write feature row to CSV file with fixed columns first, then dynamic columns.

    Args:
        output_path: Output CSV file path
        row_data: Feature dictionary (may include dynamic columns beyond FIXED_COLUMNS)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Separate fixed and dynamic columns
    fixed_data = {}
    dynamic_data = {}

    for key, value in row_data.items():
        if key in FIXED_COLUMNS:
            fixed_data[key] = value
        else:
            dynamic_data[key] = value

    # Ensure all fixed columns exist (fill with NaN if missing)
    for col in FIXED_COLUMNS:
        if col not in fixed_data:
            fixed_data[col] = np.nan

    # Combine: fixed columns (in order) + dynamic columns (sorted)
    ordered_data = {col: fixed_data[col] for col in FIXED_COLUMNS}
    ordered_data.update({k: dynamic_data[k] for k in sorted(dynamic_data)})

    # Write to CSV
    df = pd.DataFrame([ordered_data])
    df.to_csv(output_path, index=False)


class FeatureSchema:
    """Feature schema manager for V4.2 Phase A."""

    def __init__(self, close_contacts_cutoff: float = DEFAULT_CLOSE_CONTACTS_CUTOFF):
        """Initialize feature schema.

        Args:
            close_contacts_cutoff: Cutoff distance for close contacts detection (Å)
        """
        self.close_contacts_cutoff = close_contacts_cutoff

    def get_all_columns(self, features: Dict[str, Any]) -> List[str]:
        """Get all columns including dynamic ones.

        Args:
            features: Feature dictionary

        Returns:
            List of column names (fixed + dynamic sorted)
        """
        fixed = FIXED_COLUMNS
        dynamic = sorted([k for k in features if k not in FIXED_COLUMNS])
        return fixed + dynamic

    def validate_row(self, row_data: Dict[str, Any]) -> List[str]:
        """Validate feature row for schema compliance.

        Args:
            row_data: Feature dictionary

        Returns:
            List of warning messages (empty if valid)
        """
        warnings = []

        # Check required schema metadata
        if "schema_version" not in row_data or row_data["schema_version"] != SCHEMA_VERSION:
            warnings.append(f"schema_version missing or invalid: expected {SCHEMA_VERSION}")

        if "feature_status" not in row_data:
            warnings.append("feature_status missing")

        # Check thermo features (at least one energy should be present)
        thermo_keys = [
            "thermo.dG_activation", "thermo.dG_reaction",
            "thermo.dE_activation", "thermo.dE_reaction"
        ]

        def _is_valid_energy(val):
            return val is not None and not (isinstance(val, float) and np.isnan(val))

        has_energy = any(_is_valid_energy(row_data.get(k)) for k in thermo_keys)
        if not has_energy:
            warnings.append("No thermo energy values found")

        return warnings


def get_units_for_column(column_name: str) -> Optional[str]:
    """Get unit string for a given column name.

    Args:
        column_name: Column name (e.g., "thermo.dG_activation", "geom.r1")

    Returns:
        Unit string (e.g., "kcal/mol", "Å") or None if no unit
    """
    unit_map = {
        "thermo.dG_activation": "kcal/mol",
        "thermo.dG_reaction": "kcal/mol",
        "thermo.dG_activation_gibbs": "kcal/mol",
        "thermo.dG_reaction_gibbs": "kcal/mol",
        "thermo.dE_activation": "kcal/mol",
        "thermo.dE_reaction": "kcal/mol",

        "geom.r1": "Å",
        "geom.r2": "Å",
        "geom.asynch": "Å",
        "geom.rg_ts": "Å",
        "geom.min_nonbonded": "Å",

        "geom.asynch_index": None,
        "geom.close_contacts": None,
        "qc.warnings_count": None,

        "qc.has_gibbs": None,
        "qc.used_fallback_electronic": None,
        "qc.sp_report_validated": None,
        "qc.forming_bonds_valid": None,
    }
    return unit_map.get(column_name)


def generate_cache_key(
    plugin_name: str,
    input_files: Dict[str, Path],
    params: Optional[Dict[str, Any]] = None,
) -> str:
    """Generate cache key for plugin results (Contract 5).

    Phase B/C: used for caching expensive QM calculations.
    Phase A: interface is locked but caching not yet active.

    Args:
        plugin_name: Name of the plugin generating the cache
        input_files: Dictionary mapping canonical key -> Path for input files
        params: Additional parameters affecting computation

    Returns:
        Cache key string (sha1[:16] of normalized input)

    Contract specification:
        - File fingerprint includes size, mtime, and sha256[:8] of file content
        - Params are serialized with json sort_keys=True
        - Final output is sha1[:16] of combined data
    """
    import json

    key_data = {
        "plugin": plugin_name,
        "files": {},
        "params": params or {},
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
                "sha256_prefix": content_hash,
            }
        except Exception:
            continue

    key_str = json.dumps(key_data, sort_keys=True, default=str)
    full_hash = hashlib.sha1(key_str.encode('utf-8')).hexdigest()
    return full_hash[:16]


JOB_SPEC_REQUIRED_KEYS = {
    "engine",
    "command",
    "workdir",
    "input_files",
    "expected_outputs",
    "cache_key",
}

JOB_SPEC_OPTIONAL_KEYS = {
    "stdin_lines",
    "recipe_name",
    "recipe_version",
    "template_yaml_fingerprint",
    "recipe_yaml_fingerprint",
}

# MLR columns (V6.1) - default linear-model-ready columns (<=10)
DEFAULT_MLR_COLUMNS = [
    "sample_id",
    "thermo.dE_activation",
    "thermo.dE_reaction",
    # Geometry V6.1 additions
    "geom.r_avg",  # (r1+r2)/2, average forming bond distance
    "geom.dr",
    "geom.close_contacts_density",
    "ts.imag1_cm1_abs",
    "asm.distortion_total_kcal",
    "fmo.dipolar_omega_ev",
    "qc.sample_weight",
]


def write_features_raw_csv(output_path: Path, row_data: Dict[str, Any]) -> None:
    """Write feature row to features_raw.csv with fixed columns first, then dynamic columns.

    This is the full feature table including all QC columns and intermediate convenience columns.

    Args:
        output_path: Output CSV file path (features_raw.csv)
        row_data: Feature dictionary (may include dynamic columns beyond FIXED_COLUMNS)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Separate fixed and dynamic columns
    fixed_data = {}
    dynamic_data = {}

    for key, value in row_data.items():
        if key in FIXED_COLUMNS:
            fixed_data[key] = value
        else:
            dynamic_data[key] = value

    # Ensure all fixed columns exist (fill with NaN if missing)
    for col in FIXED_COLUMNS:
        if col not in fixed_data:
            fixed_data[col] = np.nan

    # Combine: fixed columns (in order) + dynamic columns (sorted)
    ordered_data = {col: fixed_data[col] for col in FIXED_COLUMNS}
    ordered_data.update({k: dynamic_data[k] for k in sorted(dynamic_data)})

    # Write to CSV
    df = pd.DataFrame([ordered_data])
    df.to_csv(output_path, index=False)


def write_features_mlr_csv(
    output_path: Path,
    row_data: Dict[str, Any],
    mlr_columns: Optional[List[str]] = None
) -> None:
    """Write MLR-ready feature row to features_mlr.csv with exactly specified columns.

    MLR output contains ONLY the columns specified, with NaN for missing values.
    No schema_version, schema_signature, or feature_status columns are included.

    Args:
        output_path: Output CSV file path (features_mlr.csv)
        row_data: Feature dictionary with all available features
        mlr_columns: List of column names to write. If None, uses DEFAULT_MLR_COLUMNS.
    """
    if mlr_columns is None:
        mlr_columns = DEFAULT_MLR_COLUMNS

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build row with exactly the specified columns, using NaN for missing values
    ordered_data = {}
    for col in mlr_columns:
        ordered_data[col] = row_data.get(col, np.nan)

    # Write to CSV
    df = pd.DataFrame([ordered_data])
    df.to_csv(output_path, index=False)


# Keep legacy write_features_csv as alias for write_features_raw_csv for compatibility
write_features_csv = write_features_raw_csv
