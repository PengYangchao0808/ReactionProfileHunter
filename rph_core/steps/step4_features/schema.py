"""
Step 4: Feature Schema (V4.2 Phase A)
====================================

Fixed columns + dynamic columns support + schema signature + CSV writer.

Author: QC Descriptors Team
Date: 2026-01-18
"""

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import pandas as pd
import numpy as np

SCHEMA_VERSION = "6.2"

# V6.2 Warning Codes
WARNING_CODES = {
    # S1 warnings
    "W_S1_MISSING_SHERMO": "S1 Shermo summary file not found",
    "W_S1_MISSING_HOAC_THERMO": "S1 HOAc thermochemistry file not found",
    "W_S1_NO_CONFORMER_ENSEMBLE": "S1 conformer ensemble not available",
    "W_S1_MISSING_S3_REACTANT_GIBBS": "S3 reactant Gibbs not available (Shermo required)",
    "W_S3_SHERMO_GIBBS_FAILED": "S3 Shermo Gibbs calculation failed",
    # CDFT warnings
    "W_CDFT_MISSING_ORBITALS": "CDFT orbital energies not available",
    "W_CDFT_ORBITAL_ORDER": "CDFT LUMO <= HOMO (orbital order invalid)",
    # P1 Enhancement: CDFT sanity check warnings
    "W_CDFT_HOMO_RANGE": "CDFT HOMO outside typical range (-30 to 0 eV)",
    "W_CDFT_LUMO_RANGE": "CDFT LUMO outside typical range (-5 to 20 eV)",
    "W_CDFT_GAP_RANGE": "CDFT gap (eta) outside typical range (0.5 to 15 eV)",
    "W_CDFT_OMEGA_RANGE": "CDFT omega outside typical range (0 to 20 eV)",
    # GEDT warnings
    "W_GEDT_NO_CHARGES": "GEDT atomic charges not available",
    "W_GEDT_FRAGMENT_CUT_FAILED": "GEDT fragment cutting failed",
    # TS validity warnings
    "W_TS_MISSING_LOG": "TS log file not found for frequency validation",
    "W_TS_IMAG_COUNT_NOT_ONE": "TS does not have exactly one imaginary frequency",
    "W_TS_NO_MODE_VECTOR": "TS imaginary frequency mode vector not available",
    # Multiwfn warnings
    "W_MULTIWFN_DISABLED": "Multiwfn is disabled in config",
    "W_MULTIWFN_FAILED": "Multiwfn execution failed",
    "W_MULTIWFN_CACHE_READ_FAILED": "Multiwfn cache read failed",
    "W_MULTIWFN_TIMEOUT": "Multiwfn execution timed out",
    "W_MULTIWFN_INVALID_OUTPUT": "Multiwfn output parsing failed",
    # Multiwfn extractor specific
    "W_MW_MISSING_INPUT": "Multiwfn missing required input files",
    "W_MW_NO_FCHK": "Multiwfn fchk file not available",
    "W_MW_NO_FORMING_BONDS": "Multiwfn forming bonds not defined",
    "W_MW_ATOM_NOT_FOUND": "Multiwfn atom index not found",
}

# MLR columns (V6.2) - versioned for backward compatibility
MLR_COLUMNS_V2 = [
    "sample_id",
    "thermo.dE_activation",
    "thermo.dE_reaction",
    "geom.r_avg",
    "geom.dr",
    "geom.close_contacts_density",
    "ts.imag1_cm1_abs",
    # V6.2 new columns
    "s1_dG_act",
    "s1_Keq_act",
    "s1_Nconf_eff",
    "s1_Sconf",
    # P1 Enhancement: Add thermodynamic average features
    "s1_E_avg_weighted",
    "s1_E_std",
    "s1_tau_CH_C_O",
    "s2_d_forming_1",
    "s2_d_forming_2",
    "s2_asynch",
    "s2_eps_homo",
    "s2_eps_lumo",
    "s2_mu",
    "s2_eta",
    "s2_omega",
    "s2_gedt_value",
    "s2_ts_validity_flag",
    "qc.sample_weight",
]

# MLR columns (deduplicated)
#
# V6.2 originally emitted Step2 TS-geometry and TS-validity keys that duplicate
# geom.* and ts.*. V3 keeps Step2 electronic (CDFT/GEDT) only.
MLR_COLUMNS_V3_DEDUP = [
    "sample_id",
    "thermo.dE_activation",
    "thermo.dE_reaction",
    "geom.r_avg",
    "geom.dr",
    "geom.close_contacts_density",
    "ts.n_imag",
    "ts.imag1_cm1_abs",
    "s1_dG_act",
    "s1_Keq_act",
    "s1_Nconf_eff",
    "s1_Sconf",
    "s1_E_avg_weighted",
    "s1_E_std",
    "s1_tau_CH_C_O",
    "s2_eps_homo",
    "s2_eps_lumo",
    "s2_mu",
    "s2_eta",
    "s2_omega",
    "s2_gedt_value",
    "qc.sample_weight",
]

# =============================================================================
# V6.4 Feature Layers (Three-Layer Separation)
# =============================================================================
# Layer 1: Reaction Data Features (from reaxys_cleaned.csv)
# Layer 2: Computational Chemistry Features (from S1/S2/S3 QC)
# Layer 3: QA Metadata Features (confidence indicators - NOT for ML training)
# =============================================================================

# Layer 1: Reaction Data Features (from experimental data source)
# These are injected from reaxys_cleaned.csv or similar
REACTION_DATA_FEATURES: List[str] = [
    "ctx_temperature",       # Reaction temperature (°C or K)
    "ctx_solvent_ep",         # Solvent dielectric constant
    "ctx_topology",           # Reaction topology (INTER/INTRA)
    "ctx_is_stepwise",        # Mechanism flag (0=concerted, 1=stepwise)
    "label_yield",            # Reaction yield (0-100)
    "label_yield_fraction",   # Yield as fraction (0-1)
]

# Layer 2: Computational Chemistry Features (from S1/S2/S3 QC outputs)
# These are the core features for ML training - NO QA metadata
COMPUTATIONAL_FEATURES: List[str] = [
    # Thermodynamic
    "thermo.dE_activation",
    "thermo.dE_reaction",
    "thermo.dG_activation",
    "thermo.dG_reaction",
    # Geometry
    "geom.natoms_ts",
    "geom.r1",
    "geom.r2",
    "geom.asynch",
    "geom.asynch_index",
    "geom.r_avg",
    "geom.dr",
    "geom.rg_ts",
    "geom.min_nonbonded",
    "geom.close_contacts",
    "geom.close_contacts_density",
    # Step1 features
    "s1_dG_act",
    "s1_Keq_act",
    "s1_Nconf_eff",
    "s1_Nconf_total",
    "s1_Sconf",
    "s1_E_span",
    "s1_E_avg_weighted",
    "s1_G_avg_weighted",
    "s1_E_variance",
    "s1_E_std",
    "s1_d_C_O_lg",
    "s1_angle_O_C_O",
    "s1_tau_lg_ring",
    "s1_d_C_H_alpha",
    "s1_tau_CH_C_O",
    # Step2 features (electronic structure)
    "s2_eps_homo",
    "s2_eps_lumo",
    "s2_mu",
    "s2_eta",
    "s2_omega",
    "s2_gedt_value",
    # Multiwfn features (local reactivity)
    "mw_fukui_fplus_atomA",
    "mw_fukui_fplus_atomB",
    "mw_fukui_fminus_atomA",
    "mw_fukui_fminus_atomB",
    "mw_fukui_f0_atomA",
    "mw_fukui_f0_atomB",
    "mw_dual_descriptor_atomA",
    "mw_dual_descriptor_atomB",
    "mw_rho_bcp_forming1",
    "mw_laplacian_bcp_forming1",
]

# Layer 3: QA Metadata Features (confidence indicators - NOT for training)
# These should be used for sample filtering, NOT as ML features
QA_METADATA_FEATURES: List[str] = [
    # TS quality
    "ts.n_imag",
    "ts.imag1_cm1_abs",
    "ts.dipole_debye",
    # QC validation
    "qc.has_gibbs",
    "qc.used_fallback_electronic",
    "qc.sp_report_validated",
    "qc.forming_bonds_valid",
    "qc.warnings_count",
    "qc.sample_weight",
    # TS validity flags
    "s2_ts_validity_flag",
    "s2_n_imag_freq",
    "s2_imag_freq_cm1",
    # Calculation level
    "calc_level_hash",
    "core_extraction_status",
    # Baseline for verification
    "thermo.dG_activation_baseline",
]

# V6.4 Deployable Features (Layer 1 + Layer 2, NO Layer 3)
# This is the SAFE feature set for ML training - no QA leakage
DEPLOYABLE_COLUMNS_V1: List[str] = [
    "sample_id",
    # Layer 1: Reaction Data (placeholder - requires ReactionContext injection)
    # "ctx_temperature",
    # "ctx_solvent_ep",
    # "ctx_topology",
    # "label_yield_fraction",
    # Layer 2: Computational Chemistry (no QA)
    "thermo.dE_activation",
    "thermo.dE_reaction",
    "geom.r_avg",
    "geom.dr",
    "geom.close_contacts_density",
    "s1_dG_act",
    "s1_Keq_act",
    "s1_Nconf_eff",
    "s1_Sconf",
    "s1_E_avg_weighted",
    "s1_E_std",
    "s1_tau_CH_C_O",
    "s2_eps_homo",
    "s2_eps_lumo",
    "s2_mu",
    "s2_eta",
    "s2_omega",
    "s2_gedt_value",
]


def validate_feature_layer(feature_name: str) -> str:
    """Validate which layer a feature belongs to.

    Args:
        feature_name: Name of the feature

    Returns:
        Layer name: "reaction_data", "computational", "qa_metadata", or "unknown"
    """
    if feature_name in REACTION_DATA_FEATURES:
        return "reaction_data"
    if feature_name in COMPUTATIONAL_FEATURES:
        return "computational"
    if feature_name in QA_METADATA_FEATURES:
        return "qa_metadata"
    # Dynamic features default to computational
    if feature_name.startswith(("s1_", "s2_", "mw_", "thermo.", "geom.")):
        return "computational"
    if feature_name.startswith(("ts.", "qc.", "calc_", "core_")):
        return "qa_metadata"
    if feature_name.startswith("ctx_") or feature_name.startswith("label_"):
        return "reaction_data"
    return "unknown"


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

        # V6.2 Step1 features
        "s1_dG_act": "kcal/mol",
        "s1_Keq_act": None,
        "s1_Sconf": "cal/mol/K",
        "s1_Nconf_eff": None,
        "s1_Nconf_total": None,
        "s1_E_span": "kcal/mol",
        # P1 Enhancement: Conformer thermodynamic features
        "s1_E_avg_weighted": "kcal/mol",
        "s1_G_avg_weighted": "kcal/mol",
        "s1_E_variance": "(kcal/mol)^2",
        "s1_E_std": "kcal/mol",
        "s1_d_C_O_lg": "Å",
        "s1_angle_O_C_O": "degrees",
        "s1_tau_lg_ring": "degrees",
        "s1_d_C_H_alpha": "Å",
        "s1_tau_CH_C_O": "degrees",

        # V6.2 Step2 features
        "s2_dGddagger": "kcal/mol",
        "s2_dHddagger": "kcal/mol",
        "s2_dSddagger": "cal/mol/K",
        "s2_TdSddagger": "kcal/mol",
        "s2_dGrxn": "kcal/mol",
        "s2_d_forming_1": "Å",
        "s2_d_forming_2": "Å",
        "s2_asynch": "Å",
        "s2_eps_homo": "eV",
        "s2_eps_lumo": "eV",
        "s2_mu": "eV",
        "s2_eta": "eV",
        "s2_omega": "eV",
        # P1 Enhancement: CDFT unit annotation
        "s2_cdft_unit": None,
        "s2_gedt_value": "e",
        "s2_n_imag_freq": None,
        "s2_imag_freq_cm1": "cm-1",
        "s2_imag_mode_overlap_forming": None,
        "s2_ts_validity_flag": None,

        # P2 Enhancement: Multiwfn features
        "mw_fukui_fplus_atomA": "e",
        "mw_fukui_fplus_atomB": "e",
        "mw_fukui_fminus_atomA": "e",
        "mw_fukui_fminus_atomB": "e",
        "mw_fukui_f0_atomA": "e",
        "mw_fukui_f0_atomB": "e",
        "mw_dual_descriptor_atomA": "e",
        "mw_dual_descriptor_atomB": "e",
        "mw_rho_bcp_forming1": "au",
        "mw_laplacian_bcp_forming1": "au",
        "mw_status": None,
        "mw_missing_reason": None,
        "mw_warnings_count": None,
        "mw_cache_hit": None,
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

    key_data: Dict[str, Any] = {
        "plugin": plugin_name,
        "files": {},
        "params": params or {},
    }
    files_data: Dict[str, Dict[str, Any]] = key_data["files"]

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

            files_data[key] = {
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

# Default MLR columns
#
# If caller does not provide an explicit column list, we default to the
# deduplicated v6.2+ set.
DEFAULT_MLR_COLUMNS = [
    # Default to the deduplicated v6.2+ column set.
    *MLR_COLUMNS_V3_DEDUP,
]


def _write_dataframe_atomic(output_path: Path, df: pd.DataFrame) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{output_path.name}.", suffix=".tmp", dir=str(output_path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        df.to_csv(tmp_path, index=False)
        tmp_path.replace(output_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)


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

    df = pd.DataFrame([ordered_data])
    _write_dataframe_atomic(output_path, df)


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

    df = pd.DataFrame([ordered_data])
    _write_dataframe_atomic(output_path, df)


# Keep legacy write_features_csv as alias for write_features_raw_csv for compatibility
write_features_csv = write_features_raw_csv
