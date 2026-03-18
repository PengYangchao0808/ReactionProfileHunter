"""
Step 4: Step2 Cyclization Extractor (V6.2 P0)
==============================================

Extract Step2 ([5+2] cycloaddition) features:
- Kinetic features (deprecated: owned by thermo.*)
- CDFT indices (HOMO, LUMO, mu, eta, omega)
- GEDT (Global Electron Density Transfer)

NOTE:
- TS geometry is emitted by the geometry extractor (geom.*)
- TS frequency/validity is emitted by the ts_quality extractor (ts.*)

Author: RPH Team
Date: 2026-02-02
"""

import re
import numpy as np
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

from .base import BaseExtractor, register_extractor
from rph_core.steps.step4_features.schema import validate_forming_bonds
from rph_core.utils.file_io import read_xyz
from rph_core.utils.geometry_tools import GeometryUtils
from rph_core.utils.fchk_reader import read_fchk_cdft_indices
from rph_core.utils.fragment_cut import cut_along_forming_bonds
from rph_core.utils.charge_reader import read_charges_with_priority


# Warning codes for Step2
W_CDFT_MISSING_ORBITALS = "W_CDFT_MISSING_ORBITALS"
W_CDFT_ORBITAL_ORDER = "W_CDFT_ORBITAL_ORDER"
# P1 Enhancement: CDFT sanity check warnings
W_CDFT_HOMO_RANGE = "W_CDFT_HOMO_RANGE"       # HOMO outside typical range (-30 to 0 eV)
W_CDFT_LUMO_RANGE = "W_CDFT_LUMO_RANGE"       # LUMO outside typical range (-5 to 20 eV)
W_CDFT_GAP_RANGE = "W_CDFT_GAP_RANGE"         # Gap (eta) outside typical range (0.5 to 15 eV)
W_CDFT_OMEGA_RANGE = "W_CDFT_OMEGA_RANGE"     # Omega outside typical range (0 to 20 eV)
W_GEDT_NO_CHARGES = "W_GEDT_NO_CHARGES"
W_GEDT_FRAGMENT_CUT_FAILED = "W_GEDT_FRAGMENT_CUT_FAILED"
W_TS_MISSING_LOG = "W_TS_MISSING_LOG"
W_TS_IMAG_COUNT_NOT_ONE = "W_TS_IMAG_COUNT_NOT_ONE"
W_TS_NO_MODE_VECTOR = "W_TS_NO_MODE_VECTOR"
W_FORMING_BONDS_INVALID = "W_FORMING_BONDS_INVALID"

# CDFT sanity check bounds (P1 Enhancement)
CDFT_HOMO_MIN = -30.0  # eV
CDFT_HOMO_MAX = 0.0    # eV
CDFT_LUMO_MIN = -5.0   # eV
CDFT_LUMO_MAX = 20.0   # eV
CDFT_GAP_MIN = 0.5     # eV (minimum reasonable HOMO-LUMO gap)
CDFT_GAP_MAX = 15.0    # eV (maximum reasonable HOMO-LUMO gap)
CDFT_OMEGA_MIN = 0.0   # eV
CDFT_OMEGA_MAX = 20.0  # eV


class Step2CyclizationExtractor(BaseExtractor):
    """Extract Step2 [5+7] cyclization features from S3 artifacts."""

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize Step2 extractor.

        Args:
            config: Configuration dict with cdft and gedt settings
        """
        super().__init__()
        self.config = config or {}
        self.cdft_config = self.config.get('stage2', {}).get('cdft', {})
        self.gedt_config = self.config.get('stage2', {}).get('gedt', {})
        self.warnings = []

    def get_plugin_name(self) -> str:
        return "step2_cyclization"

    def get_required_inputs(self) -> List[str]:
        """Return list of required context inputs.

        V6.2: All inputs are optional to allow graceful degradation.
        Extractors emit NaN values and warning codes when files are missing.
        """
        return []

    def extract(self, context) -> Dict[str, Any]:
        """Extract Step2 cyclization features.

        Args:
            context: FeatureContext with S3 path handles

        Returns:
            Dictionary of Step2 features with s2_* prefix
        """
        trace = context.get_plugin_trace(self.get_plugin_name())

        features: Dict[str, Any] = {}
        warnings: List[str] = []

        # V6.2 policy: Step2 plugin is extract-only.
        # De-dup policy: do NOT emit TS geometry/validity keys here.

        # CDFT features
        cdft_features, cdft_warnings = self._extract_cdft_features(context)
        features.update(cdft_features)
        warnings.extend(cdft_warnings)

        # GEDT features
        gedt_features, gedt_warnings = self._extract_gedt_features(context)
        features.update(gedt_features)
        warnings.extend(gedt_warnings)

        trace.warnings.extend(warnings)
        return features

    def _extract_kinetic_features(self, context) -> Dict[str, Any]:
        """Extract activation barrier and reaction energy features."""
        # Deprecated: Step2 does not own thermo/QC; kept for API compatibility.
        return {}

    def _extract_ts_geometry(self, context) -> Tuple[Dict[str, Any], List[str]]:
        """Extract TS forming bond distances and asynchronicity with validation.
        
        P0 Enhancement (V6.2): Strict validation of forming_bonds and provenance recording.
        """
        warnings: List[str] = []
        features: Dict[str, Any] = {
            's2_d_forming_1': float('nan'),
            's2_d_forming_2': float('nan'),
            's2_asynch': float('nan'),
            's2_forming_bonds_source': 'NONE',
            's2_ts_geom_source': 'NONE'
        }

        ts_xyz = context.ts_xyz
        forming_bonds = context.forming_bonds

        if not ts_xyz or not Path(ts_xyz).exists():
            return features, warnings

        features['s2_ts_geom_source'] = Path(ts_xyz).name

        try:
            coords, symbols = read_xyz(ts_xyz)
            n_atoms = len(symbols)
            
            is_valid, fb_errors = validate_forming_bonds(forming_bonds, natoms=n_atoms)

            if not is_valid:
                warnings.append(W_FORMING_BONDS_INVALID)
                return features, warnings

            features['s2_forming_bonds_source'] = 'context'

            distance_matrix = GeometryUtils.compute_distance_matrix(coords)
            bond_distances = [distance_matrix[bond[0], bond[1]].item() for bond in forming_bonds]

            features['s2_d_forming_1'] = bond_distances[0] if len(bond_distances) > 0 else float('nan')
            features['s2_d_forming_2'] = bond_distances[1] if len(bond_distances) > 1 else float('nan')

            if len(bond_distances) >= 2:
                features['s2_asynch'] = abs(bond_distances[0] - bond_distances[1])
            else:
                features['s2_asynch'] = float('nan')

        except Exception:
            pass

        return features, warnings

    def _extract_cdft_features(self, context) -> Tuple[Dict[str, Any], List[str]]:
        """Extract CDFT indices from fchk file with enhanced sanity checks.

        P1 Enhancement: Unit locking and rigorous range validation.
        All CDFT values are in eV (locked by fchk_reader.py HARTREE_TO_EV).
        """
        features = {}
        warnings = []

        ts_fchk = context.s3_ts_fchk
        intermediate_fchk = context.s3_intermediate_fchk

        fchk_path = intermediate_fchk or ts_fchk

        if fchk_path and Path(fchk_path).exists():
            try:
                cdft_indices = read_fchk_cdft_indices(fchk_path)

                features['s2_eps_homo'] = cdft_indices.get('eps_homo')
                features['s2_eps_lumo'] = cdft_indices.get('eps_lumo')
                features['s2_mu'] = cdft_indices.get('mu')
                features['s2_eta'] = cdft_indices.get('eta')
                features['s2_omega'] = cdft_indices.get('omega')

                # Sanity checks
                if features['s2_eps_homo'] is None or features['s2_eps_lumo'] is None:
                    warnings.append(W_CDFT_MISSING_ORBITALS)
                else:
                    # Basic orbital order check
                    if features['s2_eps_lumo'] <= features['s2_eps_homo']:
                        warnings.append(W_CDFT_ORBITAL_ORDER)

                    # P1 Enhancement: Range validation for all CDFT indices
                    # HOMO range check
                    if not (CDFT_HOMO_MIN <= features['s2_eps_homo'] <= CDFT_HOMO_MAX):
                        warnings.append(W_CDFT_HOMO_RANGE)

                    # LUMO range check
                    if not (CDFT_LUMO_MIN <= features['s2_eps_lumo'] <= CDFT_LUMO_MAX):
                        warnings.append(W_CDFT_LUMO_RANGE)

                    # Gap (eta) range check
                    if features['s2_eta'] is not None:
                        if not (CDFT_GAP_MIN <= features['s2_eta'] <= CDFT_GAP_MAX):
                            warnings.append(W_CDFT_GAP_RANGE)

                    # Omega (electrophilicity) range check
                    if features['s2_omega'] is not None:
                        if not (CDFT_OMEGA_MIN <= features['s2_omega'] <= CDFT_OMEGA_MAX):
                            warnings.append(W_CDFT_OMEGA_RANGE)

            except Exception as e:
                features['s2_eps_homo'] = float('nan')
                features['s2_eps_lumo'] = float('nan')
                features['s2_mu'] = float('nan')
                features['s2_eta'] = float('nan')
                features['s2_omega'] = float('nan')
                warnings.append(W_CDFT_MISSING_ORBITALS)
        else:
            features['s2_eps_homo'] = float('nan')
            features['s2_eps_lumo'] = float('nan')
            features['s2_mu'] = float('nan')
            features['s2_eta'] = float('nan')
            features['s2_omega'] = float('nan')
            warnings.append(W_CDFT_MISSING_ORBITALS)

        return features, warnings

    def _extract_gedt_features(self, context) -> Tuple[Dict[str, Any], List[str]]:
        """Extract GEDT from TS charges and fragment analysis."""
        features = {}
        warnings = []

        ts_fchk = context.s3_ts_fchk
        ts_xyz = context.ts_xyz
        forming_bonds = context.forming_bonds

        if not forming_bonds:
            features['s2_gedt_value'] = float('nan')
            features['s2_gedt_charge_type'] = 'NONE'
            features['s2_gedt_sign_convention'] = self.gedt_config.get('sign_convention', 'default')
            return features, warnings

        # Read coordinates for fragment cut
        coords = None
        atom_symbols = []
        if ts_xyz and Path(ts_xyz).exists():
            try:
                coords, atom_symbols = read_xyz(ts_xyz)
                coords = np.array(coords) if coords is not None else None
            except Exception:
                pass

        # Validate forming_bonds indices when we have TS atom count.
        if coords is not None and atom_symbols:
            is_valid, _fb_errors = validate_forming_bonds(forming_bonds, natoms=len(atom_symbols))
            if not is_valid:
                features['s2_gedt_value'] = float('nan')
                features['s2_gedt_charge_type'] = 'NONE'
                features['s2_gedt_sign_convention'] = self.gedt_config.get('sign_convention', 'default')
                warnings.append(W_FORMING_BONDS_INVALID)
                return features, warnings

        # Read charges with priority
        charge_priority = self.gedt_config.get('charge_priority', ['NBO', 'CM5', 'MULLIKEN'])
        charges, charge_type = read_charges_with_priority(
            fchk_path=ts_fchk,
            charge_priority=charge_priority
        )

        if charges is None or coords is None:
            features['s2_gedt_value'] = float('nan')
            features['s2_gedt_charge_type'] = charge_type if charges else 'NONE'
            features['s2_gedt_sign_convention'] = self.gedt_config.get('sign_convention', 'default')
            warnings.append(W_GEDT_NO_CHARGES)
            return features, warnings

        # Perform fragment cut and GEDT calculation
        try:
            gedt_result = cut_along_forming_bonds(
                coordinates=coords,
                forming_bonds=forming_bonds,
                charges=charges,
                symbols=atom_symbols if atom_symbols else None,  # Pass symbols for deterministic labeling
                config=self.gedt_config
            )

            features['s2_gedt_value'] = gedt_result.get('gedt_value', float('nan'))
            features['s2_gedt_charge_type'] = charge_type
            features['s2_gedt_sign_convention'] = gedt_result.get('gedt_sign_convention', 'default')
            features['s2_gedt_fragment_labeling'] = gedt_result.get('gedt_fragment_labeling', 'unknown')
            features['s2_q_fragment_dipole'] = gedt_result.get('q_fragment_dipole', float('nan'))
            features['s2_q_fragment_dipolarophile'] = gedt_result.get('q_fragment_dipolarophile', float('nan'))

            if gedt_result.get('fragment_a') is None:
                warnings.append(W_GEDT_FRAGMENT_CUT_FAILED)

        except Exception as e:
            features['s2_gedt_value'] = float('nan')
            features['s2_gedt_charge_type'] = charge_type
            features['s2_gedt_sign_convention'] = self.gedt_config.get('sign_convention', 'default')
            features['s2_gedt_fragment_labeling'] = 'error'
            warnings.append(W_GEDT_NO_CHARGES)

        return features, warnings

    def _extract_ts_validity(self, context) -> Tuple[Dict[str, Any], List[str]]:
        """Validate TS (imaginary frequency check)."""
        features = {}
        warnings = []

        ts_log = context.s3_ts_log

        n_imag = 0
        imag_freq_cm1 = float('nan')
        validity_flag = "ok"

        if ts_log and Path(ts_log).exists():
            try:
                with open(ts_log, 'r', encoding='utf-8') as f:
                    content = f.read()

                # Look for frequency section
                freq_match = re.search(
                    r'Frequencies\s*--\s*([-\d.]+)\s*([-\d.]+)?\s*([-\d.]+)?',
                    content
                )
                if freq_match:
                    freqs = [float(f) for f in freq_match.groups() if f]
                    imag_freqs = [f for f in freqs if f < 0]
                    n_imag = len(imag_freqs)

                    if n_imag > 0:
                        imag_freq_cm1 = abs(imag_freqs[0])

                    # Check for only one imaginary frequency
                    if n_imag != 1:
                        warnings.append(W_TS_IMAG_COUNT_NOT_ONE)
                        validity_flag = "warn"
                else:
                    warnings.append(W_TS_MISSING_LOG)
                    validity_flag = "warn"

            except Exception as e:
                warnings.append(W_TS_MISSING_LOG)
                validity_flag = "warn"
        else:
            warnings.append(W_TS_MISSING_LOG)
            validity_flag = "warn"

        features['s2_n_imag_freq'] = n_imag
        features['s2_imag_freq_cm1'] = imag_freq_cm1
        features['s2_ts_validity_flag'] = validity_flag

        return features, warnings


register_extractor(Step2CyclizationExtractor())
