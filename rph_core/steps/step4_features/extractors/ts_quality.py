"""
Step 4: TS Quality Extractor (V6.1 Phase A)
=====================================

Transition state quality features extractor.

Extracts:
- n_imag: Number of imaginary frequencies
- imag1_cm1_abs: Absolute value of most negative frequency (cm^-1)
- dipole_debye: Dipole moment magnitude (Debye)

Author: QC Descriptors Team
Date: 2026-01-27
"""

from typing import Dict, Any, List
import numpy as np
import pathlib

from .base import BaseExtractor, register_extractor


class TSQualityExtractor(BaseExtractor):
    """Extract TS quality features from Gaussian or ORCA output.

    Features (ts.* prefix):
    - ts.n_imag: Number of imaginary frequencies
    - ts.imag1_cm1_abs: Absolute value of most negative frequency (cm^-1)
    - ts.dipole_debye: Dipole moment magnitude (Debye)

    Parsing sources (in order of preference):
    1. Gaussian .log file (parsed via log_parser)
    2. ORCA .out file (parsed via log_parser)
    """

    def get_plugin_name(self) -> str:
        """Return plugin name for feature prefix."""
        return "ts_quality"

    def get_required_inputs(self) -> List[str]:
        """Return list of required context inputs."""
        return ["ts_log", "ts_orca_out", "ts_fchk"]

    def extract(self, context) -> Dict[str, Any]:
        """Extract TS quality features from context.

        Args:
            context: FeatureContext with TS output paths

        Returns:
            Dictionary with ts.* features or empty dict if parsing fails
        """
        features = {}

        # Try to parse frequencies from Gaussian log
        ts_log = context.ts_log
        if ts_log is not None and ts_log.exists():
            try:
                from ..log_parser import _extract_frequencies
                frequencies = _extract_frequencies(ts_log)
            except Exception:
                frequencies = None
        else:
            frequencies = None

        # Try to parse frequencies from ORCA .out
        if frequencies is None and context.ts_orca_out is not None and context.ts_orca_out.exists():
            try:
                from ..log_parser import _parse_ts_output
                ts_data = _parse_ts_output(context.ts_orca_out)
                if ts_data is not None and hasattr(ts_data, 'frequencies'):
                    frequencies = ts_data.frequencies
            except Exception:
                frequencies = None

        if frequencies is None or len(frequencies) == 0:
            # No frequencies available, return empty features
            return features

        # Compute n_imag (count of negative frequencies)
        imaginary_freqs = [f for f in frequencies if f < 0]
        features["ts.n_imag"] = len(imaginary_freqs)

        # Compute imag1_cm1_abs (absolute value of most negative frequency)
        if len(imaginary_freqs) > 0:
            most_negative = min(imaginary_freqs)
            features["ts.imag1_cm1_abs"] = abs(most_negative)
        else:
            features["ts.imag1_cm1_abs"] = np.nan

        # Try to parse dipole moment
        dipole_debye = None

        # Try to read from Gaussian .log
        if ts_log is not None and ts_log.exists():
            try:
                from ..log_parser import _extract_dipole_moment
                dipole_vector = _extract_dipole_moment(ts_log)
                if dipole_vector is not None and len(dipole_vector) >= 3:
                    dipole_debye = np.linalg.norm(dipole_vector)
            except Exception:
                pass

        # Try to read from ORCA .out (log_parser doesn't have ORCA dipole parsing yet)
        if dipole_debye is None and context.ts_orca_out is not None and context.ts_orca_out.exists():
            try:
                from ..log_parser import _parse_ts_output
                ts_data = _parse_ts_output(context.ts_orca_out)
                if ts_data is not None and hasattr(ts_data, 'dipole_moment'):
                    dipole_vector = ts_data.dipole_moment
                    if dipole_vector is not None and len(dipole_vector) >= 3:
                        dipole_debye = np.linalg.norm(dipole_vector)
            except Exception:
                pass

        # Set dipole feature
        if dipole_debye is not None:
            features["ts.dipole_debye"] = dipole_debye
        else:
            features["ts.dipole_debye"] = np.nan

        return features


# Register extractor
register_extractor(TSQualityExtractor())
