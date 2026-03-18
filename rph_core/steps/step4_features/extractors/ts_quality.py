"""
Step 4: TS Quality Extractor (V6.2 Phase A)
=====================================

Transition state quality features extractor.

Extracts:
- n_imag: Number of imaginary frequencies
- imag1_cm1_abs: Absolute value of most negative frequency (cm^-1)
- dipole_debye: Dipole moment magnitude (Debye)

V6.2: All inputs are optional - extractor degrades gracefully with NaN values
and warning codes when files are missing.

Author: QC Descriptors Team
Date: 2026-02-03
"""

from typing import Dict, Any, List, Optional
import re
import numpy as np
import pathlib

from .base import BaseExtractor, register_extractor


class TSQualityExtractor(BaseExtractor):
    """Extract TS quality features from Gaussian or ORCA output.

    Features (ts.* prefix):
    - ts.n_imag: Number of imaginary frequencies (0 if no frequencies found)
    - ts.imag1_cm1_abs: Absolute value of most negative frequency (cm^-1), NaN if unavailable
    - ts.dipole_debye: Dipole moment magnitude (Debye), NaN if unavailable

    V6.2: All inputs are optional. If files are missing:
    - Emit keys with NaN values
    - Record appropriate warning codes
    - Status remains OK (not SKIPPED)

    Parsing sources (in order of preference):
    1. Gaussian .log file (parsed via log_parser)
    2. ORCA .out file (parsed via log_parser)
    3. Fallback to s3_ts_log if available
    """

    def get_plugin_name(self) -> str:
        """Return plugin name for feature prefix."""
        return "ts_quality"

    def get_required_inputs(self) -> List[str]:
        """Return list of required context inputs.

        V6.2: All inputs are optional to allow graceful degradation.
        """
        return []

    def extract(self, context) -> Dict[str, Any]:
        """Extract TS quality features from context.

        Args:
            context: FeatureContext with TS output paths

        Returns:
            Dictionary with ts.* features (NaN for missing data)
        """
        features = {}
        warnings = []
        trace = context.get_plugin_trace(self.get_plugin_name())

        # Try to parse frequencies from Gaussian log
        ts_log = context.ts_log
        frequencies = None

        if ts_log is not None and pathlib.Path(ts_log).exists():
            try:
                from ..log_parser import GaussianLogParser
                with open(ts_log, 'r') as f:
                    content = f.read()
                frequencies = GaussianLogParser._extract_frequencies(content)
            except Exception as e:
                frequencies = None
        else:
            # Try fallback to s3_ts_log
            s3_ts_log = getattr(context, 's3_ts_log', None)
            if s3_ts_log is not None and pathlib.Path(s3_ts_log).exists():
                try:
                    from ..log_parser import GaussianLogParser
                    with open(s3_ts_log, 'r') as f:
                        content = f.read()
                    frequencies = GaussianLogParser._extract_frequencies(content)
                except Exception:
                    frequencies = None

        if frequencies is None or len(frequencies) == 0:
            # No frequencies available - emit NaN with warning
            warnings.append("W_TS_MISSING_LOG")
            features["ts.n_imag"] = 0
            features["ts.imag1_cm1_abs"] = np.nan
        else:
            # Compute n_imag (count of negative frequencies)
            imaginary_freqs = [f for f in frequencies if f < 0]
            n_imag = len(imaginary_freqs)
            features["ts.n_imag"] = n_imag

            # Record warning if n_imag != 1
            if n_imag != 1:
                warnings.append("W_TS_IMAG_COUNT_NOT_ONE")

            # Compute imag1_cm1_abs (absolute value of most negative frequency)
            if len(imaginary_freqs) > 0:
                most_negative = min(imaginary_freqs)
                features["ts.imag1_cm1_abs"] = abs(most_negative)
            else:
                features["ts.imag1_cm1_abs"] = np.nan

        # Try to parse dipole moment
        dipole_debye = None

        def _parse_fortran_float(token: str) -> Optional[float]:
            token = (token or "").strip()
            if not token:
                return None
            try:
                return float(token.replace("D", "E").replace("d", "e"))
            except ValueError:
                return None

        def _extract_gaussian_dipole_debye(content: str) -> Optional[float]:
            # Prefer archive Dipole= triplet when available.
            m = re.search(
                r"Dipole=\s*([+\-\d\.DEe]+)\s*,\s*([+\-\d\.DEe]+)\s*,\s*([+\-\d\.DEe]+)",
                content,
            )
            if m:
                x = _parse_fortran_float(m.group(1))
                y = _parse_fortran_float(m.group(2))
                z = _parse_fortran_float(m.group(3))
                if x is not None and y is not None and z is not None:
                    return float(np.linalg.norm([x, y, z]))

            # Fallback: field-independent basis block.
            m = re.search(
                r"Dipole moment \(field-independent basis, Debye\):\s*\n\s*X=\s*([+\-\d\.]+)\s+Y=\s*([+\-\d\.]+)\s+Z=\s*([+\-\d\.]+)(?:\s+Tot=\s*([+\-\d\.]+))?",
                content,
            )
            if m:
                if m.group(4) is not None:
                    try:
                        return float(m.group(4))
                    except ValueError:
                        pass
                try:
                    x = float(m.group(1))
                    y = float(m.group(2))
                    z = float(m.group(3))
                    return float(np.linalg.norm([x, y, z]))
                except ValueError:
                    return None

            return None

        # Try to read from Gaussian .log
        if ts_log is not None and pathlib.Path(ts_log).exists():
            try:
                from ..log_parser import GaussianLogParser
                with open(ts_log, 'r') as f:
                    content = f.read()
                dipole_debye = _extract_gaussian_dipole_debye(content)
            except Exception:
                pass

        # ORCA dipole parsing not implemented - leave as NaN

        # Try to read from ORCA .out
        if dipole_debye is None and context.ts_orca_out is not None and pathlib.Path(context.ts_orca_out).exists():
            try:
                with open(context.ts_orca_out, 'r') as f:
                    content = f.read()
                m = re.search(
                    r"Total Dipole Moment\s*:\s*([+\-\d\.]+)\s+([+\-\d\.]+)\s+([+\-\d\.]+)",
                    content,
                )
                if m:
                    dipole_debye = float(np.linalg.norm([float(m.group(1)), float(m.group(2)), float(m.group(3))]))
            except Exception:
                pass

        # Set dipole feature
        if dipole_debye is not None:
            features["ts.dipole_debye"] = dipole_debye
        else:
            features["ts.dipole_debye"] = np.nan

        # Record warnings in trace
        trace.warnings.extend(warnings)

        return features


# Register extractor
register_extractor(TSQualityExtractor())
