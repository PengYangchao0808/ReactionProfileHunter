"""
FMO/CDFT Dipolar Parser
=========================

Lightweight FMO/CDFT parser for dipolar intermediate outputs from S3.

Author: QCcalc Team
Date: 2026-01-27
"""

import logging
import re
import hashlib
from pathlib import Path
from typing import Optional, Dict, Any

from .base import BaseExtractor
from ..status import FeatureResultStatus

logger = logging.getLogger(__name__)


class FmoCdftDipolarParser(BaseExtractor):
    """
    FMO/CDFT dipolar parser that reads HOMO/LUMO from Gaussian/ORCA outputs.

    Parses dipolar intermediate files to extract:
    - HOMO energy (eV)
    - LUMO energy (eV)
    - Gap (eV)
    - Electrophilicity (omega, eV) = (HOMO² + LUMO²)/2
    """

    HARTREE_TO_EV = 27.2114

    def get_plugin_name(self) -> str:
        """Return unique plugin name."""
        return "fmo_cdft_dipolar"

    def can_submit_jobs(self) -> bool:
        """Return whether this extractor can submit QC jobs."""
        return False

    def get_required_inputs(self) -> list[str]:
        return []

    def get_required_inputs_for_context(self, context) -> list[str]:
        """Dynamic hook to determine required inputs based on available context."""
        required = []

        if hasattr(context, 's3_dir'):
            required.append('s3_dir')

        if hasattr(context, 'artifacts_index'):
            required.append('artifacts_index')

        if not required:
            return super().get_required_inputs_for_context(context)

        return required

    def extract(self, context) -> Dict[str, Any]:
        """
        Extract dipolar FMO features from context.

        Args:
            context: FeatureContext with all inputs

        Returns:
            Dictionary with fmo_cdft_dipolar_* features
        """
        trace = context.get_plugin_trace(self.get_plugin_name())
        trace.missing_fields = []
        trace.missing_paths = []
        trace.errors = []
        trace.warnings = []
        extracted_features = {}

        s3_dir = getattr(context, 's3_dir', None)
        artifacts_index = getattr(context, 'artifacts_index', None)

        if s3_dir is None:
            trace.missing_fields.append('s3_dir')
            trace.status = FeatureResultStatus.SKIPPED
            trace.errors.append("s3_dir not provided")
            return trace._extracted_features

        if artifacts_index is None:
            logger.warning("artifacts_index not provided, attempting fallback scan")
            dipolar_output = self._scan_for_dipolar_output(s3_dir)
        else:
            dipolar_output = self._resolve_dipolar_from_index(artifacts_index, s3_dir)

        if dipolar_output is None:
            trace.missing_fields.append('dipolar_output')
            trace.status = FeatureResultStatus.SKIPPED
            trace.warnings.append("No dipolar output found")
            extracted_features['fmo_cdft_dipolar.status'] = "skipped"
            extracted_features['fmo_cdft_dipolar.missing_reason'] = "dipolar_output_not_found"
            return trace._extracted_features

        try:
            homo_ev, lumo_ev, gap_ev, omega_ev = self._parse_dipolar_output(dipolar_output)

            if homo_ev is not None and lumo_ev is not None:
                trace._extracted_features['fmo_cdft_dipolar.homo_ev'] = homo_ev
                trace._extracted_features['fmo_cdft_dipolar.lumo_ev'] = lumo_ev
                trace._extracted_features['fmo_cdft_dipolar.gap_ev'] = gap_ev
                trace._extracted_features['fmo_cdft_dipolar.omega_ev'] = omega_ev
                trace._extracted_features['fmo_cdft_dipolar.status'] = "ok"
                trace._extracted_features['fmo_cdft_dipolar.missing_reason'] = None

                gap_valid = gap_ev > 0 if gap_ev is not None else False
                omega_valid = omega_ev > 0 if omega_ev is not None else False

                if not gap_valid:
                    trace.warnings.append("Gap is not positive")
                    trace._extracted_features['fmo_cdft_dipolar.gap_ev_is_invalid'] = True
                    trace._extracted_features['fmo_cdft_dipolar.gap_is_invalid_reason'] = "gap_not_positive"

                if not omega_valid:
                    trace.warnings.append("Omega is not positive")
                    trace._extracted_features['fmo_cdft_dipolar.omega_ev_is_invalid'] = True
                    trace._extracted_features['fmo_cdft_dipolar.omega_is_invalid_reason'] = "omega_not_positive"
            else:
                trace.status = FeatureResultStatus.SKIPPED
                trace.missing_fields.append('fmo_cdft_dipolar features')
                trace.errors.append("Failed to parse HOMO/LUMO")

        except Exception as e:
            logger.error(f"Failed to parse dipolar output: {e}", exc_info=True)
            trace.status = FeatureResultStatus.FAILED
            trace.errors.append(f"Exception: {str(e)}")

        return trace._extracted_features

    def _scan_for_dipolar_output(self, s3_dir: Path) -> Optional[Path]:
        """
        Scan S3 directory for dipolar output files.

        Args:
            s3_dir: S3 output directory

        Returns:
            Path to dipolar output, or None if not found
        """
        dipolar_patterns = ['*dipolar*.log', '*dipolar*.out']

        for pattern in dipolar_patterns:
            for output_path in sorted(s3_dir.glob(pattern)):
                return output_path

        return None

    def _resolve_dipolar_from_index(self, artifacts_index: Dict[str, Any], s3_dir: Path) -> Optional[Path]:
        """
        Resolve dipolar output path from artifacts_index.

        Args:
            artifacts_index: S3 artifacts_index.json contents
            s3_dir: S3 output directory

        Returns:
            Path to dipolar output, or None if not found
        """
        dipolar_info = artifacts_index.get('dipolar', {})

        if not dipolar_info:
            logger.warning("No dipolar info in artifacts_index")
            return None

        path_rel = dipolar_info.get('path_rel')
        if not path_rel:
            logger.warning("No dipolar path_rel in artifacts_index")
            return None

        dipolar_path = s3_dir / path_rel

        if dipolar_path.exists():
            sha256 = dipolar_info.get('sha256')
            if sha256:
                computed_sha = self._compute_file_hash(dipolar_path)
                if computed_sha != sha256:
                    logger.warning(f"SHA256 mismatch for {dipolar_path}: index={sha256}, file={computed_sha[:8]}...")
            return dipolar_path

        return None

    def _parse_dipolar_output(self, dipolar_path: Path) -> tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
        """
        Parse dipolar output file to extract HOMO/LUMO.

        Args:
            dipolar_path: Path to dipolar log/out file

        Returns:
            (homo_ev, lumo_ev, gap_ev, omega_ev)
        """
        suffix = dipolar_path.suffix.lower()

        content = dipolar_path.read_text()

        if suffix == '.log':
            homo_ev, lumo_ev = self._parse_gaussian_log(content)
        elif suffix == '.out':
            homo_ev, lumo_ev = self._parse_orca_out(content)
        else:
            logger.warning(f"Unsupported dipolar output format: {suffix}")
            return None, None, None, None

        if homo_ev is not None and lumo_ev is not None:
            gap_ev = lumo_ev - homo_ev
            omega_ev = (homo_ev**2 + lumo_ev**2) / 2

            return homo_ev, lumo_ev, gap_ev, omega_ev

        return None, None, None, None

    def _parse_gaussian_log(self, content: str) -> tuple[Optional[float], Optional[float]]:
        """
        Parse Gaussian .log file for HOMO/LUMO.

        Args:
            content: Gaussian log file content

        Returns:
            (homo_ev, lumo_ev)
        """
        H_TO_EV = 27.2114

        occ_matches = re.findall(r"Alpha\s+occ\.\s+eigenvalues\s+--\s+(.+)", content)
        vir_matches = re.findall(r"Alpha\s+virt\.\s+eigenvalues\s+--\s+(.+)", content)

        if not occ_matches or not vir_matches:
            return None, None

        last_occ = occ_matches[-1].split()
        homo = float(last_occ[-1]) * self.HARTREE_TO_EV

        first_virt = vir_matches[0].split()
        lumo = float(first_virt[0]) * self.HARTREE_TO_EV

        return homo, lumo

    def _parse_orca_out(self, content: str) -> tuple[Optional[float], Optional[float]]:
        """
        Parse ORCA .out file for HOMO/LUMO.

        Args:
            content: ORCA out file content

        Returns:
            (homo_ev, lumo_ev)
        """
        occ_pattern = r"ALPHA\s+EIGENVALUES\s*\(\s*(-?\d+\.\d+)\s*\)"

        occ_match = re.search(occ_pattern, content, re.IGNORECASE)

        if not occ_match:
            return None, None

        eigenvalues_str = occ_match.group(1)
        eigenvalues = list(map(float, eigenvalues_str.split()))

        if not eigenvalues:
            return None, None

        homo = max(eigenvalues)
        lumo = min([e for e in eigenvalues if e > 0])

        return homo * self.HARTREE_TO_EV, lumo * self.HARTREE_TO_EV

    @staticmethod
    def _compute_file_hash(file_path: Path) -> str:
        """
        Compute SHA256 hash of a file.

        Args:
            file_path: Path to file

        Returns:
            SHA256 hash hex string
        """
        sha256_hash = hashlib.sha256()
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                sha256_hash.update(chunk)
        return sha256_hash.hexdigest()
