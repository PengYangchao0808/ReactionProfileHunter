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
from rph_core.utils.constants import HARTREE_TO_EV

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
        # Only require s3_dir; artifacts_index is optional (fallback scan is supported).
        if hasattr(context, 's3_dir'):
            return ['s3_dir']
        return super().get_required_inputs_for_context(context)

    def extract(self, context) -> Dict[str, Any]:
        """
        Extract dipolar FMO features from context.

        Args:
            context: FeatureContext with all inputs

        Returns:
            Dictionary with fmo_cdft_dipolar_* features
        """
        trace = context.get_plugin_trace(self.get_plugin_name())
        features: Dict[str, Any] = {}

        s3_dir = getattr(context, 's3_dir', None)
        artifacts_index = getattr(context, 'artifacts_index', None)

        if s3_dir is None:
            # Should be gated by BaseExtractor.validate_inputs(), but keep a safe fallback.
            trace.warnings.append("W_FMO_MISSING_S3_DIR")
            features['fmo_cdft_dipolar.status'] = "skipped"
            features['fmo_cdft_dipolar.missing_reason'] = "s3_dir_missing"
            return features

        dipolar_output: Optional[Path] = None
        if isinstance(artifacts_index, dict):
            dipolar_output = self._resolve_dipolar_from_index(artifacts_index, s3_dir)
        elif artifacts_index is not None:
            logger.warning("artifacts_index is not a dict; attempting fallback scan")

        if dipolar_output is None:
            dipolar_output = self._scan_for_dipolar_output(s3_dir)

        if dipolar_output is None:
            trace.warnings.append("W_FMO_DIPOLAR_OUTPUT_NOT_FOUND")
            features['fmo_cdft_dipolar.status'] = "skipped"
            features['fmo_cdft_dipolar.missing_reason'] = "dipolar_output_not_found"
            return features

        try:
            homo_ev, lumo_ev, gap_ev, omega_ev = self._parse_dipolar_output(dipolar_output)
        except Exception as e:
            logger.error(f"Failed to parse dipolar output: {e}", exc_info=True)
            trace.errors.append(f"Exception: {str(e)}")
            features['fmo_cdft_dipolar.status'] = "failed"
            features['fmo_cdft_dipolar.missing_reason'] = "exception"
            return features

        if homo_ev is None or lumo_ev is None:
            trace.warnings.append("W_FMO_PARSE_HOMO_LUMO_FAILED")
            features['fmo_cdft_dipolar.status'] = "skipped"
            features['fmo_cdft_dipolar.missing_reason'] = "homo_lumo_unavailable"
            return features

        features['fmo_cdft_dipolar.homo_ev'] = homo_ev
        features['fmo_cdft_dipolar.lumo_ev'] = lumo_ev
        features['fmo_cdft_dipolar.gap_ev'] = gap_ev
        features['fmo_cdft_dipolar.omega_ev'] = omega_ev
        features['fmo_cdft_dipolar.status'] = "ok"
        features['fmo_cdft_dipolar.missing_reason'] = None

        gap_valid = gap_ev is not None and gap_ev > 0
        omega_valid = omega_ev is not None and omega_ev > 0

        if not gap_valid:
            trace.warnings.append("W_FMO_GAP_NOT_POSITIVE")
            features['fmo_cdft_dipolar.gap_ev_is_invalid'] = True
            features['fmo_cdft_dipolar.gap_is_invalid_reason'] = "gap_not_positive"

        if not omega_valid:
            trace.warnings.append("W_FMO_OMEGA_NOT_POSITIVE")
            features['fmo_cdft_dipolar.omega_ev_is_invalid'] = True
            features['fmo_cdft_dipolar.omega_is_invalid_reason'] = "omega_not_positive"

        return features

    def _scan_for_dipolar_output(self, s3_dir: Path) -> Optional[Path]:
        """
        Scan S3 directory for dipolar output files.

        Args:
            s3_dir: S3 output directory

        Returns:
            Path to dipolar output, or None if not found
        """
        # Prefer ORCA single-point outputs under L2_SP, and prefer TS side over reactant.
        candidates = []

        # Direct dipolar-named files (legacy)
        candidates.extend(s3_dir.rglob("*dipolar*.out"))
        candidates.extend(s3_dir.rglob("*dipolar*.log"))

        # Generic ORCA outputs (common in V6.x layouts)
        candidates.extend(s3_dir.rglob("*.out"))

        uniq = []
        seen = set()
        for p in candidates:
            if p in seen:
                continue
            seen.add(p)
            if p.is_file():
                uniq.append(p)

        def _score(path: Path) -> tuple[int, float, str]:
            parts = [x.lower() for x in path.parts]
            name = path.name.lower()
            score = 0
            if "l2_sp" in parts:
                score += 100
            if "ts_opt" in parts or "ts" in name:
                score += 50
            if "ts" in parts:
                score += 25
            if "reactant" in parts or "reactant" in name:
                score -= 10
            if "dipolar" in name:
                score += 30
            if name.endswith(".smd.out"):
                score -= 50
            try:
                mtime = path.stat().st_mtime
            except OSError:
                mtime = 0.0
            return (score, mtime, str(path))

        uniq.sort(key=_score, reverse=True)

        for p in uniq:
            suffix = p.suffix.lower()
            if suffix == ".out":
                try:
                    txt = p.read_text(encoding="utf-8", errors="ignore")
                except Exception:
                    continue
                if "ORBITAL ENERGIES" in txt:
                    return p
                # Still allow dipolar-named ORCA outputs even without the section.
                if "dipolar" in p.name.lower():
                    return p
            elif suffix == ".log" and "dipolar" in p.name.lower():
                return p

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
        Parse Gaussian .log file for HOMO/LUMO with enhanced robustness.

        V6.2 Enhancement: More flexible regex patterns and better error diagnostics.

        Args:
            content: Gaussian log file content

        Returns:
            (homo_ev, lumo_ev)
        """
        # More flexible patterns that handle spacing variations
        occ_pattern = r"Alpha\s+occ\.?\s+eigenvalues\s*[-—–:]+\s*(.+)"
        vir_pattern = r"Alpha\s+virt\.?\s+eigenvalues\s*[-—–:]+\s*(.+)"

        occ_matches = re.findall(occ_pattern, content, re.IGNORECASE)
        vir_matches = re.findall(vir_pattern, content, re.IGNORECASE)

        if not occ_matches:
            logger.debug("Gaussian occ eigenvalues not found in log")
            return None, None
        if not vir_matches:
            logger.debug("Gaussian virt eigenvalues not found in log")
            return None, None

        # Collect all occupied eigenvalues (handle multi-line)
        all_occ_values = []
        for line in occ_matches:
            values = line.split()
            for v in values:
                try:
                    all_occ_values.append(float(v))
                except ValueError:
                    continue

        # Collect all virtual eigenvalues
        all_vir_values = []
        for line in vir_matches:
            values = line.split()
            for v in values:
                try:
                    all_vir_values.append(float(v))
                except ValueError:
                    continue

        if not all_occ_values:
            logger.debug("No valid occ eigenvalues parsed from Gaussian log")
            return None, None
        if not all_vir_values:
            logger.debug("No valid virt eigenvalues parsed from Gaussian log")
            return None, None

        # HOMO = last occupied, LUMO = first virtual (both in Hartree, convert to eV)
        homo = all_occ_values[-1] * HARTREE_TO_EV
        lumo = all_vir_values[0] * HARTREE_TO_EV

        return homo, lumo

    def _parse_orca_out(self, content: str) -> tuple[Optional[float], Optional[float]]:
        """Parse ORCA .out file for HOMO/LUMO from ORBITAL ENERGIES section."""
        import re

        orbital_section = re.search(
            r'ORBITAL ENERGIES\s*-+\s*(.*?)(?:\n\s*-+\s*\n|\Z)',
            content,
            re.DOTALL | re.IGNORECASE
        )

        if not orbital_section:
            return None, None

        section_content = orbital_section.group(1)
        occ_values = []
        virt_values = []

        for line in section_content.split('\n'):
            line = line.strip()
            if not line or line.startswith('NO'):
                continue

            match = re.match(
                r'\s*(\d+)\s+([\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)',
                line
            )
            if match:
                occ = float(match.group(2))
                energy_ev = float(match.group(4))

                if occ > 0.5:
                    occ_values.append(energy_ev)
                else:
                    virt_values.append(energy_ev)

        if not occ_values or not virt_values:
            return None, None

        return max(occ_values), min(virt_values)

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
