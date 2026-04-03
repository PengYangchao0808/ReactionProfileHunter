"""
Step 4: Multiwfn Features Extractor (V6.2 P2)
=============================================

Extract Multiwfn-derived features:
- Fukui functions (f+, f-, f0) for forming bond atoms
- Dual descriptor values
- QTAIM BCP properties (rho, laplacian)
- NCI statistics

P2 Feature: Tier-1 Multiwfn analysis with caching and fail-open.

Author: RPH Team
Date: 2026-02-02
"""

import logging
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

from .base import BaseExtractor, register_extractor
from rph_core.utils.multiwfn_runner import (
    MultiwfnRunner,
    W_MULTIWFN_DISABLED,
    W_MULTIWFN_FAILED,
    W_MULTIWFN_CACHE_READ_FAILED,
    W_MULTIWFN_TIMEOUT,
    W_MULTIWFN_INVALID_OUTPUT
)


logger = logging.getLogger(__name__)


class MultiwfnFeaturesExtractor(BaseExtractor):
    """Extract Multiwfn-derived features from fchk/wfn files."""

    # Warning codes
    W_MW_MISSING_INPUT = "W_MW_MISSING_INPUT"
    W_MW_NO_FCHK = "W_MW_NO_FCHK"
    W_MW_NO_FORMING_BONDS = "W_MW_NO_FORMING_BONDS"
    W_MW_ATOM_NOT_FOUND = "W_MW_ATOM_NOT_FOUND"

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize Multiwfn extractor.

        Args:
            config: Configuration dict with multiwfn settings
        """
        super().__init__()
        self.config = config or {}
        self.mw_config = self.config.get('multiwfn', {})

        # Extract config parameters
        self.enabled = self.mw_config.get('enabled', False)
        self.multiwfn_path = self.mw_config.get('multiwfn_path', 'Multiwfn')
        self.cache_dir = Path(self.mw_config.get('cache_dir', '.cache/step4_multiwfn'))
        self.timeout_sec = self.mw_config.get('timeout_sec', 120)
        self.modules = self.mw_config.get('modules', {
            'fukui': True,
            'dual_descriptor': True,
            'qtaim_bcp': False,
            'nci': False
        })

        # Initialize runner if enabled
        self.runner: Optional[MultiwfnRunner] = None
        if self.enabled:
            try:
                self.runner = MultiwfnRunner(
                    multiwfn_path=self.multiwfn_path,
                    cache_dir=self.cache_dir,
                    timeout_sec=self.timeout_sec,
                    enabled_modules=self.modules
                )
            except Exception as e:
                logger.warning(f"Failed to initialize Multiwfn runner: {e}")
                self.runner = None

    def get_plugin_name(self) -> str:
        return "multiwfn_features"

    def get_required_inputs(self) -> List[str]:
        return [
            "s3_ts_fchk",
            "forming_bonds",
        ]

    def extract(self, context) -> Dict[str, Any]:
        """Extract Multiwfn features.

        Args:
            context: FeatureContext with S3 path handles

        Returns:
            Dictionary of Multiwfn features with mw_* prefix
        """
        features = {}
        warnings = []

        # Check if Multiwfn is enabled
        if not self.enabled:
            features['mw_status'] = 'disabled'
            features['mw_missing_reason'] = W_MULTIWFN_DISABLED
            return features

        if self.runner is None:
            features['mw_status'] = 'failed'
            features['mw_missing_reason'] = 'runner_init_failed'
            return features

        # Get required inputs
        ts_fchk = context.s3_ts_fchk
        forming_bonds = context.forming_bonds

        # Validate inputs
        if ts_fchk is None or not Path(ts_fchk).exists():
            features['mw_status'] = 'failed'
            features['mw_missing_reason'] = W_MW_NO_FCHK
            return features

        if not forming_bonds or len(forming_bonds) == 0:
            features['mw_status'] = 'failed'
            features['mw_missing_reason'] = W_MW_NO_FORMING_BONDS
            return features

        try:
            # Get atom indices for forming bonds (convert 0-based to 1-based)
            atom_a = forming_bonds[0][0] + 1
            atom_b = forming_bonds[0][1] + 1

            # Run complete analysis
            result = self.runner.run_complete_analysis(
                Path(ts_fchk), atom_a, atom_b
            )

            # Build feature dict
            features = self._result_to_features(result, atom_a, atom_b)

            # Add warnings
            if result.warnings:
                features['mw_warnings'] = result.warnings
                features['mw_warnings_count'] = len(result.warnings)

        except Exception as e:
            logger.error(f"Multiwfn extraction failed: {e}", exc_info=True)
            features['mw_status'] = 'failed'
            features['mw_missing_reason'] = W_MULTIWFN_FAILED
            features['mw_error_message'] = str(e)

        return features

    def _result_to_features(
        self,
        result,
        atom_a: int,
        atom_b: int
    ) -> Dict[str, Any]:
        """Convert MultiwfnResult to feature dict.

        Args:
            result: MultiwfnResult object
            atom_a: First atom index (1-based)
            atom_b: Second atom index (1-based)

        Returns:
            Feature dictionary
        """
        features = {}

        # Status
        if result.success:
            features['mw_status'] = 'ok'
        else:
            features['mw_status'] = 'failed'
            if result.warnings:
                features['mw_missing_reason'] = result.warnings[0]
            elif result.error_message:
                features['mw_missing_reason'] = result.error_message
            else:
                features['mw_missing_reason'] = 'unknown'

        # Fukui functions
        if result.fukui_fplus_a is not None:
            features['mw_fukui_fplus_atomA'] = result.fukui_fplus_a
        if result.fukui_fplus_b is not None:
            features['mw_fukui_fplus_atomB'] = result.fukui_fplus_b
        if result.fukui_fminus_a is not None:
            features['mw_fukui_fminus_atomA'] = result.fukui_fminus_a
        if result.fukui_fminus_b is not None:
            features['mw_fukui_fminus_atomB'] = result.fukui_fminus_b
        if result.fukui_f0_a is not None:
            features['mw_fukui_f0_atomA'] = result.fukui_f0_a
        if result.fukui_f0_b is not None:
            features['mw_fukui_f0_atomB'] = result.fukui_f0_b

        # Dual descriptor
        if result.dual_descriptor_a is not None:
            features['mw_dual_descriptor_atomA'] = result.dual_descriptor_a
        if result.dual_descriptor_b is not None:
            features['mw_dual_descriptor_atomB'] = result.dual_descriptor_b

        # QTAIM BCP
        if result.rho_bcp_forming1 is not None:
            features['mw_rho_bcp_forming1'] = result.rho_bcp_forming1
        if result.laplacian_bcp_forming1 is not None:
            features['mw_laplacian_bcp_forming1'] = result.laplacian_bcp_forming1

        # Atom indices for reference
        features['mw_atomA'] = atom_a
        features['mw_atomB'] = atom_b

        # Cache hit indicator
        if W_MULTIWFN_CACHE_READ_FAILED in (result.warnings or []):
            features['mw_cache_hit'] = True

        return features


register_extractor(MultiwfnFeaturesExtractor())
