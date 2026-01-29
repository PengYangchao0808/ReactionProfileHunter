"""
Step 4: QC Checks Extractor (V4.2 Phase A)
================================================

QC validation features extractor.

Author: QC Descriptors Team
Date: 2026-01-18
"""

from typing import Dict, Any, List

from .base import BaseExtractor, register_extractor


class QCChecksExtractor(BaseExtractor):
    """Extract QC validation features.

    Features (qc.* prefix):
    - has_gibbs: 1 if g_* present in sp_report, else 0
    - used_fallback_electronic: 1 if Gibbs unavailable, else 0
    - sp_report_validated: 1 if sp_report.validate() succeeded, else 0
    - forming_bonds_valid: 1 if forming_bonds provided and valid, else 0
    - warnings_count: Total warnings across all plugins (aggregated later)

    This plugin does not generate warnings_count directly (it's aggregated globally).
    """

    def get_plugin_name(self) -> str:
        return "qc_checks"

    def get_required_inputs(self) -> List[str]:
        return ["sp_report", "forming_bonds"]

    def extract(self, context) -> Dict[str, Any]:
        sp_report = context.sp_report
        forming_bonds = context.forming_bonds

        features = {}

        g_ts = getattr(sp_report, "g_ts", None) if sp_report else None
        g_reactant = getattr(sp_report, "g_reactant", None) if sp_report else None
        g_product = getattr(sp_report, "g_product", None) if sp_report else None

        features["qc.has_gibbs"] = 1 if (
            g_ts is not None and g_reactant is not None and g_product is not None
        ) else 0

        used_fallback_electronic = 0
        if features["qc.has_gibbs"] == 0:
            used_fallback_electronic = 1
        features["qc.used_fallback_electronic"] = used_fallback_electronic

        sp_report_validated = 0
        if sp_report is not None and hasattr(sp_report, "validate"):
            try:
                sp_report_validated = 1 if sp_report.validate() else 0
            except Exception:
                sp_report_validated = 0
        features["qc.sp_report_validated"] = sp_report_validated

        features["qc.forming_bonds_valid"] = 1 if forming_bonds else 0

        features["qc.warnings_count"] = 0

        return features


register_extractor(QCChecksExtractor())
