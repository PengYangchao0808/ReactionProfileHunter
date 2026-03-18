"""
Step 4: Feature Extractors (V4.2 Phase A)
========================================

Base extractor class and plugin registry.

Author: QC Descriptors Team
Date: 2026-01-18
"""

from .base import (
    BaseExtractor,
    EXTRACTORS,
    register_extractor,
    get_extractor,
    list_extractors,
)

from . import thermo, geometry, qc_checks, interaction_analysis, nics, nbo_e2, ts_quality, asm_enrichment
from .fmo_cdft_dipolar import FmoCdftDipolarParser

# V6.2 new extractors
from . import step1_activation, step2_cyclization

# P2 Enhancement: Multiwfn features
from . import multiwfn_features

register_extractor(FmoCdftDipolarParser())
