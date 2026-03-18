"""
Step 4: ASM Enrichment Extractor (V6.2)
==========================================

ASM/DIAS features extractor from S3 enrichment contract.

Reads enrichment.json written by Step3 post-QC enrichment and derives:
- asm.distortion_total_kcal: Total distortion energy
- asm.interaction_kcal: Interaction energy (if ts_TS available)

Author: QC Descriptors Team
Date: 2026-01-27
"""

from typing import Dict, Any, List, Optional
import json
import pathlib

from .base import BaseExtractor, register_extractor
from rph_core.utils.constants import HARTREE_TO_KCAL


class ASMEnrichmentExtractor(BaseExtractor):
    """Extract ASM/DIAS features from S3 enrichment contract.

    Features (asm.* prefix):
    - distortion_total_kcal: Total distortion energy from fragment SPs
    - interaction_kcal: Interaction energy (if ts_TS available)

    Reads from <S3_DIR>/S3_PostQCEnrichment/enrichment.json.
    """

    def get_plugin_name(self) -> str:
        return "asm_enrichment"

    def get_required_inputs(self) -> List[str]:
        return ["s3_dir"]

    def extract(self, context) -> Dict[str, Any]:
        features = {}

        s3_dir = context.s3_dir
        if s3_dir is None:
            features["asm.distortion_total_kcal"] = None
            features["asm.interaction_kcal"] = None
            return features

        enrichment_dir = s3_dir / "S3_PostQCEnrichment"
        if not enrichment_dir.exists():
            features["asm.distortion_total_kcal"] = None
            features["asm.interaction_kcal"] = None
            return features

        enrichment_json = enrichment_dir / "enrichment.json"
        if not enrichment_json.exists():
            features["asm.distortion_total_kcal"] = None
            features["asm.interaction_kcal"] = None
            return features

        try:
            with open(enrichment_json, 'r', encoding='utf-8') as f:
                enrichment = json.load(f)
        except Exception:
            features["asm.distortion_total_kcal"] = None
            features["asm.interaction_kcal"] = None
            return features

        sp_results = enrichment.get('sp_results', {})
        fragA_R = sp_results.get('fragA_R')
        fragB_R = sp_results.get('fragB_R')
        fragA_TS = sp_results.get('fragA_TS')
        fragB_TS = sp_results.get('fragB_TS')
        ts_TS = sp_results.get('ts_TS')

        distortion_total_kcal = None
        interaction_kcal = None

        if fragA_R is not None and fragA_TS is not None:
            distortion_A = (fragA_TS - fragA_R) * HARTREE_TO_KCAL
        else:
            distortion_A = None

        if fragB_R is not None and fragB_TS is not None:
            distortion_B = (fragB_TS - fragB_R) * HARTREE_TO_KCAL
        else:
            distortion_B = None

        if distortion_A is not None and distortion_B is not None:
            distortion_total_kcal = distortion_A + distortion_B

        features["asm.distortion_total_kcal"] = distortion_total_kcal

        if ts_TS is not None and fragA_TS is not None and fragB_TS is not None:
            interaction_kcal = (ts_TS - (fragA_TS + fragB_TS)) * HARTREE_TO_KCAL

        features["asm.interaction_kcal"] = interaction_kcal

        return features


register_extractor(ASMEnrichmentExtractor())
