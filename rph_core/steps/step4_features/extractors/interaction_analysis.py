"""
Step 4: Interaction Analysis Extractor (V4.2 Phase B/C)
================================================================

Interaction/fragment distortion analysis extractor.

Phase A: Interface stub only - SKIPPED if fragment_indices missing.
Phase B: Generate job_specs only (dry-run, no eda.*_aux features)
Phase C: Execute fragment SP and output eda.*_aux features (not implemented).

Author: QC Descriptors Team
Date: 2026-01-18
"""

from pathlib import Path
from typing import Dict, Any, List

from .base import BaseExtractor, register_extractor
from ..status import FeatureResultStatus
from ..schema import generate_cache_key


class InteractionAnalysisExtractor(BaseExtractor):
    """Extract interaction/fragment analysis features.

    Features (eda.* prefix, Phase C only, Contract 4):
    - eda.E_orb_aux: Orbital interaction energy proxy
    - eda.Q_CT_aux: Charge transfer proxy

    Requirements (Contract 4):
    - fragment_indices must be provided: Tuple[List[int], List[int]]
    - reactant_fchk must be provided for fragment analysis

    Phase A behavior:
    - Always SKIPPED with missing_fields/missing_paths if inputs missing
    - Does NOT run any external QM calculations (dry-run)

    Phase B behavior:
    - Generate job_specs for fragment SP calculations (dry-run only)
    - No eda.*_aux features written (job_specs only in trace)

    Phase C behavior (not implemented):
    - Execute fragment SP calculations
    - Parse results and write eda.E_orb_aux and eda.Q_CT_aux
    """

    def get_plugin_name(self) -> str:
        return "interaction"

    def get_required_inputs(self) -> List[str]:
        return ["reactant_fchk", "fragment_indices"]

    def extract(self, context) -> Dict[str, Any]:
        trace = context.get_plugin_trace(self.get_plugin_name())
        fragment_indices = context.fragment_indices
        reactant_fchk = context.reactant_fchk
        ts_fchk = context.ts_fchk
        job_run_policy = context.job_run_policy

        assert reactant_fchk is not None

        reactant_path = reactant_fchk

        input_files_dict = {"reactant_fchk": reactant_path}

        if ts_fchk is not None:
            input_files_dict["ts_fchk"] = ts_fchk

        params = {
            "fragment_indices": list(fragment_indices) if fragment_indices else None,
            "job_run_policy": job_run_policy,
        }

        cache_key = generate_cache_key(
            plugin_name=self.get_plugin_name(),
            input_files=input_files_dict,
            params=params,
        )

        job_spec = {
            "engine": "eda_aux",
            "recipe_name": "fragment_distortion",
            "input_files": {k: str(v) for k, v in input_files_dict.items()},
            "workdir": str(reactant_path.parent),
            "command": "fragment_sp",
            "expected_outputs": ["fragment1_energies.txt", "fragment2_energies.txt"],
            "cache_key": cache_key,
        }

        trace.job_specs = [job_spec]
        trace.add_warning("Phase B placeholder: EDA aux job_specs generated (dry-run only, eda.*_aux features not written)")

        return {}


register_extractor(InteractionAnalysisExtractor())
