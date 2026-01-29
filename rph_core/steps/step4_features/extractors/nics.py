"""
Step 4: NICS Extractor (V4.2 Phase B/C)
=============================================

Multiwfn NICS descriptor extractor (dry-run only).

Phase A: SKIPPED (no NICS job specs generated)
Phase B: Generate job_specs only (dry-run, no execution)
Phase C: Execute Multiwfn and parse NICS outputs (not implemented)

Author: QC Descriptors Team
Date: 2026-01-18
"""

from typing import Dict, Any, List
from pathlib import Path
import hashlib

from .base import BaseExtractor, register_extractor
from ..status import FeatureResultStatus, NICSTriggerPolicy
from ..schema import generate_cache_key
from ..multiwfn.recipes import load_multiwfn_recipes, render_recipe


class NICSExtractor(BaseExtractor):
    """Extract NICS aromaticity descriptors using Multiwfn.

    Features (arom.* prefix, Phase C only):
    - arom.nics_1zz
    - arom.nics_0zz
    - Additional aromaticity descriptors

    Phase A/B behavior:
    - Always SKIPPED with missing_fields/missing_paths if ts_fchk missing
    - If nics_trigger_policy == SKIP: SKIPPED with status
    - If nics_trigger_policy == GENERATE_ONLY: Generate job_specs, no features
    - If nics_trigger_policy == AUTO_RUN and job_run_policy != allow: Same as GENERATE_ONLY

    Phase C behavior (not implemented):
    - Execute Multiwfn using recipe from config
    - Parse output files for NICS values
    - Write arom.* features
    """

    def get_plugin_name(self) -> str:
        return "nics"

    def get_required_inputs(self) -> List[str]:
        return ["ts_fchk"]

    def extract(self, context) -> Dict[str, Any]:
        ts_fchk = context.ts_fchk
        nics_policy = context.nics_trigger_policy
        job_run_policy = context.job_run_policy

        trace = context.get_plugin_trace(self.get_plugin_name())

        if nics_policy == NICSTriggerPolicy.SKIP:
            trace.status = FeatureResultStatus.SKIPPED
            trace.add_warning("NICS generation skipped by policy")
            return {}

        if nics_policy == NICSTriggerPolicy.AUTO_RUN and job_run_policy != "allow":
            trace.add_warning("NICS AUTO_RUN policy set but job_run_policy=disallow, generating job_specs only")

        recipe_path = context.multiwfn_recipe_path or Path(__file__).parent.parent / "config" / "multiwfn_recipes.yaml"

        try:
            recipes = load_multiwfn_recipes(recipe_path)
        except FileNotFoundError:
            trace.add_error(f"Multiwfn recipes file not found: {recipe_path}")
            trace.status = FeatureResultStatus.FAILED
            return {}
        except Exception as e:
            trace.add_error(f"Failed to load Multiwfn recipes: {e}")
            trace.status = FeatureResultStatus.FAILED
            return {}

        recipe_name = "nics_1zz"
        if recipe_name not in recipes:
            trace.add_error(f"Recipe '{recipe_name}' not found in {recipe_path}")
            trace.status = FeatureResultStatus.FAILED
            return {}

        recipe = recipes[recipe_name]

        try:
            stat = recipe_path.stat()
            recipe_yaml_fingerprint = {
                "hash": hashlib.sha256(recipe_path.read_bytes()).hexdigest()[:16],
                "size": stat.st_size,
                "mtime": stat.st_mtime,
            }
        except Exception:
            recipe_yaml_fingerprint = None

        stdin_lines = render_recipe(
            recipe,
            input_fchk=ts_fchk,
            workdir=ts_fchk.parent,
            output_prefix="nics_1zz",
        )

        input_files = {"ts_fchk": ts_fchk}
        params = {
            "policy": nics_policy.value,
            "job_run_policy": job_run_policy,
            "recipe_name": recipe_name,
            "recipe_version": recipe.version,
        }

        cache_key = generate_cache_key(
            plugin_name=self.get_plugin_name(),
            input_files=input_files,
            params=params,
        )

        job_spec = {
            "engine": "multiwfn",
            "recipe_name": recipe_name,
            "recipe_version": recipe.version,
            "input_files": {"ts_fchk": str(ts_fchk)},
            "workdir": str(ts_fchk.parent),
            "command": "multiwfn",
            "stdin_lines": stdin_lines,
            "expected_outputs": recipe.outputs,
            "cache_key": cache_key,
            "recipe_yaml_fingerprint": recipe_yaml_fingerprint,
        }

        trace.job_specs = [job_spec]
        trace.add_warning("Phase B placeholder: NICS job_specs generated (dry-run only, not executed)")

        return {}


register_extractor(NICSExtractor())
