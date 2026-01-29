"""
Step 4: NBO E(2) Extractor (V4.2 Phase B/C)
==================================================

NBO E(2) template-based matching extractor (dry-run only).

Phase A: SKIPPED (no NBO file, template matching not applicable)
Phase B: Generate job_specs only (dry-run, no parsing)
Phase C: Parse NBO output and match against templates

Author: QC Descriptors Team
Date: 2026-01-18
"""

from typing import Dict, Any, List
from pathlib import Path
import hashlib
import yaml

from .base import BaseExtractor, register_extractor
from ..status import FeatureResultStatus
from ..schema import generate_cache_key
from ..nbo.e2_parser import find_nbo_e2_section, parse_nbo_e2_table, match_templates


class NBOE2Extractor(BaseExtractor):
    """Extract NBO E(2) interaction energies using template whitelist.

    Features (nbo.e2.* prefix, Phase C only):
    - nbo.e2.<template_name>.e2_kcal for each matched template

    Phase B behavior:
    - Generate job_specs for NBO computation
    - No features written (dry-run)

    Contract 3 requirements:
    - Donor and acceptor must match template exactly
    - 0 or >1 matches: result is None + warning
    """

    def get_plugin_name(self) -> str:
        return "nbo_e2"

    def get_required_inputs(self) -> List[str]:
        return ["ts_fchk"]

    def get_required_inputs_for_context(self, context) -> List[str]:
        if context.ts_orca_out is not None:
            return ["ts_orca_out"]
        return ["ts_fchk"]

    def extract(self, context) -> Dict[str, Any]:
        trace = context.get_plugin_trace(self.get_plugin_name())
        features = {}

        if context.ts_orca_out and context.ts_orca_out.exists():
            try:
                text = context.ts_orca_out.read_text(encoding="utf-8", errors="replace")

                section = find_nbo_e2_section(text)
                if not section:
                    trace.errors.append("NBO E(2) section not found in ORCA output")
                    trace.status = FeatureResultStatus.FAILED
                    return {}

                interactions = parse_nbo_e2_table(section)
                if not interactions:
                    trace.warnings.append("Found E(2) section but failed to parse any interactions")
                    return {}

                tmpl_path = getattr(context, "nbo_template_path", None)
                if not tmpl_path:
                    base_dir = Path(__file__).parent.parent
                    tmpl_path = base_dir / "config" / "nbo_templates.yaml"

                if not tmpl_path.exists():
                    trace.errors.append(f"NBO templates file not found: {tmpl_path}")
                    trace.status = FeatureResultStatus.FAILED
                    return {}

                with open(tmpl_path, "r") as f:
                    yaml_data = yaml.safe_load(f)

                if not isinstance(yaml_data, dict):
                    trace.errors.append("Invalid NBO templates file structure")
                    trace.status = FeatureResultStatus.FAILED
                    return {}

                templates_dict = {}
                for tmpl in yaml_data.get("nbo_e2_templates", []):
                    if isinstance(tmpl, dict) and "name" in tmpl:
                        templates_dict[tmpl["name"]] = tmpl

                matched_feats, warnings = match_templates(interactions, templates_dict)

                features.update(matched_feats)
                trace.warnings.extend(warnings)

                return features

            except Exception as e:
                trace.errors.append(f"Error parsing NBO output: {str(e)}")
                trace.status = FeatureResultStatus.FAILED
                return {}

        trace.warnings.append("Missing ts_orca_out; generating NBO job spec (Phase B)")

        template_path = context.nbo_template_path or Path(__file__).parent.parent / "config" / "nbo_templates.yaml"

        if not template_path.exists():
            trace.errors.append(f"NBO templates file not found: {template_path}")
            trace.status = FeatureResultStatus.FAILED
            return {}

        try:
            stat = template_path.stat()
            template_yaml_fingerprint = {
                "hash": hashlib.sha256(template_path.read_bytes()).hexdigest()[:16],
                "size": stat.st_size,
                "mtime": stat.st_mtime,
            }
        except Exception:
            template_yaml_fingerprint = None

        ts_fchk = context.ts_fchk
        if ts_fchk is None:
            trace.errors.append("ts_fchk is None, cannot generate job spec")
            trace.status = FeatureResultStatus.FAILED
            return {}

        assert ts_fchk is not None
        input_files = {"ts_fchk": ts_fchk}
        params = {
            "template_count": 0,
            "template_path": str(template_path),
        }

        try:
            with open(template_path, "r") as f:
                yaml_data = yaml.safe_load(f)
                if isinstance(yaml_data, dict):
                    params["template_count"] = len(yaml_data.get("nbo_e2_templates", []))
        except Exception:
            pass

        cache_key = generate_cache_key(
            plugin_name=self.get_plugin_name(),
            input_files=input_files,
            params=params,
        )

        job_spec = {
            "engine": "orca",
            "recipe_name": "nbo_e2",
            "input_files": {"ts_fchk": str(ts_fchk)},
            "workdir": str(ts_fchk.parent),
            "command": "nbo_e2_from_fchk",
            "expected_outputs": ["nbo_e2.out"],
            "cache_key": cache_key,
            "template_yaml_fingerprint": template_yaml_fingerprint,
        }

        trace.job_specs = [job_spec]

        return {}


register_extractor(NBOE2Extractor())
