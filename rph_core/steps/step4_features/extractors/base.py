"""
Step 4: Extractor Base (V4.2 Phase A/B/C)
==========================================

Base extractor class and plugin registry with dynamic dependency support.

Author: QC Descriptors Team
Date: 2026-01-18
"""

import logging
import time
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List

from ..context import FeatureContext, PluginTrace
from ..status import FeatureResultStatus
from ..path_accessor import PathAccessor
from ..schema import JOB_SPEC_REQUIRED_KEYS, JOB_SPEC_OPTIONAL_KEYS

logger = logging.getLogger(__name__)

PATH_KEYS = {
    "ts_xyz", "reactant_xyz", "product_xyz",
    "ts_fchk", "reactant_fchk", "product_fchk",
    "ts_orca_out",
}


def validate_job_spec(job_spec: Dict[str, Any]) -> List[str]:
    """Lightweight validation of job spec against schema contract."""
    errors = []
    keys = set(job_spec.keys())
    missing = JOB_SPEC_REQUIRED_KEYS - keys
    if missing:
        errors.append(f"Missing required keys: {missing}")
    return errors


class BaseExtractor(ABC):
    """Base class for all feature extraction plugins.

    All extractors must implement:
    - `get_required_inputs()` - List of required context keys
    - `extract()` - Main extraction logic
    - `get_plugin_name()` - Unique plugin identifier
    - `can_submit_jobs` - Whether plugin can trigger QC jobs (default: False)

    V6.1: Plugins with can_submit_jobs=True are only allowed in job_generate_mode
    V6.1: If job_run_policy=disallow and can_submit_jobs=True, plugin may only generate job_specs.

    Plugins are discovered and registered in the EXTRACTORS dict.
    """

    @abstractmethod
    def get_plugin_name(self) -> str:
        """Return unique plugin name.

        Used for registration and trace identification.

        Returns:
            Plugin name string (e.g., "thermo", "geometry", "qc_checks")
        """
        pass

    def can_submit_jobs(self) -> bool:
        """Return whether this extractor can submit QC jobs.

        V6.1: Plugins with can_submit_jobs=True are only allowed in job_generate_mode.
        V6.1: If job_run_policy=disallow and can_submit_jobs=True, plugin may only generate job_specs.

        Default: False (extract-only plugins).
        """
        return False

    @abstractmethod
    def get_required_inputs(self) -> List[str]:
        """Return list of required context keys (static default).

        Used by PathAccessor to validate inputs before extraction.
        Override get_required_inputs_for_context() for dynamic behavior.

        Returns:
            List of canonical context key names (e.g., ["sp_report", "ts_xyz"])
        """
        pass

    def get_required_inputs_for_context(self, context: FeatureContext) -> List[str]:
        """Dynamic hook to determine required inputs based on available context.

        Defaults to static get_required_inputs().
        Override this to support alternative dependency paths (e.g. Phase C vs Phase B).

        Args:
            context: FeatureContext containing all inputs

        Returns:
            List of required context key names
        """
        return self.get_required_inputs()

    @abstractmethod
    def extract(self, context: FeatureContext) -> Dict[str, Any]:
        """Extract features from context.

        Args:
            context: FeatureContext containing all inputs

        Returns:
            Dictionary of extracted features (with prefixed keys like "thermo.dG_activation")

        Raises:
            Exceptions should be caught by caller and recorded in PluginTrace
        """
        pass

    def validate_inputs(self, context: FeatureContext, trace: PluginTrace) -> bool:
        """Validate that all required inputs are available.

        Args:
            context: FeatureContext containing all inputs
            trace: PluginTrace to populate with missing field/path errors

        Returns:
            True if all required inputs present, False otherwise
        """
        required = self.get_required_inputs_for_context(context)
        missing = []

        for key in required:
            if not hasattr(context, key) or getattr(context, key, None) is None:
                missing.append(key)

        if not missing:
            return True

        for key in missing:
            if key in PATH_KEYS:
                trace.missing_paths.append(key)
            else:
                trace.missing_fields.append(key)

        trace.status = FeatureResultStatus.SKIPPED
        trace.errors.append(f"Missing required inputs: {missing}")
        return False

    def run(self, context: FeatureContext) -> PluginTrace:
        """Run extractor and return trace.

        Args:
            context: FeatureContext containing all inputs

        Returns:
            PluginTrace with execution results
        """
        name = self.get_plugin_name()
        trace = context.get_plugin_trace(name)

        trace.missing_fields = []
        trace.missing_paths = []
        trace.errors = []
        trace.warnings = []
        trace.job_specs = []
        trace._extracted_features = {}
        trace.runtime_ms = 0.0

        start_time = time.time()

        try:
            if self.validate_inputs(context, trace):
                features = self.extract(context)
                trace._extracted_features = features
                trace.status = FeatureResultStatus.OK

            if trace.job_specs:
                for idx, spec in enumerate(trace.job_specs):
                    spec_errors = validate_job_spec(spec)
                    for err in spec_errors:
                        trace.errors.append(f"JobSpec[{idx}] Invalid: {err}")

                if any("JobSpec" in e for e in trace.errors):
                    trace.status = FeatureResultStatus.FAILED

        except Exception as e:
            logger.error(f"{name} extractor failed: {e}", exc_info=True)
            trace.status = FeatureResultStatus.FAILED
            trace.errors.append(f"Uncaught exception: {str(e)}")

        finally:
            trace.runtime_ms = (time.time() - start_time) * 1000.0

        if trace.status == FeatureResultStatus.OK:
            has_features = bool(trace._extracted_features)
            has_jobs = bool(trace.job_specs)

            if not has_features and has_jobs:
                trace.status = FeatureResultStatus.SKIPPED

        return trace

    def get_extracted_features(self, trace: PluginTrace) -> Dict[str, Any]:
        """Get features from a PluginTrace.

        Args:
            trace: PluginTrace from run()

        Returns:
            Feature dictionary (empty if extraction failed)
        """
        return getattr(trace, "_extracted_features", {})


# Plugin registry
EXTRACTORS: Dict[str, BaseExtractor] = {}


def register_extractor(extractor: BaseExtractor) -> None:
    """Register an extractor plugin.

    Args:
        extractor: BaseExtractor instance
    """
    name = extractor.get_plugin_name()
    if name in EXTRACTORS:
        logger.warning(f"Extractor '{name}' already registered, overwriting")

    EXTRACTORS[name] = extractor
    logger.info(f"Registered extractor: {name}")


def get_extractor(name: str) -> Optional[BaseExtractor]:
    """Get a registered extractor by name.

    Args:
        name: Plugin name

    Returns:
        BaseExtractor instance or None if not found
    """
    return EXTRACTORS.get(name, None)


def list_extractors() -> List[str]:
    """List all registered extractor names.

    Returns:
        List of plugin names
    """
    return list(EXTRACTORS.keys())
