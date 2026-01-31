"""
Step 4: Feature Context and Result Types (V4.2 Phase A)
======================================================

FeatureContext, PluginTrace, and FeatureResult dataclasses.

Author: QC Descriptors Team
Date: 2026-01-18
"""

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
import time
import json

from .status import FeatureResultStatus, FeatureStatus, NICSTriggerPolicy


@dataclass
class PluginTrace:
    """Execution trace for a single plugin.

    Written to `feature_meta.json` under `trace.plugins.<plugin_name>`.
    """
    plugin_name: str
    status: FeatureResultStatus = FeatureResultStatus.SKIPPED
    runtime_ms: float = 0.0
    missing_fields: List[str] = field(default_factory=list)
    missing_paths: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    job_specs: List[Dict[str, Any]] = field(default_factory=list)
    _extracted_features: Dict[str, Any] = field(default_factory=dict)

    def record_start(self) -> float:
        """Record plugin start time.

        Returns:
            Start timestamp for later end measurement
        """
        self.runtime_ms = 0.0
        return time.time()

    def record_end(self, start_time: float) -> None:
        """Record plugin end time and calculate runtime.

        Args:
            start_time: Start timestamp from record_start()
        """
        self.runtime_ms = (time.time() - start_time) * 1000.0

    def add_error(self, error: str) -> None:
        """Add an error message to the trace.

        Args:
            error: Error message string
        """
        self.errors.append(error)

    def add_warning(self, warning: str) -> None:
        """Add a warning message to the trace.

        Args:
            warning: Warning message string
        """
        self.warnings.append(warning)


@dataclass
class FeatureContext:
    """Context object containing all inputs and configuration for feature extraction.

    Plugins access data through this context rather than raw parameters.
    All paths are validated by PathAccessor before being passed to plugins.
    """

    ts_xyz: Optional[Path] = None
    reactant_xyz: Optional[Path] = None
    product_xyz: Optional[Path] = None
    ts_fchk: Optional[Path] = None
    reactant_fchk: Optional[Path] = None
    product_fchk: Optional[Path] = None
    ts_orca_out: Optional[Path] = None
    ts_log: Optional[Path] = None
    reactant_log: Optional[Path] = None
    product_log: Optional[Path] = None
    ts_qm_output: Optional[Path] = None
    reactant_qm_output: Optional[Path] = None
    product_qm_output: Optional[Path] = None

    forming_bonds: Optional[Tuple[Tuple[int, int], ...]] = None
    fragment_indices: Optional[Tuple[List[int], List[int]]] = None

    sp_report: Optional[Any] = None

    # V6.2: S3 directory for accessing enrichment contract and dipolar artifacts
    # Resolved in feature_miner before plugin execution
    s3_dir: Optional[Path] = None
    artifacts_index: Optional[Dict[str, Any]] = None

    close_contacts_cutoff: float = 2.2
    temperature_K: float = 298.15

    ts_fingerprint: Optional[Dict[str, Any]] = None
    reactant_fingerprint: Optional[Dict[str, Any]] = None
    product_fingerprint: Optional[Dict[str, Any]] = None
    ts_fchk_fingerprint: Optional[Dict[str, Any]] = None
    reactant_fchk_fingerprint: Optional[Dict[str, Any]] = None
    product_fchk_fingerprint: Optional[Dict[str, Any]] = None
    ts_orca_out_fingerprint: Optional[Dict[str, Any]] = None
    sp_report_fingerprint: Optional[str] = None

    multiwfn_recipe_path: Optional[Path] = None
    nbo_template_path: Optional[Path] = None

    nics_trigger_policy: NICSTriggerPolicy = NICSTriggerPolicy.GENERATE_ONLY
    job_run_policy: str = "disallow"

    plugin_traces: Dict[str, PluginTrace] = field(default_factory=dict)

    def get_path(self, key: str) -> Optional[Path]:
        """Get a path from context by canonical key name.

        Args:
            key: Canonical context key (e.g., "ts_xyz", "reactant_xyz")

        Returns:
            Path object or None if not set
        """
        return getattr(self, key, None)

    def has_path(self, key: str) -> bool:
        """Check if a path key exists and is not None.

        Args:
            key: Canonical context key name

        Returns:
            True if path exists and is not None
        """
        path_obj = self.get_path(key)
        return path_obj is not None

    def get_plugin_trace(self, plugin_name: str) -> PluginTrace:
        """Get or create a PluginTrace for a plugin.

        Args:
            plugin_name: Name of the plugin

        Returns:
            PluginTrace instance
        """
        if plugin_name not in self.plugin_traces:
            self.plugin_traces[plugin_name] = PluginTrace(plugin_name=plugin_name)
        return self.plugin_traces[plugin_name]


@dataclass
class FeatureResult:
    """Complete feature extraction result with metadata.

    Contains the flattened features dictionary and reproducibility trace.
    """

    features: Dict[str, Any]
    schema_version: str = "4.2-tag"
    schema_signature: str = ""
    feature_status: FeatureStatus = FeatureStatus.OK

    # Method/solvent information
    method: str = ""
    solvent: str = ""
    temperature_K: float = 298.15

    # Enabled plugins
    enabled_plugins: List[str] = field(default_factory=list)

    ts_fingerprint: Optional[Dict[str, Any]] = None
    reactant_fingerprint: Optional[Dict[str, Any]] = None
    ts_fchk_fingerprint: Optional[Dict[str, Any]] = None
    reactant_fchk_fingerprint: Optional[Dict[str, Any]] = None
    product_fchk_fingerprint: Optional[Dict[str, Any]] = None
    product_fingerprint: Optional[Dict[str, Any]] = None
    ts_orca_out_fingerprint: Optional[Dict[str, Any]] = None
    sp_report_fingerprint: Optional[str] = None

    # Global warnings list (aggregated from all plugins)
    warnings: List[str] = field(default_factory=list)

    # Plugin traces
    plugin_traces: Dict[str, PluginTrace] = field(default_factory=dict)

    # Job policies (for reproducibility documentation)
    job_run_policy: str = "disallow"
    nics_trigger_policy: NICSTriggerPolicy = NICSTriggerPolicy.GENERATE_ONLY

    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary (excluding None/NaN from features).

        Returns:
            Dictionary representation suitable for JSON serialization
        """
        import numpy as np

        filtered_features = {
            k: v for k, v in self.features.items()
            if v is not None and not (isinstance(v, float) and np.isnan(v))
        }

        result = {
            "meta": {
                "schema_version": self.schema_version,
                "schema_signature": self.schema_signature,
                "feature_status": self.feature_status.value,
                "method": self.method,
                "solvent": self.solvent,
                "temperature_K": self.temperature_K,
                "enabled_plugins": self.enabled_plugins,
                "policies": {
                    "job_run_policy": self.job_run_policy,
                    "nics_trigger_policy": self.nics_trigger_policy.value,
                },
            },
            "trace": {
                "inputs_fingerprint": {},
                "plugins": {},
            },
            "warnings": self.warnings,
        }

        fingerprint_mapping = {
            "ts_xyz": self.ts_fingerprint,
            "reactant_xyz": self.reactant_fingerprint,
            "ts_fchk": self.ts_fchk_fingerprint,
            "reactant_fchk": self.reactant_fchk_fingerprint,
            "product_fchk": self.product_fchk_fingerprint,
            "product_xyz": self.product_fingerprint,
            "ts_orca_out": self.ts_orca_out_fingerprint,
        }

        for key, fp in fingerprint_mapping.items():
            if fp:
                result["trace"]["inputs_fingerprint"][key] = fp

        if self.sp_report_fingerprint:
            result["trace"]["inputs_fingerprint"]["sp_report_summary"] = self.sp_report_fingerprint

        for plugin_name, trace in self.plugin_traces.items():
            result["trace"]["plugins"][plugin_name] = {
                "status": trace.status.value,
                "runtime_ms": trace.runtime_ms,
                "missing_fields": trace.missing_fields,
                "missing_paths": trace.missing_paths,
                "errors": trace.errors,
                "warnings": trace.warnings,
                "cache_hit": False,
                "job_specs": trace.job_specs,
            }

        return result

    def to_json(self, output_path: Path) -> None:
        """Write feature result metadata to JSON file.

        Args:
            output_path: Path to feature_meta.json
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        result_dict = self.to_dict()
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(result_dict, f, indent=2, default=str)
