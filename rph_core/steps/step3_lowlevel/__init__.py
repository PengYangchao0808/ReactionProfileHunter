"""S3 low-level QC stage (backward-compat alias for RefinementEngine).

Phase 3+: ``LowLevelEngine`` is just a thin alias layer over
``RefinementEngine``. It accepts both the legacy constructor shape
``(config, run_id=None, parent_manifest_paths=None)`` and the inherited
``RefinementEngine`` shape ``(config, profile, run_id=None, ...)``.
"""

from collections.abc import Mapping
from pathlib import Path
from typing import Any, Optional

from rph_core.steps.fidelity_profile import FidelityProfile
from rph_core.steps.refinement.engine import RefinementEngine


class LowLevelEngine(RefinementEngine):
    """Backward-compat alias for S3 refinement execution."""

    stage_name = "S3"

    def __init__(
        self,
        config: Any,
        profile: FidelityProfile | str | Mapping[str, Path] | None = None,
        run_id: Optional[str] | Mapping[str, Path] = None,
        parent_manifest_paths: Optional[Mapping[str, Path]] = None,
    ):
        resolved_profile: FidelityProfile
        resolved_run_id: Optional[str]
        resolved_parent_manifest_paths = parent_manifest_paths

        if isinstance(profile, FidelityProfile):
            resolved_profile = profile
            if isinstance(run_id, Mapping) and parent_manifest_paths is None:
                resolved_run_id = None
                resolved_parent_manifest_paths = run_id
            else:
                resolved_run_id = run_id if isinstance(run_id, str) else None
        else:
            resolved_profile = FidelityProfile.from_config(config, "S3")
            resolved_run_id = str(profile) if isinstance(profile, str) else None
            if isinstance(profile, Mapping) and parent_manifest_paths is None and run_id is None:
                resolved_parent_manifest_paths = profile
            elif isinstance(run_id, Mapping) and parent_manifest_paths is None:
                resolved_parent_manifest_paths = run_id
            elif isinstance(run_id, str):
                resolved_run_id = run_id

        super().__init__(
            config=config,
            profile=resolved_profile,
            run_id=resolved_run_id,
            parent_manifest_paths=resolved_parent_manifest_paths,
        )


__all__ = ["LowLevelEngine"]
