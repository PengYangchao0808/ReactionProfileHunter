"""Compatibility adapter for the unified V4 checkpoint store."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from rph_core.utils.v4_checkpoint import V4Checkpoint


class VariantCheckpoint:
    """Deprecated adapter that now writes only to ``pipeline.state``.

    New code should use :class:`V4Checkpoint` directly.  The adapter remains
    temporarily available for callers outside the orchestrator and for loading
    existing integrations without recreating the former dual-checkpoint state.
    """

    def __init__(self, work_dir: Path):
        self._checkpoint = V4Checkpoint(work_dir)
        self.work_dir = Path(work_dir)
        self.path = self._checkpoint.path

    def set_reaction_id(self, reaction_id: str) -> None:
        self._checkpoint.set_reaction_id(reaction_id)

    def reusable(
        self, scope: str, stage: str, signature: str, manifest: Path
    ) -> bool:
        return self._checkpoint.reusable_scope(scope, stage, signature, manifest)

    def mark(
        self,
        scope: str,
        stage: str,
        signature: str,
        manifest: Path,
        status: str = "complete",
    ) -> None:
        self._checkpoint.mark_scope(scope, stage, signature, manifest, status)

    def mark_failed(
        self, scope: str, stage: str, signature: str, error: str
    ) -> None:
        self._checkpoint.mark_failed(scope, stage, signature, error)

    @staticmethod
    def signature(payload: Mapping[str, object]) -> str:
        return V4Checkpoint.signature(payload)
