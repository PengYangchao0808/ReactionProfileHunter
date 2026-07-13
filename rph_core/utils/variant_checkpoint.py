"""Variant-aware checkpoint for multi-variant V4 pipeline."""

from __future__ import annotations

from collections.abc import Mapping
import hashlib
import json
from pathlib import Path
from typing import cast


class VariantCheckpoint:
    """Per-variant checkpoint state stored in a single .rph/checkpoint.json file."""

    def __init__(self, work_dir: Path):
        self.work_dir: Path = Path(work_dir)
        self.rph_dir: Path = self.work_dir / ".rph"
        _ = self.rph_dir.mkdir(parents=True, exist_ok=True)
        self.path: Path = self.rph_dir / "checkpoint.json"
        self._state: dict[str, object] = self._load()

    def _load(self) -> dict[str, object]:
        if self.path.exists():
            payload = cast(object, json.loads(self.path.read_text(encoding="utf-8")))
            if isinstance(payload, dict):
                return cast(dict[str, object], payload)
        return {"schema_version": "rph_variant_checkpoint_v1", "variants": {}}

    def _save(self) -> None:
        _ = self.path.write_text(json.dumps(self._state, indent=2), encoding="utf-8")

    def set_reaction_id(self, reaction_id: str) -> None:
        self._state["reaction_id"] = reaction_id
        self._save()

    def reusable(self, scope: str, stage: str, signature: str, manifest: Path) -> bool:
        """Check whether a stage within a scope can be reused."""

        entry = self._get_entry(scope, stage)
        if not entry:
            return False
        if entry.get("signature") != signature:
            return False
        if entry.get("status") != "complete":
            return False
        recorded_manifest = entry.get("manifest")
        if isinstance(recorded_manifest, str):
            if Path(recorded_manifest).resolve() != Path(manifest).resolve():
                return False
        return Path(manifest).exists()

    def mark(self, scope: str, stage: str, signature: str, manifest: Path, status: str = "complete") -> None:
        """Mark a stage within a scope as complete."""

        self._set_entry(scope, stage, {
            "signature": signature,
            "manifest": str(Path(manifest).resolve()),
            "status": status,
        })
        self._save()

    def mark_failed(self, scope: str, stage: str, signature: str, error: str) -> None:
        """Mark a stage as failed."""

        self._set_entry(scope, stage, {
            "signature": signature,
            "status": "failed",
            "error": error,
        })
        self._save()

    def _get_entry(self, scope: str, stage: str) -> dict[str, object] | None:
        if scope in ("s0", "s3", "s4", "precursor_s1"):
            entry = self._state.get(scope)
            return cast(dict[str, object], entry) if isinstance(entry, dict) else None
        variants = self._state.get("variants", {})
        if not isinstance(variants, dict):
            return None
        variants_map = cast(dict[str, object], variants)
        scoped = variants_map.get(scope, {})
        if not isinstance(scoped, dict):
            return None
        scoped_map = cast(dict[str, object], scoped)
        entry = scoped_map.get(stage)
        return cast(dict[str, object], entry) if isinstance(entry, dict) else None

    def _set_entry(self, scope: str, stage: str, value: Mapping[str, object]) -> None:
        payload = dict(value)
        if scope in ("s0", "s3", "s4", "precursor_s1"):
            self._state[scope] = payload
            return
        variants = self._state.get("variants")
        if not isinstance(variants, dict):
            variants = {}
            self._state["variants"] = variants
        variants_map = cast(dict[str, object], variants)
        scoped = variants_map.get(scope)
        if not isinstance(scoped, dict):
            scoped = {}
            variants_map[scope] = scoped
        scoped_map = cast(dict[str, object], scoped)
        scoped_map[stage] = payload

    @staticmethod
    def signature(payload: Mapping[str, object]) -> str:
        """SHA-256 signature of a payload dict."""

        return hashlib.sha256(
            json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
        ).hexdigest()
