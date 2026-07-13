"""Small manifest-based checkpoint store for the V4 pipeline."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Optional


class V4Checkpoint:
    filename = "pipeline.state"
    stage_order = ("s0", "s1", "s2", "s3", "s4")

    def __init__(self, work_dir: Path):
        self.work_dir = Path(work_dir)
        self.path = self.work_dir / self.filename

    def load(self) -> Dict[str, Any]:
        if not self.path.exists():
            return {"schema_version": "rph_v4_s0_s4_v1", "stages": {}}
        data = json.loads(self.path.read_text(encoding="utf-8"))
        if data.get("schema_version") != "rph_v4_s0_s4_v1":
            raise RuntimeError("legacy pipeline state is unsupported by RPH V4")
        return data

    def save(self, data: Dict[str, Any]) -> None:
        self.work_dir.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(".tmp")
        temporary.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")
        temporary.replace(self.path)

    @staticmethod
    def signature(payload: Dict[str, Any]) -> str:
        raw = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()

    @staticmethod
    def file_signature(path: Path) -> str:
        """Hash a file for stage input identity without reading it into JSON."""
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def reusable(self, stage: str, signature: str, manifest: Path) -> bool:
        state = self.load()
        record = (state.get("stages", {}) or {}).get(stage, {})
        if record.get("status") != "complete":
            return False
        recorded_manifest = Path(record.get("manifest", ""))
        requested_manifest = Path(manifest)
        return (
            record.get("signature") == signature
            and recorded_manifest.resolve() == requested_manifest.resolve()
            and requested_manifest.exists()
        )

    def invalidate_from(self, stage: str) -> None:
        """Remove a stage and all downstream stages from the resume state."""
        if stage not in self.stage_order:
            raise ValueError(f"Unknown V4 stage: {stage}")
        state = self.load()
        stages = state.setdefault("stages", {})
        start = self.stage_order.index(stage)
        for downstream in self.stage_order[start:]:
            stages.pop(downstream, None)
        self.save(state)

    def mark(self, stage: str, signature: str, manifest: Path, status: str = "complete") -> None:
        if stage not in self.stage_order:
            raise ValueError(f"Unknown V4 stage: {stage}")
        state = self.load()
        stages = state.setdefault("stages", {})
        stages[stage] = {
            "signature": signature,
            "manifest": str(Path(manifest).resolve()),
            "status": status,
        }
        start = self.stage_order.index(stage) + 1
        for downstream in self.stage_order[start:]:
            stages.pop(downstream, None)
        self.save(state)

    def mark_s0(self, signature: str, manifest: Path, status: str = "complete") -> None:
        """Mark record/template S0 without invalidating independent S1 output.

        V4 S0 validates reaction metadata while S1 produces product conformers;
        both are inputs to S2.  A changed S0 result therefore invalidates S2
        and later stages, but must not discard the expensive CENSO-LITE S1
        manifest.
        """

        state = self.load()
        stages = state.setdefault("stages", {})
        stages["s0"] = {
            "signature": signature,
            "manifest": str(Path(manifest).resolve()),
            "status": status,
        }
        for downstream in ("s2", "s3", "s4"):
            stages.pop(downstream, None)
        self.save(state)
