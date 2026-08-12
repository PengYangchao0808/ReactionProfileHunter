"""Manifest-based checkpoint store for the V4 S0-S4 pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

from rph_core.utils.json_io import write_json_atomic
from rph_core.utils.run_id import RUN_ID_FIELD, is_valid_run_id


logger = logging.getLogger(__name__)


STATE_SCHEMA_V1 = "rph_v4_s0_s4_v1"
STATE_SCHEMA_V2 = "rph_v4_checkpoint_v2"
VARIANT_SCHEMA_V1 = "rph_variant_checkpoint_v1"


@dataclass(frozen=True)
class ReuseDecision:
    """Explain whether one checkpoint entry can be reused."""

    reusable: bool
    reason: str
    stage: str
    scope: Optional[str]
    recorded_signature: Optional[str] = None
    requested_signature: Optional[str] = None
    differences: tuple[dict[str, Any], ...] = ()

    def describe(self) -> str:
        label = f"{self.stage}:{self.scope}" if self.scope else self.stage
        if self.reason != "signature_mismatch":
            return f"Checkpoint {label} is not reusable: {self.reason}"
        lines = [
            f"Checkpoint {label} signature mismatch",
            f"recorded={self.recorded_signature}",
            f"current={self.requested_signature}",
        ]
        for row in self.differences[:20]:
            lines.append(
                f"{row['path']}: {row.get('recorded')!r} -> {row.get('current')!r}"
            )
        if not self.differences:
            lines.append("Recorded signature payload is unavailable (legacy checkpoint).")
        return "; ".join(lines)


class V4Checkpoint:
    """Single stage- and variant-aware checkpoint stored in ``pipeline.state``."""

    filename = "pipeline.state"
    stage_order = ("s0", "s1", "s2", "s3", "s4")

    def __init__(self, work_dir: Path):
        self.work_dir = Path(work_dir)
        self.path = self.work_dir / self.filename
        self.legacy_variant_path = self.work_dir / ".rph" / "checkpoint.json"
        self.run_id: Optional[str] = None

    @staticmethod
    def _empty_state() -> Dict[str, Any]:
        return {"schema_version": STATE_SCHEMA_V2, RUN_ID_FIELD: None, "stages": {}}

    def load(self) -> Dict[str, Any]:
        changed = False
        if not self.path.exists():
            data = self._empty_state()
        else:
            data = json.loads(self.path.read_text(encoding="utf-8"))
            schema = data.get("schema_version")
            if schema == STATE_SCHEMA_V1:
                data = {
                    **data,
                    "schema_version": STATE_SCHEMA_V2,
                    "stages": dict(data.get("stages", {}) or {}),
                }
                changed = True
            elif schema != STATE_SCHEMA_V2:
                raise RuntimeError(
                    f"Unsupported RPH V4 pipeline state schema: {schema!r}"
                )

        if RUN_ID_FIELD not in data:
            data[RUN_ID_FIELD] = None
            changed = True
        elif data.get(RUN_ID_FIELD) not in (None, "legacy") and not is_valid_run_id(
            data.get(RUN_ID_FIELD)
        ):
            logger.warning(
                "Invalid checkpoint run_id %r in %s; normalizing to None",
                data.get(RUN_ID_FIELD),
                self.path,
            )
            data[RUN_ID_FIELD] = None
            changed = True

        if not data.get("legacy_variant_checkpoint_imported"):
            changed = self._import_legacy_variant_state(data) or changed
            data["legacy_variant_checkpoint_imported"] = True
            changed = True

        self.run_id = data.get(RUN_ID_FIELD)
        if changed:
            self.save(data)
        return data

    def _import_legacy_variant_state(self, state: Dict[str, Any]) -> bool:
        """Import the old .rph checkpoint without modifying or deleting it."""

        if not self.legacy_variant_path.exists():
            return False
        legacy = json.loads(self.legacy_variant_path.read_text(encoding="utf-8"))
        if legacy.get("schema_version") != VARIANT_SCHEMA_V1:
            return False

        changed = False
        if legacy.get("reaction_id") and not state.get("reaction_id"):
            state["reaction_id"] = legacy["reaction_id"]
            changed = True

        for scope in ("s0", "s3", "s4"):
            entry = legacy.get(scope)
            if isinstance(entry, dict):
                changed = self._merge_stage_entry(state, scope, entry) or changed

        precursor = legacy.get("precursor_s1")
        if isinstance(precursor, dict):
            changed = self._merge_scope_entry(
                state, "s1", "precursor_s1", precursor
            ) or changed

        variants = legacy.get("variants", {}) or {}
        if isinstance(variants, dict):
            for scope, stage_map in variants.items():
                if not isinstance(stage_map, dict):
                    continue
                for stage, entry in stage_map.items():
                    if stage in self.stage_order and isinstance(entry, dict):
                        changed = self._merge_scope_entry(
                            state, stage, str(scope), entry
                        ) or changed
        return changed

    @staticmethod
    def _merge_stage_entry(
        state: Dict[str, Any], stage: str, entry: Mapping[str, Any]
    ) -> bool:
        stages = state.setdefault("stages", {})
        current = stages.setdefault(stage, {})
        if current.get("signature"):
            return False
        current.update(dict(entry))
        return True

    @staticmethod
    def _merge_scope_entry(
        state: Dict[str, Any], stage: str, scope: str, entry: Mapping[str, Any]
    ) -> bool:
        stages = state.setdefault("stages", {})
        stage_record = stages.setdefault(stage, {})
        scopes = stage_record.setdefault("scopes", {})
        if scope in scopes:
            return False
        scopes[scope] = dict(entry)
        return True

    def save(self, data: Dict[str, Any]) -> None:
        self.work_dir.mkdir(parents=True, exist_ok=True)
        data["schema_version"] = STATE_SCHEMA_V2
        data[RUN_ID_FIELD] = self.run_id
        write_json_atomic(self.path, data)

    def write_config_snapshot(self, config: Mapping[str, Any]) -> Path:
        """Write the immutable configuration used to start a result directory."""

        path = self.work_dir / "run.config.json"
        if path.exists():
            return path
        self.work_dir.mkdir(parents=True, exist_ok=True)
        write_json_atomic(
            path,
            {
                "schema_version": "rph_v4_run_config_v1",
                "config": dict(config),
            },
        )
        return path

    def set_run_id(self, run_id: str) -> None:
        if not is_valid_run_id(run_id):
            raise ValueError(f"Invalid run_id: {run_id!r}")
        state = self.load()
        if state.get(RUN_ID_FIELD) == run_id:
            self.run_id = run_id
            return
        if state.get(RUN_ID_FIELD) in (None, "legacy") and (
            state.get("reaction_id") or state.get("stages")
        ):
            logger.warning(
                "Resuming across runs (old run_id=%s, new run_id=%s)",
                state.get(RUN_ID_FIELD),
                run_id,
            )
        state[RUN_ID_FIELD] = run_id
        self.run_id = run_id
        self.save(state)

    def set_reaction_id(self, reaction_id: str) -> None:
        state = self.load()
        if state.get("reaction_id") == reaction_id:
            return
        state["reaction_id"] = reaction_id
        self.save(state)

    @staticmethod
    def signature(payload: Mapping[str, Any]) -> str:
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

    def check_reuse(
        self,
        stage: str,
        signature: str,
        manifest: Path,
        *,
        scope: Optional[str] = None,
        signature_payload: Optional[Mapping[str, Any]] = None,
        accept_signature_mismatch: bool = False,
    ) -> ReuseDecision:
        entry = self._get_entry(self.load(), stage, scope)
        if not entry:
            return ReuseDecision(False, "missing_entry", stage, scope)
        if entry.get("status") != "complete":
            return ReuseDecision(
                False,
                f"status_{entry.get('status', 'unknown')}",
                stage,
                scope,
            )

        requested_manifest = Path(manifest)
        recorded_manifest = entry.get("manifest")
        if not isinstance(recorded_manifest, str):
            return ReuseDecision(False, "manifest_unrecorded", stage, scope)
        if Path(recorded_manifest).resolve() != requested_manifest.resolve():
            return ReuseDecision(False, "manifest_mismatch", stage, scope)
        if not requested_manifest.exists():
            return ReuseDecision(False, "manifest_missing", stage, scope)

        recorded_signature = entry.get("signature")
        if recorded_signature != signature:
            differences = tuple(
                self.payload_differences(
                    entry.get("signature_payload"), signature_payload
                )
            )
            if accept_signature_mismatch:
                return ReuseDecision(
                    True,
                    "signature_mismatch_accepted",
                    stage,
                    scope,
                    str(recorded_signature) if recorded_signature else None,
                    signature,
                    differences,
                )
            return ReuseDecision(
                False,
                "signature_mismatch",
                stage,
                scope,
                str(recorded_signature) if recorded_signature else None,
                signature,
                differences,
            )
        return ReuseDecision(
            True,
            "reusable",
            stage,
            scope,
            str(recorded_signature) if recorded_signature else None,
            signature,
        )

    def reusable(self, stage: str, signature: str, manifest: Path) -> bool:
        return self.check_reuse(stage, signature, manifest).reusable

    def reusable_scope(
        self, scope: str, stage: str, signature: str, manifest: Path
    ) -> bool:
        return self.check_reuse(
            stage, signature, manifest, scope=self._normalized_scope(scope, stage)
        ).reusable

    def has_stage(self, stage_name: str) -> bool:
        state = self.load()
        stages = state.get("stages", {}) or {}
        return stage_name in stages

    @staticmethod
    def payload_differences(
        recorded: Any, current: Any
    ) -> Sequence[dict[str, Any]]:
        if not isinstance(recorded, Mapping) or not isinstance(current, Mapping):
            return ()

        rows: list[dict[str, Any]] = []

        def walk(old: Any, new: Any, prefix: str) -> None:
            if isinstance(old, Mapping) and isinstance(new, Mapping):
                for key in sorted(set(old) | set(new), key=str):
                    path = f"{prefix}.{key}" if prefix else str(key)
                    walk(old.get(key), new.get(key), path)
                return
            if old != new:
                rows.append({"path": prefix, "recorded": old, "current": new})

        walk(recorded, current, "")
        return rows

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

    def mark_recovered(self, stage_name: str, *, status: str) -> None:
        state = self.load()
        stages = state.setdefault("stages", {})
        entry = stages.get(stage_name)
        if not isinstance(entry, dict):
            logger.warning(
                "Cannot mark recovered stage %s in %s because no checkpoint entry exists",
                stage_name,
                self.path,
            )
            return
        if (
            entry.get("status") == status
            and entry.get("recovery_status") == status
            and isinstance(entry.get("recovered_at"), str)
            and entry.get("recovered_at")
        ):
            return
        entry["status"] = status
        entry["recovery_status"] = status
        entry["recovered_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        self.save(state)

    def mark(
        self,
        stage: str,
        signature: str,
        manifest: Path,
        status: str = "complete",
        *,
        signature_payload: Optional[Mapping[str, Any]] = None,
    ) -> None:
        if stage not in self.stage_order:
            raise ValueError(f"Unknown V4 stage: {stage}")
        state = self.load()
        stages = state.setdefault("stages", {})
        stage_record = stages.setdefault(stage, {})
        scopes = stage_record.get("scopes")
        stage_record.clear()
        if isinstance(scopes, dict):
            stage_record["scopes"] = scopes
        stage_record.update(
            self._entry(signature, manifest, status, signature_payload)
        )
        self._invalidate_downstream(stages, stage)
        self.save(state)

    def mark_scope(
        self,
        scope: str,
        stage: str,
        signature: str,
        manifest: Path,
        status: str = "complete",
        *,
        signature_payload: Optional[Mapping[str, Any]] = None,
    ) -> None:
        normalized_scope = self._normalized_scope(scope, stage)
        if normalized_scope is None:
            if stage == "s0":
                self.mark_s0(
                    signature,
                    manifest,
                    status,
                    signature_payload=signature_payload,
                )
            else:
                self.mark(
                    stage,
                    signature,
                    manifest,
                    status,
                    signature_payload=signature_payload,
                )
            return

        state = self.load()
        stages = state.setdefault("stages", {})
        stage_record = stages.setdefault(stage, {})
        for key in ("signature", "manifest", "status", "signature_payload"):
            stage_record.pop(key, None)
        scopes = stage_record.setdefault("scopes", {})
        scopes[normalized_scope] = self._entry(
            signature, manifest, status, signature_payload
        )
        self._invalidate_downstream(stages, stage)
        self.save(state)

    def mark_failed(
        self,
        scope: str,
        stage: str,
        signature: str,
        error: str,
        *,
        signature_payload: Optional[Mapping[str, Any]] = None,
    ) -> None:
        normalized_scope = self._normalized_scope(scope, stage)
        state = self.load()
        stages = state.setdefault("stages", {})
        stage_record = stages.setdefault(stage, {})
        for key in ("signature", "manifest", "status", "signature_payload"):
            stage_record.pop(key, None)
        value: Dict[str, Any] = {
            "signature": signature,
            "status": "failed",
            "error": error,
        }
        payload = self._signature_payload_with_run_id(signature_payload)
        if payload is not None:
            value["signature_payload"] = payload
        if normalized_scope is None:
            stage_record.update(value)
        else:
            stage_record.setdefault("scopes", {})[normalized_scope] = value
        self._invalidate_downstream(stages, stage)
        self.save(state)

    def mark_s0(
        self,
        signature: str,
        manifest: Path,
        status: str = "complete",
        *,
        signature_payload: Optional[Mapping[str, Any]] = None,
    ) -> None:
        """Mark S0 while preserving independent S1 and invalidating S2+."""

        state = self.load()
        stages = state.setdefault("stages", {})
        stages["s0"] = self._entry(
            signature, manifest, status, signature_payload
        )
        for downstream in ("s2", "s3", "s4"):
            stages.pop(downstream, None)
        self.save(state)

    def _entry(
        self,
        signature: str,
        manifest: Path,
        status: str,
        signature_payload: Optional[Mapping[str, Any]],
    ) -> Dict[str, Any]:
        value: Dict[str, Any] = {
            "signature": signature,
            "manifest": str(Path(manifest).resolve()),
            "status": status,
        }
        payload = self._signature_payload_with_run_id(signature_payload)
        if payload is not None:
            value["signature_payload"] = payload
        return value

    def _signature_payload_with_run_id(
        self,
        signature_payload: Optional[Mapping[str, Any]],
    ) -> Optional[Dict[str, Any]]:
        if signature_payload is None and self.run_id is None:
            return None
        payload = dict(signature_payload or {})
        payload[RUN_ID_FIELD] = self.run_id
        return payload

    def _get_entry(
        self, state: Mapping[str, Any], stage: str, scope: Optional[str]
    ) -> Optional[Dict[str, Any]]:
        stage_record = (state.get("stages", {}) or {}).get(stage, {})
        if not isinstance(stage_record, dict):
            return None
        if scope is None:
            return stage_record if stage_record.get("signature") else None
        scopes = stage_record.get("scopes", {}) or {}
        entry = scopes.get(scope) if isinstance(scopes, dict) else None
        return entry if isinstance(entry, dict) else None

    @staticmethod
    def _normalized_scope(scope: str, stage: str) -> Optional[str]:
        if scope in ("s0", "s3", "s4") and scope == stage:
            return None
        return scope

    def _invalidate_downstream(self, stages: Dict[str, Any], stage: str) -> None:
        start = self.stage_order.index(stage) + 1
        for downstream in self.stage_order[start:]:
            stages.pop(downstream, None)
