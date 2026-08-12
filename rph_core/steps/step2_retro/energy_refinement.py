"""Optional B97-3c single-point refinement for an xTB scan trajectory."""

from __future__ import annotations

import hashlib
import json
import logging
import threading
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.qc_jobs import run_single_point
from rph_core.utils.qc_models import QCJobSpec
from rph_core.utils.resource_utils import mem_to_mb
from rph_core.utils.stage_scheduler import StageSchedule, run_structure_queue


logger = logging.getLogger(__name__)


class ScanEnergyRefiner(LoggerMixin):
    """Run resumable ORCA SP jobs without changing scan geometries."""

    CACHE_SCHEMA = "s2_geometry_sp_cache_v2"

    def __init__(
        self,
        config: Dict[str, Any],
        output_dir: Path,
        *,
        event_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
        variant: Optional[str] = None,
    ):
        self.config = config
        self.output_dir = Path(output_dir)
        self.cfg = dict((config.get("step2", {}) or {}).get("energy_refinement", {}) or {})
        self.event_callback = event_callback
        self.variant = variant or "product"
        self._legacy_cache_index: Dict[tuple[str, str], Dict[str, Any]] = {}
        self._legacy_index_spec_signature = ""
        self._progress_lock = threading.Lock()
        self._progress_done = 0
        self._progress_failed = 0
        self._progress_total = 0
        self._progress_workers = 1
        self._batch_started_monotonic = 0.0
        self._batch_id = f"{self.variant}:b973c_sp"

    def _emit(self, event: str, **fields: Any) -> None:
        if self.event_callback is None:
            return
        try:
            self.event_callback(event, {"variant": self.variant, **fields})
        except Exception as exc:  # pragma: no cover - UI isolation
            self.logger.warning("[S2] Ignoring B97-3c UI callback failure for %s: %s", event, exc)

    @property
    def enabled(self) -> bool:
        return bool(self.cfg.get("enabled", False))

    def _spec(self) -> QCJobSpec:
        resources = dict(self.config.get("resources", {}) or {})
        cores_per_job = max(1, int(self.cfg.get("cores_per_job", 1)))
        memory_per_job = str(self.cfg.get("memory_per_job", "2GB"))
        configured_maxcore = self.cfg.get("maxcore")
        maxcore = (
            int(configured_maxcore)
            if configured_maxcore is not None
            else max(
                256,
                int(
                    mem_to_mb(memory_per_job)
                    * float(resources.get("orca_maxcore_safety", 0.65))
                    / cores_per_job
                ),
            )
        )
        return QCJobSpec(
            engine=str(self.cfg.get("engine", "orca")),
            task="sp",
            method=str(self.cfg.get("method", "B97-3c")),
            basis=str(self.cfg.get("basis", "")),
            aux_basis=str(self.cfg.get("aux_basis", "")),
            solvent=str(self.cfg.get("solvent", "acetone")),
            solvent_model=str(self.cfg.get("solvent_model", "CPCM")),
            route_extras=str(self.cfg.get("route_extras", "")),
            charge=int(self.cfg.get("charge", 0)),
            multiplicity=int(self.cfg.get("multiplicity", 1)),
            nproc=min(cores_per_job, max(1, int(resources.get("nproc", 1)))),
            memory=memory_per_job,
            maxcore=maxcore,
            timeout=int(self.cfg["timeout"]) if self.cfg.get("timeout") is not None else None,
        )

    @staticmethod
    def _file_hash(path: Path) -> str:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    @staticmethod
    def _signature(spec: QCJobSpec) -> str:
        payload = json.dumps(asdict(spec), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    @staticmethod
    def _write_json_atomic(path: Path, payload: Dict[str, Any]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(f".{path.name}.tmp")
        temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        temporary.replace(path)

    def _index_legacy_frame_cache(self, spec_signature: str) -> None:
        """Index v1 frame-number caches by scientific identity for migration.

        Rebuilt only when the QC spec signature changes; repeated ``refine()``
        calls for the same spec reuse the existing index instead of re-globbing.
        """

        if (
            self._legacy_index_spec_signature == spec_signature
            and self._legacy_cache_index
        ):
            return
        self._legacy_cache_index = {}
        for cache_path in sorted(self.output_dir.glob("frame_*/result.json")):
            try:
                record = json.loads(cache_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            xyz_hash = record.get("xyz_sha256")
            cached_spec = record.get("spec_sha256")
            if (
                xyz_hash
                and cached_spec == spec_signature
                and record.get("status") == "complete"
                and record.get("energy_hartree") is not None
            ):
                record["legacy_cache_path"] = str(cache_path)
                self._legacy_cache_index[(str(xyz_hash), str(cached_spec))] = record
        self._legacy_index_spec_signature = spec_signature

    def _run_one_impl(
        self, item: Dict[str, Any], spec: QCJobSpec, spec_signature: str
    ) -> Dict[str, Any]:
        index = int(item["index"])
        frame = Path(item["frame"])
        xyz_hash_value = item.get("xyz_sha256")
        xyz_hash = str(xyz_hash_value) if xyz_hash_value is not None else self._file_hash(frame)
        job_dir = self.output_dir / "cache" / f"{xyz_hash[:20]}_{spec_signature[:12]}"
        cache_path = job_dir / "result.json"
        if cache_path.is_file():
            try:
                cached = json.loads(cache_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                cached = {}
            if (
                cached.get("schema_version") == self.CACHE_SCHEMA
                and cached.get("xyz_sha256") == xyz_hash
                and cached.get("spec_sha256") == spec_signature
                and cached.get("status") == "complete"
                and cached.get("energy_hartree") is not None
            ):
                cached["reused"] = True
                cached["frame_index"] = index
                cached["point_id"] = item.get("point_id")
                cached["input_xyz"] = str(frame)
                return cached

        legacy = self._legacy_cache_index.get((xyz_hash, spec_signature))
        if legacy:
            migrated = dict(legacy)
            migrated.update(
                {
                    "schema_version": self.CACHE_SCHEMA,
                    "frame_index": index,
                    "point_id": item.get("point_id"),
                    "input_xyz": str(frame),
                    "reused": True,
                    "migrated_from": legacy.get("legacy_cache_path"),
                }
            )
            self._write_json_atomic(cache_path, migrated)
            return migrated

        result = run_single_point(spec, frame, job_dir, self.config)
        record: Dict[str, Any] = {
            "schema_version": self.CACHE_SCHEMA,
            "frame_index": index,
            "point_id": item.get("point_id"),
            "input_xyz": str(frame),
            "xyz_sha256": xyz_hash,
            "spec_sha256": spec_signature,
            "method": spec.method,
            "solvent": spec.solvent,
            "solvent_model": spec.solvent_model,
            "status": result.status,
            "energy_hartree": result.energy_hartree,
            "output_file": str(result.output_file) if result.output_file else None,
            "error": result.error,
            "reused": False,
        }
        self._write_json_atomic(cache_path, record)
        return record

    def _run_one(
        self, item: Dict[str, Any], spec: QCJobSpec, spec_signature: str
    ) -> Dict[str, Any]:
        point_id = str(item.get("point_id") or f"frame_{int(item['index']):04d}")
        job_id = f"{self.variant}:{point_id}"
        xyz_hash_value = item.get("xyz_sha256")
        xyz_hash = (
            str(xyz_hash_value)
            if xyz_hash_value is not None
            else self._file_hash(Path(item["frame"]))
        )
        job_dir = self.output_dir / "cache" / f"{xyz_hash[:20]}_{spec_signature[:12]}"
        started = time.monotonic()
        self._emit(
            "batch_job_started",
            batch=self._batch_id,
            job_id=job_id,
            point_id=point_id,
            engine=spec.engine,
            method=spec.method,
            nprocs=int(spec.nproc or 1),
            output=str(job_dir),
            started_at=time.time(),
        )
        try:
            record = self._run_one_impl(item, spec, spec_signature)
        except Exception as exc:
            elapsed = time.monotonic() - started
            self._emit(
                "batch_job_failed",
                batch=self._batch_id,
                job_id=job_id,
                point_id=point_id,
                engine=spec.engine,
                method=spec.method,
                nprocs=int(spec.nproc or 1),
                output=str(job_dir),
                elapsed_seconds=elapsed,
                error=str(exc),
            )
            self._update_progress(point_id=point_id, failed=True)
            raise

        failed = record.get("status") != "complete" or record.get("energy_hartree") is None
        elapsed = time.monotonic() - started
        self._emit(
            "batch_job_failed" if failed else "batch_job_finished",
            batch=self._batch_id,
            job_id=job_id,
            point_id=point_id,
            engine=spec.engine,
            method=spec.method,
            nprocs=int(spec.nproc or 1),
            output=record.get("output_file") or str(job_dir),
            elapsed_seconds=elapsed,
            status="failed" if failed else "cached" if record.get("reused") else "complete",
            energy_hartree=record.get("energy_hartree"),
            error=record.get("error"),
            reused=bool(record.get("reused")),
        )
        self._update_progress(point_id=point_id, failed=failed)
        return record

    def _update_progress(self, *, point_id: str, failed: bool) -> None:
        with self._progress_lock:
            if failed:
                self._progress_failed += 1
            else:
                self._progress_done += 1
            done = self._progress_done
            failed_count = self._progress_failed
            processed = done + failed_count
            total = self._progress_total
            elapsed = max(0.0, time.monotonic() - self._batch_started_monotonic)
            rate = processed * 60.0 / elapsed if elapsed > 0.0 else None
            eta = (
                max(0.0, (total - processed) * elapsed / processed)
                if processed > 0 and processed < total
                else 0.0
            )
            self._emit(
                "batch_progress",
                batch=self._batch_id,
                label="ORCA B97-3c SP refinement",
                phase="b973c_sp",
                total=total,
                done=done,
                failed=failed_count,
                running=max(0, min(total - processed, self._progress_workers)),
                current=point_id,
                elapsed_seconds=elapsed,
                rate_per_minute=rate,
                eta_seconds=eta,
            )

    def refine(
        self,
        frame_paths: Sequence[Path],
        *,
        point_ids: Optional[Sequence[str]] = None,
    ) -> Dict[str, Any]:
        if point_ids is not None and len(point_ids) != len(frame_paths):
            raise ValueError("S2 energy refinement requires one point_id per frame")
        if not self.enabled:
            return {
                "enabled": False,
                "status": "not_requested",
                "energies_hartree": [None for _ in frame_paths],
                "records": [],
            }

        spec = self._spec()
        total_cores = max(1, int((self.config.get("resources", {}) or {}).get("nproc", 1)))
        total_memory_mb = mem_to_mb(
            str((self.config.get("resources", {}) or {}).get("mem", "1GB"))
        )
        requested_jobs = max(1, int(self.cfg.get("parallel_jobs", total_cores)))
        memory_capacity = max(1, total_memory_mb // mem_to_mb(str(spec.memory or "1GB")))
        max_workers = max(
            1,
            min(
                requested_jobs,
                total_cores // int(spec.nproc or 1),
                memory_capacity,
            ),
        )
        schedule = StageSchedule(
            enabled=max_workers > 1,
            max_workers=max_workers,
            nproc_per_job=int(spec.nproc or 1),
            memory_per_job=str(spec.memory or "1GB"),
            priority_roles=("frame",),
        )
        spec_signature = self._signature(spec)
        self._index_legacy_frame_cache(spec_signature)
        items = [
            {
                "index": index,
                "point_id": point_ids[index] if point_ids is not None else f"frame_{index:04d}",
                "frame": str(Path(frame)),
                "xyz_sha256": self._file_hash(Path(frame)),
                "role": "frame",
            }
            for index, frame in enumerate(frame_paths)
        ]
        self._progress_done = 0
        self._progress_failed = 0
        self._progress_total = len(items)
        self._progress_workers = max_workers
        self._batch_started_monotonic = time.monotonic()
        self._emit(
            "batch_started",
            batch=self._batch_id,
            label="ORCA B97-3c SP refinement",
            phase="b973c_sp",
            total=len(items),
            done=0,
            failed=0,
            running=min(max_workers, len(items)),
            started_at=time.time(),
            parallel_jobs=max_workers,
            cores_per_job=int(spec.nproc or 1),
            engine=spec.engine,
            method=spec.method,
        )
        try:
            records = run_structure_queue(
                items,
                lambda item: self._run_one(item, spec, spec_signature),
                schedule,
                thread_name_prefix="rph-s2-sp",
            )
        except Exception as exc:
            self._emit(
                "batch_finished",
                batch=self._batch_id,
                label="ORCA B97-3c SP refinement",
                phase="b973c_sp",
                total=len(items),
                done=self._progress_done,
                failed=max(1, self._progress_failed),
                running=0,
                status="failed",
                elapsed_seconds=time.monotonic() - self._batch_started_monotonic,
                error=str(exc),
            )
            raise
        energies: List[Optional[float]] = [None] * len(frame_paths)
        for record in records:
            if record.get("status") == "complete" and record.get("energy_hartree") is not None:
                energies[int(record["frame_index"])] = float(record["energy_hartree"])
        completed = sum(value is not None for value in energies)
        final_status = "complete" if completed == len(frame_paths) else "partial"
        self._emit(
            "batch_finished",
            batch=self._batch_id,
            label="ORCA B97-3c SP refinement",
            phase="b973c_sp",
            total=len(items),
            done=completed,
            failed=len(items) - completed,
            running=0,
            status=final_status,
            elapsed_seconds=time.monotonic() - self._batch_started_monotonic,
            method=spec.method,
        )
        return {
            "enabled": True,
            "status": final_status,
            "method": spec.method,
            "solvent": spec.solvent,
            "solvent_model": spec.solvent_model,
            "completed": completed,
            "total": len(frame_paths),
            "parallel_jobs": max_workers,
            "cores_per_job": int(spec.nproc or 1),
            "energies_hartree": energies,
            "records": records,
        }
