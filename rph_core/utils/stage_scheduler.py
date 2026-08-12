"""Bounded, work-conserving structure scheduling for V4 QC stages."""

from __future__ import annotations

import logging
import threading
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor, as_completed
from copy import deepcopy
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Callable, Dict, Iterable, List, TypeVar

from rph_core.utils.resource_utils import mem_to_mb


T = TypeVar("T")
logger = logging.getLogger(__name__)


def _empty_role_limits() -> Mapping[str, int]:
    return MappingProxyType({})


def _structure_role(structure: Mapping[str, Any]) -> str:
    return str(structure.get("role") or structure.get("kind", "product")).lower()


@dataclass(frozen=True)
class StageSchedule:
    enabled: bool
    max_workers: int
    nproc_per_job: int
    memory_per_job: str
    priority_roles: tuple[str, ...]
    max_concurrent_by_role: Mapping[str, int] = field(default_factory=_empty_role_limits)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "max_concurrent_by_role",
            MappingProxyType(
                {
                    str(role): int(limit)
                    for role, limit in self.max_concurrent_by_role.items()
                }
            ),
        )

    def manifest_record(self) -> Dict[str, Any]:
        return {
            "enabled": self.enabled,
            "max_workers": self.max_workers,
            "nproc_per_job": self.nproc_per_job,
            "memory_per_job": self.memory_per_job,
            "priority_roles": list(self.priority_roles),
            "max_concurrent_by_role": dict(self.max_concurrent_by_role),
        }


def resolve_stage_schedule(config: Dict[str, Any], stage_key: str) -> StageSchedule:
    """Validate one stage's worker allocation against the global resource pool."""

    resources = dict(config.get("resources", {}) or {})
    stage_cfg = dict(config.get(stage_key, {}) or {})
    scheduling = dict(stage_cfg.get("scheduling", {}) or {})
    enabled = bool(scheduling.get("enabled", False))
    max_workers = int(scheduling.get("max_workers", 1)) if enabled else 1
    nproc_per_job = int(scheduling.get("nproc_per_job", resources.get("nproc", 1)))
    memory_per_job = str(scheduling.get("memory_per_job", resources.get("mem", "1GB")))
    priority_roles = tuple(
        str(role).strip().lower()
        for role in scheduling.get(
            "priority_roles",
            ("ts", "intermediate", "precursor", "product"),
        )
    )
    raw_role_limits = scheduling.get("max_concurrent_by_role", {}) or {}
    if max_workers <= 0 or nproc_per_job <= 0:
        raise ValueError(f"{stage_key}.scheduling worker and CPU counts must be positive")
    if not isinstance(raw_role_limits, Mapping):
        raise ValueError(
            f"{stage_key}.scheduling.max_concurrent_by_role must be a mapping"
        )
    total_nproc = int(resources.get("nproc", 1))
    if max_workers * nproc_per_job > total_nproc:
        raise ValueError(
            f"{stage_key}.scheduling requests {max_workers * nproc_per_job} cores "
            f"from a {total_nproc}-core resource pool"
        )
    total_memory_mb = mem_to_mb(str(resources.get("mem", "1GB")))
    requested_memory_mb = max_workers * mem_to_mb(memory_per_job)
    if requested_memory_mb > total_memory_mb:
        raise ValueError(
            f"{stage_key}.scheduling requests {requested_memory_mb} MB from a "
            f"{total_memory_mb} MB resource pool"
        )
    max_concurrent_by_role: Dict[str, int] = {}
    for raw_role, raw_limit in raw_role_limits.items():
        role = str(raw_role).strip()
        if not role or role != role.lower():
            raise ValueError(
                f"{stage_key}.scheduling.max_concurrent_by_role keys must be lowercase"
            )
        if not isinstance(raw_limit, int) or isinstance(raw_limit, bool) or raw_limit < 1:
            raise ValueError(
                f"{stage_key}.scheduling.max_concurrent_by_role[{role!r}] must be an integer >= 1"
            )
        limit = raw_limit
        if limit > max_workers:
            logger.warning(
                "%s.scheduling.max_concurrent_by_role[%r]=%d exceeds max_workers=%d; clamping",
                stage_key,
                role,
                limit,
                max_workers,
            )
            limit = max_workers
        if nproc_per_job * limit > total_nproc:
            raise ValueError(
                f"{stage_key}.scheduling.max_concurrent_by_role[{role!r}] requests "
                f"{nproc_per_job * limit} cores from a {total_nproc}-core resource pool"
            )
        max_concurrent_by_role[role] = limit
    return StageSchedule(
        enabled=enabled,
        max_workers=max_workers,
        nproc_per_job=nproc_per_job,
        memory_per_job=memory_per_job,
        priority_roles=priority_roles,
        max_concurrent_by_role=max_concurrent_by_role,
    )


def worker_config(config: Dict[str, Any], schedule: StageSchedule) -> Dict[str, Any]:
    """Return an isolated config carrying one worker's CPU and memory budget."""

    derived = deepcopy(config)
    resources = dict(derived.get("resources", {}) or {})
    resources.update(
        {
            "nproc": schedule.nproc_per_job,
            "mem": schedule.memory_per_job,
        }
    )
    derived["resources"] = resources
    return derived


def worker_theory(theory: Dict[str, Any], schedule: StageSchedule) -> Dict[str, Any]:
    """Apply the worker budget to every QC job in a stage theory block."""

    derived = deepcopy(theory)
    for section_name in ("optimization", "single_point"):
        section = dict(derived.get(section_name, {}) or {})
        section["nproc"] = schedule.nproc_per_job
        section["mem"] = schedule.memory_per_job
        derived[section_name] = section
    return derived


def run_structure_queue(
    structures: Iterable[Dict[str, Any]],
    worker: Callable[[Dict[str, Any]], T],
    schedule: StageSchedule,
    *,
    thread_name_prefix: str,
) -> List[T]:
    """Run a priority-ordered queue and return results in original input order."""

    materialized = list(structures)
    if not materialized:
        return []
    role_priority = {role: index for index, role in enumerate(schedule.priority_roles)}

    def priority(item: tuple[int, Dict[str, Any]]) -> tuple[int, int]:
        index, structure = item
        role = _structure_role(structure)
        return role_priority.get(role, len(role_priority)), index

    queued = sorted(enumerate(materialized), key=priority)
    results: List[Any] = [None] * len(materialized)
    if schedule.max_workers == 1:
        for index, structure in queued:
            results[index] = worker(structure)
        return results
    role_semaphores = {
        role: threading.BoundedSemaphore(limit)
        for role, limit in schedule.max_concurrent_by_role.items()
    }
    max_workers = min(schedule.max_workers, len(materialized))
    if not role_semaphores:
        with ThreadPoolExecutor(
            max_workers=max_workers,
            thread_name_prefix=thread_name_prefix,
        ) as executor:
            futures = {
                executor.submit(worker, structure): index
                for index, structure in queued
            }
            for future in as_completed(futures):
                results[futures[future]] = future.result()
        return results

    def run_with_role_limit(
        structure: Dict[str, Any],
        role_semaphore: threading.BoundedSemaphore | None,
    ) -> T:
        try:
            return worker(structure)
        finally:
            if role_semaphore is not None:
                role_semaphore.release()

    with ThreadPoolExecutor(
        max_workers=max_workers,
        thread_name_prefix=thread_name_prefix,
    ) as executor:
        pending = list(queued)
        futures: Dict[Any, int] = {}

        def submit_ready() -> None:
            while len(futures) < max_workers and pending:
                pending_position: int | None = None
                acquired_semaphore: threading.BoundedSemaphore | None = None
                for position, (_, structure) in enumerate(pending):
                    role = _structure_role(structure)
                    semaphore = role_semaphores.get(role)
                    if semaphore is not None and not semaphore.acquire(blocking=False):
                        continue
                    pending_position = position
                    acquired_semaphore = semaphore
                    break
                if pending_position is None:
                    return
                index, structure = pending.pop(pending_position)
                try:
                    future = executor.submit(run_with_role_limit, structure, acquired_semaphore)
                except Exception:
                    if acquired_semaphore is not None:
                        acquired_semaphore.release()
                    raise
                futures[future] = index

        submit_ready()
        while futures:
            future = next(as_completed(tuple(futures)))
            results[futures.pop(future)] = future.result()
            submit_ready()
        if pending:
            raise RuntimeError("Role-aware scheduler stalled with pending structures")
    return results
