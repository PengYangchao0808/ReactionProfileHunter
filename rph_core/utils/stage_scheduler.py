"""Bounded, work-conserving structure scheduling for V4 QC stages."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from copy import deepcopy
from dataclasses import asdict, dataclass
from typing import Any, Callable, Dict, Iterable, List, Sequence, TypeVar

from rph_core.utils.resource_utils import mem_to_mb


T = TypeVar("T")


@dataclass(frozen=True)
class StageSchedule:
    enabled: bool
    max_workers: int
    nproc_per_job: int
    memory_per_job: str
    priority_roles: tuple[str, ...]

    def manifest_record(self) -> Dict[str, Any]:
        return asdict(self)


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
    if max_workers <= 0 or nproc_per_job <= 0:
        raise ValueError(f"{stage_key}.scheduling worker and CPU counts must be positive")
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
    return StageSchedule(
        enabled=enabled,
        max_workers=max_workers,
        nproc_per_job=nproc_per_job,
        memory_per_job=memory_per_job,
        priority_roles=priority_roles,
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
        role = str(structure.get("role") or structure.get("kind", "product")).lower()
        return role_priority.get(role, len(role_priority)), index

    queued = sorted(enumerate(materialized), key=priority)
    results: List[Any] = [None] * len(materialized)
    if schedule.max_workers == 1:
        for index, structure in queued:
            results[index] = worker(structure)
        return results

    with ThreadPoolExecutor(
        max_workers=min(schedule.max_workers, len(materialized)),
        thread_name_prefix=thread_name_prefix,
    ) as executor:
        futures = {
            executor.submit(worker, structure): index
            for index, structure in queued
        }
        for future in as_completed(futures):
            results[futures[future]] = future.result()
    return results
