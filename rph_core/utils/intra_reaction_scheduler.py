"""
Intra-Reaction Parallel Scheduler (单反应内部并行调度器)
======================================================

Manages concurrent QC jobs within a single reaction, respecting
resource constraints (CPU cores, memory) per lane.

Config section: intra_reaction_parallel (in defaults.yaml)
Default: enabled=false — requires explicit user opt-in.

Author: QCcalc Team
Date: 2026-01-20
"""

from __future__ import annotations

import logging
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Callable, cast

logger = logging.getLogger(__name__)

ConfigDict = dict[str, object]
JobCallable = Callable[..., object]


@dataclass
class ResourceLane:
    """A resource allocation lane for a parallel QC job.

    Attributes:
        nproc: Number of CPU cores allocated
        mem: Memory string (e.g., "24GB")
        label: Human-readable label for logging
    """

    nproc: int = 8
    mem: str = "24GB"
    label: str = "lane_default"


@dataclass
class ParallelQCJob:
    """A single QC job to be submitted to the scheduler.

    Attributes:
        job_id: Unique identifier for logging
        lane: Resource allocation for this job
        func: Callable to execute (typically a QCTaskRunner method call)
        args: Positional arguments for func
        kwargs: Keyword arguments for func
        description: Human-readable description for logging
    """

    job_id: str
    lane: ResourceLane
    func: JobCallable
    args: tuple[object, ...] = ()
    kwargs: ConfigDict = field(default_factory=dict)
    description: str = ""


@dataclass
class ParallelResult:
    """Result from a parallel QC job execution.

    Attributes:
        job_id: Matching job identifier
        success: Whether the job completed without exception
        result: Return value from the job function (if success)
        error: Exception message (if not success)
        lane_label: Which resource lane was used
    """

    job_id: str
    success: bool
    result: object | None = None
    error: str | None = None
    lane_label: str = ""


class IntraReactionScheduler:
    """Parallel QC job scheduler for single-reaction multi-task execution.

    Reads configuration from `intra_reaction_parallel` section.
    When disabled (default), jobs execute sequentially as fallback.
    When enabled, jobs run concurrently respecting resource lane constraints.

    Usage:
        scheduler = IntraReactionScheduler(config)

        # Submit jobs
        jobs = [
            ParallelQCJob(job_id="ts_opt", lane=ResourceLane(nproc=12, mem="36GB"), ...),
            ParallelQCJob(job_id="intermediate_opt", lane=ResourceLane(nproc=4, mem="12GB"), ...),
        ]

        # Execute all
        results = scheduler.run_parallel(jobs)
    """

    def __init__(self, config: ConfigDict):
        """Initialize scheduler from config.

        Args:
            config: Full RPH configuration dict
        """

        self.config: ConfigDict = config
        parallel_cfg = self._config_section("intra_reaction_parallel")
        self.enabled: bool = self._bool_value(parallel_cfg, "enabled", False)
        self.total_cores: int = self._int_value(parallel_cfg, "total_cores", 16)
        self.total_mem_str: str = self._str_value(parallel_cfg, "total_mem", "52GB")
        self.total_mem_gb: float = self._parse_mem_gb(self.total_mem_str)

        if self.enabled:
            logger.info(
                "IntraReactionScheduler: ENABLED (total_cores=%s, total_mem=%s)",
                self.total_cores,
                self.total_mem_str,
            )
        else:
            logger.debug("IntraReactionScheduler: DISABLED (sequential fallback)")

    @staticmethod
    def _as_config_dict(value: object) -> ConfigDict:
        if not isinstance(value, dict):
            return {}
        typed_value = cast(dict[object, object], value)
        result: ConfigDict = {}
        for key, item in typed_value.items():
            result[str(key)] = item
        return result

    def _config_section(self, *keys: str) -> ConfigDict:
        current: object = self.config
        for key in keys:
            current = self._as_config_dict(current).get(key)
        return self._as_config_dict(current)

    @staticmethod
    def _bool_value(section: ConfigDict, key: str, default: bool) -> bool:
        value = section.get(key)
        return value if isinstance(value, bool) else default

    @staticmethod
    def _int_value(section: ConfigDict, key: str, default: int) -> int:
        value = section.get(key)
        if isinstance(value, bool):
            return default
        if isinstance(value, int):
            return value
        if isinstance(value, float):
            return int(value)
        return default

    @staticmethod
    def _str_value(section: ConfigDict, key: str, default: str) -> str:
        value = section.get(key)
        return value if isinstance(value, str) else default

    @staticmethod
    def _parse_mem_gb(mem_str: object) -> float:
        """Parse memory string like '24GB' to float GB value."""

        if isinstance(mem_str, bool):
            return 48.0
        if isinstance(mem_str, (int, float)):
            return float(mem_str)

        s = str(mem_str).upper().strip()
        if s.endswith("GB"):
            try:
                return float(s[:-2])
            except ValueError:
                return 48.0
        if s.endswith("MB"):
            try:
                return float(s[:-2]) / 1024.0
            except ValueError:
                return 48.0
        return 48.0

    def validate_resource_budget(self, jobs: list[ParallelQCJob]) -> bool:
        """Check that all jobs fit within total resource budget.

        Args:
            jobs: List of parallel jobs to validate

        Returns:
            True if resources are within budget
        """

        total_nproc = sum(job.lane.nproc for job in jobs)
        total_mem_gb = sum(self._parse_mem_gb(job.lane.mem) for job in jobs)

        ok = True
        if total_nproc > self.total_cores:
            logger.warning(
                "Resource OVERCOMMIT: %s cores requested vs %s available",
                total_nproc,
                self.total_cores,
            )
            ok = False
        if total_mem_gb > self.total_mem_gb:
            logger.warning(
                "Resource OVERCOMMIT: %.1fGB mem requested vs %.1fGB available",
                total_mem_gb,
                self.total_mem_gb,
            )
            ok = False
        return ok

    def run_parallel(self, jobs: list[ParallelQCJob]) -> list[ParallelResult]:
        """Execute QC jobs in parallel (if enabled) or sequentially (fallback).

        Args:
            jobs: List of ParallelQCJob to execute

        Returns:
            List of ParallelResult, one per job, in job submission order
        """

        if not jobs:
            return []

        if not self.enabled:
            return self._run_sequential(jobs)

        resource_ok = self.validate_resource_budget(jobs)
        if not resource_ok:
            logger.debug("IntraReactionScheduler: proceeding despite resource overcommit")

        max_workers = len(jobs)
        logger.info(
            "IntraReactionScheduler: Running %s jobs in parallel (max_workers=%s)",
            len(jobs),
            max_workers,
        )

        results: dict[str, ParallelResult] = {}

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_job: dict[Future[ParallelResult], ParallelQCJob] = {}

            for job in jobs:
                future = executor.submit(self._execute_job, job)
                future_to_job[future] = job

            for future in as_completed(future_to_job):
                job = future_to_job[future]
                try:
                    result = future.result()
                    results[job.job_id] = result
                    status = "✅" if result.success else "❌"
                    logger.info(
                        "  %s Parallel job '%s' (%s): %s",
                        status,
                        job.job_id,
                        job.lane.label,
                        "OK" if result.success else result.error,
                    )
                except Exception as exc:
                    results[job.job_id] = ParallelResult(
                        job_id=job.job_id,
                        success=False,
                        error=str(exc),
                        lane_label=job.lane.label,
                    )
                    logger.error("  ❌ Parallel job '%s' exception: %s", job.job_id, exc)

        return [
            results.get(
                job.job_id,
                ParallelResult(
                    job_id=job.job_id,
                    success=False,
                    error="job_not_found",
                ),
            )
            for job in jobs
        ]

    def _run_sequential(self, jobs: list[ParallelQCJob]) -> list[ParallelResult]:
        """Fallback: execute jobs sequentially."""

        logger.info(
            "IntraReactionScheduler: Running %s jobs sequentially (parallel disabled)",
            len(jobs),
        )
        results: list[ParallelResult] = []
        for job in jobs:
            result = self._execute_job(job)
            results.append(result)
            status = "✅" if result.success else "❌"
            logger.info(
                "  %s Sequential job '%s': %s",
                status,
                job.job_id,
                "OK" if result.success else result.error,
            )
        return results

    def _execute_job(self, job: ParallelQCJob) -> ParallelResult:
        """Execute a single QC job."""

        logger.debug(
            "Executing job '%s' on %s (%s cores, %s)",
            job.job_id,
            job.lane.label,
            job.lane.nproc,
            job.lane.mem,
        )
        try:
            result = job.func(*job.args, **job.kwargs)
            return ParallelResult(
                job_id=job.job_id,
                success=True,
                result=result,
                lane_label=job.lane.label,
            )
        except Exception as exc:
            logger.error("Job '%s' failed: %s", job.job_id, exc)
            return ParallelResult(
                job_id=job.job_id,
                success=False,
                error=str(exc),
                lane_label=job.lane.label,
            )

    def get_s1_conformer_lane(self) -> ResourceLane:
        """Get resource lane for S1 conformer parallel optimization."""

        s1_cfg = self._config_section("intra_reaction_parallel", "s1")
        return ResourceLane(
            nproc=self._int_value(s1_cfg, "opt_cores_per_job", 8),
            mem=self._str_value(s1_cfg, "opt_mem_per_job", "24GB"),
            label="s1_conformer",
        )

    def get_s1_molecule_lane(self) -> ResourceLane:
        """Get resource lane for one molecule-level S1 anchor job."""

        s1_cfg = self._config_section("intra_reaction_parallel", "s1")
        return ResourceLane(
            nproc=self._int_value(s1_cfg, "crest_cores", 8),
            mem=self._str_value(s1_cfg, "molecule_mem_per_job", "24GB"),
            label="s1_molecule",
        )

    def get_s3_ts_lane(self) -> ResourceLane:
        """Get resource lane for S3 TS optimization."""

        s3_cfg = self._config_section("intra_reaction_parallel", "s3")
        return ResourceLane(
            nproc=self._int_value(s3_cfg, "ts_opt_cores", 12),
            mem=self._str_value(s3_cfg, "ts_opt_mem", "36GB"),
            label="s3_ts_opt",
        )

    def get_s3_intermediate_lane(self) -> ResourceLane:
        """Get resource lane for S3 intermediate optimization."""

        s3_cfg = self._config_section("intra_reaction_parallel", "s3")
        return ResourceLane(
            nproc=self._int_value(s3_cfg, "intermediate_opt_cores", 4),
            mem=self._str_value(s3_cfg, "intermediate_opt_mem", "12GB"),
            label="s3_intermediate_opt",
        )

    def get_s3_sp_lane(self) -> ResourceLane:
        """Get resource lane for S3 parallel single-point."""

        s3_cfg = self._config_section("intra_reaction_parallel", "s3")
        sp_maxcore_per_job = self._int_value(s3_cfg, "sp_maxcore_per_job", 2600)
        sp_cores_per_job = self._int_value(s3_cfg, "sp_cores_per_job", 8)
        return ResourceLane(
            nproc=sp_cores_per_job,
            mem=f"{sp_maxcore_per_job * sp_cores_per_job // 1024}GB",
            label="s3_sp",
        )

    @property
    def s1_conformer_parallel(self) -> bool:
        """Whether S1 conformer parallel optimization is enabled."""

        s1_cfg = self._config_section("intra_reaction_parallel", "s1")
        return self.enabled and self._bool_value(s1_cfg, "conformer_parallel", False)

    @property
    def s1_molecule_parallel(self) -> bool:
        """Whether S1 molecule-level anchor parallelization is enabled."""

        s1_cfg = self._config_section("intra_reaction_parallel", "s1")
        return self.enabled and self._bool_value(s1_cfg, "molecule_parallel", False)

    @property
    def s1_max_molecule_workers(self) -> int:
        """Max number of concurrent molecule-level S1 anchor workers."""

        s1_cfg = self._config_section("intra_reaction_parallel", "s1")
        return self._int_value(s1_cfg, "max_molecule_workers", 1)

    @property
    def s1_max_conformer_workers(self) -> int:
        """Max number of concurrent conformer workers."""

        s1_cfg = self._config_section("intra_reaction_parallel", "s1")
        return self._int_value(s1_cfg, "max_conformer_workers", 2)

    @property
    def s3_parallel_intermediate_ts(self) -> bool:
        """Whether S3 intermediate+TS parallel is enabled."""

        s3_cfg = self._config_section("intra_reaction_parallel", "s3")
        return self.enabled and self._bool_value(s3_cfg, "parallel_intermediate_ts", False)

    @property
    def s3_sp_parallel(self) -> bool:
        """Whether S3 parallel SP is enabled."""

        s3_cfg = self._config_section("intra_reaction_parallel", "s3")
        return self.enabled and self._bool_value(s3_cfg, "sp_parallel", False)
