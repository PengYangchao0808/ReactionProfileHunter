import threading
import time

import pytest

from rph_core.utils.stage_scheduler import StageSchedule, resolve_stage_schedule, run_structure_queue


def _config(
    *,
    max_workers=2,
    nproc=16,
    nproc_per_job=8,
    memory_per_job="16GB",
    mem="32GB",
    max_concurrent_by_role=None,
    priority_roles=None,
):
    return {
        "resources": {"nproc": nproc, "mem": mem},
        "step3": {
            "scheduling": {
                "enabled": True,
                "max_workers": max_workers,
                "nproc_per_job": nproc_per_job,
                "memory_per_job": memory_per_job,
                "priority_roles": priority_roles or ["ts", "product"],
                "max_concurrent_by_role": max_concurrent_by_role or {},
            }
        },
    }


def test_empty_role_limits_keep_parallel_behavior():
    schedule = resolve_stage_schedule(
        _config(max_concurrent_by_role={}, priority_roles=["product"]),
        "step3",
    )
    first_pair = threading.Barrier(2)
    lock = threading.Lock()
    active = 0
    max_active = 0
    structures = [{"id": str(index), "role": "product"} for index in range(4)]

    def worker(structure):
        nonlocal active, max_active
        with lock:
            active += 1
            max_active = max(max_active, active)
        try:
            if structure["id"] in {"0", "1"}:
                first_pair.wait(timeout=1.0)
            return structure["id"]
        finally:
            with lock:
                active -= 1

    results = run_structure_queue(
        structures,
        worker,
        schedule,
        thread_name_prefix="test-role-free",
    )

    assert dict(schedule.max_concurrent_by_role) == {}
    assert max_active == 2
    assert results == ["0", "1", "2", "3"]


def test_ts_role_limit_one_runs_serially():
    schedule = resolve_stage_schedule(
        _config(max_concurrent_by_role={"ts": 1}, priority_roles=["ts"]),
        "step3",
    )
    lock = threading.Lock()
    active_ts = 0
    max_active_ts = 0
    structures = [{"id": f"ts-{index}", "role": "ts"} for index in range(4)]

    def worker(structure):
        nonlocal active_ts, max_active_ts
        with lock:
            active_ts += 1
            max_active_ts = max(max_active_ts, active_ts)
        try:
            time.sleep(0.03)
            return structure["id"]
        finally:
            with lock:
                active_ts -= 1

    results = run_structure_queue(
        structures,
        worker,
        schedule,
        thread_name_prefix="test-ts-serial",
    )

    assert max_active_ts == 1
    assert results == [f"ts-{index}" for index in range(4)]


def test_mixed_role_limits_remain_work_conserving():
    schedule = resolve_stage_schedule(
        _config(
            max_workers=4,
            nproc_per_job=4,
            memory_per_job="8GB",
            max_concurrent_by_role={"ts": 1, "product": 2},
            priority_roles=["ts", "product"],
        ),
        "step3",
    )
    ts_started = threading.Event()
    release_all = threading.Event()
    concurrency_target = threading.Event()
    lock = threading.Lock()
    active_by_role = {"ts": 0, "product": 0}
    max_by_role = {"ts": 0, "product": 0}
    structures = [
        {"id": "ts-0", "role": "ts"},
        {"id": "ts-1", "role": "ts"},
        {"id": "ts-2", "role": "ts"},
        {"id": "product-0", "role": "product"},
        {"id": "product-1", "role": "product"},
    ]

    def worker(structure):
        role = structure["role"]
        with lock:
            active_by_role[role] += 1
            max_by_role[role] = max(max_by_role[role], active_by_role[role])
            if active_by_role["ts"] == 1 and active_by_role["product"] == 2:
                concurrency_target.set()
                release_all.set()
        try:
            if role == "ts":
                ts_started.set()
                assert release_all.wait(timeout=2.0)
            else:
                assert ts_started.wait(timeout=2.0)
                assert release_all.wait(timeout=2.0)
            return structure["id"]
        finally:
            with lock:
                active_by_role[role] -= 1

    run_structure_queue(
        structures,
        worker,
        schedule,
        thread_name_prefix="test-mixed-roles",
    )

    assert concurrency_target.is_set()
    assert max_by_role == {"ts": 1, "product": 2}


@pytest.mark.parametrize("limit", [0, -1])
def test_invalid_role_limit_must_be_positive(limit):
    with pytest.raises(ValueError, match="integer >= 1"):
        resolve_stage_schedule(_config(max_concurrent_by_role={"ts": limit}), "step3")


def test_role_limit_keys_must_be_lowercase():
    with pytest.raises(ValueError, match="lowercase"):
        resolve_stage_schedule(_config(max_concurrent_by_role={"TS": 1}), "step3")


def test_role_limit_above_max_workers_is_clamped(caplog):
    caplog.set_level("WARNING")

    schedule = resolve_stage_schedule(
        _config(max_workers=2, nproc_per_job=4, max_concurrent_by_role={"product": 3}),
        "step3",
    )

    assert dict(schedule.max_concurrent_by_role) == {"product": 2}
    assert "clamping" in caplog.text


def test_role_limit_respects_nproc_pool_validation():
    with pytest.raises(ValueError, match="resource pool"):
        resolve_stage_schedule(
            _config(max_workers=2, nproc=12, nproc_per_job=8, max_concurrent_by_role={"ts": 2}),
            "step3",
        )


def test_manifest_record_includes_role_limits():
    schedule = StageSchedule(
        enabled=True,
        max_workers=2,
        nproc_per_job=8,
        memory_per_job="16GB",
        priority_roles=("ts", "product"),
        max_concurrent_by_role={"ts": 1, "product": 2},
    )

    assert schedule.manifest_record()["max_concurrent_by_role"] == {"ts": 1, "product": 2}
