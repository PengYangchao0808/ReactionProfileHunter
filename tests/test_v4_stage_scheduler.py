import json
import threading
from pathlib import Path

import rph_core.steps.refinement.engine as refinement_engine_module
from rph_core.steps.step3_lowlevel import LowLevelEngine
from rph_core.utils.config_loader import load_config
from rph_core.utils.qc_models import QCJobResult
from rph_core.utils.stage_scheduler import (
    resolve_stage_schedule,
    run_structure_queue,
    worker_config,
    worker_theory,
)


def _config():
    return {
        "resources": {"nproc": 16, "mem": "32GB"},
        "step3": {
            "scheduling": {
                "enabled": True,
                "max_workers": 2,
                "nproc_per_job": 8,
                "memory_per_job": "16GB",
                "priority_roles": ["ts", "intermediate", "product"],
            }
        },
    }


def test_stage_schedule_splits_cpu_and_memory_for_every_qc_job():
    schedule = resolve_stage_schedule(_config(), "step3")
    derived_config = worker_config(_config(), schedule)
    derived_theory = worker_theory(
        {"optimization": {"method": "B97-3c"}, "single_point": {"method": "r2SCAN-3c"}},
        schedule,
    )

    assert schedule.max_workers == 2
    assert derived_config["resources"]["nproc"] == 8
    assert derived_config["resources"]["mem"] == "16GB"
    assert derived_theory["optimization"]["nproc"] == 8
    assert derived_theory["optimization"]["mem"] == "16GB"
    assert derived_theory["single_point"]["nproc"] == 8
    assert derived_theory["single_point"]["mem"] == "16GB"


def test_dynamic_queue_starts_next_job_without_waiting_for_long_partner():
    schedule = resolve_stage_schedule(_config(), "step3")
    ts_started = threading.Event()
    release_ts = threading.Event()
    product_started = threading.Event()
    structures = [
        {"id": "product", "role": "product"},
        {"id": "ts", "role": "ts"},
        {"id": "intermediate", "role": "intermediate"},
    ]

    def worker(structure):
        if structure["role"] == "ts":
            ts_started.set()
            assert release_ts.wait(timeout=2.0)
        elif structure["role"] == "intermediate":
            assert ts_started.wait(timeout=2.0)
        else:
            product_started.set()
            release_ts.set()
        return structure["id"]

    results = run_structure_queue(
        structures,
        worker,
        schedule,
        thread_name_prefix="test-s3",
    )

    assert product_started.is_set()
    assert results == ["product", "ts", "intermediate"]


def test_s3_engine_applies_worker_budget_and_records_schedule(monkeypatch, tmp_path):
    xyz = tmp_path / "seed.xyz"
    xyz.write_text("1\nseed\nH 0.0 0.0 0.0\n", encoding="utf-8")
    calls = []

    def fake_opt(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del subprocess_callback
        calls.append(("opt", spec.nproc, spec.memory, config["resources"]["nproc"], config["resources"]["mem"]))
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_xyz=out_xyz,
            output_file=output_dir / "opt.out",
            energy_hartree=-1.0,
        )

    def fake_frequency(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / "freq.out"
        out_file.write_text("mock freq\n", encoding="utf-8")
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=out_file,
            energy_hartree=-0.9,
            frequencies_cm1=(-321.4, 120.0, 311.2) if "ts" in str(output_dir) else (25.0, 125.0, 325.0),
        )

    def fake_sp(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del subprocess_callback
        calls.append(("sp", spec.nproc, spec.memory, config["resources"]["nproc"], config["resources"]["mem"]))
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=output_dir / "sp.out",
            energy_hartree=-1.1,
        )

    monkeypatch.setattr(refinement_engine_module, "run_optimization", fake_opt)
    monkeypatch.setattr(refinement_engine_module, "run_frequency", fake_frequency)
    monkeypatch.setattr(refinement_engine_module, "run_single_point", fake_sp)
    config = load_config()
    config["resources"] = _config()["resources"]
    config["step3"] = _config()["step3"]
    structures = [
        {"id": "product", "kind": "minimum", "role": "product", "input_xyz": str(xyz)},
        {"id": "ts", "kind": "ts", "role": "ts", "input_xyz": str(xyz)},
    ]

    manifest_path = LowLevelEngine(config).run(structures, tmp_path / "S3_LowLevel")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert manifest["schema_version"] == "refinement_manifest_v1"
    assert manifest["stage"] == "S3"
    assert manifest["scheduling"]["max_workers"] == 2
    assert [row["id"] for row in manifest["structures"]] == ["product", "ts"]
    assert any(kind == "opt" for kind, *_rest in calls)
    assert any(kind == "sp" for kind, *_rest in calls)
