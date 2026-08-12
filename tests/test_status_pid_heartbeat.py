from __future__ import annotations

import json
import os
import threading
import time
from pathlib import Path
from typing import Any

import rph_core.steps.refinement.engine as refinement_engine_module
from rph_core.steps.step3_lowlevel import LowLevelEngine
from rph_core.utils import orca_interface as orca_module
from rph_core.utils.config_loader import load_config
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.utils.qc_models import QCJobResult
from rph_core.utils.s4_progress import S4ProgressReporter
from rph_core.utils.stage_progress import StageProgressReporter


def test_stage_progress_snapshot_includes_pid_and_heartbeat(tmp_path: Path) -> None:
    reporter = StageProgressReporter(
        tmp_path / "S3_LowLevel",
        "S3",
        config={"ui": {"heartbeat_seconds": 17}},
    )

    status = json.loads((tmp_path / "S3_LowLevel" / "status.json").read_text(encoding="utf-8"))
    assert status["pid"] == os.getpid()
    assert status["heartbeat_interval_seconds"] == 17
    assert time.time() - float(status["heartbeat_at"]) < 2.0
    assert reporter.heartbeat_interval_seconds == 17


def test_stage_progress_reporter_tracks_queued_structure_subprocess_and_heartbeat(tmp_path: Path) -> None:
    reporter = StageProgressReporter(tmp_path / "S3_LowLevel", "S3")
    sandbox_dir = tmp_path / "S3_LowLevel" / "product_ts" / "diagnostics" / "attempt_001_opt_ts"
    sandbox_dir.mkdir(parents=True)

    reporter.register_structure("product_ts", role="ts", kind="ts")
    reporter.start_structure("product_ts", role="ts", kind="ts")
    reporter.report_structure_subprocess("product_ts", pid=4321, sandbox_path=str(sandbox_dir))
    initial = json.loads((tmp_path / "S3_LowLevel" / "status.json").read_text(encoding="utf-8"))
    initial_heartbeat = float(initial["structures"]["product_ts"].get("heartbeat_at") or 0)
    time.sleep(0.01)
    reporter.touch_structure_heartbeat("product_ts")

    status = json.loads((tmp_path / "S3_LowLevel" / "status.json").read_text(encoding="utf-8"))
    structure = status["structures"]["product_ts"]
    assert status["total_structures"] == 1
    assert structure["status"] == "running"
    assert structure["pid"] == 4321
    assert structure["sandbox_path"] == str(sandbox_dir)
    assert float(structure["heartbeat_at"]) > initial_heartbeat


def test_s4_progress_reporter_report_structure_subprocess_sets_pid_and_sandbox_path(
    tmp_path: Path,
) -> None:
    reporter = S4ProgressReporter(
        tmp_path / "S4_HighLevel",
        [{"id": "product_ts", "kind": "ts", "fallback_xyz": "seed.xyz"}],
        config={"ui": {"heartbeat_seconds": 9}},
    )
    sandbox_dir = tmp_path / "S4_HighLevel" / "product_ts" / "opt"
    sandbox_dir.mkdir(parents=True, exist_ok=True)

    reporter.report_structure_subprocess(
        "product_ts",
        pid=4321,
        sandbox_path=str(sandbox_dir),
    )

    status = json.loads((tmp_path / "S4_HighLevel" / "status.json").read_text(encoding="utf-8"))
    structure = status["structures"][0]
    assert status["pid"] == os.getpid()
    assert status["heartbeat_interval_seconds"] == 9
    assert structure["pid"] == 4321
    assert structure["sandbox_path"] == str(sandbox_dir)


def test_s4_progress_reporter_touch_structure_heartbeat_updates_fields(tmp_path: Path) -> None:
    reporter = S4ProgressReporter(
        tmp_path / "S4_HighLevel",
        [{"id": "product_ts", "kind": "ts", "fallback_xyz": "seed.xyz"}],
    )
    sandbox_dir = tmp_path / "S4_HighLevel" / "product_ts" / "opt"
    sandbox_dir.mkdir(parents=True, exist_ok=True)
    out_file = sandbox_dir / "orca.out"
    out_file.write_text("first\n", encoding="utf-8")
    reporter.report_structure_subprocess(
        "product_ts",
        pid=4321,
        sandbox_path=str(sandbox_dir),
    )

    initial = json.loads((tmp_path / "S4_HighLevel" / "status.json").read_text(encoding="utf-8"))["structures"][0]
    initial_heartbeat = float(initial["heartbeat_at"])
    initial_output_update = float(initial["last_output_update"])

    time.sleep(0.05)
    out_file.write_text("second\n", encoding="utf-8")
    os.utime(out_file, None)
    reporter.touch_structure_heartbeat("product_ts")

    updated = json.loads((tmp_path / "S4_HighLevel" / "status.json").read_text(encoding="utf-8"))["structures"][0]
    assert float(updated["heartbeat_at"]) > initial_heartbeat
    assert float(updated["last_output_update"]) >= initial_output_update
    assert float(updated["last_output_update"]) >= out_file.stat().st_mtime - 1e-6


def test_refinement_engine_passes_subprocess_callback_and_callback_fires(
    monkeypatch: Any,
    tmp_path: Path,
) -> None:
    input_xyz = tmp_path / "seed.xyz"
    input_xyz.write_text("1\nseed\nH 0.0 0.0 0.0\n", encoding="utf-8")
    callback_flags: list[bool] = []

    class FakeReporter:
        heartbeat_interval_seconds = 0.01

        def __init__(self) -> None:
            self.subprocess_calls: list[tuple[str, int | None, str | None]] = []
            self.heartbeat_calls: list[str] = []

        def report_structure_subprocess(
            self,
            structure_id: str,
            *,
            pid: int | None,
            sandbox_path: str | None,
        ) -> None:
            self.subprocess_calls.append((structure_id, pid, sandbox_path))

        def touch_structure_heartbeat(self, structure_id: str) -> None:
            self.heartbeat_calls.append(structure_id)

    class FakeProcess:
        def __init__(self, pid: int) -> None:
            self.pid = pid

        def poll(self) -> int:
            return 0

    def fake_opt(spec, input_path, output_dir, _config, subprocess_callback=None):
        del spec
        callback_flags.append(subprocess_callback is not None)
        assert subprocess_callback is not None
        subprocess_callback(FakeProcess(3001))
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_path).read_text(encoding="utf-8"), encoding="utf-8")
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_path),
            output_xyz=out_xyz,
            output_file=output_dir / "opt.out",
            energy_hartree=-1.0,
        )

    def fake_frequency(spec, input_path, output_dir, _config, subprocess_callback=None):
        del spec
        callback_flags.append(subprocess_callback is not None)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_path),
            output_file=output_dir / "freq.out",
            energy_hartree=-0.9,
            frequencies_cm1=(25.0, 125.0, 325.0),
        )

    def fake_sp(spec, input_path, output_dir, _config, subprocess_callback=None):
        del spec
        callback_flags.append(subprocess_callback is not None)
        assert subprocess_callback is not None
        subprocess_callback(FakeProcess(3002))
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_path),
            output_file=Path(output_dir) / "sp.out",
            energy_hartree=-1.1,
        )

    monkeypatch.setattr(refinement_engine_module, "run_optimization", fake_opt)
    monkeypatch.setattr(refinement_engine_module, "run_frequency", fake_frequency)
    monkeypatch.setattr(refinement_engine_module, "run_single_point", fake_sp)
    reporter = FakeReporter()
    engine = LowLevelEngine(load_config())
    engine.set_progress_reporter(reporter)

    result = json.loads(
        engine.run(
            [
                {
                    "id": "product",
                    "role": "product",
                    "kind": "minimum",
                    "input_xyz": str(input_xyz),
                    "s1_thermochemistry_status": "complete",
                }
            ],
            tmp_path / "stage",
        ).read_text(encoding="utf-8")
    )["structures"][0]

    assert callback_flags == [True, True, True]
    assert reporter.subprocess_calls == [
        ("product", 3001, str(tmp_path / "stage" / "product" / "diagnostics" / "attempt_001_opt")),
        ("product", 3002, str(tmp_path / "stage" / "product" / "diagnostics" / "attempt_003_sp")),
    ]
    assert reporter.heartbeat_calls
    assert result["status"] == "complete"


def test_refinement_engine_heartbeat_thread_starts_and_stops_cleanly(
    monkeypatch: Any,
    tmp_path: Path,
) -> None:
    input_xyz = tmp_path / "seed.xyz"
    input_xyz.write_text("1\nseed\nH 0.0 0.0 0.0\n", encoding="utf-8")

    class FakeReporter:
        heartbeat_interval_seconds = 0.01

        def __init__(self) -> None:
            self.heartbeat_seen = threading.Event()
            self.subprocess_calls: list[tuple[str, int | None, str | None]] = []

        def report_structure_subprocess(
            self,
            structure_id: str,
            *,
            pid: int | None,
            sandbox_path: str | None,
        ) -> None:
            self.subprocess_calls.append((structure_id, pid, sandbox_path))

        def touch_structure_heartbeat(self, structure_id: str) -> None:
            self.heartbeat_seen.set()

    class FakeProcess:
        def __init__(self, pid: int) -> None:
            self.pid = pid
            self.done = False

        def poll(self) -> int | None:
            return 0 if self.done else None

    reporter = FakeReporter()

    def fake_opt(spec, input_path, output_dir, _config, subprocess_callback=None):
        del spec
        assert subprocess_callback is not None
        process = FakeProcess(4001)
        subprocess_callback(process)
        assert reporter.heartbeat_seen.wait(timeout=1.0)
        process.done = True
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_path).read_text(encoding="utf-8"), encoding="utf-8")
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_path),
            output_xyz=out_xyz,
            output_file=output_dir / "opt.out",
            energy_hartree=-1.0,
        )

    def fake_frequency(spec, input_path, output_dir, _config, subprocess_callback=None):
        del spec, subprocess_callback
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_path),
            output_file=Path(output_dir) / "freq.out",
            energy_hartree=-0.9,
            frequencies_cm1=(25.0, 125.0, 325.0),
        )

    def fake_sp(spec, input_path, output_dir, _config, subprocess_callback=None):
        del spec, subprocess_callback
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_path),
            output_file=Path(output_dir) / "sp.out",
            energy_hartree=-1.1,
        )

    monkeypatch.setattr(refinement_engine_module, "run_optimization", fake_opt)
    monkeypatch.setattr(refinement_engine_module, "run_frequency", fake_frequency)
    monkeypatch.setattr(refinement_engine_module, "run_single_point", fake_sp)
    engine = LowLevelEngine(load_config())
    engine.set_progress_reporter(reporter)
    engine.run(
        [
            {
                "id": "product",
                "role": "product",
                "kind": "minimum",
                "input_xyz": str(input_xyz),
                "s1_thermochemistry_status": "complete",
            }
        ],
        tmp_path / "stage",
    )

    assert reporter.subprocess_calls == [
        ("product", 4001, str(tmp_path / "stage" / "product" / "diagnostics" / "attempt_001_opt"))
    ]
    assert not [
        thread.name
        for thread in threading.enumerate()
        if thread.name.startswith("rph-structure-heartbeat-")
    ]


def test_orca_interface_run_orca_accepts_subprocess_callback_without_breaking_callers(
    monkeypatch: Any,
    tmp_path: Path,
) -> None:
    class FakePopen:
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            self.pid = 2468
            self.returncode = 0

        def communicate(self, timeout: int | None = None) -> tuple[str, str]:
            return "", ""

        def kill(self) -> None:
            self.returncode = -9

    monkeypatch.setattr(orca_module.subprocess, "Popen", FakePopen)

    interface = ORCAInterface(
        method="B97-3c",
        basis="",
        aux_basis="",
        nprocs=1,
        solvent="",
        orca_binary_path="/bin/true",
        config={},
    )
    inp_file = tmp_path / "job.inp"
    inp_file.write_text("! B97-3c\n", encoding="utf-8")
    callback_pids: list[int] = []

    out_file = interface._run_orca(
        inp_file,
        tmp_path,
        timeout=1,
        subprocess_callback=lambda process: callback_pids.append(process.pid),
    )
    assert out_file == inp_file.with_suffix(".out")
    assert callback_pids == [2468]
    assert interface._last_subprocess is None

    second_out_file = interface._run_orca(inp_file, tmp_path, timeout=1)
    assert second_out_file.exists()
    assert interface._last_subprocess is None
