from pathlib import Path
from typing import Any, Dict

from rph_core.utils.data_types import ScanResult
from rph_core.utils.qc_interface import XTBInterface
from rph_core.utils.xtb_runner import XTBRunner


def _write_min_xyz(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "2",
                "min",
                "H 0.000 0.000 0.000",
                "H 0.750 0.000 0.000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_xtb_runner_scan_command_includes_gfn_alpb_and_etemp(tmp_path: Path, monkeypatch) -> None:
    input_xyz = tmp_path / "input.xyz"
    _write_min_xyz(input_xyz)

    runner = XTBRunner(config={"executables": {"xtb": {"path": "python3"}}, "resources": {"nproc": 2}}, work_dir=tmp_path)
    captured = {"cmd": []}

    def _fake_run_command(cmd, log_file=None):
        captured["cmd"] = cmd
        if log_file is not None:
            log_file.write_text("", encoding="utf-8")
        return None

    def _fake_parse_scan_log(_):
        return ScanResult(success=True, energies=[-1.0, -0.9], max_energy_index=1, ts_guess_xyz=tmp_path / "guess.xyz")

    monkeypatch.setattr(runner, "_run_command", _fake_run_command)
    monkeypatch.setattr(runner, "_parse_scan_log", _fake_parse_scan_log)

    runner.run_scan(
        input_xyz=input_xyz,
        constraints={"0 1": 2.0},
        scan_range=(3.5, 1.8),
        scan_steps=8,
        solvent="acetone",
        gfn_level=1,
        etemp=450.0,
    )

    cmd = captured["cmd"]
    assert "--gfn" in cmd and cmd[cmd.index("--gfn") + 1] == "1"
    assert "--alpb" in cmd and cmd[cmd.index("--alpb") + 1] == "acetone"
    assert "--etemp" in cmd and cmd[cmd.index("--etemp") + 1] == "450.0"


def test_xtb_interface_scan_reads_step2_xtb_settings(tmp_path: Path, monkeypatch) -> None:
    input_xyz = tmp_path / "input.xyz"
    _write_min_xyz(input_xyz)

    captured: Dict[str, Any] = {"kwargs": None}

    def _fake_run_scan(self, **kwargs):
        captured["kwargs"] = kwargs
        return ScanResult(success=True, energies=[-1.0, -0.9], max_energy_index=1)

    monkeypatch.setattr(XTBRunner, "run_scan", _fake_run_scan)

    interface = XTBInterface(
        gfn_level=2,
        solvent="water",
        nproc=1,
        config={
            "executables": {"xtb": {"path": "python3"}},
            "step2": {
                "xtb_settings": {
                    "gfn_level": 0,
                    "solvent": "acetone",
                    "etemp": 300,
                }
            },
        },
    )

    interface.scan(
        xyz_file=input_xyz,
        output_dir=tmp_path / "scan_out",
        constraints={"0 1": 2.0},
        scan_range=(3.5, 1.8),
        scan_steps=10,
    )

    assert captured["kwargs"] is not None
    assert captured["kwargs"]["gfn_level"] == 0
    assert captured["kwargs"]["solvent"] == "acetone"
    assert captured["kwargs"]["etemp"] == 300.0


def test_xtb_interface_scan_sandbox_remaps_frame_paths_with_subdirs(tmp_path: Path, monkeypatch) -> None:
    input_xyz = tmp_path / "input.xyz"
    _write_min_xyz(input_xyz)

    def _fake_run_scan(self, **kwargs):
        frame_dir = self.work_dir / "scan_frames"
        frame_dir.mkdir(parents=True, exist_ok=True)
        frame = frame_dir / "frame_000.xyz"
        frame.write_text("2\nframe\nH 0.0 0.0 0.0\nH 0.7 0.0 0.0\n", encoding="utf-8")
        scan_log = self.work_dir / "xtb_scan.log"
        scan_log.write_text("", encoding="utf-8")
        ts_guess = self.work_dir / "xtb_scan_ts_guess.xyz"
        ts_guess.write_text("2\nts\nH 0.0 0.0 0.0\nH 0.8 0.0 0.0\n", encoding="utf-8")
        return ScanResult(
            success=True,
            energies=[-1.0, -0.9],
            geometries=[frame],
            max_energy_index=1,
            ts_guess_xyz=ts_guess,
            scan_log=scan_log,
        )

    monkeypatch.setattr(XTBRunner, "run_scan", _fake_run_scan)
    monkeypatch.setattr("rph_core.utils.qc_interface.is_path_toxic", lambda _path: True)

    interface = XTBInterface(
        config={
            "executables": {"xtb": {"path": "python3"}},
            "step2": {"xtb_settings": {"gfn_level": 2, "solvent": "acetone"}},
        }
    )
    result = interface.scan(
        xyz_file=input_xyz,
        output_dir=tmp_path / "scan_out",
        constraints={"0 1": 2.0},
        scan_range=(3.5, 1.8),
        scan_steps=10,
    )

    assert result.success is True
    assert isinstance(result.geometries, list)
    assert isinstance(result.geometries[0], Path)
    assert result.geometries[0].exists()
    assert "scan_frames" in str(result.geometries[0])
