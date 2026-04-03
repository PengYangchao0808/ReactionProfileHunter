from pathlib import Path
import tempfile

from rph_core.utils.gau_xtb_interface import GauXTBInterface


def _write_min_xyz(path: Path) -> None:
    path.write_text("2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.8\n")


def test_optimize_ts_uses_sandbox_for_toxic_output_path(monkeypatch, tmp_path: Path) -> None:
    xyz = tmp_path / "ts_guess.xyz"
    _write_min_xyz(xyz)

    toxic_out = tmp_path / "gau xtb [toxic]"
    toxic_out.mkdir(parents=True, exist_ok=True)

    sandbox_dir = tmp_path / "gau_xtb_sandbox"
    sandbox_dir.mkdir(parents=True, exist_ok=True)

    observed = {"cwd": None}

    def _fake_mkdtemp(prefix: str, dir: str) -> str:
        _ = (prefix, dir)
        return str(sandbox_dir)

    def _fake_run(cmd, cwd, env, capture_output, text, timeout):
        _ = (cmd, env, capture_output, text, timeout)
        observed["cwd"] = cwd
        log_file = Path(cwd) / "input.log"
        log_file.write_text("Normal termination of Gaussian\n")

        class _Result:
            returncode = 0

        return _Result()

    def _fake_extract(self, log_file, engine_type):
        _ = (self, log_file, engine_type)
        return ([[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]], ["H", "H"], None)

    monkeypatch.setattr(tempfile, "mkdtemp", _fake_mkdtemp)
    monkeypatch.setattr("rph_core.utils.gau_xtb_interface.subprocess.run", _fake_run)
    monkeypatch.setattr(
        "rph_core.utils.gau_xtb_interface.LogParser.extract_last_converged_coords",
        _fake_extract,
    )

    config = {
        "executables": {
            "gaussian": {
                "wrapper_path": "./scripts/run_g16_worker.sh",
            }
        }
    }
    interface = GauXTBInterface(config=config, nproc=1)
    result = interface.optimize_ts(xyz_file=xyz, output_dir=toxic_out, task_name="TS-TEST-001")

    assert result.success is True
    assert observed["cwd"] == str(sandbox_dir)
    assert (toxic_out / "input.log").exists()
    assert (toxic_out / "ts_final.xyz").exists()


def test_optimize_ts_sets_xtb_env_from_config(monkeypatch, tmp_path: Path) -> None:
    xyz = tmp_path / "ts_guess.xyz"
    _write_min_xyz(xyz)

    out_dir = tmp_path / "gau_xtb_ok"
    out_dir.mkdir(parents=True, exist_ok=True)

    observed = {"env": None}

    def _fake_run(cmd, cwd, env, capture_output, text, timeout):
        _ = (cmd, cwd, capture_output, text, timeout)
        observed["env"] = env
        (Path(cwd) / "input.log").write_text("Normal termination of Gaussian\n")

        class _Result:
            returncode = 0

        return _Result()

    def _fake_extract(self, log_file, engine_type):
        _ = (self, log_file, engine_type)
        return ([[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]], ["H", "H"], None)

    monkeypatch.setattr("rph_core.utils.gau_xtb_interface.subprocess.run", _fake_run)
    monkeypatch.setattr(
        "rph_core.utils.gau_xtb_interface.LogParser.extract_last_converged_coords",
        _fake_extract,
    )

    xtb_path = "/opt/software/xtb/bin/xtb"
    config = {
        "executables": {
            "gaussian": {
                "wrapper_path": "./scripts/run_g16_worker.sh",
            },
            "xtb": {
                "path": xtb_path,
            },
        }
    }

    interface = GauXTBInterface(config=config, nproc=1)
    result = interface.optimize_ts(xyz_file=xyz, output_dir=out_dir, task_name="TS-TEST-002")

    assert result.success is True
    assert observed["env"] is not None
    assert observed["env"]["XTB_PATH"] == xtb_path
    assert observed["env"]["PATH"].startswith("/opt/software/xtb/bin:")


def test_optimize_ts_sandboxes_relative_output_under_toxic_cwd(monkeypatch, tmp_path: Path) -> None:
    xyz = tmp_path / "ts_guess.xyz"
    _write_min_xyz(xyz)

    toxic_cwd = tmp_path / "repo [toxic]"
    toxic_cwd.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(toxic_cwd)

    sandbox_dir = tmp_path / "gau_xtb_rel_sandbox"
    sandbox_dir.mkdir(parents=True, exist_ok=True)

    observed = {"cwd": None}

    def _fake_mkdtemp(prefix: str, dir: str) -> str:
        _ = (prefix, dir)
        return str(sandbox_dir)

    def _fake_run(cmd, cwd, env, capture_output, text, timeout):
        _ = (cmd, env, capture_output, text, timeout)
        observed["cwd"] = cwd
        (Path(cwd) / "input.log").write_text("Normal termination of Gaussian\n")

        class _Result:
            returncode = 0

        return _Result()

    def _fake_extract(self, log_file, engine_type):
        _ = (self, log_file, engine_type)
        return ([[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]], ["H", "H"], None)

    monkeypatch.setattr(tempfile, "mkdtemp", _fake_mkdtemp)
    monkeypatch.setattr("rph_core.utils.gau_xtb_interface.subprocess.run", _fake_run)
    monkeypatch.setattr(
        "rph_core.utils.gau_xtb_interface.LogParser.extract_last_converged_coords",
        _fake_extract,
    )

    config = {
        "executables": {
            "gaussian": {
                "wrapper_path": "./scripts/run_g16_worker.sh",
            }
        }
    }
    interface = GauXTBInterface(config=config, nproc=1)
    result = interface.optimize_ts(
        xyz_file=xyz,
        output_dir=Path("relative_output"),
        task_name="TS-TEST-003",
    )

    assert result.success is True
    assert observed["cwd"] == str(sandbox_dir)
    assert (toxic_cwd / "relative_output" / "input.log").exists()
