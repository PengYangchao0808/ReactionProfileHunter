from pathlib import Path
import tempfile

from rph_core.utils.qc_interface import GaussianInterface, LogParser


def _orientation_block(header: str, xyz: float) -> str:
    return (
        f" {header}:\n"
        " ---------------------------------------------------------------------\n"
        " Center     Atomic      Atomic             Coordinates (Angstroms)\n"
        " Number     Number       Type             X           Y           Z\n"
        " ---------------------------------------------------------------------\n"
        f"      1          6           0        {xyz: .6f}    0.000000    0.000000\n"
        "      2          1           0        0.000000    0.000000    1.000000\n"
        " ---------------------------------------------------------------------\n"
    )


def test_extract_final_geometry_parses_fixture_orientation_block() -> None:
    fixture = Path(
        "tests/tmp_v2_2_test/da_reaction/S1_Anchor/ethylene/dft/conf_000.log"
    )
    content = fixture.read_text(errors="ignore")

    atoms = LogParser.extract_final_geometry(content)

    assert len(atoms) == 6
    assert atoms[0]["symbol"] == "C"


def test_extract_final_geometry_uses_last_orientation_block() -> None:
    log_content = (
        _orientation_block("Input orientation", 1.111111)
        + " Some unrelated line\n"
        + _orientation_block("Standard orientation", 2.222222)
        + " Normal termination of Gaussian 16.\n"
    )

    atoms = LogParser.extract_final_geometry(log_content)

    assert len(atoms) == 2
    assert abs(atoms[0]["x"] - 2.222222) < 1e-6


def test_extract_final_geometry_supports_input_orientation_only() -> None:
    log_content = _orientation_block("Input orientation", 3.333333)

    atoms = LogParser.extract_final_geometry(log_content)

    assert len(atoms) == 2
    assert abs(atoms[0]["x"] - 3.333333) < 1e-6


def test_constrained_optimize_reports_geometry_parse_failure(monkeypatch, tmp_path: Path) -> None:
    xyz = tmp_path / "m.xyz"
    xyz.write_text("1\ncomment\nH 0.0 0.0 0.0\n")

    out = tmp_path / "gauss"
    out.mkdir(parents=True, exist_ok=True)

    def _fake_run(cmd, cwd, capture_output, text, timeout):
        _ = (cmd, capture_output, text, timeout)
        log_file = Path(cwd) / f"{xyz.stem}_constrained.log"
        log_file.write_text(
            " Normal termination of Gaussian 16 at Sat Jan 31 21:00:33 2026.\n"
        )

        class _Result:
            returncode = 0
            stderr = ""

        return _Result()

    monkeypatch.setattr("rph_core.utils.qc_interface.subprocess.run", _fake_run)

    gi = GaussianInterface(config={"executables": {"gaussian": {"use_wrapper": False}}})
    result = gi.constrained_optimize(
        xyz_file=xyz,
        output_dir=out,
        frozen_indices=[0],
        charge=0,
        spin=1,
    )

    assert result.success is False
    assert result.converged is True
    assert result.coordinates is None
    assert result.error_message is not None
    assert "Gaussian geometry parse failed" in result.error_message


def test_constrained_optimize_reports_execution_failure_with_returncode(
    monkeypatch, tmp_path: Path
) -> None:
    xyz = tmp_path / "m.xyz"
    xyz.write_text("1\ncomment\nH 0.0 0.0 0.0\n")

    out = tmp_path / "gauss"
    out.mkdir(parents=True, exist_ok=True)

    def _fake_run(cmd, cwd, capture_output, text, timeout):
        _ = (cmd, capture_output, text, timeout)
        log_file = Path(cwd) / f"{xyz.stem}_constrained.log"
        log_file.write_text(
            " Error termination via Lnk1e in /opt/g16/l502.exe at Sat Jan 31 21:00:31 2026.\n"
        )

        class _Result:
            returncode = 1
            stderr = "SCF failed"

        return _Result()

    monkeypatch.setattr("rph_core.utils.qc_interface.subprocess.run", _fake_run)

    gi = GaussianInterface(config={"executables": {"gaussian": {"use_wrapper": False}}})
    result = gi.constrained_optimize(
        xyz_file=xyz,
        output_dir=out,
        frozen_indices=[0],
        charge=0,
        spin=1,
    )

    assert result.success is False
    assert result.converged is False
    assert result.error_message is not None
    assert "Gaussian execution failed (returncode=1)" in result.error_message
    assert "Error termination" in result.error_message


def test_constrained_optimize_returncode_zero_with_error_termination_is_nonconverged(
    monkeypatch, tmp_path: Path
) -> None:
    xyz = tmp_path / "m.xyz"
    xyz.write_text("1\ncomment\nH 0.0 0.0 0.0\n")

    out = tmp_path / "gauss"
    out.mkdir(parents=True, exist_ok=True)

    def _fake_run(cmd, cwd, capture_output, text, timeout):
        _ = (cmd, capture_output, text, timeout)
        log_file = Path(cwd) / f"{xyz.stem}_constrained.log"
        log_file.write_text(
            " Error termination via Lnk1e in /opt/g16/l123.exe at Sat Jan 31 21:00:31 2026.\n"
        )

        class _Result:
            returncode = 0
            stderr = ""

        return _Result()

    monkeypatch.setattr("rph_core.utils.qc_interface.subprocess.run", _fake_run)

    gi = GaussianInterface(config={"executables": {"gaussian": {"use_wrapper": False}}})
    result = gi.constrained_optimize(
        xyz_file=xyz,
        output_dir=out,
        frozen_indices=[0],
        charge=0,
        spin=1,
    )

    assert result.success is False
    assert result.converged is False
    assert result.error_message is not None
    assert "Gaussian did not converge" in result.error_message
    assert "Error termination" in result.error_message


def test_constrained_optimize_rejects_out_of_range_constraints_without_subprocess(
    monkeypatch, tmp_path: Path
) -> None:
    xyz = tmp_path / "m.xyz"
    xyz.write_text("2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

    called = {"run": False}

    def _fake_run(*args, **kwargs):
        _ = (args, kwargs)
        called["run"] = True
        raise AssertionError("subprocess.run should not be called for invalid constraints")

    monkeypatch.setattr("rph_core.utils.qc_interface.subprocess.run", _fake_run)

    gi = GaussianInterface(config={"executables": {"gaussian": {"use_wrapper": False}}})
    result = gi.constrained_optimize(
        xyz_file=xyz,
        output_dir=tmp_path / "gauss",
        distance_constraints=[(0, 5)],
        charge=0,
        spin=1,
    )

    assert called["run"] is False
    assert result.success is False
    assert result.converged is False
    assert result.error_message is not None
    assert "out of range" in result.error_message


def test_constrained_optimize_uses_sandbox_for_toxic_output_path(
    monkeypatch, tmp_path: Path
) -> None:
    xyz = tmp_path / "m.xyz"
    xyz.write_text("1\ncomment\nH 0.0 0.0 0.0\n")

    sandbox_dir = tmp_path / "sandbox_run"
    sandbox_dir.mkdir(parents=True, exist_ok=True)

    toxic_out = tmp_path / "gauss [toxic]"
    toxic_out.mkdir(parents=True, exist_ok=True)

    def _fake_mkdtemp(prefix: str, dir: str) -> str:
        _ = (prefix, dir)
        return str(sandbox_dir)

    def _fake_run(cmd, cwd, capture_output, text, timeout):
        _ = (cmd, capture_output, text, timeout)
        assert Path(cwd) == sandbox_dir
        log_file = sandbox_dir / f"{xyz.stem}_constrained.log"
        log_file.write_text(
            " Error termination via Lnk1e in /opt/g16/l502.exe at Sat Jan 31 21:00:31 2026.\n"
        )

        class _Result:
            returncode = 1
            stderr = "SCF failed"

        return _Result()

    monkeypatch.setattr(tempfile, "mkdtemp", _fake_mkdtemp)
    monkeypatch.setattr("rph_core.utils.qc_interface.subprocess.run", _fake_run)

    gi = GaussianInterface(config={"executables": {"gaussian": {"use_wrapper": False}}})
    result = gi.constrained_optimize(
        xyz_file=xyz,
        output_dir=toxic_out,
        frozen_indices=[0],
        charge=0,
        spin=1,
    )

    copied_log = toxic_out / f"{xyz.stem}_constrained.log"
    assert copied_log.exists()
    assert result.success is False
    assert result.converged is False
    assert result.output_file == copied_log
    assert result.error_message is not None
    assert "returncode=1" in result.error_message


def test_constrained_optimize_default_route_normalizes_def2_basis(
    monkeypatch, tmp_path: Path
) -> None:
    xyz = tmp_path / "m.xyz"
    xyz.write_text("2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.8\n")

    out = tmp_path / "gauss"
    out.mkdir(parents=True, exist_ok=True)

    def _fake_run(cmd, cwd, capture_output, text, timeout):
        _ = (cmd, capture_output, text, timeout)
        log_file = Path(cwd) / f"{xyz.stem}_constrained.log"
        log_file.write_text("QPErr --- A syntax error was detected in the input line.\n")

        class _Result:
            returncode = 1
            stderr = ""

        return _Result()

    monkeypatch.setattr("rph_core.utils.qc_interface.subprocess.run", _fake_run)

    config = {
        "executables": {"gaussian": {"use_wrapper": False}},
        "theory": {
            "optimization": {
                "method": "B3LYP",
                "basis": "def2-SVP",
                "dispersion": "GD3BJ",
                "solvent": "acetone",
            }
        },
    }
    gi = GaussianInterface(config=config)
    _ = gi.constrained_optimize(
        xyz_file=xyz,
        output_dir=out,
        distance_constraints=[(0, 1)],
        charge=0,
        spin=1,
    )

    gjf_content = (out / f"{xyz.stem}_constrained.gjf").read_text()
    route_line = next((line.strip() for line in gjf_content.splitlines() if line.strip().startswith("#p")), "")

    assert "def2-SVP" not in route_line
    assert "def2SVP" in route_line
    assert "Opt=(CalcFC,ModRedundant)" in route_line


def test_constrained_optimize_injects_modredundant_into_opt_clause(
    monkeypatch, tmp_path: Path
) -> None:
    xyz = tmp_path / "m.xyz"
    xyz.write_text("2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.8\n")

    out = tmp_path / "gauss"
    out.mkdir(parents=True, exist_ok=True)

    def _fake_run(cmd, cwd, capture_output, text, timeout):
        _ = (cmd, capture_output, text, timeout)
        log_file = Path(cwd) / f"{xyz.stem}_constrained.log"
        log_file.write_text("QPErr --- A syntax error was detected in the input line.\n")

        class _Result:
            returncode = 1
            stderr = ""

        return _Result()

    monkeypatch.setattr("rph_core.utils.qc_interface.subprocess.run", _fake_run)

    gi = GaussianInterface(config={"executables": {"gaussian": {"use_wrapper": False}}})
    _ = gi.constrained_optimize(
        xyz_file=xyz,
        output_dir=out,
        distance_constraints=[(0, 1)],
        charge=0,
        spin=1,
        route="B3LYP/def2-SVP Opt=CalcFC Freq",
    )

    gjf_content = (out / f"{xyz.stem}_constrained.gjf").read_text()
    route_line = next((line.strip() for line in gjf_content.splitlines() if line.strip().startswith("#p")), "")

    assert "Opt=CalcFC Freq ModRedundant" not in route_line
    assert "Opt=(CalcFC,ModRedundant)" in route_line
