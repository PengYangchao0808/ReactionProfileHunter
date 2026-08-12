from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from rph_core.utils import qc_jobs
from rph_core.utils.method_registry import MethodRegistry
from rph_core.utils.orca_input_renderer import OrcaInputRenderer
from rph_core.utils.orca_interface import ORCAInterface, _render_irc_block
from rph_core.utils.qc_interface import TaskKind
from rph_core.utils.qc_models import IRCJobSpec, QCJobResult, QCJobSpec


def _xyz_frame(x: float) -> str:
    return f"2\nframe\nH 0.0 0.0 0.0\nH {x:.6f} 0.0 0.0\n"


def _renderer() -> OrcaInputRenderer:
    spec = MethodRegistry.normalize_spec(
        {"engine": "orca", "method": "B97-3c", "basis": ""},
        default_engine="orca",
        default_basis="",
    )
    return OrcaInputRenderer(spec)


def test_render_irc_block_default() -> None:
    block = _render_irc_block(IRCJobSpec())

    assert "MaxIter 5" in block
    assert "Direction both" in block
    assert "InitHess read" in block
    assert "HessMode 0" in block
    assert "DE_INIT_DISPL 2.0" in block
    assert "SCALE_INIT_DISPL 0.1" in block
    assert "SCALE_DISPL_SD 0.15" in block
    assert "ADAPT_SCALE_DISPL true" in block


def test_render_irc_block_with_hessian_filename() -> None:
    block = _render_irc_block(IRCJobSpec(hessian_filename="ts.hess"))

    assert 'Hess_Filename "ts.hess"' in block


def test_render_irc_block_forward_direction() -> None:
    block = _render_irc_block(IRCJobSpec(direction="forward"))

    assert "Direction forward" in block


def test_render_irc_block_length_displacement_mode() -> None:
    block = _render_irc_block(IRCJobSpec(init_displacement_mode="length"))

    assert "DE_INIT_DISPL" not in block


def test_render_irc_block_no_adaptive_step() -> None:
    block = _render_irc_block(IRCJobSpec(adaptive_step=False))

    assert "ADAPT_SCALE_DISPL" not in block


def test_render_irc_block_has_end_token() -> None:
    assert _render_irc_block(IRCJobSpec()).splitlines()[-1] == "end"


def test_render_irc_block_max_iter_custom() -> None:
    assert "MaxIter 100" in _render_irc_block(IRCJobSpec(max_iter=100))


def test_taskkind_irc_exists() -> None:
    assert TaskKind.IRC.value == "irc"


def test_orca_input_renderer_irc_keyword() -> None:
    assert "IRC" in _renderer().render_simple_keywords("irc")


def test_orca_interface_run_irc_method_exists() -> None:
    assert hasattr(ORCAInterface, "run_irc")


def test_orca_interface_run_irc_parses_endpoints(monkeypatch, tmp_path: Path) -> None:
    input_xyz = tmp_path / "ts.xyz"
    input_xyz.write_text(_xyz_frame(0.9), encoding="utf-8")
    interface = ORCAInterface(
        method="B97-3c",
        basis="",
        aux_basis="",
        nprocs=1,
        solvent="",
        orca_binary_path="/bin/true",
        config={},
    )

    def fake_run_orca(self, inp_file: Path, output_dir: Path, timeout=None, subprocess_callback=None) -> Path:
        content = inp_file.read_text(encoding="utf-8")
        assert "! IRC B97-3c" in content
        assert "%irc" in content
        assert "Direction both" in content
        out_file = inp_file.with_suffix(".out")
        out_file.write_text(
            "FINAL SINGLE POINT ENERGY -123.456789\nORCA TERMINATED NORMALLY\n",
            encoding="utf-8",
        )
        (output_dir / "irc_forward_path.xyz").write_text(
            _xyz_frame(1.1) + _xyz_frame(1.5),
            encoding="utf-8",
        )
        (output_dir / "irc_backward_path.xyz").write_text(
            _xyz_frame(0.8) + _xyz_frame(0.6),
            encoding="utf-8",
        )
        return out_file

    monkeypatch.setattr(ORCAInterface, "_run_orca", fake_run_orca)
    result = interface.run_irc(
        spec=IRCJobSpec(method="B97-3c"),
        input_xyz=input_xyz,
        output_dir=tmp_path / "irc",
        charge=0,
        spin=1,
    )

    assert result.error is None, result.error
    assert result.status == "complete"
    assert result.output_xyz == input_xyz
    assert result.extra is not None
    endpoint_a = Path(result.extra["endpoint_a_xyz"])
    endpoint_b = Path(result.extra["endpoint_b_xyz"])
    assert endpoint_a.exists()
    assert endpoint_b.exists()
    assert "1.500000" in endpoint_a.read_text(encoding="utf-8")
    assert "0.600000" in endpoint_b.read_text(encoding="utf-8")
    assert result.extra["irc_trajectory_dir"] == tmp_path / "irc"


def test_orca_interface_run_irc_fails_without_endpoints(monkeypatch, tmp_path: Path) -> None:
    input_xyz = tmp_path / "ts.xyz"
    input_xyz.write_text(_xyz_frame(0.9), encoding="utf-8")
    interface = ORCAInterface(
        method="B97-3c",
        basis="",
        aux_basis="",
        nprocs=1,
        solvent="",
        orca_binary_path="/bin/true",
        config={},
    )

    def fake_run_orca(self, inp_file: Path, output_dir: Path, timeout=None, subprocess_callback=None) -> Path:
        out_file = inp_file.with_suffix(".out")
        out_file.write_text(
            "FINAL SINGLE POINT ENERGY -123.456789\nORCA TERMINATED NORMALLY\n",
            encoding="utf-8",
        )
        return out_file

    monkeypatch.setattr(ORCAInterface, "_run_orca", fake_run_orca)
    result = interface.run_irc(
        spec=IRCJobSpec(method="B97-3c"),
        input_xyz=input_xyz,
        output_dir=tmp_path / "irc",
        charge=0,
        spin=1,
    )

    assert result.status == "failed"
    assert result.error is not None
    assert "Unable to locate ORCA IRC endpoint geometries" in result.error


def test_neb_input_uses_sandbox_local_xyz_names(monkeypatch, tmp_path: Path) -> None:
    start_xyz = tmp_path / "start.xyz"
    end_xyz = tmp_path / "end.xyz"
    guess_xyz = tmp_path / "guess.xyz"
    for xyz_path in (start_xyz, end_xyz, guess_xyz):
        xyz_path.write_text(_xyz_frame(1.0), encoding="utf-8")
    interface = ORCAInterface(
        method="B97-3c", basis="", aux_basis="", nprocs=1, solvent="",
        orca_binary_path="/bin/true", config={},
    )

    def fake_run_orca(self, inp_file: Path, output_dir: Path, timeout=None, subprocess_callback=None) -> Path:
        content = inp_file.read_text(encoding="utf-8")
        assert str(output_dir) not in content
        assert 'NEB_END_XYZFILE "end.xyz"' not in content
        assert "NEB_END_XYZFILE \"start_neb_ts_" in content
        assert "_end.xyz\"" in content
        assert "NEB_TS_XYZFILE \"start_neb_ts_" in content
        assert "_guess.xyz\"" in content
        assert "* xyzfile 0 1 start_neb_ts_" in content
        assert '* xyzfile 0 1 "' not in content
        candidate = output_dir / "mock_NEB-CI_converged.xyz"
        candidate.write_text(_xyz_frame(1.1), encoding="utf-8")
        out_file = inp_file.with_suffix(".out")
        out_file.write_text("ORCA TERMINATED NORMALLY\n", encoding="utf-8")
        return out_file

    monkeypatch.setattr(ORCAInterface, "_run_orca", fake_run_orca)
    result = interface.run_neb_ts(
        QCJobSpec(
            engine="orca", task="neb_ts", method="B97-3c", charge=0, multiplicity=1,
            neb_end_xyz=end_xyz, neb_ts_guess_xyz=guess_xyz,
        ),
        start_xyz,
        tmp_path / "work dir with spaces",
    )

    assert result.error is None, result.error
    assert result.status == "complete"
    assert result.output_xyz is not None
    assert result.output_xyz.name.endswith("_NEB-CI_converged.xyz")
    assert result.extra["candidate_kind"] == "neb_ci_converged"


def test_neb_sandbox_copies_start_end_and_guess_xyz(monkeypatch, tmp_path: Path) -> None:
    output_dir = tmp_path / "work dir with spaces"
    output_dir.mkdir()
    for name in ("start.xyz", "end.xyz", "guess.xyz"):
        (output_dir / name).write_text(_xyz_frame(1.0), encoding="utf-8")
    inp_file = output_dir / "neb.inp"
    inp_file.write_text(
        '%neb\n  NEB_END_XYZFILE "end.xyz"\n  NEB_TS_XYZFILE "guess.xyz"\nend\n'
        '* xyzfile 0 1 start.xyz\n',
        encoding="utf-8",
    )
    interface = ORCAInterface(
        method="B97-3c", basis="", aux_basis="", nprocs=1, solvent="",
        orca_binary_path="/bin/true", config={},
    )
    interface.orca_binary = Path("fake-orca")

    def fake_execute(*, cmd, cwd, out_file, env, timeout, subprocess_callback):
        assert Path(cwd) != output_dir
        for name in ("start.xyz", "end.xyz", "guess.xyz"):
            assert (Path(cwd) / name).is_file()
        Path(out_file).write_text("ORCA TERMINATED NORMALLY\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, timed_out=False, stderr_text="")

    monkeypatch.setattr(interface, "_build_orca_runtime_env", lambda: {})
    monkeypatch.setattr(interface, "_execute_orca_process", fake_execute)
    output = interface._run_orca(inp_file, output_dir)

    assert output.is_file()


def test_qc_jobs_run_irc_dispatches_to_interface(tmp_path: Path) -> None:
    output_dir = tmp_path / "dispatch"
    expected = QCJobResult(status="complete", input_xyz=Path("/tmp/ts.xyz"), output_xyz=Path("/tmp/ts.xyz"))

    with patch.object(qc_jobs.ORCAInterface, "run_irc", return_value=expected) as mocked_run_irc:
        result = qc_jobs.run_irc(IRCJobSpec(method="B97-3c"), Path("/tmp/ts.xyz"), output_dir, {})

    assert result is expected
    mocked_run_irc.assert_called_once()
    _, kwargs = mocked_run_irc.call_args
    assert kwargs["spec"].method == "B97-3c"
    assert kwargs["input_xyz"] == Path("/tmp/ts.xyz").resolve()
    assert kwargs["output_dir"] == output_dir
    assert kwargs["charge"] == 0
    assert kwargs["spin"] == 1


def test_qc_jobs_run_irc_propagates_failure(tmp_path: Path) -> None:
    failed = QCJobResult(status="failed", input_xyz=Path("/tmp/ts.xyz"), error="boom")

    with patch.object(qc_jobs.ORCAInterface, "run_irc", return_value=failed):
        result = qc_jobs.run_irc(IRCJobSpec(method="B97-3c"), Path("/tmp/ts.xyz"), tmp_path / "failure", {})

    assert result.status == "failed"
    assert result.error == "boom"


def test_qc_jobs_run_irc_creates_output_dir(tmp_path: Path) -> None:
    output_dir = tmp_path / "new-output-dir"
    expected = QCJobResult(status="complete", input_xyz=Path("/tmp/ts.xyz"), output_xyz=Path("/tmp/ts.xyz"))

    with patch.object(qc_jobs.ORCAInterface, "run_irc", return_value=expected):
        qc_jobs.run_irc(IRCJobSpec(method="B97-3c"), Path("/tmp/ts.xyz"), output_dir, {})

    assert output_dir.is_dir()
