import subprocess
from pathlib import Path

from rph_core.steps.conformer_search import xtb_thermo


def test_xtb_thermo_uses_configured_bhess_and_removes_oversized_stack(
    monkeypatch, tmp_path: Path
):
    coord = tmp_path / "input.coord"
    coord.write_text("$coord\n0.0 0.0 0.0 h\n$end\n", encoding="utf-8")
    captured = {}

    def fake_run(cmd, *, cwd, capture_output, text, timeout, env):
        captured["cmd"] = cmd
        captured["env"] = env
        (Path(cwd) / "xtb_enso.json").write_text('{"G(T)": 0.125}', encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setenv("OMP_STACKSIZE", "4G")
    monkeypatch.setattr(xtb_thermo.subprocess, "run", fake_run)

    result = xtb_thermo.run_xtb_enso(
        xtb_bin=Path("xtb"),
        coord_file=coord,
        output_dir=tmp_path / "run",
        nproc=4,
        gfn_level=1,
        bhess_level="normal",
        omp_stacksize=None,
    )

    assert result.success is True
    assert result.g_rrho_correction_hartree == 0.125
    assert result.g_total == 0.125
    assert "--bhess" in captured["cmd"]
    assert captured["cmd"][captured["cmd"].index("--bhess") + 1] == "normal"
    assert captured["cmd"][captured["cmd"].index("--parallel") + 1] == "4"
    assert captured["env"]["OMP_NUM_THREADS"] == "4"
    assert "OMP_STACKSIZE" not in captured["env"]
