import subprocess
from pathlib import Path

import pytest

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
        (Path(cwd) / "xtb_enso.json").write_text(
            '{"temperature": 298.15, "energy": -10.0, '
            '"free energy": -9.875, "G(T)": 0.125}',
            encoding="utf-8",
        )
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
        max_scc_iterations=1000,
    )

    assert result.success is True
    assert result.g_rrho_correction_hartree == 0.125
    assert result.xtb_electronic_energy_hartree == -10.0
    assert result.xtb_total_free_energy_hartree == -9.875
    assert result.energy_ledger_residual_hartree == 0.0
    assert result.xtb_solvation_energy_hartree is None
    assert result.explicit_solvation_correction_added is False
    assert result.g_total == 0.125
    assert "--bhess" in captured["cmd"]
    assert captured["cmd"][captured["cmd"].index("--bhess") + 1] == "normal"
    assert captured["cmd"][captured["cmd"].index("--parallel") + 1] == "4"
    assert captured["cmd"][captured["cmd"].index("--iterations") + 1] == "1000"
    assert captured["env"]["OMP_NUM_THREADS"] == "4"
    assert "OMP_STACKSIZE" not in captured["env"]


def test_xtb_thermo_rejects_inconsistent_energy_ledger(monkeypatch, tmp_path: Path):
    coord = tmp_path / "input.coord"
    coord.write_text("$coord\n0.0 0.0 0.0 h\n$end\n", encoding="utf-8")

    def fake_run(cmd, *, cwd, capture_output, text, timeout, env):
        del capture_output, text, timeout, env
        (Path(cwd) / "xtb_enso.json").write_text(
            '{"energy": -10.0, "free energy": -9.8, "G(T)": 0.125}',
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(xtb_thermo.subprocess, "run", fake_run)
    result = xtb_thermo.run_xtb_enso(
        Path("xtb"), coord, tmp_path / "run", ledger_tolerance_hartree=1.0e-7
    )

    assert result.success is False
    assert result.g_rrho_correction_hartree is None
    assert result.energy_ledger_residual_hartree == pytest.approx(0.075)
    assert "Inconsistent xTB energy ledger" in str(result.error)


def test_xtb_thermo_rejects_non_finite_energy_ledger(monkeypatch, tmp_path: Path):
    coord = tmp_path / "input.coord"
    coord.write_text("$coord\n0.0 0.0 0.0 h\n$end\n", encoding="utf-8")

    def fake_run(cmd, *, cwd, capture_output, text, timeout, env):
        del capture_output, text, timeout, env
        (Path(cwd) / "xtb_enso.json").write_text(
            '{"energy": -10.0, "free energy": NaN, "G(T)": NaN}',
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(xtb_thermo.subprocess, "run", fake_run)
    result = xtb_thermo.run_xtb_enso(Path("xtb"), coord, tmp_path / "run")

    assert result.success is False
    assert result.g_rrho_correction_hartree is None
    assert result.error == "Non-finite value in xTB energy ledger"
