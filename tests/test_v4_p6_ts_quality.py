"""P6 TS quality closure tests."""
from __future__ import annotations

from pathlib import Path

from rph_core.steps import stage_calculator as stage_module
from rph_core.steps.stage_calculator import StageCalculator
from rph_core.utils.qc_models import QCJobResult
from rph_core.utils.ts_quality import analyze_ts_quality


def test_frequency_count_valid_with_one_significant_imaginary():
    result = analyze_ts_quality(
        frequency_output=None,
        optimized_xyz=None,
        forming_bonds=[(0, 1)],
        frequencies_cm1=[-150.0, 100.0, 200.0],
    )
    assert result["frequency_count_valid"] is True
    assert result["significant_imaginary_count"] == 1


def test_frequency_count_invalid_with_no_imaginary():
    result = analyze_ts_quality(
        frequency_output=None,
        optimized_xyz=None,
        forming_bonds=[(0, 1)],
        frequencies_cm1=[100.0, 200.0, 300.0],
    )
    assert result["frequency_count_valid"] is False
    assert result["quality_summary"] == "frequency_count_failed"


def test_frequency_count_invalid_with_two_significant_imaginary():
    result = analyze_ts_quality(
        frequency_output=None,
        optimized_xyz=None,
        forming_bonds=[(0, 1)],
        frequencies_cm1=[-200.0, -150.0, 100.0],
    )
    assert result["frequency_count_valid"] is False


def test_mode_displacement_none_when_no_frequency_output():
    result = analyze_ts_quality(
        frequency_output=None,
        optimized_xyz=None,
        forming_bonds=[(0, 1)],
        frequencies_cm1=[-150.0, 100.0],
    )
    assert result["mode_displacement_valid"] is None
    assert result["quality_summary"] == "frequency_valid_mode_not_checked"


def test_irc_valid_always_none():
    result = analyze_ts_quality(
        frequency_output=None,
        optimized_xyz=None,
        forming_bonds=[(0, 1)],
        frequencies_cm1=[-150.0, 100.0],
    )
    assert result["irc_valid"] is None


def test_mode_displacement_checked_with_geometry(tmp_path: Path):
    xyz = tmp_path / "opt.xyz"
    xyz.write_text(
        "3\ntest\n"
        "C 0.0 0.0 0.0\n"
        "O 1.2 0.0 0.0\n"
        "H 2.0 0.0 0.0\n",
        encoding="utf-8",
    )
    freq_out = tmp_path / "freq.log"
    freq_out.write_text(
        "Harmonic frequencies\n"
        "Frequencies -- -500.0\n"
        "\n"
        "Atom  1   X  0.50  Y  0.0  Z  0.0\n"
        "Atom  2   X  -0.50  Y  0.0  Z  0.0\n"
        "Atom  3   X  0.0  Y  0.0  Z  0.0\n",
        encoding="utf-8",
    )
    result = analyze_ts_quality(
        frequency_output=freq_out,
        optimized_xyz=xyz,
        forming_bonds=[(0, 1)],
        frequencies_cm1=[-500.0, 100.0],
    )
    assert result["frequency_count_valid"] is True
    assert result["mode_displacement_valid"] is True
    assert result["quality_summary"] == "frequency_valid_mode_aligned"
    details = result["mode_displacement_details"]
    assert "bonds" in details
    assert details["bonds"][0]["aligned"] is True


def test_quality_summary_reports_mode_failure(tmp_path: Path):
    xyz = tmp_path / "opt.xyz"
    xyz.write_text(
        "2\ntest\n"
        "C 0.0 0.0 0.0\n"
        "O 1.2 0.0 0.0\n",
        encoding="utf-8",
    )
    freq_out = tmp_path / "freq.log"
    freq_out.write_text(
        "Frequencies -- -500.0\n"
        "Red. masses -- 1.0\n"
        "Atom  AN      X      Y      Z\n"
        "   1   6    0.00   0.50   0.00\n"
        "   2   8    0.00  -0.50   0.00\n",
        encoding="utf-8",
    )
    result = analyze_ts_quality(
        frequency_output=freq_out,
        optimized_xyz=xyz,
        forming_bonds=[(0, 1)],
        frequencies_cm1=[-500.0, 100.0],
    )
    assert result["mode_displacement_valid"] is False
    assert result["quality_summary"] == "frequency_valid_mode_misaligned"


def test_mode_displacement_checked_with_orca_geometry(tmp_path: Path):
    xyz = tmp_path / "opt.xyz"
    xyz.write_text(
        "2\ntest\n"
        "C 0.0 0.0 0.0\n"
        "O 1.2 0.0 0.0\n",
        encoding="utf-8",
    )
    freq_out = tmp_path / "freq.out"
    freq_out.write_text(
        "VIBRATIONAL FREQUENCIES\n"
        "CARTESIAN DISPLACEMENTS\n"
        "-----------\n"
        "Mode:   0\n"
        "Freq:   -321.40 cm**-1 (imaginary mode)\n"
        "        dx          dy          dz\n"
        "Atom 0:  0.45  0.00  0.00\n"
        "Atom 1: -0.45  0.00  0.00\n",
        encoding="utf-8",
    )
    result = analyze_ts_quality(
        frequency_output=freq_out,
        optimized_xyz=xyz,
        forming_bonds=[(0, 1)],
        frequencies_cm1=[-321.4, 100.0],
    )
    assert result["mode_displacement_valid"] is True
    assert result["quality_summary"] == "frequency_valid_mode_aligned"


def test_stage_calculator_returns_ts_quality_fields(monkeypatch, tmp_path: Path):
    input_xyz = tmp_path / "ts.xyz"
    input_xyz.write_text(
        "2\nts\n"
        "C 0.0 0.0 0.0\n"
        "O 1.2 0.0 0.0\n",
        encoding="utf-8",
    )
    freq_out = tmp_path / "freq.out"
    freq_out.write_text(
        "VIBRATIONAL FREQUENCIES\n"
        "CARTESIAN DISPLACEMENTS\n"
        "-----------\n"
        "Mode:   0\n"
        "Freq:   -321.40 cm**-1 (imaginary mode)\n"
        "        dx          dy          dz\n"
        "Atom 0:  0.45  0.00  0.00\n"
        "Atom 1: -0.45  0.00  0.00\n",
        encoding="utf-8",
    )

    def fake_opt(spec, input_path, output_dir, config):
        return QCJobResult(
            "complete",
            input_path,
            input_path,
            output_dir / "opt.out",
            -1.0,
        )

    def fake_frequency(spec, input_path, output_dir, config):
        return QCJobResult(
            "complete",
            input_path,
            output_file=freq_out,
            frequencies_cm1=(-321.4, 120.0, 311.2),
        )

    def fake_sp(spec, input_path, output_dir, config):
        return QCJobResult("complete", input_path, None, output_dir / "sp.out", -1.1)

    monkeypatch.setattr(stage_module, "run_optimization", fake_opt)
    monkeypatch.setattr(stage_module, "run_frequency", fake_frequency)
    monkeypatch.setattr(stage_module, "run_single_point", fake_sp)

    theory = {
        "optimization": {
            "engine": "orca",
            "method": "B97-3c",
            "frequency": {
                "enabled_for_ts": True,
                "task": "freq",
                "imaginary_cutoff_cm1": -50.0,
                "require_exactly_one": True,
            },
        },
        "single_point": {"engine": "orca", "method": "r2SCAN-3c"},
    }

    result = StageCalculator({}, theory).run_structure(
        {
            "id": "ts",
            "kind": "ts",
            "input_xyz": str(input_xyz),
            "forming_bonds": [(0, 1)],
        },
        tmp_path / "stage",
    )

    assert result["ts_frequency_valid"] is True
    assert result["frequency_count_valid"] is True
    assert result["mode_displacement_valid"] is True
    assert result["irc_valid"] is None
    assert result["ts_mode_displacement_verified"] is True
    assert result["ts_quality_summary"] == "frequency_valid_mode_aligned"
