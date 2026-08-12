"""P6 TS quality closure tests."""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import rph_core.steps.refinement.engine as engine_module
from rph_core.steps.fidelity_profile import FidelityProfile
from rph_core.steps.refinement import RefinementEngine
from rph_core.utils.config_loader import load_config
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


def test_refinement_engine_returns_ts_classification_fields(monkeypatch, tmp_path: Path):
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

    def fake_opt(spec, input_path, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        return QCJobResult(
            "complete",
            Path(input_path),
            output_xyz=Path(input_path),
            output_file=Path(output_dir) / "opt.out",
            energy_hartree=-1.0,
        )

    def fake_frequency(spec, input_path, output_dir, config, subprocess_callback=None):
        del spec, output_dir, config, subprocess_callback
        return QCJobResult(
            "complete",
            Path(input_path),
            output_file=freq_out,
            frequencies_cm1=(-321.4, 120.0, 311.2),
        )

    def fake_sp(spec, input_path, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        return QCJobResult(
            "complete",
            Path(input_path),
            output_file=Path(output_dir) / "sp.out",
            energy_hartree=-1.1,
        )

    monkeypatch.setattr(engine_module, "run_optimization", fake_opt)
    monkeypatch.setattr(engine_module, "run_frequency", fake_frequency)
    monkeypatch.setattr(engine_module, "run_single_point", fake_sp)
    monkeypatch.setattr(
        engine_module.identity_module,
        "classify_ts",
        lambda **kwargs: {
            "hessian_index": 1,
            "curvature_class": "strict",
            "mode_identity": "target",
            "stationary_point_class": "valid_target_ts",
            "alignment_score": 0.9,
        },
    )
    config = load_config()
    profile = replace(FidelityProfile.from_config(config, "S3"), irc_enabled=False)
    manifest_path = RefinementEngine(config, profile).run(
        [
            {
                "id": "ts",
                "role": "ts",
                "kind": "ts",
                "input_xyz": str(input_xyz),
                "forming_bonds": [(0, 1)],
            }
        ],
        tmp_path / "stage",
    )
    result = json.loads(manifest_path.read_text(encoding="utf-8"))["structures"][0]

    assert result["frequency_status"] == "complete"
    assert result["imaginary_frequencies_cm1"] == [-321.4]
    assert result["ts_classification"]["stationary_point_class"] == "valid_target_ts"
    assert result["ts_classification"]["alignment_score"] > 0.0
