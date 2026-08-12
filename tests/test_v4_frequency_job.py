from pathlib import Path
from types import SimpleNamespace

from rph_core.utils import qc_jobs
from rph_core.utils.qc_models import QCJobSpec


def test_numerical_frequency_job_passes_numfreq_to_orca(monkeypatch, tmp_path: Path):
    input_xyz = tmp_path / "ts.xyz"
    input_xyz.write_text("1\nH\nH 0 0 0\n", encoding="utf-8")
    captured = {}

    class FakeOrcaInterface:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def single_point(self, input_path, output_dir, **kwargs):
            return SimpleNamespace(
                converged=True,
                output_file=Path(output_dir) / "numfreq.out",
                energy=-1.0,
                frequencies=(-250.0, 125.0),
                error_message=None,
            )

        def extract_frequencies_from_output(self, output_file):
            return (-250.0, 125.0)

    monkeypatch.setattr(qc_jobs, "ORCAInterface", FakeOrcaInterface)
    result = qc_jobs.run_frequency(
        QCJobSpec(engine="orca", task="numfreq", method="B97-3c"),
        input_xyz,
        tmp_path / "freq",
        {"resources": {"nproc": 2}},
    )

    assert captured["route_extras"] == "NumFreq"
    assert result.status == "complete"
    assert result.frequencies_cm1 == (-250.0, 125.0)


def test_frequency_job_parses_output_when_interface_result_has_no_frequencies(monkeypatch, tmp_path: Path):
    input_xyz = tmp_path / "ts.xyz"
    input_xyz.write_text("1\nH\nH 0 0 0\n", encoding="utf-8")

    class FakeOrcaInterface:
        def __init__(self, **kwargs):
            pass

        def single_point(self, input_path, output_dir, **kwargs):
            output_file = Path(output_dir) / "freq.out"
            output_file.parent.mkdir(parents=True, exist_ok=True)
            output_file.write_text(
                "ORCA TERMINATED NORMALLY\n"
                "VIBRATIONAL FREQUENCIES\n"
                "  0:   -250.0000 cm**-1\n"
                "  1:    125.0000 cm**-1\n"
                "Zero point energy ... 0.12340000 Eh\n"
                "Total thermal energy -1.01230000 Eh\n"
                "Total Enthalpy ... -1.01130000 Eh\n"
                "Final Gibbs free energy ... -1.04560000 Eh\n"
                "G-E(el) ... -0.04560000 Eh\n",
                encoding="utf-8",
            )
            return SimpleNamespace(
                converged=True,
                output_file=output_file,
                energy=-1.0,
                frequencies=None,
                error_message=None,
            )

        def extract_frequencies_from_output(self, output_file):
            return (-250.0, 125.0)

    monkeypatch.setattr(qc_jobs, "ORCAInterface", FakeOrcaInterface)
    result = qc_jobs.run_frequency(
        QCJobSpec(engine="orca", task="freq", method="B97-3c"),
        input_xyz,
        tmp_path / "freq",
        {"resources": {"nproc": 2}},
    )

    assert result.status == "complete"
    assert result.frequencies_cm1 == (-250.0, 125.0)
    assert result.zero_point_energy_hartree == 0.1234
    assert result.thermal_energy_hartree == -1.0123
    assert result.enthalpy_hartree == -1.0113
    assert result.gibbs_free_energy_hartree == -1.0456
    assert result.gibbs_correction_hartree == -0.0456
