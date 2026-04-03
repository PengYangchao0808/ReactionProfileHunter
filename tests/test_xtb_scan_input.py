from pathlib import Path

from rph_core.utils.data_types import QCResult
from rph_core.utils.qc_interface import XTBInterface
from rph_core.utils.xtb_runner import XTBRunner


def _make_xyz_frame(atom_count: int, energy: float) -> str:
    coord_lines = [f"C {i: .6f} 0.000000 0.000000" for i in range(1, atom_count + 1)]
    return "\n".join([str(atom_count), f" energy: {energy:.7f} xtb: 6.7.1"] + coord_lines)


def test_write_scan_input_uses_configured_force_constant(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(XTBRunner, "_verify_executable", lambda self: "xtb")
    runner = XTBRunner(config={"resources": {"nproc": 1}}, work_dir=tmp_path)

    scan_file = runner._write_scan_input(
        constraints={"0 1": 2.2},
        scan_range=(2.2, 3.5),
        scan_steps=10,
        scan_mode="concerted",
        scan_force_constant=0.5,
    )

    content = scan_file.read_text(encoding="utf-8")

    assert "force constant=0.500" in content
    assert "distance: 1, 2, 2.200" in content


def test_xtb_interface_optimize_converts_zero_based_constraints(tmp_path: Path, monkeypatch) -> None:
    captured = {}

    def _mock_optimize(self, structure, constraints=None, frozen_indices=None, solvent=None, charge=0, uhf=0, **kwargs):
        captured["constraints"] = constraints
        output_file = tmp_path / "xtbopt.xyz"
        output_file.write_text("1\nxtb\nH 0.0 0.0 0.0\n", encoding="utf-8")
        return QCResult(success=True, converged=True, coordinates=output_file, output_file=output_file)

    monkeypatch.setattr(XTBRunner, "_verify_executable", lambda self: "xtb")
    monkeypatch.setattr(XTBRunner, "optimize", _mock_optimize)

    xyz_file = tmp_path / "input.xyz"
    xyz_file.write_text("1\ninput\nH 0.0 0.0 0.0\n", encoding="utf-8")

    xtb = XTBInterface(config={"resources": {"nproc": 1}})
    result = xtb.optimize(
        xyz_file=xyz_file,
        output_dir=tmp_path / "opt",
        constraints={"0 1": 2.2, "2 3": 2.8},
    )

    assert result.success is True
    assert captured["constraints"] == {"1 2": 2.2, "3 4": 2.8}


def test_parse_scan_log_reads_energy_from_nearby_total_energy_line(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(XTBRunner, "_verify_executable", lambda self: "xtb")
    runner = XTBRunner(config={"resources": {"nproc": 1}}, work_dir=tmp_path)

    scan_log = tmp_path / "xtb_scan.log"
    scan_log.write_text(
        "\n".join(
            [
                "  * total energy  :   -64.1034035 Eh     change       -0.3365475E-05 Eh",
                "    gradient norm :     0.0007554 Eh/alpha   predicted    -0.1436597E-05 ( -57.31%)",
                "  total energy gain   :        -1.5792588 Eh     -990.9998 kcal/mol",
                "39",
                " xtb: 6.7.1 (edcfbbe)",
            ]
            + [f"C {i: .6f} 0.000000 0.000000" for i in range(1, 40)]
        )
        + "\n",
        encoding="utf-8",
    )

    result = runner._parse_scan_log(scan_log)

    assert result.success is True
    assert result.energies == [-64.1034035]
    assert result.max_energy_index == 0
    assert result.ts_guess_xyz is not None
    assert Path(result.ts_guess_xyz).exists()


def test_parse_scan_log_prefers_richer_xtbscan_when_both_logs_exist(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(XTBRunner, "_verify_executable", lambda self: "xtb")
    runner = XTBRunner(config={"resources": {"nproc": 1}}, work_dir=tmp_path)

    xtb_scan_log = tmp_path / "xtb_scan.log"
    xtb_scan_log.write_text(
        "\n".join(
            [
                "===============",
                " final structure:",
                "===============",
                _make_xyz_frame(atom_count=3, energy=-64.2679222),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    xtbscan_log = tmp_path / "xtbscan.log"
    xtbscan_log.write_text(
        "\n".join(
            [
                _make_xyz_frame(atom_count=3, energy=-64.3000000),
                _make_xyz_frame(atom_count=3, energy=-64.2000000),
                _make_xyz_frame(atom_count=3, energy=-64.1500000),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = runner._parse_scan_log(xtb_scan_log)

    assert result.success is True
    assert result.scan_log == xtbscan_log
    assert result.energies == [-64.3, -64.2, -64.15]
    assert result.max_energy_index == 2
    assert result.ts_guess_xyz is not None
    assert Path(result.ts_guess_xyz).exists()
