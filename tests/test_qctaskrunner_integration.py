"""
Test QCTaskRunner Integration
==================================
验证 QCTaskRunner 统一计算中枢的集成
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))

from rph_core.utils.data_types import QCResult
from rph_core.utils.qc_task_runner import QCTaskRunner


def test_qctaskrunner_import():
    """测试 QCTaskRunner 导入"""
    assert QCTaskRunner is not None
    assert QCTaskRunner.__name__ == "QCTaskRunner"


def test_qctaskrunner_init():
    """测试 QCTaskRunner 初始化（需要配置）"""
    config = {
        'theory': {
            'optimization': {
                'method': 'B3LYP',
                'basis': 'def2-SVP',
                'dispersion': 'GD3BJ',
                'engine': 'gaussian',
                'nproc': 16,
                'mem': '32GB'
            },
            'single_point': {
                'method': 'WB97M-V',
                'basis': 'def2-TZVPP',
                'aux_basis': 'def2/J',
                'nproc': 16,
                'maxcore': 4000,
                'solvent': 'acetone'
            }
        },
        'optimization_control': {
            'hessian': {
                'initial': 'calcfc',
                'recalc_every': 10
            },
            'step': {},
            'global': {
                'log_level': 'INFO'
            }
        }
    }
    
    runner = QCTaskRunner(config=config)
    assert runner is not None
    assert hasattr(runner, 'engine_type')


def test_qctaskrunner_methods():
    """测试 QCTaskRunner 方法"""
    config = {
        'theory': {
            'optimization': {
                'method': 'B3LYP',
                'basis': 'def2-SVP',
                'dispersion': 'GD3BJ',
                'engine': 'gaussian',
                'nproc': 16,
                'mem': '32GB'
            },
            'single_point': {
                'method': 'WB97M-V',
                'basis': 'def2-TZVPP',
                'aux_basis': 'def2/J',
                'nproc': 16,
                'maxcore': 4000,
                'solvent': 'acetone'
            }
        }
    }
    
    runner = QCTaskRunner(config=config)
    assert hasattr(runner, 'run_opt_sp_cycle')
    assert hasattr(runner, 'run_ts_opt_cycle')
    assert hasattr(runner, 'run_sp_only')


def test_qctaskrunner_run_sp_only_uses_factory_single_point_interface(monkeypatch, tmp_path: Path):
    xyz = tmp_path / "mol.xyz"
    xyz.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")

    calls = {}

    class _FakeSPInterface:
        def single_point(self, xyz_file, output_dir, charge=None, spin=None):
            calls["single_point"] = {
                "xyz_file": xyz_file,
                "output_dir": output_dir,
                "charge": charge,
                "spin": spin,
            }
            return QCResult(
                success=True,
                energy=-123.456,
                converged=True,
                output_file=output_dir / "gaussian_sp.log",
                log_file=output_dir / "gaussian_sp.log",
                chk_file=output_dir / "gaussian_sp.chk",
                fchk_file=output_dir / "gaussian_sp.fchk",
                qm_output_file=output_dir / "gaussian_sp.log",
            )

    def _fake_create_interface(engine_type, **kwargs):
        calls["engine_type"] = engine_type
        calls["kwargs"] = kwargs
        return _FakeSPInterface()

    monkeypatch.setattr(
        "rph_core.utils.qc_task_runner.QCInterfaceFactory.create_interface",
        _fake_create_interface,
    )

    config = {
        'theory': {
            'optimization': {
                'method': 'B3LYP',
                'basis': 'def2-SVP',
                'dispersion': 'GD3BJ',
                'engine': 'gaussian',
            },
            'single_point': {
                'method': 'M06-2X',
                'basis': 'def2-TZVPP',
                'engine': 'gaussian',
                'route_extras': 'VeryTightSCF',
                'solvent': 'acetone',
            },
        },
        'resources': {'nproc': 8, 'mem': '16GB'},
    }

    runner = QCTaskRunner(config=config)
    result = runner.run_sp_only(xyz, tmp_path / "sp_out", charge=1, spin=2)

    assert calls["engine_type"] == "gaussian"
    assert calls["kwargs"]["route_extras"] == "VeryTightSCF"
    assert "normalized_spec" in calls["kwargs"]
    assert calls["single_point"]["charge"] == 1
    assert calls["single_point"]["spin"] == 2
    assert result.converged is True
    assert result.energy == -123.456
    assert result.log_file == tmp_path / "sp_out" / "gaussian_sp.log"
    assert result.qm_output_file == tmp_path / "sp_out" / "gaussian_sp.log"
    assert result.fchk_file == tmp_path / "sp_out" / "gaussian_sp.fchk"


def test_qctaskrunner_run_l2_sp_uses_generic_interface(monkeypatch, tmp_path: Path):
    xyz = tmp_path / "mol.xyz"
    xyz.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")

    class _FakeSPInterface:
        def single_point(self, xyz_file, output_dir, charge=None, spin=None):
            _ = (xyz_file, charge, spin)
            return QCResult(
                success=True,
                energy=-77.0,
                converged=True,
                output_file=output_dir / "orca_sp.out",
                qm_output_file=output_dir / "orca_sp.out",
            )

    monkeypatch.setattr(
        "rph_core.utils.qc_task_runner.QCInterfaceFactory.create_interface",
        lambda engine_type, **kwargs: _FakeSPInterface(),
    )

    config = {
        'theory': {
            'optimization': {
                'method': 'B3LYP',
                'basis': 'def2-SVP',
                'dispersion': 'GD3BJ',
                'engine': 'gaussian',
            },
            'single_point': {
                'method': 'DLPNO-CCSD(T)',
                'basis': 'def2-TZVPP',
                'engine': 'orca',
                'route_extras': 'TightPNO',
                'solvent': 'acetone',
            },
        },
    }

    runner = QCTaskRunner(config=config)
    result = runner._run_l2_sp(xyz, tmp_path / "l2_sp", charge=0, spin=1)

    assert result.converged is True
    assert result.energy == -77.0
    assert result.output_file == tmp_path / "l2_sp" / "orca_sp.out"
    assert result.qm_output_file == tmp_path / "l2_sp" / "orca_sp.out"


if __name__ == '__main__':
    import pytest as _pytest
    _pytest.main([__file__, '-v'])
