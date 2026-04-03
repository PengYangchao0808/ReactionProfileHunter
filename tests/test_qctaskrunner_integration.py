"""
Test QCTaskRunner Integration
==================================
验证 QCTaskRunner 统一计算中枢的集成
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))

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


if __name__ == '__main__':
    import pytest as _pytest
    _pytest.main([__file__, '-v'])
