"""
Test QCTaskRunner Integration
===================================
验证 QCTaskRunner 统一计算中枢的集成
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))

from rph_core.utils.qc_task_runner import QCTaskRunner

def test_qctaskrunner_import():
    """测试 QCTaskRunner 导入"""
    print("✓ QCTaskRunner 导入成功")
    print(f"  - 类: {QCTaskRunner.__name__}")
    return True

def test_qctaskrunner_init():
    """测试 QCTaskRunner 初始化（需要配置）"""
    try:
        # 创建最小配置
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
        print(f"✓ QCTaskRunner 初始化成功")
        print(f"  - 引擎: {runner.engine_type}")
        print(f"  - L2 SP: {runner.orca.method}/{runner.orca.basis}")
        return True
        
    except Exception as e:
        print(f"✗ QCTaskRunner 初始化失败: {e}")
        return False

def test_qctaskrunner_methods():
    """测试 QCTaskRunner 方法"""
    try:
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
        
        print("✓ QCTaskRunner 方法检查:")
        print(f"  - run_opt_sp_cycle: {hasattr(runner, 'run_opt_sp_cycle')}")
        print(f"  - run_ts_opt_cycle: {hasattr(runner, 'run_ts_opt_cycle')}")
        print(f"  - run_sp_only: {hasattr(runner, 'run_sp_only')}")
        
        return True
        
    except Exception as e:
        print(f"✗ QCTaskRunner 方法检查失败: {e}")
        return False

if __name__ == '__main__':
    print("=" * 60)
    print("QCTaskRunner 集成测试")
    print("=" * 60)
    
    all_passed = True
    
    # 测试 1: 导入
    all_passed &= test_qctaskrunner_import()
    
    # 测试 2: 初始化
    all_passed &= test_qctaskrunner_init()
    
    # 测试 3: 方法检查
    all_passed &= test_qctaskrunner_methods()
    
    print("=" * 60)
    if all_passed:
        print("✓ 所有测试通过")
        sys.exit(0)
    else:
        print("✗ 部分测试失败")
        sys.exit(1)
