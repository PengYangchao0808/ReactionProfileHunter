"""
测试几何预处理器 (GeometryPreprocessor)
========================================

用于验证 xTB 预优化逻辑是否正常工作

Author: ReactionProfileHunter Team
Date: 2026-01-31
"""

from pathlib import Path
import sys

# 添加项目路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from rph_core.utils.geometry_preprocessor import GeometryPreprocessor
from rph_core.utils.file_io import write_xyz

def create_test_xyz_normal():
    """创建正常的测试分子（无重叠）"""
    coords = [
        [0.0, 0.0, 0.0],      # C
        [1.4, 0.0, 0.0],      # C
        [2.1, 1.2, 0.0],      # H
        [2.1, -0.6, 0.9],     # H
        [2.1, -0.6, -0.9],    # H
        [-0.5, 0.9, 0.0],     # H
        [-0.5, -0.5, 0.9],    # H
        [-0.5, -0.5, -0.9],   # H
    ]
    symbols = ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
    
    test_file = project_root / "test_normal.xyz"
    write_xyz(test_file, coords, symbols, title="Normal ethane (no overlap)")
    return test_file

def create_test_xyz_overlap():
    """创建原子重叠的测试分子"""
    coords = [
        [0.0, 0.0, 0.0],      # C
        [0.5, 0.0, 0.0],      # C (故意设置太近，造成重叠)
        [0.8, 0.8, 0.0],      # H
        [0.8, -0.4, 0.7],     # H
        [0.8, -0.4, -0.7],    # H
        [-0.3, 0.7, 0.0],     # H
        [-0.3, -0.4, 0.7],    # H
        [-0.3, -0.4, -0.7],   # H
    ]
    symbols = ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
    
    test_file = project_root / "test_overlap.xyz"
    write_xyz(test_file, coords, symbols, title="Overlapping ethane (C-C = 0.5 Å)")
    return test_file

def test_preprocessor():
    """测试预处理器"""
    print("=" * 60)
    print("测试几何预处理器")
    print("=" * 60)
    
    # 加载配置（简化版）
    config = {
        'theory': {
            'preoptimization': {
                'enabled': True,
                'gfn_level': 2,
                'solvent': 'acetone',
                'nproc': 4,
                'overlap_threshold': 1.0,
                'opt_level': 'crude'
            }
        },
        'executables': {
            'xtb': {
                'path': 'xtb'  # 假设 xTB 在系统 PATH 中
            }
        },
        'resources': {
            'nproc': 4
        }
    }
    
    preprocessor = GeometryPreprocessor(config)
    
    # 测试 1: 正常分子（无重叠）
    print("\n测试 1: 正常分子（无重叠）")
    print("-" * 60)
    normal_xyz = create_test_xyz_normal()
    result1 = preprocessor.preprocess(
        xyz_file=normal_xyz,
        output_dir=project_root / "test_output" / "normal"
    )
    print(f"结果: success={result1.success}")
    if hasattr(result1, 'skip_reason'):
        print(f"跳过原因: {result1.skip_reason}")
    
    # 测试 2: 重叠分子
    print("\n测试 2: 重叠分子（需要预优化）")
    print("-" * 60)
    overlap_xyz = create_test_xyz_overlap()
    result2 = preprocessor.preprocess(
        xyz_file=overlap_xyz,
        output_dir=project_root / "test_output" / "overlap"
    )
    print(f"结果: success={result2.success}")
    if result2.success and hasattr(result2, 'preopt_performed'):
        print(f"是否执行了预优化: {result2.preopt_performed}")
        if result2.preopt_performed:
            print(f"优化后的结构: {result2.coordinates}")
    
    print("\n" + "=" * 60)
    print("测试完成！")
    print("=" * 60)

if __name__ == "__main__":
    try:
        test_preprocessor()
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
