"""
资源配置解析测试 — resolve_resources()
======================================

测试资源继承、预算分配、自动检测功能。
本文件不依赖 RDKit，可在 base Python 环境中运行。
"""

import os
import sys
import json

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rph_core.utils.resource_utils import (
    resolve_resources,
    calc_orca_maxcore,
    mem_to_mb,
)


class TestResolveResources:
    """resolve_resources() 功能测试"""

    def test_returns_config(self):
        """resolve_resources() 必须返回 config 字典（不是 None）。"""
        cfg = {'resources': {'nproc': 16, 'mem': '32GB'}}
        result = resolve_resources(cfg)
        assert result is not None, "返回值为 None"
        assert result is cfg, "返回值不是同一对象"

    def test_resources_unchanged(self):
        """resources.nproc/mem 不为 None 时不覆盖。"""
        cfg = {'resources': {'nproc': 8, 'mem': '16GB'}}
        resolve_resources(cfg)
        assert cfg['resources']['nproc'] == 8
        assert cfg['resources']['mem'] == '16GB'

    def test_inherit_theory_optimization(self):
        """theory.optimization.nproc 为 None 时继承 resources.nproc。"""
        cfg = {'resources': {'nproc': 24, 'mem': '64GB'}, 'theory': {'optimization': {'nproc': None}}}
        resolve_resources(cfg)
        assert cfg['theory']['optimization']['nproc'] == 24

    def test_inherit_theory_single_point(self):
        """theory.single_point.nproc 为 None 时继承 resources.nproc。"""
        cfg = {'resources': {'nproc': 16, 'mem': '32GB'}, 'theory': {'single_point': {'nproc': None}}}
        resolve_resources(cfg)
        assert cfg['theory']['single_point']['nproc'] == 16

    def test_override_not_overwritten(self):
        """显式 theory.optimization.nproc=8 不被 resources.nproc=16 覆盖。"""
        cfg = {'resources': {'nproc': 16, 'mem': '32GB'}, 'theory': {'optimization': {'nproc': 8}}}
        resolve_resources(cfg)
        assert cfg['theory']['optimization']['nproc'] == 8, "显式 override 被覆盖"

    def test_inherit_step2_scan(self):
        """step2.scan.nproc 为 None 时继承。"""
        cfg = {'resources': {'nproc': 16, 'mem': '32GB'}, 'step2': {'scan': {'nproc': None}}}
        resolve_resources(cfg)
        assert cfg['step2']['scan']['nproc'] == 16

    def test_inherit_step1_crest_threads(self):
        """step1.crest.threads 为 None 时继承。"""
        cfg = {
            'resources': {'nproc': 16, 'mem': '32GB'},
            'step1': {'crest': {'threads': None}},
        }
        resolve_resources(cfg)
        assert cfg['step1']['crest']['threads'] == 16


class TestIrpDefaults:
    """intra_reaction_parallel 默认值填充测试"""

    def test_irp_s1_crest_cores(self):
        """s1.crest_cores 为 None 时自动填充。"""
        cfg = {
            'resources': {'nproc': 16, 'mem': '32GB'},
            'intra_reaction_parallel': {},
        }
        resolve_resources(cfg)
        irp = cfg['intra_reaction_parallel']
        assert irp['s1']['crest_cores'] is not None
        assert irp['s1']['crest_cores'] <= 16

    def test_irp_s3_ts_inter_split(self):
        """parallel_intermediate_ts=true 时 TS 和 Inter 按比例拆分内存。"""
        cfg = {
            'resources': {'nproc': 16, 'mem': '32GB'},
            'intra_reaction_parallel': {
                's3': {'parallel_intermediate_ts': True},
            },
        }
        resolve_resources(cfg)
        s3 = cfg['intra_reaction_parallel']['s3']
        ts_mem = int(s3['ts_opt_mem'].replace('GB', ''))
        inter_mem = int(s3['intermediate_opt_mem'].replace('GB', ''))
        assert ts_mem >= inter_mem, f"TS mem {ts_mem} < Inter mem {inter_mem}"
        assert ts_mem + inter_mem <= 32, f"Sum mem {ts_mem+inter_mem} > 32GB"

    def test_irp_s3_sp_maxcore_budget(self):
        """SP maxcore 不超过 safety budget。"""
        cfg = {
            'resources': {'nproc': 16, 'mem': '32GB', 'orca_maxcore_safety': 0.65},
            'intra_reaction_parallel': {'s3': {'sp_workers': 2}},
        }
        resolve_resources(cfg)
        s3 = cfg['intra_reaction_parallel']['s3']
        sp_max = s3['sp_maxcore_per_job']
        sp_cores = s3['sp_cores_per_job']
        sp_workers = s3.get('sp_workers', 1)
        total = sp_max * sp_cores * sp_workers / 1024
        budget = 32 * 0.65
        assert total <= budget * 1.01, f"SP {total:.1f}GB > budget {budget:.1f}GB"

    def test_irp_total_cores_back_written(self):
        """自动检测后 total_cores 回写。"""
        cfg = {
            'resources': {'nproc': 16, 'mem': '32GB'},
            'intra_reaction_parallel': {},
        }
        resolve_resources(cfg)
        irp = cfg['intra_reaction_parallel']
        assert irp.get('total_cores') == 16
        assert irp.get('total_mem') == '32GB'

    def test_irp_s3_serial_no_split(self):
        """parallel_intermediate_ts=false 时 TS 和 Inter 都用 full nproc。"""
        cfg = {
            'resources': {'nproc': 16, 'mem': '32GB'},
            'intra_reaction_parallel': {
                's3': {'parallel_intermediate_ts': False},
            },
        }
        resolve_resources(cfg)
        s3 = cfg['intra_reaction_parallel']['s3']
        assert s3['ts_opt_cores'] == 16
        assert s3['intermediate_opt_cores'] == 16


class TestCalcOrcaMaxcore:
    """calc_orca_maxcore 计算测试"""

    def test_default_safety(self):
        """默认 safety_factor=0.65。"""
        maxcore = calc_orca_maxcore('32GB', 16)
        assert maxcore == 1331, f"Expected 1331, got {maxcore}"

    def test_custom_safety(self):
        """自定义 safety_factor。"""
        maxcore = calc_orca_maxcore('32GB', 16, 0.8)
        assert maxcore == 1638, f"Expected 1638, got {maxcore}"

    def test_different_mem(self):
        """不同内存配置。"""
        maxcore = calc_orca_maxcore('64GB', 16, 0.65)
        assert maxcore == 2662, f"Expected 2662, got {maxcore}"

    def test_single_core(self):
        """单核。"""
        maxcore = calc_orca_maxcore('32GB', 1, 0.65)
        assert maxcore == 21299, f"Expected 21299, got {maxcore}"


class TestMemToMb:
    """mem_to_mb 解析测试"""

    def test_gb(self):
        """32GB → 32768 MB"""
        assert mem_to_mb('32GB') == 32768

    def test_numeric_gb(self):
        """不带单位的纯数字视为 GB。"""
        assert mem_to_mb('16') == 16384

    def test_mb(self):
        """4096MB → 4096 MB"""
        assert mem_to_mb('4096MB') == 4096
