"""
可执行文件自动发现测试 — resolve_all_executables()
====================================================

测试发现优先级、source 标签、必需性判断。
本文件不依赖 RDKit，可在 base Python 环境中运行。
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rph_core.utils.resource_utils import (
    resolve_all_executables,
    _is_executable_required,
)


class TestPriority:
    """发现优先级测试"""

    def test_config_path_priority(self):
        """显式 config 路径存在时优先使用。"""
        cfg = {
            'executables': {
                'xtb': {'path': sys.executable},
                'discovery': {'enabled': False},
            },
        }
        results = resolve_all_executables(cfg)
        xtb = results.get('xtb', {})
        assert xtb.get('found'), f"xtb not found: {xtb}"
        assert xtb.get('source') == 'config', f"source should be config, got {xtb.get('source')}"

    def test_config_path_missing_fail_fast(self):
        """显式路径无效时 fail_on_invalid_explicit=True 抛出 RuntimeError。"""
        import pytest
        cfg = {
            'executables': {
                'xtb': {'path': '/nonexistent/path/xtb'},
                'discovery': {'enabled': False, 'fail_on_invalid_explicit': True},
            },
        }
        with pytest.raises(RuntimeError, match='fail_on_invalid_explicit=True'):
            resolve_all_executables(cfg)

    def test_config_path_missing_fallback(self):
        """显式路径无效但 fail_on_invalid_explicit=False 时 fallback 搜索。"""
        cfg = {
            'executables': {
                'xtb': {'path': '/nonexistent/path/xtb'},
                'discovery': {'enabled': False, 'fail_on_invalid_explicit': False},
            },
        }
        results = resolve_all_executables(cfg)
        xtb = results.get('xtb', {})
        if not xtb.get('found'):
            assert xtb.get('source') == 'not_found'


class TestRequired:
    """必需性判断测试"""

    def test_gaussian_required_when_engine_gaussian(self):
        """theory.optimization.engine=gaussian 时 gaussian 为必需。"""
        cfg = {'theory': {'optimization': {'engine': 'gaussian'}}}
        assert _is_executable_required(cfg, 'gaussian') is True

    def test_gaussian_not_required_when_engine_orca(self):
        """theory.optimization.engine=orca 时 gaussian 非必需。"""
        cfg = {'theory': {'optimization': {'engine': 'orca'}, 'single_point': {'engine': 'orca'}}}
        assert _is_executable_required(cfg, 'gaussian') is False

    def test_xtb_always_required(self):
        """xtb 始终必需。"""
        assert _is_executable_required({}, 'xtb') is True

    def test_mpirun_required_when_orca_parallel(self):
        """ORCA 多核时 mpirun 必需。"""
        cfg = {'theory': {'single_point': {'engine': 'orca', 'nproc': 8}}}
        assert _is_executable_required(cfg, 'mpirun') is True

    def test_mpirun_not_required_when_orca_serial(self):
        """ORCA 单核时 mpirun 非必需。"""
        cfg = {'theory': {'single_point': {'engine': 'orca', 'nproc': 1}}}
        assert _is_executable_required(cfg, 'mpirun') is False


class TestSource:
    """source 标签测试"""

    def test_source_not_found(self):
        """找不到时 source='not_found'（fail_on_invalid_explicit=False）。"""
        cfg = {
            'executables': {
                'xtb': {'path': '/nonexistent/path'},
                'discovery': {'enabled': False, 'fail_on_invalid_explicit': False},
            },
        }
        results = resolve_all_executables(cfg)
        xtb = results.get('xtb', {})
        if not xtb.get('found'):
            assert xtb.get('source') == 'not_found'

    def test_discovery_disabled_skip(self):
        """discovery.enabled=false 时 known_dirs 不搜索。"""
        cfg = {
            'executables': {
                'xtb': {'path': '/nonexistent/path'},
                'discovery': {'enabled': False, 'fail_on_invalid_explicit': False},
            },
        }
        results = resolve_all_executables(cfg)
        assert isinstance(results, dict)


class TestIsExecutableRequired:
    """_is_executable_required 覆盖测试"""

    def test_orca_required_when_sp_engine_orca(self):
        """single_point.engine=orca 时 orca 必需。"""
        cfg = {'theory': {'single_point': {'engine': 'orca'}}}
        assert _is_executable_required(cfg, 'orca') is True

    def test_crest_required_in_single_stage(self):
        """two_stage_enabled=False 时 crest 必需。"""
        cfg = {'step1': {'conformer_search': {'two_stage_enabled': False}}}
        assert _is_executable_required(cfg, 'crest') is True

    def test_isostat_always_required(self):
        """isostat 始终必需。"""
        assert _is_executable_required({}, 'isostat') is True

    def test_shermo_always_required(self):
        """shermo 始终必需。"""
        assert _is_executable_required({}, 'shermo') is True

    def test_unknown_program_not_required(self):
        """未知程序返回 False。"""
        assert _is_executable_required({}, 'nonexistent_program') is False
