"""
Integration tests for V5.4 QC interface features (NBO-only).

V5.4 Update: Removed all NMR and Hirshfeld tests and references.
"""

import pytest
import tempfile
from pathlib import Path
from rph_core.utils.qc_interface import TaskKind, _select_task_resources, harvest_nbo_files, NBO_WHITELIST


class TestTaskKindV52:
    """Test TaskKind enum unified for V5.4."""

    def test_nbo_enum_exists(self):
        """Test that NBO enum value exists and is unique."""
        assert hasattr(TaskKind, 'NBO')
        assert TaskKind.NBO.value == 'nbo'

    def test_nbo_whitelist_exists(self):
        """Test that NBO_WHITELIST constant exists."""
        assert NBO_WHITELIST is not None
        assert '.47' in NBO_WHITELIST
        assert '.nbo' in NBO_WHITELIST
        assert '.3' in NBO_WHITELIST
        assert '.31' in NBO_WHITELIST
        assert '.41' in NBO_WHITELIST


class TestSelectTaskResourcesV52:
    """Test _select_task_resources function."""

    @pytest.fixture
    def base_config(self):
        """Create base configuration."""
        return {
            'resources': {
                'mem': '32GB',
                'nproc': 16
            },
            'theory': {
                'optimization': {
                    'method': 'B3LYP',
                    'basis': 'def2SVP'
                },
                'single_point': {
                    'method': 'WB97M-V',
                    'basis': 'def2-TZVPP'
                }
            }
        }

    def test_select_optimization_resources(self, base_config):
        """Test selecting optimization task resources."""
        result = _select_task_resources(TaskKind.OPTIMIZATION, base_config)

        assert 'mem' in result
        assert 'nproc' in result
        assert 'method' in result
        assert 'basis' in result


class TestTaskKindValuesV54:
    """Test all TaskKind enum values for V5.4 (NMR/Hirshfeld removed)."""

    def test_all_task_kinds_exist(self):
        """Test that all expected TaskKind values exist."""
        from rph_core.utils.qc_interface import TaskKind

        expected_kinds = [
            'OPTIMIZATION',
            'SINGLE_POINT',
            'FREQUENCY',
            'TS_OPTIMIZATION',
            'IRC',
            'NBO'
        ]

        for kind in expected_kinds:
            assert hasattr(TaskKind, kind), f"TaskKind.{kind} not found"

    def test_task_kind_values_unique(self):
        """Test that all TaskKind values are unique."""
        from rph_core.utils.qc_interface import TaskKind

        values = [tk.value for tk in TaskKind]
        assert len(values) == len(set(values)), "TaskKind values are not unique"

    def test_task_kind_value_mapping(self):
        """Test TaskKind value to expected string mapping."""
        from rph_core.utils.qc_interface import TaskKind

        assert TaskKind.OPTIMIZATION.value == "optimization"
        assert TaskKind.SINGLE_POINT.value == "single_point"
        assert TaskKind.FREQUENCY.value == "frequency"
        assert TaskKind.TS_OPTIMIZATION.value == "ts_optimization"
        assert TaskKind.IRC.value == "irc"
        assert TaskKind.NBO.value == "nbo"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
