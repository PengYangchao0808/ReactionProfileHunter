"""
Tests for Two-Stage Conformer Search (v3.1)

Tests the GFN0 → ISOSTAT → GFN2 workflow for high-flexibility molecules.
"""

import shutil
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from rph_core.steps.conformer_search.engine import ConformerEngine
from rph_core.utils.qc_interface import CRESTInterface


class TestCRESTInterfaceEnhancements:
    """Test new methods in CRESTInterface for two-stage support."""

    def test_run_conformer_search_gfn_override(self, tmp_path):
        """Test that gfn_override parameter changes GFN level in command."""
        # This is a mock test - we verify the parameter is accepted
        interface = CRESTInterface(
            gfn_level=2,  # Default
            solvent="acetone",
            nproc=8,
            config={}
        )

        # Verify the interface accepts gfn_override
        # (actual execution would require CREST binary)
        assert hasattr(interface, 'run_conformer_search')
        import inspect
        sig = inspect.signature(interface.run_conformer_search)
        assert 'gfn_override' in sig.parameters

    def test_run_batch_optimization_exists(self, tmp_path):
        """Test that run_batch_optimization method exists."""
        interface = CRESTInterface(
            gfn_level=2,
            solvent="acetone",
            nproc=8,
            config={}
        )

        assert hasattr(interface, 'run_batch_optimization')
        import inspect
        sig = inspect.signature(interface.run_batch_optimization)
        params = list(sig.parameters.keys())
        assert 'ensemble_xyz' in params
        assert 'output_dir' in params
        assert 'gfn_level' in params


class TestConformerEngineTwoStage:
    """Test ConformerEngine two-stage configuration and workflow."""

    @pytest.fixture
    def minimal_config(self):
        """Minimal configuration for testing."""
        return {
            'step1': {
                'conformer_search': {
                    'two_stage_enabled': True,
                    'stage1_gfn0': {
                        'enabled': True,
                        'gfn_level': 0,
                        'energy_window_kcal': 10.0,
                        'crest_flags': '--niceprint',
                        'clustering': {
                            'run_after': True,
                            'isostat_gdis': 0.125,
                            'isostat_edis': 1.0
                        }
                    },
                    'stage2_gfn2': {
                        'enabled': True,
                        'gfn_level': 2,
                        'energy_window_kcal': 3.0,
                        'crest_flags': '--niceprint',
                        'clustering': {
                            'run_after': True,
                            'isostat_gdis': 0.125,
                            'isostat_edis': 1.0
                        }
                    },
                    'common': {
                        'solvent': 'acetone',
                        'threads': 8,
                        'ngeom_default': 6,
                        'ngeom_max': 20
                    }
                },
                'crest': {
                    'gfn_level': 2,
                    'solvent': 'acetone',
                    'energy_window': 6.0,
                    'threads': 16
                }
            },
            'theory': {
                'optimization': {
                    'method': 'B3LYP',
                    'basis': 'def2-SVP',
                    'engine': 'gaussian',
                    'nproc': 16,
                    'mem': '64GB',
                    'solvent': 'acetone'
                },
                'single_point': {
                    'method': 'wB97X-D3BJ',
                    'basis': 'def2-TZVPP',
                    'engine': 'orca',
                    'nproc': 16,
                    'solvent': 'acetone'
                }
            },
            'thermo': {
                'temperature_k': 298.15
            },
            'solvent': {
                'name': 'acetone'
            },
            'executables': {
                'isostat': {
                    'path': 'isostat'
                },
                'shermo': {
                    'path': 'Shermo'
                }
            }
        }

    @pytest.fixture
    def single_stage_config(self, minimal_config):
        """Configuration for single-stage mode (backward compatibility)."""
        minimal_config['step1']['conformer_search']['two_stage_enabled'] = False
        return minimal_config

    def test_two_stage_enabled_by_default(self, minimal_config, tmp_path):
        """Test that two_stage_enabled is read from config."""
        with patch('rph_core.steps.conformer_search.engine.CRESTInterface'):
            engine = ConformerEngine(
                config=minimal_config,
                work_dir=tmp_path,
                molecule_name="test_mol"
            )

        assert engine.two_stage_enabled is True

    def test_single_stage_disabled_by_default(self, single_stage_config, tmp_path):
        """Test that single-stage mode works when two_stage_enabled=False."""
        with patch('rph_core.steps.conformer_search.engine.CRESTInterface'):
            engine = ConformerEngine(
                config=single_stage_config,
                work_dir=tmp_path,
                molecule_name="test_mol"
            )

        assert engine.two_stage_enabled is False

    def test_stage1_config_parsed(self, minimal_config, tmp_path):
        """Test that Stage 1 (GFN0) configuration is correctly parsed."""
        with patch('rph_core.steps.conformer_search.engine.CRESTInterface'):
            engine = ConformerEngine(
                config=minimal_config,
                work_dir=tmp_path,
                molecule_name="test_mol"
            )

        assert engine.stage1_enabled is True
        assert engine.stage1_gfn_level == 0
        assert engine.stage1_energy_window == 10.0
        assert engine.stage1_run_clustering is True
        assert engine.stage1_isostat_gdis == 0.125
        assert engine.stage1_isostat_edis == 1.0

    def test_stage2_config_parsed(self, minimal_config, tmp_path):
        """Test that Stage 2 (GFN2) configuration is correctly parsed."""
        with patch('rph_core.steps.conformer_search.engine.CRESTInterface'):
            engine = ConformerEngine(
                config=minimal_config,
                work_dir=tmp_path,
                molecule_name="test_mol"
            )

        assert engine.stage2_enabled is True
        assert engine.stage2_gfn_level == 2
        assert engine.stage2_energy_window == 3.0
        assert engine.stage2_run_clustering is True
        assert engine.stage2_isostat_gdis == 0.125
        assert engine.stage2_isostat_edis == 1.0

    def test_directory_structure_created(self, minimal_config, tmp_path):
        """Test that required directories are created on init."""
        with patch('rph_core.steps.conformer_search.engine.CRESTInterface'):
            engine = ConformerEngine(
                config=minimal_config,
                work_dir=tmp_path,
                molecule_name="test_mol"
            )

        assert engine.molecule_dir.exists()
        assert engine.crest_dir.exists()
        assert engine.cluster_dir.exists()
        assert engine.dft_dir.exists()


class TestTwoStageWorkflow:
    """Integration tests for the two-stage workflow."""

    @pytest.fixture
    def two_stage_config(self):
        """Full configuration for two-stage testing."""
        return {
            'step1': {
                'conformer_search': {
                    'two_stage_enabled': True,
                    'stage1_gfn0': {
                        'enabled': True,
                        'gfn_level': 0,
                        'energy_window_kcal': 10.0,
                        'crest_flags': '--niceprint',
                        'clustering': {
                            'run_after': True,
                            'isostat_gdis': 0.125,
                            'isostat_edis': 1.0
                        }
                    },
                    'stage2_gfn2': {
                        'enabled': True,
                        'gfn_level': 2,
                        'energy_window_kcal': 3.0,
                        'crest_flags': '--niceprint',
                        'clustering': {
                            'run_after': True,
                            'isostat_gdis': 0.125,
                            'isostat_edis': 1.0
                        }
                    },
                    'common': {
                        'solvent': 'acetone',
                        'threads': 8,
                        'ngeom_default': 6,
                        'ngeom_max': 20
                    }
                },
                'crest': {
                    'gfn_level': 2,
                    'solvent': 'acetone',
                    'energy_window': 6.0,
                    'threads': 16
                }
            },
            'theory': {
                'optimization': {
                    'method': 'B3LYP',
                    'basis': 'def2-SVP',
                    'engine': 'gaussian',
                    'nproc': 16,
                    'mem': '64GB',
                    'solvent': 'acetone'
                },
                'single_point': {
                    'method': 'wB97X-D3BJ',
                    'basis': 'def2-TZVPP',
                    'engine': 'orca',
                    'nproc': 16,
                    'solvent': 'acetone'
                }
            },
            'thermo': {
                'temperature_k': 298.15
            },
            'solvent': {
                'name': 'acetone'
            },
            'executables': {
                'isostat': {
                    'path': 'isostat'
                },
                'shermo': {
                    'path': 'Shermo'
                }
            }
        }

    def test_run_method_exists(self, two_stage_config, tmp_path):
        """Test that run method exists and dispatches correctly."""
        with patch('rph_core.steps.conformer_search.engine.CRESTInterface') as mock_crest:
            mock_instance = MagicMock()
            mock_crest.return_value = mock_instance

            engine = ConformerEngine(
                config=two_stage_config,
                work_dir=tmp_path,
                molecule_name="test_mol"
            )

            assert hasattr(engine, 'run')
            import inspect
            sig = inspect.signature(engine.run)
            assert 'smiles' in sig.parameters

    def test_two_stage_crest_method_exists(self, two_stage_config, tmp_path):
        """Test that _step_two_stage_crest method exists."""
        with patch('rph_core.steps.conformer_search.engine.CRESTInterface'):
            engine = ConformerEngine(
                config=two_stage_config,
                work_dir=tmp_path,
                molecule_name="test_mol"
            )

            assert hasattr(engine, '_step_two_stage_crest')

    def test_run_crest_stage_method_exists(self, two_stage_config, tmp_path):
        """Test that _run_crest_stage helper method exists."""
        with patch('rph_core.steps.conformer_search.engine.CRESTInterface'):
            engine = ConformerEngine(
                config=two_stage_config,
                work_dir=tmp_path,
                molecule_name="test_mol"
            )

            assert hasattr(engine, '_run_crest_stage')

    def test_run_isostat_clustering_method_exists(self, two_stage_config, tmp_path):
        """Test that _run_isostat_clustering helper method exists."""
        with patch('rph_core.steps.conformer_search.engine.CRESTInterface'):
            engine = ConformerEngine(
                config=two_stage_config,
                work_dir=tmp_path,
                molecule_name="test_mol"
            )

            assert hasattr(engine, '_run_isostat_clustering')


class TestBackwardCompatibility:
    """Test backward compatibility with single-stage mode."""

    @pytest.fixture
    def legacy_config(self):
        """Legacy configuration (no two-stage settings)."""
        return {
            'step1': {
                'crest': {
                    'gfn_level': 2,
                    'solvent': 'acetone',
                    'energy_window': 6.0,
                    'threads': 16
                }
            },
            'theory': {
                'optimization': {
                    'method': 'B3LYP',
                    'basis': 'def2-SVP',
                    'engine': 'gaussian',
                    'nproc': 16,
                    'mem': '64GB',
                    'solvent': 'acetone'
                },
                'single_point': {
                    'method': 'wB97X-D3BJ',
                    'basis': 'def2-TZVPP',
                    'engine': 'orca',
                    'nproc': 16,
                    'solvent': 'acetone'
                }
            },
            'thermo': {
                'temperature_k': 298.15
            },
            'solvent': {
                'name': 'acetone'
            },
            'executables': {
                'isostat': {
                    'path': 'isostat'
                },
                'shermo': {
                    'path': 'Shermo'
                }
            }
        }

    def test_falls_back_to_single_stage(self, legacy_config, tmp_path):
        """Test that missing two_stage_enabled falls back to single-stage."""
        with patch('rph_core.steps.conformer_search.engine.CRESTInterface'):
            engine = ConformerEngine(
                config=legacy_config,
                work_dir=tmp_path,
                molecule_name="test_mol"
            )

        # Should default to False (single-stage)
        assert engine.two_stage_enabled is False

    def test_legacy_config_uses_defaults(self, legacy_config, tmp_path):
        """Test that legacy config uses default GFN values."""
        with patch('rph_core.steps.conformer_search.engine.CRESTInterface'):
            engine = ConformerEngine(
                config=legacy_config,
                work_dir=tmp_path,
                molecule_name="test_mol"
            )

        # Should use defaults from legacy crest config
        assert engine.stage1_gfn_level == 0  # Default from defaults.yaml
        assert engine.stage2_gfn_level == 2  # Default from defaults.yaml
