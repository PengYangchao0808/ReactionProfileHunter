"""
Tests for Two-Stage Conformer Search (v3.1)

Tests the GFN0 → ISOSTAT → GFN2 workflow for high-flexibility molecules.
"""

import json
import shutil
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from rph_core.steps.conformer_search.engine import ConformerEngine
from rph_core.steps.conformer_search.state_manager import ConformerStateManager
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


class TestConformerStateAndResume:
    def test_emit_structured_progress_is_parseable(self, tmp_path, caplog):
        config = {
            "step1": {"conformer_search": {"two_stage_enabled": False}, "crest": {}},
            "theory": {
                "optimization": {"nproc": 1, "mem": "1GB", "charge": 0, "multiplicity": 1},
                "single_point": {"nproc": 1},
            },
            "thermo": {"temperature_k": 298.15},
            "solvent": {"name": "acetone"},
            "executables": {"isostat": {"path": "isostat"}, "shermo": {"path": "Shermo"}},
        }
        with patch("rph_core.steps.conformer_search.engine.CRESTInterface"):
            engine = ConformerEngine(config=config, work_dir=tmp_path, molecule_name="product")

        caplog.set_level("INFO")
        engine._emit_s1_progress("unit_test", {"x": 1})

        lines = [rec.message for rec in caplog.records if rec.message.startswith("S1_PROGRESS|")]
        assert len(lines) == 1
        payload = json.loads(lines[0].split("|", 1)[1])
        assert payload["schema"] == "s1_progress_v1"
        assert payload["event"] == "unit_test"
        assert payload["molecule"] == "product"

    def test_state_manager_persists_crest_status(self, tmp_path):
        manager = ConformerStateManager(tmp_path / "product", "product")
        manager.start_run("C1=CCCCC1", two_stage_enabled=True)
        output = tmp_path / "product" / "xtb2" / "ensemble.xyz"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("3\nstate\nH 0 0 0\nH 0 0 1\nH 0 1 0\n")
        manager.mark_crest_stage("final_ensemble", "completed", output)

        reloaded = ConformerStateManager(tmp_path / "product", "product")
        assert reloaded.get_crest_output("final_ensemble") is not None
        assert reloaded.state.get("crest", {}).get("final_ensemble", {}).get("status") == "completed"

    def test_dft_loop_skips_completed_conformer(self, tmp_path):
        config = {
            "step1": {"conformer_search": {"two_stage_enabled": False}, "crest": {}},
            "theory": {
                "optimization": {"nproc": 1, "mem": "1GB", "charge": 0, "multiplicity": 1},
                "single_point": {"nproc": 1},
            },
            "thermo": {"temperature_k": 298.15},
            "solvent": {"name": "acetone"},
            "executables": {"isostat": {"path": "isostat"}, "shermo": {"path": "Shermo"}},
        }

        with patch("rph_core.steps.conformer_search.engine.CRESTInterface"):
            engine = ConformerEngine(config=config, work_dir=tmp_path, molecule_name="product")

        conf0 = engine.dft_dir / "conf_000.xyz"
        conf1 = engine.dft_dir / "conf_001.xyz"
        xyz_text = "3\nconf\nH 0 0 0\nH 0 0 1\nH 0 1 0\n"
        conf0.write_text(xyz_text)
        conf1.write_text(xyz_text)

        cached_log = engine.dft_dir / "conf_000.log"
        cached_log.write_text("cached")
        engine.state_manager.mark_conformer_completed(
            "conf_000",
            {
                "name": "conf_000",
                "g_used": -10.0,
                "g_sum": -10.0,
                "g_conc": None,
                "h_sum": -9.0,
                "u_sum": -8.0,
                "s_total": 1.0,
                "sp_energy": -10.0,
                "log_file": str(cached_log),
            },
        )

        class _FakeThermo:
            g_conc = None
            g_sum = -11.0
            h_sum = -10.0
            u_sum = -9.0
            s_total = 2.0

        class _FakeGaussian:
            def __init__(self, **kwargs):
                self.kwargs = kwargs

            def write_input_file(self, xyz_source, gjf_file, route, title):
                Path(gjf_file).write_text("# fake")

        with patch("rph_core.utils.qc_interface.GaussianInterface", _FakeGaussian), \
             patch.object(engine, "_run_gaussian_opt", return_value=(True, None)) as mock_opt, \
             patch.object(engine, "_run_orca_sp", return_value=-11.0), \
             patch(
                 "rph_core.steps.conformer_search.engine.LogParser.extract_last_converged_coords",
                 return_value=(np.array([[0.0, 0.0, 0.0]]), ["H"], None),
             ), \
             patch("rph_core.steps.conformer_search.engine.run_shermo", return_value=_FakeThermo()):
            best_log, best_sp = engine._step_dft_opt_sp_coupled([conf0, conf1])

        assert mock_opt.call_count == 1
        assert isinstance(best_log, Path)
        assert best_sp is not None
        assert engine.state_manager.is_conformer_completed("conf_001")
