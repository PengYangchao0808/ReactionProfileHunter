"""
Integration tests for V5.2 features.
"""

import pytest
import tempfile
import json
import numpy as np
from pathlib import Path

from rph_core.steps.step2_retro.retro_scanner import RetroScanner, RetroScanResultV2
from rph_core.utils.tsv_dataset import ReactionRecord
from rph_core.utils.file_io import write_xyz


class TestRetroScannerV52:
    """Test RetroScanner V5.2 features (run_with_precursor)."""

    @pytest.fixture
    def sample_config(self):
        """Create minimal configuration."""
        return {
            'step2': {
                'ts_distance': 2.2,
                'break_distance': 3.5,
                'xtb_settings': {
                    'gfn_level': 2,
                    'solvent': 'acetone',
                    'nproc': 8
                },
                'neutral_precursor': {
                    'enabled': True
                }
            },
            'resources': {
                'nproc': 16,
                'mem': '32GB'
            }
        }

    @pytest.fixture
    def sample_record(self):
        """Create a sample reaction record."""
        return ReactionRecord(
            rx_id='39717847',
            precursor_smiles='C=CC(=O)O',
            ylide_leaving_group='AcOH',
            leaving_group='AcOH',
            product_smiles_main='O=C1C=C2CCCC(=O)[C@@H]3C[C@H]1O[C@H]23',
            raw={}
        )

    @pytest.fixture
    def reactant_complex_xyz(self, tmp_path):
        """Create a sample reactant complex XYZ file."""
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        symbols = ['C', 'O']
        xyz_path = tmp_path / "reactant_complex.xyz"
        write_xyz(xyz_path, coords, symbols, title="Reactant Complex")
        return xyz_path

    def test_run_with_precursor_disabled(self, sample_config, sample_record, reactant_complex_xyz, tmp_path):
        """Test that neutral precursor is not generated when disabled."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        sample_config['step2']['neutral_precursor']['enabled'] = False
        scanner = RetroScanner(sample_config)
        result = scanner.run_with_precursor(
            reactant_complex_xyz=reactant_complex_xyz,
            record=sample_record,
            output_dir=output_dir,
            enabled=False
        )

        assert result is not None
        assert result.ts_guess_xyz is None
        assert result.reactant_xyz is None
        assert result.forming_bonds is None
        assert result.neutral_precursor_xyz is None
        assert result.meta_json_path is None

    def test_run_with_precursor_enabled_no_meta(self, sample_config, sample_record, reactant_complex_xyz, tmp_path):
        """Test that neutral precursor is generated when enabled, without meta.json."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        scanner = RetroScanner(sample_config)
        result = scanner.run_with_precursor(
            reactant_complex_xyz=reactant_complex_xyz,
            record=sample_record,
            output_dir=output_dir,
            enabled=True,
            output_meta=False
        )

        assert result is not None
        assert result.ts_guess_xyz is None
        assert result.reactant_xyz is None
        assert result.forming_bonds is None
        assert result.neutral_precursor_xyz is not None
        assert result.meta_json_path is None

        assert result.neutral_precursor_xyz.exists()
        assert result.neutral_precursor_xyz.name == "neutral_precursor.xyz"

    def test_run_with_precursor_enabled_with_meta(self, sample_config, sample_record, reactant_complex_xyz, tmp_path):
        """Test that neutral precursor and meta.json are generated when enabled."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        scanner = RetroScanner(sample_config)
        result = scanner.run_with_precursor(
            reactant_complex_xyz=reactant_complex_xyz,
            record=sample_record,
            output_dir=output_dir,
            enabled=True,
            output_meta=True
        )

        assert result is not None
        assert result.ts_guess_xyz is None
        assert result.reactant_xyz is None
        assert result.forming_bonds is None
        assert result.neutral_precursor_xyz is not None
        assert result.meta_json_path is not None

        assert result.neutral_precursor_xyz.exists()
        assert result.meta_json_path.exists()

        meta_data = json.loads(result.meta_json_path.read_text())
        assert 'precursor_smiles' in meta_data
        assert 'leaving_small_molecule_key' in meta_data
        assert 'strategy' in meta_data
        assert 'source_reactant_complex' in meta_data
        assert meta_data['precursor_smiles'] == sample_record.precursor_smiles
        assert meta_data['leaving_small_molecule_key'] == 'AcOH'
        assert meta_data['strategy'] == 'reactant_complex'

    def test_run_with_precursor_missing_source(self, sample_config, sample_record, tmp_path):
        """Test that missing reactant complex returns None."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        missing_xyz = tmp_path / "does_not_exist.xyz"

        scanner = RetroScanner(sample_config)
        result = scanner.run_with_precursor(
            reactant_complex_xyz=missing_xyz,
            record=sample_record,
            output_dir=output_dir,
            enabled=True
        )

        assert result is not None
        assert result.neutral_precursor_xyz is None
        assert result.meta_json_path is None

    def test_run_with_precursor_strategy_fallback(self, sample_config, sample_record, reactant_complex_xyz, tmp_path):
        """Test that unsupported strategy falls back to reactant_complex."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        scanner = RetroScanner(sample_config)
        result = scanner.run_with_precursor(
            reactant_complex_xyz=reactant_complex_xyz,
            record=sample_record,
            output_dir=output_dir,
            enabled=True,
            output_meta=True,
            strategy='unsupported_strategy'
        )

        assert result is not None
        assert result.neutral_precursor_xyz is not None
        assert result.neutral_precursor_xyz.exists()

        assert result.meta_json_path is not None
        meta_data = json.loads(result.meta_json_path.read_text())
        assert meta_data['strategy'] == 'reactant_complex'
