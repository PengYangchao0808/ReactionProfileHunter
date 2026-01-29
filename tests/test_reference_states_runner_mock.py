"""
Integration tests for Reference States Runner (mock mode).

Tests directory structure creation and index generation without real QC.
"""

import pytest
import tempfile
import json
from pathlib import Path

from rph_core.reference_states import (
    ReferenceStateEntry,
    ReactionReferenceState,
    ReferenceStatesRunner
)
from rph_core.utils.tsv_dataset import ReactionRecord


class TestReferenceStatesRunnerMock:
    """Test Reference States Runner in mock mode (no QC)."""

    @pytest.fixture
    def sample_config(self):
        """Create minimal configuration."""
        return {
            'reference_states': {
                'base_dirname': 'reference_states',
                'enabled': True,
                'require_thermo': False,
                'unknown_small_molecule_policy': 'error',
                'small_molecule_map': {
                    'AcOH': {
                        'smiles': 'CC(=O)O',
                        'charge': 0,
                        'multiplicity': 1
                    }
                }
            },
            'theory': {
                'optimization': {'method': 'B3LYP', 'basis': 'def2SVP'},
                'single_point': {'method': 'WB97M-V', 'basis': 'def2-TZVPP'}
            }
        }

    @pytest.fixture
    def sample_records(self):
        """Create sample reaction records."""
        return [
            ReactionRecord(
                rx_id='39717847',
                precursor_smiles='C=CC(=O)O',
                ylide_leaving_group='AcOH',
                leaving_group='AcOH',
                product_smiles_main='O=C1C=C2CCCC(=O)[C@@H]3C[C@H]1O[C@H]23',
                raw={'solvent': 'acetone'}
            ),
            ReactionRecord(
                rx_id='39717881',
                precursor_smiles='CC(=O)O',
                ylide_leaving_group='AcOH',
                leaving_group='AcOH',
                product_smiles_main='C[C@H]1[C@H]2C(=O)CCCC3=CC(=O)[C@@H]1O[C@H]32',
                raw={'solvent': 'acetone'}
            )
        ]

    def test_directory_structure_created(self, sample_config, sample_records, tmp_path):
        """Test that reference states directory structure is created."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)
        result = runner.run(sample_records, require_thermo=False)

        ref_base = tmp_path / 'reference_states'

        assert ref_base.exists()
        assert (ref_base / 'reactions').exists()
        assert (ref_base / 'small_molecules').exists()
        assert (ref_base / 'index.json').exists()

    def test_reaction_directories_created(self, sample_config, sample_records, tmp_path):
        """Test that individual reaction directories are created."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)
        result = runner.run(sample_records, require_thermo=False)

        rx_dir = tmp_path / 'reference_states' / 'reactions'

        assert (rx_dir / 'rx_39717847').exists()
        assert (rx_dir / 'rx_39717847' / 'precursor').exists()

        assert (rx_dir / 'rx_39717881').exists()
        assert (rx_dir / 'rx_39717881' / 'precursor').exists()

    def test_small_molecule_directories_created(self, sample_config, sample_records, tmp_path):
        """Test that small molecule directories are created."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)
        result = runner.run(sample_records, require_thermo=False)

        smol_dir = tmp_path / 'reference_states' / 'small_molecules'

        assert (smol_dir / 'AcOH').exists()

    def test_energy_json_created(self, sample_config, sample_records, tmp_path):
        """Test that energy.json files are created."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)
        result = runner.run(sample_records, require_thermo=False)

        precursor_energy = (
            tmp_path / 'reference_states' / 'reactions' / 'rx_39717847' / 'precursor' / 'energy.json'
        )

        assert precursor_energy.exists()

        data = json.loads(precursor_energy.read_text())
        assert 'smiles' in data
        assert 'global_min_xyz' in data

    def test_meta_json_created(self, sample_config, sample_records, tmp_path):
        """Test that meta.json files are created."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)
        result = runner.run(sample_records, require_thermo=False)

        meta_path = (
            tmp_path / 'reference_states' / 'reactions' / 'rx_39717847' / 'meta.json'
        )

        assert meta_path.exists()

        meta_data = json.loads(meta_path.read_text())
        assert 'solvent' in meta_data
        assert meta_data['solvent'] == 'acetone'

    def test_index_json_structure(self, sample_config, sample_records, tmp_path):
        """Test that index.json has expected structure."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)
        result = runner.run(sample_records, require_thermo=False)

        index_path = tmp_path / 'reference_states' / 'index.json'
        index_data = json.loads(index_path.read_text())

        assert 'version' in index_data
        assert index_data['version'] == '1'

        assert 'csv_schema' in index_data
        assert index_data['csv_schema']['id_col'] == 'rx_id'

        assert 'reactions' in index_data
        assert '39717847' in index_data['reactions']
        assert '39717881' in index_data['reactions']

        assert 'small_molecules' in index_data
        assert 'AcOH' in index_data['small_molecules']

    def test_leaving_key_mapping(self, sample_config, sample_records, tmp_path):
        """Test that leaving small molecule keys are mapped correctly."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)
        result = runner.run(sample_records, require_thermo=False)

        index_data = json.loads(
            (tmp_path / 'reference_states' / 'index.json').read_text()
        )

        assert index_data['reactions']['39717847']['leaving_small_molecule_key'] == 'AcOH'
        assert index_data['reactions']['39717881']['leaving_small_molecule_key'] == 'AcOH'

    def test_summary_returned(self, sample_config, sample_records, tmp_path):
        """Test that run returns summary dict."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)
        result = runner.run(sample_records, require_thermo=False)

        assert result['success'] == True
        assert result['total_reactions'] == 2
        assert result['total_small_molecules'] == 1  # AcOH deduplicated
        assert 'errors' in result
        assert 'index_path' in result
