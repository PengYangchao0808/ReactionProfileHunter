"""
Integration tests for V5.2 thermo validation (P2-2).
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


class TestThermoValidationV52:
    """Test thermo validation upgrade to machine-detectable (P2-2)."""

    @pytest.fixture
    def sample_config(self):
        """Create minimal configuration."""
        return {
            'reference_states': {
                'base_dirname': 'reference_states',
                'enabled': True,
                'require_thermo': True,
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
    def sample_records_with_thermo(self):
        """Create sample reaction records with thermo data."""
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
                product_smiles_main='C[C@H]1[C@H]2C(=O)CCCC3=CC(=O)C4=O',
                raw={'solvent': 'acetone'}
            )
        ]

    @pytest.fixture
    def sample_records_missing_thermo(self):
        """Create sample reaction records without thermo data."""
        return [
            ReactionRecord(
                rx_id='39717847',
                precursor_smiles='C=CC(=O)O',
                ylide_leaving_group='AcOH',
                leaving_group='AcOH',
                product_smiles_main='O=C1C=C2CCCC(=O)[C@@H]3C[C@H]1O[C@H]23',
                raw={'solvent': 'acetone'}
            )
        ]

    def test_thermo_validation_passed(self, sample_config, sample_records_with_thermo, tmp_path, monkeypatch):
        """Test that thermo validation fails when require_thermo=True but no thermo (V5.2 P2-2)."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)

        result = runner.run(sample_records_with_thermo, require_thermo=True)
        assert result['success'] == True  # run() still succeeds

        index_data = json.loads((tmp_path / 'reference_states' / 'index.json').read_text())
        assert 'thermo_validation_passed' in index_data
        # In V5.1 (no QC), thermo validation should fail
        assert index_data['thermo_validation_passed'] == False
        assert 'thermo_validation_errors' in index_data
        # Check that errors were recorded
        assert len(result['errors']['missing_thermo_reactions']) > 0

    def test_thermo_validation_missing_reactions(self, sample_config, sample_records_missing_thermo, tmp_path):
        """Test that thermo validation fails when thermo fields missing (P2-2)."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)

        result = runner.run(sample_records_missing_thermo, require_thermo=True)
        assert result['success'] == True  # run() still succeeds

        index_data = json.loads((tmp_path / 'reference_states' / 'index.json').read_text())
        assert 'thermo_validation_passed' in index_data
        assert index_data['thermo_validation_passed'] == False
        assert 'thermo_validation_errors' in index_data
        assert index_data['thermo_validation_errors']['missing_reactions'] == 1
        assert '39717847' in result['errors']['missing_thermo_reactions']

    def test_thermo_validation_disabled(self, sample_config, sample_records_missing_thermo, tmp_path):
        """Test that thermo validation is skipped when require_thermo=False (P2-2)."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)

        result = runner.run(sample_records_missing_thermo, require_thermo=False)
        assert result['success'] == True

        index_data = json.loads((tmp_path / 'reference_states' / 'index.json').read_text())
        assert 'thermo_validation_passed' in index_data
        assert index_data['thermo_validation_passed'] is None
        assert 'thermo_validation_errors' not in index_data

    def test_thermo_validation_machine_detectable(self, sample_config, sample_records_missing_thermo, tmp_path):
        """Test that thermo validation is machine-detectable (P2-2)."""
        runner = ReferenceStatesRunner(sample_config, tmp_path)

        result = runner.run(sample_records_missing_thermo, require_thermo=True)

        # Machine-detectable: check for validation_passed flag
        index_path = tmp_path / 'reference_states' / 'index.json'
        index_data = json.loads(index_path.read_text())

        has_validation_failed = (
            'thermo_validation_passed' in index_data and
            index_data['thermo_validation_passed'] == False
        )
        has_error_details = (
            'thermo_validation_errors' in index_data and
            isinstance(index_data['thermo_validation_errors'], dict)
        )

        assert has_validation_failed, "Thermo validation failure should be machine-detectable"
        assert has_error_details, "Thermo validation error details should be available"
        assert index_data['thermo_validation_errors']['missing_reactions'] > 0
