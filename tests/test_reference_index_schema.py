"""
Unit tests for reference states index schema.
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


class TestReferenceIndexSchema:
    """Test index.json structure and schema validation."""

    @pytest.fixture
    def sample_config(self):
        """Create minimal configuration."""
        return {
            'reference_states': {
                'enabled': True,
                'base_dirname': 'reference_states',
                'small_molecule_map': {
                    'AcOH': {
                        'smiles': 'CC(=O)O',
                        'charge': 0,
                        'multiplicity': 1
                    }
                }
            },
            'theory': {
                'optimization': {},
                'single_point': {}
            }
        }

    @pytest.fixture
    def sample_entry(self):
        """Create a sample reference state entry."""
        return ReferenceStateEntry(
            smiles='CC(=O)O',
            charge=0,
            multiplicity=1,
            global_min_xyz='AcOH_global_min.xyz',
            sp_energy_hartree=-229.12345678,
            g_used_hartree=-228.89876543,
            thermo={
                'temperature_k': 298.15,
                'g_sum_hartree': -228.91234567,
                'g_conc_hartree': -228.89876543,
                'g_used_hartree': -228.89876543
            },
            artifacts={
                'freq_log': 'dft/freq.log',
                'sp_out': 'dft/_SP.out',
                'shermo_sum': 'dft/_Shermo.sum'
            }
        )

    def test_reference_state_entry_to_dict(self, sample_entry):
        """Test serialization of ReferenceStateEntry."""
        d = sample_entry.to_dict()

        assert d['smiles'] == 'CC(=O)O'
        assert d['sp_energy_hartree'] == -229.12345678
        assert d['thermo'] is not None
        assert d['thermo']['g_used_hartree'] == -228.89876543

    def test_reaction_reference_state_to_dict(self):
        """Test serialization of ReactionReferenceState."""
        precursor = ReferenceStateEntry(
            smiles='C=CC(=O)O',
            global_min_xyz='precursor_global_min.xyz'
        )

        state = ReactionReferenceState(
            rx_id='39717847',
            precursor=precursor,
            leaving_small_molecule_key='AcOH',
            raw_meta={'solvent': 'acetone'}
        )

        d = state.to_dict()

        assert d['rx_id'] == '39717847'
        assert d['precursor'] is not None
        assert d['precursor']['smiles'] == 'C=CC(=O)O'
        assert d['leaving_small_molecule_key'] == 'AcOH'
        assert d['raw_meta'] is not None

    def test_index_json_structure(self, sample_config, tmp_path):
        """Test that index.json has expected structure."""
        from rph_core.utils.tsv_dataset import ReactionRecord

        records = [
            ReactionRecord(
                rx_id='39717847',
                precursor_smiles='C=CC(=O)O',
                ylide_leaving_group='AcOH',
                leaving_group='AcOH',
                raw={}
            )
        ]

        runner = ReferenceStatesRunner(sample_config, tmp_path)
        result = runner.run(records, require_thermo=False)

        index_path = tmp_path / 'reference_states' / 'index.json'
        assert index_path.exists()

        index_data = json.loads(index_path.read_text())

        assert 'version' in index_data
        assert index_data['version'] == '1'

        assert 'csv_schema' in index_data
        assert index_data['csv_schema']['id_col'] == 'rx_id'

        assert 'config_snapshot' in index_data

        assert 'reactions' in index_data
        assert '39717847' in index_data['reactions']

        assert 'small_molecules' in index_data
        assert 'AcOH' in index_data['small_molecules']

        assert 'errors' in index_data
