"""
Unit tests for Small Molecule Catalog.
"""

import pytest
from rph_core.utils.small_molecule_catalog import (
    SmallMolecule,
    SmallMoleculeCatalog,
    UnknownSmallMoleculeError
)


class TestSmallMoleculeCatalog:
    """Test small molecule catalog loading and lookup."""

    @pytest.fixture
    def sample_config(self):
        """Create sample configuration."""
        return {
            'reference_states': {
                'small_molecule_map': {
                    'AcOH': {
                        'smiles': 'CC(=O)O',
                        'charge': 0,
                        'multiplicity': 1
                    },
                    'TFE': {
                        'smiles': 'OCC(F)(F)F',
                        'charge': 0,
                        'multiplicity': 1
                    },
                    'EmptySmiles': {
                        'smiles': '',
                        'charge': 0,
                        'multiplicity': 1
                    }
                }
            }
        }

    @pytest.fixture
    def catalog(self, sample_config):
        """Initialize catalog with sample config."""
        return SmallMoleculeCatalog(sample_config)

    def test_catalog_initialization(self, catalog):
        """Test catalog loads molecules correctly."""
        assert len(catalog) == 2  # EmptySmiles should be skipped

    def test_get_existing_key(self, catalog):
        """Test getting existing molecule."""
        acoh = catalog.get('AcOH')
        assert acoh is not None
        assert acoh.key == 'AcOH'
        assert acoh.smiles == 'CC(=O)O'
        assert acoh.charge == 0
        assert acoh.multiplicity == 1

    def test_get_nonexistent_key(self, catalog):
        """Test getting nonexistent molecule returns None."""
        unknown = catalog.get('NoSuchMolecule')
        assert unknown is None

    def test_require_existing_key(self, catalog):
        """Test requiring existing molecule."""
        acoh = catalog.require('TFE')
        assert acoh.key == 'TFE'
        assert acoh.smiles == 'OCC(F)(F)F'

    def test_require_nonexistent_key_raises(self, catalog):
        """Test requiring nonexistent molecule raises error."""
        with pytest.raises(UnknownSmallMoleculeError) as exc_info:
            catalog.require('NoSuchMolecule')

        assert 'Unknown small molecule key' in str(exc_info.value)

    def test_validate_keys_all_valid(self, catalog):
        """Test validation with all valid keys."""
        unknown = catalog.validate_keys(['AcOH', 'TFE'])
        assert unknown == []

    def test_validate_keys_with_unknown(self, catalog):
        """Test validation with unknown keys."""
        unknown = catalog.validate_keys(['AcOH', 'NoSuchMolecule', 'TFE'])
        assert unknown == ['NoSuchMolecule']

    def test_list_keys(self, catalog):
        """Test listing all available keys."""
        keys = catalog.list_keys()
        assert set(keys) == {'AcOH', 'TFE'}

    def test_empty_smiles_skipped(self, catalog):
        """Test that entries with empty SMILES are skipped."""
        assert 'EmptySmiles' not in catalog.list_keys()
