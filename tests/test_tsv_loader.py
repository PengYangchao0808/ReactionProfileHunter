"""
Unit tests for TSV loader module.
"""

import pytest
import tempfile
from pathlib import Path

from rph_core.utils.tsv_dataset import (
    ReactionRecord,
    load_tsv_records,
    collect_leaving_small_molecule_keys,
    TSVLoaderError
)


class TestTSVLoader:
    """Test TSV file loading and parsing."""

    @pytest.fixture
    def sample_tsv_file(self, tmp_path):
        """Create a sample TSV file for testing."""
        tsv_path = tmp_path / "test_data.tsv"
        content = """rx_id\tprecursor_smiles\tylide_leaving_group\tleaving_group\tproduct_smiles_main
39717847\tC=CC(=O)CCCC1=CC(=O)COC1OC(C)=O\tAcOH\tAcOH\tO=C1C=C2CCCC(=O)[C@@H]3C[C@H]1O[C@H]23
39717881\tCC=CC(=O)CCCC1=CC(=O)COC1OC(C)=O\tAcOH\tAcOH\tC[C@H]1[C@H]2C(=O)CCCC3=CC(=O)[C@@H]1O[C@H]32
39717883\tCC(=O)OC1OCC(=O)C=C1CCCC(=O)C=C(C)C\tAcOH\tAcOH\tCC1(C)[C@@H]2O[C@H]3C(=CC2=O)CCCC(=O)[C@H]31
"""
        tsv_path.write_text(content)
        return tsv_path

    def test_load_all_records(self, sample_tsv_file):
        """Test loading all records from TSV."""
        records = load_tsv_records(sample_tsv_file)

        assert len(records) == 3

    def test_filter_by_ids(self, sample_tsv_file):
        """Test filtering records by rx_id."""
        records = load_tsv_records(
            sample_tsv_file,
            filter_ids=["39717847", "39717881"]
        )

        assert len(records) == 2
        rx_ids = {r.rx_id for r in records}
        assert rx_ids == {"39717847", "39717881"}

    def test_required_fields(self, sample_tsv_file):
        """Test that required fields are parsed correctly."""
        records = load_tsv_records(sample_tsv_file)

        first = records[0]
        assert first.rx_id == "39717847"
        assert first.precursor_smiles == "C=CC(=O)CCCC1=CC(=O)COC1OC(C)=O"

    def test_optional_fields(self, sample_tsv_file):
        """Test that optional fields are handled correctly."""
        records = load_tsv_records(sample_tsv_file)

        first = records[0]
        assert first.ylide_leaving_group == "AcOH"
        assert first.leaving_group == "AcOH"
        assert first.product_smiles_main == "O=C1C=C2CCCC(=O)[C@@H]3C[C@H]1O[C@H]23"

    def test_file_not_found(self, tmp_path):
        """Test error when TSV file does not exist."""
        with pytest.raises(TSVLoaderError) as exc_info:
            load_tsv_records(tmp_path / "nonexistent.tsv")

        assert "not found" in str(exc_info.value)

    def test_missing_required_column(self, tmp_path):
        """Test error when required column is missing."""
        tsv_path = tmp_path / "invalid.tsv"
        tsv_path.write_text("rx_id\tsmiles\n123\tC")

        with pytest.raises(TSVLoaderError) as exc_info:
            load_tsv_records(tsv_path, precursor_smiles_col="missing_col")

        assert "Missing required columns" in str(exc_info.value)

    def test_filter_ids_no_match(self, sample_tsv_file):
        """Test error when filter_ids match nothing."""
        with pytest.raises(TSVLoaderError) as exc_info:
            load_tsv_records(sample_tsv_file, filter_ids=["99999999"])

        assert "No records found matching filter_ids" in str(exc_info.value)


class TestLeavingSmallMoleculeKeys:
    """Test extraction of leaving small molecule keys."""

    @pytest.fixture
    def sample_records(self):
        """Create sample reaction records."""
        return [
            ReactionRecord(
                rx_id="39717847",
                precursor_smiles="CC(=O)O",
                ylide_leaving_group="AcOH",
                leaving_group=None,
                raw={}
            ),
            ReactionRecord(
                rx_id="39717881",
                precursor_smiles="CC(=O)O",
                ylide_leaving_group=None,
                leaving_group="TFE",
                raw={}
            ),
            ReactionRecord(
                rx_id="39717883",
                precursor_smiles="CC(=O)O",
                ylide_leaving_group="",
                leaving_group="",
                raw={}
            ),
        ]

    def test_collect_primary_keys(self, sample_records):
        """Test that ylide_leaving_group is preferred."""
        keys = collect_leaving_small_molecule_keys(sample_records)

        assert "AcOH" in keys
        assert "TFE" in keys

    def test_fallback_to_leaving_group(self, sample_records):
        """Test fallback to leaving_group when ylide_leaving_group is empty."""
        keys = collect_leaving_small_molecule_keys(sample_records)

        assert len(keys) == 2  # AcOH + TFE

    def test_empty_keys_excluded(self, sample_records):
        """Test that empty strings are excluded."""
        keys = collect_leaving_small_molecule_keys(sample_records)

        assert "" not in keys
