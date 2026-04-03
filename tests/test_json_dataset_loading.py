"""Tests for JSON dataset loading compatibility."""

import json
import tempfile
from pathlib import Path

import pytest

from rph_core.utils.dataset_loader import (
    DatasetLoaderError,
    load_reaction_records,
)
from rph_core.utils.cleaner_adapter import (
    CleanerAdapterError,
    convert_cleaner_row_to_record,
    match_reaction_profile_key,
    load_cleaner_records,
)


class TestMatchReactionProfileKey:

    def test_exact_match(self):
        profiles = {"[4+3]_default": {}, "[5+2]_default": {}}
        assert match_reaction_profile_key("[4+3]_default", profiles) == "[4+3]_default"
        assert match_reaction_profile_key("[5+2]_default", profiles) == "[5+2]_default"

    def test_bracket_with_space_suffix(self):
        profiles = {"[4+3]_default": {}, "[5+2]_default": {}, "_universal": {}}
        result = match_reaction_profile_key("[4+3] cycloaddition", profiles)
        assert result == "[4+3]_default"

    def test_bracket_with_description(self):
        profiles = {"[4+3]_default": {}, "[5+2]_default": {}, "_universal": {}}
        result = match_reaction_profile_key("[5+2] pericyclic", profiles)
        assert result == "[5+2]_default"

    def test_normalized_bracket_format(self):
        profiles = {"[4+3]_default": {}, "_universal": {}}
        result = match_reaction_profile_key("[4+3]", profiles)
        assert result == "[4+3]_default"

    def test_fallback_to_universal(self):
        profiles = {"[4+3]_default": {}, "_universal": {}}
        result = match_reaction_profile_key("[3+2] cycloaddition", profiles)
        assert result == "_universal"

    def test_no_match_no_universal(self):
        profiles = {"[4+3]_default": {}}
        result = match_reaction_profile_key("[3+2]", profiles)
        assert result is None


class TestConvertCleanerRowToRecord:
    def test_record_id_field(self):
        row = {
            "record_id": "X2003-Sch2-1A",
            "substrate_smiles": "C1=CC=CC=C1C=O",
            "product_smiles": "C1=CC=C2C=CC=CC2=C1",
        }
        record = convert_cleaner_row_to_record(row, row_index=1)
        assert record is not None
        assert record.rx_id == "X2003-Sch2-1A"

    def test_reaction_family_field(self):
        row = {
            "record_id": "test_001",
            "substrate_smiles": "C1=CC=CC=C1C=O",
            "product_smiles": "C1=CC=C2C=CC=CC2=C1",
            "reaction_family": "[4+3] cycloaddition",
        }
        reaction_profiles = {"[4+3]_default": {}, "_universal": {}}
        record = convert_cleaner_row_to_record(row, row_index=1, reaction_profiles=reaction_profiles)
        assert record is not None
        assert record.raw.get("reaction_profile") == "[4+3]_default"

    def test_substrate_smiles_mapping(self):
        row = {
            "record_id": "test_002",
            "substrate_smiles": "CC(=O)O",
            "product_smiles": "CCO",
        }
        record = convert_cleaner_row_to_record(row, row_index=1)
        assert record is not None
        assert record.precursor_smiles == "CC(=O)O"

    def test_product_smiles_mapping(self):
        row = {
            "record_id": "test_003",
            "substrate_smiles": "CC(=O)O",
            "product_smiles": "CCO",
        }
        record = convert_cleaner_row_to_record(row, row_index=1)
        assert record is not None
        assert record.product_smiles_main == "CCO"

    def test_missing_precursor_returns_none(self):
        row = {
            "record_id": "test_004",
            "product_smiles": "CCO",
        }
        record = convert_cleaner_row_to_record(row, row_index=1)
        assert record is None

    def test_id_fallback_chain(self):
        row1 = {"rx_id": "id1", "record_id": "id2", "substrate_smiles": "CC"}
        record1 = convert_cleaner_row_to_record(row1, row_index=1)
        assert record1 is not None
        assert record1.rx_id == "id1"

        row2 = {"record_id": "id2", "id": "id3", "substrate_smiles": "CC"}
        record2 = convert_cleaner_row_to_record(row2, row_index=1)
        assert record2 is not None
        assert record2.rx_id == "id2"

        row3 = {"id": "id3", "substrate_smiles": "CC"}
        record3 = convert_cleaner_row_to_record(row3, row_index=1)
        assert record3 is not None
        assert record3.rx_id == "id3"

        row4 = {"substrate_smiles": "CC"}
        record4 = convert_cleaner_row_to_record(row4, row_index=10)
        assert record4 is not None
        assert record4.rx_id == "cleaner_000010"


class TestLoadCleanerRecords:
    def test_json_list_format(self, tmp_path: Path):
        json_path = tmp_path / "test_list.json"
        data = [
            {
                "record_id": "rx_001",
                "substrate_smiles": "C1=CC=CC=C1",
                "product_smiles": "C1=CC=CC=C1O",
            },
            {
                "record_id": "rx_002",
                "substrate_smiles": "CC(=O)O",
                "product_smiles": "CCO",
            },
        ]
        json_path.write_text(json.dumps(data), encoding="utf-8")

        records = load_cleaner_records(json_path)
        assert len(records) == 2
        assert records[0].rx_id == "rx_001"
        assert records[1].rx_id == "rx_002"

    def test_json_records_format(self, tmp_path: Path):
        json_path = tmp_path / "test_records.json"
        data = {
            "benchmark_name": "TestBenchmark",
            "records": [
                {
                    "record_id": "rx_001",
                    "substrate_smiles": "C1=CC=CC=C1",
                    "product_smiles": "C1=CC=CC=C1O",
                },
            ],
        }
        json_path.write_text(json.dumps(data), encoding="utf-8")

        records = load_cleaner_records(json_path)
        assert len(records) == 1
        assert records[0].rx_id == "rx_001"

    def test_json_empty_records_raises_error(self, tmp_path: Path):
        json_path = tmp_path / "test_empty.json"
        data = {"records": []}
        json_path.write_text(json.dumps(data), encoding="utf-8")

        with pytest.raises(CleanerAdapterError):
            load_cleaner_records(json_path)

    def test_json_with_reaction_family(self, tmp_path: Path):
        json_path = tmp_path / "test_profile.json"
        data = {
            "records": [
                {
                    "record_id": "rx_001",
                    "substrate_smiles": "C1=CC=CC=C1C=O",
                    "product_smiles": "C1=CC=C2C=CC=CC2=C1",
                    "reaction_family": "[4+3] cycloaddition",
                },
            ],
        }
        json_path.write_text(json.dumps(data), encoding="utf-8")

        reaction_profiles = {"[4+3]_default": {}, "_universal": {}}
        records = load_cleaner_records(json_path, reaction_profiles=reaction_profiles)
        assert len(records) == 1
        assert records[0].raw.get("reaction_profile") == "[4+3]_default"


class TestLoadReactionRecordsJson:
    def test_json_file_detection_and_loading(self, tmp_path: Path):
        json_path = tmp_path / "dataset.json"
        data = {
            "records": [
                {
                    "record_id": "X2003-Sch2-1A",
                    "substrate_smiles": "CC(=O)O",
                    "product_smiles": "CCO",
                },
            ],
        }
        json_path.write_text(json.dumps(data), encoding="utf-8")

        dataset_cfg = {"path": str(json_path)}
        records = load_reaction_records(dataset_cfg)

        assert len(records) == 1
        assert records[0].rx_id == "X2003-Sch2-1A"
        assert records[0].precursor_smiles == "CC(=O)O"

    def test_json_with_filter_ids(self, tmp_path: Path):
        json_path = tmp_path / "dataset.json"
        data = {
            "records": [
                {"record_id": "rx_001", "substrate_smiles": "CC", "product_smiles": "CCO"},
                {"record_id": "rx_002", "substrate_smiles": "CCC", "product_smiles": "CCCO"},
                {"record_id": "rx_003", "substrate_smiles": "CCCC", "product_smiles": "CCCCO"},
            ],
        }
        json_path.write_text(json.dumps(data), encoding="utf-8")

        dataset_cfg = {"path": str(json_path)}
        records = load_reaction_records(dataset_cfg, filter_ids=["rx_001", "rx_003"])

        assert len(records) == 2
        assert {r.rx_id for r in records} == {"rx_001", "rx_003"}

    def test_json_with_max_tasks(self, tmp_path: Path):
        json_path = tmp_path / "dataset.json"
        data = {
            "records": [
                {"record_id": "rx_001", "substrate_smiles": "CC", "product_smiles": "CCO"},
                {"record_id": "rx_002", "substrate_smiles": "CCC", "product_smiles": "CCCO"},
                {"record_id": "rx_003", "substrate_smiles": "CCCC", "product_smiles": "CCCCO"},
            ],
        }
        json_path.write_text(json.dumps(data), encoding="utf-8")

        dataset_cfg = {"path": str(json_path)}
        records = load_reaction_records(dataset_cfg, max_tasks=2)

        assert len(records) == 2

    def test_json_with_reaction_profiles(self, tmp_path: Path):
        json_path = tmp_path / "dataset.json"
        data = {
            "records": [
                {
                    "record_id": "rx_001",
                    "substrate_smiles": "CC",
                    "product_smiles": "CCO",
                    "reaction_family": "[4+3] cycloaddition",
                },
            ],
        }
        json_path.write_text(json.dumps(data), encoding="utf-8")

        dataset_cfg = {
            "path": str(json_path),
            "reaction_profiles": {
                "[4+3]_default": {},
                "_universal": {},
            },
        }
        records = load_reaction_records(dataset_cfg)

        assert len(records) == 1
        assert records[0].raw.get("reaction_profile") == "[4+3]_default"

    def test_tsv_still_works(self, tmp_path: Path):
        tsv_path = tmp_path / "dataset.tsv"
        tsv_content = "rx_id\tprecursor_smiles\tproduct_smiles_main\n"
        tsv_content += "rx_001\tCC\tCCO\n"
        tsv_content += "rx_002\tCCC\tCCCO\n"
        tsv_path.write_text(tsv_content, encoding="utf-8")

        dataset_cfg = {"path": str(tsv_path)}
        records = load_reaction_records(dataset_cfg)

        assert len(records) == 2
        assert records[0].rx_id == "rx_001"
        assert records[1].rx_id == "rx_002"

    def test_empty_json_raises_error(self, tmp_path: Path):
        json_path = tmp_path / "dataset.json"
        data = {"records": []}
        json_path.write_text(json.dumps(data), encoding="utf-8")

        dataset_cfg = {"path": str(json_path)}
        with pytest.raises(DatasetLoaderError):
            load_reaction_records(dataset_cfg)


class TestXiongJsonCompatibility:
    def test_xiong_json_record_structure(self, tmp_path: Path):
        json_path = tmp_path / "xiong_2003_test.json"
        data = {
            "benchmark_name": "Xiong_2003_test",
            "document_id": "Xiong et al. 2003 JACS",
            "citation": "J. Am. Chem. Soc. 2003, 125, 12694",
            "scope_policy": "Test scope",
            "fields_note": "Test fields",
            "records": [
                {
                    "record_id": "X2003-Sch2-1A",
                    "reaction_family": "[4+3] cycloaddition",
                    "topology": "intramolecular",
                    "substrate_smiles": "CC1=CC=CC=C1C(=O)C(C)C",
                    "product_smiles": "CC1=CC=C2C(=C1)C(C)(C)C2=O",
                    "yield_percent": "80",
                    "solvent": "CH2Cl2",
                },
                {
                    "record_id": "X2003-T1-1A",
                    "reaction_family": "[4+3] cycloaddition",
                    "topology": "intramolecular",
                    "substrate_smiles": "CC1=CC=CC=C1C(=O)C(C)C",
                    "product_smiles": "CC1=CC=C2C(=C1)C(C)(C)C2=O",
                    "yield_percent": "75",
                    "solvent": "CH2Cl2",
                },
            ],
        }
        json_path.write_text(json.dumps(data), encoding="utf-8")

        dataset_cfg = {
            "path": str(json_path),
            "reaction_profiles": {
                "[4+3]_default": {},
                "_universal": {},
            },
        }
        records = load_reaction_records(dataset_cfg)

        assert len(records) == 2
        assert records[0].rx_id == "X2003-Sch2-1A"
        assert records[1].rx_id == "X2003-T1-1A"
        assert records[0].raw.get("reaction_profile") == "[4+3]_default"
        assert records[0].precursor_smiles == "CC1=CC=CC=C1C(=O)C(C)C"
        assert records[0].product_smiles_main == "CC1=CC=C2C(=C1)C(C)(C)C2=O"
        assert records[0].raw.get("topology") == "intramolecular"
        assert records[0].raw.get("solvent") == "CH2Cl2"

