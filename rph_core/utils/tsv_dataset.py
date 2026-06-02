"""
TSV Dataset Loader for Reference States

Handles loading and parsing TSV/CSV files containing reaction records
for reference state calculations (precursors, small molecular species, etc.).
"""

from pathlib import Path
from typing import Optional, List, Dict, Any
from dataclasses import dataclass
import csv


@dataclass
class ReactionRecord:
    """Single reaction record from TSV/CSV dataset."""
    rx_id: str
    precursor_smiles: str
    raw: Dict[str, str]
    small_molecular_primary: Optional[str] = None
    small_molecular_fallback: Optional[str] = None
    product_smiles_main: Optional[str] = None
    ylide_smiles: Optional[str] = None
    solvent: Optional[str] = None
    base: Optional[str] = None
    temp_celsius: Optional[str] = None
    yield_: Optional[str] = None
    ee: Optional[str] = None
    dr_major: Optional[str] = None
    dr_minor: Optional[str] = None
    reaction_type: Optional[str] = None
    forming_bonds: Optional[str] = None
    breaking_bonds: Optional[str] = None
    solvent_override: Optional[str] = None

    def __post_init__(self):
        object.__setattr__(self, 'raw', self.raw or {})

    def get_small_molecular_keys(self) -> list[str]:
        """
        Get small molecular species keys from dataset record.

        Priority:
        1) small_molecular_primary (comma/semicolon separated list)
        2) small_molecular_fallback (fallback column)

        Returns:
            List of small molecular keys (may be empty)
        """
        keys: list[str] = []
        if self.small_molecular_primary and self.small_molecular_primary.strip():
            for part in self.small_molecular_primary.replace(';', ',').split(','):
                key = part.strip()
                if key:
                    keys.append(key)
        if not keys and self.small_molecular_fallback and self.small_molecular_fallback.strip():
            for part in self.small_molecular_fallback.replace(';', ',').split(','):
                key = part.strip()
                if key:
                    keys.append(key)
        return keys


class TSVLoaderError(Exception):
    """Base exception for TSV loading errors."""
    pass


def load_tsv_records(
    path: Path,
    filter_ids: Optional[List[str]] = None,
    id_col: str = "rx_id",
    precursor_smiles_col: str = "precursor_smiles",
    small_molecular_primary_col: str = "",
    small_molecular_fallback_col: str = "",
    product_smiles_col: str = "product_smiles_main",
    ylide_smiles_col: str = "ylide_smiles",
    delimiter: str = "\t"
) -> List[ReactionRecord]:
    """
    Load reaction records from TSV/CSV file.

    Args:
        path: Path to TSV/CSV file
        filter_ids: Optional list of rx_ids to filter (None = load all)
        id_col: Column name for reaction ID (default: "rx_id")
        precursor_smiles_col: Column name for precursor SMILES
        small_molecular_primary_col: Column for small molecular keys (primary, comma-sep)
        small_molecular_fallback_col: Column for small molecular keys (fallback)
        product_smiles_col: Column name for product SMILES (optional)
        ylide_smiles_col: Column name for ylide SMILES (optional)
        delimiter: Field delimiter (default: tab)

    Returns:
        List of ReactionRecord objects

    Raises:
        TSVLoaderError: If file cannot be read or required columns missing
    """
    if not path.exists():
        raise TSVLoaderError(f"TSV file not found: {path}")

    records = []

    try:
        with open(path, 'r', encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f, delimiter=delimiter)

            fieldnames = reader.fieldnames or []
            required_cols = [id_col, precursor_smiles_col]

            missing_cols = [col for col in required_cols if col not in fieldnames]
            if missing_cols:
                raise TSVLoaderError(
                    f"Missing required columns in TSV: {missing_cols}. "
                    f"Available columns: {fieldnames}"
                )

            for row_num, row in enumerate(reader, start=2):
                rx_id = row.get(id_col, "").strip()
                if not rx_id:
                    continue

                if filter_ids and rx_id not in filter_ids:
                    continue

                precursor_smiles = row.get(precursor_smiles_col, "").strip()
                if not precursor_smiles:
                    continue

                raw_dict = {k: v for k, v in row.items() if v is not None and v.strip()}

                record = ReactionRecord(
                    rx_id=rx_id,
                    precursor_smiles=precursor_smiles,
                    small_molecular_primary=_clean_str(row.get(small_molecular_primary_col)),
                    small_molecular_fallback=_clean_str(row.get(small_molecular_fallback_col)),
                    product_smiles_main=_clean_str(row.get(product_smiles_col)),
                    ylide_smiles=_clean_str(row.get(ylide_smiles_col)),
                    solvent=_clean_str(row.get("solvent")),
                    base=_clean_str(row.get("base")),
                    temp_celsius=_clean_str(row.get("temp_celsius")),
                    yield_=_clean_str(row.get("yield")),
                    ee=_clean_str(row.get("ee")),
                    dr_major=_clean_str(row.get("dr_major")),
                    dr_minor=_clean_str(row.get("dr_minor")),
                    reaction_type=_clean_str(row.get("reaction_type")),
                    forming_bonds=_clean_str(row.get("forming_bonds")),
                    breaking_bonds=_clean_str(row.get("breaking_bonds")),
                    solvent_override=_clean_str(row.get("solvent_override")),
                    raw=raw_dict
                )

                records.append(record)

    except csv.Error as e:
        raise TSVLoaderError(f"CSV parsing error in {path}: {e}")
    except Exception as e:
        raise TSVLoaderError(f"Unexpected error reading {path}: {e}")

    if not records:
        if filter_ids:
            raise TSVLoaderError(
                f"No records found matching filter_ids: {filter_ids}"
            )
        else:
            raise TSVLoaderError(f"No valid records found in {path}")

    return records


def _clean_str(value: Any) -> Optional[str]:
    """Clean string value: strip whitespace, return None if empty."""
    if value is None:
        return None
    s = str(value).strip()
    return s if s else None


def collect_small_molecular_keys_from_records(
    records: List[ReactionRecord]
) -> set[str]:
    """
    Collect all unique small molecular keys from records.

    Args:
        records: List of ReactionRecord objects

    Returns:
        Set of unique keys (e.g., {"DMDO", "acetone", "AcOH"})
    """
    keys = set()
    for record in records:
        for key in record.get_small_molecular_keys():
            keys.add(key)
    return keys
