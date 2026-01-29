"""
TSV Dataset Loader for Reference States

Handles loading and parsing TSV files containing reaction records
for reference state calculations (precursors, leaving small molecules, etc.).
"""

from pathlib import Path
from typing import Optional, List, Dict, Any
from dataclasses import dataclass, asdict
import csv


@dataclass
class ReactionRecord:
    """Single reaction record from TSV dataset."""
    rx_id: str
    precursor_smiles: str
    raw: Dict[str, str]
    ylide_leaving_group: Optional[str] = None
    leaving_group: Optional[str] = None
    product_smiles_main: Optional[str] = None
    ylide_smiles: Optional[str] = None
    solvent: Optional[str] = None
    base: Optional[str] = None
    temp_celsius: Optional[str] = None
    yield_: Optional[str] = None
    ee: Optional[str] = None
    dr_major: Optional[str] = None
    dr_minor: Optional[str] = None

    def __post_init__(self):
        object.__setattr__(self, 'raw', self.raw or {})

    def get_leaving_small_molecule_key(self) -> Optional[str]:
        """
        Get leaving small molecule key with priority:
        1) ylide_leaving_group (primary)
        2) leaving_group (fallback)
        """
        if self.ylide_leaving_group and self.ylide_leaving_group.strip():
            return self.ylide_leaving_group.strip()
        if self.leaving_group and self.leaving_group.strip():
            return self.leaving_group.strip()
        return None


class TSVLoaderError(Exception):
    """Base exception for TSV loading errors."""
    pass


def load_tsv_records(
    path: Path,
    filter_ids: Optional[List[str]] = None,
    id_col: str = "rx_id",
    precursor_smiles_col: str = "precursor_smiles",
    ylide_leaving_group_col: str = "ylide_leaving_group",
    leaving_group_col: str = "leaving_group",
    product_smiles_col: str = "product_smiles_main",
    ylide_smiles_col: str = "ylide_smiles",
    delimiter: str = "\t"
) -> List[ReactionRecord]:
    """
    Load reaction records from TSV file.

    Args:
        path: Path to TSV file
        filter_ids: Optional list of rx_ids to filter (None = load all)
        id_col: Column name for reaction ID (default: "rx_id")
        precursor_smiles_col: Column name for precursor SMILES
        ylide_leaving_group_col: Column name for leaving small molecule label (primary)
        leaving_group_col: Column name for leaving group (fallback)
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

            # Validate required columns
            fieldnames = reader.fieldnames or []
            required_cols = [id_col, precursor_smiles_col]

            missing_cols = [col for col in required_cols if col not in fieldnames]
            if missing_cols:
                raise TSVLoaderError(
                    f"Missing required columns in TSV: {missing_cols}. "
                    f"Available columns: {fieldnames}"
                )

            for row_num, row in enumerate(reader, start=2):  # Start at 2 (header = 1)
                # Filter by ID if specified
                rx_id = row.get(id_col, "").strip()
                if not rx_id:
                    continue  # Skip rows without rx_id

                if filter_ids and rx_id not in filter_ids:
                    continue

                # Extract precursor SMILES (required)
                precursor_smiles = row.get(precursor_smiles_col, "").strip()
                if not precursor_smiles:
                    continue  # Skip rows without precursor_smiles

                # Build raw dict for meta.json
                raw_dict = {k: v for k, v in row.items() if v is not None and v.strip()}

                record = ReactionRecord(
                    rx_id=rx_id,
                    precursor_smiles=precursor_smiles,
                    ylide_leaving_group=_clean_str(row.get(ylide_leaving_group_col)),
                    leaving_group=_clean_str(row.get(leaving_group_col)),
                    product_smiles_main=_clean_str(row.get(product_smiles_col)),
                    ylide_smiles=_clean_str(row.get(ylide_smiles_col)),
                    solvent=_clean_str(row.get("solvent")),
                    base=_clean_str(row.get("base")),
                    temp_celsius=_clean_str(row.get("temp_celsius")),
                    yield_=_clean_str(row.get("yield")),
                    ee=_clean_str(row.get("ee")),
                    dr_major=_clean_str(row.get("dr_major")),
                    dr_minor=_clean_str(row.get("dr_minor")),
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


def collect_leaving_small_molecule_keys(
    records: List[ReactionRecord]
) -> set[str]:
    """
    Collect all unique leaving small molecule keys from records.

    Args:
        records: List of ReactionRecord objects

    Returns:
        Set of unique keys (e.g., {"AcOH", "TFE"})
    """
    keys = set()
    for record in records:
        key = record.get_leaving_small_molecule_key()
        if key:
            keys.add(key)
    return keys
