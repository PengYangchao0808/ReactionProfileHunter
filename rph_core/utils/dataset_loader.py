# pyright: ignore
"""
Dataset loader for pipeline task generation.

Wraps the existing TSV loader with path normalization and task limits.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional


from .path_compat import normalize_path
from .cleaner_adapter import (
    CleanerAdapterError,
    extract_product_smiles,
    load_cleaner_records,
    match_reaction_profile_key,
)
from .tsv_dataset import ReactionRecord, load_tsv_records, TSVLoaderError


logger = logging.getLogger(__name__)

_warned_row_ids: set[str] = set()


class DatasetLoaderError(RuntimeError):
    """Raised when dataset loading fails for pipeline task generation."""


def load_reaction_records(
    dataset_cfg: Dict[str, Any],
    filter_ids: list[str] | None = None,
    max_tasks: int | None = None,
) -> list[ReactionRecord]:
    """Load reaction records from dataset config.

    Args:
        dataset_cfg: Dataset configuration dictionary.
        filter_ids: Optional list of rx_ids to filter.
        max_tasks: Optional max number of records to return (None/0 = all).

    Returns:
        List of ReactionRecord objects.
    """
    path_value = dataset_cfg.get("path")
    if not path_value:
        raise DatasetLoaderError("dataset.path is required for run.source=dataset")

    path = normalize_path(str(path_value))

    reaction_profiles = dataset_cfg.get("reaction_profiles")
    if not isinstance(reaction_profiles, dict):
        reaction_profiles = {}

    is_json = path.suffix.lower() == ".json"

    if is_json:
        records = _load_json_records(
            path=path,
            reaction_profiles=reaction_profiles,
            filter_ids=filter_ids,
        )
    else:
        try:
            records = load_tsv_records(
                path=path,
                filter_ids=list(filter_ids) if filter_ids else None,
                id_col=str(dataset_cfg.get("id_col", "rx_id")),
                precursor_smiles_col=str(dataset_cfg.get("precursor_smiles_col", "precursor_smiles")),
                small_molecular_primary_col=str(dataset_cfg.get("small_molecular_primary_col", "")),
                small_molecular_fallback_col=str(dataset_cfg.get("small_molecular_fallback_col", "")),
                product_smiles_col=str(dataset_cfg.get("product_smiles_col", "product_smiles_main")),
                ylide_smiles_col=str(dataset_cfg.get("ylide_smiles_col", "ylide_smiles")),
                delimiter=str(dataset_cfg.get("delimiter", "\t")),
            )
        except TSVLoaderError as exc:
            raise DatasetLoaderError(str(exc)) from exc

    non_ok_count = 0
    non_ok_examples: list[str] = []

    for record in records:
        _enrich_cleaner_metadata(record, reaction_profiles)
        map_status = str((record.raw or {}).get("map_status", "")).strip().upper()
        if map_status and map_status != "OK":
            non_ok_count += 1
            if len(non_ok_examples) < 3:
                non_ok_examples.append(record.rx_id)

    if non_ok_count > 0:
        logger.warning(
            "Dataset contains %d rows with non-OK map_status (examples: %s)",
            non_ok_count,
            ", ".join(non_ok_examples),
        )

    if max_tasks and max_tasks > 0:
        return records[:max_tasks]
    return records


def _load_json_records(
    path: Path,
    reaction_profiles: Dict[str, Any],
    filter_ids: list[str] | None = None,
) -> list[ReactionRecord]:
    logger.info("Loading JSON dataset from: %s", path)

    try:
        records = load_cleaner_records(
            path=path,
            reaction_profiles=reaction_profiles,
        )
    except CleanerAdapterError as exc:
        raise DatasetLoaderError(f"JSON dataset loading failed: {exc}") from exc

    if not records:
        raise DatasetLoaderError(f"No valid records found in JSON file: {path}")

    if filter_ids:
        filter_set = set(filter_ids)
        filtered = [rec for rec in records if rec.rx_id in filter_set]
        if not filtered:
            logger.warning(
                "No JSON records matched filter_ids (tried %d IDs from %d records)",
                len(filter_ids),
                len(records),
            )
        records = filtered

    logger.info("Loaded %d records from JSON dataset", len(records))
    return records


def _enrich_cleaner_metadata(record: ReactionRecord, reaction_profiles: Dict[str, Dict[str, Any]]) -> None:
    raw = dict(record.raw or {})

    map_status = str(raw.get("map_status", "")).strip().upper()
    raw["map_status"] = map_status

    map_confidence = _to_float(
        _first_nonempty(raw, "map_confidence", "mapping_confidence", "map_score", "mapping_score")
    )
    if map_confidence is not None:
        raw["map_confidence"] = f"{map_confidence:.6g}"

    mapped_product_smiles = _first_nonempty(
        raw,
        "mapped_product_smiles",
        "product_smiles_mapped",
        "mapped_product",
    )
    mapped_smiles = _first_nonempty(
        raw,
        "mapped_product_smiles",
        "product_smiles_mapped",
        "mapped_product",
    )

    mapped_reaction = _first_nonempty(raw, "rxn_smiles_mapped", "reaction_smiles_mapped")

    if not mapped_product_smiles:
        mapped_product_smiles = extract_product_smiles(mapped_reaction)

    if not mapped_smiles:
        mapped_smiles = extract_product_smiles(mapped_reaction)

    if not mapped_smiles:
        mapped_smiles = _first_nonempty(raw, "product_smiles_main", "product_smiles") or record.product_smiles_main

    if not mapped_smiles:
        mapped_smiles = record.precursor_smiles

    if not mapped_product_smiles:
        mapped_product_smiles = _first_nonempty(raw, "product_smiles_main", "product_smiles") or record.product_smiles_main
    if mapped_product_smiles:
        raw["mapped_product_smiles"] = mapped_product_smiles

    reaction_type = _first_nonempty(raw, "reaction_type", "rxn_type", "reaction_family")
    profile_key = match_reaction_profile_key(reaction_type, reaction_profiles)
    if profile_key:
        raw["reaction_profile"] = profile_key

    record.raw = raw


def _first_nonempty(row: Dict[str, Any], *keys: str) -> Optional[str]:
    for key in keys:
        value = row.get(key)
        if value is None:
            continue
        text = str(value).strip()
        if text:
            return text
    return None

def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
