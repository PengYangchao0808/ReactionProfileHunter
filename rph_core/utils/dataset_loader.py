# pyright: ignore
"""
Dataset loader for pipeline task generation.

Wraps the existing TSV loader with path normalization and task limits.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple


from .path_compat import normalize_path
from .cleaner_adapter import (
    extract_product_smiles,
    extract_reactant_smiles,
    extract_formed_pairs_from_reaction_smarts,
    map_pairs_to_internal_indices,
    match_reaction_profile_key,
    parse_pairs_text,
    parse_formed_pairs_from_core_bond_changes,
)
from .tsv_dataset import ReactionRecord, load_tsv_records, TSVLoaderError


logger = logging.getLogger(__name__)


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

    try:
        records = load_tsv_records(
            path=path,
            filter_ids=list(filter_ids) if filter_ids else None,
            id_col=str(dataset_cfg.get("id_col", "rx_id")),
            precursor_smiles_col=str(dataset_cfg.get("precursor_smiles_col", "precursor_smiles")),
            ylide_leaving_group_col=str(dataset_cfg.get("ylide_leaving_group_col", "ylide_leaving_group")),
            leaving_group_col=str(dataset_cfg.get("leaving_group_col_fallback", "leaving_group")),
            product_smiles_col=str(dataset_cfg.get("product_smiles_col", "product_smiles_main")),
            ylide_smiles_col=str(dataset_cfg.get("ylide_smiles_col", "ylide_smiles")),
            delimiter=str(dataset_cfg.get("delimiter", "\t")),
        )
    except TSVLoaderError as exc:
        raise DatasetLoaderError(str(exc)) from exc

    for record in records:
        _enrich_cleaner_metadata(record, reaction_profiles)

    if max_tasks and max_tasks > 0:
        return records[:max_tasks]
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

    core_changes = _first_nonempty(raw, "core_bond_changes")
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
        mapped_smiles = extract_reactant_smiles(mapped_reaction)

    if not mapped_smiles:
        mapped_smiles = record.precursor_smiles

    if not mapped_product_smiles:
        mapped_product_smiles = _first_nonempty(raw, "product_smiles_main", "product_smiles") or record.product_smiles_main

    formed_map_pairs = parse_formed_pairs_from_core_bond_changes(core_changes)
    bond_source = "core_bond_changes"

    if not formed_map_pairs:
        formed_map_pairs = parse_pairs_text(_first_nonempty(raw, "formed_bond_map_pairs"))
        if formed_map_pairs:
            bond_source = "formed_bond_map_pairs"

    if not formed_map_pairs:
        formed_map_pairs = extract_formed_pairs_from_reaction_smarts(
            _first_nonempty(raw, "reaction_smarts", "smarts", "fallback_smarts")
        )
        if formed_map_pairs:
            bond_source = "smarts_fallback"
        else:
            bond_source = "none"

    index_source_smiles = mapped_product_smiles or mapped_smiles
    formed_index_pairs = map_pairs_to_internal_indices(index_source_smiles, formed_map_pairs)
    if formed_index_pairs:
        logger.debug(
            "Dataset row %s: formed_index_pairs=%s from SMILES=%s",
            record.rx_id,
            formed_index_pairs,
            index_source_smiles[:60] if index_source_smiles else None,
        )

    if map_status and map_status != "OK":
        logger.warning(
            "Dataset row %s uses non-OK map_status=%s; continuing with formed bond metadata",
            record.rx_id,
            map_status,
        )

    if formed_map_pairs and not formed_index_pairs:
        logger.warning(
            "Dataset row %s has formed map pairs but failed map->index conversion",
            record.rx_id,
        )

    raw["bond_change_source"] = bond_source
    if mapped_product_smiles:
        raw["mapped_product_smiles"] = mapped_product_smiles
    raw["formed_bond_map_pairs"] = _pairs_to_str(formed_map_pairs)
    raw["formed_bond_index_pairs"] = _pairs_to_str(formed_index_pairs)
    raw["forming_bonds"] = raw["formed_bond_index_pairs"]

    if formed_index_pairs:
        raw["forming_bonds_index_base"] = "0"
        raw["index_base"] = "0"

    reaction_type = _first_nonempty(raw, "reaction_type", "rxn_type")
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


def _pairs_to_str(pairs: List[Tuple[int, int]]) -> str:
    if not pairs:
        return ""
    return ";".join(f"{a}-{b}" for a, b in pairs)


def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
