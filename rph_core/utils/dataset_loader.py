# pyright: ignore
"""
Dataset loader for pipeline task generation.

Wraps the existing TSV loader with path normalization and task limits.
"""

from __future__ import annotations


from .path_compat import normalize_path
from .tsv_dataset import ReactionRecord, load_tsv_records, TSVLoaderError


class DatasetLoaderError(RuntimeError):
    """Raised when dataset loading fails for pipeline task generation."""


def load_reaction_records(
    dataset_cfg: dict,
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

    if max_tasks and max_tasks > 0:
        return records[:max_tasks]
    return records
