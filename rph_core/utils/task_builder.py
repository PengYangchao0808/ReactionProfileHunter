from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from rph_core.utils.dataset_loader import load_reaction_records
from rph_core.utils.reaction_identity import normalize_record_identity


_RX_ID_SANITIZE_RE = re.compile(r"[^A-Za-z0-9_-]+")


@dataclass
class TaskSpec:
    rx_id: str
    row_id: str
    reaction_id: str
    reaction_cache_key: str = ""
    condition_id: str = ""
    product_smiles: str = ""
    meta: dict[str, Any] = field(default_factory=dict)


@dataclass
class BranchTaskSpec:
    parent: TaskSpec
    branch_id: str
    pathway_id: str
    product_smiles: str
    generation_policy: str
    flipped_map_numbers: list[int] = field(default_factory=list)
    fixed_stereocenters: list[int] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def sanitize_rx_id(rx_id: str) -> str:
    """Sanitize rx_id for filesystem paths."""
    cleaned = _RX_ID_SANITIZE_RE.sub("_", rx_id.strip())
    return cleaned or "unknown"


def build_tasks_from_run_config(run_cfg: dict[str, Any], theory_signature: str = "") -> list[TaskSpec]:
    source = str(run_cfg.get("source", "dataset"))

    if source == "single":
        raise RuntimeError(
            "run.source='single' is no longer supported. "
            "Use 'dataset' mode instead. "
            "See config/defaults.yaml for dataset configuration."
        )

    if source == "dataset":
        dataset_cfg_value = run_cfg.get("dataset") or {}
        dataset_cfg = dict(dataset_cfg_value) if isinstance(dataset_cfg_value, dict) else {}

        filter_ids_value = run_cfg.get("filter_ids")
        filter_ids: list[str] | None
        if isinstance(filter_ids_value, list):
            filter_ids = [s for s in (str(x).strip() for x in filter_ids_value) if s]
        else:
            filter_ids = None

        max_tasks_value = run_cfg.get("max_tasks")
        max_tasks = max_tasks_value if isinstance(max_tasks_value, int) else None

        records = load_reaction_records(
            dataset_cfg=dataset_cfg,
            filter_ids=filter_ids,
            max_tasks=max_tasks,
        )

        tasks: list[TaskSpec] = []
        for record in records:
            product_smiles = (record.product_smiles_main or "").strip()
            if not product_smiles:
                raise ValueError(
                    f"Missing product_smiles for rx_id={record.rx_id}. Check dataset.product_smiles_col."
                )
            meta = {
                "precursor_smiles": record.precursor_smiles,
                "ylide_smiles": record.ylide_smiles,
                "small_molecular_keys": record.get_small_molecular_keys(),
                "reaction_type": record.reaction_type,
                "reaction_profile": (record.raw or {}).get("reaction_profile"),
                "cleaner_data": record.raw or {},
            }

            ident = normalize_record_identity(
                {
                    "rx_id": record.rx_id,
                    "precursor_smiles": record.precursor_smiles,
                    "product_smiles_main": record.product_smiles_main,
                    "reaction_type": record.reaction_type,
                    "cyclo_mode": (record.raw or {}).get("cyclo_mode") or "concerted",
                },
                theory_signature=theory_signature,
            )

            meta.update(
                {
                    "row_id": ident.row_id,
                    "reaction_id": ident.reaction_id,
                    "reaction_cache_key": ident.reaction_cache_key,
                    "condition_id": ident.condition_id,
                    "reactant_smiles_canon": ident.reactant_smiles_canon,
                    "product_smiles_canon": ident.product_smiles_canon,
                }
            )

            tasks.append(
                TaskSpec(
                    rx_id=record.rx_id,
                    row_id=ident.row_id,
                    reaction_id=ident.reaction_id,
                    reaction_cache_key=ident.reaction_cache_key,
                    condition_id=ident.condition_id,
                    product_smiles=product_smiles,
                    meta=meta,
                )
            )
        return tasks

    raise ValueError(f"Unsupported run.source: {source}")
