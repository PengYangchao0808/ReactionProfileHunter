# pyright: ignore
"""
Task builder for pipeline runs.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from .dataset_loader import load_reaction_records


_RX_ID_SANITIZE_RE = re.compile(r"[^A-Za-z0-9_-]+")


@dataclass
class TaskSpec:
    rx_id: str
    product_smiles: str
    meta: dict[str, Any]


def sanitize_rx_id(rx_id: str) -> str:
    """Sanitize rx_id for filesystem paths."""
    cleaned = _RX_ID_SANITIZE_RE.sub("_", rx_id.strip())
    return cleaned or "unknown"


def build_tasks_from_run_config(run_cfg: dict[str, Any]) -> list[TaskSpec]:
    source = str(run_cfg.get("source", "single"))

    if source == "single":
        single_cfg_value = run_cfg.get("single") or {}
        single_cfg = dict(single_cfg_value) if isinstance(single_cfg_value, dict) else {}
        product_smiles_value = single_cfg.get("product_smiles")
        product_smiles = str(product_smiles_value or "").strip()
        if not product_smiles:
            raise ValueError("run.single.product_smiles is required when run.source=single")
        rx_id_value = single_cfg.get("rx_id")
        rx_id = str(rx_id_value or "manual").strip() or "manual"
        return [TaskSpec(rx_id=rx_id, product_smiles=product_smiles, meta={})]

    if source == "dataset":
        dataset_cfg_value = run_cfg.get("dataset") or {}
        dataset_cfg = dict(dataset_cfg_value) if isinstance(dataset_cfg_value, dict) else {}

        filter_ids_value = run_cfg.get("filter_ids")
        filter_ids = filter_ids_value if isinstance(filter_ids_value, list) else None

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
                    f"Missing product_smiles for rx_id={record.rx_id}. "
                    f"Check dataset.product_smiles_col."
                )
            meta = {
                "precursor_smiles": record.precursor_smiles,
                "ylide_smiles": record.ylide_smiles,
                "leaving_small_molecule_key": record.get_leaving_small_molecule_key(),
                "reaction_type": record.reaction_type,
                "reaction_profile": (record.raw or {}).get("reaction_profile"),
                "cleaner_data": record.raw or {},
            }
            tasks.append(TaskSpec(rx_id=record.rx_id, product_smiles=product_smiles, meta=meta))
        return tasks

    raise ValueError(f"Unsupported run.source: {source}")
