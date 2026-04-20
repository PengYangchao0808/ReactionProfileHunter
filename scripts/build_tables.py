from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from rph_core.utils.dataset_loader import load_reaction_records
from rph_core.utils.path_compat import normalize_path
from rph_core.utils.reaction_identity import normalize_record_identity


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in fieldnames})


def build_tables(*, config_path: Path, output_dir: Path) -> tuple[Path, Path]:
    cfg = _read_json(config_path) if config_path.suffix.lower() == ".json" else None
    if cfg is None:
        from rph_core.utils.config_loader import load_config

        cfg = load_config(config_path)

    run_cfg_obj = cfg.get("run", {})
    if not isinstance(run_cfg_obj, dict):
        raise ValueError("config.run must be a mapping")
    run_cfg = dict(run_cfg_obj)
    if str(run_cfg.get("source")) != "dataset":
        raise ValueError("build_tables requires run.source=dataset")

    dataset_cfg_obj = run_cfg.get("dataset", {})
    if not isinstance(dataset_cfg_obj, dict):
        raise ValueError("config.run.dataset must be a mapping")
    dataset_cfg = dict(dataset_cfg_obj)
    dataset_cfg["reaction_profiles"] = cfg.get("reaction_profiles", {}) or {}

    records = load_reaction_records(dataset_cfg)

    reaction_rows: dict[str, dict[str, Any]] = {}
    condition_rows: list[dict[str, Any]] = []
    for rec in records:
        ident = normalize_record_identity(
            {
                "rx_id": rec.rx_id,
                "precursor_smiles": rec.precursor_smiles,
                "product_smiles_main": rec.product_smiles_main,
                "reaction_type": rec.reaction_type or "",
                "cyclo_mode": (rec.raw or {}).get("cyclo_mode") or "concerted",
            }
        )

        raw = rec.raw or {}
        topology = str(raw.get("topology") or "")
        reaction_rows.setdefault(
            ident.reaction_id,
            {
                "reaction_id": ident.reaction_id,
                "reactant_smiles_canon": ident.reactant_smiles_canon,
                "product_smiles_canon": ident.product_smiles_canon,
                "reaction_type": ident.reaction_type,
                "cyclo_mode": ident.cyclo_mode,
                "topology": topology,
                "representative_row_id": ident.row_id,
            },
        )

        def _float_or_blank(value: Any) -> str:
            if value is None:
                return ""
            try:
                return f"{float(str(value).strip()):.6g}"
            except Exception:
                return ""

        t_c = raw.get("temperature_c") or raw.get("temp_celsius")
        t_k = raw.get("temperature_K") or raw.get("temperature_k")
        temperature_c = _float_or_blank(t_c)
        temperature_k = _float_or_blank(t_k)
        if not temperature_k and temperature_c:
            try:
                temperature_k = f"{float(temperature_c) + 273.15:.6g}"
            except Exception:
                temperature_k = ""

        condition_rows.append(
            {
                "condition_id": ident.condition_id,
                "reaction_id": ident.reaction_id,
                "row_id": ident.row_id,
                "temperature_c": temperature_c,
                "temperature_K": temperature_k,
                "solvent": "DCM",
                "catalyst": raw.get("catalyst") or "",
                "additive": raw.get("additive") or "",
                "oxidant": raw.get("oxidant") or "",
                "yield_pct": _float_or_blank(raw.get("yield_pct") or raw.get("yield")),
                "dr_major": _float_or_blank(raw.get("dr_major")),
                "dr_minor": _float_or_blank(raw.get("dr_minor")),
                "ee_pct": _float_or_blank(raw.get("ee_pct") or raw.get("ee")),
                "source_ref": raw.get("source_ref") or "",
            }
        )

    reaction_table = output_dir / "reaction_table.csv"
    condition_table = output_dir / "condition_table.csv"

    reaction_fieldnames = [
        "reaction_id",
        "reactant_smiles_canon",
        "product_smiles_canon",
        "reaction_type",
        "cyclo_mode",
        "topology",
        "representative_row_id",
    ]
    condition_fieldnames = [
        "condition_id",
        "reaction_id",
        "row_id",
        "temperature_c",
        "temperature_K",
        "solvent",
        "catalyst",
        "additive",
        "oxidant",
        "yield_pct",
        "dr_major",
        "dr_minor",
        "ee_pct",
        "source_ref",
    ]

    _write_csv(reaction_table, list(reaction_rows.values()), reaction_fieldnames)
    _write_csv(condition_table, condition_rows, condition_fieldnames)
    return reaction_table, condition_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Build reaction_table.csv and condition_table.csv for V2.1.1")
    parser.add_argument("--config", required=True, help="Path to config/defaults.yaml")
    parser.add_argument("--output", required=True, help="Output directory for tables")
    args = parser.parse_args()

    config_path = normalize_path(args.config)
    output_dir = normalize_path(args.output)
    reaction_table, condition_table = build_tables(config_path=config_path, output_dir=output_dir)
    print(str(reaction_table))
    print(str(condition_table))


if __name__ == "__main__":
    main()
