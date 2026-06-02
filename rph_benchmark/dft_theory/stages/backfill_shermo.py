#!/usr/bin/env python3
"""Shermo backfill for existing GEO benchmark results.

Reads completed geo_result.json files, re-runs Shermo on the existing
freq_log outputs, and updates G-related fields.
No QC re-computation — operates on existing artifacts only.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import yaml
from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.utils.shermo_runner import run_shermo


def _read_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid JSON payload: {path}")
    return payload


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _load_yaml(path: Path) -> Dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid YAML payload: {path}")
    return payload


def _safe_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _metric_kcal(upper: Optional[float], lower: Optional[float]) -> Optional[float]:
    if upper is None or lower is None:
        return None
    return (upper - lower) * HARTREE_TO_KCAL


def run_shermo_for_point(
    shermo_bin: Path,
    freq_log_path: Optional[str],
    l2_energy: Optional[float],
    output_dir: Path,
    thermo_cfg: Dict[str, Any],
    label: str,
) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "g_hartree": None,
        "thermo_path": None,
        "thermo_error": None,
    }

    if not freq_log_path:
        result["thermo_error"] = f"missing_{label}_freq_log"
        return result

    freq_path = Path(freq_log_path).expanduser().resolve()
    if not freq_path.exists():
        result["thermo_error"] = f"{label}_freq_log_missing:{freq_path}"
        return result

    if l2_energy is None:
        result["thermo_error"] = f"missing_{label}_l2_energy"
        return result

    if not shermo_bin.exists():
        result["thermo_error"] = f"shermo_binary_missing:{shermo_bin}"
        return result

    sum_file = output_dir / f"thermo_{label}.sum"
    try:
        thermo = run_shermo(
            shermo_bin=shermo_bin,
            freq_output=freq_path,
            sp_energy=float(l2_energy),
            output_file=sum_file,
            temperature_k=thermo_cfg.get("temperature_k"),
            pressure_atm=thermo_cfg.get("pressure_atm"),
            scl_zpe=thermo_cfg.get("scl_zpe"),
            ilowfreq=thermo_cfg.get("ilowfreq"),
            imagreal=thermo_cfg.get("imagreal"),
        )
        g_value = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
        result["g_hartree"] = float(g_value)
        result["thermo_path"] = str(sum_file.resolve())
    except Exception as exc:
        result["thermo_error"] = str(exc)

    return result


def backfill_session(session_dir: Path, config: Dict[str, Any]) -> None:
    executables = config.get("executables", {}) or {}
    shermo_cfg = executables.get("shermo", {}) or {}
    shermo_path = shermo_cfg.get("path", "Shermo") if isinstance(shermo_cfg, dict) else str(shermo_cfg)
    shermo_bin = Path(shermo_path).expanduser()

    thermo_cfg = config.get("thermo", {}) or {}

    rx_ids = ["1", "3", "8", "15"]
    geo_methods = ["GEO-1", "GEO-2", "GEO-3", "GEO-4", "GEO-5"]

    updated_count = 0
    skipped_count = 0
    failed_count = 0

    for rx_id in rx_ids:
        for method_id in geo_methods:
            result_path = session_dir / f"rx{rx_id}" / "geo" / method_id / "geo_result.json"
            if not result_path.exists():
                print(f"[SKIP] rx{rx_id}/{method_id}: no geo_result.json")
                skipped_count += 1
                continue

            data = _read_json(result_path)

            # Only backfill if status is "completed" and G values are missing
            if data.get("status") != "completed":
                print(f"[SKIP] rx{rx_id}/{method_id}: status={data.get('status')}")
                skipped_count += 1
                continue

            g_complex = data.get("g_complex_hartree")
            g_product = data.get("g_product_hartree")
            g_ts = data.get("g_ts_hartree")

            if g_complex is not None and g_product is not None and g_ts is not None:
                print(f"[SKIP] rx{rx_id}/{method_id}: G values already present")
                skipped_count += 1
                continue

            print(f"[RUN ] rx{rx_id}/{method_id}: backfilling Shermo...")

            run_dir = result_path.parent
            complex_cfg = data.get("complex", {})
            product_cfg = data.get("product", {})
            ts_cfg = data.get("ts", {})

            complex_freq = complex_cfg.get("freq_log") if isinstance(complex_cfg, dict) else None
            product_freq = product_cfg.get("freq_log") if isinstance(product_cfg, dict) else None
            ts_freq = ts_cfg.get("freq_log") if isinstance(ts_cfg, dict) else None

            complex_l2 = _safe_float(complex_cfg.get("l2_energy_hartree") if isinstance(complex_cfg, dict) else None)
            product_l2 = _safe_float(product_cfg.get("l2_energy_hartree") if isinstance(product_cfg, dict) else None)
            ts_l2 = _safe_float(ts_cfg.get("l2_energy_hartree") if isinstance(ts_cfg, dict) else None)

            complex_result = run_shermo_for_point(
                shermo_bin, complex_freq, complex_l2, run_dir / "freq", thermo_cfg, "complex"
            )
            product_result = run_shermo_for_point(
                shermo_bin, product_freq, product_l2, run_dir / "freq", thermo_cfg, "product"
            )
            ts_result = run_shermo_for_point(
                shermo_bin, ts_freq, ts_l2, run_dir / "freq", thermo_cfg, "ts"
            )

            # Only overwrite G-related fields
            data["g_complex_hartree"] = complex_result["g_hartree"]
            data["g_product_hartree"] = product_result["g_hartree"]
            data["g_ts_hartree"] = ts_result["g_hartree"]

            data["dG_activation_kcal"] = _metric_kcal(data["g_ts_hartree"], data["g_complex_hartree"])
            data["dG_reaction_kcal"] = _metric_kcal(data["g_product_hartree"], data["g_complex_hartree"])

            data.setdefault("thermo_outputs", {})
            thermo_outputs = data["thermo_outputs"]
            if isinstance(thermo_outputs, dict):
                thermo_outputs["complex"] = complex_result["thermo_path"]
                thermo_outputs["product"] = product_result["thermo_path"]
                thermo_outputs["ts"] = ts_result["thermo_path"]

            data.setdefault("thermo_errors", {})
            thermo_errors = data["thermo_errors"]
            if isinstance(thermo_errors, dict):
                thermo_errors["complex"] = complex_result["thermo_error"]
                thermo_errors["product"] = product_result["thermo_error"]
                thermo_errors["ts"] = ts_result["thermo_error"]

            _write_json(result_path, data)

            all_ok = all(
                r["g_hartree"] is not None
                for r in [complex_result, product_result, ts_result]
            )
            if all_ok:
                print(f"  -> OK  dG‡={data['dG_activation_kcal']:.2f} kcal/mol  dG_rxn={data['dG_reaction_kcal']:.2f} kcal/mol")
                updated_count += 1
            else:
                errors = {
                    k: v["thermo_error"]
                    for k, v in [("complex", complex_result), ("product", product_result), ("ts", ts_result)]
                    if v["thermo_error"]
                }
                print(f"  -> PARTIAL  errors={errors}")
                updated_count += 1

    print(f"\nDone: {updated_count} updated, {skipped_count} skipped, {failed_count} failed")


def main() -> int:
    parser = argparse.ArgumentParser(description="Backfill Shermo G values for GEO benchmark results")
    parser.add_argument("--session-dir", required=True, help="Benchmark session directory")
    parser.add_argument(
        "--defaults-config",
        default=str(PROJECT_ROOT / "config" / "defaults.yaml"),
        help="Defaults YAML config",
    )
    args = parser.parse_args()

    session_dir = Path(args.session_dir).expanduser().resolve()
    config = _load_yaml(Path(args.defaults_config).expanduser().resolve())

    backfill_session(session_dir, config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
