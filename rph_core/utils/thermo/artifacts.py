import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional

from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.utils.thermo.options import ShermoOptions
from rph_core.utils.thermo.schema import ThermoRecord, write_thermo_json
from rph_core.utils.thermo.shermo import parse_shermo_sum

logger = logging.getLogger(__name__)


def derive_thermo_json_from_sum(
    sum_file: Path,
    output_json: Path,
    options: Optional[ShermoOptions] = None,
) -> Path:
    record = parse_shermo_sum(sum_file)
    write_thermo_json(record, output_json)
    return output_json


def ensure_thermo_json_from_entry(cache_entry: Path) -> Optional[Path]:
    thermo_path = cache_entry / "thermo.json"

    if thermo_path.exists():
        try:
            data = json.loads(thermo_path.read_text(encoding="utf-8"))
            if isinstance(data, dict) and data.get("g_kcal") is not None:
                return thermo_path
            legacy_g = data.get("G") or data.get("g")
            if legacy_g is not None:
                return thermo_path
        except Exception:
            pass

    for dft_dir_name in ("finalDFT", "dft"):
        dft_dir = cache_entry / dft_dir_name
        if not dft_dir.is_dir():
            continue
        sum_candidates = sorted(
            list(dft_dir.glob("*_Shermo.sum")) + list(dft_dir.glob("*Shermo*.sum")),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        for sum_file in sum_candidates:
            try:
                derive_thermo_json_from_sum(sum_file, thermo_path)
                logger.info(f"Derived thermo.json for {cache_entry.name} from {sum_file.name}")
                _clear_missing_shermo_warning(cache_entry)
                return thermo_path
            except Exception as exc:
                logger.debug(f"Failed to derive thermo.json from {sum_file}: {exc}")
                continue

    logger.warning(f"No usable Shermo .sum found for {cache_entry.name}; thermo.json not created")
    return None


def _clear_missing_shermo_warning(cache_entry: Path) -> None:
    meta_path = cache_entry / "cache_meta.json"
    if not meta_path.exists():
        return
    try:
        with open(meta_path, "r", encoding="utf-8") as f:
            meta = json.load(f)
        warnings_list = meta.get("warnings", [])
        if "W_MISSING_SHERMO_SUM_FOR_THERMO_JSON" in warnings_list:
            warnings_list.remove("W_MISSING_SHERMO_SUM_FOR_THERMO_JSON")
            meta["warnings"] = warnings_list
            with open(meta_path, "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2)
    except Exception:
        pass


def derive_shermo_summary(
    sum_file: Path,
    output_json: Path,
    molecule_type: str = "precursor",
    options: Optional[ShermoOptions] = None,
) -> Dict[str, Any]:
    record = parse_shermo_sum(sum_file)

    g_value_au = record.g_conc_hartree if record.g_conc_hartree is not None else record.g_sum_hartree

    summary = {
        "unit": "kcal/mol",
        "temperature_K": record.temperature_K,
        f"g_{molecule_type}": g_value_au * HARTREE_TO_KCAL,
        "g_sum": record.g_sum_hartree * HARTREE_TO_KCAL,
        "g_conc": (record.g_conc_hartree * HARTREE_TO_KCAL) if record.g_conc_hartree is not None else None,
        "h_sum": record.h_sum_hartree * HARTREE_TO_KCAL,
        "u_sum": record.u_sum_hartree * HARTREE_TO_KCAL,
        "s_total": record.s_cal_mol_k,
        "derived_from_sum": str(sum_file),
        "derived_artifacts": True,
    }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    return summary
