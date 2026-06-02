from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional


@dataclass(frozen=True)
class ThermoRecord:
    g_kcal: float
    h_kcal: float
    u_kcal: Optional[float]
    s_cal_mol_k: Optional[float]

    g_sum_hartree: float
    h_sum_hartree: float
    u_sum_hartree: float
    g_conc_hartree: Optional[float]

    temperature_K: float
    pressure_atm: Optional[float]
    concentration: Optional[str]

    source: str
    source_file: Optional[Path]


THERMO_JSON_SCHEMA_VERSION = "rph.thermo.v1"


def write_thermo_json(record: ThermoRecord, output_path: Path) -> dict[str, Any]:
    data = {
        "schema_version": THERMO_JSON_SCHEMA_VERSION,
        "unit": "kcal/mol",
        "temperature_K": record.temperature_K,
        "pressure_atm": record.pressure_atm,
        "concentration": record.concentration,
        "g_kcal": record.g_kcal,
        "h_kcal": record.h_kcal,
        "u_kcal": record.u_kcal,
        "s_cal_mol_k": record.s_cal_mol_k,
        "g_sum_hartree": record.g_sum_hartree,
        "h_sum_hartree": record.h_sum_hartree,
        "u_sum_hartree": record.u_sum_hartree,
        "g_conc_hartree": record.g_conc_hartree,
        "G": record.g_kcal,
        "g": record.g_kcal,
        "source": record.source,
        "source_file": str(record.source_file) if record.source_file else None,
        "derived_artifacts": True,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    import json
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    return data


def _to_float(value) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_thermo_json(path: Path) -> Optional[ThermoRecord]:
    import json
    if not path.exists():
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        return None
    if not isinstance(data, dict):
        return None

    g_kcal = _to_float(data.get("g_kcal"))
    h_kcal = _to_float(data.get("h_kcal"))
    u_kcal = _to_float(data.get("u_kcal"))
    s_cal_mol_k = _to_float(data.get("s_cal_mol_k"))

    if g_kcal is None:
        legacy_g = _to_float(data.get("G") or data.get("g"))
        if legacy_g is None:
            return None
        g_kcal = legacy_g
    if h_kcal is None:
        h_kcal = 0.0
    if s_cal_mol_k is None:
        s_cal_mol_k = None

    g_sum_hartree = _to_float(data.get("g_sum_hartree"))
    h_sum_hartree = _to_float(data.get("h_sum_hartree"))
    u_sum_hartree = _to_float(data.get("u_sum_hartree"))
    g_conc_hartree = _to_float(data.get("g_conc_hartree"))

    src_file = data.get("source_file")
    return ThermoRecord(
        g_kcal=g_kcal,
        h_kcal=h_kcal,
        u_kcal=u_kcal,
        s_cal_mol_k=s_cal_mol_k,
        g_sum_hartree=g_sum_hartree if g_sum_hartree is not None else 0.0,
        h_sum_hartree=h_sum_hartree if h_sum_hartree is not None else 0.0,
        u_sum_hartree=u_sum_hartree if u_sum_hartree is not None else 0.0,
        g_conc_hartree=g_conc_hartree,
        temperature_K=float(data.get("temperature_K", 298.15)),
        pressure_atm=_to_float(data.get("pressure_atm")),
        concentration=data.get("concentration"),
        source=data.get("source", "thermo_json"),
        source_file=Path(src_file) if src_file else None,
    )
