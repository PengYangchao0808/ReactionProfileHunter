from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.utils.shermo_runner import ThermoResult, run_shermo

logger = logging.getLogger(__name__)


def _read_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def _safe_path_from_json(value: Any, base_dir: Path) -> Optional[Path]:
    if value is None:
        return None
    try:
        text = str(value).strip()
    except Exception:
        return None
    if not text:
        return None

    p = Path(text)
    if p.exists():
        return p
    if not p.is_absolute():
        candidate = base_dir / p
        if candidate.exists():
            return candidate
    return None


def _latest_log_file(search_dir: Path) -> Optional[Path]:
    if not search_dir.exists():
        return None
    candidates = list(search_dir.rglob("*.log")) + list(search_dir.rglob("*.out"))
    candidates = [p for p in candidates if p.is_file()]
    if not candidates:
        return None
    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


def _pick_g_h_s(thermo: ThermoResult) -> Tuple[float, float, Optional[float]]:
    g_h = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
    g_kcal = g_h * HARTREE_TO_KCAL
    h_kcal = thermo.h_sum * HARTREE_TO_KCAL
    s_cal_mol_k = thermo.s_total
    return g_kcal, h_kcal, s_cal_mol_k


@dataclass(frozen=True)
class ConditionThermoCalculator:
    config: Dict[str, Any]
    reaction_root: Path
    condition_root: Path

    def run(self) -> Path:
        output_path = self.condition_root / "thermo_calculation.json"

        temperature_k = self._load_condition_temperature_k()
        if output_path.exists():
            try:
                existing = _read_json(output_path)
                existing_t = existing.get("temperature_K")
                if existing_t is not None and abs(float(existing_t) - float(temperature_k)) < 1e-6:
                    return output_path
            except Exception:
                pass

        s3_dir = self.reaction_root / "S3_TS"
        sp_meta_path = s3_dir / "sp_matrix_metadata.json"
        resume_path = s3_dir / "s3_resume.json"

        payload: Dict[str, Any] = {
            "temperature_K": float(temperature_k),
            "source": "shermo_recalc",
            "ts": {"g_kcal": None, "h_kcal": None, "s_cal_mol_k": None},
            "reactant": {"g_kcal": None, "h_kcal": None, "s_cal_mol_k": None},
            "delta": {"dG_activation": None, "dH_activation": None, "dS_activation": None},
            "metadata": {
                "method": "",
                "solvent": "",
                "recalc_mode": "rerun_shermo",
                "s3_dir": str(s3_dir),
            },
            "artifacts": {
                "sp_matrix_metadata": str(sp_meta_path),
                "s3_resume": str(resume_path),
                "ts_freq_log": None,
                "reactant_freq_log": None,
            },
            "warnings": [],
        }

        if not sp_meta_path.exists():
            payload["warnings"].append(f"W_MISSING_SP_MATRIX_METADATA:{sp_meta_path}")
            _write_json(output_path, payload)
            return output_path

        sp_meta = _read_json(sp_meta_path)
        payload["metadata"]["method"] = str(sp_meta.get("method") or "")
        payload["metadata"]["solvent"] = str(sp_meta.get("solvent") or "")

        ts_sp_energy = sp_meta.get("e_ts")
        reactant_sp_energy = sp_meta.get("e_reactant")
        if ts_sp_energy is None or reactant_sp_energy is None:
            payload["warnings"].append("W_MISSING_SP_ENERGIES")
            _write_json(output_path, payload)
            return output_path

        ts_freq_log, reactant_freq_log = self._resolve_freq_logs(s3_dir=s3_dir, resume_path=resume_path)
        if ts_freq_log is None:
            payload["warnings"].append("W_MISSING_TS_FREQ_LOG")
        if reactant_freq_log is None:
            payload["warnings"].append("W_MISSING_REACTANT_FREQ_LOG")

        payload["artifacts"]["ts_freq_log"] = str(ts_freq_log) if ts_freq_log else None
        payload["artifacts"]["reactant_freq_log"] = str(reactant_freq_log) if reactant_freq_log else None

        shermo_bin = Path(
            (self.config.get("executables", {}) or {}).get("shermo", {}).get("path", "Shermo")
        )
        thermo_cfg = self.config.get("thermo", {}) or {}

        def _to_float(value: Any) -> Optional[float]:
            if value is None:
                return None
            try:
                return float(str(value).strip())
            except (TypeError, ValueError):
                return None

        def _to_int(value: Any) -> Optional[int]:
            if value is None:
                return None
            try:
                return int(str(value).strip())
            except (TypeError, ValueError):
                return None

        pressure_atm = _to_float(thermo_cfg.get("pressure_atm"))
        scl_zpe = _to_float(thermo_cfg.get("scl_zpe"))
        ilowfreq = _to_int(thermo_cfg.get("ilowfreq"))
        imagreal = _to_float(thermo_cfg.get("imagreal"))
        conc = thermo_cfg.get("conc")

        try:
            if ts_freq_log is not None:
                ts_out = self.condition_root / "ts_Shermo.sum"
                ts_thermo = run_shermo(
                    shermo_bin=shermo_bin,
                    freq_output=ts_freq_log,
                    sp_energy=float(ts_sp_energy),
                    output_file=ts_out,
                    temperature_k=float(temperature_k),
                    pressure_atm=pressure_atm,
                    scl_zpe=scl_zpe,
                    ilowfreq=ilowfreq,
                    imagreal=imagreal,
                    conc=str(conc) if conc is not None else None,
                )
                g_kcal, h_kcal, s_cal = _pick_g_h_s(ts_thermo)
                payload["ts"] = {"g_kcal": g_kcal, "h_kcal": h_kcal, "s_cal_mol_k": s_cal}
            if reactant_freq_log is not None:
                reactant_out = self.condition_root / "reactant_Shermo.sum"
                reactant_thermo = run_shermo(
                    shermo_bin=shermo_bin,
                    freq_output=reactant_freq_log,
                    sp_energy=float(reactant_sp_energy),
                    output_file=reactant_out,
                    temperature_k=float(temperature_k),
                    pressure_atm=pressure_atm,
                    scl_zpe=scl_zpe,
                    ilowfreq=ilowfreq,
                    imagreal=imagreal,
                    conc=str(conc) if conc is not None else None,
                )
                g_kcal, h_kcal, s_cal = _pick_g_h_s(reactant_thermo)
                payload["reactant"] = {"g_kcal": g_kcal, "h_kcal": h_kcal, "s_cal_mol_k": s_cal}
        except Exception as exc:
            payload["warnings"].append(f"W_SHERMO_RECALC_FAILED:{exc}")

        ts_g = payload["ts"].get("g_kcal")
        ts_h = payload["ts"].get("h_kcal")
        ts_s = payload["ts"].get("s_cal_mol_k")
        r_g = payload["reactant"].get("g_kcal")
        r_h = payload["reactant"].get("h_kcal")
        r_s = payload["reactant"].get("s_cal_mol_k")

        if ts_g is not None and r_g is not None:
            payload["delta"]["dG_activation"] = float(ts_g) - float(r_g)
        if ts_h is not None and r_h is not None:
            payload["delta"]["dH_activation"] = float(ts_h) - float(r_h)
        if ts_s is not None and r_s is not None:
            payload["delta"]["dS_activation"] = float(ts_s) - float(r_s)

        _write_json(output_path, payload)
        return output_path

    def _load_condition_temperature_k(self) -> float:
        manifest_path = self.condition_root / "condition_manifest.json"
        if not manifest_path.exists():
            return float((self.config.get("thermo", {}) or {}).get("temperature_k", 298.15) or 298.15)
        try:
            data = _read_json(manifest_path)
        except Exception:
            return float((self.config.get("thermo", {}) or {}).get("temperature_k", 298.15) or 298.15)

        t = data.get("temperature_K") or data.get("temperature_k")
        if t is not None:
            try:
                return float(t)
            except (TypeError, ValueError):
                pass
        t_c = data.get("temperature_c")
        if t_c is not None:
            try:
                return float(t_c) + 273.15
            except (TypeError, ValueError):
                pass
        return float((self.config.get("thermo", {}) or {}).get("temperature_k", 298.15) or 298.15)

    def _resolve_freq_logs(self, s3_dir: Path, resume_path: Path) -> Tuple[Optional[Path], Optional[Path]]:
        resume_state: Dict[str, Any] = {}
        if resume_path.exists():
            try:
                resume_state = _read_json(resume_path)
            except Exception as exc:
                logger.warning(f"Failed to read s3_resume.json: {exc}")

        ts_from_resume = _safe_path_from_json((resume_state.get("ts_opt") or {}).get("freq_log"), s3_dir)
        if ts_from_resume is None:
            ts_from_resume = _safe_path_from_json((resume_state.get("ts_opt") or {}).get("log_file"), s3_dir)

        reactant_from_resume = _safe_path_from_json((resume_state.get("reactant_opt") or {}).get("freq_log"), s3_dir)
        if reactant_from_resume is None:
            reactant_from_resume = _safe_path_from_json((resume_state.get("reactant_opt") or {}).get("log_file"), s3_dir)

        ts_log = ts_from_resume
        if ts_log is None:
            ts_log = _latest_log_file(s3_dir / "ts_opt")

        reactant_log = reactant_from_resume
        if reactant_log is None:
            reactant_log = _latest_log_file(s3_dir / "reactant_opt")

        return ts_log, reactant_log
