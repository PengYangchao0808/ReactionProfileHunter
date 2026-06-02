from __future__ import annotations

import csv
import logging
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from rph_core.utils import shermo_runner
from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.utils.json_io import read_json, read_json_dict, write_json
from rph_core.utils.small_molecule_cache import SmallMoleculeCache
from rph_core.utils.small_molecule_catalog import SmallMoleculeCatalog
from rph_core.utils.thermo import ensure_thermo_json_from_entry, read_thermo_json as read_thermo_record

logger = logging.getLogger(__name__)

R_KCAL_MOL_K = 1.987203e-3
R_CAL_MOL_K = 1.987203
SCF_DONE_PATTERN = re.compile(r"SCF Done:\s+E\([^)]+\)\s*=\s*([-+]?\d+\.\d+)")
CONDITION_FEATURE_MLR_COLUMNS = [
    "sample_id",
    "reaction_id",
    "condition_id",
    "branch_id",
    "thermo.temperature_K",
    "thermo.dE_activation",
    "thermo.dE_reaction",
    "thermo.dG_activation",
    "thermo.dG_reaction",
    "thermo.dG_activation_gibbs",
    "thermo.dG_reaction_gibbs",
    "thermo.dG_act_step1",
    "thermo.dG_act_full",
    "thermo.Keq_act",
    "s1_dG_act",
    "s1_Keq_act",
    "thermo.ts_g_kcal",
    "thermo.ts_h_kcal",
    "thermo.ts_s_cal_mol_k",
    "thermo.intermediate_g_kcal",
    "thermo.intermediate_h_kcal",
    "thermo.product_g_kcal",
    "thermo.precursor_g_kcal",
    "thermo.precursor_boltzmann_g_kcal",
    "s1_Nconf_eff",
    "s1_Sconf",
    "s1_ff_conformer_data_source",
    "s1_ff_Nconf_eff",
    "s1_ff_Sconf",
    "s1_ff_E_span",
    "s1_ff_E_avg_weighted",
    "s1_ff_E_std",
    "s1_ff_W_max",
    "s1_ff_Hpop",
    "s1_ff_N_total",
    "s1_dft_conformer_data_source",
    "s1_dft_Nconf_eff",
    "s1_dft_Sconf",
    "s1_dft_E_span",
    "s1_dft_E_avg_weighted",
    "s1_dft_E_std",
    "s1_dft_W_max",
    "s1_dft_Hpop",
    "s1_dft_N_total",
]


_S3_TS_TAIL = re.compile(r"(?:^|.*/)S3_TS/(.+)$")

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
    base_resolved = base_dir.resolve(strict=False)
    if p.is_absolute():
        candidates = [p]
    else:
        candidates = [base_dir / p]

    tail_match = _S3_TS_TAIL.match(text)
    if tail_match:
        candidates.append(base_dir / tail_match.group(1))

    for candidate in candidates:
        try:
            candidate.resolve(strict=False).relative_to(base_resolved)
        except ValueError:
            continue
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


def _pick_g_h_s(thermo: shermo_runner.ThermoResult) -> Tuple[float, float, Optional[float]]:
    g_h = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
    g_kcal = g_h * HARTREE_TO_KCAL
    h_kcal = thermo.h_sum * HARTREE_TO_KCAL
    s_cal_mol_k = thermo.s_total
    return g_kcal, h_kcal, s_cal_mol_k


def _as_dict(value: Any) -> Dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _to_optional_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _resolve_existing_path(value: Any, *base_dirs: Path) -> Optional[Path]:
    if value is None:
        return None
    try:
        text = str(value).strip()
    except Exception:
        return None
    if not text:
        return None

    raw_path = Path(text)
    candidates = []
    if raw_path.is_absolute():
        candidates.append(raw_path)
    else:
        for base_dir in base_dirs:
            candidates.append(base_dir / raw_path)
        candidates.append(raw_path)

    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _compute_keq(delta_g_kcal: Any, temperature_k: Any) -> Optional[float]:
    dg = _to_optional_float(delta_g_kcal)
    temperature = _to_optional_float(temperature_k)
    if dg is None or temperature is None or temperature <= 0.0:
        return None
    exponent = -dg / (R_KCAL_MOL_K * temperature)
    if exponent >= 700.0:
        return float("inf")
    if exponent <= -700.0:
        return 0.0
    return math.exp(exponent)


@dataclass(frozen=True)
class ConditionThermoEngine:
    config: Dict[str, Any]
    reaction_root: Path
    condition_root: Path
    artifact_root: Optional[Path] = None
    branch_id: Optional[str] = None

    def run(self) -> Path:
        output_path = self.condition_root / "thermo_calculation.json"
        reaction_pairs = self._build_reaction_pairs()

        temperature_k = self._load_condition_temperature_k()
        if output_path.exists():
            try:
                existing = read_json_dict(output_path)
                existing_t = existing.get("temperature_K")
                warnings = [str(w) for w in existing.get("warnings", []) or []]
                has_blocking_warning = any(
                    warning.startswith("W_MISSING_SP_MATRIX_METADATA")
                    or warning.startswith("W_MISSING_SP_ENERGIES")
                    for warning in warnings
                )
                steps = _as_dict(existing.get("steps"))
                species = _as_dict(existing.get("species"))
                delta = _as_dict(existing.get("delta"))
                conformer = _as_dict(existing.get("conformer"))
                has_new_fields = (
                    bool(steps)
                    and "precursor" in species
                    and "product" in species
                    and bool(conformer)
                    and all(
                        key in delta
                        for key in (
                            "dE_activation",
                            "dE_reaction",
                            "dG_reaction",
                            "dG_activation_gibbs",
                            "dG_reaction_gibbs",
                            "dG_act_step1",
                            "dG_act_full",
                            "Keq_act",
                        )
                    )
                )
                has_energy_data = any(
                    isinstance(step, dict) and step.get("dG") is not None
                    for step in steps.values()
                ) or any(
                    isinstance(entry, dict)
                    and any(entry.get(key) is not None for key in ("g_kcal", "h_kcal", "s_cal_mol_k"))
                    for entry in species.values()
                )
                if (
                    existing_t is not None
                    and abs(float(existing_t) - float(temperature_k)) < 1e-6
                    and has_new_fields
                    and has_energy_data
                    and not has_blocking_warning
                ):
                    _ = self.write_condition_features_mlr(self.condition_root)
                    return output_path
            except Exception:
                pass

        s3_dir = (self.artifact_root or self.reaction_root) / "S3_TS"
        sp_meta_path = s3_dir / "sp_matrix_metadata.json"
        resume_path = s3_dir / "s3_resume.json"

        payload: Dict[str, Any] = {
            "temperature_K": float(temperature_k),
            "source": "shermo_recalc",
            "ts": {"g_kcal": None, "h_kcal": None, "s_cal_mol_k": None},
            "intermediate": {"g_kcal": None, "h_kcal": None, "s_cal_mol_k": None},
            "delta": {
                "dE_activation": None,
                "dE_reaction": None,
                "dG_activation": None,
                "dH_activation": None,
                "dS_activation": None,
                "dG_reaction": None,
                "dG_activation_gibbs": None,
                "dG_reaction_gibbs": None,
                "dG_act_step1": None,
                "dG_act_full": None,
                "Keq_act": None,
            },
            "species": {
                "ts": {"g_kcal": None, "h_kcal": None, "s_cal_mol_k": None},
                "intermediate": {"g_kcal": None, "h_kcal": None, "s_cal_mol_k": None},
                "product": {"g_kcal": None, "h_kcal": None, "s_cal_mol_k": None},
                "precursor": {"g_kcal": None, "h_kcal": None, "s_cal_mol_k": None},
            },
            "conformer": {
                "g_precursor_boltzmann_avg_kcal": float("nan"),
                "conformer_N_eff": float("nan"),
                "conformer_S_conf_cal_mol_k": float("nan"),
                "conformer_top_k_energies": [],
                "n_conformers": 0,
                "source_dir": None,
                "ff_data_source": None,
                "ff_Nconf_eff": float("nan"),
                "ff_Sconf": float("nan"),
                "ff_E_span": float("nan"),
                "ff_E_avg_weighted": float("nan"),
                "ff_E_std": float("nan"),
                "ff_W_max": float("nan"),
                "ff_Hpop": float("nan"),
                "ff_N_total": 0,
                "dft_data_source": None,
                "dft_Nconf_eff": float("nan"),
                "dft_Sconf": float("nan"),
                "dft_E_span": float("nan"),
                "dft_E_avg_weighted": float("nan"),
                "dft_E_std": float("nan"),
                "dft_W_max": float("nan"),
                "dft_Hpop": float("nan"),
                "dft_N_total": 0,
            },
            "steps": {
                "oxidation": {
                    "reaction": self._format_reaction_label(*reaction_pairs["oxidation"]),
                    "dG": None,
                },
                "cycloaddition_activation": {
                    "reaction": self._format_reaction_label(*reaction_pairs["cycloaddition_activation"]),
                    "dG": None,
                },
                "overall_activation": {
                    "reaction": self._format_reaction_label(*reaction_pairs["overall_activation"]),
                    "dG": None,
                },
                "overall_reaction": {
                    "reaction": self._format_reaction_label(*reaction_pairs["overall_reaction"]),
                    "dG": None,
                },
            },
            "metadata": {
                "method": "",
                "solvent": "",
                "recalc_mode": "rerun_shermo",
                "s3_dir": str(s3_dir),
                "branch_id": self._resolve_branch_id(self.condition_root),
            },
            "artifacts": {
                "sp_matrix_metadata": str(sp_meta_path),
                "s3_resume": str(resume_path),
                "ts_freq_log": None,
                "intermediate_freq_log": None,
                "product_freq_log": None,
                "precursor_freq_log": None,
            },
            "warnings": [],
        }

        for key in self._reference_small_molecule_keys():
            if key not in payload["species"]:
                payload["species"][key] = {"g_kcal": None, "h_kcal": None, "s_cal_mol_k": None}

        if not sp_meta_path.exists():
            payload["warnings"].append(f"W_MISSING_SP_MATRIX_METADATA:{sp_meta_path}")
            write_json(output_path, payload)
            _ = self.write_condition_features_mlr(self.condition_root)
            return output_path

        sp_meta = read_json_dict(sp_meta_path)
        payload["metadata"]["method"] = str(sp_meta.get("method") or "")
        payload["metadata"]["solvent"] = str(sp_meta.get("solvent") or "")

        ts_sp_energy = _to_optional_float(sp_meta.get("e_ts"))
        intermediate_sp_energy = _to_optional_float(sp_meta.get("e_reactant"))
        product_sp_energy = _to_optional_float(sp_meta.get("e_product"))

        if ts_sp_energy is not None and intermediate_sp_energy is not None:
            payload["delta"]["dE_activation"] = (float(ts_sp_energy) - float(intermediate_sp_energy)) * HARTREE_TO_KCAL

        if ts_sp_energy is None or intermediate_sp_energy is None:
            payload["warnings"].append("W_MISSING_SP_ENERGIES")
            write_json(output_path, payload)
            _ = self.write_condition_features_mlr(self.condition_root)
            return output_path

        ts_freq_log, intermediate_freq_log = self._resolve_freq_logs(s3_dir=s3_dir, resume_path=resume_path)
        if ts_freq_log is None:
            payload["warnings"].append("W_MISSING_TS_FREQ_LOG")
        if intermediate_freq_log is None:
            payload["warnings"].append("W_MISSING_INTERMEDIATE_FREQ_LOG")

        payload["artifacts"]["ts_freq_log"] = str(ts_freq_log) if ts_freq_log else None
        payload["artifacts"]["intermediate_freq_log"] = str(intermediate_freq_log) if intermediate_freq_log else None

        shermo_bin = Path(
            (self.config.get("executables", {}) or {}).get("shermo", {}).get("path", "Shermo")
        )
        thermo_cfg = self.config.get("thermo", {}) or {}

        def _to_int(value: Any) -> Optional[int]:
            if value is None:
                return None
            try:
                return int(str(value).strip())
            except (TypeError, ValueError):
                return None

        pressure_atm = _to_optional_float(thermo_cfg.get("pressure_atm"))
        scl_zpe = _to_optional_float(thermo_cfg.get("scl_zpe"))
        ilowfreq = _to_int(thermo_cfg.get("ilowfreq"))
        imagreal = _to_optional_float(thermo_cfg.get("imagreal"))
        conc = thermo_cfg.get("conc")

        try:
            if ts_freq_log is not None:
                ts_thermo = self._run_shermo_for_log(
                    shermo_bin=shermo_bin,
                    freq_log=ts_freq_log,
                    sp_energy=float(ts_sp_energy),
                    output_file=self.condition_root / "ts_Shermo.sum",
                    temperature_k=float(temperature_k),
                    pressure_atm=pressure_atm,
                    scl_zpe=scl_zpe,
                    ilowfreq=ilowfreq,
                    imagreal=imagreal,
                    conc=conc,
                )
                if ts_thermo is not None:
                    g_kcal, h_kcal, s_cal = _pick_g_h_s(ts_thermo)
                    payload["ts"] = {"g_kcal": g_kcal, "h_kcal": h_kcal, "s_cal_mol_k": s_cal}
            if intermediate_freq_log is not None:
                intermediate_thermo = self._run_shermo_for_log(
                    shermo_bin=shermo_bin,
                    freq_log=intermediate_freq_log,
                    sp_energy=float(intermediate_sp_energy),
                    output_file=self.condition_root / "intermediate_Shermo.sum",
                    temperature_k=float(temperature_k),
                    pressure_atm=pressure_atm,
                    scl_zpe=scl_zpe,
                    ilowfreq=ilowfreq,
                    imagreal=imagreal,
                    conc=conc,
                )
                if intermediate_thermo is not None:
                    g_kcal, h_kcal, s_cal = _pick_g_h_s(intermediate_thermo)
                    payload["intermediate"] = {"g_kcal": g_kcal, "h_kcal": h_kcal, "s_cal_mol_k": s_cal}
        except Exception as exc:
            payload["warnings"].append(f"W_SHERMO_RECALC_FAILED:{exc}")

        conformer_payload = self._compute_conformer_boltzmann(float(temperature_k))
        payload["conformer"] = conformer_payload

        product_reference = self._resolve_s1_species_reference("product")
        precursor_reference = self._resolve_s1_species_reference("precursor")

        if product_sp_energy is None:
            product_sp_energy = _to_optional_float(product_reference.get("sp_energy"))
        precursor_sp_energy = _to_optional_float(precursor_reference.get("sp_energy"))

        product_freq_log = product_reference.get("freq_log")
        precursor_freq_log = precursor_reference.get("freq_log")
        payload["artifacts"]["product_freq_log"] = str(product_freq_log) if isinstance(product_freq_log, Path) else None
        payload["artifacts"]["precursor_freq_log"] = str(precursor_freq_log) if isinstance(precursor_freq_log, Path) else None

        if product_sp_energy is not None and intermediate_sp_energy is not None:
            payload["delta"]["dE_reaction"] = (float(product_sp_energy) - float(intermediate_sp_energy)) * HARTREE_TO_KCAL

        if product_freq_log is None:
            payload["warnings"].append("W_MISSING_PRODUCT_FREQ_LOG")
        elif product_sp_energy is None:
            payload["warnings"].append("W_MISSING_PRODUCT_SP_ENERGY")
        else:
            try:
                product_thermo = self._run_shermo_for_log(
                    shermo_bin=shermo_bin,
                    freq_log=product_freq_log,
                    sp_energy=float(product_sp_energy),
                    output_file=self.condition_root / "product_Shermo.sum",
                    temperature_k=float(temperature_k),
                    pressure_atm=pressure_atm,
                    scl_zpe=scl_zpe,
                    ilowfreq=ilowfreq,
                    imagreal=imagreal,
                    conc=conc,
                )
                if product_thermo is not None:
                    g_kcal, h_kcal, s_cal = _pick_g_h_s(product_thermo)
                    payload["species"]["product"] = {"g_kcal": g_kcal, "h_kcal": h_kcal, "s_cal_mol_k": s_cal}
            except Exception as exc:
                payload["warnings"].append(f"W_PRODUCT_SHERMO_RECALC_FAILED:{exc}")

        if precursor_freq_log is None:
            payload["warnings"].append("W_MISSING_PRECURSOR_FREQ_LOG")
        elif precursor_sp_energy is None:
            payload["warnings"].append("W_MISSING_PRECURSOR_SP_ENERGY")
        else:
            try:
                precursor_thermo = self._run_shermo_for_log(
                    shermo_bin=shermo_bin,
                    freq_log=precursor_freq_log,
                    sp_energy=float(precursor_sp_energy),
                    output_file=self.condition_root / "precursor_Shermo.sum",
                    temperature_k=float(temperature_k),
                    pressure_atm=pressure_atm,
                    scl_zpe=scl_zpe,
                    ilowfreq=ilowfreq,
                    imagreal=imagreal,
                    conc=conc,
                )
                if precursor_thermo is not None:
                    g_kcal, h_kcal, s_cal = _pick_g_h_s(precursor_thermo)
                    payload["species"]["precursor"] = {"g_kcal": g_kcal, "h_kcal": h_kcal, "s_cal_mol_k": s_cal}
            except Exception as exc:
                payload["warnings"].append(f"W_PRECURSOR_SHERMO_RECALC_FAILED:{exc}")

        ts_g = _to_optional_float(payload["ts"].get("g_kcal"))
        ts_h = _to_optional_float(payload["ts"].get("h_kcal"))
        ts_s = _to_optional_float(payload["ts"].get("s_cal_mol_k"))
        intermediate_g = _to_optional_float(payload["intermediate"].get("g_kcal"))
        intermediate_h = _to_optional_float(payload["intermediate"].get("h_kcal"))
        intermediate_s = _to_optional_float(payload["intermediate"].get("s_cal_mol_k"))
        product_g = _to_optional_float(_as_dict(payload["species"].get("product")).get("g_kcal"))
        precursor_g = _to_optional_float(_as_dict(payload["species"].get("precursor")).get("g_kcal"))
        precursor_boltzmann_g = _to_optional_float(conformer_payload.get("g_precursor_boltzmann_avg_kcal"))

        if ts_g is not None and intermediate_g is not None:
            activation_g = float(ts_g) - float(intermediate_g)
            payload["delta"]["dG_activation"] = activation_g
            payload["delta"]["dG_activation_gibbs"] = activation_g
            payload["delta"]["Keq_act"] = _compute_keq(activation_g, temperature_k)
        if ts_h is not None and intermediate_h is not None:
            payload["delta"]["dH_activation"] = float(ts_h) - float(intermediate_h)
        if ts_s is not None and intermediate_s is not None:
            payload["delta"]["dS_activation"] = float(ts_s) - float(intermediate_s)
        if product_g is not None and intermediate_g is not None:
            reaction_g = float(product_g) - float(intermediate_g)
            payload["delta"]["dG_reaction"] = reaction_g
            payload["delta"]["dG_reaction_gibbs"] = reaction_g
        if precursor_boltzmann_g is not None and precursor_g is not None:
            payload["delta"]["dG_act_step1"] = float(precursor_boltzmann_g) - float(precursor_g)

        species: Dict[str, Dict[str, Optional[float]]] = {
            "ts": dict(payload["ts"]),
            "intermediate": dict(payload["intermediate"]),
            "product": dict(_as_dict(payload["species"].get("product"))),
            "precursor": dict(_as_dict(payload["species"].get("precursor"))),
        }
        for key in self._reference_small_molecule_keys():
            if key not in species:
                species[key] = {"g_kcal": None, "h_kcal": None, "s_cal_mol_k": None}
        self._add_reference_species(species)
        payload["species"] = species

        if precursor_boltzmann_g is not None:
            species_for_step1 = {key: dict(value) for key, value in payload["species"].items()}
            species_for_step1["precursor"] = {
                "g_kcal": float(precursor_boltzmann_g),
                "h_kcal": _to_optional_float(species_for_step1.get("precursor", {}).get("h_kcal")),
                "s_cal_mol_k": _to_optional_float(species_for_step1.get("precursor", {}).get("s_cal_mol_k")),
            }
            payload["delta"]["dG_act_full"] = self._compute_delta_g(
                *reaction_pairs["oxidation"],
                species_map=species_for_step1,
            )

        payload["steps"] = {
            "oxidation": {
                "reaction": self._format_reaction_label(*reaction_pairs["oxidation"]),
                "dG": self._compute_delta_g(*reaction_pairs["oxidation"], species_map=payload["species"]),
            },
            "cycloaddition_activation": {
                "reaction": self._format_reaction_label(*reaction_pairs["cycloaddition_activation"]),
                "dG": self._compute_delta_g(
                    *reaction_pairs["cycloaddition_activation"],
                    species_map=payload["species"],
                ),
            },
            "overall_activation": {
                "reaction": self._format_reaction_label(*reaction_pairs["overall_activation"]),
                "dG": self._compute_delta_g(*reaction_pairs["overall_activation"], species_map=payload["species"]),
            },
            "overall_reaction": {
                "reaction": self._format_reaction_label(*reaction_pairs["overall_reaction"]),
                "dG": self._compute_delta_g(*reaction_pairs["overall_reaction"], species_map=payload["species"]),
            },
        }

        write_json(output_path, payload)
        _ = self.write_condition_features_mlr(self.condition_root)
        return output_path

    def write_condition_features_mlr(self, output_dir: Path) -> Path:
        output_dir = Path(output_dir)
        thermo_path = output_dir / "thermo_calculation.json"
        if not thermo_path.exists():
            _ = self.run()

        thermo = read_json_dict(thermo_path)
        manifest = read_json_dict(output_dir / "condition_manifest.json")
        delta = _as_dict(thermo.get("delta"))
        species = _as_dict(thermo.get("species"))
        conformer = _as_dict(thermo.get("conformer"))

        ts = _as_dict(thermo.get("ts")) or _as_dict(species.get("ts"))
        intermediate = _as_dict(thermo.get("intermediate")) or _as_dict(species.get("intermediate"))
        product = _as_dict(species.get("product"))
        precursor = _as_dict(species.get("precursor"))

        reaction_id = manifest.get("reaction_id") or self.reaction_root.name
        condition_id = manifest.get("condition_id") or output_dir.name
        branch_id = self._resolve_branch_id(output_dir)
        s1_dg_act = delta.get("dG_act_full")

        row = {
            "sample_id": self._build_sample_id(manifest=manifest, output_dir=output_dir),
            "reaction_id": reaction_id,
            "condition_id": condition_id,
            "branch_id": branch_id,
            "thermo.temperature_K": thermo.get("temperature_K"),
            "thermo.dE_activation": delta.get("dE_activation"),
            "thermo.dE_reaction": delta.get("dE_reaction"),
            "thermo.dG_activation": delta.get("dG_activation"),
            "thermo.dG_reaction": delta.get("dG_reaction"),
            "thermo.dG_activation_gibbs": delta.get("dG_activation_gibbs"),
            "thermo.dG_reaction_gibbs": delta.get("dG_reaction_gibbs"),
            "thermo.dG_act_step1": delta.get("dG_act_step1"),
            "thermo.dG_act_full": delta.get("dG_act_full"),
            "thermo.Keq_act": delta.get("Keq_act"),
            "s1_dG_act": s1_dg_act,
            "s1_Keq_act": _compute_keq(s1_dg_act, thermo.get("temperature_K")),
            "thermo.ts_g_kcal": ts.get("g_kcal"),
            "thermo.ts_h_kcal": ts.get("h_kcal"),
            "thermo.ts_s_cal_mol_k": ts.get("s_cal_mol_k"),
            "thermo.intermediate_g_kcal": intermediate.get("g_kcal"),
            "thermo.intermediate_h_kcal": intermediate.get("h_kcal"),
            "thermo.product_g_kcal": product.get("g_kcal"),
            "thermo.precursor_g_kcal": precursor.get("g_kcal"),
            "thermo.precursor_boltzmann_g_kcal": conformer.get("g_precursor_boltzmann_avg_kcal"),
            "s1_Nconf_eff": conformer.get("conformer_N_eff"),
            "s1_Sconf": conformer.get("conformer_S_conf_cal_mol_k"),
            "s1_ff_conformer_data_source": conformer.get("ff_data_source"),
            "s1_ff_Nconf_eff": conformer.get("ff_Nconf_eff"),
            "s1_ff_Sconf": conformer.get("ff_Sconf"),
            "s1_ff_E_span": conformer.get("ff_E_span"),
            "s1_ff_E_avg_weighted": conformer.get("ff_E_avg_weighted"),
            "s1_ff_E_std": conformer.get("ff_E_std"),
            "s1_ff_W_max": conformer.get("ff_W_max"),
            "s1_ff_Hpop": conformer.get("ff_Hpop"),
            "s1_ff_N_total": conformer.get("ff_N_total"),
            "s1_dft_conformer_data_source": conformer.get("dft_data_source"),
            "s1_dft_Nconf_eff": conformer.get("dft_Nconf_eff"),
            "s1_dft_Sconf": conformer.get("dft_Sconf"),
            "s1_dft_E_span": conformer.get("dft_E_span"),
            "s1_dft_E_avg_weighted": conformer.get("dft_E_avg_weighted"),
            "s1_dft_E_std": conformer.get("dft_E_std"),
            "s1_dft_W_max": conformer.get("dft_W_max"),
            "s1_dft_Hpop": conformer.get("dft_Hpop"),
            "s1_dft_N_total": conformer.get("dft_N_total"),
        }

        output_path = output_dir / "condition_features_mlr.csv"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=CONDITION_FEATURE_MLR_COLUMNS)
            writer.writeheader()
            writer.writerow({column: row.get(column) for column in CONDITION_FEATURE_MLR_COLUMNS})
        return output_path

    def _build_sample_id(self, manifest: Dict[str, Any], output_dir: Path) -> str:
        reaction_id = str(manifest.get("reaction_id") or self.reaction_root.name)
        condition_id = str(manifest.get("condition_id") or output_dir.name)
        branch_id = self._resolve_branch_id(output_dir)
        parts = [reaction_id, condition_id]
        if branch_id:
            parts.append(branch_id)
        return "__".join(part for part in parts if part)

    def _resolve_branch_id(self, output_dir: Optional[Path] = None) -> Optional[str]:
        if self.branch_id is not None:
            branch_id = str(self.branch_id).strip()
            if branch_id:
                return branch_id

        target_dir = output_dir or self.condition_root
        if target_dir.parent.name == "branches":
            branch_id = str(target_dir.name).strip()
            return branch_id or None
        return None

    def _run_shermo_for_log(
        self,
        *,
        shermo_bin: Path,
        freq_log: Path,
        sp_energy: float,
        output_file: Path,
        temperature_k: float,
        pressure_atm: Optional[float],
        scl_zpe: Optional[float],
        ilowfreq: Optional[int],
        imagreal: Optional[float],
        conc: Any,
    ) -> Optional[shermo_runner.ThermoResult]:
        return shermo_runner.run_shermo(
            shermo_bin=shermo_bin,
            freq_output=freq_log,
            sp_energy=sp_energy,
            output_file=output_file,
            temperature_k=temperature_k,
            pressure_atm=pressure_atm,
            scl_zpe=scl_zpe,
            ilowfreq=ilowfreq,
            imagreal=imagreal,
            conc=str(conc) if conc is not None else None,
        )

    def _compute_boltzmann_stats(
        self,
        energies_kcal: list[float],
        temperature_k: float,
    ) -> dict[str, float]:
        if not energies_kcal or temperature_k <= 0.0:
            return {
                "Nconf_eff": float("nan"),
                "Sconf": float("nan"),
                "E_span": float("nan"),
                "E_avg_weighted": float("nan"),
                "E_std": float("nan"),
                "W_max": float("nan"),
                "Hpop": float("nan"),
                "N_total": len(energies_kcal) if energies_kcal else 0,
            }

        min_energy = min(energies_kcal)
        relative_energies = [e - min_energy for e in energies_kcal]
        raw_weights = [
            math.exp(-energy / (R_KCAL_MOL_K * temperature_k))
            for energy in relative_energies
        ]
        total_weight = sum(raw_weights)
        if total_weight <= 0.0:
            return {
                "Nconf_eff": float("nan"),
                "Sconf": float("nan"),
                "E_span": max(relative_energies) - min(relative_energies) if relative_energies else float("nan"),
                "E_avg_weighted": float("nan"),
                "E_std": float("nan"),
                "W_max": float("nan"),
                "Hpop": float("nan"),
                "N_total": len(energies_kcal),
            }

        weights = [w / total_weight for w in raw_weights]
        entropy_cal = -R_CAL_MOL_K * sum(
            w * math.log(max(w, 1e-300)) for w in weights if w > 0.0
        )
        avg_energy = sum(w * e for w, e in zip(weights, energies_kcal))
        e_span = max(energies_kcal) - min(energies_kcal) if energies_kcal else 0.0
        e_std = math.sqrt(
            sum(w * (e - avg_energy) ** 2 for w, e in zip(weights, energies_kcal))
        )
        w_max = max(weights) if weights else 0.0

        return {
            "Nconf_eff": 1.0 / sum(w * w for w in weights) if weights else 0.0,
            "Sconf": entropy_cal,
            "E_span": e_span,
            "E_avg_weighted": avg_energy,
            "E_std": e_std,
            "W_max": w_max,
            "Hpop": -sum(w * math.log2(max(w, 1e-300)) for w in weights) if weights else float("nan"),
            "N_total": len(energies_kcal),
        }

    def _compute_conformer_boltzmann(self, temperature_k: float) -> Dict[str, Any]:
        result: Dict[str, Any] = {
            "g_precursor_boltzmann_avg_kcal": float("nan"),
            "conformer_N_eff": float("nan"),
            "conformer_S_conf_cal_mol_k": float("nan"),
            "conformer_top_k_energies": [],
            "n_conformers": 0,
            "source_dir": None,
            "ff_data_source": None,
            "ff_Nconf_eff": float("nan"),
            "ff_Sconf": float("nan"),
            "ff_E_span": float("nan"),
            "ff_E_avg_weighted": float("nan"),
            "ff_E_std": float("nan"),
            "ff_W_max": float("nan"),
            "ff_Hpop": float("nan"),
            "ff_N_total": 0,
            "dft_data_source": None,
            "dft_Nconf_eff": float("nan"),
            "dft_Sconf": float("nan"),
            "dft_E_span": float("nan"),
            "dft_E_avg_weighted": float("nan"),
            "dft_E_std": float("nan"),
            "dft_W_max": float("nan"),
            "dft_Hpop": float("nan"),
            "dft_N_total": 0,
        }

        if temperature_k <= 0.0:
            return result

        crest_dir = self._resolve_conformer_ensemble_dir(source="crest")
        if crest_dir is not None:
            crest_energies, crest_source = self._load_conformer_energy_values(crest_dir)
            if crest_energies:
                stats = self._compute_boltzmann_stats(crest_energies, temperature_k)
                result["source_dir"] = str(crest_dir)
                result["n_conformers"] = stats["N_total"]
                result["conformer_N_eff"] = stats["Nconf_eff"]
                result["conformer_S_conf_cal_mol_k"] = stats["Sconf"]
                result["conformer_top_k_energies"] = [
                    float(energy - min(crest_energies))
                    for energy in sorted(crest_energies)[:5]
                ]
                result["g_precursor_boltzmann_avg_kcal"] = (
                    stats["E_avg_weighted"]
                    - temperature_k * stats["Sconf"] / 1000.0
                )
                result["ff_data_source"] = crest_source
                result["ff_Nconf_eff"] = stats["Nconf_eff"]
                result["ff_Sconf"] = stats["Sconf"]
                result["ff_E_span"] = stats["E_span"]
                result["ff_E_avg_weighted"] = stats["E_avg_weighted"]
                result["ff_E_std"] = stats["E_std"]
                result["ff_W_max"] = stats["W_max"]
                result["ff_Hpop"] = stats["Hpop"]
                result["ff_N_total"] = stats["N_total"]

        dft_dir = self._resolve_conformer_ensemble_dir(source="finaldft")
        if dft_dir is not None:
            dft_energies, dft_source = self._load_conformer_energy_values(dft_dir)
            if dft_energies:
                stats = self._compute_boltzmann_stats(dft_energies, temperature_k)
                result["dft_data_source"] = dft_source
                result["dft_Nconf_eff"] = stats["Nconf_eff"]
                result["dft_Sconf"] = stats["Sconf"]
                result["dft_E_span"] = stats["E_span"]
                result["dft_E_avg_weighted"] = stats["E_avg_weighted"]
                result["dft_E_std"] = stats["E_std"]
                result["dft_W_max"] = stats["W_max"]
                result["dft_Hpop"] = stats["Hpop"]
                result["dft_N_total"] = stats["N_total"]

                if result["source_dir"] is None:
                    result["source_dir"] = str(dft_dir)
                    result["n_conformers"] = stats["N_total"]
                    result["conformer_N_eff"] = stats["Nconf_eff"]
                    result["conformer_S_conf_cal_mol_k"] = stats["Sconf"]
                    result["conformer_top_k_energies"] = [
                        float(energy - min(dft_energies))
                        for energy in sorted(dft_energies)[:5]
                    ]
                    result["g_precursor_boltzmann_avg_kcal"] = (
                        stats["E_avg_weighted"]
                        - temperature_k * stats["Sconf"] / 1000.0
                    )

        return result

    def _resolve_conformer_ensemble_dir(self, source: str = "any") -> Optional[Path]:
        for species_name in ("precursor", "product"):
            if source == "crest":
                candidate = self._s1_root / "S1_ConfGeneration" / species_name / "crest"
                if candidate.is_dir():
                    return candidate
            elif source == "fastsp":
                candidate = self._s1_root / "S1_ConfGeneration" / species_name / "fastsp"
                if candidate.is_dir():
                    return candidate
            elif source == "finaldft":
                return self._resolve_s1_species_dft_dir(species_name)
            else:
                for src in ("crest", "fastsp", "finaldft"):
                    result = self._resolve_conformer_ensemble_dir(source=src)
                    if result is not None:
                        return result
                return None
        return None

    @property
    def _s1_root(self) -> Path:
        """Root for S1 artifacts: artifact_root (branch dir) takes priority over reaction_root.

        In the existing data layout, S1/S3 artifacts live under branches/BR_*/,
        which is passed as artifact_root. In greenfield runs they may be directly
        under reaction_root.
        """
        return self.artifact_root or self.reaction_root

    def _resolve_s1_species_dft_dir(self, species_name: str) -> Optional[Path]:
        species_dir = self._s1_root / "S1_ConfGeneration" / species_name
        for dft_dir_name in ("finalDFT", "dft"):
            candidate = species_dir / dft_dir_name
            if candidate.is_dir():
                return candidate
        return None

    def _load_conformer_energy_values(self, search_dir: Path) -> tuple[list[float], str]:
        json_path = search_dir / "conformer_ensemble_energies.json"
        if json_path.exists():
            raw = read_json(json_path, default=None)
            if isinstance(raw, dict) and isinstance(raw.get("energies"), list):
                unit = str(raw.get("energy_unit") or "kcal/mol").strip().lower()
                source = str(raw.get("source") or search_dir.name)
                parsed = self._coerce_energy_list(raw["energies"], unit=unit)
                if parsed:
                    return parsed, source
            if isinstance(raw, list):
                parsed = [value for value in (_to_optional_float(item) for item in raw) if value is not None]
                if parsed:
                    return parsed, search_dir.name

        json_path = search_dir / "conformer_energies.json"
        energies = self._read_conformer_energy_json(json_path)
        if energies:
            return energies, "finaldft"

        return (
            [energy * HARTREE_TO_KCAL for _, energy in self._scan_log_energies(search_dir)],
            "log_scan",
        )

    def _read_conformer_energy_json(self, json_path: Path) -> list[float]:
        if not json_path.exists():
            return []
        raw = read_json(json_path, default=None)
        if isinstance(raw, list):
            return [value for value in (_to_optional_float(item) for item in raw) if value is not None]
        if not isinstance(raw, dict):
            return []

        unit = str(raw.get("unit") or "kcal/mol").strip().lower()
        if isinstance(raw.get("energies"), list):
            return self._coerce_energy_list(raw.get("energies") or [], unit=unit)

        if isinstance(raw.get("records"), list):
            values: list[float] = []
            for record in raw.get("records") or []:
                if not isinstance(record, dict):
                    continue
                for key in ("energy_kcal", "g_used_kcal", "sp_energy_kcal"):
                    value = _to_optional_float(record.get(key))
                    if value is not None:
                        values.append(value)
                        break
                else:
                    for key in ("sp_energy", "g_used", "energy"):
                        value = _to_optional_float(record.get(key))
                        if value is not None:
                            if unit in {"hartree", "au", "a.u."}:
                                values.append(value * HARTREE_TO_KCAL)
                            else:
                                values.append(value)
                            break
            return values
        return []

    def _coerce_energy_list(self, values: list[Any], unit: str) -> list[float]:
        parsed = [value for value in (_to_optional_float(item) for item in values) if value is not None]
        if unit in {"hartree", "au", "a.u."}:
            return [value * HARTREE_TO_KCAL for value in parsed]
        return parsed

    def _scan_log_energies(self, search_dir: Path) -> list[tuple[Path, float]]:
        if not search_dir.exists():
            return []

        candidates = sorted(
            [path for path in search_dir.rglob("*.log") if path.is_file()]
            + [path for path in search_dir.rglob("*.out") if path.is_file()]
        )
        energies: list[tuple[Path, float]] = []
        for path in candidates:
            energy = self._extract_last_scf_energy(path)
            if energy is not None:
                energies.append((path, energy))
        return energies

    def _extract_last_scf_energy(self, log_path: Path) -> Optional[float]:
        try:
            content = log_path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            return None
        matches = SCF_DONE_PATTERN.findall(content)
        if not matches:
            return None
        return _to_optional_float(matches[-1])

    def _resolve_s1_species_reference(self, species_name: str) -> Dict[str, Any]:
        species_dir = self._s1_root / "S1_ConfGeneration" / species_name
        dft_dir = self._resolve_s1_species_dft_dir(species_name)
        resolved: Dict[str, Any] = {"freq_log": None, "sp_energy": None, "dft_dir": dft_dir}
        if dft_dir is None:
            return resolved

        state = read_json_dict(species_dir / "conformer_state.json")
        summary = _as_dict(state.get("summary"))
        conformers = _as_dict(state.get("conformers"))
        best_name = summary.get("best_conformer")
        if isinstance(best_name, str) and best_name:
            best_entry = _as_dict(conformers.get(best_name))
            record = _as_dict(best_entry.get("record"))
            freq_log = _resolve_existing_path(record.get("log_file"), species_dir, dft_dir)
            if freq_log is not None:
                resolved["freq_log"] = freq_log
            resolved["sp_energy"] = _to_optional_float(record.get("sp_energy")) or _to_optional_float(record.get("g_used"))

        if resolved["sp_energy"] is None:
            resolved["sp_energy"] = _to_optional_float(summary.get("global_min_energy"))

        if resolved["freq_log"] is None or resolved["sp_energy"] is None:
            csv_log, csv_energy = self._resolve_s1_species_reference_from_csv(
                dft_dir / "conformer_thermo.csv",
                species_dir=species_dir,
                dft_dir=dft_dir,
            )
            if resolved["freq_log"] is None:
                resolved["freq_log"] = csv_log
            if resolved["sp_energy"] is None:
                resolved["sp_energy"] = csv_energy

        if resolved["freq_log"] is None or resolved["sp_energy"] is None:
            scanned = self._scan_log_energies(dft_dir)
            if scanned:
                best_log, best_energy = min(scanned, key=lambda item: item[1])
                if resolved["freq_log"] is None:
                    resolved["freq_log"] = best_log
                if resolved["sp_energy"] is None:
                    resolved["sp_energy"] = best_energy

        if resolved["freq_log"] is None:
            resolved["freq_log"] = _latest_log_file(dft_dir)
        if resolved["sp_energy"] is None and isinstance(resolved["freq_log"], Path):
            resolved["sp_energy"] = self._extract_last_scf_energy(resolved["freq_log"])
        return resolved

    def _resolve_s1_species_reference_from_csv(
        self,
        csv_path: Path,
        *,
        species_dir: Path,
        dft_dir: Path,
    ) -> Tuple[Optional[Path], Optional[float]]:
        if not csv_path.exists():
            return None, None

        best_log: Optional[Path] = None
        best_energy: Optional[float] = None
        try:
            with open(csv_path, "r", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    energy = _to_optional_float(row.get("sp_energy"))
                    if energy is None:
                        continue
                    if best_energy is None or energy < best_energy:
                        best_energy = energy
                        best_log = _resolve_existing_path(row.get("log_file"), species_dir, dft_dir)
        except Exception as exc:
            logger.warning(f"Failed to read conformer_thermo.csv for {species_dir.name}: {exc}")
        return best_log, best_energy

    def _compute_delta_g(
        self,
        reactants: Dict[str, int],
        products: Dict[str, int],
        *,
        species_map: Dict[str, Dict[str, Optional[float]]],
    ) -> Optional[float]:
        delta_g = 0.0
        for side, sign in ((reactants, -1.0), (products, 1.0)):
            for key, coefficient in side.items():
                value = _to_optional_float(_as_dict(species_map.get(key)).get("g_kcal"))
                if value is None:
                    return None
                delta_g += sign * float(value) * float(coefficient)
        return delta_g

    def _build_reaction_pairs(self) -> Dict[str, Tuple[Dict[str, int], Dict[str, int]]]:
        oxidation_reactants, oxidation_products = self._resolve_precursor_to_intermediate_stoichiometry()
        return {
            "oxidation": (oxidation_reactants, oxidation_products),
            "cycloaddition_activation": ({"intermediate": 1}, {"ts": 1}),
            "overall_activation": (
                dict(oxidation_reactants),
                self._replace_species_key(oxidation_products, "intermediate", "ts"),
            ),
            "overall_reaction": (
                self._replace_species_key(oxidation_reactants, "precursor", "intermediate"),
                self._replace_species_key(oxidation_products, "intermediate", "product"),
            ),
        }

    def _resolve_precursor_to_intermediate_stoichiometry(self) -> Tuple[Dict[str, int], Dict[str, int]]:
        reaction_profile = (self.config.get("run", {}) or {}).get("reaction_type", "[4+3]_default")
        ref_terms = (
            ((self.config.get("reaction_reference_terms", {}) or {}).get(reaction_profile, {}) or {}).get(
                "precursor_to_intermediate", {}
            )
            or {}
        )
        if not isinstance(ref_terms, dict):
            return {"precursor": 1}, {"intermediate": 1}

        reactants = self._normalize_stoichiometry_side(ref_terms.get("reactants"), fallback={"precursor": 1})
        products = self._normalize_stoichiometry_side(ref_terms.get("products"), fallback={"intermediate": 1})
        return reactants, products

    def _normalize_stoichiometry_side(self, side_terms: Any, fallback: Dict[str, int]) -> Dict[str, int]:
        if not isinstance(side_terms, dict):
            return dict(fallback)

        normalized: Dict[str, int] = {}
        for key, coefficient in side_terms.items():
            key_text = str(key)
            try:
                coefficient_int = int(coefficient)
            except (TypeError, ValueError):
                continue
            if coefficient_int <= 0:
                continue
            normalized[key_text] = coefficient_int
        return normalized or dict(fallback)

    def _replace_species_key(self, stoichiometry: Dict[str, int], source_key: str, target_key: str) -> Dict[str, int]:
        replaced: Dict[str, int] = {}
        for key, coefficient in stoichiometry.items():
            new_key = target_key if key == source_key else key
            replaced[new_key] = replaced.get(new_key, 0) + coefficient
        return replaced

    def _format_reaction_label(self, reactants: Dict[str, int], products: Dict[str, int]) -> str:
        return f"{self._format_reaction_side(reactants)} -> {self._format_reaction_side(products)}"

    def _format_reaction_side(self, species: Dict[str, int]) -> str:
        terms = []
        for key, coefficient in species.items():
            label = "TS" if key == "ts" else key
            terms.append(f"{coefficient} {label}" if coefficient > 1 else label)
        return " + ".join(terms)

    def _load_condition_temperature_k(self) -> float:
        manifest_path = self.condition_root / "condition_manifest.json"
        if not manifest_path.exists():
            return float((self.config.get("thermo", {}) or {}).get("temperature_k", 298.15) or 298.15)
        try:
            data = read_json_dict(manifest_path)
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
                resume_state = read_json_dict(resume_path)
            except Exception as exc:
                logger.warning(f"Failed to read s3_resume.json: {exc}")

        ts_from_resume = _safe_path_from_json((resume_state.get("ts_opt") or {}).get("freq_log"), s3_dir)
        if ts_from_resume is None:
            ts_from_resume = _safe_path_from_json((resume_state.get("ts_opt") or {}).get("log_file"), s3_dir)

        intermediate_from_resume = _safe_path_from_json((resume_state.get("reactant_opt") or {}).get("freq_log"), s3_dir)
        if intermediate_from_resume is None:
            intermediate_from_resume = _safe_path_from_json((resume_state.get("reactant_opt") or {}).get("log_file"), s3_dir)

        ts_log = ts_from_resume
        if ts_log is None:
            ts_log = _latest_log_file(s3_dir / "ts_opt")

        intermediate_log = intermediate_from_resume
        if intermediate_log is None:
            intermediate_log = _latest_log_file(s3_dir / "reactant_opt")

        return ts_log, intermediate_log

    def _reference_small_molecule_keys(self) -> list[str]:
        reactants, products = self._resolve_precursor_to_intermediate_stoichiometry()
        keys: list[str] = []
        seen = set()
        for side_terms in (reactants, products):
            for key in side_terms:
                if key in {"precursor", "intermediate"} or key in seen:
                    continue
                seen.add(key)
                keys.append(str(key))
        return keys

    def _add_reference_species(self, species: Dict[str, Dict[str, Optional[float]]]) -> None:
        small_mol_keys = self._reference_small_molecule_keys()
        if not small_mol_keys:
            return

        cache_root = (self.config.get("global", {}) or {}).get("small_molecule_cache_dir")
        if not cache_root:
            logger.warning(
                "small_molecule_cache_dir not configured; cannot resolve reference small-molecule thermo"
            )
            return

        catalog = SmallMoleculeCatalog(self.config)
        cache = SmallMoleculeCache(Path(cache_root))

        for key in small_mol_keys:
            if key not in species:
                species[key] = {"g_kcal": None, "h_kcal": None, "s_cal_mol_k": None}

            molecule = catalog.get(key)
            if molecule is None:
                logger.warning(f"Missing small-molecule catalog entry for {key}")
                continue

            cache_dir = cache.get_path(molecule.smiles)
            if cache_dir is None:
                continue

            thermo_path = ensure_thermo_json_from_entry(cache_dir)
            if thermo_path is None:
                continue

            thermo_record = read_thermo_record(thermo_path)
            if thermo_record is None:
                continue
            thermo_data = {
                "g_kcal": thermo_record.g_kcal,
                "h_kcal": thermo_record.h_kcal,
                "s_cal_mol_k": thermo_record.s_cal_mol_k,
            }

            species[key] = {
                "g_kcal": _to_optional_float(thermo_data.get("g_kcal")),
                "h_kcal": _to_optional_float(thermo_data.get("h_kcal")),
                "s_cal_mol_k": _to_optional_float(thermo_data.get("s_cal_mol_k")),
            }


ConditionThermoCalculator = ConditionThermoEngine
