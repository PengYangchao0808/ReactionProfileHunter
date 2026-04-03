"""
Step 4: Thermo Extractor (V4.2 Phase A)
============================================

Thermo features extractor using SPMatrixReport.

Author: QC Descriptors Team
Date: 2026-01-18
"""

from typing import Dict, Any, List, Optional, Tuple, cast
import numpy as np
import re
from pathlib import Path
from dataclasses import dataclass

from .base import BaseExtractor, register_extractor
from rph_core.utils.constants import HARTREE_TO_KCAL


@dataclass
class QCSignature:
    method: str
    basis: str
    solvent: str
    
    def matches(self, other: "QCSignature", strict: bool = True) -> bool:
        if strict:
            return (self.method.lower() == other.method.lower() and
                    self.basis.lower() == other.basis.lower() and
                    self.solvent.lower() == other.solvent.lower())
        return self.method.lower() == other.method.lower() and self.basis.lower() == other.basis.lower()


def _parse_qc_signature_from_log(log_path: Path) -> Optional[QCSignature]:
    try:
        text = Path(log_path).read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None

    method = ""
    basis = ""
    solvent = ""

    route_match = re.search(r"#[pPnNtT]?\s+(.+?)(?:\n-{4,}|\n\s*\n)", text, re.DOTALL)
    if route_match:
        route = route_match.group(1).replace("\n", " ").strip()
        
        method_patterns = [
            r"\b(RHF|UHF|ROHF|HF)\b",
            r"\b(B3LYP|B3PW91|BLYP|BP86|PBE|PBE0|M06|M06-2X|M06-L|wB97X-D|wB97X-D3|wB97X-D3BJ|wB97M-D3BJ|CAM-B3LYP)\b",
            r"\b(MP2|MP3|MP4|CCSD|CCSD\(T\)|CIS|CISD)\b",
        ]
        for pat in method_patterns:
            m = re.search(pat, route, re.IGNORECASE)
            if m:
                method = m.group(1).upper()
                break
        
        basis_patterns = [
            r"\b(def2-TZVPP?|def2-SVP|def2-QZVP)\b",
            r"\b(cc-pVDZ|cc-pVTZ|cc-pVQZ|aug-cc-pVDZ|aug-cc-pVTZ)\b",
            r"\b(6-31G\*?|6-31\+G\*?|6-31G\(d,p\)|6-311G\*?|6-311\+G\*?|6-311\+\+G\*?\*?)\b",
            r"\b(STO-3G|3-21G|4-31G)\b",
        ]
        for pat in basis_patterns:
            m = re.search(pat, route, re.IGNORECASE)
            if m:
                basis = m.group(1)
                break
        
        solvent_match = re.search(r"SCRF\s*=\s*\([^)]*solvent\s*=\s*(\w+)", route, re.IGNORECASE)
        if solvent_match:
            solvent = solvent_match.group(1).lower()
        elif "scrf" in route.lower() or "pcm" in route.lower() or "smd" in route.lower():
            solvent = "implicit"

    if not method:
        scf_match = re.search(r"SCF Done:\s+E\(([RU]?)(HF|B3LYP|B3PW91|MP2|M06|PBE0?)\)", text)
        if scf_match:
            method = scf_match.group(2).upper()
            if scf_match.group(1):
                method = scf_match.group(1).upper() + method
    
    if not basis:
        basis_match = re.search(r"Standard basis:\s+([\w\-\+\*\(\),]+)", text)
        if basis_match:
            basis = basis_match.group(1).strip()
    
    if not solvent:
        if re.search(r"Solvent\s*:\s*(\w+)", text, re.IGNORECASE):
            solvent_match = re.search(r"Solvent\s*:\s*(\w+)", text, re.IGNORECASE)
            if solvent_match:
                solvent = solvent_match.group(1).lower()
        elif "SMD" in text or "PCM" in text or "CPCM" in text:
            solvent = "implicit"

    return QCSignature(method=method, basis=basis, solvent=solvent)


class ThermoExtractor(BaseExtractor):
    """Extract thermodynamic features from SPMatrixReport.

    Features (thermo.* prefix):
    - dG_activation, dG_reaction: Gibbs free energy differences (kcal/mol)
    - dE_activation, dE_reaction: Electronic energy differences (kcal/mol)
    - energy_source_activation, energy_source_reaction: "gibbs" or "electronic"
    - method: Computational method string
    - solvent: Solvent string

    Uses Gibbs priority strategy: prefer g_* if available, fallback to e_*.
    """

    def get_plugin_name(self) -> str:
        return "thermo"

    def get_required_inputs(self) -> List[str]:
        # V6.2 extract-only: allow thermo backfill from existing logs even when
        # SPMatrixReport is absent.
        return []

    def extract(self, context) -> Dict[str, Any]:
        sp_report = context.sp_report

        # If upstream did not provide an SPMatrixReport, create a lightweight
        # placeholder so downstream QCChecks can still record qc.has_gibbs.
        if sp_report is None:
            class _FallbackSPReport:
                g_ts: Optional[float] = None
                g_reactant: Optional[float] = None
                g_product: Optional[float] = None
                e_ts_final: Optional[float] = None
                e_ts: Optional[float] = None
                e_reactant: Optional[float] = None
                e_product: Optional[float] = None
                method: str = ""
                solvent: str = ""

                def validate(self) -> bool:
                    return True

            sp_report = _FallbackSPReport()

        features = {}

        # Extract Gibbs energies (kcal/mol)
        g_ts = getattr(sp_report, "g_ts", None)
        g_reactant = getattr(sp_report, "g_reactant", None)
        g_product = getattr(sp_report, "g_product", None)

        # V6.2: Gibbs should come from Shermo (SP energy + freq log).
        # Do not backfill from Gaussian free energies when Shermo is missing.
        def _set_sp_report_attr(name: str, value: Optional[float]) -> None:
            if value is None:
                return
            try:
                setattr(sp_report, name, value)
            except Exception:
                return

        def _best_product_shermo_sum(s1_dir: Path) -> Optional[Path]:
            patterns = [
                s1_dir / "product" / "dft" / "*_Shermo.sum",
                s1_dir / "product" / "dft" / "*Shermo*.sum",
                s1_dir / "**" / "product" / "**" / "*_Shermo.sum",
            ]
            import glob

            candidates: List[Path] = []
            for pat in patterns:
                for m in glob.glob(str(pat), recursive=True):
                    candidates.append(Path(m))

            uniq: List[Path] = []
            seen = set()
            for p in candidates:
                p = Path(p)
                if p in seen or not p.exists() or not p.is_file():
                    continue
                seen.add(p)
                uniq.append(p)
            if not uniq:
                return None

            try:
                from rph_core.utils.shermo_runner import _parse_sum_file
            except Exception:
                return None

            best: Optional[Tuple[float, Path]] = None
            for p in uniq:
                try:
                    thermo = _parse_sum_file(p)
                    g_val = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
                    if g_val is None:
                        continue
                    if best is None or g_val < best[0]:
                        best = (g_val, p)
                except Exception:
                    continue
            return best[1] if best else None

        # TS Gibbs (no fallback)
        # Reactant Gibbs (no fallback)

        # Product Gibbs (allow Shermo .sum from S1)
        if g_product is None:
            s1_dir = getattr(context, "s1_dir", None)
            if s1_dir is not None and Path(s1_dir).exists():
                best_sum = _best_product_shermo_sum(Path(s1_dir))
                if best_sum is not None:
                    try:
                        from rph_core.utils.shermo_runner import _parse_sum_file
                        thermo = _parse_sum_file(best_sum)
                        g_val = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
                        if g_val is not None:
                            g_product = g_val * HARTREE_TO_KCAL
                            _set_sp_report_attr("g_product", g_product)
                    except Exception as exc:
                        trace = context.get_plugin_trace(self.get_plugin_name())
                        trace.warnings.append(f"W_THERMO_PRODUCT_SUM_PARSE_FAILED:{exc}")

        # Extract electronic energies (Hartree)
        e_ts = getattr(sp_report, "e_ts_final", None) or getattr(sp_report, "e_ts", None)
        e_reactant = getattr(sp_report, "e_reactant", None)
        e_product = getattr(sp_report, "e_product", None)

        has_gibbs_activation = g_ts is not None and g_reactant is not None
        has_gibbs_reaction = g_product is not None and g_reactant is not None
        
        ts_log_path = getattr(context, "ts_log", None) or getattr(context, "s3_ts_log", None)
        intermediate_log_path = getattr(context, "intermediate_log", None) or getattr(context, "s3_intermediate_log", None)
        
        ts_sig: Optional[QCSignature] = None
        intermediate_sig: Optional[QCSignature] = None
        lot_mismatch_activation = False
        
        if ts_log_path and Path(ts_log_path).exists():
            ts_sig = _parse_qc_signature_from_log(Path(ts_log_path))
        if intermediate_log_path and Path(intermediate_log_path).exists():
            intermediate_sig = _parse_qc_signature_from_log(Path(intermediate_log_path))
        
        if ts_sig and intermediate_sig:
            if not ts_sig.matches(intermediate_sig, strict=False):
                lot_mismatch_activation = True
                features["thermo.lot_mismatch_activation"] = True
                features["thermo.ts_method"] = f"{ts_sig.method}/{ts_sig.basis}"
                features["thermo.reactant_method"] = f"{intermediate_sig.method}/{intermediate_sig.basis}"
        
        # Activation energy
        if has_gibbs_activation:
            g_ts_val = cast(float, g_ts)
            g_reactant_val = cast(float, g_reactant)
            dG_act = g_ts_val - g_reactant_val
            
            if lot_mismatch_activation:
                features["thermo.dG_activation"] = np.nan
                features["thermo.dG_activation_gibbs"] = np.nan
                features["thermo.energy_source_activation"] = "lot_mismatch"
            elif abs(dG_act) > 200.0:
                features["thermo.dG_activation"] = np.nan
                features["thermo.dG_activation_gibbs"] = np.nan
                features["thermo.energy_source_activation"] = "out_of_range"
            else:
                features["thermo.dG_activation"] = dG_act
                features["thermo.dG_activation_gibbs"] = dG_act
                features["thermo.energy_source_activation"] = "gibbs"
        elif e_ts is not None and e_reactant is not None:
            features["thermo.dG_activation"] = (e_ts - e_reactant) * 627.509
            features["thermo.dG_activation_gibbs"] = np.nan
            features["thermo.energy_source_activation"] = "electronic"
        else:
            features["thermo.dG_activation"] = np.nan
            features["thermo.dG_activation_gibbs"] = np.nan
            features["thermo.energy_source_activation"] = "none"
        
        # Reaction energy
        if has_gibbs_reaction:
            g_product_val = cast(float, g_product)
            g_reactant_val = cast(float, g_reactant)
            features["thermo.dG_reaction"] = g_product_val - g_reactant_val
            features["thermo.dG_reaction_gibbs"] = features["thermo.dG_reaction"]
            features["thermo.energy_source_reaction"] = "gibbs"
        elif e_product is not None and e_reactant is not None:
            features["thermo.dG_reaction"] = (e_product - e_reactant) * 627.509
            features["thermo.dG_reaction_gibbs"] = np.nan
            features["thermo.energy_source_reaction"] = "electronic"
        else:
            features["thermo.dG_reaction"] = np.nan
            features["thermo.dG_reaction_gibbs"] = np.nan
            features["thermo.energy_source_reaction"] = "none"

        # Electronic energy differences (always compute if available)
        if e_ts is not None and e_reactant is not None:
            features["thermo.dE_activation"] = (e_ts - e_reactant) * 627.509
        else:
            features["thermo.dE_activation"] = np.nan

        if e_product is not None and e_reactant is not None:
            features["thermo.dE_reaction"] = (e_product - e_reactant) * 627.509
        else:
            features["thermo.dE_reaction"] = np.nan

        # Method and solvent
        features["thermo.method"] = getattr(sp_report, "method", "")
        features["thermo.solvent"] = getattr(sp_report, "solvent", "")

        return features


register_extractor(ThermoExtractor())
