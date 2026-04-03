# pyright: ignore
"""
Step 4: Step1 Activation Extractor (V6.4)
============================================

Extract Step1 (precursor activation) features:
- Thermodynamic features from Shermo
- Conformer entropy and ensemble properties
- V6.4: Hardcoded geometry features deprecated

Author: RPH Team
Date: 2026-02-02
Updated: 2026-03-17 (V6.4: Remove hardcoded geometry)
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, Any, List, Optional

from .base import BaseExtractor, register_extractor
from ..schema import WARNING_CODES
from rph_core.utils.file_io import read_xyz
from rph_core.utils.constants import HARTREE_TO_KCAL


class Step1ActivationExtractor(BaseExtractor):
    """Extract Step1 activation features from S1 artifacts."""

    def get_plugin_name(self) -> str:
        return "step1_activation"

    def get_required_inputs(self) -> List[str]:
        """Return list of required context inputs.

        V6.2: All inputs are optional to allow graceful degradation.
        Extractors emit NaN values and warning codes when files are missing.
        """
        return []

    def extract(self, context) -> Dict[str, Any]:
        """Extract Step1 activation features.

        Args:
            context: FeatureContext with S1 path handles

        Returns:
            Dictionary of Step1 features with s1_* prefix
        """
        features = {}
        trace = context.get_plugin_trace(self.get_plugin_name())

        # Thermodynamic features
        thermo_features = self._extract_thermo_features(context, trace)
        features.update(thermo_features)

        # Conformer ensemble features
        conformer_features = self._extract_conformer_features(context, trace)
        features.update(conformer_features)

        # Leaving group geometry
        lg_geometry_features = self._extract_lg_geometry(context)
        features.update(lg_geometry_features)

        # Alpha-H features
        alpha_h_features = self._extract_alpha_h_features(context)
        features.update(alpha_h_features)

        return features

    def _extract_thermo_features(self, context, trace) -> Dict[str, Any]:
        """Extract thermodynamic activation features.
        
        Physical picture: Precursor → S3_Reactant + Leaving_Group
        dG_act = G(S3_reactant) + G(leaving_group) - G(precursor)
        
        Note: The "ylide" species IS the S3_reactant (reactive intermediate).
        S1 only provides: precursor, leaving_group (HOAc), product.
        """
        features = {}
        s1_dir = getattr(context, "s1_dir", None)
        
        g_leaving_group = self._get_leaving_group_gibbs(context, trace, s1_dir)
        features['s1_leaving_group_G'] = g_leaving_group
        
        g_precursor = self._get_precursor_gibbs(context, trace, s1_dir)
        features['s1_precursor_G'] = g_precursor
        
        g_s3_intermediate = getattr(context, "s3_intermediate_g_kcal", None)
        if g_s3_intermediate is None:
            sp_report = getattr(context, "sp_report", None)
            g_s3_intermediate = getattr(sp_report, "g_intermediate", None) if sp_report else None
            if g_s3_intermediate is None:
                g_source = getattr(sp_report, "g_intermediate_source", None) if sp_report else None
                if g_source in ("shermo_failed", "missing_freq_log"):
                    trace.warnings.append("W_S3_SHERMO_GIBBS_FAILED")
                else:
                    trace.warnings.append("W_S1_MISSING_S3_INTERMEDIATE_GIBBS")
        
        if g_s3_intermediate is not None and g_leaving_group is not None and g_precursor is not None:
            if not np.isnan(g_s3_intermediate) and not np.isnan(g_leaving_group) and not np.isnan(g_precursor):
                dG_act = g_s3_intermediate + g_leaving_group - g_precursor
                
                if abs(dG_act) > 200.0:
                    trace.warnings.append("W_S1_DG_ACT_OUT_OF_RANGE")
                    features['s1_dG_act'] = float('nan')
                else:
                    features['s1_dG_act'] = dG_act
            else:
                trace.warnings.append("W_S1_MISSING_THERMO_COMPONENT")
                features['s1_dG_act'] = float('nan')
        else:
            trace.warnings.append("W_S1_MISSING_THERMO_COMPONENT")
            features['s1_dG_act'] = float('nan')

        dG_act = features.get('s1_dG_act', float('nan'))
        if not np.isnan(dG_act):
            R = 1.987e-3
            T = context.temperature_K
            features['s1_Keq_act'] = np.exp(-dG_act / (R * T))
        else:
            features['s1_Keq_act'] = float('nan')

        return features

    def _get_leaving_group_gibbs(self, context, trace, s1_dir) -> Optional[float]:
        """Get leaving group (HOAc) Gibbs energy from S1."""
        hoac_path = context.s1_hoac_thermo_file

        if (hoac_path is None or not Path(hoac_path).exists()) and s1_dir is not None:
            s1_dir = Path(s1_dir)
            canonical = s1_dir / "small_molecules" / "HOAc" / "thermo.json"
            if canonical.exists():
                hoac_path = canonical
            else:
                try:
                    from rph_core.utils.shermo_runner import find_shermo_sum_files, derive_hoac_thermo_from_sum
                    sum_files = find_shermo_sum_files(s1_dir)
                    hoac_sum = sum_files.get("hoac") or sum_files.get("leaving_group")
                    if hoac_sum is not None and Path(hoac_sum).exists():
                        derive_hoac_thermo_from_sum(Path(hoac_sum), canonical)
                        hoac_path = canonical if canonical.exists() else hoac_path
                except Exception as exc:
                    trace.warnings.append(f"W_S1_HOAC_DERIVE_FAILED:{exc}")

        if hoac_path is not None and Path(hoac_path).exists():
            try:
                with open(hoac_path, 'r') as f:
                    hoac_data = json.load(f)
                return hoac_data.get('G', hoac_data.get('g', float('nan')))
            except Exception:
                trace.warnings.append("W_S1_MISSING_HOAC_THERMO")
                return float('nan')
        else:
            trace.warnings.append("W_S1_MISSING_HOAC_THERMO")
            return float('nan')

    def _get_precursor_gibbs(self, context, trace, s1_dir) -> Optional[float]:
        """Get precursor Gibbs energy from S1 Shermo summary or .sum files."""
        shermo_path = context.s1_shermo_summary_file
        if (shermo_path is None or not Path(shermo_path).exists()) and s1_dir is not None:
            s1_dir = Path(s1_dir)
            candidate = s1_dir / "shermo_summary.json"
            if candidate.exists():
                shermo_path = candidate

        if shermo_path and Path(shermo_path).exists():
            try:
                with open(shermo_path, 'r') as f:
                    shermo_data = json.load(f)
                g_precursor = shermo_data.get('g_precursor') or shermo_data.get('G_precursor')
                if g_precursor is not None:
                    return g_precursor
            except Exception as exc:
                trace.warnings.append(f"W_S1_SHERMO_SUMMARY_READ_FAILED:{exc}")

        if s1_dir is not None:
            try:
                from rph_core.utils.shermo_runner import find_shermo_sum_files, _parse_sum_file
                sum_files = find_shermo_sum_files(Path(s1_dir))
                precursor_sum = sum_files.get("precursor")
                if precursor_sum is not None and Path(precursor_sum).exists():
                    thermo = _parse_sum_file(Path(precursor_sum))
                    g_val = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
                    return g_val * HARTREE_TO_KCAL
            except Exception as exc:
                trace.warnings.append(f"W_S1_PRECURSOR_SUM_PARSE_FAILED:{exc}")

        trace.warnings.append("W_S1_MISSING_PRECURSOR_THERMO")
        return float('nan')

    def _parse_s3_intermediate_gibbs(self, context, trace) -> Optional[float]:
        """Parse S3 intermediate Gibbs energy from Gaussian log."""
        import re
        intermediate_log = getattr(context, "s3_intermediate_log", None) or getattr(context, "intermediate_log", None)
        if intermediate_log is None:
            s3_dir = getattr(context, "s3_dir", None)
            if s3_dir is not None:
                candidates = [
                    Path(s3_dir) / "S3_intermediate_opt" / "standard" / "intermediate.log",
                    Path(s3_dir) / "S3_intermediate_opt" / "rescue" / "intermediate.log",
                ]
                for cand in candidates:
                    if cand.exists():
                        reactant_log = cand
                        break

        if reactant_log is None or not Path(reactant_log).exists():
            trace.warnings.append("W_S1_MISSING_S3_REACTANT_LOG")
            return None

        try:
            text = Path(reactant_log).read_text(encoding="utf-8", errors="ignore")
            matches = re.findall(
                r"Sum of electronic and thermal Free Energies=\s*([+\-]?\d+\.\d+)",
                text,
            )
            if matches:
                return float(matches[-1]) * HARTREE_TO_KCAL
        except Exception as exc:
            trace.warnings.append(f"W_S1_REACTANT_GIBBS_PARSE_FAILED:{exc}")

        trace.warnings.append("W_S1_MISSING_S3_REACTANT_GIBBS")
        return None

    def _extract_conformer_features(self, context, trace) -> Dict[str, Any]:
        """Extract conformer ensemble features."""
        features = {}

        conformer_path = context.s1_conformer_energies_file
        if conformer_path and Path(conformer_path).exists():
            try:
                with open(conformer_path, 'r') as f:
                    conformer_data = json.load(f)

                # Parse conformer energies
                energies = []
                if isinstance(conformer_data, dict):
                    energies = conformer_data.get('energies', [])
                elif isinstance(conformer_data, list):
                    energies = conformer_data

                if energies:
                    energies = [float(e) for e in energies]

                    # Calculate Boltzmann-weighted properties
                    R = 1.987e-3  # kcal/mol/K
                    T = context.temperature_K

                    # Boltzmann weights
                    min_e = min(energies)
                    weights = [np.exp(-(e - min_e) / (R * T)) for e in energies]
                    total_weight = sum(weights)

                    if total_weight > 0:
                        weights = [w / total_weight for w in weights]

                        # Effective number of conformers (inverse participation ratio)
                        N_eff = 1.0 / sum(w**2 for w in weights)
                        features['s1_Nconf_eff'] = N_eff

                        # Total conformer count
                        features['s1_Nconf_total'] = len(energies)

                        # Energy span (max - min)
                        features['s1_E_span'] = max(energies) - min(energies)

                        # Boltzmann-weighted entropy approximation
                        # S = -k * sum(p_i * ln(p_i))
                        S_conf = -R * sum(w * np.log(w + 1e-10) for w in weights)
                        features['s1_Sconf'] = S_conf

                        # P1 Enhancement: Boltzmann-weighted average energy (thermodynamic average)
                        # E_avg = sum(p_i * E_i)
                        E_avg_weighted = sum(w * e for w, e in zip(weights, energies))
                        features['s1_E_avg_weighted'] = E_avg_weighted

                        # P1 Enhancement: Gibbs free energy average (includes entropy term)
                        # G_avg = E_avg - T * S_conf
                        G_avg_weighted = E_avg_weighted - T * S_conf / 1000.0  # S in cal/mol/K, convert to kcal/mol/K
                        features['s1_G_avg_weighted'] = G_avg_weighted

                        # P1 Enhancement: Energy variance (conformational flexibility)
                        # Var(E) = sum(p_i * (E_i - E_avg)^2)
                        E_variance = sum(w * (e - E_avg_weighted)**2 for w, e in zip(weights, energies))
                        features['s1_E_variance'] = E_variance

                        # P1 Enhancement: Standard deviation of energies
                        features['s1_E_std'] = np.sqrt(E_variance) if E_variance > 0 else 0.0
                    else:
                        trace.warnings.append("W_S1_NO_CONFORMER_ENSEMBLE")
                        features['s1_Nconf_eff'] = float('nan')
                        features['s1_Nconf_total'] = len(energies)
                        features['s1_E_span'] = float('nan')
                        features['s1_Sconf'] = float('nan')
                        # P1 fallbacks
                        features['s1_E_avg_weighted'] = float('nan')
                        features['s1_G_avg_weighted'] = float('nan')
                        features['s1_E_variance'] = float('nan')
                        features['s1_E_std'] = float('nan')
                else:
                    trace.warnings.append("W_S1_NO_CONFORMER_ENSEMBLE")
                    features['s1_Nconf_eff'] = float('nan')
                    features['s1_Nconf_total'] = 0
                    features['s1_E_span'] = float('nan')
                    features['s1_Sconf'] = float('nan')
                    # P1 fallbacks
                    features['s1_E_avg_weighted'] = float('nan')
                    features['s1_G_avg_weighted'] = float('nan')
                    features['s1_E_variance'] = float('nan')
                    features['s1_E_std'] = float('nan')

            except Exception as e:
                trace.warnings.append("W_S1_NO_CONFORMER_ENSEMBLE")
                features['s1_Nconf_eff'] = float('nan')
                features['s1_Nconf_total'] = float('nan')
                features['s1_E_span'] = float('nan')
                features['s1_Sconf'] = float('nan')
                # P1 fallbacks
                features['s1_E_avg_weighted'] = float('nan')
                features['s1_G_avg_weighted'] = float('nan')
                features['s1_E_variance'] = float('nan')
                features['s1_E_std'] = float('nan')
        else:
            trace.warnings.append("W_S1_NO_CONFORMER_ENSEMBLE")
            features['s1_Nconf_eff'] = float('nan')
            features['s1_Nconf_total'] = float('nan')
            features['s1_E_span'] = float('nan')
            features['s1_Sconf'] = float('nan')
            # P1 fallbacks
            features['s1_E_avg_weighted'] = float('nan')
            features['s1_G_avg_weighted'] = float('nan')
            features['s1_E_variance'] = float('nan')
            features['s1_E_std'] = float('nan')

        return features

    def _extract_lg_geometry(self, context) -> Dict[str, Any]:
        """Extract leaving group geometry features.

        V6.4: DEPRECATED - Hardcoded geometry thresholds removed.
        These features are chemically fragile and require manual validation.
        Returns NaN by default.
        """
        features = {
            's1_d_C_O_lg': float('nan'),
            's1_angle_O_C_O': float('nan'),
            's1_tau_lg_ring': float('nan'),
        }

        trace = context.get_plugin_trace(self.get_plugin_name())
        trace.warnings.append("W_S1_GEOMETRY_DEPRECATED: Hardcoded geometry extraction disabled in V6.4")

        return features

    def _extract_alpha_h_features(self, context) -> Dict[str, Any]:
        """Extract alpha-H features relevant to enolization.

        V6.4: DEPRECATED - Hardcoded geometry thresholds removed.
        These features are chemically fragile and require manual validation.
        Returns NaN by default.
        """
        features = {
            's1_d_C_H_alpha': float('nan'),
            's1_tau_CH_C_O': float('nan'),
        }

        trace = context.get_plugin_trace(self.get_plugin_name())
        trace.warnings.append("W_S1_ALPHA_H_DEPRECATED: Hardcoded alpha-H geometry extraction disabled in V6.4")

        return features


register_extractor(Step1ActivationExtractor())
