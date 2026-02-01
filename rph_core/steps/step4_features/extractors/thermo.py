"""
Step 4: Thermo Extractor (V4.2 Phase A)
============================================

Thermo features extractor using SPMatrixReport.

Author: QC Descriptors Team
Date: 2026-01-18
"""

from typing import Dict, Any, List, cast
import numpy as np

from .base import BaseExtractor, register_extractor


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
        return ["sp_report"]

    def extract(self, context) -> Dict[str, Any]:
        sp_report = context.sp_report

        if sp_report is None:
            raise ValueError("sp_report is None")

        features = {}

        # Extract Gibbs energies (kcal/mol)
        g_ts = getattr(sp_report, "g_ts", None)
        g_reactant = getattr(sp_report, "g_reactant", None)
        g_product = getattr(sp_report, "g_product", None)

        # Extract electronic energies (Hartree)
        e_ts = getattr(sp_report, "e_ts_final", None) or getattr(sp_report, "e_ts", None)
        e_reactant = getattr(sp_report, "e_reactant", None)
        e_product = getattr(sp_report, "e_product", None)

        has_gibbs_activation = g_ts is not None and g_reactant is not None
        has_gibbs_reaction = g_product is not None and g_reactant is not None
        
        # Activation energy
        if has_gibbs_activation:
            g_ts_val = cast(float, g_ts)
            g_reactant_val = cast(float, g_reactant)
            features["thermo.dG_activation"] = g_ts_val - g_reactant_val
            features["thermo.dG_activation_gibbs"] = features["thermo.dG_activation"]
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
