"""
Step 4: Precursor Geometry Extractor (V6.2 Future Extension)
==============================================================

Extract geometry features from precursor molecule.

Features (prec_geom.* prefix):
- natoms: Number of atoms in precursor
- r_gyration: Radius of gyration (Å)
- max_bond_length: Maximum bond length (Å)
- min_bond_length: Minimum bond length (Å)
- avg_bond_length: Average bond length (Å)

This is an optional extension - extraction only runs if precursor_xyz is available.

Author: RPH Team
Date: 2026-02-02
"""

from typing import Dict, Any, List, Optional
import numpy as np

from .base import BaseExtractor, register_extractor
from rph_core.utils.file_io import read_xyz
from rph_core.utils.geometry_tools import GeometryUtils


class PrecursorGeometryExtractor(BaseExtractor):
    """Extract geometric features from precursor XYZ structure.

    V6.2 Extension: Provides precursor geometry context for reaction analysis.
    """

    def get_plugin_name(self) -> str:
        return "precursor_geometry"

    def get_required_inputs(self) -> List[str]:
        return []  # Optional - runs if precursor_xyz available

    def get_required_inputs_for_context(self, context) -> List[str]:
        """Dynamic: only require if precursor_xyz exists."""
        if hasattr(context, 's1_precursor_xyz') and context.s1_precursor_xyz:
            return ['s1_precursor_xyz']
        return []

    def extract(self, context) -> Dict[str, Any]:
        precursor_xyz = getattr(context, 's1_precursor_xyz', None)

        if precursor_xyz is None:
            # Return empty features if no precursor available
            return {}

        try:
            coords, symbols = read_xyz(precursor_xyz)
            if coords is None or symbols is None:
                return {}

            n_atoms = len(symbols)
            distance_matrix = GeometryUtils.compute_distance_matrix(coords)

            features = {}
            features["prec_geom.natoms"] = n_atoms

            # Radius of gyration
            center_of_mass = np.mean(coords, axis=0)
            squared_distances = np.sum((coords - center_of_mass) ** 2, axis=1)
            features["prec_geom.r_gyration"] = np.sqrt(np.mean(squared_distances))

            # Bond statistics (non-bonded cutoff = 2.0 Å)
            bond_mask = (distance_matrix > 0) & (distance_matrix < 2.5)
            if np.any(bond_mask):
                bond_lengths = distance_matrix[bond_mask]
                features["prec_geom.max_bond_length"] = float(np.max(bond_lengths))
                features["prec_geom.min_bond_length"] = float(np.min(bond_lengths))
                features["prec_geom.avg_bond_length"] = float(np.mean(bond_lengths))
            else:
                features["prec_geom.max_bond_length"] = float("nan")
                features["prec_geom.min_bond_length"] = float("nan")
                features["prec_geom.avg_bond_length"] = float("nan")

            return features

        except Exception:
            return {}


register_extractor(PrecursorGeometryExtractor())
