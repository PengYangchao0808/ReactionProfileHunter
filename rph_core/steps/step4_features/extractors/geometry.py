"""
Step 4: Geometry Extractor (V4.2 Phase A)
==============================================

Geometry features extractor from TS structure.

Author: QC Descriptors Team
Date: 2026-01-18
"""

from typing import Dict, Any, List
import numpy as np

from .base import BaseExtractor, register_extractor
from rph_core.utils.file_io import read_xyz
from rph_core.utils.geometry_tools import GeometryUtils


def calculate_radius_of_gyration(coordinates):
    """Calculate radius of gyration of a structure.

    Args:
        coordinates: Coordinate array (N, 3)

    Returns:
        Radius of gyration in Angstroms
    """
    center_of_mass = np.mean(coordinates, axis=0)
    squared_distances = np.sum((coordinates - center_of_mass) ** 2, axis=1)
    return np.sqrt(np.mean(squared_distances))


class GeometryExtractor(BaseExtractor):
    """Extract geometric features from TS XYZ structure.

    Features (geom.* prefix):
    - natoms_ts: Number of atoms
    - r1, r2: Forming bond distances (Å)
    - asynch: |r1 - r2| (Å)
    - asynch_index: asynch / (r1 + r2)
    - rg_ts: Radius of gyration (Å)
    - min_nonbonded: Minimum non-bonded distance (Å)
    - close_contacts: Number of close contacts (< cutoff)

    Requires forming_bonds to identify forming bond atom pairs.
    """

    def get_plugin_name(self) -> str:
        return "geometry"

    def get_required_inputs(self) -> List[str]:
        return ["ts_xyz", "forming_bonds"]

    def extract(self, context) -> Dict[str, Any]:
        ts_xyz = context.ts_xyz

        if ts_xyz is None:
            raise ValueError("ts_xyz is None")

        forming_bonds = context.forming_bonds

        # Read XYZ file directly
        coords, symbols = read_xyz(ts_xyz)
        if coords is None or symbols is None:
            raise ValueError("Failed to read ts_xyz file")

        n_atoms = len(symbols)
        distance_matrix = GeometryUtils.compute_distance_matrix(coords)

        features = {}
        features["geom.natoms_ts"] = n_atoms

        r1: float = float("nan")
        r2: float = float("nan")

        # Forming bond distances
        if forming_bonds:
            bond_pairs = list(forming_bonds)

            if len(bond_pairs) >= 1:
                i, j = bond_pairs[0]
                if i < n_atoms and j < n_atoms:
                    r1 = distance_matrix[i, j]

            if len(bond_pairs) >= 2:
                i, j = bond_pairs[1]
                if i < n_atoms and j < n_atoms:
                    r2 = distance_matrix[i, j]

            features["geom.r1"] = r1
            features["geom.r2"] = r2

            # Asynchronicity
            if not np.isnan(r1) and not np.isnan(r2):
                features["geom.asynch"] = abs(r1 - r2)
                sum_r = r1 + r2
                if sum_r > 0:
                    features["geom.asynch_index"] = features["geom.asynch"] / sum_r
                else:
                    features["geom.asynch_index"] = np.nan
            elif not np.isnan(r1):
                features["geom.asynch"] = r1
                features["geom.asynch_index"] = 1.0
            elif not np.isnan(r2):
                features["geom.asynch"] = r2
                features["geom.asynch_index"] = 1.0
            else:
                features["geom.asynch"] = np.nan
                features["geom.asynch_index"] = np.nan
        else:
            features["geom.r1"] = np.nan
            features["geom.r2"] = np.nan
            features["geom.asynch"] = np.nan
            features["geom.asynch_index"] = np.nan

        features["geom.rg_ts"] = calculate_radius_of_gyration(coords)

        cutoff = context.close_contacts_cutoff
        min_nonbonded = np.inf
        close_contact_count = 0

        bonded_pairs = set()
        if forming_bonds:
            for bond in forming_bonds:
                i, j = bond
                if i < n_atoms and j < n_atoms:
                    bonded_pairs.add(tuple(sorted((i, j))))

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                pair = tuple(sorted((i, j)))
                if pair in bonded_pairs:
                    continue

                dist = distance_matrix[i, j]
                if dist < 1.25:
                    continue

                if dist < min_nonbonded:
                    min_nonbonded = dist
                if dist < cutoff:
                    close_contact_count += 1

        features["geom.min_nonbonded"] = min_nonbonded if min_nonbonded != np.inf else np.nan
        features["geom.close_contacts"] = close_contact_count

        # V6.1: Add new geometry columns for linear models
        if not np.isnan(r1) and not np.isnan(r2):
            r1_val = float(r1)
            r2_val = float(r2)
            features["geom.r_avg"] = (r1_val + r2_val) / 2.0
            features["geom.dr"] = r1_val - r2_val
        else:
            features["geom.r_avg"] = np.nan
            features["geom.dr"] = np.nan

        # V6.1: Close contacts density
        if n_atoms > 0:
            features["geom.close_contacts_density"] = close_contact_count / n_atoms
        else:
            features["geom.close_contacts_density"] = np.nan


        return features


register_extractor(GeometryExtractor())
