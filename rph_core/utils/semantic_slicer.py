"""
SemanticSlicer: ASM 2.0 Semantic-Aware Fragmentation
====================================================

Houk Lab ASM 2.0 implementation for robust fragment splitting in
intramolecular [5+2] cycloaddition systems.

Features:
- No bond order dependence (Connectivity only)
- Strict BFS pathfinding (No back-flow)
- Robust Core Identification (Scoring + Strict Expansion)
- Anti-collision H-Capping vectors

Author: Houk Lab / QCcalc Team
Date: 2026-01-31
"""

from __future__ import annotations
from typing import List, Tuple, Optional, Set, Dict, Any, TYPE_CHECKING
import logging

import numpy as np

RDKIT_AVAILABLE: bool = False

try:
    from rdkit import Chem
    from rdkit.Chem import rdmolops
    RDKIT_AVAILABLE = True
except ImportError:
    pass

if TYPE_CHECKING:
    from rdkit import Chem
    from rdkit.Chem import rdmolops

logger = logging.getLogger(__name__)


class SemanticSlicerResult:
    """Result container for SemanticSlicer operations."""
    
    def __init__(
        self,
        success: bool,
        cut_bond: Optional[Tuple[int, int]] = None,
        core_indices: Optional[List[int]] = None,
        target_indices: Optional[List[int]] = None,
        h_cap_vectors: Optional[Tuple[np.ndarray, np.ndarray]] = None,
        debug_info: Optional[Dict[str, Any]] = None,
        reason: str = ""
    ):
        self.success = success
        self.cut_bond = cut_bond
        self.core_indices = core_indices or []
        self.target_indices = target_indices or []
        self.h_cap_vectors = h_cap_vectors
        self.debug_info = debug_info or {}
        self.reason = reason


class SemanticSlicer:
    """
    Houk Lab ASM 2.0: Semantic-Aware Fragmentation (Production Ready)
    
    For intramolecular [5+2] cycloaddition systems, identifies the optimal
    tether cut bond using ring-based core detection rather than geometric
    shortest paths.
    
    Features:
    - No bond order dependence (Connectivity only)
    - Strict BFS pathfinding (No back-flow)
    - Robust Core Identification (Scoring + Strict Expansion)
    - Anti-collision H-Capping
    """
    
    # Covalent radii for bond detection (Angstrom)
    COV_RADII = {
        'C': 0.77, 'H': 0.38, 'O': 0.73, 'N': 0.75, 
        'S': 1.02, 'P': 1.06, 'F': 0.64, 'Cl': 0.99,
        'Br': 1.14, 'I': 1.33
    }
    
    def __init__(
        self,
        atoms: List[str],
        coords: np.ndarray,
        forming_bonds: Tuple[Tuple[int, int], Tuple[int, int]],
        verbose: bool = True
    ):
        """
        Initialize SemanticSlicer.
        
        Args:
            atoms: List of element symbols
            coords: Atomic coordinates (N, 3) in Angstrom
            forming_bonds: ((u1, v1), (u2, v2)) forming bond atom pairs
            verbose: Enable debug logging
        """
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit is required for SemanticSlicer")
            
        self.atoms = atoms
        self.coords = np.array(coords)
        self.verbose = verbose
        
        self.forming_pairs: Set[Tuple[int, int]] = set()
        for b in forming_bonds:
            self.forming_pairs.add((min(b[0], b[1]), max(b[0], b[1])))
        
        self.forming_atoms_flat: Set[int] = set()
        for b in forming_bonds:
            self.forming_atoms_flat.add(b[0])
            self.forming_atoms_flat.add(b[1])
            
        self.mol: Any = None
        self.debug_info: Dict[str, Any] = {}

    def build_prereaction_graph(self) -> None:
        """
        Build pre-reaction molecular graph excluding forming bonds.
        
        Uses RWMol for explicit bond control and applies strict
        geometric bond detection with forming bond exclusion.
        """
        mol = Chem.RWMol()  # pyright: ignore
        
        for i, symbol in enumerate(self.atoms):
            a = Chem.Atom(symbol)  # pyright: ignore
            a.SetIntProp('original_idx', i)
            mol.AddAtom(a)

        n = len(self.atoms)
        
        for i in range(n):
            for j in range(i + 1, n):
                if tuple(sorted((i, j))) in self.forming_pairs:
                    continue
                
                dist = np.linalg.norm(self.coords[i] - self.coords[j])
                r_sum = (
                    self.COV_RADII.get(self.atoms[i], 0.8) + 
                    self.COV_RADII.get(self.atoms[j], 0.8)
                )
                
                is_bond = False
                if dist < (r_sum + 0.30):
                    is_bond = True
                elif {self.atoms[i], self.atoms[j]} == {'C', 'O'} and dist < 1.35:
                    is_bond = True

                if is_bond:
                    if mol.GetBondBetweenAtoms(i, j) is None:
                        mol.AddBond(i, j, Chem.BondType.SINGLE)  # type: ignore[union-attr]

        self.mol = mol.GetMol()
        
        try:
            rdmolops.GetSymmSSSR(self.mol)  # type: ignore[union-attr]
        except Exception:
            if self.verbose:
                logger.warning("SSSR calculation failed on pre-reaction graph")

    def _score_ring_candidate(self, ring_idxs: Tuple[int, ...]) -> int:
        """
        Score a ring for core identification.
        
        Scoring rules:
        - +2: Contains heteroatom (O, N)
        - +3: Has carbonyl-like feature (C-O < 1.32 A)
        - +2: Overlaps with reaction center
        
        Args:
            ring_idxs: Atom indices in the ring
            
        Returns:
            Integer score (higher = better core candidate)
        """
        score = 0
        ring_atoms = [self.atoms[idx] for idx in ring_idxs]
        
        # Rule 1: Heteroatom presence
        if 'O' in ring_atoms or 'N' in ring_atoms:
            score += 2
            
        # Rule 2: Carbonyl-like feature (Internal or External)
        has_carbonyl = False
        for r_idx in ring_idxs:
            if self.atoms[r_idx] == 'C':
                atom_obj = self.mol.GetAtomWithIdx(r_idx)
                for nbr in atom_obj.GetNeighbors():
                    n_idx = nbr.GetIdx()
                    if self.atoms[n_idx] == 'O':
                        d = np.linalg.norm(self.coords[r_idx] - self.coords[n_idx])
                        if d < 1.32:
                            has_carbonyl = True
                            break
            if has_carbonyl:
                break
            
        if has_carbonyl:
            score += 3
            
        # Rule 3: Proximity to reaction center
        overlap = set(ring_idxs).intersection(self.forming_atoms_flat)
        if overlap:
            score += 2

        return score

    def identify_core(self) -> Optional[Set[int]]:
        """
        Identify the dipole core using ring detection and scoring.
        
        Returns:
            Set of atom indices belonging to the core, or None if failed
        """
        if self.mol is None:
            self.build_prereaction_graph()
            
        ri = self.mol.GetRingInfo()
        atom_rings = ri.AtomRings()
        
        if not atom_rings:
            if self.verbose:
                logger.debug("No rings found in pre-reaction graph")
            return None

        best_ring: Optional[Set[int]] = None
        max_score = -1
        
        for ring in atom_rings:
            if len(ring) not in [5, 6, 7]:
                continue
            score = self._score_ring_candidate(ring)
            if score > max_score:
                max_score = score
                best_ring = set(ring)
        
        if max_score < 3:  # Minimum threshold
            if self.verbose:
                logger.debug(f"Core score too low: {max_score}")
            return None

        if best_ring is None:
            return None
            
        final_core = set(best_ring)
        for r_idx in best_ring:
            atom_sym = self.atoms[r_idx]
            if atom_sym != 'C':
                continue  # Only expand from C
            
            atom_obj = self.mol.GetAtomWithIdx(r_idx)
            for nbr in atom_obj.GetNeighbors():
                n_idx = nbr.GetIdx()
                if n_idx in final_core:
                    continue
                
                n_sym = self.atoms[n_idx]
                # Strict Element Pair + Distance
                if n_sym in ['O', 'N']:
                    d = np.linalg.norm(self.coords[r_idx] - self.coords[n_idx])
                    if d < 1.32:
                        final_core.add(n_idx)
                        
        self.debug_info['core_local_indices'] = list(final_core)
        return final_core

    def _resolve_targets(self, core_indices: Set[int]) -> Optional[Set[int]]:
        """
        Resolve target atoms (alkene end) from forming atoms not in core.
        
        Ensures exactly 2 atoms are identified as targets.
        
        Args:
            core_indices: Set of core atom indices
            
        Returns:
            Set of 2 target atom indices, or None if resolution fails
        """
        # Raw difference: forming atoms NOT in core
        raw_targets = list(self.forming_atoms_flat - core_indices)
        
        if len(raw_targets) == 2:
            return set(raw_targets)
        
        if self.verbose:
            logger.debug(f"Ambiguous targets found: {raw_targets}. Attempting resolution.")
            
        # Strategy 1: Find bonded pair within raw_targets
        pairs = []
        n_raw = len(raw_targets)
        for i in range(n_raw):
            for j in range(i + 1, n_raw):
                pairs.append((raw_targets[i], raw_targets[j]))
        
        # Check if bonded in graph
        bonded_pairs = []
        for u, v in pairs:
            if self.mol.GetBondBetweenAtoms(u, v):
                bonded_pairs.append((u, v))
                
        if len(bonded_pairs) == 1:
            return set(bonded_pairs[0])
        elif len(bonded_pairs) > 1:
            # Pick shortest bond
            best_pair = min(
                bonded_pairs, 
                key=lambda p: np.linalg.norm(self.coords[p[0]] - self.coords[p[1]])
            )
            return set(best_pair)
            
        # Strategy 2: Geometric fallback (Shortest distance among all pairs)
        if pairs:
            best_pair = min(
                pairs, 
                key=lambda p: np.linalg.norm(self.coords[p[0]] - self.coords[p[1]])
            )
            return set(best_pair)
            
        return None  # Failed to resolve

    def _bfs_check(
        self,
        start_node: int,
        target_nodes: Set[int],
        forbidden_nodes: Set[int],
        blocked_edge: Tuple[int, int]
    ) -> Tuple[bool, int]:
        """
        Strict BFS to check reachability from start to targets.
        
        Args:
            start_node: Starting atom index
            target_nodes: Set of target atom indices
            forbidden_nodes: Set of atoms to avoid (core)
            blocked_edge: (u, v) edge to block
            
        Returns:
            (can_reach, steps) tuple
        """
        if start_node in target_nodes:
            return True, 0
        
        queue = [(start_node, 0)]
        visited = {start_node}
        visited.update(forbidden_nodes)
        
        blocked_edge_sorted = tuple(sorted(blocked_edge))
        
        while queue:
            curr, steps = queue.pop(0)
            
            if curr in target_nodes:
                return True, steps
            
            atom = self.mol.GetAtomWithIdx(curr)
            for nbr in atom.GetNeighbors():
                n_idx = nbr.GetIdx()
                
                # Check blocked edge
                if tuple(sorted((curr, n_idx))) == blocked_edge_sorted:
                    continue
                
                if n_idx not in visited:
                    visited.add(n_idx)
                    queue.append((n_idx, steps + 1))
                    
        return False, 999

    def find_cut_bond(self) -> Optional[Tuple[int, int]]:
        """
        Find the optimal tether-core cut bond.
        
        Main entry point for semantic slicing. Uses ring-based core
        identification and BFS reachability to find the cut bond.
        
        Returns:
            (u, v) tuple of original atom indices, or None if failed
        """
        core_indices = self.identify_core()
        if not core_indices:
            self.debug_info['failure_reason'] = 'core_identification_failed'
            return None

        targets = self._resolve_targets(core_indices)
        if not targets:
            self.debug_info['failure_reason'] = 'target_resolution_failed'
            return None
            
        self.debug_info['targets_local'] = list(targets)

        candidates = []
        
        # Iterate boundary bonds (one atom in core, one outside)
        for bond in self.mol.GetBonds():
            u = bond.GetBeginAtomIdx()
            v = bond.GetEndAtomIdx()
            
            u_in = u in core_indices
            v_in = v in core_indices
            
            if u_in != v_in:  # Boundary bond
                core_atom = u if u_in else v
                tether_atom = v if u_in else u
                
                # Check if tether side can reach targets without going through core
                can_reach, steps = self._bfs_check(
                    start_node=tether_atom, 
                    target_nodes=targets, 
                    forbidden_nodes=core_indices,
                    blocked_edge=(core_atom, tether_atom)
                )
                
                if can_reach:
                    # Calculate geom distance to targets for tie-breaking
                    dists = [
                        np.linalg.norm(self.coords[tether_atom] - self.coords[t]) 
                        for t in targets
                    ]
                    min_geom_dist = min(dists)
                    
                    candidates.append({
                        'bond': (core_atom, tether_atom),
                        'steps': steps,
                        'dist': min_geom_dist
                    })

        if not candidates:
            self.debug_info['failure_reason'] = 'no_valid_cut_candidates'
            return None
            
        # Sort by Steps (primary) and Geometry (secondary)
        candidates.sort(key=lambda x: (x['steps'], x['dist']))
        best_cand = candidates[0]
        
        # Map back to Original Indices
        u_local, v_local = best_cand['bond']
        u_real = self.mol.GetAtomWithIdx(u_local).GetIntProp('original_idx')
        v_real = self.mol.GetAtomWithIdx(v_local).GetIntProp('original_idx')
        
        self.debug_info['chosen_cut'] = (u_real, v_real)
        self.debug_info['all_candidates'] = [
            {'bond': c['bond'], 'steps': c['steps'], 'dist': round(c['dist'], 3)}
            for c in candidates[:5]  # Top 5 for debugging
        ]
        
        return (u_real, v_real)

    def get_h_cap_vectors(
        self, 
        u_real: int, 
        v_real: int
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate H-cap placement vectors for cut bond atoms.
        
        Vectors point AWAY from the bond (correct direction for H placement).
        
        Args:
            u_real: First atom index of cut bond
            v_real: Second atom index of cut bond
            
        Returns:
            (vec_for_u, vec_for_v) normalized direction vectors
        """
        coord_u = self.coords[u_real]
        coord_v = self.coords[v_real]
        
        # H on u points away from v
        vec_u = coord_u - coord_v
        # H on v points away from u
        vec_v = coord_v - coord_u
        
        norm_u = np.linalg.norm(vec_u)
        norm_v = np.linalg.norm(vec_v)
        
        if norm_u < 1e-8 or norm_v < 1e-8:
            raise ValueError("Cut bond atoms are coincident")
        
        return (vec_u / norm_u), (vec_v / norm_v)

    def run(self) -> SemanticSlicerResult:
        """
        Execute full semantic slicing workflow.
        
        Returns:
            SemanticSlicerResult with cut bond and metadata
        """
        try:
            cut_bond = self.find_cut_bond()
            
            if cut_bond is None:
                return SemanticSlicerResult(
                    success=False,
                    debug_info=self.debug_info,
                    reason=self.debug_info.get('failure_reason', 'unknown')
                )
            
            u, v = cut_bond
            h_cap_vecs = self.get_h_cap_vectors(u, v)
            
            return SemanticSlicerResult(
                success=True,
                cut_bond=cut_bond,
                core_indices=self.debug_info.get('core_local_indices', []),
                target_indices=self.debug_info.get('targets_local', []),
                h_cap_vectors=h_cap_vecs,
                debug_info=self.debug_info,
                reason=""
            )
            
        except Exception as e:
            logger.error(f"SemanticSlicer failed: {e}", exc_info=True)
            return SemanticSlicerResult(
                success=False,
                debug_info=self.debug_info,
                reason=str(e)
            )


def semantic_slice(
    atoms: List[str],
    coords: np.ndarray,
    forming_bonds: Tuple[Tuple[int, int], Tuple[int, int]],
    verbose: bool = False
) -> SemanticSlicerResult:
    """
    Convenience function for semantic slicing.
    
    Args:
        atoms: List of element symbols
        coords: Atomic coordinates (N, 3)
        forming_bonds: ((u1, v1), (u2, v2)) forming bond pairs
        verbose: Enable debug logging
        
    Returns:
        SemanticSlicerResult
    """
    if not RDKIT_AVAILABLE:
        return SemanticSlicerResult(
            success=False,
            reason="RDKit not available"
        )
    
    slicer = SemanticSlicer(atoms, coords, forming_bonds, verbose=verbose)
    return slicer.run()
