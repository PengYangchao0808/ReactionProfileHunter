"""Torsion-aware candidate deduplication used by the fixed CENSO-LITE funnel."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

from rph_core.steps.conformer_search.torsion_signature import TorsionSignature, build_signature, signatures_equivalent
from rph_core.utils.file_io import read_xyz


@dataclass(frozen=True)
class DedupCandidate:
    path: Path
    score: float
    signature: TorsionSignature
    metadata: Dict[str, Any]


def _heavy_atom_rmsd(mol: Chem.Mol, left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> float:
    probe = Chem.Mol(mol)
    ref = Chem.Mol(mol)
    for target, coords in ((probe, left), (ref, right)):
        conf = Chem.Conformer(target.GetNumAtoms())
        for idx, xyz in enumerate(coords):
            conf.SetAtomPosition(idx, tuple(float(value) for value in xyz))
        target.RemoveAllConformers()
        target.AddConformer(conf)
    atom_map = [(idx, idx) for idx, atom in enumerate(mol.GetAtoms()) if atom.GetAtomicNum() > 1]
    if not atom_map:
        return 0.0
    return float(rdMolAlign.AlignMol(probe, ref, atomMap=atom_map))


class TorsionAwareDeduplicator:
    """Preserve chemically distinct rotamers even when Cartesian RMSD is small."""

    def __init__(self, config: Dict[str, Any]):
        cfg = dict(config or {})
        self.bin_width_deg = float(cfg.get("torsion_bin_deg", 20.0))
        self.torsion_tolerance_deg = float(cfg.get("torsion_rmsd_deg", 25.0))
        self.rmsd_prefilter = float(cfg.get("heavy_atom_rmsd_prefilter_A", 0.25))

    def deduplicate(self, mol: Chem.Mol, candidates: Iterable[DedupCandidate]) -> List[DedupCandidate]:
        kept: List[DedupCandidate] = []
        for candidate in sorted(candidates, key=lambda item: (float(item.score), str(item.path))):
            duplicate = False
            for existing in kept:
                if not signatures_equivalent(candidate.signature, existing.signature, self.torsion_tolerance_deg):
                    continue
                try:
                    _, left = read_xyz(candidate.path)
                    _, right = read_xyz(existing.path)
                    if _heavy_atom_rmsd(mol, left, right) <= self.rmsd_prefilter:
                        duplicate = True
                        break
                except (OSError, ValueError, IndexError):
                    continue
            if not duplicate:
                kept.append(candidate)
        return kept

    def annotate(self, mol: Chem.Mol, path: Path, score: float, metadata: Dict[str, Any] | None = None) -> DedupCandidate:
        coordinates, _ = read_xyz(path)
        signature = build_signature(mol, coordinates, self.bin_width_deg)
        payload = dict(metadata or {})
        payload["torsion_signature"] = signature.key()
        return DedupCandidate(Path(path), float(score), signature, payload)
