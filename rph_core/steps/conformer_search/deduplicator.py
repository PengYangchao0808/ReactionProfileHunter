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
            duplicate_index = None
            for index, existing in enumerate(kept):
                if not signatures_equivalent(candidate.signature, existing.signature, self.torsion_tolerance_deg):
                    continue
                try:
                    left, _ = read_xyz(candidate.path)
                    right, _ = read_xyz(existing.path)
                    if _heavy_atom_rmsd(mol, left, right) <= self.rmsd_prefilter:
                        duplicate_index = index
                        break
                except (OSError, ValueError, IndexError):
                    continue
            if duplicate_index is None:
                kept.append(self._with_merge_provenance(candidate))
                continue

            # CREST sampling frequency is search provenance, not a physical
            # statistical degeneracy. Merge the provenance while deliberately
            # retaining the representative's independently justified d_i.
            existing = kept[duplicate_index]
            metadata = dict(existing.metadata)
            merged_from = list(metadata.get("merged_from") or self._source_ids(existing))
            for source_id in self._source_ids(candidate):
                if source_id not in merged_from:
                    merged_from.append(source_id)
            metadata.update(
                {
                    "merged_from": merged_from,
                    "merge_count": len(merged_from),
                    "degeneracy": int(metadata.get("degeneracy", 1)),
                    "degeneracy_source": metadata.get(
                        "degeneracy_source", "default_unique_minimum"
                    ),
                }
            )
            kept[duplicate_index] = DedupCandidate(
                existing.path, existing.score, existing.signature, metadata
            )
        return kept

    @staticmethod
    def _source_ids(candidate: DedupCandidate) -> List[str]:
        values = candidate.metadata.get("merged_from")
        if isinstance(values, Sequence) and not isinstance(values, (str, bytes)) and values:
            return [str(value) for value in values]
        return [str(candidate.metadata.get("source_conformer_id") or candidate.path.stem)]

    def _with_merge_provenance(self, candidate: DedupCandidate) -> DedupCandidate:
        metadata = dict(candidate.metadata)
        merged_from = self._source_ids(candidate)
        metadata.update(
            {
                "source_conformer_id": str(
                    metadata.get("source_conformer_id") or candidate.path.stem
                ),
                "merged_from": merged_from,
                "merge_count": len(merged_from),
                "degeneracy": int(metadata.get("degeneracy", 1)),
                "degeneracy_source": metadata.get(
                    "degeneracy_source", "default_unique_minimum"
                ),
            }
        )
        return DedupCandidate(candidate.path, candidate.score, candidate.signature, metadata)

    def annotate(self, mol: Chem.Mol, path: Path, score: float, metadata: Dict[str, Any] | None = None) -> DedupCandidate:
        coordinates, _ = read_xyz(path)
        signature = build_signature(mol, coordinates, self.bin_width_deg)
        payload = dict(metadata or {})
        payload["torsion_signature"] = signature.key()
        payload.setdefault("source_conformer_id", Path(path).stem)
        payload.setdefault("merged_from", [str(payload["source_conformer_id"])])
        payload.setdefault("merge_count", len(payload["merged_from"]))
        payload.setdefault("degeneracy", 1)
        payload.setdefault("degeneracy_source", "default_unique_minimum")
        return DedupCandidate(Path(path), float(score), signature, payload)
