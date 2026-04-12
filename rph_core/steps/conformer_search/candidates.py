from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from rph_core.utils.constants import HARTREE_TO_KCAL


@dataclass
class ConformerCandidate:
    candidate_id: str
    source_xyz: Path
    source_stage: str
    rank_initial: int
    xtb_energy: Optional[float]
    xtb_free_energy: Optional[float]
    prescreen_energy: Optional[float]
    screening_energy: Optional[float]
    rerank_energy: Optional[float]
    rerank_method: Optional[str]
    ranking_history: Dict[str, float]
    included_by: List[str]
    within_window: bool
    cluster_id: Optional[int]
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class CandidateSet:
    candidates: List[ConformerCandidate]
    source_protocol: str
    ranking_basis: str
    ranking_unit: str
    total_count: int
    selected_count_before_handoff: int
    metadata: Dict[str, Any] = field(default_factory=dict)


def candidate_set_from_paths(
    paths: Sequence[Path],
    *,
    source_protocol: str,
    ranking_basis: str = "input_order",
    ranking_unit: str = "hartree",
    source_stage: str = "ensemble",
) -> CandidateSet:
    candidates: List[ConformerCandidate] = []
    for idx, path in enumerate(paths):
        resolved = Path(path).resolve()
        candidates.append(
            ConformerCandidate(
                candidate_id=f"cand_{idx:03d}",
                source_xyz=resolved,
                source_stage=source_stage,
                rank_initial=idx,
                xtb_energy=None,
                xtb_free_energy=None,
                prescreen_energy=None,
                screening_energy=None,
                rerank_energy=None,
                rerank_method=None,
                ranking_history={},
                included_by=["input"],
                within_window=True,
                cluster_id=None,
                metadata={},
            )
        )
    return CandidateSet(
        candidates=candidates,
        source_protocol=source_protocol,
        ranking_basis=ranking_basis,
        ranking_unit=ranking_unit,
        total_count=len(candidates),
        selected_count_before_handoff=len(candidates),
        metadata={},
    )


def clone_candidate_set(candidate_set: CandidateSet, candidates: List[ConformerCandidate]) -> CandidateSet:
    return CandidateSet(
        candidates=candidates,
        source_protocol=candidate_set.source_protocol,
        ranking_basis=candidate_set.ranking_basis,
        ranking_unit=candidate_set.ranking_unit,
        total_count=candidate_set.total_count,
        selected_count_before_handoff=len(candidates),
        metadata=dict(candidate_set.metadata),
    )


def sort_candidates(candidate_set: CandidateSet, basis: str) -> CandidateSet:
    key_map = {
        "rerank_energy": lambda c: c.rerank_energy,
        "screening_energy": lambda c: c.screening_energy,
        "prescreen_energy": lambda c: c.prescreen_energy,
        "xtb_free_energy": lambda c: c.xtb_free_energy,
        "xtb_energy": lambda c: c.xtb_energy,
    }
    extractor = key_map.get(basis, key_map["rerank_energy"])
    ordered = sorted(
        candidate_set.candidates,
        key=lambda c: (float(extractor(c)) if extractor(c) is not None else float(c.rank_initial), c.rank_initial),
    )
    ranked: List[ConformerCandidate] = []
    for idx, candidate in enumerate(ordered):
        updated = ConformerCandidate(
            candidate_id=candidate.candidate_id,
            source_xyz=candidate.source_xyz,
            source_stage=candidate.source_stage,
            rank_initial=idx,
            xtb_energy=candidate.xtb_energy,
            xtb_free_energy=candidate.xtb_free_energy,
            prescreen_energy=candidate.prescreen_energy,
            screening_energy=candidate.screening_energy,
            rerank_energy=candidate.rerank_energy,
            rerank_method=candidate.rerank_method,
            ranking_history=dict(candidate.ranking_history),
            included_by=list(candidate.included_by),
            within_window=candidate.within_window,
            cluster_id=candidate.cluster_id,
            metadata=dict(candidate.metadata),
        )
        ranked.append(updated)
    return CandidateSet(
        candidates=ranked,
        source_protocol=candidate_set.source_protocol,
        ranking_basis=basis,
        ranking_unit=candidate_set.ranking_unit,
        total_count=candidate_set.total_count,
        selected_count_before_handoff=len(ranked),
        metadata=dict(candidate_set.metadata),
    )


def mark_candidates_within_window(candidate_set: CandidateSet, basis: str, window_kcal: float) -> CandidateSet:
    sorted_set = sort_candidates(candidate_set, basis)
    if not sorted_set.candidates:
        return sorted_set
    threshold = float(window_kcal)
    if str(sorted_set.ranking_unit).lower() in {"hartree", "ha"}:
        threshold = float(window_kcal) / float(HARTREE_TO_KCAL)
    energy_values: List[float] = []
    for candidate in sorted_set.candidates:
        value = candidate.rerank_energy
        if basis == "screening_energy":
            value = candidate.screening_energy
        elif basis == "prescreen_energy":
            value = candidate.prescreen_energy
        elif basis == "xtb_free_energy":
            value = candidate.xtb_free_energy
        elif basis == "xtb_energy":
            value = candidate.xtb_energy
        energy_values.append(float(value) if value is not None else float(candidate.rank_initial))
    anchor = energy_values[0]
    marked: List[ConformerCandidate] = []
    for candidate, energy in zip(sorted_set.candidates, energy_values):
        within = (energy - anchor) <= threshold
        included_by = list(candidate.included_by)
        if within and "window" not in included_by:
            included_by.append("window")
        marked.append(
            ConformerCandidate(
                candidate_id=candidate.candidate_id,
                source_xyz=candidate.source_xyz,
                source_stage=candidate.source_stage,
                rank_initial=candidate.rank_initial,
                xtb_energy=candidate.xtb_energy,
                xtb_free_energy=candidate.xtb_free_energy,
                prescreen_energy=candidate.prescreen_energy,
                screening_energy=candidate.screening_energy,
                rerank_energy=candidate.rerank_energy,
                rerank_method=candidate.rerank_method,
                ranking_history=dict(candidate.ranking_history),
                included_by=included_by,
                within_window=within,
                cluster_id=candidate.cluster_id,
                metadata=dict(candidate.metadata),
            )
        )
    return CandidateSet(
        candidates=marked,
        source_protocol=sorted_set.source_protocol,
        ranking_basis=sorted_set.ranking_basis,
        ranking_unit=sorted_set.ranking_unit,
        total_count=sorted_set.total_count,
        selected_count_before_handoff=sum(1 for c in marked if c.within_window),
        metadata=dict(sorted_set.metadata),
    )


def select_candidates_by_window(candidate_set: CandidateSet, basis: str, window_kcal: float) -> CandidateSet:
    marked = mark_candidates_within_window(candidate_set, basis, window_kcal)
    selected = [candidate for candidate in marked.candidates if candidate.within_window]
    if not selected and marked.candidates:
        selected = [marked.candidates[0]]
    return clone_candidate_set(marked, selected)
