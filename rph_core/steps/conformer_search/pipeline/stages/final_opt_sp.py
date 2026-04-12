from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

from rph_core.steps.conformer_search.candidates import CandidateSet, ConformerCandidate, candidate_set_from_paths
from rph_core.steps.conformer_search.protocols import ProtocolSpec
from rph_core.utils.constants import HARTREE_TO_KCAL


@dataclass(frozen=True)
class FinalOptSPResult:
    best_log: Path
    best_result_path: Path
    best_energy: float
    stage_meta: Dict[str, Any]


def _normalize_candidates(candidates: Sequence[Path]) -> List[Path]:
    normalized: List[Path] = []
    for item in candidates:
        normalized.append(Path(item).resolve())
    return normalized


def _candidate_energy(candidate: ConformerCandidate) -> Optional[float]:
    if candidate.rerank_energy is not None:
        return float(candidate.rerank_energy)
    if candidate.xtb_free_energy is not None:
        return float(candidate.xtb_free_energy)
    if candidate.screening_energy is not None:
        return float(candidate.screening_energy)
    if candidate.prescreen_energy is not None:
        return float(candidate.prescreen_energy)
    if candidate.xtb_energy is not None:
        return float(candidate.xtb_energy)
    return None


def _delta_to_kcal(delta: float, ranking_unit: str) -> float:
    if str(ranking_unit).lower() in {"hartree", "ha"}:
        return float(delta) * float(HARTREE_TO_KCAL)
    return float(delta)


def _sort_candidates_for_selection(candidate_set: CandidateSet) -> List[ConformerCandidate]:
    def _energy_key(candidate: ConformerCandidate) -> Tuple[float, int]:
        energy = _candidate_energy(candidate)
        return (float(energy) if energy is not None else float(candidate.rank_initial), candidate.rank_initial)

    return sorted(candidate_set.candidates, key=_energy_key)


def _select_all_candidates(candidate_set: CandidateSet) -> List[ConformerCandidate]:
    return _sort_candidates_for_selection(candidate_set)


def _select_survivors_within_window(candidate_set: CandidateSet, window_kcal: float) -> List[ConformerCandidate]:
    ordered = _sort_candidates_for_selection(candidate_set)
    if not ordered:
        return []

    reference = None
    values: List[Tuple[ConformerCandidate, float]] = []
    for candidate in ordered:
        val = _candidate_energy(candidate)
        if val is None:
            continue
        value = float(val)
        if reference is None:
            reference = value
        values.append((candidate, value))

    if not values:
        return ordered

    threshold = float(window_kcal)
    if str(candidate_set.ranking_unit).lower() in {"hartree", "ha"}:
        threshold = float(window_kcal) / float(HARTREE_TO_KCAL)
    selected = [cand for cand, value in values if (value - values[0][1]) <= threshold]
    return selected if selected else [values[0][0]]


def _select_rank1(candidate_set: CandidateSet) -> List[ConformerCandidate]:
    ordered = _sort_candidates_for_selection(candidate_set)
    return ordered[:1]


def _select_top2_if_gap_small(candidate_set: CandidateSet, gap_kcal: float) -> Tuple[List[ConformerCandidate], bool]:
    ordered = _sort_candidates_for_selection(candidate_set)
    if len(ordered) < 2:
        return ordered[:1], False

    e1 = _candidate_energy(ordered[0])
    e2 = _candidate_energy(ordered[1])
    if e1 is None or e2 is None:
        return ordered[:1], False

    gap_value = _delta_to_kcal(e2 - e1, candidate_set.ranking_unit)
    if gap_value <= float(gap_kcal):
        return ordered[:2], True
    return ordered[:1], False


def _select_all_within_window(candidate_set: CandidateSet, window_kcal: float) -> List[ConformerCandidate]:
    return _select_survivors_within_window(candidate_set, window_kcal)


def _select_candidates(candidate_set: CandidateSet, protocol_spec: Optional[ProtocolSpec]) -> Tuple[List[ConformerCandidate], Dict[str, Any]]:
    if protocol_spec is None:
        selected = _select_all_candidates(candidate_set)
        return selected, {
            "handoff_mode": "optimize_all_candidates",
            "fallback_mode": None,
            "fallback_triggered": False,
            "energy_window_used": None,
            "small_gap_kcal": None,
            "selection_reason": "all_candidates_default",
        }

    mode = protocol_spec.handoff_policy.mode
    mode_requested = mode
    mode_effective = mode
    fallback_mode = protocol_spec.handoff_policy.fallback_mode
    fallback_triggered = False
    energy_window_used = None
    small_gap_kcal = protocol_spec.handoff_policy.small_gap_kcal
    ordered = _sort_candidates_for_selection(candidate_set)
    window_count = len(candidate_set.candidates)
    gap_rank1_rank2 = None
    if len(ordered) >= 2:
        e1 = _candidate_energy(ordered[0])
        e2 = _candidate_energy(ordered[1])
        if e1 is not None and e2 is not None:
            gap_rank1_rank2 = _delta_to_kcal(e2 - e1, candidate_set.ranking_unit)

    if mode == "optimize_all_candidates":
        selected = _select_all_candidates(candidate_set)
        selection_reason = "all_candidates"
        mode_effective = "optimize_all_candidates"
    elif mode == "optimize_all_survivors_within_window":
        energy_window_used = protocol_spec.funnel_policy.survivor_window_kcal or 3.0
        selected = _select_survivors_within_window(candidate_set, energy_window_used)
        selection_reason = "survivors_within_window"
        mode_effective = "optimize_all_survivors_within_window"
    elif mode == "optimize_rank1":
        selected = _select_rank1(candidate_set)
        selection_reason = "rank1"
        mode_effective = "optimize_rank1"
        if fallback_mode == "optimize_top2_if_gap_small":
            threshold = protocol_spec.handoff_policy.small_gap_kcal or 1.0
            fallback_candidates, triggered = _select_top2_if_gap_small(candidate_set, threshold)
            if triggered:
                selected = fallback_candidates
                fallback_triggered = True
                small_gap_kcal = threshold
                selection_reason = "top2_fallback_gap_small"
                mode_effective = "optimize_top2_if_gap_small"
        elif fallback_mode == "optimize_all_within_0p5_kcal":
            selected = _select_all_within_window(candidate_set, 0.5)
            fallback_triggered = len(selected) > 1
            energy_window_used = 0.5
            selection_reason = "within_0p5_kcal_window"
            mode_effective = "optimize_all_within_0p5_kcal" if len(selected) > 1 else "optimize_rank1"
    elif mode == "optimize_top2_if_gap_small":
        threshold = protocol_spec.handoff_policy.small_gap_kcal or 1.0
        selected, fallback_triggered = _select_top2_if_gap_small(candidate_set, threshold)
        small_gap_kcal = threshold
        selection_reason = "explicit_top2_if_gap_small"
        mode_effective = "optimize_top2_if_gap_small" if fallback_triggered else "optimize_rank1"
    elif mode == "optimize_all_within_0p5_kcal":
        selected = _select_all_within_window(candidate_set, 0.5)
        energy_window_used = 0.5
        fallback_triggered = len(selected) > 1
        selection_reason = "explicit_within_0p5_kcal"
        mode_effective = "optimize_all_within_0p5_kcal" if len(selected) > 1 else "optimize_rank1"
    else:
        selected = _select_all_candidates(candidate_set)
        selection_reason = "fallback_all_candidates"
        mode_effective = "optimize_all_candidates"

    meta = {
        "handoff_mode": mode,
        "mode_requested": mode_requested,
        "mode_effective": mode_effective,
        "fallback_mode": fallback_mode,
        "fallback_triggered": fallback_triggered,
        "fallback_trigger": fallback_triggered,
        "energy_window_used": energy_window_used,
        "small_gap_kcal": small_gap_kcal,
        "selection_reason": selection_reason,
        "window_count": window_count,
        "gap_rank1_rank2": gap_rank1_rank2,
    }
    return selected, meta


def execute_final_opt_sp(
    *,
    candidates: Optional[List[Path]] = None,
    candidate_set: Optional[CandidateSet] = None,
    runner: Callable[[List[ConformerCandidate]], Union[Tuple[Path, float], Tuple[Path, float, Dict[str, Any]]]],
    protocol_spec: Optional[ProtocolSpec] = None,
) -> FinalOptSPResult:
    if candidate_set is None:
        normalized_candidates = _normalize_candidates(candidates or [])
        candidate_set = candidate_set_from_paths(
            normalized_candidates,
            source_protocol=protocol_spec.name if protocol_spec else "ext",
        )

    if not candidate_set.candidates:
        raise ValueError("final_opt_sp requires at least one candidate")

    protocol = protocol_spec.name if protocol_spec else "ext"
    requested_selection_mode = protocol_spec.handoff_policy.mode if protocol_spec else "optimize_all_candidates"
    selected_candidates, selection_meta = _select_candidates(candidate_set, protocol_spec)

    freq_requested = bool(protocol_spec.freq_enabled) if protocol_spec else False
    final_sp_requested = bool(protocol_spec.final_sp_enabled) if protocol_spec else True
    runner_result = runner(selected_candidates)
    if len(runner_result) == 2:
        best_log, best_energy = runner_result
        runner_meta: Dict[str, Any] = {}
    else:
        best_log, best_energy, runner_meta = runner_result
    best_result_path = Path(best_log).resolve()

    stage_meta: Dict[str, Any] = {
        "stage": "final_opt_sp",
        "stage_status": "completed",
        "protocol": protocol,
        "selection_mode": requested_selection_mode,
        "input_candidate_count": len(candidate_set.candidates),
        "selected_candidate_count": len(selected_candidates),
        "input_candidate_ids": [candidate.candidate_id for candidate in candidate_set.candidates],
        "selected_candidate_ids": [candidate.candidate_id for candidate in selected_candidates],
        "ranking_basis_before_handoff": candidate_set.ranking_basis,
        "candidate_scores_before_handoff": {
            candidate.candidate_id: (
                candidate.rerank_energy
                if candidate.rerank_energy is not None
                else (candidate.xtb_free_energy if candidate.xtb_free_energy is not None else candidate.xtb_energy)
            )
            for candidate in candidate_set.candidates
        },
        "freq_requested": freq_requested,
        "final_sp_requested": final_sp_requested,
        **selection_meta,
        **runner_meta,
    }

    return FinalOptSPResult(
        best_log=best_result_path,
        best_result_path=best_result_path,
        best_energy=float(best_energy),
        stage_meta=stage_meta,
    )
