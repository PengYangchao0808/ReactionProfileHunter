import math
from pathlib import Path
from typing import Any, List, Optional

from rph_core.steps.conformer_search.candidates import (
    CandidateSet,
    ConformerCandidate,
    candidate_set_from_paths,
    clone_candidate_set,
    select_candidates_by_window,
    sort_candidates,
)
from rph_core.steps.conformer_search.protocols import ProtocolSpec, is_ext_protocol, is_full_protocol, is_lite_protocol, is_zero_protocol
from rph_core.utils.constants import HARTREE_TO_KCAL


class FunnelRunner:
    def __init__(self, engine: Any, protocol_spec: ProtocolSpec):
        self.engine = engine
        self.protocol_spec = protocol_spec

    def run(self, initial_xyz: Path) -> CandidateSet:
        if is_ext_protocol(self.protocol_spec):
            return self.run_ext(initial_xyz)
        if is_full_protocol(self.protocol_spec):
            return self.run_full(initial_xyz)
        if is_lite_protocol(self.protocol_spec):
            return self.run_lite(initial_xyz)
        if is_zero_protocol(self.protocol_spec):
            return self.run_zero(initial_xyz)
        return self.run_ext(initial_xyz)

    def run_ext(self, initial_xyz: Path) -> CandidateSet:
        self._log("[S1]   [2/5] EXT funnel · two-stage baseline (GFN0→GFN2)")
        ensemble = self._run_ensemble(initial_xyz)
        paths = self.engine._step_process_ensemble(ensemble)
        candidate_set = self._build_candidate_set(paths, ranking_basis="xtb_energy", source_stage="ext_funnel")
        return self._finalize_candidate_set(
            candidate_set,
            funnel_mode="ext",
            stages_executed=["rdkit_embed", "crest_two_stage", "ensemble_processing"],
        )

    def run_full(self, initial_xyz: Path) -> CandidateSet:
        self._log("[S1]   [2/5] FULL funnel · prescreen fast-SP → screening fast-SP → survivor window")
        ensemble = self.engine._step_crest_search(initial_xyz)
        paths = self.engine._step_process_ensemble(ensemble)
        candidate_set = self._build_candidate_set(paths, ranking_basis="xtb_energy", source_stage="full_funnel")
        candidate_set = self._apply_fast_sp_profile(
            candidate_set,
            profile_name="prescreen_sp",
            target_field="prescreen_energy",
            ranking_label="prescreen_energy",
        )
        prescreen_window = self.protocol_spec.funnel_policy.prescreen_window_kcal or 4.0
        candidate_set = select_candidates_by_window(candidate_set, "prescreen_energy", prescreen_window)
        candidate_set = self._apply_fast_sp_profile(
            candidate_set,
            profile_name="screening_sp",
            target_field="screening_energy",
            ranking_label="screening_energy",
        )
        approx_thermo_applied = False
        if self.protocol_spec.funnel_policy.use_mrrho_like_correction:
            candidate_set = self._apply_mrrho_like_correction(candidate_set)
            approx_thermo_applied = bool(candidate_set.metadata.get("approx_thermo_applied", False))
        screening_window = self.protocol_spec.funnel_policy.screening_window_kcal or 3.5
        ranking_basis = "xtb_free_energy" if approx_thermo_applied else "screening_energy"
        candidate_set = select_candidates_by_window(candidate_set, ranking_basis, screening_window)
        return self._finalize_candidate_set(
            candidate_set,
            funnel_mode="full",
            stages_executed=["rdkit_embed", "crest_gfn2", "ensemble_processing", "prescreen_sp", "screening_sp"],
            energy_window_kcal=screening_window,
            approx_thermo_applied=approx_thermo_applied,
        )

    def run_lite(self, initial_xyz: Path) -> CandidateSet:
        self._log("[S1]   [2/5] LITE funnel · screening-level fast-SP → cutoff")
        ensemble = self.engine._step_crest_search(initial_xyz)
        paths = self.engine._step_process_ensemble(ensemble)
        candidate_set = self._build_candidate_set(paths, ranking_basis="screening_energy", source_stage="lite_funnel")
        ranked = self._apply_fast_sp_profile(
            candidate_set,
            profile_name="screening_sp",
            target_field="screening_energy",
            ranking_label="screening_energy",
        )
        approx_thermo_applied = False
        if self.protocol_spec.funnel_policy.use_mrrho_like_correction:
            ranked = self._apply_mrrho_like_correction(ranked)
            approx_thermo_applied = bool(ranked.metadata.get("approx_thermo_applied", False))
        ranking_basis = "xtb_free_energy" if approx_thermo_applied else "screening_energy"
        ranked = sort_candidates(ranked, ranking_basis)
        boltzmann_cutoff = self.protocol_spec.funnel_policy.boltzmann_cutoff
        if boltzmann_cutoff is not None and ranked.candidates:
            energies: List[float] = []
            for candidate in ranked.candidates:
                value = candidate.xtb_free_energy if ranking_basis == "xtb_free_energy" else candidate.screening_energy
                if value is None:
                    value = candidate.rerank_energy if candidate.rerank_energy is not None else candidate.xtb_energy
                energies.append(float(value if value is not None else candidate.rank_initial))
            min_energy = min(energies)
            thermo_config = getattr(self.engine, "thermo_config", {}) or {}
            rt = 0.0019872041 * float(thermo_config.get("temperature_k", 298.15))
            weights = [math.exp(-((energy - min_energy) * HARTREE_TO_KCAL) / rt) for energy in energies]
            total = sum(weights) or 1.0
            cumulative = 0.0
            kept: List[ConformerCandidate] = []
            for candidate, weight in zip(ranked.candidates, weights):
                cumulative += weight / total
                kept.append(candidate)
                if cumulative >= float(boltzmann_cutoff):
                    break
            ranked = clone_candidate_set(ranked, kept)
            ranked.metadata["boltzmann_cutoff"] = float(boltzmann_cutoff)
        return self._finalize_candidate_set(
            ranked,
            funnel_mode="lite",
            stages_executed=["rdkit_embed", "crest_gfn2", "ensemble_processing", "screening_sp"],
            boltzmann_cutoff=boltzmann_cutoff,
            approx_thermo_applied=approx_thermo_applied,
        )

    def run_zero(self, initial_xyz: Path) -> CandidateSet:
        self._log("[S1]   [2/5] ZERO funnel · GFN2 ranking → narrow energy window")
        ensemble = self.engine._step_crest_search(initial_xyz)
        paths = self.engine._step_process_ensemble(ensemble)
        candidate_set = self._build_candidate_set(paths, ranking_basis="xtb_energy", source_stage="zero_funnel")
        candidate_set = sort_candidates(candidate_set, "xtb_energy")
        window = self.protocol_spec.funnel_policy.narrow_window_kcal or 0.5
        marked = select_candidates_by_window(candidate_set, "xtb_energy", window)
        return self._finalize_candidate_set(
            marked,
            funnel_mode="zero",
            stages_executed=["rdkit_embed", "crest_gfn2", "ensemble_processing", "narrow_window"],
            energy_window_kcal=window,
        )

    def _run_ensemble(self, initial_xyz: Path) -> Path:
        if self.engine.two_stage_enabled and self.engine.stage1_enabled and self.engine.stage2_enabled:
            return self.engine._step_two_stage_crest(initial_xyz)
        return self.engine._step_crest_search(initial_xyz)

    def _build_candidate_set(self, paths: List[Path], ranking_basis: str, source_stage: str) -> CandidateSet:
        candidate_set = candidate_set_from_paths(
            paths,
            source_protocol=self.protocol_spec.name,
            ranking_basis=ranking_basis,
            ranking_unit="hartree",
            source_stage=source_stage,
        )
        enriched: List[ConformerCandidate] = []
        for candidate in candidate_set.candidates:
            energy = self._safe_extract_energy(candidate.source_xyz)
            ranking_history = dict(candidate.ranking_history)
            if energy is not None:
                ranking_history["xtb_energy"] = energy
            enriched.append(
                ConformerCandidate(
                    candidate_id=candidate.candidate_id,
                    source_xyz=candidate.source_xyz,
                    source_stage=candidate.source_stage,
                    rank_initial=candidate.rank_initial,
                    xtb_energy=energy,
                    xtb_free_energy=None,
                    prescreen_energy=None,
                    screening_energy=None,
                    rerank_energy=None,
                    rerank_method="xtb_energy",
                    ranking_history=ranking_history,
                    included_by=list(candidate.included_by),
                    within_window=candidate.within_window,
                    cluster_id=candidate.cluster_id,
                    metadata=dict(candidate.metadata),
                )
            )
        candidate_set.candidates = enriched
        candidate_set.selected_count_before_handoff = len(enriched)
        return candidate_set

    def _apply_fast_sp_profile(
        self,
        candidate_set: CandidateSet,
        *,
        profile_name: str,
        target_field: str,
        ranking_label: str,
    ) -> CandidateSet:
        reranked: List[ConformerCandidate] = []
        for candidate in candidate_set.candidates:
            base = candidate.xtb_energy if candidate.xtb_energy is not None else float(candidate.rank_initial)
            fast_energy = None
            if hasattr(self.engine, "_run_fast_sp_profile"):
                fast_energy = self.engine._run_fast_sp_profile(
                    candidate.source_xyz,
                    profile_name=profile_name,
                    output_subdir=profile_name,
                )
            rerank_energy = float(fast_energy if fast_energy is not None else base)
            ranking_history = dict(candidate.ranking_history)
            ranking_history[ranking_label] = rerank_energy
            prescreen_energy = candidate.prescreen_energy
            screening_energy = candidate.screening_energy
            rerank_method = candidate.rerank_method
            if target_field == "prescreen_energy":
                prescreen_energy = rerank_energy
            elif target_field == "screening_energy":
                screening_energy = rerank_energy
                rerank_method = ranking_label
            reranked.append(
                ConformerCandidate(
                    candidate_id=candidate.candidate_id,
                    source_xyz=candidate.source_xyz,
                    source_stage=candidate.source_stage,
                    rank_initial=candidate.rank_initial,
                    xtb_energy=candidate.xtb_energy,
                    xtb_free_energy=candidate.xtb_free_energy,
                    prescreen_energy=prescreen_energy,
                    screening_energy=screening_energy,
                    rerank_energy=rerank_energy,
                    rerank_method=rerank_method,
                    ranking_history=ranking_history,
                    included_by=list(candidate.included_by),
                    within_window=candidate.within_window,
                    cluster_id=candidate.cluster_id,
                    metadata=dict(candidate.metadata),
                )
            )
        candidate_set.candidates = reranked
        candidate_set.ranking_basis = target_field
        return sort_candidates(candidate_set, target_field)

    def _apply_mrrho_like_correction(self, candidate_set: CandidateSet) -> CandidateSet:
        if not any(candidate.xtb_free_energy is not None for candidate in candidate_set.candidates):
            candidate_set.metadata["approx_thermo_applied"] = False
            return candidate_set
        corrected: List[ConformerCandidate] = []
        for candidate in candidate_set.candidates:
            base = candidate.rerank_energy if candidate.rerank_energy is not None else (
                candidate.xtb_energy if candidate.xtb_energy is not None else float(candidate.rank_initial)
            )
            corrected_energy = float(base)
            ranking_history = dict(candidate.ranking_history)
            ranking_history["mrrho_like"] = corrected_energy
            corrected.append(
                ConformerCandidate(
                    candidate_id=candidate.candidate_id,
                    source_xyz=candidate.source_xyz,
                    source_stage=candidate.source_stage,
                    rank_initial=candidate.rank_initial,
                    xtb_energy=candidate.xtb_energy,
                    xtb_free_energy=corrected_energy,
                    prescreen_energy=candidate.prescreen_energy,
                    screening_energy=candidate.screening_energy,
                    rerank_energy=corrected_energy,
                    rerank_method="mrrho_like",
                    ranking_history=ranking_history,
                    included_by=list(candidate.included_by),
                    within_window=candidate.within_window,
                    cluster_id=candidate.cluster_id,
                    metadata=dict(candidate.metadata),
                )
            )
        candidate_set.candidates = corrected
        candidate_set.ranking_basis = "xtb_free_energy"
        candidate_set.metadata["approx_thermo_applied"] = True
        return sort_candidates(candidate_set, "xtb_free_energy")

    def _finalize_candidate_set(
        self,
        candidate_set: CandidateSet,
        *,
        funnel_mode: str,
        stages_executed: List[str],
        energy_window_kcal: Optional[float] = None,
        boltzmann_cutoff: Optional[float] = None,
        approx_thermo_applied: bool = False,
    ) -> CandidateSet:
        ordered = sort_candidates(candidate_set, candidate_set.ranking_basis)
        gap_rank1_rank2 = None
        if len(ordered.candidates) >= 2:
            first = self._candidate_energy(ordered.candidates[0])
            second = self._candidate_energy(ordered.candidates[1])
            if first is not None and second is not None:
                gap_rank1_rank2 = (second - first) * HARTREE_TO_KCAL
        ordered.metadata["funnel_mode"] = funnel_mode
        ordered.metadata["stages_executed"] = list(stages_executed)
        ordered.metadata["energy_window_kcal"] = energy_window_kcal
        ordered.metadata["window_count"] = len(ordered.candidates)
        ordered.metadata["gap_rank1_rank2"] = gap_rank1_rank2
        ordered.metadata["approx_thermo_applied"] = approx_thermo_applied
        if boltzmann_cutoff is not None:
            ordered.metadata["boltzmann_cutoff"] = boltzmann_cutoff
        return ordered

    def _candidate_energy(self, candidate: ConformerCandidate) -> Optional[float]:
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

    def _log(self, message: str) -> None:
        logger = getattr(self.engine, "logger", None)
        if logger is not None:
            logger.info(message)

    def _safe_extract_energy(self, xyz_path: Path) -> Optional[float]:
        try:
            return float(self.engine._extract_energy_from_xyz(xyz_path))
        except Exception:
            return None
