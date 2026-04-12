from pathlib import Path
from typing import List, Optional

import pytest

from rph_core.steps.conformer_search.candidates import CandidateSet, ConformerCandidate
from rph_core.steps.conformer_search.pipeline.stages.final_opt_sp import execute_final_opt_sp
from rph_core.steps.conformer_search.protocols import FunnelPolicy, HandoffPolicy, ProtocolSpec


def _spec(name: str, *, freq: bool, handoff_mode: str = "optimize_rank1", fallback_mode: Optional[str] = None) -> ProtocolSpec:
    return ProtocolSpec(
        name=name,
        two_stage_enabled=False,
        ngeom_default=1,
        ngeom_max=1,
        funnel_policy=FunnelPolicy(
            search_mode="crest_gfn2",
            clustering_mode="isostat",
            prescreen_mode="none",
            rerank_mode="none",
            use_mrrho_like_correction=False,
            survivor_window_kcal=None,
            narrow_window_kcal=0.5 if name == "zero" else None,
            prescreen_window_kcal=None,
            screening_window_kcal=None,
            optimize_limit=1,
            top2_fallback_enabled=False,
            boltzmann_cutoff=None,
        ),
        handoff_policy=HandoffPolicy(
            enabled=True,
            mode=handoff_mode,
            fallback_mode=fallback_mode,
            small_gap_kcal=1.0 if fallback_mode == "optimize_top2_if_gap_small" else None,
            ranking_after_handoff="final_sp_minimum",
        ),
        final_opt_sp_enabled=True,
        freq_enabled=freq,
        final_sp_enabled=True,
        selection_mode=handoff_mode,
        provenance_flags={},
        raw_config={},
    )


def _candidate(path: Path, idx: int, energy: Optional[float] = None) -> ConformerCandidate:
    return ConformerCandidate(
        candidate_id=f"cand_{idx:03d}",
        source_xyz=path.resolve(),
        source_stage="ensemble",
        rank_initial=idx,
        xtb_energy=energy,
        xtb_free_energy=None,
        prescreen_energy=None,
        screening_energy=None,
        rerank_energy=energy,
        rerank_method="test",
        ranking_history={},
        included_by=["input"],
        within_window=True,
        cluster_id=None,
        metadata={},
    )


def _candidate_set(paths: List[Path], energies: List[float], protocol: str) -> CandidateSet:
    candidates = [_candidate(path, idx, energies[idx]) for idx, path in enumerate(paths)]
    return CandidateSet(
        candidates=candidates,
        source_protocol=protocol,
        ranking_basis="rerank_energy",
        ranking_unit="kcal/mol",
        total_count=len(candidates),
        selected_count_before_handoff=len(candidates),
        metadata={},
    )


def test_execute_final_opt_sp_rejects_empty_candidates() -> None:
    calls = {"count": 0}

    def runner(candidates: List[ConformerCandidate]):
        calls["count"] += 1
        return candidates[0].source_xyz, -1.0

    with pytest.raises(ValueError, match="at least one candidate"):
        execute_final_opt_sp(candidates=[], runner=runner, protocol_spec=_spec("lite", freq=True))

    assert calls["count"] == 0


def test_execute_final_opt_sp_single_candidate_returns_stage_meta(tmp_path: Path) -> None:
    input_candidate = tmp_path / "candidate.xyz"
    input_candidate.write_text("1\nC\nC 0 0 0\n", encoding="utf-8")
    call_count = 0
    last_candidates: List[ConformerCandidate] = []

    def runner(candidates: List[ConformerCandidate]):
        nonlocal call_count, last_candidates
        call_count += 1
        last_candidates = candidates
        return candidates[0].source_xyz.with_suffix(".log"), -123.456

    result = execute_final_opt_sp(
        candidates=[input_candidate],
        runner=runner,
        protocol_spec=_spec("lite", freq=True),
    )

    assert call_count == 1
    assert len(last_candidates) == 1
    assert last_candidates[0].source_xyz == input_candidate.resolve()
    assert result.best_log == input_candidate.with_suffix(".log").resolve()
    assert result.best_result_path == input_candidate.with_suffix(".log").resolve()
    assert result.best_energy == -123.456
    assert result.stage_meta["stage"] == "final_opt_sp"
    assert result.stage_meta["stage_status"] == "completed"
    assert result.stage_meta["protocol"] == "lite"
    assert result.stage_meta["selected_candidate_count"] == 1
    assert result.stage_meta["freq_requested"] is True
    assert result.stage_meta["final_sp_requested"] is True
    assert result.stage_meta["selection_mode"] == "optimize_rank1"
    assert result.stage_meta["mode_requested"] == "optimize_rank1"
    assert result.stage_meta["mode_effective"] == "optimize_rank1"


def test_execute_final_opt_sp_rank1_uses_first_candidate(tmp_path: Path) -> None:
    c1 = tmp_path / "c1.xyz"
    c2 = tmp_path / "c2.xyz"
    c1.write_text("1\nC\nC 0 0 0\n", encoding="utf-8")
    c2.write_text("1\nC\nC 1 0 0\n", encoding="utf-8")
    calls = {"last_candidates": []}

    def runner(candidates: List[ConformerCandidate]):
        calls["last_candidates"] = candidates
        return candidates[0].source_xyz.with_suffix(".log"), -10.0

    result = execute_final_opt_sp(
        candidate_set=_candidate_set([c1, c2], [0.0, 0.8], "zero"),
        runner=runner,
        protocol_spec=_spec("zero", freq=True, handoff_mode="optimize_rank1"),
    )

    assert len(calls["last_candidates"]) == 1
    assert calls["last_candidates"][0].source_xyz == c1.resolve()
    assert result.stage_meta["input_candidate_count"] == 2
    assert result.stage_meta["selected_candidate_count"] == 1
    assert result.stage_meta["selection_mode"] == "optimize_rank1"
    assert result.stage_meta["mode_effective"] == "optimize_rank1"


@pytest.mark.parametrize(
    ("protocol", "freq_requested"),
    [("lite", True), ("zero", True)],
)
def test_execute_final_opt_sp_freq_requested_follows_protocol(protocol: str, freq_requested: bool, tmp_path: Path) -> None:
    input_candidate = tmp_path / "candidate.xyz"
    input_candidate.write_text("1\nC\nC 0 0 0\n", encoding="utf-8")

    def runner(candidates: List[ConformerCandidate]):
        return candidates[0].source_xyz.with_suffix(".log"), -7.5

    result = execute_final_opt_sp(
        candidates=[input_candidate],
        runner=runner,
        protocol_spec=_spec(protocol, freq=freq_requested),
    )

    assert result.stage_meta["protocol"] == protocol
    assert result.stage_meta["freq_requested"] is freq_requested
    assert result.stage_meta["final_sp_requested"] is True


def test_execute_final_opt_sp_top2_fallback_when_gap_small(tmp_path: Path) -> None:
    c1 = tmp_path / "c1.xyz"
    c2 = tmp_path / "c2.xyz"
    c3 = tmp_path / "c3.xyz"
    for p in [c1, c2, c3]:
        p.write_text("1\nC\nC 0 0 0\n", encoding="utf-8")

    def runner(candidates: List[ConformerCandidate]):
        return candidates[0].source_xyz.with_suffix(".log"), -5.0

    result = execute_final_opt_sp(
        candidate_set=_candidate_set([c1, c2, c3], [0.0, 0.4, 1.8], "lite"),
        runner=runner,
        protocol_spec=_spec(
            "lite",
            freq=True,
            handoff_mode="optimize_rank1",
            fallback_mode="optimize_top2_if_gap_small",
        ),
    )

    assert result.stage_meta["selected_candidate_count"] == 2
    assert result.stage_meta["fallback_triggered"] is True
    assert result.stage_meta["fallback_trigger"] is True
    assert result.stage_meta["fallback_mode"] == "optimize_top2_if_gap_small"
    assert result.stage_meta["mode_requested"] == "optimize_rank1"
    assert result.stage_meta["mode_effective"] == "optimize_top2_if_gap_small"
