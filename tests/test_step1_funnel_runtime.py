from pathlib import Path

from rph_core.steps.conformer_search.funnel import FunnelRunner
from rph_core.steps.conformer_search.protocols import resolve_protocol_spec


def _base_config(protocol: str):
    return {
        "step1": {
            "protocol": protocol,
            "shared_handoff": {"enabled": True, "dft_optfreq": True, "final_sp": True},
            "protocol_stack": {
                "ext": {"funnel": {"survivor_window_kcal": 3.0}, "handoff": {"mode": "optimize_all_candidates"}},
                "full": {
                    "funnel": {
                        "search_mode": "crest_gfn2",
                        "survivor_window_kcal": 3.0,
                        "prescreen_mode": "low_cost_dft_sp",
                    },
                    "handoff": {"mode": "optimize_all_survivors_within_window"},
                },
                "lite": {
                    "funnel": {"search_mode": "crest_gfn2", "use_mrrho_like_correction": True, "boltzmann_cutoff": 0.90},
                    "handoff": {"mode": "optimize_rank1", "fallback_mode": "optimize_top2_if_gap_small", "small_gap_kcal": 1.0},
                },
                "zero": {
                    "funnel": {"search_mode": "crest_gfn2_or_skip", "narrow_window_kcal": 0.5},
                    "handoff": {"mode": "optimize_rank1", "fallback_mode": "optimize_all_within_0p5_kcal"},
                },
            },
        }
    }


class _FakeEngine:
    def __init__(self, tmp_path: Path):
        self.tmp_path = tmp_path
        self.two_stage_enabled = True
        self.stage1_enabled = True
        self.stage2_enabled = True

    def _step_two_stage_crest(self, input_xyz: Path) -> Path:
        return input_xyz

    def _step_crest_search(self, input_xyz: Path) -> Path:
        return input_xyz

    def _step_process_ensemble(self, ensemble_file: Path):
        paths = []
        for idx in range(3):
            p = self.tmp_path / f"conf_{idx:03d}.xyz"
            p.write_text(f"1\nE={float(idx):.6f}\nH 0 0 {idx}\n", encoding="utf-8")
            paths.append(p)
        return paths

    def _extract_energy_from_xyz(self, xyz_path: Path) -> float:
        line = xyz_path.read_text(encoding="utf-8").splitlines()[1]
        return float(line.split("=")[1])


def test_ext_funnel_keeps_all_candidates(tmp_path: Path) -> None:
    engine = _FakeEngine(tmp_path)
    spec = resolve_protocol_spec(_base_config("ext"), "ext")
    runner = FunnelRunner(engine, spec)

    initial = tmp_path / "init.xyz"
    initial.write_text("1\nE=0.0\nH 0 0 0\n", encoding="utf-8")
    candidate_set = runner.run(initial)

    assert candidate_set.source_protocol == "ext"
    assert candidate_set.total_count == 3
    assert len(candidate_set.candidates) == 3


def test_full_funnel_applies_survivor_window(tmp_path: Path) -> None:
    engine = _FakeEngine(tmp_path)
    spec = resolve_protocol_spec(_base_config("full"), "full")
    runner = FunnelRunner(engine, spec)

    initial = tmp_path / "init.xyz"
    initial.write_text("1\nE=0.0\nH 0 0 0\n", encoding="utf-8")
    candidate_set = runner.run(initial)

    assert candidate_set.source_protocol == "full"
    assert len(candidate_set.candidates) >= 1
    assert candidate_set.metadata["funnel_mode"] == "full"
    assert "prescreen_sp" in candidate_set.metadata["stages_executed"]
    assert candidate_set.metadata["energy_window_kcal"] == 3.5


def test_lite_funnel_sets_rerank_basis(tmp_path: Path) -> None:
    engine = _FakeEngine(tmp_path)
    spec = resolve_protocol_spec(_base_config("lite"), "lite")
    runner = FunnelRunner(engine, spec)

    initial = tmp_path / "init.xyz"
    initial.write_text("1\nE=0.0\nH 0 0 0\n", encoding="utf-8")
    candidate_set = runner.run(initial)

    assert candidate_set.source_protocol == "lite"
    assert candidate_set.ranking_basis in {"xtb_free_energy", "screening_energy"}
    assert candidate_set.metadata["funnel_mode"] == "lite"
    assert candidate_set.metadata["boltzmann_cutoff"] == 0.9


def test_zero_funnel_marks_narrow_window(tmp_path: Path) -> None:
    engine = _FakeEngine(tmp_path)
    spec = resolve_protocol_spec(_base_config("zero"), "zero")
    runner = FunnelRunner(engine, spec)

    initial = tmp_path / "init.xyz"
    initial.write_text("1\nE=0.0\nH 0 0 0\n", encoding="utf-8")
    candidate_set = runner.run(initial)

    assert candidate_set.source_protocol == "zero"
    assert candidate_set.metadata["funnel_mode"] == "zero"
    assert candidate_set.metadata["energy_window_kcal"] == 0.5
    assert "narrow_window" in candidate_set.metadata["stages_executed"]
