import json
from pathlib import Path

import pytest

from rph_core.steps.anchor.engine_factory import create_s1_engine
from rph_core.steps.anchor.handler import AnchorPhase
from rph_core.steps.conformer_search.protocols import resolve_protocol_spec


def _minimal_config(protocol: str = "ext"):
    return {
        "step1": {
            "protocol": protocol,
            "shared_handoff": {
                "enabled": True,
                "dft_optfreq": True,
                "final_sp": True,
                "ranking_after_handoff": "final_sp_minimum",
            },
            "protocol_stack": {
                "ext": {
                    "handoff": {"mode": "optimize_all_candidates"},
                },
                "full": {
                    "handoff": {"mode": "optimize_all_survivors_within_window"},
                },
                "lite": {
                    "handoff": {
                        "mode": "optimize_rank1",
                        "fallback_mode": "optimize_top2_if_gap_small",
                        "small_gap_kcal": 1.0,
                    },
                },
                "zero": {
                    "handoff": {
                        "mode": "optimize_rank1",
                        "fallback_mode": "optimize_all_within_0p5_kcal",
                    },
                },
            },
            "conformer_search": {
                "two_stage_enabled": True,
            },
        },
        "theory": {
            "optimization": {"method": "B3LYP", "basis": "def2-SVP", "engine": "gaussian"},
            "single_point": {"method": "wB97X-D4", "basis": "def2-TZVPP", "engine": "orca"},
        },
        "solvent": {"name": "acetone"},
    }


def test_create_s1_engine_accepts_lite_phase2(tmp_path: Path):
    engine = create_s1_engine(
        protocol="lite",
        config=_minimal_config("lite"),
        work_dir=tmp_path,
        molecule_name="product",
    )
    assert engine is not None
    assert engine.two_stage_enabled is False
    assert engine.ngeom_default == 1
    assert engine.max_conformers == 1
    assert engine.protocol_spec.name == "lite"
    assert engine.protocol_spec.freq_enabled is True
    assert engine.protocol_spec.final_opt_sp_enabled is True


def test_create_s1_engine_accepts_zero_phase3a(tmp_path: Path):
    engine = create_s1_engine(
        protocol="zero",
        config=_minimal_config("zero"),
        work_dir=tmp_path,
        molecule_name="product",
    )
    assert engine is not None
    assert engine.two_stage_enabled is False
    assert engine.ngeom_default == 1
    assert engine.max_conformers == 1
    assert engine.protocol_spec.name == "zero"
    assert engine.protocol_spec.freq_enabled is True
    assert engine.protocol_spec.final_opt_sp_enabled is True


def test_create_s1_engine_accepts_full_protocol(tmp_path: Path):
    engine = create_s1_engine(
        protocol="full",
        config=_minimal_config("full"),
        work_dir=tmp_path,
        molecule_name="product",
    )
    assert engine is not None
    assert engine.protocol_spec.name == "full"
    assert engine.protocol_spec.handoff_policy.mode == "optimize_all_survivors_within_window"


def test_anchor_phase_writes_phase1_provenance(monkeypatch, tmp_path: Path):
    class FakeEngine:
        def __init__(self, work_dir: Path, molecule_name: str):
            self.work_dir = work_dir
            self.molecule_name = molecule_name
            self.last_final_opt_sp_meta = {
                "stage": "disabled",
                "stage_status": "completed",
                "enabled": False,
                "protocol": "ext",
                "selection_mode": "default",
                "selected_candidate_count": 0,
                "freq_requested": False,
                "final_sp_requested": False,
            }

        def run(self, smiles: str):
            mol_dir = self.work_dir / self.molecule_name
            dft_dir = mol_dir / "finalDFT"
            dft_dir.mkdir(parents=True, exist_ok=True)
            xyz = mol_dir / f"{self.molecule_name}_global_min.xyz"
            xyz.write_text("3\nE=-1.234567\nC 0 0 0\nH 0 0 1\nH 0 1 0\n", encoding="utf-8")
            return xyz, -1.234567

        def run_optimization_only(self, smiles: str):
            return self.run(smiles)

    def fake_create_s1_engine(*, protocol, config, work_dir, molecule_name):
        assert protocol == "ext"
        return FakeEngine(work_dir=work_dir, molecule_name=molecule_name)

    monkeypatch.setattr("rph_core.steps.anchor.handler.is_small_molecule", lambda smiles, threshold: False)
    monkeypatch.setattr("rph_core.steps.anchor.handler.create_s1_engine", fake_create_s1_engine)

    anchor = AnchorPhase(config=_minimal_config("ext"), base_work_dir=tmp_path)
    result = anchor.run({"product": "C=C"})

    assert result.success is True
    assert "product" in result.anchored_molecules
    assert result.anchored_molecules["product"]["protocol"] == "ext"

    provenance_file = tmp_path / "provenance.json"
    assert provenance_file.exists()

    data = json.loads(provenance_file.read_text(encoding="utf-8"))
    assert data["schema_version"] == "s1_provenance_v1"
    assert data["protocol_spec_version"] == "protocol_spec_v1"
    assert data["protocol"] == "ext"
    assert data["has_geometry_optimization"] is True
    assert data["final_opt_sp"]["enabled"] is False
    assert data["final_opt_sp"]["stage"] == "disabled"
    assert data["final_opt_sp"]["stage_status"] == "completed"
    assert data["final_opt_sp"]["freq_requested"] is False
    assert data["final_opt_sp"]["final_sp_requested"] is False
    assert data["final_opt_sp"]["selection_mode"] == "default"
    assert data["final_opt_sp"]["stage_meta"]["stage_status"] == "completed"
    assert data["artifact_contract"]["product_geometry_level"] == "dft_optimized"
    assert data["protocol_spec"]["name"] == "ext"
    assert data["protocol_spec"]["funnel_policy"]["search_mode"] == "crest_two_stage_gfn0_to_gfn2"
    assert data["protocol_spec"]["handoff_policy"]["mode"] == "optimize_all_candidates"
    assert data["funnel"]["search_mode"] == "crest_two_stage_gfn0_to_gfn2"
    assert data["handoff"]["mode"] == "optimize_all_candidates"
    assert data["handoff"]["mode_requested"] == "optimize_all_candidates"
    assert data["handoff"]["mode_effective"] == "optimize_all_candidates"
    assert "product" in data["molecules"]
    assert data["molecules"]["product"]["status"] == "completed"


def test_anchor_phase_protocol_lite_succeeds_in_phase2(monkeypatch, tmp_path: Path):
    class FakeEngine:
        def __init__(self, work_dir: Path, molecule_name: str):
            self.work_dir = work_dir
            self.molecule_name = molecule_name
            self.last_final_opt_sp_meta = {
                "stage": "final_opt_sp",
                "stage_status": "completed",
                "enabled": True,
                "protocol": "lite",
                "selection_mode": "single_lowest_energy",
                "selected_candidate_count": 1,
                "freq_requested": True,
                "final_sp_requested": True,
            }

        def run(self, smiles: str):
            mol_dir = self.work_dir / self.molecule_name
            dft_dir = mol_dir / "finalDFT"
            dft_dir.mkdir(parents=True, exist_ok=True)
            xyz = mol_dir / f"{self.molecule_name}_global_min.xyz"
            xyz.write_text("3\nE=-2.000001\nC 0 0 0\nH 0 0 1\nH 0 1 0\n", encoding="utf-8")
            return xyz, -2.000001

        def run_optimization_only(self, smiles: str):
            return self.run(smiles)

    def fake_create_s1_engine(*, protocol, config, work_dir, molecule_name):
        assert protocol == "lite"
        return FakeEngine(work_dir=work_dir, molecule_name=molecule_name)

    monkeypatch.setattr("rph_core.steps.anchor.handler.is_small_molecule", lambda smiles, threshold: False)
    monkeypatch.setattr("rph_core.steps.anchor.handler.create_s1_engine", fake_create_s1_engine)

    anchor = AnchorPhase(config=_minimal_config("lite"), base_work_dir=tmp_path)
    result = anchor.run({"product": "C=C"})

    assert result.success is True
    assert result.anchored_molecules["product"]["protocol"] == "lite"

    data = json.loads((tmp_path / "provenance.json").read_text(encoding="utf-8"))
    assert data["protocol"] == "lite"
    assert data["final_opt_sp"]["enabled"] is True
    assert data["final_opt_sp"]["stage"] == "final_opt_sp"
    assert data["final_opt_sp"]["stage_status"] == "completed"
    assert data["final_opt_sp"]["freq_requested"] is True
    assert data["final_opt_sp"]["final_sp_requested"] is True
    assert data["final_opt_sp"]["selection_mode"] == "single_lowest_energy"
    assert data["final_opt_sp"]["selected_candidate_count"] == 1
    assert data["final_opt_sp"]["stage_meta"]["protocol"] == "lite"
    assert data["artifact_contract"]["product_geometry_level"] == "dft_optimized"
    assert data["protocol_spec"]["name"] == "lite"
    assert data["protocol_spec"]["handoff_policy"]["mode"] == "optimize_rank1"
    assert data["handoff"]["mode"] == "optimize_rank1"
    assert data["handoff"]["mode_requested"] == "optimize_rank1"
    assert data["handoff"]["mode_effective"] == "optimize_rank1"
    assert data["conformer_search"]["ngeom_max"] == 1


def test_anchor_phase_protocol_zero_succeeds_in_phase3a(monkeypatch, tmp_path: Path):
    class FakeEngine:
        def __init__(self, work_dir: Path, molecule_name: str):
            self.work_dir = work_dir
            self.molecule_name = molecule_name
            self.last_final_opt_sp_meta = {
                "stage": "final_opt_sp",
                "stage_status": "completed",
                "enabled": True,
                "protocol": "zero",
                "selection_mode": "single_lowest_energy",
                "selected_candidate_count": 1,
                "freq_requested": False,
                "final_sp_requested": True,
            }

        def run(self, smiles: str):
            mol_dir = self.work_dir / self.molecule_name
            dft_dir = mol_dir / "finalDFT"
            dft_dir.mkdir(parents=True, exist_ok=True)
            xyz = mol_dir / f"{self.molecule_name}_global_min.xyz"
            xyz.write_text("3\nE=-3.000001\nC 0 0 0\nH 0 0 1\nH 0 1 0\n", encoding="utf-8")
            return xyz, -3.000001

        def run_optimization_only(self, smiles: str):
            return self.run(smiles)

    def fake_create_s1_engine(*, protocol, config, work_dir, molecule_name):
        assert protocol == "zero"
        return FakeEngine(work_dir=work_dir, molecule_name=molecule_name)

    monkeypatch.setattr("rph_core.steps.anchor.handler.is_small_molecule", lambda smiles, threshold: False)
    monkeypatch.setattr("rph_core.steps.anchor.handler.create_s1_engine", fake_create_s1_engine)

    anchor = AnchorPhase(config=_minimal_config("zero"), base_work_dir=tmp_path)
    result = anchor.run({"product": "C=C"})

    assert result.success is True
    assert result.anchored_molecules["product"]["protocol"] == "zero"

    data = json.loads((tmp_path / "provenance.json").read_text(encoding="utf-8"))
    assert data["protocol"] == "zero"
    assert data["final_opt_sp"]["enabled"] is True
    assert data["final_opt_sp"]["stage"] == "final_opt_sp"
    assert data["final_opt_sp"]["stage_status"] == "completed"
    assert data["final_opt_sp"]["freq_requested"] is False
    assert data["final_opt_sp"]["final_sp_requested"] is True
    assert data["final_opt_sp"]["selection_mode"] == "single_lowest_energy"
    assert data["final_opt_sp"]["selected_candidate_count"] == 1
    assert data["final_opt_sp"]["stage_meta"]["protocol"] == "zero"
    assert data["artifact_contract"]["product_geometry_level"] == "dft_optimized"
    assert data["protocol_spec"]["name"] == "zero"
    assert data["protocol_spec"]["handoff_policy"]["fallback_mode"] == "optimize_all_within_0p5_kcal"
    assert data["handoff"]["fallback_mode"] == "optimize_all_within_0p5_kcal"
    assert data["handoff"]["mode_requested"] == "optimize_rank1"
    assert data["handoff"]["mode_effective"] == "optimize_rank1"


def test_resolve_protocol_spec_ext_and_default_use_compatible_semantics() -> None:
    config = _minimal_config("ext")

    ext_spec = resolve_protocol_spec(config, "ext")
    default_spec = resolve_protocol_spec(config, "default")

    assert ext_spec.name == "ext"
    assert default_spec.name == "default"
    assert ext_spec.final_opt_sp_enabled is True
    assert default_spec.final_opt_sp_enabled is True
    assert ext_spec.two_stage_enabled is True
    assert default_spec.two_stage_enabled is True
