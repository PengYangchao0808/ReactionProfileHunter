from __future__ import annotations

import inspect
import json
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
from dataclasses import dataclass

from rph_core.scheduling.models import ReactionJob, BranchJob, ConditionJob
from rph_core.scheduling.v3_scheduler import V3Scheduler, _load_dr_branch_plan, _to_float
from rph_core.utils.task_builder import TaskSpec
from rph_core.orchestrator import ReactionProfileHunter


@pytest.fixture
def tmp_output(tmp_path: Path) -> Path:
    return tmp_path / "rph_output"


@pytest.fixture
def mock_hunter(tmp_output: Path):
    hunter = MagicMock()
    hunter.config = {
        "scheduler": {"mode": "v3", "v3": {
            "precompute_small_molecules": False,
            "precursor_s1_once_per_reaction": True,
            "branch_product_only_s1": True,
            "condition_branch_thermo": True,
        }},
        "theory": {
            "optimization": {"method": "B3LYP", "basis": "def2-SVP", "solvent": "acetone"},
            "single_point": {"method": "wB97M-V", "basis": "def2-TZVPP"},
        },
        "global": {},
        "s0": {"dr_completion": {"enabled": True}},
        "thermo": {"temperature_k": 298.15},
        "reaction_profiles": {},
    }
    hunter.s4_engine = MagicMock()
    hunter.s4_engine.run = MagicMock()
    hunter.logger = MagicMock()
    return hunter


@pytest.fixture
def mock_run_cfg(tmp_output: Path) -> dict[str, object]:
    return {
        "source": "dataset",
        "output_root": str(tmp_output),
        "skip_steps": ["s4"],
        "dry_run": False,
        "resume": False,
    }


def _make_task(rx_id="rx_001", reaction_id="rx_001", product_smiles="C=C(C)C(=O)O",
               condition_id="COND_001", **meta_kw):
    meta = {
        "precursor_smiles": "CC(=O)OC",
        "leaving_small_molecule_key": "AcOH",
        "small_molecular_keys": ["AcOH"],
        "reaction_profile": "[4+3]_default",
        "cleaner_data": {"temperature_c": 25},
        **meta_kw,
    }
    return TaskSpec(
        rx_id=rx_id,
        row_id=rx_id,
        reaction_id=reaction_id,
        product_smiles=product_smiles,
        condition_id=condition_id,
        meta=meta,
    )


class TestV3SchedulerPlan:
    @patch("rph_core.scheduling.v3_scheduler.build_tasks_from_run_config")
    def test_plan_creates_reaction_jobs(self, mock_build, mock_hunter, mock_run_cfg, tmp_output):
        mock_build.return_value = [
            _make_task("rx_001", "RXN_A", "SMILES_A", "COND_001"),
            _make_task("rx_002", "RXN_A", "SMILES_A", "COND_002"),
            _make_task("rx_003", "RXN_B", "SMILES_B", "COND_003"),
        ]
        scheduler = V3Scheduler(hunter=mock_hunter, run_cfg=mock_run_cfg)
        jobs = scheduler.plan()
        assert len(jobs) == 2
        rxn_a = [j for j in jobs if j.reaction_id == "RXN_A"][0]
        rxn_b = [j for j in jobs if j.reaction_id == "RXN_B"][0]
        assert len(rxn_a.conditions) == 2
        assert len(rxn_b.conditions) == 1
        assert (tmp_output / "RXN_A").exists()
        assert (tmp_output / "RXN_B").exists()

    @patch("rph_core.scheduling.v3_scheduler.build_tasks_from_run_config")
    def test_plan_reads_dr_branch_plan(self, mock_build, mock_hunter, mock_run_cfg, tmp_output):
        mock_build.return_value = [_make_task("rx_001", "RXN_A")]
        scheduler = V3Scheduler(hunter=mock_hunter, run_cfg=mock_run_cfg)
        jobs_init = scheduler.plan()
        rxn_root = jobs_init[0].reaction_root
        s0_dir = rxn_root / "S0_Mechanism"
        s0_dir.mkdir(parents=True, exist_ok=True)
        dr_plan = {
            "branches": [
                {"branch_id": "BR_MAJOR"},
                {"branch_id": "BR_DR_001", "pathway_id": "dr1", "product_smiles": "C=C(C)C(=O)O",
                 "generation_policy": "flip", "flipped_map_numbers": [3, 4], "notes": ["test"]},
            ]
        }
        with open(s0_dir / "dr_branch_plan.json", "w") as f:
            json.dump(dr_plan, f)

        jobs = scheduler.plan()
        assert len(jobs) == 1
        assert [branch.branch_id for branch in jobs[0].branches] == ["BR_MAJOR", "BR_DR_001"]


class TestV3SchedulerHelpers:
    def test_to_float(self):
        assert _to_float("25.0") == 25.0
        assert _to_float(None) is None
        assert _to_float("abc") is None
        assert _to_float("298.15") == 298.15

    def test_load_dr_branch_plan_missing(self, tmp_path):
        assert _load_dr_branch_plan(tmp_path) is None

    def test_load_dr_branch_plan_valid(self, tmp_path):
        s0 = tmp_path / "S0_Mechanism"
        s0.mkdir()
        plan = {"branches": [{"branch_id": "BR_MAJOR"}]}
        with open(s0 / "dr_branch_plan.json", "w") as f:
            json.dump(plan, f)
        result = _load_dr_branch_plan(tmp_path)
        assert result is not None
        assert len(result["branches"]) == 1


class TestV3SchedulerDryRun:
    @patch("rph_core.scheduling.v3_scheduler.build_tasks_from_run_config")
    def test_dry_run_returns_empty(self, mock_build, mock_hunter, tmp_output):
        mock_build.return_value = [_make_task()]
        mock_run_cfg = {
            "source": "dataset",
            "output_root": str(tmp_output),
            "skip_steps": [],
            "dry_run": True,
        }
        scheduler = V3Scheduler(hunter=mock_hunter, run_cfg=mock_run_cfg)
        results = scheduler.run()
        assert results == []


class TestV3BranchExecutionContract:
    def test_orchestrator_accepts_reference_small_molecule_flag(self):
        signature = inspect.signature(ReactionProfileHunter.run_pipeline)
        parameter = signature.parameters["include_reference_small_molecules"]

        assert parameter.default is True

    def test_shared_precursor_reuses_existing_global_min(self, mock_hunter, mock_run_cfg, tmp_output):
        from rph_core.scheduling.models import ReactionJob

        scheduler = V3Scheduler(hunter=mock_hunter, run_cfg=mock_run_cfg)
        task = _make_task("rx_001", "RXN_A")
        job = ReactionJob(
            reaction_id="RXN_A",
            reaction_root=tmp_output / "RXN_A",
            representative=task,
        )
        precursor_dir = job.reaction_root / "precursor" / "S1_ConfGeneration" / "precursor"
        precursor_dir.mkdir(parents=True)
        (precursor_dir / "precursor_global_min.xyz").write_text("1\nmin\nH 0 0 0\n")

        scheduler._run_precursor_s1(job)

        assert (precursor_dir / "precursor_min.xyz").exists()
        mock_hunter.s1_engine.run.assert_not_called()

    def test_branch_run_does_not_skip_s4_or_reference_small_molecules(self, mock_hunter, mock_run_cfg, tmp_output):
        from rph_core.scheduling.models import BranchJob, ReactionJob

        mock_run_cfg.pop("skip_steps", None)
        mock_result = MagicMock(success=True, features_csv=tmp_output / "features.csv")
        mock_hunter.run_pipeline.return_value = mock_result
        scheduler = V3Scheduler(hunter=mock_hunter, run_cfg=mock_run_cfg)
        task = _make_task("rx_001", "RXN_A")
        job = ReactionJob(
            reaction_id="RXN_A",
            reaction_root=tmp_output / "RXN_A",
            representative=task,
        )
        branch = BranchJob(
            parent_reaction_id="RXN_A",
            branch_id="BR_DR_001",
            pathway_id="dr1",
            product_smiles="C=C",
            branch_root=tmp_output / "RXN_A" / "branches" / "BR_DR_001",
        )

        scheduler._run_branch(job, branch)

        kwargs = mock_hunter.run_pipeline.call_args.kwargs
        assert "s4" not in [str(s).lower() for s in kwargs["skip_steps"]]
        assert kwargs["precursor_smiles"] is None
        assert kwargs["small_molecular_keys"] == []
        assert kwargs["include_reference_small_molecules"] is False

    @patch("rph_core.scheduling.v3_scheduler.build_tasks_from_run_config")
    def test_plan_rejects_unsafe_branch_id(self, mock_build, mock_hunter, mock_run_cfg, tmp_output):
        mock_build.return_value = [_make_task("rx_001", "RXN_A")]
        scheduler = V3Scheduler(hunter=mock_hunter, run_cfg=mock_run_cfg)
        jobs_init = scheduler.plan()
        rxn_root = jobs_init[0].reaction_root
        s0_dir = rxn_root / "S0_Mechanism"
        s0_dir.mkdir(parents=True, exist_ok=True)
        with open(s0_dir / "dr_branch_plan.json", "w") as f:
            json.dump({"branches": [{"branch_id": "../escape"}]}, f)

        with pytest.raises(ValueError, match="Unsafe branch_id"):
            scheduler.plan()
