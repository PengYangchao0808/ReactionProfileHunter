import json
from pathlib import Path

from rph_core.utils.checkpoint_manager import CheckpointManager
from rph_core.utils.layout_contract import resolve_required_files


def test_mark_step_in_progress_persists_phase(tmp_path: Path) -> None:
    mgr = CheckpointManager(tmp_path)
    mgr.initialize_state(product_smiles="C", config={})

    mgr.mark_step_in_progress("s2", phase="scan", output_files={"ts_guess_xyz": "foo.xyz"})

    state = mgr.load_state()
    assert state is not None
    step = state.steps["step_s2"]
    assert step.completed is False
    assert step.metadata is not None
    assert step.metadata.get("phase") == "scan"
    assert step.output_files.get("ts_guess_xyz") == "foo.xyz"


def test_mark_step_failed_partial_records_error(tmp_path: Path) -> None:
    mgr = CheckpointManager(tmp_path)
    mgr.initialize_state(product_smiles="C", config={})

    mgr.mark_step_failed_partial("s3", phase="optimize", error_message="boom")

    state = mgr.load_state()
    assert state is not None
    step = state.steps["step_s3"]
    assert step.completed is False
    assert step.metadata is not None
    assert step.metadata.get("phase") == "optimize"
    assert step.metadata.get("error_stage") == "optimize"
    assert step.metadata.get("error_message") == "boom"


def test_s0_contract_requires_graph_and_summary(tmp_path: Path) -> None:
    required = resolve_required_files(tmp_path, "s0")
    assert set(required.keys()) == {
        "mechanism_graph_json",
        "mechanism_summary_json",
        "dr_branch_plan_json",
    }
    assert required["mechanism_graph_json"].name == "mechanism_graph.json"
    assert required["mechanism_summary_json"].name == "mechanism_summary.json"
    assert required["dr_branch_plan_json"].name == "dr_branch_plan.json"


def test_rehydrate_step2_uses_executable_forming_bonds_in_metadata_and_signature(tmp_path: Path) -> None:
    mgr = CheckpointManager(tmp_path)

    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir(parents=True, exist_ok=True)
    (s1_dir / "product_min.xyz").write_text(
        "4\nproduct\nC 0 0 0\nC 0 0 3\nC 2 0 0\nC 2 0 3\n",
        encoding="utf-8",
    )

    s2_dir = tmp_path / "S2_Retro"
    s2_dir.mkdir(parents=True, exist_ok=True)
    (s2_dir / "ts_guess.xyz").write_text("2\nts\nC 0 0 0\nC 0 0 1\n", encoding="utf-8")
    (s2_dir / "intermediate.xyz").write_text("2\nint\nC 0 0 0\nC 0 0 1\n", encoding="utf-8")
    (s2_dir / "scan_profile.json").write_text(
        json.dumps({"forming_bonds": [[0, 2], [1, 3]]}, indent=2),
        encoding="utf-8",
    )

    state = mgr.rehydrate_state_from_artifacts(product_smiles="C", config={})

    assert state is not None
    step_s2 = state.steps["step_s2"]
    assert step_s2.metadata is not None
    assert step_s2.metadata["forming_bonds"] == [[0, 2], [1, 3]]
    step2_signature = step_s2.metadata.get("step2_signature")
    assert isinstance(step2_signature, dict)
    assert step2_signature["forming_bonds"] == [[0, 2], [1, 3]]
