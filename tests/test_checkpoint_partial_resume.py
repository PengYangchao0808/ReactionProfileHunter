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
    assert set(required.keys()) == {"mechanism_graph_json", "mechanism_summary_json"}
    assert required["mechanism_graph_json"].name == "mechanism_graph.json"
    assert required["mechanism_summary_json"].name == "mechanism_summary.json"
