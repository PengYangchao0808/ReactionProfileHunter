import pytest
import json
from rph_core.utils.checkpoint_manager import CheckpointManager
from rph_core.steps.step4_features.feature_miner import FeatureMiner
from pathlib import Path

@pytest.fixture
def pipeline_dirs(tmp_path: Path) -> dict[str, Path]:
    s1 = tmp_path / "S1_ConfGeneration"
    s2 = tmp_path / "S2_Retro"
    s3 = tmp_path / "S3_TS"
    s4 = tmp_path / "S4_Data"
    for d in (s1, s2, s3, s4):
        d.mkdir(parents=True, exist_ok=True)
    return {"S1": s1, "S2": s2, "S3": s3, "S4": s4}

def _write_minimal_xyz(path: Path) -> None:
    path.write_text(
        """3
comment
C 0.0 0.0 0.0
H 0.0 1.0 0.0
H 0.0 0.0 1.0
"""
    )

def test_is_step4_complete_requires_mech_index_when_enabled(tmp_path: Path, pipeline_dirs):
    mgr = CheckpointManager(tmp_path)
    mgr.initialize_state(product_smiles="C", config={})
    mgr.mark_step_completed("s4", output_files={})

    s4_dir = pipeline_dirs["S4"]
    (s4_dir / "features_raw.csv").write_text("reaction_id,dG_activation\n")

    config = {"step4": {"mechanism_packaging": {"enabled": False}}}
    assert mgr.is_step4_complete(s4_dir, config) is True
    assert not (s4_dir / "mech_index.json").exists()
