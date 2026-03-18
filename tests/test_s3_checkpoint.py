import pytest
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

from rph_core.utils.checkpoint_manager import CheckpointManager
from rph_core.steps.step3_opt.ts_optimizer import SPMatrixReport


@pytest.fixture
def mock_config():
    return {
        "theory": {
            "optimization": {
                "engine": "gaussian",
                "method": "B3LYP",
                "basis": "def2-SVP",
                "dispersion": "GD3BJ",
                "route": "",
                "rescue_route": "",
                "nproc": 16,
                "mem": "32GB",
            },
            "single_point": {
                "engine": "orca",
                "method": "WB97M-V",
                "basis": "def2-TZVPP",
                "aux_basis": "def2/J",
                "solvent": "acetone",
                "nproc": 16,
                "maxcore": 4000,
            }
        },
        "step3": {
            "reactant_opt": {
                "charge": 0,
                "multiplicity": 1,
                "enable_nbo": False,
            }
        }
    }


def _write_valid_xyz(path: Path, energy: float = -100.0) -> None:
    path.write_text(
        f"""3
E={energy:.6f}
C 0.0 0.0 0.0
H 0.0 1.0 0.0
H 0.0 0.0 1.0
"""
    )


def _write_sp_metadata(path: Path, e_ts: float = -100.0, e_reactant: float = -99.0) -> None:
    path.write_text(
        json.dumps({
            "e_ts": e_ts,
            "e_reactant": e_reactant,
            "e_product": -98.0,
            "activation_energy_kcal": (e_ts - e_reactant) * 627.509,
            "method": "B3LYP/def2-SVP",
            "solvent": "acetone",
            "timestamp": "2026-01-01T00:00:00"
        })
    )


class TestS3Checkpoint:
    def test_is_step3_complete_no_checkpoint(self, tmp_path: Path, mock_config):
        mgr = CheckpointManager(tmp_path)
        mgr.initialize_state(product_smiles="C", config=mock_config)
        
        s3_dir = tmp_path / "S3_TransitionAnalysis"
        s3_dir.mkdir(parents=True, exist_ok=True)
        
        assert mgr.is_step3_complete(s3_dir, mock_config) is False

    def test_is_step3_complete_with_valid_checkpoint(self, tmp_path: Path, mock_config):
        mgr = CheckpointManager(tmp_path)
        mgr.initialize_state(product_smiles="C", config=mock_config)
        
        s3_dir = tmp_path / "S3_TransitionAnalysis"
        s3_dir.mkdir(parents=True, exist_ok=True)
        
        ts_final = s3_dir / "ts_final.xyz"
        _write_valid_xyz(ts_final)
        
        sp_meta = s3_dir / "sp_matrix_metadata.json"
        _write_sp_metadata(sp_meta)
        
        step3_sig = mgr._compute_step3_signature(mock_config)
        
        mgr.mark_step_completed(
            "s3",
            output_files={
                "ts_final_xyz": str(ts_final),
                "sp_matrix_metadata_json": str(sp_meta),
            },
            metadata={
                "step3_signature": step3_sig,
                "input_hashes": {
                    "ts_guess": "abc123",
                    "reactant": "def456",
                    "product": "ghi789",
                }
            }
        )
        
        assert mgr.is_step3_complete(
            s3_dir, 
            mock_config, 
            check_signature=True,
            input_hashes={
                "ts_guess": "abc123",
                "reactant": "def456",
                "product": "ghi789",
            }
        ) is True

    def test_is_step3_complete_signature_mismatch(self, tmp_path: Path, mock_config):
        mgr = CheckpointManager(tmp_path)
        mgr.initialize_state(product_smiles="C", config=mock_config)
        
        s3_dir = tmp_path / "S3_TransitionAnalysis"
        s3_dir.mkdir(parents=True, exist_ok=True)
        
        ts_final = s3_dir / "ts_final.xyz"
        _write_valid_xyz(ts_final)
        
        sp_meta = s3_dir / "sp_matrix_metadata.json"
        _write_sp_metadata(sp_meta)
        
        old_sig = mgr._compute_step3_signature(mock_config)
        
        mgr.mark_step_completed(
            "s3",
            output_files={
                "ts_final_xyz": str(ts_final),
                "sp_matrix_metadata_json": str(sp_meta),
            },
            metadata={
                "step3_signature": old_sig,
            }
        )
        
        modified_config = mock_config.copy()
        modified_config["theory"]["optimization"]["method"] = "M06-2X"
        
        assert mgr.is_step3_complete(s3_dir, modified_config, check_signature=True) is False

    def test_is_step3_complete_input_hash_mismatch(self, tmp_path: Path, mock_config):
        mgr = CheckpointManager(tmp_path)
        mgr.initialize_state(product_smiles="C", config=mock_config)
        
        s3_dir = tmp_path / "S3_TransitionAnalysis"
        s3_dir.mkdir(parents=True, exist_ok=True)
        
        ts_final = s3_dir / "ts_final.xyz"
        _write_valid_xyz(ts_final)
        
        sp_meta = s3_dir / "sp_matrix_metadata.json"
        _write_sp_metadata(sp_meta)
        
        step3_sig = mgr._compute_step3_signature(mock_config)
        
        mgr.mark_step_completed(
            "s3",
            output_files={
                "ts_final_xyz": str(ts_final),
                "sp_matrix_metadata_json": str(sp_meta),
            },
            metadata={
                "step3_signature": step3_sig,
                "input_hashes": {
                    "ts_guess": "abc123",
                    "reactant": "def456",
                    "product": "ghi789",
                }
            }
        )
        
        assert mgr.is_step3_complete(
            s3_dir, 
            mock_config, 
            check_signature=True,
            input_hashes={
                "ts_guess": "different_hash",
                "reactant": "def456",
                "product": "ghi789",
            }
        ) is False

    def test_is_step3_complete_missing_ts_final(self, tmp_path: Path, mock_config):
        mgr = CheckpointManager(tmp_path)
        mgr.initialize_state(product_smiles="C", config=mock_config)
        
        s3_dir = tmp_path / "S3_TransitionAnalysis"
        s3_dir.mkdir(parents=True, exist_ok=True)
        
        sp_meta = s3_dir / "sp_matrix_metadata.json"
        _write_sp_metadata(sp_meta)
        
        mgr.mark_step_completed(
            "s3",
            output_files={
                "ts_final_xyz": str(s3_dir / "nonexistent.xyz"),
                "sp_matrix_metadata_json": str(sp_meta),
            },
            metadata={}
        )
        
        assert mgr.is_step3_complete(s3_dir, mock_config) is False

    def test_is_step3_complete_missing_sp_metadata(self, tmp_path: Path, mock_config):
        mgr = CheckpointManager(tmp_path)
        mgr.initialize_state(product_smiles="C", config=mock_config)
        
        s3_dir = tmp_path / "S3_TransitionAnalysis"
        s3_dir.mkdir(parents=True, exist_ok=True)
        
        ts_final = s3_dir / "ts_final.xyz"
        _write_valid_xyz(ts_final)
        
        mgr.mark_step_completed(
            "s3",
            output_files={
                "ts_final_xyz": str(ts_final),
                "sp_matrix_metadata_json": str(s3_dir / "nonexistent.json"),
            },
            metadata={}
        )
        
        assert mgr.is_step3_complete(s3_dir, mock_config) is False

    def test_compute_step3_signature(self, tmp_path: Path, mock_config):
        mgr = CheckpointManager(tmp_path)
        
        sig1 = mgr._compute_step3_signature(mock_config)
        sig2 = mgr._compute_step3_signature(mock_config)
        
        assert sig1 == sig2
        assert "version" in sig1
        assert "theory_optimization" in sig1
        assert "theory_single_point" in sig1
        assert "step3_reactant_opt" in sig1

    def test_compute_file_hash(self, tmp_path: Path):
        mgr = CheckpointManager(tmp_path)
        
        test_file = tmp_path / "test.txt"
        test_file.write_text("Hello, World!")
        
        hash1 = mgr.compute_file_hash(test_file)
        hash2 = mgr.compute_file_hash(test_file)
        
        assert hash1 is not None
        assert hash1 == hash2
        assert len(hash1) == 16

    def test_compute_file_hash_nonexistent(self, tmp_path: Path):
        mgr = CheckpointManager(tmp_path)
        
        hash_val = mgr.compute_file_hash(tmp_path / "nonexistent.txt")

        assert hash_val is None

    def test_compute_step2_signature_stable(self, tmp_path: Path, mock_config):
        mgr = CheckpointManager(tmp_path)
        product_xyz = tmp_path / "product.xyz"
        product_xyz.write_text(
            "4\nproduct\nC 0 0 0\nC 1 0 0\nC 0 1 0\nC 0 0 1\n",
            encoding="utf-8",
        )

        sig1 = mgr.compute_step2_signature(
            config=mock_config,
            product_xyz=product_xyz,
            forming_bonds=((0, 1), (2, 3)),
            reaction_profile="[4+3]_default",
            scan_config={"scan_start_distance": 3.2, "scan_end_distance": 1.8, "scan_steps": 20},
        )
        sig2 = mgr.compute_step2_signature(
            config=mock_config,
            product_xyz=product_xyz,
            forming_bonds=((2, 3), (1, 0)),
            reaction_profile="[4+3]_default",
            scan_config={"scan_start_distance": 3.2, "scan_end_distance": 1.8, "scan_steps": 20},
        )

        assert sig1 == sig2
        assert sig1["product_xyz_hash"]
        assert sig1["forming_bonds"] == [[0, 1], [2, 3]]

    def test_is_step3_complete_upstream_step2_signature_mismatch(self, tmp_path: Path, mock_config):
        mgr = CheckpointManager(tmp_path)
        mgr.initialize_state(product_smiles="C", config=mock_config)

        s3_dir = tmp_path / "S3_TransitionAnalysis"
        s3_dir.mkdir(parents=True, exist_ok=True)

        ts_final = s3_dir / "ts_final.xyz"
        _write_valid_xyz(ts_final)

        sp_meta = s3_dir / "sp_matrix_metadata.json"
        _write_sp_metadata(sp_meta)

        step3_sig = mgr._compute_step3_signature(mock_config)
        mgr.mark_step_completed(
            "s3",
            output_files={
                "ts_final_xyz": str(ts_final),
                "sp_matrix_metadata_json": str(sp_meta),
            },
            metadata={
                "step3_signature": step3_sig,
                "input_hashes": {
                    "ts_guess": "abc123",
                    "reactant": "def456",
                    "product": "ghi789",
                },
                "upstream_step2_signature": {"k": "v1"},
            },
        )

        assert mgr.is_step3_complete(
            s3_dir,
            mock_config,
            check_signature=True,
            input_hashes={"ts_guess": "abc123", "reactant": "def456", "product": "ghi789"},
            upstream_step2_signature={"k": "v2"},
        ) is False


class TestS3ResumeJson:
    def test_load_s3_resume_state_empty(self, tmp_path: Path):
        from rph_core.steps.step3_opt.ts_optimizer import TSOptimizer
        
        optimizer = TSOptimizer(config={})
        state = optimizer._load_s3_resume_state(tmp_path)
        
        assert state == {}

    def test_save_and_load_s3_resume_state(self, tmp_path: Path):
        from rph_core.steps.step3_opt.ts_optimizer import TSOptimizer
        
        optimizer = TSOptimizer(config={})
        
        test_state = {
            "ts_opt": {
                "completed": True,
                "optimized_xyz": "/path/to/ts.xyz",
                "l2_energy": -100.0,
            },
            "reactant_opt": {
                "completed": True,
                "optimized_xyz": "/path/to/reactant.xyz",
                "l2_energy": -99.0,
            }
        }
        
        optimizer._save_s3_resume_state(tmp_path, test_state)
        loaded_state = optimizer._load_s3_resume_state(tmp_path)
        
        assert loaded_state == test_state

    def test_verify_ts_result_valid(self, tmp_path: Path):
        from rph_core.steps.step3_opt.ts_optimizer import TSOptimizer
        from rph_core.utils.qc_task_runner import QCOptimizationResult
        
        optimizer = TSOptimizer(config={})
        
        ts_opt_dir = tmp_path / "ts_opt"
        ts_opt_dir.mkdir(parents=True, exist_ok=True)
        l2_sp_dir = ts_opt_dir / "L2_SP"
        l2_sp_dir.mkdir(parents=True, exist_ok=True)
        
        optimized_xyz = ts_opt_dir / "ts.xyz"
        _write_valid_xyz(optimized_xyz)
        
        out_file = l2_sp_dir / "output.out"
        out_file.write_text("ORCA output")
        
        mock_result = QCOptimizationResult(
            optimized_xyz=optimized_xyz,
            converged=True,
            imaginary_count=1,
        )
        
        assert optimizer._verify_ts_result(mock_result, ts_opt_dir) is True

    def test_verify_ts_result_missing_l2_sp(self, tmp_path: Path):
        from rph_core.steps.step3_opt.ts_optimizer import TSOptimizer
        from rph_core.utils.qc_task_runner import QCOptimizationResult
        
        optimizer = TSOptimizer(config={})
        
        ts_opt_dir = tmp_path / "ts_opt"
        ts_opt_dir.mkdir(parents=True, exist_ok=True)
        
        optimized_xyz = ts_opt_dir / "ts.xyz"
        _write_valid_xyz(optimized_xyz)
        
        mock_result = QCOptimizationResult(
            optimized_xyz=optimized_xyz,
            converged=True,
            imaginary_count=1,
        )
        
        assert optimizer._verify_ts_result(mock_result, ts_opt_dir) is False

    def test_verify_ts_result_does_not_reject_name_pattern(self, tmp_path: Path):
        from rph_core.steps.step3_opt.ts_optimizer import TSOptimizer
        from rph_core.utils.qc_task_runner import QCOptimizationResult

        optimizer = TSOptimizer(config={})

        ts_opt_dir = tmp_path / "ts_opt"
        ts_opt_dir.mkdir(parents=True, exist_ok=True)
        l2_sp_dir = ts_opt_dir / "L2_SP"
        l2_sp_dir.mkdir(parents=True, exist_ok=True)

        optimized_xyz = ts_opt_dir / "ts_guess_like_name.xyz"
        _write_valid_xyz(optimized_xyz)
        (l2_sp_dir / "output.out").write_text("ORCA output")

        mock_result = QCOptimizationResult(
            optimized_xyz=optimized_xyz,
            converged=True,
            imaginary_count=1,
        )

        assert optimizer._verify_ts_result(mock_result, ts_opt_dir) is True

    def test_verify_reactant_result_valid(self, tmp_path: Path):
        from rph_core.steps.step3_opt.ts_optimizer import TSOptimizer
        from rph_core.utils.qc_task_runner import QCOptimizationResult
        
        optimizer = TSOptimizer(config={})
        
        reactant_opt_dir = tmp_path / "reactant_opt" / "standard"
        reactant_opt_dir.mkdir(parents=True, exist_ok=True)
        l2_sp_dir = reactant_opt_dir / "L2_SP"
        l2_sp_dir.mkdir(parents=True, exist_ok=True)
        
        optimized_xyz = reactant_opt_dir / "reactant.xyz"
        _write_valid_xyz(optimized_xyz)
        
        out_file = l2_sp_dir / "output.out"
        out_file.write_text("ORCA output")
        
        mock_result = QCOptimizationResult(
            optimized_xyz=optimized_xyz,
            converged=True,
            imaginary_count=0,
        )
        
        assert optimizer._verify_reactant_result(mock_result, reactant_opt_dir) is True

    def test_verify_reactant_result_with_imaginary(self, tmp_path: Path):
        from rph_core.steps.step3_opt.ts_optimizer import TSOptimizer
        from rph_core.utils.qc_task_runner import QCOptimizationResult
        
        optimizer = TSOptimizer(config={})
        
        reactant_opt_dir = tmp_path / "reactant_opt" / "standard"
        reactant_opt_dir.mkdir(parents=True, exist_ok=True)
        l2_sp_dir = reactant_opt_dir / "L2_SP"
        l2_sp_dir.mkdir(parents=True, exist_ok=True)
        
        optimized_xyz = reactant_opt_dir / "reactant.xyz"
        _write_valid_xyz(optimized_xyz)
        
        mock_result = QCOptimizationResult(
            optimized_xyz=optimized_xyz,
            converged=True,
            imaginary_count=1,
        )
        
        assert optimizer._verify_reactant_result(mock_result, reactant_opt_dir) is False
