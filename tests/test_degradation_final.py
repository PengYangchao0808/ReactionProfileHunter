"""
Simplified degradation test: Verify warning logging when fchk is missing

This test validates that:
1. Features.csv is generated even when .fchk files are missing
2. Appropriate warnings are logged for missing artifacts
3. Pipeline doesn't crash
4. Warnings are aggregated into warnings_count
"""
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock, call
import sys

# Import rph_core modules
rph_core_path = Path(__file__).parent.parent / "rph_core"
if str(rph_core_path) not in sys.path:
    sys.path.insert(0, str(rph_core_path))

# Import orchestrator
from rph_core.orchestrator import ReactionProfileHunter
from rph_core.steps.step4_features.feature_miner import FeatureMiner


@dataclass
class MockQCResult:
    """Mock QCResult with new artifact fields"""
    energy: float = 0.0
    converged: bool = False
    log_file: Path = None
    chk_file: Path = None
    fchk_file: Path = None
    qm_output_file: Path = None


@dataclass
class MockS3Result:
    """Mock S3 TransitionAnalysisResult with new artifact fields"""
    ts_final_xyz: Path = Path("/fake/ts_final.xyz")
    ts_checkpoint: Path = None
    sp_report: MagicMock()
    method_used: str = "Berny"
    ts_fchk: Path = None
    ts_log: Path = None
    ts_qm_output: Path = None
    reactant_fchk: Path = None
    reactant_log: Path = None
    reactant_qm_output: Path = None


def test_missing_fchk_logging():
    """Test that missing .chk files trigger warning logging"""
    with tempfile.TemporaryDirectory(prefix="test_missing_fchk_") as tmpdir:
        # Import
        rph_core_path = Path(__file__).parent.parent / "rph_core"
        if str(rph_core_path) not in sys.path:
            sys.path.insert(0, str(rph_core_path))

        # Mock FeatureMiner with warning tracking
        with patch.object(FeatureMiner, '__init__', autospec=False) as mock_miner:
            mock_miner = MagicMock(spec=FeatureMiner)
            mock_miner.warning = MagicMock()

        # Create temporary directory to work in
        work_dir = tmpdir / "test_missing_fchk_workdir"

        # Patch subprocess.run to simulate formchk failure
        with patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stderr="formchk failed")

        # Import and call FeatureMiner
        from rph_core.steps.step4_features.feature_miner import FeatureMiner
        from rph_core.orchestrator import ReactionProfileHunter

        from rph_core.utils.log_manager import setup_logger
        setup_logger(name=__name__, file=None)
        feature_miner = FeatureMiner(config_path=None, work_dir=work_dir)

        # Call with non-existent .chk file (should log warning)
        nonexistent_chk = work_dir / "nonexistent.chk"

        result = feature_miner.run(
            ts_final=work_dir / "ts_final.xyz",
            reactant=work_dir / "reactant.xyz",
            product=work_dir / "product_min.xyz",
            output_dir=work_dir / "features_output",
            sp_matrix_report=MagicMock(),
            forming_bonds=None
        )

        # Verify warning was called
        assert len(mock_miner.warning.call_args) > 0, "Should log warning for missing .chk"

        print("✅ test_missing_fchk_logging passed")


def test_degradation_features_csv_generation():
    """Test that features_raw.csv is generated even with missing fchk"""
    with tempfile.TemporaryDirectory(prefix="test_degradation_") as tmpdir:
        # Import and mock FeatureMiner
        rph_core_path = Path(__file__).parent.parent / "rph_core"
        if str(rph_core_path) not in sys.path:
            sys.path.insert(0, str(rph_core_path))

        from rph_core.steps.step4_features.feature_miner import FeatureMiner
        from rph_core.orchestrator import ReactionProfileHunter
        from rph_core.utils.log_manager import setup_logger
        setup_logger(name=__name__, file=None)

        # Create mock FeatureMiner
        mock_miner = MagicMock(spec=FeatureMiner)
        mock_miner.warning = MagicMock()

        # Create work directory
        work_dir = tmpdir / "test_degradation_workdir"

        # Create mock S3Result with all new artifact fields
        mock_s3_result = TransitionAnalysisResult(
            ts_final_xyz=work_dir / "ts_final.xyz",
            ts_fchk=work_dir / "ts.fchk",
            ts_log=work_dir / "ts.log",
            ts_qm_output=work_dir / "ts.out",
            reactant_fchk=work_dir / "reactant.fchk",
            reactant_log=work_dir / "reactant.log",
            reactant_qm_output=work_dir / "reactant.out",
            sp_report=MagicMock()
        )

        # Create mock PipelineResult with all new fields
        mock_hunter_result = PipelineResult(
            success=True,
            product_smiles="C=C(O)C",
            product_xyz=Path("S1_Product/product_min.xyz"),
            product_fchk=Path("S1_Product/product.fchk"),
            product_log=Path("S1_Product/product.log"),
            product_qm_output=Path("S1_Product/product.out"),
            ts_guess_xyz=Path("S2_Retro/ts_guess.xyz"),
            reactant_xyz=Path("S2_Retro/reactant.xyz"),
            ts_final_xyz=mock_s3_result["ts_final_xyz"],
            ts_fchk=mock_s3_result["ts_fchk"],
            ts_log=mock_s3_result["ts_log"],
            ts_qm_output=mock_s3_result["ts_qm_output"],
            reactant_fchk=mock_s3_result["reactant_fchk"],
            reactant_log=mock_s3_result["reactant_log"],
            reactant_qm_output=mock_s3_result["reactant_qm_output"],
            forming_bonds=((1, 2), (3, 4)),
            features_csv=Path("S4_Data/features_raw.csv")
        )

        # Mock orchestrator
        with patch.object(ReactionProfileHunter, '__init__', autospec=False) as mock_hunter:
            # Create instance without calling __init__
            mock_hunter = MagicMock(spec=ReactionProfileHunter)
            ReactionProfileHunter.__init__(mock_hunter)

            # Mock s4_engine.run() to capture parameters and return features_raw.csv
            captured_s4_params = {}
            def capture_s4_run(*args, **kwargs):
                captured_params.update(kwargs)
                return Path("S4_Data/features_raw.csv")

            mock_hunter.s4_engine.run = capture_s4_run

            # Mock Step 1-3 results with new artifact fields
            mock_s1_result = {
                "xyz": Path("S1_Product/product_min.xyz"),
                "e_sp": 0.0,
                "log": Path("S1_Product/product.log"),
                "chk": Path("S1_Product/product.fchk"),
                "qm_output": Path("S1_Product/product.out")
            }
            mock_s3_result = {
                "xyz": Path("S2_Retro/ts_guess.xyz"),
                "e_sp": -0.5,
                "log": Path("S2_Retro/ts.log"),
                "chk": Path("S2_Retro/ts.fchk"),
                "qm_output": Path("S2_Retro/ts.out"),
                }
            mock_hunter.run_pipeline.return_value = MagicMock(
                success=True,
                product_smiles="C=C(O)C",
                product_xyz=mock_s1_result["xyz"],
                product_fchk=mock_s1_result["fchk"],
                product_log=mock_s1_result["log"],
                product_qm_output=mock_s1_result["qm_output"],
                ts_guess_xyz=mock_s3_result["xyz"],
                ts_final_xyz=mock_s3_result["ts_final_xyz"],
                ts_fchk=mock_s3_result["ts_fchk"],
                ts_log=mock_s3_result["ts_log"],
                ts_qm_output=mock_s3_result["ts_qm_output"],
                reactant_fchk=mock_s3_result["reactant_fchk"],
                reactant_log=mock_s3_result["reactant_log"],
                reactant_qm_output=mock_s3_result["reactant_qm_output"],
                sp_matrix_report=mock_s3_result["sp_report"],
                forming_bonds=((1, 2), (3, 4)),
                features_csv=Path("S4_Data/features_raw.csv")
            )

        # Call orchestrator with skip steps=['s2', 's3', 's4'] to isolate Step 4
        # Mock ReactionProfileHunter.__init__ to return dummy pipeline
        mock_pipeline_result = MagicMock(
            product_smiles="C=C(O)C",
            product_xyz=Path("S1_Product/product_min.xyz"),
            ts_guess_xyz=Path("S2_Retro/ts_guess.xyz"),
            reactant_xyz=Path("S2_Retro/reactant.xyz"),
            ts_final_xyz=mock_s3_result["ts_final_xyz"],
            ts_fchk=mock_s3_result["ts_fchk"],
            ts_log=mock_s3_result["ts_log"],
            ts_qm_output=mock_s3_result["ts_qm_output"],
            reactant_fchk=mock_s3_result["reactant_fchk"],
            reactant_log=mock_s3_result["reactant_log"],
            reactant_qm_output=mock_s3_result["reactant_qm_output"],
            sp_matrix_report=mock_s3_result["sp_report"],
            forming_bonds=((1, 2), (3, 4)),
            features_csv=Path("S4_Data/features_raw.csv")
        )

        # Call orchestrator
        result = mock_hunter.run_pipeline(
            work_dir=Path.cwd() / "test_integration"
        )

        # Verify features_raw.csv exists
        features_csv = Path.cwd() / "test_integration/S4_Data/features_raw.csv"
        assert features_csv.exists(), "features_raw.csv should be created"

        # Extract warning count from features_raw.csv
        with features_csv.open('r') as f:
            content = f.read()
            # Count warning strings
            warning_count = content.count('WARNING') + content.count('⚠️')

        assert warning_count >= 3, f"Expected at least 3 warnings logged"
        print(f"✅ test_degradation_features_csv_generation passed (warnings logged: {warning_count})")


if __name__ == "__main__":
    test_missing_fchk_logging()
    test_degradation_features_csv_generation()
    print("\n✅ All degradation tests passed")
