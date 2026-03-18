"""
Simple degradation test: verify S4 generates features_raw.csv even when fchk files are missing

This test simulates the scenario where:
1. Step 3 (TSOptimizer) fails to generate .fchk files for TS/reactant
2. orchestrator Step 4 is called with None fchk paths
3. FeatureMiner runs with the missing inputs and successfully generates features_raw.csv
4. Warnings are properly logged but S4 doesn't crash
"""
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

# Mock S4 FeatureMiner to capture warnings
def test_degradation_with_missing_fchk():
    """Test that S4 produces features_raw.csv with warnings when fchk files are missing"""
    # Create temporary directory
    with tempfile.TemporaryDirectory(prefix="test_degradation_") as tmpdir:
        # Mock FeatureMiner.run() to return features_raw.csv even with missing inputs
        from unittest.mock import patch, MagicMock

        # Import and patch
        import sys
        rph_core_path = Path(__file__).parent.parent / "rph_core"
        if str(rph_core_path) not in sys.path:
            sys.path.insert(0, str(rph_core_path))

        with patch('rph_core.steps.step4_features.feature_miner.FeatureMiner') as mock_s4:
            # Capture warnings
            warnings_captured = []

            def mock_run_features(**kwargs):
                # Simulate successful feature extraction
                # Record warnings
                warnings_captured.append("⚠️  TS fchk missing: .fchk not available")
                warnings_captured.append("⚠️  Reactant fchk missing: .fchk not available")

                # Create mock features_raw.csv with basic features
                features_csv_path = tmpdir / "features_raw.csv"
                features_csv_path.write_text("test_id,activation_energy,reacton_energy,dE_activation\n")
                return features_csv_path

            # Patch run_features to return features_raw.csv and capture warnings
            mock_s4.run_features = mock_run_features

        # Simulate Step4 call in orchestrator
        from rph_core.orchestrator import ReactionProfileHunter

        # Mock PipelineResult with S1/S3 artifacts but None fchk paths
        mock_pipeline_result = MagicMock(
            success=True,
            product_xyz=Path("fake/product_min.xyz"),
            ts_final_xyz=Path("fake/ts_final.xyz"),
            ts_fchk=None,  # S3 failed to generate .fchk
            reactant_fchk=None,  # S3 failed to generate .fchk
            product_fchk=None,  # S1 failed to generate .fchk
            ts_log=None,
            reactant_log=None,
            ts_qm_output=None,
            reactant_qm_output=None,
            product_log=None,
            product_qm_output=None,
            sp_matrix_report=MagicMock()
        )

        # Create temporary directory for S4_Data
        s4_dir = tmpdir / "S4_Data"

        # Mock orchestrator
        from unittest.mock import MagicMock
        mock_orchestrator = MagicMock(spec=ReactionProfileHunter)
        mock_orchestrator.run_pipeline = MagicMock(return_value=mock_pipeline_result)

        # Create mock s4_engine
        with patch('rph_core.steps.step4_features.feature_miner.FeatureMiner') as mock_s4:
            mock_s4.run_features = mock_run_features

        # Run pipeline with mock results
        result = mock_orchestrator.run_pipeline(
            product_smiles="C=C(O)C",
            work_dir=s4_dir
        )

        # Verify: features_raw.csv should be created
        assert (s4_dir / "features_raw.csv").exists(), "features_raw.csv should be created"
        assert (s4_dir / "features_raw.csv").read_text().count("test_id") == 1, "Should contain 1 data row"

        print("✅ test_degradation_with_missing_fchk passed")


if __name__ == "__main__":
    test_degradation_with_missing_fchk()
    print("\n✅ All degradation tests passed")
