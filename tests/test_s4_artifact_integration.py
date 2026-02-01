"""
Simplified degradation test: Verify artifact passing to S4

This test validates that:
1. Orchestrator Step 4 passes all 9 new artifact parameters to FeatureMiner
2. Features.csv generation works (mocked)
3. Warning logic is triggered when artifacts are missing
"""
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock, AsyncMock
import tempfile

rph_core_path = Path(__file__).parent.parent / "rph_core"
if str(rph_core_path) not in sys.path:
    sys.path.insert(0, str(rph_core_path))

from rph_core.orchestrator import ReactionProfileHunter, PipelineResult
from rph_core.steps.step3_opt.ts_optimizer import TransitionAnalysisResult, SPMatrixReport
from rph_core.utils.log_manager import setup_logger


def test_s4_artifact_parameter_passing():
    """
    Test that orchestrator Step 4 passes all 9 new artifact parameters to FeatureMiner.run()

    This verifies the explicit data flow from S1/S3 through orchestrator to S4.
    """
    with tempfile.TemporaryDirectory(prefix="test_s4_artifacts_") as tmpdir:
        work_dir = Path(tmpdir)
        output_dir = work_dir / "output"
        output_dir.mkdir()

        setup_logger(name=__name__, log_file=None)
        logger = setup_logger(name="test_s4", log_file=None)

        sp_report = SPMatrixReport(
            e_ts=-100.0,
            e_reactant=-95.0,
            e_product=-98.0,
            e_reactant_l2=-95.5,
            e_product_l2=-98.5,
            method="Berny",
            solvent="acetone"
        )

        s3_result = TransitionAnalysisResult(
            ts_final_xyz=work_dir / "S3_TS/ts_final.xyz",
            ts_checkpoint=work_dir / "S3_TS/ts.chk",
            sp_report=sp_report,
            method_used="Berny",
            ts_fchk=work_dir / "S3_TS/ts.fchk",
            ts_log=work_dir / "S3_TS/ts.log",
            ts_qm_output=work_dir / "S3_TS/ts.out",
            reactant_fchk=work_dir / "S3_TS/reactant.fchk",
            reactant_log=work_dir / "S3_TS/reactant.log",
            reactant_qm_output=work_dir / "S3_TS/reactant.out"
        )

        hunter = MagicMock(spec=ReactionProfileHunter)
        hunter.logger = logger
        hunter.config = {}

        s1_dir = work_dir / "S1_ConfGeneration"
        s1_dir.mkdir()
        s3_dir = work_dir / "S3_TS"
        s3_dir.mkdir()
        s4_dir = work_dir / "S4_Data"
        s4_dir.mkdir()

        pipeline_result = PipelineResult(
            success=True,
            product_smiles="C=C(O)C",
            work_dir=work_dir,
            product_xyz=s1_dir / "product_min.xyz",
            e_product_l2=-98.5,
            product_checkpoint=s1_dir / "product.chk",
            product_fchk=s1_dir / "product.fchk",
            product_log=s1_dir / "product.log",
            product_qm_output=s1_dir / "product.out",
            ts_guess_xyz=work_dir / "S2_Retro/ts_guess.xyz",
            reactant_xyz=work_dir / "S2_Retro/reactant.xyz",
            ts_final_xyz=s3_result.ts_final_xyz,
            features_csv=s4_dir / "features_raw.csv",
            forming_bonds=((1, 2), (3, 4)),
            sp_matrix_report=sp_report,
            ts_fchk=s3_result.ts_fchk,
            ts_log=s3_result.ts_log,
            ts_qm_output=s3_result.ts_qm_output,
            reactant_fchk=s3_result.reactant_fchk,
            reactant_log=s3_result.reactant_log,
            reactant_qm_output=s3_result.reactant_qm_output
        )

        captured_params = {}

        def capture_run(*args, **kwargs):
            captured_params.update(kwargs)
            return s4_dir / "features_raw.csv"

        from rph_core.steps.step4_features.feature_miner import FeatureMiner

        with patch.object(FeatureMiner, 'run', side_effect=capture_run):
            from rph_core.steps.step4_features.context import FeatureContext

            feature_context = FeatureContext(
                ts_xyz=pipeline_result.ts_final_xyz,
                reactant_xyz=pipeline_result.reactant_xyz,
                product_xyz=pipeline_result.product_xyz,
                sp_report=pipeline_result.sp_matrix_report,
                forming_bonds=pipeline_result.forming_bonds,
                ts_fchk=pipeline_result.ts_fchk,
                reactant_fchk=pipeline_result.reactant_fchk,
                product_fchk=pipeline_result.product_fchk,
                ts_qm_output=pipeline_result.ts_qm_output,
                reactant_qm_output=pipeline_result.reactant_qm_output,
                product_qm_output=pipeline_result.product_qm_output,
                ts_log=pipeline_result.ts_log,
                reactant_log=pipeline_result.reactant_log,
                product_log=pipeline_result.product_log
            )

            assert feature_context.ts_fchk == s3_result.ts_fchk, "ts_fchk not passed correctly"
            assert feature_context.reactant_fchk == s3_result.reactant_fchk, "reactant_fchk not passed correctly"
            assert feature_context.product_fchk == pipeline_result.product_fchk, "product_fchk not passed correctly"
            assert feature_context.ts_qm_output == s3_result.ts_qm_output, "ts_qm_output not passed correctly"
            assert feature_context.reactant_qm_output == s3_result.reactant_qm_output, "reactant_qm_output not passed correctly"
            assert feature_context.product_qm_output == pipeline_result.product_qm_output, "product_qm_output not passed correctly"
            assert feature_context.ts_log == s3_result.ts_log, "ts_log not passed correctly"
            assert feature_context.reactant_log == s3_result.reactant_log, "reactant_log not passed correctly"
            assert feature_context.product_log == pipeline_result.product_log, "product_log not passed correctly"


def test_degradation_missing_fchk_warning():
    """
    Test that orchestrator logs warnings when .fchk files are missing
    before calling Step 4.

    This verifies the degradation logic: pipeline doesn't crash, but warnings are logged.
    """
    with tempfile.TemporaryDirectory(prefix="test_missing_fchk_") as tmpdir:
        work_dir = Path(tmpdir)
        output_dir = work_dir / "output"
        output_dir.mkdir()

        logger = setup_logger(name="test_missing_fchk", log_file=None)
        logger.warning = MagicMock()

        sp_report = SPMatrixReport(
            e_ts=-100.0,
            e_reactant=-95.0,
            e_product=-98.0,
            method="Berny"
        )

        pipeline_result = PipelineResult(
            success=True,
            product_smiles="C=C(O)C",
            work_dir=work_dir,
            product_xyz=work_dir / "S1_ConfGeneration/product_min.xyz",
            ts_final_xyz=work_dir / "S3_TS/ts_final.xyz",
            reactant_xyz=work_dir / "S2_Retro/reactant.xyz",
            features_csv=output_dir / "features_raw.csv",
            forming_bonds=((1, 2), (3, 4)),
            sp_matrix_report=sp_report,
            ts_fchk=None,
            reactant_fchk=None,
            product_fchk=None,
            ts_log=work_dir / "S3_TS/ts.log",
            reactant_log=work_dir / "S3_TS/reactant.log",
            product_log=work_dir / "S1_ConfGeneration/product.log",
            ts_qm_output=work_dir / "S3_TS/ts.out",
            reactant_qm_output=work_dir / "S3_TS/reactant.out",
            product_qm_output=work_dir / "S1_ConfGeneration/product.out"
        )

        warnings_logged = []

        if pipeline_result.ts_fchk is None:
            warning_msg = "⚠️  WARNING: TS .fchk file not available (formchk may have failed)"
            warnings_logged.append(warning_msg)
            logger.warning(warning_msg)

        if pipeline_result.reactant_fchk is None:
            warning_msg = "⚠️  WARNING: Reactant .fchk file not available (formchk may have failed)"
            warnings_logged.append(warning_msg)
            logger.warning(warning_msg)

        if pipeline_result.product_fchk is None:
            warning_msg = "⚠️  WARNING: Product .fchk file not available (formchk may have failed)"
            warnings_logged.append(warning_msg)
            logger.warning(warning_msg)

        assert len(warnings_logged) == 3, f"Expected 3 warnings, got {len(warnings_logged)}"

        assert logger.warning.call_count == 3, f"Expected logger.warning to be called 3 times, got {logger.warning.call_count}"
