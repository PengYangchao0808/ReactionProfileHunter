from pathlib import Path

from rph_core.orchestrator import PipelineResult


def test_pipeline_result_supports_artifact_fields() -> None:
    result = PipelineResult(success=True, product_smiles="C=C", work_dir=Path("/tmp"))
    result.product_fchk = Path("/tmp/product.fchk")
    result.ts_fchk = Path("/tmp/ts.fchk")
    result.reactant_fchk = Path("/tmp/reactant.fchk")
    result.product_log = Path("/tmp/product.log")
    result.ts_log = Path("/tmp/ts.log")
    result.reactant_log = Path("/tmp/reactant.log")
    result.product_qm_output = Path("/tmp/product.out")
    result.ts_qm_output = Path("/tmp/ts.out")
    result.reactant_qm_output = Path("/tmp/reactant.out")

    assert result.success is True
    assert result.product_fchk.name == "product.fchk"
    assert result.ts_qm_output.name == "ts.out"
