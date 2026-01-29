"""
Import Smoke Tests for Step 4 Features and Step 3
===================================================

Fast import validation tests to catch import chain errors early.
These tests verify that Step 4 and Step 3 modules can be imported successfully,
preventing ModuleNotFoundError during pipeline execution.

Purpose: Early detection of import chain errors before pipeline runs.
Speed: Very fast (sub-second execution).
Coverage: Step 4 feature extraction modules, Step 3 TS optimizer.

Author: QC Descriptors Team
Date: 2026-01-19
"""


def test_import_step4_feature_miner():
    """Test that FeatureMiner can be imported."""
    import rph_core.steps.step4_features.feature_miner  # noqa: F401


def test_import_step4_extractors():
    """Test that all extractor modules can be imported."""
    import rph_core.steps.step4_features.extractors.geometry  # noqa: F401
    import rph_core.steps.step4_features.extractors.thermo  # noqa: F401
    import rph_core.steps.step4_features.extractors.qc_checks  # noqa: F401
    import rph_core.steps.step4_features.extractors.interaction_analysis  # noqa: F401
    import rph_core.steps.step4_features.extractors.nics  # noqa: F401
    import rph_core.steps.step4_features.extractors.nbo_e2  # noqa: F401


def test_import_step4_context():
    """Test that context module can be imported."""
    import rph_core.steps.step4_features.context  # noqa: F401


def test_import_step4_path_accessor():
    """Test that path_accessor module can be imported."""
    import rph_core.steps.step4_features.path_accessor  # noqa: F401


def test_import_step4_base():
    """Test that base extractor module can be imported."""
    import rph_core.steps.step4_features.extractors.base  # noqa: F401


def test_import_step3_ts_optimizer():
    """Test that TSOptimizer can be imported (Step 3 was fixed for import issues)."""
    import rph_core.steps.step3_opt.ts_optimizer  # noqa: F401


def test_import_step3_sp_matrix_report():
    """Test that SPMatrixReport class can be imported from Step 3."""
    from rph_core.steps.step3_opt.ts_optimizer import SPMatrixReport  # noqa: F401
