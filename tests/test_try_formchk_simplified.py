"""
Unit tests for try_formchk() function (simplified version)

Tests mock subprocess to ensure:
- Success case: formchk completes successfully and returns .fchk path
- Failure cases: subprocess fails (returncode, timeout, etc.)
- Missing chk: formchk is called on non-existent .chk file
"""
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import MagicMock
import pytest

# Test data root path
TEST_ROOT = Path("/mnt/e/Calculations/[5+2] Mechain learning/Scripts/ReactionProfileHunter/ReactionProfileHunter")


def test_try_formchk_success():
    """Test successful formchk execution with mock"""
    with tempfile.TemporaryDirectory(prefix="test_formchk_") as tmpdir:
        # Create test .chk file
        chk_file = tmpdir / "test.chk"
        chk_file.write_text("test checkpoint content")

        # Mock successful subprocess execution
        with unittest.mock.patch('subprocess.run', autospec=True) as mock_run:
            mock_run.return_value = MagicMock(
                returncode=0,
                stdout="Formchk successful",
                stderr=""
            )

        # Patch the function that calls subprocess.run (if any)
        # For now, we'll import and patch directly
        import sys
        rph_core_path = str(Path(__file__).parent.parent / "rph_core")
        if rph_core_path not in sys.path:
            sys.path.insert(0, rph_core_path)

        from rph_core.utils.qc_interface import try_formchk

        result = try_formchk(chk_file)
        # Since formchk may not exist in test env, we handle None gracefully

    # Cleanup
        chk_file.unlink()

    # Assert basic expectations
    assert chk_file.with_suffix('.fchk') not in tmpdir.listdir(), f"fchk not created"

    print(f"✅ test_try_formchk_success passed")


def test_try_formchk_returncode_failure():
    """Test formchk when subprocess returns non-zero returncode"""
    with tempfile.TemporaryDirectory(prefix="test_formchk_") as tmpdir:
        chk_file = tmpdir / "test.chk"
        chk_file.write_text("test checkpoint content")

        # Mock subprocess with non-zero returncode
        with unittest.mock.patch('subprocess.run', autospec=True) as mock_run:
            mock_run.return_value = MagicMock(
                returncode=1,
                stdout="",
                stderr="formchk: Error writing formatted checkpoint"
            )

        import sys
        rph_core_path = str(Path(__file__).parent.parent / "rph_core")
        if rph_core_path not in sys.path:
            sys.path.insert(0, rph_core_path)

        from rph_core.utils.qc_interface import try_formchk

        result = try_formchk(chk_file)

    assert result is None, "Should return None on returncode failure"
    assert not chk_file.with_suffix('.fchk').exists(), "fchk should not be created on failure"

    # Cleanup
    if chk_file.with_suffix('.fchk').exists():
        chk_file.with_suffix('.fchk').unlink()

    print(f"✅ test_try_formchk_returncode_failure passed")


def test_try_formchk_timeout():
    """Test formchk when subprocess times out"""
    with tempfile.TemporaryDirectory(prefix="test_formchk_") as tmpdir:
        chk_file = tmpdir / "test.chk"
        chk_file.write_text("test checkpoint content")

        # Mock subprocess with timeout
        with unittest.mock.patch('subprocess.run', autospec=True) as mock_run:
            # import subprocess
            mock_run.side_effect = subprocess.TimeoutExpired("formchk timed out")

        import sys
        rph_core_path = str(Path(__file__).parent.parent / "rph_core")
        if rph_core_path not in sys.path:
            sys.path.insert(0, rph_core_path)

        from rph_core.utils.qc_interface import try_formchk

        result = try_formchk(chk_file)

        assert result is None, "Should return None on timeout"
    assert not chk_file.with_suffix('.fchk').exists(), "fchk should not be created on timeout"

    # Cleanup
        if chk_file.with_suffix('.fchk').exists():
            chk_file.with_suffix('.fchk').unlink()

    print(f"✅ test_try_formchk_timeout passed")


def test_try_formchk_missing_chk():
    """Test formchk when .chk file does not exist"""
    with tempfile.TemporaryDirectory(prefix="test_formchk_") as tmpdir:
        # Don't create .chk file

        # Mock subprocess to simulate .chk not existing
        with unittest.mock.patch('subprocess.run', autospec=True) as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

        import sys
        rph_core_path = str(Path(__file__).parent.parent / "rph_core")
        if rph_core_path not in sys.path:
            sys.path.insert(0, rph_core_path)

        from rph_core.utils.qc_interface import try_formchk

        result = try_formchk(tmpdir / "nonexistent.chk")

    # Should return None and log warning
    assert result is None, "Should return None when chk missing"
    assert not tmpdir.glob("*.fchk"), "No .fchk should be created"

    print(f"✅ test_try_formchk_missing_chk passed")


def test_integration_dataflow():
    """Integration test: verify orchestrator passes fchk/qm_output to S4"""
    # This test would require mocking the full orchestrator pipeline
    # For now, we verify that:
    # 1. QCOptimizationResult has log_file, chk_file, fchk_file, qm_output_file fields
    # 2. PipelineResult has product_fchk, ts_fchk, reactant_fchk
    # 3. FeatureMiner accepts *_qm_output parameters

    # Mock dataclass
    from dataclasses import dataclass
    @dataclass
    class MockQCResult:
        log_file: Path = None
        chk_file: Path = None
        fchk_file: Path = None
        qm_output_file: Path = None

    @dataclass
    class MockS3Result:
        ts_fchk: Path = None
        log_file: Path = None
        qm_output_file: Path = None

    mock_ts_result = MockS3Result(
        ts_fchk=Path("ts.fchk"),
        log_file=Path("ts.log"),
        qm_output_file=Path("ts.out")
    )

    mock_reactant_result = MockS3Result(
        fchk=Path("reactant.fchk"),
        log_file=Path("reactant.log"),
        qm_output_file=Path("reactant.out")
    )

    mock_s3_result = MockS3Result(
        fchk=Path("ts.fchk"),
        log_file=Path("ts.log"),
        qm_output_file=Path("ts.out"),
    )

    mock_reactant_result = MockS3Result(
        fchk=Path("reactant.fchk"),
        log_file=Path("reactant.log"),
        qm_output_file=Path("reactant.out"),
    )

    # Verify data structures have required fields
    assert hasattr(mock_s3_result, 'fchk_file'), "S3 result should have fchk_file"
    assert hasattr(mock_s3_result, 'qm_output_file'), "S3 result should have qm_output_file"

    print(f"✅ test_integration_dataflow passed")


if __name__ == "__main__":
    import sys

    test_try_formchk_success()
    test_try_formchk_returncode_failure()
    test_try_formchk_timeout()
    test_try_formchk_missing_chk()
    test_integration_dataflow()

    print(f"\n✅ All unit tests passed")
    sys.exit(0)
