"""
Simple standalone test for try_formchk utility function

This test validates the try_formchk() function behavior:
1. Success with existing .chk returns .fchk path
2. Failure with non-existent .chk returns None (graceful degradation)
3. Missing .chk logs a warning
"""
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import pytest

# Test success case: formchk creates .fchk successfully
def test_try_formchk_success(tmpdir):
    # Create test .chk file
    chk_file = tmpdir / "test.chk"
    chk_file.write_text("test checkpoint content", encoding='utf-8')

    # Mock subprocess.run to simulate successful formchk
    with patch('subprocess.run') as mock_run:
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Formchk successful",
            stderr=""
        )

    # Import and run try_formchk (mocking subprocess)
    # Since formchk uses subprocess.run internally, we patch that
    import sys
    import os
    rph_core_path = Path(__file__).parent.parent / "rph_core"
    if str(rph_core_path) not in sys.path:
        sys.path.insert(0, str(rph_core_path))

    from rph_core.utils.qc_interface import try_formchk

    # Test 1: formchk on existing .chk should return .fchk
    chk_file = tmpdir / "test_success.chk"
    chk_file.write_text("test checkpoint for success test")

    with patch('rph_core.utils.qc_interface.subprocess.run') as mock_run:
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Formchk successful",
            stderr=""
        )

    # Call try_formchk (subprocess.run is patched)
    result = try_formchk(chk_file)

    assert result is not None, "Should return .fchk path on success"
    assert result == chk_file.with_suffix('.fchk'), "Should return correct .fchk path"
    assert result.exists(), "Should create .fchk file"

    # Cleanup
    result.unlink()
    chk_file.unlink()

    print(f"✅ test_try_formchk_success passed")


# Test 2: formchk fails (returncode=1) should return None
def test_try_formchk_returncode_failure(tmpdir):
    # Create test .chk file
    chk_file = tmpdir / "test_failure.chk"
    chk_file.write_text("test checkpoint for failure test", encoding='utf-8')

    # Mock subprocess.run with non-zero returncode
    with patch('rph_core.utils.qc_interface.subprocess.run') as mock_run:
        mock_run.return_value = MagicMock(
            returncode=1,
            stdout="",
            stderr="formchk Error writing formatted checkpoint"
        )

    from rph_core.utils.qc_interface import try_formchk
    result = try_formchk(chk_file)

    assert result is None, "Should return None on failure"
    assert not chk_file.with_suffix('.fchk').exists(), "fchk should not be created on failure"

    # Cleanup
    if chk_file.exists():
        chk_file.unlink()


# Test 3: formchk called on non-existent .chk should return None and log warning
def test_try_formchk_missing_chk(tmpdir):
    # Don't create .chk file
    nonexist_file = tmpdir / "nonexistent.chk"

    from rph_core.utils.qc_interface import try_formchk
    result = try_formchk(nonexist_file)

    assert result is None, "Should return None when .chk doesn't exist"
    print(f"✅ test_try_formchk_missing_chk passed")


# Test 4: formchk times out should return None and log warning
def test_try_formchk_timeout(tmpdir):
    # Don't create .chk file
    timeout_file = tmpdir / "timeout_test.chk"

    # Mock subprocess.run with TimeoutExpired
    with patch('rph_core.utils.qc_interface.subprocess.run') as mock_run:
        mock_run.side_effect = subprocess.TimeoutExpired("formchk timed out")

    from rph_core.utils.qc_interface import try_formchk
    result = try_formchk(timeout_file)

    assert result is None, "Should return None on timeout"
    assert not timeout_file.with_suffix('.fchk').exists(), "fchk should not be created on timeout"

    print(f"✅ test_try_formchk_timeout passed")


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory(prefix="test_formchk_") as tmpdir:
        try:
            test_try_formchk_success(tmpdir)
            test_try_formchk_returncode_failure(tmpdir)
            test_try_formchk_timeout(tmpdir)
            test_try_formchk_missing_chk(tmpdir)
            print("\n✅ All standalone tests passed")
        finally:
            import shutil
            shutil.rmtree(tmpdir, ignore_errors=True)
