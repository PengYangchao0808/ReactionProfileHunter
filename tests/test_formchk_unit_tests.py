"""
Unit tests for try_formchk() function (simplified)

Tests mock subprocess to ensure:
- Success case: formchk completes successfully and returns .fchk path
- Failure cases: subprocess fails (returncode, timeout, etc.) and returns None
- Missing chk: formchk is called on non-existent .chk file
"""
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import pytest
import sys

# Test in current directory
TEST_ROOT = Path.cwd() / "test_tmpdir"
TEST_ROOT.mkdir(parents=True, exist_ok=True)


def test_try_formchk_success():
    """Test successful formchk execution"""
    # Create a mock .chk file
    chk_file = TEST_ROOT / "test.chk"
    chk_file.write_text("test checkpoint content")

    # Mock successful subprocess execution
    with patch('subprocess.run') as mock_run:
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Formchk successful",
            stderr=""
        )
        mock_run.side_effect = None

    # Import and run try_formchk
    from rph_core.utils.qc_interface import try_formchk

    result = try_formchk(chk_file)

    assert result is not None, "Should return fchk path on success"
    assert result == chk_file.with_suffix('.fchk'), f"Should return correct fchk path"
    assert result.exists(), "fchk file should exist"

    # Cleanup
    chk_file.unlink()

    print("✅ test_try_formchk_success passed")


def test_try_formchk_returncode_failure():
    """Test formchk when subprocess returns non-zero return code"""
    chk_file = TEST_ROOT / "test.chk"
    chk_file.write_text("test checkpoint content")

    # Mock subprocess with non-zero returncode
    with patch('subprocess.run') as mock_run:
        mock_run.return_value = MagicMock(
            returncode=1,
            stdout="",
            stderr="formchk: Error writing formatted checkpoint"
        )

    from rph_core.utils.qc_interface import try_formchk
    result = try_formchk(chk_file)

    assert result is None, "Should return None on failure"
    assert not chk_file.with_suffix('.fchk').exists(), "fchk should not be created on failure"

    # Cleanup
    if chk_file.with_suffix('.fchk').exists():
        chk_file.with_suffix('.fchk').unlink()


def test_try_formchk_timeout():
    """Test formchk when subprocess times out"""
    chk_file = TEST_ROOT / "test.chk"
    chk_file.write_text("test checkpoint content")

    # Mock subprocess with timeout
    with patch('subprocess.run') as mock_run:
        mock_run.side_effect = subprocess.TimeoutExpired("formchk timed out")

    from rph_core.utils.qc_interface import try_formchk
    result = try_formchk(chk_file)

    assert result is None, "Should return None on timeout"
    assert not chk_file.with_suffix('.fchk').exists(), "fchk should not be created on timeout"

    # Cleanup
    if chk_file.with_suffix('.fchk').exists():
        chk_file.with_suffix('.fchk').unlink()


def test_try_formchk_missing_chk():
    """Test formchk when .chk file does not exist"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Don't create .chk file

        # Mock subprocess to simulate .chk not existing
        with patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

        import sys
        rph_core_path = str(Path(__file__).parent.parent / "rph_core")
        if rph_core_path not in sys.path:
            sys.path.insert(0, rph_core_path)

        from rph_core.utils.qc_interface import try_formchk

        result = try_formchk(tmpdir / "nonexistent.chk")

        # Should return None and log warning
        assert result is None, "Should return None when chk missing"
        assert not list(tmpdir.glob("*.fchk")), "No .fchk should be created"

        print(f"✅ test_try_formchk_missing_chk passed")


if __name__ == "__main__":
    test_try_formchk_success()
    test_try_formchk_returncode_failure()
    test_try_formchk_timeout()
    test_try_formchk_missing_chk()
    print("\n✅ All unit tests passed")
