"""
Unit tests for try_formchk() function

Tests the formchk utility with monkeypatch to ensure:
- Success case: formchk completes successfully and returns .fchk path
- Failure cases: formchk fails (returncode, timeout, etc.) and returns None
- Missing chk: formchk is called on non-existent .chk file
"""
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def test_try_formchk_success(tmp_path):
    """Test successful formchk execution"""
    # Create a mock .chk file
    chk_file = tmp_path / "test.chk"
    chk_file.write_text("test checkpoint content")

    # Patch subprocess.run to simulate successful formchk
    with patch('subprocess.run') as mock_run:
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Formchk successful",
            stderr=""
        )
        mock_run.side_effect = None  # No file modifications needed

    # Import and run try_formchk
    from rph_core.utils.qc_interface import try_formchk

    result = try_formchk(chk_file)

    assert result is not None, "Should return fchk path on success"
    assert result == chk_file.with_suffix('.fchk'), f"Should return correct fchk path"
    assert result.exists(), "fchk file should exist"

    # Cleanup
    chk_file.unlink()


def test_try_formchk_returncode_failure(tmp_path):
    """Test formchk when subprocess returns non-zero return code"""
    chk_file = tmp_path / "test.chk"
    chk_file.write_text("test checkpoint content")

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


def test_try_formchk_timeout(tmp_path):
    """Test formchk when subprocess times out"""
    chk_file = tmp_path / "test.chk"
    chk_file.write_text("test checkpoint content")

    with patch('subprocess.run') as mock_run:
        mock_run.side_effect = subprocess.TimeoutExpired("formchk timed out")

    from rph_core.utils.qc_interface import try_formchk
    result = try_formchk(chk_file)

    assert result is None, "Should return None on timeout"
    assert not chk_file.with_suffix('.fchk').exists(), "fchk should not be created on timeout"


def test_try_formchk_missing_chk(tmp_path):
    """Test formchk when .chk file does not exist"""
    chk_file = tmp_path / "nonexistent.chk"
    # Don't create the file

    with patch('subprocess.run') as mock_run:
        mock_run.return_value = MagicMock(returncode=0)

    from rph_core.utils.qc_interface import try_formchk
    result = try_formchk(chk_file)

    # Should return None (graceful degradation) and log warning
    assert result is None, "Should return None when chk missing"


def test_try_formchk_fchk_not_created(tmp_path):
    """Test formchk when subprocess succeeds but .fchk file is not created"""
    chk_file = tmp_path / "test.chk"
    chk_file.write_text("test checkpoint content")

    with patch('subprocess.run') as mock_run:
        # Simulate successful subprocess but fchk not created
        mock_run.return_value = MagicMock(returncode=0)

    from rph_core.utils.qc_interface import try_formchk
    result = try_formchk(chk_file)

    assert result is None, "Should return None when fchk not created"


if __name__ == "__main__":
    import tempfile

    tmpdir = Path(tempfile.mkdtemp(prefix="test_formchk_"))
    try:
        test_try_formchk_success(tmpdir)
        test_try_formchk_returncode_failure(tmpdir)
        test_try_formchk_timeout(tmpdir)
        test_try_formchk_missing_chk(tmpdir)
        test_try_formchk_fchk_not_created(tmpdir)
        print("✅ All try_formchk tests passed")
    finally:
        # Cleanup
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)
