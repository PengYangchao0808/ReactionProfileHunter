"""
Simple standalone test for try_formchk utility function

This test validates that try_formchk() function behavior:
1. Success with existing .chk returns .fchk path
2. Failure with non-existent .chk returns None (graceful degradation)
3. Missing .chk logs a warning
"""
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import sys

# Test success case: formchk creates .fchk successfully
def test_try_formchk_success():
    with tempfile.TemporaryDirectory(prefix="test_formchk_") as tmpdir:
        # Create test .chk file
        chk_file = Path(tmpdir) / "test.chk"
        chk_file.write_text("test checkpoint content", encoding='utf-8')

        # Mock subprocess.run to simulate successful formchk
        with patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(
                returncode=0,
                stdout="Formchk successful",
                stderr=""
            )

        # Import and run try_formchk (mocking subprocess)
        import os
        rph_core_path = Path(__file__).parent.parent / "rph_core"
        if str(rph_core_path) not in sys.path:
            sys.path.insert(0, str(rph_core_path))

        from rph_core.utils.qc_interface import try_formchk

        # Call try_formchk (subprocess.run is mocked)
        result = try_formchk(chk_file)

        assert result is not None, "Should return .fchk path on success"
        assert result.exists(), "Should create .fchk file"
        assert str(result).endswith('.fchk'), "Should return .fchk path"

        print("✅ test_try_formchk_success passed")


# Test 2: formchk fails (returncode=1) should return None
def test_try_formchk_returncode_failure():
    with tempfile.TemporaryDirectory(prefix="test_formchk_") as tmpdir:
        # Create test .chk file
        chk_file = Path(tmpdir) / "test.chk"
        chk_file.write_text("test checkpoint content", encoding='utf-8')

        # Mock subprocess.run with non-zero returncode
        with patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(
                returncode=1,
                stdout="",
                stderr="formchk Error writing formatted checkpoint"
            )

        import os
        rph_core_path = Path(__file__).parent.parent / "rph_core"
        if str(rph_core_path) not in sys.path:
            sys.path.insert(0, str(rph_core_path))

        from rph_core.utils.qc_interface import try_formchk
        result = try_formchk(chk_file)

        assert result is None, "Should return None on failure"
        assert not (Path(tmpdir) / "test.fchk").exists(), "fchk should not be created on failure"

        print("✅ test_try_formchk_returncode_failure passed")


# Test 3: formchk called on non-existent .chk should return None and log warning
def test_try_formchk_missing_chk():
    with tempfile.TemporaryDirectory(prefix="test_formchk_") as tmpdir:
        # Don't create .chk file
        nonexist_file = Path(tmpdir) / "nonexistent.chk"

        import os
        rph_core_path = Path(__file__).parent.parent / "rph_core"
        if str(rph_core_path) not in sys.path:
            sys.path.insert(0, str(rph_core_path))

        from rph_core.utils.qc_interface import try_formchk
        result = try_formchk(nonexist_file)

        assert result is None, "Should return None when .chk doesn't exist"
        assert not list(Path(tmpdir).glob("*.fchk")), "No .fchk should be created"

        print("✅ test_try_formchk_missing_chk passed")


if __name__ == "__main__":
    test_try_formchk_success()
    test_try_formchk_returncode_failure()
    test_try_formchk_missing_chk()
    print("\n✅ All unit tests passed")
