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
    """Test that try_formchk creates .fchk when .chk exists and formchk succeeds"""
    with tempfile.TemporaryDirectory(prefix="test_formchk_") as tmpdir:
        # Create test .chk file
        chk_file = Path(tmpdir) / "test.chk"
        chk_file.write_text("test checkpoint content", encoding='utf-8')

        # Patch rph_core.utils.qc_interface.subprocess.run
        # Import the module first so we can patch it correctly
        rph_core_path = Path(__file__).parent.parent / "rph_core"
        if str(rph_core_path) not in sys.path:
            sys.path.insert(0, str(rph_core_path))

        # Patch subprocess at the module level
        with patch('rph_core.utils.qc_interface.subprocess.run', autospec=True) as mock_run:
            # Simulate successful formchk: subprocess.run returns 0 and creates .fchk file
            mock_run.return_value = MagicMock(returncode=0)

            # Also patch Path.write to create the .fchk file
            original_write = Path.write_text

            def create_fchk_side_effect(self, content, *args, **kwargs):
                # If writing .fchk, create a dummy file
                if 'fchk' in str(self):
                    Path.write_text(self, content, *args, **kwargs)
                    # Create empty .fchk file
                    fchk_path = Path(str(self).replace('.chk', '.fchk'))
                    fchk_path.touch()

            with patch.object(Path, 'write_text', side_effect=create_fchk_side_effect):
                # Now call try_formchk
                from rph_core.utils.qc_interface import try_formchk
                result = try_formchk(chk_file)

        # Verify: .fchk file should be created and returned
        assert result is not None, "Should return .fchk path on success"
        assert result.exists(), "Should create .fchk file"
        assert str(result).endswith('.fchk'), "Should return .fchk path"

        print("✅ test_try_formchk_success passed")


# Test 2: formchk fails (returncode=1) should return None
def test_try_formchk_returncode_failure():
    """Test that try_formchk returns None when formchk fails (returncode=1)"""
    with tempfile.TemporaryDirectory(prefix="test_formchk_") as tmpdir:
        # Create test .chk file
        chk_file = Path(tmpdir) / "test.chk"
        chk_file.write_text("test checkpoint content", encoding='utf-8')

        # Import and patch
        rph_core_path = Path(__file__).parent.parent / "rph_core"
        if str(rph_core_path) not in sys.path:
            sys.path.insert(0, str(rph_core_path))

        with patch('rph_core.utils.qc_interface.subprocess.run', autospec=True) as mock_run:
            # Simulate formchk failure: subprocess.run returns 1
            mock_run.return_value = MagicMock(returncode=1, stderr="formchk Error")

            from rph_core.utils.qc_interface import try_formchk
            result = try_formchk(chk_file)

        # Verify: should return None and .fchk should not be created
        assert result is None, "Should return None on failure"
        assert not (Path(tmpdir) / "test.fchk").exists(), "fchk should not be created on failure"

        print("✅ test_try_formchk_returncode_failure passed")


# Test 3: formchk called on non-existent .chk should return None and log warning
def test_try_formchk_missing_chk():
    """Test that try_formchk returns None when .chk doesn't exist and logs warning"""
    with tempfile.TemporaryDirectory(prefix="test_formchk_") as tmpdir:
        # Don't create .chk file
        nonexist_file = Path(tmpdir) / "nonexistent.chk"

        # Import
        rph_core_path = Path(__file__).parent.parent / "rph_core"
        if str(rph_core_path) not in sys.path:
            sys.path.insert(0, str(rph_core_path))

        # Patch subprocess to avoid actual formchk call
        with patch('rph_core.utils.qc_interface.subprocess.run', autospec=True) as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            from rph_core.utils.qc_interface import try_formchk
            result = try_formchk(nonexist_file)

        # Verify: should return None and log warning
        assert result is None, "Should return None when .chk doesn't exist"
        assert not list(Path(tmpdir).glob("*.fchk")), "No .fchk should be created"

        print("✅ test_try_formchk_missing_chk passed")


if __name__ == "__main__":
    test_try_formchk_success()
    test_try_formchk_returncode_failure()
    test_try_formchk_missing_chk()
    print("\n✅ All unit tests passed")
