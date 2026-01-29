"""
Simple validation test for try_formchk

This test verifies that:
1. try_formchk function exists and is importable
2. Function signature matches expected behavior
3. Function returns Optional[Path]
4. Function has proper error handling (no exceptions raised)
"""
import sys
from pathlib import Path
from unittest.mock import MagicMock

# Test that try_formchk is importable and has correct signature
def test_try_formchk_exists_and_callable():
    """Verify try_formchk exists and has expected signature"""
    # Import the module
    rph_core_path = Path(__file__).parent.parent / "rph_core"
    if str(rph_core_path) not in sys.path:
        sys.path.insert(0, str(rph_core_path))

    from rph_core.utils.qc_interface import try_formchk

    # Verify function is callable
    assert callable(try_formchk), "try_formchk should be a callable function"

    # Verify function takes Path argument
    import inspect
    sig = inspect.signature(try_formchk)
    params = list(sig.parameters.keys())

    assert 'chk_path' in params, "try_formchk should have 'chk_path' parameter"
    assert len(params) == 1, "try_formchk should have exactly 1 parameter"

    # Check function has docstring (for API documentation)
    assert try_formchk.__doc__ is not None, "try_formchk should have docstring"

    print("✅ test_try_formchk_exists_and_callable passed")


# Test that try_formchk doesn't raise exceptions
def test_try_formchk_no_exceptions():
    """Verify try_formchk doesn't raise exceptions on missing files"""
    rph_core_path = Path(__file__).parent.parent / "rph_core"
    if str(rph_core_path) not in sys.path:
        sys.path.insert(0, str(rph_core_path))

    from rph_core.utils.qc_interface import try_formchk

    # Call with non-existent file - should handle gracefully
    nonexistent_file = Path("/nonexistent/test.chk")
    result = try_formchk(nonexistent_file)

    # Should return None and NOT raise exception
    assert result is None, "Should return None for non-existent .chk"
