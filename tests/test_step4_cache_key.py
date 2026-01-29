import sys
from pathlib import Path as PathLib, Path

project_root = PathLib(__file__).parent.parent
sys.path.insert(0, str(project_root))

from rph_core.steps.step4_features.schema import generate_cache_key


def test_cache_key_deterministic_same_inputs():
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.xyz"
        test_file.write_text("C 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

        input_files = {"structure": test_file}
        params = {"method": "b3lyp", "basis": "6-31g*"}

        key1 = generate_cache_key(plugin_name="thermo", input_files=input_files, params=params)
        key2 = generate_cache_key(plugin_name="thermo", input_files=input_files, params=params)

        assert key1 == key2, "Same inputs should produce same cache key"


def test_cache_key_different_plugin_name():
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.xyz"
        test_file.write_text("C 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

        input_files = {"structure": test_file}
        params = {"method": "b3lyp"}

        key1 = generate_cache_key(plugin_name="thermo", input_files=input_files, params=params)
        key2 = generate_cache_key(plugin_name="geometry", input_files=input_files, params=params)

        assert key1 != key2, "Different plugin_name should produce different cache key"


def test_cache_key_different_file_content():
    import tempfile
    import time

    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.xyz"

        test_file.write_text("C 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        input_files = {"structure": test_file}
        params = {"method": "b3lyp"}

        key1 = generate_cache_key(plugin_name="thermo", input_files=input_files, params=params)

        time.sleep(0.01)

        test_file.write_text("C 1.0 0.0 0.0\nH 1.0 0.0 1.0\n")

        key2 = generate_cache_key(plugin_name="thermo", input_files=input_files, params=params)

        assert key1 != key2, "Different file content should produce different cache key"


def test_cache_key_different_params():
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.xyz"
        test_file.write_text("C 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

        input_files = {"structure": test_file}
        params1 = {"method": "b3lyp", "basis": "6-31g*"}
        params2 = {"method": "b3lyp", "basis": "6-311g**"}

        key1 = generate_cache_key(plugin_name="thermo", input_files=input_files, params=params1)
        key2 = generate_cache_key(plugin_name="thermo", input_files=input_files, params=params2)

        assert key1 != key2, "Different params should produce different cache key"


def test_cache_key_multiple_files():
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        file1 = Path(tmpdir) / "file1.xyz"
        file2 = Path(tmpdir) / "file2.xyz"
        file1.write_text("C 0.0 0.0 0.0\n")
        file2.write_text("H 0.0 0.0 1.0\n")

        input_files = {"ts": file1, "reactant": file2}
        params = {"method": "b3lyp"}

        key1 = generate_cache_key(plugin_name="thermo", input_files=input_files, params=params)
        key2 = generate_cache_key(plugin_name="thermo", input_files=input_files, params=params)

        assert key1 == key2, "Same multiple files should produce same cache key"


def test_cache_key_none_params():
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.xyz"
        test_file.write_text("C 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

        input_files = {"structure": test_file}

        key1 = generate_cache_key(plugin_name="thermo", input_files=input_files, params=None)
        key2 = generate_cache_key(plugin_name="thermo", input_files=input_files, params=None)

        assert key1 == key2, "None params should be handled consistently"


def test_cache_key_missing_file():
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "existing.xyz"
        missing_file = Path(tmpdir) / "missing.xyz"
        test_file.write_text("C 0.0 0.0 0.0\n")

        input_files = {"existing": test_file, "missing": missing_file}
        params = {"method": "b3lyp"}

        key = generate_cache_key(plugin_name="thermo", input_files=input_files, params=params)

        assert len(key) == 16, "Cache key should be 16 characters (sha1[:16])"


def test_cache_key_key_format():
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.xyz"
        test_file.write_text("C 0.0 0.0 0.0\n")

        input_files = {"structure": test_file}
        params = {"method": "b3lyp"}

        key = generate_cache_key(plugin_name="thermo", input_files=input_files, params=params)

        assert len(key) == 16, "Cache key should be 16 characters"
        assert all(c in "0123456789abcdef" for c in key), "Cache key should be hexadecimal"


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
