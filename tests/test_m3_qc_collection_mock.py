"""
M3-4: Mock QC Collection Tests (Simple Version)
===================================================================

Simplified tests for mock NBO collection functions.
Tests verify file existence and absolute paths without external QC calls.

Author: QC Descriptors Team
Date: 2026-01-21 (V5.4: NMR/Hirshfeld removed)
"""

import pytest
from pathlib import Path
from rph_core.utils.qc_interface import (
    NBO_WHITELIST,
    harvest_nbo_files,
)


class TestMockNBOCollection:
    """Test NBO file collection (mock only)."""
    
    def test_nbo_whitelist_exists(self):
        """NBO_WHITELIST constant should exist."""
        assert len(NBO_WHITELIST) > 0
        assert '.47' in NBO_WHITELIST
        assert '.nbo' in NBO_WHITELIST
        assert '.3' in NBO_WHITELIST
        assert '.31' in NBO_WHITELIST
        assert '.41' in NBO_WHITELIST
    
    def test_collect_nbo_files_returns_dict(self, tmp_path):
        """harvest_nbo_files should return dictionary of paths."""
        from rph_core.utils.qc_interface import harvest_nbo_files
        
        nbo_dir = tmp_path / "nbo_analysis"
        nbo_dir.mkdir(parents=True, exist_ok=True)
        
        jobname = "test_job"
        
        for ext in NBO_WHITELIST:
            (nbo_dir / f"{jobname}{ext}").write_text(f"NBO file {ext}")
        
        result = harvest_nbo_files(tmp_path, "test_job")
        
        assert isinstance(result, dict)
        assert len(result) == len(NBO_WHITELIST)
        for ext in NBO_WHITELIST:
            assert ext in result
            assert result[ext].name == f"test_job{ext}"
            assert result[ext].is_absolute()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
