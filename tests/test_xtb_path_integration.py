"""
Tests for xTB Path Search Integration
"""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
import tempfile
import json


class TestPathSearchResult:
    """Test PathSearchResult dataclass"""
    
    def test_path_search_result_creation(self):
        from rph_core.utils.data_types import PathSearchResult
        
        result = PathSearchResult(
            success=True,
            ts_guess_xyz=Path("ts_guess.xyz"),
            barrier_forward_kcal=12.5,
            barrier_backward_kcal=37.0,
            reaction_energy_kcal=-24.5,
        )
        
        assert result.success is True
        assert result.ts_guess_xyz == Path("ts_guess.xyz")
        assert result.barrier_forward_kcal == 12.5
        assert result.barrier_backward_kcal == 37.0
        assert result.reaction_energy_kcal == -24.5


class TestXTBRunnerPath:
    """Test XTBRunner.run_path method"""
    
    def test_run_path_method_exists(self):
        from rph_core.utils.xtb_runner import XTBRunner
        assert hasattr(XTBRunner, 'run_path')
    
    def test_parse_path_log_method_exists(self):
        from rph_core.utils.xtb_runner import XTBRunner
        assert hasattr(XTBRunner, '_parse_path_log')


class TestXTBInterfacePath:
    """Test XTBInterface.path method"""
    
    def test_path_method_exists(self):
        from rph_core.utils.qc_interface import XTBInterface
        assert hasattr(XTBInterface, 'path')


class TestArtifactResolver:
    """Test S3 Artifact Resolver"""
    
    def test_resolve_s3_inputs_function_exists(self):
        from rph_core.steps.step3_opt.artifact_resolver import resolve_s3_inputs
        assert callable(resolve_s3_inputs)
    
    def test_check_s2_artifacts_function_exists(self):
        from rph_core.steps.step3_opt.artifact_resolver import check_s2_artifacts
        assert callable(check_s2_artifacts)
    
    def test_check_s2_artifacts_returns_dict(self, tmp_path):
        from rph_core.steps.step3_opt.artifact_resolver import check_s2_artifacts
        
        result = check_s2_artifacts(tmp_path)
        
        assert isinstance(result, dict)
        assert "ts_guess_exists" in result
        assert "reactant_complex_exists" in result
        assert "dipolar_intermediate_exists" in result
        assert "scan_profile_exists" in result


class TestConfigPathSearch:
    """Test path_search config in defaults.yaml"""
    
    def test_path_search_config_exists(self):
        from rph_core.orchestrator import ReactionProfileHunter
        from pathlib import Path
        
        hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))
        step2_cfg = hunter.config.get("step2", {})
        path_search_cfg = step2_cfg.get("path_search", {})
        
        assert "enabled" in path_search_cfg
        assert path_search_cfg.get("enabled") is False
    
    def test_path_search_params_exist(self):
        from rph_core.orchestrator import ReactionProfileHunter
        from pathlib import Path
        
        hunter = ReactionProfileHunter(config_path=Path("config/defaults.yaml"))
        step2_cfg = hunter.config.get("step2", {})
        path_search_cfg = step2_cfg.get("path_search", {})
        
        assert "nrun" in path_search_cfg
        assert "npoint" in path_search_cfg
        assert "anopt" in path_search_cfg
        assert "kpush" in path_search_cfg
        assert "kpull" in path_search_cfg


class TestRetroScannerPathSearch:
    """Test RetroScanner.run_path_search method"""
    
    def test_run_path_search_method_exists(self):
        from rph_core.steps.step2_retro import RetroScanner
        assert hasattr(RetroScanner, 'run_path_search')


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
