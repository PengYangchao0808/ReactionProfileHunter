"""
M4-C: Mechanism Context Resolver Tests
======================================

Tests for resolve_mechanism_context() - unified S1/S2/S3 asset resolution
with directory alias support and recursive search.

Author: QC Descriptors Team
Date: 2026-01-21
"""

import pytest
from pathlib import Path
import json
import tempfile
import shutil
from typing import Dict, Optional

from rph_core.steps.step4_features.mech_packager import (
    resolve_mechanism_context,
    S3_DIR_ALIASES,
    S2_DIR_ALIASES,
    S1_DIR_ALIASES,
    MechanismContext,
)


class TestResolveMechanismContext:
    """Test M4-C: Unified mechanism context resolver."""

    def test_resolve_mechanism_context_basic(self, tmp_path):
        """Basic resolution with standard directory structure."""
        work_dir = tmp_path / "work"
        s1_dir = work_dir / "S1_ConfGeneration"
        s2_dir = work_dir / "S2_Retro"
        s3_dir = work_dir / "S3_TS"

        s1_dir.mkdir(parents=True, exist_ok=True)
        s2_dir.mkdir(parents=True, exist_ok=True)
        s3_dir.mkdir(parents=True, exist_ok=True)

        s1_product = s1_dir / "product_min.xyz"
        s1_product.write_text("S1 product")

        s2_ts_guess = s2_dir / "ts_guess.xyz"
        s2_ts_guess.write_text("S2 ts_guess")

        s2_reactant = s2_dir / "reactant_complex.xyz"
        s2_reactant.write_text("S2 reactant complex")

        s3_ts_final = s3_dir / "ts_final.xyz"
        s3_ts_final.write_text("S3 ts_final")

        s3_reactant = s3_dir / "reactant_sp.xyz"
        s3_reactant.write_text("S3 reactant sp")

        context = resolve_mechanism_context(
            s4_dir=tmp_path / "S4_Data",
            config={},
            max_recursion_depth=3
        )

        assert context.s1_dir == s1_dir
        assert context.s2_dir == s2_dir
        assert context.s3_dir == s3_dir
        assert context.s1_product == s1_product
        assert context.s2_ts_guess == s2_ts_guess
        assert context.s2_reactant_complex == s2_reactant
        assert context.s3_ts_final == s3_ts_final
        assert context.s3_reactant_sp == s3_reactant
        assert context.s1_precursor_source == "none"

    @pytest.mark.parametrize("s3_alias", S3_DIR_ALIASES)
    def test_resolve_mechanism_context_s3_aliases(self, tmp_path, s3_alias):
        """Test S3 directory alias resolution."""
        work_dir = tmp_path / "work"
        s3_dir = work_dir / s3_alias
        s3_dir.mkdir(parents=True)
        s3_ts_final = s3_dir / "ts_final.xyz"
        s3_ts_final.write_text("S3 ts final")

        context = resolve_mechanism_context(
            s4_dir=tmp_path / "S4_Data",
            config={},
            max_recursion_depth=3
        )

        assert context.s3_dir == s3_dir
        assert context.s3_ts_final == s3_ts_final
        assert context.s3_reactant_sp is None

    def test_resolve_mechanism_context_missing_s3(self, tmp_path):
        """Test graceful handling when S3 directory is missing."""
        work_dir = tmp_path / "work"
        s1_dir = work_dir / "S1_ConfGeneration"
        s2_dir = work_dir / "S2_Retro"
        s1_dir.mkdir(parents=True)
        s2_dir.mkdir(parents=True)

        context = resolve_mechanism_context(
            s4_dir=tmp_path / "S4_Data",
            config={},
            max_recursion_depth=3
        )

        assert context.s1_dir == s1_dir
        assert context.s2_dir == s2_dir
        assert context.s3_dir is None
        assert context.s3_ts_final is None
        assert context.s3_reactant_sp is None

    def test_resolve_mechanism_context_missing_s1(self, tmp_path):
        """Test graceful handling when S1 directory is missing."""
        work_dir = tmp_path / "work"
        s3_dir = work_dir / "S3_TS"
        s2_dir = work_dir / "S2_Retro"
        s3_dir.mkdir(parents=True)
        s2_dir.mkdir(parents=True)

        context = resolve_mechanism_context(
            s4_dir=tmp_path / "S4_Data",
            config={},
            max_recursion_depth=3
        )

        assert context.s1_dir is None
        assert context.s2_dir == s2_dir
        assert context.s3_dir == s3_dir
        assert context.s1_product is None
        assert context.s1_precursor is None

    def test_resolve_mechanism_context_recursion_limit(self, tmp_path):
        """Test recursion depth limiting for S3 directory search."""
        work_dir = tmp_path / "work"
        s3_dir = work_dir / "S3_TS"
        s3_dir.mkdir(parents=True)

        nested = s3_dir / "nested" / "deeper"
        nested.mkdir(parents=True)
        (nested / "S3_TS").mkdir(parents=True)
        deeper_dir = nested / "deeper" / "S3_TS"
        deeper_dir.mkdir(parents=True, exist_ok=True)
        (deeper_dir / "ts_final.xyz").write_text("Nested S3")

        context = resolve_mechanism_context(
            s4_dir=tmp_path / "S4_Data",
            config={},
            max_recursion_depth=1
        )

        assert context.s3_dir is None
        assert context.s1_dir is None
        assert context.s2_dir is None

    def test_resolve_mechanism_context_config_priority(self, tmp_path):
        """Test precursor source priority from config."""
        work_dir = tmp_path / "work"
        s1_dir = work_dir / "S1_ConfGeneration"
        s2_dir = work_dir / "S2_Retro"
        s3_dir = work_dir / "S3_TS"
        s1_dir.mkdir(parents=True, exist_ok=True)
        s2_dir.mkdir(parents=True, exist_ok=True)
        s3_dir.mkdir(parents=True, exist_ok=True)
        s1_precursor = s1_dir / "neutral_precursor.xyz"
        s1_precursor.write_text("S1 neutral precursor")

        s2_reactant = s2_dir / "neutral_precursor.xyz"
        s2_reactant.write_text("S2 neutral precursor")

        context = resolve_mechanism_context(
            s4_dir=tmp_path / "S4_Data",
            config={"precursor_source_priority": ["S2_neutral_precursor"]},
            max_recursion_depth=3
        )

        assert context.s1_precursor == s2_reactant
        assert context.s1_precursor_source == "S2_neutral_precursor"

    def test_resolve_mechanism_context_work_dir_not_found(self, tmp_path):
        """Test RuntimeError when work_dir cannot be determined."""
        context = resolve_mechanism_context(
            s4_dir=tmp_path / "S4_Data_orphaned",
            config={},
            max_recursion_depth=3
        )

        assert context.s1_dir is None
        assert context.s2_dir is None
        assert context.s3_dir is None


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
