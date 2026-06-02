"""Tests for S0 DR SMILES comparison image rendering.

Validates that render_dr_comparison_images produces the expected PNG outputs
when given a valid DRBranchPlan with original and flipped SMILES.
"""
import json
import pytest
from pathlib import Path

from pydantic import BaseModel


class _FakeBranch(BaseModel):
    branch_id: str = "BR_MAJOR"
    pathway_id: str = "primary"
    role: str = "major_or_reference"
    generation_policy: str = "reference_product"
    product_smiles: str | None = None
    flipped_map_numbers: list[int] = []
    fixed_stereocenters: list[int] = []
    notes: list[str] = []


class _FakePlan(BaseModel):
    version: str = "s0-dr-branch-plan-v1"
    reaction_id: str = "rx_test"
    status: str = "complete"
    reaction_type: str = "[4+3]"
    cyclo_mode: str = "C4_PLUS_3"
    topology: str = "INTRA_TYPE_I"
    reported_dr: dict = {}
    branches: list[_FakeBranch] = []


@pytest.fixture
def two_branch_plan() -> _FakePlan:
    return _FakePlan(
        branches=[
            _FakeBranch(
                branch_id="BR_MAJOR",
                role="major_or_reference",
                product_smiles="O=[C:1]1[CH2:2][CH2:3][C@@H:4]1O",
                fixed_stereocenters=[1, 2, 3],
            ),
            _FakeBranch(
                branch_id="BR_DR_001",
                role="diastereomer_candidate",
                product_smiles="O=[C:1]1[CH2:2][CH2:3][C@H:4]1O",
                flipped_map_numbers=[4],
                fixed_stereocenters=[1, 2, 3],
            ),
        ]
    )


@pytest.fixture
def single_branch_plan() -> _FakePlan:
    return _FakePlan(
        branches=[
            _FakeBranch(
                branch_id="BR_MAJOR",
                role="major_or_reference",
                product_smiles="C1CC1",
            )
        ]
    )


class TestRenderDRComparisonImages:
    def test_single_branch_produces_no_output(
        self, single_branch_plan, tmp_path
    ):
        from rph_core.steps.mechanism_classifier.dr_smiles_renderer import (
            render_dr_comparison_images,
        )

        output_dir = tmp_path / "viz"
        result = render_dr_comparison_images(single_branch_plan, output_dir)
        assert result is False
        assert not list(output_dir.glob("*.png"))

    def test_two_branches_produces_comparison_and_per_branch(
        self, two_branch_plan, tmp_path
    ):
        from rph_core.steps.mechanism_classifier.dr_smiles_renderer import (
            render_dr_comparison_images,
        )

        output_dir = tmp_path / "viz"
        result = render_dr_comparison_images(two_branch_plan, output_dir)
        assert result is True

        assert (output_dir / "smiles_comparison.png").exists()
        assert (output_dir / "smiles_BR_MAJOR.png").exists()
        assert (output_dir / "smiles_BR_DR_001.png").exists()

    def test_output_files_are_valid_png(self, two_branch_plan, tmp_path):
        from rph_core.steps.mechanism_classifier.dr_smiles_renderer import (
            render_dr_comparison_images,
        )

        output_dir = tmp_path / "viz"
        render_dr_comparison_images(two_branch_plan, output_dir)

        for name in ("smiles_comparison.png", "smiles_BR_MAJOR.png", "smiles_BR_DR_001.png"):
            data = (output_dir / name).read_bytes()
            assert data[:8] == b"\x89PNG\r\n\x1a\n", f"{name} is not a valid PNG"

    def test_custom_prefix(self, two_branch_plan, tmp_path):
        from rph_core.steps.mechanism_classifier.dr_smiles_renderer import (
            render_dr_comparison_images,
        )

        output_dir = tmp_path / "viz"
        render_dr_comparison_images(two_branch_plan, output_dir, prefix="dr_check")
        assert (output_dir / "dr_check_comparison.png").exists()
        assert (output_dir / "dr_check_BR_MAJOR.png").exists()

    def test_missing_smiles_is_skipped(self, tmp_path):
        from rph_core.steps.mechanism_classifier.dr_smiles_renderer import (
            render_dr_comparison_images,
        )

        plan = _FakePlan(
            branches=[
                _FakeBranch(
                    branch_id="BR_MAJOR",
                    product_smiles="C1CC1",
                ),
                _FakeBranch(
                    branch_id="BR_DR_001",
                    product_smiles=None,
                ),
            ]
        )
        output_dir = tmp_path / "viz"
        result = render_dr_comparison_images(plan, output_dir)
        assert result is True
        assert (output_dir / "smiles_BR_MAJOR.png").exists()
        assert not (output_dir / "smiles_BR_DR_001.png").exists()

    def test_zero_branches_no_crash(self, tmp_path):
        from rph_core.steps.mechanism_classifier.dr_smiles_renderer import (
            render_dr_comparison_images,
        )

        plan = _FakePlan(branches=[])
        result = render_dr_comparison_images(plan, tmp_path / "viz")
        assert result is False

    def test_rendering_non_blocking_on_bad_smiles(self, tmp_path):
        from rph_core.steps.mechanism_classifier.dr_smiles_renderer import (
            render_dr_comparison_images,
        )

        plan = _FakePlan(
            branches=[
                _FakeBranch(
                    branch_id="BR_MAJOR",
                    product_smiles="NOT_A_VALID_SMILES_ZZZZ",
                ),
                _FakeBranch(
                    branch_id="BR_DR_001",
                    product_smiles="C1CC1",
                ),
            ]
        )
        output_dir = tmp_path / "viz"
        result = render_dr_comparison_images(plan, output_dir)
        assert result is True
        assert (output_dir / "smiles_BR_DR_001.png").exists()
