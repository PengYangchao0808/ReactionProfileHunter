"""
M4-A: Template Structure Tests (NBO-only, V5.4)
================================================

Tests to ensure Gaussian input templates follow correct formatting.

V5.4 Update: Removed NMR/Hirshfeld template references - only NBO and core templates remain.
"""

import pytest
from pathlib import Path
from typing import List

# Template files that use {{charge}} {{mult}} {{coords_block}} syntax
RENDERABLE_TEMPLATES = [
    "gaussian_ts.gjf",
    "gaussian_nbo.gjf",
    "gaussian_freq.gjf",
]

# Templates that don't use these placeholders (hardcoded coordinates)
NON_RENDERABLE_TEMPLATES = [
    "gaussian_qst2.gjf",  # Uses [reactant coordinates] and [product coordinates]
    "gaussian_irc.gjf",   # Uses [coordinates]
]


class TestTemplateChargeMultLineIsolated:
    """Test that {{charge}} {{mult}} placeholder is on separate line."""

    def test_template_charge_mult_line_isolated(self):
        template_dir = Path(__file__).parent.parent / "config" / "templates"

        for template_name in RENDERABLE_TEMPLATES:
            template_path = template_dir / template_name
            assert template_path.exists(), f"Template not found: {template_name}"

            content = template_path.read_text()
            lines = content.split("\n")

            charge_mult_line_idx = None
            for idx, line in enumerate(lines):
                if "{{charge}}" in line and "{{mult}}" in line:
                    charge_mult_line_idx = idx
                    break

            assert charge_mult_line_idx is not None, (
                f"Template {template_name} does not contain {{charge}} {{mult}} placeholder"
            )

            charge_mult_line = lines[charge_mult_line_idx].strip()
            assert charge_mult_line == "{{charge}} {{mult}}", (
                f"Template {template_name}: line with {{charge}} {{mult}} must contain only these placeholders, "
                f"but found: '{charge_mult_line}'"
            )


class TestTemplateCoordsBlockNotInlineWithChargeMult:
    """
    M4-A-4: Verify that {{coords_block}} starts on NEXT line after {{charge}} {{mult}}.

    This is required by Gaussian input format - coordinates cannot be on same line
    as charge/mult.
    """

    def test_template_coords_block_not_inline_with_charge_mult(self):
        template_dir = Path(__file__).parent.parent / "config" / "templates"

        for template_name in RENDERABLE_TEMPLATES:
            template_path = template_dir / template_name
            assert template_path.exists(), f"Template not found: {template_name}"

            content = template_path.read_text()
            lines = content.split("\n")

            # Find line containing {{charge}} {{mult}}
            charge_mult_line_idx = None
            for idx, line in enumerate(lines):
                if "{{charge}}" in line and "{{mult}}" in line:
                    charge_mult_line_idx = idx
                    break

            # Check that {{coords_block}} is on NEXT line
            next_line_idx = charge_mult_line_idx + 1
            if next_line_idx < len(lines):
                next_line = lines[next_line_idx]
                assert "{{coords_block}}" in next_line, (
                    f"Template {template_name}: {{coords_block}} must be on line immediately after "
                    f"{{charge}} {{mult}}, but found '{next_line}' on next line"
                )
            else:
                pytest.fail(
                    f"Template {template_name}: {{coords_block}} not found after {{charge}} {{mult}} line"
                )


class TestTemplateStructureParameterized:
    """M4-A-5: Parameterized tests for all template files."""

    @pytest.mark.parametrize("template_name", RENDERABLE_TEMPLATES)
    def test_template_has_required_placeholders(self, template_name):
        """Verify renderable templates have all required placeholders."""
        template_dir = Path(__file__).parent.parent / "config" / "templates"
        template_path = template_dir / template_name
        assert template_path.exists()

        content = template_path.read_text()
        assert "{{charge}}" in content, f"Template {template_name} missing {{charge}} placeholder"
        assert "{{mult}}" in content, f"Template {template_name} missing {{mult}} placeholder"
        assert "{{coords_block}}" in content, f"Template {template_name} missing {{coords_block}} placeholder"

    @pytest.mark.parametrize("template_name", RENDERABLE_TEMPLATES)
    def test_template_structure_complies_with_gaussian_format(self, template_name):
        """
        Comprehensive test: verify full structure complies with Gaussian input format.
        1. {{charge}} {{mult}} on separate line
        2. {{coords_block}} on next line (not inline)
        3. No other content on those lines (except blank lines)
        """
        template_dir = Path(__file__).parent.parent / "config" / "templates"
        template_path = template_dir / template_name
        assert template_path.exists()

        content = template_path.read_text()
        lines = content.split("\n")

        # Find critical lines
        charge_mult_line_idx = None
        coords_block_line_idx = None

        for idx, line in enumerate(lines):
            if "{{charge}}" in line and "{{mult}}" in line:
                charge_mult_line_idx = idx
            elif "{{coords_block}}" in line:
                coords_block_line_idx = idx

        assert charge_mult_line_idx is not None, f"Template {template_name} missing {{charge}} {{mult}}"
        assert coords_block_line_idx is not None, f"Template {template_name} missing {{coords_block}}"

        # Validate charge/mult line
        charge_mult_line = lines[charge_mult_line_idx].strip()
        assert charge_mult_line == "{{charge}} {{mult}}", (
            f"Template {template_name}: line {charge_mult_line_idx + 1} should contain only "
            f"'{{charge}} {{mult}}', but found: '{charge_mult_line}'"
        )

        # Validate coords_block line
        assert coords_block_line_idx == charge_mult_line_idx + 1, (
            f"Template {template_name}: line {coords_block_line_idx + 1} should contain only "
            f"'{{coords_block}}', but found: '{lines[coords_block_line_idx]}'"
        )

    @pytest.mark.parametrize("template_name", NON_RENDERABLE_TEMPLATES)
    def test_non_renderable_templates_unchanged(self, template_name):
        """
        Verify non-renderable templates (with hardcoded coordinates) don't use placeholders.
        These templates use special coordinate blocks for QST2, IRC, etc.
        """
        template_dir = Path(__file__).parent.parent / "config" / "templates"
        template_path = template_dir / template_name
        assert template_path.exists()

        content = template_path.read_text()

        # These should NOT have placeholder syntax
        if template_name == "gaussian_qst2.gjf":
            assert "[reactant coordinates]" in content
            assert "[product coordinates]" in content
        elif template_name == "gaussian_irc.gjf":
            assert "[coordinates]" in content

        # Verify no placeholder syntax in non-renderable templates
        assert "{{charge}}" not in content
        assert "{{mult}}" not in content
        assert "{{coords_block}}" not in content


class TestTemplateRenderingContract:
    """Test that templates can be rendered without breaking Gaussian format."""

    @pytest.mark.parametrize("template_name", RENDERABLE_TEMPLATES)
    def test_template_rendering_preserves_format(self, template_name):
        """
        Simulate template rendering and verify Gaussian format is preserved.
        This catches potential issues like newline stripping or incorrect formatting.
        """
        template_dir = Path(__file__).parent.parent / "config" / "templates"
        template_path = template_dir / template_name
        template_content = template_path.read_text()

        # Simulate rendering with sample values
        rendered = template_content.replace("{{charge}}", "0").replace("{{mult}}", "1")
        coords_sample = "C    0.000000    0.000000    0.000000\nH    1.089000    0.000000    0.000000"
        rendered = rendered.replace("{{coords_block}}", coords_sample)

        # Basic validation: rendered content should have charge/mult line followed by coords
        lines = rendered.split("\n")

        # Find charge/mult line (now "0 1")
        charge_mult_idx = None
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if stripped == "0 1":
                charge_mult_idx = idx
                break

        assert charge_mult_idx is not None, (
            f"Template {template_name}: could not find rendered charge/mult line '0 1'"
        )

        # Verify next line has coordinates (starts with atom symbol)
        if charge_mult_idx + 1 < len(lines):
            next_line = lines[charge_mult_idx + 1]
            # Should start with an atomic symbol or have atom coordinates
            assert next_line.strip().startswith("C") or next_line.strip().startswith("H"), (
                f"Template {template_name}: line after charge/mult should contain coordinates, "
                f"but found: '{next_line}'"
            )
        else:
            pytest.fail(
                f"Template {template_name}: no line after charge/mult line '0 1'"
            )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
