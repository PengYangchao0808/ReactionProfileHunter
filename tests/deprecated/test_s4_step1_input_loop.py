"""
TDD Tests for Step1 Input Loop Closure (Task 7)
================================================

Tests to verify that:
1. .sum files can be found when .json doesn't exist
2. shermo_summary.json can be derived from .sum files
3. Derived artifacts have proper provenance

Author: RPH Team
Date: 2026-02-03
"""

import pytest
import tempfile
import json
from pathlib import Path

from rph_core.utils.constants import HARTREE_TO_KCAL


def _write_mock_sum_file(path: Path, g_sum: float, h_sum: float, u_sum: float) -> None:
    """Write mock Shermo .sum file."""
    content = f"""Sum of electronic energy and thermal correction to U = {u_sum:.6f}
Sum of electronic energy and thermal correction to H = {h_sum:.6f}
Sum of electronic energy and thermal correction to G = {g_sum:.6f}
Gibbs free energy at specified concentration = {g_sum:.6f}
Total S = 50.0
"""
    path.write_text(content)


class TestShermoSummaryDerivation:
    """Tests for deriving shermo_summary.json from .sum files."""

    def test_derive_shermo_summary_from_sum(self):
        """Should derive valid JSON from .sum file."""
        from rph_core.utils.shermo_runner import derive_shermo_summary_from_sum

        with tempfile.TemporaryDirectory(prefix="test_derive_") as tmp:
            tmpdir = Path(tmp)

            # Create mock .sum file
            sum_file = tmpdir / "precursor_Shermo.sum"
            _write_mock_sum_file(sum_file, g_sum=-500.0, h_sum=-490.0, u_sum=-495.0)

            # Derive JSON
            output_json = tmpdir / "shermo_summary.json"
            result = derive_shermo_summary_from_sum(
                sum_file,
                output_json,
                molecule_type="precursor"
            )

            # Verify JSON was created
            assert output_json.exists(), "JSON file should be created"

            # Verify content
            with open(output_json, 'r') as f:
                data = json.load(f)

            assert data["unit"] == "kcal/mol"
            expected_g = -500.0 * HARTREE_TO_KCAL
            expected_h = -490.0 * HARTREE_TO_KCAL
            expected_u = -495.0 * HARTREE_TO_KCAL
            assert abs(data["g_precursor"] - expected_g) < 0.001
            assert abs(data["g_sum"] - expected_g) < 0.001
            assert abs(data["h_sum"] - expected_h) < 0.001
            assert abs(data["u_sum"] - expected_u) < 0.001
            assert data["derived_artifacts"] == True
            assert "derived_from_sum" in data

    def test_find_shermo_sum_files(self):
        """Should find .sum files in S1 directory."""
        from rph_core.utils.shermo_runner import find_shermo_sum_files

        with tempfile.TemporaryDirectory(prefix="test_find_") as tmp:
            tmpdir = Path(tmp)

            # Create mock .sum files
            (tmpdir / "product" / "dft").mkdir(parents=True, exist_ok=True)
            sum_file = tmpdir / "product" / "dft" / "precursor_Shermo.sum"
            _write_mock_sum_file(sum_file, g_sum=-500.0, h_sum=-490.0, u_sum=-495.0)

            hoac_sum = tmpdir / "hoac_Shermo.sum"
            _write_mock_sum_file(hoac_sum, g_sum=-200.0, h_sum=-195.0, u_sum=-198.0)

            # Find .sum files
            results = find_shermo_sum_files(tmpdir)

            # Verify results
            assert results["precursor"] is not None
            assert results["hoac"] is not None

    def test_orchestrator_resolves_s1_artifacts(self):
        """Orchestrator should resolve S1 artifacts including derived shermo_summary.json."""
        from rph_core.orchestrator import ReactionProfileHunter

        with tempfile.TemporaryDirectory(prefix="test_orch_") as tmp:
            tmpdir = Path(tmp)
            work_dir = tmpdir / "work"
            work_dir.mkdir(parents=True, exist_ok=True)
            s1_dir = work_dir / "S1_ConfGeneration"
            s1_dir.mkdir(parents=True, exist_ok=True)

            # Create mock .sum file (no JSON)
            (s1_dir / "product" / "dft").mkdir(parents=True, exist_ok=True)
            sum_file = s1_dir / "product" / "dft" / "precursor_Shermo.sum"
            _write_mock_sum_file(sum_file, g_sum=-500.0, h_sum=-490.0, u_sum=-495.0)

            # Create minimal config file
            config_content = """
global:
  output_dir: {output_dir}
executables: {{}}
""".format(output_dir=tmpdir)

            config_file = tmpdir / "config.yaml"
            config_file.write_text(config_content)

            # Create orchestrator with config file
            hunter = ReactionProfileHunter(config_path=config_file)

            # Resolve S1 artifacts
            artifacts = hunter._resolve_s1_artifacts(work_dir)

            # Should find S1 directory
            assert artifacts["s1_dir"] is not None
            assert artifacts["s1_dir"].exists()

            # Should derive shermo_summary.json
            assert artifacts["s1_shermo_summary_file"] is not None
            assert artifacts["s1_shermo_summary_file"].exists()

            # Verify derived JSON content (converted to kcal/mol)
            with open(artifacts["s1_shermo_summary_file"], 'r') as f:
                data = json.load(f)
            expected_g = -500.0 * HARTREE_TO_KCAL
            assert abs(data["g_precursor"] - expected_g) < 0.001


class TestHOAcThermoDerivation:
    """Tests for HOAc thermochemistry derivation."""

    def test_derive_hoac_thermo_from_sum(self):
        """Should derive HOAc thermo.json from .sum file."""
        from rph_core.utils.shermo_runner import derive_hoac_thermo_from_sum

        with tempfile.TemporaryDirectory(prefix="test_hoac_") as tmp:
            tmpdir = Path(tmp)

            # Create mock HOAc .sum file
            sum_file = tmpdir / "HOAc_Shermo.sum"
            _write_mock_sum_file(sum_file, g_sum=-200.0, h_sum=-195.0, u_sum=-198.0)

            # Derive thermo.json
            output_json = tmpdir / "thermo.json"
            result = derive_hoac_thermo_from_sum(sum_file, output_json)

            # Verify JSON was created
            assert output_json.exists(), "HOAc thermo.json should be created"

            # Verify content has both 'g' and 'G' keys for compatibility
            with open(output_json, 'r') as f:
                data = json.load(f)

            assert data["unit"] == "kcal/mol"
            # Should be converted from Hartree to kcal/mol
            expected_g = -200.0 * HARTREE_TO_KCAL
            assert abs(data["g"] - expected_g) < 0.001
            assert abs(data["G"] - expected_g) < 0.001
            assert data["derived_artifacts"] == True
            assert "derived_from_sum" in data
