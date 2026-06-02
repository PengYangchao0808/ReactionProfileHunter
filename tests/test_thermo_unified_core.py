"""
Tests for the unified thermochemistry module (rph_core.utils.thermo).

Covers: parse_shermo_sum, ThermoRecord, write/read_thermo_json roundtrip,
derive_thermo_json_from_sum, ensure_thermo_json_from_entry, legacy fallback.
"""

import json
import os
import subprocess
import tempfile
from pathlib import Path

import pytest

from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.utils.thermo import (
    ThermoRecord,
    ShermoOptions,
    parse_shermo_sum,
    write_thermo_json,
    read_thermo_json,
    derive_thermo_json_from_sum,
    ensure_thermo_json_from_entry,
    derive_shermo_summary,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
SAMPLE_SUM = REPO_ROOT / "tests" / "tmp_v2_2_test" / "da_reaction" / "S1_Anchor" / "ethylene" / "dft" / "conf_000_Shermo.sum"

# Minimal Shermo 2.6 .sum content for unit tests (no real QC needed)
SHERMO_26_SUM = """\
 Shermo: A general code for calculating molecular thermochemistry properties
 Version 2.6

 Temperature:     298.150 K
 Pressure:          1.000 atm

 Electronic energy:        -77.0749052 a.u.
                       ========== Total ==========
 Total S:      217.494 J/mol/K      51.982 cal/mol/K    -TS:   -15.499 kcal/mol
 Thermal correction to U:    163.728 kJ/mol     39.132 kcal/mol   0.062361 a.u.
 Thermal correction to H:    166.207 kJ/mol     39.724 kcal/mol   0.063305 a.u.
 Thermal correction to G:    101.361 kJ/mol     24.226 kcal/mol   0.038606 a.u.
 Sum of electronic energy and thermal correction to U:         -77.0125445 a.u.
 Sum of electronic energy and thermal correction to H:         -77.0116003 a.u.
 Sum of electronic energy and thermal correction to G:         -77.0362988 a.u.
 Gibbs free energy at specified concentration 1.000 mol/L:   -77.0375000 a.u.
"""

# Legacy format (old Shermo versions or custom outputs)
LEGACY_SUM = """\
g(total)  -77.0362988
h(total)  -77.0116003
u(total)  -77.0125445
s(total)  51.982
"""


@pytest.fixture
def tmp_dir(tmp_path):
    """Provide a temporary directory for test outputs."""
    return tmp_path


# ---------------------------------------------------------------------------
# 1. parse_shermo_sum
# ---------------------------------------------------------------------------

class TestParseShermoSum:
    def test_parse_real_sum_file(self):
        """Parse actual ethylene Shermo .sum from test fixtures."""
        if not SAMPLE_SUM.exists():
            pytest.skip(f"Sample .sum not found: {SAMPLE_SUM}")

        record = parse_shermo_sum(SAMPLE_SUM)
        assert isinstance(record, ThermoRecord)
        assert record.g_sum_hartree == pytest.approx(-77.0362988, abs=1e-5)
        assert record.h_sum_hartree == pytest.approx(-77.0116003, abs=1e-5)
        assert record.u_sum_hartree == pytest.approx(-77.0125445, abs=1e-5)
        # Ethylene should have positive entropy in cal/mol/K
        assert record.s_cal_mol_k is not None
        assert record.s_cal_mol_k > 0
        # g_kcal should be reasonable (negative, ~-48300 kcal/mol for ethylene)
        assert record.g_kcal < 0

    def test_parse_shermo26_format(self, tmp_dir):
        """Parse Shermo 2.6 stdout format from synthetic content."""
        sum_file = tmp_dir / "test_26.sum"
        sum_file.write_text(SHERMO_26_SUM)

        record = parse_shermo_sum(sum_file)
        assert record.g_sum_hartree == pytest.approx(-77.0362988, abs=1e-5)
        assert record.h_sum_hartree == pytest.approx(-77.0116003, abs=1e-5)
        assert record.u_sum_hartree == pytest.approx(-77.0125445, abs=1e-5)
        assert record.g_conc_hartree == pytest.approx(-77.0375, abs=1e-3)
        assert record.s_cal_mol_k == pytest.approx(51.982, abs=0.01)
        # g_kcal should prefer g_conc (since it's available)
        expected_g_kcal = -77.0375 * HARTREE_TO_KCAL
        assert record.g_kcal == pytest.approx(expected_g_kcal, abs=0.1)

    def test_parse_legacy_format(self, tmp_dir):
        """Parse legacy g(total)/h(total) format as fallback."""
        sum_file = tmp_dir / "legacy.sum"
        sum_file.write_text(LEGACY_SUM)

        record = parse_shermo_sum(sum_file)
        assert record.g_sum_hartree == pytest.approx(-77.0362988, abs=1e-5)
        assert record.h_sum_hartree == pytest.approx(-77.0116003, abs=1e-5)
        assert record.u_sum_hartree == pytest.approx(-77.0125445, abs=1e-5)
        assert record.s_cal_mol_k == pytest.approx(51.982, abs=0.01)

    def test_parse_missing_fields_raises(self, tmp_dir):
        """Raise RuntimeError when required fields are missing."""
        sum_file = tmp_dir / "empty.sum"
        sum_file.write_text("Temperature: 298.15 K\n")

        with pytest.raises(RuntimeError, match="missing thermodynamic fields"):
            parse_shermo_sum(sum_file)

    def test_g_kcal_without_conc_uses_g_sum(self, tmp_dir):
        """When no concentration-based G, g_kcal uses g_sum_hartree."""
        content = """\
 Temperature:     298.150 K
 Total S:      217.494 J/mol/K      51.982 cal/mol/K    -TS:   -15.499 kcal/mol
 Sum of electronic energy and thermal correction to U:         -77.0125445 a.u.
 Sum of electronic energy and thermal correction to H:         -77.0116003 a.u.
 Sum of electronic energy and thermal correction to G:         -77.0362988 a.u.
"""
        sum_file = tmp_dir / "no_conc.sum"
        sum_file.write_text(content)

        record = parse_shermo_sum(sum_file)
        assert record.g_conc_hartree is None
        expected_g_kcal = -77.0362988 * HARTREE_TO_KCAL
        assert record.g_kcal == pytest.approx(expected_g_kcal, abs=0.1)


# ---------------------------------------------------------------------------
# 2. run_shermo_record
# ---------------------------------------------------------------------------

class TestRunShermoRecord:
    def test_non_toxic_path_feeds_auto_continue_input(self, tmp_dir, monkeypatch):
        """Shermo wrapper should send newlines instead of inheriting terminal stdin."""
        from rph_core.utils.thermo import run_shermo_record

        freq_output = tmp_dir / "freq.log"
        freq_output.write_text("freq", encoding="utf-8")
        output_file = tmp_dir / "test.sum"
        calls = {}

        def _fake_run(args, **kwargs):
            calls["args"] = args
            calls.update(kwargs)
            return subprocess.CompletedProcess(args=args, returncode=0, stdout=SHERMO_26_SUM, stderr="")

        monkeypatch.setattr("rph_core.utils.thermo.shermo.subprocess.run", _fake_run)

        record = run_shermo_record(
            shermo_bin=Path("Shermo"),
            freq_output=freq_output,
            sp_energy=-77.0,
            output_file=output_file,
            options=ShermoOptions(temperature_k=298.15),
        )

        assert calls["input"] == "\n\n"
        assert calls["capture_output"] is True
        assert calls["text"] is True
        assert calls["cwd"] is None
        assert output_file.exists()
        assert record.g_sum_hartree == pytest.approx(-77.0362988, abs=1e-5)

    def test_toxic_path_feeds_auto_continue_input(self, tmp_dir, monkeypatch):
        """Toxic-path sandboxed Shermo execution should also stay non-interactive."""
        from rph_core.utils.thermo import run_shermo_record

        freq_dir = tmp_dir / "[toxic]"
        freq_dir.mkdir(parents=True, exist_ok=True)
        freq_output = freq_dir / "freq.log"
        freq_output.write_text("freq", encoding="utf-8")
        output_file = tmp_dir / "test_toxic.sum"
        calls = {}

        def _fake_run(args, **kwargs):
            calls["args"] = args
            calls.update(kwargs)
            return subprocess.CompletedProcess(args=args, returncode=0, stdout=SHERMO_26_SUM, stderr="")

        monkeypatch.setattr("rph_core.utils.thermo.shermo.subprocess.run", _fake_run)

        record = run_shermo_record(
            shermo_bin=Path("Shermo"),
            freq_output=freq_output,
            sp_energy=-77.0,
            output_file=output_file,
            options=ShermoOptions(temperature_k=298.15),
        )

        assert calls["input"] == "\n\n"
        assert calls["capture_output"] is True
        assert calls["text"] is True
        assert calls["cwd"] is not None
        assert Path(calls["args"][1]).name == "freq.log"
        assert output_file.exists()
        assert record.g_sum_hartree == pytest.approx(-77.0362988, abs=1e-5)


# ---------------------------------------------------------------------------
# 3. ThermoRecord + write/read roundtrip
# ---------------------------------------------------------------------------

class TestThermoJsonRoundtrip:
    def test_write_read_roundtrip(self, tmp_dir):
        """Write thermo.json and read it back with all fields preserved."""
        record = ThermoRecord(
            g_kcal=-48300.0,
            h_kcal=-48290.0,
            u_kcal=-48295.0,
            s_cal_mol_k=52.0,
            g_sum_hartree=-77.0362988,
            h_sum_hartree=-77.0116003,
            u_sum_hartree=-77.0125445,
            g_conc_hartree=-77.0375,
            temperature_K=298.15,
            pressure_atm=1.0,
            concentration=None,
            source="shermo_sum",
            source_file=Path("/fake/sum.sum"),
        )
        json_path = tmp_dir / "thermo.json"
        write_thermo_json(record, json_path)

        # Verify file was created
        assert json_path.exists()

        # Read back
        loaded = read_thermo_json(json_path)
        assert loaded is not None
        assert loaded.g_kcal == pytest.approx(-48300.0)
        assert loaded.h_kcal == pytest.approx(-48290.0)
        assert loaded.s_cal_mol_k == pytest.approx(52.0)
        assert loaded.g_sum_hartree == pytest.approx(-77.0362988, abs=1e-5)
        assert loaded.g_conc_hartree == pytest.approx(-77.0375, abs=1e-3)
        assert loaded.temperature_K == 298.15

    def test_json_has_legacy_G_and_g_keys(self, tmp_dir):
        """thermo.json must contain both 'G' and 'g' for backward compat."""
        record = ThermoRecord(
            g_kcal=-12345.6,
            h_kcal=-12340.0,
            u_kcal=-12342.0,
            s_cal_mol_k=50.0,
            g_sum_hartree=-19.7,
            h_sum_hartree=-19.69,
            u_sum_hartree=-19.695,
            g_conc_hartree=None,
            temperature_K=298.15,
            pressure_atm=1.0,
            concentration=None,
            source="test",
            source_file=None,
        )
        json_path = tmp_dir / "thermo.json"
        write_thermo_json(record, json_path)

        with open(json_path) as f:
            data = json.load(f)

        assert "G" in data
        assert "g" in data
        assert data["G"] == pytest.approx(-12345.6)
        assert data["g"] == pytest.approx(-12345.6)
        assert data["g_kcal"] == pytest.approx(-12345.6)

    def test_read_v0_schema_legacy(self, tmp_dir):
        """Read legacy thermo.json with only G/g keys (no g_kcal)."""
        json_path = tmp_dir / "thermo_legacy.json"
        legacy_data = {
            "G": -48300.0,
            "g": -48300.0,
            "temperature_K": 298.15,
        }
        json_path.write_text(json.dumps(legacy_data))

        record = read_thermo_json(json_path)
        assert record is not None
        assert record.g_kcal == pytest.approx(-48300.0)

    def test_read_nonexistent_returns_none(self, tmp_dir):
        """read_thermo_json returns None for missing files."""
        result = read_thermo_json(tmp_dir / "nonexistent.json")
        assert result is None

    def test_read_invalid_json_returns_none(self, tmp_dir):
        """read_thermo_json returns None for corrupt JSON."""
        bad_file = tmp_dir / "bad.json"
        bad_file.write_text("not valid json{{{")
        result = read_thermo_json(bad_file)
        assert result is None


# ---------------------------------------------------------------------------
# 3. derive_thermo_json_from_sum
# ---------------------------------------------------------------------------

class TestDeriveThermoJsonFromSum:
    def test_derive_creates_thermo_json(self, tmp_dir):
        """derive_thermo_json_from_sum creates valid thermo.json from .sum."""
        sum_file = tmp_dir / "test.sum"
        sum_file.write_text(SHERMO_26_SUM)
        output = tmp_dir / "output_thermo.json"

        result = derive_thermo_json_from_sum(sum_file, output)

        assert result == output
        assert output.exists()

        loaded = read_thermo_json(output)
        assert loaded is not None
        assert loaded.g_kcal < 0
        assert loaded.g_sum_hartree == pytest.approx(-77.0362988, abs=1e-5)

    def test_derive_creates_parent_dirs(self, tmp_dir):
        """derive_thermo_json_from_sum creates parent directories."""
        sum_file = tmp_dir / "test.sum"
        sum_file.write_text(SHERMO_26_SUM)
        output = tmp_dir / "nested" / "deep" / "thermo.json"

        derive_thermo_json_from_sum(sum_file, output)
        assert output.exists()


# ---------------------------------------------------------------------------
# 4. ensure_thermo_json_from_entry (cache self-healing)
# ---------------------------------------------------------------------------

class TestEnsureThermoJsonFromEntry:
    def test_existing_valid_thermo_json(self, tmp_dir):
        """Return existing thermo.json when valid."""
        entry = tmp_dir / "cache_entry"
        entry.mkdir()

        # Write a valid thermo.json
        thermo_data = {
            "g_kcal": -48300.0,
            "h_kcal": -48290.0,
            "s_cal_mol_k": 52.0,
            "G": -48300.0,
            "g": -48300.0,
        }
        (entry / "thermo.json").write_text(json.dumps(thermo_data))

        result = ensure_thermo_json_from_entry(entry)
        assert result is not None
        assert result.name == "thermo.json"

    def test_self_heal_from_sum_file(self, tmp_dir):
        """Auto-generate thermo.json from .sum when missing."""
        entry = tmp_dir / "cache_entry"
        dft_dir = entry / "finalDFT"
        dft_dir.mkdir(parents=True)

        # Place a .sum file
        (dft_dir / "conf_000_Shermo.sum").write_text(SHERMO_26_SUM)

        result = ensure_thermo_json_from_entry(entry)
        assert result is not None
        assert result.name == "thermo.json"
        assert (entry / "thermo.json").exists()

        # Verify content
        loaded = read_thermo_json(result)
        assert loaded is not None
        assert loaded.g_kcal < 0

    def test_self_heal_legacy_dft_dir(self, tmp_dir):
        """Auto-heal from 'dft' (legacy) directory name."""
        entry = tmp_dir / "cache_entry"
        dft_dir = entry / "dft"
        dft_dir.mkdir(parents=True)

        (dft_dir / "mol_Shermo.sum").write_text(SHERMO_26_SUM)

        result = ensure_thermo_json_from_entry(entry)
        assert result is not None
        assert (entry / "thermo.json").exists()

    def test_no_sum_file_returns_none(self, tmp_dir):
        """Return None when no .sum file available for self-healing."""
        entry = tmp_dir / "cache_entry"
        entry.mkdir()
        # No thermo.json, no .sum file

        result = ensure_thermo_json_from_entry(entry)
        assert result is None

    def test_corrupt_existing_regenerates(self, tmp_dir):
        """Regenerate thermo.json when existing one is corrupt/empty."""
        entry = tmp_dir / "cache_entry"
        dft_dir = entry / "finalDFT"
        dft_dir.mkdir(parents=True)

        # Write corrupt thermo.json (no g_kcal, no G, no g)
        (entry / "thermo.json").write_text('{"temperature_K": 298.15}')

        # Place a good .sum file
        (dft_dir / "conf_000_Shermo.sum").write_text(SHERMO_26_SUM)

        result = ensure_thermo_json_from_entry(entry)
        assert result is not None
        loaded = read_thermo_json(result)
        assert loaded is not None
        assert loaded.g_kcal < 0


# ---------------------------------------------------------------------------
# 5. derive_shermo_summary
# ---------------------------------------------------------------------------

class TestDeriveShermoSummary:
    def test_creates_shermo_summary_json(self, tmp_dir):
        """derive_shermo_summary creates a shermo_summary.json."""
        sum_file = tmp_dir / "test.sum"
        sum_file.write_text(SHERMO_26_SUM)
        output = tmp_dir / "shermo_summary.json"

        derive_shermo_summary(sum_file, output, molecule_type="precursor")

        assert output.exists()
        with open(output) as f:
            data = json.load(f)

        assert "g_precursor" in data
        assert data["g_precursor"] < 0

    def test_molecule_type_key(self, tmp_dir):
        """Use molecule_type as key prefix in summary."""
        sum_file = tmp_dir / "test.sum"
        sum_file.write_text(SHERMO_26_SUM)
        output = tmp_dir / "shermo_summary.json"

        derive_shermo_summary(sum_file, output, molecule_type="ylide")

        with open(output) as f:
            data = json.load(f)

        assert "g_ylide" in data


# ---------------------------------------------------------------------------
# 6. ShermoOptions
# ---------------------------------------------------------------------------

class TestShermoOptions:
    def test_default_construction(self):
        """ShermoOptions with minimal args."""
        opts = ShermoOptions(temperature_k=298.15)
        assert opts.temperature_k == 298.15
        assert opts.pressure_atm is None
        assert opts.conc is None

    def test_from_config(self):
        """ShermoOptions.from_config extracts parameters from config dict."""
        config = {
            "thermo": {
                "temperature_k": 350.0,
                "pressure_atm": 2.0,
            }
        }
        # Explicit temperature_k takes priority over config
        opts = ShermoOptions.from_config(config, temperature_k=298.15)
        assert opts.temperature_k == 298.15
        assert opts.pressure_atm == 2.0

        # Without explicit temperature_k, uses config value
        opts2 = ShermoOptions.from_config(config)
        assert opts2.temperature_k == 350.0


# ---------------------------------------------------------------------------
# 7. Import smoke test
# ---------------------------------------------------------------------------

class TestImportSmoke:
    def test_all_public_symbols_importable(self):
        """Verify all public symbols from rph_core.utils.thermo are importable."""
        from rph_core.utils.thermo import (
            ThermoRecord,
            ShermoOptions,
            parse_shermo_sum,
            run_shermo_record,
            write_thermo_json,
            read_thermo_json,
            derive_thermo_json_from_sum,
            ensure_thermo_json_from_entry,
            derive_shermo_summary,
        )
        assert ThermoRecord is not None
        assert ShermoOptions is not None
        assert callable(parse_shermo_sum)
        assert callable(write_thermo_json)
        assert callable(read_thermo_json)
