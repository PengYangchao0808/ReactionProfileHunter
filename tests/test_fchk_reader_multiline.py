"""
TDD Tests for FCHK Reader Multiline Parsing
============================================

Tests for robust fchk parsing including:
- Restricted orbital energies (no Alpha/Beta prefix)
- Unrestricted Alpha/Beta orbital energies
- Multiline array parsing (R N=XXX format)
- Electron counts without '=' sign
- Charges parsing (Mulliken, CM5)

These tests should FAIL initially (RED) until fchk_reader.py is fixed.

Author: RPH Team
Date: 2026-02-03
"""

import pytest
import tempfile
from pathlib import Path
from typing import List

from rph_core.utils.fchk_reader import FCHKReader, read_fchk_cdft_indices


# =============================================================================
# Synthetic FCHK Fixtures (inline strings - no external files needed)
# =============================================================================

FCHK_RESTRICTED = """Test FCHK - Restricted Calculation
Number of atoms                       I          6
Number of electrons                   I         26
Number of alpha electrons             I         13
Number of beta electrons              I         13
Number of basis functions             I         24
Orbital Energies                      R  N=24
-1.245673D+01 -1.156723D+01 -1.023456D+01 -8.765432D+00 -7.654321D+00
-6.543210D+00 -5.432109D+00 -4.321098D+00 -3.210987D+00 -2.345678D+00
-1.234567D+00 -9.876543D-01 -5.432109D-01 -2.345678D-01  1.234567D-01
  3.456789D-01  5.678901D-01  7.890123D-01  9.876543D-01  1.234567D+00
  1.456789D+00  1.678901D+00  1.901234D+00  2.123456D+00
Mulliken Charges                       R  N=6
  1  0.123456D+00
  2 -0.234567D+00
  3  0.876543D-01
  4 -0.987654D-01
  5  0.111111D+00
  6 -0.111111D+00
CM5 Charges                            R  N=6
  1  0.098765D+00
  2 -0.187654D+00
  3  0.065432D+00
  4 -0.076543D-01
  5  0.088888D+00
  6 -0.098888D+00
Current Cartesian coordinates           R  N=18
  0.123456D+01  0.234567D+00  0.345678D+00
  0.456789D+00  0.567890D+01  0.678901D+00
  0.789012D+00  0.890123D+00  0.901234D+01
  0.111111D+01  0.111111D+01  0.000000D+00
  0.222222D+01  0.000000D+00  0.000000D+00
  0.333333D+01  0.111111D+01  0.000000D+00
Atomic numbers                         I  N=6
  6  6  6  8  1  1
"""


FCHK_UNRESTRICTED = """Test FCHK - Unrestricted Calculation
Number of atoms                       I          6
Number of electrons                   I         26
Number of alpha electrons             I         13
Number of beta electrons              I         13
Number of basis functions             I         24
Alpha Orbital Energies                R  N=24
-1.245673D+01 -1.156723D+01 -1.023456D+01 -8.765432D+00 -7.654321D+00
-6.543210D+00 -5.432109D+00 -4.321098D+00 -3.210987D+00 -2.345678D+00
-1.234567D+00 -9.876543D-01 -5.432109D-01 -2.345678D-01  1.234567D-01
  3.456789D-01  5.678901D-01  7.890123D-01  9.876543D-01  1.234567D+00
  1.456789D+00  1.678901D+00  1.901234D+00  2.123456D+00
Beta Orbital Energies                 R  N=24
-1.198765D+01 -1.109876D+01 -9.876543D+00 -7.654321D+00 -6.543210D+00
-5.432109D+00 -4.321098D+00 -3.210987D+00 -2.345678D+00 -1.567890D+00
-9.876543D-01 -5.432109D-01 -2.345678D-01  1.234567D-01  3.456789D-01
  5.678901D-01  7.890123D-01  9.876543D-01  1.234567D+00  1.456789D+00
  1.678901D+00  1.901234D+00  2.123456D+00  2.345678D+00
Mulliken Charges                       R  N=6
  1  0.123456D+00
  2 -0.234567D+00
  3  0.876543D-01
  4 -0.987654D-01
  5  0.111111D+00
  6 -0.111111D+00
Atomic numbers                         I  N=6
  6  6  6  8  1  1
"""


FCHK_ELECTRON_COUNTS_NO_EQUALS = """Test FCHK - Electron counts without '='
Number of atoms                       I          4
Number of electrons                   I         16
Number of alpha electrons             I          8
Number of beta electrons              I          8
Number of basis functions             I         12
Orbital Energies                      R  N=8
-1.234567D+01 -1.023456D+01 -8.765432D+00 -6.543210D+00
-5.432109D+00 -3.210987D+00 -1.234567D+00 -5.432109D-01
Mulliken Charges                       R  N=4
  1  0.100000D+00
  2 -0.100000D+00
  3  0.050000D+00
  4 -0.050000D+00
Atomic numbers                         I  N=4
  6  6  8  1
"""


FCHK_CHARGES_ONLY = """Test FCHK - Charges only
Number of atoms                       I          3
Mulliken Charges                       R  N=3
  1  0.500000D+00
  2 -0.300000D+00
  3 -0.200000D+00
CM5 Charges                            R  N=3
  1  0.450000D+00
  2 -0.280000D+00
  3 -0.170000D+00
Atomic numbers                         I  N=3
  6  8  1
"""


# =============================================================================
# Helper Functions
# =============================================================================

def create_temp_fchk(content: str) -> Path:
    """Create a temporary fchk file with given content."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fchk', delete=False) as f:
        f.write(content)
        return Path(f.name)


# =============================================================================
# Tests for Restricted Orbital Energies
# =============================================================================

class TestRestrictedOrbitalEnergies:
    """Tests for restricted calculation orbital energies parsing."""

    def test_parse_orbital_energies_restricted(self):
        """Should parse 'Orbital Energies' without Alpha/Beta prefix."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            # Should have orbital energies (not None)
            assert reader._raw_data.get('orbital_energies') is not None, \
                "Parser should recognize 'Orbital Energies' without Alpha/Beta prefix"
            
            energies = reader._raw_data['orbital_energies']
            assert len(energies) == 24, f"Expected 24 energies, got {len(energies)}"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_homo_lumo_from_restricted(self):
        """Should calculate HOMO/LUMO from restricted orbital energies."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            homo, lumo, unit = reader.get_homo_lumo()
            
            # HOMO should be the 13th orbital (index 12, since n_alpha=13)
            assert homo is not None, "HOMO should be extracted from restricted calculation"
            assert lumo is not None, "LUMO should be extracted from restricted calculation"
            assert unit == "eV", "Energy unit should be eV"
            assert homo < 0, "HOMO energy should be negative for stable molecule"
            assert lumo > homo, "LUMO should be higher than HOMO"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_cdft_indices_from_restricted(self):
        """Should calculate CDFT indices from restricted orbital energies."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            cdft = reader.get_cdft_indices()
            
            assert cdft['eps_homo'] is not None, "eps_homo should be numeric"
            assert cdft['eps_lumo'] is not None, "eps_lumo should be numeric"
            assert cdft['mu'] is not None, "mu (chemical potential) should be calculated"
            assert cdft['eta'] is not None, "eta (hardness) should be calculated"
            assert cdft['omega'] is not None, "omega (electrophilicity) should be calculated"
        finally:
            fchk_path.unlink(missing_ok=True)


# =============================================================================
# Tests for Unrestricted Orbital Energies
# =============================================================================

class TestUnrestrictedOrbitalEnergies:
    """Tests for unrestricted calculation orbital energies parsing."""

    def test_parse_alpha_orbital_energies(self):
        """Should parse 'Alpha Orbital Energies' section."""
        fchk_path = create_temp_fchk(FCHK_UNRESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            energies = reader._raw_data.get('orbital_energies')
            assert energies is not None, "Parser should extract Alpha orbital energies"
            assert len(energies) == 24, f"Expected 24 Alpha energies, got {len(energies)}"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_homo_from_alpha_orbital(self):
        """Should use Alpha orbital energies for HOMO index."""
        fchk_path = create_temp_fchk(FCHK_UNRESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            homo, lumo, unit = reader.get_homo_lumo()
            
            # Based on n_alpha=13, HOMO should be at index 12
            assert homo is not None, "HOMO should be extracted from Alpha orbitals"
            assert homo < 0, "HOMO should be negative"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_cdft_indices_unrestricted(self):
        """Should calculate CDFT indices from unrestricted calculation."""
        fchk_path = create_temp_fchk(FCHK_UNRESTRICTED)
        try:
            cdft = read_fchk_cdft_indices(fchk_path)
            
            assert cdft['eps_homo'] is not None, "eps_homo should be numeric"
            assert cdft['eps_lumo'] is not None, "eps_lumo should be numeric"
            assert cdft['mu'] is not None, "mu should be calculated"
            assert cdft['eta'] is not None, "eta should be calculated"
            assert cdft['omega'] is not None, "omega should be calculated"
        finally:
            fchk_path.unlink(missing_ok=True)


# =============================================================================
# Tests for Electron Counts Without Equals
# =============================================================================

class TestElectronCountsWithoutEquals:
    """Tests for electron counts in fchk format (no '=' sign)."""

    def test_parse_alpha_electron_count_no_equals(self):
        """Should parse 'Number of alpha electrons' without '=' sign."""
        fchk_path = create_temp_fchk(FCHK_ELECTRON_COUNTS_NO_EQUALS)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            n_alpha = reader._raw_data.get('n_alpha')
            n_beta = reader._raw_data.get('n_beta')
            
            assert n_alpha == 8, f"Expected n_alpha=8, got {n_alpha}"
            assert n_beta == 8, f"Expected n_beta=8, got {n_beta}"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_parse_beta_electron_count_no_equals(self):
        """Should parse 'Number of beta electrons' without '=' sign."""
        fchk_path = create_temp_fchk(FCHK_ELECTRON_COUNTS_NO_EQUALS)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            n_beta = reader._raw_data.get('n_beta')
            assert n_beta == 8, f"Expected n_beta=8, got {n_beta}"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_homo_calculation_with_correct_n_alpha(self):
        """HOMO calculation should use correct electron count."""
        fchk_path = create_temp_fchk(FCHK_ELECTRON_COUNTS_NO_EQUALS)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            homo, lumo, unit = reader.get_homo_lumo()
            
            # With n_alpha=8, HOMO should be at index 7 (0-indexed)
            assert homo is not None, "HOMO should be calculated"
            energies = reader._raw_data.get('orbital_energies')
            expected_homo = energies[7] if energies else None
            if expected_homo:
                from rph_core.utils.constants import HARTREE_TO_EV
                assert abs(homo - expected_homo * HARTREE_TO_EV) < 0.01, \
                    "HOMO should match expected orbital energy"
        finally:
            fchk_path.unlink(missing_ok=True)


# =============================================================================
# Tests for Charges Parsing
# =============================================================================

class TestChargesParsing:
    """Tests for Mulliken and CM5 charges parsing."""

    def test_parse_mulliken_charges(self):
        """Should parse Mulliken charges from fchk."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            charges = reader.get_charges("MULLIKEN")
            assert charges is not None, "Mulliken charges should be parsed"
            assert len(charges) == 6, f"Expected 6 charges, got {len(charges)}"
            # Verify charges sum to approximately zero (within 0.5 for arbitrary test data)
            assert abs(sum(charges)) < 0.5, f"Charges should sum to ~0, got {sum(charges)}"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_parse_cm5_charges(self):
        """Should parse CM5 charges from fchk."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            charges = reader.get_charges("CM5")
            assert charges is not None, "CM5 charges should be parsed"
            assert len(charges) == 6, f"Expected 6 CM5 charges, got {len(charges)}"
            # Verify charges sum to approximately zero (within 0.5 for arbitrary test data)
            assert abs(sum(charges)) < 0.5, f"CM5 charges should sum to ~0, got {sum(charges)}"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_charges_default_to_mulliken(self):
        """Default charge type should be Mulliken."""
        fchk_path = create_temp_fchk(FCHK_CHARGES_ONLY)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            default_charges = reader.get_charges()
            mulliken_charges = reader.get_charges("MULLIKEN")
            
            assert default_charges == mulliken_charges, "Default should return Mulliken charges"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_charges_with_priority_cm5_first(self):
        """Should handle CM5 charges when CM5 available."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            # Test CM5 is available
            cm5 = reader.get_charges("CM5")
            mulliken = reader.get_charges("MULLIKEN")
            
            assert cm5 is not None, "CM5 charges should be available"
            assert mulliken is not None, "Mulliken charges should be available"
            # CM5 and Mulliken should typically differ
            assert cm5 != mulliken, "CM5 and Mulliken charges should differ"
        finally:
            fchk_path.unlink(missing_ok=True)


# =============================================================================
# Tests for Multiline Array Parsing
# =============================================================================

class TestMultilineArrayParsing:
    """Tests for multiline array parsing (R N=XXX format)."""

    def test_multiline_orbital_energies_count(self):
        """Should parse correct number of orbital energies from multiline array."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            energies = reader._raw_data.get('orbital_energies')
            assert energies is not None, "Orbital energies should be parsed"
            assert len(energies) == 24, \
                f"Expected 24 energies (N=24), got {len(energies)}"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_multiline_mulliken_charges_count(self):
        """Should parse correct number of Mulliken charges from multiline array."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            charges = reader.get_charges("MULLIKEN")
            assert charges is not None, "Mulliken charges should be parsed"
            assert len(charges) == 6, \
                f"Expected 6 charges (N=6), got {len(charges)}"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_multiline_cm5_charges_count(self):
        """Should parse correct number of CM5 charges from multiline array."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            charges = reader.get_charges("CM5")
            assert charges is not None, "CM5 charges should be parsed"
            assert len(charges) == 6, \
                f"Expected 6 CM5 charges (N=6), got {len(charges)}"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_coordinates_parsing(self):
        """Should parse coordinates from multiline array."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            coords = reader.get_coordinates()
            assert coords is not None, "Coordinates should be parsed"
            assert len(coords) == 6, f"Expected 6 atoms, got {len(coords)}"
            # Each coordinate should have 3 values
            for i, coord in enumerate(coords):
                assert len(coord) == 3, f"Atom {i} should have 3 coordinates"
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_atomic_numbers_parsing(self):
        """Should parse atomic numbers from multiline array."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            atomic_numbers = reader.get_atomic_numbers()
            assert atomic_numbers is not None, "Atomic numbers should be parsed"
            assert len(atomic_numbers) == 6, \
                f"Expected 6 atomic numbers, got {len(atomic_numbers)}"
            # Check expected values (C, C, C, O, H, H)
            assert sorted(atomic_numbers) == [1, 1, 6, 6, 6, 8], \
                f"Unexpected atomic numbers: {atomic_numbers}"
        finally:
            fchk_path.unlink(missing_ok=True)


# =============================================================================
# Tests for Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_missing_file_returns_none(self):
        """Should handle missing fchk file gracefully."""
        reader = FCHKReader(Path("/nonexistent/fake.fchk"))
        result = reader.parse()
        assert result is False, "parse() should return False for missing file"

    def test_empty_content_returns_none(self):
        """Should handle empty fchk content."""
        fchk_path = create_temp_fchk("")
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            assert reader._raw_data.get('orbital_energies') is None
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_fortran_d_format_conversion(self):
        """Should correctly convert Fortran D format to Python float."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            reader = FCHKReader(fchk_path)
            reader.parse()
            
            energies = reader._raw_data.get('orbital_energies')
            assert energies is not None, "Orbital energies should be parsed"
            
            # Check that D format was converted (values should be valid floats)
            for e in energies:
                assert isinstance(e, float), f"Energy should be float, got {type(e)}"
                assert not e == float('nan'), "Energy should not be NaN"
        finally:
            fchk_path.unlink(missing_ok=True)


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    """Integration tests for fchk reading workflow."""

    def test_full_fchk_read_workflow(self):
        """Should complete full fchk read workflow successfully."""
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            # Test convenience function
            cdft = read_fchk_cdft_indices(fchk_path)
            
            assert cdft['eps_homo'] is not None
            assert cdft['eps_lumo'] is not None
            assert cdft['mu'] is not None
            assert cdft['eta'] is not None
            assert cdft['omega'] is not None
            
            # Verify calculations
            if cdft['eps_homo'] is not None and cdft['eps_lumo'] is not None:
                expected_mu = (cdft['eps_homo'] + cdft['eps_lumo']) / 2
                expected_eta = cdft['eps_lumo'] - cdft['eps_homo']
                assert abs(cdft['mu'] - expected_mu) < 0.01
                assert abs(cdft['eta'] - expected_eta) < 0.01
                if expected_eta > 0:
                    expected_omega = (expected_mu ** 2) / (2 * expected_eta)
                    assert abs(cdft['omega'] - expected_omega) < 0.01
        finally:
            fchk_path.unlink(missing_ok=True)

    def test_read_charges_convenience_function(self):
        """Should work with read_fchk_charges convenience function."""
        from rph_core.utils.fchk_reader import read_fchk_charges
        
        fchk_path = create_temp_fchk(FCHK_RESTRICTED)
        try:
            mulliken = read_fchk_charges(fchk_path, "MULLIKEN")
            cm5 = read_fchk_charges(fchk_path, "CM5")
            
            assert mulliken is not None
            assert cm5 is not None
            assert len(mulliken) == 6
            assert len(cm5) == 6
        finally:
            fchk_path.unlink(missing_ok=True)
