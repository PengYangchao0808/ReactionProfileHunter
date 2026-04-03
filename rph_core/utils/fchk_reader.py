"""
Gaussian FCHK File Reader
=========================

Parse Gaussian fchk files to extract orbital energies and atomic properties.
Robust parsing for restricted and unrestricted calculations.

Author: RPH Team
Date: 2026-02-03
"""

import re
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any, Protocol
import logging

from rph_core.utils.constants import HARTREE_TO_EV

logger = logging.getLogger(__name__)


_FLOAT_TOKEN_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\d*\.?\d+)(?:[DdEe][+-]?\d+)?$")


def _parse_fortran_float(token: str) -> Optional[float]:
    token = (token or "").strip()
    if not token:
        return None
    if not _FLOAT_TOKEN_RE.match(token):
        return None
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


class _HeaderMatcher(Protocol):
    def end(self) -> int: ...


class FCHKReader:
    """Parse Gaussian fchk files for orbital energies and charges."""

    def __init__(self, fchk_path: Path):
        """Initialize FCHK reader.

        Args:
            fchk_path: Path to fchk file
        """
        self.fchk_path = Path(fchk_path)
        self._raw_data: Dict[str, Any] = {}

    def parse(self) -> bool:
        """Parse fchk file.

        Returns:
            True if parsing successful, False otherwise
        """
        if not self.fchk_path.exists():
            logger.warning(f"FCHK file not found: {self.fchk_path}")
            return False

        try:
            with open(self.fchk_path, 'r', encoding='utf-8') as f:
                content = f.read()
            self._raw_data = self._parse_content(content)
            return True
        except Exception as e:
            logger.error(f"Failed to parse fchk {self.fchk_path}: {e}")
            return False

    def _parse_content(self, content: str) -> Dict[str, Any]:
        """Parse fchk content into structured data.

        Args:
            content: Raw fchk file content

        Returns:
            Parsed data dictionary
        """
        result = {}

        # Extract orbital energies
        result['orbital_energies'] = self._extract_orbital_energies(content)

        # Extract number of electrons
        result['n_alpha'], result['n_beta'] = self._extract_electron_counts(content)

        # Extract atomic numbers and charges
        result['atomic_numbers'] = self._extract_atomic_numbers(content)
        result['mulliken_charges'] = self._extract_mulliken_charges(content)

        # Extract CM5 charges if present
        result['cm5_charges'] = self._extract_cm5_charges(content)

        # Extract coordinates
        result['coordinates'] = self._extract_coordinates(content)

        return result

    def _collect_array_block(
        self,
        content: str,
        header_match: _HeaderMatcher,
        expected_count: int,
        kind: str,
    ) -> List[Any]:
        """Collect values following an fchk array header until expected_count.

        Gaussian fchk arrays follow a header line like:
            "Mulliken Charges  R  N=  42"
        and then a whitespace-separated stream of values that may span
        multiple lines.

        Args:
            content: Full fchk file content
            header_match: Regex match object for the header line
            expected_count: Number of values to collect
            kind: "float" or "int"

        Returns:
            List of parsed values (may be shorter than expected_count)
        """
        if expected_count <= 0:
            return []

        start = header_match.end()
        nl = content.find("\n", start)
        if nl != -1:
            start = nl + 1

        values: List[Any] = []
        for line in content[start:].splitlines():
            if len(values) >= expected_count:
                break
            if not line.strip():
                continue
            tokens = line.strip().split()
            for tok in tokens:
                if len(values) >= expected_count:
                    break
                if kind == "float":
                    # Skip plain integers to avoid accidental row indices.
                    if not any(c in tok for c in (".", "D", "d", "E", "e")):
                        continue
                    v = _parse_fortran_float(tok)
                    if v is None:
                        continue
                    values.append(v)
                else:
                    try:
                        values.append(int(tok))
                    except ValueError:
                        continue

        return values

    def _extract_orbital_energies(self, content: str) -> Optional[List[float]]:
        """Extract orbital energies from fchk (restricted or unrestricted).

        Args:
            content: Raw fchk content

        Returns:
            List of orbital energies in Hartree, or None
        """
        # Try "Alpha Orbital Energies" first (unrestricted)
        match = re.search(r"Alpha Orbital Energies\s+R\s+N=\s*(\d+)", content)

        if match:
            n_values = int(match.group(1))
            energies = self._collect_array_block(content, match, n_values, kind="float")
            if len(energies) >= n_values:
                return energies[:n_values]
            return None
        
        # Try "Orbital Energies" (restricted - no Alpha/Beta prefix)
        match = re.search(r"Orbital Energies\s+R\s+N=\s*(\d+)", content)

        if match:
            n_values = int(match.group(1))
            energies = self._collect_array_block(content, match, n_values, kind="float")
            if len(energies) >= n_values:
                return energies[:n_values]
            return None
        
        logger.debug("No Orbital Energies found in fchk")
        return None

    def _parse_fortran_array(self, values: List[str], expected_count: int) -> List[float]:
        """Backward-compatible float parsing helper.

        Note: new parsing uses _collect_array_block(); keep this for any
        legacy callers.
        """
        result: List[float] = []
        for v in values:
            parsed = _parse_fortran_float(v)
            if parsed is not None:
                result.append(parsed)
        return result

    def _extract_electron_counts(self, content: str) -> Tuple[int, int]:
        """Extract number of alpha and beta electrons.

        Supports both formats:
        - "Number of alpha electrons  I  25" (fchk style, no '=')
        - "Number of alpha electrons = 25" (with '=')

        Args:
            content: Raw fchk content

        Returns:
            Tuple of (n_alpha, n_beta)
        """
        n_alpha, n_beta = 0, 0

        # Try fchk-style format first (no '='): "Number of alpha electrons  I  25"
        # The 'I' indicates integer type in fchk format
        alpha_match = re.search(r'Number of alpha electrons\s+I\s+(\d+)', content)
        if not alpha_match:
            # Try with '=' format
            alpha_match = re.search(r'Number of alpha electrons\s*=\s*(\d+)', content)
        
        if alpha_match:
            n_alpha = int(alpha_match.group(1))

        # Same for beta
        beta_match = re.search(r'Number of beta electrons\s+I\s+(\d+)', content)
        if not beta_match:
            beta_match = re.search(r'Number of beta electrons\s*=\s*(\d+)', content)
        
        if beta_match:
            n_beta = int(beta_match.group(1))

        return n_alpha, n_beta

    def _extract_atomic_numbers(self, content: str) -> Optional[List[int]]:
        """Extract atomic numbers from fchk.

        Args:
            content: Raw fchk content

        Returns:
            List of atomic numbers, or None
        """
        match = re.search(r"Atomic numbers\s+I\s+N=\s*(\d+)", content)
        if not match:
            logger.debug("No Atomic numbers found in fchk")
            return None

        n_values = int(match.group(1))
        
        values = self._collect_array_block(content, match, n_values, kind="int")
        if len(values) >= n_values:
            return values[:n_values]
        return None

    def _extract_mulliken_charges(self, content: str) -> Optional[List[float]]:
        """Extract Mulliken charges from fchk.

        Args:
            content: Raw fchk content

        Returns:
            List of Mulliken charges, or None
        """
        match = re.search(r"Mulliken Charges\s+R\s+N=\s*(\d+)", content)
        if not match:
            logger.debug("No Mulliken Charges found in fchk")
            return None

        n_values = int(match.group(1))
        
        values = self._collect_array_block(content, match, n_values, kind="float")
        if len(values) >= n_values:
            return values[:n_values]
        return None

    def _extract_cm5_charges(self, content: str) -> Optional[List[float]]:
        """Extract CM5 charges from fchk.

        Args:
            content: Raw fchk content

        Returns:
            List of CM5 charges, or None
        """
        match = re.search(r"CM5 Charges\s+R\s+N=\s*(\d+)", content)
        if not match:
            logger.debug("No CM5 Charges found in fchk")
            return None

        n_values = int(match.group(1))
        
        values = self._collect_array_block(content, match, n_values, kind="float")
        if len(values) >= n_values:
            return values[:n_values]
        return None

    def _extract_coordinates(self, content: str) -> Optional[List[List[float]]]:
        """Extract Cartesian coordinates from fchk.

        Args:
            content: Raw fchk content

        Returns:
            List of [x, y, z] coordinates, or None
        """
        match = re.search(r"Current Cartesian coordinates\s+R\s+N=\s*(\d+)", content)
        if not match:
            logger.debug("No Current Cartesian coordinates found in fchk")
            return None

        n_values = int(match.group(1))
        expected_coords = n_values // 3  # 3 values per atom
        
        flat_values = self._collect_array_block(content, match, n_values, kind="float")
        if len(flat_values) < 3:
            return None
        coords: List[List[float]] = []
        for i in range(0, min(len(flat_values), n_values), 3):
            if i + 2 < len(flat_values):
                coords.append([flat_values[i], flat_values[i + 1], flat_values[i + 2]])
        if len(coords) >= expected_coords:
            return coords[:expected_coords]
        return coords

    def get_homo_lumo(self) -> Tuple[Optional[float], Optional[float], str]:
        """Get HOMO and LUMO energies in eV.

        Returns:
            Tuple of (homo_eV, lumo_eV, unit)
            homo_eV: HOMO energy in eV
            lumo_eV: LUMO energy in eV
            unit: Unit of energies ("eV")
        """
        energies = self._raw_data.get('orbital_energies')
        if not energies:
            return None, None, "eV"

        n_alpha = self._raw_data.get('n_alpha', 0)
        n_beta = self._raw_data.get('n_beta', 0)

        # HOMO is the highest occupied orbital
        # For restricted: use n_alpha (or n_beta if available)
        # For unrestricted: use alpha orbitals
        n_occ = n_alpha if n_alpha > 0 else n_beta
        
        if n_occ > 0 and len(energies) >= n_occ:
            homo_idx = n_occ - 1
            homo = energies[homo_idx] * HARTREE_TO_EV
        else:
            homo = None

        # LUMO is the first unoccupied orbital (next index)
        if n_occ > 0 and len(energies) > n_occ:
            lumo_idx = n_occ
            lumo = energies[lumo_idx] * HARTREE_TO_EV
        else:
            lumo = None

        return homo, lumo, "eV"

    def get_cdft_indices(self) -> Dict[str, Optional[float]]:
        """Calculate CDFT indices from orbital energies.

        Returns:
            Dictionary with keys:
            - eps_homo: HOMO energy (eV)
            - eps_lumo: LUMO energy (eV)
            - mu: Chemical potential = (HOMO + LUMO) / 2 (eV)
            - eta: Hardness = LUMO - HOMO (eV)
            - omega: Parr electrophilicity index (eV)
        """
        homo, lumo, _ = self.get_homo_lumo()

        result = {
            'eps_homo': homo,
            'eps_lumo': lumo,
            'mu': None,
            'eta': None,
            'omega': None
        }

        if homo is not None and lumo is not None:
            result['mu'] = (homo + lumo) / 2.0
            result['eta'] = lumo - homo
            if result['eta'] > 0:
                result['omega'] = (result['mu'] ** 2) / (2.0 * result['eta'])

        return result

    def get_charges(self, charge_type: str = "MULLIKEN") -> Optional[List[float]]:
        """Get atomic charges of specified type.

        Args:
            charge_type: Type of charges ("MULLIKEN", "CM5")

        Returns:
            List of charges, or None if not available
        """
        if charge_type.upper() == "CM5":
            return self._raw_data.get('cm5_charges')
        else:
            return self._raw_data.get('mulliken_charges')

    def get_atomic_numbers(self) -> Optional[List[int]]:
        """Get atomic numbers.

        Returns:
            List of atomic numbers
        """
        return self._raw_data.get('atomic_numbers')

    def get_coordinates(self) -> Optional[List[List[float]]]:
        """Get Cartesian coordinates.

        Returns:
            List of [x, y, z] coordinates in Bohr
        """
        return self._raw_data.get('coordinates')


def read_fchk_orbital_energies(fchk_path: Path) -> Tuple[Optional[float], Optional[float], str]:
    """Convenience function to read HOMO/LUMO from fchk.

    Args:
        fchk_path: Path to fchk file

    Returns:
        Tuple of (homo_eV, lumo_eV, unit)
    """
    reader = FCHKReader(fchk_path)
    if reader.parse():
        return reader.get_homo_lumo()
    return None, None, "eV"


def read_fchk_cdft_indices(fchk_path: Path) -> Dict[str, Optional[float]]:
    """Convenience function to read CDFT indices from fchk.

    Args:
        fchk_path: Path to fchk file

    Returns:
        Dictionary with cdft indices
    """
    reader = FCHKReader(fchk_path)
    if reader.parse():
        return reader.get_cdft_indices()
    return {'eps_homo': None, 'eps_lumo': None, 'mu': None, 'eta': None, 'omega': None}


def read_fchk_charges(fchk_path: Path, charge_type: str = "MULLIKEN") -> Optional[List[float]]:
    """Convenience function to read charges from fchk.

    Args:
        fchk_path: Path to fchk file
        charge_type: Type of charges ("MULLIKEN", "CM5")

    Returns:
        List of charges
    """
    reader = FCHKReader(fchk_path)
    if reader.parse():
        return reader.get_charges(charge_type)
    return None
