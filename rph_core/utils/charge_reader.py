"""
Atomic Charge Reader
====================

Unified interface for reading atomic charges from various sources:
NBO, CM5, Mulliken. Implements charge priority and fallback logic.

Author: RPH Team
Date: 2026-02-02
"""

import re
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any
import logging

logger = logging.getLogger(__name__)


class ChargeReader:
    """Read atomic charges with priority fallback."""

    def __init__(self, config: Dict[str, Any]):
        """Initialize charge reader.

        Args:
            config: Configuration dict with charge_priority list
        """
        self.priority = config.get('charge_priority', ['NBO', 'CM5', 'MULLIKEN'])

    def read_charges(
        self,
        fchk_path: Optional[Path] = None,
        nbo_path: Optional[Path] = None,
        log_path: Optional[Path] = None
    ) -> Tuple[Optional[List[float]], str]:
        """Read charges with priority fallback.

        Args:
            fchk_path: Path to fchk file
            nbo_path: Path to NBO output file
            log_path: Path to Gaussian log file

        Returns:
            Tuple of (charges_list, charge_type)
        """
        for charge_type in self.priority:
            if charge_type.upper() == 'NBO':
                charges = self._read_nbo(nbo_path or log_path)
                if charges is not None:
                    return charges, 'NBO'
            elif charge_type.upper() == 'CM5':
                charges = self._read_cm5_from_fchk(fchk_path)
                if charges is not None:
                    return charges, 'CM5'
            elif charge_type.upper() == 'MULLIKEN':
                charges = self._read_mulliken_from_fchk(fchk_path)
                if charges is not None:
                    return charges, 'MULLIKEN'

        logger.warning(f"No charges found from priority list: {self.priority}")
        return None, 'NONE'

    def _read_nbo(self, nbo_path: Optional[Path]) -> Optional[List[float]]:
        """Read charges from NBO output.

        Args:
            nbo_path: Path to NBO output file

        Returns:
            List of NPA charges, or None
        """
        if nbo_path is None or not Path(nbo_path).exists():
            return None

        try:
            with open(nbo_path, 'r', encoding='utf-8') as f:
                content = f.read()

            # Look for Natural Population Analysis section
            pattern = r'Summary of Natural Population Analysis:\s*\n\s*-?\d+\s*\n([\d\s.\-+]+)'
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                charges_str = match.group(1).strip().split()
                charges = [float(c) for c in charges_str]
                return charges

            # Alternative: look for individual atom charges
            # Pattern: "Atom No         Charge"
            lines = content.split('\n')
            charges = []
            in_section = False

            for line in lines:
                if 'Natural Population Analysis' in line or 'Summary' in line:
                    in_section = True
                    continue
                if in_section and line.strip() == '':
                    in_section = False
                if in_section:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            charge = float(parts[-1])
                            charges.append(charge)
                        except ValueError:
                            pass

            if charges:
                return charges

            logger.debug("NBO charges not found in file")
            return None

        except Exception as e:
            logger.warning(f"Failed to read NBO charges: {e}")
            return None

    def _read_cm5_from_fchk(self, fchk_path: Optional[Path]) -> Optional[List[float]]:
        """Read CM5 charges from fchk.

        Args:
            fchk_path: Path to fchk file

        Returns:
            List of CM5 charges, or None
        """
        if fchk_path is None or not Path(fchk_path).exists():
            return None

        try:
            from .fchk_reader import FCHKReader
            reader = FCHKReader(fchk_path)
            if reader.parse():
                return reader.get_charges('CM5')
            return None
        except ImportError:
            logger.warning("FCHKReader not available for CM5 parsing")
            return None
        except Exception as e:
            logger.warning(f"Failed to read CM5 from fchk: {e}")
            return None

    def _read_mulliken_from_fchk(self, fchk_path: Optional[Path]) -> Optional[List[float]]:
        """Read Mulliken charges from fchk.

        Args:
            fchk_path: Path to fchk file

        Returns:
            List of Mulliken charges, or None
        """
        if fchk_path is None or not Path(fchk_path).exists():
            return None

        try:
            from .fchk_reader import FCHKReader
            reader = FCHKReader(fchk_path)
            if reader.parse():
                return reader.get_charges('MULLIKEN')
            return None
        except ImportError:
            logger.warning("FCHKReader not available for Mulliken parsing")
            return None
        except Exception as e:
            logger.warning(f"Failed to read Mulliken from fchk: {e}")
            return None


def read_charges_with_priority(
    fchk_path: Optional[Path] = None,
    nbo_path: Optional[Path] = None,
    log_path: Optional[Path] = None,
    charge_priority: Optional[List[str]] = None
) -> Tuple[Optional[List[float]], str]:
    """Convenience function to read charges with priority.

    Args:
        fchk_path: Path to fchk file
        nbo_path: Path to NBO output file
        log_path: Path to Gaussian log file
        charge_priority: List of charge types in priority order

    Returns:
        Tuple of (charges_list, charge_type)
    """
    config = {'charge_priority': charge_priority or ['NBO', 'CM5', 'MULLIKEN']}
    reader = ChargeReader(config)
    return reader.read_charges(fchk_path, nbo_path, log_path)
