"""
Data Types Module
==================

Standardized data structures for ReactionProfileHunter.

Provides type-safe dataclasses for computational chemistry results,
replacing dictionary-based returns with structured objects.

Author: QCcalc Team
Date: 2026-01-15
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional


@dataclass
class QCResult:
    """
    Quantum chemistry calculation result.

    Standardizes return format for QC computations, providing
    type safety and clear semantics compared to dictionary returns.

    Attributes:
        success: Whether calculation completed successfully
        energy: Final energy value (None if calculation failed)
        coordinates: Optimized coordinates or file reference
        error_message: Human-readable error description (None if successful)
        converged: Whether optimization converged (False if failed)
        output_file: Output file path (if applicable)
        frequencies: Vibrational frequencies (if applicable)
        log_file: Gaussian .log or ORCA .out file path (for S4 parsing)
        chk_file: Gaussian .chk checkpoint file path
        fchk_file: Gaussian .fchk formatted checkpoint file path (for S4)
        qm_output_file: Generic QM output file (Gaussian=.log, ORCA=.out)
    """
    success: bool = False
    energy: Optional[float] = None
    coordinates: Optional[Any] = None
    error_message: Optional[str] = None
    converged: bool = False
    output_file: Optional[Path] = None
    frequencies: Optional[Any] = None
    homo: Optional[float] = None
    lumo: Optional[float] = None
    gap: Optional[float] = None
    dipole_moment: Optional[Any] = None
    charges: Optional[Dict[str, float]] = None
    log_file: Optional[Path] = None
    chk_file: Optional[Path] = None
    fchk_file: Optional[Path] = None
    qm_output_file: Optional[Path] = None
