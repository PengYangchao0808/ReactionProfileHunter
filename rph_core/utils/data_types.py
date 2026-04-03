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
from typing import Any, Dict, List, NewType, Optional, Union


MapId = NewType("MapId", int)
MolIdx = NewType("MolIdx", int)


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


@dataclass
class ScanResult:
    """
    Scan calculation result for retro scan (S2).

    Stores results from bond stretching scan calculations,
    providing structured data for TS guess selection.

    Attributes:
        success: Whether scan completed successfully
        energies: List of energy values at each scan point
        geometries: Scan point geometries (single path or list of paths)
        max_energy_index: Index of highest energy point (TS guess location)
        ts_guess_xyz: Path to TS guess structure from scan
        scan_log: Path to scan log file
    """
    success: bool = False
    energies: Optional[List[float]] = None
    geometries: Optional[Union[List[Path], Path]] = None
    max_energy_index: int = -1
    ts_guess_xyz: Optional[Path] = None
    scan_log: Optional[Path] = None


@dataclass
class PathSearchResult:
    """
    xTB Path Finder result (S2).

    Stores results from xTB meta-dynamics reaction path finder,
    providing structured data for TS guess selection.

    Attributes:
        success: Whether path search completed successfully
        path_xyz_files: List of path trajectory XYZ files
        ts_guess_xyz: Path to estimated transition state (xtbpath_ts.xyz)
        path_log: Path to path search log file
        
        # Energy information
        barrier_forward_kcal: Forward barrier height in kcal/mol
        barrier_backward_kcal: Backward barrier height in kcal/mol
        reaction_energy_kcal: Reaction energy in kcal/mol
        
        # TS quality metrics
        estimated_ts_point: Point index of estimated TS on path
        gradient_norm_at_ts: Gradient norm at estimated TS
        
        error_message: Error message if failed
    """
    success: bool = False
    path_xyz_files: Optional[List[Path]] = None
    ts_guess_xyz: Optional[Path] = None
    path_log: Optional[Path] = None
    
    # Energy information
    barrier_forward_kcal: Optional[float] = None
    barrier_backward_kcal: Optional[float] = None
    reaction_energy_kcal: Optional[float] = None
    
    # TS quality metrics
    estimated_ts_point: Optional[int] = None
    gradient_norm_at_ts: Optional[float] = None
    
    error_message: Optional[str] = None
