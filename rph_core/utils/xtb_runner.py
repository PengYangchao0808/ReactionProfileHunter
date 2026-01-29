"""
XTB Runner Module
=================

Robust execution wrapper for GFN-xTB calculations.

Provides a clean interface for running XTB optimizations with support for:
- Multi-level executable path verification
- Constrained optimizations (bond length constraints)
- Solvent effects (ALPB model)
- Charge and spin state control
- Parallel execution with proper error handling

Author: QCcalc Team
Date: 2026-01-15
"""

import shutil
import subprocess
from pathlib import Path
from typing import Dict, Optional

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.data_types import QCResult


class XTBRunner(LoggerMixin):
    """
    GFN-xTB calculation executor with robust error handling.

    Features:
    - Automatic executable discovery (config -> PATH -> fallback paths)
    - Bond length constraint support
    - Solvent model integration (GFN2-ALPB)
    - Comprehensive logging and error reporting
    - Output validation with file size checks

    Example:
        >>> runner = XTBRunner(config)
        >>> result = runner.optimize(
        ...     structure=Path("input.xyz"),
        ...     constraints={"1 5": 2.1},
        ...     solvent="acetone",
        ...     charge=0,
        ...     uhf=0
        ... )
        >>> if result.success:
        ...     print(f"Energy: {result.energy}")
    """

    # Hardcoded fallback paths for common installation locations
    FALLBACK_PATHS = [
        "/opt/xtb/bin/xtb",
        "/usr/local/bin/xtb"
    ]

    def __init__(self, config: dict, work_dir: Optional[Path] = None):
        """
        Initialize XTB runner with configuration.

        Args:
            config: Configuration dictionary containing:
                - executables.xtb.path: Preferred xtb executable path
                - resources.nproc: Number of parallel cores
            work_dir: Working directory for calculations (optional)

        Raises:
            RuntimeError: If XTB executable cannot be found
        """
        self.config = config
        self.xtb_path = self._verify_executable()

        if work_dir is None:
            self.work_dir = Path.cwd()
        else:
            self.work_dir = Path(work_dir)

        self.work_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"XTBRunner initialized: {self.xtb_path}")
        self.logger.debug(f"Working directory: {self.work_dir}")

    def _verify_executable(self) -> str:
        """
        Verify XTB executable location with multi-level fallback.

        Search strategy:
        1. Check config['executables']['xtb']['path']
        2. Try 'xtb' in system PATH
        3. Try hardcoded fallback paths
        4. Raise RuntimeError if all attempts fail

        Returns:
            Absolute path to xtb executable

        Raises:
            RuntimeError: If XTB cannot be found anywhere
        """
        # 1. Check configured path
        config_path = self.config.get('executables', {}).get('xtb', {}).get('path')
        if config_path:
            xtb_path = shutil.which(config_path)
            if xtb_path:
                self.logger.info(f"XTB found from config: {xtb_path}")
                return xtb_path
            else:
                self.logger.warning(f"Configured path not found: {config_path}")

        # 2. Try system PATH
        system_path = shutil.which('xtb')
        if system_path:
            self.logger.info(f"XTB found in system PATH: {system_path}")
            return system_path

        # 3. Try fallback paths
        for fallback in self.FALLBACK_PATHS:
            if Path(fallback).exists() and shutil.which(fallback):
                self.logger.info(f"XTB found at fallback path: {fallback}")
                return fallback

        # 4. All attempts failed
        error_msg = (
            "XTB executable not found. Searched locations:\n"
            f"  - Config path: {config_path}\n"
            f"  - System PATH\n"
            f"  - Fallback paths: {self.FALLBACK_PATHS}\n"
            "Please ensure XTB is installed and accessible."
        )
        self.logger.error(error_msg)
        raise RuntimeError(error_msg)

    def _write_constraint_input(self, constraints: Dict[str, float]) -> Path:
        """
        Write XTB constraint file from dictionary.

        Converts bond length constraints to XTB $constrain block format.

        Args:
            constraints: Dictionary with format {"atom1_idx atom2_idx": distance}
                        (e.g., {"1 5": 2.1, "10 15": 2.3})

        Returns:
            Path to constraint file

        Example:
            >>> constraints = {"1 5": 2.1, "10 15": 2.3}
            >>> constraint_file = self._write_constraint_input(constraints)
            >>> # Generates: $constrain ... distance: 1, 5, 2.100 ... $end
        """
        constraint_file = self.work_dir / "constraints.inp"

        with constraint_file.open('w') as f:
            f.write("$constrain\n")
            f.write("  force constant=1.0\n")

            for atoms, distance in constraints.items():
                # Parse atom indices from string (e.g., "1 5" -> [1, 5])
                atom1, atom2 = map(int, atoms.split())
                f.write(f"  distance: {atom1}, {atom2}, {distance:.3f}\n")

            f.write("$end\n")

        self.logger.debug(f"Wrote constraint file: {constraint_file}")
        return constraint_file

    def optimize(
        self,
        structure: Path,
        constraints: Optional[Dict[str, float]] = None,
        solvent: Optional[str] = None,
        charge: int = 0,
        uhf: int = 0
    ) -> QCResult:
        """
        Run XTB geometry optimization.

        Args:
            structure: Input XYZ file path
            constraints: Optional bond length constraints
                         Format: {"atom1_idx atom2_idx": distance}
            solvent: Solvent name for ALPB model (e.g., "acetone")
            charge: Molecular charge (default: 0)
            uhf: Number of unpaired electrons (default: 0)

        Returns:
            QCResult object with optimization results
        """
        structure = Path(structure)
        if not structure.exists():
            return QCResult(
                success=False,
                error_message=f"XTB input structure not found: {structure}"
            )

        local_structure = self.work_dir / structure.name
        if local_structure.resolve() != structure.resolve():
            local_structure.write_text(structure.read_text())

        # Build XTB command
        cmd = [self.xtb_path, local_structure.name, "--opt"]

        # Add parallel processing parameter
        nproc = self.config.get('resources', {}).get('nproc', 1)
        if nproc == 1:
            self.logger.warning("Single-core optimization (nproc=1). "
                                "Consider increasing nproc for faster performance.")
        cmd.extend(["-P", str(nproc)])

        # Add electronic parameters
        cmd.extend(["--chrg", str(charge)])
        if uhf > 0:
            cmd.extend(["--uhf", str(uhf)])

        # Add solvent model
        if solvent:
            cmd.extend(["--gfn", "2", "--alpb", str(solvent)])
            self.logger.debug(f"Using solvent model: ALPB-{solvent}")

        # Add constraints if provided
        constraint_file = None
        if constraints:
            constraint_file = self._write_constraint_input(constraints)
            cmd.extend(["--input", constraint_file.name])
            self.logger.debug(f"Applying constraints: {constraints}")

        self.logger.info(f"Running XTB optimization: {' '.join(cmd)}")

        try:
            # Execute XTB
            result = subprocess.run(
                cmd,
                cwd=self.work_dir,
                capture_output=True,
                text=True,
                check=True
            )

            # Check for output file
            output_file = self.work_dir / "xtbopt.xyz"

            if not output_file.exists():
                return QCResult(
                    success=False,
                    error_message="Optimization completed but xtbopt.xyz not found"
                )

            # Check for empty file
            if output_file.stat().st_size == 0:
                return QCResult(
                    success=False,
                    error_message="xtbopt.xyz file is empty"
                )

            # Parse output for energy
            energy = self._parse_energy(output_file)

            self.logger.info(f"XTB optimization successful. Energy: {energy:.6f} Hartree")
            self.logger.debug(f"Output file: {output_file}")

            return QCResult(
                success=True,
                energy=energy,
                coordinates=output_file,
                converged=True,
                output_file=output_file
            )

        except subprocess.CalledProcessError as e:
            # Enhanced error logging
            error_msg = (
                f"XTB optimization failed with return code {e.returncode}\n"
                f"Command: {' '.join(e.cmd)}\n"
                f"STDERR:\n{e.stderr}\n"
                f"STDOUT (last 500 chars):\n{e.stdout[-500:] if len(e.stdout) > 500 else e.stdout}"
            )
            self.logger.error(error_msg)

            return QCResult(
                success=False,
                error_message=f"XTB failed: {e.stderr}"
            )

    def _parse_energy(self, output_file: Path) -> Optional[float]:
        """
        Parse total energy from XTB output.

        XTB writes energy information to stdout, but we can also
        check the xyz file comments for energy information.

        Args:
            output_file: Path to xtbopt.xyz file

        Returns:
            Energy in Hartree, or None if parsing fails
        """
        try:
            # XTB xyz files often have energy in the comment line (line 2)
            lines = output_file.read_text().splitlines()
            if len(lines) >= 2:
                # Try to extract energy from comment line
                comment = lines[1]
                # Look for energy pattern (e.g., "energy: -123.456")
                import re
                energy_match = re.search(r'-?\d+\.\d+', comment)
                if energy_match:
                    return float(energy_match.group())

        except Exception as e:
            self.logger.warning(f"Failed to parse energy from output: {e}")

        return None
