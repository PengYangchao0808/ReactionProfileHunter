"""
Quantum Chemistry Interface Module - Sterile Pipeline Architecture
===================================================================

Replaces buggy Gaussian execution module with a robust "Sterile Pipeline"
architecture to fix I/O freezes on WSL, prevent disk overflow, and ensure
syntax compliance.

Classes:
- InputFactory: Creates sanitized Gaussian input files
- LogParser: Extracts final geometry from Gaussian logs
- LinuxSandbox: Context manager for safe execution environment
- GaussianRunner: Executes Gaussian calculations safely
- ResultHarvester: Collects and validates results
- run_gaussian_optimization: Main pipeline function

Author: QC Descriptors Team
Date: 2026-01-15
"""

import os
import re
import shutil
import tempfile
import subprocess
from pathlib import Path
from typing import List, Dict, Optional
import logging

logger = logging.getLogger(__name__)


class GaussianError(Exception):
    """Custom exception for Gaussian execution errors."""
    pass


class InputFactory:
    """
    Creates sanitized Gaussian input files.

    Enforces scientific policies and ensures syntax compliance.
    """

    @staticmethod
    def create(route, atoms, charge, mult, mem, nproc):
        """
        Generates a Gaussian input file string with strict formatting.
        """
        import re

        # 1. Basic Sanitization
        clean_route = route.encode('ascii', 'ignore').decode('ascii')
        clean_route = clean_route.replace('\r\n', '\n').replace('\r', '\n')

        # 2. Scientific Policy: Remove CalcFC (Surgical Cleanup)
        # First, remove the keyword itself (case insensitive)
        clean_route = re.sub(r'CalcFC', '', clean_route, flags=re.IGNORECASE)

        # Now clean up the artifacts left behind:

        # Case A: "(," inside parens -> "("
        clean_route = re.sub(r'\(\s*,', '(', clean_route)

        # Case B: ",)" inside parens -> ")"
        clean_route = re.sub(r',\s*\)', ')', clean_route)

        # Case C: Double commas ",," -> ","
        clean_route = re.sub(r',\s*,', ',', clean_route)

        # Case D: Empty parentheses "()" -> remove entirely
        clean_route = re.sub(r'\(\s*\)', '', clean_route)

        # Case E: Dangling "Opt=" (Critical Fix)
        # Matches "Opt=" followed by whitespace or end of string
        # Replaces with "Opt" (which implies default optimization)
        clean_route = re.sub(r'Opt=(?=\s|$)', 'Opt', clean_route)

        # 3. Collapse multiple spaces
        clean_route = re.sub(r'\s+', ' ', clean_route).strip()

        # 4. Check Policy & Log
        if "CalcFC" in route and "CalcFC" not in clean_route:
            import logging
            logging.getLogger(__name__).warning(f"Removed CalcFC keyword from route: {route}")

        # 5. Format Coordinates (Fixed Width)
        coords_lines = []
        for atom in atoms:
            line = f"{atom['symbol']:<2} {atom['x']:14.8f} {atom['y']:14.8f} {atom['z']:14.8f}"
            coords_lines.append(line)
        coords_str = "\n".join(coords_lines)

        # 6. Template Assembly
        # MUST end with exactly 2 blank lines
        template = (
            f"%chk=checkpoint.chk\n"
            f"%mem={mem}\n"
            f"%nprocshared={nproc}\n"
            f"#p {clean_route}\n"
            f"\n"
            f"RPH_Auto_Job\n"
            f"\n"
            f"{charge} {mult}\n"
            f"{coords_str}\n"
            f"\n\n"
        )

        return template


class LogParser:
    """
    Parses Gaussian log files to extract critical information.

    Focuses on robust extraction of final geometry.
    """

    @staticmethod
    def extract_final_geometry(log_content: str) -> List[Dict[str, float]]:
        """
        Extract final geometry from Gaussian log file.

        Locates the LAST "Standard orientation" or "Input orientation" block
        and returns atomic coordinates.

        Args:
            log_content: Full content of Gaussian log file

        Returns:
            List of dicts with keys 'symbol', 'x', 'y', 'z'
            Returns empty list if parsing fails
        """
        # Pattern for orientation blocks
        # Matches both "Standard orientation" and "Input orientation"
        orientation_pattern = re.compile(
            r'(?:Standard orientation|Input orientation):.*?\n'
            r'[-]+\s*\n'
            r'\s+\d+\s+\d+\s+\d+\s+\d+\s+[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+\s+.*?\n'  # header
            r'(.*?)'  # capture atoms
            r'[-]+\s*\n',
            re.DOTALL
        )

        blocks = list(orientation_pattern.finditer(log_content))

        if not blocks:
            logger.warning("No orientation blocks found in log file")
            return []

        # Get the LAST block
        last_block = blocks[-1]

        # Parse atomic coordinates
        # Format: atomic_number, atomic_type, X, Y, Z, ...
        atom_pattern = re.compile(
            r'\s+\d+\s+\d+\s+\d+\s+\d+\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+.*?\n'
        )

        atoms = []
        for match in atom_pattern.finditer(last_block.group(1)):
            x = float(match.group(1))
            y = float(match.group(2))
            z = float(match.group(3))

            # Convert atomic number to symbol (using periodic table)
            # For now, use atomic number as symbol placeholder
            # In a full implementation, you'd have a mapping
            atomic_num_line = re.search(r'\s+(\d+)\s+\d+', match.group(0))
            if atomic_num_line:
                atomic_num = int(atomic_num_line.group(1))
                symbol = LogParser._atomic_num_to_symbol(atomic_num)
            else:
                symbol = 'X'

            atoms.append({
                'symbol': symbol,
                'x': x,
                'y': y,
                'z': z
            })

        logger.debug(f"Extracted {len(atoms)} atoms from final orientation")
        return atoms

    @staticmethod
    def _atomic_num_to_symbol(atomic_num: int) -> str:
        """
        Convert atomic number to element symbol.

        Args:
            atomic_num: Atomic number (1-118)

        Returns:
            Element symbol string
        """
        periodic_table = [
            'X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
            'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
            'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
            'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
            'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
            'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
            'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
        ]

        if 1 <= atomic_num < len(periodic_table):
            return periodic_table[atomic_num]
        else:
            logger.warning(f"Unknown atomic number: {atomic_num}")
            return 'X'


class LinuxSandbox:
    """
    Context manager for safe Gaussian execution.

    Features:
    - Automatic cleanup of old sandboxes (Janitor)
    - Disk space validation
    - WSL path safety checks
    - Dedicated scratch directory
    """

    def __init__(self, base_dir: str = "/tmp", min_free_gb: float = 100.0):
        """
        Initialize sandbox context manager.

        Args:
            base_dir: Base directory for sandbox creation
            min_free_gb: Minimum required free disk space in GB

        Raises:
            RuntimeError: If base_dir is unsafe or insufficient disk space
        """
        self.base_dir = Path(base_dir)
        self.min_free_gb = min_free_gb
        self.sandbox_path: Optional[Path] = None

    def __enter__(self) -> Path:
        """
        Enter sandbox context.

        Performs safety checks, janitor cleanup, and creates sandbox.

        Returns:
            Path to sandbox directory

        Raises:
            RuntimeError: If safety checks fail
        """
        # Safety: Reject /mnt/ paths (WSL mount points cause I/O issues)
        if str(self.base_dir).startswith('/mnt/'):
            raise RuntimeError(
                f"Safety violation: base_dir '{self.base_dir}' is a WSL mount point. "
                "This causes I/O freezes on WSL. Use a native Linux path like /tmp."
            )

        # Safety: Check disk space
        free_gb = self._get_free_disk_space_gb(self.base_dir)
        if free_gb < self.min_free_gb:
            raise RuntimeError(
                f"Insufficient disk space: {free_gb:.2f} GB free, "
                f"required {self.min_free_gb} GB minimum"
            )

        # Janitor: Clean old RPH_G16_* folders (older than 6 hours)
        self._janitor_cleanup()

        # Create sandbox
        self.sandbox_path = Path(tempfile.mkdtemp(
            prefix="RPH_G16_",
            dir=str(self.base_dir)
        ))

        # CRITICAL: Create scratch directory
        scratch_path = self.sandbox_path / "scratch"
        scratch_path.mkdir(parents=True, exist_ok=True)

        logger.info(f"Sandbox created: {self.sandbox_path}")
        return self.sandbox_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exit sandbox context and cleanup.

        Args:
            exc_type: Exception type if an error occurred
            exc_val: Exception value if an error occurred
            exc_tb: Exception traceback if an error occurred
        """
        if self.sandbox_path and self.sandbox_path.exists():
            try:
                shutil.rmtree(self.sandbox_path, ignore_errors=True)
                logger.info(f"Sandbox cleaned: {self.sandbox_path}")
            except Exception as e:
                logger.warning(f"Failed to clean sandbox: {e}")

        return False  # Don't suppress exceptions

    @staticmethod
    def _get_free_disk_space_gb(path: Path) -> float:
        """
        Get free disk space in GB.

        Args:
            path: Path to check disk space

        Returns:
            Free space in GB
        """
        stat = shutil.disk_usage(path)
        return stat.free / (1024 ** 3)

    def _janitor_cleanup(self):
        """
        Clean up old RPH_G16_* sandbox directories.

        Deletes folders older than 6 hours.
        """
        import time

        if not self.base_dir.exists():
            return

        max_age_hours = 6
        current_time = time.time()
        max_age_seconds = max_age_hours * 3600

        cleaned_count = 0

        for item in self.base_dir.iterdir():
            if item.is_dir() and item.name.startswith('RPH_G16_'):
                try:
                    mtime = item.stat().st_mtime
                    age_seconds = current_time - mtime

                    if age_seconds > max_age_seconds:
                        logger.info(f"Janitor: removing old sandbox {item.name} "
                                  f"({age_seconds / 3600:.1f} hours old)")
                        shutil.rmtree(item, ignore_errors=True)
                        cleaned_count += 1
                except Exception as e:
                    logger.warning(f"Janitor: failed to clean {item.name}: {e}")

        if cleaned_count > 0:
            logger.info(f"Janitor: cleaned {cleaned_count} old sandbox(es)")


class GaussianRunner:
    """
    Executes Gaussian calculations in a controlled environment.

    Ensures proper environment setup, safe subprocess execution,
    and comprehensive error handling.
    """

    @staticmethod
    def run(sandbox_path: Path, input_content: str,
            timeout: int = 86400) -> str:
        """
        Run Gaussian calculation.

        Args:
            sandbox_path: Path to sandbox directory (must exist)
            input_content: Gaussian input file content
            timeout: Timeout in seconds (default: 24 hours)

        Returns:
            Content of Gaussian log file

        Raises:
            GaussianError: If execution fails
        """
        # Write input file
        input_file = sandbox_path / "input.gjf"
        with open(input_file, 'w', newline='\n') as f:
            f.write(input_content)

        logger.info(f"Input file written: {input_file}")

        # Prepare environment
        env = os.environ.copy()

        # Force override GAUSS_SCRDIR to sandbox scratch
        env["GAUSS_SCRDIR"] = str(sandbox_path / "scratch")

        logger.debug(f"GAUSS_SCRDIR set to: {env['GAUSS_SCRDIR']}")

        # Execute Gaussian
        try:
            logger.info(f"Starting Gaussian execution in: {sandbox_path}")

            result = subprocess.run(
                ["g16", "input.gjf"],
                cwd=str(sandbox_path),
                env=env,
                timeout=timeout,
                capture_output=True,
                text=True,
                shell=False  # Explicitly disable shell=True
            )

            # Always read log file content
            log_file = sandbox_path / "input.log"
            if log_file.exists():
                log_content = log_file.read_text()
            else:
                log_content = ""
                logger.warning(f"Log file not found: {log_file}")

            # Check return code
            if result.returncode != 0:
                error_msg = (
                    f"Gaussian failed with return code {result.returncode}\n"
                    f"STDOUT: {result.stdout}\n"
                    f"STDERR: {result.stderr}\n"
                    f"Log preview (first 500 chars): {log_content[:500]}"
                )
                logger.error(error_msg)
                raise GaussianError(log_content)

            logger.info("Gaussian execution completed successfully")
            return log_content

        except subprocess.TimeoutExpired:
            error_msg = f"Gaussian execution timed out after {timeout} seconds"
            logger.error(error_msg)

            # Try to read log file even on timeout
            log_file = sandbox_path / "input.log"
            if log_file.exists():
                log_content = log_file.read_text()
                raise GaussianError(log_content)
            else:
                raise GaussianError(error_msg)


class ResultHarvester:
    """
    Collects and validates Gaussian calculation results.

    Ensures only essential files are copied, avoiding disk overflow.
    """

    @staticmethod
    def harvest(sandbox_path: Path,
                destination_dir: Path) -> Dict[str, Optional[Path]]:
        """
        Harvest Gaussian results from sandbox.

        Args:
            sandbox_path: Path to sandbox directory
            destination_dir: Path to destination directory

        Returns:
            Dict with keys 'log' (Path) and 'chk' (Path or None)

        Raises:
            RuntimeError: If log file is missing
        """
        destination_dir = Path(destination_dir)
        destination_dir.mkdir(parents=True, exist_ok=True)

        results: Dict[str, Optional[Path]] = {}

        # Copy log file (REQUIRED)
        log_source = sandbox_path / "input.log"
        if not log_source.exists():
            raise RuntimeError(f"Log file missing: {log_source}")

        log_dest = destination_dir / "final_output.log"
        shutil.copy(log_source, log_dest)
        results['log'] = log_dest
        logger.info(f"Log file harvested: {log_dest}")

        # Copy checkpoint file (OPTIONAL)
        chk_source = sandbox_path / "checkpoint.chk"
        if chk_source.exists():
            chk_dest = destination_dir / "final_output.chk"
            shutil.copy(chk_source, chk_dest)
            results['chk'] = chk_dest
            logger.info(f"Checkpoint file harvested: {chk_dest}")
        else:
            logger.warning("Checkpoint file not found (this may be expected)")

        # Forbidden: Do NOT copy .rwf, .int, or scratch/
        forbidden_patterns = ['*.rwf', '*.int', 'scratch']
        for pattern in forbidden_patterns:
            for item in sandbox_path.glob(pattern):
                logger.debug(f"Skipped forbidden file/directory: {item}")

        return results


def run_gaussian_optimization(route, atoms, charge, mult, output_dir, config):
    """
    Main entry point for Gaussian execution.
    
    Args:
        route (str): Gaussian route line.
        atoms (list): List of atom dicts.
        charge (int): System charge.
        mult (int): Spin multiplicity.
        output_dir (Path or str): Directory to save final files.
        config (dict): Configuration dict containing 'mem', 'nproc', 'timeout'.
    """
    # 1. Unpack Config (Adapter Layer)
    mem = config.get('mem', '48GB')
    nproc = config.get('nproc', 16)
    timeout = config.get('timeout', 86400)
    
    try:
        # 2. Input Generation
        gjf_content = InputFactory.create(route, atoms, charge, mult, mem, nproc)
        
        # 3. Sandbox Execution (Safety: Check for 100GB free space)
        with LinuxSandbox(min_free_gb=100.0) as sandbox:
            runner = GaussianRunner()
            # Pass unpacked args to runner
            log_content = runner.run(sandbox, gjf_content, timeout=timeout)
            
            # 4. Result Harvesting
            # FIX: Ensure we use the 'output_dir' argument passed in
            results = ResultHarvester.harvest(sandbox, output_dir)
            
            # 5. Parsing
            final_atoms = LogParser.extract_final_geometry(log_content)
            
            return {
                "success": True,
                "final_geometry": final_atoms,
                "log_path": results['log'],
                "chk_path": results['chk']
            }

    except Exception as e:
        return {"success": False, "error": str(e)}
