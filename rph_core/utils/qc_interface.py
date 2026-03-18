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
import warnings
from pathlib import Path
from typing import Any, Dict, List, Match, Optional, Tuple
from enum import Enum
import logging

import numpy as np
from rph_core.utils.keyword_translator import KeywordTranslator

logger = logging.getLogger(__name__)


def _format_gaussian_route_block(route_line: str, max_cols: int = 80) -> str:
    """Format a Gaussian route card respecting the 80-column input limit.

    Gaussian input is parsed with a hard 80-column record limit. If we emit a
    long single route line, Gaussian truncates it and can split keywords (e.g.
    'Pop' -> 'P'), causing QPErr syntax errors.

    This formatter wraps only on whitespace and uses Gaussian's continuation
    convention (subsequent route lines begin with a single leading space).
    """

    route_line = re.sub(r"\s+", " ", (route_line or "").strip())
    if not route_line:
        return "#p"

    tokens = route_line.split(" ")
    first_prefix = "#p "
    cont_prefix = " "

    lines: List[str] = []
    current = first_prefix
    current_cap = max_cols

    for token in tokens:
        prefix = first_prefix if not lines and current == first_prefix else cont_prefix
        # current may already be a continuation line; ensure it starts with the right prefix
        if not lines and current == first_prefix:
            prefix = first_prefix
        elif lines and current.startswith(cont_prefix):
            prefix = cont_prefix

        if current == prefix:
            candidate = f"{current}{token}"
        else:
            candidate = f"{current} {token}"

        if len(candidate) <= current_cap:
            current = candidate
            continue

        # Start a new continuation line (never split tokens)
        lines.append(current.rstrip())
        current = f"{cont_prefix}{token}"
        current_cap = max_cols

    lines.append(current.rstrip())
    return "\n".join(lines)


def _normalize_nbo_keylist_block(raw: str) -> str:
    """Normalize $NBO keylist block for Gaussian.

    Accepts either:
    - full block containing $NBO ... $END
    - a bare keylist body (will be wrapped)
    """

    text = (raw or "").strip("\n")
    if not text.strip():
        return ""

    upper = text.upper()
    if "$NBO" in upper:
        # Assume user provides a full block; ensure it ends with newline.
        return text + "\n"

    return "$NBO\n" + text + "\n$END\n"


def _normalize_def2_basis_in_route(route_line: str) -> str:
    def _repl(match: Match[str]) -> str:
        token = match.group(0)
        return token.replace("-", "").replace("_", "")

    return re.sub(r"def2[-_]?[A-Za-z0-9]+", _repl, route_line, flags=re.IGNORECASE)


def _append_modredundant_to_opt(route_line: str) -> str:
    if "modredundant" in route_line.lower():
        return route_line

    def _inject_parenthesized(match: Match[str]) -> str:
        body = match.group(1)
        entries = [item.strip() for item in body.split(",") if item.strip()]
        if any(item.lower() == "modredundant" for item in entries):
            return f"Opt=({body})"
        entries.append("ModRedundant")
        return f"Opt=({','.join(entries)})"

    if re.search(r"\bopt\s*=\s*\([^)]*\)", route_line, flags=re.IGNORECASE):
        return re.sub(
            r"\bopt\s*=\s*\(([^)]*)\)",
            _inject_parenthesized,
            route_line,
            count=1,
            flags=re.IGNORECASE,
        )

    if re.search(r"\bopt\s*=\s*[^\s]+", route_line, flags=re.IGNORECASE):
        return re.sub(
            r"\bopt\s*=\s*([^\s]+)",
            lambda match: f"Opt=({match.group(1)},ModRedundant)",
            route_line,
            count=1,
            flags=re.IGNORECASE,
        )

    if re.search(r"\bopt\b", route_line, flags=re.IGNORECASE):
        return re.sub(
            r"\bopt\b",
            "Opt=ModRedundant",
            route_line,
            count=1,
            flags=re.IGNORECASE,
        )

    return f"Opt=ModRedundant {route_line}".strip()


# =============================================================================
# M3: QC Task Kind Enumeration
# =============================================================================

class TaskKind(Enum):
    """Enumeration of QC task types for unified task framework and resource selection."""
    OPTIMIZATION = "optimization"
    SINGLE_POINT = "single_point"
    FREQUENCY = "frequency"
    TS_OPTIMIZATION = "ts_optimization"
    IRC = "irc"
    NBO = "nbo"
    SCAN = "scan"


# M3-4-2: NBO output file whitelist (prevents picking up old files)
NBO_WHITELIST = {
    '.47': True,
    '.nbo': True,
    '.3': True,
    '.31': True,
    '.41': True,
    '.nbo7': True
}


def _sanitize_jobname(jobname: str) -> str:
    """
    Sanitize jobname to prevent path injection and filesystem issues.

    M3-1: Jobname sanitization for safe file naming.

    Args:
        jobname: Raw job name string

    Returns:
        Sanitized job name safe for use in filenames
    """
    # Remove or replace dangerous characters
    sanitized = re.sub(r'[^\w\-.]', '_', jobname)
    # Collapse multiple underscores
    sanitized = re.sub(r'_+', '_', sanitized)
    # Strip leading/trailing underscores and dots
    sanitized = sanitized.strip('_')
    # Limit length to avoid filesystem limits
    max_len = 128
    if len(sanitized) > max_len:
        sanitized = sanitized[:max_len]
    # Ensure non-empty
    if not sanitized or sanitized.startswith('.'):
        sanitized = 'rph_job'
    return sanitized


def _select_task_resources(
    task_kind: TaskKind,
    config: Dict[str, Any],
    fallback_to_global: bool = True
) -> Dict[str, Any]:
    """
    Select resources and theory parameters for a specific task kind.

    M3-2: Task-specific resource selection with global fallback.

    Args:
        task_kind: Type of QC task (TaskKind enum value)
        config: Full configuration dictionary
        fallback_to_global: If True, fall back to global resources when task-specific ones missing

    Returns:
        Dict with keys: 'mem', 'nproc', 'method', 'basis'
    """
    task_resources = config.get('task_resources', {}).get(task_kind.value, {})
    global_resources = config.get('resources', {})

    mem = task_resources.get('mem') or (global_resources.get('mem') if fallback_to_global else None)
    nproc = task_resources.get('nproc') or (global_resources.get('nproc') if fallback_to_global else None)

    # Determine theory key based on task kind
    theory_key_map = {
        TaskKind.OPTIMIZATION: 'optimization',
        TaskKind.TS_OPTIMIZATION: 'optimization',
        TaskKind.SINGLE_POINT: 'single_point',
        TaskKind.FREQUENCY: 'optimization',  # Freq uses optimization method/basis
        TaskKind.IRC: 'optimization',
        TaskKind.NBO: 'single_point',  # NBO typically runs on SP
    }
    theory_key = theory_key_map.get(task_kind, 'optimization')

    theory_config = config.get('theory', {}).get(theory_key, {})

    method = task_resources.get('method') or (theory_config.get('method') if fallback_to_global else None)
    basis = task_resources.get('basis') or (theory_config.get('basis') if fallback_to_global else None)

    # Provide defaults if nothing found
    if method is None:
        method = 'B3LYP'
    if basis is None:
        basis = 'def2-SVP'

    return {
        'mem': mem or '32GB',
        'nproc': nproc or 16,
        'method': method,
        'basis': basis
    }



def harvest_nbo_files(output_dir: Path, jobname: str, sub_dir: str = "nbo_analysis") -> Dict[str, Path]:
    """
    Collect NBO output files using whitelist-based pattern matching.

    M3-4-2: NBO file collection with whitelist and stale file protection.

    This function robustly searches for NBO files using a whitelist of
    standard extensions. It uses defensive globbing to avoid picking up
    old files from previous runs.

    Args:
        output_dir: Base output directory
        jobname: Job name used to search for NBO files
        sub_dir: Subdirectory within output_dir to search (default: "nbo_analysis")

    Returns:
        Dict mapping file extensions to absolute paths
        Only returns files whose extensions are in NBO_WHITELIST.
    """
    output_dir = Path(output_dir).resolve()
    search_dir = output_dir / sub_dir

    if not search_dir.exists():
        logger.debug(f"NBO directory not found: {search_dir}")
        return {}

    nbo_files: Dict[str, Path] = {}

    for ext in NBO_WHITELIST:
        glob_pattern = f"{jobname}{ext}"
        matching_files = list(search_dir.glob(glob_pattern))

        if matching_files:
            # Sort by modification time, use most recent
            if len(matching_files) > 1:
                matching_files.sort(key=lambda f: f.stat().st_mtime, reverse=True)
                logger.debug(
                    f"Multiple files found for {jobname}{ext}: "
                    f"using most recent: {matching_files[0].name}"
                )

            nbo_files[ext] = matching_files[0].resolve()
            logger.debug(f"Found NBO file: {nbo_files[ext]}")

    if not nbo_files:
        logger.debug(f"No NBO files found in {search_dir} with whitelist")

    return nbo_files


def _is_float(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


class QCTaskResult:
    """Result container for QC task execution."""

    def __init__(
        self,
        success: bool,
        log_file: Optional[Path] = None,
        chk_file: Optional[Path] = None,
        fchk_file: Optional[Path] = None,
        error: Optional[str] = None,
        nbo_files: Optional[Dict[str, Path]] = None,
        **kwargs
    ):
        self.success = success
        self.log_file = log_file
        self.chk_file = chk_file
        self.fchk_file = fchk_file
        self.error = error
        self.nbo_files = nbo_files or {}
        # Store any additional attributes
        for key, value in kwargs.items():
            setattr(self, key, value)


def run_gaussian_task(
    task_kind: TaskKind,
    xyz_file: Path,
    output_dir: Path,
    config: Dict[str, Any],
    jobname: Optional[str] = None,
    extra_route: str = "",
    charge: int = 0,
    spin: int = 1,
    collect_nbo: bool = False
) -> QCTaskResult:
    # M3-1: Sanitize jobname
    raw_jobname = jobname or xyz_file.stem
    safe_jobname = _sanitize_jobname(raw_jobname)

    # M3-2: Select task-specific resources
    resources = _select_task_resources(task_kind, config)
    mem = resources.get('mem', '32GB')
    nproc = resources.get('nproc', 16)
    method = resources.get('method', 'B3LYP')
    basis = resources.get('basis', 'def2-SVP')

    # Build route based on task kind
    route_map = {
        TaskKind.OPTIMIZATION: f"Opt",
        TaskKind.SINGLE_POINT: f"SP",
        TaskKind.FREQUENCY: f"Freq",
        TaskKind.TS_OPTIMIZATION: f"Opt=(TS, CalcFC, NoEigenTest)",
        TaskKind.IRC: f"IRC=(CalcFC)",
        TaskKind.NBO: f"SP Pop=NBO",
    }
    base_route = route_map.get(task_kind, "SP")
    if extra_route:
        base_route = f"{base_route} {extra_route}"
    route = f"{method}/{basis} {base_route}"

    try:
        # Read atoms from input file
        if xyz_file.suffix == '.xyz':
            from rph_core.utils.file_io import read_xyz
            coords, symbols = read_xyz(xyz_file)
            atoms = [{'symbol': s, 'x': x, 'y': y, 'z': z}
                    for s, (x, y, z) in zip(symbols, coords)]
        else:
            # For .gjf files, we'd need to parse differently
            # For now, require .xyz input
            raise ValueError(f"Unsupported input format: {xyz_file.suffix}")

        # Create input file
        gjf_content = InputFactory.create(route, atoms, charge, spin, mem, nproc)

        # Ensure output directory exists
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # M3-1: Sandbox execution with toxic path check
        # LinuxSandbox returns a Path (sandbox dir). It does not expose a `.path` attribute.
        if is_path_toxic(output_dir):
            logger.warning(f"Output directory is toxic, using sandbox: {output_dir}")
            base_dir = "/tmp"
        else:
            base_dir = str(output_dir)

        with LinuxSandbox(min_free_gb=50.0, base_dir=base_dir) as sandbox_dir:
            # Write input file
            input_file = Path(sandbox_dir) / f"{safe_jobname}.gjf"
            input_file.write_text(gjf_content)

            # Execute Gaussian
            runner = GaussianRunner()
            log_content = runner.run(Path(sandbox_dir), gjf_content)

            # Harvest results
            results = ResultHarvester.harvest(Path(sandbox_dir), output_dir)

        log_file = results.get('log')
        chk_file = results.get('chk')

        # Try to generate fchk
        fchk_file = None
        if chk_file and chk_file.exists():
            fchk_file = try_formchk(chk_file)

        # M3-4-2: Collect NBO files if requested
        nbo_files: Dict[str, Path] = {}
        if collect_nbo and task_kind == TaskKind.NBO:
            nbo_files = harvest_nbo_files(output_dir, safe_jobname)

        return QCTaskResult(
            success=True,
            log_file=log_file,
            chk_file=chk_file,
            fchk_file=fchk_file,
            nbo_files=nbo_files
        )

    except Exception as e:
        logger.error(f"Gaussian task {task_kind.value} failed: {e}")
        return QCTaskResult(success=False, error=str(e))


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
        # NOTE: Gaussian has a hard 80-column input record limit (route card in particular).
        # If we emit a single long route line, Gaussian truncates at 80 chars and can turn
        # e.g. "Pop" into "P" (and drop the rest), leading to QPErr syntax errors.
        clean_route = route.encode('ascii', 'ignore').decode('ascii')
        clean_route = clean_route.replace('\r\n', '\n').replace('\r', '\n')

        # Normalize any leading route header fragments; we will re-add '#p' ourselves.
        clean_route = clean_route.strip()
        if clean_route.startswith('#'):
            clean_route = clean_route.lstrip('#').strip()
        if clean_route.lower().startswith('p '):
            clean_route = clean_route[1:].strip()

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
        route_block = _format_gaussian_route_block(clean_route)

        # Optional: NBO keylist when using Pop=(NBORead)
        nbo_block = ""
        if "NBOREAD" in clean_route.upper():
            # NOTE: InputFactory path does not currently take config.
            # If caller wants NBORead, they must include a $NBO...$END block
            # as part of atoms/route handling elsewhere.
            pass

        template = (
            f"%chk=checkpoint.chk\n"
            f"%mem={mem}\n"
            f"%nprocshared={nproc}\n"
            f"{route_block}\n"
            f"\n"
            f"RPH_Auto_Job\n"
            f"\n"
            f"{charge} {mult}\n"
            f"{coords_str}\n"
            f"\n"
            f"{nbo_block}"
            f"\n"
        )

        return template


class LogParser:
    """
    Parses Gaussian log files to extract critical information.

    Focuses on robust extraction of final geometry.
    """

    @staticmethod
    def extract_final_geometry(log_content: str) -> List[Dict[str, Any]]:
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
        lines = log_content.splitlines()

        def _find_last_orientation_start() -> Optional[int]:
            last_idx = None
            orientation_header = re.compile(r"^\s*(Standard|Input)\s+orientation:\s*$", re.IGNORECASE)
            for idx, line in enumerate(lines):
                if orientation_header.match(line):
                    last_idx = idx
            return last_idx

        def _to_float(token: str) -> Optional[float]:
            try:
                return float(token.replace("D", "E").replace("d", "e"))
            except ValueError:
                return None

        start_idx = _find_last_orientation_start()
        if start_idx is None:
            logger.warning("No orientation blocks found in log file")
            return []

        idx = start_idx + 1
        while idx < len(lines) and not re.match(r"\s*-{5,}\s*$", lines[idx]):
            idx += 1
        if idx >= len(lines):
            logger.warning("Orientation block header separator not found in log file")
            return []

        idx += 1
        while idx < len(lines) and not re.match(r"\s*-{5,}\s*$", lines[idx]):
            idx += 1
        if idx >= len(lines):
            logger.warning("Orientation block column header separator not found in log file")
            return []

        idx += 1
        atoms: List[Dict[str, Any]] = []
        while idx < len(lines):
            line = lines[idx]
            if re.match(r"\s*-{5,}\s*$", line):
                break

            parts = line.strip().split()
            if len(parts) >= 6 and parts[0].isdigit():
                atomic_num = None
                try:
                    atomic_num = int(parts[1])
                except ValueError:
                    atomic_num = None

                x = _to_float(parts[3])
                y = _to_float(parts[4])
                z = _to_float(parts[5])
                if x is not None and y is not None and z is not None:
                    symbol = LogParser._atomic_num_to_symbol(atomic_num) if atomic_num is not None else "X"
                    atoms.append({"symbol": symbol, "x": x, "y": y, "z": z})

            idx += 1

        if not atoms:
            logger.warning("Orientation block found but no atomic coordinates parsed")
            return []

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


# Re-export QCResult for backward compatibility
from rph_core.utils.data_types import QCResult, ScanResult, PathSearchResult  # noqa: F401

# Additional imports for interface classes
from rph_core.utils.file_io import read_energy_from_gaussian, read_xyz
from rph_core.utils.resource_utils import resolve_executable_config
from rph_core.utils.xtb_runner import XTBRunner

logger = logging.getLogger(__name__)


class XTBInterface:
    def __init__(
        self,
        gfn_level: int = 2,
        solvent: Optional[str] = None,
        nproc: int = 1,
        config: Optional[Dict[str, Any]] = None
    ):
        self.gfn_level = gfn_level
        self.solvent = solvent
        self.nproc = nproc
        self.config = config or {}

    def optimize(
        self,
        xyz_file: Path,
        output_dir: Path,
        constraints: Optional[Dict[str, float]] = None,
        charge: int = 0,
        spin: int = 1,
        constraint_force_constant: float = 1.0,
        keepaway_constraints: Optional[Dict[str, float]] = None
    ) -> QCResult:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        config = dict(self.config)
        config.setdefault('resources', {})
        config['resources']['nproc'] = self.nproc

        runner = XTBRunner(config=config, work_dir=output_dir)

        normalized_constraints: Optional[Dict[str, float]] = None
        if constraints:
            normalized_constraints = {}
            for atoms, distance in constraints.items():
                atom1_raw, atom2_raw = map(int, str(atoms).split())
                normalized_constraints[f"{atom1_raw + 1} {atom2_raw + 1}"] = float(distance)

        normalized_keepaway: Optional[Dict[str, float]] = None
        if keepaway_constraints:
            normalized_keepaway = {}
            for atoms, distance in keepaway_constraints.items():
                atom1_raw, atom2_raw = map(int, str(atoms).split())
                normalized_keepaway[f"{atom1_raw + 1} {atom2_raw + 1}"] = float(distance)

        result = runner.optimize(
            structure=xyz_file,
            constraints=normalized_constraints,
            solvent=self.solvent,
            charge=charge,
            uhf=max(spin - 1, 0),
            constraint_force_constant=constraint_force_constant,
            keepaway_constraints=normalized_keepaway
        )
        if isinstance(result.coordinates, Path):
            result.output_file = result.coordinates
            try:
                coords, _ = read_xyz(result.output_file)
                result.coordinates = coords
            except Exception as exc:
                result.error_message = f"Unable to read XTB coordinates: {exc}"
                result.converged = False
                result.success = False
        return result

    def scan(
        self,
        xyz_file: Path,
        output_dir: Path,
        constraints: Dict[str, float],
        scan_range: Tuple[float, float],
        scan_steps: int,
        scan_mode: str = "concerted",
        scan_force_constant: float = 1.0,
        charge: int = 0,
        spin: int = 1
    ) -> ScanResult:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        config = dict(self.config)
        config.setdefault('resources', {})
        config['resources']['nproc'] = self.nproc

        use_sandbox = is_path_toxic(output_dir)
        run_dir = output_dir
        if use_sandbox:
            logger.warning("XTB scan output path is toxic, using sandbox: %s", output_dir)
            run_dir = Path(tempfile.mkdtemp(prefix="RPH_XTB_SCAN_", dir="/tmp"))

        try:
            step2_cfg = config.get("step2", {}) or {}
            xtb_settings = step2_cfg.get("xtb_settings", {}) or {}

            resolved_solvent = xtb_settings.get("solvent", self.solvent)
            resolved_gfn_level = int(xtb_settings.get("gfn_level", self.gfn_level))
            raw_etemp = xtb_settings.get("etemp")
            resolved_etemp = float(raw_etemp) if raw_etemp is not None else None

            runner = XTBRunner(config=config, work_dir=run_dir)
            result = runner.run_scan(
                input_xyz=xyz_file,
                constraints=constraints,
                scan_range=scan_range,
                scan_steps=scan_steps,
                scan_mode=scan_mode,
                scan_force_constant=scan_force_constant,
                solvent=resolved_solvent,
                gfn_level=resolved_gfn_level,
                etemp=resolved_etemp,
                charge=charge,
                uhf=max(spin - 1, 0)
            )

            if use_sandbox:
                for item in run_dir.iterdir():
                    target = output_dir / item.name
                    if item.is_dir():
                        if target.exists():
                            shutil.rmtree(target)
                        shutil.copytree(item, target)
                    else:
                        shutil.copy2(item, target)

                def _remap_output_path(path_obj: Path) -> Path:
                    path_obj = Path(path_obj)
                    try:
                        rel = path_obj.relative_to(run_dir)
                        return output_dir / rel
                    except ValueError:
                        return output_dir / path_obj.name

                if result.scan_log:
                    result.scan_log = _remap_output_path(result.scan_log)
                if result.ts_guess_xyz:
                    result.ts_guess_xyz = _remap_output_path(result.ts_guess_xyz)
                if isinstance(result.geometries, Path):
                    result.geometries = _remap_output_path(result.geometries)
                elif isinstance(result.geometries, list):
                    result.geometries = [
                        _remap_output_path(geom) if isinstance(geom, Path) else geom
                        for geom in result.geometries
                    ]

            if result.success:
                logger.info("XTB scan completed successfully: %s", output_dir)
            else:
                logger.warning("XTB scan failed: %s", output_dir)

            return result
        finally:
            if use_sandbox and run_dir.exists():
                shutil.rmtree(run_dir, ignore_errors=True)

    def path(
        self,
        start_xyz: Path,
        end_xyz: Path,
        output_dir: Path,
        nrun: int = 1,
        npoint: int = 25,
        anopt: int = 10,
        kpush: float = 0.003,
        kpull: float = -0.015,
        ppull: float = 0.05,
        alp: float = 1.2,
        charge: int = 0,
        spin: int = 1
    ) -> PathSearchResult:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        config = dict(self.config)
        config.setdefault('resources', {})
        config['resources']['nproc'] = self.nproc

        use_sandbox = is_path_toxic(output_dir)
        run_dir = output_dir
        if use_sandbox:
            logger.warning("XTB path search output path is toxic, using sandbox: %s", output_dir)
            run_dir = Path(tempfile.mkdtemp(prefix="RPH_XTB_PATH_", dir="/tmp"))

        try:
            step2_cfg = config.get("step2", {}) or {}
            xtb_settings = step2_cfg.get("xtb_settings", {}) or {}

            resolved_solvent = xtb_settings.get("solvent", self.solvent)
            resolved_gfn_level = int(xtb_settings.get("gfn_level", self.gfn_level))

            runner = XTBRunner(config=config, work_dir=run_dir)
            result = runner.run_path(
                start_xyz=start_xyz,
                end_xyz=end_xyz,
                nrun=nrun,
                npoint=npoint,
                anopt=anopt,
                kpush=kpush,
                kpull=kpull,
                ppull=ppull,
                alp=alp,
                gfn_level=resolved_gfn_level,
                solvent=resolved_solvent,
                charge=charge,
                uhf=max(spin - 1, 0),
            )

            if use_sandbox:
                for item in run_dir.iterdir():
                    target = output_dir / item.name
                    if item.is_dir():
                        if target.exists():
                            shutil.rmtree(target)
                        shutil.copytree(item, target)
                    else:
                        shutil.copy2(item, target)

                def _remap_output_path(path_obj: Path) -> Path:
                    path_obj = Path(path_obj)
                    try:
                        rel = path_obj.relative_to(run_dir)
                        return output_dir / rel
                    except ValueError:
                        return output_dir / path_obj.name

                if result.path_log:
                    result.path_log = _remap_output_path(result.path_log)
                if result.ts_guess_xyz:
                    result.ts_guess_xyz = _remap_output_path(result.ts_guess_xyz)
                if result.path_xyz_files:
                    result.path_xyz_files = [
                        _remap_output_path(f) if isinstance(f, Path) else f
                        for f in result.path_xyz_files
                    ]

            if result.success:
                logger.info("XTB path search completed successfully: %s", output_dir)
            else:
                logger.warning("XTB path search failed: %s", output_dir)

            return result
        finally:
            if use_sandbox and run_dir.exists():
                shutil.rmtree(run_dir, ignore_errors=True)


class CRESTInterface:
    """CREST conformer search interface with support for two-stage workflow.

    Supports both single-stage (GFN2 only) and two-stage (GFN0→GFN2) conformer search.
    """

    def __init__(
        self,
        gfn_level: int = 2,
        solvent: Optional[str] = None,
        nproc: int = 1,
        config: Optional[Dict[str, Any]] = None,
        additional_flags: Optional[str] = None
    ):
        self.gfn_level = gfn_level
        self.solvent = solvent
        self.nproc = nproc
        self.additional_flags = additional_flags
        self.config = config or {}

        exe_config = resolve_executable_config(
            self.config,
            'crest',
            env_vars=['CREST_PATH', 'CREST_BIN']
        )
        crest_path = exe_config.get('path')
        if crest_path:
            self.crest_cmd = str(crest_path)
        else:
            self.crest_cmd = None

    def run_conformer_search(
        self,
        xyz_file: Path,
        output_dir: Path,
        gfn_override: Optional[int] = None,
        additional_flags: Optional[str] = None
    ) -> Path:
        """Run CREST conformer search on a single structure.

        Args:
            xyz_file: Input XYZ structure file
            output_dir: Output directory for CREST artifacts
            gfn_override: Override GFN level (0, 1, or 2). If None, uses instance gfn_level
            additional_flags: Additional command-line flags for CREST

        Returns:
            Path to best structure or ensemble file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if not self.crest_cmd:
            raise RuntimeError("CREST executable not found")

        gfn_use = gfn_override if gfn_override is not None else self.gfn_level
        cmd = [self.crest_cmd, "input.xyz", f"-gfn{gfn_use}", "-T", str(self.nproc)]

        if self.solvent:
            cmd.extend(["-alpb", self.solvent])

        # Handle additional flags
        flags = additional_flags or self.additional_flags
        if flags:
            cmd.extend(flags.split())

        shutil.copy(xyz_file, output_dir / "input.xyz")
        result = subprocess.run(
            cmd,
            cwd=output_dir,
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"CREST conformer search failed: {result.stderr}")

        best_path = output_dir / "crest_best.xyz"
        if best_path.exists():
            return best_path
        ensemble_path = output_dir / "crest_conformers.xyz"
        if ensemble_path.exists():
            return ensemble_path
        fallback_path = output_dir / xyz_file.name
        shutil.copy(xyz_file, fallback_path)
        return fallback_path

    def run_batch_optimization(
        self,
        ensemble_xyz: Path,
        output_dir: Path,
        gfn_level: int = 2,
        solvent: Optional[str] = None,
        additional_flags: Optional[str] = None
    ) -> Path:
        """Batch optimize multiple structures from an ensemble file.

        This method is used in Stage 2 of the two-stage workflow to refine
        structures that passed through GFN0 screening and ISOSTAT clustering.

        Args:
            ensemble_xyz: Input ensemble XYZ file (from GFN0 stage)
            output_dir: Output directory for optimization artifacts
            gfn_level: GFN level for optimization (0, 1, or 2)
            solvent: Solvent model (ALPB)
            additional_flags: Additional command-line flags for CREST

        Returns:
            Path to optimized ensemble file (crest_ensemble.xyz)
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if not self.crest_cmd:
            raise RuntimeError("CREST executable not found")

        # Use -mdopt mode for batch ensemble optimization
        # This takes an ensemble and optimizes all structures within it
        cmd = [
            self.crest_cmd,
            "-mdopt",
            ensemble_xyz.name,
            f"-gfn{gfn_level}",
            "-T", str(self.nproc)
        ]

        sol = solvent or self.solvent
        if sol:
            cmd.extend(["--alpb", sol])

        flags = additional_flags or self.additional_flags
        if flags:
            cmd.extend(flags.split())

        # Copy ensemble file to output directory
        local_ensemble = output_dir / ensemble_xyz.name
        shutil.copy(ensemble_xyz, local_ensemble)

        result = subprocess.run(
            cmd,
            cwd=output_dir,
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            # Log the error but don't raise - may have partial results
            logger.warning(f"CREST batch optimization returned non-zero: {result.stderr}")

        # Return the optimized ensemble
        optimized_ensemble = output_dir / "crest_ensemble.xyz"
        if optimized_ensemble.exists():
            return optimized_ensemble

        # Fallback: return input if no optimization output
        logger.warning("No optimized ensemble found, returning input ensemble")
        return ensemble_xyz


toxic_chars = {' ', '[', ']', '(', ')', '{', '}'}


def is_path_toxic(path: Path) -> bool:
    """
    Check if path contains characters that break Gaussian/QC calculations.

    This is a HARD CONSTRAINT: paths containing spaces or special characters
    will cause Gaussian to fail on WSL/Windows. The function rejects paths
    with: space, brackets (), [], and braces {}.

    Usage:
        - Used before QC calculations to decide whether to use sandbox
        - Returns True if path is unsafe and requires sandbox execution

    Args:
        path: Path to check

    Returns:
        True if path contains toxic characters, False otherwise
    """
    path_str = str(path)
    return any(char in path_str for char in toxic_chars)


def try_formchk(chk_path: Path) -> Optional[Path]:
    """
    Attempt to convert Gaussian checkpoint (.chk) to formatted checkpoint (.fchk).

    This function is called after every Gaussian job completion to enable
    downstream feature extraction (Step 4) to read wavefunction data.

    Args:
        chk_path: Path to Gaussian .chk file

    Returns:
        Path to .fchk file if conversion succeeds, None otherwise
    """
    if not chk_path.exists():
        logger.warning(f"Checkpoint file not found, skipping formchk: {chk_path}")
        return None

    fchk_path = chk_path.with_suffix(".fchk")

    try:
        result = subprocess.run(
            ["formchk", str(chk_path), str(fchk_path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            logger.warning(
                f"formchk failed for {chk_path.name}: "
                f"returncode={result.returncode}, stderr={result.stderr[:200]}"
            )
            return None

        if not fchk_path.exists():
            logger.warning(f"formchk completed but .fchk not found: {fchk_path}")
            return None

        logger.info(f"Successfully generated formatted checkpoint: {fchk_path}")
        return fchk_path

    except subprocess.TimeoutExpired:
        logger.warning(f"formchk timed out for {chk_path.name}")
        return None
    except Exception as e:
        logger.warning(f"formchk raised exception for {chk_path.name}: {e}")
        return None


class GaussianInterface:
    def __init__(
        self,
        charge: int = 0,
        multiplicity: int = 1,
        nprocshared: int = 16,
        mem: str = "32GB",
        config: Optional[Dict[str, Any]] = None
    ):
        self.charge = charge
        self.multiplicity = multiplicity
        self.nprocshared = nprocshared
        self.mem = mem
        self.config = config or {}

        exe_config = self.config.get('executables', {}).get('gaussian', {})
        self.use_wrapper = exe_config.get('use_wrapper', True)
        if self.use_wrapper:
            wrapper_path = exe_config.get('wrapper_path', './scripts/run_g16_worker.sh')
            wrapper_path = Path(wrapper_path)
            if not wrapper_path.is_absolute():
                wrapper_path = (Path.cwd() / wrapper_path).resolve()
            self.gaussian_cmd = str(wrapper_path)
        else:
            resolved = resolve_executable_config(
                self.config,
                'gaussian',
                env_vars=['GAUSS_PATH', 'GAUSSIAN_PATH']
            )
            self.gaussian_cmd = str(resolved.get('path') or 'g16')

    def write_input_file(
        self,
        xyz_file: Path,
        gjf_file: Path,
        route: str,
        title: str,
        old_checkpoint: Optional[Path] = None
    ) -> Path:
        coords, symbols = read_xyz(xyz_file)
        route_line = route.strip()
        if route_line.startswith('#'):
            route_line = route_line.lstrip('#').strip()
        if route_line.lower().startswith('p '):
            route_line = route_line[1:].strip()
        if old_checkpoint and old_checkpoint.exists() and "Guess=Read" not in route_line:
            route_line = f"{route_line} Guess=Read"
        if "SCRF=" in route_line:
            route_line = route_line.replace("SCRF=(Solvent=", "SCRF=(PCM,Solvent=")

        # Optional NBO keylist injection when using NBORead.
        nbo_keylist = (
            self.config.get("theory", {}).get("optimization", {}).get("nbo_keylist")
            or self.config.get("nbo", {}).get("keylist")
        )
        nbo_block = ""
        if "NBOREAD" in route_line.upper() and nbo_keylist:
            nbo_block = _normalize_nbo_keylist_block(str(nbo_keylist))

        chk_name = f"{gjf_file.stem}.chk"
        lines = [
            f"%chk={chk_name}\n",
            f"%mem={self.mem}\n",
            f"%nprocshared={self.nprocshared}\n"
        ]
        if old_checkpoint and old_checkpoint.exists():
            lines.append(f"%oldchk={old_checkpoint.name}\n")
        lines.append(f"{_format_gaussian_route_block(route_line)}\n\n")
        lines.append(f"{title}\n\n")
        lines.append(f"{self.charge} {self.multiplicity}\n")
        for symbol, coord in zip(symbols, coords):
            lines.append(f"{symbol:<2} {coord[0]:14.8f} {coord[1]:14.8f} {coord[2]:14.8f}\n")
        lines.append("\n")
        if nbo_block:
            lines.append(nbo_block)
            lines.append("\n")
        gjf_file.parent.mkdir(parents=True, exist_ok=True)
        gjf_file.write_text("".join(lines))
        return gjf_file

    def _build_default_constrained_route(self) -> str:
        theory_opt = self.config.get("theory", {}).get("optimization", {})
        method = str(theory_opt.get("method", "B3LYP")).strip() or "B3LYP"
        basis_raw = str(theory_opt.get("basis", "def2-SVP")).strip() or "def2-SVP"
        basis = KeywordTranslator.to_gaussian_basis(basis_raw)

        route_parts: List[str] = [f"{method}/{basis}"]

        dispersion = KeywordTranslator.to_gaussian_dispersion(str(theory_opt.get("dispersion", "")))
        if dispersion:
            route_parts.append(dispersion)

        solvent = KeywordTranslator.to_gaussian_solvent(str(theory_opt.get("solvent", "")))
        if solvent:
            route_parts.append(solvent.strip())

        route_parts.append("Opt=CalcFC")
        return " ".join(route_parts)

    def optimize(
        self,
        xyz_file: Path,
        output_dir: Path,
        route: str,
        old_checkpoint: Optional[Path] = None,
        timeout: Optional[int] = None
    ) -> QCResult:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        gjf_file = output_dir / f"{xyz_file.stem}.gjf"
        log_file = output_dir / f"{xyz_file.stem}.log"

        local_old_chk = None
        if old_checkpoint and Path(old_checkpoint).exists():
            local_old_chk_path = output_dir / "previous.chk"
            try:
                shutil.copy2(old_checkpoint, local_old_chk_path)
                local_old_chk = Path("previous.chk")
            except Exception as e:
                import warnings
                warnings.warn(f"Cannot localize checkpoint: {e}")
                local_old_chk = Path(old_checkpoint).absolute()

        self.write_input_file(
            xyz_file=xyz_file,
            gjf_file=gjf_file,
            route=route,
            title=xyz_file.stem,
            old_checkpoint=local_old_chk
        )

        cmd = [self.gaussian_cmd, gjf_file.name, log_file.name] if self.use_wrapper else [self.gaussian_cmd, gjf_file.name]

        result = subprocess.run(
            cmd,
            cwd=output_dir,
            capture_output=True,
            text=True,
            timeout=timeout
        )

        if result.returncode != 0 or not log_file.exists():
            error_snippet = "Unknown error"
            if log_file.exists():
                try:
                    with open(log_file, 'r', errors='replace') as f:
                        lines = f.readlines()
                        important_lines = [l.strip() for l in lines[-50:] if "Error" in l or "termination" in l or "galloc" in l or "allocation" in l]
                        if important_lines:
                            error_snippet = " | ".join(important_lines)
                        else:
                            error_snippet = "".join(lines[-10:])
                except Exception:
                    error_snippet = f"Cannot read log file: {log_file}"
            elif result.stderr:
                error_snippet = result.stderr.strip()

            return QCResult(
                success=False,
                converged=False,
                error_message=error_snippet
            )

        log_content = log_file.read_text()
        converged = "Normal termination" in log_content
        energy = None
        try:
            energy = read_energy_from_gaussian(log_file)
        except Exception:
            energy = None

        atoms = LogParser.extract_final_geometry(log_content)
        coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms]) if atoms else np.array([])
        frequencies = self._parse_frequencies(log_content)

        chk_file = gjf_file.with_suffix(".chk")
        fchk_file = try_formchk(chk_file) if chk_file.exists() else None

        error_message = None
        if not converged:
            error_lines = [line.strip() for line in log_content.splitlines() if "Error termination" in line]
            if error_lines:
                error_message = error_lines[-1]
            else:
                error_message = "Gaussian did not terminate normally"

        return QCResult(
            success=converged,
            energy=energy,
            coordinates=coords,
            converged=converged,
            frequencies=frequencies,
            output_file=log_file,
            log_file=log_file,
            chk_file=chk_file if chk_file.exists() else None,
            fchk_file=fchk_file,
            qm_output_file=log_file,
            error_message=error_message
        )

    @staticmethod
    def _parse_frequencies(log_content: str) -> Optional[np.ndarray]:
        freqs: List[float] = []
        for line in log_content.splitlines():
            stripped = line.strip()
            if stripped.startswith("Frequencies --"):
                parts = stripped.split()[2:]
                for part in parts:
                    try:
                        freqs.append(float(part))
                    except ValueError:
                        continue
        if not freqs:
            return None
        return np.array(freqs)

    def constrained_optimize(
        self,
        xyz_file: Path,
        output_dir: Path,
        frozen_indices: Optional[List[int]] = None,
        distance_constraints: Optional[List[Tuple[int, int]]] = None,
        charge: int = 0,
        spin: int = 1,
        route: Optional[str] = None,
        timeout: Optional[int] = None
    ) -> QCResult:
        from rph_core.utils.file_io import read_xyz

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        coords, symbols = read_xyz(xyz_file)
        atom_count = len(symbols)

        invalid_reasons: List[str] = []
        normalized_frozen_indices: List[int] = []
        for raw_idx in frozen_indices or []:
            idx = int(raw_idx)
            if idx < 0 or idx >= atom_count:
                invalid_reasons.append(f"frozen index out of range: {idx} (atom_count={atom_count})")
                continue
            normalized_frozen_indices.append(idx)

        normalized_distance_constraints: List[Tuple[int, int]] = []
        for raw_bond in distance_constraints or []:
            if len(raw_bond) != 2:
                invalid_reasons.append(f"invalid bond constraint shape: {raw_bond}")
                continue
            i, j = int(raw_bond[0]), int(raw_bond[1])
            if i == j:
                invalid_reasons.append(f"invalid bond constraint with identical atoms: ({i}, {j})")
                continue
            if i < 0 or i >= atom_count or j < 0 or j >= atom_count:
                invalid_reasons.append(f"bond constraint out of range: ({i}, {j}) (atom_count={atom_count})")
                continue
            normalized_distance_constraints.append((i, j))

        if invalid_reasons:
            return QCResult(
                success=False,
                converged=False,
                error_message="; ".join(invalid_reasons),
            )

        if route is None:
            route = self._build_default_constrained_route()

        mod_redundant_lines: List[str] = []
        for idx in normalized_frozen_indices:
            mod_redundant_lines.append(f"X  {idx + 1}  F")
        for bond in normalized_distance_constraints:
            i, j = int(bond[0]), int(bond[1])
            mod_redundant_lines.append(f"B  {i + 1}  {j + 1}  F")

        if not mod_redundant_lines:
            return QCResult(success=False, converged=False, error_message="No constrained terms provided")

        mod_redundant_block = "\n".join(mod_redundant_lines)

        route_line = route.strip().lstrip('#').strip()
        if route_line.lower().startswith('p '):
            route_line = route_line[1:].strip()
        route_line = re.sub(r"\s+", " ", route_line).strip()
        route_line = _normalize_def2_basis_in_route(route_line)

        # Ensure ModRedundant is in route when constraints are present
        if mod_redundant_lines and 'modredundant' not in route_line.lower():
            route_line = _append_modredundant_to_opt(route_line)

        use_sandbox = is_path_toxic(output_dir)
        run_dir = output_dir
        if use_sandbox:
            logger.warning("Constrained Gaussian output path is toxic, using sandbox: %s", output_dir)
            run_dir = Path(tempfile.mkdtemp(prefix="RPH_G16_CONSTR_", dir="/tmp"))

        chk_name = f"{xyz_file.stem}_constrained.chk"
        log_name = f"{xyz_file.stem}_constrained.log"
        gjf_name = f"{xyz_file.stem}_constrained.gjf"
        gjf_file = run_dir / gjf_name
        log_file = run_dir / log_name

        lines = [
            f"%chk={chk_name}\n",
            f"%mem={self.mem}\n",
            f"%nprocshared={self.nprocshared}\n",
            f"{_format_gaussian_route_block(route_line)}\n\n",
            f"Constrained Optimization\n\n",
            f"{charge} {spin}\n"
        ]
        for symbol, coord in zip(symbols, coords):
            lines.append(f"{symbol:<2} {coord[0]:14.8f} {coord[1]:14.8f} {coord[2]:14.8f}\n")
        lines.append("\n")
        lines.append(mod_redundant_block)
        lines.append("\n\n")

        gjf_file.write_text("".join(lines))

        try:
            cmd = [self.gaussian_cmd, gjf_file.name, log_file.name] if self.use_wrapper else [self.gaussian_cmd, gjf_file.name]
            result = subprocess.run(
                cmd,
                cwd=str(run_dir),
                capture_output=True,
                text=True,
                timeout=timeout
            )

            if use_sandbox:
                for item in run_dir.iterdir():
                    target = output_dir / item.name
                    if item.is_dir():
                        if target.exists():
                            shutil.rmtree(target)
                        shutil.copytree(item, target)
                    else:
                        shutil.copy2(item, target)

            final_log_file = (output_dir / log_name) if use_sandbox else log_file

            if not final_log_file.exists():
                return QCResult(
                    success=False,
                    converged=False,
                    error_message=(
                        f"Gaussian log not created (returncode={result.returncode})"
                    ),
                )

            log_content = final_log_file.read_text(errors="replace")
            has_normal_termination = "Normal termination" in log_content
            has_error_termination = "Error termination" in log_content

            if result.returncode != 0:
                error_lines = [
                    line.strip()
                    for line in log_content.splitlines()
                    if (
                        "Error termination" in line
                        or "QPErr" in line
                        or "Lnk1e" in line
                        or "Cannot open" in line
                        or "No such file" in line
                        or "Error imposing constraints" in line
                    )
                ]
                detail = error_lines[0] if error_lines else (result.stderr or "").strip()
                detail = detail.splitlines()[-1] if detail else ""
                message = f"Gaussian execution failed (returncode={result.returncode})"
                if detail:
                    message = f"{message}: {detail}"
                return QCResult(
                    success=False,
                    converged=False,
                    coordinates=None,
                    output_file=final_log_file,
                    error_message=message,
                )

            atoms = LogParser.extract_final_geometry(log_content)
            result_coords = np.array([[a['x'], a['y'], a['z']] for a in atoms]) if atoms else None
            error_lines = [
                line.strip() for line in log_content.splitlines() if "Error termination" in line
            ]

            converged = has_normal_termination

            if not converged:
                if has_error_termination and error_lines:
                    message = f"Gaussian did not converge: {error_lines[-1]}"
                elif has_error_termination:
                    message = "Gaussian did not converge: error termination detected"
                else:
                    message = "Gaussian did not converge: no normal termination marker"
                return QCResult(
                    success=False,
                    converged=False,
                    coordinates=result_coords,
                    output_file=final_log_file,
                    error_message=message,
                )

            if result_coords is None:
                has_orientation_header = (
                    "Standard orientation" in log_content or "Input orientation" in log_content
                )
                if has_orientation_header:
                    message = (
                        "Gaussian geometry parse failed: orientation sections found but no coordinates parsed"
                    )
                else:
                    message = "Gaussian geometry parse failed: no orientation sections found in log"
                return QCResult(
                    success=False,
                    converged=True,
                    coordinates=None,
                    output_file=final_log_file,
                    error_message=message,
                )

            return QCResult(
                success=True,
                converged=True,
                coordinates=result_coords,
                output_file=final_log_file,
                error_message=None,
            )

        except subprocess.TimeoutExpired:
            return QCResult(success=False, converged=False, error_message="Gaussian timed out")
        except Exception as e:
            return QCResult(success=False, converged=False, error_message=str(e))
        finally:
            if use_sandbox and run_dir.exists():
                shutil.rmtree(run_dir, ignore_errors=True)


class QCInterfaceFactory:
    @staticmethod
    def create_interface(engine_type: str, **kwargs):
        """
        Factory method to create QC interface instances by engine type.

        Args:
            engine_type: Type of QC engine ('gaussian', 'orca', 'xtb')
            **kwargs: Additional arguments passed to the interface constructor

        Returns:
            Configured interface instance

        Raises:
            ValueError: If engine type is unsupported
        """
        engine = engine_type.lower().strip()
        if engine == 'gaussian':
            return GaussianInterface(
                charge=kwargs.get('charge', 0),
                multiplicity=kwargs.get('multiplicity', 1),
                nprocshared=kwargs.get('nprocshared', 16),
                mem=kwargs.get('mem', '32GB'),
                config=kwargs.get('config')
            )
        if engine == 'orca':
            from rph_core.utils.orca_interface import ORCAInterface
            return ORCAInterface(
                method=kwargs.get('method', 'M062X'),
                basis=kwargs.get('basis', 'def2-TZVPP'),
                aux_basis=kwargs.get('aux_basis', 'def2/J'),
                nprocs=kwargs.get('nprocshared', kwargs.get('nprocs', 16)),
                maxcore=kwargs.get('maxcore'),
                solvent=kwargs.get('solvent', 'acetone'),
                config=kwargs.get('config')
            )
        if engine == 'xtb':
            return XTBInterface(
                gfn_level=kwargs.get('gfn_level', 2),
                solvent=kwargs.get('solvent'),
                nproc=kwargs.get('nproc', 1),
                config=kwargs.get('config')
            )
        raise ValueError(f"Unsupported QC engine: {engine_type}")
