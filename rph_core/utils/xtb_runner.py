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
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.data_types import QCResult, ScanResult, PathSearchResult


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

    FALLBACK_PATHS = ["/opt/xtb/bin/xtb", "/usr/local/bin/xtb"]
    _ENERGY_VALUE_PATTERN = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
    _ENERGY_LINE_PATTERNS = (
        re.compile(r"\*?\s*total\s+energy(?!\s+gain)\s*:?\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", re.IGNORECASE),
        re.compile(r"\benergy\s*[=:]\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", re.IGNORECASE),
        re.compile(r"\bscf\s+done\s*:\s*e\([^)]*\)\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", re.IGNORECASE),
    )

    def __init__(self, config: Dict[str, Any], work_dir: Optional[Path] = None):
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
        self._xtb_path = self.xtb_path

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
        executables_cfg: Dict[str, Any] = self.config.get('executables', {}) or {}
        xtb_cfg: Dict[str, Any] = executables_cfg.get('xtb', {}) or {}
        config_path = xtb_cfg.get('path')
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
        configured_fallbacks = xtb_cfg.get('fallback_paths', [])
        if not isinstance(configured_fallbacks, list):
            configured_fallbacks = []
        fallback_paths = configured_fallbacks + [p for p in self.FALLBACK_PATHS if p not in configured_fallbacks]
        for fallback in fallback_paths:
            if Path(fallback).exists() and shutil.which(fallback):
                self.logger.info(f"XTB found at fallback path: {fallback}")
                return fallback

        # 4. All attempts failed
        error_msg = (
            "XTB executable not found. Searched locations:\n"
            f"  - Config path: {config_path}\n"
            f"  - System PATH\n"
            f"  - Fallback paths: {fallback_paths}\n"
            "Please ensure XTB is installed and accessible."
        )
        self.logger.error(error_msg)
        raise RuntimeError(error_msg)

    def _get_timeout_seconds(self) -> Optional[float]:
        timeout_cfg = self.config.get('optimization_control', {}).get('timeout', {})
        if not timeout_cfg.get('enabled', False):
            return None
        timeout_value = timeout_cfg.get('default_seconds')
        if timeout_value is None:
            return None
        try:
            return float(timeout_value)
        except (TypeError, ValueError):
            self.logger.warning(f"Invalid timeout value for XTB optimize: {timeout_value}")
            return None

    def _write_constraint_input(
        self,
        constraints: Dict[str, float],
        force_constant: float = 1.0,
        keepaway_constraints: Optional[Dict[str, float]] = None
    ) -> Path:
        """
        Write XTB constraint file from dictionary.

        Converts bond length constraints to XTB $constrain block format.
        Supports both forming-bond constraints and keep-away constraints.

        Args:
            constraints: Dictionary with format {"atom1_idx atom2_idx": distance}
                        (e.g., {"1 5": 2.1, "10 15": 2.3})
            force_constant: Harmonic constraint force constant (default: 1.0)
            keepaway_constraints: Optional additional constraints to prevent
                                  unintended bond formation (same format)

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
            f.write(f"  force constant={force_constant:.3f}\n")

            # Write primary constraints (forming bonds)
            for atoms, distance in constraints.items():
                atom1, atom2 = map(int, atoms.split())
                f.write(f"  distance: {atom1}, {atom2}, {distance:.3f}\n")

            # Write keep-away constraints (to prevent unintended bonds)
            if keepaway_constraints:
                for atoms, distance in keepaway_constraints.items():
                    atom1, atom2 = map(int, atoms.split())
                    f.write(f"  distance: {atom1}, {atom2}, {distance:.3f}\n")

            f.write("$end\n")

        total_constraints = len(constraints) + (len(keepaway_constraints) if keepaway_constraints else 0)
        self.logger.debug(f"Wrote constraint file: {constraint_file} "
                         f"({len(constraints)} forming, "
                         f"{len(keepaway_constraints) if keepaway_constraints else 0} keep-away, "
                         f"force_constant={force_constant:.3f})")
        return constraint_file

    def _write_fix_input(self, frozen_indices: List[int]) -> Path:
        fix_file = self.work_dir / "xcontrol.inp"
        with fix_file.open('w') as f:
            f.write("$fix\n")
            atoms_str = ",".join(str(idx + 1) for idx in frozen_indices)
            f.write(f"   atoms: {atoms_str}\n")
            f.write("$end\n")
        self.logger.debug(f"Wrote $fix file: {fix_file} (atoms: {len(frozen_indices)})")
        return fix_file

    def optimize(
        self,
        structure: Path,
        constraints: Optional[Dict[str, float]] = None,
        frozen_indices: Optional[List[int]] = None,
        solvent: Optional[str] = None,
        charge: int = 0,
        uhf: int = 0,
        constraint_force_constant: float = 1.0,
        keepaway_constraints: Optional[Dict[str, float]] = None
    ) -> QCResult:
        """
        Run XTB geometry optimization.

        Args:
            structure: Input XYZ file path
            constraints: Optional bond length constraints
                         Format: {"atom1_idx atom2_idx": distance}
            frozen_indices: Optional list of atom indices (0-based) to freeze
            solvent: Solvent name for ALPB model (e.g., "acetone")
            charge: Molecular charge (default: 0)
            uhf: Number of unpaired electrons (default: 0)
            constraint_force_constant: Force constant for constraints (default: 1.0)
            keepaway_constraints: Optional keep-away constraints to prevent
                                  unintended bond formation

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
        if constraints or keepaway_constraints:
            constraint_file = self._write_constraint_input(
                constraints or {},
                force_constant=constraint_force_constant,
                keepaway_constraints=keepaway_constraints
            )
            cmd.extend(["--input", constraint_file.name])
            self.logger.debug(f"Applying constraints: {constraints}, keepaway: {keepaway_constraints}")

        if frozen_indices:
            fix_file = self._write_fix_input(frozen_indices)
            cmd.extend(["--input", fix_file.name])
            self.logger.debug(f"Freezing {len(frozen_indices)} atoms")

        self.logger.info(f"Running XTB optimization: {' '.join(cmd)}")

        try:
            # Execute XTB
            result = subprocess.run(
                cmd,
                cwd=self.work_dir,
                capture_output=True,
                text=True,
                check=True,
                timeout=self._get_timeout_seconds(),
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

    def _run_command(self, cmd: List[str], log_file: Optional[Path] = None) -> subprocess.CompletedProcess[str]:
        try:
            result = subprocess.run(
                cmd,
                cwd=self.work_dir,
                capture_output=True,
                text=True,
                check=True
            )

            if log_file:
                log_file.write_text(result.stdout + ("\n" + result.stderr if result.stderr else ""))

            return result

        except subprocess.CalledProcessError as e:
            if log_file:
                stdout = e.stdout or ""
                stderr = e.stderr or ""
                log_file.write_text(stdout + ("\n" + stderr if stderr else ""))

            self.logger.error(
                "XTB command failed with return code %s: %s",
                e.returncode,
                " ".join(e.cmd) if isinstance(e.cmd, list) else str(e.cmd)
            )
            raise

    def run_scan(
        self,
        input_xyz: Path,
        constraints: Dict[str, float],
        scan_range: Tuple[float, float],
        scan_steps: int,
        scan_mode: str = "concerted",
        scan_force_constant: float = 1.0,
        solvent: Optional[str] = None,
        gfn_level: int = 2,
        etemp: Optional[float] = None,
        charge: int = 0,
        uhf: int = 0
    ) -> ScanResult:
        input_xyz = Path(input_xyz)
        if not input_xyz.exists():
            self.logger.error("XTB scan input structure not found: %s", input_xyz)
            return ScanResult(success=False)

        local_structure = self.work_dir / input_xyz.name
        if local_structure.resolve() != input_xyz.resolve():
            local_structure.write_text(input_xyz.read_text())

        scan_input = self._write_scan_input(
            constraints=constraints,
            scan_range=scan_range,
            scan_steps=scan_steps,
            scan_mode=scan_mode,
            scan_force_constant=scan_force_constant,
        )

        cmd = [self._xtb_path, local_structure.name, "--opt", "--input", scan_input.name]

        nproc = self.config.get('resources', {}).get('nproc', 1)
        cmd.extend(["-P", str(nproc)])
        cmd.extend(["--chrg", str(charge)])
        if uhf > 0:
            cmd.extend(["--uhf", str(uhf)])

        if solvent:
            cmd.extend(["--gfn", str(int(gfn_level)), "--alpb", str(solvent)])
        elif gfn_level != 2:
            cmd.extend(["--gfn", str(int(gfn_level))])

        if etemp is not None:
            cmd.extend(["--etemp", str(float(etemp))])

        scan_log = self.work_dir / "xtb_scan.log"
        self.logger.info("Running XTB scan: %s", " ".join(cmd))

        try:
            self._run_command(cmd, log_file=scan_log)
            parsed = self._parse_scan_log(scan_log)

            if isinstance(parsed, ScanResult):
                if parsed.scan_log is None:
                    parsed.scan_log = scan_log
                return parsed

            if isinstance(parsed, dict):
                parsed.setdefault("scan_log", scan_log)
                return ScanResult(**parsed)

            self.logger.error("_parse_scan_log returned unsupported type: %s", type(parsed).__name__)
            return ScanResult(success=False, scan_log=scan_log)

        except subprocess.CalledProcessError:
            return ScanResult(success=False, scan_log=scan_log)

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

    def _write_scan_input(
        self,
        constraints: Dict[str, float],
        scan_range: Tuple[float, float],
        scan_steps: int,
        scan_mode: str = "concerted",
        scan_force_constant: float = 1.0,
    ) -> Path:
        if not constraints:
            raise ValueError("scan constraints must not be empty")

        if scan_steps <= 0:
            raise ValueError(f"scan_steps must be positive, got {scan_steps}")

        if scan_force_constant <= 0:
            raise ValueError(f"scan_force_constant must be positive, got {scan_force_constant}")

        mode = scan_mode.strip().lower()
        if mode not in {"concerted", "sequential"}:
            raise ValueError(f"unsupported scan mode: {scan_mode}")

        start, end = scan_range
        maxcycle = self.config.get("optimization_control", {}).get("max_cycles", 100)
        scan_file = self.work_dir / "scan.inp"

        parsed_constraints: List[Tuple[int, int, float]] = []
        for atoms, distance in constraints.items():
            atom1_raw, atom2_raw = map(int, atoms.split())
            atom1 = atom1_raw + 1
            atom2 = atom2_raw + 1
            parsed_constraints.append((atom1, atom2, distance))

        with scan_file.open("w") as f:
            f.write("$constrain\n")
            f.write(f"  force constant={scan_force_constant:.3f}\n")
            for atom1, atom2, distance in parsed_constraints:
                f.write(f"  distance: {atom1}, {atom2}, {distance:.3f}\n")

            f.write("$scan\n")
            f.write(f"  mode={mode}\n")
            for idx in range(1, len(parsed_constraints) + 1):
                f.write(f"  {idx}: {start:.3f}, {end:.3f}, {scan_steps}\n")

            f.write("$opt\n")
            f.write(f"  maxcycle={maxcycle}\n")
            f.write("$end\n")

        self.logger.debug(
            "Wrote scan input file: %s (constraints=%d, mode=%s, steps=%d)",
            scan_file,
            len(parsed_constraints),
            mode,
            scan_steps,
        )
        return scan_file

    def _parse_scan_log(self, scan_log: Path) -> ScanResult:
        primary_log = Path(scan_log)
        xtbscan_log = primary_log.parent / "xtbscan.log"
        candidate_logs: List[Path] = [xtbscan_log, primary_log]

        parsed_candidates: List[Tuple[Path, List[float], List[List[str]]]] = []
        checked_logs: List[str] = []

        for candidate in candidate_logs:
            if candidate in [item[0] for item in parsed_candidates]:
                continue

            checked_logs.append(str(candidate))
            if not candidate.exists() or not candidate.is_file() or candidate.stat().st_size == 0:
                continue

            parsed_payload = self._parse_scan_frames(candidate)
            if parsed_payload is None:
                continue

            energies, frames = parsed_payload
            parsed_candidates.append((candidate, energies, frames))

        if not parsed_candidates:
            self.logger.error("No readable scan log found. Checked: %s", [str(p) for p in candidate_logs])
            return ScanResult(success=False, scan_log=Path(scan_log))

        log_path, energies, frames = max(parsed_candidates, key=lambda item: len(item[1]))

        if len(parsed_candidates) > 1:
            details = ", ".join(f"{item[0].name}:{len(item[1])}" for item in parsed_candidates)
            self.logger.info(
                "Multiple scan logs parsed (%s). Selected %s with %d points",
                details,
                log_path.name,
                len(energies),
            )
        elif checked_logs and log_path.name != Path(checked_logs[0]).name:
            self.logger.info(
                "Selected fallback scan log %s (parsed points=%d)",
                log_path.name,
                len(energies),
            )

        max_energy_index = max(range(len(energies)), key=energies.__getitem__)
        ts_guess_xyz = self.work_dir / "xtb_scan_ts_guess.xyz"
        frame_dir = self.work_dir / "scan_frames"
        frame_dir.mkdir(parents=True, exist_ok=True)
        frame_paths: List[Path] = []

        for idx, frame_lines in enumerate(frames):
            frame_path = frame_dir / f"frame_{idx:03d}.xyz"
            frame_path.write_text("\n".join(frame_lines) + "\n")
            frame_paths.append(frame_path)

        try:
            _ = ts_guess_xyz.write_text("\n".join(frames[max_energy_index]) + "\n")
        except Exception as exc:
            self.logger.error("Parsed scan but failed to write TS guess xyz: %s", exc)
            return ScanResult(
                success=False,
                energies=energies,
                geometries=frame_paths,
                max_energy_index=max_energy_index,
                scan_log=log_path,
            )

        self.logger.info(
            "Parsed %d scan points from %s; max energy index=%d, energy=%.8f",
            len(energies),
            log_path,
            max_energy_index,
            energies[max_energy_index],
        )
        return ScanResult(
            success=True,
            energies=energies,
            geometries=frame_paths,
            max_energy_index=max_energy_index,
            ts_guess_xyz=ts_guess_xyz,
            scan_log=log_path,
        )

    def _parse_scan_frames(self, log_path: Path) -> Optional[Tuple[List[float], List[List[str]]]]:
        try:
            lines = log_path.read_text(errors="replace").splitlines()
        except Exception as exc:
            self.logger.error("Failed to read scan log %s: %s", log_path, exc)
            return None

        energies: List[float] = []
        frames: List[List[str]] = []

        i = 0
        total_lines = len(lines)
        while i < total_lines:
            raw_line = lines[i].strip()
            if not raw_line:
                i += 1
                continue

            try:
                atom_count = int(raw_line)
            except ValueError:
                i += 1
                continue

            frame_end = i + 2 + atom_count
            if atom_count <= 0 or frame_end > total_lines:
                self.logger.warning(
                    "Skip malformed frame at line %d in %s (atom_count=%s)",
                    i + 1,
                    log_path,
                    atom_count,
                )
                i += 1
                continue

            comment_line = lines[i + 1]
            energy = self._extract_energy_from_line(comment_line)
            if energy is None:
                energy = self._find_nearby_energy(lines=lines, frame_start_index=i)

            if energy is None:
                self.logger.warning(
                    "Cannot parse energy for frame at line %d in %s, frame skipped",
                    i + 1,
                    log_path,
                )
                i = frame_end
                continue

            energies.append(energy)
            frame_lines = [lines[i], lines[i + 1]] + lines[i + 2:frame_end]
            frames.append(frame_lines)
            i = frame_end

        if not energies:
            self.logger.warning("No valid scan frame with energy parsed from %s", log_path)
            return None

        return energies, frames

    def _extract_energy_from_line(self, line: str) -> Optional[float]:
        for pattern in self._ENERGY_LINE_PATTERNS:
            match = pattern.search(line)
            if match is None:
                continue
            try:
                return float(match.group(1))
            except (TypeError, ValueError):
                continue

        lowered_line = line.lower()
        skip_markers = ("energy gain", "kcal/mol", "kw/mol")
        if any(marker in lowered_line for marker in skip_markers):
            return None

        preferred_tokens = ("scf", "done", "energy", "e=", "e(")
        if not any(token in lowered_line for token in preferred_tokens):
            return None

        tokens = self._ENERGY_VALUE_PATTERN.findall(line)
        if not tokens:
            return None

        try:
            return float(tokens[-1])
        except ValueError:
            return None

    def _find_nearby_energy(self, lines: List[str], frame_start_index: int, lookback: int = 120) -> Optional[float]:
        start = max(0, frame_start_index - lookback)
        for idx in range(frame_start_index - 1, start - 1, -1):
            energy = self._extract_energy_from_line(lines[idx])
            if energy is not None:
                return energy
        return None

    def run_path(
        self,
        start_xyz: Path,
        end_xyz: Path,
        nrun: int = 1,
        npoint: int = 25,
        anopt: int = 10,
        kpush: float = 0.003,
        kpull: float = -0.015,
        ppull: float = 0.05,
        alp: float = 1.2,
        gfn_level: int = 2,
        solvent: Optional[str] = None,
        charge: int = 0,
        uhf: int = 0,
        etemp: Optional[float] = None,
    ) -> PathSearchResult:
        start_xyz = Path(start_xyz)
        end_xyz = Path(end_xyz)
        if not start_xyz.exists():
            return PathSearchResult(success=False, error_message=f"start_xyz not found: {start_xyz}")
        if not end_xyz.exists():
            return PathSearchResult(success=False, error_message=f"end_xyz not found: {end_xyz}")

        local_start = self.work_dir / start_xyz.name
        local_end = self.work_dir / end_xyz.name
        if local_start.resolve() != start_xyz.resolve():
            local_start.write_text(start_xyz.read_text())
        if local_end.resolve() != end_xyz.resolve():
            local_end.write_text(end_xyz.read_text())

        path_input = self.work_dir / "path.inp"
        with path_input.open("w") as f:
            f.write("$path\n")
            f.write(f"   nrun={nrun}\n")
            f.write(f"   npoint={npoint}\n")
            f.write(f"   anopt={anopt}\n")
            f.write(f"   kpush={kpush}\n")
            f.write(f"   kpull={kpull}\n")
            f.write(f"   ppull={ppull}\n")
            f.write(f"   alp={alp}\n")
            f.write("$end\n")

        cmd = [self._xtb_path, local_start.name, "--path", local_end.name, "--input", path_input.name]
        nproc = self.config.get("resources", {}).get("nproc", 1)
        cmd.extend(["-P", str(nproc)])
        cmd.extend(["--chrg", str(charge)])
        if uhf > 0:
            cmd.extend(["--uhf", str(uhf)])
        if solvent:
            cmd.extend(["--gfn", str(gfn_level), "--alpb", str(solvent)])
        elif gfn_level != 2:
            cmd.extend(["--gfn", str(gfn_level)])
        if etemp is not None:
            cmd.extend(["--etemp", str(etemp)])

        path_log = self.work_dir / "xtb_path.log"
        self.logger.info("Running xTB path search: %s", " ".join(cmd))

        try:
            self._run_command(cmd, log_file=path_log)
            parsed = self._parse_path_log(path_log)

            path_xyz_files = sorted(self.work_dir.glob("xtbpath_*.xyz"))

            ts_guess_path = self.work_dir / "xtbpath_ts.xyz"
            if not ts_guess_path.exists():
                return PathSearchResult(
                    success=False,
                    path_xyz_files=path_xyz_files,
                    path_log=path_log,
                    error_message="xtbpath_ts.xyz not found in output",
                )

            return PathSearchResult(
                success=True,
                path_xyz_files=path_xyz_files,
                ts_guess_xyz=ts_guess_path,
                path_log=path_log,
                barrier_forward_kcal=parsed.get("barrier_forward_kcal"),
                barrier_backward_kcal=parsed.get("barrier_backward_kcal"),
                reaction_energy_kcal=parsed.get("reaction_energy_kcal"),
                estimated_ts_point=parsed.get("estimated_ts_point"),
                gradient_norm_at_ts=parsed.get("gradient_norm_at_ts"),
            )

        except Exception as e:
            self.logger.error("xTB path search failed: %s", e)
            return PathSearchResult(success=False, error_message=str(e))

    def _parse_path_log(self, log_path: Path) -> Dict[str, Any]:
        result: Dict[str, Any] = {}
        if not log_path.exists():
            return result

        content = log_path.read_text(errors="replace")

        barrier_fw_match = re.search(r"forward\s+barrier.*?:\s*([-+]?\d*\.?\d+)", content, re.IGNORECASE)
        if barrier_fw_match:
            result["barrier_forward_kcal"] = float(barrier_fw_match.group(1))

        barrier_bw_match = re.search(r"backward\s+barrier.*?:\s*([-+]?\d*\.?\d+)", content, re.IGNORECASE)
        if barrier_bw_match:
            result["barrier_backward_kcal"] = float(barrier_bw_match.group(1))

        rxn_energy_match = re.search(r"reaction\s+energy.*?:\s*([-+]?\d*\.?\d+)", content, re.IGNORECASE)
        if rxn_energy_match:
            result["reaction_energy_kcal"] = float(rxn_energy_match.group(1))

        ts_point_match = re.search(r"estimated\s+TS\s+on\s+file.*?point:\s*(\d+)", content, re.IGNORECASE)
        if ts_point_match:
            result["estimated_ts_point"] = int(ts_point_match.group(1))

        grad_norm_match = re.search(r"norm\(g\)\s+at\s+est\.\s+TS.*?:\s*([-+]?\d*\.?\d+)", content, re.IGNORECASE)
        if grad_norm_match:
            result["gradient_norm_at_ts"] = float(grad_norm_match.group(1))

        return result
