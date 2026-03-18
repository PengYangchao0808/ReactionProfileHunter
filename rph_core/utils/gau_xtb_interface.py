"""
Gau_XTB Interface - Gaussian calling XTB for TS optimization
============================================================

Implements the Gau_XTB workflow:
- Gaussian with external='./xtb.sh' calls XTB for gradients
- Supports TS optimization, IRC, and regular optimization
- Auto-detects task type from filename

Author: ReactionProfileHunter Team
Date: 2026-03-16
"""

import logging
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import read_xyz, write_xyz
from rph_core.utils.geometry_tools import LogParser
from rph_core.utils.data_types import QCResult
from rph_core.utils.qc_interface import is_path_toxic


class GauXTBInterface(LoggerMixin):
    """
    Gaussian + XTB hybrid TS optimizer.
    
    Uses Gaussian's external keyword to call XTB for gradient computation
    during TS optimization. This provides more accurate TS structures than
    XTB alone while being faster than pure DFT.
    
    Usage:
        >>> gau_xtb = GauXTBInterface(config)
        >>> result = gau_xtb.optimize_ts(
        ...     xyz_file=Path("ts_guess.xyz"),
        ...     output_dir=Path("./gau_xtb_ts_opt"),
        ...     task_name="TS-001"
        ... )
    """
    
    # Source files required for Gau_XTB
    REQUIRED_FILES = ["xtb.sh", "genxyz", "extderi"]
    
    def __init__(
        self,
        config: Dict[str, Any],
        nproc: int = 1,
        gau_xtb_dir: Optional[Path] = None
    ):
        """
        Initialize Gau_XTB interface.
        
        Args:
            config: Configuration dictionary
            nproc: Number of processors for Gaussian
            gau_xtb_dir: Path to Gau_XTB scripts directory
                        Default: scripts/Gau_XTB/
        """
        self.config = config
        self.nproc = nproc
        
        # Find Gau_XTB scripts directory
        if gau_xtb_dir is None:
            # Default to scripts/Gau_XTB in project root
            import rph_core
            project_root = Path(rph_core.__file__).parent.parent
            gau_xtb_dir = project_root / "scripts" / "Gau_XTB"
        
        self.gau_xtb_dir = Path(gau_xtb_dir)
        
        # Verify required files exist
        self._verify_dependencies()
        
        # Gaussian settings
        gaussian_cfg = config.get("executables", {}).get("gaussian", {})
        self.gaussian_bin = gaussian_cfg.get("path", "g16")
        
    def _verify_dependencies(self) -> None:
        """Verify Gau_XTB scripts exist."""
        missing = []
        for fname in self.REQUIRED_FILES:
            fpath = self.gau_xtb_dir / fname
            if not fpath.exists():
                missing.append(fname)
        
        if missing:
            raise RuntimeError(
                f"Gau_XTB required files not found in {self.gau_xtb_dir}: {missing}"
            )
        self.logger.info(f"Gau_XTB interface initialized with scripts from {self.gau_xtb_dir}")
    
    def _copy_xtb_scripts(self, work_dir: Path) -> None:
        """Copy required XTB scripts to work directory."""
        for fname in self.REQUIRED_FILES:
            src = self.gau_xtb_dir / fname
            dst = work_dir / fname
            shutil.copy2(src, dst)
            # Make executable
            dst.chmod(dst.stat().st_mode | 0o111)
        self.logger.debug(f"Copied Gau_XTB scripts to {work_dir}")
    
    def _generate_gjf(
        self,
        xyz_file: Path,
        output_dir: Path,
        task_name: str,
        route_options: Optional[str] = None
    ) -> Path:
        """
        Generate Gaussian input file with appropriate route.
        
        Auto-detects task type from task_name:
        - TS-xxx → TS optimization
        - TS-xxx-IRC → IRC calculation  
        - INT-xxx, S-xxx, P-xxx → Regular optimization
        
        Args:
            xyz_file: Input XYZ file
            output_dir: Output directory
            task_name: Task name (determines route)
            route_options: Additional route options
            
        Returns:
            Path to generated .gjf file
        """
        # Read XYZ file
        coords, symbols = read_xyz(xyz_file)
        
        # Build route based on task name
        route = self._select_route(task_name, route_options)
        
        # Generate Gaussian input
        gjf_name = f"{task_name}.gjf"
        gjf_path = output_dir / gjf_name
        
        # Build charge/multiplicity (default: neutral singlet)
        charge = 0
        mult = 1
        
        with open(gjf_path, 'w') as f:
            f.write("%nprocshared=1\n")
            f.write(f"{route}\n\n")
            f.write(f"{task_name}\n\n")
            f.write(f"{charge} {mult}\n")
            
            for symbol, (x, y, z) in zip(symbols, coords):
                f.write(f"{symbol:2s} {x:15.8f} {y:15.8f} {z:15.8f}\n")
            
            f.write("\n")
        
        self.logger.debug(f"Generated Gaussian input: {gjf_path}")
        return gjf_path
    
    def _select_route(
        self,
        task_name: str,
        extra_options: Optional[str] = None
    ) -> str:
        """
        Select Gaussian route based on task name.
        
        Auto-detection rules:
        - TS-xxx → opt=(calcall,TS,noeigen,nomicro) external='./xtb.sh'
        - TS-xxx-IRC → IRC(calcfc,stepsize=5) external='./xtb.sh'
        - INT/S/P-xxx → opt=(calcall,noeigen,nomicro) external='./xtb.sh'
        """
        name_upper = task_name.upper()
        
        # IRC task
        if "IRC" in name_upper:
            route = "#p IRC(calcfc,stepsize=5)"
        # TS task
        elif name_upper.startswith("TS"):
            route = "#p opt=(calcall,TS,noeigen,nomicro)"
        # Regular optimization (INT, S, P)
        else:
            route = "#p opt=(calcall,noeigen,nomicro)"
        
        # Add external keyword
        route += " external='./xtb.sh' nosymm"
        
        # Add extra options
        if extra_options:
            route += f" {extra_options}"
        
        return route
    
    def optimize_ts(
        self,
        xyz_file: Path,
        output_dir: Path,
        task_name: str = "TS-001",
        max_cycles: int = 50,
        extra_options: Optional[str] = None,
        charge: int = 0,
        spin: int = 1
    ) -> QCResult:
        """
        Run Gaussian+XTB TS optimization.
        
        Args:
            xyz_file: Input TS guess XYZ file
            output_dir: Output directory
            task_name: Task name (determines optimization type)
            max_cycles: Maximum optimization cycles
            extra_options: Additional Gaussian route options
            charge: Molecular charge
            spin: Spin multiplicity
            
        Returns:
            QCResult with optimized geometry
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        resolved_output_dir = output_dir.resolve()

        execution_dir = output_dir
        cleanup_execution_dir = False
        if is_path_toxic(resolved_output_dir):
            execution_dir = Path(tempfile.mkdtemp(prefix="rph_gau_xtb_", dir="/tmp"))
            cleanup_execution_dir = True
            self.logger.warning(
                "Gau_XTB output path is toxic; using sandbox execution dir: %s",
                execution_dir,
            )
        
        # Copy input XYZ to work dir
        work_xyz = execution_dir / "input.xyz"
        shutil.copy2(xyz_file, work_xyz)
        
        # Copy Gau_XTB scripts
        self._copy_xtb_scripts(execution_dir)
        
        # Generate Gaussian input
        gjf_path = self._generate_gjf(
            xyz_file=work_xyz,
            output_dir=execution_dir,
            task_name=task_name,
            route_options=extra_options
        )
        
        # Update charge/multiplicity if non-default
        if charge != 0 or spin != 1:
            content = gjf_path.read_text()
            lines = content.split('\n')
            # Find charge/multiplicity line (after route block)
            for i, line in enumerate(lines):
                if line.strip() and not line.startswith('%') and not line.startswith('#'):
                    # This should be charge/mult line
                    try:
                        lines[i] = f"{charge} {spin}"
                        gjf_path.write_text('\n'.join(lines))
                    except Exception:
                        pass
                    break
        
        output_log = output_dir / "input.log"
        
        self.logger.info(f"Running Gau_XTB TS optimization: {task_name}")
        
        try:
            # Use wrapper script for Gaussian execution
            # Write input file as input.gjf for wrapper compatibility
            input_file = execution_dir / "input.gjf"
            if gjf_path.name != "input.gjf":
                shutil.copy2(gjf_path, input_file)

            output_log = execution_dir / "input.log"
            
            # Get wrapper path from config
            gaussian_cfg = self.config.get("executables", {}).get("gaussian", {})
            wrapper_path = gaussian_cfg.get("wrapper_path", "scripts/run_g16_worker.sh")
            
            # Resolve wrapper path relative to project root
            import rph_core
            project_root = Path(rph_core.__file__).parent.parent
            wrapper_full_path = project_root / wrapper_path
            
            if not wrapper_full_path.exists():
                raise RuntimeError(f"Gaussian wrapper not found: {wrapper_full_path}")
            
            # Set environment
            env = os.environ.copy()
            xtb_path = self.config.get("executables", {}).get("xtb", {}).get("path")
            if xtb_path:
                env["XTB_PATH"] = str(xtb_path)
                xtb_bin_dir = str(Path(xtb_path).parent)
                env["PATH"] = f"{xtb_bin_dir}:{env.get('PATH', '')}"
            
            # Run Gaussian via wrapper
            self.logger.info(f"Executing Gaussian via wrapper: {wrapper_full_path}")
            
            result = subprocess.run(
                [str(wrapper_full_path), "input.gjf", "input.log"],
                cwd=str(execution_dir),
                env=env,
                capture_output=True,
                text=True,
                timeout=86400
            )

            if execution_dir != output_dir:
                for artifact in ["input.log", "input.gjf", "xtbout", "gradient", "hessian"]:
                    src = execution_dir / artifact
                    if src.exists():
                        shutil.copy2(src, output_dir / artifact)
                output_log = output_dir / "input.log"
            
            if output_log.exists() and "Normal termination" in output_log.read_text():
                parser = LogParser()
                opt_coords, opt_symbols, error_msg = parser.extract_last_converged_coords(
                    output_log, engine_type="gaussian"
                )
                
                if opt_coords is not None and opt_symbols is not None:
                    opt_xyz = output_dir / "ts_final.xyz"
                    write_xyz(opt_xyz, opt_coords, opt_symbols)
                    
                    energy = self._parse_energy(output_log)
                    imag_freq = self._check_imaginary_frequency(output_log)
                    
                    self.logger.info(f"Gau_XTB TS optimization successful: {opt_xyz}")
                    
                    return QCResult(
                        success=True,
                        energy=energy,
                        coordinates=opt_coords,
                        converged=True,
                        output_file=opt_xyz
                    )
            
            # Failed
            error_msg = "Gau_XTB optimization failed"
            if output_log.exists():
                error_msg += f": {output_log.read_text()[-500:]}"
            self.logger.error(error_msg)
            return QCResult(success=False, error_message=error_msg)

        except Exception as e:
            self.logger.error(f"Gau_XTB error: {e}")
            return QCResult(success=False, error_message=str(e))
        finally:
            if cleanup_execution_dir and execution_dir.exists():
                shutil.rmtree(execution_dir, ignore_errors=True)
    
    def _parse_energy(self, log_file: Path) -> Optional[float]:
        """Parse final energy from Gaussian log."""
        if not log_file.exists():
            return None
        
        content = log_file.read_text()
        
        patterns = [
            r"Final energy\s*=\s*([-+]?\d+\.\d+)",
            r"SCF Done:\s*E\([^)]+\)\s*=\s*([-+]?\d+\.\d+)",
            r"(?m)^\s*Energy=\s*([-+]?\d+\.\d+)",
        ]
        
        for pattern in patterns:
            matches = re.findall(pattern, content)
            if matches:
                return float(matches[-1])
        
        return None
    
    def _check_imaginary_frequency(self, log_file: Path) -> Optional[float]:
        """
        Check for imaginary frequency in Gaussian frequency output.
        
        Returns:
            Imaginary frequency in cm-1, or None if not found
        """
        if not log_file.exists():
            return None
        
        content = log_file.read_text()
        
        # Look for negative frequencies (imaginary)
        # Gaussian outputs "Frequencies" and "Freq" lines
        # Imaginary frequencies are printed as negative values
        
        # Pattern for negative frequency
        pattern = r"(?:Frequencies|Freq)\s+--\s+([-+]?\d+\.\d+)"
        matches = re.findall(pattern, content)
        
        for freq_str in matches:
            freq = float(freq_str)
            if freq < 0:
                return freq  # Return negative (imaginary) frequency
        
        return None


class GauXTBOptimizer:
    """High-level wrapper for running Gau_XTB TS optimization with retries."""
    
    def __init__(self, config: Dict[str, Any], nproc: int = 1):
        self.config = config
        self.nproc = nproc
        self.interface = None
        self.logger = logging.getLogger(__name__)
    
    def _get_gau_xtb_config(self) -> Dict[str, Any]:
        """Get Gau_XTB configuration from config."""
        step2_cfg = self.config.get("step2", {}) or {}
        return step2_cfg.get("gau_xtb", {}) or {}
    
    def optimize(
        self,
        ts_guess_xyz: Path,
        output_dir: Path,
        task_name: str,
        forming_bonds: Optional[List[Tuple[int, int]]] = None,
        max_attempts: int = 3
    ) -> Tuple[Path, str, Optional[float]]:
        """
        Run Gau_XTB TS optimization with validation and retries.
        
        Args:
            ts_guess_xyz: Input TS guess XYZ
            output_dir: Output directory
            task_name: Task name for Gaussian
            forming_bonds: Forming bonds (for metadata)
            max_attempts: Maximum optimization attempts
            
        Returns:
            Tuple of (optimized_xyz, confidence, imaginary_freq)
                - optimized_xyz: Path to optimized TS structure
                - confidence: "high", "medium", or "low"
                - imaginary_freq: Imaginary frequency if found
        """
        gau_xtb_cfg = self._get_gau_xtb_config()
        
        if not self.interface:
            self.interface = GauXTBInterface(
                config=self.config,
                nproc=self.nproc,
                gau_xtb_dir=gau_xtb_cfg.get("scripts_dir")
            )
        
        max_cycles = gau_xtb_cfg.get("max_cycles", 50)
        
        best_result = None
        best_confidence = "low"
        best_imag_freq = None
        
        for attempt in range(1, max_attempts + 1):
            attempt_dir = output_dir / f"attempt_{attempt}"
            attempt_name = f"{task_name}-attempt{attempt}"
            
            self.logger.info(
                f"[Gau_XTB] Attempt {attempt}/{max_attempts}: {attempt_name}"
            )
            
            result = self.interface.optimize_ts(
                xyz_file=ts_guess_xyz,
                output_dir=attempt_dir,
                task_name=attempt_name,
                max_cycles=max_cycles
            )
            
            if result.success:
                if result.output_file is None:
                    self.logger.warning(f"[Gau_XTB] Attempt {attempt} has no output file")
                    continue
                
                # Get the log file path - it's in the same directory as the XYZ output
                log_file = result.output_file.parent / "input.log"
                imag_freq = self.interface._check_imaginary_frequency(log_file)
                
                if imag_freq is not None and imag_freq < 0:
                    confidence = "high"
                    self.logger.info(
                        f"[Gau_XTB] Success! Found imaginary frequency: "
                        f"{imag_freq:.1f} cm-1 (confidence: {confidence})"
                    )
                    return (
                        result.output_file,
                        confidence,
                        imag_freq
                    )
                elif imag_freq is None:
                    confidence = "medium"
                    self.logger.info(
                        f"[Gau_XTB] Optimization converged, no frequency data. "
                        f"(confidence: {confidence})"
                    )
                    if best_confidence != "high":
                        best_result = result.output_file
                        best_confidence = confidence
                        best_imag_freq = None
                else:
                    self.logger.warning(
                        f"[Gau_XTB] Multiple or no imaginary frequencies found"
                    )
            else:
                self.logger.warning(
                    f"[Gau_XTB] Attempt {attempt} failed: {result.error_message}"
                )
        
        if best_result:
            return best_result, best_confidence, best_imag_freq
        
        raise RuntimeError(
            f"Gau_XTB TS optimization failed after {max_attempts} attempts"
        )
