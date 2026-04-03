"""
Multiwfn Runner
===============

Non-interactive batch runner for Multiwfn analysis.
Supports Fukui functions, dual descriptor, QTAIM, and NCI analysis.

P2 Feature: Multiwfn Tier-1 implementation with caching and fail-open.

Author: RPH Team
Date: 2026-02-02
"""

import hashlib
import json
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = logging.getLogger(__name__)

# Multiwfn warning codes
W_MULTIWFN_DISABLED = "W_MULTIWFN_DISABLED"
W_MULTIWFN_FAILED = "W_MULTIWFN_FAILED"
W_MULTIWFN_CACHE_READ_FAILED = "W_MULTIWFN_CACHE_READ_FAILED"
W_MULTIWFN_TIMEOUT = "W_MULTIWFN_TIMEOUT"
W_MULTIWFN_INVALID_OUTPUT = "W_MULTIWFN_INVALID_OUTPUT"


@dataclass
class MultiwfnResult:
    """Result from Multiwfn analysis."""
    success: bool
    fukui_fplus_a: Optional[float] = None
    fukui_fplus_b: Optional[float] = None
    fukui_fminus_a: Optional[float] = None
    fukui_fminus_b: Optional[float] = None
    fukui_f0_a: Optional[float] = None
    fukui_f0_b: Optional[float] = None
    dual_descriptor_a: Optional[float] = None
    dual_descriptor_b: Optional[float] = None
    rho_bcp_forming1: Optional[float] = None
    laplacian_bcp_forming1: Optional[float] = None
    nci_stats_ts: Optional[Dict[str, float]] = None
    nci_stats_reactant: Optional[Dict[str, float]] = None
    warnings: List[str] = field(default_factory=list)
    cache_hit: bool = False
    error_message: Optional[str] = None


class MultiwfnRunner:
    """Non-interactive Multiwfn batch runner with caching and fail-open."""

    def __init__(
        self,
        multiwfn_path: str = "Multiwfn",
        cache_dir: Path = Path(".cache/step4_multiwfn"),
        timeout_sec: int = 120,
        enabled_modules: Optional[Dict[str, bool]] = None
    ):
        """Initialize Multiwfn runner.

        Args:
            multiwfn_path: Path to Multiwfn executable
            cache_dir: Directory for cache files
            timeout_sec: Timeout for each Multiwfn execution
            enabled_modules: Dict of module_name -> enabled
        """
        self.multiwfn_path = multiwfn_path
        self.cache_dir = Path(cache_dir)
        self.timeout_sec = timeout_sec
        self.enabled_modules = enabled_modules or {
            "fukui": True,
            "dual_descriptor": True,
            "qtaim_bcp": False,
            "nci": False
        }

        # Ensure cache directory exists
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _generate_cache_key(
        self,
        input_file: Path,
        module: str,
        atom_indices: Optional[Tuple[int, int]] = None
    ) -> str:
        """Generate cache key for Multiwfn result.

        Args:
            input_file: Input file path (fchk, wfn, or wfx)
            module: Multiwfn module being run
            atom_indices: Optional atom indices for atom-specific analysis

        Returns:
            Cache key string (sha1[:16])
        """
        if not input_file.exists():
            return ""

        try:
            with open(input_file, "rb") as f:
                content = f.read()
            file_hash = hashlib.sha256(content).hexdigest()[:16]

            cache_data = {
                "file_hash": file_hash,
                "module": module,
                "atom_indices": atom_indices,
                "multiwfn_version": self._get_multiwfn_version()
            }

            key_str = json.dumps(cache_data, sort_keys=True, default=str)
            return hashlib.sha1(key_str.encode('utf-8')).hexdigest()[:16]
        except (OSError, TypeError, ValueError):
            return ""

    def _get_multiwfn_version(self) -> str:
        """Get Multiwfn version string."""
        try:
            result = subprocess.run(
                [self.multiwfn_path, "-version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                # Parse version from output
                for line in result.stdout.split('\n'):
                    if 'version' in line.lower():
                        return line.strip()
            return "unknown"
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            return "unknown"

    def _get_cache_path(self, cache_key: str) -> Path:
        """Get path for cached result."""
        return self.cache_dir / f"{cache_key}.json"

    def _run_multiwfn_script(
        self,
        input_file: Path,
        script: str
    ) -> Tuple[bool, str]:
        """Run Multiwfn with a script.

        Args:
            input_file: Input file for Multiwfn
            script: Multiwfn commands (newline-separated)

        Returns:
            Tuple of (success, output)
        """
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                script_file = tmpdir / "multiwfn_script.txt"
                output_file = tmpdir / "multiwfn_output.txt"

                # Write script
                script_file.write_text(script)

                # Run Multiwfn
                env = os.environ.copy()
                env['MULTIWFN_TMPDIR'] = str(tmpdir)

                result = subprocess.run(
                    [
                        self.multiwfn_path,
                        str(input_file),
                        "-script", str(script_file)
                    ],
                    capture_output=True,
                    text=True,
                    timeout=self.timeout_sec,
                    env=env
                )

                # Multiwfn outputs to stdout in script mode
                output = result.stdout + result.stderr

                # Check for success indicators
                success = result.returncode == 0 and "Normal termination" in output

                return success, output

        except subprocess.TimeoutExpired:
            return False, "Timeout"
        except Exception as e:
            return False, str(e)

    def _parse_fukui_from_output(self, output: str) -> Dict[str, float]:
        """Parse Fukui function values from Multiwfn output.

        Args:
            output: Multiwfn stdout/stderr

        Returns:
            Dict of fukui values
        """
        result = {}

        # Patterns for Fukui function values
        # Multiwfn typically outputs f+ (electrophilic), f- (nucleophilic), f0 (frontier)
        patterns = {
            'fukui_fplus': [r'f\+[:\s]+([-\d.]+)', r'Electrophilic.*?([-\d.]+)'],
            'fukui_fminus': [r'f\-[:\s]+([-\d.]+)', r'Nucleophilic.*?([-\d.]+)'],
            'fukui_f0': [r'f0[:\s]+([-\d.]+)', r'Frontier.*?([-\d.]+)']
        }

        for key, pattern_list in patterns.items():
            for pattern in pattern_list:
                import re
                match = re.search(pattern, output, re.IGNORECASE)
                if match:
                    try:
                        result[key] = float(match.group(1))
                        break
                    except ValueError:
                        continue

        return result

    def _parse_dual_descriptor_from_output(self, output: str) -> Dict[str, float]:
        """Parse dual descriptor values from Multiwfn output.

        Args:
            output: Multiwfn stdout/stderr

        Returns:
            Dict of dual descriptor values
        """
        result = {}

        import re
        # Dual descriptor = f+ - f-
        patterns = [
            r'Dual descriptor[:\s]+([-\d.]+)',
            r'Δf\(r\)[:\s]+([-\d.]+)'
        ]

        for pattern in patterns:
            match = re.search(pattern, output, re.IGNORECASE)
            if match:
                try:
                    result['dual_descriptor'] = float(match.group(1))
                    break
                except ValueError:
                    continue

        return result

    def _parse_qtaim_bcp_from_output(self, output: str) -> Dict[str, float]:
        """Parse QTAIM BCP values from Multiwfn output.

        Args:
            output: Multiwfn stdout/stderr

        Returns:
            Dict of QTAIM BCP values (rho, laplacian, etc.)
        """
        result = {}

        import re
        # Look for critical point analysis section
        patterns = {
            'rho_bcp': r'rho\(BCP\)[:\s]+([-\d.]+)',
            'laplacian_bcp': r'∇²rho\(BCP\)[:\s]+([-\d.]+)',
            'energy_bcp': r'Energy density\(BCP\)[:\s]+([-\d.]+)'
        }

        for key, pattern in patterns.items():
            match = re.search(pattern, output, re.IGNORECASE)
            if match:
                try:
                    result[key] = float(match.group(1))
                except ValueError:
                    continue

        return result

    def _generate_fukui_script(self, atom_a: int, atom_b: int) -> str:
        """Generate Multiwfn script for Fukui function analysis.

        Args:
            atom_a: First atom index (1-based)
            atom_b: Second atom index (1-based)

        Returns:
            Multiwfn script commands
        """
        # Multiwfn command sequence for Fukui analysis
        # 7 = FMO analysis
        # Then select Fukui function submenu
        return f"""7
1
{atom_a}
{atom_b}
0
0
"""

    def _generate_dual_descriptor_script(self, atom_a: int, atom_b: int) -> str:
        """Generate Multiwfn script for dual descriptor analysis.

        Args:
            atom_a: First atom index (1-based)
            atom_b: Second atom index (1-based)

        Returns:
            Multiwfn script commands
        """
        return f"""7
2
{atom_a}
{atom_b}
0
0
"""

    def _generate_qtaim_script(self, atom_a: int, atom_b: int) -> str:
        """Generate Multiwfn script for QTAIM BCP analysis.

        Args:
            atom_a: First atom index (1-based)
            atom_b: Second atom index (1-based)

        Returns:
            Multiwfn script commands
        """
        return f"""5
1
{atom_a}
{atom_b}
0
0
"""

    def run_fukui_analysis(
        self,
        input_file: Path,
        atom_a: int,
        atom_b: int
    ) -> Tuple[Optional[Dict[str, float]], List[str]]:
        """Run Fukui function analysis for two atoms.

        Args:
            input_file: Input file (.fchk, .wfn, or .wfx)
            atom_a: First atom index (1-based)
            atom_b: Second atom index (1-based)

        Returns:
            Tuple of (fukui_values dict, warnings list)
        """
        warnings = []
        cache_key = self._generate_cache_key(input_file, "fukui", (atom_a, atom_b))
        cache_path = self._get_cache_path(cache_key)

        # Check cache
        if cache_path.exists():
            try:
                with open(cache_path, 'r') as f:
                    cache_data = json.load(f)
                return cache_data.get('result', {}), [W_MULTIWFN_CACHE_READ_FAILED]
            except (json.JSONDecodeError, OSError, UnicodeDecodeError):
                pass

        # Generate script and run
        script = self._generate_fukui_script(atom_a, atom_b)
        success, output = self._run_multiwfn_script(input_file, script)

        if not success:
            warnings.append(W_MULTIWFN_FAILED)
            if output == "Timeout":
                warnings.append(W_MULTIWFN_TIMEOUT)
            return None, warnings

        # Parse results
        result = self._parse_fukui_from_output(output)
        result.update(self._parse_dual_descriptor_from_output(output))

        if not result:
            warnings.append(W_MULTIWFN_INVALID_OUTPUT)

        # Cache result
        if cache_key and result:
            try:
                cache_data = {
                    'result': result,
                    'cache_key': cache_key,
                    'input_file': str(input_file)
                }
                with open(cache_path, 'w') as f:
                    json.dump(cache_data, f)
            except (OSError, TypeError):
                pass

        return result if result else None, warnings

    def run_dual_descriptor_analysis(
        self,
        input_file: Path,
        atom_a: int,
        atom_b: int
    ) -> Tuple[Optional[float], List[str]]:
        """Run dual descriptor analysis for two atoms.

        Args:
            input_file: Input file (.fchk, .wfn, or .wfx)
            atom_a: First atom index (1-based)
            atom_b: Second atom index (1-based)

        Returns:
            Tuple of (dual_descriptor value, warnings list)
        """
        warnings = []
        cache_key = self._generate_cache_key(input_file, "dual_descriptor", (atom_a, atom_b))
        cache_path = self._get_cache_path(cache_key)

        # Check cache
        if cache_path.exists():
            try:
                with open(cache_path, 'r') as f:
                    cache_data = json.load(f)
                    dd_val = cache_data.get('result', {}).get('dual_descriptor')
                    return dd_val, [W_MULTIWFN_CACHE_READ_FAILED]
            except (json.JSONDecodeError, OSError, UnicodeDecodeError):
                pass

        # Generate script and run
        script = self._generate_dual_descriptor_script(atom_a, atom_b)
        success, output = self._run_multiwfn_script(input_file, script)

        if not success:
            warnings.append(W_MULTIWFN_FAILED)
            if output == "Timeout":
                warnings.append(W_MULTIWFN_TIMEOUT)
            return None, warnings

        # Parse result
        result = self._parse_dual_descriptor_from_output(output)
        dd_val = result.get('dual_descriptor')

        if dd_val is None:
            warnings.append(W_MULTIWFN_INVALID_OUTPUT)

        # Cache result
        if cache_key and dd_val is not None:
            try:
                cache_data = {
                    'result': result,
                    'cache_key': cache_key,
                    'input_file': str(input_file)
                }
                with open(cache_path, 'w') as f:
                    json.dump(cache_data, f)
            except (OSError, TypeError):
                pass

        return dd_val, warnings

    def run_qtaim_bcp_analysis(
        self,
        input_file: Path,
        atom_a: int,
        atom_b: int
    ) -> Tuple[Optional[Dict[str, float]], List[str]]:
        """Run QTAIM BCP analysis for bond between two atoms.

        Args:
            input_file: Input file (.fchk, .wfn, or .wfx)
            atom_a: First atom index (1-based)
            atom_b: Second atom index (1-based)

        Returns:
            Tuple of (qtaim_values dict, warnings list)
        """
        warnings = []
        cache_key = self._generate_cache_key(input_file, "qtaim_bcp", (atom_a, atom_b))
        cache_path = self._get_cache_path(cache_key)

        # Check cache
        if cache_path.exists():
            try:
                with open(cache_path, 'r') as f:
                    cache_data = json.load(f)
                return cache_data.get('result', {}), [W_MULTIWFN_CACHE_READ_FAILED]
            except (json.JSONDecodeError, OSError, UnicodeDecodeError):
                pass

        # Generate script and run
        script = self._generate_qtaim_script(atom_a, atom_b)
        success, output = self._run_multiwfn_script(input_file, script)

        if not success:
            warnings.append(W_MULTIWFN_FAILED)
            if output == "Timeout":
                warnings.append(W_MULTIWFN_TIMEOUT)
            return None, warnings

        # Parse results
        result = self._parse_qtaim_bcp_from_output(output)

        if not result:
            warnings.append(W_MULTIWFN_INVALID_OUTPUT)

        # Cache result
        if cache_key and result:
            try:
                cache_data = {
                    'result': result,
                    'cache_key': cache_key,
                    'input_file': str(input_file)
                }
                with open(cache_path, 'w') as f:
                    json.dump(cache_data, f)
            except (OSError, TypeError):
                pass

        return result if result else None, warnings

    def run_complete_analysis(
        self,
        input_file: Path,
        atom_a: int,
        atom_b: int
    ) -> MultiwfnResult:
        """Run complete Multiwfn analysis based on enabled modules.

        Args:
            input_file: Input file (.fchk, .wfn, or .wfx)
            atom_a: First atom index (1-based)
            atom_b: Second atom index (1-based)

        Returns:
            MultiwfnResult with all analyzed values
        """
        result = MultiwfnResult(success=False, warnings=[])

        # Check if Multiwfn is enabled
        if not any(self.enabled_modules.values()):
            result.warnings.append(W_MULTIWFN_DISABLED)
            return result

        # Run enabled modules
        try:
            # Fukui analysis
            if self.enabled_modules.get("fukui", False):
                fukui_result, fukui_warnings = self.run_fukui_analysis(
                    input_file, atom_a, atom_b
                )
                if fukui_result:
                    result.fukui_fplus_a = fukui_result.get('fukui_fplus')
                    result.fukui_fplus_b = fukui_result.get('fukui_fplus')
                    result.fukui_fminus_a = fukui_result.get('fukui_fminus')
                    result.fukui_fminus_b = fukui_result.get('fukui_fminus')
                    result.fukui_f0_a = fukui_result.get('fukui_f0')
                    result.fukui_f0_b = fukui_result.get('fukui_f0')
                result.warnings.extend(fukui_warnings)

            # Dual descriptor analysis
            if self.enabled_modules.get("dual_descriptor", False):
                dd_result, dd_warnings = self.run_dual_descriptor_analysis(
                    input_file, atom_a, atom_b
                )
                if dd_result is not None:
                    result.dual_descriptor_a = dd_result
                    result.dual_descriptor_b = dd_result
                result.warnings.extend(dd_warnings)

            # QTAIM BCP analysis
            if self.enabled_modules.get("qtaim_bcp", False):
                qtaim_result, qtaim_warnings = self.run_qtaim_bcp_analysis(
                    input_file, atom_a, atom_b
                )
                if qtaim_result:
                    result.rho_bcp_forming1 = qtaim_result.get('rho_bcp')
                    result.laplacian_bcp_forming1 = qtaim_result.get('laplacian_bcp')
                result.warnings.extend(qtaim_warnings)

            # Mark success if any module produced results
            result.success = any([
                result.fukui_fplus_a is not None,
                result.dual_descriptor_a is not None,
                result.rho_bcp_forming1 is not None
            ])

        except Exception as e:
            result.error_message = str(e)
            result.warnings.append(W_MULTIWFN_FAILED)

        return result


def run_multiwfn_feature_extraction(
    input_file: Path,
    atom_a: int,
    atom_b: int,
    multiwfn_config: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """Convenience function for Multiwfn feature extraction.

    Args:
        input_file: Input file path
        atom_a: First atom index (1-based)
        atom_b: Second atom index (1-based)
        multiwfn_config: Configuration dict with keys:
            - multiwfn_path: str
            - cache_dir: str
            - timeout_sec: int
            - enabled_modules: dict

    Returns:
        Dict of Multiwfn features with mw_* prefix
    """
    config = multiwfn_config or {}

    runner = MultiwfnRunner(
        multiwfn_path=config.get('multiwfn_path', 'Multiwfn'),
        cache_dir=Path(config.get('cache_dir', '.cache/step4_multiwfn')),
        timeout_sec=config.get('timeout_sec', 120),
        enabled_modules=config.get('enabled_modules', {
            'fukui': True,
            'dual_descriptor': True,
            'qtaim_bcp': False,
            'nci': False
        })
    )

    result = runner.run_complete_analysis(input_file, atom_a, atom_b)

    # Convert to feature dict
    features = {
        'mw_fukui_fplus_atomA': result.fukui_fplus_a,
        'mw_fukui_fplus_atomB': result.fukui_fplus_b,
        'mw_fukui_fminus_atomA': result.fukui_fminus_a,
        'mw_fukui_fminus_atomB': result.fukui_fminus_b,
        'mw_fukui_f0_atomA': result.fukui_f0_a,
        'mw_fukui_f0_atomB': result.fukui_f0_b,
        'mw_dual_descriptor_atomA': result.dual_descriptor_a,
        'mw_dual_descriptor_atomB': result.dual_descriptor_b,
        'mw_rho_bcp_forming1': result.rho_bcp_forming1,
        'mw_laplacian_bcp_forming1': result.laplacian_bcp_forming1,
        'mw_status': 'ok' if result.success else 'failed',
        'mw_missing_reason': None if result.success else (
            result.error_message or
            (result.warnings[0] if result.warnings else 'unknown')
        ),
        'mw_warnings_count': len(result.warnings),
    }

    # Remove None values for cleaner output
    features = {k: v for k, v in features.items() if v is not None}

    return features
