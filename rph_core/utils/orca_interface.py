"""
ORCA Interface Module
=====================

ORCA 量子化学软件接口 - 用于高精度单点能计算

Author: QC Descriptors Team
Date: 2026-01-10
Session: #1 - ORCAInterface._generate_input()
Session: #2 - ORCAInterface._parse_output()
Session: #3 - ORCAInterface._find_orca_binary() + _run_orca()
Session: #4 - ORCAInterface.single_point()
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import subprocess
import shutil
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Mapping, Optional, Tuple, Union

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from rph_core.utils.optimization_config import OptimizationConfig

from rph_core.utils.resource_utils import (
    find_executable,
    setup_ld_library_path,
    mem_to_mb,
    calc_orca_maxcore,
    resolve_executable_config
)
from rph_core.utils.geometry_tools import CoordinateExtractor
from rph_core.utils.data_types import QCResult
from rph_core.utils.keyword_translator import KeywordTranslator
from rph_core.utils.method_registry import MethodRegistry, NormalizedMethodSpec
from rph_core.utils.qc_models import (
    IRCJobSpec,
    NEBJobSpec,
    QCJobResult,
    QCJobSpec,
    SurfaceScanSpec,
)
from rph_core.utils.capability_validator import CapabilityValidator
from rph_core.utils.orca_failure_classifier import (
    OrcaFailureClassification,
    OrcaFailureType,
    classify_orca_failure_for_config,
    has_normal_termination,
)
from rph_core.utils.orca_input_renderer import OrcaInputRenderer

logger = logging.getLogger(__name__)

_ORCA_VERSION_RE = re.compile(r"Program\s+Version\s+([^\n\r]+)", re.IGNORECASE)


def _parse_orca_major_version(version_str: Optional[str]) -> Optional[int]:
    """Extract the ORCA major version from a version string like ``6.1.1`` -> 6."""
    if not version_str:
        return None
    match = re.search(r"(\d+)(?:\.\d+)*", str(version_str).strip())
    if not match:
        return None
    try:
        return int(match.group(1))
    except (ValueError, IndexError):
        return None


def _resolve_orca_version_from_path(path: Optional[Path]) -> Optional[int]:
    """Fallback version detection from the ORCA directory name (e.g. ``/opt/orca_6_1_1``)."""
    if path is None:
        return None
    parent_name = path.parent.name
    compact = parent_name.replace("_", ".").replace("-", ".")
    return _parse_orca_major_version(compact)


@dataclass(frozen=True)
class _OrcaRunMetadata:
    returncode: Optional[int] = None
    stderr_text: str = ""
    timed_out: bool = False


def _parse_orca_angstrom_lines(coord_lines: List[str]) -> Optional[NDArray[np.float64]]:
    """Parse ORCA CARTESIAN COORDINATES (ANGSTROEM) block lines into (N,3) array.

    ORCA angstrom coordinate lines have format: ``[element] [x] [y] [z]``.
    Returns ``None`` if no valid coordinate lines are found.
    """
    coords: List[List[float]] = []
    for line in coord_lines:
        parts = line.strip().split()
        if len(parts) >= 4:
            try:
                coords.append([float(parts[-3]), float(parts[-2]), float(parts[-1])])
            except (ValueError, IndexError):
                continue
    if coords:
        return np.array(coords)
    return None


def _write_orca_xyz_file(
    xyz_path: Path,
    coords: NDArray[np.float64],
    input_xyz: Optional[Path] = None,
) -> None:
    """Write coordinates to a plain XYZ file, inferring element symbols from input_xyz."""
    symbols: List[str] = ["X"] * coords.shape[0]
    if input_xyz is not None and input_xyz.exists():
        try:
            lines = input_xyz.read_text(encoding="utf-8").splitlines()
            if len(lines) > 2:
                symbol_lines = lines[2:2 + coords.shape[0]]
                for i, line in enumerate(symbol_lines):
                    parts = line.split()
                    if parts:
                        symbols[i] = parts[0]
        except Exception:
            pass
    n = coords.shape[0]
    xyz_lines = [f"{n}", "ORCA fallback coordinates"]
    for i in range(n):
        xyz_lines.append(f"{symbols[i]:2s}  {coords[i,0]:12.6f}  {coords[i,1]:12.6f}  {coords[i,2]:12.6f}")
    xyz_path.write_text("\n".join(xyz_lines) + "\n", encoding="utf-8")


def _orca_xyz_text_from_coords(
    coords: NDArray[np.float64],
    input_xyz: Optional[Path] = None,
) -> str:
    symbols: List[str] = ["X"] * coords.shape[0]
    if input_xyz is not None and input_xyz.exists():
        try:
            lines = input_xyz.read_text(encoding="utf-8").splitlines()
            if len(lines) > 2:
                for index, line in enumerate(lines[2:2 + coords.shape[0]]):
                    parts = line.split()
                    if parts:
                        symbols[index] = parts[0]
        except Exception:
            pass
    xyz_lines = [f"{coords.shape[0]}", "ORCA fallback coordinates"]
    for index in range(coords.shape[0]):
        xyz_lines.append(
            f"{symbols[index]:2s}  {coords[index,0]:12.6f}  {coords[index,1]:12.6f}  {coords[index,2]:12.6f}"
        )
    return "\n".join(xyz_lines) + "\n"


def _expected_atom_count(input_xyz: Optional[Path]) -> int | None:
    if input_xyz is None or not input_xyz.exists():
        return None
    try:
        lines = input_xyz.read_text(encoding="utf-8").splitlines()
        return int(lines[0].strip()) if lines else None
    except (OSError, ValueError, IndexError, UnicodeDecodeError):
        return None


def extract_orca_version_from_output(output_file: Optional[Path]) -> Optional[str]:
    if output_file is None:
        return None
    try:
        content = Path(output_file).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    match = _ORCA_VERSION_RE.search(content)
    return match.group(1).strip() if match else None


def _render_irc_block(spec: IRCJobSpec) -> str:
    direction = str(spec.direction or "both").strip().lower()
    if direction not in {"both", "forward", "backward"}:
        raise ValueError(f"Unsupported ORCA IRC direction: {spec.direction!r}")

    init_hessian = str(spec.init_hessian or "read").strip().lower()
    if init_hessian not in {"read", "calc_anfreq", "calc_numfreq"}:
        raise ValueError(f"Unsupported ORCA IRC InitHess mode: {spec.init_hessian!r}")

    init_displacement_mode = str(spec.init_displacement_mode or "energy").strip().lower()
    if init_displacement_mode not in {"energy", "length"}:
        raise ValueError(
            f"Unsupported ORCA IRC initial displacement mode: {spec.init_displacement_mode!r}"
        )

    lines = [
        "%irc",
        f"  MaxIter {int(spec.max_iter)}",
        f"  Direction {direction}",
        f"  InitHess {init_hessian}",
    ]
    if spec.hessian_filename:
        lines.append(f'  Hess_Filename "{spec.hessian_filename}"')
    lines.append(f"  HessMode {int(spec.hessian_mode)}")
    if init_displacement_mode == "energy":
        lines.append(f"  DE_INIT_DISPL {float(spec.initial_delta_energy_mEh)}")
    lines.append(f"  SCALE_INIT_DISPL {float(spec.scale_initial_displacement)}")
    lines.append(f"  SCALE_DISPL_SD {float(spec.scale_steepest_descent)}")
    if spec.adaptive_step:
        lines.append("  ADAPT_SCALE_DISPL true")
    lines.append("end")
    return "\n".join(lines)


def _irc_direction_for_filename(path: Path) -> Optional[str]:
    lowered = path.name.lower()
    if "forward" in lowered or re.search(r"(^|[_\-.])fwd([_\-.]|$)", lowered):
        return "forward"
    if (
        "backward" in lowered
        or "reverse" in lowered
        or re.search(r"(^|[_\-.])bwd([_\-.]|$)", lowered)
    ):
        return "backward"
    return None


def _extract_last_xyz_frame_text(xyz_like_file: Path) -> Optional[str]:
    try:
        lines = xyz_like_file.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        return None

    last_frame: Optional[List[str]] = None
    index = 0
    while index < len(lines):
        stripped = lines[index].strip()
        if not stripped or stripped.startswith("---"):
            index += 1
            continue
        try:
            n_atoms = int(stripped)
        except ValueError:
            index += 1
            continue

        frame_end = index + n_atoms + 2
        if frame_end > len(lines):
            break
        last_frame = lines[index:frame_end]
        index = frame_end

    if not last_frame:
        return None
    return "\n".join(last_frame).rstrip() + "\n"


def _materialize_irc_endpoint(source_file: Path, destination: Path) -> Optional[Path]:
    frame_text = _extract_last_xyz_frame_text(source_file)
    if frame_text is None:
        return None
    destination.write_text(frame_text, encoding="utf-8")
    return destination


def _collect_irc_candidate_files(output_dir: Path, *, exclude_names: Tuple[str, ...]) -> List[Path]:
    patterns = ("*.xyz", "*.irc", "*.path")
    seen: set[Path] = set()
    candidates: List[Path] = []
    for pattern in patterns:
        for candidate in sorted(output_dir.glob(pattern)):
            if (
                not candidate.is_file()
                or candidate.name in exclude_names
                or candidate in seen
            ):
                continue
            lowered = candidate.name.lower()
            if not any(token in lowered for token in ("irc", "forward", "backward", "path")):
                continue
            seen.add(candidate)
            candidates.append(candidate)
    return candidates


def _resolve_irc_endpoint_sources(
    output_dir: Path,
    *,
    direction: str,
    exclude_names: Tuple[str, ...],
) -> List[Path]:
    candidates = _collect_irc_candidate_files(output_dir, exclude_names=exclude_names)
    forward = [path for path in candidates if _irc_direction_for_filename(path) == "forward"]
    backward = [path for path in candidates if _irc_direction_for_filename(path) == "backward"]
    other = [path for path in candidates if _irc_direction_for_filename(path) is None]

    selected: List[Path] = []
    if direction in {"both", "forward"} and forward:
        selected.append(forward[0])
    if direction in {"both", "backward"} and backward:
        selected.append(backward[0])

    fallback_pool = [path for path in other + forward + backward if path not in selected]
    required = 2 if direction == "both" else 1
    while len(selected) < required and fallback_pool:
        selected.append(fallback_pool.pop(0))
    return selected


def extract_last_orca_geometry_xyz_text(
    output_file: Optional[Path],
    input_xyz: Optional[Path] = None,
) -> Optional[str]:
    if output_file is None:
        return None

    output_path = Path(output_file)
    sibling_xyz = output_path.with_suffix(".xyz")
    if sibling_xyz.is_file():
        try:
            text = sibling_xyz.read_text(encoding="utf-8")
            lines = text.splitlines()
            if lines:
                parsed_atoms = int(lines[0].strip())
                expected_atoms = _expected_atom_count(input_xyz)
                if expected_atoms is None or parsed_atoms == expected_atoms:
                    return text if text.endswith("\n") else f"{text}\n"
        except (OSError, ValueError, IndexError, UnicodeDecodeError):
            pass

    try:
        content = output_path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None

    coord_blocks = re.findall(
        r"CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n-+\n((?:\s*[A-Za-z]{1,3}\s+[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+(?:\n|$))+)",
        content,
        re.DOTALL,
    )
    if not coord_blocks:
        return None
    coord_lines = [line for line in coord_blocks[-1].split("\n") if line.strip()]
    coords = _parse_orca_angstrom_lines(coord_lines)
    if coords is None:
        return None
    expected_atoms = _expected_atom_count(input_xyz)
    if expected_atoms is not None and coords.shape[0] != expected_atoms:
        return None
    return _orca_xyz_text_from_coords(coords, input_xyz)


def _parse_orca_displacement_vectors(content: str) -> Dict[int, NDArray[np.float64]]:
    """Parse imaginary-frequency displacement vectors from an ORCA output file.

    ORCA prints a ``CARTESIAN DISPLACEMENTS`` section after frequency
    analysis.  Each mode is a block::

        CARTESIAN DISPLACEMENTS
        -----------
        Mode:   0
        Freq:   -123.45 cm**-1 (imaginary mode)
                dx          dy          dz
        Atom 0:  0.123456   0.234567   0.345678
        Atom 1: -0.012345   0.023456  -0.034567

    Only modes with negative (imaginary) frequencies are collected.

    Returns
    -------
    dict[int, np.ndarray]
        *mode_index* (0 for the first imaginary mode) → (*N_atoms*, 3) array.
        Empty dict on parse errors or when no imaginary modes exist.
    """
    result: Dict[int, NDArray[np.float64]] = {}

    section_match = re.search(
        r"CARTESIAN DISPLACEMENTS\s*\n\s*-+\s*\n",
        content,
    )
    if not section_match:
        return result

    section = content[section_match.end():]

    imaginary_count = 0
    for mode_block in re.split(r"\n\s*(?=Mode:\s*\d+\s*\n)", section):
        mode_match = re.search(r"Mode:\s*(\d+)", mode_block)
        freq_match = re.search(r"Freq:\s*([\d\.\-+Ee]+)", mode_block)
        if not mode_match or not freq_match:
            continue

        try:
            freq = float(freq_match.group(1))
        except ValueError:
            continue

        if freq >= 0:
            continue

        disp_lines: List[Tuple[float, float, float]] = []
        for line in mode_block.split("\n"):
            atom_match = re.match(r"\s*Atom\s+\d+\s*:\s+(.*)", line)
            if not atom_match:
                continue
            parts = atom_match.group(1).split()
            if len(parts) < 3:
                continue
            try:
                dx, dy, dz = float(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue
            disp_lines.append((dx, dy, dz))

        if not disp_lines:
            continue

        n_atoms = len(disp_lines)
        array = np.zeros((n_atoms, 3))
        for atom_j, (dx, dy, dz) in enumerate(disp_lines):
            array[atom_j, 0] = dx
            array[atom_j, 1] = dy
            array[atom_j, 2] = dz

        result[imaginary_count] = array
        imaginary_count += 1

    if result:
        logger.debug(
            "Extracted %d imaginary-mode displacement vector(s) from ORCA output",
            len(result),
        )
    return result


_ELEMENT_MASS: Dict[str, float] = {
    "H": 1.008, "He": 4.003, "Li": 6.941, "Be": 9.012, "B": 10.811, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305,
    "Al": 26.982, "Si": 28.086, "P": 30.974, "S": 32.065, "Cl": 35.453, "Ar": 39.948,
    "K": 39.098, "Ca": 40.078, "Fe": 55.845, "Cu": 63.546, "Zn": 65.380,
    "Br": 79.904, "I": 126.904,
}


def parse_orca6_normal_mode_vectors(
    content: str,
    element_symbols: Sequence[str],
) -> Dict[int, NDArray[np.float64]]:
    """Parse ORCA 6 ``NORMAL MODES`` blocks into mass-deweighted displacement vectors.

    ORCA 6 writes the mass-weighted mode matrix::

        NORMAL MODES
        ------------
        These modes are the Cartesian displacements weighted by the diagonal
        matrix M(i,i)=1/sqrt(m[i]) ...
        <mode-index header row>
        <coordinate-row> <value>...

    Each block repeats a header row of up to six mode indices, followed by one
    row per Cartesian coordinate (3 per atom).  Modes are mass-weighted, so the
    returned vectors are deweighted by ``sqrt(m[atom])`` and renormalised to
    unit norm, matching the semantics of ``_parse_orca_displacement_vectors``
    (Cartesian displacement mode vectors).

    Returns
    -------
    dict[int, np.ndarray]
        *mode_index* → (*N_atoms*, 3) deweighted displacement vectors.
        Empty dict on parse errors.
    """
    header_pattern = re.compile(r"^\s*\d+(\s+\d+)*\s*$")
    lines = content.splitlines()
    start = next((i for i, line in enumerate(lines) if line.strip() == "NORMAL MODES"), None)
    if start is None:
        return {}

    masses = [
        _ELEMENT_MASS.get(str(symbol).strip().capitalize(), 12.011)
        for symbol in element_symbols
    ]

    modes: Dict[int, List[float]] = {}
    active_modes: List[int] = []
    for line in lines[start + 1:]:
        if header_pattern.match(line):
            active_modes = [int(token) for token in line.split()]
            continue
        tokens = line.split()
        if not tokens or not tokens[0].isdigit():
            continue
        row = int(tokens[0])
        if active_modes and row == 0 and not modes:
            modes = {mode: [] for mode in active_modes}
            for offset, value in enumerate(tokens[1:]):
                if offset < len(active_modes):
                    modes[active_modes[offset]].append(float(value))
            continue
        if not active_modes:
            continue
        for offset, value in enumerate(tokens[1:]):
            if offset < len(active_modes):
                modes.setdefault(active_modes[offset], []).append(float(value))

    if not modes:
        return {}

    n_atoms = len(masses)
    result: Dict[int, NDArray[np.float64]] = {}
    for mode_index, values in modes.items():
        if len(values) != 3 * n_atoms:
            continue
        raw = np.asarray(values, dtype=np.float64).reshape((n_atoms, 3))
        deweighted = raw * np.sqrt(np.asarray(masses, dtype=np.float64))[:, None]
        norm = float(np.linalg.norm(deweighted))
        if norm <= 1e-12:
            continue
        result[mode_index] = deweighted / norm
    if result:
        logger.debug(
            "Extracted %d ORCA6 normal-mode displacement vector(s)",
            len(result),
        )
    return result



class ORCAInterface:
    """ORCA 接口 - 高精度单点能计算"""

    def __init__(
        self,
        method: str = "M062X",
        basis: str = "def2-TZVPP",
        aux_basis: str = "def2/J",
        nprocs: int = 16,
        maxcore: Optional[int] = None,
        solvent: str = "acetone",
        route_extras: str = "",
        scf_maxiter: Optional[int] = None,
        solvent_model: Optional[str] = None,
        orca_binary_path: Optional[str] = None,
        config: Optional[Dict[str, Any]] = None,
        method_spec: Optional[Dict[str, Any]] = None,
    ):
        """
        初始化 ORCA 接口

        Args:
            method: DFT 方法 (如 M062X, B3LYP, PWPB95)
            basis: 基组 (如 def2-TZVPP, def2-SVP)
            aux_basis: 辅助基组 (如 def2/J)
            nprocs: 并行进程数
            maxcore: 每核最大内存 (MB，None 则从配置计算）
            solvent: 溶剂名称 (用于 SMD 模型)
            orca_binary_path: ORCA 可执行文件路径 (可选)
            config: 配置字典 (可选，用于派生路径和内存）
        """
        raw_method_spec: Dict[str, Any] = dict(method_spec or {})
        if not raw_method_spec:
            raw_method_spec = {
                "engine": "orca",
                "method": method,
                "basis": basis,
                "aux_basis": aux_basis,
                "route_extras": route_extras,
                "solvent": solvent,
            }
            if maxcore is not None:
                raw_method_spec["maxcore"] = maxcore
        if solvent_model:
            raw_method_spec["solvent_model"] = solvent_model

        self.method_alias = str(raw_method_spec.get("method", method) or method)
        self.method_profile = MethodRegistry.get_profile(self.method_alias, engine="orca")
        self.method_spec: NormalizedMethodSpec = MethodRegistry.normalize_spec(
            raw_method_spec,
            default_engine="orca",
            default_basis=basis,
        )

        self.method = self.method_spec.method
        self.dispersion = self.method_spec.dispersion.keyword
        self.basis = self.method_spec.basis
        self.aux_basis = self.method_spec.aux_basis or aux_basis
        self.nprocs = nprocs
        self.solvent = self.method_spec.solvent or solvent
        self.solvent_model = self._resolve_solvent_model(config, raw_method_spec)
        if self.method_spec.family == "semiempirical_xtb" and self.solvent_model != "ALPB":
            # ORCA's GFN-xTB family accepts only ALPB (CPCM/SMD abort the run).
            _logger = logging.getLogger(f"{__name__}.{method}/{basis}")
            _logger.warning(
                "[ORCA] %s supports only ALPB solvation; translating solvent_model=%s -> ALPB",
                self.method,
                self.solvent_model,
            )
            self.solvent_model = "ALPB"
        self.route_extras = self.method_spec.route_extras
        self.scf_maxiter = int(scf_maxiter) if scf_maxiter is not None else None

        # 处理 maxcore：如果未提供，从配置派生
        if maxcore is None and config:
            res_cfg = config.get('resources', {})
            mem = res_cfg.get('mem', '32GB')
            safety_factor = res_cfg.get('orca_maxcore_safety', 0.65)
            self.maxcore = calc_orca_maxcore(mem, nprocs, safety_factor)
        else:
            self.maxcore = maxcore if maxcore is not None else 4000

        # 设置日志（必须在使用 self.logger 的任何辅助方法之前初始化）
        self.logger = logging.getLogger(f"{__name__}.{method}/{basis}")

        validation_errors = CapabilityValidator.validate(self.method_spec, task_type="sp")
        if validation_errors:
            self.logger.warning("ORCA method spec validation issues: %s", "; ".join(validation_errors))

        # 查找 ORCA 二进制文件（集成新的配置系统）
        self.orca_binary = self._find_orca_binary(orca_binary_path, config)
        self._orca_major_version = self._detect_orca_major_version(config)

        self.config = config

        # SP cache (in-memory, per-process)
        self._sp_cache = {}
        self._sp_cache_hits = 0
        self._sp_cache_misses = 0
        self._last_resolved_final_xyz: Optional[Path] = None
        self._last_orca_run_metadata = _OrcaRunMetadata()
        self._last_subprocess: subprocess.Popen[Any] | None = None

    def _pal_block(self, nprocs: Optional[int] = None) -> str:
        effective_nprocs = self.nprocs if nprocs is None else nprocs
        return f"%pal nprocs {effective_nprocs} end\n" if effective_nprocs > 1 else ""

    def _compute_sp_cache_key(
        self,
        xyz_file: Path,
        charge: Optional[int] = None,
        spin: Optional[int] = None
    ) -> Optional[str]:
        try:
            content = xyz_file.read_bytes()
        except Exception as e:
            self.logger.debug(f"SP cache skipped: failed to read {xyz_file}: {e}")
            return None

        geometry_hash = hashlib.sha256(content).hexdigest()

        if charge is None or spin is None:
            inferred_charge, inferred_spin = 0, 1
            if xyz_file.suffix == '.out':
                try:
                    inferred_charge, inferred_spin = CoordinateExtractor.get_charge_spin_from_orca_out(xyz_file)
                except Exception as e:
                    self.logger.debug(f"SP cache charge/spin fallback for {xyz_file}: {e}")
            charge = inferred_charge if charge is None else charge
            spin = inferred_spin if spin is None else spin

        payload = {
            "geometry_hash": geometry_hash,
            "method": self.method,
            "basis": self.basis,
            "aux_basis": self.aux_basis,
            "solvent": self.solvent,
            "solvent_model": self.solvent_model,
            "route_extras": self.route_extras,
            "charge": charge,
            "spin": spin
        }

        payload_str = json.dumps(payload, sort_keys=True, default=str)
        return hashlib.sha256(payload_str.encode("utf-8")).hexdigest()

    def _log_sp_cache_stats(self) -> None:
        total = self._sp_cache_hits + self._sp_cache_misses
        if total == 0:
            return
        hit_rate = (self._sp_cache_hits / total) * 100.0
        self.logger.debug(
            f"SP cache hit rate: {self._sp_cache_hits}/{total} ({hit_rate:.1f}%)"
        )

    def _get_renderer(self) -> OrcaInputRenderer:
        return OrcaInputRenderer(self.method_spec)

    def render_simple_keywords(self, task_type: str = "sp") -> str:
        return self._get_renderer().render_simple_keywords(task_type=task_type)

    def render_blocks(self) -> Dict[str, str]:
        blocks = dict(self._get_renderer().render_blocks())
        if self.scf_maxiter is not None:
            blocks["scf"] = "\n".join(
                [
                    "%scf",
                    f"  MaxIter {int(self.scf_maxiter)}",
                    "end",
                ]
            )
        return blocks

    @staticmethod
    def _resolve_solvent_model(config: Optional[Dict[str, Any]], method_spec: Mapping[str, Any]) -> str:
        candidates: List[Any] = [method_spec.get("solvent_model")]
        raw_spec = method_spec.get("raw")
        if isinstance(raw_spec, Mapping):
            candidates.append(raw_spec.get("solvent_model"))
        if isinstance(config, Mapping):
            theory = config.get("theory")
            if isinstance(theory, Mapping):
                theory_solvent = theory.get("solvent")
                if isinstance(theory_solvent, Mapping):
                    candidates.append(theory_solvent.get("model"))
            top_solvent = config.get("solvent")
            if isinstance(top_solvent, Mapping):
                candidates.append(top_solvent.get("model"))

        for candidate in candidates:
            if isinstance(candidate, str) and candidate.strip():
                return candidate.strip().upper()
        return "CPCM"

    def _solvent_name_for_orca(self) -> str:
        from rph_core.utils.solvent_map import orca_smd_solvent

        return orca_smd_solvent(self.solvent)

    def _render_route(self, task_type: str = "sp") -> str:
        route = self.render_simple_keywords(task_type=task_type)
        if not self.solvent or str(self.solvent).upper() == "NONE":
            return route
        if self.solvent_model in {"CPCM", "PCM", "ALPB"}:
            model_keyword = "ALPB" if self.solvent_model == "ALPB" else "CPCM"
            solvent_keyword = f"{model_keyword}({self._solvent_name_for_orca()})"
            if solvent_keyword.lower() not in route.lower():
                return f"{route} {solvent_keyword}"
        return route

    def _render_cpcm_block(self) -> str:
        if not self.solvent or str(self.solvent).upper() == "NONE":
            return ""
        if self.solvent_model not in {"SMD", "CPCM-SMD"}:
            return ""

        solvent_name = self._solvent_name_for_orca()
        return (
            "\n %cpcm\n"
            "    smd true\n"
            f"    SMDsolvent \"{solvent_name}\"\n"
            " end\n "
        )

    @staticmethod
    def _render_named_blocks(blocks: Dict[str, str]) -> str:
        if not blocks:
            return ""
        ordered_keys = ("basis", "method", "scf", "mdci")
        chunk: List[str] = []
        for key in ordered_keys:
            value = blocks.get(key)
            if value:
                chunk.append(value.rstrip())
        for key, value in blocks.items():
            if key in ordered_keys:
                continue
            if value:
                chunk.append(value.rstrip())
        if not chunk:
            return ""
        return "\n" + "\n".join(chunk) + "\n"

    def _generate_input(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: Optional[int] = None,
        spin: Optional[int] = None,
        extra_input_block: Optional[str] = None,
        job_tag: Optional[str] = None
    ) -> Path:
        """
        生成 ORCA 输入文件

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录

        Returns:
            生成的 .inp 文件路径
        """
        route = self._render_route(task_type="sp")
        cpcm_block = self._render_cpcm_block()
        rendered_blocks = self._render_named_blocks(self.render_blocks())
 
        if charge is None or spin is None:
            inferred_charge, inferred_spin = 0, 1
            if xyz_file.suffix == '.out':
                try:
                    inferred_charge, inferred_spin = CoordinateExtractor.get_charge_spin_from_orca_out(xyz_file)
                    self.logger.debug(
                        f"从 {xyz_file.name} 提取电荷/自旋: {inferred_charge}/{inferred_spin}"
                    )
                except Exception as e:
                    self.logger.warning(f"无法从 {xyz_file.name} 提取电荷/自旋: {e}")
            charge = inferred_charge if charge is None else charge
            spin = inferred_spin if spin is None else spin

        if job_tag is None:
            job_tag = uuid.uuid4().hex[:8]

        base_name = f"{xyz_file.stem}_{job_tag}"

        # 组装完整输入文件内容 - 直接嵌入xyz坐标
        import shutil
        xyz_copy = output_dir / f"{base_name}{xyz_file.suffix}"
        shutil.copy(xyz_file, xyz_copy)

        # 读取xyz文件内容，跳过前两行（原子数和标题）
        xyz_lines = xyz_copy.read_text().split('\n')
        if len(xyz_lines) > 2:
            xyz_content = '\n'.join(xyz_lines[2:])
        else:
            xyz_content = xyz_copy.read_text()

        extra_block = ""
        if extra_input_block:
            extra_block = f"\n{extra_input_block}\n"

        pal_block = self._pal_block()

        inp_content = f"""{route}
 %maxcore {self.maxcore}
 {pal_block}{cpcm_block}{rendered_blocks}
{extra_block}
  * xyz {charge} {spin}
{xyz_content}
  *
"""

        # 写入文件
        inp_file = output_dir / f"{base_name}.inp"
        inp_file.write_text(inp_content)

        return inp_file

    def constrained_optimize(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int,
        spin: int,
        constraints_block: str,
        timeout: Optional[int] = None
    ) -> QCResult:
        """
        ORCA constrained optimization with custom %geom block.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        route = self._render_route(task_type="opt")
        cpcm_block = self._render_cpcm_block()
        rendered_blocks = self._render_named_blocks(self.render_blocks())

        import shutil
        job_tag = uuid.uuid4().hex[:8]
        base_name = f"{xyz_file.stem}_constrained_opt_{job_tag}"
        xyz_copy = output_dir / f"{base_name}{xyz_file.suffix}"
        shutil.copy(xyz_file, xyz_copy)

        xyz_lines = xyz_copy.read_text().split('\n')
        if len(xyz_lines) > 2:
            xyz_content = '\n'.join(xyz_lines[2:])
        else:
            xyz_content = xyz_copy.read_text()

        pal_block = self._pal_block()

        inp_content = f"""{route}
 %maxcore {self.maxcore}
 {pal_block}{cpcm_block}{rendered_blocks}
{constraints_block}
  * xyz {charge} {spin}
{xyz_content}
  *
"""

        inp_file = output_dir / f"{base_name}.inp"
        inp_file.write_text(inp_content)

        try:
            out_file = self._run_orca(inp_file, output_dir, timeout=timeout)
            result = self._parse_output(
                out_file,
                require_optimization_convergence=True,
                returncode=self._last_orca_run_metadata.returncode,
                stderr_text=self._last_orca_run_metadata.stderr_text,
                timed_out=self._last_orca_run_metadata.timed_out,
            )

            result.coordinates = self._extract_final_orca_coordinates(out_file, xyz_copy)

            return result

        except Exception as e:
            self.logger.error(f"ORCA constrained optimization failed: {e}")
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=str(e)
            )

    def optimize_with_task(
        self,
        xyz_file: Path,
        output_dir: Path,
        task_type: str,
        charge: int,
        spin: int,
        geom_block: Optional[str] = None,
        max_cycles_opt: Optional[int] = None,
        timeout: Optional[int] = None,
        job_tag: Optional[str] = None,
        subprocess_callback: Callable[[subprocess.Popen[Any]], None] | None = None,
    ) -> QCResult:
        """Run an ORCA optimization with an explicit ORCA task type."""
        normalized_task = str(task_type or "").strip().lower()
        if normalized_task not in {"opt", "opt_freq", "ts", "ts_freq"}:
            raise ValueError(
                f"Unsupported ORCA task_type={task_type!r}; "
                "expected one of: opt, opt_freq, ts, ts_freq"
            )

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        route = self._render_route(task_type=normalized_task)
        cpcm_block = self._render_cpcm_block()
        rendered_blocks = self._render_named_blocks(self.render_blocks())

        if job_tag is None:
            job_tag = uuid.uuid4().hex[:8]

        base_name = f"{xyz_file.stem}_{normalized_task}_{job_tag}"
        xyz_copy = output_dir / f"{base_name}{xyz_file.suffix}"
        shutil.copy(xyz_file, xyz_copy)

        xyz_lines = xyz_copy.read_text().split('\n')
        if len(xyz_lines) > 2:
            xyz_content = '\n'.join(xyz_lines[2:])
        else:
            xyz_content = xyz_copy.read_text()

        pal_block = self._pal_block()
        rendered_geom_block = ""
        if geom_block and geom_block.strip():
            rendered_geom_block = f"{geom_block.rstrip()}\n"

        inp_content = f"""{route}
 %maxcore {self.maxcore}
 {pal_block}{cpcm_block}{rendered_blocks}
{rendered_geom_block}  * xyz {charge} {spin}
{xyz_content}
  *
"""

        inp_file = output_dir / f"{base_name}.inp"
        inp_file.write_text(inp_content)

        try:
            out_file = self._run_orca(
                inp_file,
                output_dir,
                timeout=timeout,
                subprocess_callback=subprocess_callback,
            )
            result = self._parse_output(
                out_file,
                require_optimization_convergence=True,
                returncode=self._last_orca_run_metadata.returncode,
                stderr_text=self._last_orca_run_metadata.stderr_text,
                timed_out=self._last_orca_run_metadata.timed_out,
                max_cycles_opt=max_cycles_opt,
            )
            result.coordinates = self._extract_final_orca_coordinates(out_file, xyz_copy)

            if normalized_task in {"ts", "ts_freq", "opt_freq"}:
                result.frequencies = self._extract_frequencies_from_output(out_file)

            return result

        except Exception as e:
            self.logger.error("ORCA task optimization failed (%s): %s", normalized_task, e)
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=str(e)
            )

    def run_irc(
        self,
        spec: IRCJobSpec,
        input_xyz: Path,
        output_dir: Path,
        charge: int = 0,
        spin: int = 1,
        timeout: Optional[int] = None,
        subprocess_callback: Optional[Callable[[Any], None]] = None,
    ) -> QCJobResult:
        """Execute an ORCA IRC calculation for TS endpoint discovery."""

        input_xyz = Path(input_xyz)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if not str(spec.method or "").strip():
            raise ValueError("IRCJobSpec.method is required for ORCA IRC jobs")

        route = self._render_route(task_type="irc")
        cpcm_block = self._render_cpcm_block()
        rendered_blocks = self._render_named_blocks(self.render_blocks())
        irc_block = _render_irc_block(spec)

        job_tag = uuid.uuid4().hex[:8]
        base_name = f"{input_xyz.stem}_irc_{job_tag}"
        xyz_copy = output_dir / f"{base_name}{input_xyz.suffix}"
        shutil.copy(input_xyz, xyz_copy)

        xyz_lines = xyz_copy.read_text(encoding="utf-8").splitlines()
        xyz_content = "\n".join(xyz_lines[2:]) if len(xyz_lines) > 2 else xyz_copy.read_text(encoding="utf-8")
        pal_block = self._pal_block()
        inp_content = (
            f"{route}\n"
            f" %maxcore {self.maxcore}\n"
            f" {pal_block}{cpcm_block}{rendered_blocks}\n"
            f"{irc_block}\n"
            f"  * xyz {int(charge)} {int(spin)}\n"
            f"{xyz_content}\n"
            "  *\n"
        )
        inp_file = output_dir / f"{base_name}.inp"
        inp_file.write_text(inp_content, encoding="utf-8")

        try:
            out_file = self._run_orca(
                inp_file,
                output_dir,
                timeout=timeout,
                subprocess_callback=subprocess_callback,
            )
            parsed = self._parse_output(
                out_file,
                require_optimization_convergence=False,
                returncode=self._last_orca_run_metadata.returncode,
                stderr_text=self._last_orca_run_metadata.stderr_text,
                timed_out=self._last_orca_run_metadata.timed_out,
            )
        except Exception as exc:
            self.logger.error("ORCA IRC execution failed: %s", exc)
            return QCJobResult(
                status="failed",
                input_xyz=input_xyz,
                output_xyz=input_xyz,
                error=str(exc),
            )

        if not parsed.converged:
            return QCJobResult(
                status="failed",
                input_xyz=input_xyz,
                output_xyz=input_xyz,
                output_file=getattr(parsed, "output_file", out_file),
                energy_hartree=getattr(parsed, "energy", None),
                error=getattr(parsed, "error_message", None),
                failure_type=getattr(parsed, "failure_type", None),
                failure_evidence_lines=getattr(parsed, "failure_evidence_lines", None),
                returncode=self._last_orca_run_metadata.returncode,
                stderr_text=self._last_orca_run_metadata.stderr_text,
                timed_out=self._last_orca_run_metadata.timed_out,
                orca_version=extract_orca_version_from_output(getattr(parsed, "output_file", out_file)),
            )

        normalized_direction = str(spec.direction or "both").strip().lower()
        endpoint_sources = _resolve_irc_endpoint_sources(
            output_dir,
            direction=normalized_direction,
            exclude_names=(xyz_copy.name, inp_file.name, out_file.name),
        )
        required_endpoints = 2 if normalized_direction == "both" else 1
        if len(endpoint_sources) < required_endpoints:
            error_message = (
                f"Unable to locate ORCA IRC endpoint geometries in {output_dir}; "
                f"inspect {out_file} manually"
            )
            self.logger.warning(error_message)
            return QCJobResult(
                status="failed",
                input_xyz=input_xyz,
                output_xyz=input_xyz,
                output_file=out_file,
                energy_hartree=getattr(parsed, "energy", None),
                error=error_message,
                returncode=self._last_orca_run_metadata.returncode,
                stderr_text=self._last_orca_run_metadata.stderr_text,
                timed_out=self._last_orca_run_metadata.timed_out,
                orca_version=extract_orca_version_from_output(out_file),
            )

        endpoint_a_xyz = output_dir / "irc_endpoint_a.xyz"
        materialized_a = _materialize_irc_endpoint(endpoint_sources[0], endpoint_a_xyz)
        endpoint_b_xyz: Optional[Path] = None
        materialized_b: Optional[Path] = None
        if required_endpoints == 2:
            endpoint_b_xyz = output_dir / "irc_endpoint_b.xyz"
            materialized_b = _materialize_irc_endpoint(endpoint_sources[1], endpoint_b_xyz)

        if materialized_a is None or (required_endpoints == 2 and materialized_b is None):
            error_message = (
                f"Failed to parse ORCA IRC endpoint geometries from {output_dir}; "
                f"inspect {out_file} manually"
            )
            self.logger.warning(error_message)
            return QCJobResult(
                status="failed",
                input_xyz=input_xyz,
                output_xyz=input_xyz,
                output_file=out_file,
                energy_hartree=getattr(parsed, "energy", None),
                error=error_message,
                returncode=self._last_orca_run_metadata.returncode,
                stderr_text=self._last_orca_run_metadata.stderr_text,
                timed_out=self._last_orca_run_metadata.timed_out,
                orca_version=extract_orca_version_from_output(out_file),
            )

        extra: Dict[str, Any] = {
            "endpoint_a_xyz": materialized_a,
            "irc_trajectory_dir": output_dir,
        }
        if materialized_b is not None:
            extra["endpoint_b_xyz"] = materialized_b

        return QCJobResult(
            status="complete",
            input_xyz=input_xyz,
            output_xyz=input_xyz,
            output_file=out_file,
            energy_hartree=getattr(parsed, "energy", None),
            error=getattr(parsed, "error_message", None),
            failure_type=getattr(parsed, "failure_type", None),
            failure_evidence_lines=getattr(parsed, "failure_evidence_lines", None),
            returncode=self._last_orca_run_metadata.returncode,
            stderr_text=self._last_orca_run_metadata.stderr_text,
            timed_out=self._last_orca_run_metadata.timed_out,
            orca_version=extract_orca_version_from_output(out_file),
            extra=extra,
        )

    def run_surface_scan(
        self,
        spec: SurfaceScanSpec,
        input_xyz: Path,
        output_dir: Path,
        subprocess_callback: Optional[Callable[[Any], None]] = None,
    ) -> QCJobResult:
        """Run an ORCA relaxed scan and preserve every optimized scan frame.

        This is deliberately an interface primitive.  S2 owns when a scan is
        warranted and how a peak is interpreted; the interface only renders
        ORCA input, executes it, and returns recoverable scan artefacts.
        """

        input_xyz = Path(input_xyz)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        if not input_xyz.exists():
            return QCJobResult(status="failed", input_xyz=input_xyz, error="Surface-scan input XYZ is unavailable")
        if not spec.coordinates or len(spec.coordinates) > 3:
            return QCJobResult(status="failed", input_xyz=input_xyz, error="Surface scan requires one to three coordinates")
        if spec.scan_ts and len(spec.coordinates) != 1:
            return QCJobResult(status="failed", input_xyz=input_xyz, error="ScanTS supports one local coordinate only")

        job_tag = uuid.uuid4().hex[:8]
        base_name = f"{input_xyz.stem}_surface_scan_{job_tag}"
        input_copy = output_dir / f"{base_name}_start.xyz"
        shutil.copy(input_xyz, input_copy)
        scan_lines: List[str] = []
        for coordinate in spec.coordinates:
            kind = str(coordinate.kind).upper()
            expected = {"B": 2, "A": 3, "D": 4}.get(kind)
            atoms = tuple(int(atom) for atom in coordinate.atoms)
            if expected is None or len(atoms) != expected or min(atoms, default=-1) < 0:
                return QCJobResult(status="failed", input_xyz=input_xyz, error=f"Invalid scan coordinate: {coordinate!r}")
            if int(coordinate.steps) < 3:
                return QCJobResult(status="failed", input_xyz=input_xyz, error="Surface scan requires at least three points")
            scan_lines.append(
                f"    {kind} {' '.join(str(atom) for atom in atoms)} = "
                f"{float(coordinate.start):.8f}, {float(coordinate.end):.8f}, {int(coordinate.steps)}"
            )

        route = self._render_route(task_type="opt").replace(" noautostart", "")
        if spec.scan_ts:
            route = re.sub(r"\\bOpt\\b", "ScanTS", route, count=1, flags=re.IGNORECASE)
            if "scants" not in route.lower():
                route = f"{route} ScanTS"
        geom_lines = ["%geom", "  Scan", *scan_lines, "  end"]
        if spec.simultaneous and len(spec.coordinates) > 1:
            geom_lines.append("  Simul_Scan true")
        if spec.scan_ts and spec.full_scan:
            geom_lines.append("  FullScan true")
        if spec.max_cycles is not None:
            geom_lines.append(f"  MaxIter {int(spec.max_cycles)}")
        geom_lines.append("end")
        cpcm_block = self._render_cpcm_block()
        rendered_blocks = self._render_named_blocks(self.render_blocks())
        pal_block = self._pal_block()
        inp_file = output_dir / f"{base_name}.inp"
        geom_block = "\n".join(geom_lines)
        inp_file.write_text(
            f"{route}\n%maxcore {self.maxcore}\n{pal_block}{cpcm_block}{rendered_blocks}\n"
            f"{geom_block}\n"
            f"* xyzfile {int(spec.charge)} {int(spec.multiplicity)} {input_copy.name}\n",
            encoding="utf-8",
        )
        if self.method_spec.family == "semiempirical_xtb" and self.orca_binary is not None:
            orca_dir = self.orca_binary.parent
            if not (orca_dir / "otool_xtb").is_file() and not (orca_dir / "xtb").is_file():
                return QCJobResult(
                    status="failed",
                    input_xyz=input_xyz,
                    error=(
                        "ORCA GFN-xTB requires the xtb binary in the ORCA directory; "
                        f"neither {orca_dir / 'otool_xtb'} nor {orca_dir / 'xtb'} exists. "
                        f"Fix: ln -s <path-to-xtb> {orca_dir / 'xtb'}"
                    ),
                )
        try:
            out_file = self._run_orca(
                inp_file, output_dir, timeout=spec.timeout, subprocess_callback=subprocess_callback
            )
            out_text = out_file.read_text(encoding="utf-8", errors="ignore")
        except Exception as exc:
            return QCJobResult(status="failed", input_xyz=input_xyz, error=str(exc))

        allxyz_files = sorted(output_dir.glob(f"{base_name}*.allxyz"))
        frames: List[Path] = []
        if allxyz_files:
            lines = allxyz_files[-1].read_text(encoding="utf-8", errors="ignore").splitlines()
            cursor = 0
            frame_dir = output_dir / "scan_frames"
            frame_dir.mkdir(parents=True, exist_ok=True)
            while cursor < len(lines):
                try:
                    atom_count = int(lines[cursor].strip())
                except ValueError:
                    cursor += 1
                    continue
                end = cursor + atom_count + 2
                if atom_count <= 0 or end > len(lines):
                    break
                frame = frame_dir / f"frame_{len(frames):04d}.xyz"
                frame.write_text("\n".join(lines[cursor:end]) + "\n", encoding="utf-8")
                frames.append(frame)
                cursor = end
        # ORCA 5 writes the optimized energy of each relaxed-scan point to
        # ``*.relaxscanact.dat``.  The main output contains many inner OPT
        # cycles and may omit/compact some ``FINAL SINGLE POINT ENERGY`` lines
        # under ``miniprint``; it is therefore neither a complete nor a
        # one-to-one frame-energy ledger.  Prefer the dedicated scan ledger.
        scan_energy_file: Optional[Path] = None
        scan_energies: List[float] = []
        for candidate_file in sorted(output_dir.glob(f"{base_name}*.relaxscanact.dat")):
            parsed = self._parse_relaxed_scan_energy_ledger(
                candidate_file,
                coordinate_count=len(spec.coordinates),
            )
            if len(parsed) >= len(scan_energies):
                scan_energy_file = candidate_file
                scan_energies = parsed
        energy_values = [
            float(value)
            for value in re.findall(
                r"FINAL\\s+SINGLE\\s+POINT\\s+ENERGY\\s+([-+]?\\d+(?:\\.\\d+)?(?:[Ee][-+]?\\d+)?)",
                out_text,
                flags=re.IGNORECASE,
            )
        ]
        energies: List[Optional[float]] = [None] * len(frames)
        if len(scan_energies) >= len(frames):
            energies = [float(value) for value in scan_energies[:len(frames)]]
        elif len(energy_values) >= len(frames):
            energies = [float(value) for value in energy_values[:len(frames)]]
        ts_candidates = sorted(
            [*output_dir.glob(f"{base_name}*TS*.xyz"), *output_dir.glob(f"{base_name}*opt*.xyz")]
        )
        candidate = ts_candidates[-1] if spec.scan_ts and ts_candidates else None
        normal = has_normal_termination(out_text)
        status = "complete" if normal and frames else "failed"
        error = None if status == "complete" else "ORCA relaxed scan did not yield optimized frames"
        if status == "failed":
            classification = classify_orca_failure_for_config(
                output_text=out_text,
                stderr_text=self._last_orca_run_metadata.stderr_text,
                returncode=self._last_orca_run_metadata.returncode,
                timed_out=self._last_orca_run_metadata.timed_out,
                config=self.config,
            )
            if classification is not None:
                error = classification.summary
                if classification.failure_type == OrcaFailureType.XTB_BINARY_UNAVAILABLE:
                    orca_dir = self.orca_binary.parent if self.orca_binary is not None else Path(".")
                    error = (
                        f"{classification.summary}: {orca_dir}. "
                        f"Fix: ln -s <path-to-xtb> {orca_dir / 'xtb'}"
                    )
        return QCJobResult(
            status=status,
            input_xyz=input_xyz,
            output_xyz=candidate,
            output_file=out_file,
            energy_hartree=energy_values[-1] if energy_values else None,
            error=error,
            returncode=self._last_orca_run_metadata.returncode,
            stderr_text=self._last_orca_run_metadata.stderr_text,
            timed_out=self._last_orca_run_metadata.timed_out,
            orca_version=extract_orca_version_from_output(out_file),
            extra={
                "frames": [str(frame) for frame in frames],
                "energies_hartree": energies,
                "energy_source": (
                    str(scan_energy_file)
                    if scan_energy_file is not None and len(scan_energies) >= len(frames)
                    else "output_final_single_point_energy"
                ),
                "scan_ts_candidate_xyz": str(candidate) if candidate else None,
                "scan_coordinate_count": len(spec.coordinates),
            },
        )

    @staticmethod
    def _parse_relaxed_scan_energy_ledger(
        path: Path,
        *,
        coordinate_count: int,
    ) -> List[float]:
        """Read ORCA's per-point optimized-energy ledger for a relaxed scan."""

        energies: List[float] = []
        try:
            for raw_line in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
                values = raw_line.split()
                if len(values) < coordinate_count + 1:
                    continue
                try:
                    energies.append(float(values[coordinate_count]))
                except ValueError:
                    continue
        except OSError:
            return []
        return energies

    def run_neb_ts(
        self,
        spec: QCJobSpec,
        input_xyz: Path,
        output_dir: Path,
        subprocess_callback: Optional[Callable[[Any], None]] = None,
    ) -> QCJobResult:
        """Run ORCA's double-ended NEB-TS and preserve the best TS seed.

        This intentionally lives in the low-level interface so refinement code
        never invokes ORCA directly.  A partially converged NEB is returned as
        ``partial`` when a final highest-energy image can be materialized.
        """
        input_xyz = Path(input_xyz)
        end_xyz = Path(spec.neb_end_xyz) if spec.neb_end_xyz else None
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        if end_xyz is None or not end_xyz.exists():
            return QCJobResult(status="failed", input_xyz=input_xyz, error="NEB end-point XYZ is unavailable")

        route = self._render_route(task_type="sp").replace(" noautostart", "")
        route = route.replace("! ", "! ", 1)
        route = route.replace(self.method, f"{self.method} {spec.route or 'LOOSE-NEB-TS'}", 1)
        cpcm_block = self._render_cpcm_block()
        rendered_blocks = self._render_named_blocks(self.render_blocks())
        job_tag = uuid.uuid4().hex[:8]
        base_name = f"{input_xyz.stem}_neb_ts_{job_tag}"
        start_copy = output_dir / f"{base_name}_start.xyz"
        end_copy = output_dir / f"{base_name}_end.xyz"
        shutil.copy(input_xyz, start_copy)
        shutil.copy(end_xyz, end_copy)
        # ORCA's NEB parser does not reliably honour quoted absolute paths
        # containing whitespace.  The job runner may move this input into a
        # clean sandbox, so all NEB references must be local sandbox names.
        start_ref = start_copy.name
        end_ref = end_copy.name
        ts_guess_line = ""
        if spec.neb_ts_guess_xyz and Path(spec.neb_ts_guess_xyz).exists():
            guess_copy = output_dir / f"{base_name}_guess.xyz"
            shutil.copy(Path(spec.neb_ts_guess_xyz), guess_copy)
            ts_guess_line = f'  NEB_TS_XYZFILE "{guess_copy.name}"\n'
        neb_block = (
            "%neb\n"
            f'  NEB_END_XYZFILE "{end_ref}"\n'
            f"{ts_guess_line}"
            "end\n"
        )
        pal_block = self._pal_block()
        inp_content = (
            f"{route}\n %maxcore {self.maxcore}\n {pal_block}{cpcm_block}{rendered_blocks}\n"
            # ORCA 5 accepts quoted filenames in the %neb assignments but
            # treats quotes as literal characters in `* xyzfile`.  The
            # sandbox-local basename has no whitespace, so it must be emitted
            # unquoted here.
            f"{neb_block}* xyzfile {int(spec.charge)} {int(spec.multiplicity)} {start_ref}\n"
        )
        inp_file = output_dir / f"{base_name}.inp"
        inp_file.write_text(inp_content, encoding="utf-8")
        try:
            out_file = self._run_orca(inp_file, output_dir, timeout=spec.timeout, subprocess_callback=subprocess_callback)
            text = out_file.read_text(encoding="utf-8", errors="ignore")
        except Exception as exc:
            return QCJobResult(status="failed", input_xyz=input_xyz, error=str(exc))

        # ORCA writes the climbing-image result as ``*_NEB-CI_converged.xyz``
        # for several NEB-TS variants.  That file is the intended TS seed even
        # when the output does not contain the ordinary OptTS convergence
        # banner, so do not silently discard it in favour of the warm-up
        # geometry.  A converged CI is a *candidate*, not yet a validated TS;
        # refinement performs the subsequent OptTS/FREQ validation.
        converged_candidates = sorted(
            [
                *output_dir.glob("*_NEB-TS_converged.xyz"),
                *output_dir.glob("*_NEB-CI_converged.xyz"),
            ]
        )
        partial_candidates = sorted(
            [
                *output_dir.glob("*_NEB-TS*.xyz"),
                *output_dir.glob("*_NEB-CI*.xyz"),
            ]
        )
        candidate = (converged_candidates or partial_candidates or [None])[-1]
        status = "complete" if converged_candidates else ("partial" if candidate else "failed")
        return QCJobResult(
            status=status,
            input_xyz=input_xyz,
            output_xyz=candidate,
            output_file=out_file,
            error=None if status == "complete" else "ORCA NEB-TS did not fully converge",
            extra={
                "candidate_xyz": str(candidate) if candidate else None,
                "candidate_kind": (
                    "neb_ci_converged" if candidate and "_NEB-CI_converged.xyz" in candidate.name
                    else "neb_ts_converged" if candidate and "_NEB-TS_converged.xyz" in candidate.name
                    else "partial_neb_image" if candidate else None
                ),
                "start_xyz": str(start_copy),
                "end_xyz": str(end_copy),
            },
        )

    def _parse_output(
        self,
        out_file: Path,
        *,
        require_optimization_convergence: bool = False,
        stderr_text: Optional[str] = None,
        returncode: Optional[int] = None,
        timed_out: Optional[bool] = None,
        max_cycles_opt: Optional[int] = None,
        scf_converged: Optional[bool] = None,
    ) -> QCResult:
        """
        解析 ORCA 输出文件

        Args:
            out_file: ORCA 输出文件路径

        Returns:
            QCResult 对象
        """
        run_metadata = getattr(self, "_last_orca_run_metadata", _OrcaRunMetadata())
        effective_stderr_text = run_metadata.stderr_text if stderr_text is None else stderr_text
        effective_returncode = run_metadata.returncode if returncode is None else returncode
        effective_timed_out = run_metadata.timed_out if timed_out is None else timed_out

        def _failure_result(
            message: str,
            *,
            classification: Optional[OrcaFailureClassification] = None,
        ) -> QCResult:
            result = QCResult(
                energy=0.0,
                converged=False,
                output_file=out_file,
                error_message=message,
            )
            if classification is not None:
                result.failure_type = classification.failure_type.value
                result.failure_evidence_lines = classification.evidence_lines
                setattr(result, "failure_retry_hint", classification.retry_hint)
            return result

        # 读取输出文件内容
        try:
            content = out_file.read_text(encoding="utf-8", errors="ignore")
        except Exception as e:
            return QCResult(
                energy=0.0,
                converged=False,
                output_file=out_file,
                error_message=f"Unable to read ORCA output file: {e}",
            )

        optimization_converged = "THE OPTIMIZATION HAS CONVERGED" in content
        classification = classify_orca_failure_for_config(
            output_text=content,
            stderr_text=effective_stderr_text,
            returncode=effective_returncode,
            timed_out=effective_timed_out,
            max_cycles_opt=max_cycles_opt,
            scf_converged=scf_converged,
            config=getattr(self, "config", None),
        )
        normal_termination = has_normal_termination(content)
        geometry_partial_success = False

        if classification is not None:
            geometry_partial_success = (
                normal_termination
                and require_optimization_convergence
                and not optimization_converged
                and classification.failure_type == OrcaFailureType.GEOMETRY_OPTIMIZATION_NOT_CONVERGED
            )
            if not geometry_partial_success:
                return _failure_result(classification.summary, classification=classification)

        # 检查是否正常终止
        if not normal_termination:
            if classification is not None:
                return _failure_result(classification.summary, classification=classification)
            return _failure_result("ORCA did not terminate normally")

        # 提取最终单点能
        energy_match = re.search(
            r"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            content,
        )
        if not energy_match:
            return QCResult(
                energy=0.0,
                converged=False,
                output_file=out_file,
                error_message="Could not find ORCA final single-point energy",
            )

        # 解析能量
        try:
            energy = float(energy_match.group(1))
        except ValueError:
            return QCResult(
                energy=0.0,
                converged=False,
                output_file=out_file,
                error_message=f"Invalid ORCA energy value: {energy_match.group(1)}",
            )

        coordinates = None
        try:
            coord_block = re.search(
                r'CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n-+\n((?:\s*[A-Za-z]{1,3}\s+[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+(?:\n|$))+)',
                content,
                re.DOTALL
            )
            if coord_block:
                coord_lines = [l for l in coord_block.group(1).split('\n') if l.strip()]
                coords = _parse_orca_angstrom_lines(coord_lines)
                if coords is not None:
                    coordinates = coords
                else:
                    self.logger.warning("Failed to parse ORCA coordinates from %s", out_file)
            else:
                self.logger.warning("No ORCA coordinate block found in %s", out_file)
        except Exception as e:
            self.logger.warning("Failed to extract ORCA coordinates from %s: %s", out_file, e)

        # A normal ORCA process exit is not equivalent to an optimized
        # stationary point. Preserve the terminal MaxIter classification so
        # downstream S3 rescue can distinguish it from a generic failure.
        result = QCResult(
            energy=energy,
            converged=(
                not require_optimization_convergence
                or optimization_converged
            ),
            coordinates=coordinates,
            output_file=out_file,
            error_message=(
                None
                if not require_optimization_convergence
                or optimization_converged
                else "ORCA terminated normally, but the geometry optimization did not converge"
            ),
            failure_type=(
                classification.failure_type.value
                if classification is not None and not optimization_converged
                else None
            ),
            failure_evidence_lines=(
                classification.evidence_lines
                if classification is not None and not optimization_converged
                else None
            ),
            failure_retry_hint=(
                classification.retry_hint
                if classification is not None and not optimization_converged
                else None
            ),
            optimization_converged=(
                optimization_converged if require_optimization_convergence else None
            ),
            stop_reason=(
                "native_maxiter_checkpoint"
                if geometry_partial_success
                else None
            ),
        )
        return result

    def extract_frequencies_from_output(self, out_file: Path) -> Optional[NDArray[np.float64]]:
        """Extract vibrational frequencies from a completed ORCA output."""

        return self._extract_frequencies_from_output(out_file)

    def _extract_frequencies_from_output(self, out_file: Path) -> Optional[NDArray[np.float64]]:
        """Extract the last ORCA vibrational-frequency block from an output file."""
        try:
            content = out_file.read_text(encoding="utf-8", errors="ignore")
        except OSError as exc:
            self.logger.warning("Failed to read ORCA output for frequency extraction %s: %s", out_file, exc)
            return None

        last_freq_pos = content.rfind("VIBRATIONAL FREQUENCIES")
        if last_freq_pos < 0:
            return None

        block_content = content[last_freq_pos:]
        freqs = re.findall(r'[-]?\d+\.\d+(?=\s+cm\*\*-1)', block_content)
        if not freqs:
            return None

        try:
            return np.array([float(freq) for freq in freqs], dtype=float)
        except ValueError:
            return None

    def _resolve_orca_final_xyz(
        self,
        out_file: Path,
        input_xyz: Optional[Path] = None,
    ) -> Optional[Path]:
        """Resolve ORCA final geometry: prefer sibling .xyz, fallback to .out parsing.

        ORCA writes the final optimized geometry to a ``.xyz`` file with the
        same stem as the output file (e.g. ``ts_guess_ts_opt_6693bdc1.xyz``).
        This method returns that sibling file when it exists and passes
        atom-count validation.  Only when no sibling ``.xyz`` is found does it
        fall back to parsing ``CARTESIAN COORDINATES (ANGSTROEM)`` blocks from
        the ``.out`` file.

        Returns:
            Path to the resolved final XYZ file, or ``None`` if unavailable.
        """
        sibling_xyz = out_file.with_suffix(".xyz")
        if sibling_xyz.is_file():
            try:
                lines = sibling_xyz.read_text(encoding="utf-8").splitlines()
                n_atoms = int(lines[0].strip()) if lines else 0
            except (ValueError, IndexError, UnicodeDecodeError) as exc:
                self.logger.debug("Cannot read sibling xyz %s: %s", sibling_xyz, exc)
            else:
                # If we have an input_xyz for validation, check atom count.
                if input_xyz is not None:
                    try:
                        ref_lines = input_xyz.read_text(encoding="utf-8").splitlines()
                        ref_atoms = int(ref_lines[0].strip()) if ref_lines else 0
                    except Exception:
                        ref_atoms = 0
                    if n_atoms > 0 and ref_atoms > 0 and n_atoms != ref_atoms:
                        self.logger.warning(
                            "ORCA sibling xyz atom count mismatch in %s: %d vs %d",
                            sibling_xyz, n_atoms, ref_atoms,
                        )
                        # Do not return fallback — prefer to return nothing unambiguous.
                        return None
                return sibling_xyz

        # Fallback: parse final coordinate block from .out
        self.logger.debug("No ORCA final sibling xyz for %s; parsing .out coordinates", out_file)
        coords = self._extract_orca_coordinates_from_out(out_file, input_xyz)
        if coords is not None:
            # Write a canonical final xyz file so that downstream code can use a
            # real file path instead of trusting an in-memory array.
            fallback_xyz = out_file.with_suffix(".rph_fallback.xyz")
            try:
                _write_orca_xyz_file(fallback_xyz, coords, input_xyz)
                return fallback_xyz
            except Exception:
                pass
        return None

    def _extract_orca_coordinates_from_out(
        self,
        out_file: Path,
        input_xyz: Optional[Path] = None,
    ) -> Optional[NDArray[np.float64]]:
        """Extract the LAST Cartesian coordinate block from an ORCA .out file.

        Prefer :meth:`_resolve_orca_final_xyz` for canonical geometry —
        this method is a fallback for parsing.
        """
        try:
            content = out_file.read_text(encoding="utf-8")
        except Exception as e:
            self.logger.warning("Failed to read ORCA output for coordinate extraction %s: %s", out_file, e)
            return None

        # Match through the separator line to coordinate lines, stop at blank line.
        # ORCA coordinate lines have format: element x y z (uppercase element starters
        # like C, H, N would fool a naive \n[A-Z] terminator).
        coord_blocks = re.findall(
            r'CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n-+\n((?:\s*[A-Za-z]{1,3}\s+[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+(?:\n|$))+)',
            content,
            re.DOTALL
        )
        if not coord_blocks:
            self.logger.warning("No ORCA coordinate blocks found in %s", out_file)
            return None

        coord_lines = [line for line in coord_blocks[-1].split('\n') if line.strip()]
        coordinates = _parse_orca_angstrom_lines(coord_lines)
        if coordinates is None:
            self.logger.warning("Failed to parse final ORCA coordinate block in %s", out_file)
            return None

        if input_xyz is not None:
            try:
                xyz_lines = [line for line in input_xyz.read_text(encoding="utf-8").splitlines() if line.strip()]
                expected_atoms = int(xyz_lines[0]) if xyz_lines else 0
            except Exception:
                expected_atoms = 0
            if expected_atoms > 0 and coordinates.shape[0] != expected_atoms:
                self.logger.warning(
                    "ORCA coordinate atom count mismatch in %s: parsed %d, expected %d from %s",
                    out_file, coordinates.shape[0], expected_atoms, input_xyz,
                )
                return None

        return coordinates

    def _extract_final_orca_coordinates(
        self,
        out_file: Path,
        input_xyz: Path
    ) -> Optional[NDArray[np.float64]]:
        """Legacy wrapper — prefer :meth:`_resolve_orca_final_xyz`.

        Returns the coordinates as a numpy array when a final .xyz file is
        available; otherwise falls back to parsing from the .out file.
        """
        final_xyz = self._resolve_orca_final_xyz(out_file, input_xyz)
        if final_xyz is not None:
            return self._read_xyz_coordinates(final_xyz)
        return self._extract_orca_coordinates_from_out(out_file, input_xyz)

    def _read_xyz_coordinates(self, xyz_path: Path) -> Optional[NDArray[np.float64]]:
        """Read coordinates from a plain XYZ file, returning (N,3) array."""
        try:
            from rph_core.utils.file_io import read_xyz as _io_read_xyz
            coords, _symbols = _io_read_xyz(xyz_path)
            return coords.astype(np.float64)
        except Exception as e:
            self.logger.warning("Failed to read XYZ coordinates from %s: %s", xyz_path, e)
            return None

    def _detect_orca_major_version(self, config: Optional[Dict[str, Any]] = None) -> Optional[int]:
        """Resolve the active ORCA major version from config or the binary path."""
        if config:
            orca_cfg = config.get("executables", {}).get("orca", {}) or {}
            if isinstance(orca_cfg, Mapping) and orca_cfg.get("version"):
                version = _parse_orca_major_version(str(orca_cfg.get("version")))
                if version is not None:
                    self.logger.debug("ORCA major version from config: %s", version)
                    return version
        path_version = _resolve_orca_version_from_path(self.orca_binary)
        if path_version is not None:
            self.logger.debug("ORCA major version from binary path: %s", path_version)
        return path_version

    def _find_orca_binary(self, provided_path: Optional[str] = None, config: Optional[Dict[str, Any]] = None) -> Optional[Path]:
        """
        查找 ORCA 可执行文件（集成新的配置系统）

        查找顺序:
        1. 提供的路径
        2. 配置文件中的路径
        3. 环境变量 $ORCA_PATH, $ORCA_BIN
        4. PATH 中的 'orca' 命令

        Args:
            provided_path: 用户提供的 ORCA 路径
            config: 配置字典

        Returns:
            ORCA 可执行文件的 Path 对象，如果找不到返回 None
        """
        if provided_path:
            path = Path(provided_path)
            if path.exists() and path.is_file():
                logger.info(f"使用提供的 ORCA 路径: {path}")
                return path
            else:
                logger.warning(f"提供的 ORCA 路径不存在: {provided_path}")

        exe_config = resolve_executable_config(
            config or {}, 'orca', env_vars=['ORCA_PATH', 'ORCA_BIN']
        )
        if exe_config.get('found'):
            return exe_config['path']

        logger.error("未找到 ORCA 可执行文件")
        return None

    def _build_orca_runtime_env(self) -> Dict[str, str]:
        env = os.environ.copy()

        # ORCA parallelism is controlled explicitly by %pal.  Inheriting a
        # larger OpenMP/BLAS thread count makes each MPI rank spawn additional
        # threads and invalidates the scheduler's core budget.  One math thread
        # per rank is the safe default and remains configurable for dedicated
        # benchmarks.
        resources = dict((self.config or {}).get("resources", {}) or {})
        math_threads = max(1, int(resources.get("orca_math_threads_per_rank", 1)))
        for key in (
            "OMP_NUM_THREADS",
            "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            env[key] = str(math_threads)
        env["OMP_MAX_ACTIVE_LEVELS"] = "1"
        env["OMP_THREAD_LIMIT"] = str(math_threads)

        if self.orca_binary is None:
            return env

        orca_dir = self.orca_binary.parent

        mpi_bin_from_config: Optional[str] = None
        mpi_lib_from_config: Optional[str] = None
        if self.config:
            orca_cfg = self.config.get("executables", {}).get("orca", {})
            mpi_bin_from_config = orca_cfg.get("mpi_bin_dir")
            mpi_lib_from_config = orca_cfg.get("mpi_lib_dir")

        mpi_bin_candidates: List[Path] = []
        if mpi_bin_from_config:
            mpi_bin_candidates.append(Path(mpi_bin_from_config))
            self.logger.debug(f"ORCA MPI bin dir from config: {mpi_bin_from_config}")
        mpi_bin_candidates += [
            orca_dir / "openmpi418" / "bin",
            orca_dir / "openmpi" / "bin",
            orca_dir / "mpi" / "bin",
            orca_dir / "bin",
            orca_dir,
        ]

        # Also derive MPI bin candidates from LD_LIBRARY_PATH entries
        # e.g. /opt/openmpi418/lib -> /opt/openmpi418/bin
        current_ld = env.get("LD_LIBRARY_PATH", "")
        for ld_entry in current_ld.split(":"):
            ld_entry = ld_entry.strip()
            if not ld_entry:
                continue
            ld_path = Path(ld_entry)
            if ld_path.name == "lib" or ld_path.name.startswith("lib"):
                sibling_bin = ld_path.parent / "bin"
                if sibling_bin not in mpi_bin_candidates:
                    mpi_bin_candidates.append(sibling_bin)

        selected_mpi_bin: Optional[Path] = None
        for candidate in mpi_bin_candidates:
            mpirun_path = candidate / "mpirun"
            if mpirun_path.exists() and mpirun_path.is_file():
                selected_mpi_bin = candidate
                break

        # Warn if configured MPI bin dir is missing or has no mpirun
        if mpi_bin_from_config:
            config_path = Path(mpi_bin_from_config)
            if not config_path.exists() or not (config_path / "mpirun").exists():
                self.logger.warning(
                    f"Configured MPI bin dir not found or missing mpirun: {mpi_bin_from_config}"
                )

        current_path = env.get("PATH", "")
        path_parts = [p for p in current_path.split(":") if p]
        preferred_path_parts = [str(orca_dir)]
        if selected_mpi_bin is not None:
            preferred_path_parts.append(str(selected_mpi_bin))
            self.logger.debug(f"ORCA runtime PATH prefixed with MPI bin: {selected_mpi_bin}")
        env["PATH"] = ":".join(
            preferred_path_parts + [
                p for p in path_parts
                if p not in preferred_path_parts and not ("openmpi" in p.lower() and p != str(selected_mpi_bin))
            ]
        )

        mpi_lib_candidates: List[Path] = []
        if mpi_lib_from_config:
            mpi_lib_candidates.append(Path(mpi_lib_from_config))
            self.logger.debug(f"ORCA MPI lib dir from config: {mpi_lib_from_config}")
        mpi_lib_candidates += [
            orca_dir / "openmpi418" / "lib",
            orca_dir / "openmpi" / "lib",
            orca_dir / "mpi" / "lib",
            orca_dir / "lib",
            orca_dir,
        ]

        selected_lib_paths = [str(p) for p in mpi_lib_candidates if p.exists() and p.is_dir()]
        if selected_mpi_bin is not None:
            selected_mpi_lib = selected_mpi_bin.parent / "lib"
            selected_mpi_lib_str = str(selected_mpi_lib)
            if selected_mpi_lib.exists() and selected_mpi_lib.is_dir() and selected_mpi_lib_str not in selected_lib_paths:
                selected_lib_paths.insert(0, selected_mpi_lib_str)
        if str(orca_dir) not in selected_lib_paths:
            selected_lib_paths.insert(0, str(orca_dir))

        if selected_lib_paths:
            current_ld = env.get("LD_LIBRARY_PATH", "")
            ld_parts = [p for p in current_ld.split(":") if p]
            merged = selected_lib_paths + [
                p for p in ld_parts
                if p not in selected_lib_paths and not ("openmpi" in p.lower() and p not in selected_lib_paths)
            ]
            env["LD_LIBRARY_PATH"] = ":".join(merged)

        env.setdefault("ORCA_PATH", str(orca_dir))
        env.setdefault("ORCA_DIR", str(orca_dir))

        # ORCA >= 6 launches its parallel modules from the serial driver itself;
        # the OMPI_ALLOW_RUN_AS_ROOT injection is only meaningful for the 5.x
        # external-mpirun path.  WSL2/container shared-memory tuning is handled
        # by OMPI_MCA_btl_vader_single_copy_mechanism=none at the environment
        # level instead.
        major_version = getattr(self, "_orca_major_version", None)
        self.logger.debug("Detected ORCA major version for runtime env: %s", major_version)
        if os.getuid() == 0 and (major_version is None or major_version < 6):
            env.setdefault("OMPI_ALLOW_RUN_AS_ROOT", "1")
            env.setdefault("OMPI_ALLOW_RUN_AS_ROOT_CONFIRM", "1")

        return env

    def _execute_orca_process(
        self,
        *,
        cmd: List[str],
        cwd: Path,
        out_file: Path,
        env: Dict[str, str],
        timeout: Optional[int],
        subprocess_callback: Callable[[subprocess.Popen[Any]], None] | None = None,
    ) -> _OrcaRunMetadata:
        stderr_text = ""
        timed_out = False

        with open(out_file, 'w') as out_f:
            process = subprocess.Popen(
                cmd,
                stdout=out_f,
                stderr=subprocess.PIPE,
                cwd=str(cwd),
                text=True,
                env=env,
                # A live S3 monitor must be able to terminate ORCA together
                # with its MPI children without touching unrelated jobs.
                start_new_session=(os.name == "posix"),
            )
            self._last_subprocess = process
            # The refinement monitor needs the active file when a toxic-path
            # sandbox is in use.  Keep this private, process-local metadata so
            # the public callback signature remains backward compatible.
            setattr(process, "_rph_orca_output_path", str(out_file))
            if subprocess_callback is not None:
                try:
                    subprocess_callback(process)
                except Exception as exc:
                    self.logger.warning("Ignoring ORCA subprocess callback failure: %s", exc)

            try:
                try:
                    _, stderr_text = process.communicate(timeout=timeout)
                except subprocess.TimeoutExpired:
                    process.kill()
                    timed_out = True
                    _, stderr_text = process.communicate()
            finally:
                self._last_subprocess = None

        return _OrcaRunMetadata(
            returncode=process.returncode,
            stderr_text=stderr_text or "",
            timed_out=timed_out,
        )

    def _log_orca_failure_diagnostics(
        self,
        out_file: Path,
        *,
        stderr_text: str,
        returncode: Optional[int],
        timed_out: bool,
        sandbox_dir: Optional[Path] = None,
    ) -> None:
        if stderr_text.strip():
            stderr_tail = stderr_text.splitlines()[-50:]
            self.logger.error("ORCA stderr:\n%s", "\n".join(stderr_tail))

        try:
            if out_file.exists():
                tail_lines = out_file.read_text(errors='ignore').splitlines()[-50:]
                if tail_lines:
                    label = "ORCA sandbox run" if sandbox_dir is not None else "ORCA run"
                    detail = "timed out" if timed_out else f"failed (rc={returncode})"
                    self.logger.error(
                        "%s %s; last 50 lines of output:\n%s",
                        label,
                        detail,
                        "\n".join(tail_lines),
                    )
        except Exception as tail_e:
            self.logger.warning("Failed to read ORCA output tail: %s", tail_e)

        if sandbox_dir is not None:
            self.logger.error("ORCA sandbox dir preserved for debugging: %s", sandbox_dir)

    def _run_orca(
        self,
        inp_file: Path,
        output_dir: Path,
        timeout: Optional[int] = 3600,
        subprocess_callback: Callable[[subprocess.Popen[Any]], None] | None = None,
    ) -> Path:
        """
        运行 ORCA 计算

        Args:
            inp_file: ORCA 输入文件
            output_dir: 输出目录
            timeout: 超时时间（秒），默认 1 小时

        Returns:
            ORCA 输出文件路径

        Raises:
            RuntimeError: ORCA 二进制文件未找到
            RuntimeError: ORCA 进程启动失败
        """
        # 检查 ORCA 是否可用
        if self.orca_binary is None:
            raise RuntimeError(
                "ORCA 二进制文件未找到。请通过以下方式之一指定:\n"
                "  1. 设置环境变量 $ORCA_PATH\n"
                "  2. 确保 'orca' 在 PATH 中\n"
                "  3. 在初始化时提供 orca_binary_path 参数"
            )

        out_file = inp_file.with_suffix('.out')
        self._last_orca_run_metadata = _OrcaRunMetadata()

        # ORCA / MPI 对含空格或特殊字符的路径支持较差；检测到 toxic path 时使用临时目录运行
        from rph_core.utils.path_compat import is_toxic_path

        if is_toxic_path(output_dir.resolve()) or is_toxic_path(inp_file.resolve()):
            import tempfile
            import shutil

            self.logger.debug("检测到路径包含空格/特殊字符，使用临时目录运行 ORCA")

            # Use a unique temp directory per ORCA job to avoid MPI/IO collisions.
            job_tag = inp_file.stem.split('_')[-1]
            if not (len(job_tag) == 8 and all(c in '0123456789abcdef' for c in job_tag.lower())):
                job_tag = uuid.uuid4().hex[:8]
            temp_dir = Path(tempfile.mkdtemp(prefix=f"orca_run_{job_tag}_"))
            temp_out_file: Optional[Path] = None
            succeeded = False
            try:
                temp_inp_file = temp_dir / inp_file.name
                temp_out_file = temp_inp_file.with_suffix('.out')
                shutil.copy(inp_file, temp_inp_file)

                # Copy every locally referenced XYZ file into the sandbox.
                # NEB needs a start, end and optional TS-guess geometry.  Do
                # not preserve external paths here: the input renderer emits
                # basenames specifically so ORCA only ever sees sandbox-local
                # files, including when the original worktree has spaces.
                inp_text = temp_inp_file.read_text()
                referenced_xyz: List[str] = []
                for line in inp_text.splitlines():
                    stripped = line.strip()
                    if stripped.startswith("* xyzfile "):
                        parts = stripped.split()
                        if len(parts) >= 4:
                            referenced_xyz.append(parts[-1].strip('"'))
                        continue
                    neb_xyz_match = re.match(
                        r'NEB_(?:END|TS)_XYZFILE\s+"?([^"\s]+)"?',
                        stripped,
                        re.IGNORECASE,
                    )
                    if neb_xyz_match:
                        referenced_xyz.append(neb_xyz_match.group(1))
                for xyz_ref in referenced_xyz:
                    xyz_name = Path(xyz_ref).name
                    xyz_src = output_dir / xyz_name
                    if xyz_src.is_file():
                        shutil.copy(xyz_src, temp_dir / xyz_name)

                # A mode-directed rescue starts from the Hessian generated by
                # its local mode-analysis job.  Its source directory may have
                # spaces, so make that Hessian sandbox-local too and rewrite
                # the input reference just as we do for NEB XYZ endpoints.
                # This preserves ``InHess Read`` + ``TS_Mode`` on Windows/WSL
                # toxic paths instead of silently losing the selected mode.
                inp_text = temp_inp_file.read_text(encoding="utf-8", errors="ignore")
                hessian_refs = re.findall(
                    r'(?im)^\s*InHessName\s+"([^"]+)"\s*$', inp_text
                )
                for hessian_ref in hessian_refs:
                    hessian_src = Path(hessian_ref)
                    if not hessian_src.is_file():
                        hessian_src = output_dir / hessian_src.name
                    if not hessian_src.is_file():
                        continue
                    hessian_name = hessian_src.name
                    shutil.copy(hessian_src, temp_dir / hessian_name)
                    inp_text = inp_text.replace(
                        f'InHessName "{hessian_ref}"',
                        f'InHessName "{hessian_name}"',
                    )
                temp_inp_file.write_text(inp_text, encoding="utf-8")

                cmd = [str(self.orca_binary), str(temp_inp_file.resolve())]

                env = self._build_orca_runtime_env()
                env['ORCA_TEMP_DIR'] = str(temp_dir)

                metadata = self._execute_orca_process(
                    cmd=cmd,
                    cwd=temp_dir,
                    out_file=temp_out_file,
                    env=env,
                    timeout=timeout,
                    subprocess_callback=subprocess_callback,
                )
                self._last_orca_run_metadata = metadata

                if temp_out_file.exists():
                    shutil.copy(temp_out_file, out_file)
                else:
                    out_file.touch(exist_ok=True)

                for f in temp_dir.glob('*'):
                    if f.is_file() and f not in [temp_inp_file, temp_out_file]:
                        shutil.copy(f, output_dir / f.name)

                if metadata.timed_out or metadata.returncode not in (None, 0):
                    self._log_orca_failure_diagnostics(
                        out_file,
                        stderr_text=metadata.stderr_text,
                        returncode=metadata.returncode,
                        timed_out=metadata.timed_out,
                        sandbox_dir=temp_dir,
                    )
                    return out_file

                succeeded = True
                return out_file

            except Exception as e:
                # Preserve temp dir for postmortem; show tail of .out for faster diagnosis.
                try:
                    if temp_out_file is not None and temp_out_file.exists():
                        tail_lines = temp_out_file.read_text(errors='ignore').splitlines()[-50:]
                        self.logger.error(
                            "ORCA sandbox run failed; last 50 lines of output:\n" + "\n".join(tail_lines)
                        )
                except Exception as tail_e:
                    self.logger.warning(f"Failed to read ORCA sandbox output tail: {tail_e}")

                self.logger.error(f"ORCA sandbox dir preserved for debugging: {temp_dir}")
                raise RuntimeError(f"ORCA 运行出错 (sandbox): {e}")

            finally:
                # Only cleanup on success; keep failed sandbox directory for postmortem.
                if succeeded:
                    shutil.rmtree(temp_dir, ignore_errors=True)

        inp_file_abs = inp_file.resolve()
        cmd = [str(self.orca_binary), str(inp_file_abs)]

        try:
            metadata = self._execute_orca_process(
                cmd=cmd,
                cwd=output_dir,
                out_file=out_file,
                env=self._build_orca_runtime_env(),
                timeout=timeout,
                subprocess_callback=subprocess_callback,
            )
            self._last_orca_run_metadata = metadata
            if not out_file.exists():
                out_file.touch(exist_ok=True)

            if metadata.timed_out or metadata.returncode not in (None, 0):
                self._log_orca_failure_diagnostics(
                    out_file,
                    stderr_text=metadata.stderr_text,
                    returncode=metadata.returncode,
                    timed_out=metadata.timed_out,
                )

        except Exception as e:
            # Provide output tail for faster diagnosis
            try:
                if out_file.exists():
                    tail_lines = out_file.read_text(errors='ignore').splitlines()[-50:]
                    self.logger.error(
                        "ORCA run failed; last 50 lines of output:\n" + "\n".join(tail_lines)
                    )
            except Exception as tail_e:
                self.logger.warning(f"Failed to read ORCA output tail: {tail_e}")

            raise RuntimeError(f"ORCA 运行出错: {e}")

        return out_file

    def render_neb_input(
        self,
        reactant_xyz: Path,
        product_xyz: Path,
        output_dir: Path,
        *,
        charge: int = 0,
        spin: int = 1,
        n_images: int = 10,
        interpolation: str = "XTB2",
        method_keyword: str = "XTB2",
        job_basename: str = "neb_ts",
    ) -> Path:
        """
        渲染 ORCA NEB-TS 输入文件（默认 XTB2 级别）

        起止几何以确定性文件名复制到输出目录，使 NEB 产物
        （<basename>_NEB-TS_converged.xyz 等）路径可预测。

        Returns:
            生成的 .inp 文件路径
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        import shutil
        start_copy = output_dir / f"{job_basename}_start.xyz"
        end_copy = output_dir / f"{job_basename}_end.xyz"
        shutil.copy(Path(reactant_xyz), start_copy)
        shutil.copy(Path(product_xyz), end_copy)

        route = f"! {method_keyword} NEB-TS noautostart miniprint nopop"
        pal_block = self._pal_block()
        neb_block = (
            "%neb\n"
            f'  NEB_END_XYZFILE "{end_copy.name}"\n'
            f"  NImages {int(n_images)}\n"
            f"  Interpolation {interpolation}\n"
            "end\n"
        )
        inp_content = (
            f"{route}\n"
            f" %maxcore {self.maxcore}\n"
            f" {pal_block}"
            f"{neb_block}\n"
            f"* xyzfile {int(charge)} {int(spin)} {start_copy.name}\n"
        )
        inp_file = output_dir / f"{job_basename}.inp"
        inp_file.write_text(inp_content)
        return inp_file

    def run_neb_xtb2(
        self,
        reactant_xyz: Path,
        product_xyz: Path,
        output_dir: Path,
        *,
        charge: int = 0,
        spin: int = 1,
        n_images: int = 10,
        interpolation: str = "XTB2",
        job_basename: str = "neb_ts",
        timeout: Optional[int] = None,
        subprocess_callback: Callable[[subprocess.Popen[Any]], None] | None = None,
    ) -> Path:
        """
        执行 XTB2 NEB-TS 计算并返回 .out 文件路径

        仅负责输入渲染与执行；NEB 产物解析由 orca_neb_runner 完成。
        """
        inp_file = self.render_neb_input(
            reactant_xyz,
            product_xyz,
            output_dir,
            charge=charge,
            spin=spin,
            n_images=n_images,
            interpolation=interpolation,
            job_basename=job_basename,
        )
        return self._run_orca(
            inp_file,
            Path(output_dir),
            timeout=timeout,
            subprocess_callback=subprocess_callback,
        )

    def single_point(
        self,
        xyz_file: Path,
        output_dir: Path,
        timeout: Optional[int] = None,
        charge: Optional[int] = None,
        spin: Optional[int] = None,
        subprocess_callback: Callable[[subprocess.Popen[Any]], None] | None = None,
    ) -> QCResult:
        """
        执行单点能计算（端到端流程）

        整合流程:
        1. 生成 ORCA 输入文件 (_generate_input)
        2. 运行 ORCA 计算 (_run_orca)
        3. 解析输出文件 (_parse_output)

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录
            timeout: 超时时间（秒），None = 无限制

        Returns:
            QCResult 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        cache_key = self._compute_sp_cache_key(xyz_file, charge=charge, spin=spin)
        if cache_key and cache_key in self._sp_cache:
            self._sp_cache_hits += 1
            cached_energy = self._sp_cache[cache_key]
            self.logger.debug(f"SP cache hit: {xyz_file.name}")
            self._log_sp_cache_stats()
            return QCResult(
                energy=cached_energy,
                converged=True,
                output_file=None,
                error_message=None
            )

        if cache_key:
            self._sp_cache_misses += 1

        display_dir: object = output_dir
        try:
            display_dir = output_dir.resolve().relative_to(Path.cwd().resolve())
        except Exception:
            parts = output_dir.parts
            if len(parts) >= 3:
                display_dir = Path("...") / parts[-2] / parts[-1]
        self.logger.debug(f"开始 ORCA 单点能计算: {xyz_file.name}")
        self.logger.debug(f"  方法: {self.method}/{self.basis}")
        self.logger.debug(f"  输出目录: {display_dir}")

        try:
            # 步骤 1: 生成输入文件
            self.logger.debug("  生成 ORCA 输入文件...")
            inp_file = self._generate_input(
                xyz_file,
                output_dir,
                charge=charge,
                spin=spin,
                job_tag=uuid.uuid4().hex[:8]
            )
            self.logger.debug(f"  输入文件: {inp_file}")

            # 步骤 2: 运行 ORCA
            self.logger.debug("  运行 ORCA 计算...")
            out_file = self._run_orca(
                inp_file,
                output_dir,
                timeout=timeout,
                subprocess_callback=subprocess_callback,
            )
            self.logger.debug(f"  输出文件: {out_file}")

            # 步骤 3: 解析输出
            self.logger.debug("  解析 ORCA 输出...")
            result = self._parse_output(
                out_file,
                returncode=self._last_orca_run_metadata.returncode,
                stderr_text=self._last_orca_run_metadata.stderr_text,
                timed_out=self._last_orca_run_metadata.timed_out,
            )
            # Frequency jobs are routed through this entrypoint by qc_jobs.
            # _parse_output() owns common energy parsing, so explicitly
            # populate frequencies when the output contains a Hessian block.
            if result.converged and result.frequencies is None:
                result.frequencies = self._extract_frequencies_from_output(out_file)

            if result.converged:
                self.logger.debug(f"  ✓ 计算成功: 能量 = {result.energy:.8f} Hartree")
                if cache_key and result.energy is not None:
                    self._sp_cache[cache_key] = result.energy
                    self._log_sp_cache_stats()
            else:
                self.logger.error(f"  ✗ 计算失败: {result.error_message}")

            return result

        except Exception as e:
            self.logger.error(f"ORCA 单点能计算失败: {e}")
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=str(e)
            )

    def optimize(
        self,
        xyz_file: Path,
        output_dir: Path,
        route: Optional[str] = None,
        constraints: Optional[str] = None,
        old_checkpoint: Optional[Path] = None,
        timeout: Optional[int] = None,
        charge: int = 0,
        spin: int = 1
    ) -> QCResult:
        """
        ORCA 几何优化（统一接口，兼容 Gaussian 风格）

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录
            route: Gaussian 风格的 route card（会被转换为 ORCA 关键词）
            constraints: 约束（暂时忽略，ORCA 不支持相同格式）
            old_checkpoint: checkpoint 文件（暂时忽略）
            timeout: 超时时间（秒）
            charge: 分子电荷（默认 0）
            spin: 自旋多重度（默认 1，闭壳层）

        Returns:
            QCResult 对象
        """
        from rph_core.utils.optimization_config import OptimizationConfig

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 解析 route card，提取方法信息
        if route:
            # 简化处理：使用已配置的方法和基组
            # 实际项目中应该更精确地解析 route card
            if 'TS' in route:
                # TS 优化
                return self.ts_optimization(
                    xyz_file, output_dir,
                    opt_config=OptimizationConfig(),
                    timeout=timeout, charge=charge, spin=spin
                )
            else:
                # 基态优化 - 使用 Opt 关键词
                self.logger.info(f"ORCA 基态优化: {xyz_file.name}")
                return self._run_normal_optimization(
                    xyz_file, output_dir, timeout=timeout, charge=charge, spin=spin
                )

        # 默认基态优化
        return self._run_normal_optimization(
            xyz_file, output_dir, timeout=timeout, charge=charge, spin=spin
        )

    def _run_normal_optimization(
        self,
        xyz_file: Path,
        output_dir: Path,
        timeout: Optional[int] = None,
        charge: int = 0,
        spin: int = 1
    ) -> QCResult:
        """
        ORCA 基态优化（使用 Opt 关键词）

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录
            timeout: 超时时间
            charge: 分子电荷
            spin: 自旋多重度

        Returns:
            QCResult 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 生成优化输入文件
        import shutil
        job_tag = uuid.uuid4().hex[:8]
        base_name = f"{xyz_file.stem}_opt_{job_tag}"
        xyz_copy = output_dir / f"{base_name}{xyz_file.suffix}"
        shutil.copy(xyz_file, xyz_copy)

        # 读取 XYZ 内容
        xyz_lines = xyz_copy.read_text().split('\n')
        if len(xyz_lines) > 2:
            xyz_content = '\n'.join(xyz_lines[2:])
        else:
            xyz_content = xyz_copy.read_text()

        route = self._render_route(task_type="opt_freq")
        cpcm_block = self._render_cpcm_block()
        rendered_blocks = self._render_named_blocks(self.render_blocks())

        # 组装完整输入文件
        inp_content = f"""{route}
%maxcore {self.maxcore}
{self._pal_block()}
{cpcm_block}{rendered_blocks}
 * xyzfile {charge} {spin} {xyz_copy.name}
"""

        inp_file = output_dir / f"{base_name}.inp"
        inp_file.write_text(inp_content)

        # 运行 ORCA
        try:
            out_file = self._run_orca(inp_file, output_dir, timeout=timeout)
            result = self._parse_output(
                out_file,
                require_optimization_convergence=True,
                returncode=self._last_orca_run_metadata.returncode,
                stderr_text=self._last_orca_run_metadata.stderr_text,
                timed_out=self._last_orca_run_metadata.timed_out,
            )

            # 提取频率
            freq_block = re.search(r'VIBRATIONAL FREQUENCIES\s*\n((?:\s*[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+(?:\n|$))+)', out_file.read_text(), re.DOTALL)
            if freq_block:
                freqs = re.findall(r'[\-]?\d+\.\d+', freq_block.group(1))
                result.frequencies = np.array([float(f) for f in freqs]) if freqs else None

            # 提取最终优化坐标
            result.coordinates = self._extract_final_orca_coordinates(out_file, xyz_copy)

            return result

        except Exception as e:
            self.logger.error(f"ORCA 优化失败: {e}")
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=str(e)
            )

    def ts_optimization(
        self,
        xyz_file: Path,
        output_dir: Path,
        opt_config: Optional['OptimizationConfig'] = None,
        timeout: Optional[int] = None,
        charge: int = 0,
        spin: int = 1
    ) -> QCResult:
        """
        ORCA 过渡态优化

        使用 OptTS 关键词进行 TS 优化，支持 Hessian 控制

        Args:
            xyz_file: 输入 TS 猜想 XYZ 文件
            output_dir: 输出目录
            opt_config: 优化配置对象（来自 OptimizationConfig）
            timeout: 超时时间（秒），None = 无限制（QC 计算需要长时间）
            charge: 分子电荷
            spin: 自旋多重度

        Returns:
            QCResult 对象，包含优化后的坐标、能量和频率信息

        Raises:
            RuntimeError: ORCA 二进制文件未找到
            TimeoutError: 计算超时（如果设置了 timeout）
        """
        from rph_core.utils.optimization_config import OptimizationConfig

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"开始 ORCA TS 优化: {xyz_file.name}")
        self.logger.info(f"  方法: {self.method}/{self.basis}")
        self.logger.info(f"  电荷/自旋: {charge}/{spin}")

        # 如果没有提供 opt_config，使用默认值
        if opt_config is None:
            opt_config = OptimizationConfig()
        assert opt_config is not None

        # 检查超时设置
        if timeout is None:
            if opt_config.timeout_enabled and opt_config.timeout_seconds:
                timeout = opt_config.timeout_seconds
                self.logger.info(f"  超时: {timeout} 秒")
            else:
                timeout = None
                self.logger.info("  超时: 禁用（无限制）")

        try:
            # 生成 TS 优化输入文件
            self.logger.debug("  生成 ORCA TS 优化输入文件...")
            inp_file = self._generate_ts_input(
                xyz_file, output_dir, opt_config, charge, spin
            )
            self.logger.debug(f"  输入文件: {inp_file}")

            # 运行 ORCA（无超时或自定义超时）
            self.logger.debug("  运行 ORCA TS 优化...")
            out_file = self._run_orca(inp_file, output_dir, timeout=timeout)
            self.logger.debug(f"  输出文件: {out_file}")

            # 解析输出
            self.logger.debug("  解析 ORCA TS 优化输出...")
            result = self._parse_ts_output(out_file, xyz_file)

            if result.converged:
                self.logger.info(
                    f"  ✓ TS 优化成功: 能量 = {result.energy:.8f} Hartree"
                )
                if result.frequencies is not None:
                    imaginary = [f for f in result.frequencies if f < 0]
                    if imaginary:
                        self.logger.info(
                            f"  虚频: {len(imaginary)} 个, "
                            f"最小 = {min(imaginary):.1f} cm⁻¹"
                        )
            else:
                self.logger.error(f"  ✗ TS 优化失败: {result.error_message}")

            return result

        except Exception as e:
            self.logger.error(f"ORCA TS 优化失败: {e}")
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=str(e)
            )

    def _generate_ts_input(
        self,
        xyz_file: Path,
        output_dir: Path,
        opt_config: 'OptimizationConfig',
        charge: int,
        spin: int
    ) -> Path:
        """
        生成 ORCA TS 优化输入文件

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录
            opt_config: 优化配置
            charge: 分子电荷
            spin: 自旋多重度

        Returns:
            生成的 .inp 文件路径
        """
        import shutil

        route = self._render_route(task_type="ts_freq")
        cpcm_block = self._render_cpcm_block()
        rendered_blocks = self._render_named_blocks(self.render_blocks())
 
        # 如果 xyz_file 是 .out 文件，尝试从中提取电荷和自旋

        if xyz_file.suffix == '.out':
            try:
                charge, spin = CoordinateExtractor.get_charge_spin_from_orca_out(xyz_file)
                self.logger.debug(f"从 {xyz_file.name} 提取电荷/自旋: {charge}/{spin}")
            except Exception as e:
                self.logger.warning(f"无法从 {xyz_file.name} 提取电荷/自旋: {e}")

        # 读取 XYZ 内容
        job_tag = uuid.uuid4().hex[:8]
        base_name = f"{xyz_file.stem}_ts_opt_{job_tag}"
        xyz_copy = output_dir / f"{base_name}{xyz_file.suffix}"
        shutil.copy(xyz_file, xyz_copy)
        xyz_lines = xyz_copy.read_text().split('\n')
        if len(xyz_lines) > 2:
            xyz_content = '\n'.join(xyz_lines[2:])
        else:
            xyz_content = xyz_copy.read_text()

        # 生成 %geom 块
        geom_block = opt_config.to_orca_geom_block(is_ts=True)

        # 组装完整输入文件
        inp_content = f"""{route}
%maxcore {self.maxcore}
{self._pal_block()}
{cpcm_block}{rendered_blocks}
{geom_block}
 * xyzfile {charge} {spin} {xyz_copy.name}
"""

        # 写入文件
        inp_file = output_dir / f"{base_name}.inp"
        inp_file.write_text(inp_content)

        return inp_file

    def _run_orca_no_timeout(
        self,
        inp_file: Path,
        output_dir: Path
    ) -> Path:
        """
        运行 ORCA 计算（无超时限制）

        Args:
            inp_file: ORCA 输入文件
            output_dir: 输出目录

        Returns:
            ORCA 输出文件路径

        Raises:
            RuntimeError: ORCA 二进制文件未找到或运行失败
        """
        if self.orca_binary is None:
            raise RuntimeError(
                "ORCA 二进制文件未找到。请通过以下方式之一指定:\n"
                "  1. 设置环境变量 $ORCA_PATH\n"
                "  2. 确保 'orca' 在 PATH 中\n"
                "  3. 在初始化时提供 orca_binary_path 参数"
            )

        # 使用绝对路径确保 ORCA 能找到输入文件
        inp_file_abs = inp_file.resolve()
        out_file = inp_file.with_suffix('.out')
        cmd = [str(self.orca_binary), str(inp_file_abs)]
        env = self._build_orca_runtime_env()

        try:
            with open(out_file, 'w') as out_f:
                process = subprocess.Popen(
                    cmd,
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                    cwd=str(output_dir),
                    text=True,
                    env=env
                )

                # 等待进程完成（无超时限制）
                _, stderr = process.communicate()

                # 检查返回码
                if process.returncode != 0:
                    raise RuntimeError(
                        f"ORCA 运行失败 (返回码 {process.returncode})\n"
                        f"错误信息: {stderr}"
                    )

        except Exception as e:
            raise RuntimeError(f"ORCA 运行出错: {e}")

        return out_file

    def _parse_ts_output(self, out_file: Path, input_xyz: Optional[Path] = None) -> QCResult:
        """
        解析 ORCA TS 优化输出文件

        提取能量、收敛状态、坐标和频率信息

        Args:
            out_file: ORCA 输出文件路径

        Returns:
            QCResult 对象
        """
        import numpy as np

        # 读取输出文件内容
        try:
            content = out_file.read_text()
        except Exception as e:
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=f"无法读取输出文件: {e}"
            )

        # 检查是否正常终止
        if "ORCA TERMINATED NORMALLY" not in content:
            return QCResult(
                energy=0.0,
                converged=False,
                error_message="ORCA 未正常终止"
            )

        # 检查收敛状态
        converged = "THE OPTIMIZATION HAS CONVERGED" in content

        # 提取最终能量
        energy_match = re.search(
            r'FINAL\s+SINGLE\s+POINT\s+ENERGY\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)',
            content,
        )
        energy = float(energy_match.group(1)) if energy_match else 0.0

        # 提取频率 — 使用最后一个 VIBRATIONAL FREQUENCIES 块（最终收敛频率）
        # re.search() 会匹配第一个块（优化中期的中间频率），导致虚频数错误
        frequencies = self._extract_frequencies_from_output(out_file)

        # Resolve final geometry via sibling .xyz (preferred) or .out fallback
        final_xyz_path = self._resolve_orca_final_xyz(out_file, input_xyz)
        coordinates = None
        if final_xyz_path is not None:
            coordinates = self._read_xyz_coordinates(final_xyz_path)
        # Store the resolved xyz path so downstream job handling can
            # benchmark) can reference the canonical geometry file directly.
            self._last_resolved_final_xyz = final_xyz_path
        elif input_xyz is not None:
            coordinates = self._extract_orca_coordinates_from_out(out_file, input_xyz)
            self._last_resolved_final_xyz = None
        else:
            self.logger.warning("Cannot validate ORCA TS coordinates without input XYZ for %s", out_file)
            self._last_resolved_final_xyz = None

        error_message = None if converged else "优化未收敛"

        result = QCResult(
            energy=energy,
            converged=converged,
            coordinates=coordinates,
            frequencies=frequencies,
            output_file=out_file,
            error_message=error_message
        )
        # Attach resolved final xyz path to result for downstream canonical use
        if final_xyz_path is not None:
            result.final_xyz_path = str(final_xyz_path.resolve())
        return result
