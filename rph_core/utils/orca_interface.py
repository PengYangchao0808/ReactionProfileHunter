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

from pathlib import Path
from dataclasses import dataclass
from typing import Mapping, Optional, Union, TYPE_CHECKING, Dict, Any, List, Tuple
import uuid
import hashlib
import json
import re
import subprocess
import shutil
import os
import logging
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
from rph_core.utils.capability_validator import CapabilityValidator
from rph_core.utils.orca_input_renderer import OrcaInputRenderer

logger = logging.getLogger(__name__)


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
        self.route_extras = self.method_spec.route_extras

        # 处理 maxcore：如果未提供，从配置派生
        if maxcore is None and config:
            res_cfg = config.get('resources', {})
            mem = res_cfg.get('mem', '32GB')
            safety_factor = res_cfg.get('orca_maxcore_safety', 0.65)
            self.maxcore = calc_orca_maxcore(mem, nprocs, safety_factor)
        else:
            self.maxcore = maxcore if maxcore is not None else 4000

        validation_errors = CapabilityValidator.validate(self.method_spec, task_type="sp")
        if validation_errors:
            self.logger = logging.getLogger(f"{__name__}.{method}/{basis}")
            self.logger.warning("ORCA method spec validation issues: %s", "; ".join(validation_errors))

        # 查找 ORCA 二进制文件（集成新的配置系统）
        self.orca_binary = self._find_orca_binary(orca_binary_path, config)

        # 设置日志
        self.logger = logging.getLogger(f"{__name__}.{method}/{basis}")
        self.config = config

        # SP cache (in-memory, per-process)
        self._sp_cache = {}
        self._sp_cache_hits = 0
        self._sp_cache_misses = 0
        self._last_resolved_final_xyz: Optional[Path] = None

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
        return self._get_renderer().render_blocks()

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
        if self.solvent_model in {"CPCM", "PCM"}:
            solvent_keyword = f"CPCM({self._solvent_name_for_orca()})"
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
            result = self._parse_output(out_file, require_optimization_convergence=True)

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
        timeout: Optional[int] = None,
        job_tag: Optional[str] = None,
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
            out_file = self._run_orca(inp_file, output_dir, timeout=timeout)
            result = self._parse_output(out_file, require_optimization_convergence=True)
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

    def _parse_output(
        self,
        out_file: Path,
        *,
        require_optimization_convergence: bool = False,
    ) -> QCResult:
        """
        解析 ORCA 输出文件

        Args:
            out_file: ORCA 输出文件路径

        Returns:
            QCResult 对象
        """
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

        # 提取最终单点能
        energy_match = re.search(r"FINAL SINGLE POINT ENERGY\s+([\-\d\.]+)", content)
        if not energy_match:
            return QCResult(
                energy=0.0,
                converged=False,
                error_message="无法找到能量信息"
            )

        # 解析能量
        try:
            energy = float(energy_match.group(1))
        except ValueError:
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=f"能量格式错误: {energy_match.group(1)}"
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

        # 成功解析
        return QCResult(
            energy=energy,
            converged=(
                not require_optimization_convergence
                or "THE OPTIMIZATION HAS CONVERGED" in content
            ),
            coordinates=coordinates,
            output_file=out_file,
            error_message=(
                None
                if not require_optimization_convergence
                or "THE OPTIMIZATION HAS CONVERGED" in content
                else "ORCA terminated normally, but the geometry optimization did not converge"
            )
        )

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

        if os.getuid() == 0:
            env.setdefault("OMPI_ALLOW_RUN_AS_ROOT", "1")
            env.setdefault("OMPI_ALLOW_RUN_AS_ROOT_CONFIRM", "1")

        return env

    def _run_orca(
        self,
        inp_file: Path,
        output_dir: Path,
        timeout: Optional[int] = 3600
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
            RuntimeError: ORCA 运行失败
            TimeoutError: 计算超时
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

                # Copy referenced XYZ files into sandbox temp dir
                inp_text = temp_inp_file.read_text()
                for line in inp_text.splitlines():
                    stripped = line.strip()
                    if stripped.startswith("* xyzfile "):
                        parts = stripped.split()
                        if len(parts) >= 4:
                            xyz_ref = parts[-1]
                            xyz_src = output_dir / xyz_ref
                            if xyz_src.is_file():
                                shutil.copy(xyz_src, temp_dir / xyz_ref)

                cmd = [str(self.orca_binary), str(temp_inp_file.resolve())]

                env = self._build_orca_runtime_env()
                env['ORCA_TEMP_DIR'] = str(temp_dir)

                with open(temp_out_file, 'w') as out_f:
                    process = subprocess.Popen(
                        cmd,
                        stdout=out_f,
                        stderr=subprocess.PIPE,
                        cwd=str(temp_dir),
                        text=True,
                        env=env
                    )

                    try:
                        _, stderr = process.communicate(timeout=timeout)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        raise TimeoutError(
                            f"ORCA 计算超时 (>{timeout}秒): {inp_file.name}"
                        )

                    if process.returncode != 0:
                        if stderr.strip():
                            self.logger.error("ORCA stderr:\n%s", stderr[-4000:])
                        raise RuntimeError(
                            f"ORCA 运行失败 (返回码 {process.returncode})\n"
                            f"错误信息: {stderr}"
                        )

                shutil.copy(temp_out_file, out_file)
                succeeded = True

                for f in temp_dir.glob('*'):
                    if f.is_file() and f not in [temp_inp_file, temp_out_file]:
                        shutil.copy(f, output_dir / f.name)

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
            with open(out_file, 'w') as out_f:
                process = subprocess.Popen(
                    cmd,
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                    cwd=str(output_dir),
                    text=True,
                    env=self._build_orca_runtime_env()
                )

                # 等待进程完成或超时
                try:
                    _, stderr = process.communicate(timeout=timeout)
                except subprocess.TimeoutExpired:
                    process.kill()
                    raise TimeoutError(
                        f"ORCA 计算超时 (>{timeout}秒): {inp_file.name}"
                    )

                # 检查返回码
                if process.returncode != 0:
                    if stderr.strip():
                        self.logger.error("ORCA stderr:\n%s", stderr[-4000:])
                    raise RuntimeError(
                        f"ORCA 运行失败 (返回码 {process.returncode})\n"
                        f"错误信息: {stderr}"
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

    def single_point(
        self,
        xyz_file: Path,
        output_dir: Path,
        timeout: Optional[int] = None,
        charge: Optional[int] = None,
        spin: Optional[int] = None
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
            out_file = self._run_orca(inp_file, output_dir, timeout=timeout)
            self.logger.debug(f"  输出文件: {out_file}")

            # 步骤 3: 解析输出
            self.logger.debug("  解析 ORCA 输出...")
            result = self._parse_output(out_file)
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
            result = self._parse_output(out_file, require_optimization_convergence=True)

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
        energy_match = re.search(r'FINAL SINGLE POINT ENERGY\s+([\-\d\.]+)', content)
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
