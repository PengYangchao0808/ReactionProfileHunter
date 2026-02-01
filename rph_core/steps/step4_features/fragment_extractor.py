"""
Fragment Extractor - TS结构切分工具
======================================

从过渡态结构中提取独立的反应片段，用于畸变能计算

Author: QCcalc Team
Date: 2026-01-10
"""

import json
import logging
import os
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional, cast, Any, Dict
import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import read_xyz, write_xyz
from rph_core.utils.molecular_graph import build_bond_graph, get_connected_components
from rph_core.utils.fragment_manipulation import (
    get_fragment_charges,
    get_fragment_multiplicities
)

logger = logging.getLogger(__name__)


@dataclass
class FragmentSplitResult:
    fragment_indices: Tuple[List[int], List[int]]
    split_reason: str
    split_source: str
    fragA_charge: Optional[int] = None
    fragA_mult: Optional[int] = None
    fragB_charge: Optional[int] = None
    fragB_mult: Optional[int] = None
    debug: Optional[Dict[str, Any]] = None


class FragmentExtractor(LoggerMixin):
    """
    片段提取器 - 用于双片段畸变能计算

    功能:
    1. 从TS结构中切分出两个独立片段
    2. 生成片段XYZ文件
    3. 调用DFT进行单点能计算
    4. 支持BSSE校正（可选）
    """

    def __init__(self, config: Any, sp_engine: Any):
        """
        初始化片段提取器

        Args:
            config: 配置字典，包含DFT方法、基组等
        """
        self.config = config
        self.solvent = config.get('solvent', 'acetone')

        # DI: fragment single-point engine (must be ORCAInterface or compatible)
        self.sp_engine = sp_engine

        self.logger.info("FragmentExtractor 初始化: SP engine injected")

    def extract_and_calculate(
        self,
        ts_xyz: Path,
        fragment_indices: Optional[Tuple[List[int], List[int]]],
        output_dir: Path,
        old_checkpoint: Optional[Path] = None,
        apply_bsse: bool = False,
        forming_bonds: Optional[Tuple[Tuple[int, int], ...]] = None,
        reactant_xyz: Optional[Path] = None,
        split_config: Optional[Any] = None,
        system_charge: Optional[int] = None,
        system_mult: Optional[int] = None
    ) -> dict[str, Any]:
        """
        提取片段并计算DFT单点能

        Args:
            ts_xyz: TS结构XYZ文件
            fragment_indices: (fragment_A_atoms, fragment_B_atoms); if None, derive from forming_bonds
            output_dir: 输出目录
            old_checkpoint: 可选的checkpoint文件（复用TS的轨道）
            apply_bsse: 是否应用BSSE校正
            forming_bonds: 形成键索引 (用于自动分片)
            reactant_xyz: Reactant complex XYZ (用于 fragmenter 分片)
            split_config: 分片策略配置（use_fragmenter/fragmenter）

        Returns:
            {
                'e_fragment_a_ts': float,  # A片段在TS几何下的能量
                'e_fragment_b_ts': float,  # B片段在TS几何下的能量
                'e_fragment_a_relaxed': float,  # A片段松弛后的能量
                'e_fragment_b_relaxed': float,  # B片段松弛后的能量
            }
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("开始片段能量计算...")

        split_result = self._resolve_fragment_split(
            ts_xyz=ts_xyz,
            fragment_indices=fragment_indices,
            forming_bonds=forming_bonds,
            reactant_xyz=reactant_xyz,
            split_config=split_config,
            debug_dir=output_dir,
            system_charge=system_charge,
            system_mult=system_mult
        )
        frag_A_indices, frag_B_indices = split_result.fragment_indices

        # 读取TS结构
        ts_coords, ts_symbols = read_xyz(ts_xyz)

        if split_result.split_source == "intramolecular_fragmenter":
            return self._run_vrm_workflow(
                ts_xyz=ts_xyz,
                reactant_xyz=reactant_xyz,
                forming_bonds=forming_bonds,
                split_config=split_config,
                output_dir=output_dir,
                ts_coords=ts_coords,
                ts_symbols=ts_symbols,
                system_charge=system_charge,
                system_mult=system_mult
            )

        # 1. 生成片段结构（在TS几何下）
        frag_A_xyz = output_dir / "fragment_A_at_ts.xyz"
        frag_B_xyz = output_dir / "fragment_B_at_ts.xyz"

        self._extract_fragment_xyz(
            ts_coords, ts_symbols, frag_A_indices, frag_A_xyz,
            title="Fragment A (at TS geometry)"
        )
        self._extract_fragment_xyz(
            ts_coords, ts_symbols, frag_B_indices, frag_B_xyz,
            title="Fragment B (at TS geometry)"
        )

        # 2. 计算TS几何下的单点能（统一使用 ORCA；不复用 chk/fchk）
        fragA_charge, fragA_mult, fragB_charge, fragB_mult = self._resolve_fragment_charge_mult(
            split_result=split_result,
            split_config=split_config,
            system_charge=system_charge,
            system_mult=system_mult,
            output_dir=output_dir
        )

        e_frag_A_ts = self._single_point_orca(
            frag_A_xyz,
            output_dir / "fragment_A_ts_sp",
            charge=fragA_charge,
            mult=fragA_mult
        )

        e_frag_B_ts = self._single_point_orca(
            frag_B_xyz,
            output_dir / "fragment_B_ts_sp",
            charge=fragB_charge,
            mult=fragB_mult
        )

        # 3. 几何优化片段（得到松弛能量）
        frag_A_relaxed_xyz = output_dir / "fragment_A_relaxed.xyz"
        frag_B_relaxed_xyz = output_dir / "fragment_B_relaxed.xyz"

        # 使用快速方法进行片段优化（GFN-xTB或低级别DFT）
        e_frag_A_relaxed = self._optimize_fragment(
            frag_A_xyz, frag_A_relaxed_xyz, output_dir / "fragment_A_opt"
        )

        e_frag_B_relaxed = self._optimize_fragment(
            frag_B_xyz, frag_B_relaxed_xyz, output_dir / "fragment_B_opt"
        )

        # 4. 可选：BSSE校正
        if apply_bsse:
            self.logger.info("应用BSSE校正...")
            # TODO: 实现counterpoise校正
            # 这需要计算ghost原子体系

        results = {
            'e_fragment_a_ts': e_frag_A_ts,
            'e_fragment_b_ts': e_frag_B_ts,
            'e_fragment_a_relaxed': e_frag_A_relaxed,
            'e_fragment_b_relaxed': e_frag_B_relaxed,
            'fragment_split_reason': split_result.split_reason,
            'fragment_split_source': split_result.split_source,
            'fragment_indices': split_result.fragment_indices,
            'fragment_charge_mult': {
                'fragA': {'charge': fragA_charge, 'multiplicity': fragA_mult},
                'fragB': {'charge': fragB_charge, 'multiplicity': fragB_mult}
            }
        }

        self.logger.info("片段能量计算完成:")
        self.logger.info(f"  E(A_TS) = {e_frag_A_ts:.6f} Hartree")
        self.logger.info(f"  E(B_TS) = {e_frag_B_ts:.6f} Hartree")
        self.logger.info(f"  E(A_relaxed) = {e_frag_A_relaxed:.6f} Hartree")
        self.logger.info(f"  E(B_relaxed) = {e_frag_B_relaxed:.6f} Hartree")

        return results

    def _resolve_fragment_split(
        self,
        ts_xyz: Path,
        fragment_indices: Optional[Tuple[List[int], List[int]]],
        forming_bonds: Optional[Tuple[Tuple[int, int], ...]],
        reactant_xyz: Optional[Path],
        split_config: Optional[Any],
        debug_dir: Path,
        system_charge: Optional[int],
        system_mult: Optional[int]
    ) -> FragmentSplitResult:
        split_config_dict = dict(split_config) if split_config else {}
        allow_unsafe_index_split = bool(split_config_dict.get("allow_unsafe_index_split", False))
        fail_hard = split_config_dict.get("fail_hard_on_split_error", True)
        use_fragmenter = split_config_dict.get("use_fragmenter", False)

        if allow_unsafe_index_split:
            raise RuntimeError("Unsafe index slicing is disabled by default and must not be used")

        if reactant_xyz is None or not reactant_xyz.exists():
            raise RuntimeError("reactant_xyz is required for fragment splitting")

        if forming_bonds is None:
            raise RuntimeError("forming_bonds is required for fragment splitting")

        if len(forming_bonds) != 2:
            raise RuntimeError(f"forming_bonds must have 2 pairs for [5+2], got {len(forming_bonds)}")

        # A1: reactant_xyz and ts_xyz must match atom identities (count + symbol order)
        reactant_coords, reactant_symbols = read_xyz(reactant_xyz)
        ts_coords, ts_symbols = read_xyz(ts_xyz)
        self._validate_atom_identity(reactant_symbols, ts_symbols)

        # A3: build reactant bond graph (truth source)
        graph = build_bond_graph(
            reactant_coords,
            reactant_symbols,
            scale=float(split_config_dict.get("connectivity_scale", 1.25)),
            min_dist=float(split_config_dict.get("bond_min_dist_angstrom", 0.6))
        )
        components = get_connected_components(graph)
        components_sorted = [sorted(c) for c in components]

        # Strategy A: provided indices (explicit only)
        if fragment_indices is not None:
            candidate = FragmentSplitResult(
                fragment_indices=fragment_indices,
                split_reason="provided",
                split_source="provided_indices"
            )
            try:
                self._validate_split(graph, forming_bonds, candidate.fragment_indices)
                return candidate
            except Exception as e:
                self._write_split_debug_json(
                    output_dir=debug_dir,
                    split_source=candidate.split_source,
                    split_reason=candidate.split_reason,
                    forming_bonds=forming_bonds,
                    components=components_sorted,
                    error=str(e),
                    extra={"fragment_indices": candidate.fragment_indices}
                )
                raise

        # Strategy selector: intramolecular vs bimolecular
        is_intramolecular = len(components_sorted) == 1
        is_bimolecular = len(components_sorted) == 2

        if is_intramolecular:
            if not use_fragmenter:
                self.logger.warning(
                    "Intramolecular reactant detected with use_fragmenter disabled; forcing Strategy B"
                )
                split_config_dict["use_fragmenter"] = True
                use_fragmenter = True

            fragmenter_cfg = split_config_dict.get("fragmenter")
            if not isinstance(fragmenter_cfg, dict):
                fragmenter_cfg = {}
                split_config_dict["fragmenter"] = fragmenter_cfg

            if "charge_multiplicity" not in fragmenter_cfg:
                total_charge = int(system_charge) if system_charge is not None else 0
                total_mult = int(system_mult) if system_mult is not None else 1
                fragA_charge, fragB_charge = get_fragment_charges(
                    total_charge=total_charge,
                    n_fragA=0,
                    n_fragB=0,
                    dipole_in_fragA=True
                )
                fragA_mult, fragB_mult = get_fragment_multiplicities(
                    total_multiplicity=total_mult,
                    n_fragA=0,
                    n_fragB=0,
                    dipole_in_fragA=True
                )
                fragmenter_cfg["charge_multiplicity"] = {
                    "fragA": {"charge": fragA_charge, "multiplicity": fragA_mult},
                    "fragB": {"charge": fragB_charge, "multiplicity": fragB_mult}
                }
                split_config_dict["fragmenter"] = fragmenter_cfg

            from rph_core.steps.step3_opt.intramolecular_fragmenter import IntramolecularFragmenter

            fragmenter = IntramolecularFragmenter()
            fragmenter_cfg = split_config_dict.get("fragmenter", {})
            if not isinstance(fragmenter_cfg, dict):
                fragmenter_cfg = {}
            forming_bonds_pair = cast(Tuple[Tuple[int, int], Tuple[int, int]], forming_bonds)
            fragment_result = fragmenter.fragment(
                reactant_coords=reactant_coords,
                reactant_symbols=reactant_symbols,
                ts_coords=ts_coords,
                ts_symbols=ts_symbols,
                forming_bonds=forming_bonds_pair,
                config=fragmenter_cfg
            )

            if fragment_result.status != "ok":
                self._write_split_debug_json(
                    output_dir=debug_dir,
                    split_source="intramolecular_fragmenter",
                    split_reason="fragmenter_failed",
                    forming_bonds=forming_bonds,
                    components=components_sorted,
                    error=fragment_result.reason,
                    extra={
                        "cut_bond": getattr(fragment_result, "cut_bond_indices", None),
                        "dipole_core_path": getattr(fragment_result, "dipole_core_indices", None),
                        "alkene_endpoints": getattr(fragment_result, "alkene_end_indices", None),
                        "dipole_endpoints": getattr(fragment_result, "dipole_end_indices", None)
                    }
                )
                raise RuntimeError(f"Fragmenter failed (intramolecular): {fragment_result.reason}")

            candidate = FragmentSplitResult(
                fragment_indices=(fragment_result.fragA_indices, fragment_result.fragB_indices),
                split_reason="fragmenter",
                split_source="intramolecular_fragmenter",
                fragA_charge=getattr(fragment_result, "fragA_charge", None),
                fragA_mult=getattr(fragment_result, "fragA_mult", None),
                fragB_charge=getattr(fragment_result, "fragB_charge", None),
                fragB_mult=getattr(fragment_result, "fragB_mult", None),
                debug={
                    "cut_bond": fragment_result.cut_bond_indices,
                    "dipole_core_path": fragment_result.dipole_core_indices,
                    "alkene_endpoints": fragment_result.alkene_end_indices,
                    "dipole_endpoints": fragment_result.dipole_end_indices
                }
            )
            try:
                self._validate_split(graph, forming_bonds, candidate.fragment_indices)
                return candidate
            except Exception as e:
                self._write_split_debug_json(
                    output_dir=debug_dir,
                    split_source=candidate.split_source,
                    split_reason="validate_failed",
                    forming_bonds=forming_bonds,
                    components=components_sorted,
                    error=str(e),
                    extra=candidate.debug
                )
                raise

        if is_bimolecular:
            comp_a, comp_b = components_sorted
            # Stable ordering: by (size, min_index)
            comps_ordered = sorted([comp_a, comp_b], key=lambda c: (len(c), min(c)))
            frag_a_indices, frag_b_indices = comps_ordered[0], comps_ordered[1]
            candidate = FragmentSplitResult(
                fragment_indices=(frag_a_indices, frag_b_indices),
                split_reason="topology_components",
                split_source="topology_components"
            )
            try:
                self._validate_split(graph, forming_bonds, candidate.fragment_indices)
                return candidate
            except Exception as e:
                self._write_split_debug_json(
                    output_dir=debug_dir,
                    split_source=candidate.split_source,
                    split_reason="validate_failed",
                    forming_bonds=forming_bonds,
                    components=components_sorted,
                    error=str(e),
                    extra={"fragment_indices": candidate.fragment_indices}
                )
                if fail_hard:
                    raise

        self._write_split_debug_json(
            output_dir=debug_dir,
            split_source="strategy_selector",
            split_reason="unsupported_components",
            forming_bonds=forming_bonds,
            components=components_sorted,
            error=f"Unsupported reactant component count: {len(components_sorted)}"
        )
        raise RuntimeError(f"Unsupported reactant topology: {len(components_sorted)} components")

    def _validate_atom_identity(self, reactant_symbols: List[str], ts_symbols: List[str]) -> None:
        if len(reactant_symbols) != len(ts_symbols):
            raise RuntimeError(
                f"Atom count mismatch: reactant={len(reactant_symbols)} ts={len(ts_symbols)}"
            )
        if reactant_symbols != ts_symbols:
            # Provide a small diff without dumping everything
            mismatch = []
            for idx, (a, b) in enumerate(zip(reactant_symbols, ts_symbols)):
                if a != b:
                    mismatch.append((idx, a, b))
                    if len(mismatch) >= 10:
                        break
            raise RuntimeError(
                f"Element sequence mismatch between reactant and TS (first mismatches: {mismatch})"
            )

    def _validate_split(
        self,
        graph: Dict[int, List[int]],
        forming_bonds: Tuple[Tuple[int, int], ...],
        fragment_indices: Tuple[List[int], List[int]]
    ) -> None:
        frag_a, frag_b = fragment_indices
        set_a = set(frag_a)
        set_b = set(frag_b)

        n_atoms = len(graph)
        all_atoms = set(range(n_atoms))

        if set_a & set_b:
            raise RuntimeError("Fragment split invalid: A and B overlap")
        if (set_a | set_b) != all_atoms:
            missing = sorted(all_atoms - (set_a | set_b))[:20]
            extra = sorted((set_a | set_b) - all_atoms)[:20]
            raise RuntimeError(
                f"Fragment split invalid: coverage mismatch (missing={missing}, extra={extra})"
            )

        for (i, j) in forming_bonds:
            in_a = i in set_a
            in_b = i in set_b
            j_in_a = j in set_a
            j_in_b = j in set_b
            if not ((in_a and j_in_b) or (in_b and j_in_a)):
                raise RuntimeError(
                    f"Forming bond ({i},{j}) does not cross fragments"
                )

        if not self._is_connected_subgraph(graph, set_a):
            raise RuntimeError("Fragment A is not connected in reactant graph")
        if not self._is_connected_subgraph(graph, set_b):
            raise RuntimeError("Fragment B is not connected in reactant graph")

    def _is_connected_subgraph(self, graph: Dict[int, List[int]], nodes: set[int]) -> bool:
        if not nodes:
            return False
        start = next(iter(nodes))
        visited = set([start])
        queue = [start]
        while queue:
            cur = queue.pop()
            for nb in graph.get(cur, []):
                if nb in nodes and nb not in visited:
                    visited.add(nb)
                    queue.append(nb)
        return visited == nodes

    def _write_split_debug_json(
        self,
        output_dir: Path,
        split_source: str,
        split_reason: str,
        forming_bonds: Optional[Tuple[Tuple[int, int], ...]],
        components: Optional[List[List[int]]],
        error: str,
        extra: Optional[Dict[str, Any]] = None
    ) -> None:
        try:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            path = output_dir / "split_debug.json"
            payload: Dict[str, Any] = {
                "split_source": split_source,
                "split_reason": split_reason,
                "forming_bonds": [list(b) for b in forming_bonds] if forming_bonds else None,
                "components": components,
                "error": error
            }
            if extra:
                payload["extra"] = extra
            path.write_text(json.dumps(payload, indent=2, sort_keys=True))
        except Exception as e:
            self.logger.warning(f"Failed to write split_debug.json: {e}")

    def _resolve_fragment_charge_mult(
        self,
        split_result: FragmentSplitResult,
        split_config: Optional[Any],
        system_charge: Optional[int],
        system_mult: Optional[int],
        output_dir: Path
    ) -> Tuple[int, int, int, int]:
        if split_result.fragA_charge is not None and split_result.fragA_mult is not None and \
           split_result.fragB_charge is not None and split_result.fragB_mult is not None:
            return (
                int(split_result.fragA_charge),
                int(split_result.fragA_mult),
                int(split_result.fragB_charge),
                int(split_result.fragB_mult)
            )

        cfg = dict(split_config) if split_config else {}
        cm = cfg.get("charge_multiplicity", None)
        if cm is None:
            self._write_split_debug_json(
                output_dir=output_dir,
                split_source=split_result.split_source,
                split_reason="missing_charge_multiplicity",
                forming_bonds=None,
                components=None,
                error="Fragment charge/multiplicity missing; provide via fragmenter output or split_config.charge_multiplicity"
            )
            raise RuntimeError("Fragment charge/multiplicity missing (required)")

        try:
            fragA = cm["fragA"]
            fragB = cm["fragB"]
            return (
                int(fragA["charge"]),
                int(fragA["multiplicity"]),
                int(fragB["charge"]),
                int(fragB["multiplicity"])
            )
        except Exception:
            raise RuntimeError("Invalid split_config.charge_multiplicity format")

    def _extract_fragment_xyz(
        self,
        coords: np.ndarray,
        symbols: List[str],
        indices: List[int],
        output_xyz: Path,
        title: str = "Fragment"
    ):
        """
        从完整结构中提取片段XYZ

        Args:
            coords: 完整坐标 (N, 3)
            symbols: 原子符号列表
            indices: 片段原子索引列表
            output_xyz: 输出XYZ文件
            title: 标题
        """
        # 提取片段坐标和符号
        frag_coords = coords[indices]
        frag_symbols = [symbols[i] for i in indices]

        # 保存XYZ
        write_xyz(output_xyz, frag_coords, frag_symbols, title=title)

        self.logger.info(f"✓ 提取片段: {output_xyz} ({len(indices)} 原子)")

    def _single_point_orca(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int,
        mult: int
    ) -> float:
        """
        DFT单点能计算

        Args:
            xyz_file: 输入XYZ文件
            output_dir: 输出目录
            old_checkpoint: 复用的checkpoint

        Returns:
            能量 (Hartree)
        """
        try:
            result = self.sp_engine.single_point(
                xyz_file=xyz_file,
                output_dir=output_dir,
                charge=charge,
                spin=mult
            )
            if not result.converged:
                raise RuntimeError(f"ORCA fragment SP not converged: {result.error_message}")

            energy = result.energy
            if energy is None:
                raise RuntimeError("ORCA fragment SP returned no energy")

            try:
                geom_hash = __import__('hashlib').sha256(Path(xyz_file).read_bytes()).hexdigest()
            except Exception:
                geom_hash = None

            provenance = {
                "engine": "orca",
                "theory": {
                    "method": getattr(self.sp_engine, "method", None),
                    "basis": getattr(self.sp_engine, "basis", None),
                    "aux_basis": getattr(self.sp_engine, "aux_basis", None),
                    "solvent": getattr(self.sp_engine, "solvent", None)
                },
                "charge": int(charge),
                "multiplicity": int(mult),
                "geom_hash": geom_hash,
                "source_dir": str(output_dir),
                "output_file": str(getattr(result, "output_file", "")) if getattr(result, "output_file", None) else None
            }
            try:
                (Path(output_dir) / "provenance.json").write_text(
                    json.dumps(provenance, indent=2, sort_keys=True)
                )
            except Exception as e:
                self.logger.warning(f"Failed to write SP provenance.json: {e}")

            return float(energy)

        except Exception as e:
            self.logger.error(f"单点计算失败: {e}")
            raise RuntimeError(f"片段单点计算失败: {e}")

    def _run_vrm_workflow(
        self,
        ts_xyz: Path,
        reactant_xyz: Optional[Path],
        forming_bonds: Optional[Tuple[Tuple[int, int], ...]],
        split_config: Optional[Any],
        output_dir: Path,
        ts_coords: np.ndarray,
        ts_symbols: List[str],
        system_charge: Optional[int],
        system_mult: Optional[int]
    ) -> dict[str, Any]:
        if reactant_xyz is None or not reactant_xyz.exists():
            raise RuntimeError("reactant_xyz is required for VRM workflow")
        if forming_bonds is None:
            raise RuntimeError("forming_bonds is required for VRM workflow")

        reactant_coords, reactant_symbols = read_xyz(reactant_xyz)

        split_cfg = dict(split_config) if split_config else {}
        vrm_sp_nproc = split_cfg.get("vrm_sp_nproc")
        fragmenter_cfg = split_cfg.get("fragmenter", {})
        if not isinstance(fragmenter_cfg, dict):
            fragmenter_cfg = {}

        if "charge_multiplicity" not in fragmenter_cfg:
            fragmenter_cfg["charge_multiplicity"] = {
                "fragA": {"charge": 0, "multiplicity": 1},
                "fragB": {"charge": 0, "multiplicity": 1}
            }

        from rph_core.steps.step3_opt.intramolecular_fragmenter import IntramolecularFragmenter

        fragmenter = IntramolecularFragmenter()
        forming_bonds_pair = cast(Tuple[Tuple[int, int], Tuple[int, int]], forming_bonds)
        fragment_result = fragmenter.fragment(
            reactant_coords=reactant_coords,
            reactant_symbols=reactant_symbols,
            ts_coords=ts_coords,
            ts_symbols=ts_symbols,
            forming_bonds=forming_bonds_pair,
            config=fragmenter_cfg
        )

        if fragment_result.status != "ok":
            raise RuntimeError(f"Fragmenter failed (VRM): {fragment_result.reason}")

        coords_a = (
            fragment_result.fragA_coords_capped
            if fragment_result.fragA_coords_capped is not None
            else fragment_result.fragA_coords_TS
        )
        coords_b = (
            fragment_result.fragB_coords_capped
            if fragment_result.fragB_coords_capped is not None
            else fragment_result.fragB_coords_TS
        )
        syms_a = (
            fragment_result.fragA_symbols_capped
            if fragment_result.fragA_symbols_capped is not None
            else fragment_result.fragA_symbols_TS
        )
        syms_b = (
            fragment_result.fragB_symbols_capped
            if fragment_result.fragB_symbols_capped is not None
            else fragment_result.fragB_symbols_TS
        )

        if system_charge is None:
            system_charge = 0
        if system_mult is None:
            system_mult = 1

        # Houk 修正: 继承主线电荷，而非智能推导
        # VRM 模型电荷 = 系统电荷（添加 H 不改变电荷）
        # 片段电荷直接继承 fragmenter config 中的设定
        frag_a_charge = fragment_result.fragA_charge if fragment_result.fragA_charge is not None else 0
        frag_b_charge = fragment_result.fragB_charge if fragment_result.fragB_charge is not None else 0

        frag_a_mult = self._enforce_physical_spin(list(syms_a), frag_a_charge, 1)
        frag_b_mult = self._enforce_physical_spin(list(syms_b), frag_b_charge, 1)

        self.logger.info(
            f"H-capped fragment charges (inherited): FragA={frag_a_charge} (mult={frag_a_mult}), "
            f"FragB={frag_b_charge} (mult={frag_b_mult})"
        )

        n_orig_a = len(fragment_result.fragA_indices)
        n_orig_b = len(fragment_result.fragB_indices)

        combined_coords = np.vstack([coords_a, coords_b])
        combined_symbols = list(syms_a) + list(syms_b)

        frozen_indices = list(range(n_orig_a)) + list(
            range(len(syms_a), len(syms_a) + n_orig_b)
        )

        self.logger.info("Running constrained optimization for VRM H-caps")
        # P0-1: 强制独立沙盒 - 每次 VRM 运行使用唯一目录，避免文件冲突
        run_id = f"vrm_{uuid.uuid4().hex[:8]}"
        vrm_dir = Path(output_dir) / run_id
        vrm_dir.mkdir(parents=True, exist_ok=True)
        self.logger.debug(f"VRM sandbox created: {vrm_dir}")

        old_nprocs = getattr(self.sp_engine, "nprocs", None)
        try:
            if vrm_sp_nproc is not None:
                self.logger.info(f"VRM SP 固定 nprocs={vrm_sp_nproc}")
                setattr(self.sp_engine, "nprocs", int(vrm_sp_nproc))

            # Houk 修正: VRM 模型电荷直接继承系统电荷
            vrm_model_charge = system_charge
            vrm_model_mult = self._enforce_physical_spin(combined_symbols, vrm_model_charge, system_mult)
            self.logger.info(
                f"VRM combined model: charge={vrm_model_charge}, mult={vrm_model_mult}"
            )

            optimized_coords, optimized_symbols, optimized_xyz = self._run_constrained_optimization(
                combined_coords=combined_coords,
                combined_symbols=combined_symbols,
                frozen_indices=frozen_indices,
                work_dir=vrm_dir,
                system_charge=vrm_model_charge,  # 使用校正后的电荷
                system_mult=vrm_model_mult        # 使用校正后的多重度
            )

            n_a = len(syms_a)
            coords_a_opt = optimized_coords[:n_a]
            coords_b_opt = optimized_coords[n_a:]
            syms_a_opt = optimized_symbols[:n_a]
            syms_b_opt = optimized_symbols[n_a:]

            frag_a_vrm = vrm_dir / "fragment_A_vrm.xyz"
            frag_b_vrm = vrm_dir / "fragment_B_vrm.xyz"
            write_xyz(frag_a_vrm, coords_a_opt, syms_a_opt, title="Fragment A (VRM)")
            write_xyz(frag_b_vrm, coords_b_opt, syms_b_opt, title="Fragment B (VRM)")

            # 【简化】直接使用已计算的电荷/多重度，不再调用 _enforce_physical_spin
            # 因为 get_hcapped_fragment_charge 已确保闭壳层
            e_frag_a_ts = self._single_point_orca(
                frag_a_vrm,
                vrm_dir / "fragment_A_vrm_sp",
                charge=frag_a_charge,
                mult=frag_a_mult
            )
            e_frag_b_ts = self._single_point_orca(
                frag_b_vrm,
                vrm_dir / "fragment_B_vrm_sp",
                charge=frag_b_charge,
                mult=frag_b_mult
            )

            e_ts_model = self._single_point_orca(
                optimized_xyz,
                vrm_dir / "ts_model_vrm_sp",
                charge=vrm_model_charge,
                mult=vrm_model_mult
            )

            frag_a_relaxed = vrm_dir / "fragment_A_vrm_relaxed.xyz"
            frag_b_relaxed = vrm_dir / "fragment_B_vrm_relaxed.xyz"

            e_frag_a_relaxed = self._optimize_fragment(
                frag_a_vrm, frag_a_relaxed, vrm_dir / "fragment_A_vrm_opt"
            )
            e_frag_b_relaxed = self._optimize_fragment(
                frag_b_vrm, frag_b_relaxed, vrm_dir / "fragment_B_vrm_opt"
            )
        finally:
            if vrm_sp_nproc is not None and old_nprocs is not None:
                setattr(self.sp_engine, "nprocs", old_nprocs)

        results = {
            'e_fragment_a_ts': e_frag_a_ts,
            'e_fragment_b_ts': e_frag_b_ts,
            'e_fragment_a_relaxed': e_frag_a_relaxed,
            'e_fragment_b_relaxed': e_frag_b_relaxed,
            'e_ts_model_vrm': e_ts_model,
            'fragment_split_reason': 'vrm',
            'fragment_split_source': 'intramolecular_fragmenter',
            'fragment_indices': (fragment_result.fragA_indices, fragment_result.fragB_indices),
            'fragment_charge_mult': {
                'fragA': {'charge': frag_a_charge, 'multiplicity': frag_a_mult},
                'fragB': {'charge': frag_b_charge, 'multiplicity': frag_b_mult}
            }
        }

        self.logger.info("VRM fragment energy calculation completed")

        return results

    def _run_constrained_optimization(
        self,
        combined_coords: np.ndarray,
        combined_symbols: List[str],
        frozen_indices: List[int],
        work_dir: Path,
        system_charge: int,
        system_mult: int
    ) -> Tuple[np.ndarray, List[str], Path]:
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        model_xyz = work_dir / "model_ts_h_relax.xyz"
        write_xyz(model_xyz, combined_coords, combined_symbols, title="H-cap combined model")

        run_mult = self._enforce_physical_spin(combined_symbols, system_charge, system_mult)
        current_coords = combined_coords
        current_xyz = model_xyz

        self.logger.info("Step 1/3: xTB buffer optimization (消除 Steric Clash)")
        try:
            xtb_coords, xtb_symbols, xtb_xyz = self._run_xtb_preopt(
                model_xyz=model_xyz,
                frozen_indices=frozen_indices,
                work_dir=work_dir,
                charge=system_charge,
                mult=run_mult
            )
            current_coords = xtb_coords
            current_xyz = xtb_xyz
            self.logger.info("xTB preopt completed successfully")
        except Exception as e:
            self.logger.warning(f"xTB preopt failed: {e}; continuing with original geometry")

        self.logger.info("Step 2/3: Gaussian constrained optimization (Berny 精修)")
        try:
            from rph_core.utils.qc_interface import GaussianInterface
            gauss = GaussianInterface(
                charge=system_charge,
                multiplicity=run_mult,
                nprocshared=self.config.get('resources', {}).get('nproc', 8),
                mem=self.config.get('resources', {}).get('mem', '16GB'),
                config=self.config
            )
            gauss_result = gauss.constrained_optimize(
                xyz_file=current_xyz,
                output_dir=work_dir / "gaussian_opt",
                frozen_indices=frozen_indices,
                charge=system_charge,
                spin=run_mult
            )
            if gauss_result.converged and gauss_result.coordinates is not None:
                current_coords = cast(np.ndarray, gauss_result.coordinates)
                gauss_xyz = work_dir / "model_gaussian_opt.xyz"
                write_xyz(gauss_xyz, current_coords, combined_symbols, title="Gaussian opt")
                current_xyz = gauss_xyz
                self.logger.info("Gaussian constrained opt completed successfully")
            else:
                self.logger.warning("Gaussian opt failed; trying ORCA fallback")
        except Exception as e:
            self.logger.warning(f"Gaussian opt failed: {e}; trying ORCA fallback")

        self.logger.info("Step 3/3: ORCA constrained optimization (fallback/refinement)")
        try:
            constraints_lines = ["%geom", "  Constraints"]
            for idx in frozen_indices:
                constraints_lines.append(f"    {{ C {idx} C }}")
            constraints_lines.append("  end")
            constraints_lines.append("end")
            constraints_block = "\n".join(constraints_lines)

            result = self.sp_engine.constrained_optimize(
                xyz_file=current_xyz,
                output_dir=work_dir / "orca_opt",
                charge=system_charge,
                spin=run_mult,
                constraints_block=constraints_block
            )

            if result.converged and result.coordinates is not None:
                coords = cast(np.ndarray, result.coordinates)
                optimized_xyz = work_dir / "model_ts_h_relax_opt.xyz"
                write_xyz(optimized_xyz, coords, combined_symbols, title="H-cap relaxed model (opt)")
                return coords, list(combined_symbols), optimized_xyz
        except Exception as e:
            self.logger.warning(f"ORCA constrained opt failed: {e}")

        self.logger.warning("All constrained optimization attempts failed; using best available geometry")
        final_xyz = work_dir / "model_ts_h_relax_opt.xyz"
        write_xyz(final_xyz, current_coords, combined_symbols, title="H-cap relaxed model (best)")
        return current_coords, list(combined_symbols), final_xyz

    def _enforce_physical_spin(self, symbols: List[str], charge: int, mult: int) -> int:
        atomic_numbers = {
            'H': 1, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Si': 14,
            'P': 15, 'S': 16, 'Cl': 17, 'Br': 35, 'I': 53
        }

        total_electrons = 0
        for sym in symbols:
            z = atomic_numbers.get(sym)
            if z is None:
                raise RuntimeError(f"Unknown element for spin check: {sym}")
            total_electrons += z

        total_electrons -= int(charge)
        is_odd = total_electrons % 2 == 1
        mult_is_even = int(mult) % 2 == 0

        if is_odd and mult_is_even:
            return int(mult)
        if (not is_odd) and (not mult_is_even):
            return int(mult)

        corrected = 2 if is_odd else 1
        self.logger.warning(
            f"Multiplicity parity mismatch (Ne={total_electrons}, mult={mult}); "
            f"using {corrected}"
        )
        return corrected

    def _run_xtb_preopt(
        self,
        model_xyz: Path,
        frozen_indices: List[int],
        work_dir: Path,
        charge: int,
        mult: int
    ) -> Tuple[np.ndarray, List[str], Path]:
        from rph_core.utils.xtb_runner import XTBRunner

        work_dir = Path(work_dir)
        xtb_dir = work_dir / "xtb_preopt"
        xtb_dir.mkdir(parents=True, exist_ok=True)

        try:
            xtb = XTBRunner(self.config, work_dir=xtb_dir)
            uhf = 0 if mult == 1 else mult - 1
            result = xtb.optimize(
                structure=model_xyz,
                frozen_indices=frozen_indices,
                solvent=self.solvent,
                charge=charge,
                uhf=uhf
            )

            if result.success and result.output_file:
                opt_xyz = Path(result.output_file)
                coords, symbols = read_xyz(opt_xyz)
                self.logger.info(f"xTB preopt succeeded: {len(coords)} atoms")
                return coords, symbols, opt_xyz

        except Exception as e:
            self.logger.warning(f"xTB preopt failed: {e}; using original geometry")

        coords, symbols = read_xyz(model_xyz)
        return coords, symbols, model_xyz

    def _optimize_fragment(
        self,
        frag_xyz: Path,
        output_xyz: Path,
        output_dir: Path
    ) -> float:
        """
        片段几何优化（使用快速方法）

        策略:
        1. 首先尝试GFN2-xTB优化
        2. 可选：低级别DFT精修（B3LYP/def2-SVP）

        Args:
            frag_xyz: 片段XYZ文件
            output_xyz: 输出优化后的XYZ
            output_dir: 输出目录

        Returns:
            优化后的能量 (Hartree)
        """
        from rph_core.utils.qc_interface import XTBInterface

        # 使用XTB进行快速优化
        xtb = XTBInterface(
            gfn_level=2,
            solvent=self.solvent,
            nproc=8
        )

        self.logger.info(f"优化片段: {frag_xyz.name} (GFN2-xTB)")
        result = xtb.optimize(frag_xyz, output_dir)

        if result.converged and result.coordinates is not None and result.energy is not None:
            coordinates = cast(np.ndarray, result.coordinates)
            energy = cast(float, result.energy)
            # 保存优化后的结构
            write_xyz(
                output_xyz,
                coordinates,
                self._read_symbols(frag_xyz),
                title=f"{frag_xyz.stem} relaxed",
                energy=energy
            )
            return energy
        else:
            self.logger.warning(f"XTB优化未收敛，使用原始能量")
            # 读取原始能量
            coords, energy = self._read_xyz_with_energy(frag_xyz)
            return energy

    def _read_symbols(self, xyz_file: Path) -> List[str]:
        """从XYZ文件读取原子符号"""
        coords, symbols = read_xyz(xyz_file)
        return symbols

    def _read_xyz_with_energy(self, xyz_file: Path) -> Tuple[np.ndarray, float]:
        """从XYZ文件读取坐标和能量"""
        with open(xyz_file, 'r') as f:
            lines = f.readlines()

        n_atoms = int(lines[0].strip())
        title_line = lines[1].strip()

        # 尝试从标题提取能量
        energy = 0.0
        if "energy:" in title_line or "E =" in title_line:
            import re
            match = re.search(r'[\-]?\d+\.\d+', title_line)
            if match:
                energy = float(match.group())

        # 读取坐标
        coords = []
        for line in lines[2:2+n_atoms]:
            parts = line.strip().split()
            if len(parts) >= 4:
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

        return np.array(coords), energy
