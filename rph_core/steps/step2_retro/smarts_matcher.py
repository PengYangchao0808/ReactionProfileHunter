"""
SMARTS Matcher - Product Bond Identification
==============================================

识别 [5+2] 环加成产物中的新形成 C-C 键
v2.1: 拓扑路径分析版 (Topological Path Analysis) - 修正了仅依赖键长的缺陷

Author: QCcalc Team / HY-Houk
Date: 2026-01-13
"""

import logging
from typing import Any, Dict, List, Optional, Tuple
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Geometry import Point3D

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import read_xyz

logger = logging.getLogger(__name__)


@dataclass
class BondInfo:
    atom_idx_1: int
    atom_idx_2: int
    current_length: float
    bond_type: str
    atom_types: Tuple[str, str]


@dataclass
class SMARTSMatchResult:
    matched: bool
    pattern_name: str
    bond_1: Optional[BondInfo]
    bond_2: Optional[BondInfo]
    match_atoms: Tuple[int, ...]
    confidence: float
    error_message: Optional[str] = None


@dataclass(frozen=True)
class SMARTSTemplate:
    reaction_type: str
    name: str
    smarts_pattern: Optional[str]
    bond_position_indices: Tuple[Tuple[int, int], Tuple[int, int]]
    core_atom_count: int
    identify_func: str


_TEMPLATES: Dict[str, SMARTSTemplate] = {
    "[5+2]": SMARTSTemplate(
        reaction_type="[5+2]",
        name="topology_path_321",
        smarts_pattern="[O]-[C]-[C]-[C]-[C]-[C]-[C]",
        bond_position_indices=((0, 1), (0, 2)),
        core_atom_count=5,
        identify_func="_topological_core_identification",
    ),
    "5+2": SMARTSTemplate(
        reaction_type="[5+2]",
        name="topology_path_321",
        smarts_pattern="[O]-[C]-[C]-[C]-[C]-[C]-[C]",
        bond_position_indices=((0, 1), (0, 2)),
        core_atom_count=5,
        identify_func="_topological_core_identification",
    ),
    "[4+3]": SMARTSTemplate(
        reaction_type="[4+3]",
        name="topology_ring_43",
        smarts_pattern="[O]1~[C]~[C]~[C]~[C]~[C]~[C]1",
        bond_position_indices=((0, 1), (0, 2)),
        core_atom_count=3,
        identify_func="_topological_ring_identification_43",
    ),
    "4+3": SMARTSTemplate(
        reaction_type="[4+3]",
        name="topology_ring_43",
        smarts_pattern="[O]1~[C]~[C]~[C]~[C]~[C]~[C]1",
        bond_position_indices=((0, 1), (0, 2)),
        core_atom_count=3,
        identify_func="_topological_ring_identification_43",
    ),
    "[4+2]": SMARTSTemplate(
        reaction_type="[4+2]",
        name="topology_ring_42",
        smarts_pattern="[C]1~[C]~[C]~[C]~[C]~[C]1",
        bond_position_indices=((0, 1), (0, 2)),
        core_atom_count=4,
        identify_func="_topological_ring_identification_42",
    ),
    "4+2": SMARTSTemplate(
        reaction_type="[4+2]",
        name="topology_ring_42",
        smarts_pattern="[C]1~[C]~[C]~[C]~[C]~[C]1",
        bond_position_indices=((0, 1), (0, 2)),
        core_atom_count=4,
        identify_func="_topological_ring_identification_42",
    ),
    "[3+2]": SMARTSTemplate(
        reaction_type="[3+2]",
        name="topology_ring_32",
        smarts_pattern="[#6,#7,#8]1~[#6,#7,#8]~[#6,#7,#8]~[#6,#7,#8]~[#6,#7,#8]1",
        bond_position_indices=((0, 1), (0, 2)),
        core_atom_count=3,
        identify_func="_topological_ring_identification_32",
    ),
    "3+2": SMARTSTemplate(
        reaction_type="[3+2]",
        name="topology_ring_32",
        smarts_pattern="[#6,#7,#8]1~[#6,#7,#8]~[#6,#7,#8]~[#6,#7,#8]~[#6,#7,#8]1",
        bond_position_indices=((0, 1), (0, 2)),
        core_atom_count=3,
        identify_func="_topological_ring_identification_32",
    ),
}


class SMARTSMatcher(LoggerMixin):
    """
    SMARTS 匹配引擎 (Step 2 组件)
    
    核心策略 v2.1:
    1. 几何感知构建分子图 (XYZ -> Mol)
    2. 锁定桥氧与桥头碳
    3. 路径分析: 区分 [3.2.1] 体系中的 "3碳桥" (旧键) 和 "2碳桥" (新键)
    """

    def __init__(self, config_path: Optional[Path] = None):
        self.logger.info("SMARTSMatcher (Topology Aware) 初始化完成")

    def find_reactive_bonds(
        self,
        product_xyz: Path,
        cleaner_data: Optional[Dict[str, Any]] = None,
    ) -> SMARTSMatchResult:
        try:
            coords, symbols = read_xyz(product_xyz)
        except Exception as e:
            return SMARTSMatchResult(False, "io_error", None, None, (), 0.0, str(e))

        # 1. 构建分子图 (基于几何距离)
        mol = self._xyz_to_mol_with_connectivity(coords, symbols)
        if mol is None:
             return SMARTSMatchResult(False, "mol_build_error", None, None, (), 0.0, "无法构建分子连接性")

        cleaner_result = self._match_from_cleaner_data(cleaner_data, mol, coords)
        if cleaner_result is not None and cleaner_result.matched:
            self.logger.info("✓ 使用 cleaner 数据识别形成键")
            return cleaner_result

        reaction_type = self._normalize_reaction_type(cleaner_data)
        template = self._select_template(reaction_type)
        if template is None:
            return SMARTSMatchResult(
                matched=False,
                pattern_name="template_not_found",
                bond_1=None,
                bond_2=None,
                match_atoms=(),
                confidence=0.0,
                error_message=f"未找到反应类型模板: {reaction_type}",
            )

        identify = getattr(self, template.identify_func, None)
        if not callable(identify):
            return SMARTSMatchResult(
                matched=False,
                pattern_name="template_error",
                bond_1=None,
                bond_2=None,
                match_atoms=(),
                confidence=0.0,
                error_message=f"模板识别函数不存在: {template.identify_func}",
            )

        self.logger.info("执行模板识别 (Template Registry)...")
        topo_candidate = identify(mol, coords)
        if not isinstance(topo_candidate, SMARTSMatchResult):
            return SMARTSMatchResult(
                matched=False,
                pattern_name="template_error",
                bond_1=None,
                bond_2=None,
                match_atoms=(),
                confidence=0.0,
                error_message=f"模板返回类型错误: {template.identify_func}",
            )
        topo_result = topo_candidate
        
        if topo_result.matched:
            self.logger.info(f"✓ 拓扑识别成功: {topo_result.pattern_name}")
            return topo_result

        return SMARTSMatchResult(
            matched=False, 
            pattern_name="failed", 
            bond_1=None, 
            bond_2=None, 
            match_atoms=(), 
            confidence=0.0, 
            error_message=f"未能识别 {template.reaction_type} 对应拓扑特征"
        )

    def _normalize_reaction_type(self, cleaner_data: Optional[Dict[str, Any]]) -> str:
        if not cleaner_data:
            return "[5+2]"

        candidate = cleaner_data.get("reaction_type") or cleaner_data.get("rxn_type")
        if not candidate:
            return "[5+2]"

        text = str(candidate).strip()
        return text if text else "[5+2]"

    def _select_template(self, reaction_type: str) -> Optional[SMARTSTemplate]:
        if reaction_type in _TEMPLATES:
            return _TEMPLATES[reaction_type]

        normalized = reaction_type.replace(" ", "")
        return _TEMPLATES.get(normalized)

    def _match_from_cleaner_data(
        self,
        cleaner_data: Optional[Dict[str, Any]],
        mol: Chem.Mol,
        coords: np.ndarray,
    ) -> Optional[SMARTSMatchResult]:
        if not cleaner_data:
            return None

        pairs = self._extract_cleaner_pairs(cleaner_data)
        if len(pairs) < 2:
            return None

        (a1, a2), (b1, b2) = pairs[0], pairs[1]
        atom_count = len(coords)
        if min(a1, a2, b1, b2) < 0 or max(a1, a2, b1, b2) >= atom_count:
            self.logger.warning("cleaner formed bond indices 超出范围，回退模板识别")
            return None

        return SMARTSMatchResult(
            matched=True,
            pattern_name="cleaner_formed_bond_indices",
            bond_1=self._get_bond_info(mol, a1, a2, coords),
            bond_2=self._get_bond_info(mol, b1, b2, coords),
            match_atoms=(a1, a2, b1, b2),
            confidence=1.0,
        )

    def _extract_cleaner_pairs(self, cleaner_data: Dict[str, Any]) -> List[Tuple[int, int]]:
        if "formed_bond_index_pairs" in cleaner_data:
            return self._parse_pair_payload(cleaner_data.get("formed_bond_index_pairs"))

        raw = cleaner_data.get("raw")
        if isinstance(raw, dict) and "formed_bond_index_pairs" in raw:
            return self._parse_pair_payload(raw.get("formed_bond_index_pairs"))

        return []

    def _parse_pair_payload(self, payload: Any) -> List[Tuple[int, int]]:
        if payload is None:
            return []

        parsed: List[Tuple[int, int]] = []

        if isinstance(payload, str):
            for chunk in payload.split(";"):
                pair = self._parse_pair_chunk(chunk)
                if pair is not None:
                    parsed.append(pair)
            return parsed

        if isinstance(payload, (list, tuple)):
            for item in payload:
                if isinstance(item, (list, tuple)) and len(item) == 2:
                    try:
                        parsed.append((int(item[0]), int(item[1])))
                    except (TypeError, ValueError):
                        continue
            return parsed

        return []

    def _parse_pair_chunk(self, chunk: str) -> Optional[Tuple[int, int]]:
        piece = chunk.strip()
        if not piece or "-" not in piece:
            return None

        left, right = piece.split("-", 1)
        try:
            return int(left.strip()), int(right.strip())
        except ValueError:
            return None

    def _xyz_to_mol_with_connectivity(self, coords: np.ndarray, symbols: List[str]) -> Optional[Chem.Mol]:
        mol = self._build_base_mol(coords, symbols)
        if mol is None:
            return None

        if self._try_determine_bonds(mol):
            return mol

        self.logger.warning("rdDetermineBonds 失败，回退到距离阈值连通性推断")
        return self._build_mol_with_distance_heuristic(coords, symbols)

    def _build_base_mol(self, coords: np.ndarray, symbols: List[str]) -> Optional[Chem.Mol]:
        mol = Chem.RWMol()
        num_atoms = len(symbols)
        for s in symbols:
            mol.AddAtom(Chem.Atom(s))

        conf = Chem.Conformer(num_atoms)
        for i, coord in enumerate(coords):
            conf.SetAtomPosition(i, Point3D(float(coord[0]), float(coord[1]), float(coord[2])))
        mol.AddConformer(conf)

        return mol.GetMol()

    def _try_determine_bonds(self, mol: Chem.Mol) -> bool:
        try:
            rdDetermineBonds.DetermineBonds(mol, charge=0)
        except Exception as exc:
            self.logger.debug(f"DetermineBonds failed: {exc}")
            try:
                rdDetermineBonds.DetermineConnectivity(mol)
            except Exception as connectivity_exc:
                self.logger.debug(f"DetermineConnectivity failed: {connectivity_exc}")
                return False

        if mol.GetNumBonds() == 0:
            self.logger.debug("rdDetermineBonds produced zero bonds")
            return False

        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:
            self.logger.debug(f"RDKit sanitize warning for inferred connectivity: {exc}")
        return True

    def _build_mol_with_distance_heuristic(self, coords: np.ndarray, symbols: List[str]) -> Optional[Chem.Mol]:
        mol = Chem.RWMol()
        num_atoms = len(symbols)
        for s in symbols:
            mol.AddAtom(Chem.Atom(s))

        conf = Chem.Conformer(num_atoms)
        for i, coord in enumerate(coords):
            conf.SetAtomPosition(i, Point3D(float(coord[0]), float(coord[1]), float(coord[2])))
        mol.AddConformer(conf)

        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                dist = np.linalg.norm(coords[i] - coords[j])
                s1, s2 = symbols[i], symbols[j]

                is_connected = False
                if s1 == 'H' or s2 == 'H':
                    if dist < 1.2:
                        is_connected = True
                else:
                    if dist < 1.75:
                        is_connected = True

                if is_connected:
                    mol.AddBond(i, j, Chem.BondType.SINGLE)

        result = mol.GetMol()
        try:
            Chem.SanitizeMol(result)
        except Exception as exc:
            self.logger.debug(f"RDKit sanitize warning for inferred connectivity: {exc}")
        return result

    def _topological_core_identification(self, mol: Chem.Mol, coords: np.ndarray) -> SMARTSMatchResult:
        """
        通过路径分析识别核心
        逻辑:
        1. 找桥氧 (连接且仅连接2个碳)
        2. 找桥头碳 (H1, H2)
        3. 寻找 H1 到 H2 的所有无环路径 (除去桥氧本身)
        4. 在 [5+2] 产物中，应该有一条 3-atom path (C-C-C) 和一条 2-atom path (C-C)
        5. 2-atom path 连接的就是新形成的键
        """
        
        # 1. 寻找桥氧和桥头碳
        bridge_oxygen_idx = -1
        bridgeheads = []

        for atom in mol.GetAtoms():
            if atom.GetSymbol() != 'O': continue
            
            carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']
            
            # 严格特征: 桥氧只连 2 个碳，且这两个碳之间通常没有直接键连 (在[3.2.1]中)
            if len(carbon_neighbors) == 2:
                # 排除羰基氧 (距离判断辅助)
                d1 = np.linalg.norm(coords[atom.GetIdx()] - coords[carbon_neighbors[0].GetIdx()])
                d2 = np.linalg.norm(coords[atom.GetIdx()] - coords[carbon_neighbors[1].GetIdx()])
                
                if d1 > 1.30 and d2 > 1.30:
                    bridge_oxygen_idx = atom.GetIdx()
                    bridgeheads = [n.GetIdx() for n in carbon_neighbors]
                    break
        
        if bridge_oxygen_idx == -1:
            return SMARTSMatchResult(False, "no_bridge_oxygen", None, None, (), 0.0)

        h1, h2 = bridgeheads
        self.logger.info(f"锁定桥氧: {bridge_oxygen_idx}, 桥头碳: {h1}, {h2}")

        # 2. 路径分析
        # 我们需要找到 H1 和 H2 之间的路径，但必须“屏蔽”掉桥氧，否则最短路径就是 H1-O-H2
        
        # 我们可以通过获取 H1 的邻居（除了O），和 H2 的邻居（除了O），看它们如何相连
        h1_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(h1).GetNeighbors() if n.GetIdx() != bridge_oxygen_idx]
        h2_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(h2).GetNeighbors() if n.GetIdx() != bridge_oxygen_idx]

        # 我们寻找两个独立的“桥”：
        # Bridge A (3-C): H1 - C - C - C - H2 (3 internal atoms? No, path length 4 bonds, 3 carbons inside)
        # Bridge B (2-C): H1 - C - C - H2 (path length 3 bonds, 2 carbons inside)
        # 实际上，[3.2.1] 中：
        # 3-atom bridge: H1 - C2 - C3 - C4 - H2 (3 atoms in between)
        # 2-atom bridge: H1 - C6 - C7 - H2 (2 atoms in between)
        
        # 使用 RDKit 的 GetShortestPath 可能会混淆，我们手动遍历 H1 的邻居出发的路径
        
        two_carbon_bridge_start_end = None # (start_neighbor_idx, end_neighbor_idx)
        
        # 遍历 H1 的所有非氧邻居
        for n1 in h1_neighbors:
            # 寻找从 n1 到 H2 的路径，且不经过 H1 和 O
            # 这是一个简单的寻路
            path = self._find_path_to_target(mol, start_idx=n1, target_idx=h2, 
                                           exclude_idxs={bridge_oxygen_idx, h1})
            
            if path:
                # path 是 [n1, ..., target]
                # 路径中间的原子数 = len(path) - 1 (target is h2)
                # 实际上我们关心的是这根桥上有几个碳
                # n1 是桥的第一个碳。
                # path 包含了直到 h2 之前的所有原子 + h2
                
                atoms_in_bridge = len(path) # 包括 n1 和 h2, 不包括 h1
                # Bridge atoms = atoms_in_bridge - 1 (subtract h2)
                bridge_atom_count = atoms_in_bridge - 1
                
                self.logger.info(f"发现路径: H1 -> {n1} ... -> H2, 桥原子数: {bridge_atom_count}")
                
                if bridge_atom_count == 2:
                    # 找到了 2-碳桥！这是来自烯烃的部分
                    # 形成键就是 H1-n1 和 (path倒数第二个原子)-H2
                    n_last = path[-2] # 连接 H2 的那个原子
                    two_carbon_bridge_start_end = (n1, n_last)
                    break
        
        if not two_carbon_bridge_start_end:
             return SMARTSMatchResult(False, "topology_mismatch", None, None, (), 0.0, "未找到 2-碳桥结构 ([3.2.1] 特征)")
        
        # 3. 提取键信息
        c_start, c_end = two_carbon_bridge_start_end
        
        bond_1 = self._get_bond_info(mol, h1, c_start, coords)
        bond_2 = self._get_bond_info(mol, h2, c_end, coords)
        
        return SMARTSMatchResult(
            matched=True,
            pattern_name="topology_path_321",
            bond_1=bond_1,
            bond_2=bond_2,
            match_atoms=tuple([bridge_oxygen_idx, h1, h2, c_start, c_end]),
            confidence=0.95
        )

    def _find_path_to_target(self, mol, start_idx, target_idx, exclude_idxs):
        """简单的 BFS 寻路，寻找从 start 到 target 的路径"""
        queue = [[start_idx]]
        visited = set(exclude_idxs)
        visited.add(start_idx)
        
        while queue:
            path = queue.pop(0)
            node = path[-1]
            
            if node == target_idx:
                return path
            
            for neighbor in mol.GetAtomWithIdx(node).GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx == target_idx:
                    return path + [n_idx]
                
                if n_idx not in visited and neighbor.GetSymbol() == 'C': # 桥上只应该是碳
                    visited.add(n_idx)
                    new_path = list(path)
                    new_path.append(n_idx)
                    queue.append(new_path)
        return None

    def _topological_ring_identification_43(self, mol: Chem.Mol, coords: np.ndarray) -> SMARTSMatchResult:
        return self._topological_ring_core_identification(
            mol=mol,
            coords=coords,
            ring_size=7,
            core_atom_count=3,
            pattern_name="topology_ring_43",
            require_ring_oxygen_count=1,
            require_oxygen_in_core=True,
        )

    def _topological_ring_identification_42(self, mol: Chem.Mol, coords: np.ndarray) -> SMARTSMatchResult:
        return self._topological_ring_core_identification(
            mol=mol,
            coords=coords,
            ring_size=6,
            core_atom_count=4,
            pattern_name="topology_ring_42",
        )

    def _topological_ring_identification_32(self, mol: Chem.Mol, coords: np.ndarray) -> SMARTSMatchResult:
        return self._topological_ring_core_identification(
            mol=mol,
            coords=coords,
            ring_size=5,
            core_atom_count=3,
            pattern_name="topology_ring_32",
        )

    def _topological_ring_core_identification(
        self,
        mol: Chem.Mol,
        coords: np.ndarray,
        ring_size: int,
        core_atom_count: int,
        pattern_name: str,
        require_ring_oxygen_count: Optional[int] = None,
        require_oxygen_in_core: bool = False,
    ) -> SMARTSMatchResult:
        ring_info = mol.GetRingInfo()
        atom_rings = [list(ring) for ring in ring_info.AtomRings() if len(ring) == ring_size]
        if not atom_rings:
            return SMARTSMatchResult(False, f"no_{ring_size}_membered_ring", None, None, (), 0.0)

        for ring in atom_rings:
            ring_symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
            if require_ring_oxygen_count is not None and ring_symbols.count("O") != require_ring_oxygen_count:
                continue

            cuts = self._find_cycloaddition_cut_edges(ring, core_atom_count, ring_symbols, require_oxygen_in_core)
            if cuts is None:
                continue

            (edge_i_a, edge_i_b), (edge_j_a, edge_j_b), core_atoms = cuts
            bond_1 = self._get_bond_info(mol, edge_i_a, edge_i_b, coords)
            bond_2 = self._get_bond_info(mol, edge_j_a, edge_j_b, coords)

            confidence = 0.90
            if require_oxygen_in_core:
                confidence = 0.93

            return SMARTSMatchResult(
                matched=True,
                pattern_name=pattern_name,
                bond_1=bond_1,
                bond_2=bond_2,
                match_atoms=tuple(core_atoms + [edge_i_a, edge_i_b, edge_j_a, edge_j_b]),
                confidence=confidence,
            )

        return SMARTSMatchResult(False, "topology_mismatch", None, None, (), 0.0)

    def _find_cycloaddition_cut_edges(
        self,
        ring: List[int],
        core_atom_count: int,
        ring_symbols: List[str],
        require_oxygen_in_core: bool,
    ) -> Optional[Tuple[Tuple[int, int], Tuple[int, int], List[int]]]:
        ring_size = len(ring)
        for i in range(ring_size):
            for j in range(i + 1, ring_size):
                if j == i + 1 or (i == 0 and j == ring_size - 1):
                    continue

                distance = j - i
                if distance not in (core_atom_count, ring_size - core_atom_count):
                    continue

                segment_a = [ring[(i + step) % ring_size] for step in range(1, distance + 1)]
                segment_b = [ring[(j + step) % ring_size] for step in range(1, ring_size - distance + 1)]

                core_atoms = segment_a if len(segment_a) == core_atom_count else segment_b
                if require_oxygen_in_core:
                    oxygen_in_core = any(
                        ring_symbols[ring.index(atom_idx)] == "O" for atom_idx in core_atoms
                    )
                    if not oxygen_in_core:
                        continue

                edge_1 = (ring[i], ring[(i + 1) % ring_size])
                edge_2 = (ring[j], ring[(j + 1) % ring_size])
                return edge_1, edge_2, core_atoms

        return None

    def _get_bond_info(self, mol, idx1, idx2, coords):
        """辅助函数：获取键信息"""
        length = float(np.linalg.norm(coords[idx1] - coords[idx2]))
        atom1 = mol.GetAtomWithIdx(idx1).GetSymbol()
        atom2 = mol.GetAtomWithIdx(idx2).GetSymbol()
        return BondInfo(idx1, idx2, length, "single", (atom1, atom2))
