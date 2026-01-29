"""
SMARTS Matcher - Product Bond Identification
==============================================

识别 [5+2] 环加成产物中的新形成 C-C 键
v2.1: 拓扑路径分析版 (Topological Path Analysis) - 修正了仅依赖键长的缺陷

Author: QCcalc Team / HY-Houk
Date: 2026-01-13
"""

import logging
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from rdkit import Chem
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

    def find_reactive_bonds(self, product_xyz: Path) -> SMARTSMatchResult:
        try:
            coords, symbols = read_xyz(product_xyz)
        except Exception as e:
            return SMARTSMatchResult(False, "io_error", None, None, (), 0.0, str(e))

        # 1. 构建分子图 (基于几何距离)
        mol = self._xyz_to_mol_with_connectivity(coords, symbols)
        if mol is None:
             return SMARTSMatchResult(False, "mol_build_error", None, None, (), 0.0, "无法构建分子连接性")

        # 2. 执行拓扑路径分析
        self.logger.info("执行拓扑路径分析 (Topological Path Analysis)...")
        topo_result = self._topological_core_identification(mol, coords)
        
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
            error_message="未能识别 [3.2.1] 氧杂桥环拓扑特征"
        )

    def _xyz_to_mol_with_connectivity(self, coords: np.ndarray, symbols: List[str]) -> Optional[Chem.Mol]:
        """将 XYZ 转换为带键的 RDKit Mol (保持不变)"""
        mol = Chem.RWMol()
        num_atoms = len(symbols)
        for s in symbols:
            mol.AddAtom(Chem.Atom(s))
        
        conf = Chem.Conformer(num_atoms)
        for i, coord in enumerate(coords):
            conf.SetAtomPosition(i, Point3D(float(coord[0]), float(coord[1]), float(coord[2])))
        mol.AddConformer(conf)

        dist_mat = np.zeros((num_atoms, num_atoms))
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                dist = np.linalg.norm(coords[i] - coords[j])
                dist_mat[i, j] = dist_mat[j, i] = dist
                s1, s2 = symbols[i], symbols[j]
                
                is_connected = False
                if s1 == 'H' or s2 == 'H':
                    if dist < 1.2: is_connected = True
                else:
                    if dist < 1.75: is_connected = True # 稍微放宽以容纳张力键
                
                if is_connected:
                    mol.AddBond(i, j, Chem.BondType.SINGLE)
        
        try:
            Chem.SanitizeMol(mol)
        except:
            pass
        return mol

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

    def _get_bond_info(self, mol, idx1, idx2, coords):
        """辅助函数：获取键信息"""
        length = float(np.linalg.norm(coords[idx1] - coords[idx2]))
        atom1 = mol.GetAtomWithIdx(idx1).GetSymbol()
        atom2 = mol.GetAtomWithIdx(idx2).GetSymbol()
        return BondInfo(idx1, idx2, length, "single", (atom1, atom2))