"""
Clean Program Output Adapter
===========================

适配 Clean 程序输出的 CSV 文件，转换为标准化格式。

Clean 程序输出的 CSV 包含以下关键字段：
- reaction_type: 反应类型 (4+3, 5+2, etc.)
- precursor_smiles: 前体 SMILES
- product_smiles_main: 产物 SMILES
- precursor_type: 前体类型 (allenamide, etc.)
- topology: 拓扑类型 (INTER, INTRA_TYPE_I, etc.)
- cyclo_mode: 环加成模式 ([4+3], [5+2], etc.)
- core_atom_map: 核心原子映射 (JSON)
- core_bond_changes: 成键/断键变化 (如 "6-7:formed;10-11:formed")

Author: RPH Team
Date: 2026-03-18
"""

from __future__ import annotations

import json
import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class CleanRecord:
    """
    Clean 程序输出的标准化记录
    
    Attributes:
        reaction_id: 反应唯一 ID (rxn_key_hash)
        reaction_type: 反应类型
        precursor_smiles: 前体 SMILES
        product_smiles: 产物 SMILES
        precursor_type: 前体类型
        precursor_subtype: 前体子类型
        topology: 拓扑类型
        topology_evidence: 拓扑证据
        cyclo_mode: 环加成模式
        cyclo_mode_reason: 环加成证据
        core_atom_map: 核心原子映射 {idx: idx}
        core_bond_changes: 成键/断键变化
        new_ring_size: 新环大小
        raw: 原始数据
    """
    reaction_id: str
    reaction_type: str
    precursor_smiles: str
    product_smiles: str
    precursor_type: Optional[str] = None
    precursor_subtype: Optional[str] = None
    topology: str = "UNKNOWN"
    topology_evidence: str = ""
    cyclo_mode: str = "UNKNOWN"
    cyclo_mode_reason: str = ""
    core_atom_map: Optional[Dict[int, int]] = None
    core_bond_changes: Optional[Dict[str, List[Tuple[int, int]]]] = None
    new_ring_size: Optional[int] = None
    raw: Optional[Dict[str, str]] = None
    
    def __post_init__(self):
        if self.core_atom_map is None:
            self.core_atom_map = {}
        if self.core_bond_changes is None:
            self.core_bond_changes = {'forming': [], 'breaking': []}
        if self.raw is None:
            self.raw = {}


class CleanAdapter:
    """
    Clean 程序 CSV 输出适配器
    
    负责解析 Clean 程序输出的 CSV 文件，转换为标准化的 CleanRecord。
    """
    
    # CSV 字段名映射 (Clean 输出 -> 内部字段)
    FIELD_MAPPING = {
        'rxn_key_hash': 'reaction_id',
        'precursor_smiles': 'precursor_smiles',
        'product_smiles_main': 'product_smiles',
        'reaction_type': 'reaction_type',
        'reaction_type_confidence': 'reaction_type_confidence',
        'precursor_type': 'precursor_type',
        'precursor_subtype': 'precursor_subtype',
        'topology': 'topology',
        'topology_evidence': 'topology_evidence',
        'cyclo_mode': 'cyclo_mode',
        'cyclo_mode_reason': 'cyclo_mode_reason',
        'core_atom_map': 'core_atom_map',
        'core_bond_changes': 'core_bond_changes',
        'new_ring_size': 'new_ring_size',
    }
    
    def __init__(self):
        """初始化适配器"""
        pass
    
    def parse_csv(self, csv_path: Path) -> List[CleanRecord]:
        """
        解析 Clean 程序输出的 CSV 文件
        
        Args:
            csv_path: CSV 文件路径
            
        Returns:
            CleanRecord 列表
        """
        records = []
        
        if not csv_path.exists():
            logger.error(f"CSV file not found: {csv_path}")
            return records
        
        try:
            with open(csv_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                
                for row in reader:
                    try:
                        record = self._parse_row(row)
                        if record:
                            records.append(record)
                    except Exception as e:
                        logger.warning(f"Failed to parse row: {e}")
                        continue
                        
        except Exception as e:
            logger.error(f"Failed to read CSV: {e}")
            return []
        
        logger.info(f"Parsed {len(records)} records from {csv_path}")
        return records
    
    def parse_row(self, row: Dict[str, str]) -> Optional[CleanRecord]:
        """
        解析单行 CSV
        
        Args:
            row: CSV 行字典
            
        Returns:
            CleanRecord 或 None
        """
        try:
            return self._parse_row(row)
        except Exception as e:
            logger.warning(f"Failed to parse row: {e}")
            return None
    
    def _parse_row(self, row: Dict[str, str]) -> Optional[CleanRecord]:
        """内部解析方法"""
        # 提取基本字段
        reaction_id = row.get('rxn_key_hash', '')
        if not reaction_id:
            return None
        
        # 解析 core_atom_map (JSON 字符串)
        core_atom_map = self._parse_atom_map(row.get('core_atom_map', '{}'))
        
        # 解析 core_bond_changes，同时进行 MapId -> MolIdx 转换
        bond_changes = self._parse_bond_changes(
            row.get('core_bond_changes', ''),
            atom_map=core_atom_map
        )
        
        if not bond_changes.get('forming'):
            alt_forming = self._extract_forming_from_alt_sources(row)
            if alt_forming:
                bond_changes['forming'] = alt_forming
                logger.debug(f"[CleanAdapter] Using forming_bonds from alt sources: {alt_forming}")
        
        # 解析 new_ring_size
        new_ring_size = None
        if row.get('new_ring_size'):
            try:
                new_ring_size = int(row['new_ring_size'])
            except ValueError:
                pass
        
        # 构建原始数据字典
        raw = {k: v for k, v in row.items() if v}
        
        return CleanRecord(
            reaction_id=reaction_id,
            reaction_type=row.get('reaction_type', 'unknown'),
            precursor_smiles=row.get('precursor_smiles', ''),
            product_smiles=row.get('product_smiles_main', ''),
            precursor_type=row.get('precursor_type'),
            precursor_subtype=row.get('precursor_subtype'),
            topology=row.get('topology', 'UNKNOWN'),
            topology_evidence=row.get('topology_evidence', ''),
            cyclo_mode=row.get('cyclo_mode', 'UNKNOWN'),
            cyclo_mode_reason=row.get('cyclo_mode_reason', ''),
            core_atom_map=core_atom_map,
            core_bond_changes=bond_changes,
            new_ring_size=new_ring_size,
            raw=raw
        )
    
    def _parse_atom_map(self, atom_map_str: str) -> Dict[int, int]:
        """
        解析 core_atom_map 字段
        
        格式: JSON，如 '{"10": 10, "11": 11, ...}'
        """
        if not atom_map_str:
            return {}
        
        try:
            atom_map = json.loads(atom_map_str)
            # 转换为 int key
            return {int(k): v for k, v in atom_map.items()}
        except (json.JSONDecodeError, ValueError) as e:
            logger.warning(f"Failed to parse core_atom_map: {e}")
            return {}
    
    def _parse_bond_changes(
        self, 
        changes_str: str,
        atom_map: Optional[Dict[int, int]] = None
    ) -> Dict[str, List[Tuple[int, int]]]:
        """
        解析 core_bond_changes 字段
        
        格式: "6-7:formed;10-11:formed;6-9:broken"
        
        Args:
            changes_str: core_bond_changes 字符串
            atom_map: 可选的 MapId -> MolIdx 映射，用于坐标系转换
        """
        forming = []
        breaking = []
        
        if not changes_str:
            return {'forming': forming, 'breaking': breaking}
        
        # 准备 MapId -> MolIdx 转换
        # atom_map 格式: {MapId: MolIdx}
        # 需要反向映射: {MolIdx: MapId} 来将 MapId 转换为 MolIdx
        reverse_map = {}
        if atom_map:
            reverse_map = {v: k for k, v in atom_map.items()}
        
        for item in changes_str.split(';'):
            item = item.strip()
            if not item:
                continue
            
            # 解析 "6-7:formed" 格式
            if ':' not in item:
                continue
            
            try:
                bond_part, change_type = item.split(':', 1)
                change_type = change_type.strip()
                
                if '-' not in bond_part:
                    continue
                
                a, b = bond_part.split('-')
                a = int(a.strip())
                b = int(b.strip())
                
                # 如果有 atom_map，将 MapId 转换为 MolIdx
                # core_bond_changes 中的索引是 MapId
                # 需要找到对应的 MolIdx
                if atom_map and reverse_map:
                    # a 和 b 是 MapId，转换为 MolIdx
                    # atom_map 是 {MapId: MolIdx}
                    # 我们需要: 如果 atom_map[a] 存在，用它
                    a_molidx = atom_map.get(a, a)
                    b_molidx = atom_map.get(b, b)
                else:
                    a_molidx = a
                    b_molidx = b
                
                if change_type == 'formed':
                    forming.append((a_molidx, b_molidx))
                elif change_type == 'broken':
                    breaking.append((a_molidx, b_molidx))
                    
            except (ValueError, IndexError):
                logger.warning(f"Failed to parse bond change: {item}")
                continue
        
        return {'forming': forming, 'breaking': breaking}
    
    def _extract_forming_from_alt_sources(self, row: Dict[str, str]) -> List[Tuple[int, int]]:
        """从备用字段提取 forming_bonds (当 core_bond_changes 为空时)."""
        for key in ['formed_bond_index_pairs', 'forming_bonds', 'formed_bonds']:
            if key in row and row[key]:
                bonds = self._normalize_bond_list(row[key])
                if bonds:
                    return bonds
        return []
    
    def _normalize_bond_list(self, value: Any) -> List[Tuple[int, int]]:
        """归一化各种格式的 bond 列表."""
        if isinstance(value, list):
            result = []
            for pair in value:
                if isinstance(pair, (list, tuple)) and len(pair) == 2:
                    try:
                        result.append((int(pair[0]), int(pair[1])))
                    except (ValueError, TypeError):
                        continue
                elif isinstance(pair, str) and '-' in pair:
                    try:
                        a, b = pair.split('-', 1)
                        result.append((int(a.strip()), int(b.strip())))
                    except (ValueError, IndexError):
                        continue
            return result
        if isinstance(value, str):
            result = []
            for item in value.split(';'):
                item = item.strip()
                if not item:
                    continue
                for delim in ['-', ',', ' ']:
                    if delim in item:
                        try:
                            parts = item.split(delim)
                            if len(parts) >= 2:
                                a = int(parts[0].strip())
                                b = int(parts[1].strip())
                                result.append((a, b))
                                break
                        except (ValueError, IndexError):
                            continue
            return result
        return []
    
    def filter_by_reaction_type(
        self, 
        records: List[CleanRecord], 
        reaction_type: str
    ) -> List[CleanRecord]:
        """按反应类型过滤"""
        return [r for r in records if r.reaction_type == reaction_type]
    
    def filter_by_cyclo_mode(
        self, 
        records: List[CleanRecord], 
        cyclo_mode: str
    ) -> List[CleanRecord]:
        """按环加成模式过滤"""
        return [r for r in records if r.cyclo_mode == cyclo_mode]


class CleanAdapterV2(CleanAdapter):
    """
    增强版适配器，支持更多字段
    
    额外解析反应条件字段：
    - temperature, yield, de, ee, dr
    - solvent, base, leaving_group
    """
    
    EXTRA_FIELDS = [
        'yield', 'yield_individual', 'yield_overall',
        'de', 'ee', 'dr_major', 'dr_minor',
        'temp_celsius', 'temp_kelvin',
        'solvent', 'base', 'leaving_group',
        'mass_balance_status', 'stereo_consistency'
    ]
    
    def parse_with_conditions(self, row: Dict[str, str]) -> Dict[str, Any]:
        """解析包含反应条件的完整记录"""
        base = self._parse_row(row)
        if base is None:
            return {}
        
        # 提取条件字段
        result = {
            'reaction_id': base.reaction_id,
            'reaction_type': base.reaction_type,
            'cyclo_mode': base.cyclo_mode,
            'topology': base.topology,
            'precursor_smiles': base.precursor_smiles,
            'product_smiles': base.product_smiles,
            'core_atom_map': base.core_atom_map,
            'core_bond_changes': base.core_bond_changes,
            
            # 条件字段
            'conditions': {
                'temperature_celsius': self._parse_numeric(row.get('temp_celsius')),
                'temperature_kelvin': self._parse_numeric(row.get('temp_kelvin')),
                'yield': self._parse_numeric(row.get('yield')),
                'yield_individual': self._parse_numeric(row.get('yield_individual')),
                'de': self._parse_numeric(row.get('de')),
                'ee': self._parse_numeric(row.get('ee')),
                'dr_major': self._parse_numeric(row.get('dr_major')),
                'dr_minor': self._parse_numeric(row.get('dr_minor')),
                'solvent': row.get('solvent'),
                'base': row.get('base'),
                'leaving_group': row.get('leaving_group'),
                'mass_balance_status': row.get('mass_balance_status'),
                'stereo_consistency': row.get('stereo_consistency'),
            }
        }
        
        return result
    
    def _parse_numeric(self, value: Optional[str]) -> Optional[float]:
        """解析数值"""
        if not value:
            return None
        try:
            return float(value)
        except ValueError:
            return None