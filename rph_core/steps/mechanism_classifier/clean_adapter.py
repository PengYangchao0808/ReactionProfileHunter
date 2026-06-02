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
    map_confidence: Optional[float] = None
    low_confidence: bool = False
    raw: Optional[Dict[str, Any]] = None
    
    def __post_init__(self):
        if self.core_atom_map is None:
            self.core_atom_map = {}
        if self.core_bond_changes is None:
            self.core_bond_changes = {'forming': [], 'breaking': [], 'order_changed': []}
        elif 'order_changed' not in self.core_bond_changes:
            self.core_bond_changes['order_changed'] = []
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
    
    def parse_row(self, row: Dict[str, Any]) -> Optional[CleanRecord]:
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
    
    def _parse_row(self, row: Dict[str, Any]) -> Optional[CleanRecord]:
        """内部解析方法"""
        # 提取基本字段
        reaction_id = row.get('rxn_key_hash', '')
        if not reaction_id:
            return None

        raw_obj = row.get('raw')
        row_raw = raw_obj if isinstance(raw_obj, dict) else {}

        # 解析 core_atom_map (JSON 字符串)
        core_atom_map = self._parse_atom_map(row.get('core_atom_map', '{}'))
        
        # 解析 core_bond_changes，同时进行 MapId -> MolIdx 转换
        bond_changes = self._parse_bond_changes(
            row.get('core_bond_changes', ''),
            atom_map=core_atom_map
        )
        
        # 解析 new_ring_size
        new_ring_size = None
        if row.get('new_ring_size'):
            try:
                new_ring_size = int(row['new_ring_size'])
            except ValueError:
                pass

        map_confidence = self._parse_float(
            row.get('map_confidence')
            or row.get('map_score')
            or row.get('mapping_score')
            or row_raw.get('map_confidence')
            or row_raw.get('map_score')
            or row_raw.get('mapping_score')
        )
        low_confidence = self._parse_bool(row.get('low_confidence') or row_raw.get('low_confidence'))
        if map_confidence is not None and map_confidence < 0.8:
            low_confidence = True
            logger.warning(
                f"[CleanAdapter] reaction {reaction_id} map_confidence={map_confidence:.3f} below 0.8; "
                "retaining record as low-confidence"
            )

        # 构建原始数据字典
        raw = {k: v for k, v in row.items() if k != 'raw' and v not in (None, '')}
        raw.update({k: v for k, v in row_raw.items() if v not in (None, '')})
        if map_confidence is not None:
            raw['map_confidence'] = map_confidence
        if low_confidence:
            raw['low_confidence'] = True

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
            map_confidence=map_confidence,
            low_confidence=low_confidence,
            raw=raw
        )

    @staticmethod
    def _parse_float(value: Any) -> Optional[float]:
        if value is None:
            return None
        text = str(value).strip()
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            return None

    @staticmethod
    def _parse_bool(value: Any) -> bool:
        if isinstance(value, bool):
            return value
        if value is None:
            return False
        return str(value).strip().lower() in {'1', 'true', 'yes', 'y', 'on'}
    
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
        
        格式: "6-7:formed;10-11:formed;6-9:broken;12-19:order_changed(AROMATIC→SINGLE)"
        
        Args:
            changes_str: core_bond_changes 字符串
            atom_map: 可选的 MapId -> MolIdx 映射，用于坐标系转换
            
        Returns:
            dict with keys 'forming', 'breaking', 'order_changed'.
            'order_changed' stores (map_i, map_j) tuples for bonds whose
            order changed between precursor and product (e.g. AROMATIC→SINGLE).
        """
        forming: List[Tuple[int, int]] = []
        breaking: List[Tuple[int, int]] = []
        order_changed: List[Tuple[int, int]] = []
        
        if not changes_str:
            return {'forming': forming, 'breaking': breaking, 'order_changed': order_changed}
        
        for item in changes_str.split(';'):
            item = item.strip()
            if not item:
                continue
            
            # 解析 "6-7:formed" 或 "12-19:order_changed(AROMATIC→SINGLE)" 格式
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
                
                # core_bond_changes 中的索引是 MapId，直接保留
                if change_type == 'formed':
                    forming.append((a, b))
                elif change_type == 'broken':
                    breaking.append((a, b))
                elif change_type.startswith('order_changed'):
                    order_changed.append((a, b))
                    
            except (ValueError, IndexError):
                logger.warning(f"Failed to parse bond change: {item}")
                continue
        
        return {'forming': forming, 'breaking': breaking, 'order_changed': order_changed}
    
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
