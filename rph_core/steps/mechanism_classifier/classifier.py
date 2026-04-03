"""
Mechanism Classifier - S0 Stage Main Entry
==========================================

S0 机理分类模块主入口。

功能：
1. 解析 Clean 程序输出 CSV
2. 构建 DAG 图模型
3. 提供统一的分类接口
4. 输出边属性给 S2

Author: RPH Team
Date: 2026-03-18
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Dict, Optional, Any

from rph_core.steps.mechanism_classifier.models import (
    MechanismGraph,
    GraphEdge,
    CycloMode,
)
from rph_core.steps.mechanism_classifier.clean_adapter import (
    CleanAdapter,
    CleanRecord,
)
from rph_core.steps.mechanism_classifier.graph_builder import (
    GraphBuilder,
)

logger = logging.getLogger(__name__)


class MechanismClassifier:
    """
    S0 机理分类器主入口
    
    整合 Clean 适配器、图构建器，提供统一的分类接口。
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        初始化分类器
        
        Args:
            config: 配置字典，支持以下键：
                - generate_intermediates: 是否生成中间体 SMILES
                - enable_dr_paths: 是否启用 dr 路径枚举
        """
        self.config = config or {}
        
        # 初始化组件
        self.adapter = CleanAdapter()
        self.builder = GraphBuilder()
        
        # 配置选项
        self.generate_intermediates = self.config.get(
            "generate_intermediates", False
        )
        self.enable_dr_paths = self.config.get(
            "enable_dr_paths", False
        )
    
    def classify_from_csv(
        self,
        csv_path: Path,
        reaction_ids: Optional[List[str]] = None,
        reaction_type: Optional[str] = None,
    ) -> List[MechanismGraph]:
        """
        从 CSV 文件分类
        
        Args:
            csv_path: Clean 程序输出的 CSV 文件路径
            reaction_ids: 可选，要处理的具体反应 ID 列表
            reaction_type: 可选，按反应类型过滤 (如 "4+3")
            
        Returns:
            MechanismGraph 列表
        """
        logger.info(f"Loading records from {csv_path}")
        
        # 1. 解析 CSV
        records = self.adapter.parse_csv(csv_path)
        
        if not records:
            logger.warning(f"No records parsed from {csv_path}")
            return []
        
        # 2. 过滤 (如果指定)
        if reaction_ids:
            records = [r for r in records if r.reaction_id in reaction_ids]
            logger.info(f"Filtered to {len(records)} records by reaction_ids")
        
        if reaction_type:
            records = [r for r in records if r.reaction_type == reaction_type]
            logger.info(f"Filtered to {len(records)} records by reaction_type={reaction_type}")
        
        # 3. 构建图
        graphs = self.builder.build_batch(records)
        
        # 4. 可选：添加 dr 路径
        if self.enable_dr_paths:
            graphs = self._add_dr_paths(graphs)
        
        logger.info(f"Classified {len(graphs)} reactions")
        return graphs
    
    def classify_single(
        self,
        record: CleanRecord
    ) -> MechanismGraph:
        """
        单条记录分类
        
        Args:
            record: CleanRecord
            
        Returns:
            MechanismGraph
        """
        return self.builder.build(record)
    
    def classify_by_reaction_id(
        self,
        csv_path: Path,
        reaction_id: str
    ) -> Optional[MechanismGraph]:
        """
        按 reaction_id 分类单条反应
        
        Args:
            csv_path: Clean CSV 路径
            reaction_id: 要分类的反应 ID
            
        Returns:
            MechanismGraph 或 None（如果未找到）
        """
        records = self.adapter.parse_csv(csv_path)
        
        for record in records:
            if record.reaction_id == reaction_id:
                return self.builder.build(record)
        
        logger.warning(f"Reaction {reaction_id} not found in {csv_path}")
        return None
    
    def classify_from_dict(
        self,
        data: Dict[str, str]
    ) -> Optional[MechanismGraph]:
        """
        从字典分类
        
        Args:
            data: CSV 行字典
            
        Returns:
            MechanismGraph 或 None
        """
        record = self.adapter.parse_row(data)
        if record is None:
            return None
        
        return self.builder.build(record)
    
    def _add_dr_paths(
        self,
        graphs: List[MechanismGraph]
    ) -> List[MechanismGraph]:
        """
        为图添加 dr (diastereomeric ratio) 路径
        
        简单的 placeholder 实现。
        后续可以基于立体化学信息扩展。
        
        Args:
            graphs: MechanismGraph 列表
            
        Returns:
            添加 dr 路径后的图列表
        """
        # 当前实现：仅为 [4+3] 反应添加 endo/exo 路径
        # 这是一个简化版本
        
        result = []
        
        for graph in graphs:
            if graph.cyclo_mode == CycloMode.C4_PLUS_3:
                # TODO: 基于立体化学数据添加 endo/exo 路径
                # 当前跳过，因为需要实际的立体化学数据
                result.append(graph)
            else:
                result.append(graph)
        
        return result
    
    def get_edges_for_pathway(
        self, 
        graph: MechanismGraph, 
        pathway: str = "primary"
    ) -> List[Dict[str, Any]]:
        """
        获取指定路径的边属性
        
        供 S2 模块使用
        
        Args:
            graph: MechanismGraph
            pathway: 路径 ID
            
        Returns:
            边属性字典列表
        """
        edges = graph.get_edges_for_pathway(pathway)
        
        return [
            {
                "source": edge.source,
                "target": edge.target,
                "forming_bonds": edge.forming_bonds,
                "breaking_bonds": edge.breaking_bonds,
                "ts_type": edge.ts_type,
                "pathway_id": edge.pathway_id
            }
            for edge in edges
        ]
    
    def export_to_json(
        self,
        graphs: List[MechanismGraph],
        output_path: Path
    ) -> None:
        """
        导出到 JSON 文件
        
        Args:
            graphs: MechanismGraph 列表
            output_path: 输出路径
        """
        data = []
        
        for graph in graphs:
            cyclo = graph.cyclo_mode if isinstance(graph.cyclo_mode, str) else graph.cyclo_mode.value
            topo = graph.topology if isinstance(graph.topology, str) else graph.topology.value
            
            graph_dict = {
                "reaction_id": graph.reaction_id,
                "reaction_type": graph.reaction_type,
                "cyclo_mode": cyclo,
                "topology": topo,
                "precursor_type": graph.precursor_type,
                "nodes": [
                    {
                        "node_id": n.node_id,
                        "smiles": n.smiles,
                        "state_type": n.state_type.value if hasattr(n.state_type, 'value') else n.state_type,
                        "role": n.role
                    }
                    for n in graph.nodes
                ],
                "edges": self.get_edges_for_pathway(graph),
                "pathways": [
                    {
                        "pathway_id": p.pathway_id,
                        "description": p.description,
                        "stereochemistry": p.stereochemistry
                    }
                    for p in graph.pathways
                ]
            }
            data.append(graph_dict)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"Exported {len(graphs)} graphs to {output_path}")
    
    def export_summary(
        self,
        graphs: List[MechanismGraph],
        output_path: Path
    ) -> None:
        """
        导出摘要信息
        
        Args:
            graphs: MechanismGraph 列表
            output_path: 输出路径
        """
        summary = {
            "total_graphs": len(graphs),
            "by_cyclo_mode": {},
            "by_topology": {},
            "by_reaction_type": {}
        }
        
        # 统计
        for graph in graphs:
            # 按环加成模式统计
            mode = graph.cyclo_mode if isinstance(graph.cyclo_mode, str) else graph.cyclo_mode.value
            summary["by_cyclo_mode"][mode] = summary["by_cyclo_mode"].get(mode, 0) + 1
            
            # 按拓扑统计
            topo = graph.topology if isinstance(graph.topology, str) else graph.topology.value
            summary["by_topology"][topo] = summary["by_topology"].get(topo, 0) + 1
            
            # 按反应类型统计
            rtype = graph.reaction_type
            summary["by_reaction_type"][rtype] = summary["by_reaction_type"].get(rtype, 0) + 1
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Exported summary to {output_path}")


# ============ 便捷函数 ============

def classify_reaction(
    csv_path: Path,
    reaction_id: str = None,
    reaction_type: str = None
) -> List[MechanismGraph]:
    """
    便捷函数：分类反应
    
    Args:
        csv_path: Clean CSV 路径
        reaction_id: 可选的反应 ID
        reaction_type: 可选的反应类型
        
    Returns:
        MechanismGraph 列表
    """
    classifier = MechanismClassifier()
    return classifier.classify_from_csv(
        csv_path, 
        reaction_ids=[reaction_id] if reaction_id else None,
        reaction_type=reaction_type
    )


def get_mechanism_info(graph: MechanismGraph) -> Dict[str, Any]:
    """
    获取机理信息的便捷函数
    
    Args:
        graph: MechanismGraph
        
    Returns:
        信息字典
    """
    return {
        "reaction_id": graph.reaction_id,
        "reaction_type": graph.reaction_type,
        "cyclo_mode": graph.cyclo_mode.value,
        "topology": graph.topology.value,
        "precursor_type": graph.precursor_type,
        "n_steps": len(graph.edges),
        "is_two_step": len(graph.edges) == 2,
        "pathways": graph.get_pathways()
    }
