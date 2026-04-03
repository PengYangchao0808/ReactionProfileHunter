"""
Mechanism Graph Builder
=======================

基于 Clean 程序输出的记录构建 DAG (有向无环图) 图模型。

主要功能：
1. 将 CleanRecord 转换为 MechanismGraph
2. 解析 core_bond_changes 生成 forming_bonds/breaking_bonds
3. 支持两步反应 (如 [4+3], [5+2]) 的中间体节点
4. 支持 dr (diastereomeric ratio) 路径枚举

Author: RPH Team
Date: 2026-03-18
"""

from __future__ import annotations

import logging
from typing import List, Dict, Optional, Any, Tuple

from rph_core.steps.mechanism_classifier.models import (
    MechanismGraph,
    GraphNode,
    GraphEdge,
    PathwayInfo,
    NodeState,
    CycloMode,
    TopologyType,
)
from rph_core.steps.mechanism_classifier.clean_adapter import CleanRecord

logger = logging.getLogger(__name__)


class GraphBuilder:
    """
    机理图构建器
    
    将 CleanRecord 转换为 MechanismGraph (DAG)
    """
    
    # 环加成模式到节点/边配置的映射
    CYCLO_MODE_CONFIG = {
        "[4+3]": {
            "steps": 2,
            "ts_types": ["stepwise_first_C_C", "concerted"],
        },
        "[5+2]": {
            "steps": 2,
            "ts_types": ["stepwise_first_C_C", "concerted"],
        },
        "[4+2]": {
            "steps": 1,
            "ts_types": ["concerted"],
        },
        "[3+2]": {
            "steps": 1,
            "ts_types": ["concerted"],
        },
    }
    
    def __init__(self, config: Optional[Dict] = None):
        """
        初始化构建器
        
        Args:
            config: 可选配置字典
        """
        self.config = config or {}
    
    def build(self, record: CleanRecord) -> MechanismGraph:
        """
        从 CleanRecord 构建 MechanismGraph
        
        Args:
            record: CleanRecord 对象
            
        Returns:
            MechanismGraph
        """
        # 1. 解析环加成模式
        cyclo_mode = self._normalize_cyclo_mode(record.cyclo_mode)
        
        # 2. 解析拓扑类型
        topology = self._normalize_topology(record.topology)
        
        # 3. 获取配置
        config = self.CYCLO_MODE_CONFIG.get(
            cyclo_mode.value, 
            self.CYCLO_MODE_CONFIG["[4+3]"]
        )
        
        # 4. 创建节点
        nodes = self._create_nodes(record, config)
        
        # 5. 创建边
        edges = self._create_edges(record, config)
        
        # 6. 创建默认路径
        pathways = self._create_default_pathways(cyclo_mode)
        
        # 7. 组装模型
        mechanism_name = f"{cyclo_mode.value}_{record.precursor_type or 'unknown'}"
        
        return MechanismGraph(
            reaction_id=record.reaction_id,
            reaction_type=record.reaction_type,
            cyclo_mode=cyclo_mode,
            topology=topology,
            precursor_type=record.precursor_type,
            nodes=nodes,
            edges=edges,
            pathways=pathways,
            source_data=record.raw
        )
    
    def _normalize_cyclo_mode(self, mode: str) -> CycloMode:
        """
        解析环加成模式字符串
        
        Args:
            mode: 模式字符串 (如 "[4+3]", "4+3", "[4+3]_default")
            
        Returns:
            CycloMode 枚举
        """
        mode = mode.upper().replace(' ', '').strip()
        
        # 直接匹配
        for cyclo in CycloMode:
            if cyclo.value.replace(' ', '') == mode:
                return cyclo
        
        # 兼容 "4+3" 格式
        mode_map = {
            "4+3": CycloMode.C4_PLUS_3,
            "5+2": CycloMode.C5_PLUS_2,
            "4+2": CycloMode.C4_PLUS_2,
            "3+2": CycloMode.C3_PLUS_2,
        }
        
        return mode_map.get(mode, CycloMode.UNKNOWN)
    
    def _normalize_topology(self, topology: str) -> TopologyType:
        """
        解析拓扑类型字符串
        
        Args:
            topology: 拓扑字符串
            
        Returns:
            TopologyType 枚举
        """
        topology = topology.upper().strip()
        
        topo_map = {
            "INTER": TopologyType.INTER,
            "INTRA_TYPE_I": TopologyType.INTRA_TYPE_I,
            "INTRA_TYPE_II": TopologyType.INTRA_TYPE_II,
            "INTRA_UNKNOWN": TopologyType.INTRA_UNKNOWN,
        }
        
        return topo_map.get(topology, TopologyType.INTER)
    
    def _create_nodes(
        self, 
        record: CleanRecord,
        config: Dict
    ) -> List[GraphNode]:
        """
        创建节点列表
        
        Args:
            record: CleanRecord
            config: 环加成配置
            
        Returns:
            GraphNode 列表
        """
        nodes = []
        
        # 1. 反应物节点 (Reactants Complex)
        nodes.append(GraphNode(
            node_id="N_reactants",
            smiles=record.precursor_smiles,
            state_type=NodeState.REACTANT,
            role="reactants",
            charge=0,
            multiplicity=1
        ))
        
        # 2. 中间体节点 (如果是两步反应)
        if config["steps"] >= 2:
            intermediate_smiles = self._infer_intermediate_smiles(
                record.precursor_smiles,
                record.core_bond_changes,
                record.core_atom_map
            )
            
            nodes.append(GraphNode(
                node_id="N_intermediate",
                smiles=intermediate_smiles,
                state_type=NodeState.INTERMEDIATE,
                role="intermediate",
                charge=0,
                multiplicity=1,
                atom_map=record.core_atom_map
            ))
        
        # 3. 产物节点
        nodes.append(GraphNode(
            node_id="N_product",
            smiles=record.product_smiles,
            state_type=NodeState.PRODUCT,
            role="product",
            charge=0,
            multiplicity=1
        ))
        
        return nodes
    
    def _infer_intermediate_smiles(
        self,
        precursor_smiles: str,
        bond_changes: Optional[Dict[str, List[Tuple[int, int]]]],
        atom_map: Optional[Dict[int, int]]
    ) -> Optional[str]:
        """
        推断中间体 SMILES
        
        当前版本返回 None，占位符由 S1/S2 阶段自行生成。
        
        Args:
            precursor_smiles: 前体 SMILES
            bond_changes: 成键/断键变化
            atom_map: 原子映射
            
        Returns:
            中间体 SMILES 或 None
        """
        return None
    
    def _create_edges(
        self,
        record: CleanRecord,
        config: Dict
    ) -> List[GraphEdge]:
        """
        创建边列表
        
        Args:
            record: CleanRecord
            config: 环加成配置
            
        Returns:
            GraphEdge 列表
        """
        edges = []
        bond_changes = record.core_bond_changes
        
        forming = bond_changes.get('forming', [])
        breaking = bond_changes.get('breaking', [])
        
        if config["steps"] == 1:
            # 单步反应：Reactants → Product
            edges.append(GraphEdge(
                source="N_reactants",
                target="N_product",
                forming_bonds=forming,
                breaking_bonds=breaking,
                ts_type=config["ts_types"][0],
                pathway_id="primary"
            ))
        
        else:
            # 两步反应：Reactants → Intermediate → Product
            
            # Step 1: Reactants → Intermediate
            # 第一个成键形成中间体
            step1_forming = [forming[0]] if forming else []
            edges.append(GraphEdge(
                source="N_reactants",
                target="N_intermediate",
                forming_bonds=step1_forming,
                breaking_bonds=breaking,
                ts_type=config["ts_types"][0],
                pathway_id="primary"
            ))
            
            # Step 2: Intermediate → Product
            # 第二个成键形成产物环
            step2_forming = forming[1:] if len(forming) > 1 else forming
            edges.append(GraphEdge(
                source="N_intermediate",
                target="N_product",
                forming_bonds=step2_forming,
                breaking_bonds=[],
                ts_type=config["ts_types"][1] if len(config["ts_types"]) > 1 else "concerted",
                pathway_id="primary"
            ))
        
        return edges
    
    def _create_default_pathways(
        self, 
        cyclo_mode: CycloMode
    ) -> List[PathwayInfo]:
        """创建默认路径"""
        return [
            PathwayInfo(
                pathway_id="primary",
                description=f"Primary {cyclo_mode.value} pathway"
            )
        ]
    
    def add_dr_pathway(
        self,
        base_graph: MechanismGraph,
        dr_type: str,
        description: str,
        stereochemistry: str = None
    ) -> MechanismGraph:
        """
        添加 dr (diastereomeric ratio) 路径
        
        用于平行枚举 endo/exo 等立体选择性路径
        
        Args:
            base_graph: 基础机理图
            dr_type: dr 类型 (如 "endo", "exo")
            description: 路径描述
            stereochemistry: 立体化学描述
            
        Returns:
            新的 MechanismGraph (包含 dr 路径)
        """
        # 复制基础图 (深拷贝)
        import copy
        new_graph = copy.deepcopy(base_graph)
        
        # 为每条边添加新的 pathway_id
        pathway_id = f"dr_{dr_type}"
        
        new_edges = []
        for edge in new_graph.edges:
            new_edge = GraphEdge(
                source=edge.source,
                target=edge.target,
                forming_bonds=edge.forming_bonds,
                breaking_bonds=edge.breaking_bonds,
                ts_type=edge.ts_type,
                pathway_id=pathway_id
            )
            new_edges.append(new_edge)
        
        new_graph.edges = new_edges
        
        # 添加路径信息
        new_graph.pathways.append(PathwayInfo(
            pathway_id=pathway_id,
            description=description,
            stereochemistry=stereochemistry
        ))
        
        return new_graph
    
    def build_batch(
        self,
        records: List[CleanRecord]
    ) -> List[MechanismGraph]:
        """
        批量构建
        
        Args:
            records: CleanRecord 列表
            
        Returns:
            MechanismGraph 列表
        """
        graphs = []
        
        for record in records:
            try:
                graph = self.build(record)
                
                # 验证 DAG
                if graph.validate_dag():
                    graphs.append(graph)
                else:
                    logger.warning(
                        f"Cycle detected in graph for {record.reaction_id}, "
                        "skipping..."
                    )
                    
            except Exception as e:
                logger.warning(
                    f"Failed to build graph for {record.reaction_id}: {e}"
                )
                continue
        
        logger.info(f"Built {len(graphs)} graphs from {len(records)} records")
        return graphs
    
    def get_edge_summary(self, graph: MechanismGraph) -> Dict[str, Any]:
        """
        获取边的摘要信息
        
        用于调试和日志
        
        Args:
            graph: MechanismGraph
            
        Returns:
            摘要字典
        """
        summary = {
            "reaction_id": graph.reaction_id,
            "cyclo_mode": graph.cyclo_mode.value,
            "topology": graph.topology.value,
            "n_nodes": len(graph.nodes),
            "n_edges": len(graph.edges),
            "edges": []
        }
        
        for edge in graph.edges:
            summary["edges"].append({
                "source": edge.source,
                "target": edge.target,
                "forming_bonds": edge.forming_bonds,
                "breaking_bonds": edge.breaking_bonds,
                "ts_type": edge.ts_type,
                "pathway": edge.pathway_id
            })
        
        return summary
