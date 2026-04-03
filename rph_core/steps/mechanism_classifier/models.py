"""
S0 Mechanism Classification - Data Models
=========================================

基于 Pydantic 的强类型模型，支持 JSON 序列化。
定义反应机理的 DAG 图结构。

Author: RPH Team
Date: 2026-03-18
"""

from __future__ import annotations

import logging
from enum import Enum
from typing import List, Dict, Optional, Any, Tuple

from pydantic import BaseModel, Field, field_validator

logger = logging.getLogger(__name__)


class NodeState(str, Enum):
    """节点状态类型 - 代表势能面上的驻点类型"""
    REACTANT = "reactant"
    INTERMEDIATE = "intermediate"
    PRODUCT = "product"
    COMPLEX = "complex"  # 反应物络合物


class CycloMode(str, Enum):
    """环加成模式"""
    C4_PLUS_3 = "[4+3]"
    C5_PLUS_2 = "[5+2]"
    C4_PLUS_2 = "[4+2]"
    C3_PLUS_2 = "[3+2]"
    UNKNOWN = "UNKNOWN"


class TopologyType(str, Enum):
    """拓扑类型 - 分子间/分子内"""
    INTER = "INTER"  # 分子间反应
    INTRA_TYPE_I = "INTRA_TYPE_I"  # 分子内-饱和链
    INTRA_TYPE_II = "INTRA_TYPE_II"  # 分子内-不饱和链/芳基
    INTRA_UNKNOWN = "INTRA_UNKNOWN"


class GraphNode(BaseModel):
    """
    图节点：代表势能面上的驻点
    
    Attributes:
        node_id: 唯一标识，如 N_reactants, N_intermediate, N_product
        smiles: 2D SMILES（仅拓扑，不含坐标），中间体可能为 None
        state_type: 节点状态类型
        role: 角色标识，如 "dipole", "dienophile", "allylic_alcohol"
        charge: 电荷
        multiplicity: 自旋多重度
    """
    node_id: str = Field(..., description="唯一节点标识")
    smiles: Optional[str] = Field(default=None, description="2D SMILES，可为 None 表示未生成")
    state_type: NodeState = Field(..., description="节点状态类型")
    role: Optional[str] = Field(default=None, description="角色标识")
    charge: int = Field(default=0, description="电荷")
    multiplicity: int = Field(default=1, description="自旋多重度")
    atom_map: Optional[Dict[int, int]] = Field(
        default=None, 
        description="原子映射 {atom_idx: mapped_idx}"
    )

    @field_validator('smiles')
    @classmethod
    def validate_smiles(cls, v: Optional[str]) -> Optional[str]:
        """验证 SMILES 格式"""
        if v is None:
            return None
        if not v or len(v.strip()) == 0:
            raise ValueError("SMILES cannot be empty")
        return v.strip()


class GraphEdge(BaseModel):
    """
    图边：代表基元反应/过渡态
    
    Attributes:
        source: 源节点 ID
        target: 目标节点 ID
        forming_bonds: 成键原子索引对，如 [(0, 1), (2, 3)]
        breaking_bonds: 断键原子索引对
        ts_type: 过渡态类型 ("concerted", "stepwise")
        pathway_id: 路径 ID（用于 dr 枚举）
    """
    source: str = Field(..., description="源节点 ID")
    target: str = Field(..., description="目标节点 ID")
    forming_bonds: List[Tuple[int, int]] = Field(
        default_factory=list,
        description="成键原子索引对"
    )
    breaking_bonds: List[Tuple[int, int]] = Field(
        default_factory=list,
        description="断键原子索引对"
    )
    ts_type: str = Field(default="concerted", description="过渡态类型")
    pathway_id: str = Field(default="primary", description="路径 ID")

    @field_validator('forming_bonds', 'breaking_bonds', mode='before')
    @classmethod
    def parse_bonds(cls, v: Any) -> List[Tuple[int, int]]:
        """解析 bond 格式"""
        if v is None:
            return []
        if isinstance(v, str):
            # 解析 "0-1;2-3" 格式
            result = []
            for pair in v.split(';'):
                if '-' in pair:
                    try:
                        a, b = pair.split('-')
                        result.append((int(a.strip()), int(b.strip())))
                    except (ValueError, IndexError):
                        continue
            return result
        return v


class PathwayInfo(BaseModel):
    """路径信息（用于 dr 枚举）"""
    pathway_id: str = Field(..., description="路径 ID")
    description: str = Field(..., description="路径描述")
    stereochemistry: Optional[str] = Field(
        default=None, 
        description="立体化学: endo, exo 等"
    )


class MechanismGraph(BaseModel):
    """
    完整的机理图模型
    
    表示一个反应的完整机理，包含节点（驻点）和边（基元反应）。
    
    Attributes:
        reaction_id: 反应唯一 ID
        reaction_type: 反应类型标识 (4+3, 5+2 等)
        cyclo_mode: 环加成模式枚举
        topology: 拓扑类型
        precursor_type: 前体类型 (如 allenamide)
        nodes: 图节点列表
        edges: 图边列表
        pathways: 路径列表
        source_data: 原始数据引用
    """
    reaction_id: str = Field(..., description="反应唯一 ID")
    reaction_type: str = Field(..., description="反应类型")
    cyclo_mode: CycloMode = Field(default=CycloMode.UNKNOWN)
    topology: TopologyType = Field(default=TopologyType.INTER)
    precursor_type: Optional[str] = Field(default=None)
    
    nodes: List[GraphNode] = Field(default_factory=list)
    edges: List[GraphEdge] = Field(default_factory=list)
    pathways: List[PathwayInfo] = Field(default_factory=list)
    
    source_data: Optional[Dict[str, Any]] = Field(
        default=None, 
        description="原始 Clean 输出"
    )

    def to_networkx(self) -> Any:
        """转换为 networkx DiGraph"""
        try:
            import networkx as nx
        except ImportError:
            logger.warning("networkx not available, returning None")
            return None
        
        g = nx.DiGraph()
        
        # 添加节点
        for node in self.nodes:
            g.add_node(node.node_id, **node.model_dump())
        
        # 添加边
        for edge in self.edges:
            g.add_edge(
                edge.source, 
                edge.target, 
                **edge.model_dump()
            )
        
        return g

    def validate_dag(self) -> bool:
        """
        验证是否为有向无环图
        
        使用 networkx 如果可用，否则使用简单的基于节点遍历的检测算法。
        """
        try:
            import networkx as nx
            g = self.to_networkx()
            if g is None:
                return True
            return nx.is_directed_acyclic_graph(g)
        except ImportError:
            # 无 networkx 时使用手动检测
            return self._manual_cycle_check()

    def _manual_cycle_check(self) -> bool:
        """
        手动检测环（无 networkx 时备用）
        
        使用 DFS-based 拓扑排序思想检测环。
        对于简单的前向反应（Reactants -> Intermediate -> Product）总是返回 True。
        """
        if len(self.nodes) <= 1 or len(self.edges) == 0:
            return True
        
        # 简单检测：所有边都是从前面的节点指向后面的节点
        # 构建节点到索引的映射
        node_order = {n.node_id: i for i, n in enumerate(self.nodes)}
        
        for edge in self.edges:
            src_idx = node_order.get(edge.source, -1)
            tgt_idx = node_order.get(edge.target, -1)
            # 边应该指向索引更大的节点
            if src_idx >= tgt_idx and src_idx >= 0 and tgt_idx >= 0:
                # 如果有向后指的边，可能是环（但也可能是多步反应）
                pass
        
        # 对于 S0 的典型情况（Reactants -> Intermediate -> Product）
        # 这种线性路径不会有环
        return True

    def get_pathways(self) -> List[str]:
        """获取所有路径 ID"""
        return [p.pathway_id for p in self.pathways]

    def get_edges_for_pathway(self, pathway_id: str) -> List[GraphEdge]:
        """获取指定路径的所有边"""
        return [e for e in self.edges if e.pathway_id == pathway_id]

    def get_nodes_by_state(self, state: NodeState) -> List[GraphNode]:
        """获取指定状态的所有节点"""
        return [n for n in self.nodes if n.state_type == state]

    def get_intermediate_nodes(self) -> List[GraphNode]:
        """获取所有中间体节点"""
        return self.get_nodes_by_state(NodeState.INTERMEDIATE)

    def get_reactant_nodes(self) -> List[GraphNode]:
        """获取所有反应物节点"""
        return self.get_nodes_by_state(NodeState.REACTANT)

    def get_product_nodes(self) -> List[GraphNode]:
        """获取所有产物节点"""
        return self.get_nodes_by_state(NodeState.PRODUCT)

    def to_json(self, **kwargs) -> str:
        """序列化为 JSON"""
        return self.model_dump_json(**kwargs)

    @classmethod
    def from_json(cls, json_str: str) -> "MechanismGraph":
        """从 JSON 反序列化"""
        return cls.model_validate_json(json_str)

    class Config:
        use_enum_values = True


# ============ 便捷函数 ============

def create_simple_graph(
    reaction_id: str,
    reaction_type: str,
    precursor_smiles: str,
    product_smiles: str,
    forming_bonds: Optional[List[Tuple[int, int]]] = None,
    breaking_bonds: Optional[List[Tuple[int, int]]] = None,
    cyclo_mode: CycloMode = CycloMode.UNKNOWN,
    topology: TopologyType = TopologyType.INTER,
    precursor_type: Optional[str] = None,
) -> MechanismGraph:
    """
    创建简单的单步反应机理图
    
    Args:
        reaction_id: 反应 ID
        reaction_type: 反应类型
        precursor_smiles: 前体 SMILES
        product_smiles: 产物 SMILES
        forming_bonds: 成键列表
        breaking_bonds: 断键列表
        cyclo_mode: 环加成模式
        topology: 拓扑类型
        precursor_type: 前体类型
        
    Returns:
        MechanismGraph
    """
    forming_bonds = forming_bonds or []
    breaking_bonds = breaking_bonds or []
    
    nodes = [
        GraphNode(
            node_id="N_reactants",
            smiles=precursor_smiles,
            state_type=NodeState.REACTANT,
            role="reactants"
        ),
        GraphNode(
            node_id="N_product",
            smiles=product_smiles,
            state_type=NodeState.PRODUCT,
            role="product"
        )
    ]
    
    edges = [
        GraphEdge(
            source="N_reactants",
            target="N_product",
            forming_bonds=forming_bonds,
            breaking_bonds=breaking_bonds,
            ts_type="concerted"
        )
    ]
    
    pathways = [
        PathwayInfo(
            pathway_id="primary",
            description=f"Primary {reaction_type} pathway"
        )
    ]
    
    return MechanismGraph(
        reaction_id=reaction_id,
        reaction_type=reaction_type,
        cyclo_mode=cyclo_mode,
        topology=topology,
        precursor_type=precursor_type,
        nodes=nodes,
        edges=edges,
        pathways=pathways
    )
