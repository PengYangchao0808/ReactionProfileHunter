"""
Mechanism Classifier - S0 Stage
================================

S0 机理分类模块，用于从 Clean 程序输出构建 DAG 图模型。

主要功能：
1. 解析 Clean CSV 输出
2. 构建有向无环图 (DAG) 表示反应机理
3. 输出边属性 (forming_bonds/breaking_bonds) 给 S2

Usage:
    from rph_core.steps.mechanism_classifier import MechanismClassifier
    
    classifier = MechanismClassifier()
    graphs = classifier.classify_from_csv(Path("cleaned/reaxys_cleaned.csv"))
"""

from rph_core.steps.mechanism_classifier.classifier import MechanismClassifier
from rph_core.steps.mechanism_classifier.models import (
    MechanismGraph,
    GraphNode,
    GraphEdge,
    NodeState,
    CycloMode,
    TopologyType,
)
from rph_core.steps.mechanism_classifier.clean_adapter import CleanAdapter, CleanRecord
from rph_core.steps.mechanism_classifier.graph_builder import GraphBuilder

__all__ = [
    "MechanismClassifier",
    "MechanismGraph",
    "GraphNode",
    "GraphEdge",
    "NodeState",
    "CycloMode",
    "TopologyType",
    "CleanAdapter",
    "CleanRecord",
    "GraphBuilder",
]
