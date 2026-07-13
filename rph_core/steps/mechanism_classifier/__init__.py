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

from __future__ import annotations

from importlib import import_module

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
    "build_dr_branch_plan",
]

_EXPORTS = {
    "MechanismClassifier": ("classifier", "MechanismClassifier"),
    "MechanismGraph": ("models", "MechanismGraph"),
    "GraphNode": ("models", "GraphNode"),
    "GraphEdge": ("models", "GraphEdge"),
    "NodeState": ("models", "NodeState"),
    "CycloMode": ("models", "CycloMode"),
    "TopologyType": ("models", "TopologyType"),
    "CleanAdapter": ("clean_adapter", "CleanAdapter"),
    "CleanRecord": ("clean_adapter", "CleanRecord"),
    "GraphBuilder": ("graph_builder", "GraphBuilder"),
    "build_dr_branch_plan": ("dr_completion", "build_dr_branch_plan"),
}


def __getattr__(name: str):
    """Load legacy graph machinery only when a caller actually needs it."""

    try:
        module_name, attribute_name = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(name) from exc
    module = import_module(f"{__name__}.{module_name}")
    value = getattr(module, attribute_name)
    globals()[name] = value
    return value
