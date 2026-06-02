"""
V3 UI 进度数据模型。

本模块仅定义 V3 批处理界面所需的纯数据结构，
用于承载规划阶段、小分子缓存阶段、反应执行阶段以及最终汇总阶段的状态。

设计原则：
- 仅包含数据模型，不包含任何业务流程逻辑。
- 显示层与状态层解耦，UI 渲染代码不应放在这里。
- 所有可变默认值均通过 ``field(default_factory=...)`` 提供。
- 路径相关字段保留为字符串，便于 JSON/状态快照序列化。

这些模型服务于 ReactionProfileHunter V3 UI 重构，
作为运行期状态快照在 orchestrator / UI 层之间传递。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto

__all__ = [
    "V3StageState",
    "V3ItemState",
    "ReactionPlanRecord",
    "PlanningProgress",
    "SmallMoleculeRecord",
    "SmallMoleculeCacheProgress",
    "BranchProgress",
    "ReactionExecutionProgress",
    "FinalSummary",
    "V3RunProgress",
]


# ============================================================================
# V3 状态枚举
# ============================================================================


class V3StageState(Enum):
    """V3 顶层阶段状态。

    该枚举用于描述较粗粒度的阶段级状态，
    例如 planning / small_molecule_cache / reaction_execution。

    约定：
    - PENDING: 尚未进入该阶段
    - RUNNING: 当前正在执行
    - COMPLETED: 阶段已成功结束
    - FAILED: 阶段执行失败
    """

    PENDING = auto()
    RUNNING = auto()
    COMPLETED = auto()
    FAILED = auto()

    @property
    def icon(self) -> str:
        """返回用于终端/UI 展示的图标。"""
        return {
            V3StageState.PENDING: "○",
            V3StageState.RUNNING: "→",
            V3StageState.COMPLETED: "✓",
            V3StageState.FAILED: "✗",
        }[self]

    @property
    def color(self) -> str:
        """返回与状态对应的颜色名称。"""
        return {
            V3StageState.PENDING: "dim",
            V3StageState.RUNNING: "cyan",
            V3StageState.COMPLETED: "green",
            V3StageState.FAILED: "red",
        }[self]


class V3ItemState(Enum):
    """V3 细粒度条目状态。

    该枚举用于描述阶段内部的单个对象状态，
    例如单个 reaction、branch、small molecule key 或 condition 记录。

    除常规状态外，还包括：
    - CACHED: 结果直接来自缓存
    - WAITING_LOCK: 等待锁，通常表示共享资源占用中
    - THERMO_DERIVED: 热力学校正由已有结果派生
    - MISSING_THERMO: 缺失热力学信息，但流程继续降级运行
    """

    PENDING = auto()
    RUNNING = auto()
    COMPLETED = auto()
    FAILED = auto()
    CACHED = auto()
    WAITING_LOCK = auto()
    THERMO_DERIVED = auto()
    MISSING_THERMO = auto()

    @property
    def icon(self) -> str:
        """返回条目状态图标。"""
        return {
            V3ItemState.PENDING: "○",
            V3ItemState.RUNNING: "→",
            V3ItemState.COMPLETED: "✓",
            V3ItemState.FAILED: "✗",
            V3ItemState.CACHED: "↻",
            V3ItemState.WAITING_LOCK: "⌛",
            V3ItemState.THERMO_DERIVED: "Δ",
            V3ItemState.MISSING_THERMO: "?",
        }[self]

    @property
    def color(self) -> str:
        """返回条目状态颜色。"""
        return {
            V3ItemState.PENDING: "dim",
            V3ItemState.RUNNING: "cyan",
            V3ItemState.COMPLETED: "green",
            V3ItemState.FAILED: "red",
            V3ItemState.CACHED: "blue",
            V3ItemState.WAITING_LOCK: "yellow",
            V3ItemState.THERMO_DERIVED: "magenta",
            V3ItemState.MISSING_THERMO: "yellow",
        }[self]

    @property
    def label(self) -> str:
        """返回稳定的机器可读标签。"""
        return {
            V3ItemState.PENDING: "pending",
            V3ItemState.RUNNING: "running",
            V3ItemState.COMPLETED: "completed",
            V3ItemState.FAILED: "failed",
            V3ItemState.CACHED: "cached",
            V3ItemState.WAITING_LOCK: "waiting_lock",
            V3ItemState.THERMO_DERIVED: "thermo_derived",
            V3ItemState.MISSING_THERMO: "missing_thermo",
        }[self]


# ============================================================================
# Planning 阶段数据模型
# ============================================================================


@dataclass
class ReactionPlanRecord:
    """单个反应在 planning 阶段的状态记录。

    该记录用于描述：
    - 当前 reaction_id 是否已进入规划流程
    - S0 机制分析是否 ready / reused / failed
    - 已形成多少 branch / condition
    - 规划阶段的即时细节文本

    ``branch_ids`` 保留稳定顺序，
    便于上层 UI 以预期顺序显示主支路与次支路。
    """

    reaction_id: str
    state: V3ItemState = V3ItemState.PENDING
    s0_status: str = ""
    branch_count: int = 0
    condition_count: int = 0
    branch_ids: list[str] = field(default_factory=list)
    detail: str = ""


@dataclass
class PlanningProgress:
    """planning 阶段总体进度。

    字段语义：
    - ``total_rows``: 输入表原始行数
    - ``total_reactions``: 实际形成的反应记录数
    - ``total_conditions``: 规划期已解析的条件总数
    - ``reaction_plans``: 逐反应的规划记录
    - ``branch_sequence_formed``: 是否已建立全局 branch 顺序
    - ``detail``: 阶段级的一行摘要

    该对象为 V3 顶层进度树中的第一个阶段节点。
    """

    state: V3StageState = V3StageState.PENDING
    total_rows: int = 0
    total_reactions: int = 0
    total_conditions: int = 0
    reaction_plans: dict[str, ReactionPlanRecord] = field(default_factory=dict)
    branch_sequence_formed: bool = False
    detail: str = ""


# ============================================================================
# Small Molecule Cache 阶段数据模型
# ============================================================================


@dataclass
class SmallMoleculeRecord:
    """单个小分子缓存键的状态记录。

    这里的 ``key`` 通常为缓存命中与共享引用的稳定标识。
    ``detail`` 可用于承载诸如：
    - ``cache hit``
    - ``computing``
    - ``optimization_only``
    - ``failed on SP``

    该对象不负责缓存逻辑，只负责承载状态快照。
    """

    key: str
    state: V3ItemState = V3ItemState.PENDING
    detail: str = ""


@dataclass
class SmallMoleculeCacheProgress:
    """小分子全局缓存阶段进度。

    该阶段独立于单个 reaction 的主执行流程，
    主要反映共享小分子 anchor/cache 的预热与复用情况。

    统计字段说明：
    - ``summary_total``: 总 key 数
    - ``summary_hit``: 缓存命中数
    - ``summary_computed``: 新完成计算数
    - ``summary_failed``: 失败 key 数
    - ``summary_pending``: 尚未结束的 key 数

    ``cache_root`` 为字符串路径，便于状态序列化。
    """

    state: V3StageState = V3StageState.PENDING
    cache_root: str = ""
    molecules: dict[str, SmallMoleculeRecord] = field(default_factory=dict)
    summary_total: int = 0
    summary_hit: int = 0
    summary_computed: int = 0
    summary_failed: int = 0
    summary_pending: int = 0
    detail: str = ""


# ============================================================================
# Reaction Execution 阶段数据模型
# ============================================================================


@dataclass
class BranchProgress:
    """单个 branch 的执行进度。

    该记录聚焦于 branch 级别的工作流状态，覆盖：
    - 是否已根据计划跳过 S0
    - precursor 是否已关联
    - 小分子缓存引用是否已物化
    - 当前 pipeline 所在步骤及其条目状态
    - S3/S4 之后可填充的 ``dg_activation``

    ``pipeline_step`` 预期值通常为 ``s1`` / ``s2`` / ``s3`` / ``s4``。
    """

    branch_id: str
    state: V3ItemState = V3ItemState.PENDING
    s0_skipped_by_plan: bool = False
    precursor_linked: bool = False
    sm_refs_materialized: bool = False
    pipeline_step: str = ""
    pipeline_step_state: V3ItemState = V3ItemState.PENDING
    pipeline_detail: str = ""
    dg_activation: "float | None" = None
    detail: str = ""


@dataclass
class ReactionExecutionProgress:
    """单个 reaction 的执行态快照。

    该对象用于承载当前反应在执行阶段的完整视图，
    既包含 branch 维度进度，也包含后处理与 condition 维度状态。

    关键字段：
    - ``current_branch_index`` / ``total_branches``: branch 遍历位置
    - ``shared_precursor_s1_state``: 共享 precursor S1 状态
    - ``sm_cache_refs_ready``: 小分子缓存引用是否就绪
    - ``branches``: branch_id -> BranchProgress
    - ``dr_aggregation_state``: 分支聚合结果状态
    - ``condition_records``: ``COND_id/BR_id`` 粒度状态
    - ``condition_dr_state``: condition 级 DR 汇总状态
    """

    reaction_id: str
    state: V3StageState = V3StageState.PENDING
    current_branch_index: int = 0
    total_branches: int = 0
    shared_precursor_s1_state: V3ItemState = V3ItemState.PENDING
    shared_precursor_s1_detail: str = ""
    sm_cache_refs_ready: bool = False
    branches: dict[str, BranchProgress] = field(default_factory=dict)
    dr_aggregation_state: V3ItemState = V3ItemState.PENDING
    dr_aggregation_detail: str = ""
    condition_post_processing: bool = False
    condition_total: int = 0
    condition_branch_total: int = 0
    condition_records: dict[str, V3ItemState] = field(default_factory=dict)
    condition_dr_state: dict[str, V3ItemState] = field(default_factory=dict)
    detail: str = ""


# ============================================================================
# 最终汇总数据模型
# ============================================================================


@dataclass
class FinalSummary:
    """V3 运行完成后的最终汇总。

    该对象聚合规划、缓存、branch、condition 等多个维度的结果，
    适合作为最终 UI 面板、日志摘要或 JSON 运行报告的上游数据源。

    记录粒度：
    - reaction 维度: success / failed / partial
    - branch 维度: success / failed
    - condition 维度: complete / failed
    - small molecule cache 维度: hit / computed / failed

    所有 ``failed_*`` 列表均保留原始标识串，
    便于上层直接显示或进一步索引明细文件。
    """

    total_reactions: int = 0
    reaction_success: int = 0
    reaction_failed: int = 0
    reaction_records: dict[str, str] = field(default_factory=dict)
    total_branches: int = 0
    branch_success: int = 0
    branch_failed: int = 0
    branch_records: dict[str, str] = field(default_factory=dict)
    total_conditions: int = 0
    condition_complete: int = 0
    condition_failed: int = 0
    sm_cache_total: int = 0
    sm_cache_hit: int = 0
    sm_cache_computed: int = 0
    sm_cache_failed: int = 0
    failed_planning: list[str] = field(default_factory=list)
    failed_cache: list[str] = field(default_factory=list)
    failed_branches: list[str] = field(default_factory=list)
    failed_conditions: list[str] = field(default_factory=list)
    output_root: str = ""
    elapsed_sec: float = 0.0


# ============================================================================
# 顶层运行进度容器
# ============================================================================


@dataclass
class V3RunProgress:
    """V3 运行进度的顶层容器。

    这是 V3 UI 重构中最外层的状态对象，
    负责将多个阶段的独立快照拼装为一个统一结构。

    结构层次：
    - ``planning``: 规划阶段
    - ``small_molecule_cache``: 小分子缓存阶段
    - ``reaction_execution``: 当前 reaction 的执行状态
    - ``summary``: 已完成或最终的统计摘要

    ``current_reaction_id`` 与 ``current_reaction_index``
    用于在 batch 运行时显示当前焦点位置。
    """

    planning: PlanningProgress = field(default_factory=PlanningProgress)
    small_molecule_cache: SmallMoleculeCacheProgress = field(
        default_factory=SmallMoleculeCacheProgress
    )
    current_reaction_index: int = 0
    total_reactions: int = 0
    current_reaction_id: str = ""
    reaction_execution: "ReactionExecutionProgress | None" = None
    summary: FinalSummary = field(default_factory=FinalSummary)
