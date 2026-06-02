# pyright: reportAny=false, reportExplicitAny=false, reportDeprecated=false, reportUnusedImport=false, reportUnknownArgumentType=false, reportUnknownVariableType=false, reportUnknownMemberType=false, reportUnnecessaryComparison=false, reportImplicitOverride=false, reportUnusedCallResult=false, reportImplicitStringConcatenation=false
"""
V3 阶段化运行进度显示渲染器。

本模块只负责把 ``V3RunProgress`` 及其子记录渲染为顺序日志式的 Rich 输出，
不承担调度、状态推进、缓存决策或任何业务逻辑。

设计目标：
- 输出严格面向日志模式（sequential log output）
- 复用共享 Console 单例，避免和其他 UI 组件冲突
- 对尚在演进中的 V3 数据模型保持宽容读取，尽量从字段名推断展示信息
- 在重复调用 update 钩子时做最小去重，避免同一行被反复刷屏
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping, Sequence
from math import isfinite
from typing import TYPE_CHECKING, Any, Optional

from rich.console import Console
from rich.markup import escape

from rph_core.utils.shared_console import get_console

if TYPE_CHECKING:
    from rph_core.utils.v3_progress import (
        BranchProgress,
        FinalSummary,
        PlanningProgress,
        ReactionExecutionProgress,
        ReactionPlanRecord,
        SmallMoleculeCacheProgress,
        SmallMoleculeRecord,
        V3ItemState,
        V3RunProgress,
        V3StageState,
    )


console: Console = get_console()

_SEPARATOR = "═" * 60
_REACTION_STATUS_WIDTH = 24
_CACHE_KEY_WIDTH = 12
_CACHE_STATUS_WIDTH = 15
_SETUP_LABEL_WIDTH = 28
_BRANCH_LABEL_WIDTH = 28
_CONDITION_LABEL_WIDTH = 26
_SUMMARY_LABEL_WIDTH = 16

_STATE_STYLE_FALLBACKS: dict[str, tuple[str, str]] = {
    "PENDING": ("○", "dim"),
    "RUNNING": ("→", "cyan"),
    "COMPLETED": ("✓", "green"),
    "FAILED": ("✗", "red"),
    "CACHED": ("✓", "green"),
    "WAITING_LOCK": ("○", "yellow"),
    "THERMO_DERIVED": ("✓", "green"),
    "MISSING_THERMO": ("✗", "red"),
}

_CACHE_STATUS_LABELS: dict[str, str] = {
    "PENDING": "pending",
    "RUNNING": "computing",
    "COMPLETED": "computed",
    "FAILED": "failed",
    "CACHED": "cache hit",
    "WAITING_LOCK": "waiting lock",
    "THERMO_DERIVED": "thermo derived",
    "MISSING_THERMO": "missing thermo",
}

_GENERIC_STATUS_LABELS: dict[str, str] = {
    "PENDING": "pending",
    "RUNNING": "running",
    "COMPLETED": "ready",
    "FAILED": "failed",
    "CACHED": "reused",
    "WAITING_LOCK": "waiting lock",
    "THERMO_DERIVED": "thermo derived",
    "MISSING_THERMO": "missing thermo",
}

_PIPELINE_STEP_ORDER: tuple[tuple[str, str], ...] = (
    ("s1", "S1 product anchor"),
    ("s2", "S2 TS guess"),
    ("s3", "S3 intermediate/TS/SP"),
    ("s4", "S4 features"),
)

_PIPELINE_STEP_ALIASES: dict[str, str] = {
    "s1": "s1",
    "s1 product anchor": "s1",
    "product anchor": "s1",
    "product_anchor": "s1",
    "product-anchor": "s1",
    "anchor": "s1",
    "precursor s1": "s1",
    "s2": "s2",
    "s2 ts guess": "s2",
    "ts guess": "s2",
    "ts_guess": "s2",
    "retro scan": "s2",
    "retro_scan": "s2",
    "scan": "s2",
    "s3": "s3",
    "s3 intermediate/ts/sp": "s3",
    "intermediate/ts/sp": "s3",
    "sp matrix": "s3",
    "sp_matrix": "s3",
    "intermediate": "s3",
    "ts/sp": "s3",
    "s4": "s4",
    "s4 features": "s4",
    "features": "s4",
    "feature extraction": "s4",
    "feature_extraction": "s4",
}


class V3StageDisplay(ABC):
    """V3 阶段化进度显示接口。"""

    @abstractmethod
    def stage0_planning_header(self, progress: V3RunProgress) -> None:
        """渲染 Stage 0 标题与数据集统计。"""

    @abstractmethod
    def stage0_planning_update(self, progress: V3RunProgress) -> None:
        """渲染 Stage 0 的逐反应规划更新。"""

    @abstractmethod
    def stage0_planning_complete(self, progress: V3RunProgress) -> None:
        """渲染 Stage 0 完成提示。"""

    @abstractmethod
    def stage1_sm_cache_header(self, progress: V3RunProgress) -> None:
        """渲染 Stage 1 标题与缓存根目录。"""

    @abstractmethod
    def stage1_sm_cache_update(self, progress: V3RunProgress) -> None:
        """渲染 Stage 1 的逐小分子缓存更新。"""

    @abstractmethod
    def stage1_sm_cache_complete(self, progress: V3RunProgress) -> None:
        """渲染 Stage 1 汇总信息。"""

    @abstractmethod
    def stage2_reaction_header(self, progress: V3RunProgress) -> None:
        """渲染 Stage 2 当前反应头部。"""

    @abstractmethod
    def stage2_reaction_setup(self, progress: V3RunProgress) -> None:
        """渲染 Stage 2 共享反应准备区块。"""

    @abstractmethod
    def stage2_branch_start(self, progress: V3RunProgress, branch_id: str) -> None:
        """渲染某个分支开始时的快照。"""

    @abstractmethod
    def stage2_branch_update(self, progress: V3RunProgress, branch_id: str) -> None:
        """渲染某个分支运行中的快照。"""

    @abstractmethod
    def stage2_branch_complete(self, progress: V3RunProgress, branch_id: str) -> None:
        """渲染某个分支完成时的快照。"""

    @abstractmethod
    def stage2_dr_aggregation(self, progress: V3RunProgress) -> None:
        """渲染分支结果汇总与 DR 聚合状态。"""

    @abstractmethod
    def stage2_condition_header(self, progress: V3RunProgress) -> None:
        """渲染条件后处理区块头部。"""

    @abstractmethod
    def stage2_condition_update(
        self, progress: V3RunProgress, cond_branch_key: str
    ) -> None:
        """渲染某个 condition-branch 组合任务的更新。"""

    @abstractmethod
    def stage2_reaction_complete(self, progress: V3RunProgress) -> None:
        """渲染反应收尾。"""

    @abstractmethod
    def stage3_summary(self, progress: V3RunProgress) -> None:
        """渲染最终运行摘要。"""


class LoggerV3StageDisplay(V3StageDisplay):
    """V3 默认日志式显示器。"""

    def __init__(self) -> None:
        """初始化顺序日志渲染器及其去重缓存。"""
        self.console: Console = console
        self._stage0_header_snapshot: str = ""
        self._stage1_header_snapshot: str = ""
        self._stage2_header_snapshot: str = ""
        self._stage3_summary_snapshot: str = ""

        self._planning_snapshots: dict[str, str] = {}
        self._cache_snapshots: dict[str, str] = {}
        self._reaction_setup_snapshots: dict[str, str] = {}
        self._branch_snapshots: dict[tuple[str, str], str] = {}
        self._dr_snapshots: dict[str, str] = {}
        self._condition_header_snapshots: dict[str, str] = {}
        self._condition_line_snapshots: dict[tuple[str, str], str] = {}
        self._condition_dr_snapshots: dict[tuple[str, str], str] = {}
        self._condition_dr_section_printed: set[str] = set()
        self._reaction_complete_seen: set[str] = set()
        self._completion_snapshots: dict[str, str] = {}

    # ------------------------------------------------------------------
    # 基础输出辅助
    # ------------------------------------------------------------------

    def _print(self, message: str = "") -> None:
        """输出单行文本。"""
        self.console.print(message)
        self._flush()

    def _flush(self) -> None:
        """强制刷新输出缓冲区，避免终端/管道阻塞。"""
        try:
            if hasattr(self.console, "file") and self.console.file:
                self.console.file.flush()
        except Exception:
            pass

    def _safe_text(self, value: Any) -> str:
        """将动态文本转义为 Rich 安全文本。"""
        return escape("" if value is None else str(value))

    def _pad(self, value: Any, width: int) -> str:
        """先按原始文本补齐宽度，再做 Rich 转义。"""
        raw = "" if value is None else str(value)
        return escape(f"{raw:<{width}}")

    def _state_name(self, state: Any) -> str:
        """尽量把枚举/字符串/对象状态归一为大写名字。"""
        if state is None:
            return ""
        if isinstance(state, str):
            return state.strip().replace("-", "_").replace(" ", "_").upper()
        name = getattr(state, "name", None)
        if isinstance(name, str) and name:
            return name.strip().replace("-", "_").replace(" ", "_").upper()
        label = getattr(state, "label", None)
        if isinstance(label, str) and label:
            return label.strip().replace("-", "_").replace(" ", "_").upper()
        return str(state).split(".")[-1].strip().replace("-", "_").replace(" ", "_").upper()

    def _state_markup(self, state: Any, *, fallback: str = "PENDING") -> str:
        """生成带 Rich 颜色的状态图标。"""
        icon = getattr(state, "icon", None)
        color = getattr(state, "color", None)
        if icon and color:
            return f"[{color}]{self._safe_text(icon)}[/]"
        state_name = self._state_name(state) or fallback
        icon_text, color_text = _STATE_STYLE_FALLBACKS.get(
            state_name, _STATE_STYLE_FALLBACKS[fallback]
        )
        return f"[{color_text}]{icon_text}[/]"

    def _status_label(self, state: Any, mapping: dict[str, str]) -> str:
        """根据状态枚举名称生成展示文案。"""
        state_name = self._state_name(state)
        if state_name in mapping:
            return mapping[state_name]
        label = getattr(state, "label", None)
        if isinstance(label, str) and label:
            return label.replace("_", " ")
        if state_name:
            return state_name.lower().replace("_", " ")
        return "pending"

    def _fingerprint(self, *parts: Any) -> str:
        """将若干值压缩为稳定快照键。"""
        return "|".join(repr(part) for part in parts)

    def _emit_if_changed(
        self, store: dict[Any, str], key: Any, block: str
    ) -> bool:
        """仅当 block 变化时输出，以减少重复日志。"""
        if not block:
            return False
        if store.get(key) == block:
            return False
        store[key] = block
        for line in block.splitlines():
            self._print(line)
        return True

    def _emit_header_if_changed(self, attr_name: str, block: str) -> bool:
        """用于少量全局头部的去重输出。"""
        if not block:
            return False
        if getattr(self, attr_name) == block:
            return False
        setattr(self, attr_name, block)
        for line in block.splitlines():
            self._print(line)
        return True

    def _get_attr(self, obj: Any, name: str, default: Any = None) -> Any:
        """同时兼容 dataclass 对象与字典式记录。"""
        if obj is None:
            return default
        if isinstance(obj, Mapping):
            return obj.get(name, default)
        return getattr(obj, name, default)

    def _get_first(self, obj: Any, names: Sequence[str], default: Any = None) -> Any:
        """按候选字段名顺序提取首个非空值。"""
        for name in names:
            value = self._get_attr(obj, name, None)
            if value is not None:
                return value
        return default

    def _iter_values(self, container: Any) -> list[Any]:
        """把 mapping/list/tuple 统一为值列表。"""
        if container is None:
            return []
        if isinstance(container, Mapping):
            return list(container.values())
        if isinstance(container, Sequence) and not isinstance(
            container, (str, bytes, bytearray)
        ):
            return list(container)
        return []

    def _format_metric_line(self, label: str, value: str) -> str:
        """格式化头部中的绿色勾选统计行。"""
        return (
            f"{self._state_markup('COMPLETED')} "
            f"{self._pad(label, 32)}{self._safe_text(value)}"
        )

    def _format_status_line(
        self,
        label: str,
        state: Any,
        detail: str,
        *,
        indent: str = "  ",
        label_width: int = _SETUP_LABEL_WIDTH,
    ) -> str:
        """格式化通用的“图标 + 标签 + 详情”行。"""
        return (
            f"{indent}{self._state_markup(state)} "
            f"{self._pad(label, label_width)}{self._safe_text(detail)}"
        ).rstrip()

    def _to_float(self, value: Any) -> Optional[float]:
        """把任意数值候选转为 float。"""
        if value is None:
            return None
        try:
            number = float(value)
        except (TypeError, ValueError):
            return None
        if not isfinite(number):
            return None
        return number

    def _format_dg(self, value: Any) -> str:
        """格式化 ΔG‡ 数值。"""
        number = self._to_float(value)
        if number is None:
            return ""
        return f"ΔG‡={number:.1f}"

    def _planning_progress(self, progress: Any) -> Any:
        """提取 Stage 0 规划进度对象。"""
        return self._get_attr(progress, "planning", None)

    def _sm_cache_progress(self, progress: Any) -> Any:
        """提取 Stage 1 小分子缓存进度对象。"""
        return self._get_attr(progress, "small_molecule_cache", None)

    def _reaction_progress(self, progress: Any) -> Any:
        """提取 Stage 2 反应执行进度对象。"""
        return self._get_attr(progress, "reaction_execution", None)

    def _summary_progress(self, progress: Any) -> Any:
        """提取最终汇总对象。"""
        return self._get_attr(progress, "summary", None)

    def _planning_records(self, planning: Any) -> list[Any]:
        """获取所有反应规划记录。"""
        return self._iter_values(self._get_attr(planning, "reaction_plans", {}))

    def _sm_records(self, cache_progress: Any) -> list[Any]:
        """获取所有小分子缓存记录。"""
        return self._iter_values(self._get_attr(cache_progress, "molecules", {}))

    def _branch_records(self, reaction_progress: Any) -> list[Any]:
        """获取所有分支记录。"""
        return self._iter_values(self._get_attr(reaction_progress, "branches", {}))

    def _get_planning_record(self, progress: Any, reaction_id: str) -> Any:
        """按 reaction_id 从规划阶段记录里查找单个反应。"""
        planning = self._planning_progress(progress)
        plans = self._get_attr(planning, "reaction_plans", None)
        if isinstance(plans, Mapping):
            if reaction_id in plans:
                return plans[reaction_id]
            for plan in plans.values():
                if self._get_first(plan, ["reaction_id", "id", "key"], "") == reaction_id:
                    return plan
        for plan in self._planning_records(planning):
            if self._get_first(plan, ["reaction_id", "id", "key"], "") == reaction_id:
                return plan
        return None

    def _get_branch_record(self, reaction_progress: Any, branch_id: str) -> Any:
        """按 branch_id 从反应执行记录里查找单个分支。"""
        branches = self._get_attr(reaction_progress, "branches", None)
        if isinstance(branches, Mapping):
            if branch_id in branches:
                return branches[branch_id]
            for branch in branches.values():
                if self._get_first(branch, ["branch_id", "id", "key"], "") == branch_id:
                    return branch
        for branch in self._branch_records(reaction_progress):
            if self._get_first(branch, ["branch_id", "id", "key"], "") == branch_id:
                return branch
        return None

    def _reaction_id(self, progress: Any, reaction_progress: Any = None) -> str:
        """优先从顶层 current_reaction_id 提取当前反应 ID。"""
        reaction_id = self._get_attr(progress, "current_reaction_id", None)
        if reaction_id:
            return str(reaction_id)
        if reaction_progress is None:
            reaction_progress = self._reaction_progress(progress)
        return str(self._get_first(reaction_progress, ["reaction_id", "id", "key"], ""))

    def _reaction_index_total(self, progress: Any) -> tuple[int, int]:
        """提取当前反应索引与总反应数。"""
        index = self._get_attr(progress, "current_reaction_index", 0) or 0
        total = self._get_attr(progress, "total_reactions", 0) or 0
        if not total:
            planning = self._planning_progress(progress)
            total = self._get_attr(planning, "total_reactions", 0) or 0
        try:
            index_int = int(index)
        except (TypeError, ValueError):
            index_int = 0
        try:
            total_int = int(total)
        except (TypeError, ValueError):
            total_int = 0
        return index_int, total_int

    def _planning_status_text(self, record: Any) -> str:
        """生成 Stage 0 单个反应的状态文案。"""
        state = self._get_attr(record, "state", None)
        state_name = self._state_name(state)
        s0_status = self._get_attr(record, "s0_status", None)
        detail = self._get_attr(record, "detail", None)
        if state_name in {"RUNNING", "FAILED"} and detail:
            return str(detail)
        if s0_status:
            return str(s0_status)
        if detail:
            return str(detail)
        if state_name == "CACHED":
            return "S0_Mechanism reused"
        if state_name == "COMPLETED":
            return "S0_Mechanism ready"
        if state_name == "RUNNING":
            return "building mechanism graph..."
        if state_name == "FAILED":
            return "mechanism planning failed"
        return "pending"

    def _cache_status_text(self, record: Any) -> str:
        """生成 Stage 1 单个小分子的状态文案。"""
        state = self._get_attr(record, "state", None)
        detail = self._get_attr(record, "status", None)
        if detail:
            return str(detail)
        return self._status_label(state, _CACHE_STATUS_LABELS)

    def _default_shared_precursor_detail(self, state: Any) -> str:
        """为共享 precursor S1 状态生成兜底说明。"""
        state_name = self._state_name(state)
        if state_name == "RUNNING":
            return "conformer search / DFT opt"
        if state_name in {"COMPLETED", "CACHED"}:
            return "ready"
        if state_name == "FAILED":
            return "shared precursor setup failed"
        return "waiting"

    def _default_sm_cache_detail(self, state: Any) -> str:
        """为共享小分子缓存引用生成兜底说明。"""
        state_name = self._state_name(state)
        if state_name in {"COMPLETED", "CACHED", "THERMO_DERIVED"}:
            return "ready"
        if state_name == "RUNNING":
            return "resolving global cache refs"
        if state_name == "FAILED":
            return "small molecule cache refs failed"
        return "waiting"

    def _planning_to_setup_detail(self, plan_record: Any) -> str:
        """把 Stage 0 规划记录转换为 Stage 2 的 mechanism graph 文案。"""
        if plan_record is None:
            return "waiting for Stage 0"
        s0_status = str(self._get_attr(plan_record, "s0_status", "") or "").lower()
        detail = self._get_attr(plan_record, "detail", None)
        if "reused" in s0_status or "cache" in s0_status:
            return "reused from Stage 0"
        if "ready" in s0_status or self._state_name(self._get_attr(plan_record, "state", None)) in {
            "COMPLETED",
            "CACHED",
        }:
            return "reused from Stage 0"
        if detail:
            return str(detail)
        return "waiting for Stage 0"

    def _branch_ids(self, progress: Any, reaction_progress: Any) -> list[str]:
        """提取当前反应的分支顺序列表。"""
        branch_ids = self._get_first(reaction_progress, ["branch_ids", "branch_sequence"], None)
        if isinstance(branch_ids, Sequence) and not isinstance(
            branch_ids, (str, bytes, bytearray)
        ):
            return [str(item) for item in branch_ids]
        plan = self._get_planning_record(progress, self._reaction_id(progress, reaction_progress))
        plan_branch_ids = self._get_attr(plan, "branch_ids", None)
        if isinstance(plan_branch_ids, Sequence) and not isinstance(
            plan_branch_ids, (str, bytes, bytearray)
        ):
            return [str(item) for item in plan_branch_ids]
        branches = self._get_attr(reaction_progress, "branches", None)
        if isinstance(branches, Mapping):
            return [str(key) for key in branches.keys()]
        return [
            str(self._get_first(branch, ["branch_id", "id", "key"], ""))
            for branch in self._branch_records(reaction_progress)
            if self._get_first(branch, ["branch_id", "id", "key"], "")
        ]

    def _branch_ready_state(self, flag: Any) -> str:
        """把布尔/空值转换为简单状态名。"""
        if flag is True:
            return "COMPLETED"
        if flag is False:
            return "PENDING"
        return "PENDING"

    def _normalize_pipeline_step(self, value: Any) -> str:
        """归一化当前分支所在的流水线步骤。"""
        if value is None:
            return ""
        text = str(value).strip().lower().replace("-", " ").replace("_", " ")
        text = " ".join(text.split())
        return _PIPELINE_STEP_ALIASES.get(text, "")

    def _branch_pipeline_lines(self, branch: Any) -> list[str]:
        """构造分支内部的 S1-S4 渲染行。"""
        overall_state = self._get_attr(branch, "state", None)
        overall_name = self._state_name(overall_state)
        current_step = self._normalize_pipeline_step(self._get_attr(branch, "pipeline_step", None))
        current_state = self._get_attr(branch, "pipeline_step_state", None) or overall_state
        current_state_name = self._state_name(current_state)
        current_detail = str(self._get_attr(branch, "detail", "") or "")
        dg_text = self._format_dg(self._get_attr(branch, "dg_activation", None))

        running_defaults = {
            "s1": "product conformer search",
            "s2": "ts guess generation",
            "s3": "intermediate / TS / SP",
            "s4": "feature extraction",
        }
        complete_defaults = {
            "s1": "product anchor ready",
            "s2": "ts guess ready",
            "s3": dg_text or "branch energetics ready",
            "s4": "features ready",
        }

        current_index: Optional[int] = None
        for index, (step_key, _) in enumerate(_PIPELINE_STEP_ORDER):
            if current_step == step_key:
                current_index = index
                break

        success_terminal = overall_name in {"COMPLETED", "CACHED", "THERMO_DERIVED"}
        failed_terminal = overall_name in {"FAILED", "MISSING_THERMO"}

        lines: list[str] = []
        for index, (step_key, label) in enumerate(_PIPELINE_STEP_ORDER):
            detail = "waiting"
            state: Any = "PENDING"

            if success_terminal:
                state = "COMPLETED"
                detail = complete_defaults.get(step_key, "ready")
                if step_key == "s4" and current_step == "s4" and current_detail:
                    detail = current_detail
            elif current_index is None:
                if overall_name == "RUNNING" and index == 0:
                    state = current_state or "RUNNING"
                    detail = current_detail or running_defaults.get(step_key, "running")
                elif failed_terminal and index == 0:
                    state = overall_state or "FAILED"
                    detail = current_detail or "failed"
            elif index < current_index:
                state = "COMPLETED"
                detail = complete_defaults.get(step_key, "ready")
            elif index == current_index:
                state = current_state or overall_state or "RUNNING"
                if current_state_name in {"COMPLETED", "CACHED"}:
                    detail = complete_defaults.get(step_key, "ready")
                elif current_state_name in {"FAILED", "MISSING_THERMO"}:
                    detail = current_detail or "failed"
                else:
                    detail = current_detail or running_defaults.get(step_key, "running")
            elif failed_terminal and current_index is not None and index > current_index:
                state = "PENDING"
                detail = "waiting"

            lines.append(
                self._format_status_line(
                    label,
                    state,
                    detail,
                    indent="  ",
                    label_width=_BRANCH_LABEL_WIDTH,
                )
            )
        return lines

    def _branch_block(self, progress: Any, branch_id: str) -> str:
        """构造单个分支的完整展示区块。"""
        reaction_progress = self._reaction_progress(progress)
        branch = self._get_branch_record(reaction_progress, branch_id)
        if branch is None:
            return ""

        s0_flag = self._get_attr(branch, "s0_skipped_by_plan", None)
        precursor_flag = self._get_attr(branch, "precursor_linked", None)
        sm_flag = self._get_attr(branch, "sm_refs_materialized", None)

        lines = [self._safe_text(branch_id)]
        lines.append(
            self._format_status_line(
                "S0 mechanism",
                self._branch_ready_state(s0_flag),
                "skipped by scheduler plan" if s0_flag else "waiting",
                label_width=_BRANCH_LABEL_WIDTH,
            )
        )
        lines.append(
            self._format_status_line(
                "precursor",
                self._branch_ready_state(precursor_flag),
                "linked from reaction shared setup" if precursor_flag else "waiting",
                label_width=_BRANCH_LABEL_WIDTH,
            )
        )
        lines.append(
            self._format_status_line(
                "small molecule refs",
                self._branch_ready_state(sm_flag),
                "materialized from global cache" if sm_flag else "waiting",
                label_width=_BRANCH_LABEL_WIDTH,
            )
        )
        lines.extend(self._branch_pipeline_lines(branch))
        return "\n".join(lines)

    def _branch_result_detail(self, branch: Any) -> str:
        """生成 DR 聚合前的单分支摘要。"""
        dg_text = self._format_dg(self._get_attr(branch, "dg_activation", None))
        if dg_text:
            return dg_text
        detail = self._get_attr(branch, "detail", None)
        if detail:
            return str(detail)
        return self._status_label(self._get_attr(branch, "state", None), _GENERIC_STATUS_LABELS)

    def _dr_aggregation_detail(self, reaction_progress: Any) -> str:
        """生成 DR 聚合阶段的一行详情。"""
        detail = self._get_first(
            reaction_progress,
            ["dr_aggregation_detail", "dr_detail", "detail"],
            None,
        )
        if detail:
            return str(detail)
        state_name = self._state_name(self._get_attr(reaction_progress, "dr_aggregation_state", None))
        if state_name in {"COMPLETED", "CACHED"}:
            return "dr_prediction_default.json ready"
        if state_name == "FAILED":
            return "dr aggregation failed"
        return "writing dr_prediction_default.json"

    def _condition_counts(self, progress: Any, reaction_progress: Any) -> tuple[int, int, int]:
        """提取条件后处理统计：conditions / branches / total jobs。"""
        reaction_id = self._reaction_id(progress, reaction_progress)
        plan = self._get_planning_record(progress, reaction_id)

        condition_count = self._get_first(
            reaction_progress,
            ["condition_count", "total_conditions", "conditions_total"],
            None,
        )
        if condition_count is None:
            condition_count = self._get_attr(plan, "condition_count", 0) or 0

        branch_count = self._get_first(
            reaction_progress,
            ["branch_count", "total_branches"],
            None,
        )
        if branch_count is None:
            branch_count = self._get_attr(plan, "branch_count", None)
        if branch_count is None:
            branch_count = len(self._branch_ids(progress, reaction_progress))

        total_jobs = self._get_first(
            reaction_progress,
            [
                "condition_branch_jobs",
                "condition_job_count",
                "total_condition_branch_jobs",
            ],
            None,
        )
        if total_jobs is None:
            try:
                total_jobs = int(condition_count) * int(branch_count)
            except (TypeError, ValueError):
                total_jobs = 0

        try:
            condition_int = int(condition_count)
        except (TypeError, ValueError):
            condition_int = 0
        try:
            branch_int = int(branch_count)
        except (TypeError, ValueError):
            branch_int = 0
        try:
            jobs_int = int(total_jobs)
        except (TypeError, ValueError):
            jobs_int = 0
        return condition_int, branch_int, jobs_int

    def _split_condition_branch_key(self, key: str) -> tuple[str, str]:
        """尽量从 cond_branch_key 中拆出 condition_id 与 branch_id。"""
        raw = str(key)
        for token in (" × ", " x ", " X ", "::", "|", "@", " -> ", "/"):
            if token in raw:
                left, right = raw.split(token, 1)
                return left.strip(), right.strip()
        return raw.strip(), ""

    def _condition_job_containers(self, reaction_progress: Any) -> list[Any]:
        """枚举可能承载 condition-branch 记录的容器。"""
        names = [
            "condition_jobs",
            "condition_branch_jobs",
            "condition_progress",
            "condition_branch_progress",
            "conditions",
        ]
        return [self._get_attr(reaction_progress, name, None) for name in names]

    def _condition_dr_containers(self, reaction_progress: Any) -> list[Any]:
        """枚举可能承载 DR-per-condition 记录的容器。"""
        names = [
            "condition_dr",
            "condition_dr_status",
            "dr_per_condition",
            "condition_dr_progress",
        ]
        return [self._get_attr(reaction_progress, name, None) for name in names]

    def _match_condition_record(
        self, container: Any, cond_branch_key: str, condition_id: str, branch_id: str
    ) -> Any:
        """在单个容器中查找指定的 condition-branch 记录。"""
        if container is None:
            return None

        if isinstance(container, Mapping):
            if cond_branch_key in container:
                return container[cond_branch_key]
            if condition_id and branch_id:
                nested = container.get(condition_id)
                if isinstance(nested, Mapping) and branch_id in nested:
                    return nested[branch_id]
            for key, value in container.items():
                if str(key) == cond_branch_key:
                    return value
                if isinstance(value, Mapping) and condition_id and branch_id and branch_id in value:
                    return value[branch_id]
                record_condition_id = self._get_first(
                    value,
                    ["condition_id", "condition_key", "cond_id", "condition"],
                    "",
                )
                record_branch_id = self._get_first(value, ["branch_id", "branch", "key"], "")
                if str(record_condition_id) == condition_id and str(record_branch_id) == branch_id:
                    return value
            return None

        if isinstance(container, Sequence) and not isinstance(
            container, (str, bytes, bytearray)
        ):
            for value in container:
                record_key = self._get_first(value, ["key", "job_key", "id"], "")
                record_condition_id = self._get_first(
                    value,
                    ["condition_id", "condition_key", "cond_id", "condition"],
                    "",
                )
                record_branch_id = self._get_first(value, ["branch_id", "branch"], "")
                if str(record_key) == cond_branch_key:
                    return value
                if str(record_condition_id) == condition_id and str(record_branch_id) == branch_id:
                    return value
        return None

    def _get_condition_record(self, reaction_progress: Any, cond_branch_key: str) -> Any:
        """查找某个 condition-branch 组合记录。"""
        condition_id, branch_id = self._split_condition_branch_key(cond_branch_key)
        for container in self._condition_job_containers(reaction_progress):
            record = self._match_condition_record(
                container, cond_branch_key, condition_id, branch_id
            )
            if record is not None:
                return record
        return None

    def _iter_condition_dr_entries(self, reaction_progress: Any) -> list[tuple[str, Any]]:
        """把 DR-per-condition 容器统一摊平成 (condition_id, record) 列表。"""
        entries: list[tuple[str, Any]] = []
        for container in self._condition_dr_containers(reaction_progress):
            if container is None:
                continue
            if isinstance(container, Mapping):
                for key, value in container.items():
                    condition_id = self._get_first(
                        value,
                        ["condition_id", "condition_key", "cond_id", "condition"],
                        key,
                    )
                    entries.append((str(condition_id), value))
            elif isinstance(container, Sequence) and not isinstance(
                container, (str, bytes, bytearray)
            ):
                for value in container:
                    condition_id = self._get_first(
                        value,
                        ["condition_id", "condition_key", "cond_id", "condition", "id", "key"],
                        "",
                    )
                    if condition_id:
                        entries.append((str(condition_id), value))
        return entries

    def _stringify_failure_items(self, items: Any) -> list[str]:
        """把失败项列表压缩成适合单行展示的字符串列表。"""
        if items is None:
            return []
        if isinstance(items, Mapping):
            values = list(items.values()) if items else []
            if values:
                return self._stringify_failure_items(values)
            return [str(key) for key in items.keys()]
        if isinstance(items, Sequence) and not isinstance(items, (str, bytes, bytearray)):
            result: list[str] = []
            for item in items:
                if isinstance(item, str):
                    result.append(item)
                    continue
                reaction_id = self._get_first(item, ["reaction_id", "rx_id"], "")
                branch_id = self._get_first(item, ["branch_id", "branch"], "")
                condition_id = self._get_first(
                    item,
                    ["condition_id", "condition_key", "cond_id", "condition"],
                    "",
                )
                key = self._get_first(item, ["key", "id"], "")
                if reaction_id and branch_id:
                    result.append(f"{reaction_id}/{branch_id}")
                elif reaction_id and condition_id:
                    result.append(f"{reaction_id}/{condition_id}")
                elif key:
                    result.append(str(key))
                else:
                    result.append(str(item))
            return result
        return [str(items)]

    def _summary_counts_from_cache(self, cache_progress: Any) -> tuple[int, int, int, int, int]:
        """从 Stage 1 进度对象中提取缓存统计。"""
        molecules = self._sm_records(cache_progress)
        total = self._get_first(cache_progress, ["summary_total", "total_molecules", "total"], None)
        hit = self._get_first(cache_progress, ["summary_hit", "summary_hits", "summary_cached"], None)
        computed = self._get_first(
            cache_progress,
            ["summary_computed", "summary_completed", "summary_done"],
            None,
        )
        failed = self._get_first(cache_progress, ["summary_failed", "failed"], None)
        pending = self._get_first(cache_progress, ["summary_pending", "pending"], None)

        if None not in (total, hit, computed, failed, pending):
            return int(total), int(hit), int(computed), int(failed), int(pending)

        total_int = len(molecules)
        hit_int = 0
        computed_int = 0
        failed_int = 0
        pending_int = 0
        for molecule in molecules:
            state_name = self._state_name(self._get_attr(molecule, "state", None))
            if state_name == "CACHED":
                hit_int += 1
            elif state_name in {"COMPLETED", "THERMO_DERIVED"}:
                computed_int += 1
            elif state_name in {"FAILED", "MISSING_THERMO"}:
                failed_int += 1
            else:
                pending_int += 1
        return total_int, hit_int, computed_int, failed_int, pending_int

    # ------------------------------------------------------------------
    # Stage 0 — 任务构建与机理图规划
    # ------------------------------------------------------------------

    def stage0_planning_header(self, progress: V3RunProgress) -> None:
        """渲染 Stage 0 顶部标题与数据集统计。"""
        planning = self._planning_progress(progress)
        total_rows = self._get_attr(planning, "total_rows", 0) or 0
        total_reactions = self._get_attr(planning, "total_reactions", 0) or 0
        total_conditions = self._get_attr(planning, "total_conditions", 0) or 0
        block = "\n".join(
            [
                f"[bold cyan]{_SEPARATOR}[/]",
                "[bold cyan][RUN STAGE 0] 任务构建与机理图统一规划[/]",
                f"[bold cyan]{_SEPARATOR}[/]",
                self._format_metric_line("Dataset loaded", f"{total_rows} rows"),
                self._format_metric_line("Reaction groups", f"{total_reactions} reactions"),
                self._format_metric_line("Condition tasks", f"{total_conditions} conditions"),
                "Mechanism Planning",
            ]
        )
        self._emit_header_if_changed("_stage0_header_snapshot", block)

    def stage0_planning_update(self, progress: V3RunProgress) -> None:
        """渲染 Stage 0 的逐反应规划更新。"""
        planning = self._planning_progress(progress)
        for record in self._planning_records(planning):
            reaction_id = str(self._get_first(record, ["reaction_id", "id", "key"], ""))
            if not reaction_id:
                continue
            state = self._get_attr(record, "state", None)
            status_text = self._planning_status_text(record)
            branch_count = self._get_attr(record, "branch_count", 0) or 0
            condition_count = self._get_attr(record, "condition_count", 0) or 0
            line = (
                f"  {self._state_markup(state)} {self._pad(reaction_id, 8)}  "
                f"{self._pad(status_text, _REACTION_STATUS_WIDTH)}  "
                f"branches={branch_count}  conditions={condition_count}"
            )
            self._emit_if_changed(self._planning_snapshots, reaction_id, line)

    def stage0_planning_complete(self, progress: V3RunProgress) -> None:
        """渲染 Stage 0 完成提示。"""
        planning = self._planning_progress(progress)
        state = self._get_attr(planning, "state", None) or "COMPLETED"
        state_name = self._state_name(state)
        text = "Planning complete" if state_name != "FAILED" else "Planning complete with failures"
        line = f"{self._state_markup(state, fallback='COMPLETED')} {self._safe_text(text)}"
        self._emit_if_changed(self._completion_snapshots, "stage0", line)

    # ------------------------------------------------------------------
    # Stage 1 — 全局小分子缓存预计算
    # ------------------------------------------------------------------

    def stage1_sm_cache_header(self, progress: V3RunProgress) -> None:
        """渲染 Stage 1 标题与缓存根目录。"""
        cache_progress = self._sm_cache_progress(progress)
        cache_root = self._get_attr(cache_progress, "cache_root", "") or ""
        block = "\n".join(
            [
                f"[bold cyan]{_SEPARATOR}[/]",
                "[bold cyan][RUN STAGE 1] 全局小分子缓存预计算[/]",
                f"[bold cyan]{_SEPARATOR}[/]",
                f"Cache root: {self._safe_text(cache_root)}",
                "Small Molecule Cache",
            ]
        )
        self._emit_header_if_changed("_stage1_header_snapshot", block)

    def stage1_sm_cache_update(self, progress: V3RunProgress) -> None:
        """渲染 Stage 1 的逐小分子缓存更新。"""
        cache_progress = self._sm_cache_progress(progress)
        for record in self._sm_records(cache_progress):
            key = str(self._get_first(record, ["key", "sm_key", "id"], ""))
            if not key:
                continue
            state = self._get_attr(record, "state", None)
            status_text = self._cache_status_text(record)
            detail = str(self._get_attr(record, "detail", "") or "")
            line = (
                f"  {self._state_markup(state)} {self._pad(key, _CACHE_KEY_WIDTH)}  "
                f"{self._pad(status_text, _CACHE_STATUS_WIDTH)}  {self._safe_text(detail)}"
            ).rstrip()
            self._emit_if_changed(self._cache_snapshots, key, line)

    def stage1_sm_cache_complete(self, progress: V3RunProgress) -> None:
        """渲染 Stage 1 汇总信息。"""
        cache_progress = self._sm_cache_progress(progress)
        total, hit, computed, failed, pending = self._summary_counts_from_cache(
            cache_progress
        )
        block = "\n".join(
            [
                "Summary",
                f"  total: {total} | hit: {hit} | computed: {computed} | failed: {failed} | pending: {pending}",
            ]
        )
        self._emit_if_changed(self._completion_snapshots, "stage1", block)

    # ------------------------------------------------------------------
    # Stage 2 — Reaction Execution
    # ------------------------------------------------------------------

    def stage2_reaction_header(self, progress: V3RunProgress) -> None:
        """渲染 Stage 2 当前反应头部。"""
        reaction_progress = self._reaction_progress(progress)
        reaction_id = self._reaction_id(progress, reaction_progress)
        current_index, total_reactions = self._reaction_index_total(progress)
        block = "\n".join(
            [
                f"[bold cyan]{_SEPARATOR}[/]",
                "[bold cyan][RUN STAGE 2] Reaction Execution[/]",
                f"Current: {self._safe_text(reaction_id)}  [{current_index}/{total_reactions}]",
                f"[bold cyan]{_SEPARATOR}[/]",
            ]
        )
        self._emit_header_if_changed("_stage2_header_snapshot", block)

    def stage2_reaction_setup(self, progress: V3RunProgress) -> None:
        """渲染 Stage 2 共享反应准备区块。"""
        reaction_progress = self._reaction_progress(progress)
        reaction_id = self._reaction_id(progress, reaction_progress)
        plan = self._get_planning_record(progress, reaction_id)

        mechanism_state = self._get_first(
            reaction_progress,
            ["mechanism_graph_state", "shared_mechanism_state"],
            self._get_attr(plan, "state", "PENDING"),
        )
        mechanism_detail = str(
            self._get_first(
                reaction_progress,
                ["mechanism_graph_detail", "shared_mechanism_detail"],
                self._planning_to_setup_detail(plan),
            )
        )

        branch_ids = self._branch_ids(progress, reaction_progress)
        branch_state: Any = "COMPLETED" if branch_ids else "PENDING"
        branch_detail = ", ".join(branch_ids) if branch_ids else "pending"

        precursor_state = self._get_attr(reaction_progress, "shared_precursor_s1_state", None) or "PENDING"
        precursor_detail = str(
            self._get_first(
                reaction_progress,
                ["shared_precursor_s1_detail", "shared_precursor_detail"],
                self._default_shared_precursor_detail(precursor_state),
            )
        )

        sm_state = self._get_first(
            reaction_progress,
            ["small_molecule_refs_state", "shared_sm_cache_state", "small_molecule_cache_state"],
            "PENDING",
        )
        sm_detail = str(
            self._get_first(
                reaction_progress,
                ["small_molecule_refs_detail", "shared_sm_cache_detail", "small_molecule_cache_detail"],
                self._default_sm_cache_detail(sm_state),
            )
        )

        block = "\n".join(
            [
                f"[{self._safe_text(reaction_id)}] Shared Reaction Setup",
                self._format_status_line(
                    "mechanism graph",
                    mechanism_state,
                    mechanism_detail,
                    label_width=_SETUP_LABEL_WIDTH,
                ),
                self._format_status_line(
                    "branch sequence",
                    branch_state,
                    branch_detail,
                    label_width=_SETUP_LABEL_WIDTH,
                ),
                self._format_status_line(
                    "precursor S1 once",
                    precursor_state,
                    precursor_detail,
                    label_width=_SETUP_LABEL_WIDTH,
                ),
                self._format_status_line(
                    "small molecule cache refs",
                    sm_state,
                    sm_detail,
                    label_width=_SETUP_LABEL_WIDTH,
                ),
            ]
        )
        self._emit_if_changed(self._reaction_setup_snapshots, reaction_id, block)

    def stage2_branch_start(self, progress: V3RunProgress, branch_id: str) -> None:
        """渲染某个分支开始时的快照。"""
        self._render_branch(progress, branch_id)

    def stage2_branch_update(self, progress: V3RunProgress, branch_id: str) -> None:
        """渲染某个分支运行中的快照。"""
        self._render_branch(progress, branch_id)

    def stage2_branch_complete(self, progress: V3RunProgress, branch_id: str) -> None:
        """渲染某个分支完成时的快照。"""
        self._render_branch(progress, branch_id)

    def _render_branch(self, progress: Any, branch_id: str) -> None:
        """统一处理分支块渲染。"""
        reaction_id = self._reaction_id(progress)
        block = self._branch_block(progress, branch_id)
        self._emit_if_changed(self._branch_snapshots, (reaction_id, branch_id), block)

    def stage2_dr_aggregation(self, progress: V3RunProgress) -> None:
        """渲染分支结果汇总与 DR 聚合状态。"""
        reaction_progress = self._reaction_progress(progress)
        reaction_id = self._reaction_id(progress, reaction_progress)
        lines = [f"[{self._safe_text(reaction_id)}] Branch Results"]
        for branch in self._branch_records(reaction_progress):
            branch_id = str(self._get_first(branch, ["branch_id", "id", "key"], ""))
            if not branch_id:
                continue
            detail = self._branch_result_detail(branch)
            lines.append(
                f"  {self._state_markup(self._get_attr(branch, 'state', None))} "
                f"{self._pad(branch_id, 14)}{self._safe_text(detail)}"
            )
        lines.append(f"[{self._safe_text(reaction_id)}] DR Aggregation")
        lines.append(
            self._format_status_line(
                "",
                self._get_attr(reaction_progress, "dr_aggregation_state", None) or "RUNNING",
                self._dr_aggregation_detail(reaction_progress),
                indent="  ",
                label_width=0,
            )
        )
        block = "\n".join(lines)
        self._emit_if_changed(self._dr_snapshots, reaction_id, block)

    def stage2_condition_header(self, progress: V3RunProgress) -> None:
        """渲染条件后处理区块头部。"""
        reaction_progress = self._reaction_progress(progress)
        reaction_id = self._reaction_id(progress, reaction_progress)
        condition_count, branch_count, total_jobs = self._condition_counts(
            progress, reaction_progress
        )
        block = "\n".join(
            [
                f"[{self._safe_text(reaction_id)}] Condition Post-processing",
                f"  conditions: {condition_count}",
                f"  branches: {branch_count}",
                f"  total condition-branch jobs: {total_jobs}",
            ]
        )
        self._emit_if_changed(self._condition_header_snapshots, reaction_id, block)

    def stage2_condition_update(
        self, progress: V3RunProgress, cond_branch_key: str
    ) -> None:
        """渲染某个 condition-branch 组合任务的更新。"""
        reaction_progress = self._reaction_progress(progress)
        reaction_id = self._reaction_id(progress, reaction_progress)
        record = self._get_condition_record(reaction_progress, cond_branch_key)
        if record is not None:
            condition_id = str(
                self._get_first(
                    record,
                    ["condition_id", "condition_key", "cond_id", "condition"],
                    self._split_condition_branch_key(cond_branch_key)[0],
                )
            )
            branch_id = str(
                self._get_first(
                    record,
                    ["branch_id", "branch"],
                    self._split_condition_branch_key(cond_branch_key)[1],
                )
            )
            label = f"{condition_id} × {branch_id}" if branch_id else condition_id
            detail = str(
                self._get_attr(record, "detail", None)
                or self._status_label(self._get_attr(record, "state", None), _GENERIC_STATUS_LABELS)
            )
            line = (
                f"  {self._state_markup(self._get_attr(record, 'state', None))} "
                f"{self._pad(label, _CONDITION_LABEL_WIDTH)}{self._safe_text(detail)}"
            )
            self._emit_if_changed(
                self._condition_line_snapshots,
                (reaction_id, cond_branch_key),
                line,
            )
        self._render_condition_dr_updates(reaction_id, reaction_progress)

    def _render_condition_dr_updates(self, reaction_id: str, reaction_progress: Any) -> None:
        """渲染 DR per condition 子区块中的变化项。"""
        entries = self._iter_condition_dr_entries(reaction_progress)
        if not entries:
            return
        if reaction_id not in self._condition_dr_section_printed:
            self._condition_dr_section_printed.add(reaction_id)
            self._print("  DR per condition:")
        for condition_id, record in entries:
            detail = str(
                self._get_attr(record, "detail", None)
                or self._status_label(self._get_attr(record, "state", None), _GENERIC_STATUS_LABELS)
            )
            line = (
                f"    {self._state_markup(self._get_attr(record, 'state', None))} "
                f"{self._pad(condition_id, 10)}{self._safe_text(detail)}"
            )
            self._emit_if_changed(
                self._condition_dr_snapshots,
                (reaction_id, condition_id),
                line,
            )

    def stage2_reaction_complete(self, progress: V3RunProgress) -> None:
        """渲染反应收尾。"""
        reaction_progress = self._reaction_progress(progress)
        reaction_id = self._reaction_id(progress, reaction_progress)
        self._render_condition_dr_updates(reaction_id, reaction_progress)
        if reaction_id and reaction_id not in self._reaction_complete_seen:
            self._reaction_complete_seen.add(reaction_id)

    # ------------------------------------------------------------------
    # Stage 3 — Final Summary
    # ------------------------------------------------------------------

    def stage3_summary(self, progress: V3RunProgress) -> None:
        """渲染最终运行摘要。"""
        summary = self._summary_progress(progress)
        planning = self._planning_progress(progress)
        cache_progress = self._sm_cache_progress(progress)

        total_reactions = int(self._get_attr(summary, "total_reactions", 0) or 0)
        reaction_success = int(self._get_attr(summary, "reaction_success", 0) or 0)
        reaction_failed = int(self._get_attr(summary, "reaction_failed", 0) or 0)

        total_branches = int(self._get_attr(summary, "total_branches", 0) or 0)
        branch_success = int(self._get_attr(summary, "branch_success", 0) or 0)
        branch_failed = int(self._get_attr(summary, "branch_failed", 0) or 0)

        condition_total = self._get_first(
            summary,
            ["total_conditions", "condition_total", "conditions_total"],
            self._get_attr(planning, "total_conditions", 0) or 0,
        )
        condition_complete = self._get_first(
            summary,
            ["condition_complete", "condition_success", "conditions_complete"],
            None,
        )
        condition_failed = self._get_first(
            summary,
            ["condition_failed", "conditions_failed"],
            None,
        )
        if condition_failed is None:
            condition_failed = len(
                self._stringify_failure_items(self._get_attr(summary, "failed_conditions", []))
            )
        if condition_complete is None:
            try:
                condition_complete = int(condition_total) - int(condition_failed)
            except (TypeError, ValueError):
                condition_complete = 0

        sm_cached = self._get_first(
            summary,
            ["sm_cache_hit", "sm_cache_cached", "sm_cache_hits"],
            None,
        )
        sm_computed = self._get_first(
            summary,
            ["sm_cache_computed", "sm_cache_completed", "sm_cache_done"],
            None,
        )
        sm_failed = self._get_first(
            summary,
            ["sm_cache_failed", "sm_cache_failures"],
            None,
        )
        if None in (sm_cached, sm_computed, sm_failed):
            _, sm_cached, sm_computed, sm_failed, _ = self._summary_counts_from_cache(
                cache_progress
            )

        failed_planning = self._stringify_failure_items(
            self._get_attr(summary, "failed_planning", [])
        )
        failed_cache = self._stringify_failure_items(
            self._get_attr(summary, "failed_cache", [])
        )
        failed_branches = self._stringify_failure_items(
            self._get_attr(summary, "failed_branches", [])
        )
        failed_conditions = self._stringify_failure_items(
            self._get_attr(summary, "failed_conditions", [])
        )

        lines = [
            f"[bold cyan]{_SEPARATOR}[/]",
            "[bold cyan][RUN STAGE 3] Final Run Summary[/]",
            f"[bold cyan]{_SEPARATOR}[/]",
            (
                f"{self._pad('Reactions', _SUMMARY_LABEL_WIDTH)}"
                f"{reaction_success}/{total_reactions} succeeded   {reaction_failed} failed"
            ),
            (
                f"{self._pad('Branches', _SUMMARY_LABEL_WIDTH)}"
                f"{branch_success}/{total_branches} succeeded   {branch_failed} failed"
            ),
            (
                f"{self._pad('Conditions', _SUMMARY_LABEL_WIDTH)}"
                f"{int(condition_complete)}/{int(condition_total)} complete  {int(condition_failed)} failed"
            ),
            (
                f"{self._pad('Small Molecules', _SUMMARY_LABEL_WIDTH)}"
                f"{int(sm_cached)} cached | {int(sm_computed)} computed | {int(sm_failed)} failed"
            ),
        ]

        if failed_planning:
            lines.append(
                f"{self._pad('Failed Planning:', 18)}{self._safe_text(', '.join(failed_planning))}"
            )
        if failed_cache:
            lines.append(
                f"{self._pad('Failed Cache:', 18)}{self._safe_text(', '.join(failed_cache))}"
            )
        if failed_branches:
            lines.append(
                f"{self._pad('Failed Branches:', 18)}{self._safe_text(', '.join(failed_branches))}"
            )
        if failed_conditions:
            lines.append(
                f"{self._pad('Failed Conditions:', 18)}{self._safe_text(', '.join(failed_conditions))}"
            )

        block = "\n".join(lines)
        self._emit_header_if_changed("_stage3_summary_snapshot", block)


class SilentV3StageDisplay(V3StageDisplay):
    """V3 静默显示器：所有渲染钩子均为空操作。"""

    def stage0_planning_header(self, progress: V3RunProgress) -> None:
        """静默模式下不输出 Stage 0 头部。"""
        pass

    def stage0_planning_update(self, progress: V3RunProgress) -> None:
        """静默模式下不输出 Stage 0 更新。"""
        pass

    def stage0_planning_complete(self, progress: V3RunProgress) -> None:
        """静默模式下不输出 Stage 0 完成提示。"""
        pass

    def stage1_sm_cache_header(self, progress: V3RunProgress) -> None:
        """静默模式下不输出 Stage 1 头部。"""
        pass

    def stage1_sm_cache_update(self, progress: V3RunProgress) -> None:
        """静默模式下不输出 Stage 1 更新。"""
        pass

    def stage1_sm_cache_complete(self, progress: V3RunProgress) -> None:
        """静默模式下不输出 Stage 1 汇总。"""
        pass

    def stage2_reaction_header(self, progress: V3RunProgress) -> None:
        """静默模式下不输出 Stage 2 头部。"""
        pass

    def stage2_reaction_setup(self, progress: V3RunProgress) -> None:
        """静默模式下不输出共享反应准备区块。"""
        pass

    def stage2_branch_start(self, progress: V3RunProgress, branch_id: str) -> None:
        """静默模式下不输出分支开始快照。"""
        pass

    def stage2_branch_update(self, progress: V3RunProgress, branch_id: str) -> None:
        """静默模式下不输出分支运行快照。"""
        pass

    def stage2_branch_complete(self, progress: V3RunProgress, branch_id: str) -> None:
        """静默模式下不输出分支完成快照。"""
        pass

    def stage2_dr_aggregation(self, progress: V3RunProgress) -> None:
        """静默模式下不输出 DR 聚合区块。"""
        pass

    def stage2_condition_header(self, progress: V3RunProgress) -> None:
        """静默模式下不输出条件后处理头部。"""
        pass

    def stage2_condition_update(
        self, progress: V3RunProgress, cond_branch_key: str
    ) -> None:
        """静默模式下不输出条件任务更新。"""
        pass

    def stage2_reaction_complete(self, progress: V3RunProgress) -> None:
        """静默模式下不输出反应收尾。"""
        pass

    def stage3_summary(self, progress: V3RunProgress) -> None:
        """静默模式下不输出最终摘要。"""
        pass


__all__ = [
    "V3StageDisplay",
    "LoggerV3StageDisplay",
    "SilentV3StageDisplay",
]
