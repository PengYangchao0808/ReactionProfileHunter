"""
Task Progress System — V4 任务导向进度追踪

替代传统的 S0-S4 线性阶段模型，使用任务 DAG 追踪每个原子计算操作的状态和进度。

核心概念:
- TaskSpec: 任务定义 (id, name, description, dependencies, weight)
- TaskState: 任务状态枚举
- TaskProgressTracker: 任务状态机 + 进度追踪
- TaskRegistry: 预定义 V4 任务清单
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum, auto
from typing import Any, Dict, List, Optional, Set, Tuple


class TaskState(Enum):
    """任务状态机。"""

    PENDING = auto()      # 等待依赖完成
    QUEUED = auto()       # 依赖已满足，等待资源
    RUNNING = auto()      # 正在执行
    COMPLETED = auto()    # 成功完成
    FAILED = auto()       # 执行失败
    CACHED = auto()       # 从 checkpoint 加载
    SKIPPED = auto()      # 被配置跳过

    @property
    def icon(self) -> str:
        return {
            TaskState.PENDING: " ",
            TaskState.QUEUED: "○",
            TaskState.RUNNING: "→",
            TaskState.COMPLETED: "✓",
            TaskState.FAILED: "✗",
            TaskState.CACHED: "↻",
            TaskState.SKIPPED: "⊘",
        }[self]

    @property
    def color(self) -> str:
        return {
            TaskState.PENDING: "dim",
            TaskState.QUEUED: "yellow",
            TaskState.RUNNING: "cyan",
            TaskState.COMPLETED: "green",
            TaskState.FAILED: "red",
            TaskState.CACHED: "blue",
            TaskState.SKIPPED: "dim",
        }[self]

    @property
    def label(self) -> str:
        return {
            TaskState.PENDING: "pending",
            TaskState.QUEUED: "queued",
            TaskState.RUNNING: "running",
            TaskState.COMPLETED: "completed",
            TaskState.FAILED: "failed",
            TaskState.CACHED: "cached",
            TaskState.SKIPPED: "skipped",
        }[self]


@dataclass
class TaskSpec:
    """任务定义。"""

    task_id: str
    name: str
    description: str = ""
    phase: str = ""           # 所属阶段: analysis/anchor/scan/ts/energy/features/post
    depends_on: list[str] = field(default_factory=list)
    weight: int = 10          # 相对权重 (影响总进度计算)
    cacheable: bool = True
    optional: bool = False


@dataclass
class TaskRecord:
    """运行时任务记录。"""

    spec: TaskSpec
    state: TaskState = TaskState.PENDING
    progress_pct: int = 0     # 0-100
    detail: str = ""          # 当前详细状态 (如 "Berny iter 12/20")
    elapsed_sec: float = 0.0
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    result_summary: str = ""  # 完成后的一行摘要 (如 "E_sp = -456.12 Ha")
    error_message: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def is_terminal(self) -> bool:
        return self.state in (
            TaskState.COMPLETED,
            TaskState.FAILED,
            TaskState.CACHED,
            TaskState.SKIPPED,
        )

    @property
    def display_name(self) -> str:
        return f"{self.spec.name}"


# ============================================================================
# V4 预定义任务清单
# ============================================================================

V4_TASK_REGISTRY: dict[str, TaskSpec] = {
    # === Analysis Phase ===
    "mechanism": TaskSpec(
        task_id="mechanism",
        name="Mechanism Analysis",
        description="Reaction type classification and topology analysis",
        phase="analysis",
        weight=1,
        cacheable=True,
    ),

    # === Anchor Phase ===
    "product_anchor": TaskSpec(
        task_id="product_anchor",
        name="Product Anchor",
        description="Conformer search + DFT opt + SP for product",
        phase="anchor",
        depends_on=["mechanism"],
        weight=15,
        cacheable=True,
    ),
    "precursor_anchor": TaskSpec(
        task_id="precursor_anchor",
        name="Precursor Anchor",
        description="Conformer search + DFT opt + SP for precursor",
        phase="anchor",
        depends_on=["mechanism"],
        weight=15,
        cacheable=True,
        optional=True,
    ),
    "smallmol_anchor": TaskSpec(
        task_id="smallmol_anchor",
        name="Small Molecule Anchor",
        description="Conformer search + DFT opt + SP for small molecules",
        phase="anchor",
        depends_on=["mechanism"],
        weight=10,
        cacheable=True,
        optional=True,
    ),

    # === Scan Phase ===
    "retro_scan": TaskSpec(
        task_id="retro_scan",
        name="Retro Scan",
        description="xTB inward scan from product to TS guess",
        phase="scan",
        depends_on=["product_anchor"],
        weight=8,
        cacheable=True,
    ),

    # === TS Phase ===
    "ts_opt": TaskSpec(
        task_id="ts_opt",
        name="TS Optimization",
        description="Transition state optimization (Berny/QST2)",
        phase="ts",
        depends_on=["retro_scan"],
        weight=20,
        cacheable=True,
    ),
    "irc_verify": TaskSpec(
        task_id="irc_verify",
        name="IRC Verification",
        description="IRC path verification",
        phase="ts",
        depends_on=["ts_opt"],
        weight=12,
        cacheable=True,
    ),
    "reactant_opt": TaskSpec(
        task_id="reactant_opt",
        name="Reactant Optimization",
        description="Reactant geometry optimization",
        phase="ts",
        depends_on=["ts_opt"],
        weight=15,
        cacheable=True,
    ),

    # === Energy Phase ===
    "sp_matrix": TaskSpec(
        task_id="sp_matrix",
        name="SP Matrix",
        description="Single point energy matrix (TS, reactant, product)",
        phase="energy",
        depends_on=["ts_opt", "reactant_opt"],
        weight=18,
        cacheable=True,
    ),
    "thermochemistry": TaskSpec(
        task_id="thermochemistry",
        name="Thermochemistry",
        description="Frequency analysis and thermochemistry",
        phase="energy",
        depends_on=["reactant_opt"],
        weight=12,
        cacheable=True,
    ),

    # === Features Phase ===
    "geom_features": TaskSpec(
        task_id="geom_features",
        name="Geometry Features",
        description="Geometric parameter extraction",
        phase="features",
        depends_on=["product_anchor", "ts_opt"],
        weight=2,
        cacheable=True,
    ),
    "elec_features": TaskSpec(
        task_id="elec_features",
        name="Electronic Features",
        description="Electronic descriptor extraction (CDFT, NBO, etc.)",
        phase="features",
        depends_on=["product_anchor", "ts_opt", "sp_matrix"],
        weight=3,
        cacheable=True,
    ),
    "thermo_features": TaskSpec(
        task_id="thermo_features",
        name="Thermo Features",
        description="Thermochemistry feature extraction",
        phase="features",
        depends_on=["thermochemistry"],
        weight=2,
        cacheable=True,
    ),
    "nbo_features": TaskSpec(
        task_id="nbo_features",
        name="NBO Analysis",
        description="NBO population and interaction analysis",
        phase="features",
        depends_on=["reactant_opt"],
        weight=5,
        cacheable=True,
        optional=True,
    ),

    # === Post Phase ===
    "dr_aggregate": TaskSpec(
        task_id="dr_aggregate",
        name="DR Aggregation",
        description="Diastereomer ratio aggregation",
        phase="post",
        depends_on=["geom_features", "elec_features", "thermo_features"],
        weight=1,
        cacheable=False,
    ),
    "condition_thermo": TaskSpec(
        task_id="condition_thermo",
        name="Condition Thermo",
        description="Condition-specific thermochemistry calculations",
        phase="post",
        depends_on=["dr_aggregate"],
        weight=1,
        cacheable=False,
        optional=True,
    ),
}


# ============================================================================
# TaskProgressTracker — 任务状态机
# ============================================================================

class TaskProgressTracker:
    """追踪一组任务的状态和进度。"""

    def __init__(self, task_specs: Optional[dict[str, TaskSpec]] = None) -> None:
        self.tasks: dict[str, TaskRecord] = {}
        self._task_order: list[str] = []
        self._start_time: Optional[datetime] = None

        specs = task_specs or V4_TASK_REGISTRY
        for task_id, spec in specs.items():
            self.tasks[task_id] = TaskRecord(spec=spec)
            self._task_order.append(task_id)

    def start(self) -> None:
        self._start_time = datetime.now()

    @property
    def elapsed_sec(self) -> float:
        if self._start_time is None:
            return 0.0
        return (datetime.now() - self._start_time).total_seconds()

    def set_state(self, task_id: str, state: TaskState) -> None:
        if task_id not in self.tasks:
            return
        task = self.tasks[task_id]
        now = datetime.now()
        if state == TaskState.RUNNING and task.start_time is None:
            task.start_time = now
        if state in (
            TaskState.COMPLETED,
            TaskState.FAILED,
            TaskState.CACHED,
            TaskState.SKIPPED,
        ):
            task.end_time = now
            if task.start_time:
                task.elapsed_sec = (now - task.start_time).total_seconds()
        task.state = state

    def set_progress(self, task_id: str, pct: int, detail: str = "") -> None:
        if task_id not in self.tasks:
            return
        self.tasks[task_id].progress_pct = max(0, min(100, pct))
        if detail:
            self.tasks[task_id].detail = detail

    def set_result_summary(self, task_id: str, summary: str) -> None:
        if task_id not in self.tasks:
            return
        self.tasks[task_id].result_summary = summary

    def set_detail(self, task_id: str, detail: str) -> None:
        if task_id not in self.tasks:
            return
        self.tasks[task_id].detail = detail

    def set_error(self, task_id: str, message: str) -> None:
        if task_id not in self.tasks:
            return
        self.tasks[task_id].error_message = message
        self.tasks[task_id].state = TaskState.FAILED

    def get_ready_tasks(self) -> list[str]:
        """返回所有依赖已满足且状态为 PENDING 的任务。"""
        ready: list[str] = []
        for task_id, task in self.tasks.items():
            if task.state != TaskState.PENDING:
                continue
            deps_satisfied = all(
                self.tasks.get(dep, TaskRecord(spec=TaskSpec("", ""))).is_terminal
                for dep in task.spec.depends_on
            )
            if deps_satisfied:
                ready.append(task_id)
        return ready

    @property
    def overall_progress(self) -> int:
        """基于任务权重计算总体进度。"""
        total_weight = sum(t.spec.weight for t in self.tasks.values() if not t.spec.optional)
        if total_weight == 0:
            return 0
        completed_weight = 0
        for task in self.tasks.values():
            if task.spec.optional and task.state == TaskState.PENDING:
                continue
            if task.is_terminal:
                completed_weight += task.spec.weight
            elif task.state == TaskState.RUNNING:
                completed_weight += task.spec.weight * task.progress_pct / 100
        return int((completed_weight / total_weight) * 100)

    @property
    def completed_count(self) -> int:
        return sum(1 for t in self.tasks.values() if t.state == TaskState.COMPLETED)

    @property
    def cached_count(self) -> int:
        return sum(1 for t in self.tasks.values() if t.state == TaskState.CACHED)

    @property
    def failed_count(self) -> int:
        return sum(1 for t in self.tasks.values() if t.state == TaskState.FAILED)

    @property
    def running_count(self) -> int:
        return sum(1 for t in self.tasks.values() if t.state == TaskState.RUNNING)

    @property
    def total_count(self) -> int:
        return len(self.tasks)

    @property
    def active_tasks(self) -> list[TaskRecord]:
        """返回当前活跃（RUNNING 或 QUEUED）的任务。"""
        return [
            t for t in self.tasks.values()
            if t.state in (TaskState.RUNNING, TaskState.QUEUED)
        ]

    @property
    def failed_tasks(self) -> list[TaskRecord]:
        return [t for t in self.tasks.values() if t.state == TaskState.FAILED]

    def get_task_order(self) -> list[str]:
        """按拓扑顺序返回任务 ID 列表。"""
        visited: set[str] = set()
        order: list[str] = []

        def visit(task_id: str) -> None:
            if task_id in visited or task_id not in self.tasks:
                return
            visited.add(task_id)
            for dep in self.tasks[task_id].spec.depends_on:
                visit(dep)
            order.append(task_id)

        for task_id in self._task_order:
            visit(task_id)
        return order

    def get_phase_tasks(self, phase: str) -> list[TaskRecord]:
        """返回指定阶段的任务。"""
        return [t for t in self.tasks.values() if t.spec.phase == phase]

    def get_summary(self) -> str:
        """生成一行摘要。"""
        parts: list[str] = []
        total = self.total_count
        completed = self.completed_count + self.cached_count
        failed = self.failed_count
        running = self.running_count

        parts.append(f"{completed}/{total} tasks")
        if running > 0:
            parts.append(f"{running} running")
        if failed > 0:
            parts.append(f"{failed} failed")
        if self.elapsed_sec > 0:
            parts.append(f"⏱ {self.elapsed_sec:.1f}s")
        return " | ".join(parts)
