"""
UI管理模块 - 支持两种模式：
1. 日志模式 (默认): 顺序打印关键流程，清晰可读
2. Rich模式: 保留动态进度面板 (适合需要实时监控的场景)

通过环境变量 RPH_UI_MODE 控制:
- RPH_UI_MODE=log    (默认) - 简洁日志模式
- RPH_UI_MODE=rich   - 动态进度模式
- RPH_UI_MODE=none   - 完全静默，只保留文件日志

其他环境变量:
- RPH_NO_PROGRESS=1  - 禁用任何进度显示
- RPH_ALT_SCREEN=1   - Rich模式下使用alternate screen (默认关闭)
"""

import os
import platform
import re
import time
from collections import deque
from contextlib import contextmanager
from datetime import datetime
from typing import Any, Deque, Dict, Optional, cast

# Rich imports - only needed for rich mode
from rich.console import Console
from rich.layout import Layout
from rich.live import Live
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskID,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.rule import Rule
from rich.table import Table

from rph_core.utils.shared_console import get_console

# TYPE_CHECKING guard to avoid circular import at runtime
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from rph_core.utils.task_progress import TaskProgressTracker

console: Console = get_console()

# Global singleton
_shared_ui_manager: Optional['BaseUIManager'] = None


def get_progress_manager() -> 'BaseUIManager':
    """Factory function that returns appropriate UI manager based on config."""
    global _shared_ui_manager
    if _shared_ui_manager is None:
        mode = os.environ.get("RPH_UI_MODE", "log").lower().strip()
        if mode == "rich":
            _shared_ui_manager = RichProgressManager()
        elif mode == "none":
            _shared_ui_manager = SilentProgressManager()
        else:  # default: log
            _shared_ui_manager = LoggerProgressManager()
    return _shared_ui_manager


# ============================================================================
# ErrorFormatter — 将原始 subprocess 错误压缩为结构化摘要
# ============================================================================

class ErrorFormatter:
    """将原始异常 + subprocess stderr 压缩为结构化摘要。

    过滤低价值噪音（SIGSEGV backtrace、空行），
    提取关键错误信息并限制输出长度。
    """

    SIGSEGV_PATTERN: re.Pattern = re.compile(
        r"Program received signal SIGSEGV.*?"
        r"Backtrace for this error:\s*",
        re.DOTALL,
    )

    # 匹配 CREST/xTB 输出中的有用错误行
    USEFUL_ERROR_PATTERNS: list[re.Pattern] = [
        re.compile(r"(error|ERROR|Error|FAILED|failed|FATAL|abort|abnormal).*", re.I),
        re.compile(r"recieved SIGINT.*"),
        re.compile(r"normally terminated.*"),
    ]

    MAX_STDERR_LINES: int = 5
    MAX_LINE_LENGTH: int = 200

    @classmethod
    def compact(cls, exc: Exception, stderr: str = "", stdout: str = "") -> str:
        """将异常 + stderr 压缩为简洁的错误摘要。

        Args:
            exc: 原始异常对象
            stderr: 子进程 stderr 输出
            stdout: 子进程 stdout 输出（可选，用于提取最后的正常输出）

        Returns:
            压缩后的错误字符串（通常 3-6 行）
        """
        parts: list[str] = []

        # 1. 异常类型和消息
        exc_msg = str(exc)
        parts.append(f"{type(exc).__name__}: {exc_msg[:300]}")

        # 2. 清理 stderr — 去掉 SIGSEGV backtrace
        clean_stderr = stderr if stderr else ""
        clean_stderr = cls.SIGSEGV_PATTERN.sub("", clean_stderr)

        # 提取有用的错误行
        useful_lines: list[str] = []
        for line in clean_stderr.split("\n"):
            stripped = line.strip()
            if not stripped:
                continue
            if any(p.search(stripped) for p in cls.USEFUL_ERROR_PATTERNS):
                useful_lines.append(stripped[: cls.MAX_LINE_LENGTH])

        if useful_lines:
            for line in useful_lines[: cls.MAX_STDERR_LINES]:
                parts.append(f"  └ {line}")

        # 3. 如果没有 stderr 信息，检查 stdout 最后几行
        if not useful_lines and stdout:
            stdout_lines = stdout.strip().split("\n")
            last_lines = [l.strip() for l in stdout_lines[-3:] if l.strip()]
            if last_lines:
                parts.append("  stdout tail:")
                for line in last_lines:
                    parts.append(f"    {line[: cls.MAX_LINE_LENGTH]}")

        return "\n".join(parts)


# ============================================================================
# BaseUIManager — 统一接口
# ============================================================================

class BaseUIManager:
    """Abstract base class defining UI manager interface."""

    def __init__(self) -> None:
        self._batch_mode: bool = False
        self._reaction_id: str = ""
        self._branch_id: str = ""

    def set_context(
        self,
        batch_mode: bool = False,
        reaction_id: str = "",
        branch_id: str = "",
    ) -> None:
        """设置当前运行上下文（用于 batch 模式下的层级显示）。

        Args:
            batch_mode: 是否处于批量运行模式
            reaction_id: 当前反应 ID
            branch_id: 当前分支 ID
        """
        self._batch_mode = batch_mode
        self._reaction_id = reaction_id
        self._branch_id = branch_id

    def start(self, title: str = "Processing...") -> None:
        """Initialize the UI."""
        raise NotImplementedError

    def stop(self) -> None:
        """Clean up and finalize the UI."""
        raise NotImplementedError

    def add_step(self, step_id: str, description: str, total: int = 100) -> None:
        """Register a new pipeline step."""
        raise NotImplementedError

    def update_step(
        self,
        step_id: str,
        completed: Optional[int] = None,
        description: Optional[str] = None,
        advance: Optional[int] = None,
    ) -> None:
        """Update step progress."""
        raise NotImplementedError

    def enter_phase(self, step_id: str, phase_name: str) -> None:
        """Mark entry into a phase within a step."""
        raise NotImplementedError

    def set_subtask(
        self, step_id: str, name: str, current: int, total: int
    ) -> None:
        """Update subtask progress."""
        raise NotImplementedError

    def log_event(self, step_id: str, message: str) -> None:
        """Log an event/message."""
        raise NotImplementedError

    def section_header(
        self, step_id: str, title: str, description: str = ""
    ) -> None:
        """统一的步骤标题渲染 — 替代 print_step_header()."""
        raise NotImplementedError

    def pipeline_banner(
        self, version: str, reaction_id: str = "", branch_id: str = ""
    ) -> None:
        """统一的 pipeline 启动横幅 — 替代 print_pipeline_header() + pm.start()."""
        raise NotImplementedError

    def error_block(
        self,
        step_id: str,
        error: Exception,
        context: Optional[Dict[str, Any]] = None,
    ) -> None:
        """结构化错误输出 — 替代裸 traceback."""
        raise NotImplementedError

    def result_summary(self, result: Any) -> None:
        """统一结果摘要 — 替代 print_result_summary()."""
        raise NotImplementedError

    def render_tasks(self, tracker: "TaskProgressTracker") -> None:
        """渲染任务树 — V4 任务导向显示."""
        raise NotImplementedError

    @contextmanager
    def manage(self, title: str = "Pipeline Progress"):
        """Context manager for easy usage."""
        try:
            self.start(title)
            yield self
        finally:
            self.stop()


# ============================================================================
# LoggerProgressManager — 默认日志模式
# ============================================================================

class LoggerProgressManager(BaseUIManager):
    """
    简洁日志式UI - 默认模式

    特点:
    - 顺序打印，不刷新，历史可追溯
    - 每个步骤有明确的 [START] / [DONE] 标记
    - 进度更新内联显示，不换行
    - 重要信息用颜色高亮
    - batch 模式下自动简化输出
    """

    def __init__(self) -> None:
        super().__init__()
        self._is_wsl: bool = self._detect_wsl()
        self._disabled: bool = self._env_truthy("RPH_NO_PROGRESS")
        self._title: str = "Processing..."
        self._steps: Dict[str, Dict[str, Any]] = {}
        self._step_order: list[str] = []
        self._current_step: Optional[str] = None
        self._start_time: Optional[datetime] = None
        # 跟踪所有已注册步骤的最终状态（用于 batch 摘要）
        self._step_statuses: Dict[str, str] = {}

    # ---- helpers ----

    def _detect_wsl(self) -> bool:
        if os.environ.get("WSL_DISTRO_NAME"):
            return True
        release = platform.release().lower()
        return "microsoft" in release or "wsl" in release

    def _env_truthy(self, name: str) -> bool:
        value = os.environ.get(name, "").strip().lower()
        return value in {"1", "true", "yes", "on"}

    def _timestamp(self) -> str:
        """统一时间戳格式: [HH:MM:SS]"""
        return datetime.now().strftime("%H:%M:%S")

    def _print(self, message: str, style: Optional[str] = None) -> None:
        """Print with optional rich styling."""
        if self._disabled:
            return
        if style:
            console.print(f"[{style}]{message}[/{style}]")
        else:
            console.print(message)

    def _progress_bar(self, pct: int, width: int = 20) -> str:
        """Generate a simple ASCII progress bar (unicode-safe fallback)."""
        try:
            filled = int(width * pct / 100)
            bar = "█" * filled + "░" * (width - filled)
        except UnicodeEncodeError:
            filled = int(width * pct / 100)
            bar = "#" * filled + "-" * (width - filled)
        return f"[{bar}]"

    def _indent(self) -> str:
        """batch 模式下的缩进前缀。"""
        if self._batch_mode:
            if self._branch_id:
                return f"  │  [{self._branch_id}] "
            if self._reaction_id:
                return f"  │  "
            return "  "
        return ""

    # ---- BaseUIManager interface ----

    def set_context(
        self,
        batch_mode: bool = False,
        reaction_id: str = "",
        branch_id: str = "",
    ) -> None:
        super().set_context(batch_mode, reaction_id, branch_id)

    def start(self, title: str = "Processing...") -> None:
        """Start the pipeline UI."""
        if self._disabled:
            return
        self._title = title
        self._start_time = datetime.now()

        if self._batch_mode:
            # batch 模式：只显示一行简洁标识
            ctx_parts: list[str] = []
            if self._reaction_id:
                ctx_parts.append(self._reaction_id)
            if self._branch_id:
                ctx_parts.append(self._branch_id)
            ctx = "/".join(ctx_parts) if ctx_parts else ""
            self._print(f"\n▶ {'RPH Pipeline' if not ctx else ctx}", "bold blue")
            self._print(f"  {'='*50}", "dim")
        else:
            self._print(f"\n{'='*70}", "bold blue")
            self._print(f"▶ {title}", "bold cyan")
            self._print(f"{'='*70}\n", "bold blue")

    def stop(self) -> None:
        """Finalize the pipeline UI."""
        if self._disabled or not self._start_time:
            return
        elapsed = (datetime.now() - self._start_time).total_seconds()

        if self._batch_mode:
            self._print(f"  {'='*50}", "dim")
            self._print(f"  ⏱ {elapsed:.1f}s", "dim")
        else:
            self._print(f"\n{'='*70}", "bold blue")
            self._print(f"⏱  总耗时: {elapsed:.1f}s", "dim")
            self._print(f"{'='*70}\n", "bold blue")

    def add_step(self, step_id: str, description: str, total: int = 100) -> None:
        """Register a new step."""
        if self._disabled:
            return
        self._steps[step_id] = {
            "description": description,
            "total": total,
            "completed": 0,
            "phase": None,
            "subtask": None,
        }
        if step_id not in self._step_order:
            self._step_order.append(step_id)
        self._current_step = step_id

        indent = self._indent()
        self._print(
            f"{indent}[{self._timestamp()}] [START] {step_id.upper()}: {description}",
            "bold yellow",
        )

    def update_step(
        self,
        step_id: str,
        completed: Optional[int] = None,
        description: Optional[str] = None,
        advance: Optional[int] = None,
    ) -> None:
        """Update step progress."""
        if self._disabled or step_id not in self._steps:
            return

        step = self._steps[step_id]
        if completed is not None:
            step["completed"] = completed
        elif advance is not None:
            step["completed"] += advance

        if description:
            step["description"] = description
            # 从描述中提取最终状态标记
            for tag in ("[OK]", "[SKIPPED]", "[REUSED]", "[DEGRADED]", "[FAILED]"):
                if tag in description:
                    self._step_statuses[step_id] = tag.strip("[]")
                    break

        pct = int((step["completed"] / step["total"]) * 100)

        # 只在重要进度点打印 (0%, 25%, 50%, 75%, 100% 或描述变化)
        if pct in [0, 25, 50, 75, 100] or description:
            bar = self._progress_bar(pct)
            indent = self._indent()
            self._print(
                f"{indent}[{self._timestamp()}] {step_id.upper()} {bar} {pct:3d}% | {step['description']}",
                "cyan",
            )

    def enter_phase(self, step_id: str, phase_name: str) -> None:
        """Mark entry into a phase."""
        if self._disabled or step_id not in self._steps:
            return
        self._steps[step_id]["phase"] = phase_name
        indent = self._indent()
        self._print(f"{indent}  → 阶段: {phase_name}", "dim")

    def set_subtask(
        self, step_id: str, name: str, current: int, total: int
    ) -> None:
        """Update subtask progress."""
        if self._disabled or step_id not in self._steps:
            return
        self._steps[step_id]["subtask"] = {
            "name": name,
            "current": current,
            "total": total,
        }
        # 只在特定间隔打印子任务更新
        if total > 0 and (
            current == 1
            or current == total
            or current % max(1, total // 4) == 0
        ):
            indent = self._indent()
            self._print(f"{indent}  • {name}: {current}/{total}", "dim")

    def log_event(self, step_id: str, message: str) -> None:
        """Log an event."""
        if self._disabled:
            return
        # 过滤掉重复或低价值日志
        skip_patterns = ["Step registered", "Phase:"]
        if any(p in message for p in skip_patterns):
            return
        indent = self._indent()
        self._print(f"{indent}    ℹ {message}", "dim")

    def section_header(
        self, step_id: str, title: str, description: str = ""
    ) -> None:
        """统一的步骤标题渲染。batch 模式压缩，single 模式保留 Rule 样式。"""
        if self._disabled:
            return

        if self._batch_mode:
            # batch: 简洁一行
            header = f"── {title}"
            if description:
                header += f" ({description})"
            indent = self._indent()
            self._print(f"{indent}{header}", "bold magenta")
        else:
            # single: 完整 Rich Rule 分隔线
            console.print()
            console.print(
                Rule(
                    f"[step.header]{step_id}[/] : [step.title]{title}[/]",
                    style="magenta",
                )
            )
            if description:
                console.print(f"[dim italic center]{description}[/]")
            console.print()

    def pipeline_banner(
        self, version: str, reaction_id: str = "", branch_id: str = ""
    ) -> None:
        """统一的 pipeline 启动横幅。batch 模式只显示标识行。"""
        if self._disabled:
            return

        if self._batch_mode or (reaction_id and branch_id):
            # batch / 子分支: 简洁一行
            ctx = f"{reaction_id}/{branch_id}" if branch_id else (reaction_id or "")
            label = f"▶ RPH Pipeline: {ctx}" if ctx else f"▶ RPH Pipeline v{version}"
            self._print(label, "bold cyan")
            return

        # single: 完整 Panel
        grid = Table.grid(expand=True)
        grid.add_column(justify="center", ratio=1)
        grid.add_row(f"[bold cyan]Reaction Profile Hunter[/] [dim]v{version}[/]")
        grid.add_row(
            "[dim]Automated Reaction Mechanism Discovery & Analysis Pipeline[/]"
        )
        panel = Panel(
            grid,
            style="cyan",
            border_style="dim cyan",
            padding=(1, 2),
        )
        console.print(panel)
        console.print()

    def error_block(
        self,
        step_id: str,
        error: Exception,
        context: Optional[Dict[str, Any]] = None,
    ) -> None:
        """结构化错误输出。batch 模式下压缩为 2 行。"""
        if self._disabled:
            return

        indent = self._indent()
        error_msg = ErrorFormatter.compact(error)

        if self._batch_mode:
            # batch: 2 行摘要
            self._print(
                f"{indent}[{self._timestamp()}] ✗ {step_id.upper()} FAILED",
                "bold red",
            )
            for line in error_msg.split("\n")[:3]:
                self._print(f"{indent}    {line}", "red")
        else:
            # single: 带边框的错误块
            self._print("")
            if context:
                ctx_str = ", ".join(f"{k}={v}" for k, v in context.items())
                self._print(f"  Context: {ctx_str}", "dim")
            self._print(f"  [{step_id.upper()}] {error_msg}", "red")
            self._print("")

    def result_summary(self, result: Any) -> None:
        """统一结果摘要渲染。"""
        if self._disabled:
            return

        if self._batch_mode:
            # batch: 一行简洁状态
            status = "✓ OK" if getattr(result, "success", False) else "✗ FAILED"
            style = "green" if getattr(result, "success", False) else "red"
            indent = self._indent()
            csv_path = getattr(result, "features_csv", None)
            detail = f" → {csv_path}" if csv_path else ""
            self._print(f"{indent}[{status}]{detail}", f"bold {style}")
            return

        # single: 完整表格 (委托给 print_result_summary)
        print_result_summary(result)

    def render_tasks(self, tracker: "TaskProgressTracker") -> None:
        """渲染任务树 — V4 任务导向显示。"""
        if self._disabled:
            return
        from rph_core.utils.task_progress import TaskState

        order = tracker.get_task_order()
        indent = self._indent()

        for task_id in order:
            task = tracker.tasks.get(task_id)
            if not task:
                continue
            state = task.state
            icon = state.icon
            color = state.color
            name = task.spec.name
            pct = task.progress_pct
            detail = task.detail
            summary = task.result_summary
            elapsed = task.elapsed_sec

            if state == TaskState.RUNNING and pct > 0:
                bar = self._progress_bar(pct)
                line = f"{indent}│  [{icon}] {name:<26} {bar} {pct:3d}%"
            elif state == TaskState.PENDING:
                deps = task.spec.depends_on
                dep_str = f"(waiting for {', '.join(deps)})" if deps else "(queued)"
                line = f"{indent}│  [{icon}] {name:<26} {dep_str}"
            elif state in (TaskState.COMPLETED, TaskState.CACHED):
                time_str = f"{elapsed:.1f}s" if elapsed > 0 else ""
                sum_str = f"  {summary}" if summary else ""
                line = f"{indent}│  [{color}]{icon}[/{color}] {name:<26}{time_str}{sum_str}"
            elif state == TaskState.FAILED:
                err = task.error_message[:60] if task.error_message else "failed"
                line = f"{indent}│  [red]{icon} {name:<26} {err}[/red]"
            elif state == TaskState.SKIPPED:
                line = f"{indent}│  [dim]{icon} {name:<26} [skipped][/dim]"
            else:
                line = f"{indent}│  [{color}]{icon}[/{color}] {name:<26}"

            self._print(line)

            if detail and state == TaskState.RUNNING:
                self._print(f"{indent}│        └─ {detail}", "dim")

    def get_step_status(self, step_id: str) -> str:
        """获取某步骤的最终状态标记。"""
        return self._step_statuses.get(step_id, "")

    def all_steps_cached(self) -> bool:
        """检查所有注册步骤是否都已被缓存/跳过/复用。"""
        for step_id in self._step_order:
            status = self._step_statuses.get(step_id, "")
            if status not in ("SKIPPED", "REUSED"):
                return False
        return len(self._step_order) > 0


# ============================================================================
# SilentProgressManager — 完全静默
# ============================================================================

class SilentProgressManager(BaseUIManager):
    """完全静默模式 - 只保留文件日志，不输出到终端."""

    def start(self, title: str = "Processing...") -> None:
        pass

    def stop(self) -> None:
        pass

    def add_step(self, step_id: str, description: str, total: int = 100) -> None:
        pass

    def update_step(
        self,
        step_id: str,
        completed: Optional[int] = None,
        description: Optional[str] = None,
        advance: Optional[int] = None,
    ) -> None:
        pass

    def enter_phase(self, step_id: str, phase_name: str) -> None:
        pass

    def set_subtask(
        self, step_id: str, name: str, current: int, total: int
    ) -> None:
        pass

    def log_event(self, step_id: str, message: str) -> None:
        pass

    def section_header(
        self, step_id: str, title: str, description: str = ""
    ) -> None:
        pass

    def pipeline_banner(
        self, version: str, reaction_id: str = "", branch_id: str = ""
    ) -> None:
        pass

    def error_block(
        self,
        step_id: str,
        error: Exception,
        context: Optional[Dict[str, Any]] = None,
    ) -> None:
        pass

    def result_summary(self, result: Any) -> None:
        pass

    def render_tasks(self, tracker: "TaskProgressTracker") -> None:
        pass


# ============================================================================
# RichProgressManager — 富文本动态UI
# ============================================================================

class RichProgressManager(BaseUIManager):
    """
    富文本动态UI - 适合需要实时监控的场景
    保留原有Rich Live功能，但默认禁用
    """

    def __init__(self) -> None:
        super().__init__()
        self.console: Console = console
        self._is_wsl: bool = self._detect_wsl()
        self._live_disabled: bool = self._env_truthy(
            "RPH_NO_PROGRESS"
        ) or self._env_truthy("RPH_NO_LIVE")
        self._use_alt_screen: bool = self._env_truthy("RPH_ALT_SCREEN")
        self._refresh_interval: float = 0.40 if self._is_wsl else 0.10
        self._last_refresh: float = 0.0
        self._pending_refresh: bool = False

        self.progress: Progress = Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            console=self.console,
            transient=True,
        )
        self.live: Optional[Live] = None
        self._tasks: Dict[str, Any] = {}
        self._step_order: list[str] = []
        self._step_descriptions: Dict[str, str] = {}
        self._step_totals: Dict[str, int] = {}
        self._step_completed: Dict[str, int] = {}
        self._active_phase: Dict[str, str] = {}
        self._active_subtasks: Dict[str, Dict[str, Any]] = {}
        self._events: Deque[str] = deque(maxlen=15)
        self._title: str = "Processing..."

    def _detect_wsl(self) -> bool:
        if os.environ.get("WSL_DISTRO_NAME"):
            return True
        release = platform.release().lower()
        return "microsoft" in release or "wsl" in release

    def _env_truthy(self, name: str) -> bool:
        value = os.environ.get(name, "").strip().lower()
        return value in {"1", "true", "yes", "on"}

    def _render_overview_panel(self) -> Panel:
        table = Table.grid(expand=True)
        table.add_column(ratio=1)
        if not self._step_order:
            table.add_row("[dim]No pipeline steps registered.[/]")
        else:
            for step_id in self._step_order:
                description = self._step_descriptions.get(step_id, step_id)
                completed = self._step_completed.get(step_id, 0)
                total = self._step_totals.get(step_id, 100)
                pct = int((completed / total) * 100) if total > 0 else 0
                table.add_row(f"{description} [dim]({pct:>3d}%)[/]")
        return Panel(table, title=f"Overview | {self._title}", border_style="cyan")

    def _render_detail_panel(self) -> Panel:
        table = Table(
            show_header=True, header_style="bold cyan", box=None, expand=True
        )
        table.add_column("Step", style="bold")
        table.add_column("Detail")
        has_detail = False
        for step_id in self._step_order:
            details: list[str] = []
            phase_name = self._active_phase.get(step_id)
            if phase_name:
                details.append(f"Phase: {phase_name}")
            subtask = self._active_subtasks.get(step_id)
            if subtask:
                details.append(
                    f"{subtask['name']}: {subtask['current']:03d} / {subtask['total']:03d}"
                )
            if details:
                has_detail = True
                table.add_row(step_id.upper(), " | ".join(details))
        if not has_detail:
            table.add_row("-", "[dim]No active sub-tasks.[/]")
        return Panel(table, title="Details", border_style="magenta")

    def _render_event_panel(self) -> Panel:
        event_table = Table.grid(expand=True)
        event_table.add_column(ratio=1)
        if self._events:
            for event in self._events:
                event_table.add_row(event)
        else:
            event_table.add_row("[dim]No events yet.[/]")
        return Panel(event_table, title="Event Log", border_style="yellow")

    def _render_layout(self) -> Layout:
        layout = Layout(name="root")
        layout.split_column(
            Layout(name="overview", ratio=2),
            Layout(name="details", ratio=2),
            Layout(name="events", ratio=3),
        )
        layout["overview"].update(self._render_overview_panel())
        layout["details"].update(self._render_detail_panel())
        layout["events"].update(self._render_event_panel())
        return layout

    def _refresh(self, force: bool = False) -> None:
        if self.live:
            now = time.monotonic()
            if force or (now - self._last_refresh) >= self._refresh_interval:
                self.live.update(self._render_layout(), refresh=True)
                self._last_refresh = now
                self._pending_refresh = False
            else:
                self._pending_refresh = True

    def start(self, title: str = "Processing...") -> None:
        """Start the live display."""
        if self._live_disabled:
            return
        if self.live:
            self.stop()
        self._title = title
        self._last_refresh = 0.0
        self._pending_refresh = False
        use_alt_screen = self._use_alt_screen and bool(
            getattr(self.console, "is_terminal", False)
        )
        self.live = Live(
            self._render_layout(),
            console=self.console,
            refresh_per_second=2 if self._is_wsl else 8,
            auto_refresh=False,
            screen=use_alt_screen,
            transient=True,
        )
        self.live.start()
        self._refresh(force=True)

    def stop(self) -> None:
        """Stop the live display."""
        if self.live:
            if self._pending_refresh:
                self._refresh(force=True)
            self.live.stop()
            self.live = None

    def add_step(self, step_id: str, description: str, total: int = 100) -> None:
        """Add a progress task for a specific step."""
        task_id = self.progress.add_task(description, total=total)
        self._tasks[step_id] = task_id
        if step_id not in self._step_order:
            self._step_order.append(step_id)
        self._step_descriptions[step_id] = description
        self._step_totals[step_id] = total
        self._step_completed[step_id] = 0
        self.log_event(step_id, f"Step registered: {description}")
        self._refresh()

    def update_step(
        self,
        step_id: str,
        completed: Optional[int] = None,
        description: Optional[str] = None,
        advance: Optional[int] = None,
    ) -> None:
        """Update a progress task."""
        if step_id not in self._tasks:
            return
        if description is not None:
            self._step_descriptions[step_id] = description
        prior_completed = self._step_completed.get(step_id, 0)
        if completed is not None:
            self._step_completed[step_id] = completed
        elif advance is not None:
            self._step_completed[step_id] = prior_completed + advance
        self.progress.update(
            cast(TaskID, self._tasks[step_id]),
            completed=completed,
            description=description,
            advance=advance,
        )
        self._refresh()

    def enter_phase(self, step_id: str, phase_name: str) -> None:
        if step_id not in self._tasks:
            return
        self._active_phase[step_id] = phase_name
        self.log_event(step_id, f"Phase: {phase_name}")
        self._refresh()

    def set_subtask(
        self, step_id: str, name: str, current: int, total: int
    ) -> None:
        if step_id not in self._tasks:
            return
        safe_total = max(total, 1)
        safe_current = min(max(current, 0), safe_total)
        self._active_subtasks[step_id] = {
            "name": name,
            "current": safe_current,
            "total": safe_total,
        }
        self._refresh()

    def log_event(self, step_id: str, message: str) -> None:
        timestamp = datetime.now().strftime("%H:%M:%S")
        self._events.append(
            f"[dim]{timestamp}[/] [{step_id.upper()}] {message}"
        )
        self._refresh()

    # RichProgressManager 不重写 section_header/pipeline_banner/error_block —
    # Rich Live 模式直接使用面板渲染，不需要额外的控制台输出。

    def section_header(
        self, step_id: str, title: str, description: str = ""
    ) -> None:
        """Rich 模式：section header 作为 log event 显示。"""
        self.log_event(step_id, f"--- {title} ---")

    def pipeline_banner(
        self, version: str, reaction_id: str = "", branch_id: str = ""
    ) -> None:
        """Rich 模式：无额外控制台输出，Live 面板已包含标题。"""
        pass

    def error_block(
        self,
        step_id: str,
        error: Exception,
        context: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Rich 模式：错误显示为 log event。"""
        self.log_event(step_id, f"ERROR: {ErrorFormatter.compact(error)}")

    def result_summary(self, result: Any) -> None:
        """Rich 模式：结果摘要作为 log event。"""
        status = "OK" if getattr(result, "success", False) else "FAILED"
        self.log_event("DONE", f"Pipeline {status}")

    def render_tasks(self, tracker: "TaskProgressTracker") -> None:
        """Rich 模式：任务列表作为 log event。"""
        from rph_core.utils.task_progress import TaskState
        for task_id in tracker.get_task_order():
            task = tracker.tasks.get(task_id)
            if not task:
                continue
            icon = task.state.icon
            name = task.spec.name
            if task.state == TaskState.RUNNING:
                detail = f" ({task.progress_pct}%)" if task.progress_pct > 0 else ""
                self.log_event(task_id, f"{icon} {name}{detail}")
            elif task.state == TaskState.PENDING:
                deps = task.spec.depends_on
                dep_str = f" waiting for {', '.join(deps)}" if deps else ""
                self.log_event(task_id, f"{icon} {name}{dep_str}")
            else:
                self.log_event(task_id, f"{icon} {name}")


# ============================================================================
# 辅助函数 - 保持向后兼容
# ============================================================================

def print_pipeline_header(version: str) -> None:
    """打印Pipeline头部信息 (向后兼容 — 委托给 UIManager)."""
    pm = get_progress_manager()
    pm.pipeline_banner(version)


def print_step_header(
    step_id: str, title: str, description: str = ""
) -> None:
    """打印步骤头部 (向后兼容 — 委托给 UIManager)."""
    pm = get_progress_manager()
    pm.section_header(step_id, title, description)


@contextmanager
def status(status_text: str, spinner: str = "dots"):
    """上下文管理器用于显示临时状态."""
    with console.status(
        f"[bold cyan]{status_text}[/]", spinner=spinner
    ) as status_ctx:
        yield status_ctx


def print_result_summary(result: Any) -> None:
    """打印执行结果摘要 (向后兼容)."""
    console.print()

    if result.success:
        status_text = "[bold green]✓ PIPELINE SUCCESS[/]"
        border_style = "green"
    else:
        status_text = (
            f"[bold red]✗ PIPELINE FAILED[/] (Step: {result.error_step})"
        )
        border_style = "red"

    table = Table(
        title="Execution Summary",
        box=None,
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("Metric", style="dim")
    table.add_column("Value", style="bold")

    table.add_row("SMILES", str(result.product_smiles))
    table.add_row("Work Dir", str(result.work_dir))

    if result.e_product_l2:
        table.add_row("Product E (L2)", f"{result.e_product_l2:.6f} Ha")

    if result.sp_matrix_report:
        try:
            dg_act = result.sp_matrix_report.get_activation_energy()
            if dg_act is not None:
                table.add_row(
                    "ΔG‡ (Act)", f"[yellow]{dg_act:.2f} kcal/mol[/]"
                )
            dg_rxn = result.sp_matrix_report.get_reaction_energy()
            if dg_rxn is not None:
                table.add_row("ΔG (Rxn)", f"{dg_rxn:.2f} kcal/mol")
        except Exception:
            pass

    table.add_section()
    table.add_row("Product XYZ", _format_path(result.product_xyz))
    table.add_row("TS Final XYZ", _format_path(result.ts_final_xyz))
    table.add_row("Features", _format_path(result.features_csv))

    if not result.success and result.error_message:
        table.add_section()
        table.add_row("Error", f"[red]{result.error_message}[/]")

    panel = Panel(
        table,
        title=status_text,
        border_style=border_style,
        padding=(1, 2),
    )
    console.print(panel)


def _format_path(path: Optional[Any]) -> str:
    """格式化路径显示."""
    if not path:
        return "[dim]-[/]"
    p = str(path)
    if len(p) > 60:
        return f"...{p[-57:]}"
    return p
