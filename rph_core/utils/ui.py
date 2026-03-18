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

from collections import deque
from datetime import datetime
import os
import platform
import time
from typing import Optional, Any, Dict, Deque, cast
from contextlib import contextmanager

# Rich imports - only needed for rich mode
from rich.layout import Layout
from rich.live import Live
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn, TaskID
from rich.rule import Rule
from rich.table import Table

from rph_core.utils.shared_console import get_console

console = get_console()

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


class BaseUIManager:
    """Abstract base class defining UI manager interface."""
    
    def start(self, title: str = "Processing..."):
        """Initialize the UI."""
        raise NotImplementedError
    
    def stop(self):
        """Clean up and finalize the UI."""
        raise NotImplementedError
    
    def add_step(self, step_id: str, description: str, total: int = 100):
        """Register a new pipeline step."""
        raise NotImplementedError
    
    def update_step(self, step_id: str, completed: Optional[int] = None, 
                    description: Optional[str] = None, advance: Optional[int] = None):
        """Update step progress."""
        raise NotImplementedError
    
    def enter_phase(self, step_id: str, phase_name: str):
        """Mark entry into a phase within a step."""
        raise NotImplementedError
    
    def set_subtask(self, step_id: str, name: str, current: int, total: int):
        """Update subtask progress."""
        raise NotImplementedError
    
    def log_event(self, step_id: str, message: str):
        """Log an event/message."""
        raise NotImplementedError
    
    @contextmanager
    def manage(self, title: str = "Pipeline Progress"):
        """Context manager for easy usage."""
        try:
            self.start(title)
            yield self
        finally:
            self.stop()


class LoggerProgressManager(BaseUIManager):
    """
    简洁日志式UI - 默认模式
    
    特点:
    - 顺序打印，不刷新，历史可追溯
    - 每个步骤有明确的 [START] / [DONE] 标记
    - 进度更新内联显示，不换行
    - 重要信息用颜色高亮
    """
    
    def __init__(self):
        self._is_wsl = self._detect_wsl()
        self._disabled = self._env_truthy("RPH_NO_PROGRESS")
        self._title = "Processing..."
        self._steps: Dict[str, Dict[str, Any]] = {}
        self._step_order: list = []
        self._current_step: Optional[str] = None
        self._start_time: Optional[datetime] = None
        
    def _detect_wsl(self) -> bool:
        if os.environ.get("WSL_DISTRO_NAME"):
            return True
        release = platform.release().lower()
        return "microsoft" in release or "wsl" in release
    
    def _env_truthy(self, name: str) -> bool:
        value = os.environ.get(name, "").strip().lower()
        return value in {"1", "true", "yes", "on"}
    
    def _timestamp(self) -> str:
        return datetime.now().strftime("%H:%M:%S")
    
    def _print(self, message: str, style: Optional[str] = None):
        """Print with optional rich styling."""
        if self._disabled:
            return
        if style:
            console.print(f"[{style}]{message}[/{style}]")
        else:
            console.print(message)
    
    def start(self, title: str = "Processing..."):
        """Start the pipeline UI."""
        if self._disabled:
            return
        self._title = title
        self._start_time = datetime.now()
        self._print(f"\n{'='*70}", "bold blue")
        self._print(f"▶ {title}", "bold cyan")
        self._print(f"{'='*70}\n", "bold blue")
    
    def stop(self):
        """Finalize the pipeline UI."""
        if self._disabled or not self._start_time:
            return
        elapsed = (datetime.now() - self._start_time).total_seconds()
        self._print(f"\n{'='*70}", "bold blue")
        self._print(f"⏱  总耗时: {elapsed:.1f}s", "dim")
        self._print(f"{'='*70}\n", "bold blue")
    
    def add_step(self, step_id: str, description: str, total: int = 100):
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
        self._print(f"[{self._timestamp()}] [START] {step_id.upper()}: {description}", "bold yellow")
    
    def update_step(self, step_id: str, completed: Optional[int] = None, 
                    description: Optional[str] = None, advance: Optional[int] = None):
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
        
        pct = int((step["completed"] / step["total"]) * 100)
        
        # 只在重要进度点打印 (0%, 25%, 50%, 75%, 100% 或描述变化)
        if pct in [0, 25, 50, 75, 100] or description:
            bar = self._progress_bar(pct)
            self._print(f"  [{self._timestamp()}] {step_id.upper()} {bar} {pct:3d}% | {step['description']}", "cyan")
    
    def _progress_bar(self, pct: int, width: int = 20) -> str:
        """Generate a simple ASCII progress bar."""
        filled = int(width * pct / 100)
        bar = "█" * filled + "░" * (width - filled)
        return f"[{bar}]"
    
    def enter_phase(self, step_id: str, phase_name: str):
        """Mark entry into a phase."""
        if self._disabled or step_id not in self._steps:
            return
        self._steps[step_id]["phase"] = phase_name
        self._print(f"    → 阶段: {phase_name}", "dim")
    
    def set_subtask(self, step_id: str, name: str, current: int, total: int):
        """Update subtask progress."""
        if self._disabled or step_id not in self._steps:
            return
        self._steps[step_id]["subtask"] = {"name": name, "current": current, "total": total}
        # 只在特定间隔打印子任务更新
        if total > 0 and (current == 1 or current == total or current % max(1, total // 4) == 0):
            self._print(f"    • {name}: {current}/{total}", "dim")
    
    def log_event(self, step_id: str, message: str):
        """Log an event."""
        if self._disabled:
            return
        # 过滤掉重复或低价值日志
        skip_patterns = ["Step registered", "Phase:"]
        if any(p in message for p in skip_patterns):
            return
        self._print(f"      ℹ {message}", "dim")
    
    def complete_step(self, step_id: str, status: str = "OK"):
        """Mark step as complete with status."""
        if self._disabled or step_id not in self._steps:
            return
        step = self._steps[step_id]
        color = "green" if status == "OK" else "yellow" if status == "SKIPPED" else "red"
        self._print(f"[{self._timestamp()}] [DONE] {step_id.upper()}: {step['description']} [{status}]\n", f"bold {color}")


class SilentProgressManager(BaseUIManager):
    """完全静默模式 - 只保留文件日志，不输出到终端."""
    
    def start(self, title: str = "Processing..."): pass
    def stop(self): pass
    def add_step(self, step_id: str, description: str, total: int = 100): pass
    def update_step(self, step_id: str, completed: Optional[int] = None, 
                    description: Optional[str] = None, advance: Optional[int] = None): pass
    def enter_phase(self, step_id: str, phase_name: str): pass
    def set_subtask(self, step_id: str, name: str, current: int, total: int): pass
    def log_event(self, step_id: str, message: str): pass


class RichProgressManager(BaseUIManager):
    """
    富文本动态UI - 适合需要实时监控的场景
    保留原有Rich Live功能，但默认禁用
    """
    
    def __init__(self):
        self.console = console
        self._is_wsl = self._detect_wsl()
        self._live_disabled = self._env_truthy("RPH_NO_PROGRESS") or self._env_truthy("RPH_NO_LIVE")
        self._use_alt_screen = self._env_truthy("RPH_ALT_SCREEN")
        self._refresh_interval = 0.40 if self._is_wsl else 0.10
        self._last_refresh = 0.0
        self._pending_refresh = False
        
        self.progress = Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            console=self.console,
            transient=True
        )
        self.live: Optional[Live] = None
        self._tasks: Dict[str, Any] = {}
        self._step_order: list = []
        self._step_descriptions: Dict[str, str] = {}
        self._step_totals: Dict[str, int] = {}
        self._step_completed: Dict[str, int] = {}
        self._active_phase: Dict[str, str] = {}
        self._active_subtasks: Dict[str, Dict[str, Any]] = {}
        self._events: Deque[str] = deque(maxlen=15)
        self._title = "Processing..."

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
        table = Table(show_header=True, header_style="bold cyan", box=None, expand=True)
        table.add_column("Step", style="bold")
        table.add_column("Detail")

        has_detail = False
        for step_id in self._step_order:
            details = []
            phase_name = self._active_phase.get(step_id)
            if phase_name:
                details.append(f"Phase: {phase_name}")

            subtask = self._active_subtasks.get(step_id)
            if subtask:
                details.append(f"{subtask['name']}: {subtask['current']:03d} / {subtask['total']:03d}")

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

    def _refresh(self, force: bool = False):
        if self.live:
            now = time.monotonic()
            if force or (now - self._last_refresh) >= self._refresh_interval:
                self.live.update(self._render_layout(), refresh=True)
                self._last_refresh = now
                self._pending_refresh = False
            else:
                self._pending_refresh = True

    def start(self, title: str = "Processing..."):
        """Start the live display."""
        if self._live_disabled:
            return
        if self.live:
            self.stop()
        self._title = title
        self._last_refresh = 0.0
        self._pending_refresh = False
        use_alt_screen = self._use_alt_screen and bool(getattr(self.console, "is_terminal", False))
        self.live = Live(
            self._render_layout(),
            console=self.console,
            refresh_per_second=2 if self._is_wsl else 8,
            auto_refresh=False,
            screen=use_alt_screen,
            transient=True
        )
        self.live.start()
        self._refresh(force=True)

    def stop(self):
        """Stop the live display."""
        if self.live:
            if self._pending_refresh:
                self._refresh(force=True)
            self.live.stop()
            self.live = None

    def add_step(self, step_id: str, description: str, total: int = 100):
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

    def update_step(self, step_id: str, completed: Optional[int] = None, 
                    description: Optional[str] = None, advance: Optional[int] = None):
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
            advance=advance
        )
        self._refresh()

    def enter_phase(self, step_id: str, phase_name: str):
        if step_id not in self._tasks:
            return
        self._active_phase[step_id] = phase_name
        self.log_event(step_id, f"Phase: {phase_name}")
        self._refresh()

    def set_subtask(self, step_id: str, name: str, current: int, total: int):
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

    def log_event(self, step_id: str, message: str):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self._events.append(f"[dim]{timestamp}[/] [{step_id.upper()}] {message}")
        self._refresh()


# ============================================================================
# 辅助函数 - 保持向后兼容
# ============================================================================

def print_pipeline_header(version: str):
    """打印Pipeline头部信息."""
    grid = Table.grid(expand=True)
    grid.add_column(justify="center", ratio=1)
    grid.add_row(f"[bold cyan]Reaction Profile Hunter[/] [dim]v{version}[/]")
    grid.add_row("[dim]Automated Reaction Mechanism Discovery & Analysis Pipeline[/]")
    
    panel = Panel(
        grid,
        style="cyan",
        border_style="dim cyan",
        padding=(1, 2),
    )
    console.print(panel)
    console.print()


def print_step_header(step_id: str, title: str, description: str = ""):
    """打印步骤头部."""
    console.print()
    console.print(Rule(f"[step.header]{step_id}[/] : [step.title]{title}[/]", style="magenta"))
    if description:
        console.print(f"[dim italic center]{description}[/]")
    console.print()


@contextmanager
def status(status_text: str, spinner: str = "dots"):
    """上下文管理器用于显示临时状态."""
    with console.status(f"[bold cyan]{status_text}[/]", spinner=spinner) as status_ctx:
        yield status_ctx


def print_result_summary(result: Any):
    """打印执行结果摘要."""
    console.print()
    
    if result.success:
        status_text = "[bold green]✓ PIPELINE SUCCESS[/]"
        border_style = "green"
    else:
        status_text = f"[bold red]✗ PIPELINE FAILED[/] (Step: {result.error_step})"
        border_style = "red"
        
    table = Table(title="Execution Summary", box=None, show_header=True, header_style="bold cyan")
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
                table.add_row("ΔG‡ (Act)", f"[yellow]{dg_act:.2f} kcal/mol[/]")
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
        padding=(1, 2)
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