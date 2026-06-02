"""批量运行结果总览 — Rich Table 渲染。

在所有 reaction/branch 处理完成后，输出一个汇总表格，
显示每个 pipeline run 的关键指标和状态。
"""

from datetime import datetime
from pathlib import Path
from typing import Any, List, Optional

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from rph_core.utils.shared_console import get_console

console: Console = get_console()


def render_dashboard(
    results: List[Any],
    total_time: Optional[float] = None,
    title: str = "RPH Batch Run Summary",
) -> None:
    """渲染批量运行结果总览表格。

    Args:
        results: PipelineResult 对象列表
        total_time: 总耗时（秒），可选
        title: 表格标题
    """
    if not results:
        console.print("[dim]No results to display.[/]")
        return

    table = Table(
        title=title,
        box=None,
        show_header=True,
        header_style="bold cyan",
        expand=True,
    )
    table.add_column("Reaction", style="cyan", no_wrap=True)
    table.add_column("Branch", style="dim", no_wrap=True)
    table.add_column("Type", style="dim")
    table.add_column("ΔG‡ (kcal/mol)", justify="right")
    table.add_column("Status", justify="center")
    table.add_column("Output", style="dim")

    success_count = 0
    failed_count = 0

    for r in results:
        reaction_id = _extract_reaction_id(r)
        branch_id = _extract_branch_id(r)
        reaction_type = getattr(r, "reaction_type", None) or "-"
        dg_act = _extract_activation_energy(r)
        success = getattr(r, "success", False)

        if success:
            success_count += 1
            status = "[green]✓[/]"
        else:
            failed_count += 1
            status = f"[red]✗[/]"

        dg_str = f"{dg_act:.2f}" if dg_act is not None else "-"
        features_path = _format_output_path(getattr(r, "features_csv", None))

        table.add_row(
            reaction_id,
            branch_id,
            str(reaction_type),
            dg_str,
            status,
            features_path,
        )

    table.add_section()
    summary_style = "green" if failed_count == 0 else "yellow" if success_count > 0 else "red"
    total = len(results)
    time_str = f"⏱ {total_time:.1f}s" if total_time is not None else ""
    table.add_row(
        f"Total: {total}",
        "",
        "",
        "",
        f"[{summary_style}]{success_count} ✓ / {failed_count} ✗[/]",
        time_str,
    )

    if failed_count > 0:
        table.add_section()
        table.add_row(
            "[bold red]Failed Runs:[/]",
            "", "", "", "", ""
        )
        for r in results:
            if not getattr(r, "success", False):
                rid = _extract_reaction_id(r)
                bid = _extract_branch_id(r)
                err_step = getattr(r, "error_step", "Unknown")
                err_msg = str(getattr(r, "error_message", "") or "")[:80]
                table.add_row(
                    f"  {rid}", bid, "", "",
                    f"[red]{err_step}[/]",
                    f"[dim red]{err_msg}[/]",
                )

    panel = Panel(
        table,
        border_style="cyan",
        padding=(1, 2),
    )
    console.print()
    console.print(panel)
    console.print()


def _extract_reaction_id(result: Any) -> str:
    work_dir = getattr(result, "work_dir", None)
    if work_dir is None:
        return "-"
    path = Path(str(work_dir))
    parts = path.parts
    for i, part in enumerate(parts):
        if part.startswith("RXN_") or part.startswith("rx_"):
            return part
    parent = path.parent
    if parent.name.startswith("RXN_") or parent.name.startswith("rx_"):
        return parent.name
    return path.name[:20]


def _extract_branch_id(result: Any) -> str:
    work_dir = getattr(result, "work_dir", None)
    if work_dir is None:
        return "-"
    path = Path(str(work_dir))
    if "branches" in path.parts:
        idx = path.parts.index("branches")
        if idx + 1 < len(path.parts):
            return path.parts[idx + 1]
    return path.name if path.name.startswith("BR_") else "-"


def _extract_activation_energy(result: Any) -> Optional[float]:
    sp_report = getattr(result, "sp_matrix_report", None)
    if sp_report is None:
        return None
    try:
        return sp_report.get_activation_energy()
    except Exception:
        return None


def _format_output_path(path: Any) -> str:
    if not path:
        return "[dim]-[/]"
    s = str(path)
    if len(s) > 50:
        return f"[dim]...{s[-47:]}[/]"
    return f"[dim]{s}[/]"


def render_timing_summary(
    start_time: datetime,
    end_time: Optional[datetime] = None,
    reaction_count: int = 0,
    branch_count: int = 0,
) -> None:
    """渲染计时摘要。

    Args:
        start_time: 开始时间
        end_time: 结束时间（默认 now）
        reaction_count: 反应数量
        branch_count: 分支数量
    """
    end = end_time or datetime.now()
    elapsed = (end - start_time).total_seconds()

    table = Table(box=None, show_header=False, expand=False)
    table.add_column(style="dim")
    table.add_column(style="bold")

    table.add_row("Reactions", str(reaction_count))
    table.add_row("Branches", str(branch_count))
    table.add_row("Total time", f"{elapsed:.1f}s")
    if reaction_count > 0:
        table.add_row("Avg / reaction", f"{elapsed / reaction_count:.1f}s")

    panel = Panel(table, title="[bold]Timing Summary[/]", border_style="dim cyan")
    console.print()
    console.print(panel)
    console.print()
