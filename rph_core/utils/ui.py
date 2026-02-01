from rich.console import Console
from rich.panel import Panel
from rich.rule import Rule
from rich.table import Table
from rich.theme import Theme
from typing import Optional, Any
from contextlib import contextmanager

custom_theme = Theme({
    "info": "dim cyan",
    "warning": "yellow",
    "error": "bold red",
    "success": "bold green",
    "step.header": "bold magenta",
    "step.title": "bold white",
})

console = Console(theme=custom_theme)

def print_pipeline_header(version: str):
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
    console.print()
    console.print(Rule(f"[step.header]{step_id}[/] : [step.title]{title}[/]", style="magenta"))
    if description:
        console.print(f"[dim italic center]{description}[/]")
    console.print()

@contextmanager
def status(status_text: str, spinner: str = "dots"):
    with console.status(f"[bold cyan]{status_text}[/]", spinner=spinner) as status_ctx:
        yield status_ctx

def print_result_summary(result: Any):
    console.print()
    
    if result.success:
        status_text = "[bold green]PIPELINE SUCCESS[/]"
        border_style = "green"
    else:
        status_text = f"[bold red]PIPELINE FAILED[/] (Step: {result.error_step})"
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
        except:
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
    if not path:
        return "[dim]-[/]"
    p = str(path)
    if len(p) > 60:
        return f"...{p[-57:]}"
    return p
