"""
Shared Console instance for Rich integration.
Ensures logging and UI components use the same console for Live/Progress compatibility.
"""

from rich.console import Console
from rich.theme import Theme

# Define common theme to be shared across UI and Logging
RPH_THEME = Theme({
    "info": "dim cyan",
    "warning": "yellow",
    "error": "bold red",
    "success": "bold green",
    "step.header": "bold magenta",
    "step.title": "bold white",
})

# Singleton console instance
_shared_console = Console(theme=RPH_THEME)

def get_console() -> Console:
    """Get the global shared Rich console."""
    return _shared_console
