"""Shared console helpers for Rich and plain-text UI output."""

from __future__ import annotations

import sys
from typing import Any, TextIO

try:
    from rich.console import Console as _RichConsole
    from rich.theme import Theme as _RichTheme

    _has_rich = True
except ImportError:  # pragma: no cover - exercised in rich-less environments
    _RichConsole = None  # type: ignore[assignment]
    _RichTheme = None  # type: ignore[assignment]
    _has_rich = False

HAS_RICH = _has_rich


class _PlainConsole:
    """Minimal console fallback when Rich is unavailable."""

    def __init__(self) -> None:
        self.file: TextIO = sys.stdout

    def print(self, *objects: Any, sep: str = " ", end: str = "\n", **kwargs: Any) -> None:
        text = sep.join("" if item is None else str(item) for item in objects)
        _written = self.file.write(text + end)
        _flushed = self.file.flush()


_THEME_STYLES = {
    "info": "dim cyan",
    "warning": "yellow",
    "error": "bold red",
    "success": "bold green",
    "step.header": "bold magenta",
    "step.title": "bold white",
    "status.pending": "dim",
    "status.running": "cyan",
    "status.complete": "green",
    "status.failed": "red",
    "status.cached": "blue",
    "status.degraded": "yellow",
    "task.opt": "bright_blue",
    "task.freq": "bright_magenta",
    "task.sp": "bright_cyan",
    "energy": "green",
}

RPH_THEME = _RichTheme(_THEME_STYLES) if HAS_RICH and _RichTheme is not None else None
_shared_console = _RichConsole(theme=RPH_THEME) if HAS_RICH and _RichConsole is not None else _PlainConsole()


def get_console() -> Any:
    """Get the global shared console, Rich when available."""

    return _shared_console
