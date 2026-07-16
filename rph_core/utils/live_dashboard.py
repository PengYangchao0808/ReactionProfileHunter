# pyright: reportMissingImports=false
"""Lifecycle management for the non-full-screen embedded Rich dashboard."""

from __future__ import annotations

import threading
from typing import Any

from rph_core.utils.dashboard_renderers import HAS_RICH, render_dashboard
from rph_core.utils.dashboard_state import DashboardStateReducer

try:
    from rich.live import Live
except ImportError:  # pragma: no cover
    Live = None  # type: ignore[assignment,misc]


class LiveDashboardManager:
    """Render reducer snapshots on one UI thread without touching QC threads."""

    def __init__(
        self,
        console: Any,
        *,
        log_path: str = "",
        refresh_per_second: float = 4.0,
    ) -> None:
        self.console = console
        self.reducer = DashboardStateReducer(log_path=log_path)
        self.refresh_per_second = max(1.0, min(10.0, float(refresh_per_second)))
        self._live: Any = None
        self._thread: threading.Thread | None = None
        self._stop = threading.Event()
        self._started = False
        self._closed = False
        self._failed = False
        self._last_revision = -1
        self._failure_message = ""

    @property
    def failed(self) -> bool:
        return self._failed

    @property
    def failure_message(self) -> str:
        return self._failure_message

    def start(self) -> bool:
        if self._started:
            return not self._failed
        self._started = True
        if not HAS_RICH or Live is None:
            self._fail("Rich Live is unavailable")
            return False
        try:
            initial = render_dashboard(self.reducer.snapshot(), width=self._width())
            self._live = Live(
                initial,
                console=self.console,
                screen=False,
                auto_refresh=False,
                transient=False,
                vertical_overflow="crop",
            )
            self._live.start(refresh=True)
            self._thread = threading.Thread(
                target=self._refresh_loop,
                name="rph-live-dashboard",
                daemon=True,
            )
            self._thread.start()
            return True
        except Exception as exc:  # pragma: no cover - depends on terminal backend
            self._fail(f"could not start Rich Live: {exc}")
            return False

    def event_callback(self, event: str, record: dict[str, Any]) -> None:
        if not self._failed:
            self.reducer.apply(event, record)

    def finish(self, status: str = "complete") -> None:
        if self._closed:
            return
        self.reducer.finish(status)
        self._stop.set()
        if self._thread is not None and self._thread is not threading.current_thread():
            self._thread.join(timeout=2.0)
        if self._live is not None and not self._failed:
            try:
                self._live.update(
                    render_dashboard(self.reducer.snapshot(), width=self._width(), final=True),
                    refresh=True,
                )
            except Exception as exc:  # pragma: no cover
                self._fail(f"final dashboard refresh failed: {exc}")

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self._stop.set()
        if self._thread is not None and self._thread is not threading.current_thread():
            self._thread.join(timeout=2.0)
        if self._live is not None:
            try:
                self._live.stop()
            except Exception:  # pragma: no cover - terminal cleanup is best effort
                pass

    def _refresh_loop(self) -> None:
        interval = 1.0 / self.refresh_per_second
        while not self._stop.wait(interval):
            revision = self.reducer.revision
            if revision == self._last_revision:
                continue
            try:
                snapshot = self.reducer.snapshot()
                self._live.update(
                    render_dashboard(snapshot, width=self._width()),
                    refresh=True,
                )
                self._last_revision = revision
            except Exception as exc:  # pragma: no cover - terminal backend dependent
                self._fail(f"dashboard refresh failed: {exc}")
                return

    def _width(self) -> int:
        return int(getattr(self.console, "width", 120) or 120)

    def _fail(self, message: str) -> None:
        self._failed = True
        self._failure_message = str(message)
        self._stop.set()
        if self._live is not None:
            try:
                self._live.stop()
            except Exception:
                pass
