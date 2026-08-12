# pyright: reportMissingImports=false
"""Lifecycle management for the non-full-screen embedded Rich dashboard."""

from __future__ import annotations

import logging
import threading
import time
from typing import Any

from rph_core.utils.dashboard_renderers import HAS_RICH, render_dashboard
from rph_core.utils.dashboard_state import DashboardStateReducer
from rph_core.utils.run_id import RUN_ID_FIELD

try:
    from rich.live import Live
except ImportError:  # pragma: no cover
    Live = None  # type: ignore[assignment,misc]


logger = logging.getLogger(__name__)


class LiveDashboardManager:
    """Render reducer snapshots on one UI thread without touching QC threads."""

    def __init__(
        self,
        console: Any,
        *,
        log_path: str = "",
        refresh_per_second: float = 4.0,
        run_id: str | None = None,
    ) -> None:
        self.console = console
        self.reducer = DashboardStateReducer(log_path=log_path)
        self.run_id = str(run_id).strip() if run_id is not None else None
        if self.run_id == "":
            self.run_id = None
        self.refresh_per_second = max(1.0, min(10.0, float(refresh_per_second)))
        self._live: Any = None
        self._thread: threading.Thread | None = None
        self._stop = threading.Event()
        self._started = False
        self._closed = False
        self._failed = False
        self._last_revision = -1
        self._last_elapsed_second = -1
        self._failure_message = ""

    @property
    def failed(self) -> bool:
        return self._failed

    @property
    def failure_message(self) -> str:
        return self._failure_message

    def set_run_id(self, run_id: str | None) -> None:
        self.run_id = str(run_id).strip() if run_id is not None else None
        if self.run_id == "":
            self.run_id = None

    def start(self) -> bool:
        if self._started:
            return not self._failed
        self._started = True
        if not HAS_RICH or Live is None:
            self._fail("Rich Live is unavailable")
            return False
        try:
            initial_snapshot = self.reducer.snapshot()
            initial = render_dashboard(initial_snapshot, width=self._width())
            self._live = Live(
                initial,
                console=self.console,
                screen=False,
                auto_refresh=False,
                transient=False,
                vertical_overflow="crop",
            )
            self._live.start(refresh=True)
            self._last_revision = initial_snapshot.revision
            self._last_elapsed_second = self._elapsed_second(initial_snapshot)
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
        if not self._failed and not self._should_ignore_run_id(record):
            self.reducer.apply(event, record)

    def _should_ignore_run_id(self, record: dict[str, Any]) -> bool:
        try:
            if self.run_id is None:
                return False
            payload_run_id = record.get(RUN_ID_FIELD)
            if payload_run_id in (None, ""):
                return False
            observed_run_id = str(payload_run_id)
            if observed_run_id != self.run_id:
                logger.debug(
                    "Ignoring event from stale run_id %s (expected %s)",
                    observed_run_id,
                    self.run_id,
                )
                return True
        except Exception:  # pragma: no cover - defensive legacy compatibility
            logger.debug("Run-id event filtering failed; accepting dashboard event", exc_info=True)
        return False

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
            try:
                snapshot = self.reducer.snapshot()
                elapsed_second = self._elapsed_second(snapshot)
                # Rich Live redraws by erasing and repainting the complete
                # dashboard.  Avoid doing that on every polling tick: a state
                # change needs a prompt render, while the only idle value that
                # changes is the elapsed-time display, at one-second cadence.
                if (
                    revision == self._last_revision
                    and elapsed_second == self._last_elapsed_second
                ):
                    continue
                self._live.update(
                    render_dashboard(snapshot, width=self._width()),
                    refresh=True,
                )
                self._last_revision = revision
                self._last_elapsed_second = elapsed_second
            except Exception as exc:  # pragma: no cover - terminal backend dependent
                self._fail(f"dashboard refresh failed: {exc}")
                return

    def _width(self) -> int:
        return int(getattr(self.console, "width", 120) or 120)

    @staticmethod
    def _elapsed_second(snapshot: Any) -> int:
        return max(0, int(time.time() - float(snapshot.started_at)))

    def _fail(self, message: str) -> None:
        self._failed = True
        self._failure_message = str(message)
        self._stop.set()
        if self._live is not None:
            try:
                self._live.stop()
            except Exception:
                pass
