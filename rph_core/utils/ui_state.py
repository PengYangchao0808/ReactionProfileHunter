"""Unified UI status normalization for V4 terminal reporting."""

from __future__ import annotations

from enum import Enum


class UiStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETE = "complete"
    FAILED = "failed"
    DEGRADED = "degraded"
    CACHED = "cached"


_STATUS_MAP: dict[str, UiStatus] = {
    "pending": UiStatus.PENDING,
    "running": UiStatus.RUNNING,
    "complete": UiStatus.COMPLETE,
    "completed": UiStatus.COMPLETE,
    "completed_with_failures": UiStatus.COMPLETE,
    "failed": UiStatus.FAILED,
    "degraded": UiStatus.DEGRADED,
    "opt_failed_sp_complete": UiStatus.DEGRADED,
    "ts_frequency_unverified": UiStatus.DEGRADED,
    "cached": UiStatus.CACHED,
    "reused": UiStatus.CACHED,
    "skipped": UiStatus.CACHED,
}

_STATUS_STYLE: dict[UiStatus, tuple[str, str]] = {
    UiStatus.PENDING: ("○", "dim"),
    UiStatus.RUNNING: ("→", "cyan"),
    UiStatus.COMPLETE: ("✓", "green"),
    UiStatus.FAILED: ("✗", "red"),
    UiStatus.DEGRADED: ("⚠", "yellow"),
    UiStatus.CACHED: ("↻", "blue"),
}


def normalize_status(raw: str | None) -> UiStatus:
    return _STATUS_MAP.get(str(raw or "").strip().lower(), UiStatus.PENDING)


def status_markup(status: UiStatus) -> str:
    icon, color = _STATUS_STYLE[status]
    return f"[{color}]{icon}[/]"


def status_style(status: UiStatus) -> tuple[str, str]:
    return _STATUS_STYLE[status]
