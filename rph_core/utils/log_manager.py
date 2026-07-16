"""
Logging Manager
================

统一的日志管理系统 (Rich Enhanced)
"""

from __future__ import annotations

import importlib
import logging
import sys
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from rich.console import Console
    from rich.logging import RichHandler as RichHandlerType

try:
    from rich.logging import RichHandler as ImportedRichHandler

    _rich_handler_cls: type[RichHandlerType] | None = ImportedRichHandler
    _has_rich = True
except ImportError:
    _has_rich = False
    _rich_handler_cls = None

HAS_RICH = _has_rich
_V4_LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
_V4_DETAIL_LOGGERS = (
    "rph_core.utils.orca_interface",
    "rph_core.utils.resource_utils",
    "rph_core.steps.conformer_search.censo_lite",
    "rph_core.steps.conformer_search.censo_lite_runtime",
    "rph_core.steps.conformer_search.xtb_thermo",
    "rph_core.utils.s4_progress",
)


def _build_v4_formatter() -> logging.Formatter:
    return logging.Formatter(_V4_LOG_FORMAT)


def _get_console() -> Console:
    shared_console = importlib.import_module("rph_core.utils.shared_console")
    console_factory = cast("Callable[[], Console]", getattr(shared_console, "get_console"))
    return console_factory()


def setup_v4_logging(
    log_file: Path | None,
    level: int = logging.INFO,
    rich_console: bool = True,
) -> None:
    """Configure root logging for the V4 pipeline without disturbing host handlers."""

    root = logging.getLogger()
    root.setLevel(level)

    for handler in list(root.handlers):
        if getattr(handler, "_rph_v4_console", False):
            root.removeHandler(handler)
            handler.close()

    if rich_console and HAS_RICH:
        if _rich_handler_cls is None:
            raise RuntimeError("Rich logging unavailable while HAS_RICH is True")
        console_handler: logging.Handler = _rich_handler_cls(
            console=_get_console(),
            rich_tracebacks=True,
            show_time=True,
            show_path=False,
            markup=True,
        )
    else:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(_build_v4_formatter())

    console_handler.setLevel(level)
    setattr(console_handler, "_rph_v4_console", True)
    root.addHandler(console_handler)

    resolved_log_file = Path(log_file).resolve() if log_file is not None else None
    active_file_handler: logging.Handler | None = None

    for handler in list(root.handlers):
        owned_log_path = getattr(handler, "_rph_v4_log_path", None)
        if owned_log_path is None:
            continue
        if resolved_log_file is None or owned_log_path != str(resolved_log_file):
            root.removeHandler(handler)
            handler.close()
            continue
        if active_file_handler is None:
            handler.setLevel(logging.DEBUG)
            handler.setFormatter(_build_v4_formatter())
            active_file_handler = handler
            continue
        root.removeHandler(handler)
        handler.close()

    if resolved_log_file is None or active_file_handler is not None:
        if resolved_log_file is not None:
            for logger_name in _V4_DETAIL_LOGGERS:
                logging.getLogger(logger_name).setLevel(logging.DEBUG)
        return

    resolved_log_file.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(resolved_log_file, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(_build_v4_formatter())
    setattr(file_handler, "_rph_v4_log_path", str(resolved_log_file))
    root.addHandler(file_handler)
    for logger_name in _V4_DETAIL_LOGGERS:
        logging.getLogger(logger_name).setLevel(logging.DEBUG)


def setup_logger(
    name: str = "ReactionProfileHunter",
    log_file: Path | None = None,
    level: int = logging.INFO,
    format_string: str | None = None,
) -> logging.Logger:
    """
    设置日志系统

    Args:
        name: Logger 名称
        log_file: 日志文件路径（可选）
        level: 日志级别
        format_string: 日志格式字符串 (Rich 模式下忽略)

    Returns:
        配置好的 Logger 对象
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # 清除已有的 handlers
    logger.handlers.clear()

    # 控制台输出 (优先使用 Rich)
    if HAS_RICH:
        if _rich_handler_cls is None:
            raise RuntimeError("RichHandler unavailable while HAS_RICH is True")
        console_handler: logging.Handler = _rich_handler_cls(
            console=_get_console(),
            rich_tracebacks=True,
            show_time=True,
            show_path=False,
            markup=True,
        )
        # RichHandler 自带格式化，通常不需要 formatter
    else:
        if format_string is None:
            format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        formatter = logging.Formatter(format_string)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)

    console_handler.setLevel(level)
    logger.addHandler(console_handler)
    logger.propagate = False

    # 文件输出（如果指定）- 始终使用标准格式
    if log_file:
        if format_string is None:
            format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_formatter = logging.Formatter(format_string)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    return logger


class LoggerMixin:
    """
    Logger 混入类

    为类提供便捷的 logger 属性
    """

    @property
    def logger(self) -> logging.Logger:
        """获取该类的 logger"""
        return logging.getLogger(f"{self.__class__.__module__}.{self.__class__.__name__}")
