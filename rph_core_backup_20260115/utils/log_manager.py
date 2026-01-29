"""
Logging Manager
================

统一的日志管理系统 (Rich Enhanced)
"""

import logging
import sys
from pathlib import Path
from typing import Optional

try:
    from rich.logging import RichHandler
    HAS_RICH = True
except ImportError:
    HAS_RICH = False

def setup_logger(
    name: str = "ReactionProfileHunter",
    log_file: Optional[Path] = None,
    level: int = logging.INFO,
    format_string: Optional[str] = None
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
        console_handler = RichHandler(
            rich_tracebacks=True,
            show_time=True,
            show_path=False,
            markup=True
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
