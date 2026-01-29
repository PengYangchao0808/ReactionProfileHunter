"""
QC Runner Module
================
量子化学计算执行器：错误处理、重试逻辑、超时管理

Author: QC Descriptors Team
Date: 2026-01-13
Purpose: 替代旧版 Bash 脚本中的错误处理逻辑
"""

import logging
import time
import signal
import subprocess
from pathlib import Path
from typing import Callable, Optional, Dict, Any, TypeVar, ParamSpec
from functools import wraps
from enum import Enum

logger = logging.getLogger(__name__)

P = ParamSpec('P')
T = TypeVar('T')


class QCFailureType(Enum):
    """QC 计算失败类型"""
    TRANSIENT = "transient"
    PERMANENT = "permanent"
    TIMEOUT = "timeout"
    OUT_OF_MEMORY = "out_of_memory"
    CONVERGENCE_FAILURE = "convergence_failure"


class QCRunnerError(Exception):
    """QC Runner 基础异常"""
    def __init__(self, message: str, failure_type: QCFailureType):
        self.message = message
        self.failure_type = failure_type
        super().__init__(f"[{failure_type.value}] {message}")


class QCTimeoutError(QCRunnerError):
    """超时异常"""
    def __init__(self, message: str = "QC calculation timed out"):
        super().__init__(message, QCFailureType.TIMEOUT)


class QCMemoryError(QCRunnerError):
    """内存不足异常"""
    def __init__(self, message: str = "QC calculation ran out of memory"):
        super().__init__(message, QCFailureType.OUT_OF_MEMORY)


class QCConvergenceError(QCRunnerError):
    """收敛失败异常"""
    def __init__(self, message: str = "QC calculation failed to converge"):
        super().__init__(message, QCFailureType.CONVERGENCE_FAILURE)


def analyze_log_for_errors(log_content: str) -> Optional[QCFailureType]:
    """
    分析 QC 日志，识别错误类型

    Args:
        log_content: 日志文件内容

    Returns:
        错误类型，如果无法识别返回 None
    """
    log_content_lower = log_content.lower()

    if "out of memory" in log_content_lower or "memory allocation failed" in log_content_lower:
        return QCFailureType.OUT_OF_MEMORY

    if "convergence failure" in log_content_lower or "failed to converge" in log_content_lower:
        return QCFailureType.CONVERGENCE_FAILURE

    if "scf not converged" in log_content_lower:
        return QCFailureType.CONVERGENCE_FAILURE

    return None


def is_retryable(error: QCRunnerError) -> bool:
    """
    判断错误是否可重试

    Args:
        error: QC Runner 错误

    Returns:
        是否可重试
    """
    retryable_types = {
        QCFailureType.TRANSIENT,
        QCFailureType.CONVERGENCE_FAILURE
    }
    return error.failure_type in retryable_types


class RetryConfig:
    """重试配置"""
    def __init__(
        self,
        max_retries: int = 3,
        initial_backoff: float = 5.0,
        max_backoff: float = 300.0,
        backoff_multiplier: float = 2.0
    ):
        self.max_retries = max_retries
        self.initial_backoff = initial_backoff
        self.max_backoff = max_backoff
        self.backoff_multiplier = backoff_multiplier

    def get_backoff(self, attempt: int) -> float:
        """计算退避时间（指数退避）"""
        backoff = self.initial_backoff * (self.backoff_multiplier ** (attempt - 1))
        return min(backoff, self.max_backoff)


def with_retry(config: Optional[RetryConfig] = None, cleanup: Optional[Callable] = None):
    """
    重试装饰器

    Args:
        config: 重试配置（默认：3次重试，指数退避）
        cleanup: 失败时的清理函数
    """
    if config is None:
        config = RetryConfig()

    def decorator(func: Callable[P, T]) -> Callable[P, T]:
        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            last_exception: Optional[Exception] = None

            for attempt in range(1, config.max_retries + 1):
                try:
                    logger.info(f"执行 {func.__name__} (尝试 {attempt}/{config.max_retries})")
                    result = func(*args, **kwargs)
                    if attempt > 1:
                        logger.info(f"✓ {func.__name__} 在第 {attempt} 次尝试成功")
                    return result

                except QCRunnerError as e:
                    last_exception = e
                    logger.warning(f"第 {attempt} 次尝试失败: {e}")

                    if not is_retryable(e):
                        logger.error(f"错误类型 {e.failure_type.value} 不可重试")
                        raise

                    if attempt < config.max_retries:
                        backoff = config.get_backoff(attempt)
                        logger.info(f"等待 {backoff:.1f} 秒后重试...")
                        time.sleep(backoff)

                except Exception as e:
                    last_exception = e
                    logger.error(f"未捕获的异常: {type(e).__name__}: {e}")
                    raise

            if cleanup:
                try:
                    cleanup()
                except Exception as e:
                    logger.error(f"清理失败: {e}")

            raise QCRunnerError(
                f"{func.__name__} 在 {config.max_retries} 次尝试后失败: {last_exception}",
                QCFailureType.PERMANENT
            ) from last_exception

        return wrapper
    return decorator


def run_with_timeout(
    cmd: list,
    timeout: int,
    cwd: Optional[Path] = None,
    env: Optional[Dict[str, str]] = None,
    log_file: Optional[Path] = None
) -> subprocess.CompletedProcess:
    """
    带超时的子进程执行

    Args:
        cmd: 命令列表
        timeout: 超时时间（秒）
        cwd: 工作目录
        env: 环境变量
        log_file: 日志文件路径

    Returns:
        完成的子进程对象

    Raises:
        QCTimeoutError: 超时
        subprocess.CalledProcessError: 命令执行失败
    """
    logger.debug(f"执行命令: {' '.join(cmd)}")

    with subprocess.Popen(
        cmd,
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1
    ) as process:
        try:
            stdout, _ = process.communicate(timeout=timeout)

            if log_file:
                log_file.parent.mkdir(parents=True, exist_ok=True)
                with open(log_file, 'w') as f:
                    f.write(stdout)

            if process.returncode != 0:
                error_type = analyze_log_for_errors(stdout)
                if error_type:
                    if error_type == QCFailureType.OUT_OF_MEMORY:
                        raise QCMemoryError()
                    elif error_type == QCFailureType.CONVERGENCE_FAILURE:
                        raise QCConvergenceError()

                raise subprocess.CalledProcessError(
                    process.returncode, cmd, output=stdout
                )

            return subprocess.CompletedProcess(
                args=cmd,
                returncode=process.returncode,
                stdout=stdout,
                stderr=""
            )

        except subprocess.TimeoutExpired:
            process.kill()
            stdout, _ = process.communicate()

            if log_file:
                log_file.parent.mkdir(parents=True, exist_ok=True)
                with open(log_file, 'w') as f:
                    f.write(stdout)

            raise QCTimeoutError(
                f"命令在 {timeout} 秒后超时: {' '.join(cmd)}"
            )


class QCJobStatus(Enum):
    """QC 任务状态"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    TIMEOUT = "timeout"


class QCJob:
    """QC 计算任务"""
    def __init__(
        self,
        job_id: str,
        cmd: list,
        work_dir: Path,
        timeout: int = 3600,
        config: Optional[Dict[str, Any]] = None
    ):
        self.job_id = job_id
        self.cmd = cmd
        self.work_dir = Path(work_dir)
        self.timeout = timeout
        self.config = config or {}
        self.status = QCJobStatus.PENDING
        self.result: Optional[subprocess.CompletedProcess] = None
        self.error: Optional[Exception] = None
        self.log_file = self.work_dir / f"{job_id}.log"

    def run(self) -> subprocess.CompletedProcess:
        """执行任务"""
        self.status = QCJobStatus.RUNNING
        logger.info(f"开始任务 {self.job_id}")

        try:
            self.result = run_with_timeout(
                self.cmd,
                timeout=self.timeout,
                cwd=self.work_dir,
                log_file=self.log_file
            )
            self.status = QCJobStatus.COMPLETED
            logger.info(f"任务 {self.job_id} 完成")
            return self.result

        except QCTimeoutError as e:
            self.status = QCJobStatus.TIMEOUT
            self.error = e
            logger.error(f"任务 {self.job_id} 超时: {e}")
            raise

        except Exception as e:
            self.status = QCJobStatus.FAILED
            self.error = e
            logger.error(f"任务 {self.job_id} 失败: {e}")
            raise

    def cleanup(self):
        """清理临时文件"""
        if self.log_file.exists():
            self.log_file.unlink()
        logger.debug(f"任务 {self.job_id} 清理完成")
