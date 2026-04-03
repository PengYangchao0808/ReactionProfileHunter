"""
Resource Utilities Module
=======================
资源管理工具：内存换算、路径解析、环境变量处理

Author: QC Descriptors Team
Date: 2026-01-13
Purpose: 替代旧版 Bash 脚本中的资源管理逻辑
"""

import re
import os
import shutil
import logging
from pathlib import Path
from typing import Optional, Dict, Tuple

logger = logging.getLogger(__name__)


def get_project_root() -> Path:
    """
    获取项目根目录的绝对路径 (本文件 rph_core/utils/resource_utils.py 的上级两级)

    Returns:
        项目根目录的绝对路径
    """
    current_file = Path(__file__).resolve()
    project_root = current_file.parent.parent.parent
    return project_root


# ========== 内存换算函数 ==========
def mem_to_mb(mem_str: str) -> int:
    """
    将内存字符串转换为 MB

    支持格式:
    - "32GB", "16GB", "32G" -> 32 * 1024 MB
    - "4096MB", "2048MB", "4096M" -> 直接返回 MB 值
    - 纯数字 (如 "16") -> 视为 GB

    Args:
        mem_str: 内存字符串

    Returns:
        内存大小（MB）
    """
    if not mem_str:
        logger.warning("内存字符串为空，返回默认值 4000 MB")
        return 4000

    m = mem_str.strip().upper()
    match = re.match(r'^(\d+)\s*(GB?|MB?)?$', m)

    if match:
        val = int(match.group(1))
        unit = match.group(2) or 'GB'

        if unit.startswith('G'):
            return val * 1024
        else:
            return val

    logger.warning(f"无法解析内存字符串: {mem_str}，返回默认值 4000 MB")
    return 4000


def calc_orca_maxcore(mem: str, nproc: int, safety_factor: float = 0.8) -> int:
    """
    计算 ORCA 的 maxcore 参数（MB per core）

    旧脚本逻辑（Config.sh 第90-96行）:
    ```bash
    if [[ "$MEM" =~ ^([0-9]+)[[:space:]]*GB$ ]]; then
      ORCA_MAXCORE_MB="$(( ${BASH_REMATCH[1]} * 1024 * 2 ))"
    fi
    ```

    Args:
        mem: 内存字符串（如 "32GB"）
        nproc: 并行核数
        safety_factor: 安全系数（避免内存溢出）

    Returns:
        maxcore 值（MB per core）
    """
    total_mb = mem_to_mb(mem)
    maxcore = int(total_mb * safety_factor / nproc)

    logger.info(f"ORCA maxcore 计算: {mem} / {nproc} cores * {safety_factor} = {maxcore} MB/core")
    return maxcore


# ========== 路径解析函数 ==========
def find_executable(
    program_name: str,
    config_path: Optional[str] = None,
    env_vars: Optional[list] = None,
    allow_path_search: bool = True
) -> Optional[Path]:
    """
    查找可执行文件路径（多级回退策略）

    查找顺序（参考旧脚本 Auto_Calc_20251227.sh 的 ensure_exec）:
    1. 配置文件中的绝对路径
    2. 环境变量
    3. 系统 PATH

    Args:
        program_name: 程序名称（如 "orca"）
        config_path: 配置文件中的路径
        env_vars: 要检查的环境变量列表（如 ["ORCA_PATH", "ORCA_BIN"]）
        allow_path_search: 是否允许在系统 PATH 中查找

    Returns:
        可执行文件的 Path 对象，如果找不到返回 None
    """
    # 1. 配置文件路径
    if config_path:
        path = Path(config_path)
        if path.exists() and path.is_file() and os.access(path, os.X_OK):
            logger.info(f"找到可执行文件（配置）: {path}")
            return path
        else:
            logger.warning(f"配置路径不存在或不可执行: {config_path}")

    # 2. 环境变量
    if env_vars:
        for env_var in env_vars:
            env_path = os.environ.get(env_var)
            if env_path:
                path = Path(env_path)
                if path.exists() and path.is_file() and os.access(path, os.X_OK):
                    logger.info(f"找到可执行文件（环境变量 {env_var}）: {path}")
                    return path
                else:
                    logger.warning(f"环境变量 {env_var}={env_path} 指向的文件不存在或不可执行")

    # 3. 系统 PATH
    if allow_path_search:
        which_result = shutil.which(program_name)
        if which_result:
            path = Path(which_result)
            logger.info(f"找到可执行文件（PATH）: {path}")
            return path

    logger.error(f"未找到可执行文件: {program_name}")
    return None


def setup_ld_library_path(ld_library_paths: list):
    """
    设置 LD_LIBRARY_PATH 环境变量（ORCA MPI 支持）

    Args:
        ld_library_paths: 要添加的路径列表
    """
    if not ld_library_paths:
        return

    current_path = os.environ.get('LD_LIBRARY_PATH', '')
    new_paths = []

    for path_str in ld_library_paths:
        # Support colon-separated strings (e.g. "/opt/openmpi418/lib:/opt/software/orca")
        individual_paths = path_str.split(':') if ':' in path_str else [path_str]
        for single_path in individual_paths:
            single_path = single_path.strip()
            if not single_path:
                continue
            path = Path(single_path)
            if path.exists():
                new_paths.append(str(path))
            else:
                logger.warning(f"库路径不存在: {single_path}")

    if new_paths:
        if current_path:
            new_path = ":".join(new_paths) + ":" + current_path
        else:
            new_path = ":".join(new_paths)

        os.environ['LD_LIBRARY_PATH'] = new_path
        logger.info(f"已更新 LD_LIBRARY_PATH: {new_path}")


# ========== 配置派生函数 ==========
def resolve_executable_config(
    config: dict,
    program_key: str,
    env_vars: Optional[list] = None,
    allow_path_search: bool = True
) -> Dict:
    exe_cfg = config.get('executables', {}).get(program_key, {})
    res_cfg = config.get('resources', {})

    exe_path = find_executable(
        program_key,
        config_path=exe_cfg.get('path'),
        env_vars=env_vars,
        allow_path_search=allow_path_search
    )

    source = "not_found"
    if exe_path:
        config_path = exe_cfg.get('path')
        if config_path and exe_path == Path(config_path):
            source = "config"
        elif env_vars:
            for env_var in env_vars:
                env_path = os.environ.get(env_var)
                if env_path and exe_path == Path(env_path):
                    source = "env"
                    break
        else:
            source = "path"

    ld_path_str = exe_cfg.get('ld_library_path')

    if ld_path_str:
        setup_ld_library_path([ld_path_str])

    return {
        'path': exe_path,
        'ld_library_path': ld_path_str,
        'found': exe_path is not None,
        'source': source
    }


# ========== 验证函数 ==========
def validate_executable(path: Path, program_name: str = "") -> bool:
    """
    验证可执行文件是否可用

    Args:
        path: 可执行文件路径
        program_name: 程序名称（用于错误提示）

    Returns:
        是否可用
    """
    if not path:
        logger.error(f"{program_name}: 路径为空")
        return False

    if not path.exists():
        logger.error(f"{program_name}: 文件不存在: {path}")
        return False

    if not path.is_file():
        logger.error(f"{program_name}: 不是文件: {path}")
        return False

    if not os.access(path, os.X_OK):
        logger.error(f"{program_name}: 文件不可执行: {path}")
        return False

    return True
