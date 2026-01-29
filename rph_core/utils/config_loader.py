"""
Configuration Loader
====================

配置文件加载工具
"""

import yaml
from pathlib import Path
from typing import Dict, Any
import logging

logger = logging.getLogger(__name__)


def load_config(config_path: Path) -> Dict[str, Any]:
    """
    加载 YAML 配置文件

    Args:
        config_path: 配置文件路径

    Returns:
        配置字典
    """
    if not config_path.exists():
        raise FileNotFoundError(f"配置文件不存在: {config_path}")

    with open(config_path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)

    logger.info(f"已加载配置文件: {config_path}")
    return config


def save_config(config: Dict[str, Any], output_path: Path):
    """
    保存配置到 YAML 文件

    Args:
        config: 配置字典
        output_path: 输出文件路径
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', encoding='utf-8') as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True)

    logger.info(f"已保存配置文件: {output_path}")


def merge_configs(*configs: Dict[str, Any]) -> Dict[str, Any]:
    """
    合并多个配置字典（后面的覆盖前面的）

    Args:
        *configs: 配置字典

    Returns:
        合并后的配置字典
    """
    result = {}
    for config in configs:
        result.update(config)
    return result


def get_executable_config(
    config: Dict[str, Any],
    program_key: str,
    default_path: str = ""
) -> Dict[str, str]:
    """
    从配置中提取可执行文件配置

    Args:
        config: 配置字典
        program_key: 程序键名（如 "orca", "gaussian"）
        default_path: 默认路径

    Returns:
        {"path": "...", "ld_library_path": "..."}
    """
    exe_cfg = config.get('executables', {}).get(program_key, {})

    return {
        'path': exe_cfg.get('path', default_path),
        'ld_library_path': exe_cfg.get('ld_library_path', None)
    }
