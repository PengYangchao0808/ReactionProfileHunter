#!/usr/bin/env python3
import sys
import logging
import yaml
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from rph_core.utils.resource_utils import (
    resolve_executable_config,
    validate_executable,
    mem_to_mb,
    calc_orca_maxcore
)

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def test_config_loading():
    config_path = Path(__file__).parent / "config" / "defaults.yaml"
    logger.info(f"加载配置文件: {config_path}")
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def test_executable_resolution(config):
    logger.info("\n" + "="*60)
    logger.info("测试可执行文件解析")
    logger.info("="*60)

    programs = ['orca', 'gaussian', 'xtb', 'crest']

    for program in programs:
        logger.info(f"\n检查 {program.upper()}:")
        exe_config = resolve_executable_config(config, program, allow_path_search=False)
        logger.info(f"  路径: {exe_config.get('path')}")
        logger.info(f"  找到: {exe_config.get('found')}")
        logger.info(f"  来源: {exe_config.get('source')}")

        if exe_config.get('found'):
            is_valid = validate_executable(exe_config['path'], program)
            logger.info(f"  验证: {'✓ 通过' if is_valid else '✗ 失败'}")
        else:
            logger.warning(f"  未找到 {program} 可执行文件")


def test_resource_config(config):
    logger.info("\n" + "="*60)
    logger.info("测试资源配置")
    logger.info("="*60)

    res_cfg = config.get('resources', {})
    mem = res_cfg.get('mem', '32GB')
    nproc = res_cfg.get('nproc', 16)
    orca_safety = res_cfg.get('orca_maxcore_safety', 0.8)

    logger.info(f"\n内存: {mem}")
    logger.info(f"核数: {nproc}")
    logger.info(f"ORCA maxcore 安全系数: {orca_safety}")

    mem_mb = mem_to_mb(mem)
    logger.info(f"内存 (MB): {mem_mb}")

    maxcore = calc_orca_maxcore(mem, nproc, orca_safety)
    logger.info(f"ORCA maxcore: {maxcore} MB/core")


def test_first_smiles():
    logger.info("\n" + "="*60)
    logger.info("测试第一个 SMILES 示例")
    logger.info("="*60)

    csv_path = Path(__file__).parent / "reaxys_cleaned.csv"
    logger.info(f"加载测试数据: {csv_path}")

    import csv
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        first_row = next(reader)

    smiles = first_row['product_smiles_main']
    rx_id = first_row['rx_id']

    logger.info(f"\n反应 ID: {rx_id}")
    logger.info(f"产物 SMILES: {smiles}")
    logger.info("✓ 成功读取第一个示例")


def main():
    logger.info("开始路径配置验证")
    logger.info("="*60)

    config = test_config_loading()
    test_executable_resolution(config)
    test_resource_config(config)
    test_first_smiles()

    logger.info("\n" + "="*60)
    logger.info("验证完成")
    logger.info("="*60)


if __name__ == "__main__":
    main()
