"""
File I/O Utilities
===================

文件读写辅助工具
"""

from pathlib import Path
from typing import Tuple, List
import numpy as np
import logging

logger = logging.getLogger(__name__)


def read_xyz(xyz_file: Path) -> Tuple[np.ndarray, List[str]]:
    """
    读取 XYZ 文件

    Args:
        xyz_file: XYZ 文件路径

    Returns:
        (coordinates, symbols) 元组
        - coordinates: (N, 3) numpy 数组
        - symbols: 元素符号列表
    """
    if not xyz_file.exists():
        raise FileNotFoundError(f"XYZ 文件不存在: {xyz_file}")

    with open(xyz_file, 'r') as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())
    title = lines[1].strip()

    coords = []
    symbols = []

    for line in lines[2:2+n_atoms]:
        parts = line.strip().split()
        if len(parts) >= 4:
            symbols.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return np.array(coords), symbols


def write_xyz(xyz_file: Path, coordinates: np.ndarray, symbols: List[str],
              title: str = "", energy: float = None):
    """
    写入 XYZ 文件

    Args:
        xyz_file: 输出文件路径
        coordinates: (N, 3) 坐标数组
        symbols: 元素符号列表
        title: 标题行
        energy: 能量值（可选，会写入第二行）
    """
    n_atoms = len(symbols)

    with open(xyz_file, 'w') as f:
        f.write(f"{n_atoms}\n")

        # 写入标题行（可能包含能量）
        if energy is not None:
            f.write(f"{title} energy: {energy:.10f}\n")
        else:
            f.write(f"{title}\n")

        # 写入坐标
        for i, (symbol, coord) in enumerate(zip(symbols, coordinates)):
            f.write(f"{symbol:2s} {coord[0]:15.10f} {coord[1]:15.10f} {coord[2]:15.10f}\n")

    logger.debug(f"✓ XYZ 文件已写入: {xyz_file}")


def read_gjf(gjf_file: Path) -> Tuple[np.ndarray, List[str], int, int]:
    """
    读取 Gaussian 输入文件

    Args:
        gjf_file: GJF 文件路径

    Returns:
        (coordinates, symbols, charge, multiplicity) 元组
    """
    with open(gjf_file, 'r') as f:
        lines = f.readlines()

    # 跳过路由卡
    coord_start = 0
    for i, line in enumerate(lines):
        if line.strip().isdigit():
            charge_mult = line.strip().split()
            charge = int(charge_mult[0])
            multiplicity = int(charge_mult[1])
            coord_start = i + 1
            break

    # 读取坐标
    coords = []
    symbols = []

    for line in lines[coord_start:]:
        line = line.strip()
        if not line or line.startswith('--'):
            break
        parts = line.split()
        if len(parts) >= 4:
            symbols.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return np.array(coords), symbols, charge, multiplicity


def read_energy_from_gaussian(log_file: Path) -> float:
    """
    从 Gaussian 输出文件读取能量

    Args:
        log_file: Gaussian .log 文件路径

    Returns:
        能量值 (Hartree)
    """
    with open(log_file, 'r') as f:
        content = f.read()

    # 查找 SCF Done 行
    import re
    scf_match = re.search(r'SCF Done:\s+E\([^)]+\)\s*=\s*([\-\d.]+)', content)

    if scf_match:
        return float(scf_match.group(1))

    raise ValueError(f"无法从 {log_file} 提取能量")
