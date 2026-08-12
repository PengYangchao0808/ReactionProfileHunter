"""
File I/O Utilities
===================

文件读写辅助工具
"""

from functools import lru_cache
from pathlib import Path
from typing import List, Optional, Sequence, Tuple
import numpy as np
import logging

logger = logging.getLogger(__name__)


def _xyz_cache_key(path: Path) -> Tuple[str, int, int]:
    """构造 XYZ 解析缓存键；文件修改后自动失效。"""
    stat = Path(path).stat()
    return str(Path(path)), int(stat.st_mtime_ns), int(stat.st_size)


@lru_cache(maxsize=256)
def _read_xyz_uncached(path_str: str, mtime_ns: int, size: int) -> Tuple[np.ndarray, Tuple[str, ...]]:
    """Cached XYZ parse.  Keyed by (path, mtime_ns, size) so edits invalidate."""
    del mtime_ns, size
    lines = Path(path_str).read_text().splitlines()
    n_atoms = int(lines[0].strip())
    coords = []
    symbols = []
    for line in lines[2:2 + n_atoms]:
        parts = line.strip().split()
        if len(parts) >= 4:
            symbols.append(parts[0])
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
    return np.array(coords, dtype=float, copy=True), tuple(symbols)


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
    path_key = _xyz_cache_key(Path(xyz_file))
    coords, symbols = _read_xyz_uncached(*path_key)
    return np.array(coords, copy=True), list(symbols)


def read_xyz_many(xyz_files: Sequence[Path]) -> List[Tuple[np.ndarray, List[str]]]:
    """Batch XYZ read; each file is resolved through the shared cache."""
    return [read_xyz(Path(path)) for path in xyz_files]


def clear_xyz_cache() -> None:
    """Drop the in-memory XYZ parse cache (used by tests and long-running hosts)."""
    _read_xyz_uncached.cache_clear()


def write_xyz(xyz_file: Path, coordinates: np.ndarray, symbols: List[str],
              title: str = "", energy: Optional[float] = None):
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
        for symbol, coord in zip(symbols, coordinates):
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


def split_multixyz(xyz_file: Path, output_dir: Path, prefix: str = "frame") -> List[Path]:
    """
    将 multi-frame XYZ 文件拆分为单个 frame 文件。

    xTB PATH 输出的 xtbpath.xyz 包含多个 XYZ 帧（由 --- 或连续 natoms 行分隔）。
    此函数将其逐帧拆分为 ``output_dir/{prefix}_{idx:04d}.xyz`` 格式的文件。

    Args:
        xyz_file: multi-frame XYZ 文件路径。
        output_dir: 输出目录（会被创建）。
        prefix: 帧文件前缀（默认 "frame"）。

    Returns:
        List[Path]: 按帧序排列的单个帧 XYZ 文件路径列表。
        如果文件为空或无法解析则返回空列表。
    """
    if not xyz_file.exists():
        logger.error(f"Multi-XYZ 文件不存在: {xyz_file}")
        return []

    raw = xyz_file.read_text()
    lines = raw.splitlines()
    if not lines:
        logger.warning(f"multi-XYZ 文件为空: {xyz_file}")
        return []

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    frame_paths: List[Path] = []
    i = 0
    frame_idx = 0
    while i < len(lines):
        line = lines[i].strip()
        # 跳过空行和分隔线
        if not line or line.startswith("---"):
            i += 1
            continue
        try:
            natoms = int(line)
        except ValueError:
            i += 1
            continue

        if i + natoms + 2 > len(lines):
            logger.warning(
                "multi-XYZ 在第 %d 帧处截断（期望 %d 行，剩余 %d 行）",
                frame_idx, natoms + 2, len(lines) - i,
            )
            break

        frame_lines = lines[i:i + natoms + 2]
        frame_path = output_dir / f"{prefix}_{frame_idx:04d}.xyz"
        frame_path.write_text("\n".join(frame_lines))
        frame_paths.append(frame_path)

        frame_idx += 1
        i += natoms + 2

    logger.info(
        "拆分 multi-XYZ: %s → %d frames → %s",
        xyz_file.name, len(frame_paths), output_dir,
    )
    return frame_paths
