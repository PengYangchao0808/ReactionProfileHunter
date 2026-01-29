"""
Geometry Tools - Universal Geometric Operations
===============================================

通用几何操作工具（非业务逻辑）

Author: QCcalc Team
Date: 2026-01-09
"""

import numpy as np
from typing import Tuple, List, Optional
import logging
from pathlib import Path
import re

from rph_core.utils.file_io import read_xyz

logger = logging.getLogger(__name__)


class GeometryUtils:
    """
    通用几何工具类

    提供各种几何计算和转换功能
    这些是底层工具，不包含特定的化学反应逻辑
    """

    @staticmethod
    def calculate_distance(
        coords: np.ndarray,
        atom_i: int,
        atom_j: int
    ) -> float:
        """
        计算两个原子之间的距离

        Args:
            coords: 坐标数组 (N, 3)
            atom_i: 第一个原子索引
            atom_j: 第二个原子索引

        Returns:
            距离 (Å)
        """
        if atom_i >= len(coords) or atom_j >= len(coords):
            raise IndexError(f"原子索引超出范围: {atom_i}, {atom_j}")

        vec = coords[atom_i] - coords[atom_j]
        return float(np.linalg.norm(vec))

    @staticmethod
    def calculate_angle(
        coords: np.ndarray,
        atom_i: int,
        atom_j: int,
        atom_k: int
    ) -> float:
        """
        计算键角 (i-j-k)

        Args:
            coords: 坐标数组 (N, 3)
            atom_i: 第一个原子索引
            atom_j: 中心原子索引
            atom_k: 第三个原子索引

        Returns:
            角度 (度)
        """
        # 向量 ji 和 jk
        vec_ji = coords[atom_i] - coords[atom_j]
        vec_jk = coords[atom_k] - coords[atom_j]

        # 计算角度
        cosine = np.dot(vec_ji, vec_jk) / (
            np.linalg.norm(vec_ji) * np.linalg.norm(vec_jk)
        )

        # 防止数值误差
        cosine = np.clip(cosine, -1.0, 1.0)

        angle_rad = np.arccos(cosine)
        angle_deg = np.degrees(angle_rad)

        return float(angle_deg)

    @staticmethod
    def calculate_dihedral(
        coords: np.ndarray,
        atom_i: int,
        atom_j: int,
        atom_k: int,
        atom_l: int
    ) -> float:
        """
        计算二面角 (i-j-k-l)

        Args:
            coords: 坐标数组 (N, 3)
            atom_i, atom_j, atom_k, atom_l: 四个原子索引

        Returns:
            二面角 (度)
        """
        # 参考: https://en.wikipedia.org/wiki/Dihedral_angle

        b1 = coords[atom_j] - coords[atom_i]
        b2 = coords[atom_k] - coords[atom_j]
        b3 = coords[atom_l] - coords[atom_k]

        # 法向量
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)

        # 计算二面角
        cosine = np.dot(n1, n2) / (
            np.linalg.norm(n1) * np.linalg.norm(n2)
        )
        cosine = np.clip(cosine, -1.0, 1.0)

        angle_rad = np.arccos(cosine)
        angle_deg = np.degrees(angle_rad)

        # 确定符号
        # 使用叉积来判断方向
        m = np.cross(n1, b2)
        x = np.dot(m, n2)

        if x < 0:
            angle_deg = -angle_deg

        return float(angle_deg)

    @staticmethod
    def get_center_of_mass(
        coords: np.ndarray,
        symbols: Optional[List[str]] = None,
        masses: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        计算质心

        Args:
            coords: 坐标数组 (N, 3)
            symbols: 原子符号（可选，用于查找质量）
            masses: 原子质量数组（可选，覆盖symbols）

        Returns:
            质心坐标 (3,)
        """
        n_atoms = len(coords)

        if masses is None:
            if symbols is None:
                # 假设所有质量为1
                masses = np.ones(n_atoms)
            else:
                # 从原子符号获取质量
                masses = GeometryUtils._get_atomic_masses(symbols)

        # 计算质心
        total_mass = np.sum(masses)
        center_of_mass = np.average(coords, axis=0, weights=masses)

        return center_of_mass

    @staticmethod
    def _get_atomic_masses(symbols: List[str]) -> np.ndarray:
        """
        获取原子质量

        Args:
            symbols: 原子符号列表

        Returns:
            质量数组 (原子质量单位)
        """
        # 常见原子质量（简化版）
        atomic_masses = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
            'F': 18.998, 'P': 30.974, 'S': 32.065, 'Cl': 35.453,
            'Br': 79.904, 'I': 126.904
        }

        masses = []
        for symbol in symbols:
            # 处理特殊情况（如电荷标记）
            symbol_clean = symbol.split('.')[0]
            mass = atomic_masses.get(symbol_clean, 12.011)  # 默认碳质量
            masses.append(mass)

        return np.array(masses)

    @staticmethod
    def rotate_coords(
        coords: np.ndarray,
        axis: np.ndarray,
        angle_deg: float,
        center: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        绕轴旋转坐标

        Args:
            coords: 原始坐标 (N, 3)
            axis: 旋转轴向量 (3,)
            angle_deg: 旋转角度（度）
            center: 旋转中心（可选，默认为质心）

        Returns:
            旋转后的坐标 (N, 3)
        """
        coords = coords.copy()

        if center is None:
            center = np.mean(coords, axis=0)

        # 将坐标平移到以center为原点
        coords_shifted = coords - center

        # 归一化旋转轴
        axis = axis / np.linalg.norm(axis)

        # 旋转矩阵（Rodrigues公式）
        angle_rad = np.radians(angle_deg)
        cos_theta = np.cos(angle_rad)
        sin_theta = np.sin(angle_rad)

        # 叉积矩阵
        K = np.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
        ])

        # 旋转矩阵
        R = np.eye(3) + sin_theta * K + (1 - cos_theta) * np.dot(K, K)

        # 应用旋转
        coords_rotated = np.dot(coords_shifted, R.T)

        # 平移回原位置
        coords_final = coords_rotated + center

        return coords_final

    @staticmethod
    def translate_coords(
        coords: np.ndarray,
        vector: np.ndarray
    ) -> np.ndarray:
        """
        平移坐标

        Args:
            coords: 原始坐标 (N, 3)
            vector: 平移向量 (3,)

        Returns:
            平移后的坐标 (N, 3)
        """
        return coords + vector

    @staticmethod
    def calculate_rmsd(
        coords1: np.ndarray,
        coords2: np.ndarray
    ) -> float:
        """
        计算RMSD（均方根偏差）

        Args:
            coords1: 第一个坐标集 (N, 3)
            coords2: 第二个坐标集 (N, 3)

        Returns:
            RMSD (Å)
        """
        if coords1.shape != coords2.shape:
            raise ValueError("坐标数组形状不匹配")

        diff = coords1 - coords2
        mse = np.mean(np.sum(diff**2, axis=1))
        rmsd = np.sqrt(mse)

        return float(rmsd)


class BondOperations:
    """
    键操作类（通用工具）

    提供键拉伸、断裂等操作
    """

    @staticmethod
    def stretch_bond(
        coords: np.ndarray,
        atom_i: int,
        atom_j: int,
        target_length: float,
        fix_atoms: Optional[List[int]] = None
    ) -> np.ndarray:
        """
        拉伸单个键到目标长度

        Args:
            coords: 原始坐标 (N, 3)
            atom_i: 第一个原子索引
            atom_j: 第二个原子索引
            target_length: 目标键长 (Å)
            fix_atoms: 固定不动的原子索引列表

        Returns:
            拉伸后的坐标 (N, 3)
        """
        coords = coords.copy()

        # 当前向量
        vec = coords[atom_j] - coords[atom_i]
        current_length = np.linalg.norm(vec)

        if current_length == 0:
            raise ValueError(f"原子 {atom_i} 和 {atom_j} 重合")

        # 拉伸因子
        scale_factor = target_length / current_length

        # 新向量
        new_vec = vec * scale_factor

        # 更新位置（沿键方向移动）
        if fix_atoms is None or atom_j not in fix_atoms:
            coords[atom_j] = coords[atom_i] + new_vec

        return coords


class CoordinateExtractor:
    """
    从各种量子化学软件输出文件中提取坐标

    参考旧脚本 SP-20251226.sh 的多级回退策略
    """

    @staticmethod
    def extract_coords_from_orca_out(out_file: Path) -> Optional[np.ndarray]:
        """
        从 ORCA .out 文件中提取坐标（多级回退策略）

        参考旧脚本 SP-20251226.sh 第305-345行:
        1. CARTESIAN COORDINATES (ANGSTROEM/ANGSTROM)
        2. CARTESIAN COORDINATES (A.U.) + Bohr→Å转换
        3. * xyz 块

        Args:
            out_file: ORCA 输出文件

        Returns:
            坐标数组 (N, 3)，如果提取失败返回 None
        """
        content = out_file.read_text()

        # 方案1: CARTESIAN COORDINATES (ANGSTROEM)
        coords = CoordinateExtractor._parse_angstrom_block(content)
        if coords is not None:
            logger.info(f"从 CARTESIAN COORDINATES (ANGSTROEM) 提取坐标: {out_file.name}")
            return coords

        # 方案2: CARTESIAN COORDINATES (A.U.) + Bohr→Å转换
        coords = CoordinateExtractor._parse_au_block(content)
        if coords is not None:
            logger.info(f"从 CARTESIAN COORDINATES (A.U.) 提取坐标并转换单位: {out_file.name}")
            return coords

        # 方案3: * xyz 块
        coords = CoordinateExtractor._parse_xyz_block(content)
        if coords is not None:
            logger.info(f"从 * xyz 块提取坐标: {out_file.name}")
            return coords

        logger.warning(f"无法从 ORCA .out 提取坐标: {out_file.name}")
        return None

    @staticmethod
    def _parse_angstrom_block(content: str) -> Optional[np.ndarray]:
        r"""
        解析 CARTESIAN COORDINATES (ANGSTROEM) 块

        旧脚本 awk 逻辑:
        ```
        awk -v RS='' '/CARTESIAN COORDINATES[[:space:]]*\((ANGSTROEM|ANGSTROM)\)/{block=$0} END{print block}'
        ```
        """
        pattern = r'CARTESIAN COORDINATES\s*\(ANGSTROEM\).*?\n((?:[A-Za-z]+\s+[\-]?\d+\.\d+\s+[\-]?\d+\.\d+\s+[\-]?\d+\.\d+\s+[\-]?\d+\.\d+\n)+)'

        match = re.search(pattern, content, re.IGNORECASE | re.DOTALL)
        if not match:
            return None

        lines = match.group(1).strip().split('\n')
        coords = []

        for line in lines:
            parts = line.strip().split()
            # 格式：atom x y z (可能带电荷）
            if len(parts) >= 4:
                try:
                    x, y, z = float(parts[-3]), float(parts[-2]), float(parts[-1])
                    coords.append([x, y, z])
                except ValueError:
                    continue

        return np.array(coords) if coords else None

    @staticmethod
    def _parse_au_block(content: str) -> Optional[np.ndarray]:
        """
        解析 CARTESIAN COORDINATES (A.U.) 块并转换单位

        Bohr→Å 转换因子：0.529177210903

        旧脚本 awk 逻辑:
        ```
        awk -v k=0.529177210903 '...'
        ```
        """
        pattern = r'CARTESIAN COORDINATES\s*\(A\.U\.\).*?\n((?:[A-Za-z]+\s+[\-]?\d+\.\d+\s+[\-]?\d+\.\d+\s+[\-]?\d+\.\d+\s+[\-]?\d+\.\d+\n)+)'

        match = re.search(pattern, content, re.IGNORECASE | re.DOTALL)
        if not match:
            return None

        lines = match.group(1).strip().split('\n')
        coords = []

        bohr_to_angstrom = 0.529177210903

        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    x = float(parts[-3]) * bohr_to_angstrom
                    y = float(parts[-2]) * bohr_to_angstrom
                    z = float(parts[-1]) * bohr_to_angstrom
                    coords.append([x, y, z])
                except ValueError:
                    continue

        return np.array(coords) if coords else None

    @staticmethod
    def _parse_xyz_block(content: str) -> Optional[np.ndarray]:
        r"""
        解析 * xyz 块

        旧脚本 awk 逻辑:
        ```
        awk '/^\* *xyz/i {flag=1; next}
            /^\* / && flag {exit}
            flag && NF>=4 {print}' "$src"
        ```
        """
        # 查找 * xyz ... * 之间的坐标块
        pattern = r'\* xyz\s+(\d+)\s+(\d+)((?:[^\*]|\n)+?)(?=\*|$)'

        match = re.search(pattern, content, re.IGNORECASE | re.DOTALL)
        if not match:
            return None

        coords_text = match.group(3).strip()
        coords = []

        for line in coords_text.split('\n'):
            parts = line.strip().split()
            if len(parts) >= 4:
                # 格式：element x y z
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    coords.append([x, y, z])
                except ValueError:
                    continue

        return np.array(coords) if coords else None

    @staticmethod
    def get_charge_spin_from_orca_out(out_file: Path) -> Tuple[int, int]:
        """
        从 ORCA .out 文件中提取电荷和自旋

        参考旧脚本 SP-20251226.sh 第348-365行:
        1. 优先 * xyz 行
        2. 回退：Total Charge / Multiplicity
        3. 兜底：0 1

        Args:
            out_file: ORCA 输出文件

        Returns:
            (charge, spin) 元组
        """
        content = out_file.read_text()
        charge, spin = 0, 1

        # 方案1: * xyz 行
        match = re.search(r'\*\s*xyz\s+(\-?\d+)\s+(\d+)', content, re.IGNORECASE)
        if match:
            charge = int(match.group(1))
            spin = int(match.group(2))
            logger.debug(f"从 * xyz 行提取电荷/自旋: {charge}/{spin}")
            return charge, spin

        # 方案2: Total Charge / Multiplicity
        charge_match = re.search(r'Total\s+Charge\s*:\s*(\-?\d+)', content, re.IGNORECASE)
        spin_match = re.search(r'Multiplicity\s*:\s*(\d+)', content, re.IGNORECASE)

        if charge_match:
            charge = int(charge_match.group(1))
        if spin_match:
            spin = int(spin_match.group(1))

        logger.debug(f"提取电荷/自旋: {charge}/{spin}")
        return charge, spin


class LogParser:
    """
    鲁棒的量子化学日志解析器

    从 Gaussian (.log) 和 ORCA (.out) 文件中提取最后一次收敛的几何结构

    参考旧脚本 SP-20251226.sh 的多级回退策略
    """

    @staticmethod
    def extract_last_converged_coords(
        log_file: Path,
        engine_type: str = 'auto'
    ) -> Tuple[Optional[np.ndarray], Optional[List[str]], Optional[str]]:
        """
        从日志文件中提取最后一次收敛的几何结构

        Args:
            log_file: 日志文件路径 (.log 或 .out)
            engine_type: 引擎类型 ('gaussian', 'orca', 'auto')

        Returns:
            (coordinates, symbols, error_message) 元组
            - coordinates: (N, 3) 坐标数组
            - symbols: 元素符号列表
            - error_message: 如果失败，返回错误信息
        """
        if not log_file.exists():
            return None, None, f"文件不存在: {log_file}"

        suffix = log_file.suffix.lower()

        # 自动检测引擎类型
        if engine_type == 'auto':
            if suffix == '.log':
                engine_type = 'gaussian'
            elif suffix == '.out':
                engine_type = 'orca'
            elif suffix == '.xyz':
                engine_type = 'xyz'
            else:
                return None, None, f"无法识别的文件后缀: {suffix}"


        if engine_type == 'xyz':
            try:
                coords, symbols = read_xyz(log_file)
            except Exception as e:
                return None, None, f"读取 XYZ 失败: {e}"
            return coords, symbols, None

        if engine_type == 'xyz':
            try:
                coords, symbols = read_xyz(log_file)
            except Exception as e:
                return None, None, f"读取 XYZ 失败: {e}"
            return coords, symbols, None

        try:
            content = log_file.read_text()
        except Exception as e:
            return None, None, f"读取文件失败: {e}"

        if engine_type == 'gaussian':
            return LogParser._parse_gaussian_log(content)
        elif engine_type == 'orca':
            coords, symbols, error = LogParser._parse_orca_out(content)
            if coords is not None:
                return coords, symbols, error
            inp_path = log_file.with_suffix('.inp')
            if inp_path.exists():
                inp_content = inp_path.read_text(errors='ignore')
                coords, symbols, inp_error = LogParser._parse_orca_input(inp_content)
                if coords is not None:
                    logger.info("从 ORCA .inp 输入块提取坐标")
                    return coords, symbols, None
                if inp_error:
                    return None, None, inp_error
            return coords, symbols, error
        else:
            return None, None, f"不支持的引擎类型: {engine_type}"

    @staticmethod
    def _parse_gaussian_log(content: str) -> Tuple[Optional[np.ndarray], Optional[List[str]], Optional[str]]:
        """
        解析 Gaussian .log 文件

        逻辑：
        1. 搜索最后一个 "Standard orientation" 块
        2. 回退到 "Input orientation"（如果 Standard 不存在）
        3. 提取坐标和元素符号
        """
        lines = content.splitlines()
        atomic_to_symbol = {
            1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
            9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
            16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 26: 'Fe', 35: 'Br', 53: 'I'
        }

        def _find_last_index(keyword: str) -> Optional[int]:
            keyword_lower = keyword.lower()
            indices = [index for index, line in enumerate(lines) if keyword_lower in line.lower()]
            return indices[-1] if indices else None

        def _extract_from_index(start_index: int) -> Tuple[Optional[np.ndarray], Optional[List[str]]]:
            index = start_index + 1
            while index < len(lines) and not re.match(r'\s*-{5,}\s*$', lines[index]):
                index += 1
            if index >= len(lines):
                return None, None

            index += 1
            while index < len(lines) and not re.match(r'\s*-{5,}\s*$', lines[index]):
                index += 1
            if index >= len(lines):
                return None, None

            index += 1
            coords: List[List[float]] = []
            symbols: List[str] = []

            for line in lines[index:]:
                if re.match(r'\s*-{5,}\s*$', line):
                    break
                parts = line.strip().split()
                if len(parts) < 6 or not parts[0].isdigit():
                    continue
                try:
                    atomic_number = int(parts[1])
                    x = float(parts[3])
                    y = float(parts[4])
                    z = float(parts[5])
                except ValueError:
                    continue

                symbol = atomic_to_symbol.get(atomic_number, f'X{atomic_number}')
                coords.append([x, y, z])
                symbols.append(symbol)

            if not coords:
                return None, None

            return np.array(coords), symbols

        standard_index = _find_last_index("Standard orientation")
        if standard_index is not None:
            coords, symbols = _extract_from_index(standard_index)
            if coords is not None:
                logger.info("使用最后一个 Standard Orientation 块")
                return coords, symbols, None

        input_index = _find_last_index("Input orientation")
        if input_index is not None:
            coords, symbols = _extract_from_index(input_index)
            if coords is not None:
                logger.info("使用最后一个 Input Orientation 块")
                return coords, symbols, None

        return None, None, "未找到 Standard Orientation 或 Input Orientation 块"

    @staticmethod
    def _parse_orca_out(content: str) -> Tuple[Optional[np.ndarray], Optional[List[str]], Optional[str]]:
        """
        解析 ORCA .out 文件

        逻辑（多级回退）：
        1. CARTESIAN COORDINATES (ANGSTROEM) - 使用最后一次出现
        2. CARTESIAN COORDINATES (A.U.) + Bohr→Å转换 - 使用最后一次出现
        3. * xyz 块 - 使用最后一次出现

        参考旧脚本 awk 逻辑
        """
        lines = content.splitlines()

        def _scan_block(header_pattern: str) -> Optional[Tuple[np.ndarray, List[str]]]:
            header_indices = [
                index for index, line in enumerate(lines)
                if re.search(header_pattern, line, re.IGNORECASE)
            ]
            if not header_indices:
                return None

            start = header_indices[-1] + 1
            while start < len(lines) and not lines[start].strip():
                start += 1
            if start < len(lines) and re.match(r"\s*-{5,}\s*$", lines[start]):
                start += 1

            coords: List[List[float]] = []
            symbols: List[str] = []
            for line in lines[start:]:
                if not line.strip() or re.match(r"\s*-{5,}\s*$", line):
                    break
                match = re.match(r"\s*([A-Za-z]+)\s+([\-]?\d+\.\d+)\s+([\-]?\d+\.\d+)\s+([\-]?\d+\.\d+)\s*$", line)
                if not match:
                    continue
                symbol = match.group(1)
                coords.append([float(match.group(2)), float(match.group(3)), float(match.group(4))])
                symbols.append(symbol)

            if coords:
                return np.array(coords), symbols
            return None

        angstrom_block = _scan_block(r"CARTESIAN COORDINATES\s*\((?:ANGSTROEM|ANGSTROM)\)")
        if angstrom_block:
            coords, symbols = angstrom_block
            logger.info("从最后一个 CARTESIAN COORDINATES (ANGSTROEM) 提取坐标")
            return coords, symbols, None

        au_block = _scan_block(r"CARTESIAN COORDINATES\s*\(A\.U\.\)")
        if au_block:
            coords, symbols = au_block
            bohr_to_angstrom = 0.529177210903
            coords = coords * bohr_to_angstrom
            logger.info("从最后一个 CARTESIAN COORDINATES (A.U.) 提取坐标并转换单位")
            return coords, symbols, None

        xyz_pattern = r'\*\s*xyz\s+([\d\-]+)\s+([\d]+)((?:[^\*]|\n)+?)(?=\*|$)'
        xyz_matches = list(re.finditer(xyz_pattern, content, re.IGNORECASE | re.DOTALL))
        if xyz_matches:
            coords_text = xyz_matches[-1].group(3)
            coords = []
            symbols = []

            for line in coords_text.strip().split('\n'):
                parts = line.strip().split()
                if len(parts) >= 4:
                    try:
                        symbol = parts[0]
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])

                        coords.append([x, y, z])
                        symbols.append(symbol)
                    except (ValueError, IndexError):
                        continue

            if coords:
                logger.info("从最后一个 * xyz 块提取坐标")
                return np.array(coords), symbols, None

        return None, None, "无法从 ORCA .out 提取坐标"

    @staticmethod
    def _parse_coord_block(block_text: str) -> Tuple[Optional[np.ndarray], Optional[List[str]]]:
        """
        通用坐标块解析器

        格式：atom x y z (可能带电荷或其他列）

        Args:
            block_text: 坐标块文本

        Returns:
            (coords, symbols) 元组
        """
        coords = []
        symbols = []

        for line in block_text.strip().split('\n'):
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    symbol = parts[0]
                    x = float(parts[-3])
                    y = float(parts[-2])
                    z = float(parts[-1])

                    coords.append([x, y, z])
                    symbols.append(symbol)
                except (ValueError, IndexError):
                    continue

        return (np.array(coords), symbols) if coords else (None, None)

    @staticmethod
    def _parse_orca_input(content: str) -> Tuple[Optional[np.ndarray], Optional[List[str]], Optional[str]]:
        xyz_pattern = r'\*\s*xyz\s+[\d\-]+\s+[\d]+\s*\n((?:[^\*]|\n)+?)(?=\*|$)'
        matches = list(re.finditer(xyz_pattern, content, re.IGNORECASE | re.DOTALL))
        if not matches:
            return None, None, "无法从 ORCA .inp 提取坐标"

        coords_text = matches[-1].group(1)
        coords = []
        symbols = []
        for line in coords_text.strip().split('\n'):
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    symbol = parts[0]
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                except (ValueError, IndexError):
                    continue
                coords.append([x, y, z])
                symbols.append(symbol)

        if coords:
            return np.array(coords), symbols, None
        return None, None, "无法从 ORCA .inp 提取坐标"



