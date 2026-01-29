"""
Fragment Extractor - TS结构切分工具
======================================

从过渡态结构中提取独立的反应片段，用于畸变能计算

Author: QCcalc Team
Date: 2026-01-10
"""

import logging
from pathlib import Path
from typing import List, Tuple, Dict
import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import read_xyz, write_xyz
from rph_core.utils.qc_interface import GaussianInterface

logger = logging.getLogger(__name__)


class FragmentExtractor(LoggerMixin):
    """
    片段提取器 - 用于双片段畸变能计算

    功能:
    1. 从TS结构中切分出两个独立片段
    2. 生成片段XYZ文件
    3. 调用DFT进行单点能计算
    4. 支持BSSE校正（可选）
    """

    def __init__(self, config: dict):
        """
        初始化片段提取器

        Args:
            config: 配置字典，包含DFT方法、基组等
        """
        self.config = config
        self.method = config.get('method', 'B3LYP')
        self.basis = config.get('basis', 'def2-SVP')
        self.dispersion = config.get('dispersion', 'GD3BJ')
        self.nprocshared = config.get('nprocshared', 16)
        self.mem = config.get('mem', '32GB')
        self.solvent = config.get('solvent', 'acetone')

        # 初始化Gaussian接口
        self.gaussian = GaussianInterface(
            charge=0,
            multiplicity=1,
            nprocshared=self.nprocshared,
            mem=self.mem
        )

        self.logger.info(f"FragmentExtractor 初始化: {self.method}/{self.basis}")

    def extract_and_calculate(
        self,
        ts_xyz: Path,
        fragment_indices: Tuple[List[int], List[int]],
        output_dir: Path,
        old_checkpoint: Path = None,
        apply_bsse: bool = False
    ) -> Dict[str, float]:
        """
        提取片段并计算DFT单点能

        Args:
            ts_xyz: TS结构XYZ文件
            fragment_indices: (fragment_A_atoms, fragment_B_atoms)
            output_dir: 输出目录
            old_checkpoint: 可选的checkpoint文件（复用TS的轨道）
            apply_bsse: 是否应用BSSE校正

        Returns:
            {
                'e_fragment_a_ts': float,  # A片段在TS几何下的能量
                'e_fragment_b_ts': float,  # B片段在TS几何下的能量
                'e_fragment_a_relaxed': float,  # A片段松弛后的能量
                'e_fragment_b_relaxed': float,  # B片段松弛后的能量
            }
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("开始片段能量计算...")

        # 读取TS结构
        ts_coords, ts_symbols = read_xyz(ts_xyz)

        # 提取片段
        frag_A_indices, frag_B_indices = fragment_indices

        # 1. 生成片段结构（在TS几何下）
        frag_A_xyz = output_dir / "fragment_A_at_ts.xyz"
        frag_B_xyz = output_dir / "fragment_B_at_ts.xyz"

        self._extract_fragment_xyz(
            ts_coords, ts_symbols, frag_A_indices, frag_A_xyz,
            title="Fragment A (at TS geometry)"
        )
        self._extract_fragment_xyz(
            ts_coords, ts_symbols, frag_B_indices, frag_B_xyz,
            title="Fragment B (at TS geometry)"
        )

        # 2. 计算TS几何下的单点能（复用TS checkpoint）
        ts_chk = ts_xyz.with_suffix('.chk')
        if old_checkpoint and old_checkpoint.exists():
            ts_chk = old_checkpoint

        e_frag_A_ts = self._single_point(
            frag_A_xyz,
            output_dir / "fragment_A_ts_sp",
            old_checkpoint=ts_chk
        )

        e_frag_B_ts = self._single_point(
            frag_B_xyz,
            output_dir / "fragment_B_ts_sp",
            old_checkpoint=ts_chk
        )

        # 3. 几何优化片段（得到松弛能量）
        frag_A_relaxed_xyz = output_dir / "fragment_A_relaxed.xyz"
        frag_B_relaxed_xyz = output_dir / "fragment_B_relaxed.xyz"

        # 使用快速方法进行片段优化（GFN-xTB或低级别DFT）
        e_frag_A_relaxed = self._optimize_fragment(
            frag_A_xyz, frag_A_relaxed_xyz, output_dir / "fragment_A_opt"
        )

        e_frag_B_relaxed = self._optimize_fragment(
            frag_B_xyz, frag_B_relaxed_xyz, output_dir / "fragment_B_opt"
        )

        # 4. 可选：BSSE校正
        if apply_bsse:
            self.logger.info("应用BSSE校正...")
            # TODO: 实现counterpoise校正
            # 这需要计算ghost原子体系

        results = {
            'e_fragment_a_ts': e_frag_A_ts,
            'e_fragment_b_ts': e_frag_B_ts,
            'e_fragment_a_relaxed': e_frag_A_relaxed,
            'e_fragment_b_relaxed': e_frag_B_relaxed
        }

        self.logger.info("片段能量计算完成:")
        self.logger.info(f"  E(A_TS) = {e_frag_A_ts:.6f} Hartree")
        self.logger.info(f"  E(B_TS) = {e_frag_B_ts:.6f} Hartree")
        self.logger.info(f"  E(A_relaxed) = {e_frag_A_relaxed:.6f} Hartree")
        self.logger.info(f"  E(B_relaxed) = {e_frag_B_relaxed:.6f} Hartree")

        return results

    def _extract_fragment_xyz(
        self,
        coords: np.ndarray,
        symbols: List[str],
        indices: List[int],
        output_xyz: Path,
        title: str = "Fragment"
    ):
        """
        从完整结构中提取片段XYZ

        Args:
            coords: 完整坐标 (N, 3)
            symbols: 原子符号列表
            indices: 片段原子索引列表
            output_xyz: 输出XYZ文件
            title: 标题
        """
        # 提取片段坐标和符号
        frag_coords = coords[indices]
        frag_symbols = [symbols[i] for i in indices]

        # 保存XYZ
        write_xyz(output_xyz, frag_coords, frag_symbols, title=title)

        self.logger.info(f"✓ 提取片段: {output_xyz} ({len(indices)} 原子)")

    def _single_point(
        self,
        xyz_file: Path,
        output_dir: Path,
        old_checkpoint: Path = None
    ) -> float:
        """
        DFT单点能计算

        Args:
            xyz_file: 输入XYZ文件
            output_dir: 输出目录
            old_checkpoint: 复用的checkpoint

        Returns:
            能量 (Hartree)
        """
        # 构建route card（单点能）
        route = (
            f"# {self.method}/{self.basis} "
            f"EmpiricalDispersion={self.dispersion} "
            f"SCRF=(SMD,solvent={self.solvent}) SP"
        )

        try:
            # 调用Gaussian单点计算
            result = self.gaussian.optimize(
                xyz_file=xyz_file,
                output_dir=output_dir,
                route=route,
                old_checkpoint=old_checkpoint
            )

            if not result.converged:
                self.logger.warning(f"单点计算未收敛: {xyz_file.name}")

            return result.energy

        except Exception as e:
            self.logger.error(f"单点计算失败: {e}")
            raise RuntimeError(f"片段单点计算失败: {e}")

    def _optimize_fragment(
        self,
        frag_xyz: Path,
        output_xyz: Path,
        output_dir: Path
    ) -> float:
        """
        片段几何优化（使用快速方法）

        策略:
        1. 首先尝试GFN2-xTB优化
        2. 可选：低级别DFT精修（B3LYP/def2-SVP）

        Args:
            frag_xyz: 片段XYZ文件
            output_xyz: 输出优化后的XYZ
            output_dir: 输出目录

        Returns:
            优化后的能量 (Hartree)
        """
        from rph_core.utils.qc_interface import XTBInterface

        # 使用XTB进行快速优化
        xtb = XTBInterface(
            gfn_level=2,
            solvent=self.solvent,
            nproc=8
        )

        self.logger.info(f"优化片段: {frag_xyz.name} (GFN2-xTB)")
        result = xtb.optimize(frag_xyz, output_dir)

        if result.converged:
            # 保存优化后的结构
            write_xyz(
                output_xyz,
                result.coordinates,
                self._read_symbols(frag_xyz),
                title=f"{frag_xyz.stem} relaxed",
                energy=result.energy
            )
            return result.energy
        else:
            self.logger.warning(f"XTB优化未收敛，使用原始能量")
            # 读取原始能量
            coords, energy = self._read_xyz_with_energy(frag_xyz)
            return energy

    def _read_symbols(self, xyz_file: Path) -> List[str]:
        """从XYZ文件读取原子符号"""
        coords, symbols = read_xyz(xyz_file)
        return symbols

    def _read_xyz_with_energy(self, xyz_file: Path) -> Tuple[np.ndarray, float]:
        """从XYZ文件读取坐标和能量"""
        with open(xyz_file, 'r') as f:
            lines = f.readlines()

        n_atoms = int(lines[0].strip())
        title_line = lines[1].strip()

        # 尝试从标题提取能量
        energy = 0.0
        if "energy:" in title_line or "E =" in title_line:
            import re
            match = re.search(r'[\-]?\d+\.\d+', title_line)
            if match:
                energy = float(match.group())

        # 读取坐标
        coords = []
        for line in lines[2:2+n_atoms]:
            parts = line.strip().split()
            if len(parts) >= 4:
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

        return np.array(coords), energy
