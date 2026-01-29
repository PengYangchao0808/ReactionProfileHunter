"""
Gaussian IRC Driver - Intrinsic Reaction Coordinate Calculation
================================================================

IRC计算驱动器 - 验证TS连接正确的reactant和product

Author: QCcalc Team
Date: 2026-01-10
"""

import subprocess
import re
import logging
from pathlib import Path
from typing import Optional, Tuple
from dataclasses import dataclass

import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.qc_interface import GaussianInterface
from rph_core.utils.file_io import read_xyz
from rph_core.utils.resource_utils import get_project_root


logger = logging.getLogger(__name__)


@dataclass
class IRCResult:
    """IRC计算结果"""
    converged: bool
    forward_endpoint: np.ndarray
    reverse_endpoint: np.ndarray
    forward_path: Optional[np.ndarray] = None  # Shape: (n_points, n_atoms, 3)
    reverse_path: Optional[np.ndarray] = None
    forward_n_points: int = 0
    reverse_n_points: int = 0
    output_file: Optional[Path] = None
    ts_energy: float = 0.0
    reactant_energy: float = 0.0
    product_energy: float = 0.0


class IRCDriver(LoggerMixin):
    """
    IRC 计算驱动器 (Step 3 - TS验证)

    功能:
    - 生成 Gaussian IRC 输入文件
    - 提交 IRC 计算 (正向和反向)
    - 解析 IRC 路径和端点
    - 验证 TS 连接正确的底物和产物
    """

    def __init__(self, method: str = "B3LYP", basis: str = "def2-SVP",
                 dispersion: str = "GD3BJ", nprocshared: int = 16,
                 mem: str = "32GB", max_points: int = 50,
                 step_size: int = 10, config: dict = {}):
        """
        初始化 IRC 驱动器

        Args:
            method: DFT 方法
            basis: 基组
            dispersion: 色散校正
            nprocshared: 共享处理器数
            mem: 内存大小
            max_points: IRC最大点数
            step_size: IRC步长 (0.001 amu^1/2 bohr)
            config: 配置字典
        """
        self.method = method
        self.basis = basis
        self.dispersion = dispersion
        self.nprocshared = nprocshared
        self.mem = mem
        self.max_points = max_points
        self.step_size = step_size
        self.config = config

        # 构建 Gaussian route card
        # IRC=(CalcFC, MaxPoints=N, StepSize=S) 会同时计算正向和反向
        self.route = (
            f"# {method}/{basis} EmpiricalDispersion={dispersion} "
            f"IRC=(CalcFC, MaxPoints={max_points}, StepSize={step_size})"
        )

        self.gaussian = GaussianInterface(
            charge=0,
            multiplicity=1,
            nprocshared=nprocshared,
            mem=mem,
            config=self.config
        )

        self.logger.info(f"IRCDriver 初始化: {self.route}")

    def run_irc(self, ts_xyz: Path, output_dir: Path,
                direction: str = "both") -> IRCResult:
        """
        执行 IRC 计算

        Args:
            ts_xyz: TS XYZ 文件
            output_dir: 输出目录
            direction: IRC 方向 ("forward", "reverse", "both")

        Returns:
            IRCResult 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"开始 IRC 计算: {ts_xyz.name}")
        self.logger.info(f"方向: {direction}")

        # 1. 生成 Gaussian IRC 输入文件
        gjf_file = output_dir / f"{ts_xyz.stem}_irc.gjf"
        self._write_irc_input(ts_xyz, gjf_file)

        # 2. 提交 Gaussian IRC 计算
        self._submit_irc_job(gjf_file, output_dir)

        # 3. 解析 IRC 结果
        log_file = output_dir / f"{ts_xyz.stem}_irc.log"

        if not log_file.exists():
            raise RuntimeError(f"IRC 输出文件不存在: {log_file}")

        result = self._parse_irc_result(log_file, ts_xyz)

        self.logger.info(f"✓ IRC 计算完成")
        self.logger.info(f"  正向点数: {result.forward_n_points}")
        self.logger.info(f"  反向点数: {result.reverse_n_points}")
        self.logger.info(f"  收敛: {'是' if result.converged else '否'}")

        return result

    def _write_irc_input(self, ts_xyz: Path, gjf_file: Path):
        """
        生成 Gaussian IRC 输入文件

        Args:
            ts_xyz: TS XYZ 文件
            gjf_file: 输出 GJF 文件路径
        """
        # 读取 XYZ 文件
        coords, symbols = read_xyz(ts_xyz)

        # 构建输入内容
        content = f"""%chk={gjf_file.stem}.chk
%mem={self.mem}
%nproc={self.nprocshared}
{self.route}

IRC from TS

0 1
"""

        # 添加坐标
        for symbol, (x, y, z) in zip(symbols, coords):
            content += f"{symbol} {x:12.6f} {y:12.6f} {z:12.6f}\n"

        content += "\n\n"  # Gaussian 需要空行

        # 写入文件
        gjf_file.write_text(content, encoding='utf-8')
        self.logger.info(f"✓ IRC 输入文件已生成: {gjf_file}")

    def _submit_irc_job(self, gjf_file: Path, output_dir: Path):
        """
        提交 Gaussian IRC 计算
        """
        self.logger.info("提交 Gaussian IRC 计算...")
        log_file = gjf_file.with_suffix('.log')

        exe_config = self.config.get('executables', {}).get('gaussian', {})
        use_wrapper = exe_config.get('use_wrapper', True)

        if use_wrapper:
            wrapper_path = exe_config.get('wrapper_path', './scripts/run_g16_worker.sh')
            if not Path(wrapper_path).is_absolute():
                wrapper_path = (get_project_root() / wrapper_path).resolve()
            g16_cmd = str(wrapper_path)
        else:
            g16_cmd = 'g16'

        try:
            cmd = f'"{g16_cmd}" {gjf_file.name} {log_file.name}'
            subprocess.run(cmd, shell=True, cwd=str(output_dir), check=True)
            self.logger.info(f"✓ Gaussian IRC 完成: {log_file}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Gaussian IRC 执行失败: {e}")
            raise RuntimeError(f"Gaussian IRC 失败")
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"Gaussian IRC 计算超时: {gjf_file.name}")
        except FileNotFoundError:
            raise RuntimeError("未找到 Gaussian (g16)，请确保已安装并在 PATH 中")

    def _parse_irc_result(self, log_file: Path, ts_xyz: Path) -> IRCResult:
        """
        解析 Gaussian IRC 输出文件

        Args:
            log_file: Gaussian .log 文件
            ts_xyz: 原始 TS XYZ 文件

        Returns:
            IRCResult 对象
        """
        with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        # 1. 检查收敛
        converged = self._check_irc_convergence(content)

        # 2. 提取能量
        energies = self._extract_irc_energies(content)
        ts_energy = energies.get('ts', 0.0)
        reactant_energy = energies.get('reactant', 0.0)
        product_energy = energies.get('product', 0.0)

        # 3. 提取 IRC 路径
        forward_path, reverse_path = self._extract_irc_paths(content)

        # 4. 提取端点坐标
        if forward_path is not None and len(forward_path) > 0:
            forward_endpoint = forward_path[-1]  # 最后一个点
            forward_n_points = len(forward_path)
        else:
            # 读取原始 TS 坐标作为后备
            coords_ts, _ = read_xyz(ts_xyz)
            forward_endpoint = coords_ts
            forward_n_points = 0

        if reverse_path is not None and len(reverse_path) > 0:
            reverse_endpoint = reverse_path[-1]
            reverse_n_points = len(reverse_path)
        else:
            coords_ts, _ = read_xyz(ts_xyz)
            reverse_endpoint = coords_ts
            reverse_n_points = 0

        return IRCResult(
            converged=converged,
            forward_endpoint=forward_endpoint,
            reverse_endpoint=reverse_endpoint,
            forward_path=forward_path,
            reverse_path=reverse_path,
            forward_n_points=forward_n_points,
            reverse_n_points=reverse_n_points,
            output_file=log_file,
            ts_energy=ts_energy,
            reactant_energy=reactant_energy,
            product_energy=product_energy
        )

    def _check_irc_convergence(self, content: str) -> bool:
        """
        检查 IRC 是否收敛

        Args:
            content: Gaussian 输出文件内容

        Returns:
            是否收敛
        """
        # 查找收敛标记
        # Gaussian IRC 会在每个方向完成后报告 "IRC converged"
        if "IRC converged" in content:
            return True

        # 检查是否正常终止
        if "Normal termination" in content:
            return True

        return False

    def _extract_irc_energies(self, content: str) -> dict:
        """
        提取 IRC 各点的能量

        Args:
            content: Gaussian 输出文件内容

        Returns:
            能量字典 {'ts': float, 'reactant': float, 'product': float}
        """
        energies = {}

        # 提取 TS 能量 (第一个点)
        ts_match = re.search(
            r'IRC.*Point\s+1.*?\n.*?SCF Done:\s+E\([^)]+\)\s*=\s*([\-\\d.]+)',
            content,
            re.DOTALL
        )
        if ts_match:
            energies['ts'] = float(ts_match.group(1))

        # 提取 reactant 能量 (反向终点)
        # 查找最后一个 "Reactant" 或 "Reverse" 方向的能量
        reverse_matches = re.findall(
            r'IRC.*Reverse.*?\n.*?SCF Done:\s+E\([^)]+\)\s*=\s*([\-\\d.]+)',
            content,
            re.DOTALL
        )
        if reverse_matches:
            energies['reactant'] = float(reverse_matches[-1])

        # 提取 product 能量 (正向终点)
        forward_matches = re.findall(
            r'IRC.*Forward.*?\n.*?SCF Done:\s+E\([^)]+\)\s*=\s*([\-\\d.]+)',
            content,
            re.DOTALL
        )
        if forward_matches:
            energies['product'] = float(forward_matches[-1])

        return energies

    def _extract_irc_paths(self, content: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        提取 IRC 路径坐标

        Args:
            content: Gaussian 输出文件内容

        Returns:
            (forward_path, reverse_path) 元组
            每个 path 的形状: (n_points, n_atoms, 3)
        """
        # 查找所有 IRC 点的坐标
        # Gaussian 会在每个 IRC 点输出 "Input orientation"

        # 提取所有坐标块
        coord_blocks = re.findall(
            r'Input orientation:.*?\n'
            r'-*\s*\n'
            r'((?:\s*\d+\s+\d+\s+\d+\s+[\-\\d\.]+\s+[\-\\d\.]+\s+[\-\\d\.]+\s+[\-\\d\.]+\n)+)',
            content,
            re.DOTALL
        )

        if not coord_blocks:
            return None, None

        # 解析所有坐标
        all_coords = []
        for block in coord_blocks:
            coords = []
            for line in block.strip().split('\n'):
                parts = line.split()
                if len(parts) >= 6:
                    x = float(parts[3])
                    y = float(parts[4])
                    z = float(parts[5])
                    coords.append([x, y, z])
            if coords:
                all_coords.append(np.array(coords))

        if not all_coords:
            return None, None

        all_coords = np.array(all_coords)  # (n_total_points, n_atoms, 3)

        # 分离正向和反向路径
        # Gaussian IRC 输出: 先正向 (forward)，再反向 (reverse)
        # 通常在中间会有 "IRC REVERSE" 这样的标记

        reverse_marker = content.find("IRC REVERSE")
        if reverse_marker == -1:
            reverse_marker = content.find("Reverse direction")

        if reverse_marker != -1:
            # 估算反向路径开始的点索引
            # 计算在 reverse_marker 之前有多少个坐标块
            content_before_reverse = content[:reverse_marker]
            n_forward_blocks = content_before_reverse.count("Input orientation:")

            if n_forward_blocks > 0 and n_forward_blocks < len(all_coords):
                forward_path = all_coords[:n_forward_blocks]
                reverse_path = all_coords[n_forward_blocks:]
            else:
                # 无法分离，返回全部作为正向
                forward_path = all_coords
                reverse_path = None
        else:
            # 没有反向标记，假设全部是正向
            forward_path = all_coords
            reverse_path = None

        return forward_path, reverse_path

    def calculate_reaction_profile(
        self,
        irc_result: IRCResult
    ) -> dict:
        """
        计算反应能垒和反应能

        Args:
            irc_result: IRC 计算结果

        Returns:
            反应剖面参数字典
        """
        # 转换为 kcal/mol (1 Hartree = 627.509 kcal/mol)
        hartree_to_kcal = 627.509

        delta_E_forward = (irc_result.ts_energy - irc_result.reactant_energy) * hartree_to_kcal
        delta_E_reverse = (irc_result.ts_energy - irc_result.product_energy) * hartree_to_kcal
        delta_E_reaction = (irc_result.product_energy - irc_result.reactant_energy) * hartree_to_kcal

        profile = {
            'delta_E_forward': delta_E_forward,  # 正向能垒 (kcal/mol)
            'delta_E_reverse': delta_E_reverse,  # 反向能垒 (kcal/mol)
            'delta_E_reaction': delta_E_reaction,  # 反应能 (kcal/mol)
            'ts_energy_hartree': irc_result.ts_energy,
            'reactant_energy_hartree': irc_result.reactant_energy,
            'product_energy_hartree': irc_result.product_energy,
        }

        self.logger.info(f"反应剖面:")
        self.logger.info(f"  ΔG‡ (正向) = {delta_E_forward:.2f} kcal/mol")
        self.logger.info(f"  ΔG‡ (反向) = {delta_E_reverse:.2f} kcal/mol")
        self.logger.info(f"  ΔGᵣ (反应) = {delta_E_reaction:.2f} kcal/mol")

        return profile
