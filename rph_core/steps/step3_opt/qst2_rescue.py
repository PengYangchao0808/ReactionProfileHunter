"""
Gaussian QST2 Rescue Driver
============================

QST2救援策略驱动器

Author: QCcalc Team
Date: 2026-01-09
"""

import subprocess
import re
import logging
from pathlib import Path
from typing import Tuple, Any, Optional
from dataclasses import dataclass
import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.qc_interface import GaussianInterface
from rph_core.utils.qc_runner import run_with_timeout, QCTimeoutError
from rph_core.utils.file_io import read_xyz


logger = logging.getLogger(__name__)


@dataclass
class QST2Result:
    """QST2优化结果"""
    converged: bool
    coordinates: Optional[np.ndarray]
    energy: float
    imaginary_freq: float
    frequencies: Optional[np.ndarray]
    output_file: Path
    method_used: str = "QST2"


class QST2RescueDriver(LoggerMixin):
    """
    QST2 救援驱动器 (Step 3 - 救援策略)

    功能:
    - 使用 reactant 和 product 两个端点进行 QST2 优化
    - 在 Berny TS 失败时启用
    """

    def __init__(self, method: str = "B3LYP", basis: str = "def2-SVP",
                 dispersion: str = "GD3BJ", nprocshared: int = 16,
                 mem: str = "32GB", config: Optional[dict[str, Any]] = None):
        """
        初始化 QST2 驱动器

        Args:
            method: DFT 方法
            basis: 基组
            dispersion: 色散校正
            nprocshared: 共享处理器数
            mem: 内存大小
            config: 配置字典 (用于 wrapper script 支持)
        """
        self.method = method
        self.basis = basis
        self.dispersion = dispersion
        self.nprocshared = nprocshared
        self.mem = mem

        self.config = config or {}
        timeout_cfg = self.config.get('optimization_control', {}).get('timeout', {})
        timeout_default = timeout_cfg.get('default_seconds', 21600)
        step3_cfg = self.config.get('step3', {})
        timeout_value = step3_cfg.get('qst2_timeout_seconds', timeout_default)
        try:
            self.timeout_seconds = max(1, int(timeout_value))
        except (TypeError, ValueError):
            self.timeout_seconds = 21600

        exe_config = self.config.get('executables', {}).get('gaussian', {})
        use_wrapper = exe_config.get('use_wrapper', True)

        if use_wrapper:
            wrapper_path = exe_config.get('wrapper_path', './scripts/run_g16_worker.sh')
            if not Path(wrapper_path).is_absolute():
                wrapper_path = (Path.cwd() / wrapper_path).resolve()
            self.gaussian_cmd = str(wrapper_path)
        else:
            self.gaussian_cmd = 'g16'

        # QST2 route card
        self.route = (
            f"# {method}/{basis} EmpiricalDispersion={dispersion} "
            f"Opt=(QST2, CalcFC) Freq"
        )

        self.logger.info(f"QST2RescueDriver 初始化: {self.route}")
        
        self.gaussian = GaussianInterface(
            charge=0, multiplicity=1, # 默认
            nprocshared=nprocshared,
            mem=mem
        )

    def optimize(self, reactant: Path, product: Path,
                 output_dir: Path) -> QST2Result:
        """
        执行 QST2 优化

        Args:
            reactant: 底物复合物 XYZ 文件
            product: 产物 XYZ 文件
            output_dir: 输出目录

        Returns:
            QST2Result 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("开始 QST2 优化（救援策略）")
        
        gjf_file = output_dir / "qst2_input.gjf"
        log_file = output_dir / "qst2_input.log"
        
        # 1. 生成 QST2 输入文件
        self._write_qst2_input(reactant, product, gjf_file)

        # 2. 提交计算
        try:
            cmd = [self.gaussian_cmd, gjf_file.name, log_file.name]
            run_with_timeout(cmd=cmd, timeout=self.timeout_seconds, cwd=output_dir)
        except (subprocess.CalledProcessError, FileNotFoundError, OSError) as e:
            self.logger.error(f"QST2 运行失败: {e}")
            raise RuntimeError(f"QST2 运行失败: {e}")
        except QCTimeoutError:
            raise RuntimeError(f"QST2 计算超时: {gjf_file.name}")
            
        # 3. 解析结果
        from rph_core.utils.gaussian_log_parser import GaussianLogParser
        res = GaussianLogParser.parse_log(log_file)
        
        imag_freq = 0.0
        if res is None:
            raise RuntimeError(f"QST2 日志解析失败: {log_file}")

        if res.frequencies is not None:
            imag_freqs = [f for f in res.frequencies if f < 0]
            if imag_freqs:
                imag_freq = imag_freqs[0]

        return QST2Result(
            converged=res.converged,
            coordinates=res.coordinates,
            energy=float(res.energy) if res.energy is not None else 0.0,
            imaginary_freq=imag_freq,
            frequencies=res.frequencies,
            output_file=Path(res.output_file) if res.output_file is not None else log_file,
            method_used="QST2"
        )

    def _write_qst2_input(self, reactant: Path, product: Path, gjf_file: Path):
        """
        生成 QST2 输入文件

        Args:
            reactant: 底物 XYZ 文件
            product: 产物 XYZ 文件
            gjf_file: 输出的 GJF 文件
        """
        # 读取两个结构
        coords_r, symbols_r = read_xyz(reactant)
        coords_p, symbols_p = read_xyz(product)

        # 验证原子数相同
        n_atoms_r = len(symbols_r)
        n_atoms_p = len(symbols_p)

        if n_atoms_r != n_atoms_p:
            raise ValueError(
                f"Reactant ({n_atoms_r} atoms) 和 Product ({n_atoms_p} atoms) "
                f"原子数不同，无法进行 QST2"
            )

        # 构建 GJF 文件
        gjf_content = f"%chk=qst2.chk\n"
        gjf_content += f"%nprocshared={self.nprocshared}\n"
        gjf_content += f"%mem={self.mem}\n"
        gjf_content += f"{self.route}\n\n"
        gjf_content += f"QST2 TS Search\n\n"

        # 第一个结构: Reactant
        gjf_content += f"0 1\n"
        for i, (symbol, coord) in enumerate(zip(symbols_r, coords_r)):
            gjf_content += f"{symbol:2s} {coord[0]:15.10f} {coord[1]:15.10f} {coord[2]:15.10f}\n"

        gjf_content += "\n"

        # 第二个结构: Product
        gjf_content += f"0 1\n"
        for i, (symbol, coord) in enumerate(zip(symbols_p, coords_p)):
            gjf_content += f"{symbol:2s} {coord[0]:15.10f} {coord[1]:15.10f} {coord[2]:15.10f}\n"

        gjf_content += "\n"

        # 写入文件
        with open(gjf_file, 'w') as f:
            f.write(gjf_content)

        self.logger.info(f"✓ QST2 输入文件已生成: {gjf_file}")

    def _submit_gaussian_job(self, gjf_file: Path, output_dir: Path):
        """提交 Gaussian 计算"""
        self.logger.info("提交 QST2 Gaussian 计算...")

        cmd = ["g16", str(gjf_file)]

        try:
            result = subprocess.run(
                cmd,
                cwd=str(output_dir),
                capture_output=True,
                text=True,
                timeout=3600
            )

            if result.returncode != 0:
                raise RuntimeError(f"QST2 Gaussian 失败: {result.stderr}")

            self.logger.info("✓ QST2 Gaussian 计算完成")

        except subprocess.TimeoutExpired:
            raise RuntimeError(f"QST2 Gaussian 计算超时")
        except FileNotFoundError:
            raise RuntimeError("未找到 Gaussian (g16)")

    def _parse_qst2_result(self, log_file: Path) -> QST2Result:
        """解析 QST2 结果（与 Berny 类似）"""
        import numpy as np

        with open(log_file, 'r') as f:
            content = f.read()

        # 1. 提取能量
        energy_match = re.search(r'SCF Done:\s+E\([^)]+\)\s*=\s*([\-\d.]+)', content)
        if not energy_match:
            raise RuntimeError("无法提取能量")
        energy = float(energy_match.group(1))

        # 2. 提取频率
        frequencies = self._extract_frequencies(content)
        imaginary_freqs = [f for f in frequencies if f < 0]

        if len(imaginary_freqs) != 1:
            raise RuntimeError(f"期望1个虚频,实际{len(imaginary_freqs)}个")

        # 3. 提取坐标
        coords = self._extract_optimized_coordinates(content)

        return QST2Result(
            converged=True,
            coordinates=coords,
            energy=energy,
            imaginary_freq=imaginary_freqs[0],
            frequencies=frequencies,
            output_file=log_file,
            method_used="QST2"
        )

    def _extract_frequencies(self, content: str) -> 'np.ndarray':
        """提取频率（与Berny相同）"""
        import numpy as np

        freq_pattern = r'Frequencies --\s*([\d\-.]+)\s*([\d\-.]+)\s*([\d\-.]+)'
        frequencies = []

        for match in re.finditer(freq_pattern, content):
            for i in range(1, 4):
                try:
                    freq = float(match.group(i))
                    frequencies.append(freq)
                except (ValueError, IndexError):
                    continue

        return np.array(frequencies) if frequencies else np.array([])

    def _extract_optimized_coordinates(self, content: str) -> 'np.ndarray':
        """提取优化后的坐标（与Berny相同）"""
        import numpy as np

        pattern = r'Standard orientation:\s*-*\s*((?:\s*\d+\s+\d+\s+\d+\s+[\-\d\.]+\s+[\-\d\.]+\s+[\-\d\.]+\s+[\-\d\.]+\n)+)'

        matches = list(re.finditer(pattern, content, re.MULTILINE))

        if matches:
            last_match = matches[-1]
            coord_text = last_match.group(1)

            coords = []
            for line in coord_text.strip().split('\n'):
                parts = line.split()
                if len(parts) >= 6:
                    x = float(parts[3])
                    y = float(parts[4])
                    z = float(parts[5])
                    coords.append([x, y, z])

            return np.array(coords)

        # 如果无法提取，返回空数组
        return np.array([])
