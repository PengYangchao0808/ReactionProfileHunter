"""
ORCA Interface Module
=====================

ORCA 量子化学软件接口 - 用于高精度单点能计算

Author: QC Descriptors Team
Date: 2026-01-10
Session: #1 - ORCAInterface._generate_input()
Session: #2 - ORCAInterface._parse_output()
Session: #3 - ORCAInterface._find_orca_binary() + _run_orca()
Session: #4 - ORCAInterface.single_point()
"""

from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Union, TYPE_CHECKING
import re
import subprocess
import shutil
import os
import logging
import numpy as np

if TYPE_CHECKING:
    from rph_core.utils.optimization_config import OptimizationConfig

from rph_core.utils.resource_utils import (
    find_executable,
    setup_ld_library_path,
    mem_to_mb,
    calc_orca_maxcore,
    resolve_executable_config
)
from rph_core.utils.geometry_tools import CoordinateExtractor
from rph_core.utils.qc_interface import QCResult

logger = logging.getLogger(__name__)


class ORCAInterface:
    """ORCA 接口 - 高精度单点能计算"""

    def __init__(
        self,
        method: str = "M062X",
        basis: str = "def2-TZVPP",
        aux_basis: str = "def2/J",
        nprocs: int = 16,
        maxcore: Optional[int] = None,
        solvent: str = "acetone",
        orca_binary_path: Optional[str] = None,
        config: Optional[dict] = None
    ):
        """
        初始化 ORCA 接口

        Args:
            method: DFT 方法 (如 M062X, B3LYP, PWPB95)
            basis: 基组 (如 def2-TZVPP, def2-SVP)
            aux_basis: 辅助基组 (如 def2/J)
            nprocs: 并行进程数
            maxcore: 每核最大内存 (MB，None 则从配置计算）
            solvent: 溶剂名称 (用于 SMD 模型)
            orca_binary_path: ORCA 可执行文件路径 (可选)
            config: 配置字典 (可选，用于派生路径和内存）
        """
        self.method = method
        self.basis = basis
        self.aux_basis = aux_basis
        self.nprocs = nprocs
        self.solvent = solvent

        # 处理 maxcore：如果未提供，从配置派生
        if maxcore is None and config:
            res_cfg = config.get('resources', {})
            mem = res_cfg.get('mem', '32GB')
            safety_factor = res_cfg.get('orca_maxcore_safety', 0.8)
            self.maxcore = calc_orca_maxcore(mem, nprocs, safety_factor)
        else:
            self.maxcore = maxcore if maxcore is not None else 4000

        # 查找 ORCA 二进制文件（集成新的配置系统）
        self.orca_binary = self._find_orca_binary(orca_binary_path, config)

        # 设置日志
        self.logger = logging.getLogger(f"{__name__}.{method}/{basis}")

    def _is_double_hybrid(self) -> bool:
        """
        检查当前方法是否为双杂化泛函

        双杂化泛函需要额外的 /C 辅助基组 (如 def2-TZVPP/C)

        Returns:
            是否为双杂化泛函
        """
        double_hybrid_functionals = [
            "PWPB95", "DSD-PBEP86", "DSD-PBEP95",
            "B2PLYP", "B2GPPLYP", "DSD-BLYP"
        ]

        return self.method.upper() in [fh.upper() for fh in double_hybrid_functionals]

    def _generate_input(self, xyz_file: Path, output_dir: Path) -> Path:
        """
        生成 ORCA 输入文件

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录

        Returns:
            生成的 .inp 文件路径
        """
        # 构建路由行
        route = f"! {self.method} {self.basis} {self.aux_basis} RIJCOSX tightSCF"
        route += " noautostart miniprint nopop"

        # 双杂化泛函需要额外的 /C 辅助基组
        if self._is_double_hybrid():
            # 从主基组名称推导 /C 基组
            c_basis = self.basis + "/C"
            route += f" {c_basis}"

        # 构建溶剂块
        cpcm_block = ""
        if self.solvent and self.solvent.upper() != "NONE":
            cpcm_block = f"""
%cpcm
   smd true
   SMDsolvent "{self.solvent}"
end
"""

        # 如果 xyz_file 是 .out 文件，尝试从中提取电荷和自旋
        charge, spin = 0, 1
        if xyz_file.suffix == '.out':
            try:
                charge, spin = CoordinateExtractor.get_charge_spin_from_orca_out(xyz_file)
                self.logger.debug(f"从 {xyz_file.name} 提取电荷/自旋: {charge}/{spin}")
            except Exception as e:
                self.logger.warning(f"无法从 {xyz_file.name} 提取电荷/自旋: {e}")

        # 组装完整输入文件内容 - 直接嵌入xyz坐标
        import shutil
        xyz_copy = output_dir / xyz_file.name
        shutil.copy(xyz_file, xyz_copy)

        # 读取xyz文件内容，跳过前两行（原子数和标题）
        xyz_lines = xyz_copy.read_text().split('\n')
        if len(xyz_lines) > 2:
            xyz_content = '\n'.join(xyz_lines[2:])
        else:
            xyz_content = xyz_copy.read_text()

        inp_content = f"""{route}
%maxcore {self.maxcore}
%pal nprocs {self.nprocs} end
{cpcm_block}
 * xyz {charge} {spin}
{xyz_content}
 *
"""

        # 写入文件
        inp_file = output_dir / f"{xyz_file.stem}.inp"
        inp_file.write_text(inp_content)

        return inp_file

    def _parse_output(self, out_file: Path) -> QCResult:
        """
        解析 ORCA 输出文件

        Args:
            out_file: ORCA 输出文件路径

        Returns:
            QCResult 对象
        """
        # 读取输出文件内容
        try:
            content = out_file.read_text()
        except Exception as e:
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=f"无法读取输出文件: {e}"
            )

        # 检查是否正常终止
        if "ORCA TERMINATED NORMALLY" not in content:
            return QCResult(
                energy=0.0,
                converged=False,
                error_message="ORCA 未正常终止"
            )

        # 提取最终单点能
        energy_match = re.search(r"FINAL SINGLE POINT ENERGY\s+([\-\d\.]+)", content)
        if not energy_match:
            return QCResult(
                energy=0.0,
                converged=False,
                error_message="无法找到能量信息"
            )

        # 解析能量
        try:
            energy = float(energy_match.group(1))
        except ValueError:
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=f"能量格式错误: {energy_match.group(1)}"
            )

        # 成功解析
        return QCResult(
            energy=energy,
            converged=True,
            output_file=out_file,
            error_message=None
        )

    def _find_orca_binary(self, provided_path: Optional[str] = None, config: Optional[dict] = None) -> Optional[Path]:
        """
        查找 ORCA 可执行文件（集成新的配置系统）

        查找顺序:
        1. 提供的路径
        2. 配置文件中的路径
        3. 环境变量 $ORCA_PATH, $ORCA_BIN
        4. PATH 中的 'orca' 命令

        Args:
            provided_path: 用户提供的 ORCA 路径
            config: 配置字典

        Returns:
            ORCA 可执行文件的 Path 对象，如果找不到返回 None
        """
        if provided_path:
            path = Path(provided_path)
            if path.exists() and path.is_file():
                logger.info(f"使用提供的 ORCA 路径: {path}")
                return path
            else:
                logger.warning(f"提供的 ORCA 路径不存在: {provided_path}")

        if config:
            exe_config = resolve_executable_config(config, 'orca', env_vars=['ORCA_PATH', 'ORCA_BIN'])
            logger.info(f"从配置获取 ORCA: {exe_config}")
            if exe_config.get('found'):
                return exe_config['path']

        orca_cmd = shutil.which('orca')
        if orca_cmd:
            logger.info(f"从系统 PATH 找到 ORCA: {orca_cmd}")
            return Path(orca_cmd)

        logger.error("未找到 ORCA 可执行文件")
        return None

    def _run_orca(
        self,
        inp_file: Path,
        output_dir: Path,
        timeout: Optional[int] = 3600
    ) -> Path:
        """
        运行 ORCA 计算

        Args:
            inp_file: ORCA 输入文件
            output_dir: 输出目录
            timeout: 超时时间（秒），默认 1 小时

        Returns:
            ORCA 输出文件路径

        Raises:
            RuntimeError: ORCA 二进制文件未找到
            RuntimeError: ORCA 运行失败
            TimeoutError: 计算超时
        """
        # 检查 ORCA 是否可用
        if self.orca_binary is None:
            raise RuntimeError(
                "ORCA 二进制文件未找到。请通过以下方式之一指定:\n"
                "  1. 设置环境变量 $ORCA_PATH\n"
                "  2. 确保 'orca' 在 PATH 中\n"
                "  3. 在初始化时提供 orca_binary_path 参数"
            )

        out_file = inp_file.with_suffix('.out')

        # ORCA MPI 路径解析不能处理空格，使用临时目录
        if ' ' in str(output_dir.resolve()) or ' ' in str(inp_file.resolve()):
            import tempfile
            import shutil

            self.logger.warning("检测到路径包含空格，使用临时目录运行 ORCA")

            temp_dir = Path(tempfile.mkdtemp(prefix='orca_run_'))
            try:
                temp_inp_file = temp_dir / "orca_input.inp"
                temp_out_file = temp_dir / "orca_output.out"
                shutil.copy(inp_file, temp_inp_file)

                cmd = [str(self.orca_binary), str(temp_inp_file.name)]

                env = os.environ.copy()
                env['ORCA_TEMP_DIR'] = str(temp_dir)

                with open(temp_out_file, 'w') as out_f:
                    process = subprocess.Popen(
                        cmd,
                        stdout=out_f,
                        stderr=subprocess.PIPE,
                        cwd=str(temp_dir),
                        text=True,
                        env=env
                    )

                    try:
                        _, stderr = process.communicate(timeout=timeout)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        raise TimeoutError(
                            f"ORCA 计算超时 (>{timeout}秒): {inp_file.name}"
                        )

                    if process.returncode != 0:
                        raise RuntimeError(
                            f"ORCA 运行失败 (返回码 {process.returncode})\n"
                            f"错误信息: {stderr}"
                        )

                shutil.copy(temp_out_file, out_file)

                for f in temp_dir.glob('*'):
                    if f.is_file() and f not in [temp_inp_file, temp_out_file]:
                        shutil.copy(f, output_dir / f.name)

                return out_file

            finally:
                shutil.rmtree(temp_dir, ignore_errors=True)

        inp_file_abs = inp_file.resolve()
        cmd = [str(self.orca_binary), str(inp_file_abs)]

        try:
            with open(out_file, 'w') as out_f:
                process = subprocess.Popen(
                    cmd,
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                    cwd=str(output_dir),
                    text=True
                )

                # 等待进程完成或超时
                try:
                    _, stderr = process.communicate(timeout=timeout)
                except subprocess.TimeoutExpired:
                    process.kill()
                    raise TimeoutError(
                        f"ORCA 计算超时 (>{timeout}秒): {inp_file.name}"
                    )

                # 检查返回码
                if process.returncode != 0:
                    raise RuntimeError(
                        f"ORCA 运行失败 (返回码 {process.returncode})\n"
                        f"错误信息: {stderr}"
                    )

        except Exception as e:
            raise RuntimeError(f"ORCA 运行出错: {e}")

        return out_file

    def single_point(
        self,
        xyz_file: Path,
        output_dir: Path,
        timeout: int = 3600
    ) -> QCResult:
        """
        执行单点能计算（端到端流程）

        整合流程:
        1. 生成 ORCA 输入文件 (_generate_input)
        2. 运行 ORCA 计算 (_run_orca)
        3. 解析输出文件 (_parse_output)

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录
            timeout: 超时时间（秒），默认 1 小时

        Returns:
            QCResult 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"开始 ORCA 单点能计算: {xyz_file.name}")
        self.logger.info(f"  方法: {self.method}/{self.basis}")
        self.logger.info(f"  输出目录: {output_dir}")

        try:
            # 步骤 1: 生成输入文件
            self.logger.debug("  生成 ORCA 输入文件...")
            inp_file = self._generate_input(xyz_file, output_dir)
            self.logger.debug(f"  输入文件: {inp_file}")

            # 步骤 2: 运行 ORCA
            self.logger.debug("  运行 ORCA 计算...")
            out_file = self._run_orca(inp_file, output_dir, timeout=timeout)
            self.logger.debug(f"  输出文件: {out_file}")

            # 步骤 3: 解析输出
            self.logger.debug("  解析 ORCA 输出...")
            result = self._parse_output(out_file)

            if result.converged:
                self.logger.info(f"  ✓ 计算成功: 能量 = {result.energy:.8f} Hartree")
            else:
                self.logger.error(f"  ✗ 计算失败: {result.error_message}")

            return result

        except Exception as e:
            self.logger.error(f"ORCA 单点能计算失败: {e}")
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=str(e)
            )

    def optimize(
        self,
        xyz_file: Path,
        output_dir: Path,
        route: Optional[str] = None,
        constraints: Optional[str] = None,
        old_checkpoint: Optional[Path] = None,
        timeout: Optional[int] = None
    ) -> QCResult:
        """
        ORCA 几何优化（统一接口，兼容 Gaussian 风格）

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录
            route: Gaussian 风格的 route card（会被转换为 ORCA 关键词）
            constraints: 约束（暂时忽略，ORCA 不支持相同格式）
            old_checkpoint: checkpoint 文件（暂时忽略）
            timeout: 超时时间（秒）

        Returns:
            QCResult 对象
        """
        from rph_core.utils.optimization_config import OptimizationConfig

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 解析 route card，提取方法信息
        if route:
            # 简化处理：使用已配置的方法和基组
            # 实际项目中应该更精确地解析 route card
            if 'TS' in route:
                # TS 优化
                return self.ts_optimization(
                    xyz_file, output_dir,
                    opt_config=OptimizationConfig(),
                    timeout=timeout
                )
            else:
                # 基态优化 - 使用 Opt 关键词
                self.logger.info(f"ORCA 基态优化: {xyz_file.name}")
                return self._run_normal_optimization(
                    xyz_file, output_dir, timeout=timeout
                )

        # 默认基态优化
        return self._run_normal_optimization(
            xyz_file, output_dir, timeout=timeout
        )

    def _run_normal_optimization(
        self,
        xyz_file: Path,
        output_dir: Path,
        timeout: Optional[int] = None
    ) -> QCResult:
        """
        ORCA 基态优化（使用 Opt 关键词）

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录
            timeout: 超时时间

        Returns:
            QCResult 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 生成优化输入文件
        import shutil
        xyz_copy = output_dir / xyz_file.name
        shutil.copy(xyz_file, xyz_copy)

        # 读取 XYZ 内容
        xyz_lines = xyz_copy.read_text().split('\n')
        if len(xyz_lines) > 2:
            xyz_content = '\n'.join(xyz_lines[2:])
        else:
            xyz_content = xyz_copy.read_text()

        # 构建路由行
        route = f"! {self.method} {self.basis} {self.aux_basis} Opt tightSCF Freq"
        route += " noautostart miniprint nopop"

        # 双杂化泛函需要额外的 /C 辅助基组
        if self._is_double_hybrid():
            c_basis = self.basis + "/C"
            route += f" {c_basis}"

        # 构建溶剂块
        cpcm_block = ""
        if self.solvent and self.solvent.upper() != "NONE":
            cpcm_block = f"""
%cpcm
   smd true
   SMDsolvent "{self.solvent}"
end
"""

        # 组装完整输入文件
        inp_content = f"""{route}
%maxcore {self.maxcore}
%pal nprocs {self.nprocs} end
{cpcm_block}
 * xyzfile 0 1 {xyz_copy.name}
 *
"""

        inp_file = output_dir / f"{xyz_file.stem}_opt.inp"
        inp_file.write_text(inp_content)

        # 运行 ORCA
        try:
            out_file = self._run_orca(inp_file, output_dir, timeout=timeout)
            result = self._parse_output(out_file)

            # 提取频率
            freq_block = re.search(r'VIBRATIONAL FREQUENCIES\s*\n(.*?)(?=\n\n|\n[A-Z])', out_file.read_text(), re.DOTALL)
            if freq_block:
                freqs = re.findall(r'[\-]?\d+\.\d+', freq_block.group(1))
                result.frequencies = np.array([float(f) for f in freqs]) if freqs else None

            # 提取坐标
            coord_block = re.search(r'CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n(.*?)(?=\n\n|\n[A-Z])', out_file.read_text(), re.DOTALL)
            if coord_block:
                coord_lines = [l for l in coord_block.group(1).split('\n') if l.strip()]
                if len(coord_lines) >= 3:
                    result.coordinates = np.array([
                        [float(coord_lines[i+1].split()[1]),
                         float(coord_lines[i+1].split()[2]),
                         float(coord_lines[i+1].split()[3])]
                        for i in range(0, len(coord_lines)-1, 3)
                    ])

            return result

        except Exception as e:
            self.logger.error(f"ORCA 优化失败: {e}")
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=str(e)
            )

    def ts_optimization(
        self,
        xyz_file: Path,
        output_dir: Path,
        opt_config: Optional['OptimizationConfig'] = None,
        timeout: Optional[int] = None,
        charge: int = 0,
        spin: int = 1
    ) -> QCResult:
        """
        ORCA 过渡态优化

        使用 OptTS 关键词进行 TS 优化，支持 Hessian 控制

        Args:
            xyz_file: 输入 TS 猜想 XYZ 文件
            output_dir: 输出目录
            opt_config: 优化配置对象（来自 OptimizationConfig）
            timeout: 超时时间（秒），None = 无限制（QC 计算需要长时间）
            charge: 分子电荷
            spin: 自旋多重度

        Returns:
            QCResult 对象，包含优化后的坐标、能量和频率信息

        Raises:
            RuntimeError: ORCA 二进制文件未找到
            TimeoutError: 计算超时（如果设置了 timeout）
        """
        from rph_core.utils.optimization_config import OptimizationConfig

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"开始 ORCA TS 优化: {xyz_file.name}")
        self.logger.info(f"  方法: {self.method}/{self.basis}")
        self.logger.info(f"  电荷/自旋: {charge}/{spin}")

        # 如果没有提供 opt_config，使用默认值
        if opt_config is None:
            opt_config = OptimizationConfig()

        # 检查超时设置
        if timeout is None:
            if opt_config.timeout_enabled and opt_config.timeout_seconds:
                timeout = opt_config.timeout_seconds
                self.logger.info(f"  超时: {timeout} 秒")
            else:
                timeout = None
                self.logger.info("  超时: 禁用（无限制）")

        try:
            # 生成 TS 优化输入文件
            self.logger.debug("  生成 ORCA TS 优化输入文件...")
            inp_file = self._generate_ts_input(
                xyz_file, output_dir, opt_config, charge, spin
            )
            self.logger.debug(f"  输入文件: {inp_file}")

            # 运行 ORCA（无超时或自定义超时）
            self.logger.debug("  运行 ORCA TS 优化...")
            if timeout is None:
                # 无超时限制（QC 计算通常需要长时间）
                out_file = self._run_orca_no_timeout(inp_file, output_dir)
            else:
                out_file = self._run_orca(inp_file, output_dir, timeout=timeout)
            self.logger.debug(f"  输出文件: {out_file}")

            # 解析输出
            self.logger.debug("  解析 ORCA TS 优化输出...")
            result = self._parse_ts_output(out_file)

            if result.converged:
                self.logger.info(
                    f"  ✓ TS 优化成功: 能量 = {result.energy:.8f} Hartree"
                )
                if result.frequencies is not None:
                    imaginary = [f for f in result.frequencies if f < 0]
                    if imaginary:
                        self.logger.info(
                            f"  虚频: {len(imaginary)} 个, "
                            f"最小 = {min(imaginary):.1f} cm⁻¹"
                        )
            else:
                self.logger.error(f"  ✗ TS 优化失败: {result.error_message}")

            return result

        except Exception as e:
            self.logger.error(f"ORCA TS 优化失败: {e}")
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=str(e)
            )

    def _generate_ts_input(
        self,
        xyz_file: Path,
        output_dir: Path,
        opt_config: 'OptimizationConfig',
        charge: int,
        spin: int
    ) -> Path:
        """
        生成 ORCA TS 优化输入文件

        Args:
            xyz_file: 输入 XYZ 文件
            output_dir: 输出目录
            opt_config: 优化配置
            charge: 分子电荷
            spin: 自旋多重度

        Returns:
            生成的 .inp 文件路径
        """
        import shutil

        # 构建路由行
        route = f"! {self.method} {self.basis} {self.aux_basis} OptTS Freq tightSCF"
        route += " noautostart miniprint nopop"

        # 双杂化泛函需要额外的 /C 辅助基组
        if self._is_double_hybrid():
            c_basis = self.basis + "/C"
            route += f" {c_basis}"

        # 构建溶剂块
        cpcm_block = ""
        if self.solvent and self.solvent.upper() != "NONE":
            cpcm_block = f"""
%cpcm
   smd true
   SMDsolvent "{self.solvent}"
end
"""

        # 如果 xyz_file 是 .out 文件，尝试从中提取电荷和自旋
        if xyz_file.suffix == '.out':
            try:
                charge, spin = CoordinateExtractor.get_charge_spin_from_orca_out(xyz_file)
                self.logger.debug(f"从 {xyz_file.name} 提取电荷/自旋: {charge}/{spin}")
            except Exception as e:
                self.logger.warning(f"无法从 {xyz_file.name} 提取电荷/自旋: {e}")

        # 读取 XYZ 内容
        xyz_copy = output_dir / xyz_file.name
        shutil.copy(xyz_file, xyz_copy)
        xyz_lines = xyz_copy.read_text().split('\n')
        if len(xyz_lines) > 2:
            xyz_content = '\n'.join(xyz_lines[2:])
        else:
            xyz_content = xyz_copy.read_text()

        # 生成 %geom 块
        geom_block = opt_config.to_orca_geom_block(is_ts=True)

        # 组装完整输入文件
        inp_content = f"""{route}
%maxcore {self.maxcore}
%pal nprocs {self.nprocs} end
{cpcm_block}
{geom_block}
 * xyzfile {charge} {spin} {xyz_copy.name}
 *
"""

        # 写入文件
        inp_file = output_dir / f"{xyz_file.stem}_ts_opt.inp"
        inp_file.write_text(inp_content)

        return inp_file

    def _run_orca_no_timeout(
        self,
        inp_file: Path,
        output_dir: Path
    ) -> Path:
        """
        运行 ORCA 计算（无超时限制）

        Args:
            inp_file: ORCA 输入文件
            output_dir: 输出目录

        Returns:
            ORCA 输出文件路径

        Raises:
            RuntimeError: ORCA 二进制文件未找到或运行失败
        """
        if self.orca_binary is None:
            raise RuntimeError(
                "ORCA 二进制文件未找到。请通过以下方式之一指定:\n"
                "  1. 设置环境变量 $ORCA_PATH\n"
                "  2. 确保 'orca' 在 PATH 中\n"
                "  3. 在初始化时提供 orca_binary_path 参数"
            )

        # 使用绝对路径确保 ORCA 能找到输入文件
        inp_file_abs = inp_file.resolve()
        out_file = inp_file.with_suffix('.out')
        cmd = [str(self.orca_binary), str(inp_file_abs)]

        try:
            with open(out_file, 'w') as out_f:
                process = subprocess.Popen(
                    cmd,
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                    cwd=str(output_dir),
                    text=True
                )

                # 等待进程完成（无超时限制）
                _, stderr = process.communicate()

                # 检查返回码
                if process.returncode != 0:
                    raise RuntimeError(
                        f"ORCA 运行失败 (返回码 {process.returncode})\n"
                        f"错误信息: {stderr}"
                    )

        except Exception as e:
            raise RuntimeError(f"ORCA 运行出错: {e}")

        return out_file

    def _parse_ts_output(self, out_file: Path) -> QCResult:
        """
        解析 ORCA TS 优化输出文件

        提取能量、收敛状态、坐标和频率信息

        Args:
            out_file: ORCA 输出文件路径

        Returns:
            QCResult 对象
        """
        import numpy as np

        # 读取输出文件内容
        try:
            content = out_file.read_text()
        except Exception as e:
            return QCResult(
                energy=0.0,
                converged=False,
                error_message=f"无法读取输出文件: {e}"
            )

        # 检查是否正常终止
        if "ORCA TERMINATED NORMALLY" not in content:
            return QCResult(
                energy=0.0,
                converged=False,
                error_message="ORCA 未正常终止"
            )

        # 检查收敛状态
        converged = "THE OPTIMIZATION HAS CONVERGED" in content

        # 提取最终能量
        energy_match = re.search(r'FINAL SINGLE POINT ENERGY\s+([\-\d\.]+)', content)
        energy = float(energy_match.group(1)) if energy_match else 0.0

        # 提取频率
        frequencies = None
        freq_block = re.search(r'VIBRATIONAL FREQUENCIES\s*\n(.*?)(?=\n\n|\n[A-Z])', content, re.DOTALL)
        if freq_block:
            freqs = re.findall(r'[\-]?\d+\.\d+', freq_block.group(1))
            frequencies = np.array([float(f) for f in freqs]) if freqs else None

        # 提取坐标
        coordinates = None
        coord_block = re.search(r'CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n(.*?)(?=\n\n|\n[A-Z])', content, re.DOTALL)
        if coord_block:
            coord_lines = [l for l in coord_block.group(1).split('\n') if l.strip()]
            if len(coord_lines) >= 3:
                coordinates = np.array([
                    [float(coord_lines[i+1].split()[1]),
                     float(coord_lines[i+1].split()[2]),
                     float(coord_lines[i+1].split()[3])]
                    for i in range(0, len(coord_lines)-1, 3)
                ])

        error_message = None if converged else "优化未收敛"

        return QCResult(
            energy=energy,
            converged=converged,
            coordinates=coordinates if coordinates is not None else np.array([]),
            frequencies=frequencies,
            output_file=out_file,
            error_message=error_message
        )
