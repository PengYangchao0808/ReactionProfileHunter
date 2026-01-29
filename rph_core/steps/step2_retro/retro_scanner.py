"""
Step 2: Retro Scanner
======================

逆向扫描模块 - 从产物逆向生成TS初猜和底物

Author: QCcalc Team
Date: 2026-01-09
"""

import json
import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Optional, List, Any

import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import read_xyz, write_xyz
from rph_core.utils.geometry_tools import LogParser
from rph_core.utils.data_types import QCResult
from rph_core.utils.xtb_runner import XTBRunner
from rph_core.utils.tsv_dataset import ReactionRecord
from .smarts_matcher import SMARTSMatcher
from .bond_stretcher import BondStretcher

logger = logging.getLogger(__name__)


@dataclass
class RetroScanResultV2:
    """Result object for retro scanning with extended contract."""
    ts_guess_xyz: Optional[Path]
    reactant_xyz: Optional[Path]
    forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]]
    neutral_precursor_xyz: Optional[Path]
    meta_json_path: Optional[Path]


class RetroScanner(LoggerMixin):
    """
    逆向扫描引擎 (Step 2) - v3.0

    v2.1 核心创新: Product-First 策略
    v3.0 更新: 适配分子自治目录结构

    职责:
    1. 识别产物中的形成键 (SMARTS)
    2. 路径 A: 生成 TS 初猜 (拉伸至 2.2Å + 受限优化)
    3. 路径 B: 生成底物 (拉伸至 3.5Å + 无限制优化)

    输入:
    - product_xyz: 产物全局最低构象 (来自 Step 1，可能是文件或目录)

    输出:
    - ts_guess_xyz: TS 初猜结构
    - reactant_xyz: 底物复合物结构
    """

    # 默认参数
    DEFAULT_TS_DISTANCE = 2.2  # Å - TS 典型距离 (✅ PROMOTE.md 标准)
    DEFAULT_BREAK_DISTANCE = 3.50  # Å - 断裂距离

    def __init__(self, config: dict, molecule_name: Optional[str] = None):
        """
        初始化逆向扫描引擎 - v3.0

        Args:
            config: 配置字典
            molecule_name: 分子名称（用于定位 v3.0 目录结构）
        """
        self.config = config
        self.ts_distance = config.get('ts_distance', self.DEFAULT_TS_DISTANCE)
        self.break_distance = config.get('break_distance', self.DEFAULT_BREAK_DISTANCE)
        self.molecule_name = molecule_name

        self.xtb_runner = XTBRunner(config)

        # 初始化 SMARTS 匹配器
        self.smarts_matcher = SMARTSMatcher()

        # 初始化键拉伸器
        self.bond_stretcher = BondStretcher()

        self.logger.info(f"RetroScanner v3.0 初始化: TS距离={self.ts_distance}Å, 断裂距离={self.break_distance}Å")
        if molecule_name:
            self.logger.info(f"  目标分子: {molecule_name}")

    def run(self, product_xyz: Path, output_dir: Path, molecule_name: Optional[str] = None) -> Tuple[Path, Path, Tuple[Tuple[int, int], Tuple[int, int]]]:
        """
        执行逆向扫描 - v3.0 (Legacy Path)

        NOTE: This is the legacy run() method maintained for backward compatibility.
        For new code requiring neutral precursor support, use run_with_precursor() instead.

        Args:
            product_xyz: 产物 XYZ 文件路径（v2.1）或 S1_Product 目录（v3.0）
            output_dir: 输出目录
            molecule_name: 分子名称（v3.0，用于定位分子自治目录）

        Returns:
            (ts_guess_xyz, reactant_xyz, forming_bonds) 元组
            - forming_bonds: ((atom_idx1, atom_idx2), (atom_idx3, atom_idx4)) 形成键索引
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # v3.0 目录适配逻辑
        if product_xyz.is_dir():
            # v3.0: product_xyz 是 S1_Product 目录
            self.logger.info("检测到 v3.0 目录结构")

            # 尝试定位分子目录
            if molecule_name:
                molecule_dir = product_xyz / molecule_name
            else:
                # 自动查找分子目录
                molecule_dirs = [d for d in product_xyz.iterdir() if d.is_dir() and not d.name.startswith('.')]
                if len(molecule_dirs) == 1:
                    molecule_dir = molecule_dirs[0]
                    self.logger.info(f"自动检测到分子目录: {molecule_dir.name}")
                elif len(molecule_dirs) > 1:
                    raise RuntimeError(
                        f"发现多个分子目录: {[d.name for d in molecule_dirs]}。"
                        f"请显式指定 molecule_name 参数。"
                    )
                else:
                    raise RuntimeError(
                        f"在 {product_xyz} 中未找到分子目录。"
                    )

            # 查找全局最低构象文件
            global_min_patterns = [
                f"{molecule_dir.name}_global_min.xyz",
                "global_min.xyz",
            ]

            product_min_xyz = None
            for pattern in global_min_patterns:
                candidate = molecule_dir / pattern
                if candidate.exists():
                    product_min_xyz = candidate
                    self.logger.info(f"找到产物文件: {product_min_xyz}")
                    break

            if product_min_xyz is None:
                # 回退：查找 dft 目录中的 SP 输出
                dft_dir = molecule_dir / "dft"
                if dft_dir.exists():
                    sp_files = list(dft_dir.glob("*_SP.out"))
                    if sp_files:
                        product_min_xyz = sp_files[0]
                        self.logger.info(f"使用 SP 输出文件: {product_min_xyz}")

            if product_min_xyz is None:
                raise RuntimeError(
                    f"无法在 {molecule_dir} 中找到产物结构文件。"
                    f"查找的文件: {global_min_patterns}"
                )
        else:
            # v2.1: product_xyz 直接是文件路径
            self.logger.info("使用 v2.1 直接文件路径")
            product_min_xyz = product_xyz

        self.logger.info("Step 2 开始: 逆向扫描")
        self.logger.info(f"输入产物: {product_min_xyz}")

        # 使用 LogParser 提取坐标（确保兼容性）
        coords, symbols, error = LogParser.extract_last_converged_coords(
            product_min_xyz,
            engine_type='auto'
        )

        if coords is None:
            self.logger.warning(f"LogParser 失败: {error}，回退到常规 XYZ 读取")
            coords, symbols = read_xyz(product_min_xyz)
        else:
            self.logger.info(f"使用 LogParser 成功提取 {len(coords)} 个原子坐标")
            if symbols is None:
                self.logger.warning("未提取到元素符号，回退到常规 XYZ 读取")
                _, symbols = read_xyz(product_min_xyz)

        # 1. 识别反应键位点 (SMARTS 匹配)
        self.logger.info("识别反应键位点...")
        match_result = self.smarts_matcher.find_reactive_bonds(product_min_xyz)

        if not match_result.matched:
            raise ValueError(f"SMARTS 匹配失败: {match_result.error_message}")

        if match_result.bond_1 is None or match_result.bond_2 is None:
            raise ValueError("SMARTS 匹配失败: 缺少形成键信息")

        bond_1 = (match_result.bond_1.atom_idx_1, match_result.bond_1.atom_idx_2)
        bond_2 = (match_result.bond_2.atom_idx_1, match_result.bond_2.atom_idx_2)
        bond_1_length = match_result.bond_1.current_length
        bond_2_length = match_result.bond_2.current_length

        self.logger.info(f"  识别到形成键1: {bond_1} (当前 {bond_1_length:.2f}Å)")
        self.logger.info(f"  识别到形成键2: {bond_2} (当前 {bond_2_length:.2f}Å)")
        self.logger.info(f"  匹配模式: {match_result.pattern_name} (置信度 {match_result.confidence:.2f})")

        # ==========================================
        # 路径 A: 生成 TS 初猜 (拉伸至 TS 距离)
        # ==========================================
        self.logger.info(f"路径 A: 拉伸键至 {self.ts_distance}Å (TS 初猜)...")

        ts_raw_coords = self.bond_stretcher.stretch_two_bonds(
            coords.copy(), bond_1, bond_2, target_length=self.ts_distance
        )

        # 保存拉伸后的结构
        ts_raw_xyz = output_dir / "ts_raw_stretched.xyz"
        write_xyz(ts_raw_xyz, ts_raw_coords, symbols, title="TS Raw Stretched")

        # 受限优化：保持键长固定
        ts_constraints = {
            f"{bond_1[0]+1} {bond_1[1]+1}": self.ts_distance,
            f"{bond_2[0]+1} {bond_2[1]+1}": self.ts_distance
        }

        xtb_settings = self.config.get('step2', {}).get('xtb_settings', {})
        solvent = xtb_settings.get('solvent', 'acetone')
        charge = 0
        uhf = 0

        ts_opt_dir = output_dir / "ts_opt"
        ts_opt_dir.mkdir(parents=True, exist_ok=True)
        self.xtb_runner.work_dir = ts_opt_dir
        self.logger.info(f"运行 XTB 受限优化 (溶剂={solvent})...")
        ts_result = self.xtb_runner.optimize(
            structure=ts_raw_xyz,
            constraints=ts_constraints,
            solvent=solvent,
            charge=charge,
            uhf=uhf
        )

        if ts_result.success:
            ts_guess_xyz = output_dir / "ts_guess.xyz"
            if ts_result.output_file is None:
                ts_guess_xyz = ts_raw_xyz
                self.logger.warning("TS 优化成功但缺少输出文件，回退到拉伸结构")
            else:
                ts_coordinates = Path(ts_result.output_file)
                try:
                    shutil.copy(ts_coordinates, ts_guess_xyz)
                except shutil.SameFileError:
                    pass
                self.logger.info(f"TS 受限优化成功: {ts_guess_xyz}")
        else:
            ts_guess_xyz = ts_raw_xyz
            self.logger.warning(f"TS 受限优化失败: {ts_result.error_message}, 使用拉伸后结构")


        # ==========================================
        # 路径 B: 生成底物 (拉伸至断裂距离 + 松弛)
        # ==========================================
        self.logger.info(f"路径 B: 拉伸键至 {self.break_distance}Å (底物)...")

        reactant_raw_coords = self.bond_stretcher.stretch_two_bonds(
            coords.copy(), bond_1, bond_2, target_length=self.break_distance
        )

        reactant_raw_xyz = output_dir / "reactant_raw_stretched.xyz"
        write_xyz(reactant_raw_xyz, reactant_raw_coords, symbols, title="Reactant Raw Stretched")

        reactant_opt_dir = output_dir / "reactant_opt"
        reactant_opt_dir.mkdir(parents=True, exist_ok=True)
        self.xtb_runner.work_dir = reactant_opt_dir
        self.logger.info(f"运行 XTB 无约束松弛优化 (溶剂={solvent})...")
        reactant_result = self.xtb_runner.optimize(
            structure=reactant_raw_xyz,
            constraints=None,
            solvent=solvent,
            charge=charge,
            uhf=uhf
        )

        if reactant_result.success:
            reactant_xyz = output_dir / "reactant_complex.xyz"
            if reactant_result.output_file is None:
                reactant_xyz = reactant_raw_xyz
                self.logger.warning("底物优化成功但缺少输出文件，回退到拉伸结构")
            else:
                reactant_coordinates = Path(reactant_result.output_file)
                try:
                    shutil.copy(reactant_coordinates, reactant_xyz)
                except shutil.SameFileError:
                    pass
                self.logger.info(f"底物松弛优化成功: {reactant_xyz}")
        else:
            reactant_xyz = reactant_raw_xyz
            self.logger.warning(f"底物松弛优化失败: {reactant_result.error_message}, 使用拉伸后结构")


        self.logger.info(f"Step 2 完成:")
        self.logger.info(f"  TS 初猜: {ts_guess_xyz}")
        self.logger.info(f"  底物: {reactant_xyz}")

        # 返回初猜路径和识别到的键指数
        return (ts_guess_xyz, reactant_xyz, (bond_1, bond_2))

    def run_with_precursor(
        self,
        reactant_complex_xyz: Path,
        record: ReactionRecord,
        output_dir: Path,
        enabled: Optional[bool] = None,
        strategy: str = "reactant_complex",
        output_meta: bool = False
    ) -> RetroScanResultV2:
        """
        Run retro scan with neutral precursor support (V5.2-M1 contract).

        Args:
            reactant_complex_xyz: Path to reactant_complex.xyz (from run() method)
            record: Reaction record containing precursor_smiles
            output_dir: Output directory for neutral_precursor.xyz and meta.json
            enabled: Override neutral_precursor.enabled (None = use config default)
            strategy: Strategy for neutral precursor generation ("reactant_complex" only supported)
            output_meta: Whether to write meta.json

        Returns:
            RetroScanResultV2 with:
            - ts_guess_xyz, reactant_xyz, forming_bonds: None (not generated by this method)
            - neutral_precursor_xyz: Path to neutral precursor (or None if disabled/error)
            - meta_json_path: Path to meta.json (or None if not written)

        Note:
            This method should be called after run() generates reactant_complex.xyz.
            For "reactant_complex" strategy, neutral precursor = reactant complex (direct copy).
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        check_enabled = enabled
        if check_enabled is None:
            check_enabled = self.config.get('step2', {}).get('neutral_precursor', {}).get('enabled', False)

        if not check_enabled:
            self.logger.info("Neutral precursor generation disabled (enabled=False)")
            return RetroScanResultV2(
                ts_guess_xyz=None,
                reactant_xyz=None,
                forming_bonds=None,
                neutral_precursor_xyz=None,
                meta_json_path=None
            )

        self.logger.info(f"Neutral precursor generation enabled with strategy: {strategy}")

        if strategy != "reactant_complex":
            self.logger.warning(f"Unsupported strategy: {strategy}, defaulting to reactant_complex")
            strategy = "reactant_complex"

        reactant_complex_xyz = Path(reactant_complex_xyz)
        if not reactant_complex_xyz.exists():
            self.logger.error(f"Reactant complex file not found: {reactant_complex_xyz}")
            return RetroScanResultV2(
                ts_guess_xyz=None,
                reactant_xyz=None,
                forming_bonds=None,
                neutral_precursor_xyz=None,
                meta_json_path=None
            )

        neutral_precursor_xyz = output_dir / "neutral_precursor.xyz"
        try:
            shutil.copy(reactant_complex_xyz, neutral_precursor_xyz)
            self.logger.info(f"Neutral precursor created: {neutral_precursor_xyz}")
        except Exception as e:
            self.logger.error(f"Failed to copy reactant complex to neutral precursor: {e}")
            return RetroScanResultV2(
                ts_guess_xyz=None,
                reactant_xyz=None,
                forming_bonds=None,
                neutral_precursor_xyz=None,
                meta_json_path=None
            )

        meta_json_path = None
        if output_meta:
            meta_json_path = output_dir / "meta.json"
            meta_data = {
                'precursor_smiles': record.precursor_smiles,
                'leaving_small_molecule_key': record.get_leaving_small_molecule_key(),
                'strategy': strategy,
                'source_reactant_complex': str(reactant_complex_xyz)
            }
            try:
                import json
                meta_json_path.write_text(json.dumps(meta_data, indent=2))
                self.logger.info(f"Meta JSON written: {meta_json_path}")
            except Exception as e:
                self.logger.warning(f"Failed to write meta.json: {e}")
                meta_json_path = None

        return RetroScanResultV2(
            ts_guess_xyz=None,
            reactant_xyz=None,
            forming_bonds=None,
            neutral_precursor_xyz=neutral_precursor_xyz,
            meta_json_path=meta_json_path
        )

    def _generate_constraints(self, bond_1: Tuple[int, int], bond_2: Tuple[int, int], target_dist: float) -> str:
        """
        生成 XTB 约束文件内容

        Args:
            bond_1: 第一个键 (i, j)
            bond_2: 第二个键 (k, l)
            target_dist: 目标距离

        Returns:
            约束文件内容字符串
        """
        # XTB 使用 1-indexed
        return f"""$constrain
  force constant=0.5
  distance: {bond_1[0]+1}, {bond_1[1]+1}, {target_dist:.3f}
  distance: {bond_2[0]+1}, {bond_2[1]+1}, {target_dist:.3f}
$end
"""
