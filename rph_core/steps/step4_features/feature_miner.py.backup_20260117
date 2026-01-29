"""
Step 4: Feature Miner
======================

特征提取模块 - 畸变能、几何、电子特征提取

Author: QCcalc Team
Date: 2026-01-09
Updated: 2026-01-11 (Session #14: 集成 L2 SP 矩阵)
"""

import logging
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Optional, Tuple, List

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import read_xyz
from .log_parser import GaussianLogParser
from .distortion_calculator import DistortionCalculator
from .electronic_extractor import ElectronicExtractor
from .fragment_extractor import FragmentExtractor

logger = logging.getLogger(__name__)


# 常量
HARTREE_TO_KCAL = 627.5094740631  # Hartree to kcal/mol


class FeatureMiner(LoggerMixin):
    """
    特征提取器 (Step 4)

    提取特征:
    - 能量特征 (dG‡, dGr, Edist, Eint)
    - 几何特征 (r1, r2, asynchronicity)
    - 电子特征 (charges, HOMO/LUMO)
    - 位阻特征 (buried volume)

    v2.1 新增:
    - 双片段畸变能 (E_dist,A + E_dist,B)
    - Checkpoint复用支持
    - [Session #14] L2 高精度能量特征 (优先使用 SP 矩阵)

    输入:
    - ts_final: TS 结构 (来自 Step 3)
    - reactant: 底物 (来自 Step 2)
    - product: 产物 (来自 Step 1)
    - sp_matrix_report: S3.5 SP 矩阵报告 [NEW in Session #14]

    输出:
    - features.csv: 特征表 (包含 L1 和 L2 精度特征)
    """

    def __init__(self, config: dict):
        """
        初始化特征提取器

        Args:
            config: 配置字典
        """
        self.config = config
        self.enable_nbo = config.get('enable_nbo', False)

        self.log_parser = GaussianLogParser()
        self.dist_calc = DistortionCalculator()
        self.elec_extractor = ElectronicExtractor()

        # P0: 真实的片段能量计算
        dft_config = {
            'method': config.get('method', 'B3LYP'),
            'basis': config.get('basis', 'def2-SVP'),
            'dispersion': config.get('dispersion', 'GD3BJ'),
            'nprocshared': config.get('nprocshared', 16),
            'mem': config.get('mem', '32GB'),
            'solvent': config.get('solvent', 'acetone')
        }
        self.fragment_extractor = FragmentExtractor(dft_config)

        self.logger.info("FeatureMiner 初始化完成 (支持 L2 SP 矩阵)")

    def run(
        self,
        ts_final: Path,
        reactant: Path,
        product: Path,
        output_dir: Path,
        forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]] = None,
        fragment_indices: Optional[Tuple[List[int], List[int]]] = None,
        sp_matrix_report: Optional[object] = None  # [NEW in Session #14]
    ) -> Path:
        """
        执行特征提取 (支持 L2 高精度能量)

        Args:
            ts_final: TS 结构 (XYZ 或 .log)
            reactant: 底物结构 (XYZ 或 .log)
            product: 产物结构 (XYZ 或 .log)
            output_dir: 输出目录
            forming_bonds: 形成键的原子索引 ((i, j), (k, l))
            fragment_indices: 双片段索引 (fragment_A_atoms, fragment_B_atoms)
            sp_matrix_report: S3.5 SP 矩阵报告 [NEW in Session #14]

        Returns:
            features.csv: 特征表路径
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("Step 4 开始: 特征提取 (Distortion/Interaction)")

        # [NEW in Session #14] 优先使用 L2 SP 矩阵
        if sp_matrix_report is not None:
            self.logger.info("  使用 S3.5 L2 SP 矩阵计算特征 (高精度)")
            features = self._extract_features_from_sp_matrix(
                sp_matrix_report,
                ts_final,
                forming_bonds
            )
        else:
            self.logger.warning("  未提供 SP 矩阵，使用 L1 能量 (精度较低)")
            # 确定输入文件类型（优先使用.log，否则使用.xyz）
            ts_input = self._resolve_input_file(ts_final)
            reactant_input = self._resolve_input_file(reactant)
            product_input = self._resolve_input_file(product)

            # 提取特征 (L1 精度)
            features = self._extract_features(ts_input, reactant_input, product_input,
                                             forming_bonds, fragment_indices)

        # 保存为 CSV
        features_csv = output_dir / "features.csv"
        self._save_to_csv(features, features_csv)

        self.logger.info(f"  ✓ 特征提取完成: {features_csv}")

        return features_csv

    def _resolve_input_file(self, input_path: Path) -> Path:
        """
        解析输入文件：优先使用.log，如果不存在则使用原始文件

        Args:
            input_path: 原始路径

        Returns:
            实际存在的文件路径
        """
        # 首先尝试.log文件
        log_path = input_path.with_suffix('.log')
        if log_path.exists():
            return log_path

        # 否则使用原始文件
        if input_path.exists():
            return input_path

        # 都不存在则报错
        raise FileNotFoundError(f"无法找到输入文件: {input_path} 或 {log_path}")

    def _extract_features(
        self,
        ts_input: Path,
        reactant_input: Path,
        product_input: Path,
        forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]] = None,
        fragment_indices: Optional[Tuple[List[int], List[int]]] = None
    ) -> Dict[str, float]:
        """提取真实的物理特征"""
        features = {}

        # 1. 解析所有日志（或能量文件）
        res_ts = self.log_parser.parse_log(ts_input)
        res_rea = self.log_parser.parse_log(reactant_input)
        res_pro = self.log_parser.parse_log(product_input)

        if not all([res_ts, res_rea, res_pro]):
            self.logger.error("解析日志失败，部分特征将缺失")
            return {}

        # 2. 能量特征 (kcal/mol)
        features["dG_activation"] = self.dist_calc.calculate_activation_energy(
            res_ts.energy, res_rea.energy
        )
        features["dG_reaction"] = self.dist_calc.calculate_reaction_energy(
            res_pro.energy, res_rea.energy
        )

        # 3. 双片段畸变能分析 (P0 - [5+2]反应核心)
        if fragment_indices and res_ts.coordinates.size > 0:
            # 检查是否有TS checkpoint可用于片段计算
            ts_checkpoint = ts_input.with_suffix('.chk')

            distortion_features = self._calculate_dual_fragment_distortion(
                res_ts=res_ts,
                res_rea=res_rea,
                fragment_indices=fragment_indices,
                ts_xyz=ts_input,
                ts_checkpoint=ts_checkpoint if ts_checkpoint.exists() else None,
                output_dir=output_dir.parent / "fragment_calculation"
            )
            features.update(distortion_features)
        else:
            # 后备方案：简单畸变能（单片段近似）
            self.logger.warning("未提供fragment_indices，使用单片段近似")
            features["E_distortion_sub"] = 0.0  # 占位符
            features["E_distortion_rea"] = 0.0  # 占位符
            features["E_interaction"] = features["dG_activation"]  # 简化处理
            features["E_distortion_total"] = 0.0

        # 4. 几何特征
        if res_ts.coordinates.size > 0 and forming_bonds:
            (i, j), (k, l) = forming_bonds

            # 计算 TS 下的键长
            def dist(idx1, idx2):
                p1 = res_ts.coordinates[idx1]
                p2 = res_ts.coordinates[idx2]
                return np.linalg.norm(p1 - p2)

            r1 = dist(i, j)
            r2 = dist(k, l)

            features["r1_forming"] = r1
            features["r2_forming"] = r2
            features["asynchronicity"] = self.dist_calc.calculate_asynchronicity(r1, r2)

        # 5. 电子特征
        elec_features = self.elec_extractor.extract(res_ts)
        features.update(elec_features)

        # 6. 反应活性指数
        if res_ts.homo and res_ts.lumo:
            indices = self.elec_extractor.calculate_reactivity_indices(
                res_ts.homo, res_ts.lumo
            )
            features.update(indices)

        return features

    def _calculate_dual_fragment_distortion(
        self,
        res_ts: 'QCResult',
        res_rea: 'QCResult',
        fragment_indices: Tuple[List[int], List[int]],
        ts_xyz: Path,
        ts_checkpoint: Optional[Path] = None,
        output_dir: Optional[Path] = None
    ) -> Dict[str, float]:
        """
        计算双片段畸变能 (P0优先级) - 真实DFT计算

        公式 (PROMOTE.md Section 7):
        - E_dist,A = E(A_TS) - E(A_relaxed)
        - E_dist,B = E(B_TS) - E(B_relaxed)
        - E_dist,total = E_dist,A + E_dist,B
        - E_int = ΔG‡ - E_dist,total

        实现策略:
        1. 使用FragmentExtractor从TS结构切分片段
        2. 对每个片段进行DFT单点计算（在TS几何下）
        3. 对每个片段进行几何优化（得到松弛能量）
        4. 计算畸变能

        Args:
            res_ts: TS计算结果
            res_rea: 底物计算结果
            fragment_indices: (fragment_A_atom_indices, fragment_B_atom_indices)
            ts_xyz: TS XYZ文件路径
            ts_checkpoint: TS checkpoint文件（复用轨道）
            output_dir: 片段计算输出目录

        Returns:
            畸变能特征字典
        """
        if output_dir is None:
            output_dir = ts_xyz.parent / "fragment_calculation"

        # 使用真实的片段提取和DFT计算
        self.logger.info("使用真实DFT计算片段能量...")
        fragment_energies = self.fragment_extractor.extract_and_calculate(
            ts_xyz=ts_xyz,
            fragment_indices=fragment_indices,
            output_dir=output_dir,
            old_checkpoint=ts_checkpoint,
            apply_bsse=False  # TODO: 可选的BSSE校正
        )

        e_fragment_A_at_ts = fragment_energies['e_fragment_a_ts']
        e_fragment_B_at_ts = fragment_energies['e_fragment_b_ts']
        e_fragment_A_relaxed = fragment_energies['e_fragment_a_relaxed']
        e_fragment_B_relaxed = fragment_energies['e_fragment_b_relaxed']

        # 使用DistortionCalculator计算
        distortion_results = self.dist_calc.calculate_distortion_interaction(
            e_ts=res_ts.energy,
            e_fragment_a_at_ts=e_fragment_A_at_ts,
            e_fragment_b_at_ts=e_fragment_B_at_ts,
            e_fragment_a_relaxed=e_fragment_A_relaxed,
            e_fragment_b_relaxed=e_fragment_B_relaxed
        )

        self.logger.info("✓ 双片段畸变能计算完成 (真实DFT):")
        self.logger.info(f"  E_dist,A = {distortion_results['e_distortion_a']:.2f} kcal/mol")
        self.logger.info(f"  E_dist,B = {distortion_results['e_distortion_b']:.2f} kcal/mol")
        self.logger.info(f"  E_dist,total = {distortion_results['e_distortion_total']:.2f} kcal/mol")
        self.logger.info(f"  E_int = {distortion_results['e_interaction']:.2f} kcal/mol")

        # 映射到feature名称
        return {
            "E_distortion_sub": distortion_results["e_distortion_a"],  # 底物片段畸变
            "E_distortion_rea": distortion_results["e_distortion_b"],   # 试剂片段畸变
            "E_interaction": distortion_results["e_interaction"],
            "E_distortion_total": distortion_results["e_distortion_total"]
        }

    def _extract_features_from_sp_matrix(
        self,
        sp_matrix_report: object,
        ts_final: Path,
        forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]] = None
    ) -> Dict[str, float]:
        """
        [NEW in Session #14] 基于 S3.5 SP 矩阵提取 L2 高精度能量特征

        使用 SPMatrixReport 中的 L2 能量计算高精度特征:
        - dG_activation_L2: L2 活化能
        - dG_reaction_L2: L2 反应能
        - E_distortion_A_L2, E_distortion_B_L2: L2 畸变能
        - E_interaction_L2: L2 相互作用能

        Args:
            sp_matrix_report: S3.5 SPMatrixReport 对象
            ts_final: TS XYZ 文件路径
            forming_bonds: 形成键原子索引

        Returns:
            特征字典 (包含 L2 能量特征)
        """
        features = {}

        # 1. L2 能量特征 (kcal/mol)
        features["dG_activation_L2"] = sp_matrix_report.get_activation_energy()
        features["dG_reaction_L2"] = sp_matrix_report.get_reaction_energy()

        # 2. L2 畸变能分析
        features["E_distortion_A_L2"] = sp_matrix_report.get_distortion_energy_a()
        features["E_distortion_B_L2"] = sp_matrix_report.get_distortion_energy_b()

        # 如果没有 relaxed 能量，畸变能为 0，总畸变能计算为 TS - Reactant
        if features["E_distortion_A_L2"] == 0.0 and features["E_distortion_B_L2"] == 0.0:
            # 降级：使用总畸变能近似
            # E_dist,total ≈ ΔG‡ (忽略相互作用)
            features["E_distortion_total_L2"] = features["dG_activation_L2"]
        else:
            features["E_distortion_total_L2"] = (
                features["E_distortion_A_L2"] + features["E_distortion_B_L2"]
            )

        # 3. L2 相互作用能
        # E_int = ΔG‡ - E_dist,total
        features["E_interaction_L2"] = (
            features["dG_activation_L2"] - features["E_distortion_total_L2"]
        )

        # 4. 几何特征 (从 TS 结构读取)
        if forming_bonds and ts_final.exists():
            try:
                coords, symbols = read_xyz(ts_final)

                def dist(i, j):
                    p1 = coords[i-1]  # 原子索引从1开始
                    p2 = coords[j-1]
                    return np.linalg.norm(p1 - p2)

                (i, j), (k, l) = forming_bonds
                r1 = dist(i, j)
                r2 = dist(k, l)

                features["r1_forming"] = r1
                features["r2_forming"] = r2
                features["asynchronicity"] = abs(r1 - r2)

            except Exception as e:
                self.logger.warning(f"  几何特征提取失败: {e}")

        # 5. 元数据特征
        features["L2_method"] = sp_matrix_report.method
        features["L2_solvent"] = sp_matrix_report.solvent
        features["L2_available"] = 1.0  # 标记 L2 数据可用

        self.logger.info(f"  ✓ L2 特征提取完成:")
        self.logger.info(f"      ΔG‡_L2 = {features['dG_activation_L2']:.3f} kcal/mol")
        self.logger.info(f"      ΔG_rxn_L2 = {features['dG_reaction_L2']:.3f} kcal/mol")
        self.logger.info(f"      E_dist,A_L2 = {features['E_distortion_A_L2']:.3f} kcal/mol")
        self.logger.info(f"      E_dist,B_L2 = {features['E_distortion_B_L2']:.3f} kcal/mol")

        return features

    def _save_to_csv(self, features: Dict[str, float], output_path: Path):
        """
        保存特征到 CSV 文件

        Args:
            features: 特征字典
            output_path: 输出文件路径
        """
        # 转换为 DataFrame
        df = pd.DataFrame([features])

        # 保存
        df.to_csv(output_path, index=False)

        self.logger.info(f"✓ 特征已保存: {output_path}")
        self.logger.info(f"  特征数量: {len(features)}")

        # 打印特征摘要
        self.logger.info("特征摘要:")
        for key, value in features.items():
            if value is not None:
                self.logger.info(f"  {key}: {value}")
