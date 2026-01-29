"""
Anchor Phase Handler - v3.0 (Adapted for new ConformerEngine)
=======================================================
分子锚定阶段（Step 1）

职责:
1. 统一处理底物池（Substrate Pool）和产物（Product）
2. 使用 ConformerEngine 进行系综搜索 + OPT-SP 耦合
3. 输出：每个分子的全局最低能结构及其 SP 能量

Author: QCcalc Team
Date: 2026-01-14
Version: v3.0 (Molecule-Autonomous Architecture)
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

from rph_core.utils.log_manager import LoggerMixin
from rph_core.steps.conformer_search.engine import ConformerEngine

logger = logging.getLogger(__name__)


@dataclass
class AnchorPhaseResult:
    """锚定阶段结果（新结构）"""
    success: bool
    anchored_molecules: Dict[str, Dict[str, any]]
    # 格式: {"reactant_A": {"xyz": Path, "e_sp": float}, ...}
    error_message: Optional[str] = None


class AnchorPhase(LoggerMixin):
    """
    统一锚定阶段处理器 - v3.0

    职责：
    1. 接收分子列表（底物 + 产物）
    2. 对每个分子执行系综搜索 + DFT OPT-SP 耦合
    3. 返回全局最低能结构和 SP 能量

    输入示例:
        {
            "reactant_A": "CC(=O)C",
            "reactant_B": "C=C=C",
            "product": "CC1CC(=O)C(=O)O"
        }

    输出示例:
        {
            "reactant_A": {"xyz": path/to/A_SP.out, "e_sp": -234.5678 Hartree},
            "reactant_B": {"xyz": path/to/B_SP.out, "e_sp": -123.4567 Hartree},
            "product": {"xyz": path/to/P_SP.out, "e_sp": -238.9012 Hartree}
        }
    """

    def __init__(self, config: dict, base_work_dir: Path):
        """
        初始化 AnchorPhase - v3.0

        Args:
            config: 配置字典
            base_work_dir: 基础工作目录（例如 S1_Product/）
        """
        self.config = config
        self.base_work_dir = Path(base_work_dir)
        self.base_work_dir.mkdir(parents=True, exist_ok=True)

        # 初始化 ConformerEngine（但不在 run 方法中创建实例）
        # ConformerEngine 将在内部创建分子自治目录结构
        # 例如：S1_Product/[Molecule_Name]/crest 和 .../dft

        self.logger.info("AnchorPhase v3.0 初始化完成（分子自治架构）")

    def run(
        self,
        molecules: Dict[str, str]
    ) -> AnchorPhaseResult:
        """
        执行统一锚定流程（新结构）

        Args:
            molecules: 分子字典 {名称: SMILES}

        Returns:
            AnchorPhaseResult 对象
        """
        self.logger.info("=" * 60)
        self.logger.info("Anchor Phase: 统一锚定底物池和产物（v3.0）")
        self.logger.info("=" * 60)
        self.logger.info(f"分子数量: {len(molecules)}")

        anchored_molecules: Dict[str, Dict[str, any]] = {}
        error_messages = []

        # 对每个分子执行系综搜索 + OPT-SP 耦合
        for name, smiles in molecules.items():
            self.logger.info(f"\n>>> 锚定 {name}...")
            self.logger.info(f"    SMILES: {smiles}")

            try:
                # ConformerEngine 内部会创建：base_work_dir/[name]/crest 和 .../dft
                self.logger.info(f"    步骤 1/2: 系综搜索 + DFT OPT-SP 耦合...")

                # ConformerEngine 内部管理目录，不需要传递 work_dir
                # 每次调用 ConformerEngine.run 会创建该分子的自治目录结构
                temp_engine = ConformerEngine(
                    config=self.config,
                    work_dir=self.base_work_dir,
                    molecule_name=name
                )

                # 调用 ConformerEngine.run (OPT-SP 已耦合在内部）
                best_sp_out, sp_energy = temp_engine.run(smiles=smiles)

                self.logger.info(f"    ✓ OPT-SP 完成: {best_sp_out}")
                self.logger.info(f"    ✓ SP 能量: {sp_energy:.8f} Hartree")

                # 保存结果（新结构）
                anchored_molecules[name] = {
                    "xyz": best_sp_out,  # SP 输出文件路径
                    "e_sp": sp_energy  # SP 能量
                }

                self.logger.info(f"  ✓ {name} 锚定完成")

            except Exception as e:
                error_msg = f"{name} 锚定失败: {e}"
                self.logger.error(error_msg, exc_info=True)
                error_messages.append(error_msg)

                # 失败时仍然记录部分结果
                mol_dir = self.base_work_dir / name
                potential_xyz = list((mol_dir / "dft").glob("*_SP.out"))
                if potential_xyz:
                    anchored_molecules[name] = {
                        "xyz": potential_xyz[0],
                        "e_sp": None,
                        "failed": True
                    }

        # 检查是否至少有一个成功
        success = len(error_messages) < len(molecules)

        result = AnchorPhaseResult(
            success=success,
            anchored_molecules=anchored_molecules,
            error_message="; ".join(error_messages) if error_messages else None
        )

        self.logger.info("=" * 60)
        if success:
            self.logger.info(f"✓ Anchor Phase 完成: {len(anchored_molecules)}/{len(molecules)} 个分子成功")
        else:
            self.logger.error(f"✗ Anchor Phase 失败: {result.error_message}")
        self.logger.info("=" * 60)

        return result

    def get_sp_energy(self, molecule_name: str) -> Optional[float]:
        """
        获取分子的 SP 能量

        Args:
            molecule_name: 分子名称

        Returns:
            SP 能量，如果未找到返回 None
        """
        return self.anchored_molecules.get(molecule_name, {}).get("e_sp")

    def get_sp_output_xyz(self, molecule_name: str) -> Optional[Path]:
        """
        获取分子的 SP 输出文件路径

        Args:
            molecule_name: 分子名称

        Returns:
            XYZ 文件路径（SP 输出），如果未找到返回 None
        """
        return self.anchored_molecules.get(molecule_name, {}).get("xyz")
