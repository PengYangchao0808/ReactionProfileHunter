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
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass

from rph_core.utils.log_manager import LoggerMixin
from rph_core.steps.conformer_search.engine import ConformerEngine
from rph_core.utils.molecule_utils import is_small_molecule
from rph_core.utils.small_molecule_cache import SmallMoleculeCache

logger = logging.getLogger(__name__)


@dataclass
class AnchorPhaseResult:
    """锚定阶段结果（新结构）"""
    success: bool
    anchored_molecules: Dict[str, Dict[str, Any]]
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

    def __init__(self, config: Dict[str, Any], base_work_dir: Optional[Path] = None):
        """
        初始化 AnchorPhase - v3.0

        Args:
            config: 配置字典
            base_work_dir: 基础工作目录（例如 S1_ConfGeneration/）
        """
        self.config = config
        self.base_work_dir = Path(base_work_dir) if base_work_dir else Path.cwd()
        self.base_work_dir.mkdir(parents=True, exist_ok=True)

        # Initialize SmallMoleculeCache
        cache_dir = self.config.get("global", {}).get("small_molecule_cache_dir")
        if cache_dir:
            self.cache_root = Path(cache_dir).resolve()
        else:
            self.cache_root = self.base_work_dir.parent.resolve() / "SmallMolecules"
        
        self.small_mol_cache = SmallMoleculeCache(self.cache_root)

        # 结果缓存
        self.anchored_molecules: Dict[str, Dict[str, Any]] = {}

        # 初始化 ConformerEngine（但不在 run 方法中创建实例）
        # ConformerEngine 将在内部创建分子自治目录结构
        # 例如：S1_ConfGeneration/[Molecule_Name]/crest 和 .../dft

        self.logger.info("AnchorPhase v3.0 初始化完成（分子自治架构）")

    def run(
        self,
        molecules: Dict[str, str],
        base_work_dir: Optional[Path] = None
    ) -> AnchorPhaseResult:
        """
        执行统一锚定流程（新结构）

        Args:
            molecules: 分子字典 {名称: SMILES}
            base_work_dir: 可选，覆盖初始化的基础工作目录

        Returns:
            AnchorPhaseResult 对象
        """
        if base_work_dir:
            self.base_work_dir = Path(base_work_dir)
            self.base_work_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("=" * 60)
        self.logger.info("Anchor Phase: 统一锚定底物池和产物（v3.0）")
        self.logger.info("=" * 60)
        self.logger.info(f"分子数量: {len(molecules)}")

        anchored_molecules: Dict[str, Dict[str, Any]] = {}
        error_messages = []

        for name, smiles in molecules.items():
            self.logger.info(f"\n>>> 锚定 {name}...")
            self.logger.info(f"    SMILES: {smiles}")

            best_sp_out: Optional[Path] = None
            sp_energy: Optional[float] = None

            try:
                is_small = is_small_molecule(smiles)
                cache_hit = False

                if is_small and self.small_mol_cache.exists(smiles):
                    cache_dir = self.small_mol_cache.get_path(smiles)
                    if isinstance(cache_dir, Path):
                        self.logger.info(f"    ✓ Small molecule {name} found in cache, skipping")
                        local_mol_dir = self.base_work_dir / name
                        local_mol_dir.mkdir(parents=True, exist_ok=True)
                        local_dft_dir = local_mol_dir / "dft"
                        local_dft_dir.mkdir(parents=True, exist_ok=True)

                        cached_min_xyz = cache_dir / "molecule_min.xyz"
                        local_min_xyz = local_mol_dir / f"{name}_global_min.xyz"
                        shutil.copy(cached_min_xyz, local_min_xyz)

                        cached_dft = cache_dir / "dft"
                        if cached_dft.exists():
                            for f in cached_dft.glob("*"):
                                if f.is_file():
                                    shutil.copy(f, local_dft_dir / f.name)
                        
                        best_sp_out = local_min_xyz
                        with open(best_sp_out, "r") as f:
                            lines = f.readlines()
                            comment = lines[1] if len(lines) > 1 else ""
                            import re
                            match = re.search(r"E=([-+]?\d*\.\d+|\d+)", comment)
                            if match:
                                sp_energy = float(match.group(1))
                            else:
                                try:
                                    sp_energy = float(comment.split()[0])
                                except (ValueError, IndexError):
                                    sp_energy = 0.0
                        
                        cache_hit = True

                if not cache_hit:
                    self.logger.info("    步骤 1/2: 系综搜索 + DFT OPT-SP 耦合...")
                    temp_engine = ConformerEngine(
                        config=self.config,
                        work_dir=self.base_work_dir,
                        molecule_name=name
                    )
                    best_sp_out, sp_energy = temp_engine.run(smiles=smiles)

                    if is_small:
                        cache_dir = self.small_mol_cache.get_or_create(smiles, name=name)
                        shutil.copy(best_sp_out, cache_dir / "molecule_min.xyz")
                        mol_dft_dir = self.base_work_dir / name / "dft"
                        if mol_dft_dir.exists():
                            dest_dft = cache_dir / "dft"
                            if dest_dft.exists():
                                shutil.rmtree(dest_dft)
                            shutil.copytree(mol_dft_dir, dest_dft)
                        self.logger.info(f"    ✓ Small molecule {name} saved to cache")

                self.logger.info(f"    ✓ OPT-SP 完成: {best_sp_out}")
                self.logger.info(f"    ✓ SP 能量: {sp_energy:.8f} Hartree")

                mol_dir = self.base_work_dir / name / "dft"
                log_file = None
                chk_file = None
                fchk_file = None

                if best_sp_out and best_sp_out.suffix == ".log":
                    log_file = best_sp_out
                elif best_sp_out and best_sp_out.suffix == ".out":
                    log_file = best_sp_out
                else:
                    potential_logs = list(mol_dir.glob("*.log")) + list(mol_dir.glob("*.out"))
                    if potential_logs:
                        log_file = potential_logs[0]

                if best_sp_out and best_sp_out.stem:
                    potential_chk = (
                        list(mol_dir.glob(f"{best_sp_out.stem}.chk"))
                        + list(mol_dir.glob("*.chk"))
                    )
                    if potential_chk:
                        chk_file = potential_chk[0]
                        from rph_core.utils.qc_interface import try_formchk
                        fchk_file = try_formchk(chk_file)

                anchored_molecules[name] = {
                    "xyz": best_sp_out,
                    "e_sp": sp_energy,
                    "log": log_file,
                    "chk": chk_file,
                    "fchk": fchk_file,
                    "qm_output": best_sp_out
                }

                self.logger.info(f"  ✓ {name} 锚定完成")

            except Exception as e:
                error_msg = f"{name} 锚定失败: {e}"
                self.logger.error(error_msg, exc_info=True)
                error_messages.append(error_msg)

                mol_dir = self.base_work_dir / name
                potential_xyz = list((mol_dir / "dft").glob("*_SP.out"))
                if potential_xyz:
                    anchored_molecules[name] = {
                        "xyz": potential_xyz[0],
                        "e_sp": None,
                        "failed": True
                    }

        success = len(error_messages) < len(molecules)

        self.anchored_molecules = anchored_molecules

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
