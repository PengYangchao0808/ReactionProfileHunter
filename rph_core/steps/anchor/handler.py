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
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass

from rph_core.utils.log_manager import LoggerMixin
from rph_core.steps.conformer_search.engine import ConformerEngine
from rph_core.utils.molecule_utils import is_small_molecule
from rph_core.utils.small_molecule_cache import SmallMoleculeCache
from rph_core.utils.ui import get_progress_manager

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

        # Small molecule threshold from config (V6.2 Extension)
        self.small_mol_threshold = self.config.get("s1", {}).get(
            "small_molecule_threshold", 10
        )
        self.logger.info(f"Small molecule threshold (heavy atoms): {self.small_mol_threshold}")

        # Removed: SmallMoleculeCache initialization in __init__ (moved to run())
        # This ensures cache is created relative to the actual run directory
        # cache_dir = self.config.get("global", {}).get("small_molecule_cache_dir") ...
        
        self.small_mol_cache = None # Will be initialized in run()

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

        # Dynamic cache initialization relative to current run directory
        cache_dir_cfg = self.config.get("global", {}).get("small_molecule_cache_dir")
        if cache_dir_cfg:
            cache_root = Path(cache_dir_cfg).resolve()
            self.logger.info(f"Using configured small molecule cache: {cache_root}")
        else:
            cache_root = self.base_work_dir / "SmallMolecules"
            self.logger.info(f"Using run-local small molecule cache: {cache_root}")
        self.small_mol_cache = SmallMoleculeCache(cache_root)

        self.logger.info("=" * 60)
        self.logger.info("[S1] Anchor Phase: 统一锚定底物池和产物（v3.0）")
        self.logger.info("=" * 60)
        self.logger.info(f"[S1] 分子数量: {len(molecules)}")
        self._emit_anchor_progress(
            "anchor_started",
            {
                "total": len(molecules),
                "molecules": list(molecules.keys()),
                "work_dir": str(self.base_work_dir),
            },
        )

        pm = get_progress_manager()
        total_mols = len(molecules)
        anchored_molecules: Dict[str, Dict[str, Any]] = {}
        error_messages = []
        molecule_statuses: Dict[str, Dict[str, Any]] = {}

        for idx, (name, smiles) in enumerate(molecules.items()):
            if pm:
                start_progress = int((idx / max(total_mols, 1)) * 100)
                pm.update_step(
                    "s1",
                    completed=min(start_progress, 99),
                    description=f"S1: Anchoring [{idx+1}/{total_mols}] '{name}'"
                )
                pm.set_subtask("S1", "Anchor Molecules", idx + 1, total_mols)
            
            self.logger.info(f"\n[S1] >>> ⚓ 锚定任务启动: {name} ...")
            self.logger.info(f"[S1]     SMILES: {smiles}")

            # Write status file
            status_file = self.base_work_dir / ".rph_step_status.json"
            status_data = {
                "step": "s1",
                "molecule": name,
                "index": idx + 1,
                "total": total_mols,
                "smiles": smiles,
                "status": "running",
                "phase": "anchor_start",
                "molecule_statuses": molecule_statuses,
            }
            with open(status_file, "w") as f:
                json.dump(status_data, f, indent=4)

            molecule_statuses[name] = {
                "status": "running",
                "index": idx + 1,
                "total": total_mols,
                "smiles": smiles,
            }
            self._emit_anchor_progress(
                "molecule_started",
                {
                    "molecule": name,
                    "index": idx + 1,
                    "total": total_mols,
                },
            )

            best_sp_out: Optional[Path] = None
            sp_energy: Optional[float] = None

            try:
                is_small = is_small_molecule(smiles, threshold=self.small_mol_threshold)
                cache_hit = False
                
                theory_signature = self._build_theory_signature()

                if is_small and self.small_mol_cache.exists(smiles, theory_signature):
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
                    lock_file = None
                    if is_small:
                        lock_file = self.small_mol_cache.acquire_compute_lock(smiles)
                        if lock_file is None:
                            self.logger.info(f"    ⏳ Another process is computing {name}, waiting...")
                            import time
                            time.sleep(1)
                            if self.small_mol_cache.exists(smiles, theory_signature):
                                cache_dir = self.small_mol_cache.get_path(smiles)
                                if isinstance(cache_dir, Path):
                                    self.logger.info(f"    ✓ Small molecule {name} found in cache after wait")
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
                        temp_engine = ConformerEngine(
                            config=self.config,
                            work_dir=self.base_work_dir,
                            molecule_name=name
                        )
                        
                        if is_small:
                            self.logger.info(f"    🧪 任务类型: 刚性小分子优化 ({name})")
                            self.logger.info("    ℹ️  说明: 分子较小，跳过构象搜索，直接进行单构象优化。")
                            best_sp_out, sp_energy = temp_engine.run_optimization_only(smiles=smiles)
                        else:
                            self.logger.info(f"    🧬 任务类型: 柔性分子构象搜索与优化 ({name})")
                            self.logger.info("    ℹ️  说明: 执行系综搜索 (CREST) + DFT OPT-SP 耦合循环。")
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
                            self.small_mol_cache.write_cache_meta(smiles, theory_signature)
                            self.logger.info(f"    ✓ Small molecule {name} saved to cache")
                    
                    if lock_file:
                        self.small_mol_cache.release_compute_lock(lock_file)

                self.logger.info(f"    ✓ OPT-SP 完成: {best_sp_out}")
                self.logger.info(f"    ✓ SP 能量: {sp_energy:.8f} Hartree")

                conformer_state = self.base_work_dir / name / "conformer_state.json"

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
                    "qm_output": best_sp_out,
                    "conformer_state": conformer_state if conformer_state.exists() else None,
                }

                molecule_statuses[name] = {
                    "status": "completed",
                    "index": idx + 1,
                    "total": total_mols,
                    "smiles": smiles,
                    "e_sp": sp_energy,
                    "conformer_state": str(conformer_state) if conformer_state.exists() else "",
                }
                self._emit_anchor_progress(
                    "molecule_completed",
                    {
                        "molecule": name,
                        "index": idx + 1,
                        "total": total_mols,
                        "energy_hartree": sp_energy,
                    },
                )
                with open(status_file, "w") as f:
                    json.dump(
                        {
                            "step": "s1",
                            "molecule": name,
                            "index": idx + 1,
                            "total": total_mols,
                            "smiles": smiles,
                            "status": "running",
                            "phase": "anchor_running",
                            "molecule_statuses": molecule_statuses,
                        },
                        f,
                        indent=4,
                    )

                self.logger.info(f"  ✓ {name} 锚定完成")

            except Exception as e:
                error_msg = f"{name} 锚定失败: {e}"
                self.logger.error(error_msg, exc_info=True)
                error_messages.append(error_msg)
                molecule_statuses[name] = {
                    "status": "failed",
                    "index": idx + 1,
                    "total": total_mols,
                    "smiles": smiles,
                    "error": str(e),
                }
                self._emit_anchor_progress(
                    "molecule_failed",
                    {
                        "molecule": name,
                        "index": idx + 1,
                        "total": total_mols,
                        "error": str(e),
                    },
                )

                mol_dir = self.base_work_dir / name
                potential_xyz = list((mol_dir / "dft").glob("*_SP.out"))
                if potential_xyz:
                    anchored_molecules[name] = {
                        "xyz": potential_xyz[0],
                        "e_sp": None,
                        "failed": True
                    }

                with open(status_file, "w") as f:
                    json.dump(
                        {
                            "step": "s1",
                            "molecule": name,
                            "index": idx + 1,
                            "total": total_mols,
                            "smiles": smiles,
                            "status": "running",
                            "phase": "anchor_running",
                            "molecule_statuses": molecule_statuses,
                        },
                        f,
                        indent=4,
                    )

            if pm:
                done_progress = int(((idx + 1) / max(total_mols, 1)) * 100)
                pm.update_step(
                    "s1",
                    completed=min(done_progress, 99),
                    description=f"S1: Anchoring [{idx+1}/{total_mols}] '{name}'"
                )

        success = len(error_messages) < len(molecules)

        self.anchored_molecules = anchored_molecules

        result = AnchorPhaseResult(
            success=success,
            anchored_molecules=anchored_molecules,
            error_message="; ".join(error_messages) if error_messages else None
        )

        self.logger.info("=" * 60)
        self._emit_anchor_progress(
            "anchor_finished",
            {
                "status": "completed" if success else "failed",
                "total": total_mols,
                "completed": sum(1 for item in molecule_statuses.values() if item.get("status") == "completed"),
                "failed": sum(1 for item in molecule_statuses.values() if item.get("status") == "failed"),
            },
        )
        if success:
            self.logger.info(f"[S1] ✓ Anchor Phase 完成: {len(anchored_molecules)}/{len(molecules)} 个分子成功")
        else:
            self.logger.error(f"[S1] ✗ Anchor Phase 失败: {result.error_message}")
        self.logger.info("=" * 60)

        status_file = self.base_work_dir / ".rph_step_status.json"
        with open(status_file, "w") as f:
            json.dump(
                {
                    "step": "s1",
                    "status": "completed" if success else "failed",
                    "total": total_mols,
                    "completed": sum(1 for item in molecule_statuses.values() if item.get("status") == "completed"),
                    "failed": sum(1 for item in molecule_statuses.values() if item.get("status") == "failed"),
                    "molecule_statuses": molecule_statuses,
                },
                f,
                indent=4,
            )

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

    def _build_theory_signature(self) -> Dict[str, Any]:
        theory_opt = self.config.get("theory", {}).get("optimization", {})
        theory_sp = self.config.get("theory", {}).get("single_point", {})
        solvent_config = self.config.get("solvent", {})

        return {
            "opt_method": theory_opt.get("method", "B3LYP"),
            "opt_basis": theory_opt.get("basis", "def2-SVP"),
            "opt_engine": theory_opt.get("engine", "gaussian"),
            "sp_method": theory_sp.get("method", "wB97X-D3BJ"),
            "sp_basis": theory_sp.get("basis", "def2-TZVPP"),
            "sp_engine": theory_sp.get("engine", "orca"),
            "solvent": solvent_config.get("name", "acetone"),
        }

    def _emit_anchor_progress(self, event: str, payload: Dict[str, Any]) -> None:
        message = {
            "schema": "s1_anchor_progress_v1",
            "event": event,
            "payload": payload,
        }
        self.logger.info(
            "S1_ANCHOR_PROGRESS|"
            + json.dumps(self._json_safe(message), sort_keys=True, ensure_ascii=True, default=str, allow_nan=False)
        )

    def _json_safe(self, value: Any) -> Any:
        if isinstance(value, float) and not math.isfinite(value):
            return None
        if isinstance(value, dict):
            return {str(k): self._json_safe(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [self._json_safe(item) for item in value]
        return value
