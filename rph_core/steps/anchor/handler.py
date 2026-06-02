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
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import asdict, dataclass

from rph_core.utils.log_manager import LoggerMixin
from rph_core.steps.anchor.engine_factory import create_s1_engine
from rph_core.steps.conformer_search.protocols import resolve_protocol_spec
from rph_core.utils.intra_reaction_scheduler import IntraReactionScheduler, ParallelQCJob
from rph_core.utils.molecule_utils import is_small_molecule
from rph_core.utils.small_molecule_cache import SmallMoleculeCache
from rph_core.utils.ui import get_progress_manager
from rph_core.version import __version__

logger = logging.getLogger(__name__)


@dataclass
class AnchorPhaseResult:
    """锚定阶段结果（新结构）"""
    success: bool
    anchored_molecules: Dict[str, Dict[str, Any]]
    # 格式: {"reactant_A": {"xyz": Path, "e_sp": float}, ...}
    error_message: Optional[str] = None


@dataclass
class MoleculeAnchorOutcome:
    """Thread-safe return payload for one molecule-level anchor job."""

    name: str
    smiles: str
    index: int
    total: int
    status: str
    anchored_data: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


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
        step1_cfg = self.config.get("step1", {})
        s1_legacy_cfg = self.config.get("s1", {})
        self.protocol = str(step1_cfg.get("protocol", s1_legacy_cfg.get("protocol", "ext"))).lower()
        self.protocol_spec = resolve_protocol_spec(self.config, self.protocol)
        self.small_mol_threshold = step1_cfg.get(
            "small_molecule_threshold", s1_legacy_cfg.get("small_molecule_threshold", 10)
        )
        self.logger.info(f"[S1] Protocol tier: {self.protocol}")
        self.logger.info(f"[S1] Small-molecule cutoff: {self.small_mol_threshold} heavy atoms")

        # Removed: SmallMoleculeCache initialization in __init__ (moved to run())
        # This ensures cache is created relative to the actual run directory
        # cache_dir = self.config.get("global", {}).get("small_molecule_cache_dir") ...
        
        self.small_mol_cache = None # Will be initialized in run()

        # 结果缓存
        self.anchored_molecules: Dict[str, Dict[str, Any]] = {}

        # 初始化 ConformerEngine（但不在 run 方法中创建实例）
        # ConformerEngine 将在内部创建分子自治目录结构
        # 例如：S1_ConfGeneration/[Molecule_Name]/crest 和 .../dft

        self.logger.info("[S1] AnchorPhase ready (V4.0 CENSO-like model)")

    def _display_path(self, path: Path) -> str:
        path = Path(path)
        for root in (self.base_work_dir, Path.cwd()):
            try:
                rel = path.resolve().relative_to(root.resolve())
                return str(rel)
            except Exception:
                continue
        parts = path.parts
        if len(parts) >= 3:
            return str(Path("...") / parts[-2] / parts[-1])
        return path.name or str(path)

    def _protocol_tier_label(self) -> str:
        mapping = {
            "ext": "RPH baseline",
            "default": "RPH baseline",
            "full": "CENSO-like full funnel",
            "lite": "CENSO-like lite funnel",
            "zero": "CENSO-like zero funnel",
        }
        return mapping.get(self.protocol_spec.name, self.protocol_spec.name)

    def _protocol_funnel_summary(self) -> str:
        mapping = {
            "ext": "GFN0→GFN2 two-stage ensemble → all candidates",
            "default": "GFN0→GFN2 two-stage ensemble → all candidates",
            "full": "GFN2 ensemble → prescreen fast-SP → screening fast-SP → survivor window",
            "lite": "GFN2 ensemble → screening-level fast-SP → cutoff",
            "zero": "GFN2 ensemble → narrow energy window",
        }
        return mapping.get(self.protocol_spec.name, "protocol-defined funnel")

    def _protocol_handoff_summary(self) -> str:
        return "shared DFT handoff → Gaussian OPT/FREQ → ORCA final SP"

    def run_single_molecule(
        self,
        name: str,
        smiles: str,
        base_work_dir: Optional[Path] = None,
    ) -> AnchorPhaseResult:
        return self.run({name: smiles}, base_work_dir=base_work_dir)

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
            self.logger.info(f"[S1] Cache root   : {self._display_path(cache_root)}")
        else:
            cache_root = self.base_work_dir / "SmallMolecules"
            self.logger.info(f"[S1] Cache root   : {self._display_path(cache_root)}")
        self.small_mol_cache = SmallMoleculeCache(cache_root)

        self.logger.info("=" * 72)
        self.logger.info("[S1] Conformer Search & Anchor Phase — V4.0 CENSO-like model")
        self.logger.info(f"[S1] Tier         : {self.protocol_spec.name.upper()} · {self._protocol_tier_label()}")
        self.logger.info(f"[S1] Funnel       : {self._protocol_funnel_summary()}")
        self.logger.info(f"[S1] Handoff      : {self._protocol_handoff_summary()}")
        self.logger.info(f"[S1] Molecules    : {len(molecules)}")
        self.logger.info("=" * 72)
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
        error_messages: List[str] = []
        molecule_statuses: Dict[str, Dict[str, Any]] = {
            name: {
                "status": "queued",
                "index": idx + 1,
                "total": total_mols,
                "smiles": smiles,
            }
            for idx, (name, smiles) in enumerate(molecules.items())
        }
        status_file = self.base_work_dir / ".rph_step_status.json"
        self._write_anchor_status_file(
            status_file,
            status="running",
            phase="anchor_queued",
            molecule_statuses=molecule_statuses,
            total=total_mols,
        )

        outcomes = self._run_molecule_anchor_jobs(
            molecules=molecules,
            molecule_statuses=molecule_statuses,
            status_file=status_file,
            pm=pm,
        )

        for outcome in outcomes:
            name = outcome.name
            if outcome.status == "completed" and outcome.anchored_data is not None:
                anchored_molecules[name] = outcome.anchored_data
                conformer_state = outcome.anchored_data.get("conformer_state")
                molecule_statuses[name] = {
                    "status": "completed",
                    "index": outcome.index,
                    "total": outcome.total,
                    "smiles": outcome.smiles,
                    "e_sp": outcome.anchored_data.get("e_sp"),
                    "conformer_state": str(conformer_state) if isinstance(conformer_state, Path) else "",
                }
            else:
                error_msg = outcome.error or f"{name} 锚定失败: unknown error"
                error_messages.append(error_msg)
                molecule_statuses[name] = {
                    "status": "failed",
                    "index": outcome.index,
                    "total": outcome.total,
                    "smiles": outcome.smiles,
                    "error": error_msg,
                }
                if outcome.anchored_data is not None:
                    anchored_molecules[name] = outcome.anchored_data

            completed_count = sum(1 for item in molecule_statuses.values() if item.get("status") == "completed")
            failed_count = sum(1 for item in molecule_statuses.values() if item.get("status") == "failed")
            self._write_anchor_status_file(
                status_file,
                status="running",
                phase="anchor_running",
                molecule_statuses=molecule_statuses,
                current=outcome,
                total=total_mols,
            )
            if pm:
                done_progress = int(((completed_count + failed_count) / max(total_mols, 1)) * 100)
                pm.update_step(
                    "s1",
                    completed=min(done_progress, 99),
                    description=f"S1: Anchored {completed_count}/{total_mols} molecules ({failed_count} failed)",
                )

        success = len(error_messages) < len(molecules)

        self.anchored_molecules = anchored_molecules
        self._write_s1_provenance(anchored_molecules=anchored_molecules, molecule_statuses=molecule_statuses)

        result = AnchorPhaseResult(
            success=success,
            anchored_molecules=anchored_molecules,
            error_message="; ".join(error_messages) if error_messages else None
        )

        self.logger.info("=" * 72)
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
        self.logger.info("=" * 72)

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

    def _run_molecule_anchor_jobs(
        self,
        *,
        molecules: Dict[str, str],
        molecule_statuses: Dict[str, Dict[str, Any]],
        status_file: Path,
        pm: Any,
    ) -> List[MoleculeAnchorOutcome]:
        total_mols = len(molecules)
        items = [(idx, name, smiles) for idx, (name, smiles) in enumerate(molecules.items())]
        if not items:
            return []

        scheduler = IntraReactionScheduler(self.config)
        max_workers = min(total_mols, max(1, scheduler.s1_max_molecule_workers))
        parallel_enabled = scheduler.s1_molecule_parallel and max_workers > 1

        if not parallel_enabled:
            self.logger.info("[S1] Molecule-level parallel: disabled; running anchors sequentially")
            outcomes: List[MoleculeAnchorOutcome] = []
            for idx, name, smiles in items:
                molecule_statuses[name] = {
                    "status": "running",
                    "index": idx + 1,
                    "total": total_mols,
                    "smiles": smiles,
                }
                self._write_anchor_status_file(
                    status_file,
                    status="running",
                    phase="anchor_running",
                    molecule_statuses=molecule_statuses,
                    total=total_mols,
                    molecule=name,
                    index=idx + 1,
                    smiles=smiles,
                )
                if pm:
                    start_progress = int((idx / max(total_mols, 1)) * 100)
                    pm.update_step(
                        "s1",
                        completed=min(start_progress, 99),
                        description=f"S1: Anchoring [{idx+1}/{total_mols}] '{name}'",
                    )
                    pm.set_subtask("S1", "Anchor Molecules", idx + 1, total_mols)
                outcomes.append(self._anchor_single_molecule(name, smiles, idx, total_mols))
            return outcomes

        lane = scheduler.get_s1_molecule_lane()
        self.logger.info(
            "[S1] Molecule-level parallel: enabled (%s molecules, max_workers=%s, crest_cores/job=%s)",
            total_mols,
            max_workers,
            lane.nproc,
        )
        self._emit_anchor_progress(
            "molecule_parallel_started",
            {"total": total_mols, "max_workers": max_workers, "lane_nproc": lane.nproc},
        )
        if pm:
            pm.update_step(
                "s1",
                completed=0,
                description=f"S1: Anchoring {total_mols} molecules in parallel ({max_workers} workers)",
            )

        outcomes = []
        for batch_start in range(0, total_mols, max_workers):
            batch = items[batch_start:batch_start + max_workers]
            batch_ids = [name for _, name, _ in batch]
            self.logger.info(
                "[S1] Molecule parallel batch %s-%s/%s: %s",
                batch_start + 1,
                batch_start + len(batch),
                total_mols,
                ", ".join(batch_ids),
            )
            for idx, name, smiles in batch:
                molecule_statuses[name] = {
                    "status": "running",
                    "index": idx + 1,
                    "total": total_mols,
                    "smiles": smiles,
                }
            self._write_anchor_status_file(
                status_file,
                status="running",
                phase="anchor_parallel_batch",
                molecule_statuses=molecule_statuses,
                total=total_mols,
            )

            jobs = [
                ParallelQCJob(
                    job_id=name,
                    lane=lane,
                    func=self._anchor_single_molecule,
                    args=(name, smiles, idx, total_mols),
                    description=f"S1 anchor for {name}",
                )
                for idx, name, smiles in batch
            ]
            lookup = {name: (idx, smiles) for idx, name, smiles in batch}
            for result in scheduler.run_parallel(jobs):
                if result.success and isinstance(result.result, MoleculeAnchorOutcome):
                    outcomes.append(result.result)
                    continue
                idx, smiles = lookup.get(result.job_id, (-1, ""))
                error = result.error or "molecule_parallel_job_failed"
                outcomes.append(
                    MoleculeAnchorOutcome(
                        name=result.job_id,
                        smiles=smiles,
                        index=idx + 1 if idx >= 0 else 0,
                        total=total_mols,
                        status="failed",
                        anchored_data=self._find_failed_molecule_output(result.job_id),
                        error=f"{result.job_id} 锚定失败: {error}",
                    )
                )
        return outcomes

    def _anchor_single_molecule(self, name: str, smiles: str, idx: int, total_mols: int) -> MoleculeAnchorOutcome:
        self.logger.info(f"\n[S1] >>> Molecule {idx + 1}/{total_mols} · {name}")
        self.logger.info(f"[S1]     SMILES      : {smiles}")
        self.logger.info(f"[S1]     Search tier : {self._protocol_tier_label()}")
        self._emit_anchor_progress(
            "molecule_started",
            {"molecule": name, "index": idx + 1, "total": total_mols},
        )

        best_sp_out: Optional[Path] = None
        sp_energy: Optional[float] = None
        final_opt_sp_stage_meta: Dict[str, Any] = self._default_final_opt_sp_meta(stage_status="not_run")
        lock_file: Optional[Path] = None

        try:
            if self.small_mol_cache is None:
                raise RuntimeError("SmallMoleculeCache was not initialized")

            is_small = is_small_molecule(smiles, threshold=self.small_mol_threshold)
            cache_hit = False
            theory_signature = self._build_theory_signature()

            if is_small and self.small_mol_cache.exists(smiles, theory_signature):
                best_sp_out, sp_energy = self._reuse_small_molecule_cache(name, smiles)
                cache_hit = True
                final_opt_sp_stage_meta = self._default_final_opt_sp_meta(stage_status="cache_reused")

            if not cache_hit:
                if is_small:
                    lock_file = self.small_mol_cache.acquire_compute_lock(smiles)
                    if lock_file is None:
                        self.logger.info(f"    ⏳ Small-molecule cache busy for {name}; checking completed cache")
                        if self.small_mol_cache.exists(smiles, theory_signature):
                            best_sp_out, sp_energy = self._reuse_small_molecule_cache(name, smiles)
                            cache_hit = True

                if not cache_hit:
                    temp_engine = create_s1_engine(
                        protocol=self.protocol,
                        config=self.config,
                        work_dir=self.base_work_dir,
                        molecule_name=name,
                    )

                    if is_small:
                        self.logger.info(f"    🧪 任务类型: 刚性小分子优化 ({name})")
                        self.logger.info("    ℹ️  说明: 分子较小，跳过构象搜索，直接进行单构象优化。")
                        best_sp_out, sp_energy = temp_engine.run_optimization_only(smiles=smiles)
                    else:
                        self.logger.info(f"    🧬 任务类型: 柔性分子构象搜索与优化 ({name})")
                        self.logger.info("    ℹ️  说明: 执行系综搜索 (CREST) + DFT OPT-SP 耦合循环。")
                        best_sp_out, sp_energy = temp_engine.run(smiles=smiles)

                    stage_meta_from_engine = getattr(temp_engine, "last_final_opt_sp_meta", None)
                    if isinstance(stage_meta_from_engine, dict) and stage_meta_from_engine:
                        final_opt_sp_stage_meta = {
                            **self._default_final_opt_sp_meta(stage_status="completed"),
                            **stage_meta_from_engine,
                        }
                    else:
                        final_opt_sp_stage_meta = self._default_final_opt_sp_meta(stage_status="completed")

                    if is_small:
                        self._store_small_molecule_cache(name, smiles, best_sp_out, theory_signature)

            if best_sp_out is None or sp_energy is None:
                raise RuntimeError("S1 engine did not produce a final structure and SP energy")

            self.logger.info(f"    ✓ Finalized output: {self._display_path(best_sp_out)}")
            self.logger.info(f"    ✓ SP 能量: {sp_energy:.8f} Hartree")
            anchored_data = self._build_anchored_molecule_record(
                name=name,
                best_sp_out=best_sp_out,
                sp_energy=sp_energy,
                final_opt_sp_stage_meta=final_opt_sp_stage_meta,
            )
            self._emit_anchor_progress(
                "molecule_completed",
                {
                    "molecule": name,
                    "index": idx + 1,
                    "total": total_mols,
                    "energy_hartree": sp_energy,
                },
            )
            self.logger.info(f"  ✓ {name} anchor complete")
            return MoleculeAnchorOutcome(
                name=name,
                smiles=smiles,
                index=idx + 1,
                total=total_mols,
                status="completed",
                anchored_data=anchored_data,
            )
        except Exception as exc:
            error_msg = f"{name} 锚定失败: {exc}"
            self.logger.error(error_msg, exc_info=True)
            self._emit_anchor_progress(
                "molecule_failed",
                {"molecule": name, "index": idx + 1, "total": total_mols, "error": str(exc)},
            )
            return MoleculeAnchorOutcome(
                name=name,
                smiles=smiles,
                index=idx + 1,
                total=total_mols,
                status="failed",
                anchored_data=self._find_failed_molecule_output(name),
                error=error_msg,
            )
        finally:
            if lock_file is not None and self.small_mol_cache is not None:
                self.small_mol_cache.release_compute_lock(lock_file)

    def _reuse_small_molecule_cache(self, name: str, smiles: str) -> Tuple[Path, float]:
        if self.small_mol_cache is None:
            raise RuntimeError("SmallMoleculeCache was not initialized")
        cache_dir = self.small_mol_cache.get_path(smiles)
        if not isinstance(cache_dir, Path):
            raise RuntimeError(f"Invalid small-molecule cache path for {name}")
        self.logger.info(f"    ✓ Small molecule {name} found in cache, skipping")
        local_mol_dir = self.base_work_dir / name
        local_mol_dir.mkdir(parents=True, exist_ok=True)
        local_dft_dir = local_mol_dir / "finalDFT"
        local_dft_dir.mkdir(parents=True, exist_ok=True)

        cached_min_xyz = cache_dir / "molecule_min.xyz"
        local_min_xyz = local_mol_dir / f"{name}_global_min.xyz"
        shutil.copy(cached_min_xyz, local_min_xyz)

        cached_dft = cache_dir / "finalDFT"
        if not cached_dft.exists():
            cached_dft = cache_dir / "dft"
        if cached_dft.exists():
            for item in cached_dft.glob("*"):
                if item.is_file():
                    shutil.copy(item, local_dft_dir / item.name)
        return local_min_xyz, self._extract_xyz_energy(local_min_xyz)

    def _store_small_molecule_cache(
        self,
        name: str,
        smiles: str,
        best_sp_out: Path,
        theory_signature: Dict[str, Any],
    ) -> None:
        if self.small_mol_cache is None:
            raise RuntimeError("SmallMoleculeCache was not initialized")
        cache_dir = self.small_mol_cache.get_or_create(smiles, name=name)
        shutil.copy(best_sp_out, cache_dir / "molecule_min.xyz")
        mol_dft_dir = self.base_work_dir / name / "finalDFT"
        if not mol_dft_dir.exists():
            mol_dft_dir = self.base_work_dir / name / "dft"
        if mol_dft_dir.exists():
            dest_dft = cache_dir / "finalDFT"
            if dest_dft.exists():
                shutil.rmtree(dest_dft)
            shutil.copytree(mol_dft_dir, dest_dft)
        self.small_mol_cache.write_cache_meta(smiles, theory_signature)
        self._ensure_cache_thermo_json(cache_dir)
        self.logger.info(f"    ✓ Cached small molecule: {name}")

    def _ensure_cache_thermo_json(self, cache_dir: Path) -> None:
        from rph_core.utils.thermo import ensure_thermo_json_from_entry

        _ = ensure_thermo_json_from_entry(cache_dir)

    def _extract_xyz_energy(self, xyz_path: Path) -> float:
        import re

        with open(xyz_path, "r") as handle:
            lines = handle.readlines()
        comment = lines[1] if len(lines) > 1 else ""
        match = re.search(r"E=([-+]?\d*\.\d+|\d+)", comment)
        if match:
            return float(match.group(1))
        try:
            return float(comment.split()[0])
        except (ValueError, IndexError):
            return 0.0

    def _build_anchored_molecule_record(
        self,
        *,
        name: str,
        best_sp_out: Path,
        sp_energy: float,
        final_opt_sp_stage_meta: Dict[str, Any],
    ) -> Dict[str, Any]:
        conformer_state = self.base_work_dir / name / "conformer_state.json"
        mol_dir = self.base_work_dir / name / "finalDFT"
        if not mol_dir.exists():
            mol_dir = self.base_work_dir / name / "dft"

        log_file = None
        chk_file = None
        fchk_file = None
        if best_sp_out.suffix in {".log", ".out"}:
            log_file = best_sp_out
        else:
            potential_logs = list(mol_dir.glob("*.log")) + list(mol_dir.glob("*.out"))
            if potential_logs:
                log_file = potential_logs[0]

        if best_sp_out.stem:
            potential_chk = list(mol_dir.glob(f"{best_sp_out.stem}.chk")) + list(mol_dir.glob("*.chk"))
            if potential_chk:
                chk_file = potential_chk[0]
                from rph_core.utils.qc_interface import try_formchk

                fchk_file = try_formchk(chk_file)

        return {
            "xyz": best_sp_out,
            "e_sp": sp_energy,
            "log": log_file,
            "chk": chk_file,
            "fchk": fchk_file,
            "qm_output": best_sp_out,
            "conformer_state": conformer_state if conformer_state.exists() else None,
            "protocol": self.protocol,
            "engine": "ConformerEngine",
            "final_opt_sp_stage_meta": final_opt_sp_stage_meta,
        }

    def _find_failed_molecule_output(self, name: str) -> Optional[Dict[str, Any]]:
        mol_dir = self.base_work_dir / name
        potential_xyz = list((mol_dir / "finalDFT").glob("*_SP.out"))
        if not potential_xyz:
            potential_xyz = list((mol_dir / "dft").glob("*_SP.out"))
        if not potential_xyz:
            return None
        return {"xyz": potential_xyz[0], "e_sp": None, "failed": True}

    def _write_anchor_status_file(
        self,
        status_file: Path,
        *,
        status: str,
        phase: str,
        molecule_statuses: Dict[str, Dict[str, Any]],
        total: int,
        current: Optional[MoleculeAnchorOutcome] = None,
        molecule: Optional[str] = None,
        index: Optional[int] = None,
        smiles: Optional[str] = None,
    ) -> None:
        payload: Dict[str, Any] = {
            "step": "s1",
            "status": status,
            "phase": phase,
            "total": total,
            "completed": sum(1 for item in molecule_statuses.values() if item.get("status") == "completed"),
            "failed": sum(1 for item in molecule_statuses.values() if item.get("status") == "failed"),
            "molecule_statuses": molecule_statuses,
        }
        if current is not None:
            payload.update(
                {
                    "molecule": current.name,
                    "index": current.index,
                    "smiles": current.smiles,
                }
            )
        elif molecule is not None:
            payload.update({"molecule": molecule, "index": index, "smiles": smiles})
        status_file.write_text(json.dumps(self._json_safe(payload), indent=4), encoding="utf-8")

    def _write_s1_provenance(
        self,
        *,
        anchored_molecules: Dict[str, Dict[str, Any]],
        molecule_statuses: Dict[str, Dict[str, Any]],
    ) -> None:
        theory_signature = self._build_theory_signature()
        final_opt_sp_meta = self._resolve_final_opt_sp_meta(anchored_molecules)
        final_opt_enabled = bool(final_opt_sp_meta.get("enabled", self.protocol_spec.final_opt_sp_enabled))
        funnel_payload = {
            "search_mode": self.protocol_spec.funnel_policy.search_mode,
            "clustering_mode": self.protocol_spec.funnel_policy.clustering_mode,
            "prescreen_mode": self.protocol_spec.funnel_policy.prescreen_mode,
            "rerank_mode": self.protocol_spec.funnel_policy.rerank_mode,
            "use_mrrho_like_correction": self.protocol_spec.funnel_policy.use_mrrho_like_correction,
            "survivor_window_kcal": self.protocol_spec.funnel_policy.survivor_window_kcal,
            "narrow_window_kcal": self.protocol_spec.funnel_policy.narrow_window_kcal,
            "prescreen_window_kcal": self.protocol_spec.funnel_policy.prescreen_window_kcal,
            "screening_window_kcal": self.protocol_spec.funnel_policy.screening_window_kcal,
            "optimize_limit": self.protocol_spec.funnel_policy.optimize_limit,
            "top2_fallback_enabled": self.protocol_spec.funnel_policy.top2_fallback_enabled,
            "boltzmann_cutoff": self.protocol_spec.funnel_policy.boltzmann_cutoff,
            "ranking_basis": str(final_opt_sp_meta.get("ranking_basis_before_handoff", "input_order")),
            "candidate_total": int(final_opt_sp_meta.get("input_candidate_count", 0)),
            "survivor_count": int(final_opt_sp_meta.get("selected_candidate_count", 0)),
            "candidate_ids": list(final_opt_sp_meta.get("input_candidate_ids", [])),
            "window_threshold": final_opt_sp_meta.get("energy_window_used"),
            "energy_window_kcal": final_opt_sp_meta.get("energy_window_used"),
            "window_count": int(final_opt_sp_meta.get("window_count", 0)),
            "gap_rank1_rank2": final_opt_sp_meta.get("gap_rank1_rank2"),
            "stages_executed": list(final_opt_sp_meta.get("stages_executed", [])),
            "approx_thermo_applied": bool(final_opt_sp_meta.get("approx_thermo_applied", False)),
            "fallback_triggered": bool(final_opt_sp_meta.get("fallback_triggered", False)),
        }
        handoff_payload = {
            "mode": str(final_opt_sp_meta.get("mode_effective", self.protocol_spec.handoff_policy.mode)),
            "mode_requested": str(final_opt_sp_meta.get("mode_requested", self.protocol_spec.handoff_policy.mode)),
            "mode_effective": str(final_opt_sp_meta.get("mode_effective", self.protocol_spec.handoff_policy.mode)),
            "fallback_mode": self.protocol_spec.handoff_policy.fallback_mode,
            "small_gap_kcal": self.protocol_spec.handoff_policy.small_gap_kcal,
            "ranking_after_handoff": self.protocol_spec.handoff_policy.ranking_after_handoff,
            "selected_candidate_count": int(final_opt_sp_meta.get("selected_candidate_count", 0)),
            "selected_candidate_ids": list(final_opt_sp_meta.get("selected_candidate_ids", [])),
            "selection_rule": str(final_opt_sp_meta.get("selection_mode", self.protocol_spec.handoff_policy.mode)),
            "selection_reason": str(final_opt_sp_meta.get("selection_reason", "unspecified")),
            "candidate_scores_before_handoff": final_opt_sp_meta.get("candidate_scores_before_handoff", {}),
            "window_count": int(final_opt_sp_meta.get("window_count", 0)),
            "gap_rank1_rank2": final_opt_sp_meta.get("gap_rank1_rank2"),
            "fallback_trigger": bool(final_opt_sp_meta.get("fallback_trigger", final_opt_sp_meta.get("fallback_triggered", False))),
            "fallback_triggered": bool(final_opt_sp_meta.get("fallback_triggered", False)),
        }
        payload: Dict[str, Any] = {
            "schema_version": "s1_provenance_v1",
            "protocol_spec_version": "protocol_spec_v1",
            "rph_version": __version__,
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "protocol": self.protocol_spec.name,
            "engine": "ConformerEngine",
            "has_geometry_optimization": True,
            "final_opt_sp": {
                "enabled": final_opt_enabled,
                "stage": "final_opt_sp" if final_opt_enabled else "disabled",
                "stage_status": str(final_opt_sp_meta.get("stage_status", "unknown")),
                "freq_requested": bool(final_opt_sp_meta.get("freq_requested", self.protocol_spec.freq_enabled)),
                "final_sp_requested": bool(final_opt_sp_meta.get("final_sp_requested", self.protocol_spec.final_sp_enabled)),
                "selection_mode": str(final_opt_sp_meta.get("selection_mode", self.protocol_spec.selection_mode)),
                "protocol": str(final_opt_sp_meta.get("protocol", self.protocol_spec.name)),
                "selected_candidate_count": int(final_opt_sp_meta.get("selected_candidate_count", 0)),
                "freq_enabled": bool(final_opt_sp_meta.get("freq_requested", self.protocol_spec.freq_enabled)),
                "final_sp_enabled": bool(final_opt_sp_meta.get("final_sp_requested", self.protocol_spec.final_sp_enabled)),
                "stage_meta": dict(final_opt_sp_meta),
            },
            "protocol_spec": {
                "name": self.protocol_spec.name,
                "two_stage_enabled": self.protocol_spec.two_stage_enabled,
                "ngeom_default": self.protocol_spec.ngeom_default,
                "ngeom_max": self.protocol_spec.ngeom_max,
                "funnel_policy": asdict(self.protocol_spec.funnel_policy),
                "handoff_policy": asdict(self.protocol_spec.handoff_policy),
                "provenance_flags": self.protocol_spec.provenance_flags,
            },
            "funnel": funnel_payload,
            "handoff": handoff_payload,
            "artifact_contract": {
                "product_geometry_level": "dft_optimized",
                "has_geometry_optimization": True,
            },
            "theory_signature": theory_signature,
            "conformer_search": {
                "two_stage_enabled": self.protocol_spec.two_stage_enabled,
                "ngeom_default": self.protocol_spec.ngeom_default,
                "ngeom_max": self.protocol_spec.ngeom_max,
            },
            "molecules": {},
        }

        for name, record in anchored_molecules.items():
            xyz_path = record.get("xyz")
            payload["molecules"][name] = {
                "status": molecule_statuses.get(name, {}).get("status", "unknown"),
                "smiles": molecule_statuses.get(name, {}).get("smiles", ""),
                "energy_hartree": record.get("e_sp"),
                "xyz": str(xyz_path) if isinstance(xyz_path, Path) else None,
            }

        provenance_file = self.base_work_dir / "provenance.json"
        provenance_file.write_text(
            json.dumps(self._json_safe(payload), indent=2, ensure_ascii=False, sort_keys=True),
            encoding="utf-8",
        )

    def _default_final_opt_sp_meta(self, *, stage_status: str) -> Dict[str, Any]:
        enabled = bool(self.protocol_spec.final_opt_sp_enabled)
        return {
            "stage": "final_opt_sp" if enabled else "disabled",
            "stage_status": stage_status,
            "enabled": enabled,
            "protocol": self.protocol_spec.name,
            "selection_mode": self.protocol_spec.selection_mode,
            "mode_requested": self.protocol_spec.handoff_policy.mode,
            "mode_effective": self.protocol_spec.handoff_policy.mode,
            "selected_candidate_count": 0,
            "freq_requested": bool(self.protocol_spec.freq_enabled),
            "final_sp_requested": bool(self.protocol_spec.final_sp_enabled),
            "fallback_trigger": False,
            "fallback_triggered": False,
            "window_count": 0,
            "gap_rank1_rank2": None,
        }

    def _resolve_final_opt_sp_meta(self, anchored_molecules: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        resolved = self._default_final_opt_sp_meta(stage_status="unknown")
        for record in anchored_molecules.values():
            stage_meta = record.get("final_opt_sp_stage_meta")
            if isinstance(stage_meta, dict) and stage_meta:
                resolved = {**resolved, **stage_meta}
                break
        return resolved

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
