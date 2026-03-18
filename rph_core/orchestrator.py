"""
Reaction Profile Hunter Orchestrator
======================================

串行四步走架构的总指挥 + Step 3.5 SP矩阵集成

Author: QCcalc Team
Date: 2026-01-09
Session: #13 - 集成 S3.5 SP矩阵
"""

import logging
from pathlib import Path
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple
from datetime import datetime
import json
import re
import copy

from rph_core.utils.log_manager import setup_logger
from rph_core.utils.path_compat import normalize_path, is_toxic_path
from rph_core.utils.task_builder import build_tasks_from_run_config, sanitize_rx_id
from rph_core.utils.config_loader import load_config
from rph_core.utils.small_molecule_catalog import SmallMoleculeCatalog
from rph_core.utils.checkpoint_manager import CheckpointManager
from rph_core.utils.checkpoint_manager import PipelineState
from rph_core.utils.optimization_config import normalize_qc_config
from rph_core.utils.forming_bonds_resolver import resolve_forming_bonds, write_mechanism_meta
from rph_core.utils.file_io import read_xyz
from rph_core.utils.cleaner_adapter import map_pairs_to_xyz_indices, parse_pairs_text
from rph_core.version import __version__
from rph_core.utils import ui, notify
from rph_core.utils.ui import get_progress_manager
from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.steps.runners import run_step2, run_step3, run_step4


@dataclass
class PipelineResult:
    """流水线结果"""
    success: bool
    product_smiles: Optional[str] = None
    work_dir: Optional[Path] = None

    # Step outputs
    product_xyz: Optional[Path] = None
    e_product_l2: Optional[float] = None
    product_checkpoint: Optional[Path] = None
    product_thermo: Optional[Path] = None
    product_fchk: Optional[Path] = None
    product_log: Optional[Path] = None
    product_qm_output: Optional[Path] = None
    ts_guess_xyz: Optional[Path] = None
    substrate_xyz: Optional[Path] = None
    intermediate_xyz: Optional[Path] = None
    ts_final_xyz: Optional[Path] = None
    features_csv: Optional[Path] = None

    # Metadata from steps
    forming_bonds: Optional[Tuple[Tuple[int, int], ...]] = None
    sp_matrix_report: Optional[Any] = None
    ts_fchk: Optional[Path] = None
    ts_log: Optional[Path] = None
    ts_qm_output: Optional[Path] = None
    intermediate_fchk: Optional[Path] = None
    intermediate_log: Optional[Path] = None
    intermediate_qm_output: Optional[Path] = None

    # Error tracking
    error_step: Optional[str] = None
    error_message: Optional[str] = None

    def __str__(self):
        if self.success:
            l2_info = f", L2: {self.e_product_l2:.6f} Ha" if self.e_product_l2 else ""
            return f"✅ Pipeline 成功: {self.product_smiles}\n" \
                   f"   Product: {self.product_xyz}{l2_info}\n" \
                   f"   TS Final: {self.ts_final_xyz}\n" \
                   f"   Features: {self.features_csv}"
        else:
            return f"❌ Pipeline 失败: {self.error_step}\n" \
                   f"   错误: {self.error_message}"


class ReactionProfileHunter:
    """
    Reaction Profile Hunter v3.0 (分子自治架构)

    设计模式:
    - Orchestrator 只是"调度员"，不是"工人"
    - 每个 Step 是独立的"工人"，有明确的输入输出
    - v3.0 核心改进：分子自治目录结构 + OPT-SP 耦合循环

    串行流程:
    S1: AnchorPhase (分子锚定 + CREST + DFT OPT-SP 耦合)
      → S2: RetroScanner (逆向扫描，从 S1_ConfGeneration/[Molecule]/dft 读取)
      → S3: TSOptimizer (TS优化，使用 LogParser 提取坐标)
      → S4: FeatureMiner (特征提取)

    v3.0 目录结构:
    S1_ConfGeneration/[Molecule_Name]/
        ├── crest/          # CREST 搜索结果
        ├── dft/            # DFT OPT + SP (无子目录，扁平结构)
        └── [Molecule_Name]_global_min.xyz
    """

    def __init__(self, config_path: Optional[Path] = None, log_level: Optional[str] = None):
        """
        初始化 Reaction Profile Hunter

        Args:
            config_path: 配置文件路径（可选，默认使用defaults.yaml）
        """
        # 加载配置
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "defaults.yaml"
        self.config_path: Path = Path(config_path)

        self.config = load_config(config_path)
        if log_level:
            self.config.setdefault("global", {})["log_level"] = log_level
        self.config, qc_fixes = normalize_qc_config(self.config, auto_fix=True)
        self.logger = setup_logger(
            "ReactionProfileHunter",
            level=self.config.get('global', {}).get('log_level', 'INFO')
        )
        if qc_fixes:
            self.logger.warning(f"QC config normalized: {len(qc_fixes)} change(s)")
            for fix in qc_fixes:
                self.logger.debug(
                    f"QC config fix [{fix['field']}]: {fix['original']} -> {fix['updated']}"
                )

        ui.print_pipeline_header(__version__)
        self.logger.info(f"Reaction Profile Hunter v{__version__} 初始化 (含 S3.5)")
        
        # Initialize small molecule catalog
        self.small_mol_catalog = SmallMoleculeCatalog(self.config)

        # 延迟初始化各步骤引擎（懒加载）
        self._s1_engine = None
        self._s1_engine_type = "ProductAnchor"  # 默认类型
        self._s2_engine = None
        self._s3_engine = None
        self._s4_engine = None

    @property
    def s1_engine(self):
        """Step 1 引擎（懒加载）- v3.0: 分子自治架构"""
        if self._s1_engine is None:
            # v3.0 强制使用新的 AnchorPhase
            try:
                from rph_core.steps.anchor.handler import AnchorPhase
                self._s1_engine = AnchorPhase(
                    config=self.config,
                    base_work_dir=Path.cwd()  # 默认值，将在运行时更新
                )
                self._s1_engine_type = "AnchorPhase_v3"
                self.logger.debug("✓ Step 1 使用 AnchorPhase v3.0（分子自治架构）")
                return self._s1_engine
            except ImportError as e:
                self.logger.error(f"无法导入 AnchorPhase v3.0: {e}")
                self.logger.error("v3.0 要求必须使用新的 AnchorPhase")
                raise RuntimeError(
                    "ReactionProfileHunter v3.0 要求必须使用 AnchorPhase v3.0。"
                    "请确保 rph_core/steps/anchor/handler.py 存在且可导入。"
                )
        return self._s1_engine

    @property
    def s2_engine(self):
        """Step 2 引擎 (懒加载)"""
        if self._s2_engine is None:
            from rph_core.steps.step2_retro import RetroScanner
            self._s2_engine = RetroScanner(self.config)
            self.logger.debug("✓ Step 2 (RetroScanner) 已初始化")
        return self._s2_engine

    @property
    def s3_engine(self):
        """Step 3 引擎 (懒加载)"""
        if self._s3_engine is None:
            from rph_core.steps.step3_opt import TSOptimizer
            self._s3_engine = TSOptimizer(self.config)
            self.logger.debug("✓ Step 3 (TSOptimizer) 已初始化")
        return self._s3_engine



    @property
    def s4_engine(self):
        """Step 4 引擎 (懒加载)"""
        if self._s4_engine is None:
            from rph_core.steps.step4_features import FeatureMiner
            self._s4_engine = FeatureMiner(self.config)
            self.logger.debug("✓ Step 4 (FeatureMiner) 已初始化")
        return self._s4_engine

    def _resolve_profile_key(
        self,
        reaction_profile: Optional[str] = None,
        cleaner_data: Optional[Dict[str, Any]] = None,
    ) -> Optional[str]:
        reaction_profiles = self.config.get("reaction_profiles", {}) or {}

        explicit_profile = (
            reaction_profile
            or self.config.get("reaction_profile")
            or (self.config.get("run", {}) or {}).get("reaction_profile")
            or (self.config.get("run", {}) or {}).get("reaction_type")
        )
        if explicit_profile:
            profile_key = str(explicit_profile).strip()
            if profile_key:
                return profile_key

        if cleaner_data:
            cleaner_profile = cleaner_data.get("reaction_profile")
            if cleaner_profile:
                profile_key = str(cleaner_profile).strip()
                if profile_key:
                    return profile_key

            reaction_type = cleaner_data.get("reaction_type") or cleaner_data.get("rxn_type")
            if reaction_type:
                rt = str(reaction_type).strip()
                candidates: List[str] = []

                def _add(candidate: str) -> None:
                    if candidate and candidate not in candidates:
                        candidates.append(candidate)

                normalized = rt.replace(" ", "")
                _add(rt)
                _add(f"{rt}_default")
                _add(normalized)
                _add(f"{normalized}_default")

                if normalized.startswith("[") and normalized.endswith("]") and len(normalized) > 2:
                    inner = normalized[1:-1]
                    _add(inner)
                    _add(f"{inner}_default")
                else:
                    bracketed = f"[{normalized}]"
                    _add(bracketed)
                    _add(f"{bracketed}_default")

                for key in candidates:
                    if key in reaction_profiles:
                        return key

        return None

    def _resolve_forward_scan_config(
        self,
        reaction_profile: Optional[str] = None,
        cleaner_data: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        step2_cfg = self.config.get("step2", {}) or {}
        forward_cfg = dict(step2_cfg.get("forward_scan", {}) or {})
        forward_cfg.update(dict(step2_cfg.get("scan", {}) or {}))
        reaction_profiles = self.config.get("reaction_profiles", {}) or {}

        profile_key = self._resolve_profile_key(reaction_profile=reaction_profile, cleaner_data=cleaner_data)
        profile_cfg = {}
        if profile_key and isinstance(reaction_profiles, dict):
            profile_cfg = dict(reaction_profiles.get(str(profile_key), {}) or {})
        if not profile_cfg and isinstance(reaction_profiles, dict):
            profile_cfg = dict(reaction_profiles.get("_universal", {}) or {})

        scan_cfg = dict(profile_cfg.get("scan", {}) or {})

        scan_start_distance = scan_cfg.get(
            "scan_start_distance",
            scan_cfg.get("initial_distance", scan_cfg.get("ts_distance", forward_cfg.get("scan_start_distance", 2.2))),
        )
        scan_end_distance = scan_cfg.get(
            "scan_end_distance",
            scan_cfg.get("break_distance", forward_cfg.get("scan_end_distance", 3.5)),
        )

        merged = dict(forward_cfg)
        merged["scan_start_distance"] = float(scan_start_distance)
        merged["scan_end_distance"] = float(scan_end_distance)
        merged["scan_steps"] = int(scan_cfg.get("scan_steps", merged.get("scan_steps", 10)))
        merged["scan_mode"] = str(scan_cfg.get("scan_mode", merged.get("scan_mode", "concerted")))
        merged["scan_force_constant"] = float(
            scan_cfg.get("scan_force_constant", merged.get("scan_force_constant", 1.0))
        )
        merged["min_valid_points"] = int(scan_cfg.get("min_valid_points", merged.get("min_valid_points", 5)))
        merged["reject_boundary_maximum"] = bool(
            scan_cfg.get("reject_boundary_maximum", merged.get("reject_boundary_maximum", True))
        )
        merged["require_local_peak"] = bool(scan_cfg.get("require_local_peak", merged.get("require_local_peak", False)))
        merged["boundary_retry_once"] = bool(scan_cfg.get("boundary_retry_once", merged.get("boundary_retry_once", True)))
        merged["boundary_retry_delta"] = float(scan_cfg.get("boundary_retry_delta", merged.get("boundary_retry_delta", 0.3)))
        merged["boundary_retry_extra_steps"] = int(
            scan_cfg.get("boundary_retry_extra_steps", merged.get("boundary_retry_extra_steps", 6))
        )
        merged["allow_boundary_degradation"] = bool(
            scan_cfg.get("allow_boundary_degradation", merged.get("allow_boundary_degradation", True))
        )
        return merged

    def _resolve_product_xyz_for_s2(self, product_xyz: Path) -> Path:
        product_xyz = Path(product_xyz)
        if not product_xyz.exists():
            raise FileNotFoundError(f"S2 输入产物路径不存在: {product_xyz}")

        if product_xyz.is_file():
            return product_xyz

        candidates = [
            product_xyz / "product_min.xyz",
            product_xyz / "product" / "product_global_min.xyz",
            product_xyz / "product_global_min.xyz",
            product_xyz / "product" / "global_min.xyz",
            product_xyz / "global_min.xyz",
        ]
        resolved = next((p for p in candidates if p.exists()), None)
        if resolved is None:
            raise RuntimeError(
                f"无法在 {product_xyz} 中找到产物结构文件用于 S2。"
                f"已尝试: {[str(p) for p in candidates]}"
            )
        return resolved

    def _normalize_forming_bonds(
        self,
        raw_forming_bonds: Any,
        *,
        atom_count: Optional[int] = None,
        index_base: Optional[Any] = None,
        require_exact_two: bool = True,
    ) -> Tuple[Tuple[int, int], ...]:
        if raw_forming_bonds is None:
            return tuple()

        if isinstance(raw_forming_bonds, str):
            parsed: List[Tuple[int, int]] = []
            for chunk in re.split(r"[;,]", raw_forming_bonds):
                piece = chunk.strip()
                if not piece:
                    continue
                if "-" not in piece:
                    raise RuntimeError(f"Invalid forming_bonds token '{piece}', expected 'i-j'")
                left, right = piece.split("-", 1)
                try:
                    parsed.append((int(left.strip()), int(right.strip())))
                except (TypeError, ValueError) as exc:
                    raise RuntimeError(f"Invalid forming_bonds token '{piece}', expected integer pair") from exc
            raw_pairs = parsed
        else:
            if not isinstance(raw_forming_bonds, (list, tuple)):
                return tuple()

            raw_pairs = []
            for pair in raw_forming_bonds:
                if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                    continue
                try:
                    i, j = int(pair[0]), int(pair[1])
                except (TypeError, ValueError):
                    continue
                raw_pairs.append((i, j))

        if not raw_pairs:
            return tuple()

        normalized: List[Tuple[int, int]]

        index_base_flag: Optional[int]
        if index_base is None:
            index_base_flag = None
        elif isinstance(index_base, int):
            if index_base not in (0, 1):
                raise RuntimeError(f"Invalid forming_bonds index_base={index_base}; expected 0, 1, or auto")
            index_base_flag = index_base
        else:
            index_base_text = str(index_base).strip().lower()
            if index_base_text in {"", "auto"}:
                index_base_flag = None
            elif index_base_text in {"0", "zero", "0-based", "zero_based", "zero-based"}:
                index_base_flag = 0
            elif index_base_text in {"1", "one", "1-based", "one_based", "one-based"}:
                index_base_flag = 1
            else:
                raise RuntimeError(f"Invalid forming_bonds index_base='{index_base}'; expected 0, 1, or auto")

        if index_base_flag is None:
            has_zero = any(i == 0 or j == 0 for i, j in raw_pairs)
            if has_zero:
                index_base_flag = 0
            else:
                if atom_count is None:
                    raise RuntimeError(
                        "forming_bonds index base is ambiguous; provide index_base or product atom count"
                    )
                min_idx = min(min(i, j) for i, j in raw_pairs)
                max_idx = max(max(i, j) for i, j in raw_pairs)
                zero_based_valid = min_idx >= 0 and max_idx <= atom_count - 1
                one_based_valid = min_idx >= 1 and max_idx <= atom_count
                if zero_based_valid and one_based_valid:
                    raise RuntimeError(
                        "forming_bonds index base is ambiguous (valid as both 0-based and 1-based); "
                        "set cleaner/config index_base explicitly"
                    )
                if one_based_valid:
                    index_base_flag = 1
                elif zero_based_valid:
                    index_base_flag = 0
                else:
                    raise RuntimeError(
                        f"forming_bonds indices out of range for atom_count={atom_count}: {raw_pairs}"
                    )

        if index_base_flag == 1:
            normalized = [(i - 1, j - 1) for i, j in raw_pairs]
        else:
            normalized = list(raw_pairs)

        canonical: List[Tuple[int, int]] = []
        for i, j in normalized:
            if i == j:
                raise RuntimeError(f"Invalid forming_bonds pair ({i}, {j}): self-bond is not allowed")
            if i < 0 or j < 0:
                raise RuntimeError(f"Invalid forming_bonds pair ({i}, {j}): negative index")
            if atom_count is not None and (i >= atom_count or j >= atom_count):
                raise RuntimeError(
                    f"forming_bonds pair ({i}, {j}) out of range for atom_count={atom_count}"
                )
            canonical.append((min(i, j), max(i, j)))

        canonical = sorted(set(canonical))

        if require_exact_two and len(canonical) != 2:
            raise RuntimeError(
                f"S2 requires exactly 2 forming bonds for this workflow; got {len(canonical)}: {canonical}"
            )

        unique_atoms = {idx for pair in canonical for idx in pair}
        if require_exact_two and len(unique_atoms) != 4:
            raise RuntimeError(
                f"S2 requires 4 unique forming-bond atoms; got {len(unique_atoms)} from {canonical}"
            )

        return tuple(canonical)

    def _resolve_forming_bonds_for_s2(
        self,
        cleaner_data: Optional[Dict[str, Any]] = None,
        product_xyz_file: Optional[Path] = None,
    ) -> Tuple[Tuple[int, int], ...]:
        step2_cfg = self.config.get("step2", {}) or {}
        cleaner_cfg = self.config.get("cleaner", {}) or {}

        atom_count: Optional[int] = None
        if product_xyz_file is not None and Path(product_xyz_file).exists():
            try:
                coords, _ = read_xyz(Path(product_xyz_file))
                atom_count = int(len(coords))
            except Exception as exc:
                self.logger.warning(f"[S2] Failed to read atom_count from {product_xyz_file}: {exc}")

        cleaner_raw = ((cleaner_data or {}).get("raw", {}) or {}) if cleaner_data else {}

        cleaner_xyz_pairs_raw = (
            (cleaner_data or {}).get("formed_bond_xyz_pairs")
            or cleaner_raw.get("formed_bond_xyz_pairs")
        )
        cleaner_xyz_pairs = self._normalize_forming_bonds(
            cleaner_xyz_pairs_raw,
            atom_count=atom_count,
            index_base=0,
            require_exact_two=True,
        )
        if cleaner_xyz_pairs:
            self.logger.info(f"[S2] Using cleaner-derived formed_bond_xyz_pairs: {cleaner_xyz_pairs}")
            return cleaner_xyz_pairs

        map_pairs_raw = (
            (cleaner_data or {}).get("formed_bond_map_pairs")
            or cleaner_raw.get("formed_bond_map_pairs")
        )
        map_pairs = parse_pairs_text(str(map_pairs_raw) if map_pairs_raw is not None else None)
        mapped_product_smiles = (
            (cleaner_data or {}).get("mapped_product_smiles")
            or cleaner_raw.get("mapped_product_smiles")
            or cleaner_raw.get("product_smiles_mapped")
            or cleaner_raw.get("mapped_product")
        )
        if map_pairs and mapped_product_smiles and product_xyz_file is not None and Path(product_xyz_file).exists():
            try:
                mapped_xyz_pairs = map_pairs_to_xyz_indices(
                    mapped_smiles=str(mapped_product_smiles),
                    map_pairs=map_pairs,
                    xyz_file=Path(product_xyz_file),
                )
            except Exception as exc:
                mapped_xyz_pairs = []
                self.logger.warning(f"[S2] Failed map->XYZ conversion for forming bonds: {exc}")

            mapped_xyz_normalized = self._normalize_forming_bonds(
                mapped_xyz_pairs,
                atom_count=atom_count,
                index_base=0,
                require_exact_two=True,
            )
            if mapped_xyz_normalized:
                self.logger.info(f"[S2] Using map->XYZ resolved forming_bonds: {mapped_xyz_normalized}")
                return mapped_xyz_normalized

        configured = (
            (cleaner_data or {}).get("forming_bonds")
            or (cleaner_data or {}).get("formed_bond_index_pairs")
            or ((cleaner_data or {}).get("raw", {}) or {}).get("formed_bond_index_pairs")
            or step2_cfg.get("forming_bonds")
            or (step2_cfg.get("scan", {}) or {}).get("forming_bonds")
            or (step2_cfg.get("forward_scan", {}) or {}).get("forming_bonds")
            or cleaner_cfg.get("forming_bonds")
            or self.config.get("forming_bonds")
        )
        configured_index_base = (
            (cleaner_data or {}).get("forming_bonds_index_base")
            or (cleaner_data or {}).get("index_base")
            or ((cleaner_data or {}).get("raw", {}) or {}).get("forming_bonds_index_base")
            or ((cleaner_data or {}).get("raw", {}) or {}).get("index_base")
            or step2_cfg.get("forming_bonds_index_base")
            or cleaner_cfg.get("forming_bonds_index_base")
            or self.config.get("forming_bonds_index_base")
            or "auto"
        )

        configured_bonds = self._normalize_forming_bonds(
            configured,
            atom_count=atom_count,
            index_base=configured_index_base,
            require_exact_two=True,
        )
        if configured_bonds:
            if product_xyz_file is not None and Path(product_xyz_file).exists():
                try:
                    from rph_core.utils.geometry_tools import GeometryUtils

                    coords, _symbols = read_xyz(Path(product_xyz_file))
                    for bond in configured_bonds:
                        dist = GeometryUtils.calculate_distance(coords, int(bond[0]), int(bond[1]))
                        if dist > 2.5:
                            self.logger.warning(
                                f"[S2] Forming bond {bond} has distance {dist:.3f} A in product XYZ "
                                "- possible index mapping error (expected < 2.0 A for bonded atoms)"
                            )
                except Exception as exc:
                    self.logger.debug(f"[S2] Forming bond distance check failed: {exc}")
            self.logger.info(f"[S2] Using configured forming_bonds: {configured_bonds}")
            return configured_bonds

        raise RuntimeError("S2 requires forming_bonds from cleaner/config/dataset; SMARTS fallback has been removed")

    def _build_step2_signature(
        self,
        *,
        work_dir: Path,
        product_xyz_file: Path,
        forming_bonds: Tuple[Tuple[int, int], ...],
        reaction_profile: Optional[str],
        scan_config: Dict[str, Any],
    ) -> Dict[str, Any]:
        checkpoint_mgr = CheckpointManager(work_dir)
        return checkpoint_mgr.compute_step2_signature(
            config=self.config,
            product_xyz=product_xyz_file,
            forming_bonds=forming_bonds,
            reaction_profile=reaction_profile,
            scan_config=scan_config,
        )

    def run_pipeline(
        self,
        product_smiles: str,
        work_dir: Path,
        skip_steps: Optional[list[str]] = None,
        precursor_smiles: Optional[str] = None,
        leaving_group_key: Optional[str] = None,
        reaction_profile: Optional[str] = None,
        cleaner_data: Optional[Dict[str, Any]] = None,
    ) -> PipelineResult:
        """
        执行 v2.1 串行四步走架构

        数据流:
        S1_Output (Product_Min) ──┬──> S2_Input
                                  │
        S2_Output (TS_Guess, R) ──┼──> S3_Input
                                  │
        S3_Output (TS_Final) ─────┼──> S4_Input
                                  │
        S1 + S2 + S3 ─────────────┴──> S4_Input

        Args:
            product_smiles: 产物 SMILES
            work_dir: 工作目录
            skip_steps: 要跳过的步骤列表（用于调试）

        Returns:
            PipelineResult 对象
        """
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"🚀 任务启动: {product_smiles}")
        self.logger.info(f"📁 工作目录: {work_dir}")

        result = PipelineResult(
            success=False,
            product_smiles=product_smiles,
            work_dir=work_dir
        )

        def _notify(success: bool, error_step: Optional[str] = None, error_message: Optional[str] = None) -> None:
            title = "RPH 任务完成" if success else "RPH 任务失败"
            if success:
                message = f"{product_smiles} 完成: {work_dir}"
            else:
                step = error_step or "Unknown"
                message = f"{product_smiles} 失败于 {step}: {error_message or 'unknown error'}"
            notify.notify_completion(title, message, self.config)

        skip_steps = skip_steps or []

        pm = get_progress_manager()
        pm.start(f"RPH Pipeline: {product_smiles}")
        pm.add_step("s1", "Step 1: Anchor")
        pm.add_step("s2", "Step 2: Guess Builder")
        pm.add_step("s3", "Step 3: Analyzer")
        pm.add_step("s4", "Step 4: Features")

        try:
            resume_enabled = bool((self.config.get("run", {}) or {}).get("resume", True))
            checkpoint_mgr = CheckpointManager(work_dir)
            if resume_enabled:
                state = checkpoint_mgr.load_state()
                if state is None:
                    run_cfg = self.config.get("run", {}) or {}
                    rehydrate_enabled = bool(run_cfg.get("resume_rehydrate", True))
                    if rehydrate_enabled:
                        rehydrate_policy = run_cfg.get("resume_rehydrate_policy", "best_effort")
                        state = checkpoint_mgr.rehydrate_state_from_artifacts(
                            product_smiles=product_smiles,
                            config=self.config,
                            policy=rehydrate_policy
                        )

                    if state is None:
                        state = PipelineState(
                            product_smiles=product_smiles,
                            work_dir=str(work_dir),
                            start_time=datetime.now().isoformat(),
                            last_update=datetime.now().isoformat(),
                            steps={},
                            config_snapshot={},
                        )
                    checkpoint_mgr.save_state(state)

            current_step2_signature: Optional[Dict[str, Any]] = None

            # Resume: reuse Step1/Step2 outputs to avoid repeated heavy QC on S3 failures.
            if resume_enabled and 's1' not in skip_steps and checkpoint_mgr.is_step_completed('s1'):
                product_xyz = checkpoint_mgr.get_step_output('s1', 'product_xyz')
                if product_xyz and Path(product_xyz).exists():
                    result.product_xyz = Path(product_xyz)
                    product_sp = checkpoint_mgr.get_step_metadata('s1', 'e_product_sp')
                    if product_sp is not None:
                        result.e_product_l2 = float(product_sp)
                    fchk = checkpoint_mgr.get_step_output('s1', 'product_fchk')
                    if fchk and Path(fchk).exists():
                        result.product_fchk = Path(fchk)
                    logp = checkpoint_mgr.get_step_output('s1', 'product_log')
                    if logp and Path(logp).exists():
                        result.product_log = Path(logp)
                    qmp = checkpoint_mgr.get_step_output('s1', 'product_qm_output')
                    if qmp and Path(qmp).exists():
                        result.product_qm_output = Path(qmp)
                    chk = checkpoint_mgr.get_step_output('s1', 'product_checkpoint')
                    if chk and Path(chk).exists():
                        result.product_checkpoint = Path(chk)
                    thermo = checkpoint_mgr.get_step_output('s1', 'product_thermo')
                    if thermo and Path(thermo).exists():
                        result.product_thermo = Path(thermo)

                    self.logger.info(f"✅ Resume: Step1 already complete, reuse product: {result.product_xyz}")
                    skip_steps.append('s1')
                    pm.update_step("s1", completed=100, description="Step 1: Anchor [REUSED]")

            if resume_enabled and 's2' not in skip_steps and checkpoint_mgr.is_step_completed('s2'):
                can_reuse_step2 = False
                profile_key = self._resolve_profile_key(
                    reaction_profile=reaction_profile,
                    cleaner_data=cleaner_data,
                )
                try:
                    if result.product_xyz is not None:
                        product_xyz_file = self._resolve_product_xyz_for_s2(result.product_xyz)
                        resolved_forming_bonds = self._resolve_forming_bonds_for_s2(
                            cleaner_data=cleaner_data,
                            product_xyz_file=product_xyz_file,
                        )
                        expected_scan_cfg = self._resolve_forward_scan_config(
                            reaction_profile=profile_key,
                            cleaner_data=cleaner_data,
                        )
                        expected_scan_cfg["output_dir"] = work_dir / "S2_Retro"
                        current_step2_signature = self._build_step2_signature(
                            work_dir=work_dir,
                            product_xyz_file=product_xyz_file,
                            forming_bonds=resolved_forming_bonds,
                            reaction_profile=profile_key,
                            scan_config=expected_scan_cfg,
                        )
                        cached_step2_signature = checkpoint_mgr.get_step_metadata('s2', 'step2_signature')
                        if cached_step2_signature == current_step2_signature:
                            can_reuse_step2 = True
                        else:
                            self.logger.info("S2 checkpoint: signature mismatch, recomputing Step2")
                    else:
                        self.logger.info("S2 checkpoint: missing S1 product context; recomputing Step2")
                except Exception as exc:
                    self.logger.warning(f"S2 checkpoint: failed to validate signature ({exc}), recomputing Step2")

                if can_reuse_step2:
                    ts_guess = checkpoint_mgr.get_step_output('s2', 'ts_guess_xyz')
                    substrate_xyz = checkpoint_mgr.get_step_output('s2', 'substrate_xyz')
                    if ts_guess and substrate_xyz and Path(ts_guess).exists() and Path(substrate_xyz).exists():
                        result.ts_guess_xyz = Path(ts_guess)
                        result.substrate_xyz = Path(substrate_xyz)
                        intermediate_xyz = checkpoint_mgr.get_step_output('s2', 'intermediate_xyz')
                        if intermediate_xyz and Path(intermediate_xyz).exists():
                            result.intermediate_xyz = Path(intermediate_xyz)
                        cached_forming_bonds = checkpoint_mgr.get_step_metadata('s2', 'forming_bonds')
                        try:
                            restored_forming_bonds = self._normalize_forming_bonds(
                                cached_forming_bonds,
                                index_base=0,
                                require_exact_two=True,
                            )
                        except Exception as exc:
                            self.logger.warning(f"S2 checkpoint: invalid cached forming_bonds ({exc})")
                            restored_forming_bonds = tuple()
                        result.forming_bonds = restored_forming_bonds if restored_forming_bonds else None
                        self.logger.info(
                            f"✅ Resume: Step2 already complete, reuse ts_guess/substrate: {result.ts_guess_xyz}, {result.substrate_xyz}"
                        )
                        skip_steps.append('s2')
                        pm.update_step("s2", completed=100, description="Step 2: Guess Builder [REUSED]")

            s3_dir = work_dir / "S3_TransitionAnalysis"
            if resume_enabled and 's3' not in skip_steps:
                input_hashes: Dict[str, str] | None = None
                if result.ts_guess_xyz and result.substrate_xyz and result.product_xyz:
                    input_hashes = {
                        'ts_guess': checkpoint_mgr.compute_file_hash(result.ts_guess_xyz) or '',
                        'substrate': checkpoint_mgr.compute_file_hash(result.substrate_xyz) or '',
                        'product': checkpoint_mgr.compute_file_hash(result.product_xyz) or '',
                    }

                if checkpoint_mgr.is_step3_complete(
                    s3_dir,
                    self.config,
                    check_signature=True,
                    input_hashes=input_hashes,
                    upstream_step2_signature=current_step2_signature,
                ):
                    self.logger.info("✅ Resume: Step3 checkpoint valid, attempting to restore S3 outputs...")

                    ts_final = checkpoint_mgr.get_step_output('s3', 'ts_final_xyz')
                    sp_meta_path = checkpoint_mgr.get_step_output('s3', 'sp_matrix_metadata_json') or str(s3_dir / "sp_matrix_metadata.json")

                    if ts_final and Path(ts_final).exists() and Path(sp_meta_path).exists():
                        try:
                            with open(sp_meta_path, 'r') as f:
                                sp_meta = json.load(f)

                            from rph_core.steps.step3_opt.ts_optimizer import SPMatrixReport
                            sp_report = SPMatrixReport(
                                e_ts=sp_meta.get('e_ts', 0.0),
                                e_reactant=sp_meta.get('e_reactant'),
                                e_product=sp_meta.get('e_product'),
                                e_ts_final=sp_meta.get('e_ts'),
                                g_ts=sp_meta.get('g_ts'),
                                g_reactant=sp_meta.get('g_reactant'),
                                g_product=sp_meta.get('g_product'),
                                g_ts_source=sp_meta.get('g_ts_source'),
                                g_ts_error=sp_meta.get('g_ts_error'),
                                g_reactant_source=sp_meta.get('g_reactant_source'),
                                g_reactant_error=sp_meta.get('g_reactant_error'),
                                method=sp_meta.get('method', ''),
                                solvent=sp_meta.get('solvent', ''),
                            )

                            result.ts_final_xyz = Path(ts_final)
                            result.sp_matrix_report = sp_report

                            ts_fchk_str = checkpoint_mgr.get_step_output('s3', 'ts_fchk')
                            result.ts_fchk = Path(ts_fchk_str) if ts_fchk_str else None
                            ts_log_str = checkpoint_mgr.get_step_output('s3', 'ts_log')
                            result.ts_log = Path(ts_log_str) if ts_log_str else None
                            ts_qm_str = checkpoint_mgr.get_step_output('s3', 'ts_qm_output')
                            result.ts_qm_output = Path(ts_qm_str) if ts_qm_str else None
                            intermediate_fchk_str = checkpoint_mgr.get_step_output('s3', 'intermediate_fchk')
                            result.intermediate_fchk = Path(intermediate_fchk_str) if intermediate_fchk_str else None
                            intermediate_log_str = checkpoint_mgr.get_step_output('s3', 'intermediate_log')
                            result.intermediate_log = Path(intermediate_log_str) if intermediate_log_str else None
                            intermediate_qm_str = checkpoint_mgr.get_step_output('s3', 'intermediate_qm_output')
                            result.intermediate_qm_output = Path(intermediate_qm_str) if intermediate_qm_str else None

                            skip_steps.append('s3')
                            pm.update_step("s3", completed=100, description="Step 3: Analyzer [REUSED]")
                            self.logger.info("✅ Resume: Step3 restored from checkpoint")

                            if result.product_xyz and result.ts_final_xyz:
                                forming_cfg = self.config.get('step4', {}).get('forming_bonds', {}) or {}
                                resolved = resolve_forming_bonds(
                                    product_xyz=result.product_xyz,
                                    ts_xyz=result.ts_final_xyz,
                                    s3_dir=s3_dir,
                                    s4_dir=work_dir / "S4_Data",
                                    config=forming_cfg,
                                    write_meta=forming_cfg.get('write_meta', True)
                                )
                                result.forming_bonds = resolved.forming_bonds
                                if resolved.warnings:
                                    for w in resolved.warnings:
                                        self.logger.warning(f"Forming bonds resolver: {w}")

                            dg_act = sp_report.get_activation_energy()
                            dg_rxn = sp_report.get_reaction_energy()
                            if dg_act is not None:
                                self.logger.info(f"      ΔG‡ = {dg_act:.3f} kcal/mol (from checkpoint)")
                            if dg_rxn is not None:
                                self.logger.info(f"      ΔG_rxn = {dg_rxn:.3f} kcal/mol (from checkpoint)")
                        except Exception as e:
                            self.logger.warning(f"⚠️ Failed to restore S3 from checkpoint: {e}, will recompute S3")

            # === Step 1: 产物锚定 (Product Anchor - v3.0 分子自治架构) ===
            if 's1' not in skip_steps:
                try:
                    ui.print_step_header("Step 1", "Product Anchor", "Global Minimum Search (v3.0)")
                    self.logger.info(">>> Step 1: 寻找产物全局最低构象 (v3.0 OPT-SP 耦合)...")

                    # v3.0: 使用 AnchorPhase 处理产物
                    s1_work_dir = work_dir / "S1_ConfGeneration"

                    # 设置 AnchorPhase 的工作目录
                    self.s1_engine.base_work_dir = s1_work_dir
                    self.s1_engine.base_work_dir.mkdir(parents=True, exist_ok=True)

                    # Resolve leaving group if provided
                    resolved_lg_smiles = None
                    if leaving_group_key:
                        mol_obj = self.small_mol_catalog.get(leaving_group_key)
                        if mol_obj:
                            resolved_lg_smiles = mol_obj.smiles
                        else:
                            self.logger.warning(f"Leaving group key '{leaving_group_key}' not found in catalog. Skipping.")

                    # Build molecules dictionary for S1
                    molecules = {"product": product_smiles}
                    if precursor_smiles:
                        molecules["precursor"] = precursor_smiles
                    if resolved_lg_smiles:
                        molecules["leaving_group"] = resolved_lg_smiles

                    pm.update_step("s1", description="Running AnchorPhase (CREST + DFT)...")
                    anchor_result = self.s1_engine.run(
                        molecules=molecules
                    )
                    pm.update_step("s1", completed=100, description="Step 1: Anchor [OK]")

                    # 检查执行结果
                    if not anchor_result.success:
                        raise RuntimeError(f"AnchorPhase v3.0 失败: {anchor_result.error_message}")

                    # 从 AnchorPhaseResult 中提取产物数据（v3.0 结构）
                    product_data = anchor_result.anchored_molecules.get("product", {})

                    # v3.0 新结构：xyz 是 SP 输出文件路径，e_sp 是 SP 能量
                    product_sp_out = product_data.get("xyz")
                    e_product_sp = product_data.get("e_sp")

                    if product_sp_out is None or e_product_sp is None:
                        raise RuntimeError(
                            f"AnchorPhase 未返回完整的产物数据。"
                            f"product_data = {product_data}"
                        )

                    # 使用 LogParser 从 SP 输出中提取最终坐标
                    from rph_core.utils.geometry_tools import LogParser
                    coords, symbols, error = LogParser.extract_last_converged_coords(
                        product_sp_out,
                        engine_type='auto'
                    )

                    if coords is None:
                        self.logger.warning(f"无法从 {product_sp_out} 提取坐标: {error}")
                        # 回退：直接使用 SP 输出路径
                        product_min_xyz = product_sp_out
                    else:
                        # 创建最终的产物 XYZ 文件
                        product_min_xyz = s1_work_dir / "product_min.xyz"
                        from rph_core.utils.file_io import write_xyz
                        # 确保 symbols 不为 None
                        if symbols is None:
                            self.logger.warning("未提取到符号，从 SP 输出文件读取")
                            from rph_core.utils.file_io import read_xyz
                            _, fallback_symbols = read_xyz(product_sp_out)
                            symbols = fallback_symbols
                        write_xyz(product_min_xyz, coords, symbols, title=f"Product SP E={e_product_sp:.6f}")

                    # 保存结果（v3.0 使用 e_sp 而不是 e_l2）
                    result.product_xyz = product_min_xyz
                    result.e_product_l2 = e_product_sp  # 保持向后兼容，使用 SP 能量
                    result.product_fchk = product_data.get("fchk")
                    result.product_log = product_data.get("log")
                    result.product_qm_output = product_data.get("qm_output")

                    product_thermo_file = s1_work_dir / "product" / "dft" / "conformer_thermo.csv"
                    if product_thermo_file.exists():
                        result.product_thermo = product_thermo_file

                    # 保存 checkpoint 路径（如果存在）
                    product_checkpoint = product_data.get("chk")
                    if product_checkpoint and product_checkpoint.exists():
                        result.product_checkpoint = product_checkpoint
                        self.logger.info(f"    ✓ S1 checkpoint 可用: {product_checkpoint}")

                    self.logger.info(f"    ✓ 产物锚定完成: {product_min_xyz}")
                    self.logger.info(f"    ✓ SP 能量: {e_product_sp:.8f} Hartree")
                    if product_sp_out.suffix == ".xyz":
                        self.logger.info(f"    ✓ SP 输出文件(已是global_min): {product_sp_out}")
                    else:
                        self.logger.info(f"    ✓ SP 输出文件: {product_sp_out}")

                    if resume_enabled:
                        sp_report_any: Any = result.sp_matrix_report
                        checkpoint_mgr.mark_step_completed(
                            "s1",
                            output_files={
                                "product_xyz": str(result.product_xyz),
                                "product_fchk": str(result.product_fchk) if result.product_fchk else "",
                                "product_log": str(result.product_log) if result.product_log else "",
                                "product_qm_output": str(result.product_qm_output) if result.product_qm_output else "",
                                "product_checkpoint": str(result.product_checkpoint) if result.product_checkpoint else "",
                                "product_thermo": str(result.product_thermo) if result.product_thermo else "",
                            },
                            metadata={"e_product_sp": float(e_product_sp)},
                        )

                except Exception as e:
                    result.error_step = "Step1_ProductAnchor_v3"
                    result.error_message = str(e)
                    self.logger.error(f"Step 1 v3.0 失败: {e}", exc_info=True)
                    _notify(False, result.error_step, result.error_message)
                    return result
            else:
                self.logger.warning("⚠️  跳过 Step 1")
                pm.update_step("s1", completed=100, description="Step 1: Anchor [SKIPPED]")
                s1_candidates = [
                    work_dir / "S1_ConfGeneration",
                    work_dir.parent / "S1_test" / "S1_ConfGeneration",
                    work_dir.parent / "S1_ConfGeneration",
                    work_dir / "S1_Product",
                    work_dir.parent / "S1_test" / "S1_Product",
                    work_dir.parent / "S1_Product"
                ]
                for candidate in s1_candidates:
                    if candidate.exists():
                        result.product_xyz = candidate
                        self.logger.info(f"    ✓ 复用 S1 输出目录: {candidate}")
                        break

            if 's2' not in skip_steps and result.product_xyz:
                try:
                    ui.print_step_header("Step 2", "TS Guess Builder", "Intermediate Optimization + Inward Scan")
                    self.logger.info(">>> Step 2: 生成中间体并通过 inward scan 构建 TS 初猜...")

                    pm.update_step("s2", description="Scanning bond coordinates (xTB)...")
                    step2_artifacts = run_step2(
                        hunter=self,
                        product_xyz=result.product_xyz,
                        work_dir=work_dir,
                        reaction_profile=reaction_profile,
                        cleaner_data=cleaner_data,
                    )

                    pm.update_step("s2", completed=100, description="Step 2: Guess Builder [OK]")
                    result.ts_guess_xyz = step2_artifacts.ts_guess_xyz
                    result.substrate_xyz = step2_artifacts.substrate_xyz
                    result.intermediate_xyz = step2_artifacts.intermediate_xyz
                    result.forming_bonds = step2_artifacts.forming_bonds
                    current_step2_signature = step2_artifacts.step2_signature
                    self.logger.info(f"    ✓ TS initial guess: {step2_artifacts.ts_guess_xyz}")
                    self.logger.info(f"    ✓ Substrate: {step2_artifacts.substrate_xyz}")
                    self.logger.info(f"    ✓ Intermediate: {step2_artifacts.intermediate_xyz}")
                    self.logger.info(f"    ✓ S2 status/confidence: {step2_artifacts.status}/{step2_artifacts.ts_guess_confidence}")

                    if resume_enabled:
                        checkpoint_mgr.mark_step_completed(
                            "s2",
                            output_files={
                                "ts_guess_xyz": str(result.ts_guess_xyz),
                                "substrate_xyz": str(result.substrate_xyz),
                                "intermediate_xyz": str(result.intermediate_xyz) if result.intermediate_xyz else "",
                            },
                            metadata={
                                "s2_generation_method": step2_artifacts.generation_method,
                                "scan_profile_json": str(step2_artifacts.scan_profile_json) if step2_artifacts.scan_profile_json else "",
                                "forming_bonds": [list(b) for b in result.forming_bonds] if result.forming_bonds else None,
                                "step2_signature": current_step2_signature,
                                "status": step2_artifacts.status,
                                "ts_guess_confidence": step2_artifacts.ts_guess_confidence,
                                "degraded_reasons": list(step2_artifacts.degraded_reasons),
                            },
                        )
                except Exception as e:
                    result.error_step = "Step2_TSGuessBuilder"
                    result.error_message = str(e)
                    self.logger.error(f"Step 2 失败: {e}", exc_info=True)
                    _notify(False, result.error_step, result.error_message)
                    return result
            else:
                self.logger.warning("⚠️  跳过 Step 2")
                pm.update_step("s2", completed=100, description="Step 2: Guess Builder [SKIPPED]")

            # === Step 3: 反应分析 (Transition Analyzer) ===
            if 's3' not in skip_steps and result.ts_guess_xyz:
                try:
                    ui.print_step_header("Step 3", "Transition Analyzer", "TS Optimization & Verification")
                    self.logger.info(">>> Step 3: 反应中心全分析 (TS优化 + Reactant/Fragments SP)...")

                    # 传递S1的checkpoint以复用轨道
                    old_checkpoint = result.product_checkpoint
                    if old_checkpoint:
                        self.logger.info(f"  复用S1 checkpoint: {old_checkpoint.name}")

                    if result.product_xyz and result.product_xyz.is_dir():
                        product_dir = result.product_xyz
                        candidates = [
                            product_dir / "product_min.xyz",
                            product_dir / "product" / "product_global_min.xyz",
                            product_dir / "product_global_min.xyz",
                            product_dir / "product" / "global_min.xyz",
                            product_dir / "global_min.xyz"
                        ]
                        product_file = next((p for p in candidates if p.exists()), None)
                        if product_file is None:
                            raise RuntimeError(
                                f"无法在 {product_dir} 中找到产物结构文件用于 S3。"
                                f"已尝试: {[str(p) for p in candidates]}"
                            )
                        result.product_xyz = product_file
                        self.logger.info(f"  ✓ 使用产物文件: {result.product_xyz}")

                    if result.substrate_xyz is None or result.product_xyz is None:
                        raise RuntimeError("Step3 输入缺失: substrate 或 product 为 None")

                    pm.update_step("s3", description="Optimizing Transition State (Berny/QST2)...")
                    step3_artifacts = run_step3(
                        hunter=self,
                        ts_guess_xyz=result.ts_guess_xyz,
                        intermediate_xyz=result.intermediate_xyz,
                        product_xyz=result.product_xyz,
                        work_dir=work_dir,
                        e_product_l2=result.e_product_l2,
                        product_thermo=result.product_thermo,
                        forming_bonds=result.forming_bonds,
                        old_checkpoint=old_checkpoint,
                    )
                    pm.update_step("s3", completed=100, description="Step 3: Analyzer [OK]")

                    result.ts_final_xyz = step3_artifacts.ts_final_xyz
                    result.sp_matrix_report = step3_artifacts.sp_report
                    
                    result.ts_fchk = step3_artifacts.ts_fchk
                    result.ts_log = step3_artifacts.ts_log
                    result.ts_qm_output = step3_artifacts.ts_qm_output
                    result.intermediate_fchk = step3_artifacts.intermediate_fchk
                    result.intermediate_log = step3_artifacts.intermediate_log
                    result.intermediate_qm_output = step3_artifacts.intermediate_qm_output

                    self.logger.info("    ✓ S3 完成")

                    sp_meta_path = s3_dir / "sp_matrix_metadata.json"
                    if sp_meta_path.exists() and current_step2_signature is not None:
                        try:
                            with open(sp_meta_path, "r", encoding="utf-8") as f:
                                sp_meta_data = json.load(f)
                            sp_meta_data["upstream_step2_signature"] = current_step2_signature
                            with open(sp_meta_path, "w", encoding="utf-8") as f:
                                json.dump(sp_meta_data, f, indent=2)
                        except Exception as exc:
                            self.logger.warning(f"Failed to annotate sp_matrix_metadata.json with S2 signature: {exc}")

                    if result.product_xyz and result.ts_final_xyz and result.forming_bonds is None:
                        forming_cfg = self.config.get('step4', {}).get('forming_bonds', {}) or {}
                        resolved = resolve_forming_bonds(
                            product_xyz=result.product_xyz,
                            ts_xyz=result.ts_final_xyz,
                            s3_dir=s3_dir,
                            s4_dir=work_dir / "S4_Data",
                            config=forming_cfg,
                            write_meta=forming_cfg.get('write_meta', True)
                        )
                        result.forming_bonds = resolved.forming_bonds
                        if resolved.warnings:
                            for w in resolved.warnings:
                                self.logger.warning(f"Forming bonds resolver: {w}")
                    elif result.forming_bonds is not None:
                        try:
                            write_mechanism_meta(
                                {
                                    "version": "1",
                                    "index_base": 0,
                                    "forming_bonds": [list(b) for b in result.forming_bonds],
                                    "source": {
                                        "derived_from_steps": ["S2"],
                                        "step2_signature": current_step2_signature,
                                    },
                                    "validation": {"status": "pass", "warnings": []},
                                },
                                s3_dir / "mechanism_meta.json",
                            )
                        except Exception as exc:
                            self.logger.warning(f"Failed to write S2-derived mechanism_meta.json: {exc}")
                    
                    if resume_enabled:
                        checkpoint_mgr.mark_step_completed(
                            "s3",
                            output_files={
                                "ts_final_xyz": str(result.ts_final_xyz),
                                "ts_fchk": str(result.ts_fchk) if result.ts_fchk else "",
                                "ts_log": str(result.ts_log) if result.ts_log else "",
                                "ts_qm_output": str(result.ts_qm_output) if result.ts_qm_output else "",
                                "intermediate_fchk": str(result.intermediate_fchk) if result.intermediate_fchk else "",
                                "intermediate_log": str(result.intermediate_log) if result.intermediate_log else "",
                                "intermediate_qm_output": str(result.intermediate_qm_output) if result.intermediate_qm_output else "",
                                "sp_matrix_metadata_json": str(s3_dir / "sp_matrix_metadata.json"),
                            },
                            metadata={
                                "step3_signature": checkpoint_mgr._compute_step3_signature(self.config),
                                "input_hashes": {
                                    "ts_guess": checkpoint_mgr.compute_file_hash(result.ts_guess_xyz) if result.ts_guess_xyz else "",
                                    "intermediate": checkpoint_mgr.compute_file_hash(result.intermediate_xyz) if result.intermediate_xyz else "",
                                    "product": checkpoint_mgr.compute_file_hash(result.product_xyz) if result.product_xyz else "",
                                },
                                "energies": {
                                    "e_ts": getattr(result.sp_matrix_report, "e_ts_final", None),
                                    "e_reactant": getattr(result.sp_matrix_report, "e_reactant", None),
                                    "e_product": getattr(result.sp_matrix_report, "e_product", None),
                                },
                                "upstream_step2_signature": current_step2_signature,
                            }
                        )
                    get_act = getattr(result.sp_matrix_report, "get_activation_energy", None)
                    get_rxn = getattr(result.sp_matrix_report, "get_reaction_energy", None)
                    dg_act = get_act() if callable(get_act) else None
                    dg_rxn = get_rxn() if callable(get_rxn) else None
                    if dg_act is not None:
                        self.logger.info(f"      ΔG‡ = {dg_act:.3f} kcal/mol")
                    else:
                        self.logger.info("      ΔG‡ = N/A")
                    if dg_rxn is not None:
                        self.logger.info(f"      ΔG_rxn = {dg_rxn:.3f} kcal/mol")
                    else:
                        self.logger.info("      ΔG_rxn = N/A")
                except Exception as e:
                    result.error_step = "Step3_TransitionAnalyzer"
                    result.error_message = str(e)
                    self.logger.error(f"Step 3 失败: {e}", exc_info=True)
                    _notify(False, result.error_step, result.error_message)
                    return result
            else:
                self.logger.warning("⚠️  跳过 Step 3")
                pm.update_step("s3", completed=100, description="Step 3: Analyzer [SKIPPED]")

            # === Step 4: 特征挖掘 (Feature Miner) ===
            if 's4' not in skip_steps and result.ts_final_xyz:
                try:
                    ui.print_step_header("Step 4", "Feature Miner", "Extracting Features (Extract-Only)")
                    self.logger.info(">>> Step 4: 提取物理有机特征...")

                    # Check S3 artifacts for S4 and warn if missing
                    if result.ts_fchk is None:
                        self.logger.warning("TS fchk missing: formchk failed or not produced. S4 will degrade.")
                    if result.intermediate_fchk is None:
                        self.logger.warning("Reactant fchk missing: formchk failed or not produced. S4 will degrade.")
                    if result.product_fchk is None:
                        self.logger.warning("Product fchk missing: formchk failed or not produced. S4 will degrade.")
                    if result.ts_log is None and result.ts_qm_output is None:
                        self.logger.warning("TS log/out missing: Gaussian .log or ORCA .out not available. S4 will degrade.")
                    if result.intermediate_log is None and result.intermediate_qm_output is None:
                        self.logger.warning("Reactant log/out missing: Gaussian .log or ORCA .out not available. S4 will degrade.")
                    if result.product_log is None and result.product_qm_output is None:
                        self.logger.warning("Product log/out missing: Gaussian .log or ORCA .out not available. S4 will degrade.")

                    if result.substrate_xyz is None or result.product_xyz is None:
                        raise RuntimeError("Step4 输入缺失: reactant 或 product 为 None")

                    if current_step2_signature is not None:
                        sp_meta_path = s3_dir / "sp_matrix_metadata.json"
                        if sp_meta_path.exists():
                            try:
                                with open(sp_meta_path, "r", encoding="utf-8") as f:
                                    sp_meta_data = json.load(f)
                                s3_upstream_sig = sp_meta_data.get("upstream_step2_signature")
                                if s3_upstream_sig is not None and s3_upstream_sig != current_step2_signature:
                                    raise RuntimeError(
                                        "S3/S4 contract mismatch: S3 metadata Step2 signature does not match current Step2 inputs"
                                    )
                            except RuntimeError:
                                raise
                            except Exception as exc:
                                self.logger.warning(f"Step4 contract check skipped due to metadata read issue: {exc}")

                    pm.update_step("s4", description="Extracting features (thermo, geom, qc)...")
                    step4_artifacts = run_step4(
                        hunter=self,
                        ts_final_xyz=result.ts_final_xyz,
                        substrate_xyz=result.substrate_xyz,
                        product_xyz=result.product_xyz,
                        work_dir=work_dir,
                        forming_bonds=result.forming_bonds,
                        sp_matrix_report=result.sp_matrix_report,
                        ts_fchk=result.ts_fchk,
                        intermediate_fchk=result.intermediate_fchk,
                        product_fchk=result.product_fchk,
                        ts_log=result.ts_log,
                        intermediate_log=result.intermediate_log,
                        product_log=result.product_log,
                        ts_qm_output=result.ts_qm_output,
                        intermediate_qm_output=result.intermediate_qm_output,
                        product_qm_output=result.product_qm_output,
                    )
                    pm.update_step("s4", completed=100, description="Step 4: Features [OK]")
                    result.features_csv = step4_artifacts.features_csv
                    self.logger.info(f"    ✓ 特征提取完成: {step4_artifacts.features_csv}")
                except Exception as e:
                    result.error_step = "Step4_FeatureMiner"
                    result.error_message = str(e)
                    self.logger.error(f"Step 4 失败: {e}", exc_info=True)
                    _notify(False, result.error_step, result.error_message)
                    return result
            else:
                self.logger.warning("⚠️  跳过 Step 4")
                pm.update_step("s4", completed=100, description="Step 4: Features [SKIPPED]")

            # 成功完成
            result.success = True
            ui.print_result_summary(result)
            self.logger.info(f"✅ 任务完成! 数据已保存至: {work_dir}")
            _notify(True)

            return result
        finally:
            pm.stop()

    def run_batch(
        self,
        smiles_list: list[str],
        work_dir: Path,
        max_workers: int = 1,
        skip_steps: Optional[list[str]] = None
    ):
        if max_workers <= 0:
            max_workers = int((self.config.get("run", {}) or {}).get("max_workers", 1) or 1)
        skip_steps = list(skip_steps) if skip_steps is not None else []
        self.logger.info(f"开始批量处理: {len(smiles_list)} 个分子")
        self.logger.info(f"并发数: {max_workers}")

        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        from concurrent.futures import ProcessPoolExecutor, as_completed
        from tqdm import tqdm

        results = []
        tasks = []

        for idx, smiles in enumerate(smiles_list):
            task_work_dir = work_dir / f"task_{idx}_{smiles[:10]}"
            tasks.append((idx, smiles, task_work_dir, skip_steps))

        if max_workers == 1:
            for idx, smiles, task_dir, task_skip_steps in tasks:
                try:
                    result = self.run_pipeline(smiles, task_dir, task_skip_steps)
                    results.append(result)
                except Exception as e:
                    self.logger.error(f"Task {idx} 失败: {e}", exc_info=True)
            self.logger.info(f"批量处理完成: {len(results)}/{len(smiles_list)} 成功")
            return results

        config_path_str = str(self.config_path) if getattr(self, "config_path", None) else None
        log_level = self.config.get('global', {}).get('log_level', 'INFO')

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    _run_pipeline_batch_task,
                    config_path_str,
                    log_level,
                    smiles,
                    str(task_dir),
                    list(task_skip_steps),
                ): (idx, smiles, task_dir)
                for idx, smiles, task_dir, task_skip_steps in tasks
            }

            # 使用 tqdm 显示进度
            with tqdm(total=len(tasks), desc="处理进度") as pbar:
                for future in as_completed(futures):
                    idx, smiles, _ = futures[future]
                    try:
                        result = future.result()
                        results.append(result)
                        status = "✅" if result.success else "❌"
                        pbar.set_postfix_str(f"{status} {idx}: {smiles[:15]}")
                    except Exception as e:
                        self.logger.error(f"Task {idx} 失败: {e}")
                        pbar.set_postfix_str(f"❌ {idx}: Error")

        self.logger.info(f"批量处理完成: {len(results)}/{len(smiles_list)} 成功")
        return results

    def _resolve_s1_artifacts(self, work_dir: Path) -> Dict[str, Optional[Path]]:
        """Resolve S1 artifact paths for Step4 Step1 activation features.

        Searches for and optionally derives S1 artifacts needed by S4:
        - shermo_summary.json (from .sum files if missing)
        - HOAc thermo.json (from HOAc .sum if available)
        - conformer_energies.json
        - precursor xyz

        Args:
            work_dir: Pipeline working directory

        Returns:
            Dictionary with resolved artifact paths
        """
        from rph_core.utils.shermo_runner import (
            find_shermo_sum_files,
            derive_shermo_summary_from_sum,
            derive_hoac_thermo_from_sum
        )

        artifacts: Dict[str, Optional[Path]] = {
            "s1_dir": None,
            "s1_shermo_summary_file": None,
            "s1_hoac_thermo_file": None,
            "s1_conformer_energies_file": None,
            "s1_precursor_xyz": None
        }

        # Find S1 directory
        s1_candidates = [
            work_dir / "S1_ConfGeneration",
            work_dir.parent / "S1_test" / "S1_ConfGeneration",
            work_dir.parent / "S1_ConfGeneration",
            work_dir / "S1_Product",
            work_dir.parent / "S1_test" / "S1_Product",
            work_dir.parent / "S1_Product"
        ]
        s1_dir = None
        for candidate in s1_candidates:
            if candidate.exists():
                s1_dir = candidate
                artifacts["s1_dir"] = s1_dir
                break

        if s1_dir is None:
            self.logger.warning("S1 directory not found, Step1 activation features will be NaN")
            return artifacts

        # Look for existing shermo_summary.json
        shermo_summary = s1_dir / "shermo_summary.json"
        if shermo_summary.exists():
            artifacts["s1_shermo_summary_file"] = shermo_summary
        else:
            # Try to derive from .sum files
            sum_files = find_shermo_sum_files(s1_dir)
            if sum_files.get("precursor") or sum_files.get("ylide"):
                # Derive summary for precursor/ylide
                precursor_sum = sum_files.get("precursor")
                ylide_sum = sum_files.get("ylide")

                summary_data = {
                    "unit": "kcal/mol",
                    "temperature_K": 298.15,
                    "derived_artifacts": True
                }

                if precursor_sum:
                    from rph_core.utils.shermo_runner import _parse_sum_file
                    thermo = _parse_sum_file(precursor_sum)
                    g_val = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
                    summary_data["g_precursor"] = g_val * HARTREE_TO_KCAL
                    summary_data["derived_from_precursor"] = str(precursor_sum)

                if ylide_sum:
                    from rph_core.utils.shermo_runner import _parse_sum_file
                    thermo = _parse_sum_file(ylide_sum)
                    g_val = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum
                    summary_data["g_ylide"] = g_val * HARTREE_TO_KCAL
                    summary_data["derived_from_ylide"] = str(ylide_sum)

                # Write derived summary
                shermo_summary.parent.mkdir(parents=True, exist_ok=True)
                with open(shermo_summary, 'w') as f:
                    import json
                    json.dump(summary_data, f, indent=2)
                artifacts["s1_shermo_summary_file"] = shermo_summary
                self.logger.info(f"Derived shermo_summary.json from .sum files")

        # Look for HOAc thermo
        hoac_thermo = s1_dir / "small_molecules" / "HOAc" / "thermo.json"
        if hoac_thermo.exists():
            artifacts["s1_hoac_thermo_file"] = hoac_thermo
        else:
            # Try to derive from HOAc .sum
            sum_files = find_shermo_sum_files(s1_dir)
            hoac_sum = sum_files.get("hoac")
            if hoac_sum is not None:
                hoac_thermo.parent.mkdir(parents=True, exist_ok=True)
                derive_hoac_thermo_from_sum(hoac_sum, hoac_thermo)
                artifacts["s1_hoac_thermo_file"] = hoac_thermo
                self.logger.info(f"Derived HOAc thermo.json from {hoac_sum.name}")
        
        if artifacts["s1_hoac_thermo_file"] is None:
            global_cache_dir = self.config.get("global", {}).get("small_molecule_cache_dir")
            if global_cache_dir:
                from rph_core.utils.small_molecule_cache import SmallMoleculeCache
                from rph_core.utils.molecule_utils import get_molecule_key
                
                cache = SmallMoleculeCache(Path(global_cache_dir))
                hoac_smiles = "CC(=O)O"
                hoac_key = get_molecule_key(hoac_smiles)
                
                if hoac_key:
                    global_hoac_thermo = cache.cache_root / hoac_key / "thermo.json"
                    if global_hoac_thermo.exists():
                        artifacts["s1_hoac_thermo_file"] = global_hoac_thermo
                        self.logger.info(f"Found HOAc thermo.json in global cache: {global_hoac_thermo}")
                    else:
                        global_hoac_dft = cache.cache_root / hoac_key / "dft"
                        if global_hoac_dft.exists():
                            sum_files = find_shermo_sum_files(global_hoac_dft)
                            hoac_sum = sum_files.get("hoac")
                            if hoac_sum is not None:
                                global_hoac_thermo.parent.mkdir(parents=True, exist_ok=True)
                                derive_hoac_thermo_from_sum(hoac_sum, global_hoac_thermo)
                                artifacts["s1_hoac_thermo_file"] = global_hoac_thermo
                                self.logger.info(f"Derived HOAc thermo.json from global cache .sum: {hoac_sum.name}")

        # Look for conformer_energies.json
        conformer_energies = s1_dir / "conformer_energies.json"
        if not conformer_energies.exists():
            # Check in molecule subdirectories
            for mol_dir in s1_dir.iterdir():
                if mol_dir.is_dir() and not mol_dir.name.startswith('.'):
                    candidate = mol_dir / "dft" / "conformer_energies.json"
                    if candidate.exists():
                        conformer_energies = candidate
                        break
        if conformer_energies.exists():
            artifacts["s1_conformer_energies_file"] = conformer_energies

        # Look for precursor xyz
        precursor_xyz = s1_dir / "precursor" / "precursor_min.xyz"
        if precursor_xyz.exists():
            artifacts["s1_precursor_xyz"] = precursor_xyz

        return artifacts


# =============================================================================
# 命令行接口
# =============================================================================

def _resolve_run_config(config: dict[str, Any], args) -> dict[str, Any]:
    run_cfg = copy.deepcopy(config.get("run", {}) or {})
    global_cfg = config.get("global", {}) or {}
    run_cfg.setdefault("source", "single")
    run_cfg.setdefault("output_root", global_cfg.get("work_dir_base", "./rph_output"))
    run_cfg.setdefault("workdir_naming", "rx_{rx_id}")
    run_cfg.setdefault("resume", True)
    run_cfg.setdefault("dry_run", False)
    run_cfg.setdefault("max_tasks", 0)
    run_cfg.setdefault("filter_ids", [])
    run_cfg.setdefault("filter_rx_id", None)

    if args.output:
        run_cfg["output_root"] = args.output

    if args.smiles:
        run_cfg["source"] = "single"
        run_cfg["single"] = _deep_merge_dict(run_cfg.get("single", {}), {
            "rx_id": "manual",
            "product_smiles": args.smiles,
        })

    if getattr(args, "reaction_type", None):
        reaction_profile = str(args.reaction_type).strip()
        if reaction_profile:
            run_cfg["reaction_profile"] = reaction_profile
            run_cfg["reaction_type"] = reaction_profile

    rx_id_value = getattr(args, "rx_id", None)
    if rx_id_value:
        rx_id_text = str(rx_id_value).strip()
        if rx_id_text:
            run_cfg["filter_rx_id"] = rx_id_text

    _apply_filter_rx_id_to_filter_ids(run_cfg)

    return run_cfg


def _apply_filter_rx_id_to_filter_ids(run_cfg: dict[str, Any]) -> None:
    rx_id_value = run_cfg.get("filter_rx_id")
    if not rx_id_value:
        return

    rx_id_text = str(rx_id_value).strip()
    if not rx_id_text:
        return

    canonical = rx_id_text[3:] if rx_id_text.startswith("rx_") else rx_id_text
    run_cfg["filter_rx_id"] = canonical

    existing = run_cfg.get("filter_ids")
    if isinstance(existing, list):
        filter_ids = [str(x).strip() for x in existing if str(x).strip()]
    else:
        filter_ids = []

    if canonical not in filter_ids:
        filter_ids.append(canonical)
    prefixed = f"rx_{canonical}"
    if prefixed not in filter_ids:
        filter_ids.append(prefixed)

    run_cfg["filter_ids"] = filter_ids


def _deep_merge_dict(base: Any, override: Any) -> Any:
    if not isinstance(base, dict):
        base = {}
    if not isinstance(override, dict):
        return copy.deepcopy(override)
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge_dict(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def _run_pipeline_batch_task(
    config_path: Optional[str],
    log_level: str,
    product_smiles: str,
    task_dir: str,
    skip_steps: list[str],
) -> PipelineResult:
    hunter = ReactionProfileHunter(config_path=Path(config_path) if config_path else None, log_level=log_level)
    return hunter.run_pipeline(product_smiles=product_smiles, work_dir=Path(task_dir), skip_steps=skip_steps)


def _run_tasks(hunter: ReactionProfileHunter, run_cfg: dict[str, Any]) -> list[PipelineResult]:
    run_cfg_for_build = copy.deepcopy(run_cfg)
    if str(run_cfg_for_build.get("source", "single")) == "dataset":
        dataset_cfg_value = run_cfg_for_build.get("dataset") or {}
        dataset_cfg = dict(dataset_cfg_value) if isinstance(dataset_cfg_value, dict) else {}
        dataset_cfg["reaction_profiles"] = hunter.config.get("reaction_profiles", {}) or {}
        run_cfg_for_build["dataset"] = dataset_cfg

    tasks = build_tasks_from_run_config(run_cfg_for_build)
    output_root = normalize_path(str(run_cfg.get("output_root", "./rph_output")))

    if is_toxic_path(output_root):
        hunter.logger.warning(
            f"Output path contains toxic characters: {output_root}. "
            "Consider using a safe directory without spaces/brackets."
        )

    output_root.mkdir(parents=True, exist_ok=True)

    global_cfg = hunter.config.setdefault("global", {})
    if not global_cfg.get("small_molecule_cache_dir"):
        global_cache_dir = output_root / "small_molecules"
        global_cfg["small_molecule_cache_dir"] = str(global_cache_dir)
        hunter.logger.info(f"Global small molecule cache configured: {global_cache_dir}")

    results = []
    for task in tasks:
        rx_id = sanitize_rx_id(task.rx_id)
        work_dir = output_root / run_cfg["workdir_naming"].format(rx_id=rx_id)

        if run_cfg.get("dry_run", False):
            hunter.logger.info(f"[dry-run] {task.rx_id} -> {work_dir}")
            continue

        if run_cfg.get("resume", True):
            s4_dir = work_dir / "S4_Data"
            if s4_dir.exists():
                checkpoint_mgr = CheckpointManager(work_dir)
                if checkpoint_mgr.is_step4_complete(s4_dir, hunter.config):
                    hunter.logger.info(f"Skip {task.rx_id}: Step4 already complete")
                    continue

        result = hunter.run_pipeline(
            product_smiles=task.product_smiles,
            work_dir=work_dir,
            skip_steps=[],
            precursor_smiles=task.meta.get("precursor_smiles"),
            leaving_group_key=task.meta.get("leaving_small_molecule_key"),
            reaction_profile=task.meta.get("reaction_profile") or run_cfg.get("reaction_profile"),
            cleaner_data=task.meta.get("cleaner_data") if isinstance(task.meta.get("cleaner_data"), dict) else None,
        )
        results.append(result)

    return results


def main():
    """命令行主入口"""
    import argparse

    parser = argparse.ArgumentParser(
        description="ReactionProfileHunter v6.2 - 过渡态搜索与特征提取"
    )
    parser.add_argument(
        '--smiles',
        type=str,
        default=None,
        help='产物 SMILES 字符串（可选；默认使用 config.run）'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='输出目录（覆盖 config.run.output_root）'
    )
    parser.add_argument(
        '--config',
        type=str,
        default=None,
        help='配置文件路径'
    )
    parser.add_argument(
        '--log-level',
        type=str,
        default=None,
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='日志级别（覆盖 config.global.log_level）'
    )
    parser.add_argument(
        '--reaction-type',
        type=str,
        default="[5+2]_default",
        help='反应类型，用于匹配 reaction_profiles 配置'
    )
    parser.add_argument(
        '--rx-id',
        type=str,
        default=None,
        help='仅运行指定反应ID（可传 9422028 或 rx_9422028）'
    )

    args = parser.parse_args()

    try:
        hunter = ReactionProfileHunter(
            config_path=Path(args.config) if args.config else None,
            log_level=args.log_level
        )

        run_cfg = _resolve_run_config(hunter.config, args)
        results = _run_tasks(hunter, run_cfg)

        if not results and run_cfg.get("dry_run", False):
            return 0

        success_count = sum(1 for r in results if r.success)
        hunter.logger.info(f"批量处理完成: {success_count}/{len(results)} 成功")

        return 0 if success_count == len(results) else 1

    except Exception as e:
        logging.error(f"程序异常: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
