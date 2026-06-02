# pyright: reportAttributeAccessIssue=false, reportPossiblyUnboundVariable=false, reportOperatorIssue=false, reportGeneralTypeIssues=false
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
import numpy as np

from rph_core.utils.log_manager import setup_logger
from rph_core.utils.path_compat import normalize_path, is_toxic_path
from rph_core.utils.task_builder import build_tasks_from_run_config, sanitize_rx_id
from rph_core.utils.config_loader import load_config
from rph_core.utils.small_molecule_catalog import SmallMoleculeCatalog
from rph_core.utils.checkpoint_manager import CheckpointManager
from rph_core.utils.optimization_config import normalize_qc_config
from rph_core.utils.forming_bonds_resolver import resolve_forming_bonds, write_mechanism_meta
from rph_core.utils.file_io import read_xyz
from rph_core.version import __version__
from rph_core.utils import ui, notify
from rph_core.utils.ui import get_progress_manager, SilentProgressManager
from rph_core.utils.task_progress import (
    TaskProgressTracker,
    TaskState,
    V4_TASK_REGISTRY,
)
from rph_core.steps.runners import run_step2, run_step3
from rph_core.utils.layout_contract import resolve_step_dir
from rph_core.utils.json_io import read_json, read_json_dict, write_json

# S4 FeatureMiner has been moved to RPH_Postprocess/rph_features/
# S4 is no longer part of the DFT pipeline. Use rph-features CLI externally:
#   rph-features extract --rph-run <work_dir> --output <features_dir>
FeatureMiner = None


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
    # V7.1: S3-optimized intermediate data for feature extraction
    s3_intermediate_xyz: Optional[Path] = None
    s3_intermediate_l2_energy: Optional[float] = None

    # Error tracking
    error_step: Optional[str] = None
    error_message: Optional[str] = None

    def __str__(self):
        if self.success:
            l2_info = f", L2: {self.e_product_l2:.6f} Ha" if self.e_product_l2 else ""
            return f"✅ Pipeline 成功: {self.product_smiles}\n" \
                   f"   Product: {self.product_xyz}{l2_info}\n" \
                   f"   TS Final: {self.ts_final_xyz}"
        else:
            return f"❌ Pipeline 失败: {self.error_step}\n" \
                   f"   错误: {self.error_message}"


class ReactionProfileHunter:
    _header_printed: bool = False
    _init_log_printed: bool = False

    """Reaction Profile Hunter v3.0 (分子自治架构)

    设计模式:
    - Orchestrator 只是"调度员"，不是"工人"
    - 每个 Step 是独立的"工人"，有明确的输入输出
    - v3.0 核心改进：分子自治目录结构 + OPT-SP 耦合循环

    串行流程:
    S1: AnchorPhase (分子锚定 + CREST + DFT OPT-SP 耦合)
      → S2: RetroScanner (逆向扫描，从 S1_ConfGeneration/[Molecule]/finalDFT 读取)
      → S3: TSOptimizer (TS优化，使用 LogParser 提取坐标)

    S4 feature extraction has been moved to RPH_Postprocess/rph_features/

    v3.0 目录结构:
    S1_ConfGeneration/[Molecule_Name]/
        ├── crest/          # CREST 搜索结果
        ├── xtb/            # 两阶段 xTB 子目录（ext 协议）
        ├── cluster/        # 聚类输出
        ├── prescan/        # full 协议快速 SP 预筛选
        ├── fastsp/         # full/lite 协议快速 SP 筛选
        ├── finalDFT/       # DFT OPT + SP (无子目录，扁平结构)
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

        if not ReactionProfileHunter._header_printed:
            ui.print_pipeline_header(__version__)
            ReactionProfileHunter._header_printed = True
        if not ReactionProfileHunter._init_log_printed:
            self.logger.info(f"Reaction Profile Hunter v{__version__} 初始化 (含 S3.5)")
            ReactionProfileHunter._init_log_printed = True
        
        # Initialize small molecule catalog
        self.logger.info("Initializing SmallMoleculeCatalog...")
        self.small_mol_catalog = SmallMoleculeCatalog(self.config)
        self.logger.info("Initialization complete.")

        # 延迟初始化各步骤引擎（懒加载）
        self._s0_engine = None  # S0: 机理分类
        self._s1_engine = None
        self._s1_engine_type = "ProductAnchor"  # 默认类型
        self._s2_engine = None
        self._s3_engine = None
        self._s4_engine = None

    @property
    def s0_engine(self):
        """Step 0 引擎 (懒加载) - 机理分类"""
        if self._s0_engine is None:
            try:
                from rph_core.steps.mechanism_classifier import MechanismClassifier
                s0_config = self.config.get("s0", {})
                self._s0_engine = MechanismClassifier(s0_config)
                self.logger.debug("✓ Step 0 (MechanismClassifier) 已初始化")
            except ImportError as e:
                self.logger.warning(f"S0 机理分类模块不可用: {e}")
                self._s0_engine = None
        return self._s0_engine

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
        """Step 4 has been moved to RPH_Postprocess/rph_features/"""
        raise ImportError(
            "Step 4 feature extraction has been moved to RPH_Postprocess. "
            "Use the standalone CLI:\n"
            "  rph-features extract --rph-run <work_dir> --output <features_dir>\n"
            "Or install: pip install -e /path/to/RPH_Postprocess/rph_features"
        )

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

    def _resolve_product_xyz_for_s2(self, product_xyz: Path) -> Path:
        product_xyz = Path(product_xyz)
        if not product_xyz.exists():
            raise FileNotFoundError(f"S2 输入产物路径不存在: {product_xyz}")

        if product_xyz.is_file():
            return product_xyz

        candidates = [
            product_xyz / "product_min.xyz",
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

    @staticmethod
    def _is_truthy_flag(value: Any) -> bool:
        if isinstance(value, bool):
            return value
        if value is None:
            return False
        return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}

    def _build_smarts_fallback_context(
        self,
        cleaner_data: Optional[Dict[str, Any]],
    ) -> Optional[Dict[str, Any]]:
        if not cleaner_data:
            return None

        raw_obj = cleaner_data.get("raw")
        raw: Dict[str, Any] = dict(raw_obj) if isinstance(raw_obj, dict) else {}

        reaction_type = self._extract_reaction_type_token(
            cleaner_data.get("reaction_type")
            or cleaner_data.get("rxn_type")
            or cleaner_data.get("reaction_family")
            or cleaner_data.get("reaction_profile")
            or raw.get("reaction_type")
            or raw.get("rxn_type")
            or raw.get("reaction_family")
            or raw.get("reaction_profile")
        )
        reaction_profile = (
            cleaner_data.get("reaction_profile")
            or raw.get("reaction_profile")
            or reaction_type
        )

        context: Dict[str, Any] = {}
        if reaction_type:
            context["reaction_type"] = reaction_type
        if reaction_profile:
            context["reaction_profile"] = str(reaction_profile)

        raw_context: Dict[str, Any] = {}
        if reaction_type:
            raw_context["reaction_type"] = reaction_type
        if reaction_profile:
            raw_context["reaction_profile"] = str(reaction_profile)
        if raw_context:
            context["raw"] = raw_context

        return context or None

    def _derive_forming_from_order_changed(
        self,
        cleaner_data: Optional[Dict[str, Any]],
        product_xyz_file: Path,
    ) -> Optional[Tuple[Tuple[int, int], ...]]:
        """When all 'formed' bonds are already present in the product geometry,
        derive executable forming bonds from core_bond_changes order_changed entries
        that represent actual bond-order transitions (e.g. AROMATIC→SINGLE).
        """
        if not cleaner_data:
            return None

        cbc_raw = cleaner_data.get("core_bond_changes") or (cleaner_data.get("raw") or {}).get("core_bond_changes")
        if not cbc_raw or not isinstance(cbc_raw, str):
            return None

        from rph_core.utils.file_io import read_xyz
        try:
            coords, _ = read_xyz(product_xyz_file)
        except Exception:
            return None

        import numpy as np
        candidates: list[tuple[int, int]] = []
        for item in cbc_raw.split(";"):
            item = item.strip()
            if ":" not in item:
                continue
            bond_part, change_type = item.split(":", 1)
            if not change_type.strip().startswith("order_changed"):
                continue
            if "-" not in bond_part:
                continue
            try:
                a, b = bond_part.split("-")
                a, b = int(a.strip()), int(b.strip())
            except (ValueError, IndexError):
                continue
            # Only consider bonds where the atoms are NOT already bonded in product
            # (distance > 1.6 Å means not a single bond)
            if a < len(coords) and b < len(coords):
                dist = float(np.linalg.norm(np.array(coords[a]) - np.array(coords[b])))
                if dist > 1.6:
                    candidates.append((a, b))
                    self.logger.info(
                        f"[S2] order_changed candidate ({a},{b}): dist={dist:.3f} Å — not bonded in product"
                    )
                else:
                    self.logger.debug(
                        f"[S2] order_changed ({a},{b}): dist={dist:.3f} Å — already bonded, skip"
                    )

        if len(candidates) >= 2:
            result = tuple(candidates[:2])
            self.logger.info(f"[S2] Derived forming_bonds from order_changed: {result}")
            return result

        return None

    def _validate_forming_bonds_against_product(
        self,
        forming_bonds: List[List[int]] | Tuple[Tuple[int, int], ...],
        product_xyz_path: Path,
        index_base: int,
    ) -> Tuple[Tuple[Tuple[int, int], ...], bool]:
        import numpy as np

        validated_bonds = tuple((int(pair[0]), int(pair[1])) for pair in (forming_bonds or []))

        if not product_xyz_path or not Path(product_xyz_path).exists():
            self.logger.warning("[S2] Cannot validate forming_bonds: product_xyz not found")
            return validated_bonds, False

        try:
            coords, _ = read_xyz(Path(product_xyz_path))
        except Exception as exc:
            self.logger.warning(f"[S2] Cannot read product_xyz for validation: {exc}")
            return validated_bonds, False

        if coords is None or len(coords) == 0:
            return validated_bonds, False

        suspicious_count = 0
        zero_based = int(index_base) == 0

        for pair in validated_bonds:
            i = int(pair[0]) if zero_based else int(pair[0]) - 1
            j = int(pair[1]) if zero_based else int(pair[1]) - 1

            if i < 0 or j < 0 or i >= len(coords) or j >= len(coords):
                self.logger.warning(
                    f"[S2] Forming bond {pair} indices out of range for product XYZ ({len(coords)} atoms)"
                )
                suspicious_count += 1
                continue

            dist = float(np.linalg.norm(np.array(coords[i]) - np.array(coords[j])))

            if dist < 1.2:
                self.logger.warning(
                    f"[S2] Forming bond {pair}: distance={dist:.3f} Å in product — "
                    "anomalously close, possibly duplicate atom or data error"
                )
                suspicious_count += 1
            elif dist > 3.5:
                self.logger.warning(
                    f"[S2] Forming bond {pair}: distance={dist:.3f} Å in product — "
                    "atoms very distant, possibly wrong pair / possible index mapping error"
                )
                suspicious_count += 1
            else:
                self.logger.info(f"[S2] Forming bond {pair}: distance={dist:.3f} Å in product — OK")

        # Trigger fallback if ANY forming bond is suspicious (not all)
        needs_fallback = suspicious_count > 0 and len(validated_bonds) > 0
        if needs_fallback:
            self.logger.warning(
                f"[S2] {suspicious_count}/{len(validated_bonds)} forming bonds failed "
                "product-side validation — may need SMARTS fallback"
            )

        return validated_bonds, needs_fallback

    def _finalize_forming_bonds_for_s2(
        self,
        *,
        forming_bonds: Tuple[Tuple[int, int], ...],
        product_xyz_file: Optional[Path],
        cleaner_data: Optional[Dict[str, Any]],
        source_label: str,
        allow_smarts_fallback: bool = True,
    ) -> Tuple[Tuple[int, int], ...]:
        resolved = tuple((int(pair[0]), int(pair[1])) for pair in forming_bonds)
        final_source_label = source_label

        if resolved and product_xyz_file is not None:
            resolved, needs_fallback = self._validate_forming_bonds_against_product(
                forming_bonds=resolved,
                product_xyz_path=Path(product_xyz_file),
                index_base=0,
            )
            if needs_fallback and allow_smarts_fallback:
                self.logger.warning(
                    f"[S2] {source_label} forming_bonds failed validation, attempting SMARTS fallback"
                )
                try:
                    from rph_core.steps.step2_retro.smarts_matcher import SMARTSMatcher

                    matcher = SMARTSMatcher()
                    smarts_result = matcher.find_reactive_bonds(
                        product_xyz=Path(product_xyz_file),
                        cleaner_data=self._build_smarts_fallback_context(cleaner_data),
                    )
                    if smarts_result.matched and smarts_result.bond_1 and smarts_result.bond_2:
                        smarts_bonds = self._normalize_forming_bonds(
                            (
                                (smarts_result.bond_1.atom_idx_1, smarts_result.bond_1.atom_idx_2),
                                (smarts_result.bond_2.atom_idx_1, smarts_result.bond_2.atom_idx_2),
                            ),
                            index_base=0,
                            require_exact_two=True,
                        )
                        self.logger.info(f"[S2] SMARTS fallback forming_bonds: {smarts_bonds}")
                        resolved = smarts_bonds
                        final_source_label = f"SMARTS fallback after {source_label}"
                    elif getattr(smarts_result, "error_message", None):
                        self.logger.warning(f"[S2] SMARTS fallback failed: {smarts_result.error_message}")

                    # If SMARTS failed, try deriving from core_bond_changes order_changed
                    if needs_fallback and cleaner_data:
                        oc_bonds = self._derive_forming_from_order_changed(
                            cleaner_data, product_xyz_file
                        )
                        if oc_bonds:
                            resolved = oc_bonds
                            final_source_label = f"order_changed-derived after {source_label}"
                except Exception as exc:
                    self.logger.warning(f"[S2] SMARTS fallback failed: {exc}")

        self.logger.info(f"[S2] Using {final_source_label} forming_bonds: {resolved}")
        return resolved

    def _resolve_forming_bonds_for_s2(
        self,
        cleaner_data: Optional[Dict[str, Any]] = None,
        product_xyz_file: Optional[Path] = None,
        work_dir: Optional[Path] = None,
    ) -> Tuple[Tuple[int, int], ...]:
        atom_count: Optional[int] = None
        if product_xyz_file is not None and Path(product_xyz_file).exists():
            try:
                coords, _ = read_xyz(Path(product_xyz_file))
                atom_count = int(len(coords))
            except Exception as exc:
                self.logger.warning(f"[S2] Failed to read atom_count from {product_xyz_file}: {exc}")

        deferred_s0_bonds: Tuple[Tuple[int, int], ...] = tuple()

        if work_dir is not None:
            s0_dir = resolve_step_dir(work_dir, "s0")
            s1_dir = resolve_step_dir(work_dir, "s1")
            mechanism_graph_payload = self._load_json_artifact(s0_dir / "mechanism_graph.json") or {}
            smiles_mapping_payload = self._load_json_artifact(s0_dir / "atom_map_smiles.json") or {}
            xyz_mapping_payload = self._load_json_artifact(s1_dir / "atom_map_xyz.json") or {}

            artifact_resolvers = [
                (
                    "S1 atom_map_xyz full notation",
                    lambda: self._resolve_forming_bonds_from_xyz_mapping_artifact(
                        xyz_mapping_payload,
                        atom_count=atom_count,
                    ),
                ),
                (
                    "S0/S1 atom-map annotations",
                    lambda: self._resolve_forming_bonds_from_mapping_annotations(
                        smiles_mapping_payload,
                        mechanism_graph_payload,
                        xyz_mapping_payload,
                        atom_count=atom_count,
                    ),
                ),
                (
                    "S0 mechanism_graph primary edges + S1 atom_map_xyz",
                    lambda: self._resolve_forming_bonds_from_graph_edges_and_xyz_mapping(
                        smiles_mapping_payload,
                        mechanism_graph_payload,
                        xyz_mapping_payload,
                        atom_count=atom_count,
                    ),
                ),
            ]

            for source_label, resolver in artifact_resolvers:
                try:
                    artifact_bonds = resolver()
                except Exception as exc:
                    self.logger.debug(f"[S2] Failed to resolve forming_bonds from {source_label}: {exc}")
                    continue

                if artifact_bonds:
                    return self._finalize_forming_bonds_for_s2(
                        forming_bonds=artifact_bonds,
                        product_xyz_file=product_xyz_file,
                        cleaner_data=cleaner_data,
                        source_label=source_label,
                    )

            s0_summary_payload = self._load_json_artifact(s0_dir / "mechanism_summary.json") or {}
            s0_bonds_raw = s0_summary_payload.get("forming_bonds")
            if s0_bonds_raw:
                try:
                    s0_bonds = self._normalize_forming_bonds(
                        s0_bonds_raw,
                        atom_count=atom_count,
                        index_base=0,
                        require_exact_two=True,
                    )
                except Exception as exc:
                    self.logger.debug(f"[S2] Failed to parse S0 mechanism summary forming_bonds: {exc}")
                else:
                    if s0_bonds:
                        if self._is_truthy_flag(s0_summary_payload.get("low_confidence")):
                            self.logger.warning(
                                "[S2] S0 mechanism summary marked low_confidence; deferring legacy summary fallback"
                            )
                            deferred_s0_bonds = s0_bonds
                        else:
                            return self._finalize_forming_bonds_for_s2(
                                forming_bonds=s0_bonds,
                                product_xyz_file=product_xyz_file,
                                cleaner_data=cleaner_data,
                                source_label="S0 mechanism summary legacy fallback",
                            )

        if product_xyz_file is not None and Path(product_xyz_file).exists():
            try:
                from rph_core.steps.step2_retro.smarts_matcher import SMARTSMatcher
                matcher = SMARTSMatcher()
                smarts_result = matcher.find_reactive_bonds(
                    product_xyz=Path(product_xyz_file),
                    cleaner_data=self._build_smarts_fallback_context(cleaner_data),
                )
                if smarts_result.matched and smarts_result.bond_1 and smarts_result.bond_2:
                    auto_bonds = (
                        (smarts_result.bond_1.atom_idx_1, smarts_result.bond_1.atom_idx_2),
                        (smarts_result.bond_2.atom_idx_1, smarts_result.bond_2.atom_idx_2),
                    )
                    return self._finalize_forming_bonds_for_s2(
                        forming_bonds=auto_bonds,
                        product_xyz_file=product_xyz_file,
                        cleaner_data=cleaner_data,
                        source_label="SMARTS-derived",
                        allow_smarts_fallback=False,
                    )
            except Exception as exc:
                self.logger.debug(f"[S2] SMARTS auto-detection failed: {exc}")

        if deferred_s0_bonds:
            return self._finalize_forming_bonds_for_s2(
                forming_bonds=deferred_s0_bonds,
                product_xyz_file=product_xyz_file,
                cleaner_data=cleaner_data,
                source_label="deferred low-confidence S0 mechanism summary",
                allow_smarts_fallback=False,
            )

        s0_dir_path = str(s0_dir) if work_dir is not None else "<no work_dir>"
        raise RuntimeError(
            f"S2 requires forming_bonds from S0/S1 artifacts or SMARTS fallback; "
            f"all sources failed. S0 dir checked: {s0_dir_path}, "
            f"S0 artifacts exist: {(s0_dir / 'mechanism_summary.json').exists() if work_dir else False}, "
            f"product_xyz: {product_xyz_file}"
        )

    def _resolve_small_molecular_keys(
        self,
        reaction_profile: Optional[str],
        dataset_keys: Optional[list[str]],
        include_reference_terms: bool = True,
    ) -> list[str]:
        """
        Resolve small molecular species keys for S1 conformer search.

        Priority:
        1) reaction_reference_terms config (per reaction profile): extracts
           all non-precursor/non-intermediate species from reactants+products
        2) dataset_keys (from CSV small_molecular columns): supplemental keys

        All keys are validated against SmallMoleculeCatalog.
        """
        from rph_core.utils.small_molecule_catalog import SmallMoleculeCatalog

        catalog = SmallMoleculeCatalog(self.config)
        ref_terms = self.config.get("reaction_reference_terms", {}) or {}
        if not isinstance(ref_terms, dict):
            ref_terms = {}

        profile_key = reaction_profile or "[4+3]_default"
        profile_cfg = ref_terms.get(profile_key, {})
        if not isinstance(profile_cfg, dict):
            profile_cfg = {}

        pti = profile_cfg.get("precursor_to_intermediate", {})
        if not isinstance(pti, dict):
            pti = {}

        reserved = {"precursor", "intermediate"}
        auto_keys: list[str] = []

        if include_reference_terms:
            for side in ("reactants", "products"):
                side_dict = pti.get(side, {})
                if isinstance(side_dict, dict):
                    for key in side_dict:
                        if key not in reserved and key not in auto_keys:
                            auto_keys.append(key)

        merged_keys: list[str] = list(auto_keys)
        if dataset_keys:
            for key in dataset_keys:
                if key and key not in merged_keys:
                    merged_keys.append(key)

        validated: list[str] = []
        for key in merged_keys:
            if catalog.get(key) is not None:
                validated.append(key)
            else:
                self.logger.warning(
                    f"Small molecular key '{key}' not in catalog, skipping S1 processing"
                )

        if validated:
            self.logger.info(
                f"Resolved small molecular keys for S1: {validated} "
                f"(auto={auto_keys}, dataset={dataset_keys or []})"
            )
        return validated

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

    def _extract_reaction_type_token(self, value: Any) -> Optional[str]:
        if value is None:
            return None

        text = str(value).strip()
        if not text:
            return None

        compact = text.replace(" ", "")
        if compact.endswith("_default"):
            compact = compact[:-8]

        bracket_match = re.search(r"\[(\d+\+\d+)\]", compact)
        if bracket_match:
            return f"[{bracket_match.group(1)}]"

        plain_match = re.search(r"(\d+\+\d+)", compact)
        if plain_match:
            return f"[{plain_match.group(1)}]"

        return compact or None

    def _normalize_cleaner_data_for_pipeline(
        self,
        cleaner_data: Optional[Dict[str, Any]],
        *,
        product_smiles: str,
        precursor_smiles: Optional[str] = None,
        reaction_profile: Optional[str] = None,
    ) -> Optional[Dict[str, Any]]:
        if not cleaner_data:
            return None

        row = dict(cleaner_data)
        raw_obj = row.get("raw")
        raw: Dict[str, Any] = dict(raw_obj) if isinstance(raw_obj, dict) else {}

        def _pick_nonempty(*values: Any) -> Optional[str]:
            for value in values:
                if value is None:
                    continue
                text = str(value).strip()
                if text:
                    return text
            return None

        rxn_key_hash = _pick_nonempty(
            row.get("rxn_key_hash"),
            row.get("reaction_id"),
            row.get("rx_id"),
            row.get("record_id"),
            row.get("id"),
            raw.get("rxn_key_hash"),
            raw.get("reaction_id"),
            raw.get("rx_id"),
            raw.get("record_id"),
            raw.get("id"),
        )
        if rxn_key_hash:
            row["rxn_key_hash"] = rxn_key_hash
            raw.setdefault("rxn_key_hash", rxn_key_hash)

        precursor = _pick_nonempty(
            row.get("precursor_smiles"),
            row.get("reactant_smiles"),
            row.get("substrate_smiles"),
            raw.get("precursor_smiles"),
            raw.get("reactant_smiles"),
            raw.get("substrate_smiles"),
            precursor_smiles,
        )
        if precursor:
            row["precursor_smiles"] = precursor
            raw.setdefault("precursor_smiles", precursor)

        product_main = _pick_nonempty(
            row.get("product_smiles_main"),
            row.get("product_smiles"),
            raw.get("product_smiles_main"),
            raw.get("product_smiles"),
            product_smiles,
        )
        if product_main:
            row["product_smiles_main"] = product_main
            raw.setdefault("product_smiles_main", product_main)

        reaction_type = _pick_nonempty(
            self._extract_reaction_type_token(row.get("reaction_type")),
            self._extract_reaction_type_token(row.get("rxn_type")),
            self._extract_reaction_type_token(raw.get("reaction_type")),
            self._extract_reaction_type_token(raw.get("rxn_type")),
            self._extract_reaction_type_token(row.get("reaction_family")),
            self._extract_reaction_type_token(raw.get("reaction_family")),
            self._extract_reaction_type_token(row.get("reaction_profile")),
            self._extract_reaction_type_token(raw.get("reaction_profile")),
            self._extract_reaction_type_token(reaction_profile),
        )
        if reaction_type:
            row["reaction_type"] = reaction_type
            raw.setdefault("reaction_type", reaction_type)

        profile_text = _pick_nonempty(
            row.get("reaction_profile"),
            raw.get("reaction_profile"),
            reaction_profile,
        )
        if profile_text:
            row["reaction_profile"] = profile_text
            raw.setdefault("reaction_profile", profile_text)

        row["raw"] = raw
        return row

    def _get_cleaner_value(
        self,
        cleaner_data: Optional[Dict[str, Any]],
        *keys: str,
    ) -> Optional[str]:
        if not cleaner_data:
            return None

        raw = cleaner_data.get("raw", {}) or {}
        for source in (cleaner_data, raw):
            if not isinstance(source, dict):
                continue
            for key in keys:
                value = source.get(key)
                if value is None:
                    continue
                text = str(value).strip()
                if text:
                    return text
        return None

    def _load_json_artifact(self, path: Path) -> Optional[Dict[str, Any]]:
        payload = read_json(path, default=None)
        return payload if isinstance(payload, dict) else None

    def _normalize_int_mapping(self, payload: Any) -> Dict[int, int]:
        if not isinstance(payload, dict):
            return {}

        normalized: Dict[int, int] = {}
        for key, value in payload.items():
            try:
                normalized[int(key)] = int(value)
            except (TypeError, ValueError):
                continue
        return normalized

    def _resolve_primary_graph_forming_bonds(
        self,
        mechanism_graph_payload: Optional[Dict[str, Any]],
    ) -> Tuple[Tuple[int, int], ...]:
        if not isinstance(mechanism_graph_payload, dict):
            return tuple()

        edges_payload = mechanism_graph_payload.get("edges")
        if not isinstance(edges_payload, list):
            return tuple()

        primary_pairs: List[Tuple[int, int]] = []
        fallback_pairs: List[Tuple[int, int]] = []
        for edge in edges_payload:
            if not isinstance(edge, dict):
                continue
            raw_pairs = edge.get("forming_bonds")
            if not raw_pairs:
                continue

            target_pairs = fallback_pairs
            pathway_id = str(edge.get("pathway_id") or "").strip().lower()
            if pathway_id in {"", "primary"}:
                target_pairs = primary_pairs

            for pair in raw_pairs:
                if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                    continue
                try:
                    target_pairs.append((int(pair[0]), int(pair[1])))
                except (TypeError, ValueError):
                    continue

        graph_pairs = primary_pairs if primary_pairs else fallback_pairs
        if not graph_pairs:
            return tuple()

        return self._normalize_forming_bonds(
            graph_pairs,
            index_base=0,
            require_exact_two=True,
        )

    def _resolve_forming_bonds_from_xyz_mapping_artifact(
        self,
        xyz_mapping_payload: Optional[Dict[str, Any]],
        *,
        atom_count: Optional[int],
    ) -> Tuple[Tuple[int, int], ...]:
        if not isinstance(xyz_mapping_payload, dict):
            return tuple()

        full_notation = xyz_mapping_payload.get("forming_bonds_full_notation")
        if not isinstance(full_notation, list):
            return tuple()

        xyz_pairs_1based: List[Tuple[int, int]] = []
        for entry in full_notation:
            if not isinstance(entry, dict):
                continue
            xyz_pair = entry.get("product_xyz_1based")
            if not isinstance(xyz_pair, (list, tuple)) or len(xyz_pair) != 2:
                continue
            if xyz_pair[0] is None or xyz_pair[1] is None:
                continue
            try:
                xyz_pairs_1based.append((int(xyz_pair[0]), int(xyz_pair[1])))
            except (TypeError, ValueError):
                continue

        if not xyz_pairs_1based:
            return tuple()

        return self._normalize_forming_bonds(
            xyz_pairs_1based,
            atom_count=atom_count,
            index_base=1,
            require_exact_two=True,
        )

    def _resolve_forming_bonds_from_mapping_annotations(
        self,
        smiles_mapping_payload: Optional[Dict[str, Any]],
        mechanism_graph_payload: Optional[Dict[str, Any]],
        xyz_mapping_payload: Optional[Dict[str, Any]],
        *,
        atom_count: Optional[int],
    ) -> Tuple[Tuple[int, int], ...]:
        if not isinstance(xyz_mapping_payload, dict):
            return tuple()

        map_to_product_xyz = self._normalize_int_mapping(xyz_mapping_payload.get("map_to_product_xyz_1based"))
        if not map_to_product_xyz:
            return tuple()

        annotations: Any = None
        if isinstance(smiles_mapping_payload, dict):
            annotations = smiles_mapping_payload.get("forming_bonds_annotated")
        if not isinstance(annotations, list) and isinstance(mechanism_graph_payload, dict):
            annotations = mechanism_graph_payload.get("forming_bonds_annotated")
        if not isinstance(annotations, list):
            return tuple()

        xyz_pairs_1based: List[Tuple[int, int]] = []
        for entry in annotations:
            if not isinstance(entry, dict):
                continue
            map_pair = entry.get("map_space")
            if not isinstance(map_pair, (list, tuple)) or len(map_pair) != 2:
                continue
            try:
                left_map = int(map_pair[0])
                right_map = int(map_pair[1])
            except (TypeError, ValueError):
                continue

            left_xyz = map_to_product_xyz.get(left_map)
            right_xyz = map_to_product_xyz.get(right_map)
            if left_xyz is None or right_xyz is None:
                continue
            xyz_pairs_1based.append((left_xyz, right_xyz))

        if not xyz_pairs_1based:
            return tuple()

        return self._normalize_forming_bonds(
            xyz_pairs_1based,
            atom_count=atom_count,
            index_base=1,
            require_exact_two=True,
        )

    def _resolve_forming_bonds_from_graph_edges_and_xyz_mapping(
        self,
        smiles_mapping_payload: Optional[Dict[str, Any]],
        mechanism_graph_payload: Optional[Dict[str, Any]],
        xyz_mapping_payload: Optional[Dict[str, Any]],
        *,
        atom_count: Optional[int],
    ) -> Tuple[Tuple[int, int], ...]:
        graph_pairs = self._resolve_primary_graph_forming_bonds(mechanism_graph_payload)
        if not graph_pairs or not isinstance(xyz_mapping_payload, dict):
            return tuple()

        map_to_product_xyz = self._normalize_int_mapping(xyz_mapping_payload.get("map_to_product_xyz_1based"))
        if not map_to_product_xyz:
            return tuple()

        smiles_atom_mapping: Any = None
        if isinstance(smiles_mapping_payload, dict):
            smiles_atom_mapping = smiles_mapping_payload.get("smiles_atom_mapping")
        if not isinstance(smiles_atom_mapping, dict) and isinstance(mechanism_graph_payload, dict):
            smiles_atom_mapping = mechanism_graph_payload.get("smiles_atom_mapping")
        if not isinstance(smiles_atom_mapping, dict):
            smiles_atom_mapping = {}

        product_smiles_to_map = self._normalize_int_mapping(smiles_atom_mapping.get("product_smiles_to_map"))

        via_smiles_mapping: List[Tuple[int, int]] = []
        if product_smiles_to_map:
            for left_idx, right_idx in graph_pairs:
                left_map = product_smiles_to_map.get(int(left_idx))
                right_map = product_smiles_to_map.get(int(right_idx))
                if left_map is None or right_map is None:
                    via_smiles_mapping = []
                    break

                left_xyz = map_to_product_xyz.get(left_map)
                right_xyz = map_to_product_xyz.get(right_map)
                if left_xyz is None or right_xyz is None:
                    via_smiles_mapping = []
                    break
                via_smiles_mapping.append((left_xyz, right_xyz))

        if via_smiles_mapping:
            return self._normalize_forming_bonds(
                via_smiles_mapping,
                atom_count=atom_count,
                index_base=1,
                require_exact_two=True,
            )

        via_direct_map_ids: List[Tuple[int, int]] = []
        for left_idx, right_idx in graph_pairs:
            left_xyz = map_to_product_xyz.get(int(left_idx))
            right_xyz = map_to_product_xyz.get(int(right_idx))
            if left_xyz is None or right_xyz is None:
                via_direct_map_ids = []
                break
            via_direct_map_ids.append((left_xyz, right_xyz))

        if not via_direct_map_ids:
            return tuple()

        return self._normalize_forming_bonds(
            via_direct_map_ids,
            atom_count=atom_count,
            index_base=1,
            require_exact_two=True,
        )

    def _build_smiles_atom_mapping_payload(
        self,
        precursor_smiles: Optional[str],
        product_smiles: Optional[str],
    ) -> Dict[str, Dict[int, int]]:
        mapping: Dict[str, Dict[int, int]] = {
            "precursor_smiles_to_map": {},
            "product_smiles_to_map": {},
            "map_to_precursor_smiles": {},
            "map_to_product_smiles": {},
        }

        try:
            from rdkit import Chem
        except ImportError:
            return mapping

        for smiles, forward_key, reverse_key in (
            (precursor_smiles, "precursor_smiles_to_map", "map_to_precursor_smiles"),
            (product_smiles, "product_smiles_to_map", "map_to_product_smiles"),
        ):
            if not smiles:
                continue
            try:
                mol = Chem.MolFromSmiles(smiles)
            except Exception:
                mol = None
            if mol is None:
                continue

            for atom in mol.GetAtoms():
                map_num = int(atom.GetAtomMapNum())
                if map_num <= 0:
                    continue
                idx = int(atom.GetIdx())
                mapping[forward_key][idx] = map_num
                mapping[reverse_key][map_num] = idx

        return mapping

    @staticmethod
    def _has_atom_map_numbers(smiles: str) -> bool:
        """True if SMILES contains atom map indicators like [C:1] or [N:2]."""
        return bool(re.search(r"\[[A-Z][a-z]?:\d+\]", smiles))

    @staticmethod
    def _kabsch_align(P: "np.ndarray", Q: "np.ndarray") -> "np.ndarray":
        """Optimal rotation alignment of point cloud P to Q (Kabsch algorithm)."""
        import numpy as np

        H = P.T @ Q
        U, _, Vt = np.linalg.svd(H)
        rotation = Vt.T @ U.T
        if np.linalg.det(rotation) < 0:
            Vt[-1, :] *= -1
            rotation = Vt.T @ U.T
        return P @ rotation

    @staticmethod
    def _hungarian_match(
        P: "np.ndarray",
        Q: "np.ndarray",
    ) -> "tuple[np.ndarray, np.ndarray]":
        """Optimal 1-to-1 bipartite matching minimising total Euclidean distance."""
        import numpy as np
        from scipy.optimize import linear_sum_assignment

        cost = np.zeros((len(P), len(Q)), dtype=float)
        for i in range(len(P)):
            diff = Q - P[i]
            cost[i, :] = np.sqrt(np.sum(diff * diff, axis=1))
        return linear_sum_assignment(cost)

    def _build_mol_to_xyz_index_mapping(
        self,
        mapped_smiles: Optional[str],
        xyz_path: Optional[Path],
    ) -> Dict[int, int]:
        if not mapped_smiles or xyz_path is None or not Path(xyz_path).exists():
            return {}

        try:
            import numpy as np
            from rdkit import Chem
            from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
        except ImportError:
            self.logger.warning("RDKit or numpy unavailable; XYZ atom mapping skipped")
            return {}

        try:
            mol = Chem.MolFromSmiles(mapped_smiles)
        except Exception as exc:
            self.logger.warning(f"Failed to parse mapped SMILES for XYZ mapping: {exc}")
            return {}
        if mol is None:
            return {}

        mol_with_coords = Chem.Mol(mol)
        embed_status = -1
        try:
            embed_status = rdDistGeom.EmbedMolecule(mol_with_coords, randomSeed=0xF00D)
            if embed_status != 0:
                embed_status = rdDistGeom.EmbedMolecule(
                    mol_with_coords,
                    randomSeed=0xF00D,
                    useRandomCoords=True,
                )
            if embed_status == 0:
                try:
                    rdForceFieldHelpers.UFFOptimizeMolecule(mol_with_coords, maxIters=200)
                except Exception:
                    pass
        except Exception:
            embed_status = -1

        if embed_status == 0:
            try:
                xyz_coords, xyz_symbols = read_xyz(Path(xyz_path))
                conf = mol_with_coords.GetConformer()
                rd_coords = np.array(
                    [
                        [
                            float(conf.GetAtomPosition(i).x),
                            float(conf.GetAtomPosition(i).y),
                            float(conf.GetAtomPosition(i).z),
                        ]
                        for i in range(mol_with_coords.GetNumAtoms())
                    ],
                    dtype=float,
                )
                xyz_coords = np.asarray(xyz_coords, dtype=float)
                if rd_coords.size > 0 and xyz_coords.size > 0:
                    rd_centered = rd_coords - rd_coords.mean(axis=0)
                    xyz_centered = xyz_coords - xyz_coords.mean(axis=0)
                    rd_symbols = [atom.GetSymbol() for atom in mol_with_coords.GetAtoms()]

                    rd_aligned = self._kabsch_align(rd_centered, xyz_centered)

                    mol_to_xyz: Dict[int, int] = {}
                    elements = sorted(set(rd_symbols))
                    for element in elements:
                        rd_indices = [i for i, s in enumerate(rd_symbols) if s == element]
                        xyz_indices = [j for j, s in enumerate(xyz_symbols) if s == element]
                        if len(rd_indices) > len(xyz_indices):
                            self.logger.debug(
                                f"XYZ has fewer {element} atoms ({len(xyz_indices)}) than RDKit ({len(rd_indices)}); "
                                f"falling back to greedy match for this element"
                            )
                            used_xyz = set()
                            for rd_i in rd_indices:
                                candidates = [j for j in xyz_indices if j not in used_xyz]
                                if not candidates:
                                    break
                                best = min(candidates, key=lambda j: float(np.linalg.norm(xyz_centered[j] - rd_aligned[rd_i])))
                                mol_to_xyz[rd_i] = best
                                used_xyz.add(best)
                            continue

                        P = rd_aligned[rd_indices]
                        Q = xyz_centered[xyz_indices]
                        row_ind, col_ind = self._hungarian_match(P, Q)
                        for r, c in zip(row_ind, col_ind):
                            mol_to_xyz[int(rd_indices[r])] = int(xyz_indices[c])

                    if len(mol_to_xyz) == mol_with_coords.GetNumAtoms():
                        return mol_to_xyz
            except Exception as exc:
                self.logger.debug(f"Coordinate-based XYZ mapping fallback triggered: {exc}")

        try:
            from rdkit import Chem
            from rph_core.utils.cleaner_adapter import get_map_to_xyz_dict

            fallback_map_to_xyz = get_map_to_xyz_dict(mapped_smiles, Path(xyz_path))
            if not fallback_map_to_xyz:
                return {}

            mol_for_maps = Chem.MolFromSmiles(mapped_smiles)
            if mol_for_maps is None:
                return {}

            map_to_mol_idx = {
                int(atom.GetAtomMapNum()): int(atom.GetIdx())
                for atom in mol_for_maps.GetAtoms()
                if int(atom.GetAtomMapNum()) > 0
            }
            return {
                int(mol_idx): int(fallback_map_to_xyz[map_num])
                for map_num, mol_idx in map_to_mol_idx.items()
                if map_num in fallback_map_to_xyz
            }
        except Exception as exc:
            self.logger.warning(f"Failed to build mol↔XYZ atom mapping for {xyz_path}: {exc}")
            return {}

    def _build_and_save_xyz_mapping(
        self,
        *,
        work_dir: Path,
        cleaner_data: Optional[Dict[str, Any]],
        product_xyz: Optional[Path],
    ) -> Optional[Path]:
        s1_dir = resolve_step_dir(work_dir, "s1")
        s1_dir.mkdir(parents=True, exist_ok=True)

        mapping_path = s1_dir / "atom_map_xyz.json"
        s0_dir = resolve_step_dir(work_dir, "s0")

        product_xyz_path: Optional[Path] = None
        if product_xyz is not None:
            try:
                product_xyz_path = self._resolve_product_xyz_for_s2(Path(product_xyz))
            except Exception as exc:
                self.logger.warning(f"[S1→S2] Failed to resolve product XYZ for mapping: {exc}")

        if product_xyz_path is None or not product_xyz_path.exists():
            for candidate in (s1_dir / "product_min.xyz", s1_dir / "product" / "product_min.xyz"):
                if candidate.exists():
                    product_xyz_path = candidate
                    break

        precursor_xyz_path = self._resolve_s1_artifacts(work_dir).get("s1_precursor_xyz")

        smiles_mapping_payload = self._load_json_artifact(s0_dir / "atom_map_smiles.json") or {}
        mechanism_graph_payload = self._load_json_artifact(s0_dir / "mechanism_graph.json") or {}
        graph_source_data = mechanism_graph_payload.get("source_data")
        graph_source_data = graph_source_data if isinstance(graph_source_data, dict) else {}

        cleaner_context: Dict[str, Any] = {}
        cleaner_raw: Dict[str, Any] = {}
        if graph_source_data:
            cleaner_context.update(graph_source_data)
            cleaner_raw.update(graph_source_data)
        if cleaner_data:
            cleaner_context.update(cleaner_data)
            if isinstance(cleaner_data.get("raw"), dict):
                cleaner_raw.update(cleaner_data.get("raw", {}) or {})
        cleaner_context["raw"] = cleaner_raw

        smiles_atom_mapping = smiles_mapping_payload.get("smiles_atom_mapping")
        if not isinstance(smiles_atom_mapping, dict):
            smiles_atom_mapping = mechanism_graph_payload.get("smiles_atom_mapping")
        if not isinstance(smiles_atom_mapping, dict):
            smiles_atom_mapping = self._build_smiles_atom_mapping_payload(
                self._get_cleaner_value(
                    cleaner_context,
                    "mapped_precursor_smiles",
                    "precursor_smiles_mapped",
                    "mapped_reactant_smiles",
                    "reactant_smiles_mapped",
                    "precursor_smiles",
                ),
                self._get_cleaner_value(
                    cleaner_context,
                    "mapped_product_smiles",
                    "product_smiles_mapped",
                    "mapped_product",
                    "product_smiles_main",
                    "product_smiles",
                ),
            )

        forming_bonds_annotated = smiles_mapping_payload.get("forming_bonds_annotated")
        if not isinstance(forming_bonds_annotated, list):
            forming_bonds_annotated = mechanism_graph_payload.get("forming_bonds_annotated")
        if not isinstance(forming_bonds_annotated, list):
            forming_bonds_annotated = []

        mapped_product_smiles = self._get_cleaner_value(
            cleaner_context,
            "mapped_product_smiles",
            "product_smiles_mapped",
            "mapped_product",
            "product_smiles_main",
            "product_smiles",
        )
        mapped_precursor_smiles = self._get_cleaner_value(
            cleaner_context,
            "mapped_precursor_smiles",
            "precursor_smiles_mapped",
            "mapped_reactant_smiles",
            "reactant_smiles_mapped",
            "precursor_smiles",
        )

        rxn_mapped = self._get_cleaner_value(cleaner_context, "rxn_smiles_mapped", "reaction_smiles_mapped")
        if rxn_mapped and ">>" in rxn_mapped:
            if not mapped_product_smiles or not self._has_atom_map_numbers(mapped_product_smiles):
                product_candidate = rxn_mapped.split(">>")[1].strip()
                if product_candidate and self._has_atom_map_numbers(product_candidate):
                    mapped_product_smiles = product_candidate
            if not mapped_precursor_smiles or not self._has_atom_map_numbers(mapped_precursor_smiles):
                precursor_candidate = rxn_mapped.split(">>")[0].strip()
                if precursor_candidate and self._has_atom_map_numbers(precursor_candidate):
                    mapped_precursor_smiles = precursor_candidate

        product_mol_to_xyz = self._build_mol_to_xyz_index_mapping(mapped_product_smiles, product_xyz_path)
        precursor_mol_to_xyz = self._build_mol_to_xyz_index_mapping(mapped_precursor_smiles, precursor_xyz_path)

        map_to_product_smiles = self._normalize_int_mapping(smiles_atom_mapping.get("map_to_product_smiles"))
        map_to_precursor_smiles = self._normalize_int_mapping(smiles_atom_mapping.get("map_to_precursor_smiles"))

        map_to_product_xyz_1based = {
            int(map_num): int(product_mol_to_xyz[smiles_idx]) + 1
            for map_num, smiles_idx in map_to_product_smiles.items()
            if smiles_idx in product_mol_to_xyz
        }
        map_to_precursor_xyz_1based = {
            int(map_num): int(precursor_mol_to_xyz[smiles_idx]) + 1
            for map_num, smiles_idx in map_to_precursor_smiles.items()
            if smiles_idx in precursor_mol_to_xyz
        }

        forming_bonds_full_notation = []
        for entry in forming_bonds_annotated:
            if not isinstance(entry, dict):
                continue
            map_space = entry.get("map_space")
            if not isinstance(map_space, (list, tuple)) or len(map_space) != 2:
                continue
            try:
                map_pair = [int(map_space[0]), int(map_space[1])]
            except (TypeError, ValueError):
                continue

            product_xyz_pair = [map_to_product_xyz_1based.get(map_pair[0]), map_to_product_xyz_1based.get(map_pair[1])]
            precursor_xyz_pair = [map_to_precursor_xyz_1based.get(map_pair[0]), map_to_precursor_xyz_1based.get(map_pair[1])]

            forming_bonds_full_notation.append(
                {
                    "map_space": map_pair,
                    "product_xyz_1based": product_xyz_pair if all(v is not None for v in product_xyz_pair) else None,
                    "precursor_xyz_1based": precursor_xyz_pair if all(v is not None for v in precursor_xyz_pair) else None,
                    "bond_type_product": str(entry.get("bond_type_product") or "UNKNOWN"),
                    "bond_type_precursor": str(entry.get("bond_type_precursor") or "NONE"),
                }
            )

        payload = {
            "map_to_product_xyz_1based": map_to_product_xyz_1based,
            "map_to_precursor_xyz_1based": map_to_precursor_xyz_1based,
            "forming_bonds_full_notation": forming_bonds_full_notation,
            "index_base_convention": "all XYZ indices are 1-based; map numbers are 1-based",
        }

        with open(mapping_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)

        if not map_to_product_xyz_1based:
            self.logger.warning(f"[S1→S2] Product XYZ atom mapping incomplete: {mapping_path}")
        else:
            self.logger.info(f"[S1→S2] XYZ atom mapping saved: {mapping_path}")

        return mapping_path

    def _run_s0(
        self,
        *,
        work_dir: Path,
        product_smiles: str,
        cleaner_data: Optional[Dict[str, Any]],
        checkpoint_mgr: CheckpointManager,
        resume_enabled: bool,
        pm: Any = None,
    ) -> Optional[Dict[str, Any]]:
        s0_cfg = self.config.get("s0", {}) or {}
        s0_dir = resolve_step_dir(work_dir, "s0")
        s0_dir.mkdir(parents=True, exist_ok=True)

        if pm:
            pm.section_header("Step 0", "Mechanism Classifier", "Building reaction graph from cleaner data")

        def _write_status(
            status: str,
            reason: Optional[str] = None,
            extra: Optional[Dict[str, Any]] = None,
        ) -> None:
            artifact: Dict[str, Any] = {
                "status": status,
                "reason": reason,
                "timestamp": datetime.now().isoformat(),
            }
            if extra:
                artifact.update(extra)
            status_path = s0_dir / "s0_status.json"
            try:
                with open(status_path, "w", encoding="utf-8") as f:
                    json.dump(artifact, f, indent=2)
            except Exception as exc:
                self.logger.debug(f"[S0] Failed to write s0_status.json: {exc}")

        if not bool(s0_cfg.get("enabled", True)):
            _write_status("skipped", "disabled_by_config")
            return {"status": "skipped", "reason": "disabled_by_config"}
        if pm:
            pm.update_step("s0", completed=10, description="Validating cleaner data...")
        if not cleaner_data:
            self.logger.info("[S0] cleaner_data unavailable; skip")
            _write_status("skipped", "cleaner_data_unavailable")
            return {"status": "skipped", "reason": "cleaner_data_unavailable"}
        if self.s0_engine is None:
            self.logger.warning("[S0] MechanismClassifier unavailable; skip")
            _write_status("skipped", "engine_unavailable")
            return {"status": "skipped", "reason": "engine_unavailable"}

        if resume_enabled and checkpoint_mgr.is_step_completed("s0"):
            graph_path = checkpoint_mgr.get_step_output("s0", "mechanism_graph_json")
            summary_path = checkpoint_mgr.get_step_output("s0", "mechanism_summary_json")
            dr_plan_path = checkpoint_mgr.get_step_output("s0", "dr_branch_plan_json")
            if (
                graph_path and Path(graph_path).exists()
                and summary_path and Path(summary_path).exists()
                and dr_plan_path and Path(dr_plan_path).exists()
            ):
                self.logger.info("✅ Resume: Step0 already complete")
                result = {
                    "status": "complete",
                    "resume_reused": True,
                    "mechanism_graph_json": graph_path,
                    "mechanism_summary_json": summary_path,
                    "dr_branch_plan_json": dr_plan_path,
                }
                _write_status("complete", extra=result)
                return result

        checkpoint_mgr.mark_step_in_progress("s0", phase="classifying")
        if pm:
            pm.update_step("s0", completed=30, description="Classifying reaction mechanism...")
        cleaner_row = self._normalize_cleaner_data_for_pipeline(
            cleaner_data,
            product_smiles=product_smiles,
        ) or dict(cleaner_data)
        if not cleaner_row.get("rxn_key_hash"):
            fallback_rxn_id = (
                cleaner_row.get("reaction_id")
                or cleaner_row.get("rx_id")
                or cleaner_row.get("record_id")
                or cleaner_row.get("id")
                or "manual"
            )
            cleaner_row["rxn_key_hash"] = str(fallback_rxn_id)

        graph = self.s0_engine.classify_from_dict(cleaner_row)
        if graph is None:
            from rph_core.steps.mechanism_classifier.models import create_simple_graph

            precursor = str(
                cleaner_row.get("precursor_smiles")
                or cleaner_row.get("ylide_smiles")
                or cleaner_row.get("reactant_smiles")
                or ""
            ).strip()
            reaction_type = str(cleaner_row.get("reaction_type") or cleaner_row.get("rxn_type") or "unknown")
            reaction_id = str(cleaner_row.get("rxn_key_hash") or "manual")
            raw_forming = cleaner_row.get("forming_bonds")
            if not raw_forming:
                raw = cleaner_row.get("raw")
                if isinstance(raw, dict):
                    raw_forming = raw.get("forming_bonds")
            try:
                normalized = self._normalize_forming_bonds(
                    raw_forming,
                    atom_count=None,
                    index_base="auto",
                    require_exact_two=False,
                )
            except Exception:
                normalized = tuple()

            if precursor and normalized:
                graph = create_simple_graph(
                    reaction_id=reaction_id,
                    reaction_type=reaction_type,
                    precursor_smiles=precursor,
                    product_smiles=product_smiles,
                    forming_bonds=list(normalized),
                )
            else:
                _write_status("skipped", "insufficient_cleaner_fields")
                return {"status": "skipped", "reason": "insufficient_cleaner_fields"}
        if graph is None:
            _write_status("skipped", "classification_failed")
            return {"status": "skipped", "reason": "classification_failed"}

        if pm:
            pm.update_step("s0", completed=60, description="Extracting forming bonds...")

        from rph_core.steps.mechanism_classifier.dr_completion import (
            attach_dr_completion,
            build_disabled_dr_branch_plan,
        )

        dr_cfg = (s0_cfg.get("dr_completion", {}) or {})
        dr_enabled = bool(dr_cfg.get("enabled", True))
        if dr_enabled:
            dr_plan = attach_dr_completion(graph)
        else:
            dr_plan = build_disabled_dr_branch_plan(graph)
        graph_path = s0_dir / "mechanism_graph.json"
        summary_path = s0_dir / "mechanism_summary.json"
        dr_plan_path = s0_dir / "dr_branch_plan.json"
        graph_payload = graph.model_dump(mode="json")
        graph_source_data = graph_payload.get("source_data")
        graph_source_data = graph_source_data if isinstance(graph_source_data, dict) else {}
        low_confidence = self._is_truthy_flag(graph_source_data.get("low_confidence"))
        if pm:
            pm.update_step("s0", completed=80, description="Saving mechanism graph...")

        with open(graph_path, "w", encoding="utf-8") as f:
            json.dump(graph_payload, f, indent=2)
        with open(dr_plan_path, "w", encoding="utf-8") as f:
            json.dump(dr_plan.model_dump(mode="json"), f, indent=2, ensure_ascii=False)

        try:
            import importlib

            visualizer_module = importlib.import_module("rph_core.steps.mechanism_classifier.visualizer")
            MechanismVisualizer = getattr(visualizer_module, "MechanismVisualizer")
            viz = MechanismVisualizer()
            viz.render(graph, s0_dir / "mechanism_graph.png")
        except Exception as exc:
            self.logger.debug(f"[S0] Visualization skipped: {exc}")

        if graph.smiles_atom_mapping:
            mapping_path = s0_dir / "atom_map_smiles.json"
            mapping_payload = {
                "smiles_atom_mapping": graph.smiles_atom_mapping.model_dump(mode="json"),
                "forming_bonds_annotated": (
                    [fb.model_dump(mode="json") for fb in graph.forming_bonds_annotated]
                    if graph.forming_bonds_annotated else None
                ),
                "description": {
                    "precursor_smiles_to_map": "RDKit mol.GetAtomWithIdx(idx).GetIdx() → GetAtomMapNum()",
                    "product_smiles_to_map": "Same for product",
                    "index_base_convention": "SMILES indices are 0-based (RDKit internal); Map# are 1-based",
                },
            }
            with open(mapping_path, "w", encoding="utf-8") as f:
                json.dump(mapping_payload, f, indent=2, ensure_ascii=False)
            self.logger.info(f"[S0] SMILES atom mapping saved: {mapping_path}")

        pathway_edges = graph.get_edges_for_pathway("primary")
        forming_bonds = tuple()
        if pathway_edges:
            all_pairs = []
            for edge in pathway_edges:
                all_pairs.extend(edge.forming_bonds)
            forming_bonds = tuple(
                sorted(
                    (min(int(pair[0]), int(pair[1])), max(int(pair[0]), int(pair[1])))
                    for pair in all_pairs
                )
            )

        summary = {
            "reaction_id": graph.reaction_id,
            "reaction_type": graph.reaction_type,
            "cyclo_mode": getattr(graph.cyclo_mode, "value", graph.cyclo_mode),
            "topology": getattr(graph.topology, "value", graph.topology),
            "low_confidence": low_confidence,
            "forming_bonds": [list(pair) for pair in forming_bonds],
            "atom_mapping_files": {
                "smiles_mapping": str(s0_dir / "atom_map_smiles.json"),
                "xyz_mapping": "S1_ConfGeneration/atom_map_xyz.json (generated after S1)",
            },
            "dr_branch_plan_json": str(dr_plan_path),
            "dr_completion": {
                "status": dr_plan.status,
                "branch_count": len(dr_plan.branches),
                "branch_ids": [branch.branch_id for branch in dr_plan.branches],
            },
        }
        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)

        _write_status(
            "complete",
            extra={
                "mechanism_graph_json": str(graph_path),
                "mechanism_summary_json": str(summary_path),
                "dr_branch_plan_json": str(dr_plan_path),
                "low_confidence": low_confidence,
                "forming_bonds": [list(pair) for pair in forming_bonds],
                "reaction_type": summary.get("reaction_type"),
                "dr_completion": summary.get("dr_completion"),
            },
        )
        if resume_enabled:
            checkpoint_mgr.mark_step_completed(
                "s0",
                output_files={
                    "mechanism_graph_json": str(graph_path),
                    "mechanism_summary_json": str(summary_path),
                    "dr_branch_plan_json": str(dr_plan_path),
                },
                metadata=summary,
            )
        summary["status"] = "complete"

        self.logger.info("[S0] Classification complete:")
        self.logger.info(f"      Reaction Type: {graph.reaction_type}")
        self.logger.info(f"      Cyclo Mode: {getattr(graph.cyclo_mode, 'value', graph.cyclo_mode)}")
        self.logger.info(f"      Topology: {getattr(graph.topology, 'value', graph.topology)}")
        self.logger.info(f"      Forming Bonds: {forming_bonds}")
        self.logger.info(f"      Nodes: {len(graph.nodes)} | Edges: {len(graph.edges)}")

        return summary

    def run_pipeline(
        self,
        product_smiles: str,
        work_dir: Path,
        skip_steps: Optional[list[str]] = None,
        precursor_smiles: Optional[str] = None,
        reaction_id: str = "",
        branch_id: str = "",
        leaving_group_key: Optional[str] = None,
        small_molecular_keys: Optional[list[str]] = None,
        reaction_profile: Optional[str] = None,
        cleaner_data: Optional[Dict[str, Any]] = None,
        include_reference_small_molecules: bool = True,
        quiet_resume: bool = False,
    ) -> PipelineResult:
        self.logger.info("DEBUG: Inside run_pipeline (start)")
        """
        执行 v3.0 DFT pipeline (S0-S3)

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
            include_reference_small_molecules: 是否自动加入 reaction_reference_terms 中的小分子参考项

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

        is_batch = bool(reaction_id)
        if is_batch:
            from rph_core.utils.ui import LoggerProgressManager
            pm = LoggerProgressManager()
        else:
            pm = get_progress_manager()
            
        pm.set_context(batch_mode=is_batch, reaction_id=reaction_id, branch_id=branch_id)

        task_tracker = TaskProgressTracker(V4_TASK_REGISTRY)
        task_tracker.start()

        if quiet_resume:
            resume_enabled_early = bool((self.config.get("run", {}) or {}).get("resume", True))
            checkpoint_mgr_early = CheckpointManager(work_dir)
            cached_count = 0
            if resume_enabled_early:
                if checkpoint_mgr_early.is_step_completed("s1"):
                    cached_count += 1
                if checkpoint_mgr_early.is_step_completed("s2"):
                    cached_count += 1
                s3_dir = resolve_step_dir(work_dir, "s3")
                if checkpoint_mgr_early.is_step3_complete(s3_dir, self.config):
                    cached_count += 1
            if cached_count >= 3 and "s1" not in skip_steps:
                pm = SilentProgressManager()
                pm.set_context(batch_mode=is_batch, reaction_id=reaction_id, branch_id=branch_id)

        ctx = f" | {reaction_id}/{branch_id}" if is_batch else ""
        pm.start(f"RPH v{__version__}{ctx}")

        def _set_task(task_id: str, state: TaskState, detail: str = "", summary: str = "") -> None:
            self.logger.debug(f"DEBUG: Setting task {task_id} to {state}")
            task_tracker.set_state(task_id, state)
            if detail:
                task_tracker.set_detail(task_id, detail)
            if summary:
                task_tracker.set_result_summary(task_id, summary)

        def _render_tasks() -> None:
            try:
                pm.render_tasks(task_tracker)
            except Exception as e:
                self.logger.debug(f"UI rendering suppressed due to error: {e}")

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
                        checkpoint_mgr.initialize_state(product_smiles=product_smiles, config=self.config)
                        state = checkpoint_mgr.load_state()
                    if state is not None:
                        checkpoint_mgr.save_state(state)

            effective_cleaner_data = self._normalize_cleaner_data_for_pipeline(
                cleaner_data,
                product_smiles=product_smiles,
                precursor_smiles=precursor_smiles,
                reaction_profile=reaction_profile,
            )
            if 's0' not in skip_steps:
                _set_task("mechanism", TaskState.RUNNING)
                _render_tasks()
                try:
                    s0_summary = self._run_s0(
                        work_dir=work_dir,
                        product_smiles=product_smiles,
                        cleaner_data=effective_cleaner_data,
                        checkpoint_mgr=checkpoint_mgr,
                        resume_enabled=resume_enabled,
                        pm=pm,
                    )
                    if s0_summary and effective_cleaner_data is not None:
                        if s0_summary.get("reaction_type") and not effective_cleaner_data.get("reaction_type"):
                            effective_cleaner_data["reaction_type"] = s0_summary["reaction_type"]
                    s0_status = str((s0_summary or {}).get("status") or "").lower()
                    s0_reason = str((s0_summary or {}).get("reason") or "").strip()
                    rxn_type = s0_summary.get("reaction_type", "") if s0_summary else ""
                    topology = s0_summary.get("topology", "") if s0_summary else ""
                    summary = f"[{rxn_type} {topology}]" if rxn_type else ""
                    if s0_status == "skipped":
                        _set_task("mechanism", TaskState.SKIPPED, summary=s0_reason)
                    elif s0_status == "degraded":
                        _set_task("mechanism", TaskState.COMPLETED, summary="degraded")
                    else:
                        _set_task("mechanism", TaskState.COMPLETED, summary=summary)
                except Exception as exc:
                    if resume_enabled:
                        checkpoint_mgr.mark_step_failed_partial("s0", phase="classifying", error_message=str(exc))
                    self.logger.warning(f"[S0] mechanism classification failed: {exc}")
                    _set_task("mechanism", TaskState.FAILED, summary=str(exc))
            else:
                _set_task("mechanism", TaskState.SKIPPED)
            _render_tasks()

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
                    cleaner_data=effective_cleaner_data,
                )
                try:
                    if result.product_xyz is not None:
                        product_xyz_file = self._resolve_product_xyz_for_s2(result.product_xyz)
                        resolved_forming_bonds = self._resolve_forming_bonds_for_s2(
                            cleaner_data=effective_cleaner_data,
                            product_xyz_file=product_xyz_file,
                            work_dir=work_dir,
                        )
                        reaction_profiles = self.config.get("reaction_profiles", {}) or {}
                        profile_cfg = {}
                        if isinstance(reaction_profiles, dict) and profile_key:
                            profile_cfg = dict(reaction_profiles.get(str(profile_key), {}) or {})
                        expected_scan_cfg = dict((self.config.get("step2", {}) or {}).get("scan", {}) or {})
                        if profile_cfg:
                            expected_scan_cfg.update(dict(profile_cfg.get("scan", {}) or {}))
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
                    intermediate_xyz_val = checkpoint_mgr.get_step_output('s2', 'intermediate_xyz') or checkpoint_mgr.get_step_output('s2', 'substrate_xyz')
                    if ts_guess and intermediate_xyz_val and Path(ts_guess).exists() and Path(intermediate_xyz_val).exists():
                        result.ts_guess_xyz = Path(ts_guess)
                        result.intermediate_xyz = Path(intermediate_xyz_val)
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
                            f"✅ Resume: Step2 already complete, reuse ts_guess/intermediate: {result.ts_guess_xyz}, {result.intermediate_xyz}"
                        )
                        skip_steps.append('s2')
                        _set_task("retro_scan", TaskState.CACHED)

            s3_dir = resolve_step_dir(work_dir, "s3")
            if resume_enabled and 's3' not in skip_steps:
                input_hashes: Dict[str, str] | None = None
                if result.ts_guess_xyz and result.intermediate_xyz and result.product_xyz:
                    input_hashes = {
                        'ts_guess': checkpoint_mgr.compute_file_hash(result.ts_guess_xyz) or '',
                        'intermediate': checkpoint_mgr.compute_file_hash(result.intermediate_xyz) or '',
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
                            _set_task("ts_opt", TaskState.CACHED)
                            _set_task("irc_verify", TaskState.CACHED)
                            _set_task("reactant_opt", TaskState.CACHED)
                            _set_task("sp_matrix", TaskState.CACHED)
                            _set_task("thermochemistry", TaskState.CACHED)
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

            # === Anchor Phase: product / precursor / small molecules ===
            if 's1' not in skip_steps:
                try:
                    if resume_enabled:
                        checkpoint_mgr.mark_step_in_progress("s1", phase="anchor")

                    s1_work_dir = work_dir / "S1_ConfGeneration"
                    self.s1_engine.base_work_dir = s1_work_dir
                    self.s1_engine.base_work_dir.mkdir(parents=True, exist_ok=True)

                    resolved_sm_keys = self._resolve_small_molecular_keys(
                        reaction_profile=reaction_profile,
                        dataset_keys=small_molecular_keys,
                        include_reference_terms=include_reference_small_molecules,
                    )

                    molecules: dict[str, str] = {"product": product_smiles}
                    if precursor_smiles:
                        molecules["precursor"] = precursor_smiles
                    from rph_core.utils.small_molecule_catalog import SmallMoleculeCatalog
                    sm_catalog = SmallMoleculeCatalog(self.config)
                    if leaving_group_key:
                        leaving_group = sm_catalog.get(leaving_group_key)
                        if leaving_group and leaving_group.smiles:
                            molecules["leaving_group"] = leaving_group.smiles
                        else:
                            self.logger.warning(
                                f"Small molecular key '{leaving_group_key}' not in catalog, skipping S1 processing"
                            )
                    for key in resolved_sm_keys:
                        if key == leaving_group_key:
                            continue
                        mol = sm_catalog.get(key)
                        if mol and mol.smiles:
                            molecules[key] = mol.smiles

                    _set_task("product_anchor", TaskState.RUNNING)
                    if precursor_smiles:
                        _set_task("precursor_anchor", TaskState.RUNNING)
                    if resolved_sm_keys or leaving_group_key:
                        _set_task("smallmol_anchor", TaskState.RUNNING)
                    _render_tasks()

                    anchor_result = self.s1_engine.run(
                        molecules=molecules
                    )

                    _set_task("product_anchor", TaskState.COMPLETED,
                              summary=f"E_sp = {result.e_product_l2:.6f} Ha" if result.e_product_l2 else "")
                    if precursor_smiles and "precursor" in anchor_result.anchored_molecules:
                        _set_task("precursor_anchor", TaskState.COMPLETED)
                    if resolved_sm_keys or leaving_group_key:
                        has_small = any(k in anchor_result.anchored_molecules for k in (resolved_sm_keys or []) + ([leaving_group_key] if leaving_group_key else []))
                        if has_small:
                            _set_task("smallmol_anchor", TaskState.COMPLETED)
                    _render_tasks()

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

                    thermo_candidates = [
                        s1_work_dir / "product" / "finalDFT" / "conformer_thermo.csv",
                        s1_work_dir / "product" / "dft" / "conformer_thermo.csv",
                    ]
                    for product_thermo_file in thermo_candidates:
                        if product_thermo_file.exists():
                            result.product_thermo = product_thermo_file
                            break

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
                    if resume_enabled:
                        checkpoint_mgr.mark_step_failed_partial("s1", phase="anchor", error_message=str(e))
                    result.error_step = "Step1_ProductAnchor_v3"
                    result.error_message = str(e)
                    self.logger.error(f"Step 1 v3.0 失败: {e}", exc_info=True)
                    _notify(False, result.error_step, result.error_message)
                    return result
            else:
                self.logger.warning("⚠️  跳过 Step 1")
                _set_task("product_anchor", TaskState.SKIPPED)
                _set_task("precursor_anchor", TaskState.SKIPPED)
                _set_task("smallmol_anchor", TaskState.SKIPPED)
                candidate = resolve_step_dir(work_dir, "s1")
                if candidate.exists():
                    result.product_xyz = candidate
                    self.logger.info(f"    ✓ 复用 S1 输出目录: {candidate}")

            if result.product_xyz is not None:
                try:
                    self._build_and_save_xyz_mapping(
                        work_dir=work_dir,
                        cleaner_data=effective_cleaner_data,
                        product_xyz=result.product_xyz,
                    )
                except Exception as exc:
                    self.logger.warning(f"[S1→S2] Failed to build XYZ atom mapping: {exc}")

            if 's2' not in skip_steps and result.product_xyz:
                try:
                    _set_task("retro_scan", TaskState.RUNNING)
                    _render_tasks()
                    if resume_enabled:
                        checkpoint_mgr.mark_step_in_progress("s2", phase="scan")

                    step2_artifacts = run_step2(
                        hunter=self,
                        product_xyz=result.product_xyz,
                        work_dir=work_dir,
                        reaction_profile=reaction_profile,
                        cleaner_data=effective_cleaner_data,
                    )

                    scan_summary = ""
                    if step2_artifacts.scan_profile_json and step2_artifacts.scan_profile_json.exists():
                        try:
                            data = json.loads(step2_artifacts.scan_profile_json.read_text(encoding="utf-8"))
                            scan_steps = data.get("scan_steps", "?")
                            scan_summary = f"{scan_steps} steps"
                        except Exception:
                            pass
                    _set_task("retro_scan", TaskState.COMPLETED, summary=scan_summary)
                    result.ts_guess_xyz = step2_artifacts.ts_guess_xyz
                    result.intermediate_xyz = step2_artifacts.intermediate_xyz
                    result.forming_bonds = step2_artifacts.forming_bonds
                    current_step2_signature = step2_artifacts.step2_signature
                    self.logger.info(f"    ✓ TS initial guess: {step2_artifacts.ts_guess_xyz}")
                    self.logger.info(f"    ✓ Intermediate: {step2_artifacts.intermediate_xyz}")
                    self.logger.info(f"    ✓ S2 status/confidence: {step2_artifacts.status}/{step2_artifacts.ts_guess_confidence}")

                    if resume_enabled:
                        checkpoint_mgr.mark_step_completed(
                            "s2",
                        output_files={
                            "ts_guess_xyz": str(result.ts_guess_xyz),
                            "intermediate_xyz": str(result.intermediate_xyz),
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
                    if resume_enabled:
                        checkpoint_mgr.mark_step_failed_partial("s2", phase="scan", error_message=str(e))
                    result.error_step = "Step2_TSGuessBuilder"
                    result.error_message = str(e)
                    self.logger.error(f"Step 2 失败: {e}", exc_info=True)
                    _notify(False, result.error_step, result.error_message)
                    return result
            else:
                self.logger.warning("⚠️  跳过 Step 2")
                _set_task("retro_scan", TaskState.SKIPPED)

            # === TS Phase: TS opt / IRC / reactant opt / SP matrix / thermochemistry ===
            if 's3' not in skip_steps and result.ts_guess_xyz:
                try:
                    pm.section_header("Step 3", "Transition Analyzer", "TS Optimization & Verification")
                    self.logger.info(">>> Step 3: 反应中心全分析 (TS优化 + Intermediate/Fragments SP)...")
                    if resume_enabled:
                        checkpoint_mgr.mark_step_in_progress("s3", phase="optimize")

                    # 传递S1的checkpoint以复用轨道
                    old_checkpoint = result.product_checkpoint
                    if old_checkpoint:
                        self.logger.info(f"  复用S1 checkpoint: {old_checkpoint.name}")

                    if result.product_xyz and result.product_xyz.is_dir():
                        product_dir = result.product_xyz
                        candidates = [
                            product_dir / "product_min.xyz",
                        ]
                        product_file = next((p for p in candidates if p.exists()), None)
                        if product_file is None:
                            raise RuntimeError(
                                f"无法在 {product_dir} 中找到产物结构文件用于 S3。"
                                f"已尝试: {[str(p) for p in candidates]}"
                            )
                        result.product_xyz = product_file
                        self.logger.info(f"  ✓ 使用产物文件: {result.product_xyz}")

                    if result.intermediate_xyz is None or result.product_xyz is None:
                        raise RuntimeError("Step3 输入缺失: intermediate 或 product 为 None")

                    step3_intermediate_xyz = result.intermediate_xyz
                    step3_product_xyz = result.product_xyz

                    _set_task("ts_opt", TaskState.RUNNING)
                    _set_task("reactant_opt", TaskState.RUNNING)
                    _render_tasks()
                    step3_artifacts = run_step3(
                        hunter=self,
                        ts_guess_xyz=result.ts_guess_xyz,
                        intermediate_xyz=step3_intermediate_xyz,
                        product_xyz=step3_product_xyz,
                        work_dir=work_dir,
                        e_product_l2=result.e_product_l2,
                        product_thermo=result.product_thermo,
                        forming_bonds=result.forming_bonds,
                        old_checkpoint=old_checkpoint,
                    )
                    ts_summary = ""
                    if step3_artifacts.sp_report:
                        dg = step3_artifacts.sp_report.get_activation_energy()
                        if dg is not None:
                            ts_summary = f"ΔG‡ = {dg:.2f} kcal/mol"
                    _set_task("ts_opt", TaskState.COMPLETED, summary=ts_summary)
                    _set_task("reactant_opt", TaskState.COMPLETED)
                    _set_task("irc_verify", TaskState.COMPLETED)
                    _set_task("sp_matrix", TaskState.COMPLETED)
                    _set_task("thermochemistry", TaskState.COMPLETED)
                    _render_tasks()

                    result.ts_final_xyz = step3_artifacts.ts_final_xyz
                    result.sp_matrix_report = step3_artifacts.sp_report
                    
                    result.ts_fchk = step3_artifacts.ts_fchk
                    result.ts_log = step3_artifacts.ts_log
                    result.ts_qm_output = step3_artifacts.ts_qm_output
                    result.intermediate_fchk = step3_artifacts.intermediate_fchk
                    result.intermediate_log = step3_artifacts.intermediate_log
                    result.intermediate_qm_output = step3_artifacts.intermediate_qm_output
                    result.s3_intermediate_xyz = step3_artifacts.intermediate_xyz
                    result.s3_intermediate_l2_energy = step3_artifacts.intermediate_l2_energy

                    self.logger.info("    ✓ S3 完成")

                    sp_meta_path = s3_dir / "sp_matrix_metadata.json"
                    if sp_meta_path.exists() and current_step2_signature is not None:
                        try:
                            with open(sp_meta_path, "r", encoding="utf-8") as f:
                                sp_meta_data = json.load(f)
                            sp_meta_data["upstream_step2_signature"] = current_step2_signature
                            with open(sp_meta_path, "w", encoding="utf-8") as f:
                                json.dump(sp_meta_data, f, indent=2)
                            prov_path = s3_dir / "step3_provenance.json"
                            if prov_path.exists():
                                with open(prov_path, "r", encoding="utf-8") as f:
                                    prov_data = json.load(f)
                                prov_data["upstream_step2_signature"] = current_step2_signature
                                with open(prov_path, "w", encoding="utf-8") as f:
                                    json.dump(prov_data, f, indent=2)
                        except Exception as exc:
                            self.logger.warning(f"Failed to annotate S3 metadata with S2 signature: {exc}")

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
                                    "e_intermediate": getattr(result.sp_matrix_report, "e_reactant", None),
                                    "e_reactant": getattr(result.sp_matrix_report, "e_reactant", None),  # legacy alias for e_intermediate
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
                    if resume_enabled:
                        checkpoint_mgr.mark_step_failed_partial("s3", phase="optimize", error_message=str(e))
                    result.error_step = "Step3_TransitionAnalyzer"
                    result.error_message = str(e)
                    self.logger.error(f"Step 3 失败: {e}", exc_info=True)
                    _notify(False, result.error_step, result.error_message)
                    return result
            else:
                self.logger.warning("⚠️  跳过 Step 3")
                _set_task("ts_opt", TaskState.SKIPPED)
                _set_task("irc_verify", TaskState.SKIPPED)
                _set_task("reactant_opt", TaskState.SKIPPED)
                _set_task("sp_matrix", TaskState.SKIPPED)
                _set_task("thermochemistry", TaskState.SKIPPED)

            # === Features Phase: geometry / electronic / thermo / NBO ===
            # Rehydrate fields from disk if lost during checkpoint resume or step skip
            if not result.ts_final_xyz:
                candidate = s3_dir / "ts_final.xyz"
                if candidate.exists():
                    result.ts_final_xyz = candidate
                    self.logger.info(f"    ↻ Fallback ts_final_xyz from disk: {candidate}")
            if not result.intermediate_xyz:
                s2_dir = resolve_step_dir(work_dir, "s2")
                candidate = s2_dir / "intermediate.xyz"
                if not candidate.exists():
                    candidate = s2_dir / "reactant_complex.xyz"
                if candidate.exists():
                    result.intermediate_xyz = candidate
                    self.logger.info(f"    ↻ Fallback intermediate_xyz from disk: {candidate}")
            if result.product_xyz and result.product_xyz.is_dir():
                product_file = result.product_xyz / "product_min.xyz"
                if product_file.exists():
                    result.product_xyz = product_file
                    self.logger.info(f"    ↻ Resolved product_xyz from directory to: {product_file}")
            if not result.ts_guess_xyz:
                s2_dir = resolve_step_dir(work_dir, "s2")
                candidate = s2_dir / "ts_guess.xyz"
                if candidate.exists():
                    result.ts_guess_xyz = candidate
                    self.logger.info(f"    ↻ Fallback ts_guess_xyz from disk: {candidate}")
            # Rehydrate fchk/log/qm_output paths from checkpoint or disk when S3 was skipped
            if not result.ts_fchk:
                _v = checkpoint_mgr.get_step_output('s3', 'ts_fchk')
                result.ts_fchk = Path(_v) if _v else None
            if not result.ts_log:
                _v = checkpoint_mgr.get_step_output('s3', 'ts_log')
                result.ts_log = Path(_v) if _v else None
            if not result.ts_qm_output:
                _v = checkpoint_mgr.get_step_output('s3', 'ts_qm_output')
                result.ts_qm_output = Path(_v) if _v else None
            if not result.intermediate_fchk:
                _v = checkpoint_mgr.get_step_output('s3', 'intermediate_fchk')
                result.intermediate_fchk = Path(_v) if _v else None
            if not result.intermediate_log:
                _v = checkpoint_mgr.get_step_output('s3', 'intermediate_log')
                result.intermediate_log = Path(_v) if _v else None
            if not result.intermediate_qm_output:
                _v = checkpoint_mgr.get_step_output('s3', 'intermediate_qm_output')
                result.intermediate_qm_output = Path(_v) if _v else None
            # Also try to rehydrate sp_matrix_report if missing
            if result.sp_matrix_report is None:
                sp_meta_path = s3_dir / "sp_matrix_metadata.json"
                if sp_meta_path.exists():
                    try:
                        with open(sp_meta_path, 'r') as _f:
                            sp_meta = json.load(_f)
                        from rph_core.steps.step3_opt.ts_optimizer import SPMatrixReport
                        result.sp_matrix_report = SPMatrixReport(
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
                        self.logger.info("    ↻ Rehydrated sp_matrix_report from sp_matrix_metadata.json")
                    except Exception as _exc:
                        self.logger.debug(f"    ↻ Failed to rehydrate sp_matrix_report: {_exc}")
            if 's4' not in skip_steps and result.ts_final_xyz:
                self.logger.info(
                    ">>> Step 4: Feature extraction has been moved to RPH_Postprocess. "
                    "DFT pipeline stops at S3.\n"
                    "    To extract features separately:\n"
                    "      rph-features extract --rph-run %s --output <features_dir>",
                    work_dir,
                )
                _set_task("geom_features", TaskState.SKIPPED)
                _set_task("elec_features", TaskState.SKIPPED)
                _set_task("thermo_features", TaskState.SKIPPED)
                _set_task("nbo_features", TaskState.SKIPPED)
                _render_tasks()
            elif result.ts_final_xyz is None:
                self.logger.info("Step 4 skipped (no TS final geometry)")
            else:
                self.logger.info("Step 4 skipped by --skip-steps")

            # 成功完成
            result.success = True
            _render_tasks()
            elapsed = task_tracker.elapsed_sec
            summary = task_tracker.get_summary()
            self.logger.info(f"✅ {summary}")
            self.logger.info(f"   数据已保存至: {work_dir}")
            _notify(True)

            return result
            
        except KeyboardInterrupt:
            result.error_step = "KeyboardInterrupt"
            result.error_message = "Pipeline was aborted by user."
            self.logger.error("\n[!] 进程被用户中断 (KeyboardInterrupt)")
            _notify(False, result.error_step, result.error_message)
            return result
            
        except Exception as e:
            result.error_step = "UnexpectedError"
            result.error_message = str(e)
            self.logger.error(f"\n[!] 管道执行遭遇未捕获的致命错误: {e}", exc_info=True)
            _notify(False, result.error_step, result.error_message)
            return result
            
        finally:
            try:
                pm.stop()
            except Exception as e:
                self.logger.error(f"UI资源释放失败: {e}")

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

    def _resolve_s1_artifacts(self, work_dir: Path) -> Dict[str, Any]:
        """Resolve S1 artifact paths for Step4 Step1 activation features.

        Searches for and optionally derives S1 artifacts needed by S4:
        - shermo_summary.json (from .sum files if missing)
        - conformer_energies.json
        - precursor xyz

        Args:
            work_dir: Pipeline working directory

        Returns:
            Dictionary with resolved artifact paths
        """
        from rph_core.utils.shermo_runner import find_shermo_sum_files

        artifacts: Dict[str, Any] = {
            "s1_dir": None,
            "s1_shermo_summary_file": None,
            "s1_conformer_energies_file": None,
            "s1_crest_ensemble_energies_file": None,
            "s1_ensemble_xyz_file": None,
            "s1_conformer_dir": None,
            "s1_conformer_thermo_csv": None,
            "s1_precursor_xyz": None,
            "s1_precursor_conformer_dir": None,
            "s1_precursor_conformer_energies_file": None,
            "s1_precursor_conformer_thermo_csv": None,
            "s1_precursor_ensemble_xyz_file": None,
            "s1_precursor_crest_ensemble_energies_file": None,
            "s1_atom_map_xyz_file": None,
        }

        # Find S1 directory
        s1_candidates = [resolve_step_dir(work_dir, "s1")]
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
                from rph_core.utils.thermo import derive_shermo_summary
                sum_to_use = sum_files.get("precursor") or sum_files.get("ylide")
                if sum_to_use is not None:
                    derive_shermo_summary(sum_to_use, shermo_summary, molecule_type="precursor")
                    artifacts["s1_shermo_summary_file"] = shermo_summary
                    self.logger.info(f"Derived shermo_summary.json from .sum files via unified API")

        # Resolve small molecule thermo based on reaction_reference_terms stoichiometry
        # instead of hardcoded molecule list. This ensures S4 thermodynamics is consistent
        # with S0 mechanism graph semantics across all reaction types.
        from rph_core.utils.small_molecule_catalog import SmallMoleculeCatalog
        catalog = SmallMoleculeCatalog(self.config)
        ref_terms = self.config.get("reaction_reference_terms", {}) or {}
        profile_key = self.config.get("run", {}).get("reaction_type", "[4+3]_default")
        profile_cfg = (ref_terms.get(profile_key, {}) or {}) if isinstance(ref_terms, dict) else {}
        pti = (profile_cfg.get("precursor_to_intermediate", {}) or {}) if isinstance(profile_cfg, dict) else {}
        reserved = {"precursor", "intermediate"}

        s1_small_molecule_gibbs: Dict[str, Optional[float]] = {}

        for side in ("reactants", "products"):
            side_dict = (pti.get(side, {}) or {}) if isinstance(pti, dict) else {}
            for mol_key in side_dict:
                if mol_key in reserved:
                    continue
                mol = catalog.get(mol_key)
                if mol is None:
                    self.logger.warning(f"Small molecule '{mol_key}' not in catalog, skipping thermo resolution")
                    s1_small_molecule_gibbs[mol_key] = None
                    continue
                thermo_path = self._resolve_small_molecule_thermo(mol_key, mol.smiles, s1_dir)
                if thermo_path:
                    artifacts[f"s1_{mol_key.lower()}_thermo_file"] = thermo_path
                    try:
                        import json as _json
                        with open(thermo_path, 'r') as _f:
                            _data = _json.load(_f)
                        g_val = _data.get('g_kcal') or _data.get('G') or _data.get('g')
                        s1_small_molecule_gibbs[mol_key] = float(g_val) if g_val is not None else None
                    except Exception:
                        s1_small_molecule_gibbs[mol_key] = None
                else:
                    self.logger.warning(f"No thermo data resolved for small molecule '{mol_key}'")
                    s1_small_molecule_gibbs[mol_key] = None

        artifacts["s1_small_molecule_gibbs"] = s1_small_molecule_gibbs

        # Look for conformer_energies.json
        conformer_energies = s1_dir / "conformer_energies.json"
        if not conformer_energies.exists():
            # Check in molecule subdirectories
            for mol_dir in s1_dir.iterdir():
                if mol_dir.is_dir() and not mol_dir.name.startswith('.'):
                    candidate_new = mol_dir / "finalDFT" / "conformer_energies.json"
                    if candidate_new.exists():
                        conformer_energies = candidate_new
                        break

                    candidate_legacy = mol_dir / "dft" / "conformer_energies.json"
                    if candidate_legacy.exists():
                        conformer_energies = candidate_legacy
                        break
        if conformer_energies.exists():
            artifacts["s1_conformer_energies_file"] = conformer_energies
            conformer_parent = conformer_energies.parent
            artifacts["s1_conformer_dir"] = conformer_parent
            thermo_csv = conformer_parent / "conformer_thermo.csv"
            if thermo_csv.exists():
                artifacts["s1_conformer_thermo_csv"] = thermo_csv

        # Resolve crest ensemble energies (GFN2 full conformer ensemble)
        for mol_dir in s1_dir.iterdir():
            if not mol_dir.is_dir() or mol_dir.name.startswith('.'):
                continue
            for crest_sub in ("crest", "xtb/stage2_gfn2"):
                crest_dir = mol_dir / crest_sub
                cee = crest_dir / "conformer_ensemble_energies.json"
                if cee.is_file():
                    artifacts["s1_crest_ensemble_energies_file"] = cee
                    break
            if artifacts["s1_crest_ensemble_energies_file"] is not None:
                break
        if artifacts["s1_crest_ensemble_energies_file"] is None:
            for crest_sub in ("crest",):
                crest_dir = s1_dir / crest_sub
                cee = crest_dir / "conformer_ensemble_energies.json"
                if cee.is_file():
                    artifacts["s1_crest_ensemble_energies_file"] = cee
                    break

        # Resolve CREST GFN2 full ensemble XYZ for geometry-based extractors
        for mol_dir in s1_dir.iterdir():
            if not mol_dir.is_dir() or mol_dir.name.startswith('.'):
                continue
            for search_sub in (
                "xtb/stage2_gfn2",
                "crest",
            ):
                search_dir = mol_dir / search_sub
                if not search_dir.exists():
                    continue
                for xyz_name in (
                    "crest_ensemble.xyz",
                    "crest_conformers.xyz",
                    "ensemble.xyz",
                ):
                    candidate = search_dir / xyz_name
                    if candidate.is_file() and candidate.stat().st_size > 500:
                        artifacts["s1_ensemble_xyz_file"] = candidate
                        break
                if artifacts["s1_ensemble_xyz_file"] is not None:
                    break
            if artifacts["s1_ensemble_xyz_file"] is not None:
                break

        if artifacts["s1_ensemble_xyz_file"] is None:
            for search_sub in ("crest",):
                search_dir = s1_dir / search_sub
                if search_dir.exists():
                    for xyz_name in ("crest_conformers.xyz", "ensemble.xyz"):
                        candidate = search_dir / xyz_name
                        if candidate.is_file() and candidate.stat().st_size > 500:
                            artifacts["s1_ensemble_xyz_file"] = candidate
                            break

        # Look for precursor xyz
        precursor_xyz = s1_dir / "precursor" / "precursor_min.xyz"
        if precursor_xyz.exists():
            artifacts["s1_precursor_xyz"] = precursor_xyz

        # Look for precursor conformer ensemble
        precursor_conf_dir = s1_dir / "precursor" / "finalDFT"
        if precursor_conf_dir.exists():
            artifacts["s1_precursor_conformer_dir"] = precursor_conf_dir
            prec_energies = precursor_conf_dir / "conformer_energies.json"
            if prec_energies.exists():
                artifacts["s1_precursor_conformer_energies_file"] = prec_energies
            prec_thermo = precursor_conf_dir / "conformer_thermo.csv"
            if prec_thermo.exists():
                artifacts["s1_precursor_conformer_thermo_csv"] = prec_thermo

        # Resolve precursor CREST GFN2 full ensemble data
        for precursor_sub in ("precursor",):
            prec_s1 = s1_dir / precursor_sub
            if not prec_s1.is_dir():
                continue
            for crest_sub in ("crest",):
                crest_prec_dir = prec_s1 / crest_sub
                pcce = crest_prec_dir / "conformer_ensemble_energies.json"
                if pcce.is_file():
                    artifacts["s1_precursor_crest_ensemble_energies_file"] = pcce
                    break
            for search_sub in (
                "xtb/stage2_gfn2",
                "crest",
            ):
                search_dir = prec_s1 / search_sub
                if not search_dir.exists():
                    continue
                for xyz_name in (
                    "crest_ensemble.xyz",
                    "crest_conformers.xyz",
                    "ensemble.xyz",
                ):
                    candidate = search_dir / xyz_name
                    if candidate.is_file() and candidate.stat().st_size > 500:
                        artifacts["s1_precursor_ensemble_xyz_file"] = candidate
                        break
                if artifacts["s1_precursor_ensemble_xyz_file"] is not None:
                    break
            if artifacts["s1_precursor_ensemble_xyz_file"] is not None:
                break

        # Look for atom_map_xyz.json (containing product↔precursor atom mapping)
        atom_map = s1_dir / "atom_map_xyz.json"
        if atom_map.exists():
            artifacts["s1_atom_map_xyz_file"] = atom_map

        return artifacts

    def _resolve_small_molecule_thermo(self, mol_key: str, smiles: str, s1_dir: Path) -> Optional[Path]:
        from rph_core.utils.small_molecule_cache import SmallMoleculeCache
        from rph_core.utils.thermo import derive_thermo_json_from_sum
        from rph_core.utils.molecule_utils import get_molecule_key

        # Check local per-reaction small_molecules dir
        local_thermo = s1_dir / "small_molecules" / mol_key / "thermo.json"
        if local_thermo.exists():
            return local_thermo

        # Try deriving from .sum files in the molecule's own S1 directory
        mol_s1_dir = s1_dir / mol_key
        if mol_s1_dir.is_dir():
            for dft_dir_name in ("finalDFT", "dft"):
                dft_dir = mol_s1_dir / dft_dir_name
                if dft_dir.is_dir():
                    sum_candidates = list(dft_dir.glob("*_Shermo.sum")) + list(dft_dir.glob("*Shermo*.sum"))
                    if sum_candidates:
                        mol_sum = sum_candidates[0]
                        local_thermo.parent.mkdir(parents=True, exist_ok=True)
                        derive_thermo_json_from_sum(mol_sum, local_thermo)
                        self.logger.info(f"Derived {mol_key} thermo.json from {mol_sum}")
                        return local_thermo

        # Check global small molecule cache
        global_cache_dir = self.config.get("global", {}).get("small_molecule_cache_dir")
        if not global_cache_dir:
            return None

        cache = SmallMoleculeCache(Path(global_cache_dir))
        mol_cache_key = get_molecule_key(smiles)
        if not mol_cache_key:
            return None

        cached_thermo = cache.find_thermo(smiles)
        if cached_thermo:
            self.logger.info(f"Found {mol_key} thermo in global cache: {cached_thermo}")
            return cached_thermo

        # Try deriving from .sum in cache finalDFT (primary) or dft (legacy) dir
        cache_mol_dir = cache.cache_root / mol_cache_key
        if cache_mol_dir.is_dir():
            for dft_dir_name in ("finalDFT", "dft"):
                cache_dft = cache_mol_dir / dft_dir_name
                if cache_dft.is_dir():
                    sum_candidates = list(cache_dft.glob("*_Shermo.sum")) + list(cache_dft.glob("*Shermo*.sum"))
                    if sum_candidates:
                        mol_sum = sum_candidates[0]
                        dest = cache_mol_dir / "thermo.json"
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        derive_thermo_json_from_sum(mol_sum, dest)
                        self.logger.info(f"Derived {mol_key} thermo from global cache .sum: {mol_sum}")
                        return dest

        return None


# =============================================================================
# 命令行接口
# =============================================================================

def _resolve_run_config(config: dict[str, Any], args) -> dict[str, Any]:
    run_cfg = copy.deepcopy(config.get("run", {}) or {})
    global_cfg = config.get("global", {}) or {}
    run_cfg.setdefault("source", "dataset")
    run_cfg.setdefault("output_root", global_cfg.get("work_dir_base", "./rph_output"))
    run_cfg.setdefault("workdir_naming", "rx_{rx_id}")
    run_cfg.setdefault("resume", True)
    run_cfg.setdefault("dry_run", False)
    run_cfg.setdefault("max_tasks", 0)
    run_cfg.setdefault("filter_ids", [])
    run_cfg.setdefault("filter_rx_id", None)

    if args.output:
        run_cfg["output_root"] = args.output

    if getattr(args, "reaction_type", None):
        reaction_profile = str(args.reaction_type).strip()
        if reaction_profile:
            run_cfg["reaction_profile"] = reaction_profile
            run_cfg["reaction_type"] = reaction_profile

    skip_steps_value = getattr(args, "skip_steps", None)
    if skip_steps_value:
        parsed = [s.strip().lower() for s in skip_steps_value.split(",") if s.strip()]
        run_cfg["skip_steps"] = parsed

    # S4 feature extraction is off by default in V3.0
    # --features enables it; --skip-features is backward-compat no-op
    compute_features = run_cfg.get('compute_features', config.get('run', {}).get('compute_features', False))
    if getattr(args, 'features', False):
        compute_features = True
    current_skip = run_cfg.get("skip_steps", [])
    if not compute_features and 's4' not in current_skip:
        current_skip.append('s4')
    run_cfg["skip_steps"] = current_skip

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


def _load_dr_branch_plan(reaction_root: Path) -> dict[str, Any] | None:
    plan_path = reaction_root / "S0_Mechanism" / "dr_branch_plan.json"
    if not plan_path.exists():
        return None
    try:
        with open(plan_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def _write_branch_manifest(
    branch_root: Path,
    branch_id: str,
    pathway_id: str,
    product_smiles: str,
    generation_policy: str,
    flipped_map_numbers: list[int],
    fixed_stereocenters: list[int],
    notes: list[str],
    dr_plan_path: str,
    parent_product_smiles: str,
) -> None:
    manifest = {
        "version": "rph-branch-v1",
        "branch_id": branch_id,
        "pathway_id": pathway_id,
        "source_dr_branch_plan": dr_plan_path,
        "generation_policy": generation_policy,
        "product_smiles": product_smiles,
        "flipped_map_numbers": flipped_map_numbers,
        "fixed_stereocenters": fixed_stereocenters,
        "notes": notes,
        "parent_product_smiles": parent_product_smiles,
        "status": "PENDING",
    }
    branch_root.mkdir(parents=True, exist_ok=True)
    with open(branch_root / "branch_manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)


def _run_tasks_v3(hunter: ReactionProfileHunter, run_cfg: dict[str, Any]) -> list[PipelineResult]:
    import os
    import importlib
    from rph_core.utils.v3_progress import V3RunProgress
    from rph_core.utils.v3_stage_display import LoggerV3StageDisplay, SilentV3StageDisplay

    hunter.logger.info("DEBUG: Importing v3_scheduler...")
    V3Scheduler = getattr(importlib.import_module("rph_core.scheduling.v3_scheduler"), "V3Scheduler")

    progress = V3RunProgress()
    mode = os.environ.get("RPH_UI_MODE", "log").lower().strip()
    if mode == "none":
        display = SilentV3StageDisplay()
    else:
        display = LoggerV3StageDisplay()

    hunter.logger.info("DEBUG: Creating V3Scheduler and calling run()...")
    return V3Scheduler(
        hunter=hunter,
        run_cfg=run_cfg,
        display=display,
        progress=progress,
    ).run()


def _run_tasks(hunter: ReactionProfileHunter, run_cfg: dict[str, Any]) -> list[PipelineResult]:
    run_cfg_for_build = copy.deepcopy(run_cfg)
    if str(run_cfg_for_build.get("source", "dataset")) == "dataset":
        dataset_cfg_value = run_cfg_for_build.get("dataset") or {}
        dataset_cfg = dict(dataset_cfg_value) if isinstance(dataset_cfg_value, dict) else {}
        dataset_cfg["reaction_profiles"] = hunter.config.get("reaction_profiles", {}) or {}
        run_cfg_for_build["dataset"] = dataset_cfg

    def _theory_signature_string(config: dict[str, Any]) -> str:
        theory_opt = (config.get("theory", {}) or {}).get("optimization", {}) or {}
        theory_sp = (config.get("theory", {}) or {}).get("single_point", {}) or {}
        opt_method = str(theory_opt.get("method") or "").strip()
        opt_basis = str(theory_opt.get("basis") or "").strip()
        sp_method = str(theory_sp.get("method") or "").strip()
        sp_basis = str(theory_sp.get("basis") or "").strip()
        solvent_cfg = (config.get("solvent", {}) or {})
        solvent_name = str(solvent_cfg.get("name") or "").strip().lower()
        if not solvent_name:
            solvent_name = str(theory_opt.get("solvent") or theory_sp.get("solvent") or "").strip().lower()
        solvent_token = solvent_name or "dcm"
        tokens = [t for t in [opt_method, opt_basis, sp_method, sp_basis] if t]
        base = "_".join(tokens) if tokens else ""
        if base:
            return f"{base}_SMD_{solvent_token.upper()}"
        return f"SMD_{solvent_token.upper()}"

    tasks = build_tasks_from_run_config(run_cfg_for_build, theory_signature=_theory_signature_string(hunter.config))
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

    from collections import defaultdict
    import json
    from rph_core.utils.path_manager import get_reaction_root, get_branch_root
    from rph_core.steps.condition_thermo import ConditionThermoCalculator
    from rph_core.steps.condition_feature_merger import ConditionFeatureMerger

    def _to_float(value: Any) -> float | None:
        if value is None:
            return None
        try:
            return float(str(value).strip())
        except (TypeError, ValueError):
            return None



    grouped: dict[str, list[Any]] = defaultdict(list)
    for task in tasks:
        reaction_id = getattr(task, "reaction_id", None) or sanitize_rx_id(task.rx_id)
        grouped[str(reaction_id)].append(task)

    results: list[PipelineResult] = []
    skip_steps = run_cfg.get("skip_steps", ["s4"])
    thermo_cfg = hunter.config.get("thermo", {}) or {}
    default_temperature_k = _to_float(thermo_cfg.get("temperature_k")) or 298.15

    for reaction_id, condition_tasks in grouped.items():
        representative = condition_tasks[0]
        reaction_root = get_reaction_root(output_root, reaction_id)
        reaction_root.mkdir(parents=True, exist_ok=True)

        row_ids = [str(getattr(t, "row_id", t.rx_id)) for t in condition_tasks]
        condition_ids = [str(getattr(t, "condition_id", f"COND_{t.rx_id}")) for t in condition_tasks]
        reaction_manifest = {
            "version": "rph_v3.0.0",
            "reaction_id": reaction_id,
            "reaction_cache_key": str(getattr(representative, 'reaction_cache_key', '')),
            "reactant_smiles_canon": (representative.meta or {}).get("reactant_smiles_canon", ""),
            "product_smiles_canon": (representative.meta or {}).get("product_smiles_canon", ""),
            "reaction_type": (representative.meta or {}).get("reaction_type") or "",
            "cyclo_mode": str(((representative.meta or {}).get("cleaner_data") or {}).get("cyclo_mode") or "concerted"),
            "topology": str(((representative.meta or {}).get("cleaner_data") or {}).get("topology") or ""),
            "row_ids": row_ids,
            "condition_ids": condition_ids,
            "theory_signature": _theory_signature_string(hunter.config),
            "status": "PENDING",
        }
        manifest_path = reaction_root / "reaction_manifest.json"
        if manifest_path.exists():
            try:
                existing = read_json_dict(manifest_path)
                if isinstance(existing, dict):
                    for k, v in reaction_manifest.items():
                        if k not in existing or existing.get(k) in (None, ""):
                            existing[k] = v
                    reaction_manifest = existing
            except Exception:
                pass
        write_json(manifest_path, reaction_manifest)

        reaction_features_dir = reaction_root / "reaction_features"

        if not run_cfg.get("dry_run", False):
            if run_cfg.get("resume", True):
                checkpoint_mgr = CheckpointManager(reaction_root)
                s3_dir = reaction_root / "S3_TS"
                if s3_dir.exists() and checkpoint_mgr.is_step3_complete(s3_dir, hunter.config):
                    hunter.logger.info(
                        "%s: Step3 complete. S4 feature extraction is now external — "
                        "use: rph-features extract --rph-run %s",
                        reaction_id, reaction_root,
                    )
                else:
                    skip_steps_reaction = [s for s in list(skip_steps) if str(s).strip()]
                    hunter.logger.info(f"DEBUG: Calling hunter.run_pipeline for {reaction_id}...")
                    result = hunter.run_pipeline(
                        product_smiles=representative.product_smiles,
                        work_dir=reaction_root,
                        skip_steps=skip_steps_reaction,
                        precursor_smiles=(representative.meta or {}).get("precursor_smiles"),
                        leaving_group_key=(representative.meta or {}).get("leaving_small_molecule_key"),
                        small_molecular_keys=(representative.meta or {}).get("small_molecular_keys"),
                        reaction_profile=(representative.meta or {}).get("reaction_profile") or run_cfg.get("reaction_profile"),
                        cleaner_data=(representative.meta or {}).get("cleaner_data")
                        if isinstance((representative.meta or {}).get("cleaner_data"), dict)
                        else None,
                    )
                    results.append(result)
            else:
                skip_steps_reaction = [s for s in list(skip_steps) if str(s).strip()]
                result = hunter.run_pipeline(
                    product_smiles=representative.product_smiles,
                    work_dir=reaction_root,
                    skip_steps=skip_steps_reaction,
                    precursor_smiles=(representative.meta or {}).get("precursor_smiles"),
                    leaving_group_key=(representative.meta or {}).get("leaving_small_molecule_key"),
                    small_molecular_keys=(representative.meta or {}).get("small_molecular_keys"),
                    reaction_profile=(representative.meta or {}).get("reaction_profile") or run_cfg.get("reaction_profile"),
                    cleaner_data=(representative.meta or {}).get("cleaner_data")
                    if isinstance((representative.meta or {}).get("cleaner_data"), dict)
                    else None,
                )
                results.append(result)

            try:
                checkpoint_mgr = CheckpointManager(reaction_root)
                s3_dir = reaction_root / "S3_TS"
                s3_complete = s3_dir.exists() and checkpoint_mgr.is_step3_complete(s3_dir, hunter.config)

                updated = read_json_dict(manifest_path)
                if isinstance(updated, dict):
                    updated["status"] = "COMPLETE" if s3_complete else "PARTIAL"
                    if not updated.get("theory_signature"):
                        updated["theory_signature"] = _theory_signature_string(hunter.config)
                    write_json(manifest_path, updated)
            except Exception as exc:
                hunter.logger.debug(f"Failed to update manifest status for {reaction_id}: {exc}")

            s0_cfg = hunter.config.get("s0", {}) or {}
            dr_cfg = (s0_cfg.get("dr_completion", {}) or {})
            dr_enabled = bool(dr_cfg.get("enabled", True))
            if dr_enabled:
                dr_plan = _load_dr_branch_plan(reaction_root)
                if dr_plan and isinstance(dr_plan, dict):
                    branches = dr_plan.get("branches") or []
                    if isinstance(branches, list) and len(branches) > 1:
                        from rph_core.utils.path_manager import get_branch_root
                        dr_plan_path = str(reaction_root / "S0_Mechanism" / "dr_branch_plan.json")

                        for branch_entry in branches:
                            branch_id = str(branch_entry.get("branch_id", ""))
                            if not branch_id or branch_id == "BR_MAJOR":
                                continue

                            pathway_id = str(branch_entry.get("pathway_id", ""))
                            generation_policy = str(branch_entry.get("generation_policy", ""))
                            branch_product_smiles = (
                                str(branch_entry.get("product_smiles") or "").strip()
                                or representative.product_smiles
                            )
                            flipped_map_numbers = [
                                int(m) for m in (branch_entry.get("flipped_map_numbers") or [])
                            ]
                            fixed_stereocenters = [
                                int(m) for m in (branch_entry.get("fixed_stereocenters") or [])
                            ]
                            notes = [str(n) for n in (branch_entry.get("notes") or [])]

                            branch_root = get_branch_root(reaction_root, branch_id)
                            _write_branch_manifest(
                                branch_root=branch_root,
                                branch_id=branch_id,
                                pathway_id=pathway_id,
                                product_smiles=branch_product_smiles,
                                generation_policy=generation_policy,
                                flipped_map_numbers=flipped_map_numbers,
                                fixed_stereocenters=fixed_stereocenters,
                                notes=notes,
                                dr_plan_path=dr_plan_path,
                                parent_product_smiles=representative.product_smiles,
                            )

                            branch_skip_steps = ["s0"] + [s for s in list(skip_steps) if str(s).strip()]
                            hunter.logger.info(
                                f"  DR branch {branch_id}: product_smiles={branch_product_smiles}, "
                                f"work_dir={branch_root}"
                            )
                            try:
                                branch_result = hunter.run_pipeline(
                                    product_smiles=branch_product_smiles,
                                    work_dir=branch_root,
                                    skip_steps=branch_skip_steps,
                                    precursor_smiles=(representative.meta or {}).get("precursor_smiles"),
                                    leaving_group_key=(representative.meta or {}).get("leaving_small_molecule_key"),
                                    small_molecular_keys=(representative.meta or {}).get("small_molecular_keys"),
                                    reaction_profile=(representative.meta or {}).get("reaction_profile") or run_cfg.get("reaction_profile"),
                                    cleaner_data=(representative.meta or {}).get("cleaner_data")
                                    if isinstance((representative.meta or {}).get("cleaner_data"), dict)
                                    else None,
                                )
                                results.append(branch_result)
                                hunter.logger.info(f"    ✓ DR branch {branch_id} complete → {branch_result.features_csv}")

                                with open(branch_root / "branch_manifest.json", "r", encoding="utf-8") as bf:
                                    bm = json.load(bf)
                                bm["status"] = "COMPLETE" if branch_result.success else "FAILED"
                                with open(branch_root / "branch_manifest.json", "w", encoding="utf-8") as bf:
                                    json.dump(bm, bf, indent=2, ensure_ascii=False)
                            except Exception as exc:
                                hunter.logger.warning(f"    ✗ DR branch {branch_id} failed: {exc}")

                        try:
                            from rph_core.steps.dr_aggregator import DRAggregator

                            dr_aggregator = DRAggregator()
                            dr_prediction_path = reaction_features_dir / "dr_prediction.json"
                            dr_aggregator.write(
                                reaction_root=reaction_root,
                                output_path=dr_prediction_path,
                                temperature_k=default_temperature_k,
                            )
                            hunter.logger.info(f"  DR prediction written → {dr_prediction_path}")
                        except Exception as exc:
                            hunter.logger.warning(f"DR aggregation failed for {reaction_id}: {exc}")

        for condition_task in condition_tasks:
            condition_id = str(getattr(condition_task, "condition_id", f"COND_{condition_task.rx_id}"))
            row_id = str(getattr(condition_task, "row_id", condition_task.rx_id))
            condition_root = reaction_root / "conditions" / condition_id
            condition_root.mkdir(parents=True, exist_ok=True)

            cleaner_data = (condition_task.meta or {}).get("cleaner_data")
            cleaner_data = cleaner_data if isinstance(cleaner_data, dict) else {}
            temperature_c = _to_float(cleaner_data.get("temperature_c") or cleaner_data.get("temp_celsius"))
            temperature_k = _to_float(cleaner_data.get("temperature_K") or cleaner_data.get("temperature_k"))
            if temperature_k is None and temperature_c is not None:
                temperature_k = temperature_c + 273.15
            if temperature_k is None:
                temperature_k = default_temperature_k

            condition_manifest = {
                "version": "rph_v3.0.0",
                "condition_id": condition_id,
                "reaction_id": reaction_id,
                "row_id": row_id,
                "temperature_c": temperature_c,
                "temperature_K": temperature_k,
                "solvent": cleaner_data.get("solvent") or "DCM",
                "catalyst": cleaner_data.get("catalyst"),
                "additive": cleaner_data.get("additive"),
                "oxidant": cleaner_data.get("oxidant"),
                "yield_pct": _to_float(cleaner_data.get("yield_pct") or cleaner_data.get("yield")),
                "dr_major": _to_float(cleaner_data.get("dr_major")),
                "dr_minor": _to_float(cleaner_data.get("dr_minor")),
                "ee_pct": _to_float(cleaner_data.get("ee_pct") or cleaner_data.get("ee")),
                "source_ref": cleaner_data.get("source_ref"),
            }
            write_json(condition_root / "condition_manifest.json", condition_manifest)

            if not run_cfg.get("dry_run", False):
                try:
                    ConditionThermoCalculator(
                        config=hunter.config,
                        reaction_root=reaction_root,
                        condition_root=condition_root,
                    ).run()
                except Exception as e:
                    hunter.logger.warning(f"Condition thermo failed for {condition_id}: {e}")

                try:
                    ConditionFeatureMerger(
                        reaction_root=reaction_root,
                        condition_root=condition_root,
                    ).run()
                except Exception as e:
                    hunter.logger.warning(f"Condition merge failed for {condition_id}: {e}")

    return results


def main():
    """命令行主入口"""
    import argparse

    parser = argparse.ArgumentParser(
        description=f"ReactionProfileHunter V{__version__} - 过渡态搜索与特征提取"
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
    parser.add_argument(
        '--skip-steps',
        type=str,
        default=None,
        help='跳过指定步骤，逗号分隔。例如: --skip-steps s2,s3,s4 仅运行 S1'
    )
    parser.add_argument(
        '--skip-features',
        action='store_true',
        default=False,
        help='跳过 Step 4 特征提取 (已默认关闭，此参数用于向后兼容)'
    )
    parser.add_argument(
        '--features',
        action='store_true',
        default=False,
        help='启用 Step 4 特征提取 (默认不执行 S4)'
    )

    args = parser.parse_args()

    try:
        hunter = ReactionProfileHunter(
            config_path=Path(args.config) if args.config else None,
            log_level=args.log_level
        )

        run_cfg = _resolve_run_config(hunter.config, args)
        hunter.logger.info("Run config resolved. Starting _run_tasks_v3...")
        results = _run_tasks_v3(hunter, run_cfg)
        hunter.logger.info("_run_tasks_v3 finished.")

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
