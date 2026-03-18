"""
Checkpoint Manager - 断点续传支持
====================================

管理Reaction Profile Hunter的断点续传功能

Author: QCcalc Team
Date: 2026-01-10
"""

import json
import logging
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Any, Tuple
from datetime import datetime

logger = logging.getLogger(__name__)


@dataclass
class StepCheckpoint:
    """单个步骤的checkpoint信息"""
    step_name: str
    completed: bool
    timestamp: str
    output_files: Dict[str, str]  # {"product_xyz": "path/to/file.xyz"}
    metadata: Optional[Dict[str, Any]] = None

    def to_dict(self) -> dict[str, Any]:
        """转换为字典"""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> 'StepCheckpoint':
        """从字典创建"""
        return cls(**data)


@dataclass
class PipelineState:
    """整个pipeline的状态"""
    product_smiles: str
    work_dir: str
    start_time: str
    last_update: str

    # 各步骤状态
    steps: Dict[str, StepCheckpoint]

    # 全局配置
    config_snapshot: Optional[Dict[str, Any]] = None

    def to_dict(self) -> dict[str, Any]:
        """转换为字典"""
        data = asdict(self)
        data['steps'] = {k: v.to_dict() for k, v in self.steps.items()}
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> 'PipelineState':
        """从字典创建"""
        steps_data = data.pop('steps', {})
        steps = {k: StepCheckpoint.from_dict(v) for k, v in steps_data.items()}
        return cls(steps=steps, **data)


class CheckpointManager:
    """
    Checkpoint管理器 - 断点续传核心

    功能:
    1. 保存/加载.state文件
    2. 检查步骤完成状态
    3. 恢复pipeline执行
    """

    STATE_FILENAME = "pipeline.state"

    def __init__(self, work_dir: Path):
        """
        初始化Checkpoint管理器

        Args:
            work_dir: 工作目录
        """
        self.work_dir = Path(work_dir)
        self.state_file = self.work_dir / self.STATE_FILENAME

        self.logger = logging.getLogger(f"{__name__}[{work_dir.name}]")

    def save_state(self, state: PipelineState):
        """
        保存pipeline状态

        Args:
            state: PipelineState对象
        """
        state.last_update = datetime.now().isoformat()

        with open(self.state_file, 'w', encoding='utf-8') as f:
            json.dump(state.to_dict(), f, indent=2, ensure_ascii=False)

        self.logger.info(f"✓ 状态已保存: {self.state_file}")

    def load_state(self) -> Optional[PipelineState]:
        """
        加载pipeline状态

        Returns:
            PipelineState对象，如果文件不存在返回None
        """
        if not self.state_file.exists():
            self.logger.warning(f"状态文件不存在: {self.state_file}")
            return None

        try:
            with open(self.state_file, 'r', encoding='utf-8') as f:
                data = json.load(f)

            state = PipelineState.from_dict(data)
            self.logger.info(f"✓ 状态已加载: {self.state_file}")

            # 打印已完成步骤
            completed_steps = [k for k, v in state.steps.items() if v.completed]
            self.logger.info(f"  已完成步骤: {completed_steps}")

            return state

        except Exception as e:
            self.logger.error(f"加载状态文件失败: {e}")
            return None

    def is_step_completed(self, step_name: str) -> bool:
        """
        检查步骤是否已完成

        Args:
            step_name: 步骤名称 (s1, s2, s3, s4)

        Returns:
            是否已完成
        """
        state = self.load_state()
        if state is None:
            return False

        step_key = f"step_{step_name}"
        if step_key not in state.steps:
            return False

        return state.steps[step_key].completed

    def is_step3_complete(
        self,
        s3_dir: Path,
        config: Dict[str, Any],
        check_signature: bool = True,
        input_hashes: Optional[Dict[str, str]] = None,
        upstream_step2_signature: Optional[Dict[str, Any]] = None,
    ) -> bool:
        """
        V6.3: 检查 Step3 是否完成且可安全复用
        """
        import hashlib

        if not self.is_step_completed('s3'):
            return False

        state = self.load_state()
        if state is None:
            return False

        step_s3 = state.steps.get('step_s3')
        if step_s3 is None or not step_s3.completed:
            return False

        output_files = step_s3.output_files or {}
        metadata = step_s3.metadata or {}

        ts_final_xyz = output_files.get('ts_final_xyz')
        if not ts_final_xyz or not Path(ts_final_xyz).exists():
            self.logger.debug("S3 checkpoint: ts_final_xyz missing or not exists")
            return False

        ts_final_path = Path(ts_final_xyz)

        try:
            with open(ts_final_xyz, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.isdigit():
                    self.logger.debug("S3 checkpoint: ts_final_xyz invalid format")
                    return False
        except OSError:
            return False

        sp_meta_path = output_files.get('sp_matrix_metadata_json') or str(s3_dir / "sp_matrix_metadata.json")
        if not Path(sp_meta_path).exists():
            self.logger.debug("S3 checkpoint: sp_matrix_metadata.json missing")
            return False

        try:
            with open(sp_meta_path, 'r') as f:
                sp_meta = json.load(f)
            if 'e_ts' not in sp_meta or 'e_reactant' not in sp_meta:
                self.logger.debug("S3 checkpoint: sp_matrix_metadata missing energy fields")
                return False
        except (json.JSONDecodeError, OSError) as e:
            self.logger.debug(f"S3 checkpoint: sp_matrix_metadata.json parse failed: {e}")
            return False

        if check_signature:
            cached_sig = metadata.get('step3_signature', {})
            current_sig = self._compute_step3_signature(config)

            if cached_sig != current_sig:
                self.logger.info(
                    f"S3 checkpoint: signature mismatch (cached vs current), will recompute S3"
                )
                return False

        if input_hashes is not None:
            cached_hashes = metadata.get('input_hashes', {})
            for key in ['ts_guess', 'intermediate', 'product']:
                if key in input_hashes:
                    if cached_hashes.get(key) != input_hashes[key]:
                        self.logger.info(
                            f"S3 checkpoint: input hash mismatch for {key}, will recompute S3"
                        )
                        return False

        if upstream_step2_signature is not None:
            cached_upstream = metadata.get('upstream_step2_signature')
            if cached_upstream != upstream_step2_signature:
                self.logger.info(
                    "S3 checkpoint: upstream Step2 signature mismatch, will recompute S3"
                )
                return False

        self.logger.info("S3 checkpoint: all checks passed, can reuse")
        return True

    def _compute_step3_signature(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        计算 Step3 配置签名用于复用判定

        包含字段：
        - theory.optimization: engine, method, basis, dispersion, route, rescue_route, nproc, mem
        - theory.single_point: engine, method, basis, aux_basis, solvent, nproc, maxcore
        - step3.intermediate_opt: charge, multiplicity, enable_nbo
        - version: rph_core version
        """
        from rph_core.version import __version__

        theory_opt = config.get('theory', {}).get('optimization', {})
        theory_sp = config.get('theory', {}).get('single_point', {})
        step3_cfg = config.get('step3', {})
        intermediate_opt = step3_cfg.get('intermediate_opt', {})

        return {
            'version': __version__,
            'theory_optimization': {
                'engine': theory_opt.get('engine', 'gaussian'),
                'method': theory_opt.get('method', 'B3LYP'),
                'basis': theory_opt.get('basis', 'def2-SVP'),
                'dispersion': theory_opt.get('dispersion', 'GD3BJ'),
                'route': theory_opt.get('route', ''),
                'rescue_route': theory_opt.get('rescue_route', ''),
                'nproc': theory_opt.get('nproc', 16),
                'mem': theory_opt.get('mem', '32GB'),
            },
            'theory_single_point': {
                'engine': theory_sp.get('engine', 'orca'),
                'method': theory_sp.get('method', 'WB97M-V'),
                'basis': theory_sp.get('basis', 'def2-TZVPP'),
                'aux_basis': theory_sp.get('aux_basis', 'def2/J'),
                'solvent': theory_sp.get('solvent', 'acetone'),
                'nproc': theory_sp.get('nproc', 16),
                'maxcore': theory_sp.get('maxcore', 4000),
            },
            'step3_intermediate_opt': {
                'charge': intermediate_opt.get('charge', 0),
                'multiplicity': intermediate_opt.get('multiplicity', 1),
                'enable_nbo': intermediate_opt.get('enable_nbo', False),
            }
        }

    def compute_step2_signature(
        self,
        *,
        config: Dict[str, Any],
        product_xyz: Path,
        forming_bonds: Tuple[Tuple[int, int], ...],
        reaction_profile: Optional[str],
        scan_config: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        from rph_core.version import __version__

        step2_cfg = config.get("step2", {}) or {}
        scan_cfg = dict(step2_cfg.get("scan", {}) or {})
        if scan_config:
            scan_cfg.update(scan_config)

        canonical_bonds = sorted((min(i, j), max(i, j)) for i, j in forming_bonds)
        product_hash = self.compute_file_hash(product_xyz) or ""

        return {
            "version": __version__,
            "reaction_profile": reaction_profile or "",
            "product_xyz_hash": product_hash,
            "forming_bonds": [list(pair) for pair in canonical_bonds],
            "scan": {
                "scan_start_distance": scan_cfg.get("scan_start_distance"),
                "scan_end_distance": scan_cfg.get("scan_end_distance"),
                "scan_steps": scan_cfg.get("scan_steps"),
                "scan_mode": scan_cfg.get("scan_mode"),
                "scan_force_constant": scan_cfg.get("scan_force_constant"),
                "min_valid_points": scan_cfg.get("min_valid_points"),
                "reject_boundary_maximum": scan_cfg.get("reject_boundary_maximum"),
                "require_local_peak": scan_cfg.get("require_local_peak"),
            },
        }

    def compute_file_hash(self, file_path: Path) -> Optional[str]:
        """
        计算文件的 SHA256 哈希值

        Args:
            file_path: 文件路径

        Returns:
            哈希字符串或 None
        """
        import hashlib

        if not file_path or not file_path.exists():
            return None

        try:
            with open(file_path, 'rb') as f:
                return hashlib.sha256(f.read()).hexdigest()[:16]
        except OSError as e:
            self.logger.debug(f"Failed to compute hash for {file_path}: {e}")
            return None

    def rehydrate_state_from_artifacts(
        self,
        product_smiles: str,
        config: Dict[str, Any],
        policy: str = "best_effort"
    ) -> Optional[PipelineState]:
        """
        V6.3: Reconstruct pipeline.state from existing artifact files.
        """
        self.logger.info(f"正在从物理产物回填状态 (Policy: {policy})...")

        state = PipelineState(
            product_smiles=product_smiles,
            work_dir=str(self.work_dir),
            start_time=datetime.now().isoformat(),
            last_update=datetime.now().isoformat(),
            steps={},
            config_snapshot=config
        )

        s1_dir = self.work_dir / "S1_ConfGeneration"
        product_min = s1_dir / "product_min.xyz"
        if product_min.exists():
            s1_files = {"product_xyz": str(product_min)}
            chk = s1_dir / "product" / "dft" / "conf_000.chk"
            if chk.exists(): s1_files["product_checkpoint"] = str(chk)
            fchk = s1_dir / "product" / "dft" / "conf_000.fchk"
            if fchk.exists(): s1_files["product_fchk"] = str(fchk)
            log = s1_dir / "product" / "dft" / "conf_000.log"
            if log.exists(): s1_files["product_log"] = str(log)
            thermo = s1_dir / "product" / "dft" / "conformer_thermo.csv"
            if thermo.exists(): s1_files["product_thermo"] = str(thermo)

            state.steps["step_s1"] = StepCheckpoint(
                step_name="s1", completed=True, timestamp=datetime.now().isoformat(),
                output_files=s1_files, metadata={}
            )
            self.logger.info("  ✓ 回填 Step 1")

        s2_candidates = [
            self.work_dir / "S2_GuessGeneration",
            self.work_dir / "S2_Retro",
        ]
        s2_dir = next((d for d in s2_candidates if d.exists() and d.is_dir()), self.work_dir / "S2_Retro")
        ts_guess = s2_dir / "ts_guess.xyz"
        intermediate = s2_dir / "reactant_complex.xyz"
        scan_profile = s2_dir / "scan_profile.json"
        if ts_guess.exists() and intermediate.exists():
            output_files = {
                "ts_guess_xyz": str(ts_guess),
                "reactant_xyz": str(intermediate),
                "reactant_complex_xyz": str(intermediate),
            }

            metadata: Dict[str, Any] = {
                "status": "COMPLETE",
                "ts_guess_confidence": "high",
                "degraded_reasons": [],
                "scan_profile_json": str(scan_profile) if scan_profile.exists() else "",
            }

            forming_bonds_for_sig: Tuple[Tuple[int, int], ...] = tuple()
            if scan_profile.exists():
                try:
                    with open(scan_profile, 'r') as f:
                        scan_data = json.load(f)
                    
                    generation_method = scan_data.get("generation_method", "")
                    metadata["generation_method"] = generation_method
                    
                    if generation_method == "xtb_path_search":
                        ts_quality = scan_data.get("ts_quality", {}) or {}
                        metadata["status"] = str(ts_quality.get("status", "COMPLETE"))
                        metadata["ts_guess_confidence"] = str(ts_quality.get("ts_guess_confidence", "high"))
                        metadata["degraded_reasons"] = list(ts_quality.get("degraded_reasons", []))
                    else:
                        quality = scan_data.get("scan_quality", {}) or {}
                        metadata["status"] = str(quality.get("status", "COMPLETE"))
                        metadata["ts_guess_confidence"] = str(quality.get("ts_guess_confidence", "high"))
                        metadata["degraded_reasons"] = list(quality.get("degraded_reasons", []))
                    
                    raw_bonds = scan_data.get("forming_bonds") or []
                    parsed_bonds = []
                    for pair in raw_bonds:
                        if isinstance(pair, list) and len(pair) == 2:
                            parsed_bonds.append((int(pair[0]), int(pair[1])))
                    if parsed_bonds:
                        forming_bonds_for_sig = tuple(parsed_bonds)
                        metadata["forming_bonds"] = [list(b) for b in parsed_bonds]
                except Exception as exc:
                    self.logger.debug(f"Failed to read S2 scan_profile during rehydrate: {exc}")

            if product_min.exists() and forming_bonds_for_sig:
                try:
                    metadata["step2_signature"] = self.compute_step2_signature(
                        config=config,
                        product_xyz=product_min,
                        forming_bonds=forming_bonds_for_sig,
                        reaction_profile=None,
                        scan_config=None,
                    )
                except Exception as exc:
                    self.logger.debug(f"Failed to compute Step2 signature during rehydrate: {exc}")

            state.steps["step_s2"] = StepCheckpoint(
                step_name="s2", completed=True, timestamp=datetime.now().isoformat(),
                output_files=output_files,
                metadata=metadata,
            )
            self.logger.info("  ✓ 回填 Step 2")

        s3_dir = self.work_dir / "S3_TransitionAnalysis"
        sp_meta_path = s3_dir / "sp_matrix_metadata.json"

        if sp_meta_path.exists():
            ts_final: Optional[Path] = s3_dir / "ts_final.xyz"
            if not ts_final.exists():
                resume_file = s3_dir / "s3_resume.json"
                if resume_file.exists():
                    try:
                        with open(resume_file, 'r') as f:
                            r_data = json.load(f)
                            resume_path = r_data.get("ts_opt", {}).get("optimized_xyz")
                            ts_final = Path(resume_path) if resume_path else None
                    except (json.JSONDecodeError, OSError):
                        ts_final = None

            if ts_final is not None and ts_final.exists():
                signature_file = s3_dir / "step3_signature.json"
                can_rehydrate_s3 = False
                sig_data = {}
                
                if signature_file.exists():
                    try:
                        with open(signature_file, 'r') as f:
                            sig_data = json.load(f)
                        cached_sig = sig_data.get("step3_signature")
                        if cached_sig == self._compute_step3_signature(config):
                            can_rehydrate_s3 = True
                        else:
                            self.logger.warning("  ! Step 3 签名不匹配，不执行回填")
                    except (json.JSONDecodeError, OSError, TypeError) as e:
                        self.logger.warning(f"  ! 读取签名文件失败: {e}")
                elif policy == "best_effort":
                    self.logger.info("  ! Step 3 缺失签名文件，基于 best_effort 强制回填")
                    can_rehydrate_s3 = True
                
                if can_rehydrate_s3:
                    s3_files = {
                        "ts_final_xyz": str(ts_final),
                        "sp_matrix_metadata_json": str(sp_meta_path)
                    }
                    def _find_file(patterns: List[str], base: Path) -> Optional[str]:
                        for p in patterns:
                            matches = list(base.glob(p))
                            if matches: return str(matches[0])
                        return None

                    s3_files["ts_fchk"] = _find_file(["ts_opt/**/*.fchk", "ts_opt/*.fchk"], s3_dir) or ""
                    s3_files["ts_log"] = _find_file(["ts_opt/**/*.log", "ts_opt/*.log"], s3_dir) or ""
                    s3_files["intermediate_fchk"] = _find_file(["S3_intermediate_opt/**/*.fchk", "S3_intermediate_opt/*.fchk"], s3_dir) or ""
                    s3_files["intermediate_log"] = _find_file(["S3_intermediate_opt/**/*.log", "S3_intermediate_opt/*.log"], s3_dir) or ""

                    energies = {}
                    try:
                        with open(sp_meta_path, 'r') as f:
                            m = json.load(f)
                            energies = {"e_ts": m.get("e_ts"), "e_intermediate": m.get("e_intermediate"), "e_product": m.get("e_product")}
                    except (json.JSONDecodeError, OSError) as exc:
                        self.logger.debug(f"Failed to load SP metadata during rehydrate: {exc}")

                    rehydrated_hashes = sig_data.get("input_hashes")
                    if not rehydrated_hashes:
                        rehydrated_hashes = {
                            "ts_guess": self.compute_file_hash(ts_guess) or "",
                            "intermediate": self.compute_file_hash(intermediate) or "",
                            "product": self.compute_file_hash(product_min) or ""
                        }

                    state.steps["step_s3"] = StepCheckpoint(
                        step_name="s3", completed=True, timestamp=datetime.now().isoformat(),
                        output_files=s3_files,
                        metadata={
                            "step3_signature": sig_data.get("step3_signature") or self._compute_step3_signature(config),
                            "input_hashes": rehydrated_hashes,
                            "energies": energies
                        }
                    )
                    self.logger.info("  ✓ 回填 Step 3")

        if not state.steps:
            self.logger.warning("未发现可回填的产物")
            return None

        return state

    def is_step4_complete(self, s4_dir: Path, config: Dict[str, Any]) -> bool:
        """
        M2-A: 检查 Step4 是否完成（考虑机制打包完整性）
        """

        # 如果机制打包未启用，使用原始逻辑
        mech_config = config.get('step4', {}).get('mechanism_packaging', {})
        if not mech_config.get('enabled', False):
            return self.is_step_completed('s4')

        # 检查 features_raw.csv (V6.1: 更新为3文件契约)
        features_raw_csv = s4_dir / "features_raw.csv"
        if not features_raw_csv.exists() or features_raw_csv.stat().st_size == 0:
            return False

        # 检查 mech_index.json
        mech_index_path = s4_dir / "mech_index.json"
        if not mech_index_path.exists():
            return False

        try:
            with open(mech_index_path, 'r', encoding='utf-8') as f:
                mech_index = json.load(f)
        except Exception as e:
            self.logger.warning(f"mech_index.json 读取失败: {e}")
            return False

        # 检查 schema_version
        expected_schema = mech_config.get('schema_version', 'mech_index_v1')
        actual_schema = mech_index.get('schema_version')
        if actual_schema != expected_schema:
            self.logger.warning(
                f"mech_index.schema_version 不匹配: 期望 {expected_schema}, 实际 {actual_schema}"
            )
            return False

        # 所有检查通过
        return True

    def mark_step_completed(
        self,
        step_name: str,
        output_files: Dict[str, str],
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        标记步骤为已完成

        Args:
            step_name: 步骤名称 (s1, s2, s3, s4)
            output_files: 输出文件字典
            metadata: 额外的元数据
        """
        # 加载现有状态
        state = self.load_state()

        if state is None:
            self.logger.warning("状态文件不存在，无法标记步骤")
            return

        # 更新步骤状态
        step_key = f"step_{step_name}"
        state.steps[step_key] = StepCheckpoint(
            step_name=step_name,
            completed=True,
            timestamp=datetime.now().isoformat(),
            output_files=output_files,
            metadata=metadata or {}
        )

        # 保存状态
        self.save_state(state)

    def get_step_output(self, step_name: str, output_key: str) -> Optional[str]:
        """
        获取步骤的输出文件路径

        Args:
            step_name: 步骤名称
            output_key: 输出文件键 (如 "product_xyz")

        Returns:
            文件路径，如果不存在返回None
        """
        state = self.load_state()
        if state is None:
            return None

        step_key = f"step_{step_name}"
        if step_key not in state.steps:
            return None

        return state.steps[step_key].output_files.get(output_key)

    def get_step_metadata(self, step_name: str, meta_key: str) -> Optional[Any]:
        state = self.load_state()
        if state is None:
            return None
        step_key = f"step_{step_name}"
        if step_key not in state.steps:
            return None
        metadata = state.steps[step_key].metadata or {}
        return metadata.get(meta_key)

    def initialize_state(self, product_smiles: str, config: Dict[str, Any]):
        """
        初始化pipeline状态

        Args:
            product_smiles: 产物SMILES
            config: 配置字典
        """
        state = PipelineState(
            product_smiles=product_smiles,
            work_dir=str(self.work_dir),
            start_time=datetime.now().isoformat(),
            last_update=datetime.now().isoformat(),
            steps={
                "step_s1": StepCheckpoint("s1", False, "", {}),
                "step_s2": StepCheckpoint("s2", False, "", {}),
                "step_s3": StepCheckpoint("s3", False, "", {}),
                "step_s4": StepCheckpoint("s4", False, "", {})
            },
            config_snapshot=config
        )

        self.save_state(state)
        self.logger.info("✓ Pipeline状态已初始化")


def load_checkpoint_state(work_dir: Path) -> Optional[PipelineState]:
    """
    便捷函数：加载checkpoint状态

    Args:
        work_dir: 工作目录

    Returns:
        PipelineState对象或None
    """
    manager = CheckpointManager(work_dir)
    return manager.load_state()


def save_checkpoint_state(
    work_dir: Path,
    product_smiles: str,
    config: Dict[str, Any]
):
    """
    便捷函数：初始化checkpoint状态

    Args:
        work_dir: 工作目录
        product_smiles: 产物SMILES
        config: 配置字典
    """
    manager = CheckpointManager(work_dir)
    manager.initialize_state(product_smiles, config)
