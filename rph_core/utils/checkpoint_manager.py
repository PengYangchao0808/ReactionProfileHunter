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

from rph_core.utils.layout_contract import (
    canonical_output_files,
    check_step_minimal_complete,
    iter_step_ids,
    resolve_required_files,
    resolve_step_dir,
    seed_steps_template,
)

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

    def _ensure_seeded_steps(self, state: PipelineState) -> None:
        seeded = seed_steps_template()
        for key, value in seeded.items():
            if key not in state.steps:
                state.steps[key] = StepCheckpoint(
                    step_name=value["step_name"],
                    completed=value["completed"],
                    timestamp=value["timestamp"],
                    output_files=value["output_files"],
                    metadata=value["metadata"],
                )

    def save_state(self, state: PipelineState):
        """
        保存pipeline状态

        Args:
            state: PipelineState对象
        """
        state.last_update = datetime.now().isoformat()

        with open(self.state_file, 'w', encoding='utf-8') as f:
            json.dump(state.to_dict(), f, indent=2, ensure_ascii=False)

        self.logger.debug(f"✓ 状态已保存: {self.state_file}")

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
            self._ensure_seeded_steps(state)
            self.logger.debug(f"✓ 状态已加载: {self.state_file}")

            # 打印已完成步骤
            completed_steps = [k for k, v in state.steps.items() if v.completed]
            self.logger.debug(f"  已完成步骤: {completed_steps}")

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
        - step3.reactant_opt: charge, multiplicity, enable_nbo (TSOptimizer 实际使用的配置)
        - version: rph_core version
        """
        from rph_core.version import __version__

        theory_opt = config.get('theory', {}).get('optimization', {})
        theory_sp = config.get('theory', {}).get('single_point', {})
        step3_cfg = config.get('step3', {})
        # FIX: 使用 reactant_opt 而非 intermediate_opt (TSOptimizer 实际读取的配置)
        reactant_opt = step3_cfg.get('reactant_opt', {})

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
            # FIX: 改用 reactant_opt (TSOptimizer 实际使用)
            'step3_reactant_opt': {
                'charge': reactant_opt.get('charge', 0),
                'multiplicity': reactant_opt.get('multiplicity', 1),
                'enable_nbo': reactant_opt.get('enable_nbo', False),
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

        path_search_cfg = step2_cfg.get("path_search", {})
        reaction_profiles = config.get("reaction_profiles", {}) or {}
        profile_cfg = reaction_profiles.get(str(reaction_profile), {}) if reaction_profile else {}

        canonical_bonds = sorted((min(i, j), max(i, j)) for i, j in forming_bonds)
        product_hash = self.compute_file_hash(product_xyz) or ""

        return {
            "version": __version__,
            "reaction_profile": reaction_profile or "",
            "s2_strategy": profile_cfg.get("s2_strategy", "retro_scan"),
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
            "path_search": {
                "enabled": path_search_cfg.get("enabled", False),
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
        self._ensure_seeded_steps(state)

        now = datetime.now().isoformat()
        completed_any = False
        for step_id in iter_step_ids():
            step_key = f"step_{step_id}"
            if check_step_minimal_complete(self.work_dir, step_id, config):
                outputs = canonical_output_files(self.work_dir, step_id)
                metadata: Dict[str, Any] = {}
                if step_id == "s2":
                    scan_profile = resolve_step_dir(self.work_dir, "s2") / "scan_profile.json"
                    if scan_profile.exists():
                        metadata["scan_profile_json"] = str(scan_profile)
                        try:
                            with open(scan_profile, "r", encoding="utf-8") as f:
                                scan_data = json.load(f)
                            raw_bonds = scan_data.get("forming_bonds") or []
                            parsed_bonds = []
                            for pair in raw_bonds:
                                if isinstance(pair, list) and len(pair) == 2:
                                    parsed_bonds.append((int(pair[0]), int(pair[1])))
                            if parsed_bonds:
                                metadata["forming_bonds"] = [list(pair) for pair in parsed_bonds]
                                product_xyz = resolve_required_files(self.work_dir, "s1").get("product_xyz")
                                if product_xyz and product_xyz.exists():
                                    metadata["step2_signature"] = self.compute_step2_signature(
                                        config=config,
                                        product_xyz=product_xyz,
                                        forming_bonds=tuple(parsed_bonds),
                                        reaction_profile=None,
                                        scan_config=None,
                                    )
                        except Exception as exc:
                            self.logger.debug(f"Failed reading S2 scan_profile during rehydrate: {exc}")
                if step_id == "s3":
                    metadata["step3_signature"] = self._compute_step3_signature(config)
                state.steps[step_key] = StepCheckpoint(
                    step_name=step_id,
                    completed=True,
                    timestamp=now,
                    output_files=outputs,
                    metadata=metadata,
                )
                completed_any = True
                self.logger.info(f"  ✓ 回填 Step {step_id.upper()}")

        for step_id in iter_step_ids():
            step_key = f"step_{step_id}"
            if state.steps.get(step_key, StepCheckpoint(step_id, False, "", {})).completed:
                continue
            status_file = resolve_step_dir(self.work_dir, step_id) / ".rph_step_status.json"
            if status_file.exists():
                try:
                    with open(status_file, "r", encoding="utf-8") as f:
                        status_data = json.load(f)
                    metadata = {
                        "phase": status_data.get("description") or status_data.get("phase") or "in_progress",
                        "status_file": str(status_file),
                    }
                    state.steps[step_key] = StepCheckpoint(
                        step_name=step_id,
                        completed=False,
                        timestamp=now,
                        output_files=canonical_output_files(self.work_dir, step_id),
                        metadata=metadata,
                    )
                except Exception as exc:
                    self.logger.debug(f"Failed reading partial status for {step_id}: {exc}")

        if not completed_any:
            self.logger.warning("未发现可回填的产物")
            return None

        return state

    def is_step4_complete(self, s4_dir: Path, config: Dict[str, Any]) -> bool:
        """
        M2-A: 检查 Step4 是否完成（考虑机制打包完整性）
        """

        if not check_step_minimal_complete(self.work_dir, "s4", config):
            return False

        mech_config = config.get('step4', {}).get('mechanism_packaging', {})
        if not mech_config.get('enabled', False):
            return True

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

    def mark_step_in_progress(
        self,
        step_name: str,
        phase: str,
        output_files: Optional[Dict[str, str]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        state = self.load_state()
        if state is None:
            now = datetime.now().isoformat()
            state = PipelineState(
                product_smiles="",
                work_dir=str(self.work_dir),
                start_time=now,
                last_update=now,
                steps={},
                config_snapshot={},
            )

        self._ensure_seeded_steps(state)
        step_key = f"step_{step_name}"
        base_metadata = dict(metadata or {})
        base_metadata["phase"] = phase
        state.steps[step_key] = StepCheckpoint(
            step_name=step_name,
            completed=False,
            timestamp=datetime.now().isoformat(),
            output_files=output_files or {},
            metadata=base_metadata,
        )
        self.save_state(state)

    def mark_step_failed_partial(
        self,
        step_name: str,
        phase: str,
        error_message: str,
        output_files: Optional[Dict[str, str]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        failure_metadata = dict(metadata or {})
        failure_metadata.update({"phase": phase, "error_stage": phase, "error_message": error_message})
        self.mark_step_in_progress(
            step_name=step_name,
            phase=phase,
            output_files=output_files,
            metadata=failure_metadata,
        )

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
                "step_s0": StepCheckpoint("s0", False, "", {}),
                "step_s1": StepCheckpoint("s1", False, "", {}),
                "step_s2": StepCheckpoint("s2", False, "", {}),
                "step_s3": StepCheckpoint("s3", False, "", {}),
                "step_s4": StepCheckpoint("s4", False, "", {})
            },
            config_snapshot=config
        )

        self._ensure_seeded_steps(state)

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
