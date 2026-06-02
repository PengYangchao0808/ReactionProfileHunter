"""
Checkpoint Manager - 断点续传支持
====================================

管理Reaction Profile Hunter的断点续传功能
P0/P1: atomic write, file locks, provenance rehydrate

Author: QCcalc Team
Date: 2026-01-10 / Updated 2026-05
"""

import json
import logging
import os
import time
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Any, Tuple
from datetime import datetime

from rph_core.utils.json_io import write_json as atomic_write_json
from rph_core.utils.layout_contract import (
    canonical_output_files,
    check_step_minimal_complete,
    iter_step_ids,
    resolve_required_files,
    resolve_step_dir,
    seed_steps_template,
)

logger = logging.getLogger(__name__)


class FileRunLock:
    def __init__(self, lock_path: Path, stale_after_sec: int = 21600):
        self.lock_path = Path(lock_path)
        self.stale_after_sec = stale_after_sec
        self._owned = False

    def acquire(self) -> bool:
        self.lock_path.parent.mkdir(parents=True, exist_ok=True)
        if self.lock_path.exists():
            if self._is_stale():
                stale_path = self.lock_path.with_suffix(
                    self.lock_path.suffix + f".stale.{int(time.time())}"
                )
                self.lock_path.rename(stale_path)
                logger.warning("Stale lock replaced: %s -> %s", self.lock_path.name, stale_path.name)
            else:
                return False
        try:
            fd = os.open(str(self.lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            with os.fdopen(fd, "w") as f:
                json.dump({"pid": os.getpid(), "started_at": datetime.now().isoformat()}, f)
            self._owned = True
            return True
        except FileExistsError:
            return False

    def release(self):
        if not self._owned:
            return
        try:
            if self.lock_path.exists():
                self.lock_path.unlink()
        except OSError:
            pass
        self._owned = False

    def _is_stale(self) -> bool:
        try:
            mtime = self.lock_path.stat().st_mtime
            return time.time() - mtime > self.stale_after_sec
        except OSError:
            return True

    def __enter__(self):
        if not self.acquire():
            raise RuntimeError(f"Could not acquire lock: {self.lock_path}")
        return self

    def __exit__(self, *args):
        self.release()
        return False


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

    reaction_id: Optional[str] = None

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
    PROVENANCE_S2_FILENAME = "step2_provenance.json"
    PROVENANCE_S3_FILENAME = "step3_provenance.json"

    def __init__(self, work_dir: Path):
        self.work_dir = Path(work_dir)
        self.state_file = self.work_dir / self.STATE_FILENAME
        self._state_lock_path = self.work_dir / f".{self.STATE_FILENAME}.lock"
        self._run_lock_path = self.work_dir / ".rph_run.lock"
        self._run_lock: Optional[FileRunLock] = None

        self.logger = logging.getLogger(f"{__name__}[{work_dir.name}]")

    @property
    def provenance_s2_path(self) -> Path:
        return resolve_step_dir(self.work_dir, "s2") / self.PROVENANCE_S2_FILENAME

    @property
    def provenance_s3_path(self) -> Path:
        return resolve_step_dir(self.work_dir, "s3") / self.PROVENANCE_S3_FILENAME

    def acquire_run_lock(self) -> bool:
        if self._run_lock is None:
            self._run_lock = FileRunLock(self._run_lock_path)
        return self._run_lock.acquire()

    def release_run_lock(self):
        if self._run_lock is not None:
            self._run_lock.release()
            self._run_lock = None

    def _locked_update(self, updater):
        with FileRunLock(self._state_lock_path, stale_after_sec=300):
            state = self.load_state()
            new_state = updater(state)
            self.save_state(new_state)
            return new_state

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
        state.last_update = datetime.now().isoformat()
        atomic_write_json(self.state_file, state.to_dict())
        self.logger.debug("✓ 状态已保存: %s", self.state_file)

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
        - step3 TS rescue policy and Gaussian TS keywords
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
            },
            'step3_ts_rescue': {
                'gaussian_keywords': {
                    'berny': step3_cfg.get('gaussian_keywords', {}).get('berny', ''),
                    'ts_rescue': step3_cfg.get('gaussian_keywords', {}).get('ts_rescue', ''),
                },
                'disable_qst2_rescue': step3_cfg.get('disable_qst2_rescue', True),
                'policy': step3_cfg.get('ts_rescue_policy', {}),
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

        path_search_cfg = step2_cfg.get("path_search", {}) or {}
        if path_search_cfg.get("enabled") is True:
            path_search_mode = "always"
        elif path_search_cfg.get("enabled") is False:
            path_search_mode = "never"
        else:
            path_search_mode = str(path_search_cfg.get("mode", "rescue")).strip().lower()
        reaction_profiles = config.get("reaction_profiles", {}) or {}
        profile_cfg = reaction_profiles.get(str(reaction_profile), {}) if reaction_profile else {}

        canonical_bonds = sorted((min(int(i), int(j)), max(int(i), int(j))) for i, j in forming_bonds)
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
                "mode": path_search_mode,
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

    def _rehydrate_from_provenance(
        self,
        product_smiles: str,
        config: Dict[str, Any],
    ) -> Optional[PipelineState]:
        s2_prov = self._load_provenance(self.provenance_s2_path)
        s3_prov = self._load_provenance(self.provenance_s3_path)

        now = datetime.now().isoformat()
        state = PipelineState(
            product_smiles=product_smiles,
            work_dir=str(self.work_dir),
            start_time=now,
            last_update=now,
            steps={},
            config_snapshot=config,
        )
        self._ensure_seeded_steps(state)

        s2_complete = bool(s2_prov and check_step_minimal_complete(self.work_dir, "s2", config))
        s3_complete = bool(s3_prov and check_step_minimal_complete(self.work_dir, "s3", config))
        s4_path = resolve_step_dir(self.work_dir, "s4") / "features_raw.csv"

        if s2_complete:
            outputs = canonical_output_files(self.work_dir, "s2")
            s2_prov_data = s2_prov or {}
            meta = dict(s2_prov_data.get("metadata", {}) or {})
            meta["step2_signature"] = s2_prov_data.get("step2_signature")
            meta["forming_bonds"] = s2_prov_data.get("forming_bonds")
            meta["_rehydrate_source"] = self.PROVENANCE_S2_FILENAME
            state.steps["step_s2"] = StepCheckpoint(
                step_name="s2", completed=True, timestamp=now,
                output_files=outputs, metadata=meta,
            )
            self.logger.info("  ✓ 从 provenance 回填 Step2")

        if s3_complete:
            outputs = canonical_output_files(self.work_dir, "s3")
            s3_prov_data = s3_prov or {}
            meta = dict(s3_prov_data.get("metadata", {}) or {})
            meta["step3_signature"] = s3_prov_data.get("step3_signature")
            meta["input_hashes"] = s3_prov_data.get("input_hashes")
            meta["upstream_step2_signature"] = s3_prov_data.get("upstream_step2_signature")
            meta["_rehydrate_source"] = self.PROVENANCE_S3_FILENAME
            state.steps["step_s3"] = StepCheckpoint(
                step_name="s3", completed=True, timestamp=now,
                output_files=outputs, metadata=meta,
            )
            self.logger.info("  ✓ 从 provenance 回填 Step3")

        if s4_path.exists() and s4_path.stat().st_size > 0:
            outputs = canonical_output_files(self.work_dir, "s4")
            state.steps["step_s4"] = StepCheckpoint(
                step_name="s4", completed=True, timestamp=now,
                output_files=outputs, metadata={"_rehydrate_source": "features_raw.csv"},
            )
            self.logger.info("  ✓ 回填 Step4 (features_raw.csv)")

        completed_any = s2_complete or s3_complete or (s4_path.exists() and s4_path.stat().st_size > 0)
        if not completed_any:
            return None

        self.logger.info("✓ 从 provenance 精确恢复状态完成")
        return state

    @staticmethod
    def _load_provenance(path: Path) -> Optional[Dict[str, Any]]:
        if not path.exists():
            return None
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return None

    def rehydrate_state_from_artifacts(
        self,
        product_smiles: str,
        config: Dict[str, Any],
        policy: str = "best_effort"
    ) -> Optional[PipelineState]:
        # P0-3: Provenance-first recovery
        prov_state = self._rehydrate_from_provenance(product_smiles, config)
        if prov_state is not None:
            self.save_state(prov_state)
            return prov_state

        if policy != "best_effort":
            return None

        self.logger.info("Provenance 文件缺失，降级到 best-effort artifact scan...")
        self.logger.info("正在从物理产物回填状态 (Policy: %s)...", policy)

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
        def _do(state):
            if state is None:
                self.logger.warning("状态文件不存在，无法标记步骤")
                return state
            step_key = f"step_{step_name}"
            state.steps[step_key] = StepCheckpoint(
                step_name=step_name,
                completed=True,
                timestamp=datetime.now().isoformat(),
                output_files=output_files,
                metadata=metadata or {}
            )
            return state
        self._locked_update(_do)

    def mark_step_in_progress(
        self,
        step_name: str,
        phase: str,
        output_files: Optional[Dict[str, str]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        def _do(state):
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
            return state
        self._locked_update(_do)

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
        def _do(_state):
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
            return state
        self._locked_update(_do)


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
