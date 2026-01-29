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
from typing import Dict, List, Optional, Any
from datetime import datetime

logger = logging.getLogger(__name__)


@dataclass
class StepCheckpoint:
    """单个步骤的checkpoint信息"""
    step_name: str
    completed: bool
    timestamp: str
    output_files: Dict[str, str]  # {"product_xyz": "path/to/file.xyz"}
    metadata: Dict[str, Any] = None

    def to_dict(self) -> dict:
        """转换为字典"""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> 'StepCheckpoint':
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
    config_snapshot: Dict[str, Any] = None

    def to_dict(self) -> dict:
        """转换为字典"""
        data = asdict(self)
        data['steps'] = {k: v.to_dict() for k, v in self.steps.items()}
        return data

    @classmethod
    def from_dict(cls, data: dict) -> 'PipelineState':
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

    def is_step4_complete(self, s4_dir: Path, config: Dict[str, Any]) -> bool:
        """
        M2-A: 检查 Step4 是否完成（考虑机制打包完整性）
        当 mechanism_packaging.enabled=true 时，必须同时满足：
        - features_raw.csv 存在且有效
        - mech_index.json 存在
        - mech_index.schema_version == config 中期望版本

        Args:
            s4_dir: S4_Data 目录路径
            config: 配置字典

        Returns:
            Step4 是否完成
        """
        import json

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
