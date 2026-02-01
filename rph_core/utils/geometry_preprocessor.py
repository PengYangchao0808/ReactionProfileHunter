"""
几何结构预处理模块 (Geometry Preprocessor)
==========================================

在 Gaussian/ORCA 优化前自动检测并修复原子重叠问题

策略 (Houk 建议):
1. 检测原子间距 < overlap_threshold (默认 1.0 Å)
2. 若检测到重叠，触发 xTB 预优化 (GFN2-xTB)
3. 将 xTB 优化后的结构传递给 B3LYP/def2-SVP
4. 最终在优化后的结构上做高精度单点能

Author: ReactionProfileHunter Team
Date: 2026-01-31
Purpose: 避免 "Don't use a diamond saw to cut a sandwich" 问题
"""

from pathlib import Path
from typing import Optional, Tuple
import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.xtb_runner import XTBRunner
from rph_core.utils.data_types import QCResult


class GeometryPreprocessor(LoggerMixin):
    """
    几何结构预处理器
    
    功能:
    - 检测原子重叠 (overlap detection)
    - xTB 预优化 (pre-optimization)
    - 结构修复 (geometry repair)
    
    Example:
        >>> preprocessor = GeometryPreprocessor(config)
        >>> result = preprocessor.preprocess(xyz_file, output_dir)
        >>> if result.success:
        ...     optimized_xyz = result.coordinates
    """
    
    def __init__(self, config: dict):
        """
        初始化预处理器
        
        Args:
            config: 配置字典，需包含 theory.preoptimization 节
        """
        self.config = config
        self.preopt_cfg = config.get('theory', {}).get('preoptimization', {})
        
        self.enabled = self.preopt_cfg.get('enabled', True)
        self.overlap_threshold = self.preopt_cfg.get('overlap_threshold', 1.0)
        self.gfn_level = self.preopt_cfg.get('gfn_level', 2)
        self.solvent = self.preopt_cfg.get('solvent', 'acetone')
        
        if self.enabled:
            self.logger.info(f"GeometryPreprocessor enabled (overlap_threshold={self.overlap_threshold} Å)")
        else:
            self.logger.info("GeometryPreprocessor disabled (preoptimization.enabled=false)")
    
    def preprocess(
        self,
        xyz_file: Path,
        output_dir: Path,
        charge: int = 0,
        uhf: int = 0,
        force_preopt: bool = False
    ) -> QCResult:
        """
        预处理几何结构：检测重叠 → 必要时执行 xTB 预优化
        
        Args:
            xyz_file: 输入 XYZ 文件路径
            output_dir: 输出目录
            charge: 分子电荷
            uhf: 未配对电子数
            force_preopt: 强制执行预优化（忽略重叠检测）
        
        Returns:
            QCResult 对象，包含预处理后的结构 (或原始结构若无需预优化)
            - success=True 表示成功（无论是否执行了预优化）
            - coordinates 包含最终使用的 XYZ 文件路径
            - error_message 包含状态信息：
                - "preopt_disabled": 预优化被禁用
                - "no_overlap": 无原子重叠，跳过预优化
                - "preopt_success": 预优化成功
                - None: 预优化失败但不致命
        """
        xyz_file = Path(xyz_file)
        output_dir = Path(output_dir)
        
        if not self.enabled and not force_preopt:
            self.logger.info("Preoptimization disabled. Skipping geometry check.")
            return QCResult(
                success=True,
                coordinates=xyz_file,
                error_message="preopt_disabled"
            )
        
        has_overlap, min_distance = self._check_atom_overlap(xyz_file)
        
        if not has_overlap and not force_preopt:
            self.logger.info(
                f"No atom overlap detected (min_distance={min_distance:.3f} Å > threshold={self.overlap_threshold} Å). "
                "Skipping xTB preoptimization."
            )
            return QCResult(
                success=True,
                coordinates=xyz_file,
                error_message="no_overlap"
            )
        
        # Only create output_dir when we actually need to run preoptimization
        output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger.warning(
            f"Atom overlap detected (min_distance={min_distance:.3f} Å < threshold={self.overlap_threshold} Å). "
            "Triggering xTB preoptimization to avoid DFT crashes..."
        )
        
        xtb_runner = XTBRunner(self.config, work_dir=output_dir / "xtb_preopt")
        
        preopt_result = xtb_runner.optimize(
            structure=xyz_file,
            solvent=self.solvent,
            charge=charge,
            uhf=uhf
        )
        
        if not preopt_result.success:
            self.logger.error(f"xTB preoptimization failed: {preopt_result.error_message}")
            return QCResult(
                success=True,
                error_message=f"xtb_preopt_failed: {preopt_result.error_message}",
                coordinates=xyz_file
            )
        
        preopt_xyz = preopt_result.coordinates
        self.logger.info(f"xTB preoptimization successful. Output: {preopt_xyz}")
        
        has_overlap_after, min_distance_after = self._check_atom_overlap(preopt_xyz)
        if has_overlap_after:
            self.logger.warning(
                f"Warning: Overlap still exists after xTB preopt (min_distance={min_distance_after:.3f} Å). "
                "Proceeding with caution..."
            )
        else:
            self.logger.info(
                f"Overlap resolved after xTB preopt (min_distance={min_distance_after:.3f} Å)."
            )
        
        return QCResult(
            success=True,
            coordinates=preopt_xyz,
            energy=preopt_result.energy,
            converged=True,
            error_message="preopt_success"
        )
    
    def _check_atom_overlap(self, xyz_file: Path) -> Tuple[bool, float]:
        """检测 XYZ 文件中是否存在原子重叠"""
        try:
            coords, elements = self._parse_xyz(xyz_file)
            
            n_atoms = len(coords)
            min_distance = float('inf')
            
            for i in range(n_atoms):
                for j in range(i + 1, n_atoms):
                    dist = np.linalg.norm(coords[i] - coords[j])
                    min_distance = min(min_distance, dist)
            
            has_overlap = min_distance < self.overlap_threshold
            
            return has_overlap, min_distance
        
        except Exception as e:
            self.logger.error(f"Failed to check atom overlap: {e}")
            return False, float('inf')
    
    def _parse_xyz(self, xyz_file: Path) -> Tuple[np.ndarray, list]:
        """解析 XYZ 文件获取坐标和元素"""
        lines = xyz_file.read_text().strip().splitlines()
        
        n_atoms = int(lines[0])
        
        coords = []
        elements = []
        
        for line in lines[2:2 + n_atoms]:
            parts = line.split()
            if len(parts) < 4:
                continue
            
            element = parts[0]
            x, y, z = map(float, parts[1:4])
            
            elements.append(element)
            coords.append([x, y, z])
        
        return np.array(coords), elements
