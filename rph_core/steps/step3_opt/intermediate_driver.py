"""
S3 Intermediate DFT Optimization Driver
=======================================

Performs DFT-level optimization on the intermediate structure from S2.

Author: QCcalc Team
Date: 2026-03-16
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any

from rph_core.utils.log_manager import LoggerMixin

logger = logging.getLogger(__name__)


class IntermediateDriver(LoggerMixin):
    """
    Intermediate DFT optimization driver.
    
    Wraps qc_task_runner.run_opt_sp_cycle() for parameter encapsulation.
    Reuses global theory.optimization and theory.single_point configurations.
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize intermediate driver.
        
        Args:
            config: Full configuration dictionary
        """
        self.config = config
        self.cfg = config.get('step3', {}).get('intermediate_opt', {})
        self.enabled = self.cfg.get('enabled', True)
        
        self.theory_opt = config.get('theory', {}).get('optimization', {})
        self.theory_sp = config.get('theory', {}).get('single_point', {})
        
        self.charge = self.cfg.get('charge', 0)
        self.multiplicity = self.cfg.get('multiplicity', 1)
        
        self._qc_runner = None
        
        logger.info(f"[IntermediateDriver] Initialized: enabled={self.enabled}")
        logger.info(f"  Method: {self.theory_opt.get('method')}/{self.theory_opt.get('basis')}")
        logger.info(f"  SP: {self.theory_sp.get('method')}/{self.theory_sp.get('basis')}")
    
    @property
    def qc_runner(self):
        """延迟初始化 QC Runner"""
        if self._qc_runner is None:
            from rph_core.utils.qc_task_runner import QCTaskRunner
            self._qc_runner = QCTaskRunner(config=self.config)
        return self._qc_runner
    
    def run(
        self,
        intermediate_xyz: Path,
        output_dir: Path,
    ) -> Dict[str, Any]:
        """
        对中间体执行 DFT 优化 + 频率 + L2 SP
        
        完整流程 (复用 run_opt_sp_cycle):
        1. xTB 预优化 (qc_task_runner 内置)
        2. DFT 几何优化 + 频率
        3. L2 单点能量
        
        Args:
            intermediate_xyz: S2 生成的 intermediate.xyz (反应中间体)
            output_dir: 输出目录
        
        Returns:
            dict: {
                'optimized_xyz': Path - 优化后的结构
                'l2_energy': float - L2 单点能量 (Hartree)
                'converged': bool - 是否收敛
                'opt_output': Path - 优化输出文件
                'sp_output': Path - SP 输出文件
                'error': str - 错误信息 (如有)
            }
        """
        if not self.enabled:
            logger.info("[IntermediateDriver] Disabled, skipping")
            return {'converged': False, 'skipped': True}
        
        intermediate_xyz = Path(intermediate_xyz)
        if not intermediate_xyz.exists():
            logger.error(f"[IntermediateDriver] Input not found: {intermediate_xyz}")
            return {'converged': False, 'error': f'Input not found: {intermediate_xyz}'}
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("=" * 60)
        logger.info("[IntermediateDriver] Starting DFT optimization for intermediate")
        logger.info(f"  Input: {intermediate_xyz}")
        logger.info(f"  Output: {output_dir}")
        logger.info("=" * 60)
        
        try:
            # 复用 run_opt_sp_cycle
            result = self.qc_runner.run_opt_sp_cycle(
                xyz_file=intermediate_xyz,
                output_dir=output_dir,
                charge=self.charge,
                spin=self.multiplicity,
                enable_l2_sp=True,
                enable_nbo=False,
            )
            
            # 解析结果
            if result.converged:
                logger.info(f"[IntermediateDriver] Optimization converged")
                logger.info(f"  L2 Energy: {result.l2_energy:.6f} Hartree")
            else:
                logger.warning(f"[IntermediateDriver] Optimization failed: {result.error_message}")
            
            return {
                'optimized_xyz': result.optimized_xyz,
                'l2_energy': result.l2_energy,
                'converged': result.converged,
                'opt_output': result.log_file,
                'sp_output': result.l2_sp_result.output_file if result.l2_sp_result else None,
                'freq_output': result.freq_log,
                'error': result.error_message,
            }
            
        except Exception as e:
            logger.exception(f"[IntermediateDriver] Unexpected error: {e}")
            return {'converged': False, 'error': str(e)}
