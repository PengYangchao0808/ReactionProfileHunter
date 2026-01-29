"""
统一优化配置模块 (Unified Optimization Config)
================================================

提供引擎无关的几何优化配置，支持 ORCA 和 Gaussian 的关键词自动转换

Author: QC Descriptors Team
Date: 2026-01-13
Purpose: ReactionProfileHunter v2.1 - 统一优化参数配置
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple
from enum import Enum
import logging
import re

from rph_core.utils.keyword_translator import KeywordTranslator

logger = logging.getLogger(__name__)


@dataclass
class OptimizationConfig:
    """
    统一优化配置（引擎无关）

    支持自动转换为 Gaussian route card 或 ORCA %geom block
    """
    timeout_enabled: bool = False
    timeout_seconds: Optional[int] = None

    oscillation_window_size: int = 10
    oscillation_energy_tolerance: float = 1e-4
    oscillation_max_count: int = 3

    initial_hessian: Optional[str] = None
    recalc_hess_every: int = 10

    max_step: int = 30
    trust_radius: float = 0.3
    adaptive_step: bool = True

    convergence: str = "normal"

    ts_eigentest: bool = False
    ts_follow_mode: int = 1
    ts_step_type: str = "eigenvalue"

    @classmethod
    def from_config(cls, config_dict: Dict[str, Any]) -> 'OptimizationConfig':
        """从配置字典创建 OptimizationConfig"""
        opt_ctrl = config_dict.get('optimization_control', {})

        timeout_cfg = opt_ctrl.get('timeout', {})
        osc_cfg = opt_ctrl.get('oscillation', {})
        hess_cfg = opt_ctrl.get('hessian', {})
        step_cfg = opt_ctrl.get('step', {})
        conv_cfg = opt_ctrl.get('convergence', {})
        ts_cfg = opt_ctrl.get('ts', {})

        return cls(
            timeout_enabled=timeout_cfg.get('enabled', False),
            timeout_seconds=timeout_cfg.get('default_seconds'),
            oscillation_window_size=osc_cfg.get('window_size', 10),
            oscillation_energy_tolerance=osc_cfg.get('energy_tolerance', 1e-4),
            oscillation_max_count=osc_cfg.get('max_oscillation_count', 3),
            initial_hessian=hess_cfg.get('initial', None),
            recalc_hess_every=hess_cfg.get('recalc_every', 10),
            max_step=step_cfg.get('max_step', 30),
            trust_radius=step_cfg.get('trust_radius', 0.3),
            adaptive_step=step_cfg.get('adaptive', True),
            convergence=conv_cfg.get('level', 'normal'),
            ts_eigentest=ts_cfg.get('eigentest', False),
            ts_follow_mode=ts_cfg.get('follow_mode', 1),
            ts_step_type=ts_cfg.get('step_type', 'eigenvalue')
        )

    def to_gaussian_route(
        self,
        method: str,
        basis: str,
        dispersion: Optional[str] = None,
        is_ts: bool = False
    ) -> str:
        """
        生成 Gaussian route card (增强兼容性版)

        Args:
            method: DFT 方法（如 B3LYP）
            basis: 基组（如 def2-SVP）
            dispersion: 色散校正（如 GD3BJ）
            is_ts: 是否为 TS 优化

        Returns:
            Gaussian route card 字符串
        """
        # 1. 基础 Route 头: 使用 #p 输出更详细信息
        # 防御性编程: def2-SVP -> def2SVP (防止 G09/G16 解析横杠出错)
        clean_basis = basis.replace("-", "") if basis.lower().startswith("def2") else basis
        route_parts = [f"#p {method}/{clean_basis}"]

        # 2. 色散校正 (使用简写 em=gd3bj，更安全)
        if dispersion:
            # 用户建议: em=gd3bj 比 EmpiricalDispersion=GD3BJ 更简洁且不易报错
            route_parts.append(f"em={dispersion}")

        # 3. 优化选项构建
        opt_options = []
        
        # TS 特有选项
        if is_ts:
            opt_options.append("TS")
            if self.initial_hessian is None:
                pass
            elif self.initial_hessian == "calcfc":
                opt_options.append("CalcFC")
            elif self.initial_hessian == "calcall":
                opt_options.append("CalcAll")
            elif self.initial_hessian == "read":
                opt_options.append("ReadFC")

            if not self.ts_eigentest:
                opt_options.append("NoEigenTest")

            # RecalcFC 在 TS 中很有用
            if self.recalc_hess_every > 0 and self.initial_hessian != "calcall":
                opt_options.append(f"RecalcFC={self.recalc_hess_every}")

        # 通用选项
        # MaxStep 对于震荡控制很重要
        if self.max_step != 30:
            opt_options.append(f"MaxStep={self.max_step}")
            
        # 4. 组装 Opt 关键词
        if opt_options:
            # 例如: Opt=(TS,CalcFC,MaxStep=10)
            route_parts.append(f"Opt=({','.join(opt_options)})")
        else:
            route_parts.append("Opt") # 默认优化

        # 5. 收敛标准 (合并进 Opt 或是独立写均可，为防止由括号解析问题，这里使用追加逻辑)
        if self.convergence != "normal":
            # 将 Tight/Loose 注入到 Opt 选项中是最稳妥的做法
            last_opt = route_parts.pop() # 取出 "Opt" 或 "Opt=(...)"
            
            if "(" in last_opt:
                # 已有括号: Opt=(TS,CalcFC) -> Opt=(TS,CalcFC,Tight)
                new_opt = last_opt[:-1] + f",{self.convergence.capitalize()})"
            else:
                # 无括号: Opt -> Opt=(Tight)
                new_opt = f"Opt=({self.convergence.capitalize()})"
            
            route_parts.append(new_opt)

        # 6. SCF 选项 (XQC 是收敛的最后一道防线)
        route_parts.append("SCF=(XQC,MaxCycle=256)")

        route_str = " ".join(route_parts)
        logger.debug(f"生成的 Gaussian route: {route_str}")

        return route_str

    def to_orca_geom_block(self, is_ts: bool = False) -> str:
        """
        生成 ORCA %geom 块

        Args:
            is_ts: 是否为 TS 优化

        Returns:
            ORCA %geom 块字符串
        """
        lines = ["%geom"]

        if is_ts:
            lines.append("   Constraints")
            lines.append("      {B} 2 1 2")  # 使用内坐标约束
            lines.append("   end")

        if self.initial_hessian is None:
            lines.append("   Calc_Hess false")
        elif self.initial_hessian == "calcfc":
            lines.append("   Calc_Hess true")
        elif self.initial_hessian == "calcall":
            lines.append("   Calc_Hess true")
            lines.append(f"   Recalc_Hess 1")
        else:
            lines.append("   Calc_Hess false")

        if self.recalc_hess_every > 1 and self.initial_hessian != "calcall":
            lines.append(f"   Recalc_Hess {self.recalc_hess_every}")

        if self.max_step != 30:
            lines.append(f"   MaxStep {self.max_step / 100}")  # 转换为 Bohr

        if self.trust_radius != 0.3:
            lines.append(f"   Trust {self.trust_radius}")

        if self.adaptive_step:
            lines.append("   Constraints true")

        if is_ts:
            lines.append(f"   TSMode {self.ts_follow_mode}")
            if self.ts_step_type == "eigenvalue":
                lines.append("   StepType Eigenvalue")

        lines.append("end")

        geom_block = "\n".join(lines)
        logger.debug(f"生成的 ORCA %geom 块:\n{geom_block}")

        return geom_block

    def update_for_rescue(self, rescue_params: Dict[str, Any]):
        """根据救援参数更新配置"""
        if 'recalc_hess_every' in rescue_params:
            self.recalc_hess_every = rescue_params['recalc_hess_every']
        if 'max_step' in rescue_params:
            self.max_step = rescue_params['max_step']
        if 'trust_radius' in rescue_params:
            self.trust_radius = rescue_params['trust_radius']
        if 'initial_hessian' in rescue_params:
            self.initial_hessian = rescue_params['initial_hessian']

    def to_dict(self) -> Dict[str, Any]:
        """转换为字典格式"""
        return {
            'timeout': {'enabled': self.timeout_enabled, 'default_seconds': self.timeout_seconds},
            'oscillation': {'window_size': self.oscillation_window_size, 'energy_tolerance': self.oscillation_energy_tolerance, 'max_oscillation_count': self.oscillation_max_count},
            'hessian': {'initial': self.initial_hessian, 'recalc_every': self.recalc_hess_every},
            'step': {'max_step': self.max_step, 'trust_radius': self.trust_radius, 'adaptive': self.adaptive_step},
            'convergence': {'level': self.convergence},
            'ts': {'eigentest': self.ts_eigentest, 'follow_mode': self.ts_follow_mode, 'step_type': self.ts_step_type}
        }


def build_gaussian_route_from_config(
    config: Dict[str, Any],
    rescue: bool = False,
    include_freq: bool = True
) -> str:
    theory_opt = config.get('theory', {}).get('optimization', {})
    method = theory_opt.get('method', 'B3LYP')
    basis_raw = theory_opt.get('basis', 'def2SVP')
    basis = KeywordTranslator.to_gaussian_basis(basis_raw)

    dispersion = KeywordTranslator.to_gaussian_dispersion(theory_opt.get('dispersion'))
    solvent = KeywordTranslator.to_gaussian_solvent(theory_opt.get('solvent'))

    route_parts = [f"#p {method}/{basis}"]
    if dispersion:
        route_parts.append(dispersion)
    if solvent:
        route_parts.append(solvent)

    if rescue:
        route_parts.append("Opt=(CalcFC,NoEigenTest)")
    else:
        route_parts.append("Opt=CalcFC")
    if include_freq:
        route_parts.append("Freq")

    return " ".join(route_parts)


def _normalize_def2_in_route(route: str) -> str:
    pattern = re.compile(r"def2[-_]?([A-Za-z0-9]+)", re.IGNORECASE)
    return pattern.sub(lambda match: f"def2{match.group(1)}", route)


def _normalize_noeigentest_in_route(route: str) -> str:
    if "NoEigenTest" not in route or "Opt=" not in route:
        return route
    if "Opt=(" in route:
        return route

    tokens = [token for token in route.split() if token != "NoEigenTest"]
    route_no_token = " ".join(tokens)

    def _inject_noeigentest(match: re.Match) -> str:
        opt_value = match.group(1)
        return f"Opt=({opt_value},NoEigenTest)"

    normalized = re.sub(r"Opt=([^\s]+)", _inject_noeigentest, route_no_token, count=1)
    return normalized


def normalize_qc_config(
    config: Dict[str, Any],
    auto_fix: bool = True,
    strict: bool = False
) -> Tuple[Dict[str, Any], List[Dict[str, str]]]:
    fixes: List[Dict[str, str]] = []

    theory = config.setdefault('theory', {})
    optimization = theory.setdefault('optimization', {})
    engine = optimization.get('engine', 'gaussian')

    if str(engine).lower() == 'gaussian':
        route = optimization.get('route')
        if not route:
            new_route = build_gaussian_route_from_config(config, rescue=False)
            optimization['route'] = new_route
            fixes.append({
                'field': 'theory.optimization.route',
                'original': '<missing>',
                'updated': new_route,
                'reason': 'default_from_config'
            })
        else:
            normalized = _normalize_def2_in_route(route)
            normalized = _normalize_noeigentest_in_route(normalized)
            if normalized != route:
                if strict and not auto_fix:
                    raise ValueError(f"Invalid Gaussian route: {route}")
                if auto_fix:
                    optimization['route'] = normalized
                fixes.append({
                    'field': 'theory.optimization.route',
                    'original': route,
                    'updated': normalized,
                    'reason': 'gaussian_route_normalization'
                })

        rescue_route = optimization.get('rescue_route')
        if not rescue_route:
            new_route = build_gaussian_route_from_config(config, rescue=True)
            optimization['rescue_route'] = new_route
            fixes.append({
                'field': 'theory.optimization.rescue_route',
                'original': '<missing>',
                'updated': new_route,
                'reason': 'default_from_config'
            })
        else:
            normalized = _normalize_def2_in_route(rescue_route)
            normalized = _normalize_noeigentest_in_route(normalized)
            if normalized != rescue_route:
                if strict and not auto_fix:
                    raise ValueError(f"Invalid Gaussian route: {rescue_route}")
                if auto_fix:
                    optimization['rescue_route'] = normalized
                fixes.append({
                    'field': 'theory.optimization.rescue_route',
                    'original': rescue_route,
                    'updated': normalized,
                    'reason': 'gaussian_route_normalization'
                })

    step3 = config.get('step3', {})
    gaussian_keywords = step3.get('gaussian_keywords', {})
    for key, value in gaussian_keywords.items():
        if isinstance(value, str) and value:
            normalized = _normalize_def2_in_route(value)
            normalized = _normalize_noeigentest_in_route(normalized)
            if normalized != value:
                if strict and not auto_fix:
                    raise ValueError(f"Invalid Gaussian keyword '{key}': {value}")
                if auto_fix:
                    gaussian_keywords[key] = normalized
                fixes.append({
                    'field': f"step3.gaussian_keywords.{key}",
                    'original': value,
                    'updated': normalized,
                    'reason': 'gaussian_route_normalization'
                })

    return config, fixes
