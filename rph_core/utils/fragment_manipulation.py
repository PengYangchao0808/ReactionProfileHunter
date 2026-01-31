"""
Fragment Manipulation Utilities
==============================

H-cap and fragment manipulation operations for fragmenter.

Author: QCcalc Team
Date: 2026-01-27
Updated: 2026-01-31 - Fixed electron counting for H-capped fragments
"""

from typing import List, Tuple, Dict, Optional
import logging

import numpy as np

logger = logging.getLogger(__name__)

# 原子序数表（用于电子计数）
ATOMIC_NUMBERS: Dict[str, int] = {
    'H': 1, 'He': 2,
    'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
    'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
    'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32,
    'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39,
    'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46,
    'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54
}


def count_electrons(symbols: List[str], charge: int) -> int:
    """
    计算分子的电子总数
    
    Args:
        symbols: 原子符号列表
        charge: 分子电荷
    
    Returns:
        电子总数 = Σ(原子序数) - 电荷
    """
    total_z = sum(ATOMIC_NUMBERS.get(s, 0) for s in symbols)
    return total_z - charge


def is_closed_shell(symbols: List[str], charge: int, multiplicity: int = 1) -> bool:
    """
    检查是否为闭壳层体系
    
    闭壳层要求：电子数为偶数，多重度为 1
    
    Args:
        symbols: 原子符号列表
        charge: 分子电荷
        multiplicity: 自旋多重度
    
    Returns:
        True 如果是闭壳层
    """
    n_electrons = count_electrons(symbols, charge)
    return n_electrons % 2 == 0 and multiplicity == 1


def infer_charge_for_closed_shell(symbols: List[str], preferred_charge: int = 0) -> int:
    """
    推断使体系为闭壳层的电荷
    
    逻辑：
    - 如果 Σ(Z) 为偶数，preferred_charge 若为偶数则可用，否则 ±1
    - 如果 Σ(Z) 为奇数，preferred_charge 若为奇数则可用，否则 ±1
    
    Args:
        symbols: 原子符号列表
        preferred_charge: 优先使用的电荷
    
    Returns:
        使电子数为偶数的电荷值
    """
    total_z = sum(ATOMIC_NUMBERS.get(s, 0) for s in symbols)
    
    # 电子数 = total_z - charge
    # 要求电子数为偶数
    # 即 (total_z - charge) % 2 == 0
    # 等价于 total_z % 2 == charge % 2
    
    if total_z % 2 == preferred_charge % 2:
        # preferred_charge 已经满足闭壳层条件
        return preferred_charge
    else:
        # 需要调整电荷
        # 优先选择最接近 preferred_charge 且满足条件的值
        candidates = [preferred_charge - 1, preferred_charge + 1]
        for c in candidates:
            if total_z % 2 == c % 2:
                return c
        # 理论上不会到达这里
        return preferred_charge


def h_cap_fragment(
    coords: np.ndarray,
    symbols: List[str],
    cap_positions: List[Tuple[int, np.ndarray]],
    cap_bond_lengths: Optional[Dict[str, float]] = None
) -> Tuple[np.ndarray, List[str]]:
    """
    Add hydrogen atoms to cap open valences.

    Args:
        coords: (N, 3) coordinates
        symbols: List of element symbols
        cap_positions: List of (atom_idx, cap_direction_vec)
            - atom_idx: Index of atom to cap
            - cap_direction_vec: Vector pointing from atom toward cap position
              (should be colinear with bond being cut, pointing away from cut bond)
        cap_bond_lengths: Bond length by element (Å)
            Default: {'C': 1.09, 'N': 1.01, 'O': 0.96, 'H': 0.76}

    Returns:
        (capped_coords, capped_symbols) - Extended coordinates and symbols
    """
    if cap_bond_lengths is None:
        cap_bond_lengths = {
            'C': 1.09, 'N': 1.01, 'O': 0.96, 'H': 0.76,
            'Si': 1.48, 'P': 1.42, 'S': 1.35, 'Cl': 1.27,
            'Br': 1.42, 'I': 1.54
        }

    capped_coords = coords.copy()
    capped_symbols = symbols.copy()

    for atom_idx, direction_vec in cap_positions:
        norm = np.linalg.norm(direction_vec)
        if norm < 1e-8:
            raise ValueError("Cap direction vector is near zero")
        direction = direction_vec / norm

        element = symbols[atom_idx]
        bond_length = cap_bond_lengths.get(element, 1.09)

        cap_pos = coords[atom_idx] + direction * bond_length

        capped_coords = np.vstack([capped_coords, cap_pos])
        capped_symbols.append('H')

    return capped_coords, capped_symbols


def get_fragment_charges(
    total_charge: int,
    n_fragA: int,
    n_fragB: int,
    dipole_in_fragA: bool = True
) -> Tuple[int, int]:
    """
    Distribute total charge between two fragments (legacy interface).

    For oxidopyrylium [5+2], assign +1 to dipole (fragment A).
    
    NOTE: This function does NOT consider H-capping. For H-capped fragments,
    use get_hcapped_fragment_charge() instead.

    Args:
        total_charge: Total charge of the system
        n_fragA: Number of atoms in fragment A
        n_fragB: Number of atoms in fragment B
        dipole_in_fragA: If True, assign formal +1 charge to fragment A

    Returns:
        (charge_fragA, charge_fragB) tuple
    """
    formal_dipole_charge = 1
    if dipole_in_fragA:
        charge_fragA = formal_dipole_charge
        charge_fragB = total_charge - charge_fragA
    else:
        charge_fragB = formal_dipole_charge
        charge_fragA = total_charge - charge_fragB

    return charge_fragA, charge_fragB


def get_hcapped_fragment_charge(
    symbols: List[str],
    n_h_caps: int,
    preferred_charge: int = 0,
    force_closed_shell: bool = True
) -> Tuple[int, int]:
    """
    为 H-capped 片段计算正确的电荷和多重度
    
    核心逻辑：
    1. 计算片段的总核电荷 Σ(Z)（包括 H-cap 原子）
    2. 为保证闭壳层（电子数为偶数），调整电荷使 (Σ(Z) - charge) 为偶数
    
    Args:
        symbols: H-capped 后的原子符号列表（已包含 H-cap 原子）
        n_h_caps: 添加的 H-cap 原子数量（用于日志，实际计算基于 symbols）
        preferred_charge: 优先使用的电荷（如果满足闭壳层条件）
        force_closed_shell: 是否强制闭壳层（默认 True）
    
    Returns:
        (charge, multiplicity) 元组
        - charge: 使电子数为偶数的电荷
        - multiplicity: 1 (闭壳层) 或 2 (如果无法闭壳层)
    """
    total_z = sum(ATOMIC_NUMBERS.get(s, 0) for s in symbols)
    
    if force_closed_shell:
        # 计算使电子数为偶数的电荷
        charge = infer_charge_for_closed_shell(symbols, preferred_charge)
        multiplicity = 1
        
        n_electrons = total_z - charge
        logger.debug(
            f"H-capped fragment: {len(symbols)} atoms, Σ(Z)={total_z}, "
            f"charge={charge}, electrons={n_electrons}, mult={multiplicity}"
        )
        
        if charge != preferred_charge:
            logger.warning(
                f"Adjusted H-capped fragment charge from {preferred_charge} to {charge} "
                f"to maintain closed-shell (even electrons)"
            )
    else:
        # 不强制闭壳层，使用 preferred_charge
        charge = preferred_charge
        n_electrons = total_z - charge
        
        # 根据电子奇偶性确定多重度
        if n_electrons % 2 == 0:
            multiplicity = 1  # 闭壳层
        else:
            multiplicity = 2  # 双重态（最低合理多重度）
            logger.warning(
                f"H-capped fragment has odd electrons ({n_electrons}), "
                f"using multiplicity={multiplicity}"
            )
    
    return charge, multiplicity


def get_fragment_multiplicities(
    total_multiplicity: int,
    n_fragA: int,
    n_fragB: int,
    dipole_in_fragA: bool = True
) -> Tuple[int, int]:
    """
    Distribute total multiplicity between two fragments (legacy interface).

    Default: assign multiplicity to fragment containing dipole.
    
    NOTE: For H-capped fragments, use get_hcapped_fragment_charge() instead,
    which correctly computes multiplicity based on electron count.

    Args:
        total_multiplicity: Total multiplicity of the system
        n_fragA: Number of atoms in fragment A
        n_fragB: Number of atoms in fragment B
        dipole_in_fragA: If True, assign multiplicity to fragment A

    Returns:
        (mult_fragA, mult_fragB) tuple
    """
    if dipole_in_fragA:
        return total_multiplicity, 1
    else:
        return 1, total_multiplicity


def compute_vrm_model_charge(
    combined_symbols: List[str],
    system_charge: int,
    n_h_caps: int = 2
) -> Tuple[int, int]:
    """
    计算 VRM 合并模型的正确电荷和多重度
    
    VRM 模型 = FragA(H-capped) + FragB(H-capped)
    原子数 = 原始 TS 原子数 + n_h_caps
    
    Args:
        combined_symbols: 合并后的原子符号列表（包含所有 H-cap）
        system_charge: 原始 TS 体系的电荷
        n_h_caps: H-cap 原子数量（默认 2，因为断一根键两端各加一个 H）
    
    Returns:
        (charge, multiplicity) 元组
    """
    return get_hcapped_fragment_charge(
        symbols=combined_symbols,
        n_h_caps=n_h_caps,
        preferred_charge=system_charge,
        force_closed_shell=True
    )
