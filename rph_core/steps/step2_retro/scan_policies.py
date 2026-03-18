from enum import Enum
from dataclasses import dataclass
from typing import Dict, Tuple, List, Sequence

class ScanPolicy(Enum):
    """V5.1 Scan Policies"""
    CONCERTED_FIXED = "concerted"           # Legacy behavior: strong dual constraints
    POLICY_A_SINGLE_FREE = "policy_a"       # Single-side freedom (Constrain B, Free A)
    POLICY_B_MIRROR = "policy_b"            # Mirror single-side (Constrain A, Free B)
    POLICY_C_WEAK_DUAL = "policy_c"         # Weak dual constraint (k=0.01-0.05)

@dataclass
class ConstraintConfig:
    """Constraint Configuration"""
    constrained_bonds: List[int]  # Bond indices to constrain (e.g., [0], [1], or [0,1])
    constraint_type: str         # "fixed" | "harmonic"
    force_constant: float        # Force constant (fixed=0.5, harmonic=0.03)

class ScanPolicySelector:
    """Policy Selector for Step 2 Scanning"""
    
    POLICY_CONFIGS = {
        ScanPolicy.CONCERTED_FIXED: ConstraintConfig(
            constrained_bonds=[0, 1],
            constraint_type="fixed",
            force_constant=0.5
        ),
        ScanPolicy.POLICY_A_SINGLE_FREE: ConstraintConfig(
            constrained_bonds=[1],  # Only constrain Bond B
            constraint_type="fixed",
            force_constant=0.5
        ),
        ScanPolicy.POLICY_B_MIRROR: ConstraintConfig(
            constrained_bonds=[0],  # Only constrain Bond A
            constraint_type="fixed", 
            force_constant=0.5
        ),
        ScanPolicy.POLICY_C_WEAK_DUAL: ConstraintConfig(
            constrained_bonds=[0, 1],
            constraint_type="harmonic",
            force_constant=0.03  # Very weak
        ),
    }
    
    def select_policy(self, bonds: Sequence[Tuple[int, int]], policy_name: str, target_distance: float) -> Tuple[Dict[str, float], float]:
        """
        Produce xTB constraints for the specified scan policy.
        
        Args:
            bonds: list of (atom1_idx, atom2_idx) for forming bonds.
            policy_name: string name of the ScanPolicy.
            target_distance: the scan target distance to set for the constrained bond(s).
            
        Returns:
            (constraints_dict, force_constant)
        """
        try:
            policy = ScanPolicy(policy_name)
        except ValueError:
            policy = ScanPolicy.CONCERTED_FIXED

        config = self.POLICY_CONFIGS[policy]
        constraints = {}
        for bond_idx in config.constrained_bonds:
            if bond_idx < len(bonds):
                bond = bonds[bond_idx]
                # Note: bonds holds the atom indices (0-indexed). xTB constraint usually uses 1-indexed.
                # However, XTBInterface normalizes constraints by adding 1. So we pass 0-indexed string.
                constraints[f"{bond[0]} {bond[1]}"] = target_distance
                
        return constraints, config.force_constant
