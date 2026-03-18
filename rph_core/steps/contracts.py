from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from rph_core.utils.data_types import MapId, MolIdx


@dataclass
class Step2Input:
    product_xyz: Path
    forming_bonds: Tuple[Tuple[MapId, MapId], ...]


@dataclass
class Step2Artifacts:
    ts_guess_xyz: Path
    substrate_xyz: Path
    intermediate_xyz: Path
    forming_bonds: Tuple[Tuple[MolIdx, MolIdx], ...]
    forming_bonds_map: Optional[Tuple[Tuple[MapId, MapId], ...]] = None
    generation_method: str = "unknown"
    status: str = "COMPLETE"
    ts_guess_confidence: str = "high"
    degraded_reasons: Tuple[str, ...] = tuple()
    step2_signature: Optional[Dict[str, Any]] = None
    scan_profile_json: Optional[Path] = None


@dataclass
class Step3Artifacts:
    ts_final_xyz: Path
    sp_report: Any
    ts_fchk: Optional[Path] = None
    ts_log: Optional[Path] = None
    ts_qm_output: Optional[Path] = None
    intermediate_fchk: Optional[Path] = None
    intermediate_log: Optional[Path] = None
    intermediate_qm_output: Optional[Path] = None
    intermediate_xyz: Optional[Path] = None
    intermediate_l2_energy: Optional[float] = None
    intermediate_opt_output: Optional[Path] = None
    intermediate_sp_output: Optional[Path] = None


@dataclass
class Step4Artifacts:
    features_csv: Path
