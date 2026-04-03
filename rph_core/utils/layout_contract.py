from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional


STEP_ORDER = ("s0", "s1", "s2", "s3", "s4")

STEP_DIRS: Dict[str, str] = {
    "s0": "S0_Mechanism",
    "s1": "S1_ConfGeneration",
    "s2": "S2_Retro",
    "s3": "S3_TS",
    "s4": "S4_Data",
}


@dataclass(frozen=True)
class StepContract:
    step_id: str
    required_files: Dict[str, str]


CONTRACTS: Dict[str, StepContract] = {
    "s0": StepContract(
        "s0",
        {
            "mechanism_graph_json": "mechanism_graph.json",
            "mechanism_summary_json": "mechanism_summary.json",
        },
    ),
    "s1": StepContract("s1", {"product_xyz": "product_min.xyz"}),
    "s2": StepContract(
        "s2",
        {
            "ts_guess_xyz": "ts_guess.xyz",
            "intermediate_xyz": "intermediate.xyz",
        },
    ),
    "s3": StepContract(
        "s3",
        {
            "ts_final_xyz": "ts_final.xyz",
            "sp_matrix_metadata_json": "sp_matrix_metadata.json",
        },
    ),
    "s4": StepContract("s4", {"features_raw_csv": "features_raw.csv"}),
}


def resolve_step_dir(work_dir: Path, step_id: str) -> Path:
    if step_id not in STEP_DIRS:
        raise KeyError(f"Unknown step id: {step_id}")
    return Path(work_dir) / STEP_DIRS[step_id]


def resolve_required_files(work_dir: Path, step_id: str) -> Dict[str, Path]:
    contract = CONTRACTS[step_id]
    step_dir = resolve_step_dir(work_dir, step_id)
    return {k: step_dir / v for k, v in contract.required_files.items()}


def check_step_minimal_complete(
    work_dir: Path,
    step_id: str,
    config: Optional[Mapping[str, Any]] = None,
) -> bool:
    required = resolve_required_files(work_dir, step_id)
    for path in required.values():
        if not path.exists() or path.stat().st_size == 0:
            return False

    if step_id == "s4":
        step4_cfg = ((config or {}).get("step4", {}) or {})
        mech_cfg = (step4_cfg.get("mechanism_packaging", {}) or {})
        if bool(mech_cfg.get("enabled", False)):
            mech_index = resolve_step_dir(work_dir, "s4") / "mech_index.json"
            if not mech_index.exists() or mech_index.stat().st_size == 0:
                return False
    return True


def seed_steps_template() -> Dict[str, Dict[str, Any]]:
    seeded: Dict[str, Dict[str, Any]] = {}
    for step_id in STEP_ORDER:
        seeded[f"step_{step_id}"] = {
            "step_name": step_id,
            "completed": False,
            "timestamp": "",
            "output_files": {},
            "metadata": {},
        }
    return seeded


def canonical_output_files(work_dir: Path, step_id: str) -> Dict[str, str]:
    resolved = resolve_required_files(work_dir, step_id)
    return {k: str(v) for k, v in resolved.items() if v.exists()}


def iter_step_ids() -> Iterable[str]:
    return STEP_ORDER
