import hashlib
import hashlib
import json
import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from rph_core.utils.file_io import write_xyz
from rph_core.utils.geometry_tools import LogParser

logger = logging.getLogger(__name__)


@dataclass
class InspectionResult:
    should_skip: bool
    reason: str
    files_found: List[Path]


class ResultInspector:
    def __init__(self, work_dir: Path, config: dict, context: Optional[Dict[str, str]] = None, strict_mode: bool = True):
        self.work_dir = Path(work_dir)
        self.config = config
        self.context = context or {}
        self.strict = strict_mode
        self.force_rerun = self._read_flag("FORCE_RERUN")
        self.ignore_sig = self._read_flag("IGNORE_SIG")

    def check_step(self, step_name: str) -> InspectionResult:
        if self.force_rerun:
            return InspectionResult(False, "force_rerun_enabled", [])

        current_sig = self._generate_signature(step_name)
        sig_file = self.work_dir / f"{step_name}_sig.json"
        if not self.ignore_sig and sig_file.exists():
            try:
                saved_sig = json.loads(sig_file.read_text()).get("signature")
            except json.JSONDecodeError:
                return InspectionResult(False, f"signature_corrupted:{step_name}", [])
            if saved_sig != current_sig:
                return InspectionResult(False, f"signature_mismatch:{step_name}", [])

        checker = getattr(self, f"_check_{step_name.lower()}", None)
        if not checker:
            raise ValueError(f"Unknown step: {step_name}")

        result = checker()
        if result.should_skip and not self.ignore_sig:
            if not sig_file.exists():
                self.save_signature(step_name)
        return result

    def save_signature(self, step_name: str) -> None:
        sig_file = self.work_dir / f"{step_name}_sig.json"
        payload = {
            "signature": self._generate_signature(step_name),
            "payload": self._signature_payload(step_name)
        }
        sig_file.write_text(json.dumps(payload, sort_keys=True))

    def _read_flag(self, name: str) -> bool:
        env_value = os.getenv(name)
        if env_value is not None:
            return env_value.lower() in {"1", "true", "yes"}
        return bool(self.config.get(name, False))

    def _signature_payload(self, step_name: str) -> Dict[str, Optional[str]]:
        theory_opt = self.config.get("theory", {}).get("optimization", {})
        theory_sp = self.config.get("theory", {}).get("single_point", {})
        payload = {
            "smiles": self.context.get("smiles"),
            "charge": self.context.get("charge", 0),
            "multiplicity": self.context.get("multiplicity", 1),
            "opt_method": theory_opt.get("method"),
            "opt_basis": theory_opt.get("basis"),
            "opt_dispersion": theory_opt.get("dispersion"),
            "opt_solvent": theory_opt.get("solvent"),
            "opt_engine": theory_opt.get("engine"),
            "sp_method": theory_sp.get("method"),
            "sp_basis": theory_sp.get("basis"),
            "sp_solvent": theory_sp.get("solvent"),
            "sp_engine": theory_sp.get("engine"),
            "step": step_name
        }
        return payload

    def _generate_signature(self, step_name: str) -> str:
        payload = self._signature_payload(step_name)
        payload_str = json.dumps(payload, sort_keys=True)
        return hashlib.sha256(payload_str.encode("utf-8")).hexdigest()

    def _is_fake_path(self, path: Path) -> bool:
        lower = str(path).lower()
        return "fake" in lower or "debug" in lower

    def _check_normal_termination(self, log_file: Path, engine: str) -> bool:
        if not log_file.exists():
            return False
        try:
            with open(log_file, "rb") as handle:
                handle.seek(0, 2)
                size = handle.tell()
                handle.seek(max(size - 4096, 0), 0)
                tail = handle.read().decode("utf-8", errors="ignore")
        except Exception:
            return False

        engine_lower = engine.lower()
        if engine_lower == "orca":
            return "ORCA TERMINATED NORMALLY" in tail
        if engine_lower == "gaussian":
            return "Normal termination" in tail
        return False

    def _get_orca_energy(self, out_file: Path) -> float:
        try:
            content = out_file.read_text(errors="ignore")
        except Exception:
            return float("inf")
        match = re.search(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", content)
        if match:
            return float(match.group(1))
        return float("inf")

    def _select_global_min_xyz(self) -> Optional[Path]:
        candidates = [
            self.work_dir / "S1_ConfGeneration" / "product" / "product_global_min.xyz",
            self.work_dir / "S1_ConfGeneration" / "product_global_min.xyz",
            self.work_dir / "S1_Product" / "product" / "product_global_min.xyz",
            self.work_dir / "S1_Product" / "product_global_min.xyz"
        ]
        for candidate in candidates:
            if candidate.exists() and candidate.stat().st_size > 0:
                return candidate
        return None

    def _global_min_target(self) -> Path:
        candidates = [
            self.work_dir / "S1_ConfGeneration" / "product" / "product_global_min.xyz",
            self.work_dir / "S1_ConfGeneration" / "product_global_min.xyz",
            self.work_dir / "S1_Product" / "product" / "product_global_min.xyz",
            self.work_dir / "S1_Product" / "product_global_min.xyz"
        ]
        for candidate in candidates:
            if candidate.parent.exists():
                return candidate
        return candidates[0]

    def _select_sp_dir(self) -> Optional[Path]:
        candidates = [
            self.work_dir / "S1_ConfGeneration" / "product" / "dft",
            self.work_dir / "S1_ConfGeneration" / "dft",
            self.work_dir / "S1_Product" / "product" / "dft",
            self.work_dir / "S1_Product" / "dft"
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return None

    def _check_s1(self) -> InspectionResult:
        global_min = self._select_global_min_xyz()

        sp_dir = self._select_sp_dir()
        if not sp_dir:
            files = [global_min] if global_min else []
            return InspectionResult(False, "missing_sp_dir", files)

        valid_sps: List[Path] = []
        energies: Dict[Path, float] = {}
        for out_file in sp_dir.glob("*.out"):
            if self._is_fake_path(out_file):
                continue
            if self._check_normal_termination(out_file, engine="orca"):
                energy = self._get_orca_energy(out_file)
                valid_sps.append(out_file)
                energies[out_file] = energy

        if not valid_sps:
            files = [global_min] if global_min else []
            return InspectionResult(False, "no_valid_sp", files)

        best_sp = sorted(valid_sps, key=lambda path: energies.get(path, float("inf")))[0]
        if not global_min:
            coords, symbols, _ = LogParser.extract_last_converged_coords(best_sp, engine_type="auto")
            if coords is None or symbols is None:
                return InspectionResult(False, "global_min_missing", [best_sp])
            global_min = self._global_min_target()
            global_min.parent.mkdir(parents=True, exist_ok=True)
            write_xyz(global_min, coords, symbols, title=f"Global Min SP E={energies[best_sp]:.6f}")
        return InspectionResult(True, "s1_complete", [path for path in [global_min, best_sp] if path])

    def _check_s2(self) -> InspectionResult:
        s2_dir = self.work_dir / "S2_ReactionPath"
        reactant = s2_dir / "reactant_complex.xyz"
        ts_guess = s2_dir / "ts_guess.xyz"
        if reactant.exists() and reactant.stat().st_size > 0 and ts_guess.exists() and ts_guess.stat().st_size > 0:
            return InspectionResult(True, "s2_complete", [reactant, ts_guess])
        return InspectionResult(False, "missing_s2_outputs", [])

    def _check_s3(self) -> InspectionResult:
        s3_dir = self.work_dir / "S3_TS"
        ts_xyz = s3_dir / "ts_final.xyz"
        if not ts_xyz.exists() or ts_xyz.stat().st_size == 0:
            return InspectionResult(False, "missing_ts_xyz", [])

        if self.strict:
            logs = [
                path for path in s3_dir.iterdir()
                if path.is_file() and path.suffix in (".log", ".out")
            ]
            logs = [log for log in logs if not self._is_fake_path(log) and "ts" in log.name.lower()]
            for log in logs:
                if self._check_normal_termination(log, engine="gaussian") or self._check_normal_termination(log, engine="orca"):
                    return InspectionResult(True, "s3_complete", [ts_xyz])
            return InspectionResult(False, "ts_not_terminated", [ts_xyz])

        return InspectionResult(True, "s3_complete_loose", [ts_xyz])

    def _check_s4(self) -> InspectionResult:
        """
        V6.1-A2: Check Step4 completion (using is_step4_complete logic)

        使用 checkpoint_manager.is_step4_complete() 考虑 mechanism packaging 完整性
        """
        from rph_core.utils.checkpoint_manager import CheckpointManager

        # V6.1: Use is_step4_complete to check Step4 completion
        # 而不是只检查 features.csv
        checkpoint_mgr = CheckpointManager(self.work_dir)
        s4_dir = self.work_dir / "S4_Data"

        if s4_dir.exists():
            s4_complete = checkpoint_mgr.is_step4_complete(s4_dir, self.config)
            if s4_complete:
                return InspectionResult(True, "s4_complete", [])

        # Fallback to原始逻辑（如果 is_step4_complete 失败）
        candidates = [
            self.work_dir / "features_raw.csv",
        ]
        for csv_path in candidates:
            if csv_path.exists() and csv_path.stat().st_size > 100:
                return InspectionResult(True, "s4_complete", [csv_path])
        return InspectionResult(False, "missing_features", [])
