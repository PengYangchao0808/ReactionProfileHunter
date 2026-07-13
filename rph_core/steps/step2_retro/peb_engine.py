import json
import logging
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from rph_core.utils.file_io import read_xyz, write_xyz
from rph_core.utils.geometry_tools import GeometryUtils, LogParser
from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.naming_compat import (
    INTERMEDIATE_XYZ,
)
from rph_core.utils.qc_interface import XTBInterface
from rph_core.utils.scan_profile_plotter import (
    find_ts_and_dipole_guess,
    compute_scan_distances,
)
from rph_core.utils.ui import get_progress_manager

from .bond_stretcher import BondStretcher
from .geometry_guard import (
    compare_graph_topology,
    check_scan_trajectory,
    detect_risky_contacts,
    generate_keepaway_constraints,
    TopologyGuardResult,
    RiskyContactResult,
)

logger = logging.getLogger(__name__)


class PEBScanEngine(LoggerMixin):
    DEFAULT_SCAN_START_DISTANCE = 3.5
    DEFAULT_SCAN_END_DISTANCE = 1.8
    DEFAULT_SCAN_STEPS = 20
    DEFAULT_SCAN_MODE = "concerted"
    DEFAULT_SCAN_FORCE_CONSTANT = 0.5
    DEFAULT_MIN_VALID_POINTS = 5
    DEFAULT_INTERMEDIATE_MIN_RMSD = 0.15

    def __init__(self, config: Dict[str, Any], molecule_name: Optional[str] = None):
        self.config = config
        self.step2_cfg = config.get("step2", {}) if isinstance(config, dict) else {}
        self.molecule_name = molecule_name
        self.bond_stretcher = BondStretcher()
        self._seed_guard_result: Optional[Dict[str, Any]] = None
        self.logger.info("[S2] PEB scan engine initialized")

    def _optimize_intermediate(
        self,
        seed: Path,
        out_dir: Path,
        forming_bonds: Tuple[Tuple[int, int], ...],
        scan_start_distance: float,
    ) -> Path:
        return Path(seed)

    def _update_ui_status(self, output_dir: Path, status_text: str) -> None:
        pm = get_progress_manager()
        pm.update_step("s2", description=status_text)

        status_file = output_dir / ".rph_step_status.json"
        try:
            with open(status_file, "w") as f:
                json.dump({"step": "s2", "description": status_text}, f)
        except Exception as exc:
            self.logger.warning(f"[S2] Failed to write status file: {exc}")

    def _resolve_product_file(self, product_xyz: Path) -> Path:
        product_xyz = Path(product_xyz)
        if not product_xyz.exists():
            raise FileNotFoundError(f"S2 product input not found: {product_xyz}")
        if product_xyz.is_file():
            return product_xyz

        candidates = [
            product_xyz / "product_min.xyz",
        ]
        resolved = next((p for p in candidates if p.exists()), None)
        if resolved is None:
            raise RuntimeError(
                f"Cannot resolve product structure from {product_xyz}; tried: {[str(p) for p in candidates]}"
            )
        return resolved

    def _validate_forming_bonds(
        self,
        forming_bonds: Sequence[Tuple[int, int]],
        product_xyz_path: Optional[Path] = None,
    ) -> Tuple[Tuple[int, int], ...]:
        """Validate forming bonds, with optional distance check."""
        normalized: List[Tuple[int, int]] = []
        for pair in forming_bonds:
            if not isinstance(pair, (tuple, list)) or len(pair) != 2:
                continue
            i, j = int(pair[0]), int(pair[1])
            if i == j:
                continue
            normalized.append((min(i, j), max(i, j)))

        if not normalized:
            raise RuntimeError("S2 requires non-empty forming_bonds; got empty/invalid input")

        # Distance-based sanity check
        if product_xyz_path is not None and Path(product_xyz_path).exists():
            try:
                coords, _ = read_xyz(Path(product_xyz_path))
                suspicious = []
                for pair in normalized:
                    i, j = pair
                    if i >= len(coords) or j >= len(coords):
                        self.logger.warning(
                            f"[S2] Forming bond {pair} has atom index >= total atoms "
                            f"({len(coords)}) in product XYZ — index space mismatch"
                        )
                        suspicious.append(pair)
                        continue
                    dist = GeometryUtils.calculate_distance(coords, i, j)
                    if dist > 5.0 or dist < 0.5:
                        self.logger.warning(
                            f"[S2] Forming bond {pair} distance={dist:.2f} Å looks invalid "
                            f"(expected 1.5-5.0 Å for forming bonds) — possible index space mismatch"
                        )
                        suspicious.append(pair)
                if suspicious:
                    self.logger.warning(
                        f"[S2] {len(suspicious)}/{len(normalized)} forming bonds "
                        f"failed distance validation: {suspicious}"
                    )
            except Exception as exc:
                self.logger.debug(f"[S2] Distance validation unavailable: {exc}")

        return tuple(sorted(set(normalized)))

    def _get_topology_guard_config(self) -> Dict[str, Any]:
        """Get topology guard configuration with defaults."""
        scan_cfg = dict(self.step2_cfg.get("scan", {}) or {})
        return {
            "enabled": scan_cfg.get("topology_guard_enabled", True),
            "graph_scale": scan_cfg.get("topology_graph_scale", 1.25),
            "near_bond_ratio": scan_cfg.get("risk_contact_ratio_threshold", 0.85),
            "near_bond_max": scan_cfg.get("risk_contact_abs_cutoff_A", 2.2),
            "min_shrink_ratio": scan_cfg.get("min_shrink_ratio", 0.75),
            "max_risky_pairs": scan_cfg.get("keep_apart_max_pairs", 6),
            "keep_apart_floor": scan_cfg.get("keep_apart_floor_A", 3.0),
            "constraint_force": scan_cfg.get("keep_apart_force_constant", 0.5),
            "retry_force": scan_cfg.get("keep_apart_retry_force_constant", 1.0),
            "retry_once": scan_cfg.get("topology_retry_once", True),
        }

    def _resolve_scan_params(self, scan_config: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        scan_cfg = dict((self.step2_cfg.get("scan", {}) or {}))
        if scan_config:
            scan_cfg.update(scan_config)

        start = float(scan_cfg.get("scan_start_distance", scan_cfg.get("start_distance", self.DEFAULT_SCAN_START_DISTANCE)))
        end = float(scan_cfg.get("scan_end_distance", scan_cfg.get("end_distance", self.DEFAULT_SCAN_END_DISTANCE)))
        steps = int(scan_cfg.get("scan_steps", scan_cfg.get("steps", self.DEFAULT_SCAN_STEPS)))
        mode = str(scan_cfg.get("scan_mode", self.DEFAULT_SCAN_MODE))
        force_constant = float(scan_cfg.get("scan_force_constant", self.DEFAULT_SCAN_FORCE_CONSTANT))
        min_valid_points = int(scan_cfg.get("min_valid_points", self.DEFAULT_MIN_VALID_POINTS))
        reject_boundary_maximum = bool(scan_cfg.get("reject_boundary_maximum", True))
        require_local_peak = bool(scan_cfg.get("require_local_peak", False))
        boundary_retry_once = bool(scan_cfg.get("boundary_retry_once", True))
        boundary_retry_delta = float(scan_cfg.get("boundary_retry_delta", 0.3))
        boundary_retry_extra_steps = int(scan_cfg.get("boundary_retry_extra_steps", 6))
        allow_boundary_degradation = bool(scan_cfg.get("allow_boundary_degradation", True))

        if start <= end:
            raise RuntimeError(f"S2 scan requires scan_start_distance > scan_end_distance, got scan_start={start}, scan_end={end}")
        if steps <= 1:
            raise RuntimeError(f"S2 scan_steps must be > 1, got {steps}")

        return {
            "scan_start_distance": start,
            "scan_end_distance": end,
            "scan_steps": steps,
            "scan_mode": mode,
            "scan_force_constant": force_constant,
            "min_valid_points": min_valid_points,
            "reject_boundary_maximum": reject_boundary_maximum,
            "require_local_peak": require_local_peak,
            "boundary_retry_once": boundary_retry_once,
            "boundary_retry_delta": boundary_retry_delta,
            "boundary_retry_extra_steps": boundary_retry_extra_steps,
            "allow_boundary_degradation": allow_boundary_degradation,
            "scan_policy": scan_cfg.get("scan_policy", "policy_c"),
        }

    def _execute_scan(
        self,
        start_xyz: Path,
        output_dir: Path,
        bonds: Tuple[Tuple[int, int], ...],
        params: Dict[str, Any],
        direction: str,
        charge: int = 0,
        spin: int = 1,
    ) -> Tuple[Any, List[float], int, bool, bool]:
        """核心扫描执行 (V5.1) - with ScanPolicySelector"""
        scan_policy_name = params.get("scan_policy", "policy_c")
        from rph_core.steps.step2_retro.scan_policies import ScanPolicySelector
        selector = ScanPolicySelector()

        if direction == "inward":
            start_dist = params["scan_start_distance"]
            end_dist = params["scan_end_distance"]
        elif direction == "outward":
            start_dist = params["scan_end_distance"]
            end_dist = params["scan_start_distance"]
        else:
            raise ValueError(f"Unknown direction: {direction}")

        scan_constraints, force_constant = selector.select_policy(bonds, scan_policy_name, start_dist)

        xtb_settings = self.step2_cfg.get("xtb_settings", {}) or {}
        solvent = str(xtb_settings.get("solvent", self.config.get("theory", {}).get("optimization", {}).get("solvent", "acetone")))
        scan_cfg = self.step2_cfg.get("scan", {}) or {}
        nproc = int(scan_cfg.get("nproc") or xtb_settings.get("nproc") or self.config.get("resources", {}).get("nproc", 1))

        xtb = XTBInterface(solvent=solvent, nproc=nproc, config=self.config)
        result = xtb.scan(
            xyz_file=start_xyz,
            output_dir=output_dir,
            constraints=scan_constraints,
            scan_range=(start_dist, end_dist),
            scan_steps=params["scan_steps"],
            scan_mode=params["scan_mode"],
            scan_force_constant=force_constant,
            charge=charge,
            spin=spin,
        )

        if not result.success or not result.energies:
            raise RuntimeError(f"S2 {direction} scan failed or returned no energies")

        energies_local = [float(e) for e in result.energies]
        if len(energies_local) < params["min_valid_points"]:
            raise RuntimeError(f"S2 {direction} scan returned too few valid points: {len(energies_local)} < {params['min_valid_points']}")

        max_idx_local = max(range(len(energies_local)), key=energies_local.__getitem__)
        boundary_local = max_idx_local in {0, len(energies_local) - 1}
        local_peak_local = False
        if 0 < max_idx_local < len(energies_local) - 1:
            local_peak_local = (
                energies_local[max_idx_local] > energies_local[max_idx_local - 1]
                and energies_local[max_idx_local] > energies_local[max_idx_local + 1]
            )
        return result, energies_local, max_idx_local, boundary_local, local_peak_local

    def run(
        self,
        product_xyz: Path,
        output_dir: Path,
        forming_bonds: Sequence[Tuple[int, int]],
        scan_config: Optional[Dict[str, Any]] = None,
    ) -> Tuple[Path, Path, Path, Tuple[Tuple[int, int], ...], Path, str, str, Tuple[str, ...]]:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        bonds = self._validate_forming_bonds(forming_bonds, product_xyz_path=product_xyz)
        product_file = self._resolve_product_file(product_xyz)
        params = self._resolve_scan_params(scan_config)

        xtb_settings = self.step2_cfg.get("xtb_settings", {}) or {}
        charge = int(xtb_settings.get("charge", 0))
        spin = int(xtb_settings.get("multiplicity", 1))

        coords, symbols, _ = LogParser.extract_last_converged_coords(product_file, engine_type="auto")
        if coords is None or symbols is None:
            coords, symbols = read_xyz(product_file)
        assert symbols is not None

        for pair in bonds:
            dist = GeometryUtils.calculate_distance(coords, int(pair[0]), int(pair[1]))
            if dist > 3.0:
                self.logger.warning(
                    f"[S2] possible index mapping error: Forming bond {pair} distance is {dist:.2f} Å in product geometry"
                )

        seed_targets: List[Tuple[Tuple[int, int], float]] = [
            ((int(i), int(j)), float(params["scan_start_distance"])) for (i, j) in bonds
        ]
        seed_coords = self.bond_stretcher.stretch_bonds(coords, seed_targets)
        seed_xyz = output_dir / "intermediate_seed.xyz"
        write_xyz(seed_xyz, seed_coords, symbols, title="intermediate_seed")
        intermediate_seed = self._optimize_intermediate(
            seed_xyz,
            output_dir,
            bonds,
            float(params["scan_start_distance"]),
        )

        try:
            inter_coords, inter_symbols = read_xyz(intermediate_seed)
            min_dist = float("inf")
            for i in range(len(inter_symbols)):
                for j in range(i + 1, len(inter_symbols)):
                    d = GeometryUtils.calculate_distance(inter_coords, i, j)
                    if d < min_dist:
                        min_dist = d
            if min_dist < 0.5:
                self.logger.warning(
                    f"[S2] Intermediate geometry warning: unusually short interatomic distance {min_dist:.2f} Å"
                )
        except Exception as exc:
            self.logger.warning(f"[S2] Intermediate geometry warning: failed to validate intermediate ({exc})")

        degraded_reasons: List[str] = []
        status = "COMPLETE"
        ts_guess_confidence = "high"

        scan_dir = output_dir / "retro_scan"
        scan_result, energies, max_idx, boundary_max, peak_ok = self._execute_scan(
            start_xyz=product_file,
            output_dir=scan_dir,
            bonds=bonds,
            params=params,
            direction="outward",
            charge=charge,
            spin=spin,
        )

        reject_boundary = bool(params.get("reject_boundary_maximum", True))
        retry_once = bool(params.get("boundary_retry_once", True))
        allow_degradation = bool(params.get("allow_boundary_degradation", True))

        if reject_boundary and boundary_max and retry_once:
            params_retry = dict(params)
            params_retry["scan_start_distance"] = float(params["scan_start_distance"]) + float(params.get("boundary_retry_delta", 0.3))
            params_retry["scan_steps"] = int(params["scan_steps"]) + int(params.get("boundary_retry_extra_steps", 6))
            retry_dir = output_dir / "retro_scan_retry"
            scan_result, energies, max_idx, boundary_max, peak_ok = self._execute_scan(
                start_xyz=product_file,
                output_dir=retry_dir,
                bonds=bonds,
                params=params_retry,
                direction="outward",
                charge=charge,
                spin=spin,
            )

            if boundary_max:
                if allow_degradation:
                    status = "DEGRADED"
                    ts_guess_confidence = "low"
                    degraded_reasons.append("boundary_maximum_persisted_after_retry")
                else:
                    raise RuntimeError("S2 boundary maximum persisted after retry")
            params = params_retry
        elif reject_boundary and boundary_max:
            if allow_degradation:
                status = "DEGRADED"
                ts_guess_confidence = "low"
                degraded_reasons.append("boundary_maximum_detected")
            else:
                raise RuntimeError("S2 boundary maximum detected")

        if scan_result is None or getattr(scan_result, "ts_guess_xyz", None) is None:
            raise RuntimeError("S2 scan did not provide ts_guess geometry")

        ts_guess_xyz_final = output_dir / "ts_guess.xyz"
        shutil.copy2(Path(scan_result.ts_guess_xyz), ts_guess_xyz_final)

        dipolar_xyz = output_dir / INTERMEDIATE_XYZ
        shutil.copy2(Path(intermediate_seed), dipolar_xyz)

        reactant_xyz = output_dir / "reactant_complex.xyz"
        shutil.copy2(dipolar_xyz, reactant_xyz)

        scan_profile_json = output_dir / "scan_profile.json"
        with open(scan_profile_json, "w") as f:
            json.dump(
                {
                    "generation_method": "retro_scan",
                    "product_xyz": str(product_file),
                    "intermediate_xyz": str(dipolar_xyz),
                    "forming_bonds": [list(pair) for pair in bonds],
                    "scan_parameters": params,
                    "scan_quality": {
                        "max_energy_index": int(max_idx),
                        "boundary_maximum": bool(boundary_max),
                        "local_peak_ok": bool(peak_ok),
                        "status": status,
                        "ts_guess_confidence": ts_guess_confidence,
                        "degraded_reasons": degraded_reasons,
                    },
                "energies_hartree": energies,
                },
                f,
                indent=2,
            )

        return (
            ts_guess_xyz_final,
            reactant_xyz,
            dipolar_xyz,
            bonds,
            scan_profile_json,
            status,
            ts_guess_confidence,
            tuple(degraded_reasons),
        )

