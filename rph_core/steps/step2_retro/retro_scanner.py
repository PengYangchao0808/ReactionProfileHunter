import json
import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from rph_core.utils.file_io import read_xyz, write_xyz
from rph_core.utils.geometry_tools import GeometryUtils, LogParser
from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.naming_compat import (
    INTERMEDIATE_XYZ,
    REACTANT_COMPLEX_XYZ,
    create_intermediate_alias,
)
from rph_core.utils.qc_interface import XTBInterface
from rph_core.utils.gau_xtb_interface import GauXTBOptimizer
from rph_core.utils.scan_profile_plotter import (
    find_ts_and_dipole_guess,
    compute_scan_distances,
)
from rph_core.utils.tsv_dataset import ReactionRecord
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


@dataclass
class RetroScanResultV2:
    ts_guess_xyz: Optional[Path]
    reactant_xyz: Optional[Path]
    forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]]
    neutral_precursor_xyz: Optional[Path]
    meta_json_path: Optional[Path]


class RetroScanner(LoggerMixin):
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
        self.logger.info("[S2] RetroScanner initialized (direct product scan)")

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
            product_xyz / "product" / "product_global_min.xyz",
            product_xyz / "product_global_min.xyz",
            product_xyz / "product" / "global_min.xyz",
            product_xyz / "global_min.xyz",
        ]
        resolved = next((p for p in candidates if p.exists()), None)
        if resolved is None:
            raise RuntimeError(
                f"Cannot resolve product structure from {product_xyz}; tried: {[str(p) for p in candidates]}"
            )
        return resolved

    def _validate_forming_bonds(self, forming_bonds: Sequence[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
        normalized: List[Tuple[int, int]] = []
        for pair in forming_bonds:
            if not isinstance(pair, (tuple, list)) or len(pair) != 2:
                continue
            i, j = int(pair[0]), int(pair[1])
            if i == j:
                continue
            normalized.append((i, j))

        if not normalized:
            raise RuntimeError("S2 requires non-empty forming_bonds; got empty/invalid input")
        return tuple(normalized)

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
            raise RuntimeError(f"S2 inward scan requires start_distance > end_distance, got start={start}, end={end}")
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
        nproc = int(self.config.get("resources", {}).get("nproc", 1))

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

    def _map_bonds(self, forming_bonds: Sequence[Tuple[int, int]], atom_map: Optional[Dict[int, int]]) -> Tuple[Tuple[int, int], ...]:
        if not atom_map:
            return self._validate_forming_bonds(forming_bonds)
        mapped = []
        for pair in forming_bonds:
            if int(pair[0]) in atom_map and int(pair[1]) in atom_map:
                mapped.append((atom_map[int(pair[0])], atom_map[int(pair[1])]))
            else:
                self.logger.warning(f"Could not map MapId bond {pair} using atom_map. Falling back to using it as MolIdx.")
                mapped.append((int(pair[0]), int(pair[1])))
        return self._validate_forming_bonds(mapped)

    def run_retro_scan(
        self,
        product_xyz: Path,
        output_dir: Path,
        forming_bonds: Sequence[Tuple[int, int]],
        scan_config: Optional[Dict[str, Any]] = None,
        atom_map: Optional[Dict[int, int]] = None,
    ) -> Tuple[Path, Path, Path, Tuple[Tuple[int, int], ...], Path, str, str, Tuple[str, ...], Optional[Path]]:
        """向外扫描 (V5.1 Product-Seeded Relaxed Scan)"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        bonds = self._map_bonds(forming_bonds, atom_map)
        product_file = self._resolve_product_file(product_xyz)
        params = self._resolve_scan_params(scan_config)

        self.logger.info("[S2] Retro Scan (Product -> Reactant) started")
        xtb_settings = self.step2_cfg.get("xtb_settings", {}) or {}
        charge = int(xtb_settings.get("charge", 0))
        spin = int(xtb_settings.get("multiplicity", 1))

        # Direct outward scan from Product
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

        if scan_result.ts_guess_xyz is None or not Path(scan_result.ts_guess_xyz).exists():
            raise RuntimeError("S2 retro scan did not provide ts_guess geometry")

        scan_start = params["scan_start_distance"]
        scan_end = params["scan_end_distance"]
        scan_steps = params["scan_steps"]
        distances = compute_scan_distances(scan_start, scan_end, scan_steps, direction="outward")
        
        ts_distance, dipole_distance = find_ts_and_dipole_guess(
            distances, energies, energies_in_hartree=True
        )
        
        ts_guess_idx = 0
        dipole_idx = 0
        if ts_distance is not None and len(distances) > 0:
            ts_guess_idx = min(range(len(distances)), key=lambda i: abs(distances[i] - ts_distance))
            self.logger.info(f"[S2] Knee point algorithm: TS guess at index {ts_guess_idx}, distance {distances[ts_guess_idx]:.3f} Å")
        
        if dipole_distance is not None and len(distances) > 0:
            dipole_idx = min(range(len(distances)), key=lambda i: abs(distances[i] - dipole_distance))
            self.logger.info(f"[S2] Knee point algorithm: Dipole at index {dipole_idx}, distance {distances[dipole_idx]:.3f} Å")

        if isinstance(scan_result.geometries, list) and len(scan_result.geometries) > ts_guess_idx:
            ts_guess_geom = scan_result.geometries[ts_guess_idx]
            if ts_guess_geom and Path(ts_guess_geom).exists():
                self.logger.info(f"[S2] Using knee point TS guess: {ts_guess_geom}")
            else:
                ts_guess_geom = scan_result.ts_guess_xyz
        else:
            ts_guess_geom = scan_result.ts_guess_xyz
        
        if isinstance(scan_result.geometries, list) and len(scan_result.geometries) > dipole_idx:
            dipole_geom = scan_result.geometries[dipole_idx]
            if dipole_geom and Path(dipole_geom).exists():
                self.logger.info(f"[S2] Using knee point dipole: {dipole_geom}")
            else:
                dipole_geom = ts_guess_geom
        else:
            dipole_geom = ts_guess_geom

        intermediate_xyz = output_dir / INTERMEDIATE_XYZ
        shutil.copy2(Path(dipole_geom), intermediate_xyz)
        create_intermediate_alias(Path(dipole_geom), output_dir)
        self.logger.info(f"[S2] Intermediate structure: {intermediate_xyz}")
        self.logger.info(f"[S2] Backward compatibility alias: {REACTANT_COMPLEX_XYZ}")
        
        # Backward compatibility: keep both variable names for existing code
        reactant_complex_xyz = intermediate_xyz

        degraded_reasons: List[str] = []
        status = "COMPLETE"
        ts_guess_confidence = "high"

        trajectory_check: Dict[str, Any] = {"checked": False, "off_path_count": 0}
        guard_cfg = self._get_topology_guard_config()
        if guard_cfg.get("enabled", True) and isinstance(scan_result.geometries, list):
            coords, symbols, _ = LogParser.extract_last_converged_coords(product_file, engine_type="auto")
            if coords is None or symbols is None:
                coords, symbols = read_xyz(product_file)
            assert symbols is not None, "symbols should not be None after fallback"
            trajectory_check = check_scan_trajectory(
                product_coords=coords,
                symbols=symbols,
                forming_bonds=list(bonds),
                frame_paths=[g for g in scan_result.geometries if isinstance(g, Path)],
                graph_scale=guard_cfg.get("graph_scale", 1.25),
            )
            if trajectory_check.get("off_path_count", 0) > 0:
                degraded_reasons.append("scan_topology_drift_detected")
                status = "DEGRADED"
                ts_guess_confidence = "low"

        scan_profile_json = output_dir / "scan_profile.json"
        with open(scan_profile_json, "w") as f:
            json.dump({
                "generation_method": "retro_scan_from_product",
                "product_xyz": str(product_file),
                "intermediate_xyz": str(intermediate_xyz),
                "reactant_complex_xyz": str(reactant_complex_xyz),
                "forming_bonds": [list(pair) for pair in bonds],
                "scan_parameters": params,
                "knee_point_algorithm": {
                    "ts_distance": ts_distance,
                    "dipole_distance": dipole_distance,
                    "ts_geometry_index": ts_guess_idx,
                    "dipole_geometry_index": dipole_idx,
                },
                "scan_quality": {
                    "max_energy_index": max_idx,
                    "max_energy": energies[max_idx] if energies else 0,
                    "boundary_maximum": boundary_max,
                    "local_peak_ok": peak_ok,
                    "status": status,
                    "ts_guess_confidence": ts_guess_confidence,
                    "degraded_reasons": degraded_reasons,
                    "trajectory_check": trajectory_check,
                },
                "energies_hartree": energies,
            }, f, indent=2)

        ts_guess_xyz_final = output_dir / "ts_guess.xyz"
        shutil.copy2(Path(ts_guess_geom), ts_guess_xyz_final)
        shutil.copy2(Path(ts_guess_geom), output_dir / "ts_guess_s2.1.xyz")

        # Gau_XTB TS optimization for S2.1
        ts_guess_s2_1_gau_xtb = None
        gau_xtb_energy = None
        if self.step2_cfg.get("gau_xtb", {}).get("enabled", False):
            self.logger.info("[S2.1] Running Gau_XTB TS optimization")
            ts_guess_s2_1_gau_xtb, gau_xtb_energy = self._run_gau_xtb_optimization(
                ts_guess_xyz=ts_guess_xyz_final,
                output_dir=output_dir / "gau_xtb_s2.1",
                task_name="TS-S2.1-001",
                forming_bonds=list(bonds)
            )
        
        return (
            ts_guess_xyz_final,
            reactant_complex_xyz,
            reactant_complex_xyz,
            bonds,
            scan_profile_json,
            status,
            ts_guess_confidence,
            tuple(degraded_reasons),
            ts_guess_s2_1_gau_xtb
        )

    def run_with_precursor(
        self,
        reactant_complex_xyz: Path,
        record: ReactionRecord,
        output_dir: Path,
        enabled: Optional[bool] = None,
        strategy: str = "reactant_complex",
        output_meta: bool = False,
    ) -> RetroScanResultV2:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        check_enabled = enabled
        if check_enabled is None:
            check_enabled = self.config.get("step2", {}).get("neutral_precursor", {}).get("enabled", False)

        if not check_enabled:
            self.logger.info("[S2] Neutral precursor generation disabled (enabled=False)")
            return RetroScanResultV2(
                ts_guess_xyz=None,
                reactant_xyz=None,
                forming_bonds=None,
                neutral_precursor_xyz=None,
                meta_json_path=None,
            )

        self.logger.info(f"[S2] Neutral precursor generation enabled with strategy: {strategy}")
        if strategy != "reactant_complex":
            self.logger.warning(f"[S2] Unsupported strategy: {strategy}, defaulting to reactant_complex")
            strategy = "reactant_complex"

        reactant_complex_xyz = Path(reactant_complex_xyz)
        if not reactant_complex_xyz.exists():
            self.logger.error(f"[S2] Reactant complex file not found: {reactant_complex_xyz}")
            return RetroScanResultV2(
                ts_guess_xyz=None,
                reactant_xyz=None,
                forming_bonds=None,
                neutral_precursor_xyz=None,
                meta_json_path=None,
            )

        neutral_precursor_xyz = output_dir / "neutral_precursor.xyz"
        try:
            shutil.copy(reactant_complex_xyz, neutral_precursor_xyz)
            self.logger.info(f"[S2] Neutral precursor created: {neutral_precursor_xyz}")
        except Exception as exc:
            self.logger.error(f"[S2] Failed to copy reactant complex to neutral precursor: {exc}")
            return RetroScanResultV2(
                ts_guess_xyz=None,
                reactant_xyz=None,
                forming_bonds=None,
                neutral_precursor_xyz=None,
                meta_json_path=None,
            )

        meta_json_path = None
        if output_meta:
            meta_json_path = output_dir / "meta.json"
            meta_data = {
                "precursor_smiles": record.precursor_smiles,
                "leaving_small_molecule_key": record.get_leaving_small_molecule_key(),
                "strategy": strategy,
                "source_reactant_complex": str(reactant_complex_xyz),
            }
            try:
                meta_json_path.write_text(json.dumps(meta_data, indent=2))
                self.logger.info(f"[S2] Meta JSON written: {meta_json_path}")
            except Exception as exc:
                self.logger.warning(f"[S2] Failed to write meta.json: {exc}")
                meta_json_path = None

        return RetroScanResultV2(
            ts_guess_xyz=None,
            reactant_xyz=None,
            forming_bonds=None,
            neutral_precursor_xyz=neutral_precursor_xyz,
            meta_json_path=meta_json_path,
        )

    def _generate_constraints(self, bond_1: Tuple[int, int], bond_2: Tuple[int, int], target_dist: float) -> str:
        return f"""$constrain
  force constant=0.5
  distance: {bond_1[0]+1}, {bond_1[1]+1}, {target_dist:.3f}
  distance: {bond_2[0]+1}, {bond_2[1]+1}, {target_dist:.3f}
$end
"""

    def run_path_search(
        self,
        start_xyz: Path,
        end_xyz: Path,
        output_dir: Path,
        forming_bonds: Optional[Sequence[Tuple[int, int]]] = None,
    ) -> Tuple[Path, Path, Path, Tuple[Tuple[int, int], ...], Path, str, str, Tuple[str, ...], Optional[Path]]:
        from rph_core.utils.qc_interface import XTBInterface
        from rph_core.utils.data_types import PathSearchResult

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        path_cfg = self.step2_cfg.get("path_search", {}) or {}
        if not path_cfg.get("enabled", False):
            raise ValueError("path_search is not enabled in config")

        self.logger.info("[S2] Path Search (xTB --path) started")

        xtb_settings = self.step2_cfg.get("xtb_settings", {}) or {}
        charge = int(xtb_settings.get("charge", 0))
        spin = int(xtb_settings.get("multiplicity", 1))

        path_dir = output_dir / "path_search"
        path_dir.mkdir(parents=True, exist_ok=True)

        xtb = XTBInterface(
            solvent=xtb_settings.get("solvent", "acetone"),
            nproc=int(self.config.get("resources", {}).get("nproc", 1)),
            config=self.config
        )

        result: PathSearchResult = xtb.path(
            start_xyz=start_xyz,
            end_xyz=end_xyz,
            output_dir=path_dir,
            nrun=path_cfg.get("nrun", 1),
            npoint=path_cfg.get("npoint", 25),
            anopt=path_cfg.get("anopt", 10),
            kpush=path_cfg.get("kpush", 0.003),
            kpull=path_cfg.get("kpull", -0.015),
            ppull=path_cfg.get("ppull", 0.05),
            alp=path_cfg.get("alp", 1.2),
            charge=charge,
            spin=spin,
        )

        if not result.success:
            raise RuntimeError(f"S2 path search failed: {result.error_message}")

        if result.ts_guess_xyz is None or not result.ts_guess_xyz.exists():
            raise RuntimeError("S2 path search did not produce ts_guess")

        ts_guess_xyz = output_dir / "ts_guess.xyz"
        shutil.copy2(result.ts_guess_xyz, ts_guess_xyz)

        reactant_complex_xyz = output_dir / "reactant_complex.xyz"
        if Path(start_xyz).resolve() != Path(reactant_complex_xyz).resolve():
            shutil.copy2(start_xyz, reactant_complex_xyz)

        degraded_reasons: List[str] = []
        status = "COMPLETE"
        ts_guess_confidence = "high"

        if result.gradient_norm_at_ts is not None and result.gradient_norm_at_ts > 0.05:
            degraded_reasons.append("high_gradient_norm_at_ts")
            ts_guess_confidence = "medium"

        path_profile_json = output_dir / "scan_profile.json"
        with open(path_profile_json, "w") as f:
            json.dump({
                "generation_method": "xtb_path_search",
                "start_xyz": str(start_xyz),
                "end_xyz": str(end_xyz),
                "reactant_complex_xyz": str(reactant_complex_xyz),
                "ts_guess_xyz": str(ts_guess_xyz),
                "forming_bonds": [list(pair) for pair in (forming_bonds or [])],
                "path_parameters": {
                    "nrun": path_cfg.get("nrun", 1),
                    "npoint": path_cfg.get("npoint", 25),
                    "anopt": path_cfg.get("anopt", 10),
                    "kpush": path_cfg.get("kpush", 0.003),
                    "kpull": path_cfg.get("kpull", -0.015),
                },
                "energies": {
                    "barrier_forward_kcal": result.barrier_forward_kcal,
                    "barrier_backward_kcal": result.barrier_backward_kcal,
                    "reaction_energy_kcal": result.reaction_energy_kcal,
                },
                "ts_quality": {
                    "estimated_ts_point": result.estimated_ts_point,
                    "gradient_norm_at_ts": result.gradient_norm_at_ts,
                    "status": status,
                    "ts_guess_confidence": ts_guess_confidence,
                    "degraded_reasons": degraded_reasons,
                },
            }, f, indent=2)

        bonds = tuple(forming_bonds) if forming_bonds else tuple()

        # Gau_XTB TS optimization for S2.2
        ts_guess_s2_2_gau_xtb = None
        gau_xtb2_distance = None
        if self.step2_cfg.get("gau_xtb", {}).get("enabled", False):
            self.logger.info("[S2.2] Running Gau_XTB TS optimization")
            ts_guess_s2_2_gau_xtb, gau_xtb2_energy = self._run_gau_xtb_optimization(
                ts_guess_xyz=ts_guess_xyz,
                output_dir=output_dir / "gau_xtb_s2.2",
                task_name="TS-S2.2-001",
                forming_bonds=list(forming_bonds) if forming_bonds else []
            )
            
            if ts_guess_s2_2_gau_xtb and gau_xtb2_energy is not None and forming_bonds:
                coords, symbols = read_xyz(ts_guess_s2_2_gau_xtb)
                bond_pair = forming_bonds[0]
                atom_i, atom_j = bond_pair  # type: ignore
                dist = float(np.linalg.norm(np.array(coords[atom_i]) - np.array(coords[atom_j])))
                gau_xtb2_distance = dist

        return (
            ts_guess_xyz,
            reactant_complex_xyz,
            reactant_complex_xyz,
            bonds,
            path_profile_json,
            status,
            ts_guess_confidence,
            tuple(degraded_reasons),
            ts_guess_s2_2_gau_xtb
        )

    def _run_gau_xtb_optimization(
        self,
        ts_guess_xyz: Path,
        output_dir: Path,
        task_name: str,
        forming_bonds: List[Tuple[int, int]]
    ) -> Tuple[Optional[Path], Optional[float]]:
        """Run Gau_XTB TS optimization on a TS guess.
        
        Returns:
            Tuple of (optimized_xyz_path, energy_hartree)
        """
        gau_xtb_cfg = self.step2_cfg.get("gau_xtb", {})
        
        nproc = self.config.get("resources", {}).get("nproc", 1)
        
        optimizer = GauXTBOptimizer(
            config=self.config,
            nproc=nproc
        )
        
        try:
            opt_xyz, confidence, imag_freq = optimizer.optimize(
                ts_guess_xyz=ts_guess_xyz,
                output_dir=output_dir,
                task_name=task_name,
                forming_bonds=forming_bonds,
                max_attempts=gau_xtb_cfg.get("max_attempts", 3)
            )
            
            output_name = f"ts_guess_{task_name.lower().replace('-', '_')}_gau_xtb.xyz"
            final_xyz = output_dir / output_name
            shutil.copy2(opt_xyz, final_xyz)
            
            energy_hartree = None
            if optimizer.interface is not None:
                log_path = opt_xyz.parent / "input.log"
                if log_path.exists():
                    energy_hartree = optimizer.interface._parse_energy(log_path)
            
            self.logger.info(
                f"[S2] Gau_XTB TS optimization complete: {final_xyz} "
                f"(confidence: {confidence}, imag_freq: {imag_freq})"
            )
            
            return final_xyz, energy_hartree
            
        except Exception as e:
            self.logger.error(f"[S2] Gau_XTB TS optimization failed: {e}")
            return None, None
