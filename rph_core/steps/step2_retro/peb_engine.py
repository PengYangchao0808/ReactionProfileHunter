import json
import logging
import math
import shutil
from pathlib import Path
from statistics import median
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from rph_core.utils.file_io import read_xyz, write_xyz
from rph_core.utils.geometry_tools import GeometryUtils, LogParser
from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.naming_compat import (
    INTERMEDIATE_XYZ,
)
from rph_core.utils.qc_interface import XTBInterface
from rph_core.utils.scan_profile_plotter import (
    HARTREE_TO_KCAL,
    compute_scan_distances,
    plot_scan_profile,
)
from rph_core.utils.ui import get_progress_manager

from .bond_stretcher import BondStretcher
from .geometry_guard import (
    check_scan_trajectory,
)

logger = logging.getLogger(__name__)


class PEBScanEngine(LoggerMixin):
    DEFAULT_SCAN_START_DISTANCE = 3.5
    DEFAULT_SCAN_END_DISTANCE = 1.8
    DEFAULT_SCAN_STEPS = 20
    DEFAULT_SCAN_MODE = "concerted"
    DEFAULT_SCAN_FORCE_CONSTANT = 0.5
    DEFAULT_MIN_VALID_POINTS = 5
    def __init__(self, config: Dict[str, Any], molecule_name: Optional[str] = None):
        self.config = config
        self.step2_cfg = config.get("step2", {}) if isinstance(config, dict) else {}
        self.molecule_name = molecule_name
        self.bond_stretcher = BondStretcher()
        self.logger.info("[S2] PEB scan engine initialized")

    def _optimize_intermediate(
        self,
        seed: Path,
        out_dir: Path,
        forming_bonds: Tuple[Tuple[int, int], ...],
        scan_start_distance: float,
    ) -> Path:
        """Retain the stretched endpoint only as a diagnostic artifact."""
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
        intermediate_cfg = dict(scan_cfg.get("intermediate_selection", {}) or {})

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
            "intermediate_selection": {
                "method": str(intermediate_cfg.get("method", "stable_midpoint")),
                "endpoint_exclusion_points": int(intermediate_cfg.get("endpoint_exclusion_points", 1)),
                "stable_quantile": float(intermediate_cfg.get("stable_quantile", 0.5)),
                "min_candidate_points": int(intermediate_cfg.get("min_candidate_points", 2)),
                "require_topology_valid": bool(intermediate_cfg.get("require_topology_valid", True)),
                "allow_degraded_fallback": bool(intermediate_cfg.get("allow_degraded_fallback", False)),
                "endpoint_jump_warning_ratio": float(intermediate_cfg.get("endpoint_jump_warning_ratio", 5.0)),
                "endpoint_jump_warning_min_kcal": float(
                    intermediate_cfg.get("endpoint_jump_warning_min_kcal", 5.0)
                ),
            },
        }

    @staticmethod
    def _quantile_threshold(values: Sequence[float], quantile: float) -> float:
        if not values:
            raise RuntimeError("Cannot compute a stability threshold from an empty sequence")
        ordered = sorted(float(value) for value in values)
        bounded = min(1.0, max(0.0, float(quantile)))
        rank = max(0, min(len(ordered) - 1, int(math.ceil(bounded * len(ordered))) - 1))
        return ordered[rank]

    @classmethod
    def _select_intermediate_frame(
        cls,
        energies: Sequence[float],
        frame_paths: Sequence[Path],
        peak_index: int,
        off_path_indices: Sequence[int],
        selection_config: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Select an existing, stable scan frame between the TS peak and endpoint."""
        if len(energies) != len(frame_paths):
            raise RuntimeError(
                "S2 cannot select intermediate: energy/frame count mismatch "
                f"({len(energies)} energies, {len(frame_paths)} frames)"
            )
        if not 0 <= int(peak_index) < len(energies):
            raise RuntimeError(f"S2 peak index is outside scan trajectory: {peak_index}")

        endpoint_index = len(energies) - 1
        endpoint_exclusion = max(1, int(selection_config.get("endpoint_exclusion_points", 1)))
        candidate_stop = endpoint_index - endpoint_exclusion
        candidate_indices = list(range(int(peak_index) + 1, candidate_stop + 1))
        minimum_candidates = max(1, int(selection_config.get("min_candidate_points", 2)))
        allow_fallback = bool(selection_config.get("allow_degraded_fallback", False))
        if len(candidate_indices) < minimum_candidates and not allow_fallback:
            raise RuntimeError(
                "S2 cannot select intermediate: only "
                f"{len(candidate_indices)} interior frame(s) after peak {peak_index}; "
                f"requires {minimum_candidates}"
            )

        off_path = {int(index) for index in off_path_indices}
        require_topology = bool(selection_config.get("require_topology_valid", True))
        valid_indices = [
            index for index in candidate_indices if not require_topology or index not in off_path
        ]
        degraded_fallback = False
        if not valid_indices:
            if not allow_fallback or not candidate_indices:
                raise RuntimeError("S2 cannot select intermediate: no topology-valid interior frames")
            valid_indices = list(candidate_indices)
            degraded_fallback = True

        roughness: Dict[int, float] = {}
        for index in valid_indices:
            left_delta = abs(float(energies[index]) - float(energies[index - 1]))
            right_delta = abs(float(energies[index + 1]) - float(energies[index]))
            roughness[index] = max(left_delta, right_delta) * HARTREE_TO_KCAL

        threshold = cls._quantile_threshold(
            list(roughness.values()),
            float(selection_config.get("stable_quantile", 0.5)),
        )
        stable_indices = [index for index in valid_indices if roughness[index] <= threshold]
        midpoint_index = (float(peak_index) + float(endpoint_index)) / 2.0
        selected_index = min(
            stable_indices,
            key=lambda index: (
                abs(float(index) - midpoint_index),
                float(energies[index]),
                roughness[index],
                index,
            ),
        )

        return {
            "index": int(selected_index),
            "frame_xyz": str(Path(frame_paths[selected_index])),
            "rule": "stable_midpoint",
            "midpoint_index": midpoint_index,
            "candidate_indices": candidate_indices,
            "stable_indices": stable_indices,
            "off_path_indices": sorted(off_path),
            "roughness_kcal_mol": roughness[selected_index],
            "stability_threshold_kcal_mol": threshold,
            "endpoint_index": endpoint_index,
            "endpoint_exclusion_points": endpoint_exclusion,
            "degraded_fallback": degraded_fallback,
        }

    @staticmethod
    def _forming_bond_distances_by_frame(
        frame_paths: Sequence[Path],
        forming_bonds: Sequence[Tuple[int, int]],
    ) -> List[Optional[List[float]]]:
        distances: List[Optional[List[float]]] = []
        for frame_path in frame_paths:
            try:
                frame_coords, _ = read_xyz(Path(frame_path))
                distances.append(
                    [
                        float(GeometryUtils.calculate_distance(frame_coords, int(i), int(j)))
                        for i, j in forming_bonds
                    ]
                )
            except Exception:
                distances.append(None)
        return distances

    @staticmethod
    def _endpoint_jump_detected(
        energies: Sequence[float],
        peak_index: int,
        selection_config: Dict[str, Any],
    ) -> bool:
        if len(energies) - int(peak_index) < 4:
            return False
        step_changes = [
            abs(float(energies[index]) - float(energies[index - 1])) * HARTREE_TO_KCAL
            for index in range(int(peak_index) + 1, len(energies))
        ]
        endpoint_change = step_changes[-1]
        baseline = median(step_changes[:-1])
        ratio = float(selection_config.get("endpoint_jump_warning_ratio", 5.0))
        absolute_minimum = float(selection_config.get("endpoint_jump_warning_min_kcal", 5.0))
        return endpoint_change >= absolute_minimum and endpoint_change > ratio * max(baseline, 1.0e-12)

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

        frame_paths = [Path(path) for path in (getattr(scan_result, "geometries", None) or [])]
        if len(frame_paths) != len(energies):
            raise RuntimeError(
                "S2 scan trajectory is incomplete: "
                f"{len(energies)} energies but {len(frame_paths)} geometry frames"
            )

        topology_config = self._get_topology_guard_config()
        trajectory_quality: Dict[str, Any] = {
            "checked": False,
            "total_frames": len(frame_paths),
            "off_path_indices": [],
            "off_path_count": 0,
            "frame_issues": [],
        }
        if bool(topology_config.get("enabled", True)):
            trajectory_reference_coords, trajectory_reference_symbols = read_xyz(frame_paths[0])
            trajectory_quality = check_scan_trajectory(
                product_coords=np.asarray(trajectory_reference_coords, dtype=float),
                symbols=list(trajectory_reference_symbols),
                forming_bonds=bonds,
                frame_paths=frame_paths,
                graph_scale=float(topology_config.get("graph_scale", 1.25)),
            )
            trajectory_quality["reference_frame"] = str(frame_paths[0])
        off_path_indices = {
            int(index) for index in trajectory_quality.get("off_path_indices", [])
        }
        if max_idx in off_path_indices:
            status = "DEGRADED"
            ts_guess_confidence = "low"
            degraded_reasons.append("ts_guess_off_path")

        selection_config = dict(params.get("intermediate_selection", {}) or {})
        selection_error: Optional[str] = None
        try:
            intermediate_selection = self._select_intermediate_frame(
                energies=energies,
                frame_paths=frame_paths,
                peak_index=max_idx,
                off_path_indices=sorted(off_path_indices),
                selection_config=selection_config,
            )
            intermediate_idx: Optional[int] = int(intermediate_selection["index"])
        except RuntimeError as exc:
            selection_error = str(exc)
            intermediate_idx = None
            intermediate_selection = {
                "index": None,
                "rule": "stable_midpoint",
                "status": "unavailable",
                "error": selection_error,
            }
            status = "FAILED"
            ts_guess_confidence = "low"
            degraded_reasons.append("no_topology_valid_intermediate")

        if bool(intermediate_selection.get("degraded_fallback")):
            status = "DEGRADED"
            ts_guess_confidence = "low"
            degraded_reasons.append("intermediate_used_off_path_fallback")

        ts_guess_xyz_final = output_dir / "ts_guess.xyz"
        shutil.copy2(frame_paths[max_idx], ts_guess_xyz_final)

        dipolar_xyz = output_dir / INTERMEDIATE_XYZ
        reactant_xyz = output_dir / "reactant_complex.xyz"
        if intermediate_idx is not None:
            shutil.copy2(frame_paths[intermediate_idx], dipolar_xyz)
            shutil.copy2(dipolar_xyz, reactant_xyz)

        reaction_coordinate = compute_scan_distances(
            float(params["scan_start_distance"]),
            float(params["scan_end_distance"]),
            len(energies),
            direction="outward",
        )
        relative_energies = [
            (float(energy) - float(energies[0])) * HARTREE_TO_KCAL for energy in energies
        ]
        endpoint_jump = self._endpoint_jump_detected(energies, max_idx, selection_config)
        if endpoint_jump:
            self.logger.warning(
                "[S2] Endpoint energy discontinuity detected; endpoint frame is excluded from "
                "intermediate selection"
            )

        if intermediate_idx is not None:
            intermediate_selection.update(
                {
                    "target_distance_angstrom": float(reaction_coordinate[intermediate_idx]),
                    "energy_hartree": float(energies[intermediate_idx]),
                    "relative_energy_kcal_mol": float(relative_energies[intermediate_idx]),
                    "output_xyz": str(dipolar_xyz),
                }
            )
        ts_selection = {
            "index": int(max_idx),
            "frame_xyz": str(frame_paths[max_idx]),
            "output_xyz": str(ts_guess_xyz_final),
            "rule": "maximum_energy",
            "target_distance_angstrom": float(reaction_coordinate[max_idx]),
            "energy_hartree": float(energies[max_idx]),
            "relative_energy_kcal_mol": float(relative_energies[max_idx]),
        }
        trajectory_quality.update(
            {
                "endpoint_index": len(energies) - 1,
                "endpoint_excluded": True,
                "endpoint_jump_detected": bool(endpoint_jump),
            }
        )

        scan_profile_json = output_dir / "scan_profile.json"
        profile_payload: Dict[str, Any] = {
            "profile_schema_version": "s2_scan_profile_v2",
            "generation_method": "retro_scan",
            "product_xyz": str(product_file),
            "intermediate_xyz": str(dipolar_xyz) if intermediate_idx is not None else None,
            "forming_bonds": [list(pair) for pair in bonds],
            "scan_parameters": params,
            "scan_quality": {
                "max_energy_index": int(max_idx),
                "intermediate_index": intermediate_idx,
                "boundary_maximum": bool(boundary_max),
                "local_peak_ok": bool(peak_ok),
                "status": status,
                "ts_guess_confidence": ts_guess_confidence,
                "intermediate_confidence": "high" if status == "COMPLETE" else "low",
                "degraded_reasons": degraded_reasons,
            },
            "reaction_coordinate_angstrom": reaction_coordinate,
            "forming_bond_distances_angstrom": self._forming_bond_distances_by_frame(
                frame_paths, bonds
            ),
            "energies_hartree": energies,
            "relative_energies_kcal_mol": relative_energies,
            "selections": {
                "ts_guess": ts_selection,
                "intermediate": intermediate_selection,
            },
            "trajectory_quality": trajectory_quality,
            "scan_plot": None,
        }
        scan_profile_json.write_text(json.dumps(profile_payload, indent=2), encoding="utf-8")

        try:
            plot_path = plot_scan_profile(scan_profile_json)
            if plot_path is not None:
                profile_payload["scan_plot"] = str(plot_path)
                scan_profile_json.write_text(json.dumps(profile_payload, indent=2), encoding="utf-8")
        except Exception as exc:
            self.logger.warning("[S2] Failed to render scan profile: %s", exc, exc_info=True)

        if selection_error is not None:
            raise RuntimeError(
                f"{selection_error}; S2 scan diagnostics were written to {scan_profile_json}"
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

