import json
import logging
import math
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from rph_core.steps.step2_retro.path_profile import build_orca_scan_profile
from rph_core.steps.step2_retro.path_profile import build_xtb_path_profile
from rph_core.steps.step2_retro.path_selector import policy_from_config, select_path_seeds
from rph_core.utils.file_io import read_xyz
from rph_core.utils.bond_pairs import canonicalize_bond_pairs
from rph_core.utils.geometry_tools import GeometryUtils, LogParser, kabsch_rmsd
from rph_core.utils.json_io import write_text_atomic
from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.qc_interface import XTBInterface
from rph_core.utils.naming_compat import (
    INTERMEDIATE_XYZ,
)
from rph_core.utils.scan_profile_plotter import (
    HARTREE_TO_KCAL,
    compute_scan_distances,
    plot_scan_profile,
)

from .energy_refinement import ScanEnergyRefiner
from .geometry_guard import (
    check_scan_trajectory,
)
from .relaxed_scan_rescue import B97CRelaxedScanRescuer
from .scan_trajectory import CompositeProfileBuilder, ScanAttempt, attempt_manifest

logger = logging.getLogger(__name__)


class PEBScanEngine(LoggerMixin):
    DEFAULT_SCAN_START_DISTANCE = 3.5
    DEFAULT_SCAN_END_DISTANCE = 1.8
    DEFAULT_SCAN_STEPS = 20
    DEFAULT_SCAN_MODE = "concerted"
    DEFAULT_SCAN_FORCE_CONSTANT = 0.5
    DEFAULT_MIN_VALID_POINTS = 5
    def __init__(
        self,
        config: Dict[str, Any],
        molecule_name: Optional[str] = None,
        event_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
    ):
        self.config = config
        self.step2_cfg = config.get("step2", {}) if isinstance(config, dict) else {}
        self.molecule_name = molecule_name
        self.event_callback = event_callback
        self._active_scan_params: Dict[str, Any] = {}
        self.last_profile_payload: Optional[Dict[str, Any]] = None
        self.logger.info("[S2] PEB scan engine initialized")

    def _emit_progress(self, event: str, **fields: Any) -> None:
        if self.event_callback is None:
            return
        payload = {"variant": self.molecule_name or "product", **fields}
        try:
            self.event_callback(event, payload)
        except Exception as exc:  # pragma: no cover - UI isolation
            self.logger.warning("[S2] Ignoring UI callback failure for %s: %s", event, exc)

    @staticmethod
    def _archive_previous_outputs(output_dir: Path) -> Optional[Path]:
        """Remove stale published nodes while retaining their diagnostics."""

        published_names = (
            "manifest.json",
            "scan_profile.json",
            "scan_profile.png",
            "outputs",
            "ts_guess.xyz",
            "intermediate.xyz",
            "reactant_complex.xyz",
            "intermediate_seed.xyz",
            "intermediate_atom_mapping.json",
            "ts_atom_mapping.json",
            "retro_scan",
            "retro_scan_retry",
            "retro_scan_refined",
            "attempts",
        )
        existing = [output_dir / name for name in published_names if (output_dir / name).exists()]
        diagnostics_dir = output_dir / "diagnostics"
        if diagnostics_dir.exists():
            existing.append(diagnostics_dir)
        existing.extend(sorted(output_dir.glob("retro_scan_topology_retry_*")))
        if not existing:
            return None
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        archive_dir = output_dir / "superseded" / stamp
        archive_dir.mkdir(parents=True, exist_ok=False)
        for source in existing:
            shutil.move(str(source), str(archive_dir / source.name))
        return archive_dir

    @staticmethod
    def _atomic_copy(source: Path, destination: Path) -> Path:
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(f".{destination.name}.tmp")
        shutil.copy2(source, temporary)
        temporary.replace(destination)
        return destination

    @staticmethod
    def _require_energy(value: Optional[float], context: str) -> float:
        if value is None:
            raise RuntimeError(context)
        return float(value)

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
        try:
            normalized = canonicalize_bond_pairs(forming_bonds)
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"S2 received invalid forming_bonds: {forming_bonds!r}") from exc

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

        return normalized

    def _get_topology_guard_config(self) -> Dict[str, Any]:
        """Get topology guard configuration with defaults."""
        scan_cfg = dict(self.step2_cfg.get("scan", {}) or {})
        return {
            "enabled": scan_cfg.get("topology_guard_enabled", True),
            "graph_scale": scan_cfg.get("topology_graph_scale", 1.25),
            "reclassify_isolated_suspects": bool(
                scan_cfg.get("reclassify_isolated_topology_suspects", True)
            ),
            "isolated_max_energy_jump_kcal": float(
                scan_cfg.get("isolated_topology_max_energy_jump_kcal", 5.0)
            ),
            "isolated_max_reaction_core_rmsd_A": float(
                scan_cfg.get(
                    "isolated_topology_max_reaction_core_rmsd_A", 0.15
                )
            ),
        }

    def _resolve_scan_params(self, scan_config: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        scan_cfg = dict((self.step2_cfg.get("scan", {}) or {}))
        if scan_config:
            scan_cfg.update(scan_config)

        start = float(scan_cfg.get("scan_start_distance", scan_cfg.get("start_distance", self.DEFAULT_SCAN_START_DISTANCE)))
        end = float(scan_cfg.get("scan_end_distance", scan_cfg.get("end_distance", self.DEFAULT_SCAN_END_DISTANCE)))
        raw_steps = scan_cfg.get("scan_steps", scan_cfg.get("steps", "auto"))
        configured_steps = 0 if str(raw_steps).lower() == "auto" else int(raw_steps)
        coarse_step = float(
            scan_cfg.get("coarse_step_A", scan_cfg.get("max_step_A", 0.20)) or 0.20
        )
        if coarse_step <= 0.0:
            raise RuntimeError(f"S2 coarse_step_A must be positive, got {coarse_step}")
        derived_steps = int(math.ceil((start - end) / coarse_step)) + 1
        steps = derived_steps if configured_steps <= 1 else max(configured_steps, derived_steps)
        mode = str(scan_cfg.get("scan_mode", self.DEFAULT_SCAN_MODE))
        raw_force_constant = scan_cfg.get("scan_force_constant")
        force_constant = float(raw_force_constant) if raw_force_constant is not None else None
        min_valid_points = int(scan_cfg.get("min_valid_points", self.DEFAULT_MIN_VALID_POINTS))
        refinement_cfg = dict(scan_cfg.get("candidate_refinement", {}) or {})
        endpoint_extension_cfg = dict(scan_cfg.get("endpoint_extension", {}) or {})
        corridor_cfg = dict(scan_cfg.get("refinement_corridor", {}) or {})
        selection_cfg = dict(scan_cfg.get("selection", {}) or {})
        anchor_cfg = dict(scan_cfg.get("anchor_detection", {}) or {})

        if start <= end:
            raise RuntimeError(f"S2 scan requires scan_start_distance > scan_end_distance, got scan_start={start}, scan_end={end}")
        if steps <= 1:
            raise RuntimeError(f"S2 scan_steps must be > 1, got {steps}")

        return {
            "scan_start_distance": start,
            "scan_end_distance": end,
            "scan_steps": steps,
            "configured_scan_steps": configured_steps,
            "coarse_step_A": coarse_step,
            "max_step_A": coarse_step,
            "scan_mode": mode,
            "scan_force_constant": force_constant,
            "min_valid_points": min_valid_points,
            "scan_policy": scan_cfg.get("scan_policy", "policy_c"),
            "anchor_detection": {
                "persistent_drift_points": max(
                    1, int(anchor_cfg.get("persistent_drift_points", 2))
                ),
                "plateau_onset": {
                    "slope_ratio_to_ts": float(
                        (anchor_cfg.get("plateau_onset", {}) or {}).get(
                            "slope_ratio_to_ts", 0.25
                        )
                    ),
                    "slope_change_ratio_to_ts": float(
                        (anchor_cfg.get("plateau_onset", {}) or {}).get(
                            "slope_change_ratio_to_ts", 0.15
                        )
                    ),
                    "absolute_slope_ceiling_kcal_mol_A": float(
                        (anchor_cfg.get("plateau_onset", {}) or {}).get(
                            "absolute_slope_ceiling_kcal_mol_A", 12.0
                        )
                    ),
                    "minimum_consecutive_intervals": max(
                        1,
                        int(
                            (anchor_cfg.get("plateau_onset", {}) or {}).get(
                                "minimum_consecutive_intervals", 2
                            )
                        ),
                    ),
                    "allow_segmented_regression_fallback": bool(
                        (anchor_cfg.get("plateau_onset", {}) or {}).get(
                            "allow_segmented_regression_fallback", True
                        )
                    ),
                },
            },
            "candidate_refinement": {
                "enabled": bool(refinement_cfg.get("enabled", False)),
                "min_overlap_points": int(refinement_cfg.get("min_overlap_points", 2)),
                "coordinate_tolerance_A": float(
                    refinement_cfg.get("coordinate_tolerance_A", 0.002)
                ),
                "max_overlap_rmsd_A": float(
                    refinement_cfg.get("max_overlap_rmsd_A", 0.75)
                ),
                "max_reaction_core_rmsd_A": float(
                    refinement_cfg.get("max_reaction_core_rmsd_A", 0.50)
                ),
                "max_overlap_energy_gap_kcal": float(
                    refinement_cfg.get("max_overlap_energy_gap_kcal", 25.0)
                ),
                "max_overlap_shape_residual_kcal": float(
                    refinement_cfg.get("max_overlap_shape_residual_kcal", 15.0)
                ),
                "reaction_core_depth": int(
                    refinement_cfg.get("reaction_core_depth", 2)
                ),
            },
            "refinement_corridor": {
                "enabled": bool(corridor_cfg.get("enabled", True)),
                "strategy": str(
                    corridor_cfg.get("strategy", "unified_ts_int_corridor")
                ),
                "step_A": float(corridor_cfg.get("step_A", 0.05)),
                "left_margin_A": float(corridor_cfg.get("left_margin_A", 0.25)),
                "right_margin_A": float(corridor_cfg.get("right_margin_A", 0.25)),
                "topology_boundary_buffer_A": float(
                    corridor_cfg.get("topology_boundary_buffer_A", 0.05)
                ),
                "stop_on_persistent_topology_drift": bool(
                    corridor_cfg.get("stop_on_persistent_topology_drift", True)
                ),
                "require_seamless_ts_int_overlap": bool(
                    corridor_cfg.get("require_seamless_ts_int_overlap", True)
                ),
                "min_overlap_points": max(
                    1, int(corridor_cfg.get("min_overlap_points", 2))
                ),
                "coordinate_merge_tolerance_A": float(
                    corridor_cfg.get("coordinate_merge_tolerance_A", 0.01)
                ),
                "continuity_max_energy_jump_kcal": float(
                    corridor_cfg.get("continuity_max_energy_jump_kcal", 3.0)
                ),
                "continuity_max_geometry_rmsd_A": float(
                    corridor_cfg.get("continuity_max_geometry_rmsd_A", 0.35)
                ),
            },
            "endpoint_extension": {
                "enabled": bool(endpoint_extension_cfg.get("enabled", False)),
                "increment_A": float(endpoint_extension_cfg.get("increment_A", 0.25)),
                "maximum_coordinate_A": float(
                    endpoint_extension_cfg.get("maximum_coordinate_A", max(start, 4.2))
                ),
                "min_dissociation_side_span_A": float(
                    endpoint_extension_cfg.get("min_dissociation_side_span_A", 0.50)
                ),
                "min_valid_post_ts_points": int(
                    endpoint_extension_cfg.get("min_valid_post_ts_points", 5)
                ),
                "max_extensions": int(endpoint_extension_cfg.get("max_extensions", 3)),
            },
            "selection": {
                "preferred_energy_source": str(
                    selection_cfg.get("preferred_energy_source", "b973c")
                ).lower(),
                "ts_min_prominence_kcal_mol": max(
                    0.0,
                    float(selection_cfg.get("ts_min_prominence_kcal_mol", 0.40)),
                ),
                "ts_min_clean_neighbors": max(
                    1, int(selection_cfg.get("ts_min_clean_neighbors", 1))
                ),
                "valid_corridor_weak_peak_min_prominence_kcal_mol": max(
                    0.0,
                    float(
                        selection_cfg.get(
                            "valid_corridor_weak_peak_min_prominence_kcal_mol",
                            0.10,
                        )
                    ),
                ),
                "valid_corridor_weak_peak_min_barrier_kcal_mol": max(
                    0.0,
                    float(
                        selection_cfg.get(
                            "valid_corridor_weak_peak_min_barrier_kcal_mol",
                            3.0,
                        )
                    ),
                ),
                "ts_min_reactant_barrier_kcal_mol": max(
                    0.0,
                    float(selection_cfg.get("ts_min_reactant_barrier_kcal_mol", 3.0)),
                ),
                "require_intermediate": bool(selection_cfg.get("require_intermediate", False)),
                "max_nonreactive_scaffold_rmsd_A": max(
                    0.0,
                    float(selection_cfg.get("max_nonreactive_scaffold_rmsd_A", 0.75)),
                ),
                "full_endpoint_min_clean_frames_from_boundary": max(
                    1,
                    int(
                        selection_cfg.get(
                            "full_endpoint_min_clean_frames_from_boundary", 3
                        )
                    ),
                ),
                "ts_seed_reactant_backoff_A": max(
                    0.0,
                    float(selection_cfg.get("ts_seed_reactant_backoff_A", 0.0)),
                ),
                "int_min_basin_prominence_kcal_mol": max(
                    0.0,
                    float(
                        selection_cfg.get(
                            "int_min_basin_prominence_kcal_mol", 0.50
                        )
                    ),
                ),
                "int_plateau_fallback_enabled": bool(
                    selection_cfg.get("int_plateau_fallback_enabled", True)
                ),
                "int_plateau_min_consecutive_frames": max(
                    2, int(selection_cfg.get("int_plateau_min_consecutive_frames", 3))
                ),
                "int_plateau_min_ts_separation_A": max(
                    0.0,
                    float(selection_cfg.get("int_plateau_min_ts_separation_A", 0.10)),
                ),
                "int_plateau_energy_window_kcal_mol": max(
                    0.0,
                    float(selection_cfg.get("int_plateau_energy_window_kcal_mol", 2.0)),
                ),
                "int_plateau_barrier_fraction": max(
                    0.0,
                    float(selection_cfg.get("int_plateau_barrier_fraction", 0.25)),
                ),
                "int_plateau_max_slope_kcal_mol_A": max(
                    0.0,
                    float(selection_cfg.get("int_plateau_max_slope_kcal_mol_A", 40.0)),
                ),
                "endpoint_exclusion_frames": max(
                    0, int(selection_cfg.get("endpoint_exclusion_frames", 2))
                ),
                "min_reaction_progress": max(
                    0.0,
                    float(selection_cfg.get("min_reaction_progress", 0.35)),
                ),
                "min_valid_neighbor_window": max(
                    1, int(selection_cfg.get("min_valid_neighbor_window", 1))
                ),
                "allow_monotonic_shoulder": bool(
                    selection_cfg.get("allow_monotonic_shoulder", True)
                ),
                "shoulder_max_abs_slope_kcal_mol_per_A": max(
                    0.0,
                    float(
                        selection_cfg.get(
                            "shoulder_max_abs_slope_kcal_mol_per_A", 20.0
                        )
                    ),
                ),
                "shoulder_min_curvature_signal": max(
                    0.0,
                    float(selection_cfg.get("shoulder_min_curvature_signal", 0.05)),
                ),
                "allow_shared_search_seed": bool(
                    selection_cfg.get("allow_shared_search_seed", True)
                ),
            },
            "terminate_after_consecutive_off_path": int(
                scan_cfg.get("terminate_after_consecutive_off_path", 2)
            ),
        }

    @classmethod
    def _detect_scan_anchors(
        cls,
        energies: Sequence[Optional[float]],
        reaction_coordinate: Sequence[float],
        off_path_indices: Sequence[int],
        persistent_drift_points: int = 2,
        anchor_config: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Locate P, the gradient TS candidate, plateau onset M, drift D and endpoint E."""
        size = len(reaction_coordinate)
        if len(energies) != size or size < 3:
            raise RuntimeError("S2 anchor detection requires aligned energy/coordinate data")
        off_path = {int(index) for index in off_path_indices}
        drift_index = cls._persistent_off_path_start(
            sorted(off_path), size, persistent_drift_points
        )
        if drift_index is None and off_path:
            # Isolated suspects have already been reclassified by geometry_guard;
            # any remaining single point is therefore a real topology boundary.
            drift_index = min(off_path)
        endpoint_index = max(range(size), key=lambda index: float(reaction_coordinate[index]))
        pre_drift = sorted(
            [
            index
            for index in range(size)
            if index not in off_path
            and (drift_index is None or index < drift_index)
            and energies[index] is not None
            ],
            key=lambda index: float(reaction_coordinate[index]),
        )
        if len(pre_drift) < 3:
            raise RuntimeError("S2 cannot define P-M-D-E anchors before topology drift")
        absolute_maximum_index = max(
            pre_drift,
            key=lambda index: cls._require_energy(
                energies[index], "S2 anchor detection missing pre-drift energy"
            ),
        )
        absolute_maximum_coordinate = float(
            reaction_coordinate[absolute_maximum_index]
        )
        product_domain = [
            index
            for index in pre_drift
            if float(reaction_coordinate[index])
            <= absolute_maximum_coordinate + 1.0e-10
        ]
        product_index = min(
            product_domain,
            key=lambda index: cls._require_energy(
                energies[index], "S2 anchor detection missing product-domain energy"
            ),
        )
        product_coordinate = float(reaction_coordinate[product_index])
        rising_domain = [
            index
            for index in pre_drift
            if float(reaction_coordinate[index]) >= product_coordinate - 1.0e-10
        ]
        if len(rising_domain) < 3:
            raise RuntimeError("S2 cannot define a rising domain after the product minimum")

        x_values = np.asarray(
            [float(reaction_coordinate[index]) for index in rising_domain], dtype=float
        )
        y_values = np.asarray(
            [
                cls._require_energy(
                    energies[index], "S2 anchor detection missing rising-domain energy"
                )
                * HARTREE_TO_KCAL
                for index in rising_domain
            ],
            dtype=float,
        )
        if len(np.unique(x_values)) != len(x_values):
            raise RuntimeError("S2 anchor detection requires unique reaction coordinates")
        gradients = np.gradient(y_values, x_values, edge_order=2 if len(x_values) > 2 else 1)
        positive_positions = [
            position
            for position in range(1, len(rising_domain) - 1)
            if float(gradients[position]) > 0.0
        ]
        if not positive_positions:
            raise RuntimeError("S2 cannot define TS1: no positive coarse-scan gradient")
        ts_position = max(positive_positions, key=lambda position: float(gradients[position]))
        ts_index = int(rising_domain[ts_position])
        ts_gradient = abs(float(gradients[ts_position]))

        plateau_cfg = dict((anchor_config or {}).get("plateau_onset", {}) or {})
        slope_ratio = max(0.0, float(plateau_cfg.get("slope_ratio_to_ts", 0.25)))
        change_ratio = max(
            0.0, float(plateau_cfg.get("slope_change_ratio_to_ts", 0.15))
        )
        absolute_ceiling = max(
            0.0,
            float(plateau_cfg.get("absolute_slope_ceiling_kcal_mol_A", 12.0)),
        )
        minimum_run = max(
            1, int(plateau_cfg.get("minimum_consecutive_intervals", 2))
        )
        slope_limit = max(absolute_ceiling, slope_ratio * ts_gradient)
        change_limit = max(1.0e-8, change_ratio * ts_gradient)
        stable_positions: List[int] = []
        for position in range(ts_position + 1, len(rising_domain)):
            previous = float(gradients[position - 1])
            current = float(gradients[position])
            if abs(current) <= slope_limit and abs(current - previous) <= change_limit:
                stable_positions.append(position)

        plateau_position: Optional[int] = None
        for start_position in stable_positions:
            if all(
                (start_position + offset) in stable_positions
                for offset in range(minimum_run)
            ):
                plateau_position = start_position
                break
        plateau_method = "sustained_low_slope"
        plateau_confidence = "high"

        if plateau_position is None and bool(
            plateau_cfg.get("allow_segmented_regression_fallback", True)
        ):
            segmented_candidates: List[Tuple[float, int]] = []
            for split in range(ts_position + 2, len(rising_domain) - 1):
                left_x, left_y = x_values[: split + 1], y_values[: split + 1]
                right_x, right_y = x_values[split:], y_values[split:]
                if len(left_x) < 3 or len(right_x) < 2:
                    continue
                left_fit = np.polyfit(left_x, left_y, 1)
                right_fit = np.polyfit(right_x, right_y, 1)
                if abs(float(right_fit[0])) > max(
                    absolute_ceiling, slope_ratio * abs(float(left_fit[0]))
                ):
                    continue
                residual = float(
                    np.sum(np.square(left_y - np.polyval(left_fit, left_x)))
                    + np.sum(np.square(right_y - np.polyval(right_fit, right_x)))
                )
                segmented_candidates.append((residual, split))
            if segmented_candidates:
                _, plateau_position = min(segmented_candidates)
                plateau_method = "segmented_regression_fallback"
                plateau_confidence = "medium"

        if plateau_position is None:
            fallback_positions = list(range(ts_position + 1, len(rising_domain)))
            if not fallback_positions:
                raise RuntimeError("S2 cannot define plateau onset after TS1")
            plateau_position = min(
                fallback_positions,
                key=lambda position: (
                    abs(float(gradients[position])) / max(ts_gradient, 1.0e-8)
                    + abs(float(gradients[position]) - float(gradients[position - 1]))
                    / max(ts_gradient, 1.0e-8)
                ),
            )
            plateau_method = "minimum_slope_variation_fallback"
            plateau_confidence = "low"

        plateau_index = int(rising_domain[plateau_position])
        last_valid_index = pre_drift[-1]
        boundary_index = drift_index if drift_index is not None else endpoint_index
        return {
            "product_index": int(product_index),
            "coarse_ts_index": ts_index,
            "plateau_onset_index": plateau_index,
            "absolute_energy_maximum_index": int(absolute_maximum_index),
            "topology_drift_index": int(drift_index) if drift_index is not None else None,
            "last_valid_before_drift_index": int(last_valid_index),
            "scan_endpoint_index": int(endpoint_index),
            "intermediate_boundary_index": int(boundary_index),
            "boundary_source": "topology_drift" if drift_index is not None else "scan_endpoint",
            "plateau_onset_method": plateau_method,
            "plateau_onset_confidence": plateau_confidence,
            "coarse_gradients_kcal_per_mol_A": {
                str(index): float(gradients[position])
                for position, index in enumerate(rising_domain)
            },
            "ts_gradient_kcal_per_mol_A": float(gradients[ts_position]),
            "plateau_slope_limit_kcal_per_mol_A": slope_limit,
            "plateau_slope_change_limit_kcal_per_mol_A": change_limit,
        }

    @staticmethod
    def _select_node_by_anchor_index(
        anchors: Dict[str, Any],
        anchor_key: str,
        frame_paths: Sequence[Path],
        reaction_coordinate: Sequence[float],
        energies: Sequence[Optional[float]],
        off_path_indices: Sequence[int],
    ) -> Dict[str, Any]:
        """Pick the frame at the named anchor index (e.g. plateau_onset_index).

        Topology-invalid anchor frames fall back to the nearest topology-valid
        neighbour on the plateau side. This is the only TS selection rule:
        the plateau onset IS the transition-state region entry, and S3 OptTS
        relaxes it to the precise TS once the scan coordinate is freed.
        """
        anchor_index = int(anchors[anchor_key])
        off_path = {int(index) for index in off_path_indices}
        if anchor_index in off_path:
            eligible = [
                index
                for index in range(len(reaction_coordinate))
                if index not in off_path and energies[index] is not None
            ]
            if eligible:
                anchor_index = min(
                    eligible,
                    key=lambda index: abs(float(reaction_coordinate[index]) - float(reaction_coordinate[anchor_index])),
                )
        return {
            "index": int(anchor_index),
            "frame_xyz": str(Path(frame_paths[anchor_index])),
            "rule": "anchor:" + anchor_key,
            "target_coordinate_A": float(reaction_coordinate[anchor_index]),
            "candidate_indices": [int(anchor_index)],
        }

    @staticmethod
    def _select_reactant_side_ts_seed(
        *,
        peak_index: int,
        forming_bond_distances: Sequence[Optional[Sequence[Optional[float]]]],
        invalid_indices: Sequence[int],
        backoff_A: float,
        path_arclength: Optional[Sequence[float]] = None,
    ) -> Tuple[int, Dict[str, Any]]:
        """Compatibility helper for historical callers; not used by ``run``.

        The production path now gets both TS and INT seeds from
        ``select_path_seeds``.  Keep this helper temporarily for old fixtures
        and offline comparisons, but do not use it to publish S2 decisions.

        Choose an S3 seed before an energy-located TS peak along PATH.

        The energy peak remains the TS locator. This only chooses a nearby
        PATH geometry for S3 when the peak frame is too far along the forming-
        bond coordinate.  PATH order, not the mean forming-bond distance, is
        authoritative for deciding which side is reactant-like: asynchronous
        bond formation can make the mean distance non-monotonic near product.
        """

        def _mean_distance(index: int) -> Optional[float]:
            if index < 0 or index >= len(forming_bond_distances):
                return None
            distances = forming_bond_distances[index]
            if not distances or any(value is None for value in distances):
                return None
            valid_distances = [float(value) for value in distances if value is not None]
            if len(valid_distances) != len(distances):
                return None
            return float(sum(valid_distances) / len(valid_distances))

        invalid = {int(index) for index in invalid_indices}
        usable = [
            index
            for index in range(len(forming_bond_distances))
            if index not in invalid and _mean_distance(index) is not None
        ]
        if path_arclength is not None and len(path_arclength) == len(
            forming_bond_distances
        ):
            coordinates = [float(value) for value in path_arclength]
        else:
            coordinates = [float(index) for index in range(len(forming_bond_distances))]
        ordered = sorted(usable, key=lambda index: coordinates[index])
        peak_mean = _mean_distance(int(peak_index))
        metadata: Dict[str, Any] = {
            "energy_peak_index": int(peak_index),
            "seed_index": int(peak_index),
            "seed_rule": "energy_peak_geometry",
            "seed_backoff_requested_A": float(backoff_A),
            "seed_backoff_applied_A": 0.0,
            "seed_target_mean_forming_bond_distance_A": peak_mean,
            "seed_mean_forming_bond_distance_A": peak_mean,
        }
        if (
            backoff_A <= 0.0
            or peak_mean is None
            or int(peak_index) not in ordered
        ):
            return int(peak_index), metadata

        peak_position = ordered.index(int(peak_index))
        pre_peak = ordered[:peak_position]
        if not pre_peak:
            metadata["seed_backoff_status"] = "no_clean_pre_peak_frame"
            return int(peak_index), metadata

        peak_distances = forming_bond_distances[int(peak_index)]
        if not peak_distances or any(value is None for value in peak_distances):
            metadata["seed_backoff_status"] = "peak_bond_coordinate_unavailable"
            return int(peak_index), metadata
        peak_distance_values = [float(value) for value in peak_distances if value is not None]
        target_mean = float(peak_mean + backoff_A)
        candidates: List[int] = []
        for index in pre_peak:
            candidate_distances = forming_bond_distances[index]
            if (
                not candidate_distances
                or len(candidate_distances) != len(peak_distance_values)
                or any(value is None for value in candidate_distances)
            ):
                continue
            candidate_values = [float(value) for value in candidate_distances if value is not None]
            if all(
                candidate_distance >= peak_distance - 1.0e-6
                for candidate_distance, peak_distance in zip(
                    candidate_values,
                    peak_distance_values,
                )
            ):
                candidates.append(index)
        if not candidates:
            boundary_index = int(pre_peak[0])
            boundary_mean = _mean_distance(boundary_index)
            metadata.update(
                {
                    "seed_index": boundary_index,
                    "seed_rule": "path_order_boundary_fallback",
                    "seed_backoff_status": "no_bond_consistent_pre_peak_frame",
                    "seed_mean_forming_bond_distance_A": boundary_mean,
                    "clean_path_boundary_index": boundary_index,
                }
            )
            return boundary_index, metadata

        def _distance_or_fail(index: int) -> float:
            mean_distance = _mean_distance(index)
            if mean_distance is None:
                raise RuntimeError("S2 TS seed backoff candidate lacks bond distances")
            return float(mean_distance)

        seed_index = min(
            candidates,
            key=lambda index: (
                abs(_distance_or_fail(index) - target_mean),
                -_distance_or_fail(index),
            ),
        )
        seed_mean = _distance_or_fail(seed_index)
        metadata.update(
            {
                "seed_index": int(seed_index),
                "seed_rule": "reactant_side_mean_distance_backoff",
                "seed_backoff_status": "applied",
                "seed_backoff_applied_A": float(seed_mean - peak_mean),
                "seed_target_mean_forming_bond_distance_A": target_mean,
                "seed_mean_forming_bond_distance_A": seed_mean,
            }
        )
        return int(seed_index), metadata

    def _select_refined_path_nodes(
        self,
        *,
        anchors: Dict[str, Any],
        frame_paths: Sequence[Path],
        reaction_coordinate: Sequence[float],
        method_energies: Sequence[Optional[float]],
        off_path_indices: Sequence[int],
        path_arclength: Optional[np.ndarray],
        selection_config: Mapping[str, Any],
    ) -> Tuple[Dict[str, Any], Dict[str, Any], int, Optional[int], str]:
        """Compatibility helper for historical callers; not used by ``run``.

        Production selection is centralized in ``select_path_seeds``; this
        method remains only for old fixtures and offline comparisons.

        Select a curve-consistent TS and an S3 INT-search seed.

        A PATH geometry is not itself a stationary point.  TS selection therefore
        uses only extrema on the method-consistent energy curve.  The optional
        INT structure is either a resolvable pre-TS basin or, after TS seeding,
        a late low-gradient pre-TS platform selected by the caller.
        """
        invalid = {int(index) for index in off_path_indices}
        valid_indices = [
            index
            for index, energy in enumerate(method_energies)
            if energy is not None and index not in invalid
        ]
        if not valid_indices:
            raise RuntimeError("S2 cannot select PATH nodes without valid energies")

        if path_arclength is not None and len(path_arclength) == len(frame_paths):
            coordinates = [float(path_arclength[index]) for index in range(len(frame_paths))]
        else:
            coordinates = [float(value) for value in reaction_coordinate]
        ordered = sorted(valid_indices, key=lambda index: coordinates[index])
        energies_kcal = {
            index: self._require_energy(
                method_energies[index],
                "S2 refined PATH selection missing method energy",
            )
            * HARTREE_TO_KCAL
            for index in valid_indices
        }
        ts_prominence_cutoff = float(
            selection_config.get("ts_min_prominence_kcal_mol", 0.40)
        )
        int_prominence_cutoff = float(
            selection_config.get("int_min_basin_prominence_kcal_mol", 0.50)
        )

        minimum_neighbors = max(
            1, int(selection_config.get("ts_min_clean_neighbors", 1))
        )
        maxima: List[Tuple[int, float]] = []
        weak_maxima: List[Tuple[int, float]] = []
        for position in range(minimum_neighbors, len(ordered) - minimum_neighbors):
            left, index, right = ordered[position - 1 : position + 2]
            rise = energies_kcal[index] - energies_kcal[left]
            fall = energies_kcal[index] - energies_kcal[right]
            if rise > 0.0 and fall > 0.0:
                prominence = min(rise, fall)
                weak_maxima.append((index, prominence))
                if prominence >= ts_prominence_cutoff:
                    maxima.append((index, prominence))

        if maxima:
            ts_index, ts_prominence = max(
                maxima,
                # A multi-step reaction can expose more than one maximum.  The
                # first significant maximum from the reactant-side endpoint is
                # the relevant seed for this product-forming PATH.
                key=lambda item: (coordinates[item[0]], -item[1]),
            )
            ts_rule = "refined_curve_local_maximum"
            ts_confidence = "high" if ts_prominence >= 0.50 else "medium"
            ts_candidates = [index for index, _ in maxima]
        elif weak_maxima:
            ts_index, ts_prominence = max(
                weak_maxima,
                key=lambda item: (coordinates[item[0]], -item[1]),
            )
            ts_rule = "refined_curve_weak_local_maximum"
            ts_confidence = "low"
            ts_candidates = [index for index, _ in weak_maxima]
        else:
            interior = ordered[1:-1]
            if interior:
                ts_index = max(interior, key=lambda index: energies_kcal[index])
                ts_rule = "refined_curve_internal_maximum_fallback"
                ts_confidence = "low"
                ts_candidates = list(interior)
            else:
                fallback = anchors.get("plateau_onset_index", ordered[0])
                ts_index = int(fallback) if int(fallback) in valid_indices else ordered[0]
                ts_rule = "anchor_fallback_no_curve_extremum"
                ts_confidence = "low"
                ts_candidates = [ts_index]

        ts_position = ordered.index(ts_index)
        basin_candidates: List[Tuple[int, float]] = []
        for position in range(1, ts_position):
            left, index, right = ordered[position - 1 : position + 2]
            left_wall = energies_kcal[left] - energies_kcal[index]
            right_wall = energies_kcal[right] - energies_kcal[index]
            prominence = min(left_wall, right_wall)
            if (
                left_wall > 0.0
                and right_wall > 0.0
                and prominence >= int_prominence_cutoff
            ):
                basin_candidates.append((index, prominence))

        intermediate_idx: Optional[int] = None
        if basin_candidates:
            intermediate_idx, basin_prominence = max(
                basin_candidates,
                key=lambda item: (item[1], -energies_kcal[item[0]]),
            )
            int_mode = "stable_basin_candidate"
            int_rule = "refined_curve_pre_ts_local_minimum"
            int_candidates = [index for index, _ in basin_candidates]
            int_extra = {"basin_prominence_kcal_mol": float(basin_prominence)}
        else:
            # A monotonic path does not justify an arbitrary midpoint MIN
            # calculation.  Preserve the fact that no distinct INT was found,
            # and let downstream scheduling omit the duplicate S3 structure.
            intermediate_idx = None
            int_mode = "shared_ts_fallback"
            int_rule = "no_resolved_pre_ts_basin_shared_ts_seed"
            int_candidates = [int(ts_index)]
            int_extra = {
                "reason": "no_resolved_pre_ts_basin",
                "shared_ts_index": int(ts_index),
                "s3_job_required": False,
            }

        ts_selection = {
            "index": int(ts_index),
            "frame_xyz": str(Path(frame_paths[ts_index])),
            "rule": ts_rule,
            "actual_method": "refined_path_curve",
            "confidence": ts_confidence,
            "prominence_kcal_mol": float(
                next(
                    (
                        prominence
                        for index, prominence in maxima + weak_maxima
                        if index == ts_index
                    ),
                    0.0,
                )
            ),
            "target_coordinate_A": float(reaction_coordinate[ts_index]),
            "candidate_indices": [int(index) for index in ts_candidates],
        }
        int_selection: Dict[str, Any] = {
            "index": (
                int(ts_index)
                if int_mode == "shared_ts_fallback"
                else None if intermediate_idx is None else int(intermediate_idx)
            ),
            "frame_xyz": (
                str(Path(frame_paths[ts_index]))
                if int_mode == "shared_ts_fallback"
                else None if intermediate_idx is None else str(Path(frame_paths[intermediate_idx]))
            ),
            "rule": int_rule,
            "selection_mode": int_mode,
            "selection_status": (
                "shared_with_ts"
                if int_mode == "shared_ts_fallback"
                else "selected" if intermediate_idx is not None else "unavailable"
            ),
            "stationary_point_claimed": False,
            "target_coordinate_A": (
                float(reaction_coordinate[ts_index])
                if int_mode == "shared_ts_fallback"
                else None if intermediate_idx is None else float(reaction_coordinate[intermediate_idx])
            ),
            "candidate_indices": [int(index) for index in int_candidates],
            **int_extra,
        }
        return ts_selection, int_selection, int(ts_index), intermediate_idx, (
            "shared_with_ts"
            if int_mode == "shared_ts_fallback"
            else "selected" if intermediate_idx is not None else "unavailable"
        )

    @staticmethod
    def _path_requires_relaxed_scan(
        trajectory_quality: Mapping[str, Any],
    ) -> bool:
        """Return whether PEB supplied no usable path geometry at all.

        S2 selects *search seeds*, rather than proving stationary points.  A
        weak peak, a platform, or shared TS/INT frame is therefore still a
        valid S3 input when the underlying PATH is usable.  The expensive
        B97-3c relaxed scan is reserved for the one unambiguous PEB failure:
        topology drift persists from the first frame, leaving no clean path
        segment from which either seed can be selected.
        """

        if not bool(trajectory_quality.get("checked", False)):
            return False
        persistent_start = trajectory_quality.get(
            "persistent_off_path_start",
            trajectory_quality.get("topology_drift_index"),
        )
        usable_end = trajectory_quality.get(
            "usable_end_index",
            trajectory_quality.get("last_valid_before_drift_index"),
        )
        return (
            persistent_start == 0
            and usable_end in {-1, None}
        )

    @staticmethod
    def _nonreactive_scaffold_admission(
        reference_xyz: Path,
        candidate_xyz: Path,
        forming_bonds: Sequence[Tuple[int, int]],
        maximum_rmsd: float,
    ) -> Dict[str, Any]:
        """Reject a PATH frame that distorts the mapped non-reactive scaffold."""

        try:
            reference, _ = read_xyz(Path(reference_xyz))
            candidate, _ = read_xyz(Path(candidate_xyz))
            if reference.shape != candidate.shape:
                return {"accepted": False, "reason": "atom_count_mismatch"}
            reactive_atoms = {int(atom) for pair in forming_bonds for atom in pair}
            scaffold = [index for index in range(len(reference)) if index not in reactive_atoms]
            if len(scaffold) < 3:
                return {
                    "accepted": True,
                    "nonreactive_scaffold_rmsd_A": None,
                    "reason": "insufficient_nonreactive_atoms",
                }
            rmsd = float(kabsch_rmsd(reference[scaffold], candidate[scaffold]))
            return {
                "accepted": rmsd <= maximum_rmsd,
                "nonreactive_scaffold_rmsd_A": rmsd,
                "maximum_nonreactive_scaffold_rmsd_A": maximum_rmsd,
            }
        except (OSError, ValueError, np.linalg.LinAlgError) as exc:
            return {"accepted": False, "reason": f"geometry_check_failed:{exc}"}


    @staticmethod
    def _select_late_pre_ts_platform_seed(
        *,
        frame_paths: Sequence[Path],
        reaction_coordinate: Sequence[float],
        method_energies: Sequence[Optional[float]],
        off_path_indices: Sequence[int],
        path_arclength: Optional[np.ndarray],
        ts_peak_index: int,
        ts_seed_index: int,
        selection_config: Mapping[str, Any],
    ) -> Tuple[Optional[int], Dict[str, Any]]:
        """Compatibility helper for historical callers; not used by ``run``.

        Production selection is centralized in ``select_path_seeds``; this
        method remains only for old fixtures and offline comparisons.

        Choose the latest clean, low-gradient platform before the TS seed.

        A monotonic PATH can still contain a chemically useful precursor-side
        plateau even though it has no stationary local minimum.  This routine
        deliberately returns a search seed only; it never asserts that an INT
        has been located.
        """
        if not bool(selection_config.get("int_plateau_fallback_enabled", True)):
            return None, {"reason": "int_plateau_fallback_disabled"}
        if not (0 <= int(ts_peak_index) < len(frame_paths)) or not (
            0 <= int(ts_seed_index) < len(frame_paths)
        ):
            return None, {"reason": "invalid_ts_index_for_int_platform"}

        invalid = {int(index) for index in off_path_indices}
        if path_arclength is not None and len(path_arclength) == len(frame_paths):
            coordinates = [float(path_arclength[index]) for index in range(len(frame_paths))]
        else:
            coordinates = [float(value) for value in reaction_coordinate]
        valid = [
            index
            for index, energy in enumerate(method_energies)
            if energy is not None and index not in invalid
        ]
        ordered = sorted(valid, key=lambda index: coordinates[index])
        if int(ts_seed_index) not in ordered or int(ts_peak_index) not in ordered:
            return None, {"reason": "ts_not_in_clean_path_segment"}

        seed_position = ordered.index(int(ts_seed_index))
        pre_ts = ordered[:seed_position]
        minimum_frames = max(
            2, int(selection_config.get("int_plateau_min_consecutive_frames", 3))
        )
        if len(pre_ts) < minimum_frames:
            return None, {
                "reason": "insufficient_clean_pre_ts_frames",
                "clean_pre_ts_frame_count": len(pre_ts),
            }

        energy_kcal = {
            index: PEBScanEngine._require_energy(
                method_energies[index],
                "S2 late pre-TS platform selection missing method energy",
            )
            * HARTREE_TO_KCAL
            for index in ordered
        }
        pre_ts_floor = min(energy_kcal[index] for index in pre_ts)
        barrier = max(0.0, energy_kcal[int(ts_peak_index)] - pre_ts_floor)
        energy_window = max(
            float(selection_config.get("int_plateau_energy_window_kcal_mol", 2.0)),
            float(selection_config.get("int_plateau_barrier_fraction", 0.25)) * barrier,
        )
        min_separation = float(
            selection_config.get("int_plateau_min_ts_separation_A", 0.10)
        )
        energy_limit = pre_ts_floor + energy_window
        eligible = [
            index
            for index in pre_ts
            if energy_kcal[index] <= energy_limit
            and coordinates[int(ts_seed_index)] - coordinates[index] + 1.0e-8 >= min_separation
        ]
        eligible_set = set(eligible)
        groups: List[List[int]] = []
        current: List[int] = []
        for index in pre_ts:
            if index in eligible_set:
                current.append(index)
            elif current:
                groups.append(current)
                current = []
        if current:
            groups.append(current)

        max_slope = float(
            selection_config.get("int_plateau_max_slope_kcal_mol_A", 40.0)
        )
        accepted: List[Tuple[List[int], float]] = []
        for group in groups:
            if len(group) < minimum_frames:
                continue
            slopes = []
            for left, right in zip(group, group[1:]):
                delta_s = coordinates[right] - coordinates[left]
                if delta_s <= 1.0e-8:
                    continue
                slopes.append(abs(energy_kcal[right] - energy_kcal[left]) / delta_s)
            if max_slope > 0.0 and slopes and max(slopes) > max_slope:
                continue
            accepted.append((group, max(slopes, default=0.0)))

        if not accepted:
            return None, {
                "reason": "no_continuous_low_gradient_pre_ts_platform",
                "clean_pre_ts_frame_count": len(pre_ts),
                "energy_window_kcal_mol": energy_window,
                "energy_limit_kcal_mol": energy_limit,
            }

        group, group_max_slope = max(
            accepted, key=lambda item: coordinates[item[0][-1]]
        )
        index = int(group[-1])
        return index, {
            "rule": "late_pre_ts_low_gradient_platform",
            "selection_mode": "late_pre_ts_platform_fallback",
            "selection_status": "selected",
            "stationary_point_claimed": False,
            "candidate_indices": [int(value) for value in group],
            "pre_ts_floor_kcal_mol": pre_ts_floor,
            "barrier_from_pre_ts_floor_kcal_mol": barrier,
            "energy_window_kcal_mol": energy_window,
            "energy_limit_kcal_mol": energy_limit,
            "platform_max_slope_kcal_mol_A": group_max_slope,
            "clean_pre_ts_frame_count": len(pre_ts),
            "s3_job_required": True,
        }


    def _xtb_path_config(self) -> Dict[str, Any]:
        return dict((self.step2_cfg.get("xtb_path", {}) or {}))

    def _xtb_path_settings(self) -> Tuple[int, int]:
        xtb_settings = self.step2_cfg.get("xtb_settings", {}) or {}
        charge = int(xtb_settings.get("charge", 0))
        spin = int(xtb_settings.get("multiplicity", 1))
        return charge, spin

    @staticmethod
    def _path_mean_forming_bond_distances(
        frame_paths: Sequence[Path],
        forming_bonds: Sequence[Tuple[int, int]],
    ) -> List[Optional[float]]:
        """Mean forming-bond distance per PATH frame; None when unreadable.

        PATH frames come from a meta-dynamics relaxation, so the forming-bond
        distance is the natural reaction coordinate that lets the PATH curve
        overlay the coarse PEB scan curve on the same x-axis.
        """
        means: List[Optional[float]] = []
        for frame_path in frame_paths:
            try:
                frame_coords, _ = read_xyz(Path(frame_path))
                distances = [
                    float(GeometryUtils.calculate_distance(frame_coords, int(i), int(j)))
                    for i, j in forming_bonds
                ]
                means.append(float(sum(distances) / max(1, len(distances))))
            except Exception:
                means.append(None)
        return means

    def _path_arclength(self, frame_paths: Sequence[Path]) -> np.ndarray:
        """Cumulative Kabsch-aligned RMSD path-progress coordinate s.

        For PATH meta-dynamics profiles, the forming-bond mean distance is
        insufficient as a reaction coordinate (different geometries can share
        the same mean).  This method computes the proper path-progress
        variable:

            s[0] = 0
            s[i] = s[i-1] + kabsch_rmsd(frame[i-1], frame[i])

        s is used for gradient computations in TS/INT selection, replacing
        the default forming-bond-mean axis for PATH data.
        """
        n = len(frame_paths)
        if n < 2:
            return np.zeros(n, dtype=float)
        s = np.zeros(n, dtype=float)
        prev_coords, _ = read_xyz(Path(frame_paths[0]))
        for i in range(1, n):
            curr_coords, _ = read_xyz(Path(frame_paths[i]))
            step = kabsch_rmsd(np.asarray(prev_coords, dtype=float),
                               np.asarray(curr_coords, dtype=float))
            s[i] = s[i - 1] + step
            prev_coords = curr_coords
        return s

    def _path_frame_xtb_sp(
        self,
        frame_paths: Sequence[Path],
        output_dir: Path,
        charge: int,
        spin: int,
    ) -> List[Optional[float]]:
        """Single-point GFN2-xTB energy per PATH frame.

        ``PathSearchResult`` from xTB ``--path`` carries geometries but no
        per-frame energy table, so each frame is re-evaluated at GFN2-xTB SP
        to populate the PATH energy curve. Cheap (~1-3 s per frame).
        """
        import shutil
        import tempfile

        from rph_core.utils.xtb_runner import XTBRunner
        from rph_core.utils.resource_utils import resolve_executable_config
        from rph_core.utils.qc_interface import is_path_toxic

        output_dir.mkdir(parents=True, exist_ok=True)
        xtb_settings = self.step2_cfg.get("xtb_settings", {}) or {}
        path_cfg = self._xtb_path_config()
        gfn_level = int(path_cfg.get("gfn_level", 2))
        solvent = str(xtb_settings.get("solvent", "acetone"))
        nproc = int(xtb_settings.get("nproc") or self.config.get("resources", {}).get("nproc", 1))
        uhf = max(spin - 1, 0)

        xtb_cfg = resolve_executable_config(self.config, "xtb", env_vars=["XTB_PATH", "XTB_BIN"])
        xtb_exec = xtb_cfg["path"]

        use_sandbox = is_path_toxic(output_dir)
        work_dir = output_dir
        if use_sandbox:
            work_dir = Path(tempfile.mkdtemp(prefix="RPH_XTB_SP_", dir="/tmp"))
            self.logger.info("[S2] xTB SP output path is toxic; using sandbox: %s", output_dir)

        runner = XTBRunner(config=self.config, work_dir=work_dir)

        energies: List[Optional[float]] = []
        try:
            for index, frame in enumerate(frame_paths):
                try:
                    frame = Path(frame)
                    frame_log = work_dir / f"frame_{index:03d}.log"
                    cmd = [xtb_exec, frame.name, "--sp"]
                    cmd.extend(["-P", str(nproc)])
                    cmd.extend(["--chrg", str(charge)])
                    if uhf > 0:
                        cmd.extend(["--uhf", str(uhf)])
                    cmd.extend(["--gfn", str(gfn_level), "--alpb", solvent])

                    local_frame = work_dir / frame.name
                    if not local_frame.exists() or local_frame.resolve() != frame.resolve():
                        local_frame.write_text(frame.read_text())

                    runner._run_command(cmd, log_file=frame_log)
                    energy = runner._parse_energy(frame_log)
                    energies.append(float(energy) if energy is not None else None)
                except Exception as exc:
                    self.logger.warning("[S2] xTB SP on PATH frame %d failed: %s", index, exc)
                    energies.append(None)
        finally:
            if use_sandbox and work_dir != output_dir:
                try:
                    for item in work_dir.iterdir():
                        target = output_dir / item.name
                        if item.is_dir():
                            if target.exists():
                                shutil.rmtree(target)
                            shutil.copytree(item, target)
                        else:
                            shutil.copy2(item, target)
                    shutil.rmtree(work_dir, ignore_errors=True)
                except Exception as exc:
                    self.logger.warning("[S2] xTB SP sandbox cleanup failed: %s", exc)

        valid = sum(1 for e in energies if e is not None)
        self.logger.info("[S2] xTB SP on PATH: %d/%d frames with energies", valid, len(energies))
        return energies

    def _execute_xtb_path(
        self,
        start_xyz: Path,
        end_xyz: Path,
        output_dir: Path,
        forming_bonds: Sequence[Tuple[int, int]],
        *,
        path_name: str = "valid_corridor",
        topology_config: Optional[Dict[str, Any]] = None,
        reference_coords: Optional[Sequence[Sequence[float]]] = None,
        reference_symbols: Optional[Sequence[str]] = None,
        reference_path: Optional[Path] = None,
        consecutive_off_path: int = 2,
    ) -> Tuple[ScanAttempt, Dict[str, Any]]:
        """Run xTB2 meta-dynamics PATH and wrap the result as a ScanAttempt.

        Returns (attempt, path_metadata). Raises RuntimeError on PATH failure
        so the caller can fall back to coarse-only selection.
        """
        path_cfg = self._xtb_path_config()
        if not bool(path_cfg.get("enabled", True)):
            raise RuntimeError("S2 xtb_path disabled in config")
        npoint_key = f"npoint_{path_name}"
        npoint = int(path_cfg.get(npoint_key, path_cfg.get("npoint", 28)))
        charge, spin = self._xtb_path_settings()

        xtb_settings = self.step2_cfg.get("xtb_settings", {}) or {}
        solvent = str(xtb_settings.get("solvent", "acetone"))
        nproc = int(xtb_settings.get("nproc") or self.config.get("resources", {}).get("nproc", 1))
        xtb = XTBInterface(
            gfn_level=int(path_cfg.get("gfn_level", 2)),
            solvent=solvent,
            nproc=nproc,
            config=self.config,
        )

        path_dir = output_dir / "xtb_path" / str(path_name)
        path_dir.mkdir(parents=True, exist_ok=True)
        batch_name = f"{self.molecule_name or 'product'}:xtb_path:{path_name}"

        self._emit_progress(
            "batch_started",
            batch=batch_name,
            phase="xtb_path",
            label=f"xTB2 meta-dynamics PATH search ({path_name})",
            started_at=time.time(),
            engine="xtb",
            method=f"GFN{int(path_cfg.get('gfn_level', 2))}",
            npoint=npoint,
        )
        path_started = time.monotonic()
        try:
            path_result = xtb.path(
                start_xyz=Path(start_xyz),
                end_xyz=Path(end_xyz),
                output_dir=path_dir,
                nrun=int(path_cfg.get("nrun", 1)),
                npoint=npoint,
                anopt=int(path_cfg.get("anopt", 10)),
                kpush=float(path_cfg.get("kpush", 0.003)),
                kpull=float(path_cfg.get("kpull", -0.015)),
                ppull=float(path_cfg.get("ppull", 0.05)),
                alp=float(path_cfg.get("alp", 1.2)),
                charge=charge,
                spin=spin,
            )
        except Exception as exc:
            self._emit_progress(
                "batch_finished",
                batch=batch_name,
                phase="xtb_path",
                status="failed",
                error=str(exc),
                elapsed_seconds=time.monotonic() - path_started,
            )
            raise RuntimeError(f"S2 xTB PATH execution failed: {exc}") from exc

        if not path_result.success or not path_result.path_xyz_files:
            self._emit_progress(
                "batch_finished",
                batch=batch_name,
                phase="xtb_path",
                status="failed",
                error=path_result.error_message or "no path frames",
                elapsed_seconds=time.monotonic() - path_started,
            )
            raise RuntimeError(
                f"S2 xTB PATH produced no frames: {path_result.error_message or 'unknown'}"
            )

        frame_paths = tuple(Path(p) for p in path_result.path_xyz_files)
        path_arclength = self._path_arclength(frame_paths)
        target_distances = self._path_mean_forming_bond_distances(frame_paths, forming_bonds)
        sp_dir = path_dir / "frame_sp"
        sp_energies = self._path_frame_xtb_sp(frame_paths, sp_dir, charge, spin)

        valid_indices = [i for i, e in enumerate(sp_energies) if e is not None]
        min_required = max(3, len(frame_paths) // 2)
        if len(valid_indices) < min_required:
            raise RuntimeError(
                f"S2 xTB PATH produced too few frames with energies: "
                f"{len(valid_indices)} valid, {len(frame_paths)} total "
                f"(min required: {min_required})"
            )
        if len(valid_indices) < len(frame_paths):
            self.logger.warning(
                "[S2] xTB PATH SP: %d/%d frames have valid energies; "
                "check individual frame SP logs",
                len(valid_indices), len(frame_paths),
            )

        coordinates = tuple(
            float(d) if d is not None else float("nan")
            for d in target_distances
        )
        energies = tuple(float(e) if e is not None else float("nan") for e in sp_energies)

        trajectory_quality: Dict[str, Any] = {
            "xtb_path": True,
            "path_name": str(path_name),
            "npoint_requested": npoint,
            "npoint_returned": len(frame_paths),
            "barrier_forward_kcal": path_result.barrier_forward_kcal,
            "barrier_backward_kcal": path_result.barrier_backward_kcal,
            "reaction_energy_kcal": path_result.reaction_energy_kcal,
            "gradient_norm_at_ts": path_result.gradient_norm_at_ts,
        }
        off_path_indices: Tuple[int, ...] = ()
        if (
            topology_config is not None
            and reference_coords is not None
            and reference_symbols is not None
            and reference_path is not None
        ):
            trajectory_quality, off_path, _ = self._assess_trajectory(
                frame_paths,
                forming_bonds,
                topology_config,
                consecutive_off_path,
                [float(value) for value in energies],
                reference_coords,
                reference_symbols,
                Path(reference_path),
            )
            trajectory_quality.update(
                {
                    "xtb_path": True,
                    "path_name": str(path_name),
                    "npoint_requested": npoint,
                    "npoint_returned": len(frame_paths),
                    "barrier_forward_kcal": path_result.barrier_forward_kcal,
                    "barrier_backward_kcal": path_result.barrier_backward_kcal,
                    "reaction_energy_kcal": path_result.reaction_energy_kcal,
                    "gradient_norm_at_ts": path_result.gradient_norm_at_ts,
                }
            )
            off_path_indices = tuple(sorted(off_path))

        attempt = ScanAttempt(
            attempt_id=f"NN_xtb_path_{path_name}",
            kind="xtb_path",
            directory=path_dir,
            frame_paths=frame_paths,
            target_coordinates_A=coordinates,
            xtb_energies_hartree=energies,
            off_path_indices=off_path_indices,
            trajectory_quality=trajectory_quality,
            scan_policy="xtb_path",
        )
        path_metadata = {
            "path_name": str(path_name),
            "ts_guess_xyz": str(path_result.ts_guess_xyz) if path_result.ts_guess_xyz else None,
            "path_log": str(path_result.path_log) if path_result.path_log else None,
            "estimated_ts_point": path_result.estimated_ts_point,
            "barrier_forward_kcal": path_result.barrier_forward_kcal,
            "barrier_backward_kcal": path_result.barrier_backward_kcal,
            "reaction_energy_kcal": path_result.reaction_energy_kcal,
            "gradient_norm_at_ts": path_result.gradient_norm_at_ts,
            "start_xyz": str(start_xyz),
            "end_xyz": str(end_xyz),
            "path_arclength": [float(v) for v in path_arclength],
        }
        self._emit_progress(
            "batch_finished",
            batch=batch_name,
            phase="xtb_path",
            status="complete",
            npoint_returned=len(frame_paths),
            elapsed_seconds=time.monotonic() - path_started,
        )
        return attempt, path_metadata

    @staticmethod
    def _product_connected_valid_indices(
        frame_count: int,
        off_path_indices: Sequence[int],
    ) -> List[int]:
        """Return the topology-valid PATH segment connected to the product end.

        xTB PATH frames are ordered from the stretched endpoint to the product.
        A distorted endpoint may legitimately be retained for exploration, but
        selection must never jump across an off-topology section to use an
        earlier, disconnected fragment of that PATH.
        """
        invalid = {int(index) for index in off_path_indices}
        end = int(frame_count) - 1
        if end < 0 or end in invalid:
            return []
        start = end
        while start - 1 >= 0 and start - 1 not in invalid:
            start -= 1
        return list(range(start, end + 1))

    @classmethod
    def _path_selection_anchors(
        cls,
        *,
        energies: Sequence[Optional[float]],
        reaction_coordinate: Sequence[float],
        valid_indices: Sequence[int],
        off_path_indices: Sequence[int],
    ) -> Dict[str, Any]:
        """Build safe PATH-local anchors for a product-connected valid segment."""
        if not valid_indices:
            raise RuntimeError("S2 PATH has no topology-valid segment connected to product")
        valid = [int(index) for index in valid_indices]
        energy_valid = [
            index for index in valid if energies[index] is not None
        ]
        if not energy_valid:
            raise RuntimeError("S2 PATH product-connected segment has no energies")
        maximum = max(
            energy_valid,
            key=lambda index: cls._require_energy(
                energies[index], "S2 PATH anchor energy is missing"
            ),
        )
        return {
            "product_index": valid[-1],
            "coarse_ts_index": maximum,
            "plateau_onset_index": valid[0],
            "absolute_energy_maximum_index": maximum,
            "topology_drift_index": valid[0] - 1 if valid[0] > 0 else None,
            "last_valid_before_drift_index": valid[0],
            "scan_endpoint_index": 0,
            "product_connected_valid_indices": valid,
            "excluded_path_indices": sorted(
                set(range(len(reaction_coordinate))) - set(valid)
                | {int(index) for index in off_path_indices}
            ),
        }

    @staticmethod
    def _branch_relative_energies(
        energies: Sequence[Optional[float]],
        *,
        reference_side: str = "first",
    ) -> List[Optional[float]]:
        if reference_side == "last":
            reference = next(
                (float(value) for value in reversed(energies) if value is not None),
                None,
            )
        else:
            reference = next((float(value) for value in energies if value is not None), None)
        return [
            None if value is None or reference is None
            else (float(value) - reference) * HARTREE_TO_KCAL
            for value in energies
        ]

    @staticmethod
    def _mean_frame_coordinates(
        forming_bond_distances: Sequence[Optional[Sequence[float]]],
    ) -> List[float]:
        coordinates: List[float] = []
        for index, distances in enumerate(forming_bond_distances):
            if distances and all(value is not None for value in distances):
                coordinates.append(
                    float(sum(float(value) for value in distances) / len(distances))
                )
            else:
                coordinates.append(float(index))
        return coordinates

    @staticmethod
    def _seed_frame_index(seed: Optional[Mapping[str, Any]]) -> Optional[int]:
        if not seed:
            return None
        frame_index_raw = seed.get("frame_index")
        if frame_index_raw is None:
            return None
        try:
            frame_index = int(frame_index_raw)
        except (TypeError, ValueError):
            return None
        return frame_index if frame_index >= 0 else None

    @staticmethod
    def _selector_ts_rule(seed_evidence: str, fallback: str = "unresolved") -> str:
        if seed_evidence in {"knee_shifted", "peak_knee_shifted"}:
            return "energy_knee_right_shifted"
        if seed_evidence in {"local_peak", "peak_and_basin"}:
            return "refined_curve_local_maximum"
        if seed_evidence == "monotonic_shoulder":
            return "late_pre_ts_low_gradient_platform"
        return fallback

    @staticmethod
    def _selector_int_rule(
        *,
        selection_source: str,
        seed_evidence: str,
        int_seed: Optional[Mapping[str, Any]],
        has_independent_int: bool,
        fallback: str = "no_int_search_seed",
    ) -> str:
        if not int_seed:
            return fallback
        if has_independent_int:
            return (
                "b973c_relaxed_scan_post_ts_basin"
                if selection_source == "orca_relaxed_scan"
                else "refined_curve_pre_ts_local_minimum"
            )
        if str(int_seed.get("selection_mode")) == "ts_to_effective_endpoint_midpoint":
            return "ts_to_effective_endpoint_midpoint"
        if str(int_seed.get("selection_mode")) == "stretch_plateau":
            return "stretch_side_low_energy_plateau"
        if seed_evidence in {"knee_shifted", "peak_knee_shifted"}:
            return "ts_to_effective_endpoint_midpoint"
        if seed_evidence == "monotonic_shoulder":
            return "late_pre_ts_low_gradient_platform"
        if bool(int_seed.get("shared_with_ts", False)):
            return "no_resolved_pre_ts_basin_shared_ts_seed"
        return "search_seed"

    @staticmethod
    def _selector_int_selection_mode(
        *,
        seed_evidence: str,
        int_seed: Optional[Mapping[str, Any]],
        has_independent_int: bool,
    ) -> str:
        if not int_seed:
            return "unavailable"
        if has_independent_int:
            return "stable_basin_candidate"
        if str(int_seed.get("selection_mode")) == "ts_to_effective_endpoint_midpoint":
            return "ts_to_effective_endpoint_midpoint"
        if str(int_seed.get("selection_mode")) == "stretch_plateau":
            return "stretch_side_low_energy_plateau"
        if seed_evidence in {"knee_shifted", "peak_knee_shifted"}:
            return "ts_to_effective_endpoint_midpoint"
        if seed_evidence == "monotonic_shoulder":
            return "late_pre_ts_platform_fallback"
        if bool(int_seed.get("shared_with_ts", False)):
            return "shared_ts_fallback"
        return "search_seed"


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

    def _execute_scan(
        self,
        start_xyz: Path,
        output_dir: Path,
        bonds: Tuple[Tuple[int, int], ...],
        params: Dict[str, Any],
        direction: str,
        charge: int = 0,
        spin: int = 1,
    ) -> Tuple[Any, List[float]]:
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

        scan_constraints, policy_force_constant = selector.select_policy(
            bonds, scan_policy_name, start_dist
        )
        force_constant = (
            float(params["scan_force_constant"])
            if params.get("scan_force_constant") is not None
            else float(policy_force_constant)
        )

        xtb_settings = self.step2_cfg.get("xtb_settings", {}) or {}
        solvent = str(xtb_settings.get("solvent", self.config.get("theory", {}).get("optimization", {}).get("solvent", "acetone")))
        scan_cfg = self.step2_cfg.get("scan", {}) or {}
        nproc = int(scan_cfg.get("nproc") or xtb_settings.get("nproc") or self.config.get("resources", {}).get("nproc", 1))

        xtb = XTBInterface(solvent=solvent, nproc=nproc, config=self.config)
        result = xtb.scan(
            xyz_file=start_xyz,
            output_dir=output_dir,
            constraints=scan_constraints,
            fixed_constraints=dict(params.get("fixed_constraints", {}) or {}),
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

        return result, energies_local

    @staticmethod
    def _persistent_off_path_start(
        off_path_indices: Sequence[int],
        total_frames: int,
        consecutive: int,
    ) -> Optional[int]:
        required = max(1, int(consecutive))
        off_path = {int(index) for index in off_path_indices}
        for start in range(total_frames):
            if all((start + offset) in off_path for offset in range(required)):
                return start
        return None

    def _assess_trajectory(
        self,
        frame_paths: Sequence[Path],
        bonds: Sequence[Tuple[int, int]],
        topology_config: Dict[str, Any],
        consecutive_off_path: int,
        energies: Sequence[float],
        reference_coords: Sequence[Sequence[float]],
        reference_symbols: Sequence[str],
        reference_path: Path,
    ) -> Tuple[Dict[str, Any], set[int], Optional[int]]:
        quality: Dict[str, Any] = {
            "checked": False,
            "total_frames": len(frame_paths),
            "off_path_indices": [],
            "off_path_count": 0,
            "frame_issues": [],
        }
        if bool(topology_config.get("enabled", True)):
            quality = check_scan_trajectory(
                product_coords=np.asarray(reference_coords, dtype=float),
                symbols=list(reference_symbols),
                forming_bonds=bonds,
                frame_paths=frame_paths,
                graph_scale=float(topology_config.get("graph_scale", 1.25)),
            )
            quality["reference_frame"] = str(reference_path)
        off_path = {int(index) for index in quality.get("off_path_indices", [])}
        raw_off_path = set(off_path)
        isolated_suspects = [
            index
            for index in sorted(raw_off_path)
            if 0 < index < len(frame_paths) - 1
            and index - 1 not in raw_off_path
            and index + 1 not in raw_off_path
        ]
        reclassified: List[int] = []
        isolated_diagnostics: List[Dict[str, Any]] = []
        if bool(topology_config.get("reclassify_isolated_suspects", False)):
            jump_limit = float(topology_config.get("isolated_max_energy_jump_kcal", 5.0))
            core_rmsd_limit = float(
                topology_config.get("isolated_max_reaction_core_rmsd_A", 0.15)
            )
            core_builder = CompositeProfileBuilder(
                forming_bonds=bonds,
                reaction_core_depth=2,
            )
            core_indices = core_builder._reaction_core_indices(reference_path)
            for index in isolated_suspects:
                left_jump = (
                    abs(float(energies[index]) - float(energies[index - 1]))
                    * HARTREE_TO_KCAL
                )
                right_jump = (
                    abs(float(energies[index + 1]) - float(energies[index]))
                    * HARTREE_TO_KCAL
                )
                try:
                    left_core_rmsd = core_builder._aligned_rmsd(
                        frame_paths[index - 1], frame_paths[index], core_indices
                    )
                    right_core_rmsd = core_builder._aligned_rmsd(
                        frame_paths[index], frame_paths[index + 1], core_indices
                    )
                    bridge_core_rmsd = core_builder._aligned_rmsd(
                        frame_paths[index - 1], frame_paths[index + 1], core_indices
                    )
                except (OSError, ValueError, np.linalg.LinAlgError):
                    left_core_rmsd = right_core_rmsd = bridge_core_rmsd = float("inf")
                accepted = (
                    max(left_jump, right_jump) <= jump_limit
                    and max(left_core_rmsd, right_core_rmsd, bridge_core_rmsd)
                    <= core_rmsd_limit
                )
                isolated_diagnostics.append(
                    {
                        "frame_index": index,
                        "left_energy_jump_kcal": left_jump,
                        "right_energy_jump_kcal": right_jump,
                        "left_reaction_core_rmsd_A": left_core_rmsd,
                        "right_reaction_core_rmsd_A": right_core_rmsd,
                        "bridge_reaction_core_rmsd_A": bridge_core_rmsd,
                        "reclassified": accepted,
                    }
                )
                if accepted:
                    reclassified.append(index)
            off_path -= set(reclassified)
        quality["raw_off_path_indices"] = sorted(raw_off_path)
        quality["isolated_suspect_indices"] = isolated_suspects
        quality["reclassified_isolated_indices"] = reclassified
        quality["isolated_suspect_diagnostics"] = isolated_diagnostics
        quality["off_path_indices"] = sorted(off_path)
        quality["off_path_count"] = len(off_path)
        persistent_start = self._persistent_off_path_start(
            sorted(off_path), len(frame_paths), consecutive_off_path
        )
        quality["persistent_off_path_start"] = persistent_start
        quality["usable_end_index"] = (
            persistent_start - 1 if persistent_start is not None else len(frame_paths) - 1
        )
        return quality, off_path, persistent_start

    def _make_scan_attempt(
        self,
        *,
        attempt_id: str,
        kind: str,
        directory: Path,
        scan_result: Any,
        energies: Sequence[float],
        params: Dict[str, Any],
        bonds: Sequence[Tuple[int, int]],
        topology_config: Dict[str, Any],
        consecutive_off_path: int,
        reference_coords: Sequence[Sequence[float]],
        reference_symbols: Sequence[str],
        reference_path: Path,
        parent_attempt_id: Optional[str] = None,
        seed_xyz: Optional[Path] = None,
        seed_source_attempt: Optional[str] = None,
        seed_source_frame_index: Optional[int] = None,
        seed_coordinate_A: Optional[float] = None,
        direction: str = "outward",
    ) -> Tuple[ScanAttempt, set[int], Optional[int]]:
        frame_paths = tuple(
            Path(path) for path in (getattr(scan_result, "geometries", None) or [])
        )
        if len(frame_paths) != len(energies):
            raise RuntimeError(
                "S2 scan trajectory is incomplete: "
                f"{len(energies)} energies but {len(frame_paths)} geometry frames"
            )
        coordinates = tuple(
            compute_scan_distances(
                float(params["scan_start_distance"]),
                float(params["scan_end_distance"]),
                len(energies),
                direction=direction,
            )
        )
        quality, off_path, persistent_start = self._assess_trajectory(
            frame_paths,
            bonds,
            topology_config,
            consecutive_off_path,
            energies,
            reference_coords,
            reference_symbols,
            reference_path,
        )
        quality["scan_direction"] = direction
        attempt = ScanAttempt(
            attempt_id=attempt_id,
            kind=kind,
            parent_attempt_id=parent_attempt_id,
            seed_xyz=Path(seed_xyz) if seed_xyz is not None else None,
            seed_source_attempt=seed_source_attempt,
            seed_source_frame_index=seed_source_frame_index,
            seed_coordinate_A=seed_coordinate_A,
            directory=Path(directory),
            frame_paths=frame_paths,
            target_coordinates_A=coordinates,
            xtb_energies_hartree=tuple(float(value) for value in energies),
            off_path_indices=tuple(sorted(off_path)),
            trajectory_quality=quality,
            scan_policy=str(params.get("scan_policy", "policy_c")),
            fixed_constraints=dict(params.get("fixed_constraints", {}) or {}),
        )
        return attempt, off_path, persistent_start

    @staticmethod
    def _parent_seed(
        parent: ScanAttempt,
        target_coordinate_A: float,
    ) -> Tuple[Path, int, float]:
        """Choose a topology-valid parent frame nearest the new scan start."""

        off_path = set(parent.off_path_indices)
        eligible = [
            index for index in range(len(parent.frame_paths)) if index not in off_path
        ]
        if not eligible:
            raise RuntimeError(
                f"S2 attempt {parent.attempt_id} has no topology-valid frame for warm start"
            )
        index = min(
            eligible,
            key=lambda item: abs(
                float(parent.target_coordinates_A[item]) - float(target_coordinate_A)
            ),
        )
        return (
            Path(parent.frame_paths[index]),
            int(index),
            float(parent.target_coordinates_A[index]),
        )

    @classmethod
    def _attempt_needs_endpoint_extension(
        cls,
        attempt: ScanAttempt,
        extension_config: Dict[str, Any],
        anchor_config: Optional[Dict[str, Any]] = None,
    ) -> Tuple[bool, Dict[str, Any]]:
        anchors = cls._detect_scan_anchors(
            attempt.xtb_energies_hartree,
            attempt.target_coordinates_A,
            attempt.off_path_indices,
            int((anchor_config or {}).get("persistent_drift_points", 2)),
            anchor_config,
        )
        plateau_index = int(anchors["plateau_onset_index"])
        plateau_coordinate = float(attempt.target_coordinates_A[plateau_index])
        coordinate_max = max(attempt.target_coordinates_A)
        last_valid_index = int(anchors["last_valid_before_drift_index"])
        valid_after_plateau = max(0, last_valid_index - plateau_index)
        required_span = float(extension_config.get("min_dissociation_side_span_A", 0.50))
        required_points = max(1, int(extension_config.get("min_valid_post_ts_points", 5)))
        available_span = float(attempt.target_coordinates_A[last_valid_index]) - plateau_coordinate
        reasons: List[str] = []
        if anchors["topology_drift_index"] is None:
            if plateau_index == int(anchors["scan_endpoint_index"]):
                reasons.append("plateau_onset_at_scan_endpoint")
            if available_span < required_span:
                reasons.append("insufficient_plateau_endpoint_span")
            if valid_after_plateau < required_points:
                reasons.append("insufficient_valid_plateau_endpoint_points")
        return bool(reasons), {
            "plateau_onset_coordinate_A": plateau_coordinate,
            "coordinate_max_A": coordinate_max,
            "available_plateau_boundary_span_A": available_span,
            "required_plateau_boundary_span_A": required_span,
            "valid_plateau_boundary_points": valid_after_plateau,
            "required_valid_plateau_boundary_points": required_points,
            "topology_drift_index": anchors["topology_drift_index"],
            "reasons": reasons,
        }

    @staticmethod
    def _unified_selection_records(
        *,
        profile: Any,
        selection: Any,
        method_energies: Sequence[Optional[float]],
        energy_source: str,
    ) -> Tuple[Dict[str, Any], Dict[str, Any], Optional[int], Optional[int], str]:
        """Adapt the shared selector result to the S2 manifest contract.

        This is the only adapter between ``SeedSelection`` and the historical
        ``selections.ts_guess``/``selections.intermediate`` records.  The
        adapter never makes a second chemical decision; it only enriches the
        selector's canonical frame indices with output metadata.
        """

        diagnostics = dict(getattr(selection, "diagnostics", {}) or {})
        ts_seed = dict(getattr(selection, "ts_search_seed", None) or {})
        int_seed = dict(getattr(selection, "int_search_seed", None) or {})
        ts_index_raw = ts_seed.get("frame_index")
        ts_index = None if ts_index_raw is None else int(ts_index_raw)
        peak_index_raw = diagnostics.get("energy_peak_index", ts_index)
        peak_index = None if peak_index_raw is None else int(peak_index_raw)
        source = str(profile.source)
        source_attempt = profile.source_provenance.get("attempt_id")
        state = str(getattr(selection, "s2_state", "unresolved"))
        evidence = str(getattr(selection, "seed_evidence", "none"))
        rejection_reason = getattr(selection, "rejection_reason", None)

        if ts_index is None:
            ts_selection: Dict[str, Any] = {
                "index": None,
                "energy_peak_index": peak_index,
                "seed_index": None,
                "frame_xyz": None,
                "rule": "unified_selector_unresolved",
                "actual_method": "unified_selector",
                "configured_method": "unified_selector",
                "selection_status": "unavailable",
                "reason": rejection_reason,
                "candidate_indices": [],
                "source": None,
            }
        else:
            frame = profile.frames[ts_index]
            confidence = str(ts_seed.get("confidence", "medium"))
            ts_rule = (
                "unified_selector_knee_shifted_ts"
                if evidence in {"knee_shifted", "peak_knee_shifted"}
                else
                "unified_selector_peak"
                if evidence in {"local_peak", "peak_and_basin"}
                else "unified_selector_monotonic_shoulder"
                if evidence == "monotonic_shoulder"
                else "unified_selector"
            )
            ts_selection = {
                "index": ts_index,
                "energy_peak_index": peak_index,
                "seed_index": ts_index,
                "frame_xyz": str(frame.xyz),
                "rule": ts_rule,
                "actual_method": "unified_selector",
                "configured_method": "unified_selector",
                "confidence": confidence,
                "candidate_indices": [ts_index],
                "selection_status": "selected",
                "energy_source": energy_source,
                "source": source,
                "source_attempt": source_attempt,
                "energy_hartree": (
                    None
                    if ts_index >= len(method_energies)
                    or method_energies[ts_index] is None
                    else float(method_energies[ts_index])
                ),
            }

        int_index_raw = int_seed.get("frame_index")
        int_index = None if int_index_raw is None else int(int_index_raw)
        shared = bool(int_seed.get("shared_with_ts", False))
        has_independent_int = bool(getattr(selection, "has_independent_int", False))
        if int_index is None:
            int_selection: Dict[str, Any] = {
                "index": None,
                "frame_xyz": None,
                "rule": "unified_selector_no_int_seed",
                "selection_mode": "unavailable",
                "selection_status": "unavailable",
                "stationary_point_claimed": False,
                "reason": "no_int_search_seed",
                "candidate_indices": [],
                "source": None,
            }
            int_status = "unavailable"
        else:
            int_frame = profile.frames[int_index]
            int_mode = "stable_basin_candidate" if has_independent_int else "shared_ts_fallback"
            int_selection = {
                "index": int_index,
                "frame_xyz": str(int_frame.xyz),
                "rule": "unified_selector_pre_peak_basin"
                if has_independent_int
                else "unified_selector_ts_to_effective_endpoint_midpoint"
                if str(int_seed.get("selection_mode"))
                == "ts_to_effective_endpoint_midpoint"
                else "unified_selector_shared_ts",
                "selection_mode": (
                    str(int_seed.get("selection_mode"))
                    if int_seed.get("selection_mode")
                    else int_mode
                ),
                "selection_status": "selected",
                "stationary_point_claimed": False,
                "reason": None if has_independent_int else "no_resolved_pre_ts_basin",
                "candidate_indices": [int_index],
                "shared_ts_index": int(ts_index) if shared and ts_index is not None else None,
                "energy_source": energy_source,
                "source": source,
                "source_attempt": source_attempt,
                "energy_hartree": (
                    None
                    if int_index >= len(method_energies)
                    or method_energies[int_index] is None
                    else float(method_energies[int_index])
                ),
            }
            int_status = "shared_with_ts" if shared else "selected"

        return ts_selection, int_selection, peak_index, int_index, int_status

    @staticmethod
    def _validate_unified_selection_contract(
        *,
        s2_state: str,
        selection_source: Optional[str],
        ts_selection: Mapping[str, Any],
        intermediate_selection: Mapping[str, Any],
        s3_dispatch: Mapping[str, Any],
        ts_guess_xyz: Path,
        intermediate_xyz: Path,
    ) -> None:
        """Reject mixed selector/legacy state before publishing S2 artifacts."""

        state = str(s2_state or "unresolved")
        resolution = str(s3_dispatch.get("resolution") or "unresolved")
        submit_ts = bool(s3_dispatch.get("submit_ts", False))
        submit_intermediate = bool(s3_dispatch.get("submit_intermediate", False))
        ts_index = ts_selection.get("index")
        int_index = intermediate_selection.get("index")

        if state != resolution:
            raise RuntimeError(
                f"S2 selection contract mismatch: s2_state={state!r}, "
                f"s3_dispatch.resolution={resolution!r}"
            )
        if state == "unresolved":
            if submit_ts or submit_intermediate:
                raise RuntimeError("Unresolved S2 selection cannot submit S3 jobs")
            if ts_index is not None or int_index is not None:
                raise RuntimeError("Unresolved S2 selection retains an active seed index")
            if ts_guess_xyz.exists() or intermediate_xyz.exists():
                raise RuntimeError("Unresolved S2 selection retains a canonical XYZ artifact")
            if selection_source is not None:
                raise RuntimeError("Unresolved S2 selection must not have an active source")
            return

        if state not in {"path_seeded", "rescue_seeded"}:
            raise RuntimeError(f"Unknown S2 selection state: {state!r}")
        if selection_source is None or ts_index is None or not submit_ts:
            raise RuntimeError("Seeded S2 selection lacks a canonical TS seed")
        if not ts_guess_xyz.exists():
            raise RuntimeError("Seeded S2 selection lacks ts_guess.xyz")
        if submit_intermediate:
            if int_index is None or not intermediate_xyz.exists():
                raise RuntimeError("S3 intermediate dispatch lacks intermediate.xyz")

    @staticmethod
    def _relative_scan_energies(
        energies: Sequence[Optional[float]],
        *,
        reference: str = "last",
    ) -> List[Optional[float]]:
        values = [None if value is None else float(value) for value in energies]
        valid = [index for index, value in enumerate(values) if value is not None]
        if not valid:
            return [None] * len(values)
        reference_index = valid[-1] if reference == "last" else valid[0]
        reference_energy = float(values[reference_index])
        return [
            None if value is None else (float(value) - reference_energy) * HARTREE_TO_KCAL
            for value in values
        ]

    @staticmethod
    def _scan_coordinates(
        frames: Sequence[Path],
        forming_bonds: Sequence[Tuple[int, int]],
    ) -> List[float]:
        coordinates: List[float] = []
        for index, frame in enumerate(frames):
            try:
                coords, _symbols = read_xyz(Path(frame))
                distances = [
                    float(np.linalg.norm(coords[int(atom_i)] - coords[int(atom_j)]))
                    for atom_i, atom_j in forming_bonds
                ]
                coordinates.append(float(np.mean(distances)))
            except (OSError, ValueError, IndexError):
                coordinates.append(float(index))
        return coordinates

    @staticmethod
    def _topology_rescue_decision(
        profile: Optional[Any],
        selection: Optional[Any],
        minimum_clean_frames_after_knee: int,
    ) -> Dict[str, Any]:
        """Decide whether topology distortion truncates the usable curve.

        The rescue boundary is defined by the curve itself: the first
        topology-invalid frame must occur after a supported knee, and the
        clean corridor must contain enough frames after that knee for seed
        optimization.  No absolute reaction-coordinate threshold is used.
        """

        excluded = sorted(
            int(index) for index in (getattr(profile, "excluded_frames", ()) or ())
        )
        minimum_frames = max(1, int(minimum_clean_frames_after_knee))
        diagnostics = getattr(selection, "diagnostics", {}) if selection is not None else {}
        if not isinstance(diagnostics, Mapping):
            diagnostics = {}
        knee_value = diagnostics.get("knee_frame_index")
        try:
            knee_index = None if knee_value is None else int(knee_value)
        except (TypeError, ValueError):
            knee_index = None

        if not excluded:
            return {
                "distortion_frame_index": None,
                "knee_frame_index": knee_index,
                "post_knee_frame_indices": [],
                "clean_frames_after_knee": 0,
                "minimum_clean_frames_after_knee": minimum_frames,
                "rescue_required": False,
                "decision": "retain_primary_no_distortion",
            }

        distortion_index = excluded[0]
        frames = tuple(getattr(profile, "frames", ()) or ())
        post_knee_indices = []
        if knee_index is not None and knee_index < distortion_index:
            post_knee_indices = [
                int(getattr(frame, "frame_index", index))
                for index, frame in enumerate(frames)
                if knee_index < int(getattr(frame, "frame_index", index)) < distortion_index
                and bool(getattr(frame, "topology_valid", True))
            ]

        if knee_index is None:
            decision = "rescue_no_supported_knee"
            rescue_required = True
        elif knee_index >= distortion_index:
            decision = "rescue_knee_not_before_distortion"
            rescue_required = True
        elif len(post_knee_indices) < minimum_frames:
            decision = "rescue_insufficient_post_knee_support"
            rescue_required = True
        else:
            decision = "retain_primary_knee_supported"
            rescue_required = False

        return {
            "distortion_frame_index": distortion_index,
            "knee_frame_index": knee_index,
            "post_knee_frame_indices": post_knee_indices,
            "clean_frames_after_knee": len(post_knee_indices),
            "minimum_clean_frames_after_knee": minimum_frames,
            "rescue_required": rescue_required,
            "decision": decision,
        }

    def _run_orca_gfn2_workflow(
        self,
        product_xyz: Path,
        output_dir: Path,
        forming_bonds: Sequence[Tuple[int, int]],
        scan_config: Optional[Dict[str, Any]],
    ) -> Tuple[
        Optional[Path], Optional[Path], Optional[Path], Tuple[Tuple[int, int], ...],
        Path, str, str, Tuple[str, ...]
    ]:
        """Run the canonical ORCA-GFN2 → B97-3c-SP → rescue workflow."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        archived_outputs = self._archive_previous_outputs(output_dir)
        bonds = self._validate_forming_bonds(forming_bonds, product_xyz_path=product_xyz)
        product_file = self._resolve_product_file(product_xyz)
        params = self._resolve_scan_params(scan_config)
        selection_cfg = dict(params.get("selection", {}) or {})
        policy = policy_from_config(selection_cfg)
        degraded_reasons: List[str] = []

        primary_dir = output_dir / "scan"
        primary = B97CRelaxedScanRescuer(
            self.config,
            event_callback=self.event_callback,
            variant=self.molecule_name,
            scan_method="GFN2-xTB",
            scan_role="primary",
        ).run(
            product_file,
            primary_dir,
            bonds,
            trigger_reasons=["primary_orca_gfn2_relaxed_scan"],
        )
        primary_frames = [Path(value) for value in primary.get("frames") or []]
        primary_energies = [
            None if value is None else float(value)
            for value in primary.get("energies_hartree") or []
        ]
        b973c_refinement: Dict[str, Any] = {"status": "not_run"}
        b973c_energies: List[Optional[float]] = []
        primary_profile = None
        if primary.get("status") == "complete" and len(primary_frames) >= 3:
            refiner = ScanEnergyRefiner(
                self.config,
                output_dir / "energy_refinement",
                event_callback=self.event_callback,
                variant=self.molecule_name or "product",
            )
            b973c_refinement = refiner.refine(
                primary_frames,
                point_ids=[f"gfn2_{index:04d}" for index in range(len(primary_frames))],
            )
            b973c_energies = [
                None if value is None else float(value)
                for value in b973c_refinement.get("energies_hartree") or []
            ]
            b973c_energies = (b973c_energies + [None] * len(primary_frames))[:len(primary_frames)]
            if all(value is not None for value in b973c_energies):
                primary_profile = build_orca_scan_profile(
                    frames=primary_frames,
                    energies_hartree=b973c_energies,
                    forming_bonds=bonds,
                    product_xyz=product_file,
                    energy_source="ORCA B97-3c SP",
                    source_provenance={
                        "scan_engine": "orca",
                        "scan_method": "GFN2-xTB",
                        "energy_refinement_engine": "orca",
                        "energy_refinement_method": "B97-3c",
                        "gfn2_energy_source": primary.get("energy_source"),
                    },
                )

        primary_selection = (
            select_path_seeds(primary_profile, policy)
            if primary_profile is not None
            else None
        )
        scan_or_sp_incomplete = (
            primary_profile is None
            or len(primary_frames) < 3
            or len(primary_energies) != len(primary_frames)
            or any(value is None for value in primary_energies)
        )
        topology_decision = self._topology_rescue_decision(
            primary_profile,
            primary_selection,
            policy.minimum_clean_frames_after_knee,
        )
        topology_rescue_required = bool(topology_decision["rescue_required"])
        # GFN2-only robustness mode: with the B97-3c rescue disabled, topology
        # drift no longer forces the rescue branch -- the primary GFN2
        # selection (knee-shifted over clean frames) is retained so drift is
        # tolerated instead of triggering a second, expensive scan.
        rescue_cfg = dict((self.config or {}).get("step2", {}).get("rescue", {}) or {})
        rescue_enabled = bool(rescue_cfg.get("enabled", True))
        if topology_rescue_required and not rescue_enabled:
            degraded_reasons.append("orca_gfn2_topology_drift_tolerated_no_rescue")
        rescue_required = scan_or_sp_incomplete or (topology_rescue_required and rescue_enabled)
        rescue_payload: Dict[str, Any] = {
            "status": "not_run",
            "resolution": "unresolved",
            "s2_state": "unresolved",
            "ts_xyz": None,
            "intermediate_xyz": None,
        }
        active_profile = primary_profile
        active_frames = primary_frames
        active_energies = b973c_energies
        selection_source = "orca_gfn2_b973c_sp"
        s2_state = "unresolved"
        selection: Dict[str, Any] = {}
        if rescue_required:
            degraded_reasons.append(
                "orca_gfn2_scan_or_sp_incomplete"
                if scan_or_sp_incomplete
                else "orca_gfn2_topology_before_knee_support"
            )
            active_profile = None
            rescue_payload = B97CRelaxedScanRescuer(
                self.config,
                event_callback=self.event_callback,
                variant=self.molecule_name,
                scan_method="B97-3c",
                scan_role="rescue",
            ).run(
                product_file,
                output_dir / "rescue",
                bonds,
                trigger_reasons=[
                    "orca_gfn2_topology_before_knee_support"
                    if not scan_or_sp_incomplete
                    else "orca_gfn2_scan_or_sp_incomplete"
                ],
            )
            selection_source = "orca_b973c_relaxed_scan"
            s2_state = str(rescue_payload.get("s2_state") or "unresolved")
            active_frames = [Path(value) for value in rescue_payload.get("frames") or []]
            active_energies = [
                None if value is None else float(value)
                for value in rescue_payload.get("energies_hartree") or []
            ]
            if rescue_payload.get("scan_profile"):
                rescue_profile_path = Path(str(rescue_payload["scan_profile"]))
                try:
                    rescue_profile_payload = json.loads(
                        rescue_profile_path.read_text(encoding="utf-8")
                    )
                    selection = dict(rescue_profile_payload.get("rescue") or rescue_payload)
                except (OSError, ValueError):
                    selection = dict(rescue_payload)
            else:
                selection = dict(rescue_payload)
        elif primary_profile is not None:
            selected = primary_selection or select_path_seeds(primary_profile, policy)
            selection = selected.to_dict()
            s2_state = "gfn2_seeded" if selected.s2_state != "unresolved" else "unresolved"
        else:
            selection = {"s2_state": "unresolved", "rejection_reason": "gfn2_scan_or_sp_incomplete"}

        if active_profile is None and active_frames and len(active_frames) == len(active_energies) and all(
            value is not None for value in active_energies
        ):
            active_profile = build_orca_scan_profile(
                frames=active_frames,
                energies_hartree=active_energies,
                forming_bonds=bonds,
                product_xyz=product_file,
                energy_source="ORCA B97-3c relaxed scan",
                source_provenance={
                    "scan_engine": "orca",
                    "scan_method": "B97-3c",
                    "scan_role": "rescue",
                },
            )
        if active_profile is None:
            active_profile = build_orca_scan_profile(
                frames=active_frames,
                energies_hartree=active_energies,
                forming_bonds=bonds,
                product_xyz=product_file,
                energy_source="ORCA B97-3c SP",
            )

        ts_seed = dict(selection.get("ts_search_seed") or {})
        int_seed = dict(selection.get("int_search_seed") or {})
        if s2_state == "unresolved":
            ts_seed = {}
            int_seed = {}
        def seed_path(seed: Mapping[str, Any]) -> Optional[Path]:
            value = seed.get("xyz") or seed.get("frame_xyz")
            if value:
                return Path(str(value))
            index = seed.get("frame_index")
            if index is None or not (0 <= int(index) < len(active_frames)):
                return None
            return active_frames[int(index)]

        ts_source = seed_path(ts_seed)
        int_source = seed_path(int_seed)
        dispatch = {
            "resolution": "direct_seeded" if s2_state == "gfn2_seeded" else "rescue_seeded" if s2_state == "rescue_seeded" else "unresolved",
            "submit_ts": bool(ts_source and s2_state != "unresolved"),
            "submit_intermediate": bool(int_source and s2_state != "unresolved"),
            "neb_eligible": False,
            "source": selection_source,
            "path_fully_distorted": bool(rescue_required),
        }
        ts_guess = output_dir / "ts_guess.xyz"
        intermediate = output_dir / INTERMEDIATE_XYZ
        if dispatch["submit_ts"]:
            self._atomic_copy(ts_source, ts_guess)
        if dispatch["submit_intermediate"]:
            self._atomic_copy(int_source, intermediate)

        active_gfn2 = primary_energies if active_profile is primary_profile else []
        active_b973c = active_energies
        coordinates = self._scan_coordinates(active_frames, bonds)
        relative_b973c = self._relative_scan_energies(active_b973c)
        relative_gfn2 = self._relative_scan_energies(primary_energies)
        selection_payload = {
            "authority": "unified_selector",
            "source": selection_source,
            "s2_state": s2_state,
            "seed_evidence": selection.get("seed_evidence", "none"),
            "ts_search_seed": ts_seed or None,
            "int_search_seed": int_seed or None,
            "has_independent_int": bool(selection.get("has_independent_int", False)),
            "rejection_reason": selection.get("rejection_reason"),
            "diagnostics": dict(selection.get("diagnostics") or selection.get("selection_diagnostics") or {}),
        }
        ts_index = ts_seed.get("frame_index")
        int_index = int_seed.get("frame_index")
        profile_payload: Dict[str, Any] = {
            "profile_schema_version": "s2_scan_profile_v11",
            "generation_method": "orca_gfn2_relaxed_scan_b973c_sp",
            "product_xyz": str(product_file),
            "forming_bonds": [list(pair) for pair in bonds],
            "scan_engine": "orca",
            "scan_method": "GFN2-xTB",
            "energy_refinement_engine": "orca",
            "energy_refinement_method": "B97-3c",
            "selection_source": selection_source,
            "s2_state": s2_state,
            "seed_evidence": selection.get("seed_evidence", "none"),
            "endpoint_evidence": dict(selection.get("endpoint_evidence") or {}),
            "knee_evidence": dict(selection.get("knee_evidence") or {}),
            "ts_search_seed": ts_seed or None,
            "int_search_seed": int_seed or None,
            "has_independent_int": bool(selection.get("has_independent_int", False)),
            "rejection_reason": selection.get("rejection_reason"),
            "selection_decision": {
                "rule": "knee + adaptive right shift; INT = TS-to-effective-endpoint midpoint",
                "selected_branch": selection_source,
                "ts_seed_index": ts_index,
                "int_seed_index": int_index,
            },
            "selection_policy": {
                "preferred_source": "b973c",
                "actual_source": "B97-3c",
                "ts_source": "B97-3c",
                "intermediate_source": "B97-3c",
                "algorithm": "endpoint_knee_shift_midpoint_v1",
                "topology_rescue_rule": (
                    "supported_knee_before_first_distortion_plus_clean_post_knee_frames"
                ),
                "minimum_clean_frames_after_knee": policy.minimum_clean_frames_after_knee,
            },
            "seed_selection": selection_payload,
            "s3_dispatch": dispatch,
            "scan_parameters": params,
            "trajectory_quality": {
                "checked": True,
                "topology_state": (
                    "incomplete"
                    if scan_or_sp_incomplete
                    else "distorted_before_knee_support"
                    if topology_rescue_required
                    else "tail_distorted_after_knee"
                    if primary_profile is not None and primary_profile.excluded_frames
                    else "valid"
                ),
                "off_path_indices": list(active_profile.excluded_frames),
                "excluded_frames": list(active_profile.excluded_frames),
                "topology_rescue_decision": topology_decision,
            },
            "reaction_coordinate_angstrom": coordinates,
            "energy_curves": {
                "gfn2": {
                    "status": "complete" if primary_energies else "not_run",
                    "engine": "orca",
                    "method": "GFN2-xTB",
                    "energies_hartree": primary_energies,
                    "relative_energies_kcal_mol": relative_gfn2,
                },
                "b973c": {
                    "status": "complete" if all(value is not None for value in active_b973c) else "incomplete",
                    "engine": "orca",
                    "method": "B97-3c",
                    "energies_hartree": active_b973c,
                    "relative_energies_kcal_mol": relative_b973c,
                },
            },
            "energies_hartree": active_b973c,
            "relative_energies_kcal_mol": relative_b973c,
            "frames": [str(frame) for frame in active_frames],
            "energy_refinement": b973c_refinement,
            "primary_scan": dict(primary),
            "rescue": dict(rescue_payload),
            "selections": {
                "ts_guess": {
                    **ts_seed,
                    "index": ts_index,
                    "frame_xyz": str(ts_source) if ts_source else None,
                    "rule": "knee + adaptive right shift",
                    "source": selection_source,
                },
                "intermediate": {
                    **int_seed,
                    "index": int_index,
                    "frame_xyz": str(int_source) if int_source else None,
                    "rule": "TS-to-effective-endpoint midpoint",
                    "source": selection_source,
                },
            },
            "artifact_lifecycle": {
                "stale_outputs_archived_to": str(archived_outputs) if archived_outputs else None,
                "legacy_xtb_scan_path_removed": True,
            },
            "scan_plot": None,
        }
        scan_profile_json = output_dir / "scan_profile.json"
        self.last_profile_payload = profile_payload
        try:
            plot_path = plot_scan_profile(
                scan_profile_json,
                output_dir / "scan_profile.png",
                profile_payload=profile_payload,
            )
            profile_payload["scan_plot"] = str(plot_path) if plot_path else None
        except Exception as exc:
            self.logger.warning("[S2] Failed to render ORCA scan profile: %s", exc, exc_info=True)
        write_text_atomic(scan_profile_json, json.dumps(profile_payload, indent=2), encoding="utf-8")

        status = "COMPLETE" if dispatch["submit_ts"] else "DEGRADED"
        if not dispatch["submit_ts"]:
            degraded_reasons.append(str(selection.get("rejection_reason") or "s2_unresolved"))
        confidence = str(ts_seed.get("confidence") or "low" if dispatch["submit_ts"] else "unresolved")
        return (
            ts_guess if dispatch["submit_ts"] else None,
            intermediate if dispatch["submit_intermediate"] else None,
            intermediate if dispatch["submit_intermediate"] else None,
            bonds,
            scan_profile_json,
            status,
            confidence,
            tuple(dict.fromkeys(degraded_reasons)),
        )

    def run(
        self,
        product_xyz: Path,
        output_dir: Path,
        forming_bonds: Sequence[Tuple[int, int]],
        scan_config: Optional[Dict[str, Any]] = None,
    ) -> Tuple[
        Optional[Path],
        Optional[Path],
        Optional[Path],
        Tuple[Tuple[int, int], ...],
        Path,
        str,
        str,
        Tuple[str, ...],
    ]:
        # A bare hand-built config without the V4 S2 method block is retained
        # for old offline fixtures only.  The shipped defaults always contain
        # ``orca_gfn2_scan`` and therefore cannot enter this branch.
        legacy_fixture_mode = (
            str(self.step2_cfg.get("method", "") or "").strip().lower() == "xtb_path"
            or "xtb_path" in self.step2_cfg
            or (
                "method" not in self.step2_cfg
                and "orca_gfn2_scan" not in self.step2_cfg
                and "xtb_path" not in self.step2_cfg
            )
        )
        if not legacy_fixture_mode:
            return self._run_orca_gfn2_workflow(
                product_xyz, output_dir, forming_bonds, scan_config
            )
        # The legacy xTB coarse/PATH implementation below is retained only as
        # compatibility code for old offline fixtures; it is unreachable from
        # the V4 S2 entrypoint and never launches a calculation.
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        archived_outputs = self._archive_previous_outputs(output_dir)

        bonds = self._validate_forming_bonds(forming_bonds, product_xyz_path=product_xyz)
        product_file = self._resolve_product_file(product_xyz)
        params = self._resolve_scan_params(scan_config)
        self._active_scan_params = dict(params)

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

        degraded_reasons: List[str] = []
        status = "COMPLETE"
        ts_guess_confidence = "high"

        topology_config = self._get_topology_guard_config()
        consecutive_off_path = int(params.get("terminate_after_consecutive_off_path", 2))
        attempts_root = output_dir / "attempts"
        attempts_root.mkdir(parents=True, exist_ok=True)
        attempts: List[ScanAttempt] = []
        composite_candidates: List[ScanAttempt] = []
        coverage_checks: List[Dict[str, Any]] = []

        def execute_attempt(
            kind: str,
            attempt_params: Dict[str, Any],
            *,
            parent_attempt: Optional[ScanAttempt] = None,
            direction: str = "outward",
            seed_target_coordinate_A: Optional[float] = None,
        ) -> Tuple[ScanAttempt, set[int], Optional[int]]:
            attempt_id = f"{len(attempts):02d}_{kind}"
            attempt_dir = attempts_root / attempt_id
            labels = {
                "coarse": "xTB coarse scan",
                "endpoint_extension": "xTB endpoint extension",
                "ts_refinement": "xTB TS local refinement (0.05 Å)",
                "intermediate_refinement": "xTB INT local refinement (0.05 Å)",
            }
            label = labels.get(kind, f"xTB {kind.replace('_', ' ')}")
            batch_id = f"{self.molecule_name or 'product'}:xtb:{attempt_id}"
            configured_total = int(attempt_params["scan_steps"])
            seed_xyz = product_file
            seed_source_attempt: Optional[str] = None
            seed_source_frame_index: Optional[int] = None
            seed_coordinate_A: Optional[float] = None
            if parent_attempt is not None:
                seed_target = (
                    float(seed_target_coordinate_A)
                    if seed_target_coordinate_A is not None
                    else (
                        float(attempt_params["scan_start_distance"])
                        if direction == "inward"
                        else float(attempt_params["scan_end_distance"])
                    )
                )
                (
                    seed_xyz,
                    seed_source_frame_index,
                    seed_coordinate_A,
                ) = self._parent_seed(
                    parent_attempt,
                    seed_target,
                )
                seed_source_attempt = parent_attempt.attempt_id
            attempt_started = time.monotonic()
            self._emit_progress(
                "batch_started",
                batch=batch_id,
                phase="xtb_scan",
                attempt_id=attempt_id,
                attempt_kind=kind,
                label=label,
                current=attempt_id,
                total=configured_total,
                done=0,
                running=1,
                failed=0,
                started_at=time.time(),
                engine="xtb",
                method="GFN2-xTB",
                coordinate_min_A=float(attempt_params["scan_end_distance"]),
                coordinate_max_A=float(attempt_params["scan_start_distance"]),
                seed_xyz=str(seed_xyz),
                seed_source_attempt=seed_source_attempt,
                seed_source_frame_index=seed_source_frame_index,
                seed_coordinate_A=seed_coordinate_A,
                direction=direction,
            )
            try:
                result, attempt_energies = self._execute_scan(
                    start_xyz=seed_xyz,
                    output_dir=attempt_dir,
                    bonds=bonds,
                    params=attempt_params,
                    direction=direction,
                    charge=charge,
                    spin=spin,
                )
            except Exception as exc:
                self._emit_progress(
                    "batch_finished",
                    batch=batch_id,
                    phase="xtb_scan",
                    attempt_id=attempt_id,
                    label=label,
                    current=attempt_id,
                    total=configured_total,
                    done=0,
                    running=0,
                    failed=1,
                    status="failed",
                    error=str(exc),
                    elapsed_seconds=time.monotonic() - attempt_started,
                )
                raise
            attempt, attempt_off_path, attempt_persistent = self._make_scan_attempt(
                attempt_id=attempt_id,
                kind=kind,
                directory=attempt_dir,
                scan_result=result,
                energies=attempt_energies,
                params=attempt_params,
                bonds=bonds,
                topology_config=topology_config,
                consecutive_off_path=consecutive_off_path,
                reference_coords=np.asarray(coords, dtype=float).tolist(),
                reference_symbols=symbols,
                reference_path=product_file,
                parent_attempt_id=parent_attempt.attempt_id if parent_attempt else None,
                seed_xyz=seed_xyz,
                seed_source_attempt=seed_source_attempt,
                seed_source_frame_index=seed_source_frame_index,
                seed_coordinate_A=seed_coordinate_A,
                direction=direction,
            )
            attempts.append(attempt)
            actual_total = len(attempt.frame_paths)
            self._emit_progress(
                "batch_finished",
                batch=batch_id,
                phase="xtb_scan",
                attempt_id=attempt_id,
                attempt_kind=kind,
                label=label,
                current=attempt_id,
                total=actual_total,
                done=actual_total,
                running=0,
                failed=0,
                status="complete",
                elapsed_seconds=time.monotonic() - attempt_started,
                off_path_count=len(attempt.off_path_indices),
                output=str(attempt.directory),
            )
            return attempt, attempt_off_path, attempt_persistent

        active_params = dict(params)
        active_attempt, active_off_path, _ = execute_attempt(
            "coarse", active_params
        )
        composite_candidates.append(active_attempt)
        coarse_scan_summary: Dict[str, Any] = {
            "attempt_id": active_attempt.attempt_id,
            "directory": str(active_attempt.directory),
            "energies_hartree": list(active_attempt.xtb_energies_hartree),
            "absolute_energy_maximum_index": max(
                range(len(active_attempt.xtb_energies_hartree)),
                key=active_attempt.xtb_energies_hartree.__getitem__,
            ),
            "reaction_coordinate_angstrom": list(active_attempt.target_coordinates_A),
            "frame_paths": [str(path) for path in active_attempt.frame_paths],
        }

        extension_cfg = dict(params.get("endpoint_extension", {}) or {})
        if bool(extension_cfg.get("enabled", False)):
            for _ in range(max(0, int(extension_cfg.get("max_extensions", 3)))):
                needs_extension, coverage = self._attempt_needs_endpoint_extension(
                    active_attempt,
                    extension_cfg,
                    dict(params.get("anchor_detection", {}) or {}),
                )
                coverage["attempt_id"] = active_attempt.attempt_id
                coverage_checks.append(coverage)
                if not needs_extension:
                    break
                current_max = float(active_params["scan_start_distance"])
                maximum = float(extension_cfg.get("maximum_coordinate_A", current_max))
                new_max = min(
                    maximum,
                    current_max + float(extension_cfg.get("increment_A", 0.25)),
                )
                if new_max <= current_max + 1.0e-8:
                    coverage["extension_stopped"] = "maximum_coordinate_reached"
                    break
                extension_params = dict(active_params)
                extension_params["scan_start_distance"] = new_max
                max_step = float(extension_params.get("max_step_A") or 0.0)
                if max_step > 0.0:
                    extension_params["scan_steps"] = max(
                        int(extension_params["scan_steps"]),
                        int(
                            math.ceil(
                                (new_max - float(extension_params["scan_end_distance"]))
                                / max_step
                            )
                        )
                        + 1,
                    )
                active_attempt, active_off_path, _ = execute_attempt(
                    "endpoint_extension",
                    extension_params,
                    parent_attempt=active_attempt,
                )
                active_params = extension_params
                composite_candidates.append(active_attempt)

        anchor_cfg = dict(params.get("anchor_detection", {}) or {})
        selection_cfg = dict(params.get("selection", {}) or {})
        coarse_anchors = self._detect_scan_anchors(
            active_attempt.xtb_energies_hartree,
            active_attempt.target_coordinates_A,
            sorted(active_off_path),
            int(anchor_cfg.get("persistent_drift_points", 2)),
            anchor_cfg,
        )
        coarse_ts_selection = self._select_node_by_anchor_index(
            coarse_anchors,
            "plateau_onset_index",
            active_attempt.frame_paths,
            active_attempt.target_coordinates_A,
            active_attempt.xtb_energies_hartree,
            sorted(active_off_path),
        )
        coarse_ts_coordinate = float(
            active_attempt.target_coordinates_A[int(coarse_ts_selection["index"])]
        )
        coarse_ts_selection.update(
            {
                "target_distance_angstrom": coarse_ts_coordinate,
                "source_attempt": active_attempt.attempt_id,
                "status": "coarse_ts_guess_1",
            }
        )

        path_start_index = int(coarse_anchors["last_valid_before_drift_index"])
        full_endpoint_index = int(coarse_anchors["scan_endpoint_index"])
        product_index = int(coarse_anchors["product_index"])
        if path_start_index == product_index:
            raise RuntimeError(
                "S2 xTB PATH cannot start and end on the same frame "
                f"(index={path_start_index})"
            )
        path_end_frame = active_attempt.frame_paths[product_index]

        path_cfg = self._xtb_path_config()
        endpoint_strategy = str(
            path_cfg.get("endpoint_strategy", "full_endpoint_with_valid_corridor")
        ).strip().lower()
        explore_full_endpoint = endpoint_strategy == "full_endpoint_with_valid_corridor"
        retain_valid_corridor = bool(path_cfg.get("retain_valid_corridor_path", True))
        path_attempt: Optional[ScanAttempt] = None
        path_metadata: Dict[str, Any] = {}
        valid_path_attempt: Optional[ScanAttempt] = None
        valid_path_metadata: Dict[str, Any] = {}
        full_path_attempt: Optional[ScanAttempt] = None
        full_path_metadata: Dict[str, Any] = {}
        path_context = {
            "topology_config": topology_config,
            "reference_coords": np.asarray(coords, dtype=float).tolist(),
            "reference_symbols": symbols,
            "reference_path": product_file,
            "consecutive_off_path": int(params["terminate_after_consecutive_off_path"]),
        }

        def execute_path_branch(
            name: str, start_index: int,
        ) -> Tuple[Optional[ScanAttempt], Dict[str, Any]]:
            if start_index == product_index:
                return None, {}
            try:
                attempt, metadata = self._execute_xtb_path(
                    start_xyz=active_attempt.frame_paths[start_index],
                    end_xyz=path_end_frame,
                    output_dir=output_dir,
                    forming_bonds=bonds,
                    path_name=name,
                    **path_context,
                )
                attempts.append(attempt)
                return attempt, metadata
            except RuntimeError as exc:
                self.logger.warning("[S2] xTB PATH branch %s failed: %s", name, exc)
                degraded_reasons.append(f"xtb_path_{name}_failed:{exc}")
                return None, {}

        if explore_full_endpoint and full_endpoint_index != path_start_index:
            full_path_attempt, full_path_metadata = execute_path_branch(
                "full_endpoint", full_endpoint_index
            )
        if retain_valid_corridor or full_path_attempt is None:
            valid_path_attempt, valid_path_metadata = execute_path_branch(
                "valid_corridor", path_start_index
            )

        branch_attempts = {
            name: (attempt, metadata)
            for name, attempt, metadata in (
                ("full_endpoint", full_path_attempt, full_path_metadata),
                ("valid_corridor", valid_path_attempt, valid_path_metadata),
            )
            if attempt is not None
        }
        branch_refinements: Dict[str, Dict[str, Any]] = {}
        branch_b973c_energies: Dict[str, List[Optional[float]]] = {}
        branch_evaluations: Dict[str, Dict[str, Any]] = {}
        preferred_b973c = selection_cfg.get("preferred_energy_source") == "b973c"
        peak_cutoff = float(selection_cfg.get("ts_min_prominence_kcal_mol", 0.40))
        weak_peak_cutoff = float(
            selection_cfg.get(
                "valid_corridor_weak_peak_min_prominence_kcal_mol", 0.10
            )
        )
        weak_peak_barrier = float(
            selection_cfg.get(
                "valid_corridor_weak_peak_min_barrier_kcal_mol", 3.0
            )
        )
        full_boundary_min_frames = max(
            1,
            int(
                selection_cfg.get(
                    "full_endpoint_min_clean_frames_from_boundary", 3
                )
            ),
        )

        # Both PATHs must be evaluated before a seed branch is chosen.  The
        # former implementation selected full_endpoint eagerly, which made the
        # valid-corridor B97-3c peak informational only.
        for branch_name, (attempt, metadata) in branch_attempts.items():
            branch_frames = list(attempt.frame_paths)
            refinement = ScanEnergyRefiner(
                self.config,
                output_dir / "energy_refinement",
                event_callback=self.event_callback,
                variant=self.molecule_name or "product",
            ).refine(
                branch_frames,
                point_ids=[f"{branch_name}_{index:04d}" for index in range(len(branch_frames))],
            )
            b97_energies = [
                None if value is None else float(value)
                for value in list(refinement.get("energies_hartree") or [])
            ]
            if len(b97_energies) != len(branch_frames):
                b97_energies = (b97_energies + [None] * len(branch_frames))[:len(branch_frames)]
            branch_refinements[branch_name] = dict(refinement)
            branch_b973c_energies[branch_name] = b97_energies
            b97_complete = bool(b97_energies) and all(value is not None for value in b97_energies)
            method_energies = (
                b97_energies
                if preferred_b973c and b97_complete
                else [float(value) for value in attempt.xtb_energies_hartree]
            )
            connected = self._product_connected_valid_indices(
                len(branch_frames), attempt.off_path_indices
            )
            excluded = sorted(set(range(len(branch_frames))) - set(connected))
            branch_anchors = self._path_selection_anchors(
                energies=attempt.xtb_energies_hartree,
                reaction_coordinate=attempt.target_coordinates_A,
                valid_indices=connected,
                off_path_indices=excluded,
            )
            candidate, _int_candidate, peak_index, _int_index, _int_status = (
                self._select_refined_path_nodes(
                    anchors=branch_anchors,
                    frame_paths=branch_frames,
                    reaction_coordinate=attempt.target_coordinates_A,
                    method_energies=method_energies,
                    off_path_indices=excluded,
                    path_arclength=np.asarray(metadata.get("path_arclength") or [], dtype=float),
                    selection_config=selection_cfg,
                )
            )
            significant_peak = (
                candidate.get("rule") == "refined_curve_local_maximum"
                and float(candidate.get("prominence_kcal_mol") or 0.0) >= peak_cutoff
            )
            path_arclength = list(metadata.get("path_arclength") or [])
            path_arclength_values = (
                [float(value) for value in path_arclength if value is not None]
                if len(path_arclength) == len(branch_frames)
                and all(value is not None for value in path_arclength)
                else None
            )
            if path_arclength_values is not None:
                arclength_lookup = tuple(path_arclength_values)
                ordered_connected = sorted(
                    connected, key=lambda index: arclength_lookup[index]
                )
            else:
                ordered_connected = sorted(connected)
            peak_position = (
                ordered_connected.index(int(peak_index))
                if int(peak_index) in ordered_connected
                else None
            )
            reactant_index = ordered_connected[0] if ordered_connected else None
            peak_energy = (
                None
                if peak_index is None
                else method_energies[int(peak_index)]
            )
            reactant_energy = (
                None
                if reactant_index is None
                else method_energies[int(reactant_index)]
            )
            barrier_from_reactant = (
                (float(peak_energy) - float(reactant_energy))
                * HARTREE_TO_KCAL
                if reactant_index is not None
                and peak_energy is not None
                and reactant_energy is not None
                else None
            )
            local_peak = candidate.get("rule") in {
                "refined_curve_local_maximum",
                "refined_curve_weak_local_maximum",
            }
            topology_clean = not excluded
            weak_credible = bool(
                topology_clean
                and b97_complete
                and local_peak
                and float(candidate.get("prominence_kcal_mol") or 0.0)
                >= weak_peak_cutoff
                and barrier_from_reactant is not None
                and float(barrier_from_reactant) >= weak_peak_barrier
            )
            full_boundary_safe = bool(
                topology_clean
                or (
                    peak_position is not None
                    and int(peak_position) >= full_boundary_min_frames
                )
            )
            branch_evaluations[branch_name] = {
                "eligible": bool(significant_peak),
                "strong_peak": bool(significant_peak),
                "weak_credible": weak_credible,
                "peak_index": int(peak_index),
                "peak_rule": candidate.get("rule"),
                "peak_prominence_kcal_mol": float(candidate.get("prominence_kcal_mol") or 0.0),
                "barrier_from_reactant_kcal_mol": barrier_from_reactant,
                "energy_source": "B97-3c" if preferred_b973c and b97_complete else "xtb_fallback",
                "b973c_complete": bool(b97_complete),
                "topology_clean": topology_clean,
                "clean_frames_before_peak": peak_position,
                "full_endpoint_boundary_safe": full_boundary_safe,
                "product_connected_valid_indices": connected,
                "selection_excluded_indices": excluded,
            }

        selected_branch: Optional[str] = None
        selection_rule: str
        if branch_evaluations.get("valid_corridor", {}).get("eligible"):
            selected_branch = "valid_corridor"
            selection_rule = "precise_significant_peak_preferred"
        elif branch_evaluations.get("valid_corridor", {}).get("weak_credible"):
            selected_branch = "valid_corridor"
            selection_rule = "valid_corridor_credible_weak_peak_preferred"
        elif (
            branch_evaluations.get("full_endpoint", {}).get("eligible")
            and branch_evaluations.get("full_endpoint", {}).get(
                "full_endpoint_boundary_safe"
            )
        ):
            selected_branch = "full_endpoint"
            selection_rule = "full_endpoint_significant_peak_fallback"
        elif branch_evaluations.get("full_endpoint", {}).get("weak_credible"):
            selected_branch = "full_endpoint"
            selection_rule = "full_endpoint_credible_weak_peak_fallback"
        elif "full_endpoint" in branch_attempts:
            selected_branch = "full_endpoint"
            if branch_evaluations.get("full_endpoint", {}).get("topology_clean"):
                selection_rule = "full_endpoint_low_confidence_fallback"
            else:
                selection_rule = "full_endpoint_boundary_seed_fallback"
        elif "valid_corridor" in branch_attempts:
            selected_branch = "valid_corridor"
            selection_rule = "valid_corridor_fallback_no_significant_peak"
        else:
            selection_rule = "coarse_only_no_path_branch"

        selection_decision: Dict[str, Any] = {
            "selected_branch": selected_branch,
            "rule": selection_rule,
            "significant_peak_threshold_kcal_mol": peak_cutoff,
            "weak_peak_threshold_kcal_mol": weak_peak_cutoff,
            "weak_peak_min_barrier_kcal_mol": weak_peak_barrier,
            "candidates": branch_evaluations,
        }
        if selection_rule in {
            "full_endpoint_boundary_seed_fallback",
            "full_endpoint_low_confidence_fallback",
            "valid_corridor_fallback_no_significant_peak",
        }:
            status = "DEGRADED"
            ts_guess_confidence = "low"
            degraded_reasons.append("no_credible_refined_path_peak")
        if selected_branch is not None:
            path_attempt, path_metadata = branch_attempts[selected_branch]
        else:
            self.logger.warning(
                "[S2] all xTB PATH branches failed; falling back to coarse-only selection"
            )
            status = "DEGRADED"
            ts_guess_confidence = "medium"

        if path_attempt is not None:
            frame_paths = list(path_attempt.frame_paths)
            reaction_coordinate = [
                float(value) for value in path_attempt.target_coordinates_A
            ]
            xtb_energies = [
                float(value) for value in path_attempt.xtb_energies_hartree
            ]
            product_connected_indices = self._product_connected_valid_indices(
                len(frame_paths), path_attempt.off_path_indices
            )
            off_path_indices = set(range(len(frame_paths))) - set(product_connected_indices)
            composite_points = [
                {
                    "point_id": f"path_{index:04d}",
                    "composite_index": index,
                    "source_attempt": path_attempt.attempt_id,
                    "source_frame_index": index,
                    "frame_xyz": str(frame),
                    "xyz_sha256": None,
                    "target_coordinate_A": float(coordinate),
                    "xtb_energy_hartree": float(energy),
                    "topology_valid": index not in off_path_indices,
                }
                for index, (frame, coordinate, energy) in enumerate(
                    zip(frame_paths, reaction_coordinate, xtb_energies)
                )
            ]
            composite_profile = {
                "backbone_attempt_id": path_attempt.attempt_id,
                "accepted_attempt_ids": [path_attempt.attempt_id],
                "points": composite_points,
                "coverage": {
                    "coordinate_min_A": min(reaction_coordinate),
                    "coordinate_max_A": max(reaction_coordinate),
                    "point_count": len(composite_points),
                    "topology_valid_point_count": len(product_connected_indices),
                    "complete_xtb_curve": True,
                },
                "continuity_checks": [],
            }
            trajectory_quality = dict(path_attempt.trajectory_quality)
            trajectory_quality.update(
                {
                    "total_frames": len(composite_points),
                    "off_path_indices": sorted(off_path_indices),
                    "off_path_count": len(off_path_indices),
                    "product_connected_valid_indices": product_connected_indices,
                    "coverage_checks": coverage_checks,
                    "source_attempts": [path_attempt.attempt_id],
                    "xtb_path": True,
                    "selection_branch": path_metadata.get("path_name"),
                }
            )
        else:
            frame_paths = [Path(p) for p in active_attempt.frame_paths]
            reaction_coordinate = [
                float(value) for value in active_attempt.target_coordinates_A
            ]
            xtb_energies = [
                float(value) for value in active_attempt.xtb_energies_hartree
            ]
            off_path_indices = set(int(i) for i in active_off_path)
            composite_points = [
                {
                    "point_id": f"coarse_{index:04d}",
                    "composite_index": index,
                    "source_attempt": active_attempt.attempt_id,
                    "source_frame_index": index,
                    "frame_xyz": str(frame),
                    "xyz_sha256": None,
                    "target_coordinate_A": float(coordinate),
                    "xtb_energy_hartree": float(energy),
                    "topology_valid": index not in off_path_indices,
                }
                for index, (frame, coordinate, energy) in enumerate(
                    zip(frame_paths, reaction_coordinate, xtb_energies)
                )
            ]
            composite_profile = {
                "backbone_attempt_id": active_attempt.attempt_id,
                "accepted_attempt_ids": [active_attempt.attempt_id],
                "points": composite_points,
                "coverage": {
                    "coordinate_min_A": min(reaction_coordinate),
                    "coordinate_max_A": max(reaction_coordinate),
                    "point_count": len(composite_points),
                    "topology_valid_point_count": len(
                        [p for p in composite_points if p["topology_valid"]]
                    ),
                    "complete_xtb_curve": True,
                },
                "continuity_checks": [],
            }
            trajectory_quality = {
                "checked": bool(topology_config.get("enabled", True)),
                "total_frames": len(composite_points),
                "off_path_indices": sorted(off_path_indices),
                "off_path_count": len(off_path_indices),
                "coverage_checks": coverage_checks,
                "source_attempts": [active_attempt.attempt_id],
                "xtb_path": False,
            }

        if path_attempt is not None:
            anchors = self._path_selection_anchors(
                energies=xtb_energies,
                reaction_coordinate=reaction_coordinate,
                valid_indices=trajectory_quality["product_connected_valid_indices"],
                off_path_indices=sorted(off_path_indices),
            )
        else:
            anchors = self._detect_scan_anchors(
                xtb_energies,
                reaction_coordinate,
                sorted(off_path_indices),
                int(anchor_cfg.get("persistent_drift_points", 2)),
                anchor_cfg,
            )

        # PATH 模式: 终点应为最后一帧（产物），而非坐标极值
        # _detect_scan_anchors 的 endpoint_index = argmax(coordinate)
        # 这对粗扫正确（最大距离 = 扫描终点），但对 PATH 语义相反
        # （最大距离 = PATH 起点 frame 0，最小距离 = 产物 frame N-1）。
        all_indices = list(range(len(composite_points)))
        refinement_frames = [frame_paths[index] for index in all_indices]
        refinement_point_ids = [
            str(composite_points[index]["point_id"]) for index in all_indices
        ]
        selected_path_name = str(path_metadata.get("path_name", "valid_corridor"))
        valid_refinement = dict(branch_refinements.get(selected_path_name) or {})
        if not valid_refinement:
            valid_refinement = ScanEnergyRefiner(
                self.config,
                output_dir / "energy_refinement",
                event_callback=self.event_callback,
                variant=self.molecule_name or "product",
            ).refine(refinement_frames, point_ids=refinement_point_ids)
        b973c_energies: List[Optional[float]] = [None] * len(composite_points)
        for local_index, composite_index in enumerate(all_indices):
            values = list(valid_refinement.get("energies_hartree") or [])
            if local_index < len(values) and values[local_index] is not None:
                b973c_energies[composite_index] = float(values[local_index])
        energy_refinement = dict(valid_refinement)
        energy_refinement["energies_hartree"] = b973c_energies
        energy_refinement["scope"] = "path_full_coverage"
        energy_refinement["requested_indices"] = all_indices
        energy_refinement["requested_point_ids"] = refinement_point_ids
        completed_refinement_indices = [
            index for index in all_indices if b973c_energies[index] is not None
        ]
        energy_refinement["completed_indices"] = completed_refinement_indices
        energy_refinement["completed_point_ids"] = [
            str(composite_points[index]["point_id"])
            for index in completed_refinement_indices
        ]
        energy_refinement["total_composite_points"] = len(composite_points)

        b973c_complete = bool(completed_refinement_indices) and all(
            b973c_energies[index] is not None for index in all_indices
        )
        method_energies: Sequence[Optional[float]] = (
            b973c_energies
            if preferred_b973c and b973c_complete
            else [float(value) for value in xtb_energies]
        )
        ts_selection_energy_source = (
            "B97-3c" if preferred_b973c and b973c_complete else "xtb_fallback"
        )
        intermediate_energy_source = ts_selection_energy_source
        if preferred_b973c and not b973c_complete:
            status = "DEGRADED"
            ts_guess_confidence = "medium"
            degraded_reasons.append("b973c_path_incomplete")
        b973c_coverage = {
            "selection_scope": "path_full_coverage",
            "path_complete": b973c_complete,
            "available_indices": [
                index for index, value in enumerate(b973c_energies) if value is not None
            ],
            "fallback_reason": None if b973c_complete else "b973c_path_incomplete",
        }
        energy_refinement["coverage"] = b973c_coverage

        path_arclength = (
            np.asarray(path_metadata.get("path_arclength", []), dtype=float)
            if path_metadata and path_metadata.get("path_arclength")
            else None
        )
        xtb_path_profile = build_xtb_path_profile(
            frame_paths=frame_paths,
            energies_hartree=method_energies,
            forming_bonds=bonds,
            product_xyz=product_file,
            off_path_indices=sorted(off_path_indices),
            source_provenance={
                "path_name": selected_path_name,
                "attempt_id": composite_profile.get("backbone_attempt_id"),
                "energy_source": ts_selection_energy_source,
                "point_ids": refinement_point_ids,
                "path_arclength": list(path_arclength) if path_arclength is not None else [],
            },
        )
        xtb_selector_selection = select_path_seeds(
            xtb_path_profile,
            policy_from_config(selection_cfg),
        )
        xtb_selector_payload = xtb_selector_selection.to_dict()
        profile_selection_source = None
        profile_s2_state = str(xtb_selector_payload.get("s2_state") or "unresolved")
        profile_seed_evidence = str(
            xtb_selector_payload.get("seed_evidence") or "none"
        )
        profile_ts_search_seed = (
            None
            if xtb_selector_payload.get("ts_search_seed") is None
            else dict(xtb_selector_payload["ts_search_seed"])
        )
        profile_int_search_seed = (
            None
            if xtb_selector_payload.get("int_search_seed") is None
            else dict(xtb_selector_payload["int_search_seed"])
        )
        profile_has_independent_int = bool(
            xtb_selector_payload.get("has_independent_int", False)
        )
        profile_rejection_reason = xtb_selector_payload.get("rejection_reason")
        profile_selection_diagnostics = dict(
            xtb_selector_payload.get("diagnostics") or {}
        )

        selection_started = time.monotonic()
        self._emit_progress(
            "step_started",
            step="node_selection_render",
            label="TS/INT selection and profile rendering",
            index=3,
            total_steps=3,
            point_count=len(composite_points),
        )

        (
            ts_selection,
            intermediate_selection,
            max_idx,
            intermediate_idx,
            intermediate_selection_status,
        ) = self._unified_selection_records(
            profile=xtb_path_profile,
            selection=xtb_selector_selection,
            method_energies=method_energies,
            energy_source=ts_selection_energy_source,
        )
        ts_selection.update(
            {
                "xtb_path_estimated_ts_point": path_metadata.get("estimated_ts_point"),
                "source_branch": selected_path_name,
                "target_coordinate_A": (
                    None
                    if ts_selection.get("index") is None
                    else float(reaction_coordinate[int(ts_selection["index"])])
                ),
                "point_id": (
                    None
                    if ts_selection.get("index") is None
                    else str(composite_points[int(ts_selection["index"])] ["point_id"])
                ),
            }
        )
        if intermediate_selection.get("index") is not None:
            intermediate_selection.update(
                {
                    "target_coordinate_A": float(
                        reaction_coordinate[int(intermediate_selection["index"])]
                    ),
                    "point_id": str(
                        composite_points[int(intermediate_selection["index"])] ["point_id"]
                    ),
                    "source_branch": selected_path_name,
                }
            )
        ts_guess_confidence = str(ts_selection.get("confidence", ts_guess_confidence))
        selection_decision = {
            "selected_branch": selected_path_name,
            "rule": ts_selection.get("rule"),
            "configured_method": "unified_selector",
            "actual_method": ts_selection.get("actual_method"),
        }

        forming_bond_distances = self._forming_bond_distances_by_frame(frame_paths, bonds)
        ts_peak_idx = None if max_idx is None else int(max_idx)
        ts_seed_idx = ts_selection.get("index")

        path_requires_rescue = self._path_requires_relaxed_scan(
            trajectory_quality
        )
        selection_source = (
            "xtb_peb"
            if xtb_selector_selection.s2_state == "path_seeded"
            else None
        )
        path_seed_available = (
            selection_source == "xtb_peb"
            and not path_requires_rescue
            and ts_seed_idx is not None
        )
        if not path_requires_rescue and not path_seed_available:
            status = "DEGRADED"
            degraded_reasons.append("s2_selector_unresolved")
        s3_dispatch: Dict[str, Any] = {
            "resolution": "path_seeded" if path_seed_available else "unresolved",
            "submit_ts": bool(path_seed_available),
            "submit_intermediate": bool(path_seed_available and intermediate_idx is not None),
            # S2 only provides seeds.  Whether an optimized intermediate is
            # suitable for an S3 endpoint is deliberately outside this gate.
            "neb_eligible": False,
            "source": "peb_path",
            "path_fully_distorted": path_requires_rescue,
        }
        profile_selection_source = selection_source
        rescue_payload: Dict[str, Any] = {"status": "not_required"}
        ts_seed_source = (
            Path(frame_paths[int(ts_seed_idx)])
            if path_seed_available
            else None
        )
        intermediate_seed_idx_raw = intermediate_selection.get("index")
        intermediate_seed_idx = (
            None
            if intermediate_seed_idx_raw is None
            else int(intermediate_seed_idx_raw)
        )
        intermediate_source = (
            Path(frame_paths[intermediate_seed_idx])
            if intermediate_seed_idx is not None and path_seed_available
            else None
        )
        if intermediate_source is not None:
            # A shared TS/INT frame is still a legitimate pair of independent
            # S3 search seeds on a usable path.  Keep its index so the
            # published intermediate.xyz reaches the orchestrator.
            intermediate_idx = intermediate_seed_idx
        if path_requires_rescue:
            rescue_payload = B97CRelaxedScanRescuer(
                self.config,
                event_callback=self.event_callback,
                variant=self.molecule_name,
            ).run(
                product_file,
                output_dir / "rescue",
                bonds,
                trigger_reasons=["path_fully_topology_distorted"],
            )
            profile_selection_source = "orca_relaxed_scan"
            profile_s2_state = str(
                rescue_payload.get("s2_state")
                or rescue_payload.get("resolution")
                or "unresolved"
            )
            profile_seed_evidence = str(
                rescue_payload.get("seed_evidence") or "none"
            )
            profile_ts_search_seed = (
                None
                if rescue_payload.get("ts_search_seed") is None
                else dict(rescue_payload["ts_search_seed"])
            )
            profile_int_search_seed = (
                None
                if rescue_payload.get("int_search_seed") is None
                else dict(rescue_payload["int_search_seed"])
            )
            profile_has_independent_int = bool(
                rescue_payload.get("has_independent_int", False)
            )
            profile_rejection_reason = rescue_payload.get("rejection_reason")
            profile_selection_diagnostics = dict(
                rescue_payload.get("selection_diagnostics") or {}
            )
            rescue_state = profile_s2_state
            if rescue_state == "rescue_seeded":
                selection_source = "orca_relaxed_scan"
                for key in (
                    "energy_peak_rule",
                    "energy_peak_candidate_indices",
                    "energy_peak_coordinate_A",
                    "seed_index",
                    "seed_rule",
                    "seed_backoff_requested_A",
                    "seed_backoff_applied_A",
                    "seed_target_mean_forming_bond_distance_A",
                    "seed_mean_forming_bond_distance_A",
                    "seed_backoff_status",
                    "clean_path_boundary_index",
                    "xtb_path_estimated_ts_point",
                ):
                    ts_selection.pop(key, None)
                for key in (
                    "pre_ts_floor_kcal_mol",
                    "barrier_from_pre_ts_floor_kcal_mol",
                    "energy_window_kcal_mol",
                    "energy_limit_kcal_mol",
                    "platform_max_slope_kcal_mol_A",
                    "clean_pre_ts_frame_count",
                    "shared_ts_index",
                ):
                    intermediate_selection.pop(key, None)
                s3_dispatch = {
                    "resolution": "rescue_seeded",
                    "submit_ts": bool(rescue_payload.get("ts_xyz")),
                    "submit_intermediate": bool(rescue_payload.get("intermediate_xyz")),
                    "neb_eligible": False,
                    "source": "b973c_relaxed_scan",
                    "path_fully_distorted": True,
                }
                ts_selection_energy_source = "B97-3c"
                intermediate_energy_source = "B97-3c"
                ts_seed_source = Path(str(rescue_payload["ts_xyz"]))
                intermediate_raw = rescue_payload.get("intermediate_xyz")
                intermediate_source = Path(str(intermediate_raw)) if intermediate_raw else None
                rescue_ts_seed = profile_ts_search_seed or {}
                rescue_int_seed = profile_int_search_seed or {}
                rescue_ts_index = self._seed_frame_index(rescue_ts_seed)
                rescue_int_index = self._seed_frame_index(rescue_int_seed)
                rescue_int_mode = self._selector_int_selection_mode(
                    seed_evidence=profile_seed_evidence,
                    int_seed=rescue_int_seed,
                    has_independent_int=profile_has_independent_int,
                )
                ts_selection.update(
                    {
                        "index": rescue_ts_index,
                        "frame_xyz": str(ts_seed_source),
                        "rule": self._selector_ts_rule(
                            profile_seed_evidence,
                            fallback="b973c_relaxed_scan_seed",
                        ),
                        "actual_method": (
                            "b973c_scants"
                            if rescue_payload.get("scan_ts")
                            else "b973c_simultaneous_relaxed_scan"
                        ),
                        "confidence": rescue_ts_seed.get("confidence", "medium"),
                        "source_branch": "s2_rescue",
                        "source": selection_source,
                        "candidate_indices": (
                            [] if rescue_ts_index is None else [int(rescue_ts_index)]
                        ),
                    }
                )
                ts_guess_confidence = str(
                    rescue_ts_seed.get("confidence", ts_guess_confidence)
                )
                intermediate_selection_status = (
                    "shared_with_ts"
                    if rescue_int_mode == "shared_ts_fallback"
                    else "selected"
                    if rescue_int_seed
                    else "unavailable"
                )
                intermediate_selection.update(
                    {
                        "index": rescue_int_index,
                        "frame_xyz": (
                            str(intermediate_source)
                            if intermediate_source is not None
                            else rescue_int_seed.get("xyz")
                        ),
                        "rule": self._selector_int_rule(
                            selection_source=selection_source,
                            seed_evidence=profile_seed_evidence,
                            int_seed=rescue_int_seed,
                            has_independent_int=profile_has_independent_int,
                        ),
                        "selection_mode": rescue_int_mode,
                        "selection_status": intermediate_selection_status,
                        "stationary_point_claimed": False,
                        "reason": (
                            None
                            if rescue_int_seed and rescue_int_mode == "stable_basin_candidate"
                            else "no_resolved_pre_ts_basin"
                            if rescue_int_mode == "shared_ts_fallback"
                            else "no_int_search_seed"
                            if not rescue_int_seed
                            else "search_seed"
                        ),
                        "source_branch": "s2_rescue",
                        "source": selection_source,
                        "candidate_indices": (
                            []
                            if rescue_int_index is None
                            else [int(rescue_int_index)]
                        ),
                        "shared_ts_index": (
                            int(rescue_ts_index)
                            if rescue_int_mode == "shared_ts_fallback"
                            and rescue_ts_index is not None
                            else None
                        ),
                    }
                )
                intermediate_idx = rescue_int_index
            else:
                s3_dispatch.update(
                    {
                        "resolution": "unresolved",
                        "submit_ts": False,
                        "submit_intermediate": False,
                        "source": "b973c_relaxed_scan",
                    }
                )
                selection_source = None
                profile_selection_source = None
                ts_seed_source = None
                intermediate_source = None
                intermediate_idx = None
                ts_selection = {
                    "index": None,
                    "frame_xyz": None,
                    "rule": "unresolved_rescue",
                    "selection_status": "unavailable",
                    "reason": profile_rejection_reason or "s2_rescue_unresolved",
                    "source": None,
                    "candidate_indices": [],
                }
                intermediate_selection = {
                    "index": None,
                    "frame_xyz": None,
                    "rule": "unresolved_rescue",
                    "selection_mode": "unavailable",
                    "selection_status": "unavailable",
                    "reason": "s2_rescue_unresolved",
                    "source": None,
                    "candidate_indices": [],
                }
                degraded_reasons.append("s2_rescue_unresolved")
                status = "DEGRADED"

        ts_guess_xyz_final = output_dir / "ts_guess.xyz"
        if ts_seed_source is not None and s3_dispatch["submit_ts"]:
            self._atomic_copy(ts_seed_source, ts_guess_xyz_final)

        dipolar_xyz = output_dir / INTERMEDIATE_XYZ
        if intermediate_source is not None and s3_dispatch["submit_intermediate"]:
            self._atomic_copy(intermediate_source, dipolar_xyz)
        else:
            intermediate_idx = None

        xtb_relative_energies = [
            (float(energy) - float(xtb_energies[0])) * HARTREE_TO_KCAL
            for energy in xtb_energies
        ]
        b973c_reference_index = next(
            (index for index, value in enumerate(b973c_energies) if value is not None),
            None,
        )
        b973c_relative_energies: List[Optional[float]] = [None] * len(b973c_energies)
        if b973c_reference_index is not None:
            reference_energy = self._require_energy(
                b973c_energies[b973c_reference_index],
                "S2 B97-3c reference energy is missing",
            )
            b973c_relative_energies = [
                (float(value) - reference_energy) * HARTREE_TO_KCAL
                if value is not None
                else None
                for value in b973c_energies
            ]

        valid_branch_refinement: Optional[Dict[str, Any]] = None
        if valid_path_attempt is not None and valid_path_attempt is not path_attempt:
            valid_frames = list(valid_path_attempt.frame_paths)
            valid_branch_refinement = dict(
                branch_refinements.get("valid_corridor") or {}
            )
            if not valid_branch_refinement:
                valid_branch_refinement = ScanEnergyRefiner(
                    self.config,
                    output_dir / "energy_refinement",
                    event_callback=self.event_callback,
                    variant=self.molecule_name or "product",
                ).refine(
                    valid_frames,
                    point_ids=[f"valid_corridor_{index:04d}" for index in range(len(valid_frames))],
                )
        def path_branch_payload(
            name: str,
            attempt: ScanAttempt,
            metadata: Mapping[str, Any],
            b97_energies: Sequence[Optional[float]],
            refinement: Optional[Mapping[str, Any]],
            *,
            selected_for_s3: bool,
        ) -> Dict[str, Any]:
            branch_xtb = [float(value) for value in attempt.xtb_energies_hartree]
            branch_b97 = [
                None if value is None else float(value) for value in b97_energies
            ]
            product_connected = self._product_connected_valid_indices(
                len(attempt.frame_paths), attempt.off_path_indices
            )
            excluded = sorted(
                set(range(len(attempt.frame_paths))) - set(product_connected)
            )
            return {
                "name": name,
                "role": (
                    "selection_path"
                    if selected_for_s3
                    else "exploratory_diagnostic"
                    if name == "full_endpoint"
                    else "valid_corridor_control"
                ),
                "selected_for_s3": selected_for_s3,
                "metadata": dict(metadata),
                "point_ids": [f"{name}_{index:04d}" for index in range(len(attempt.frame_paths))],
                "frame_paths": [str(frame) for frame in attempt.frame_paths],
                "reaction_coordinate_angstrom": [
                    float(value) for value in attempt.target_coordinates_A
                ],
                "path_arclength": list(metadata.get("path_arclength") or []),
                "forming_bond_distances_angstrom": self._forming_bond_distances_by_frame(
                    attempt.frame_paths, bonds
                ),
                "trajectory_quality": {
                    **dict(attempt.trajectory_quality),
                    "product_connected_valid_indices": product_connected,
                    "selection_excluded_indices": excluded,
                },
                "energy_refinement": dict(refinement or {}),
                "energy_curves": {
                    "xtb": {
                        "energies_hartree": branch_xtb,
                        "relative_energies_kcal_mol": self._branch_relative_energies(branch_xtb),
                    },
                    "b973c": {
                        "energies_hartree": branch_b97,
                        "relative_energies_kcal_mol": self._branch_relative_energies(branch_b97),
                    },
                },
            }

        path_branches: Dict[str, Any] = {}
        for branch_name, (branch_attempt, branch_metadata) in branch_attempts.items():
            if branch_name == selected_path_name:
                branch_b97 = b973c_energies
                branch_refinement = energy_refinement
            else:
                branch_b97 = branch_b973c_energies.get(branch_name, [])
                branch_refinement = branch_refinements.get(branch_name, {})
            path_branches[branch_name] = path_branch_payload(
                branch_name,
                branch_attempt,
                branch_metadata,
                branch_b97,
                branch_refinement,
                selected_for_s3=branch_name == selected_path_name,
            )
            path_branches[branch_name]["selection_candidate"] = dict(
                branch_evaluations.get(branch_name) or {}
            )

        ts_selection["source"] = selection_source
        intermediate_selection["source"] = selection_source
        published_ts_seed_idx_raw = ts_selection.get("index")
        published_ts_seed_idx = (
            None
            if published_ts_seed_idx_raw is None
            else int(published_ts_seed_idx_raw)
        )
        published_ts_peak_idx = (
            published_ts_seed_idx
            if selection_source == "orca_relaxed_scan"
            else None
            if ts_peak_idx is None
            else int(ts_peak_idx)
        )
        published_intermediate_idx = None if intermediate_idx is None else int(intermediate_idx)

        selection_coordinates = [float(value) for value in reaction_coordinate]
        selection_point_ids = [str(point["point_id"]) for point in composite_points]
        selection_source_attempts = [
            str(point["source_attempt"]) for point in composite_points
        ]
        selection_method_energies = [
            None if value is None else float(value) for value in method_energies
        ]
        selection_relative_energies = (
            b973c_relative_energies
            if ts_selection_energy_source == "B97-3c"
            else xtb_relative_energies
        )
        selection_xtb_energies: List[Optional[float]] = [
            float(value) for value in xtb_energies
        ]
        selection_source_branch = selected_path_name
        if selection_source == "orca_relaxed_scan":
            rescue_frames = [
                Path(value) for value in list(rescue_payload.get("frames") or [])
            ]
            rescue_energies = [
                None if value is None else float(value)
                for value in list(rescue_payload.get("energies_hartree") or [])
            ]
            rescue_forming_bond_distances = self._forming_bond_distances_by_frame(
                rescue_frames,
                bonds,
            )
            selection_coordinates = self._mean_frame_coordinates(
                rescue_forming_bond_distances
            )
            selection_point_ids = [
                f"rescue_{index:04d}" for index in range(len(rescue_frames))
            ]
            selection_source_attempts = [
                "b973c_relaxed_scan" for _ in range(len(rescue_frames))
            ]
            selection_method_energies = rescue_energies
            selection_relative_energies = self._branch_relative_energies(
                rescue_energies,
                reference_side="last",
            )
            selection_xtb_energies = [None] * len(rescue_frames)
            selection_source_branch = "s2_rescue"

        if published_intermediate_idx is not None:
            intermediate_energy = selection_method_energies[published_intermediate_idx]
            if intermediate_energy is None:
                raise RuntimeError(
                    "S2 selected an intermediate without a method-consistent energy"
                )
            intermediate_relative_energy = selection_relative_energies[
                published_intermediate_idx
            ]
            intermediate_xtb_energy = selection_xtb_energies[published_intermediate_idx]
            int_selection_mode = str(
                intermediate_selection.get("selection_mode", "midpoint_fallback")
            )
            int_method_label = (
                "PRE_TS_BASIN"
                if int_selection_mode == "stable_basin_candidate"
                else "PRE_TS_PLATFORM"
                if int_selection_mode == "late_pre_ts_platform_fallback"
                else "TS_ENDPOINT_MIDPOINT"
                if int_selection_mode == "ts_to_effective_endpoint_midpoint"
                else "PRE_TS_SEARCH"
            )
            intermediate_selection.update(
                {
                    "status": f"B973C_{int_method_label}_INT_SEARCH_SEED"
                    if intermediate_energy_source == "B97-3c"
                    else f"XTB_{int_method_label}_INT_SEARCH_SEED",
                    "role": "int_search_seed",
                    "stationary_point_claimed": False,
                    "confirmation_level": "guess_only",
                    "validation_required": [
                        "s3_unconstrained_optimization",
                        "s3_frequency",
                    ],
                    "target_distance_angstrom": float(
                        selection_coordinates[published_intermediate_idx]
                    ),
                    "point_id": selection_point_ids[published_intermediate_idx],
                    "source_attempt": selection_source_attempts[
                        published_intermediate_idx
                    ],
                    "energy_hartree": float(intermediate_energy),
                    "energy_source": intermediate_energy_source,
                    "target_coordinate_A": float(
                        selection_coordinates[published_intermediate_idx]
                    ),
                    "xtb_energy_hartree": (
                        None
                        if intermediate_xtb_energy is None
                        else float(intermediate_xtb_energy)
                    ),
                    "relative_energy_kcal_mol": (
                        float(intermediate_relative_energy)
                        if intermediate_relative_energy is not None
                        else None
                    ),
                    "selection_status": intermediate_selection_status,
                    "source_branch": selection_source_branch,
                    "source": selection_source,
                    "output_xyz": str(dipolar_xyz),
                }
            )
        if published_ts_seed_idx is not None and published_ts_peak_idx is not None:
            ts_method_label = "SCAN_TS" if selection_source == "orca_relaxed_scan" else "PATH_TS"
            ts_peak_xtb_energy = selection_xtb_energies[published_ts_peak_idx]
            ts_seed_xtb_energy = selection_xtb_energies[published_ts_seed_idx]
            ts_relative_energy = selection_relative_energies[published_ts_peak_idx]
            ts_selection.update(
                {
                    "status": f"B973C_LOCAL_{ts_method_label}_TS_GUESS"
                    if ts_selection_energy_source == "B97-3c"
                    else f"XTB_LOCAL_{ts_method_label}_TS_GUESS",
                    "confirmation_level": "guess_only",
                    "validation_required": ["s3_opt_ts", "s3_frequency"],
                    "output_xyz": str(ts_guess_xyz_final) if ts_guess_xyz_final.exists() else None,
                    "target_distance_angstrom": float(
                        selection_coordinates[published_ts_seed_idx]
                    ),
                    "target_coordinate_A": float(
                        selection_coordinates[published_ts_seed_idx]
                    ),
                    "point_id": selection_point_ids[published_ts_seed_idx],
                    "source_attempt": selection_source_attempts[published_ts_seed_idx],
                    "energy_hartree": self._require_energy(
                        selection_method_energies[published_ts_peak_idx],
                        "S2 selected a TS without a method-consistent energy",
                    ),
                    "seed_energy_hartree": self._require_energy(
                        selection_method_energies[published_ts_seed_idx],
                        "S2 selected a TS seed without a method-consistent energy",
                    ),
                    "energy_source": ts_selection_energy_source,
                    "xtb_energy_hartree": (
                        None
                        if ts_peak_xtb_energy is None
                        else float(ts_peak_xtb_energy)
                    ),
                    "seed_xtb_energy_hartree": (
                        None
                        if ts_seed_xtb_energy is None
                        else float(ts_seed_xtb_energy)
                    ),
                    "relative_energy_kcal_mol": (
                        None
                        if ts_relative_energy is None
                        else float(ts_relative_energy)
                    ),
                    "source_branch": selection_source_branch,
                    "source": selection_source,
                }
            )
        self._validate_unified_selection_contract(
            s2_state=profile_s2_state,
            selection_source=selection_source,
            ts_selection=ts_selection,
            intermediate_selection=intermediate_selection,
            s3_dispatch=s3_dispatch,
            ts_guess_xyz=ts_guess_xyz_final,
            intermediate_xyz=dipolar_xyz,
        )
        trajectory_quality.update(
            {
                "endpoint_index": int(anchors["scan_endpoint_index"]),
                "endpoint_retained": True,
                "endpoint_excluded": int(anchors["scan_endpoint_index"])
                in off_path_indices,
                "topology_drift_index": anchors["topology_drift_index"],
                "last_valid_before_drift_index": anchors[
                    "last_valid_before_drift_index"
                ],
            }
        )
        for point, distances in zip(composite_points, forming_bond_distances):
            point["actual_forming_bond_distances_A"] = distances
        accepted_attempt_ids = set(composite_profile.get("accepted_attempt_ids", []))
        attempt_records = attempt_manifest(attempts)
        for attempt_record in attempt_records:
            attempt_record["selected_for_composite"] = (
                attempt_record.get("attempt_id") in accepted_attempt_ids
            )
        anchor_keys = {
            "product": "product_index",
            "coarse_ts": "coarse_ts_index",
            "plateau_onset": "plateau_onset_index",
            "absolute_energy_maximum": "absolute_energy_maximum_index",
            "topology_drift": "topology_drift_index",
            "last_valid_before_drift": "last_valid_before_drift_index",
            "scan_endpoint": "scan_endpoint_index",
        }
        anchor_payload: Dict[str, Any] = {}
        for name, key in anchor_keys.items():
            anchor_index = anchors.get(key)
            anchor_payload[name] = None if anchor_index is None else {
                "index": int(anchor_index),
                "point_id": composite_points[int(anchor_index)]["point_id"],
                "frame_xyz": str(frame_paths[int(anchor_index)]),
                "target_distance_angstrom": float(reaction_coordinate[int(anchor_index)]),
                "energy_hartree": (
                    float(xtb_energies[int(anchor_index)])
                ),
                "energy_source": "xTB",
                "topology_valid": int(anchor_index) not in off_path_indices,
            }

        scan_profile_json = output_dir / "scan_profile.json"
        profile_payload: Dict[str, Any] = {
            "profile_schema_version": "s2_scan_profile_v10",
            "generation_method": "xtb_path_full_coverage",
            "product_xyz": str(product_file),
            "intermediate_xyz": str(dipolar_xyz) if intermediate_idx is not None else None,
            "selection_source": profile_selection_source,
            "s2_state": profile_s2_state,
            "seed_evidence": profile_seed_evidence,
            "endpoint_evidence": dict(
                profile_selection_diagnostics.get("endpoints") or {}
            ),
            "knee_evidence": {
                "frame_index": profile_selection_diagnostics.get("knee_frame_index"),
                "coordinate_A": profile_selection_diagnostics.get("knee_coordinate_A"),
                "anchor_type": profile_selection_diagnostics.get("knee_anchor_type"),
                "right_shift_A": profile_selection_diagnostics.get(
                    "ts_right_shift_applied_A"
                ),
            },
            "ts_search_seed": profile_ts_search_seed,
            "int_search_seed": profile_int_search_seed,
            "has_independent_int": profile_has_independent_int,
            "rejection_reason": profile_rejection_reason,
            "forming_bonds": [list(pair) for pair in bonds],
            "scan_parameters": params,
            "coarse_scan": coarse_scan_summary,
            "attempts": attempt_records,
            "composite_profile": composite_profile,
            "anchors": anchor_payload,
            "xtb_path": path_metadata if path_attempt is not None else None,
            "xtb_paths": {
                "full_endpoint": full_path_metadata or None,
                "valid_corridor": valid_path_metadata or None,
            },
            "path_branches": path_branches,
            "selection_decision": {
                **selection_decision,
                "ts_seed_index": published_ts_seed_idx,
                "energy_peak_index": published_ts_peak_idx,
                "ts_confidence": ts_selection.get("confidence"),
                "int_policy": intermediate_selection.get("selection_mode"),
            },
            "seed_selection": {
                "authority": "unified_selector",
                "source": selection_source,
                "s2_state": profile_s2_state,
                "seed_evidence": profile_seed_evidence,
                "rejection_reason": profile_rejection_reason,
                "ts_search_seed": profile_ts_search_seed,
                "int_search_seed": profile_int_search_seed,
                "has_independent_int": profile_has_independent_int,
                "diagnostics": profile_selection_diagnostics,
            },
            "s3_dispatch": s3_dispatch,
            "rescue": rescue_payload,
            "coarse_candidates": {
                "ts_guess_1": coarse_ts_selection,
            },
            "derivative_analysis": {
                "method": ts_selection_energy_source,
                "rule": selection_decision["rule"],
                "actual_method": ts_selection.get("actual_method"),
                "selected_index": published_ts_seed_idx,
                "energy_peak_index": published_ts_peak_idx,
            },
            "scan_quality": {
                "ts_index": published_ts_seed_idx,
                "energy_peak_index": published_ts_peak_idx,
                "ts_seed_index": published_ts_seed_idx,
                "plateau_onset_index": int(anchors["plateau_onset_index"]),
                "absolute_energy_maximum_index": int(
                    anchors["absolute_energy_maximum_index"]
                ),
                "selection_energy_source": ts_selection_energy_source,
                "intermediate_energy_source": intermediate_energy_source,
                "intermediate_index": published_intermediate_idx,
                "intermediate_selection_status": intermediate_selection_status,
                "status": status,
                "ts_guess_confidence": ts_guess_confidence,
                "selected_path_branch": selection_source_branch,
                "intermediate_confidence": (
                    "high"
                    if intermediate_selection.get("selection_mode") == "stable_basin_candidate"
                    else "medium"
                    if intermediate_selection.get("selection_mode") == "late_pre_ts_platform_fallback"
                    else "medium"
                    if intermediate_selection.get("selection_mode")
                    == "ts_to_effective_endpoint_midpoint"
                    else "unresolved_shared_ts"
                    if intermediate_selection.get("selection_mode") == "shared_ts_fallback"
                    else "unavailable"
                ),
                "degraded_reasons": degraded_reasons,
            },
            "reaction_coordinate_angstrom": reaction_coordinate,
            "forming_bond_distances_angstrom": forming_bond_distances,
            "energies_hartree": xtb_energies,
            "relative_energies_kcal_mol": xtb_relative_energies,
            "xtb_energies_hartree": xtb_energies,
            "xtb_relative_energies_kcal_mol": xtb_relative_energies,
            "energy_refinement": energy_refinement,
            "energy_curves": {
                "xtb": {
                    "status": "complete",
                    "point_ids": [point["point_id"] for point in composite_points],
                    "energies_hartree": xtb_energies,
                    "relative_energies_kcal_mol": xtb_relative_energies,
                },
                "b973c": {
                    "status": energy_refinement.get("status", "not_requested"),
                    "point_ids": [point["point_id"] for point in composite_points],
                    "energies_hartree": b973c_energies,
                    "relative_energies_kcal_mol": b973c_relative_energies,
                    "reference_point_index": b973c_reference_index,
                },
            },
            "selection_policy": {
                "preferred_source": selection_cfg.get("preferred_energy_source", "b973c"),
                "actual_source": ts_selection_energy_source,
                "ts_source": ts_selection_energy_source,
                "intermediate_source": intermediate_energy_source,
                "xtb_fallback_used": (
                    ts_selection_energy_source.startswith("xtb")
                    or intermediate_energy_source.startswith("xtb")
                ),
                "fallback_reason": b973c_coverage[
                    "fallback_reason"
                ],
                "b973c_coverage": b973c_coverage,
            },
            "selections": {
                "ts_guess": ts_selection,
                "intermediate": intermediate_selection,
            },
            "trajectory_quality": trajectory_quality,
            "artifact_lifecycle": {
                "stale_outputs_archived_to": str(archived_outputs) if archived_outputs else None,
                "legacy_reactant_complex_written": False,
            },
            "scan_plot": None,
            "plot_provenance": {
                "x_coordinate_type": "mean_forming_bond_distance",
                "x_min_angstrom": min(reaction_coordinate) if reaction_coordinate else None,
                "x_max_angstrom": max(reaction_coordinate) if reaction_coordinate else None,
                "n_points": len(reaction_coordinate),
                "generation_method": "xtb_path_full_coverage",
                "path_arclength_available": bool(path_metadata and path_metadata.get("path_arclength")),
            },
        }
        write_text_atomic(
            scan_profile_json,
            json.dumps(profile_payload, indent=2),
            encoding="utf-8",
        )

        try:
            diagnostic_png = scan_profile_json.parent / "scan_profile.png"

            plot_path = plot_scan_profile(
                scan_profile_json,
                output_path=diagnostic_png,
            )
            if plot_path is not None:
                profile_payload["scan_plot"] = str(plot_path)
                profile_payload["scan_figures"] = {
                    "profile": {"png": str(diagnostic_png)}
                }

            if plot_path is not None or profile_payload.get("scan_figures"):
                write_text_atomic(
                    scan_profile_json,
                    json.dumps(profile_payload, indent=2),
                    encoding="utf-8",
                )
        except Exception as exc:
            self.logger.warning("[S2] Failed to render scan profile: %s", exc, exc_info=True)

        self._emit_progress(
            "step_finished",
            step="node_selection_render",
            label="TS/INT selection and profile rendering",
            index=3,
            total_steps=3,
            status="complete" if status == "COMPLETE" else "degraded",
            point_count=len(composite_points),
            ts_index=(None if max_idx is None else int(max_idx)),
            intermediate_index=intermediate_idx,
            elapsed_seconds=time.monotonic() - selection_started,
        )

        return (
            ts_guess_xyz_final if ts_guess_xyz_final.exists() else None,
            dipolar_xyz if intermediate_idx is not None else None,
            dipolar_xyz if intermediate_idx is not None else None,
            bonds,
            scan_profile_json,
            status,
            ts_guess_confidence,
            tuple(degraded_reasons),
        )
