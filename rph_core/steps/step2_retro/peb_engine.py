import json
import logging
import math
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from rph_core.utils.file_io import read_xyz
from rph_core.utils.bond_pairs import canonicalize_bond_pairs
from rph_core.utils.geometry_tools import GeometryUtils, LogParser, kabsch_rmsd
from rph_core.utils.json_io import write_text_atomic
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

from .energy_refinement import ScanEnergyRefiner
from .geometry_guard import (
    check_scan_trajectory,
)
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
                    float(selection_cfg.get("ts_min_prominence_kcal_mol", 0.15)),
                ),
                "int_min_basin_prominence_kcal_mol": max(
                    0.0,
                    float(
                        selection_cfg.get(
                            "int_min_basin_prominence_kcal_mol", 0.50
                        )
                    ),
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
    def _select_int_stable_point(
        anchors: Dict[str, Any],
        frame_paths: Sequence[Path],
        reaction_coordinate: Sequence[float],
        energies: Sequence[Optional[float]],
        off_path_indices: Sequence[int],
    ) -> Dict[str, Any]:
        """Pick the topology-valid frame with the smallest local |gradient|
        on the M -> last_valid plateau; fall back to the plateau midpoint
        when there is insufficient derivative support.
        """
        m_index = int(anchors["plateau_onset_index"])
        last_valid = int(anchors["last_valid_before_drift_index"])
        off_path = {int(index) for index in off_path_indices}
        candidates = [
            index
            for index in range(m_index + 1, last_valid + 1)
            if index not in off_path and energies[index] is not None
        ]
        if not candidates:
            raise RuntimeError("S2 INT selection: empty M->last_valid plateau")
        if len(candidates) == 1:
            selected = candidates[0]
            rule = "single_plateau_frame"
        else:
            x_values = np.asarray(
                [float(reaction_coordinate[index]) for index in candidates],
                dtype=float,
            )
            y_values = np.asarray(
                [
                    PEBScanEngine._require_energy(
                        energies[index], "S2 INT stable-point missing energy"
                    )
                    * HARTREE_TO_KCAL
                    for index in candidates
                ],
                dtype=float,
            )
            unique_x = len(np.unique(x_values))
            if unique_x >= 3 and len(candidates) >= 3:
                gradients = np.gradient(y_values, x_values, edge_order=2)
                grad_map = {index: abs(float(grad)) for index, grad in zip(candidates, gradients)}
                selected = min(candidates, key=lambda index: (grad_map[index], abs(reaction_coordinate[index] - reaction_coordinate[m_index])))
                rule = "min_local_gradient"
            else:
                selected = candidates[len(candidates) // 2]
                rule = "plateau_midpoint_fallback"
        return {
            "index": int(selected),
            "frame_xyz": str(Path(frame_paths[selected])),
            "rule": rule,
            "candidate_indices": list(candidates),
            "target_coordinate_A": float(reaction_coordinate[selected]),
        }


    @staticmethod
    def _select_path_nodes(
        engine,
        anchors: Dict[str, Any],
        frame_paths: Sequence[Path],
        reaction_coordinate: Sequence[float],
        method_energies: Sequence[Optional[float]],
        path_metadata: Dict[str, Any],
        off_path_indices: Sequence[int],
        path_arclength: Optional[np.ndarray] = None,
        preferred_int_index: Optional[int] = None,
    ) -> Tuple[Dict[str, Any], Dict[str, Any], int, Optional[int], str]:
        """Pick TS and INT from the refined PATH profile.

        TS strategy:
          1. If xTB PATH returned ``estimated_ts_point``, use that frame index
             directly. xTB's meta-dynamics TS estimate is a converged result.
          2. Otherwise fall back to Kneedle on the refined energy curve.
          3. As a last resort use the plateau_onset anchor M.

        INT strategy:
          - Kneedle ``dipole_distance`` (minimum-gradient point between TS and
            the energy peak). If Kneedle fails, take the frame with the
            smallest local |gradient| on the post-TS plateau; if that also
            fails, return ``intermediate_idx = None`` (DEGRADED).
        """
        from rph_core.utils.scan_profile_plotter import find_ts_and_dipole_guess

        valid_energies = [
            (i, float(method_energies[i]))
            for i in range(len(method_energies))
            if method_energies[i] is not None and i not in set(off_path_indices)
        ]
        ts_index: Optional[int] = None
        ts_rule = "unknown"
        ts_actual_method = "unknown"

        # --- 优先级 1: xTB 日志报告的 estimated TS point ---
        # 正则从 "norm(g) at est. TS, point: 0.01196   5" 解析,
        # 已转换为 0-based index.
        estimated = path_metadata.get("estimated_ts_point")
        ts_guess_xyz = path_metadata.get("ts_guess_xyz")
        if estimated is not None and 0 <= int(estimated) < len(frame_paths):
            ts_index = int(estimated)
            ts_rule = "xtb_path_estimated_ts"
            ts_actual_method = "xtb_path_estimated_ts"
        elif ts_guess_xyz:
            # fallback: 通过结构 RMSD 匹配 xtbpath_ts.xyz
            ts_target = Path(ts_guess_xyz)
            if ts_target.exists():
                from rph_core.utils.file_io import read_xyz
                from rph_core.utils.geometry_tools import kabsch_rmsd
                ts_coords, _ = read_xyz(ts_target)
                best_i, best_rmsd = 0, float("inf")
                for i, frame in enumerate(frame_paths):
                    try:
                        fc, _ = read_xyz(Path(frame))
                        r = kabsch_rmsd(np.asarray(fc, dtype=float),
                                        np.asarray(ts_coords, dtype=float))
                        if r < best_rmsd:
                            best_i, best_rmsd = i, r
                    except Exception:
                        continue
                if best_rmsd < 0.01:  # 0.01 Å 阈值
                    ts_index = best_i
                    ts_rule = "xtb_path_ts_xyz_rmsd"
                    ts_actual_method = "xtb_path_ts_xyz_rmsd"

        # --- 优先级 2: B97-3c 路径最高能点验证 / 微调 ---
        # 当 method_energies 全部非 None（即 B97-3c 完整覆盖），
        # 与 xTB TS 对比并微调
        if ts_index is not None and all(e is not None for e in method_energies):
            full_method_energies = [float(value) for value in method_energies if value is not None]
            nrg_max = max(range(len(full_method_energies)), key=full_method_energies.__getitem__)
            frame_shift = abs(nrg_max - ts_index)
            if frame_shift <= 2:
                ts_index = nrg_max
                ts_rule += "+b973c_validated"
            elif frame_shift > 2:
                engine.logger.warning(
                    "[S2] Energy max (frame %d) differs from xTB TS (frame %d) by %d frames",
                    nrg_max, ts_index, frame_shift,
                )
                ts_rule += "+b973c_disagreement"

        # --- 优先级 3: xTB 能量最高点 ---
        if ts_index is None:
            xb_max = max(valid_energies, key=lambda item: item[1])[0]
            if 0 < xb_max < len(frame_paths) - 1:  # 不在端点
                ts_index = xb_max
                ts_rule = "xtb_energy_maximum"
                ts_actual_method = "xtb_energy_maximum"

        # --- 优先级 4: Kneedle on path_arclength (仅作 fallback) ---
        if ts_index is None and len(valid_energies) >= 3:
            if path_arclength is not None:
                # PATH 模式：使用弧长排序，保持帧序
                arc_x = [float(path_arclength[i]) for i, _ in valid_energies]
            else:
                # 粗扫回退：使用 forming-bond mean（此时单调）
                arc_x = [float(reaction_coordinate[i]) for i, _ in valid_energies]
            energies_kcal = [e * HARTREE_TO_KCAL for _, e in valid_energies]
            ts_distance, _ = find_ts_and_dipole_guess(
                distances=arc_x,
                energies=energies_kcal,
                energies_in_hartree=False,
            )
            if ts_distance is not None:
                nearest_x = arc_x if path_arclength is not None else reaction_coordinate
                ts_index = min(
                    range(len(nearest_x)),
                    key=lambda i: abs(float(nearest_x[i]) - float(ts_distance)),
                )
                ts_rule = "kneedle_ts"
                ts_actual_method = "kneedle_ts"

        # --- 优先级 5: plateau onset (最后备选) ---
        if ts_index is None:
            ts_index = int(anchors["plateau_onset_index"])
            ts_rule = "anchor:plateau_onset_index"
            ts_actual_method = "anchor_plateau_onset_fallback"

        ts_index = max(0, min(int(ts_index), len(frame_paths) - 1))
        ts_selection = {
            "index": int(ts_index),
            "frame_xyz": str(Path(frame_paths[ts_index])),
            "rule": ts_rule,
            "actual_method": ts_actual_method,
            "target_coordinate_A": float(reaction_coordinate[ts_index]),
            "candidate_indices": [int(ts_index)],
            "ts_guess_xyz_source": str(ts_guess_xyz) if ts_guess_xyz else None,
        }

        intermediate_idx: Optional[int] = None
        int_rule = "unavailable"
        int_status = "unavailable"
        int_selection: Dict[str, Any] = {
            "index": None,
            "rule": int_rule,
            "status": int_status,
            "selection_status": int_status,
        }

        post_ts = [
            i for i, _ in valid_energies
            if i > ts_index
        ]
        # --- INT: 安全平台中心（来自粗扫稳定平台） ---
        if preferred_int_index is not None and 0 <= preferred_int_index < len(frame_paths):
            intermediate_idx = preferred_int_index
            # 确保 INT 在 TS 之后
            if intermediate_idx <= ts_index:
                intermediate_idx = min(ts_index + 1, len(frame_paths) - 1)
            int_rule = "coarse_safe_plateau_center"
            int_status = "selected"
        else:
            valid_post_ts = [i for i in post_ts if method_energies[i] is not None]
            if valid_post_ts and len(valid_post_ts) >= 2:
                coords_post = [float(reaction_coordinate[i]) for i in valid_post_ts]
                energies_post = [
                    float(method_energies[i]) * HARTREE_TO_KCAL for i in valid_post_ts
                ]
                sorted_pts = sorted(valid_energies, key=lambda item: reaction_coordinate[item[0]])
                distances_full = [reaction_coordinate[i] for i, _ in sorted_pts]
                energies_full = [e * HARTREE_TO_KCAL for _, e in sorted_pts]
                _, dipole_distance = find_ts_and_dipole_guess(
                    distances=distances_full,
                    energies=energies_full,
                    energies_in_hartree=False,
                )
                if dipole_distance is not None:
                    intermediate_idx = min(
                        valid_post_ts,
                        key=lambda i: abs(float(reaction_coordinate[i]) - float(dipole_distance)),
                    )
                    int_rule = "kneedle_dipole"
                if intermediate_idx is None and len(coords_post) >= 3:
                    if path_arclength is not None:
                        x_for_grad = np.asarray(
                            [float(path_arclength[i]) for i in valid_post_ts], dtype=float
                        )
                    else:
                        x_for_grad = np.asarray(coords_post, dtype=float)
                    gradients = np.gradient(
                        np.asarray(energies_post, dtype=float),
                        x_for_grad,
                        edge_order=2,
                    )
                    local_idx = int(np.argmin(np.abs(gradients)))
                    intermediate_idx = valid_post_ts[local_idx]
                    int_rule = "min_gradient_post_ts"
            if intermediate_idx is not None:
                int_selection = {
                    "index": int(intermediate_idx),
                    "frame_xyz": str(Path(frame_paths[intermediate_idx])),
                    "rule": int_rule,
                    "status": "selected",
                    "selection_status": "selected",
                    "target_coordinate_A": float(reaction_coordinate[intermediate_idx]),
                    "candidate_indices": list(post_ts),
                }

        return ts_selection, int_selection, int(ts_index), intermediate_idx, (
            "selected" if intermediate_idx is not None else "unavailable"
        )

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
        """Select a curve-consistent TS and an S3 INT-search seed.

        A PATH geometry is not itself a stationary point.  TS selection therefore
        uses only extrema on the method-consistent energy curve.  The optional
        INT structure is either a resolvable pre-TS basin or a midpoint geometry
        that S3 may relax without claiming that an intermediate was found.
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
            index: float(method_energies[index]) * HARTREE_TO_KCAL
            for index in valid_indices
        }
        ts_prominence_cutoff = float(
            selection_config.get("ts_min_prominence_kcal_mol", 0.15)
        )
        int_prominence_cutoff = float(
            selection_config.get("int_min_basin_prominence_kcal_mol", 0.50)
        )

        maxima: List[Tuple[int, float]] = []
        weak_maxima: List[Tuple[int, float]] = []
        for position in range(1, len(ordered) - 1):
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
                key=lambda item: (energies_kcal[item[0]], item[1]),
            )
            ts_rule = "refined_curve_local_maximum"
            ts_confidence = "high" if ts_prominence >= 0.50 else "medium"
            ts_candidates = [index for index, _ in maxima]
        elif weak_maxima:
            ts_index, ts_prominence = max(
                weak_maxima,
                key=lambda item: (energies_kcal[item[0]], item[1]),
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
            pre_ts_interior = ordered[1:ts_position]
            if pre_ts_interior:
                left_endpoint = ordered[0]
                midpoint = 0.5 * (coordinates[left_endpoint] + coordinates[ts_index])
                intermediate_idx = min(
                    pre_ts_interior,
                    key=lambda index: abs(coordinates[index] - midpoint),
                )
                int_mode = "midpoint_fallback"
                int_rule = "pre_ts_arclength_midpoint"
                int_candidates = list(pre_ts_interior)
                int_extra = {
                    "left_endpoint_index": int(left_endpoint),
                    "target_path_coordinate": float(midpoint),
                    "reason": "no_resolved_pre_ts_basin",
                }
            else:
                int_mode = "unavailable"
                int_rule = "no_pre_ts_interior_frame"
                int_candidates = []
                int_extra = {"reason": "no_pre_ts_interior_frame"}

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
            "index": None if intermediate_idx is None else int(intermediate_idx),
            "frame_xyz": (
                None if intermediate_idx is None else str(Path(frame_paths[intermediate_idx]))
            ),
            "rule": int_rule,
            "selection_mode": int_mode,
            "selection_status": "selected" if intermediate_idx is not None else "unavailable",
            "stationary_point_claimed": False,
            "target_coordinate_A": (
                None
                if intermediate_idx is None
                else float(reaction_coordinate[intermediate_idx])
            ),
            "candidate_indices": [int(index) for index in int_candidates],
            **int_extra,
        }
        return ts_selection, int_selection, int(ts_index), intermediate_idx, (
            "selected" if intermediate_idx is not None else "unavailable"
        )


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
    ) -> Tuple[ScanAttempt, Dict[str, Any]]:
        """Run xTB2 meta-dynamics PATH and wrap the result as a ScanAttempt.

        Returns (attempt, path_metadata). Raises RuntimeError on PATH failure
        so the caller can fall back to coarse-only selection.
        """
        path_cfg = self._xtb_path_config()
        if not bool(path_cfg.get("enabled", True)):
            raise RuntimeError("S2 xtb_path disabled in config")
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

        path_dir = output_dir / "xtb_path"
        path_dir.mkdir(parents=True, exist_ok=True)

        self._emit_progress(
            "batch_started",
            batch=f"{self.molecule_name or 'product'}:xtb_path",
            phase="xtb_path",
            label="xTB2 meta-dynamics PATH search",
            started_at=time.time(),
            engine="xtb",
            method=f"GFN{int(path_cfg.get('gfn_level', 2))}",
            npoint=int(path_cfg.get("npoint", 28)),
        )
        path_started = time.monotonic()
        try:
            path_result = xtb.path(
                start_xyz=Path(start_xyz),
                end_xyz=Path(end_xyz),
                output_dir=path_dir,
                nrun=int(path_cfg.get("nrun", 1)),
                npoint=int(path_cfg.get("npoint", 28)),
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
                batch=f"{self.molecule_name or 'product'}:xtb_path",
                phase="xtb_path",
                status="failed",
                error=str(exc),
                elapsed_seconds=time.monotonic() - path_started,
            )
            raise RuntimeError(f"S2 xTB PATH execution failed: {exc}") from exc

        if not path_result.success or not path_result.path_xyz_files:
            self._emit_progress(
                "batch_finished",
                batch=f"{self.molecule_name or 'product'}:xtb_path",
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

        attempt = ScanAttempt(
            attempt_id="NN_xtb_path",
            kind="xtb_path",
            directory=path_dir,
            frame_paths=frame_paths,
            target_coordinates_A=coordinates,
            xtb_energies_hartree=energies,
            off_path_indices=(),
            trajectory_quality={
                "xtb_path": True,
                "npoint_requested": int(path_cfg.get("npoint", 28)),
                "npoint_returned": len(frame_paths),
                "barrier_forward_kcal": path_result.barrier_forward_kcal,
                "barrier_backward_kcal": path_result.barrier_backward_kcal,
                "reaction_energy_kcal": path_result.reaction_energy_kcal,
                "gradient_norm_at_ts": path_result.gradient_norm_at_ts,
            },
            scan_policy="xtb_path",
        )
        path_metadata = {
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
            batch=f"{self.molecule_name or 'product'}:xtb_path",
            phase="xtb_path",
            status="complete",
            npoint_returned=len(frame_paths),
            elapsed_seconds=time.monotonic() - path_started,
        )
        return attempt, path_metadata


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

    def run(
        self,
        product_xyz: Path,
        output_dir: Path,
        forming_bonds: Sequence[Tuple[int, int]],
        scan_config: Optional[Dict[str, Any]] = None,
    ) -> Tuple[
        Path,
        Optional[Path],
        Optional[Path],
        Tuple[Tuple[int, int], ...],
        Path,
        str,
        str,
        Tuple[str, ...],
    ]:
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

        refinement_cfg = dict(params.get("candidate_refinement", {}) or {})
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
        product_index = int(coarse_anchors["product_index"])
        if path_start_index == product_index:
            raise RuntimeError(
                "S2 xTB PATH cannot start and end on the same frame "
                f"(index={path_start_index})"
            )
        path_start_frame = active_attempt.frame_paths[path_start_index]
        path_end_frame = active_attempt.frame_paths[product_index]

        path_attempt: Optional[ScanAttempt] = None
        path_metadata: Dict[str, Any] = {}
        try:
            path_attempt, path_metadata = self._execute_xtb_path(
                start_xyz=path_start_frame,
                end_xyz=path_end_frame,
                output_dir=output_dir,
                forming_bonds=bonds,
            )
            attempts.append(path_attempt)
        except RuntimeError as exc:
            self.logger.warning(
                "[S2] xTB PATH failed (%s); falling back to coarse-only anchor selection",
                exc,
            )
            degraded_reasons.append(f"xtb_path_failed:{exc}")
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
            off_path_indices = set()
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
                    "topology_valid": True,
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
                    "topology_valid_point_count": len(composite_points),
                    "complete_xtb_curve": True,
                },
                "continuity_checks": [],
            }
            trajectory_quality = {
                "checked": False,
                "total_frames": len(composite_points),
                "off_path_indices": [],
                "off_path_count": 0,
                "coverage_checks": coverage_checks,
                "source_attempts": [path_attempt.attempt_id],
                "xtb_path": True,
            }
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
        if path_attempt is not None:
            path_end = len(reaction_coordinate) - 1
            anchors["scan_endpoint_index"] = path_end
            if "product_index" not in anchors:
                anchors["product_index"] = path_end

        energy_refinement_cfg = dict(self.step2_cfg.get("energy_refinement", {}) or {})
        all_indices = list(range(len(composite_points)))
        refinement_frames = [frame_paths[index] for index in all_indices]
        refinement_point_ids = [
            str(composite_points[index]["point_id"]) for index in all_indices
        ]
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

        preferred_b973c = selection_cfg.get("preferred_energy_source") == "b973c"
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

        selection_started = time.monotonic()
        self._emit_progress(
            "step_started",
            step="node_selection_render",
            label="TS/INT selection and profile rendering",
            index=3,
            total_steps=3,
            point_count=len(composite_points),
        )

        # --- 从粗扫稳定平台计算 INT 候选 ---
        preferred_int_index: Optional[int] = None
        if path_attempt is not None:
            p_onset = coarse_anchors.get("plateau_onset_index")
            p_drift = coarse_anchors.get("topology_drift_index") or coarse_anchors.get("scan_endpoint_index")
            if p_onset is not None and p_drift is not None:
                p_onset = int(p_onset)
                p_drift = int(p_drift)
                if p_drift > p_onset:
                    plateau_len = p_drift - p_onset + 1
                    margin = max(1, math.ceil(0.20 * plateau_len))
                    safe_start = p_onset + margin
                    safe_end = p_drift - margin
                    if safe_start <= safe_end:
                        # 粗扫安全平台中心坐标
                        coarse_rc = coarse_scan_summary.get("reaction_coordinate_angstrom", [])
                        center_idx = (safe_start + safe_end) // 2
                        if center_idx < len(coarse_rc):
                            target_coord = float(coarse_rc[center_idx])
                            # 映射到 PATH 最近帧
                            preferred_int_index = min(
                                range(len(reaction_coordinate)),
                                key=lambda i: abs(float(reaction_coordinate[i]) - target_coord),
                            )

        ts_selection, intermediate_selection, max_idx, intermediate_idx, intermediate_selection_status = (
            self._select_refined_path_nodes(
                anchors=anchors,
                frame_paths=frame_paths,
                reaction_coordinate=reaction_coordinate,
                method_energies=method_energies,
                off_path_indices=sorted(off_path_indices),
                path_arclength=(
                    np.asarray(path_metadata.get("path_arclength", []), dtype=float)
                    if path_metadata and path_metadata.get("path_arclength")
                    else None
                ),
                selection_config=selection_cfg,
            )
        )
        ts_selection.update(
            {
                "configured_method": "refined_path_curve",
                "actual_method": ts_selection.get("actual_method", "refined_path_curve"),
                "xtb_path_estimated_ts_point": path_metadata.get("estimated_ts_point"),
            }
        )
        ts_guess_confidence = str(ts_selection.get("confidence", ts_guess_confidence))
        if intermediate_idx is None:
            intermediate_selection.setdefault("reason", "no_int_search_seed")

        ts_guess_xyz_final = output_dir / "ts_guess.xyz"
        self._atomic_copy(frame_paths[max_idx], ts_guess_xyz_final)

        dipolar_xyz = output_dir / INTERMEDIATE_XYZ
        if intermediate_idx is not None:
            self._atomic_copy(frame_paths[intermediate_idx], dipolar_xyz)

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
        if intermediate_idx is not None:
            intermediate_energy = method_energies[intermediate_idx]
            if intermediate_energy is None:
                raise RuntimeError(
                    "S2 selected an intermediate without a method-consistent energy"
                )
            intermediate_relative_energy = (
                b973c_relative_energies[intermediate_idx]
                if intermediate_energy_source == "B97-3c"
                else xtb_relative_energies[intermediate_idx]
            )
            int_selection_mode = str(
                intermediate_selection.get("selection_mode", "midpoint_fallback")
            )
            int_method_label = (
                "PRE_TS_BASIN"
                if int_selection_mode == "stable_basin_candidate"
                else "PRE_TS_MIDPOINT"
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
                    "target_distance_angstrom": float(reaction_coordinate[intermediate_idx]),
                    "point_id": composite_points[intermediate_idx]["point_id"],
                    "source_attempt": composite_points[intermediate_idx]["source_attempt"],
                    "energy_hartree": float(intermediate_energy),
                    "energy_source": intermediate_energy_source,
                    "xtb_energy_hartree": float(xtb_energies[intermediate_idx]),
                    "relative_energy_kcal_mol": (
                        float(intermediate_relative_energy)
                        if intermediate_relative_energy is not None
                        else None
                    ),
                    "selection_status": "selected",
                    "output_xyz": str(dipolar_xyz),
                }
            )
        ts_method_label = "PATH_TS"
        ts_selection.update(
            {
                "status": f"B973C_LOCAL_{ts_method_label}_TS_GUESS"
                if ts_selection_energy_source == "B97-3c"
                else f"XTB_LOCAL_{ts_method_label}_TS_GUESS",
                "confirmation_level": "guess_only",
                "validation_required": ["s3_opt_ts", "s3_frequency"],
                "output_xyz": str(ts_guess_xyz_final),
                "target_distance_angstrom": float(reaction_coordinate[max_idx]),
                "point_id": composite_points[max_idx]["point_id"],
                "source_attempt": composite_points[max_idx]["source_attempt"],
                "energy_hartree": self._require_energy(
                    method_energies[max_idx],
                    "S2 selected a TS without a method-consistent energy",
                ),
                "energy_source": ts_selection_energy_source,
                "xtb_energy_hartree": float(xtb_energies[max_idx]),
                "relative_energy_kcal_mol": (
                    self._require_energy(
                        b973c_relative_energies[max_idx],
                        "S2 selected a B97-3c TS without a relative energy",
                    )
                    if ts_selection_energy_source == "B97-3c"
                    and b973c_relative_energies[max_idx] is not None
                    else float(xtb_relative_energies[max_idx])
                ),
            }
        )
        trajectory_quality.update(
            {
                "endpoint_index": int(anchors["scan_endpoint_index"]),
                "endpoint_excluded": True,
                "topology_drift_index": anchors["topology_drift_index"],
                "last_valid_before_drift_index": anchors[
                    "last_valid_before_drift_index"
                ],
            }
        )
        forming_bond_distances = self._forming_bond_distances_by_frame(frame_paths, bonds)
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
            "profile_schema_version": "s2_scan_profile_v9",
            "generation_method": "xtb_path_full_coverage",
            "product_xyz": str(product_file),
            "intermediate_xyz": str(dipolar_xyz) if intermediate_idx is not None else None,
            "forming_bonds": [list(pair) for pair in bonds],
            "scan_parameters": params,
            "coarse_scan": coarse_scan_summary,
            "attempts": attempt_records,
            "composite_profile": composite_profile,
            "anchors": anchor_payload,
            "xtb_path": path_metadata if path_attempt is not None else None,
            "coarse_candidates": {
                "ts_guess_1": coarse_ts_selection,
            },
            "derivative_analysis": {
                "method": ts_selection_energy_source,
                "rule": "xtb_path_or_kneedle",
                "actual_method": ts_selection.get("actual_method"),
                "selected_index": int(max_idx),
            },
            "scan_quality": {
                "ts_index": int(max_idx),
                "plateau_onset_index": int(anchors["plateau_onset_index"]),
                "absolute_energy_maximum_index": int(
                    anchors["absolute_energy_maximum_index"]
                ),
                "selection_energy_source": ts_selection_energy_source,
                "intermediate_energy_source": intermediate_energy_source,
                "intermediate_index": intermediate_idx,
                "intermediate_selection_status": intermediate_selection_status,
                "status": status,
                "ts_guess_confidence": ts_guess_confidence,
                "intermediate_confidence": (
                    "high"
                    if intermediate_idx is not None
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
            ts_index=int(max_idx),
            intermediate_index=intermediate_idx,
            elapsed_seconds=time.monotonic() - selection_started,
        )

        return (
            ts_guess_xyz_final,
            dipolar_xyz if intermediate_idx is not None else None,
            dipolar_xyz if intermediate_idx is not None else None,
            bonds,
            scan_profile_json,
            status,
            ts_guess_confidence,
            tuple(degraded_reasons),
        )
