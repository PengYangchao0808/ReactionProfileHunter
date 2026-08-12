"""Pure structural identity classifiers for TS, intermediates, and minima."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from rph_core.utils.file_io import read_xyz
from rph_core.utils.geometry_tools import GeometryUtils, kabsch_rmsd


_FREQUENCY_MATCH_TOLERANCE_CM1 = 5.0


def classify_ts(
    *,
    frequencies_cm1: Sequence[float],
    forming_bonds: Sequence[Tuple[int, int]],
    normal_modes: Optional[Any] = None,
    imaginary_cutoff_cm1: float = -50.0,
    soft_mode_window_cm1: Tuple[float, float] = (-50.0, -10.0),
    mode_alignment_threshold: float = 0.30,
) -> Dict[str, Any]:
    """Classify a TS candidate by Hessian order, curvature, and mode identity."""

    soft_low, soft_high = _normalize_soft_mode_window(soft_mode_window_cm1, imaginary_cutoff_cm1)
    candidate_imaginaries = [
        float(value) for value in frequencies_cm1 if float(value) < soft_high
    ]
    strict_imaginaries = [
        value for value in candidate_imaginaries if value <= float(imaginary_cutoff_cm1)
    ]
    hessian_index = len(candidate_imaginaries)

    if hessian_index == 1 and (
        strict_imaginaries or candidate_imaginaries[0] <= soft_low
    ):
        curvature_class = "strict"
    elif hessian_index == 1:
        curvature_class = "soft"
    elif hessian_index == 0:
        curvature_class = "none"
    else:
        curvature_class = "multi_imaginary"

    if normal_modes is None or hessian_index != 1:
        mode_identity = "unavailable"
    else:
        mode_identity = _check_mode_alignment(
            normal_modes=normal_modes,
            imaginary_frequency_cm1=candidate_imaginaries[0],
            forming_bonds=forming_bonds,
            threshold=mode_alignment_threshold,
        )

    if hessian_index == 1 and curvature_class == "strict" and mode_identity == "target":
        stationary_point_class = "valid_target_ts"
    elif hessian_index == 1 and curvature_class == "soft" and mode_identity == "target":
        stationary_point_class = "soft_target_ts"
    elif hessian_index == 1 and mode_identity == "unrelated":
        stationary_point_class = "first_order_wrong_mode"
    elif hessian_index >= 2:
        stationary_point_class = "higher_order_saddle"
    elif hessian_index == 0:
        stationary_point_class = "minimum_after_optts"
    else:
        stationary_point_class = "unclassifiable"

    return {
        "hessian_index": hessian_index,
        "curvature_class": curvature_class,
        "mode_identity": mode_identity,
        "stationary_point_class": stationary_point_class,
    }


def classify_int(
    *,
    opt_xyz: Optional[Path],
    forming_bonds: Sequence[Tuple[int, int]],
    precursor_ref: Optional[Path],
    product_ref: Optional[Path],
    atom_mapping: Optional[Dict[int, int]] = None,
    frequencies_cm1: Sequence[float] = (),
    imaginary_cutoff_cm1: float = -10.0,
    rmsd_threshold: float = 0.3,
    progress_window: Tuple[float, float] = (0.2, 0.8),
    coordinates: Optional[Sequence[Tuple[float, float, float]]] = None,
    precursor_coordinates: Optional[Sequence[Tuple[float, float, float]]] = None,
    product_coordinates: Optional[Sequence[Tuple[float, float, float]]] = None,
) -> Dict[str, Any]:
    """Classify an intermediate using frequency, progress, and mapped RMSD."""

    current_coords = _resolve_coordinates(coordinates=coordinates, xyz_path=opt_xyz)
    if current_coords is None:
        return {
            "identity": "opt_failed",
            "avg_progress": None,
            "per_bond_progress": {},
            "per_bond_distance": {},
            "rmsd_to_product": None,
            "rmsd_to_precursor": None,
        }

    significant_imaginaries = [
        float(value) for value in frequencies_cm1 if float(value) < imaginary_cutoff_cm1
    ]
    if significant_imaginaries:
        return {
            "identity": "imaginary_frequency",
            "avg_progress": None,
            "per_bond_progress": {},
            "per_bond_distance": {},
            "rmsd_to_product": None,
            "rmsd_to_precursor": None,
            "imag_count": len(significant_imaginaries),
        }

    precursor_coords = _resolve_coordinates(
        coordinates=precursor_coordinates,
        xyz_path=precursor_ref,
    )
    product_coords = _resolve_coordinates(
        coordinates=product_coordinates,
        xyz_path=product_ref,
    )

    per_bond_progress: Dict[Tuple[int, int], float] = {}
    per_bond_distance: Dict[Tuple[int, int], float] = {}

    for raw_atom_i, raw_atom_j in forming_bonds:
        atom_i = int(raw_atom_i)
        atom_j = int(raw_atom_j)
        bond = (atom_i, atom_j)
        try:
            current_distance = _measure_distance(current_coords, atom_i, atom_j)
        except (IndexError, TypeError, ValueError):
            continue

        per_bond_distance[bond] = current_distance

        if precursor_coords is None or product_coords is None:
            continue

        try:
            precursor_distance = _measure_distance(precursor_coords, atom_i, atom_j)
            product_distance = _measure_distance(product_coords, atom_i, atom_j)
        except (IndexError, TypeError, ValueError):
            continue

        denominator = precursor_distance - product_distance
        progress = 0.5 if abs(denominator) < 1e-12 else (precursor_distance - current_distance) / denominator
        per_bond_progress[bond] = float(progress)

    avg_progress = (
        float(sum(per_bond_progress.values()) / len(per_bond_progress))
        if per_bond_progress
        else 0.5
    )

    rmsd_to_product = _safe_compute_rmsd(current_coords, product_coords, atom_mapping)
    rmsd_to_precursor = _safe_compute_rmsd(current_coords, precursor_coords, atom_mapping)

    if rmsd_to_product is not None and rmsd_to_product < rmsd_threshold:
        identity = "collapsed_to_product"
    elif rmsd_to_precursor is not None and rmsd_to_precursor < rmsd_threshold:
        identity = "collapsed_to_precursor"
    elif float(progress_window[0]) < avg_progress < float(progress_window[1]):
        identity = "distinct_intermediate"
    else:
        identity = "topology_ambiguous"

    return {
        "identity": identity,
        "avg_progress": avg_progress,
        "per_bond_progress": per_bond_progress,
        "per_bond_distance": per_bond_distance,
        "rmsd_to_product": rmsd_to_product,
        "rmsd_to_precursor": rmsd_to_precursor,
    }


_INT_V2_DEFAULTS: Dict[str, Any] = {
    "classification_version": "int_identity_v2",
    "imaginary_cutoff_cm1": -10.0,
    "product_rmsd_collapse_ang": 0.30,
    "formed_bond_max_ang": 1.75,
    "dipolar_open_bond_min_ang": 2.40,
    "dipolar_open_bond_max_ang": 4.60,
    "ts_rmsd_merged_ang": 0.40,
    "ts_energy_merged_kcal": 0.5,
    "ts_energy_degenerate_kcal": 0.3,
    "int_below_ts_min_kcal": 0.2,
    "int_above_product_min_kcal": 2.0,
    "above_ts_invalid_kcal": 0.5,
    "product_like_energy_max_kcal": 2.0,
}

_KCAL_PER_HARTREE = 627.509


class IntIdentityResult:
    """Structured INT identity classification (v2)."""

    __slots__ = (
        "identity",
        "usable_for_ml",
        "classification_version",
        "reason_codes",
        "metrics",
        "missing_evidence",
    )

    def __init__(
        self,
        *,
        identity: str,
        usable_for_ml: bool,
        classification_version: str = "int_identity_v2",
        reason_codes: Optional[Sequence[str]] = None,
        metrics: Optional[Mapping[str, Any]] = None,
        missing_evidence: Optional[Sequence[str]] = None,
    ) -> None:
        self.identity = identity
        self.usable_for_ml = usable_for_ml
        self.classification_version = classification_version
        self.reason_codes = list(reason_codes or [])
        self.metrics = dict(metrics or {})
        self.missing_evidence = list(missing_evidence or [])

    def as_dict(self) -> Dict[str, Any]:
        return {
            "identity": self.identity,
            "usable_for_ml": self.usable_for_ml,
            "classification_version": self.classification_version,
            "reason_codes": self.reason_codes,
            "metrics": dict(self.metrics),
            "missing_evidence": self.missing_evidence,
        }


def classify_int_v2(
    *,
    forming_bonds: Sequence[Tuple[int, int]],
    opt_xyz: Optional[Path] = None,
    product_ref: Optional[Path] = None,
    ts_ref: Optional[Path] = None,
    precursor_ref: Optional[Path] = None,
    frequencies_cm1: Sequence[float] = (),
    energy_int_hartree: Optional[float] = None,
    energy_ts_hartree: Optional[float] = None,
    energy_product_hartree: Optional[float] = None,
    atom_mapping: Optional[Dict[int, int]] = None,
    converged: Optional[bool] = None,
    thresholds: Optional[Mapping[str, Any]] = None,
    coordinates: Optional[Sequence[Tuple[float, float, float]]] = None,
    product_coordinates: Optional[Sequence[Tuple[float, float, float]]] = None,
    ts_coordinates: Optional[Sequence[Tuple[float, float, float]]] = None,
    precursor_coordinates: Optional[Sequence[Tuple[float, float, float]]] = None,
) -> IntIdentityResult:
    """Evidence-based intermediate classification.

    The legacy ``classify_int`` silently defaulted to ``distinct_intermediate``
    when reference geometries were unavailable (avg_progress fell back to 0.5).
    v2 requires explicit evidence: product-side forming-bond distances, mapped
    RMSD to the sibling product, sibling-TS proximity, and energy ordering
    (TS >= INT > product).  Precursor-side reaction progress is intentionally
    NOT used: precursor and product atom orderings differ, so raw-index
    distances on the precursor are invalid.
    """
    thresholds = {**_INT_V2_DEFAULTS, **dict(thresholds or {})}
    missing_evidence: List[str] = []
    metrics: Dict[str, Any] = {}

    current = _resolve_coordinates(coordinates=coordinates, xyz_path=opt_xyz)
    product = _resolve_coordinates(
        coordinates=product_coordinates, xyz_path=product_ref
    )
    ts = _resolve_coordinates(coordinates=ts_coordinates, xyz_path=ts_ref)

    def _finish(
        identity: str,
        usable: bool,
        reason_codes: Sequence[str],
        extra_metrics: Optional[Mapping[str, Any]] = None,
    ) -> IntIdentityResult:
        payload = dict(metrics)
        if extra_metrics:
            payload.update(extra_metrics)
        return IntIdentityResult(
            identity=identity,
            usable_for_ml=usable,
            classification_version=str(thresholds["classification_version"]),
            reason_codes=reason_codes,
            metrics=payload,
            missing_evidence=missing_evidence,
        )

    if current is None:
        return _finish("unclassified", False, ["opt_geometry_missing"])
    if product is not None and product.shape != current.shape:
        return _finish(
            "unclassified",
            False,
            ["atom_count_mismatch"],
            {"current_atoms": int(current.shape[0]), "product_atoms": int(product.shape[0])},
        )
    if converged is False:
        return _finish("unclassified", False, ["opt_not_converged"])
    if energy_int_hartree is None:
        missing_evidence.append("energy_int_hartree")
    if product is None:
        missing_evidence.append("product_ref")
    if energy_product_hartree is None:
        missing_evidence.append("energy_product_hartree")
    if ts is None:
        missing_evidence.append("ts_ref")
    if energy_ts_hartree is None:
        missing_evidence.append("energy_ts_hartree")

    significant_imaginaries = [
        float(value)
        for value in frequencies_cm1
        if float(value) < float(thresholds["imaginary_cutoff_cm1"])
    ]
    if significant_imaginaries:
        return _finish(
            "imaginary_frequency",
            False,
            ["significant_imaginary"],
            {
                "imaginary_count": len(significant_imaginaries),
                "imaginary_frequencies_cm1": significant_imaginaries,
            },
        )

    metrics["current_atoms"] = int(current.shape[0])

    fb_dist_int: Dict[Tuple[int, int], float] = {}
    for raw_atom_i, raw_atom_j in forming_bonds:
        atom_i, atom_j = int(raw_atom_i), int(raw_atom_j)
        try:
            fb_dist_int[(atom_i, atom_j)] = _measure_distance(current, atom_i, atom_j)
        except (IndexError, TypeError, ValueError):
            continue
    metrics["per_bond_distance"] = {
        str(pair): round(distance, 4) for pair, distance in fb_dist_int.items()
    }

    fb_dist_product: Dict[Tuple[int, int], float] = {}
    if product is not None:
        for pair in fb_dist_int:
            atom_i, atom_j = pair
            try:
                fb_dist_product[pair] = _measure_distance(product, atom_i, atom_j)
            except (IndexError, TypeError, ValueError):
                continue
    metrics["product_fb_distances"] = {
        str(pair): round(distance, 4) for pair, distance in fb_dist_product.items()
    }

    fb_dist_ts: Dict[Tuple[int, int], float] = {}
    if ts is not None:
        for pair in fb_dist_int:
            atom_i, atom_j = pair
            try:
                fb_dist_ts[pair] = _measure_distance(ts, atom_i, atom_j)
            except (IndexError, TypeError, ValueError):
                continue
    metrics["ts_fb_distances"] = {
        str(pair): round(distance, 4) for pair, distance in fb_dist_ts.items()
    }

    if fb_dist_int:
        formed = all(
            distance <= float(thresholds["formed_bond_max_ang"])
            for distance in fb_dist_int.values()
        )
        open_dipolar = all(
            float(thresholds["dipolar_open_bond_min_ang"])
            <= distance
            <= float(thresholds["dipolar_open_bond_max_ang"])
            for distance in fb_dist_int.values()
        )
        bond_state = "formed" if formed else ("open_dipolar" if open_dipolar else "mixed_or_borderline")
    else:
        bond_state = "unmeasured"
    metrics["bond_state"] = bond_state

    rmsd_product = _safe_compute_rmsd(current, product, atom_mapping)
    rmsd_ts = _safe_compute_rmsd(current, ts, atom_mapping)
    metrics["rmsd_to_product"] = (
        round(rmsd_product, 4) if rmsd_product is not None else None
    )
    metrics["rmsd_to_ts"] = round(rmsd_ts, 4) if rmsd_ts is not None else None

    e_gap_product_kcal = None
    if energy_int_hartree is not None and energy_product_hartree is not None:
        e_gap_product_kcal = (energy_int_hartree - energy_product_hartree) * _KCAL_PER_HARTREE
        metrics["int_minus_product_kcal"] = round(e_gap_product_kcal, 3)
    e_gap_ts_kcal = None
    if energy_int_hartree is not None and energy_ts_hartree is not None:
        e_gap_ts_kcal = (energy_int_hartree - energy_ts_hartree) * _KCAL_PER_HARTREE
        metrics["int_minus_ts_kcal"] = round(e_gap_ts_kcal, 3)

    product_like_geometry = bool(
        (rmsd_product is not None and rmsd_product <= float(thresholds["product_rmsd_collapse_ang"]))
        or bond_state == "formed"
    )
    if product_like_geometry:
        if (
            e_gap_product_kcal is not None
            and e_gap_product_kcal <= float(thresholds["product_like_energy_max_kcal"])
        ):
            return _finish("collapsed_to_product", False, ["product_collapse"])
        return _finish(
            "product_like_geometry_energy_inconsistent",
            False,
            ["product_like_geometry", "energy_not_product_like"],
        )

    ts_proximity = (
        rmsd_ts is not None
        and rmsd_ts <= float(thresholds["ts_rmsd_merged_ang"])
    )
    ts_energy_degenerate = (
        e_gap_ts_kcal is not None
        and abs(e_gap_ts_kcal) <= float(thresholds["ts_energy_merged_kcal"])
    )
    if ts_proximity and ts_energy_degenerate:
        return _finish(
            "merged_with_ts",
            False,
            ["ts_proximity", "ts_energy_degenerate"],
        )
    if (
        ts is None
        and e_gap_ts_kcal is not None
        and abs(e_gap_ts_kcal) <= float(thresholds["ts_energy_degenerate_kcal"])
    ):
        return _finish("ts_energy_degenerate", False, ["ts_geometry_missing"])

    if e_gap_ts_kcal is not None and e_gap_ts_kcal > float(thresholds["above_ts_invalid_kcal"]):
        return _finish("above_ts_not_minimum", False, ["energy_order_invalid"])

    if ts is None or e_gap_ts_kcal is None:
        if (
            e_gap_product_kcal is not None
            and e_gap_product_kcal >= float(thresholds["int_above_product_min_kcal"])
            and bond_state in {"open_dipolar", "mixed_or_borderline", "unmeasured"}
        ):
            return _finish(
                "candidate_intermediate_unverified_energy",
                False,
                ["ts_evidence_missing"],
            )
        return _finish("unclassified", False, ["insufficient_evidence"])

    well_depth_kcal = -e_gap_ts_kcal
    if (
        well_depth_kcal < -float(thresholds["int_below_ts_min_kcal"])
        or e_gap_product_kcal is None
        or e_gap_product_kcal < float(thresholds["int_above_product_min_kcal"])
        or bond_state not in {"open_dipolar", "mixed_or_borderline", "unmeasured"}
    ):
        reason_codes = []
        if well_depth_kcal < -float(thresholds["int_below_ts_min_kcal"]):
            reason_codes.append("above_ts")
        if e_gap_product_kcal is not None and e_gap_product_kcal < float(
            thresholds["int_above_product_min_kcal"]
        ):
            reason_codes.append("too_close_to_product_energy")
        if bond_state not in {"open_dipolar", "mixed_or_borderline", "unmeasured"}:
            reason_codes.append(f"bond_state={bond_state}")
        return _finish("topology_ambiguous", False, reason_codes or ["ambiguous"])

    return _finish(
        "dipolar_intermediate",
        True,
        ["valid_dipolar_well"],
        {
            "well_depth_kcal": round(well_depth_kcal, 3),
            "above_product_kcal": round(e_gap_product_kcal, 3),
        },
    )


def classify_minimum(
    *,
    frequencies_cm1: Sequence[float],
    imaginary_cutoff_cm1: float = -10.0,
    expected_role: Optional[str] = None,
    rmsd_to_expected: Optional[float] = None,
    rmsd_threshold: float = 0.3,
) -> Dict[str, Any]:
    """Validate a precursor/product minimum for frequencies and identity drift."""

    significant_imaginaries = [
        float(value) for value in frequencies_cm1 if float(value) < imaginary_cutoff_cm1
    ]
    if significant_imaginaries:
        identity = "imaginary_frequency"
        identity_status = "not_checked"
    elif rmsd_to_expected is not None and rmsd_to_expected > rmsd_threshold:
        identity = "identity_drift"
        identity_status = "mismatched"
    elif expected_role is not None and rmsd_to_expected is not None:
        identity = "valid_minimum"
        identity_status = "role_matched"
    else:
        identity = "not_checked"
        identity_status = "not_checked"

    return {
        "identity": identity,
        "imaginary_count": len(significant_imaginaries),
        "significant_imaginary_frequencies_cm1": significant_imaginaries,
        "identity_status": identity_status,
    }


def _check_mode_alignment(
    *,
    normal_modes: Any,
    imaginary_frequency_cm1: float,
    forming_bonds: Sequence[Tuple[int, int]],
    threshold: float,
) -> str:
    coordinates, displacements = _extract_mode_payload(normal_modes, imaginary_frequency_cm1)
    if coordinates is None or displacements is None:
        return "unavailable"

    geometry = _as_coordinate_array(coordinates)
    mode = _as_coordinate_array(displacements)
    if geometry.shape != mode.shape or geometry.size == 0:
        return "unavailable"

    projections = []
    for raw_atom_i, raw_atom_j in forming_bonds:
        atom_i = int(raw_atom_i)
        atom_j = int(raw_atom_j)
        if atom_i < 0 or atom_j < 0 or atom_i >= len(geometry) or atom_j >= len(geometry):
            return "unavailable"

        bond_vector = geometry[atom_j] - geometry[atom_i]
        norm = float(np.linalg.norm(bond_vector))
        if norm <= 1e-12:
            return "unavailable"

        relative_displacement = mode[atom_j] - mode[atom_i]
        projection = abs(float(np.dot(relative_displacement, bond_vector / norm)))
        projections.append(projection)

    if not projections:
        return "unavailable"
    return "target" if min(projections) >= float(threshold) else "unrelated"


def _extract_mode_payload(
    normal_modes: Any,
    target_frequency_cm1: float,
) -> Tuple[Optional[Any], Optional[Any]]:
    if isinstance(normal_modes, Mapping):
        coordinates = normal_modes.get("coordinates") or normal_modes.get("geometry")
        direct_displacements = (
            normal_modes.get("displacements")
            or normal_modes.get("atomic_displacements")
        )
        if coordinates is not None and direct_displacements is not None:
            return coordinates, direct_displacements

        modes = normal_modes.get("modes")
        if isinstance(modes, Mapping):
            for raw_frequency, payload in modes.items():
                frequency = _coerce_float(raw_frequency)
                if frequency is None or not _frequency_matches(frequency, target_frequency_cm1):
                    continue
                if isinstance(payload, Mapping):
                    return (
                        payload.get("coordinates") or payload.get("geometry") or coordinates,
                        payload.get("displacements") or payload.get("atomic_displacements"),
                    )
                return coordinates, payload

        if isinstance(modes, Sequence) and not isinstance(modes, (str, bytes)):
            for payload in modes:
                if not isinstance(payload, Mapping):
                    continue
                frequency = _coerce_float(
                    payload.get("frequency_cm1", payload.get("frequency"))
                )
                if frequency is None or not _frequency_matches(frequency, target_frequency_cm1):
                    continue
                return (
                    payload.get("coordinates") or payload.get("geometry") or coordinates,
                    payload.get("displacements") or payload.get("atomic_displacements"),
                )

    return None, None


def _measure_distance(
    coords: Sequence[Tuple[float, float, float]] | np.ndarray,
    atom_i: int,
    atom_j: int,
) -> float:
    return float(GeometryUtils.calculate_distance(_as_coordinate_array(coords), atom_i, atom_j))


def _compute_mapped_rmsd(
    coords_a: Sequence[Tuple[float, float, float]] | np.ndarray,
    coords_b: Sequence[Tuple[float, float, float]] | np.ndarray,
    atom_mapping: Optional[Dict[int, int]],
) -> float:
    array_a = _as_coordinate_array(coords_a)
    array_b = _as_coordinate_array(coords_b)

    if atom_mapping is None:
        if array_a.shape != array_b.shape:
            raise ValueError("Coordinate arrays must match when no atom mapping is provided")
        return float(kabsch_rmsd(array_a, array_b))

    ordered_pairs = sorted((int(source), int(target)) for source, target in atom_mapping.items())
    if not ordered_pairs:
        raise ValueError("Atom mapping cannot be empty")

    subset_a = []
    subset_b = []
    for source_index, target_index in ordered_pairs:
        if source_index < 0 or target_index < 0:
            raise ValueError("Atom mapping indices must be non-negative")
        if source_index >= len(array_a) or target_index >= len(array_b):
            raise ValueError("Atom mapping index out of range")
        subset_a.append(array_a[source_index])
        subset_b.append(array_b[target_index])

    return float(kabsch_rmsd(np.asarray(subset_a, dtype=float), np.asarray(subset_b, dtype=float)))


def _safe_compute_rmsd(
    current_coords: np.ndarray,
    reference_coords: Optional[np.ndarray],
    atom_mapping: Optional[Dict[int, int]],
) -> Optional[float]:
    if reference_coords is None:
        return None
    try:
        return _compute_mapped_rmsd(current_coords, reference_coords, atom_mapping)
    except ValueError:
        return None


def _resolve_coordinates(
    *,
    coordinates: Optional[Sequence[Tuple[float, float, float]]],
    xyz_path: Optional[Path],
) -> Optional[np.ndarray]:
    if coordinates is not None:
        return _as_coordinate_array(coordinates)
    if xyz_path is None:
        return None
    try:
        loaded_coordinates, _symbols = read_xyz(Path(xyz_path))
    except (FileNotFoundError, OSError, ValueError):
        return None
    return np.asarray(loaded_coordinates, dtype=float)


def _as_coordinate_array(
    coordinates: Sequence[Tuple[float, float, float]] | np.ndarray,
) -> np.ndarray:
    array = np.asarray(coordinates, dtype=float)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError("Coordinates must have shape (N, 3)")
    return array


def _normalize_soft_mode_window(
    soft_mode_window_cm1: Sequence[float],
    strict_cutoff_cm1: float,
) -> Tuple[float, float]:
    raw_window = tuple(soft_mode_window_cm1) if soft_mode_window_cm1 else (strict_cutoff_cm1, -10.0)
    if len(raw_window) != 2:
        raise ValueError("soft_mode_window_cm1 must contain exactly two values")

    soft_low = float(raw_window[0])
    soft_high = float(raw_window[1])
    if soft_low > soft_high:
        soft_low, soft_high = soft_high, soft_low
    return soft_low, soft_high


def _coerce_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _frequency_matches(value: float, target_frequency_cm1: float) -> bool:
    return abs(float(value) - float(target_frequency_cm1)) <= _FREQUENCY_MATCH_TOLERANCE_CM1


__all__ = [
    "classify_ts",
    "classify_int",
    "classify_int_v2",
    "classify_minimum",
    "IntIdentityResult",
    "_measure_distance",
    "_compute_mapped_rmsd",
]
