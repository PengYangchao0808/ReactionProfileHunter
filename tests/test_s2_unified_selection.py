from pathlib import Path

import pytest

from rph_core.steps.step2_retro import path_profile as path_profile_module
from rph_core.steps.step2_retro.path_profile import (
    build_orca_scan_profile,
    build_xtb_path_profile,
    compute_forming_bond_distances_by_frame,
    compute_neighbor_rmsds,
    compute_path_arclength,
    scaffold_rmsd_admission,
)
from rph_core.steps.step2_retro.path_selector import (
    SelectionPolicy,
    policy_from_config,
    select_path_seeds,
)
from rph_core.utils.qc_jobs import resolve_surface_scan_spec
from rph_core.utils.qc_models import SurfaceScanCoordinate, SurfaceScanSpec


def _write_xyz(path: Path, distance: float) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"2\nframe\nH 0 0 0\nH {distance:.6f} 0 0\n",
        encoding="utf-8",
    )
    return path


def _hartree_profile(relative_kcal, reference_hartree=-100.0):
    return [
        None
        if value is None
        else reference_hartree + float(value) / 627.509
        for value in relative_kcal
    ]


def _patch_topology(monkeypatch, *, off_path=(), reasons=None):
    reason_map = dict(reasons or {})

    def fake_check_scan_trajectory(*, frame_paths, **_kwargs):
        frame_issues = [
            {
                "frame_index": int(index),
                "reason": reason_map.get(int(index), "topology_invalid"),
            }
            for index in off_path
        ]
        return {
            "checked": True,
            "off_path_indices": tuple(int(index) for index in off_path),
            "frame_issues": frame_issues,
        }

    monkeypatch.setattr(
        path_profile_module,
        "check_scan_trajectory",
        fake_check_scan_trajectory,
    )


def _xtb_profile(tmp_path: Path, *, relative_kcal, off_path=()):
    distances = [3.4 - 0.2 * index for index in range(9)]
    product = _write_xyz(tmp_path / "xtb_product.xyz", distances[-1])
    frames = [
        _write_xyz(tmp_path / "xtb_frames" / f"frame_{index:03d}.xyz", distance)
        for index, distance in enumerate(distances)
    ]
    profile = build_xtb_path_profile(
        frame_paths=frames,
        energies_hartree=_hartree_profile(relative_kcal),
        forming_bonds=[(0, 1)],
        product_xyz=product,
        off_path_indices=off_path,
        source_provenance={"case": "xtb"},
    )
    return profile, frames


def _orca_profile(
    tmp_path: Path,
    *,
    relative_kcal=None,
    energies_hartree=None,
    off_path=(),
):
    distances = [1.5 + 0.2 * index for index in range(11)]
    product = _write_xyz(tmp_path / "orca_product.xyz", distances[0])
    frames = [
        _write_xyz(tmp_path / "orca_frames" / f"frame_{index:03d}.xyz", distance)
        for index, distance in enumerate(distances)
    ]
    profile = build_orca_scan_profile(
        frames=frames,
        energies_hartree=(
            _hartree_profile(relative_kcal)
            if energies_hartree is None
            else list(energies_hartree)
        ),
        forming_bonds=[(0, 1)],
        product_xyz=product,
        energy_source="orca.relaxscanact.dat",
        scan_ts_candidate_xyz=frames[6],
        source_provenance={"case": "orca", "off_path": list(off_path)},
    )
    return profile, frames


def test_path_profile_geometry_helpers_follow_two_atom_distance_path(tmp_path):
    frames = [
        _write_xyz(tmp_path / f"frame_{index:03d}.xyz", distance)
        for index, distance in enumerate([1.0, 1.2, 1.6])
    ]

    distances = compute_forming_bond_distances_by_frame(frames, [(0, 1)])
    arclength = compute_path_arclength(frames)
    neighbor_rmsd = compute_neighbor_rmsds(frames)
    scaffold_gate = scaffold_rmsd_admission(frames[0], frames[1], [(0, 1)], 0.1)

    assert distances == [[1.0], [1.2], [1.6]]
    assert list(arclength) == pytest.approx([0.0, 0.1, 0.3])
    assert neighbor_rmsd == pytest.approx([0.1, 0.2, 0.2])
    assert scaffold_gate["accepted"] is True
    assert scaffold_gate["reason"] == "insufficient_nonreactive_atoms"


def test_clean_xtb_path_selects_path_seed_without_rescue(tmp_path, monkeypatch):
    _patch_topology(monkeypatch)
    profile, frames = _xtb_profile(
        tmp_path,
        relative_kcal=[0.0, 1.0, 2.0, 1.0, 0.2, 4.2, 1.5, 0.5, 0.0],
    )

    selection = select_path_seeds(profile, SelectionPolicy())

    assert selection.s2_state == "path_seeded"
    assert selection.seed_evidence == "knee_shifted"
    assert selection.ts_search_seed["frame_index"] == 6
    assert selection.ts_search_seed["knee_frame_index"] == 5
    assert selection.ts_search_seed["stationary_point_claimed"] is False
    assert selection.int_search_seed["frame_index"] == 7
    assert selection.int_search_seed["selection_mode"] == (
        "ts_to_effective_endpoint_midpoint"
    )
    assert selection.has_independent_int is False


def test_xtb_shared_search_seed_stays_path_seeded(tmp_path, monkeypatch):
    _patch_topology(monkeypatch)
    profile, frames = _xtb_profile(
        tmp_path,
        relative_kcal=[0.0, 0.2, 0.4, 0.8, 1.2, 4.0, 1.5, 0.5, 0.0],
    )

    selection = select_path_seeds(profile, SelectionPolicy())

    assert selection.s2_state == "path_seeded"
    assert selection.seed_evidence == "knee_shifted"
    assert selection.ts_search_seed["frame_index"] == 6
    assert selection.int_search_seed["frame_index"] == 7
    assert selection.int_search_seed["selection_mode"] == (
        "ts_to_effective_endpoint_midpoint"
    )
    assert selection.has_independent_int is False


def test_orca_internal_peak_rescue_seeds_peak_frame(tmp_path, monkeypatch):
    _patch_topology(monkeypatch)
    profile, frames = _orca_profile(
        tmp_path,
        relative_kcal=[1.0, 1.2, 1.5, 1.7, 1.0, 2.0, 5.5, 2.0, 0.0, 0.2, 0.0],
    )

    selection = select_path_seeds(profile, SelectionPolicy())

    assert selection.s2_state == "rescue_seeded"
    assert selection.seed_evidence == "knee_shifted"
    assert selection.ts_search_seed["frame_index"] > selection.knee_evidence["frame_index"]
    assert selection.int_search_seed["frame_index"] > selection.ts_search_seed["frame_index"]


def test_orca_monotonic_shoulder_can_seed_rescue(tmp_path, monkeypatch):
    _patch_topology(monkeypatch)
    profile, frames = _orca_profile(
        tmp_path,
        relative_kcal=[5.0, 4.8, 4.5, 4.35, 4.3, 4.2, 3.8, 1.5, 0.0, 0.1, 0.0],
    )
    policy = policy_from_config(
        {"shoulder_max_abs_slope_kcal_mol_per_A": 5.0},
        {},
    )

    selection = select_path_seeds(profile, policy)

    assert selection.s2_state == "rescue_seeded"
    assert selection.seed_evidence == "knee_shifted"
    assert selection.ts_search_seed["frame_index"] > selection.knee_evidence["frame_index"]
    assert selection.int_search_seed["frame_index"] > selection.ts_search_seed["frame_index"]


def test_orca_strictly_monotonic_profile_stays_unresolved(tmp_path, monkeypatch):
    _patch_topology(monkeypatch)
    profile, _frames = _orca_profile(
        tmp_path,
        relative_kcal=[7.0, 6.0, 5.0, 4.5, 4.0, 3.0, 2.0, 1.0, 0.0, 0.1, 0.0],
    )
    policy = policy_from_config(
        {"shoulder_max_abs_slope_kcal_mol_per_A": 5.0},
        {},
    )

    selection = select_path_seeds(profile, policy)

    assert selection.s2_state == "unresolved"
    assert selection.seed_evidence == "none"
    assert selection.ts_search_seed is None
    assert selection.int_search_seed is None
    assert selection.rejection_reason == "no_valid_knee_point"


def test_orca_endpoint_highest_energy_is_excluded_from_ts_seeding(tmp_path, monkeypatch):
    _patch_topology(monkeypatch)
    profile, _frames = _orca_profile(
        tmp_path,
        relative_kcal=[1.0, 1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 40.0, 50.0],
    )
    policy = policy_from_config(
        {"shoulder_max_abs_slope_kcal_mol_per_A": 5.0},
        {},
    )

    selection = select_path_seeds(profile, policy)

    assert selection.s2_state == "unresolved"
    assert selection.ts_search_seed is None
    assert selection.rejection_reason == "no_valid_knee_point"
    assert 9 not in selection.diagnostics["filtered_frame_indices"]
    assert 10 not in selection.diagnostics["filtered_frame_indices"]


@pytest.mark.parametrize(
    ("name", "energies_hartree", "off_path", "expected_rejection"),
    [
        (
            "incomplete_ledger",
            _hartree_profile([1.0, 1.2, 1.5, 1.7, 1.0, None, 3.8, 1.5, 0.0, 0.1, 0.0]),
            (),
            "incomplete_profile",
        ),
        (
            "frame_count_mismatch",
            _hartree_profile([1.0, 1.2, 1.5, 1.7, 1.0, 0.0]),
            (),
            "incomplete_profile",
        ),
        (
            "topology_drift",
            _hartree_profile([1.0, 1.2, 1.5, 1.7, 1.0, 2.0, 5.5, 2.0, 0.0, 0.2, 0.0]),
            (4, 5, 6, 7, 8),
            "no_valid_frames_after_filters",
        ),
    ],
)
def test_orca_unresolved_profiles_preserve_selector_diagnostics(
    tmp_path,
    monkeypatch,
    name,
    energies_hartree,
    off_path,
    expected_rejection,
):
    _patch_topology(
        monkeypatch,
        off_path=off_path,
        reasons={index: f"{name}_excluded" for index in off_path},
    )
    profile, _frames = _orca_profile(
        tmp_path / name,
        energies_hartree=energies_hartree,
        off_path=off_path,
    )
    policy = policy_from_config(
        {"shoulder_max_abs_slope_kcal_mol_per_A": 5.0},
        {},
    )

    selection = select_path_seeds(profile, policy)

    assert selection.s2_state == "unresolved"
    assert selection.rejection_reason == expected_rejection
    assert selection.diagnostics["profile_source"] == "orca_relaxed_scan"
    assert selection.diagnostics["frame_count"] == 11
    if name in {"incomplete_ledger", "frame_count_mismatch"}:
        assert profile.complete is False
        assert selection.diagnostics["complete"] is False
    if name == "topology_drift":
        assert profile.excluded_frames == (4, 5, 6, 7, 8)
        assert selection.diagnostics["excluded_frame_indices"] == [4, 5, 6, 7, 8]


def test_resolve_surface_scan_spec_inherits_shared_s3_orca_defaults():
    spec = SurfaceScanSpec(
        method="",
        charge=0,
        multiplicity=1,
        coordinates=(
            SurfaceScanCoordinate(
                kind="B",
                atoms=(0, 1),
                start=1.5,
                end=3.1,
                steps=17,
            ),
        ),
        simultaneous=False,
        scan_ts=True,
        full_scan=True,
    )

    resolved = resolve_surface_scan_spec(
        spec,
        {
            "resources": {"nproc": 6},
            "theory": {
                "s3_low_level": {
                    "optimization": {
                        "method": "B97-3c",
                        "solvent": "acetone",
                        "solvent_model": "CPCM",
                        "max_cycles": 88,
                        "timeout": 12345,
                    }
                }
            },
        },
    )

    assert resolved.method == "B97-3c"
    assert resolved.solvent == "acetone"
    assert resolved.solvent_model == "CPCM"
    assert resolved.nproc == 6
    assert resolved.max_cycles == 88
    assert resolved.timeout == 12345
    assert resolved.scan_ts is True
    assert resolved.coordinates == spec.coordinates


def test_plateau_seed_selected_when_stretch_side_flat(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A flat low-energy plateau on the stretched side wins over the midpoint."""
    from rph_core.steps.step2_retro.path_selector import (
        _detect_stretch_plateau,
        policy_from_config,
        select_path_seeds,
    )

    # 11 frames: TS peak at frame 5-6, plateau at frames 7-9 (stretched side).
    relative_kcal = [
        -8.0, -6.0, -4.0, -2.0, 0.0, 3.0,  # 0-5 rising to TS
        -0.4, -0.3, -0.2, -0.1, 0.0,        # 6-10 plateau (stretched side)
    ]
    profile, frames = _orca_profile(
        tmp_path,
        relative_kcal=relative_kcal,
        off_path=(),
    )
    policy = policy_from_config(
        {
            "selection": {
                "int_seed": {"mode": "ts_to_effective_endpoint_midpoint"},
                "int_plateau": {
                    "min_consecutive_frames": 3,
                    "energy_window_kcal_mol": 2.0,
                    "min_ts_separation_A": 0.10,
                    "max_slope_kcal_mol_A": 40.0,
                },
            }
        }
    )
    selection = select_path_seeds(profile, policy)

    assert selection is not None
    assert selection.has_independent_int is True
    int_seed = selection.int_search_seed
    assert int_seed is not None
    assert int_seed["selection_mode"] == "stretch_plateau"
    assert int_seed["frame_index"] in (6, 7, 8, 9, 10)


def test_plateau_detector_returns_none_without_plateau(tmp_path: Path) -> None:
    """No plateau -> fall back to midpoint (has_independent_int stays False)."""
    from rph_core.steps.step2_retro.path_selector import (
        policy_from_config,
        select_path_seeds,
    )

    relative_kcal = [-8.0, -5.0, -2.0, 1.0, 4.0, 3.0, 2.0, 1.5, 1.0, 0.5]
    profile, frames = _orca_profile(tmp_path, relative_kcal=relative_kcal)
    policy = policy_from_config(
        {"selection": {"int_seed": {"mode": "ts_to_effective_endpoint_midpoint"}}}
    )
    selection = select_path_seeds(profile, policy)

    assert selection is not None
    assert selection.has_independent_int is False
    int_seed = selection.int_search_seed
    assert int_seed is None or int_seed["selection_mode"] != "stretch_plateau"


def test_compare_graph_topology_grace_for_forming_bond_contacts() -> None:
    """New edges touching a forming-bond atom are expected for open intermediates."""
    import numpy as np

    from rph_core.steps.step2_retro.geometry_guard import compare_graph_topology

    # Product: two atoms forming the bond at 1.5 A.
    product_coords = np.asarray(
        [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 2.0, 0.0]]
    )
    # Candidate: forming bond open to 3.0 A; atom 0 now contacts atom 2 (new
    # edge touching a forming-bond atom) while the forming bond is lost.
    candidate_coords = np.asarray(
        [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.1, 0.0, 0.0]]
    )
    symbols = ["H", "H", "H"]

    result = compare_graph_topology(
        product_coords,
        candidate_coords,
        symbols,
        forming_bonds=[(0, 1)],
    )

    # The new 0-2 contact touches forming-bond atom 0 -> excluded; no lost
    # non-forming edges; valid.
    assert result.is_valid is True
    assert result.new_edges == []


def test_compare_graph_topology_grace_edges_tolerance() -> None:
    """Up to topology_grace_edges non-forming new edges are tolerated."""
    import numpy as np

    from rph_core.steps.step2_retro.geometry_guard import compare_graph_topology

    product_coords = np.asarray(
        [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 5.5, 0.0]]
    )
    # Candidate: forming bond open; atoms 2 and 3 approach to 1.4 A (new
    # C-C edge NOT touching a forming-bond atom) -> needs grace to pass.
    candidate_coords = np.asarray(
        [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 4.4, 0.0]]
    )
    symbols = ["H", "H", "C", "C"]

    strict = compare_graph_topology(
        product_coords, candidate_coords, symbols, forming_bonds=[(0, 1)]
    )
    graceful = compare_graph_topology(
        product_coords,
        candidate_coords,
        symbols,
        forming_bonds=[(0, 1)],
        topology_grace_edges=1,
    )

    assert strict.is_valid is False
    assert graceful.is_valid is True
