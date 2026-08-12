"""Single-artifact S2 figure for the canonical ORCA scan workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import numpy as np

HARTREE_TO_KCAL = 627.509
MM = 1 / 25.4
COLORS = {
    "b97": "#1A1A1A",
    "gfn2": "#6F8FAF",
    "ts": "#D55E00",
    "int": "#7A5195",
    "endpoint": "#E69F00",
    "knee": "#009E73",
}


def render_s2_profile_figure(
    scan_profile: Mapping[str, Any], output_path: Path
) -> Path:
    """Render the exploratory and seed-selection PATHs in one PNG."""
    import matplotlib.pyplot as plt
    from matplotlib import rc_context

    branches = _plot_branches(scan_profile)
    panel_heights = [0.34 if branch is None else 1.0 for _, branch, _, _ in branches]
    with rc_context({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 8, "legend.fontsize": 7,
        "axes.linewidth": 0.7, "savefig.transparent": False,
    }):
        figure = plt.figure(
            figsize=(180 * MM, (62 * (0.34 + sum(panel_heights))) * MM),
        )
        grid = figure.add_gridspec(
            nrows=len(branches) + 1,
            ncols=1,
            height_ratios=[0.34] + panel_heights,
            hspace=0.52,
        )
        header = figure.add_subplot(grid[0, 0])
        _render_header(header, scan_profile)
        for index, (title, branch, exploratory, unavailable_reason) in enumerate(branches):
            axis = figure.add_subplot(grid[index + 1, 0])
            if branch is None:
                _render_unavailable_panel(axis, title, unavailable_reason)
            else:
                _render_panel(axis, branch, title=title, exploratory=exploratory)
        figure.subplots_adjust(left=0.12, right=0.985, bottom=0.065, top=0.975)
        output_path = Path(output_path).with_suffix(".png")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path, dpi=360, pad_inches=0.02)
        plt.close(figure)
    return output_path


def _plot_branches(
    profile: Mapping[str, Any],
) -> list[tuple[str, Optional[Mapping[str, Any]], bool, Optional[str]]]:
    branches = dict(profile.get("path_branches") or {})
    if not branches:
        title = (
            "B97-3c relaxed scan"
            if profile.get("profile_kind") == "relaxed_scan_rescue"
            else "ORCA GFN2-xTB relaxed scan + B97-3c SP"
        )
        rendered = [(title, profile, False, None)]
        rescue_profile = profile.get("_rescue_profile_for_plot")
        if isinstance(rescue_profile, Mapping):
            rendered.append(("B97-3c relaxed-scan rescue", dict(rescue_profile), False, None))
        return rendered
    ordered = ["full_endpoint", "valid_corridor"]
    ordered.extend(name for name in branches if name not in ordered)
    selection_policy = dict(profile.get("selection_policy") or {})
    selection_decision = dict(profile.get("selection_decision") or {})
    selections = dict(profile.get("selections") or {})
    scan_quality = dict(profile.get("scan_quality") or {})
    anchors = dict(profile.get("anchors") or {})
    root_quality = dict(profile.get("trajectory_quality") or {})
    ts_search_seed = dict(profile.get("ts_search_seed") or {})
    int_search_seed = dict(profile.get("int_search_seed") or {})
    xtb_paths = dict(profile.get("xtb_paths") or {})
    rendered: list[tuple[str, Optional[Mapping[str, Any]], bool, Optional[str]]] = []
    has_new_rule_selection = isinstance(profile.get("selection_replay"), Mapping)
    for name in ordered:
        if name not in branches:
            title = (
                "Full-endpoint exploratory PATH"
                if name == "full_endpoint"
                else "Last-valid-to-product control PATH"
            )
            reason = (
                "Not independently calculated: full endpoint and valid-corridor "
                "start from the same topology-valid frame."
                if name == "full_endpoint" and xtb_paths.get("valid_corridor")
                else "PATH branch unavailable; no curve was used for selection."
            )
            rendered.append((title, None, name == "full_endpoint", reason))
            continue
        branch = dict(branches[name] or {})
        quality = dict(branch.get("trajectory_quality") or {})
        selected = bool(branch.get("selected_for_s3", False))
        branch_profile: dict[str, Any] = {
            "branch_name": name,
            "xtb_path": {"path_arclength": branch.get("path_arclength") or []},
            "reaction_coordinate_angstrom": branch.get("reaction_coordinate_angstrom") or [],
            "energy_curves": branch.get("energy_curves") or {},
            "selection_policy": selection_policy,
            "selection_decision": selection_decision,
            "scan_quality": scan_quality,
            "scan_parameters": dict(profile.get("scan_parameters") or {}),
            "anchors": anchors,
            "seed_evidence": None if has_new_rule_selection else profile.get("seed_evidence"),
            "s2_state": profile.get("s2_state"),
            "selection_source": profile.get("selection_source"),
            "has_independent_int": profile.get("has_independent_int"),
            "rejection_reason": profile.get("rejection_reason"),
            "endpoint_direction": profile.get("endpoint_direction"),
            "endpoint_exclusion_frames": profile.get("endpoint_exclusion_frames"),
            "trajectory_quality": {
                "off_path_indices": quality.get("selection_excluded_indices")
                or quality.get("off_path_indices")
                or [],
                "endpoint_index": root_quality.get("endpoint_index"),
                "endpoint_excluded": root_quality.get("endpoint_excluded"),
            },
            "selections": {} if has_new_rule_selection else selections if selected else {},
            "ts_search_seed": {} if has_new_rule_selection else ts_search_seed if selected else {},
            "int_search_seed": {} if has_new_rule_selection else int_search_seed if selected else {},
        }
        if name == "full_endpoint":
            title = "Full-endpoint exploratory PATH"
        else:
            title = "Last-valid-to-product control PATH"
        rendered.append((title, branch_profile, name == "full_endpoint", None))
    rescue_profile = profile.get("_rescue_profile_for_plot")
    if isinstance(rescue_profile, Mapping):
        rendered.append(("B97-3c relaxed-scan rescue", dict(rescue_profile), False, None))
    return rendered


def _render_header(ax: Any, profile: Mapping[str, Any]) -> None:
    replay = profile.get("selection_replay")
    if isinstance(replay, Mapping):
        _render_new_rule_header(ax, replay)
        return
    decision = dict(profile.get("selection_decision") or {})
    ts = _selection_payload(profile, "ts_guess", "ts_search_seed")
    intermediate = _selection_payload(profile, "intermediate", "int_search_seed")
    selected = str(decision.get("selected_branch") or "PATH")
    rule = _selection_rule_label(str(decision.get("rule") or "path selection"))
    int_mode = _intermediate_mode(profile, intermediate)
    int_label = {
        "stable_basin_candidate": "resolved pre-TS basin",
        "late_pre_ts_platform_fallback": "late pre-TS platform seed",
        "shared_ts_fallback": "shared with TS (no clean independent pre-TS seed)",
        "search_seed": "search seed",
        "ts_to_effective_endpoint_midpoint": "TS-to-effective-endpoint midpoint seed",
    }.get(int_mode, int_mode.replace("_", " "))
    ts_index = _display_index(ts)
    int_index = _display_index(intermediate)
    rescue_profile = profile.get("profile_kind") == "relaxed_scan_rescue"
    state_line = _selection_state_line(profile)
    ax.axis("off")
    if rescue_profile:
        rescue = dict(profile.get("rescue") or {})
        admission = dict(rescue.get("candidate_admission") or {})
        mode = "ScanTS" if rescue.get("scan_ts") else "simultaneous relaxed scan"
        title_y, mode_y, ts_y, state_y = (0.90, 0.62, 0.34, 0.06) if state_line else (0.86, 0.51, 0.16, 0.0)
        ax.text(0.0, title_y, "S2 B97-3c relaxed-scan rescue", fontsize=9.0, fontweight="bold", va="top")
        ax.text(
            0.0, mode_y,
            f"Mode: {mode}    Rule: B97-3c relaxed-scan knee + right shift    Energy reference: stretched endpoint",
            fontsize=7.5, va="top",
        )
        knee = dict(profile.get("knee_evidence") or {})
        endpoint = dict(profile.get("endpoint_evidence") or {})
        knee_index = _coerce_int(knee.get("frame_index"))
        valid_endpoint = _coerce_int(endpoint.get("effective_endpoint_index"))
        ax.text(
            0.0, ts_y,
            f"Knee frame: {knee_index if knee_index is not None else '-'}    TS frame: {ts_index}    INT frame: {int_index}    "
            f"Effective endpoint: {valid_endpoint if valid_endpoint is not None else '-'}    "
            f"Barrier diagnostic: {_format_energy(admission.get('barrier_from_reactant_kcal_mol'))} kcal mol$^{{-1}}$",
            fontsize=7.5, va="top",
        )
        if state_line:
            ax.text(0.0, state_y, state_line, fontsize=7.3, va="top")
        return

    title_y, path_y, seed_y, state_y = (0.90, 0.62, 0.34, 0.06) if state_line else (0.86, 0.51, 0.16, 0.0)
    ax.text(0.0, title_y, "S2 ORCA relaxed-scan seed selection", fontsize=9.0, fontweight="bold", va="top")
    ax.text(0.0, path_y, f"Selected scan: {selected}    Rule: {rule}", fontsize=7.5, va="top")
    ax.text(0.0, seed_y, f"TS seed frame: {ts_index}    INT: {int_label} (frame {int_index})", fontsize=7.5, va="top")
    if state_line:
        ax.text(0.0, state_y, state_line, fontsize=7.3, va="top")


def _render_new_rule_header(ax: Any, selection: Mapping[str, Any]) -> None:
    """Render only the possible seeds from the current unified selector."""
    ts = dict(selection.get("ts_search_seed") or {})
    intermediate = dict(selection.get("int_search_seed") or {})
    knee = dict(selection.get("knee_evidence") or {})
    endpoint = dict(selection.get("endpoint_evidence") or {})
    dispatch = dict(selection.get("s3_dispatch") or {})
    ts_index = _display_index(ts)
    int_index = _display_index(intermediate)
    knee_index = _coerce_int(knee.get("frame_index"))
    endpoint_index = _coerce_int(endpoint.get("effective_endpoint_index"))
    shift = _format_energy(ts.get("right_shift_A"))
    source = str(selection.get("selection_source") or "B97-3c relaxed-scan rescue")
    if source == "orca_relaxed_scan":
        source = "B97-3c relaxed-scan rescue"
    state = str(selection.get("s2_state") or "unresolved")
    assessment = "possible seeds found" if state in {"path_seeded", "gfn2_seeded", "rescue_seeded"} else "no possible seed"
    evidence = str(selection.get("seed_evidence") or "none").replace("_", " ")
    ax.axis("off")
    ax.text(0.0, 0.90, "S2 unified seed selection", fontsize=9.0, fontweight="bold", va="top")
    ax.text(
        0.0, 0.62,
        f"Source: {source}    New-rule assessment: {assessment}    Evidence: {evidence}",
        fontsize=7.5, va="top",
    )
    ax.text(
        0.0, 0.36,
        f"Knee: frame {knee_index if knee_index is not None else '-'}    "
        f"TS seed: frame {ts_index} (knee + {shift} Å)    "
        f"INT seed: frame {int_index} (TS-to-effective-endpoint midpoint)",
        fontsize=7.5, va="top",
    )
    ax.text(
        0.0, 0.10,
        f"Effective endpoint: frame {endpoint_index if endpoint_index is not None else '-'}    "
        f"Possible S3 dispatch: TS={'yes' if dispatch.get('submit_ts') else 'no'}, "
        f"INT={'yes' if dispatch.get('submit_intermediate') else 'no'}",
        fontsize=7.3, va="top",
    )
def _render_unavailable_panel(ax: Any, title: str, reason: Optional[str]) -> None:
    ax.set_title(title, loc="left", fontsize=7.4, pad=4)
    ax.text(0.5, 0.5, reason or "PATH branch unavailable.", ha="center", va="center", fontsize=7.5,
            transform=ax.transAxes, wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)


def _render_panel(
    ax: Any,
    profile: Mapping[str, Any],
    *,
    title: str,
    exploratory: bool,
) -> None:
    x, xlabel = _x_axis(profile)
    gfn2 = _relative_curve(profile, "gfn2")
    b97 = _relative_curve(profile, "b973c")
    count = max(len(gfn2), len(b97))
    if len(x) != count:
        raise ValueError("S2 profile coordinate and energy curves are misaligned")
    excluded = {
        int(index)
        for index in (profile.get("trajectory_quality") or {}).get("off_path_indices", [])
        if 0 <= int(index) < count
    }
    mask = set() if exploratory else excluded
    if len(gfn2) == count and any(value is not None for value in gfn2):
        ax.plot(x, _masked(gfn2, mask), color=COLORS["gfn2"], linestyle=(0, (4, 2)),
                linewidth=1.0, marker="o", markersize=2.6, markevery=2, label="ORCA GFN2-xTB")
    if len(b97) == count and any(value is not None for value in b97):
        ax.plot(x, _masked(b97, mask), color=COLORS["b97"], linewidth=1.2,
                marker="o", markersize=2.7, label="B97-3c")
    endpoint_range = _endpoint_exclusion_range(profile, count)
    if endpoint_range is not None:
        left, right = _index_bounds(x, *endpoint_range)
        ax.axvspan(
            left,
            right,
            color=COLORS["endpoint"],
            alpha=0.12,
            linewidth=0,
            label="endpoint-exclusion buffer",
            zorder=0,
        )
    if excluded:
        for start, end in _contiguous_ranges(sorted(excluded)):
            left, right = _index_bounds(x, start, end)
            ax.axvspan(
                left, right, color="#B53B3B", alpha=0.10, linewidth=0,
                label="topology-distorted / excluded" if start == min(excluded) else None,
                zorder=0,
            )
        boundaries = sorted({min(excluded), max(excluded)})
        ys = [_selected_energy(profile, index) for index in boundaries]
        points = [(x[index], y) for index, y in zip(boundaries, ys) if y is not None]
        if points:
            ax.scatter(*zip(*points), marker="x", color="#8B1E1E", s=22,
                       linewidths=0.9, zorder=6)
    endpoint = dict(profile.get("endpoint_evidence") or {})
    valid_endpoint_index = _index(
        {"index": endpoint.get("effective_endpoint_index")}, count
    )
    if valid_endpoint_index is not None:
        ax.axvline(
            x[valid_endpoint_index],
            color=COLORS["endpoint"],
            linestyle=":",
            linewidth=0.9,
            label="effective endpoint",
            zorder=2,
        )
    knee = dict(profile.get("knee_evidence") or {})
    knee_index = _index({"index": knee.get("frame_index")}, count)
    if knee_index is not None:
        _seed(
            ax,
            x[knee_index],
            _selected_energy(profile, knee_index),
            "o",
            COLORS["knee"],
            "energy knee",
            filled=False,
            annotation="knee",
        )
        ax.axvline(x[knee_index], color=COLORS["knee"], linestyle="-.", linewidth=0.8, zorder=1)
    ts = _selection_payload(profile, "ts_guess", "ts_search_seed")
    intermediate = _selection_payload(profile, "intermediate", "int_search_seed")
    ts_index = _index(ts, count)
    if ts_index is not None:
        ts_annotation = _seed_annotation(profile, ts, kind="ts")
        _seed(ax, x[ts_index], _selected_energy(profile, ts_index), "^", COLORS["ts"],
              "TS seed", filled=True, annotation=ts_annotation)
        ax.axvline(x[ts_index], color=COLORS["ts"], linestyle="--", linewidth=0.8, zorder=1)
    int_index = _index(intermediate, count)
    shared_with_ts = _shared_ts_seed(intermediate)
    if int_index is not None and not shared_with_ts:
        int_mode = _intermediate_mode(profile, intermediate)
        platform = int_mode == "late_pre_ts_platform_fallback"
        fallback = int_mode == "midpoint_fallback"
        endpoint_midpoint = int_mode == "ts_to_effective_endpoint_midpoint"
        int_annotation = _seed_annotation(profile, intermediate, kind="int")
        _seed(ax, x[int_index], _selected_energy(profile, int_index), "D", COLORS["int"],
              "INT platform seed"
              if platform
              else "INT endpoint midpoint seed"
              if endpoint_midpoint
              else "INT search seed"
              if fallback
              else "INT basin candidate",
              filled=not (fallback or platform or endpoint_midpoint), annotation=int_annotation)
    _style(ax)
    ax.set_title(title, loc="left", fontsize=7.4, pad=4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Relative energy, ΔE (kcal mol$^{-1}$)")
    ax.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.88,
              edgecolor="#D0D0D0", borderpad=0.45)


def _x_axis(profile: Mapping[str, Any]) -> tuple[list[float], str]:
    count = max(
        len(_relative_curve(profile, "gfn2")),
        len(_relative_curve(profile, "b973c")),
    )
    coordinate = list(profile.get("reaction_coordinate_angstrom") or [])
    if len(coordinate) == count and all(value is not None for value in coordinate):
        return [float(value) for value in coordinate], "Mean forming-bond distance (Å)"
    return list(np.arange(count, dtype=float)), "Scan frame"


def _relative_curve(profile: Mapping[str, Any], name: str) -> list[Optional[float]]:
    curve = dict((profile.get("energy_curves") or {}).get(name) or {})
    relative = list(
        curve.get("relative_energies_kcal_mol")
        or (profile.get("relative_energies_kcal_mol") if name in {"gfn2", "xtb"} else [])
        or []
    )
    if relative:
        return [None if value is None else float(value) for value in relative]
    energies = list(
        curve.get("energies_hartree")
        or (profile.get("energies_hartree") if name in {"gfn2", "xtb"} else [])
        or []
    )
    reference = next((float(value) for value in energies if value is not None), None)
    return [None if value is None or reference is None else (float(value) - reference) * HARTREE_TO_KCAL for value in energies]


def _masked(values: Sequence[Optional[float]], indices: set[int]) -> np.ndarray:
    result = np.asarray([np.nan if value is None else float(value) for value in values], dtype=float)
    for index in indices:
        result[index] = np.nan
    return result


def _selection_payload(
    profile: Mapping[str, Any],
    selection_key: str,
    fallback_key: str,
) -> dict[str, Any]:
    fallback = dict(profile.get(fallback_key) or {})
    selections = dict(profile.get("selections") or {})
    selection = dict(selections.get(selection_key) or {})
    return {**fallback, **selection}


def _selection_index(selection: Mapping[str, Any]) -> Optional[int]:
    for key in ("index", "frame_index"):
        value = selection.get(key)
        index = _coerce_int(value)
        if index is not None:
            return index
    return None


def _index(selection: Mapping[str, Any], count: int) -> Optional[int]:
    index = _selection_index(selection)
    if index is None:
        return None
    return index if 0 <= index < count else None


def _selected_energy(profile: Mapping[str, Any], index: int) -> Optional[float]:
    selection = dict(profile.get("selection_policy") or {})
    preferred = "b973c" if str(selection.get("actual_source", "")).lower().startswith("b97") else "gfn2"
    preferred_curve = _relative_curve(profile, preferred)
    if index < len(preferred_curve) and preferred_curve[index] is not None:
        return preferred_curve[index]
    gfn2_curve = _relative_curve(profile, "gfn2")
    return gfn2_curve[index] if index < len(gfn2_curve) else None


def _format_energy(value: Any) -> str:
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return "n/a"


def _seed(ax: Any, x: float, y: Optional[float], marker: str, color: str, label: str,
          *, filled: bool, annotation: Optional[str] = None) -> None:
    if y is None:
        return
    ax.scatter([x], [y], marker=marker, s=35, facecolor=color if filled else "white",
               edgecolor="white" if filled else color, linewidth=0.9, label=label, zorder=8)
    if annotation:
        y_offset = 8 if marker == "^" else -10
        ax.annotate(
            annotation,
            xy=(x, y),
            xytext=(4, y_offset),
            textcoords="offset points",
            fontsize=6.0,
            color=color,
            va="bottom" if y_offset >= 0 else "top",
            bbox={
                "boxstyle": "round,pad=0.15",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.78,
            },
            zorder=9,
        )


def _contiguous_ranges(indices: Sequence[int]) -> list[tuple[int, int]]:
    if not indices:
        return []
    ranges: list[tuple[int, int]] = []
    start = previous = int(indices[0])
    for value in indices[1:]:
        value = int(value)
        if value != previous + 1:
            ranges.append((start, previous))
            start = value
        previous = value
    ranges.append((start, previous))
    return ranges


def _display_index(selection: Mapping[str, Any]) -> str:
    index = _selection_index(selection)
    return str(index) if index is not None else "-"


def _selection_state_line(profile: Mapping[str, Any]) -> Optional[str]:
    s2_state = profile.get("s2_state")
    seed_evidence = profile.get("seed_evidence")
    parts = []
    if s2_state:
        parts.append(f"s2_state={s2_state}")
    if seed_evidence:
        parts.append(f"seed_evidence={seed_evidence}")
    return "    ".join(parts) if parts else None


def _seed_annotation(
    profile: Mapping[str, Any],
    selection: Mapping[str, Any],
    *,
    kind: str,
) -> Optional[str]:
    parts = []
    index = _selection_index(selection)
    if index is not None:
        parts.append(f"frame={index}")
    seed_evidence = profile.get("seed_evidence")
    if seed_evidence:
        parts.append(f"evidence={str(seed_evidence).replace('_', ' ')}")
    confidence = _seed_confidence(profile, selection, kind)
    if confidence:
        parts.append(f"conf={confidence}")
    return "\n".join(parts) if parts else None


def _seed_confidence(
    profile: Mapping[str, Any],
    selection: Mapping[str, Any],
    kind: str,
) -> Optional[str]:
    scan_quality = dict(profile.get("scan_quality") or {})
    fallback_key = "ts_search_seed" if kind == "ts" else "int_search_seed"
    fallback = dict(profile.get(fallback_key) or {})
    decision = dict(profile.get("selection_decision") or {})
    candidates = [selection.get("confidence"), fallback.get("confidence")]
    if kind == "ts":
        candidates.extend([
            decision.get("ts_confidence"),
            scan_quality.get("ts_guess_confidence"),
        ])
    else:
        candidates.append(scan_quality.get("intermediate_confidence"))
    for candidate in candidates:
        if candidate in (None, "", "unavailable"):
            continue
        return str(candidate)
    return None


def _shared_ts_seed(selection: Mapping[str, Any]) -> bool:
    if bool(selection.get("shared_with_ts", False)):
        return True
    return str(selection.get("selection_mode") or "") == "shared_ts_fallback"


def _intermediate_mode(profile: Mapping[str, Any], intermediate: Mapping[str, Any]) -> str:
    mode = str(intermediate.get("selection_mode") or "")
    if mode:
        return mode
    if _shared_ts_seed(intermediate):
        return "shared_ts_fallback"
    if _selection_index(intermediate) is None:
        return "not selected"
    if bool(profile.get("has_independent_int", False)):
        return "stable_basin_candidate"
    if str(profile.get("seed_evidence") or "") == "monotonic_shoulder":
        return "late_pre_ts_platform_fallback"
    return "search_seed"


def _endpoint_exclusion_range(
    profile: Mapping[str, Any],
    count: int,
) -> Optional[tuple[int, int]]:
    if count <= 0:
        return None
    frames = _endpoint_exclusion_frames(profile)
    endpoint_index = _endpoint_index(profile, count)
    endpoint_excluded = bool((profile.get("trajectory_quality") or {}).get("endpoint_excluded"))
    if frames == 0:
        return None
    if frames is None:
        if endpoint_index is None or not endpoint_excluded:
            return None
        frames = 1
    direction = _endpoint_direction(profile, count)
    if direction == "start":
        end = min(count - 1, max(0, frames - 1))
        return 0, end
    if direction == "end":
        start = max(0, count - frames)
        return start, count - 1
    if endpoint_index is None:
        return None
    return endpoint_index, endpoint_index


def _endpoint_exclusion_frames(profile: Mapping[str, Any]) -> Optional[int]:
    scan_parameters = dict(profile.get("scan_parameters") or {})
    selection_cfg = dict(scan_parameters.get("selection") or {})
    for candidate in (
        profile.get("endpoint_exclusion_frames"),
        (profile.get("selection_decision") or {}).get("endpoint_exclusion_frames"),
        (profile.get("selection_policy") or {}).get("endpoint_exclusion_frames"),
        selection_cfg.get("endpoint_exclusion_frames"),
    ):
        value = _coerce_int(candidate)
        if value is not None:
            return max(0, value)
    return None


def _endpoint_direction(profile: Mapping[str, Any], count: int) -> Optional[str]:
    for candidate in (
        profile.get("endpoint_direction"),
        (profile.get("selection_decision") or {}).get("endpoint_direction"),
        (profile.get("selection_policy") or {}).get("endpoint_direction"),
        (profile.get("trajectory_quality") or {}).get("endpoint_direction"),
    ):
        if not isinstance(candidate, str):
            continue
        normalized = candidate.strip().lower()
        if normalized in {"start", "first", "left", "begin"}:
            return "start"
        if normalized in {"end", "last", "right", "endpoint", "scan_endpoint"}:
            return "end"
    endpoint_index = _endpoint_index(profile, count)
    if endpoint_index == 0:
        return "start"
    if endpoint_index == count - 1:
        return "end"
    return None


def _endpoint_index(profile: Mapping[str, Any], count: int) -> Optional[int]:
    anchors = dict(profile.get("anchors") or {})
    anchor = anchors.get("scan_endpoint") or {}
    if not hasattr(anchor, "get"):
        anchor = {}
    for candidate in (
        anchor.get("index"),
        anchors.get("scan_endpoint_index"),
        (profile.get("trajectory_quality") or {}).get("endpoint_index"),
    ):
        index = _coerce_int(candidate)
        if index is None:
            continue
        if 0 <= index < count:
            return index
    return None


def _index_bounds(x: Sequence[float], start: int, end: int) -> tuple[float, float]:
    left = float(x[start])
    right = float(x[end])
    if start > 0:
        left = 0.5 * (float(x[start - 1]) + float(x[start]))
    if end + 1 < len(x):
        right = 0.5 * (float(x[end]) + float(x[end + 1]))
    return left, right


def _coerce_int(value: Any) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _selection_rule_label(rule: str) -> str:
    return {
        "precise_significant_peak_preferred": "valid-corridor B97-3c peak preferred",
        "valid_corridor_credible_weak_peak_preferred": "clean valid-corridor weak peak preferred",
        "full_endpoint_significant_peak_fallback": "full-endpoint B97-3c peak fallback",
        "full_endpoint_credible_weak_peak_fallback": "clean full-endpoint weak peak fallback",
        "full_endpoint_boundary_seed_fallback": "topology-boundary TS seed fallback",
        "full_endpoint_low_confidence_fallback": "low-confidence full-endpoint fallback",
        "valid_corridor_fallback_no_significant_peak": "valid-corridor fallback",
    }.get(rule, rule.replace("_", " "))


def _style(ax: Any) -> None:
    from matplotlib.ticker import MaxNLocator

    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.4, alpha=0.55)
    ax.grid(axis="x", visible=False)
    ax.tick_params(direction="out", width=0.7, length=3)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=4))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=4))
