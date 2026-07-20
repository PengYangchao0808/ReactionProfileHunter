"""Concise, single-panel S2 refined-PATH profile figure."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import numpy as np

HARTREE_TO_KCAL = 627.509
MM = 1 / 25.4
COLORS = {"b97": "#1A1A1A", "xtb": "#9A9A9A", "ts": "#D55E00", "int": "#7A5195"}


def render_s2_profile_figure(
    scan_profile: Mapping[str, Any], output_path: Path
) -> Path:
    """Render B97-3c/xTB energy curves plus TS and optional INT search seed."""
    import matplotlib.pyplot as plt
    from matplotlib import rc_context

    x, xlabel = _x_axis(scan_profile)
    xtb = _relative_curve(scan_profile, "xtb")
    b97 = _relative_curve(scan_profile, "b973c")
    count = len(xtb)
    if len(x) != count:
        raise ValueError("S2 profile coordinate and energy curves are misaligned")
    off_path = {
        int(index)
        for index in (scan_profile.get("trajectory_quality") or {}).get("off_path_indices", [])
        if 0 <= int(index) < count
    }
    selections = dict(scan_profile.get("selections") or {})
    ts = dict(selections.get("ts_guess") or {})
    intermediate = dict(selections.get("intermediate") or {})

    with rc_context({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 8, "legend.fontsize": 7,
        "axes.linewidth": 0.7, "savefig.transparent": False,
    }):
        fig, ax = plt.subplots(figsize=(180 * MM, 88 * MM))
        ax.plot(x, _masked(xtb, off_path), color=COLORS["xtb"], linestyle=(0, (4, 2)),
                linewidth=1.0, marker="o", markersize=2.6, markevery=2, label="GFN2-xTB")
        if any(value is not None for value in b97):
            ax.plot(x, _masked(b97, off_path), color=COLORS["b97"], linewidth=1.2,
                    marker="o", markersize=2.7, label="B97-3c")
        if off_path:
            ys = [_selected_energy(scan_profile, index) for index in sorted(off_path)]
            points = [(x[index], y) for index, y in zip(sorted(off_path), ys) if y is not None]
            if points:
                ax.scatter(*zip(*points), marker="x", color="#8B1E1E", s=24,
                           linewidths=0.9, label="off-path / excluded", zorder=6)
        ts_index = _index(ts, count)
        if ts_index is not None:
            _seed(ax, x[ts_index], _selected_energy(scan_profile, ts_index), "^", COLORS["ts"],
                  "TS seed", f"TS seed · frame {ts_index}", (-48, 16), filled=True)
            ax.axvline(x[ts_index], color=COLORS["ts"], linestyle="--", linewidth=0.8, zorder=1)
        int_index = _index(intermediate, count)
        if int_index is not None:
            fallback = intermediate.get("selection_mode") == "midpoint_fallback"
            _seed(ax, x[int_index], _selected_energy(scan_profile, int_index), "D", COLORS["int"],
                  "INT search seed" if fallback else "INT basin candidate",
                  f"{'INT search seed' if fallback else 'INT basin candidate'} · frame {int_index}",
                  (14, -22), filled=not fallback)
        _style(ax)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Relative energy, ΔE (kcal mol$^{-1}$)")
        ax.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.88,
                  edgecolor="#D0D0D0", borderpad=0.45)
        output_path = Path(output_path).with_suffix(".png")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=360, bbox_inches="tight", pad_inches=0.02)
        plt.close(fig)
    return output_path


def _x_axis(profile: Mapping[str, Any]) -> tuple[list[float], str]:
    count = len(_relative_curve(profile, "xtb"))
    arc = list(((profile.get("xtb_path") or {}).get("path_arclength") or []))
    if len(arc) == count and all(value is not None for value in arc):
        return [float(value) for value in arc], "Reaction-path coordinate, s (Å)"
    coordinate = list(profile.get("reaction_coordinate_angstrom") or [])
    if len(coordinate) == count and all(value is not None for value in coordinate):
        return [float(value) for value in coordinate], "Mean forming-bond distance (Å)"
    return list(np.arange(count, dtype=float)), "PATH frame"


def _relative_curve(profile: Mapping[str, Any], name: str) -> list[Optional[float]]:
    curve = dict((profile.get("energy_curves") or {}).get(name) or {})
    relative = list(
        curve.get("relative_energies_kcal_mol")
        or (profile.get("relative_energies_kcal_mol") if name == "xtb" else [])
        or []
    )
    if relative:
        return [None if value is None else float(value) for value in relative]
    energies = list(
        curve.get("energies_hartree")
        or (profile.get("energies_hartree") if name == "xtb" else [])
        or []
    )
    reference = next((float(value) for value in energies if value is not None), None)
    return [None if value is None or reference is None else (float(value) - reference) * HARTREE_TO_KCAL for value in energies]


def _masked(values: Sequence[Optional[float]], indices: set[int]) -> np.ndarray:
    result = np.asarray([np.nan if value is None else float(value) for value in values], dtype=float)
    for index in indices:
        result[index] = np.nan
    return result


def _index(selection: Mapping[str, Any], count: int) -> Optional[int]:
    try:
        index = int(selection.get("index"))
    except (TypeError, ValueError):
        return None
    return index if 0 <= index < count else None


def _selected_energy(profile: Mapping[str, Any], index: int) -> Optional[float]:
    selection = dict(profile.get("selection_policy") or {})
    preferred = "b973c" if str(selection.get("actual_source", "")).lower().startswith("b97") else "xtb"
    value = _relative_curve(profile, preferred)[index]
    return value if value is not None else _relative_curve(profile, "xtb")[index]


def _seed(ax: Any, x: float, y: Optional[float], marker: str, color: str, label: str,
          text: str, offset: tuple[int, int], *, filled: bool) -> None:
    if y is None:
        return
    ax.scatter([x], [y], marker=marker, s=35, facecolor=color if filled else "white",
               edgecolor="white" if filled else color, linewidth=0.9, label=label, zorder=8)
    ax.annotate(text, xy=(x, y), xytext=offset, textcoords="offset points", fontsize=7,
                color=color, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 1.2},
                arrowprops={"arrowstyle": "-", "color": color, "lw": 0.8})


def _style(ax: Any) -> None:
    from matplotlib.ticker import MaxNLocator

    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.4, alpha=0.55)
    ax.grid(axis="x", visible=False)
    ax.tick_params(direction="out", width=0.7, length=3)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=4))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=4))
