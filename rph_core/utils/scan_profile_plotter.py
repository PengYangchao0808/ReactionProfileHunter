"""
Scan Profile Plotter - Visualization tool for XTB scan energy profiles.

This module provides utilities to plot scan energy profiles with bond distance
on the x-axis and energy on the y-axis.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

HARTREE_TO_KCAL = 627.509  # Hartree to kcal/mol conversion factor


def compute_scan_distances(
    start_distance: float,
    end_distance: float,
    num_steps: int,
    direction: str = "outward",
) -> List[float]:
    """
    Compute bond distances for each scan point.
    
    Args:
        start_distance: Starting bond distance in Angstrom
        end_distance: Ending bond distance in Angstrom  
        num_steps: Number of scan steps
        direction: Scan direction ("outward" or "inward")
    
    Returns:
        List of bond distances for each scan point
    """
    if direction == "outward":
        start = end_distance
        end = start_distance
    else:
        start = start_distance
        end = end_distance
    
    return np.linspace(start, end, num_steps).tolist()


def plot_scan_profile(
    scan_profile_json: Path,
    output_path: Optional[Path] = None,
    energy_unit: str = "kcal",
    show_peak: bool = True,
    figsize: Tuple[int, int] = (10, 6),
    ts_distance: Optional[float] = None,
    dipole_distance: Optional[float] = None,
    path_ts_distance: Optional[float] = None,
    path_ts_energy: Optional[float] = None,
    gau_xtb_distance: Optional[float] = None,
    gau_xtb_energy: Optional[float] = None,
    gau_xtb2_distance: Optional[float] = None,
    gau_xtb2_energy: Optional[float] = None,
) -> Optional[Path]:
    """
    Plot scan energy profile from scan_profile.json.
    
    Args:
        scan_profile_json: Path to scan_profile.json
        output_path: Optional output path for the plot. If None, saves next to json.
        energy_unit: Energy unit ("hartree" or "kcal")
        show_peak: Whether to highlight the maximum energy point
        figsize: Figure size (width, height) in inches
        ts_distance: TS guess distance from knee point algorithm (Angstrom)
        dipole_distance: Dipole intermediate distance from knee point algorithm (Angstrom)
        path_ts_distance: TS guess distance from xTB path search (Angstrom)
        path_ts_energy: TS guess energy from xTB path search (Hartree)
        gau_xtb_distance: TS guess distance from Gau_XTB optimization (Angstrom)
        gau_xtb_energy: TS guess energy from Gau_XTB optimization (Hartree)
    
    Returns:
        Path to saved plot, or None if plotting failed
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available, skipping plot generation")
        return None
    
    scan_profile_json = Path(scan_profile_json)
    if not scan_profile_json.exists():
        logger.error(f"Scan profile JSON not found: {scan_profile_json}")
        return None
    
    with open(scan_profile_json) as f:
        data = json.load(f)
    
    energies = data.get("energies_hartree")
    logger.debug(f"plot_scan_profile: found energies_hartree = {type(energies)}, len={len(energies) if energies else 0}")
    if not energies:
        logger.error(f"No energies found in scan profile. Available keys: {list(data.keys())}")
        return None
    
    params = data.get("scan_parameters", {})
    scan_start = params.get("scan_start_distance", 3.5)
    scan_end = params.get("scan_end_distance", 1.8)
    scan_steps = params.get("scan_steps", 20)
    direction = "outward" if scan_start > scan_end else "inward"
    
    distances = compute_scan_distances(scan_start, scan_end, scan_steps, direction)
    
    if len(distances) != len(energies):
        logger.warning(
            f"Distance points ({len(distances)}) != energy points ({len(energies)}), "
            f"adjusting..."
        )
        distances = np.linspace(distances[0], distances[-1], len(energies)).tolist()
    
    if energy_unit == "kcal":
        energies = [e * HARTREE_TO_KCAL for e in energies]
        ylabel = "Energy (kcal/mol)"
    else:
        ylabel = "Energy (Hartree)"
    
    fig, ax = plt.subplots(figsize=figsize)
    
    ax.plot(distances, energies, "b-o", markersize=4, linewidth=1.5, label="Scan Energy")
    
    markers_added = []
    
    if show_peak and energies:
        max_idx = np.argmax(energies)
        max_dist = distances[max_idx]
        max_energy = energies[max_idx]
        ax.scatter([max_dist], [max_energy], color="red", s=100, zorder=5, marker="^", label=f"Peak (d={max_dist:.2f}Å)")
        ax.axvline(x=max_dist, color="red", linestyle="--", alpha=0.5)
        markers_added.append("peak")
    
    if dipole_distance is not None and energies:
        dipole_idx = min(range(len(distances)), key=lambda i: abs(distances[i] - dipole_distance))
        dipole_dist = distances[dipole_idx]
        dipole_energy = energies[dipole_idx]
        ax.scatter([dipole_dist], [dipole_energy], color="green", s=120, zorder=5, marker="s", label=f"Dipole (d={dipole_dist:.2f}Å)")
        ax.axvline(x=dipole_dist, color="green", linestyle=":", alpha=0.7)
        markers_added.append("dipole")
    
    if ts_distance is not None and energies:
        ts_idx = min(range(len(distances)), key=lambda i: abs(distances[i] - ts_distance))
        ts_dist = distances[ts_idx]
        ts_energy = energies[ts_idx]
        ax.scatter([ts_dist], [ts_energy], color="orange", s=150, zorder=6, marker="*", label=f"TS Guess (d={ts_dist:.2f}Å)")
        ax.axvline(x=ts_dist, color="orange", linestyle="-.", alpha=0.7)
        markers_added.append("ts")
    
    if path_ts_distance is not None and energies:
        if path_ts_energy is not None:
            path_e = path_ts_energy * HARTREE_TO_KCAL if energy_unit == "kcal" else path_ts_energy
        else:
            path_idx = min(range(len(distances)), key=lambda i: abs(distances[i] - path_ts_distance))
            path_e = energies[path_idx]
        ax.scatter([path_ts_distance], [path_e], color="purple", s=200, zorder=7, marker="D", label=f"Path TS (d={path_ts_distance:.2f}Å, E={path_e:.2f})")
        ax.axvline(x=path_ts_distance, color="purple", linestyle="--", alpha=0.8)
    
    if gau_xtb_distance is not None and energies:
        if gau_xtb_energy is not None:
            gau_e = gau_xtb_energy * HARTREE_TO_KCAL if energy_unit == "kcal" else gau_xtb_energy
        else:
            gau_idx = min(range(len(distances)), key=lambda i: abs(distances[i] - gau_xtb_distance))
            gau_e = energies[gau_idx]
        ax.scatter([gau_xtb_distance], [gau_e], color="brown", s=250, zorder=8, marker="p", label=f"Gau_XTB S2.1 (d={gau_xtb_distance:.2f}Å, E={gau_e:.2f})")
        ax.axvline(x=gau_xtb_distance, color="brown", linestyle="-", linewidth=2, alpha=0.8)
    
    if gau_xtb2_distance is not None and energies:
        if gau_xtb2_energy is not None:
            gau2_e = gau_xtb2_energy * HARTREE_TO_KCAL if energy_unit == "kcal" else gau_xtb2_energy
        else:
            gau2_idx = min(range(len(distances)), key=lambda i: abs(distances[i] - gau_xtb2_distance))
            gau2_e = energies[gau2_idx]
        ax.scatter([gau_xtb2_distance], [gau2_e], color="red", s=300, zorder=9, marker="h", label=f"Gau_XTB S2.2 (d={gau_xtb2_distance:.2f}Å, E={gau2_e:.2f})")
        ax.axvline(x=gau_xtb2_distance, color="red", linestyle="-", linewidth=2, alpha=0.8)
    
    forming_bonds = data.get("forming_bonds", [])
    if forming_bonds:
        bond_str = ", ".join([f"{a}-{b}" for a, b in forming_bonds])
        title_suffix = f"\nForming bonds: {bond_str}"
    else:
        title_suffix = ""
    
    direction_label = "outward" if direction == "outward" else "inward"
    ax.set_xlabel(f"Bond Distance (Å) - {direction_label} scan")
    ax.set_ylabel(ylabel)
    ax.set_title(f"S2 Scan Energy Profile{title_suffix}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    ax.invert_xaxis() if direction == "outward" else None
    
    plt.tight_layout()
    
    if output_path is None:
        output_path = scan_profile_json.with_suffix(".png")
    
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    logger.info(f"Scan profile plot saved to: {output_path}")
    return output_path


def plot_path_search_profile(
    path_profile_json: Path,
    output_path: Optional[Path] = None,
    figsize: Tuple[int, int] = (10, 6),
) -> Optional[Path]:
    """
    Plot path search energy profile from path_profile.json.
    
    Args:
        path_profile_json: Path to path_profile.json
        output_path: Optional output path for the plot
        figsize: Figure size (width, height) in inches
    
    Returns:
        Path to saved plot, or None if plotting failed
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available, skipping plot generation")
        return None
    
    path_profile_json = Path(path_profile_json)
    if not path_profile_json.exists():
        logger.error(f"Path profile JSON not found: {path_profile_json}")
        return None
    
    with open(path_profile_json) as f:
        data = json.load(f)
    
    energies = data.get("energies", {})
    barrier_forward = energies.get("barrier_forward_kcal")
    barrier_backward = energies.get("barrier_backward_kcal")
    reaction_energy = energies.get("reaction_energy_kcal")
    
    if not any([barrier_forward, barrier_backward, reaction_energy]):
        logger.error("No energy data found in path profile")
        return None
    
    fig, ax = plt.subplots(figsize=figsize)
    
    labels = []
    values = []
    colors = []
    
    if barrier_forward is not None:
        labels.append("Forward Barrier (‡→P)")
        values.append(barrier_forward)
        colors.append("green")
    
    if reaction_energy is not None:
        labels.append("Reaction Energy (R→P)")
        values.append(reaction_energy)
        colors.append("blue")
    
    if barrier_backward is not None:
        labels.append("Backward Barrier (‡→R)")
        values.append(barrier_backward)
        colors.append("orange")
    
    bars = ax.bar(labels, values, color=colors, alpha=0.7)
    
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.annotate(
            f"{val:.1f} kcal/mol",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center", va="bottom",
            fontsize=10
        )
    
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title("S2.2 Path Search Energy Profile")
    ax.axhline(y=0, color="black", linestyle="-", linewidth=0.5)
    ax.grid(True, axis="y", alpha=0.3)
    
    plt.tight_layout()
    
    if output_path is None:
        output_path = path_profile_json.with_suffix(".png")
    
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    logger.info(f"Path search profile plot saved to: {output_path}")
    return output_path


def find_ts_and_dipole_guess(
    distances: List[float],
    energies: List[float],
    energies_in_hartree: bool = True,
) -> Tuple[Optional[float], Optional[float]]:
    """
    自动从势能面扫描数据中寻找 TS 初猜和偶极子中间体初猜。
    
    算法原理:
    1. Kneedle 算法: 归一化后找最大偏离基准弦的点 (TS 拐点)
    2. 最小梯度法: 在 TS 和峰之间找最平缓的点 (偶极子)
    
    Args:
        distances: 键长数据 (Angstrom), 必须从小到大排序
        energies: 能量数据
        energies_in_hartree: 如果为 True, 转换为 kcal/mol
    
    Returns:
        (ts_distance, dipole_distance): TS初猜键长和偶极子初猜键长 (Angstrom)
    """
    if not distances or not energies or len(distances) != len(energies):
        logger.warning("Invalid scan data for knee point detection")
        return None, None
    
    x = np.array(distances)
    e = np.array(energies)
    
    if energies_in_hartree:
        e = e * HARTREE_TO_KCAL
    
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    e = e[sort_idx]
    
    if len(x) < 3:
        logger.warning("Insufficient scan points for knee detection")
        return None, None
    
    peak_idx = np.argmax(e)
    x_peak = x[peak_idx]
    
    if peak_idx == 0:
        logger.warning("Peak at first point, cannot detect knee")
        return None, None
    
    x_roi = x[:peak_idx + 1]
    e_roi = e[:peak_idx + 1]
    
    if len(x_roi) < 3:
        return x[peak_idx], x[peak_idx]
    
    x_norm = (x_roi - np.min(x_roi)) / (np.max(x_roi) - np.min(x_roi) + 1e-10)
    e_norm = (e_roi - np.min(e_roi)) / (np.max(e_roi) - np.min(e_roi) + 1e-10)
    
    p1 = np.array([x_norm[0], e_norm[0]])
    p2 = np.array([x_norm[-1], e_norm[-1]])
    
    line_vec = p2 - p1
    line_len = np.linalg.norm(line_vec)
    if line_len < 1e-10:
        return x_roi[-1], x_roi[-1]
    
    line_unit = line_vec / line_len
    
    distances_from_line = []
    for i in range(len(x_norm)):
        p0 = np.array([x_norm[i], e_norm[i]])
        point_vec = p0 - p1
        cross = np.abs(line_unit[0] * point_vec[1] - line_unit[1] * point_vec[0])
        dist = cross
        distances_from_line.append(dist)
    
    ts_idx = np.argmax(distances_from_line)
    x_ts = x_roi[ts_idx]
    
    if ts_idx < peak_idx - 1:
        x_dipole_region = x_roi[ts_idx:peak_idx]
        e_dipole_region = e_roi[ts_idx:peak_idx]
        
        if len(x_dipole_region) >= 2:
            gradients = np.gradient(e_dipole_region, x_dipole_region)
            min_grad_idx = np.argmin(np.abs(gradients))
            x_dipole = x_dipole_region[min_grad_idx]
        else:
            x_dipole = (x_ts + x_peak) / 2.0
    else:
        x_dipole = (x_ts + x_peak) / 2.0
    
    logger.info(f"Knee point detection: TS={x_ts:.3f}Å, Peak={x_peak:.3f}Å, Dipole={x_dipole:.3f}Å")
    
    return float(x_ts), float(x_dipole)
