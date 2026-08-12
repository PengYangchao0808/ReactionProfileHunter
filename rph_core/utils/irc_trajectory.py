"""ORCA IRC trajectory parsing and intermediate-shoulder detection."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

_KCAL_PER_HARTREE = 627.509


@dataclass(frozen=True)
class IRCFrame:
    frame_index: int
    coords: Tuple[Tuple[float, float, float], ...]
    energy_hartree: Optional[float]
    source_file: Optional[Path]


@dataclass(frozen=True)
class IRCTrajectory:
    direction: str
    frames: List[IRCFrame] = field(default_factory=list)
    source_files: List[Path] = field(default_factory=list)

    @property
    def energies_hartree(self) -> List[float]:
        return [frame.energy_hartree for frame in self.frames if frame.energy_hartree is not None]

    @property
    def frame_count(self) -> int:
        return len(self.frames)


@dataclass(frozen=True)
class IRCShoulder:
    frame_index: int
    kind: str
    energy_kcal_mol: float
    source_file: Optional[Path]
    coords: Tuple[Tuple[float, float, float], ...]


def _energy_from_comment(comment: str) -> Optional[float]:
    """Extract the energy (Hartree) from an ORCA trajectory comment line."""
    tokens = comment.split()
    for token in tokens:
        if "energy" in token.lower() or "E_h" in token:
            cleaned = token.strip("=#(),")
            if cleaned.lower().startswith("energy"):
                cleaned = cleaned.split("=", 1)[-1].strip()
            try:
                return float(cleaned)
            except ValueError:
                continue
    for token in tokens:
        cleaned = token.strip("=#(),")
        try:
            value = float(cleaned)
        except ValueError:
            continue
        return value
    return None


def parse_irc_trajectory_file(path: Path) -> Optional[IRCTrajectory]:
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        return None
    frames: List[IRCFrame] = []
    index = 0
    while index < len(lines):
        stripped = lines[index].strip()
        if not stripped or stripped.startswith("---"):
            index += 1
            continue
        try:
            n_atoms = int(stripped)
        except ValueError:
            index += 1
            continue
        frame_end = index + n_atoms + 2
        if frame_end > len(lines):
            break
        comment = lines[index + 1].strip()
        coords: List[Tuple[float, float, float]] = []
        for row in lines[index + 2:frame_end]:
            tokens = row.split()
            if len(tokens) >= 4:
                try:
                    coords.append((float(tokens[1]), float(tokens[2]), float(tokens[3])))
                except ValueError:
                    continue
        if len(coords) == n_atoms:
            frames.append(
                IRCFrame(
                    frame_index=len(frames),
                    coords=tuple(coords),
                    energy_hartree=_energy_from_comment(comment),
                    source_file=path,
                )
            )
        index = frame_end
    if not frames:
        return None
    direction = "unknown"
    if "forward" in path.name.lower() or "fwd" in path.name.lower():
        direction = "forward"
    elif "backward" in path.name.lower() or "bwd" in path.name.lower():
        direction = "backward"
    return IRCTrajectory(direction=direction, frames=frames, source_files=[path])


def parse_irc_trajectories(irc_dir: Path) -> List[IRCTrajectory]:
    trajectories: List[IRCTrajectory] = []
    for candidate in sorted(Path(irc_dir).glob("*.xyz")):
        trajectory = parse_irc_trajectory_file(candidate)
        if trajectory is not None:
            trajectories.append(trajectory)
    return trajectories


def detect_intermediate_shoulder(
    trajectory: IRCTrajectory,
    energy_window_kcal_mol: float = 2.0,
) -> Optional[IRCShoulder]:
    """Search for an energy stationary point or shoulder between TS and endpoints.

    A stationary point is a local minimum (E[i] < E[i-1] and E[i] < E[i+1])
    that sits below the trajectory ends.  A shoulder is a sign change in the
    first energy difference (a flattening) whose energy lies within the window
    of the local baseline.
    """
    energies = trajectory.energies_hartree
    if len(energies) < 5:
        return None
    first_energy = energies[0]
    last_energy = energies[-1]
    baseline = min(first_energy, last_energy)

    for i in range(2, len(energies) - 2):
        current = energies[i]
        window = [
            energies[i - 2],
            energies[i - 1],
            energies[i + 1],
            energies[i + 2],
        ]
        if current < min(window):
            delta_kcal = (current - baseline) * _KCAL_PER_HARTREE
            if delta_kcal < energy_window_kcal_mol:
                return IRCShoulder(
                    frame_index=i,
                    kind="stationary_point",
                    energy_kcal_mol=delta_kcal,
                    source_file=trajectory.frames[i].source_file,
                    coords=trajectory.frames[i].coords,
                )

    differences = [
        (energies[j + 1] - energies[j]) * _KCAL_PER_HARTREE
        for j in range(len(energies) - 1)
    ]
    for i in range(1, len(differences) - 1):
        if differences[i - 1] * differences[i] < 0.0:
            shoulder_energy = (energies[i] - baseline) * _KCAL_PER_HARTREE
            if shoulder_energy < energy_window_kcal_mol:
                return IRCShoulder(
                    frame_index=i,
                    kind="shoulder",
                    energy_kcal_mol=shoulder_energy,
                    source_file=trajectory.frames[i].source_file,
                    coords=trajectory.frames[i].coords,
                )
    return None


def frame_to_xyz(
    trajectory: IRCTrajectory,
    frame_index: int,
    destination: Path,
    symbols: Sequence[str],
) -> Optional[Path]:
    if frame_index < 0 or frame_index >= trajectory.frame_count:
        return None
    frame = trajectory.frames[frame_index]
    if len(symbols) != len(frame.coords):
        return None
    lines = [str(len(frame.coords)), f"frame {frame.frame_index} from IRC"]
    for symbol, (x, y, z) in zip(symbols, frame.coords):
        lines.append(f"{symbol:<3} {x: .8f} {y: .8f} {z: .8f}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return destination


def classify_irc_endpoints(
    trajectories: Sequence[IRCTrajectory],
    forming_bonds: Sequence[Tuple[int, int]],
) -> Optional[Dict[str, Any]]:
    """Identify the stretched-end (dipolar-intermediate) vs product-end IRC endpoints.

    The reaction coordinate is the forming-bond distance: the product side
    converges to a short bonded distance (~1.5-1.7 A) while the stretched side
    ends at a larger distance.  We classify each trajectory's terminal frame by
    its mean forming-bond distance and return the last frame of the stretched
    trajectory as the new INT seed.
    """
    if not trajectories:
        return None
    candidates: List[Dict[str, Any]] = []
    for trajectory in trajectories:
        if not trajectory.frames:
            continue
        frame = trajectory.frames[-1]
        distances = []
        for atom_i, atom_j in forming_bonds:
            if atom_i < 0 or atom_j < 0 or atom_i >= len(frame.coords) or atom_j >= len(frame.coords):
                continue
            dx = frame.coords[atom_j][0] - frame.coords[atom_i][0]
            dy = frame.coords[atom_j][1] - frame.coords[atom_i][1]
            dz = frame.coords[atom_j][2] - frame.coords[atom_i][2]
            distances.append((dx * dx + dy * dy + dz * dz) ** 0.5)
        if not distances:
            continue
        mean_distance = sum(distances) / len(distances)
        candidates.append(
            {
                "trajectory": trajectory,
                "frame": frame,
                "mean_fb_distance": mean_distance,
                "distances": distances,
            }
        )
    if not candidates:
        return None
    stretched = max(candidates, key=lambda item: item["mean_fb_distance"])
    product = min(candidates, key=lambda item: item["mean_fb_distance"])
    return {
        "stretch_frame": stretched["frame"],
        "stretch_trajectory": stretched["trajectory"],
        "stretch_mean_fb_distance": stretched["mean_fb_distance"],
        "stretch_distances": stretched["distances"],
        "product_frame": product["frame"],
        "product_mean_fb_distance": product["mean_fb_distance"],
        "product_distances": product["distances"],
    }


__all__ = [
    "IRCShoulder",
    "IRCFrame",
    "IRCTrajectory",
    "classify_irc_endpoints",
    "detect_intermediate_shoulder",
    "frame_to_xyz",
    "parse_irc_trajectories",
    "parse_irc_trajectory_file",
]
