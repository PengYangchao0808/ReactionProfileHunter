"""Single-artifact plotting entry point for the S2 refined PATH profile."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Optional

import numpy as np

logger = logging.getLogger(__name__)

HARTREE_TO_KCAL = 627.509


def compute_scan_distances(
    start_distance: float,
    end_distance: float,
    num_steps: int,
    direction: str = "outward",
) -> List[float]:
    """Return the coarse-scan coordinate convention used by S2."""
    start, end = (
        (end_distance, start_distance)
        if direction == "outward"
        else (start_distance, end_distance)
    )
    return np.linspace(start, end, num_steps).tolist()


def plot_scan_profile(
    scan_profile_json: Path,
    output_path: Optional[Path] = None,
) -> Optional[Path]:
    """Render the sole S2 figure as ``scan_profile.png`` beside its manifest."""
    try:
        import matplotlib

        matplotlib.use("Agg")
    except ImportError:
        logger.warning("matplotlib not available; skipping S2 profile rendering")
        return None

    from rph_core.utils.s2_profile_figures import render_s2_profile_figure

    scan_profile_json = Path(scan_profile_json)
    if not scan_profile_json.is_file():
        logger.error("S2 scan profile is missing: %s", scan_profile_json)
        return None
    png_path = (
        scan_profile_json.parent / "scan_profile.png"
        if output_path is None
        else Path(output_path).with_suffix(".png")
    )
    try:
        render_s2_profile_figure(
            json.loads(scan_profile_json.read_text(encoding="utf-8")),
            png_path,
        )
        _remove_legacy_profile_artifacts(png_path)
        logger.info("Saved S2 profile figure: %s", png_path)
        return png_path
    except Exception:
        logger.warning("Failed to render S2 profile from %s", scan_profile_json, exc_info=True)
        return None


def _remove_legacy_profile_artifacts(png_path: Path) -> None:
    """Remove only known redundant plot artifacts after a PNG was written."""
    if not png_path.is_file():
        return
    legacy_paths = [png_path.with_suffix(".pdf"), png_path.with_suffix(".svg")]
    outputs = png_path.parent / "outputs"
    legacy_paths.extend(
        outputs / relative for relative in (
            Path("diagnostics/s2_path_selection_diagnostic.pdf"),
            Path("diagnostics/s2_path_selection_diagnostic.png"),
            Path("publication/s2_path_profile_publication.pdf"),
            Path("publication/s2_path_profile_publication.svg"),
            Path("publication/s2_path_profile_publication.png"),
        )
    )
    for path in legacy_paths:
        try:
            path.unlink()
        except FileNotFoundError:
            pass
        except OSError as exc:
            logger.warning("Could not remove obsolete S2 plot %s: %s", path, exc)
    for directory in (outputs / "diagnostics", outputs / "publication", outputs):
        try:
            directory.rmdir()
        except (FileNotFoundError, OSError):
            pass
