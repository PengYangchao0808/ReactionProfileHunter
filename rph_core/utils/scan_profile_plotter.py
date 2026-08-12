"""Single-artifact plotting entry point for the S2 refined PATH profile."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, List, Mapping, Optional

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
    *,
    profile_payload: Optional[Mapping[str, Any]] = None,
) -> Optional[Path]:
    """Render the sole S2 figure as ``scan_profile.png`` beside its manifest.

    ``profile_payload`` skips the manifest re-read when the caller already holds
    the in-memory payload (single-write path in the S2 engine).
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
    except ImportError:
        logger.warning("matplotlib not available; skipping S2 profile rendering")
        return None

    from rph_core.utils.s2_profile_figures import render_s2_profile_figure

    scan_profile_json = Path(scan_profile_json)
    png_path = (
        scan_profile_json.parent / "scan_profile.png"
        if output_path is None
        else Path(output_path).with_suffix(".png")
    )
    try:
        if profile_payload is not None:
            profile_payload = dict(profile_payload)
        else:
            if not scan_profile_json.is_file():
                logger.error("S2 scan profile is missing: %s", scan_profile_json)
                return None
            profile_payload = json.loads(scan_profile_json.read_text(encoding="utf-8"))
        # The rescue profile is written beside the QC artifacts, but the
        # variant owns one canonical figure.  Inject the rescue payload into
        # the main profile only for rendering; the JSON schema remains
        # unchanged and no second PNG is published.
        rescue = dict(profile_payload.get("rescue") or {})
        rescue_profile_path = rescue.get("scan_profile")
        if rescue_profile_path:
            rescue_path = _resolve_profile_path(
                str(rescue_profile_path), scan_profile_json.parent
            )
            if rescue_path.is_file() and rescue_path.resolve() != scan_profile_json.resolve():
                try:
                    profile_payload["_rescue_profile_for_plot"] = json.loads(
                        rescue_path.read_text(encoding="utf-8")
                    )
                except (OSError, ValueError) as exc:
                    logger.warning("Could not load rescue profile for combined S2 figure %s: %s", rescue_path, exc)
        _attach_new_rule_selection(profile_payload)
        render_s2_profile_figure(profile_payload, png_path)
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


def _resolve_profile_path(value: str, anchor: Path) -> Path:
    """Resolve native and WSL-style artifact paths on either host."""
    path = Path(value)
    if not path.is_absolute():
        path = Path(anchor) / path
    if path.is_file():
        return path
    normalized = value.replace("\\", "/")
    if normalized.lower().startswith("/mnt/") and len(normalized) > 6:
        drive = normalized[5]
        remainder = normalized[6:].lstrip("/")
        windows_path = Path(f"{drive.upper()}:/{remainder}")
        if windows_path.is_file():
            return windows_path
    return path


def _attach_new_rule_selection(profile: dict[str, Any]) -> None:
    """Attach canonical new-rule seeds for plotting, replaying legacy scans."""
    rescue_profile = profile.get("_rescue_profile_for_plot")
    if not isinstance(rescue_profile, Mapping):
        return
    persisted = profile.get("selection_replay")
    if isinstance(persisted, Mapping):
        display = dict(rescue_profile)
        _apply_selection_to_display_profile(display, persisted)
        profile["_rescue_profile_for_plot"] = display
        return
    rescue = dict(rescue_profile.get("rescue") or {})
    existing_ts = dict(rescue.get("ts_search_seed") or {})
    existing_knee = dict(rescue_profile.get("knee_evidence") or {})
    if existing_ts and existing_knee:
        replay = {
            "mode": "selector",
            "algorithm": "endpoint_knee_shift_midpoint_v1",
            "s2_state": rescue.get("s2_state") or rescue_profile.get("s2_state"),
            "seed_evidence": rescue.get("seed_evidence") or rescue_profile.get("seed_evidence"),
            "ts_search_seed": existing_ts,
            "int_search_seed": dict(rescue.get("int_search_seed") or rescue_profile.get("int_search_seed") or {}),
            "endpoint_evidence": dict(rescue_profile.get("endpoint_evidence") or {}),
            "knee_evidence": existing_knee,
            "diagnostics": dict(rescue_profile.get("selection_diagnostics") or {}),
            "selection_source": "orca_relaxed_scan",
        }
    else:
        replay = _run_offline_rescue_replay(profile, rescue_profile)
    if not replay:
        return
    profile["selection_replay"] = replay
    display = dict(rescue_profile)
    _apply_selection_to_display_profile(display, replay)
    profile["_rescue_profile_for_plot"] = display


def _run_offline_rescue_replay(
    parent_profile: Mapping[str, Any],
    rescue_profile: Mapping[str, Any],
) -> Optional[dict[str, Any]]:
    """Replay saved rescue frames using the unified selector only."""
    from rph_core.steps.step2_retro.path_selector import replay_rescue_selection

    rescue = dict(rescue_profile.get("rescue") or {})
    raw_frames = list(rescue.get("frames") or [])
    energies = list(rescue.get("energies_hartree") or [])
    product = parent_profile.get("product_xyz")
    if not raw_frames or len(raw_frames) != len(energies) or not product:
        return None
    frames = [str(_resolve_profile_path(str(value), Path.cwd())) for value in raw_frames]
    product_path = _resolve_profile_path(str(product), Path.cwd())
    forming_bonds = [
        tuple(int(atom) for atom in pair)
        for pair in (parent_profile.get("forming_bonds") or [])
        if len(pair) == 2
    ]
    if not forming_bonds:
        forming_bonds = [
            tuple(int(atom) for atom in item.get("atoms", ()))
            for item in (rescue.get("coordinates") or [])
            if len(item.get("atoms", ())) == 2
        ]
    if not forming_bonds or not product_path.is_file() or any(
        not Path(frame).is_file() for frame in frames
    ):
        return None
    replay_payload = dict(rescue)
    replay_payload["frames"] = frames
    replay_payload["energies_hartree"] = energies
    replay = replay_rescue_selection(
        replay_payload,
        forming_bonds=forming_bonds,
        product_xyz=product_path,
    )
    if replay.get("s2_state") != "rescue_seeded":
        return None
    return replay


def _apply_selection_to_display_profile(
    display: dict[str, Any],
    selection: Mapping[str, Any],
) -> None:
    """Replace display-only selection fields with new-rule possibilities."""
    ts = dict(selection.get("ts_search_seed") or {})
    intermediate = dict(selection.get("int_search_seed") or {})
    display.update(
        {
            "s2_state": selection.get("s2_state"),
            "seed_evidence": selection.get("seed_evidence"),
            "selection_source": selection.get("selection_source") or "orca_relaxed_scan",
            "ts_search_seed": ts or None,
            "int_search_seed": intermediate or None,
            "endpoint_evidence": dict(selection.get("endpoint_evidence") or {}),
            "knee_evidence": dict(selection.get("knee_evidence") or {}),
            "selection_policy": {"actual_source": "B97-3c"},
            "selection_decision": {
                "selected_branch": "b973c_relaxed_scan",
                "rule": "knee + adaptive right shift; INT = TS-to-effective-endpoint midpoint",
            },
            "selections": {
                "ts_guess": {
                    **ts,
                    "index": ts.get("frame_index"),
                    "rule": "knee + adaptive right shift",
                },
                "intermediate": {
                    **intermediate,
                    "index": intermediate.get("frame_index"),
                    "rule": "TS-to-effective-endpoint midpoint",
                },
            },
        }
    )
