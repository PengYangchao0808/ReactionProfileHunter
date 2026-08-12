"""Shared stale-output archival helpers for resumable stage directories."""

from __future__ import annotations

import logging
import shutil
from datetime import datetime, timezone
from pathlib import Path


logger = logging.getLogger(__name__)


def archive_stale_stage_outputs(
    stage_output_dir: Path,
    *,
    stage_name: str,
    preserve: list[str] | None = None,
    run_id: str | None = None,
) -> Path | None:
    """Move stale stage outputs into a superseded archive directory.

    Returns the archive directory when entries were moved successfully, or
    ``None`` when nothing needed archiving or when archival failed.
    """

    stage_output_dir = Path(stage_output_dir)
    preserved_names = set(preserve or [])
    normalized_run_id = str(run_id).strip() if run_id is not None else ""

    try:
        if not stage_output_dir.exists():
            return None
        if not stage_output_dir.is_dir():
            logger.warning(
                "Stage %s: cannot archive stale outputs because %s is not a directory",
                stage_name,
                stage_output_dir,
            )
            return None

        stale_entries = [
            path
            for path in sorted(stage_output_dir.iterdir(), key=lambda item: item.name)
            if path.name != "superseded" and path.name not in preserved_names
        ]
        if not stale_entries:
            return None

        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        archive_name = (
            f"{stamp}__run-{normalized_run_id}" if normalized_run_id else stamp
        )
        archive_dir = stage_output_dir / "superseded" / archive_name
        archive_dir.mkdir(parents=True, exist_ok=False)

        moved_count = 0
        for source in stale_entries:
            try:
                shutil.move(str(source), str(archive_dir / source.name))
                moved_count += 1
            except Exception as exc:
                logger.warning(
                    "Stage %s: failed while archiving stale entry %s to %s after moving %d entries: %s",
                    stage_name,
                    source,
                    archive_dir,
                    moved_count,
                    exc,
                )
                return None

        logger.info(
            "Stage %s: archived %d stale entries to %s",
            stage_name,
            moved_count,
            archive_dir,
        )
        return archive_dir
    except Exception as exc:
        logger.warning(
            "Stage %s: unable to archive stale outputs in %s: %s",
            stage_name,
            stage_output_dir,
            exc,
        )
        return None
