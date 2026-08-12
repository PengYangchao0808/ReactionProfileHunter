from __future__ import annotations

import json
import logging
import os
import tempfile
from pathlib import Path
from typing import Any, Callable

logger = logging.getLogger(__name__)


def write_json_atomic(
    path: Path,
    payload: Any,
    *,
    indent: int = 2,
    compact: bool = False,
    default: Callable[[Any], Any] = str,
    encoding: str = "utf-8",
) -> Path:
    if compact:
        content = json.dumps(payload, separators=(",", ":"), default=default)
    else:
        content = json.dumps(payload, indent=indent, default=default)
    return write_text_atomic(path, content, encoding=encoding)


def write_text_atomic(path: Path, content: str, encoding: str = "utf-8") -> Path:
    path = Path(path)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding=encoding,
            dir=path.parent,
            prefix=f"{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            logger.debug(
                "Atomic write source path %s for destination %s",
                temporary_path,
                path,
            )
            handle.write(content)
        os.replace(temporary_path, path)
        return path
    except Exception:
        if temporary_path is not None:
            try:
                temporary_path.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                pass
        raise
