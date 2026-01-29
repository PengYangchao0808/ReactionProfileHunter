"""
Path compatibility helpers.

Handles Windows-style paths when running under Linux/WSL and provides
basic detection of toxic path characters.
"""

from __future__ import annotations

import os
import re
from pathlib import Path


_WINDOWS_PATH_RE = re.compile(r"^([A-Za-z]):[\\/](.*)$")


def normalize_path(path: str | Path) -> Path:
    """Normalize a path string to a pathlib.Path.

    - If running on non-Windows and the input is a Windows drive path,
      convert it to /mnt/<drive>/... for WSL-style environments.
    - Otherwise return Path(path) directly.
    """
    if isinstance(path, Path):
        return path

    if os.name != "nt":
        match = _WINDOWS_PATH_RE.match(path)
        if match:
            drive = match.group(1).lower()
            rest = match.group(2).replace("\\", "/")
            return Path(f"/mnt/{drive}/{rest}")

    return Path(path)


def is_toxic_path(path: Path) -> bool:
    """Return True if path contains characters known to break QC tooling."""
    path_str = str(path)
    for ch in (" ", "[", "]", "(", ")", "{", "}"):
        if ch in path_str:
            return True
    return False
