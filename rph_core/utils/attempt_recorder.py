"""Helpers for immutable per-attempt diagnostic directories."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rph_core.utils.json_io import write_json_atomic


@dataclass(frozen=True)
class AttemptRecord:
    attempt_id: str
    attempt_kind: str
    retry_level: int
    directory: Path
    input_inp_path: Path
    input_xyz_path: Path
    output_out_path: Path
    stderr_log_path: Path
    last_geometry_xyz_path: Path
    runtime_json_path: Path
    error_summary_path: Path | None


class AttemptRecorder:
    def __init__(self, structure_dir: Path, *, structure_id: str):
        self.structure_dir = Path(structure_dir)
        self.structure_id = str(structure_id)
        self._diagnostics_dir = self.structure_dir / "diagnostics"
        self._attempts: list[AttemptRecord] = []
        self._counter = self._existing_attempt_counter()

    def _existing_attempt_counter(self) -> int:
        """Continue attempt numbering when a failed structure is rerun in place."""

        if not self._diagnostics_dir.is_dir():
            return 0
        counters = []
        for directory in self._diagnostics_dir.iterdir():
            match = re.fullmatch(r"attempt_(\d+)_.*", directory.name)
            if directory.is_dir() and match is not None:
                counters.append(int(match.group(1)))
        return max(counters, default=0)

    def next_attempt(self, kind: str, retry_level: int = 0) -> AttemptRecord:
        self._counter += 1
        self._diagnostics_dir.mkdir(parents=True, exist_ok=True)
        attempt_id = f"attempt_{self._counter:03d}_{str(kind).strip()}"
        attempt_dir = self._diagnostics_dir / attempt_id
        attempt_dir.mkdir(parents=True, exist_ok=False)
        record = AttemptRecord(
            attempt_id=attempt_id,
            attempt_kind=str(kind).strip(),
            retry_level=int(retry_level),
            directory=attempt_dir,
            input_inp_path=attempt_dir / "input.inp",
            input_xyz_path=attempt_dir / "input.xyz",
            output_out_path=attempt_dir / "output.out",
            stderr_log_path=attempt_dir / "stderr.log",
            last_geometry_xyz_path=attempt_dir / "last_geometry.xyz",
            runtime_json_path=attempt_dir / "runtime.json",
            error_summary_path=attempt_dir / "error_summary.json",
        )
        self._attempts.append(record)
        return record

    def finalize_attempt(
        self,
        record: AttemptRecord,
        *,
        runtime: dict[str, Any],
        error_summary: dict[str, Any] | None = None,
        last_geometry_xyz: str | None = None,
    ) -> None:
        write_json_atomic(record.runtime_json_path, runtime)
        if error_summary is not None and record.error_summary_path is not None:
            write_json_atomic(record.error_summary_path, error_summary)
        if last_geometry_xyz is not None:
            payload = last_geometry_xyz if last_geometry_xyz.endswith("\n") else f"{last_geometry_xyz}\n"
            record.last_geometry_xyz_path.write_text(payload, encoding="utf-8")

    def history(self) -> list[AttemptRecord]:
        return list(self._attempts)
