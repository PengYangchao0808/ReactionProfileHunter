import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class IsostatResult:
    cluster_xyz: Path
    log_file: Path
    n_within_window: Optional[int]


logger = logging.getLogger(__name__)


def _count_within_window(log_file: Path, energy_window: float) -> Optional[int]:
    try:
        content = log_file.read_text(errors="ignore").splitlines()
    except Exception:
        return None
    count = 0
    for line in content:
        # Relaxed check for DE and = (can be separated by spaces)
        if "DE" not in line or "=" not in line:
            continue
        parts = line.replace("=", " = ").split()
        for i, token in enumerate(parts):
            if token == "DE":
                # Check next token (could be "=" or value)
                val_idx = i + 1
                if val_idx < len(parts) and parts[val_idx] == "=":
                    val_idx += 1
                
                if val_idx < len(parts):
                    try:
                        val = float(parts[val_idx])
                        if val <= energy_window:
                            count += 1
                    except ValueError:
                        continue
    return count if count > 0 else None


def _resolve_isostat_path(isostat_bin: Path) -> Optional[Path]:
    if isostat_bin.exists() and isostat_bin.is_file() and os.access(str(isostat_bin), os.X_OK):
        return isostat_bin
    which_path = shutil.which(str(isostat_bin))
    if which_path:
        return Path(which_path)
    return None


def run_isostat(
    isostat_bin: Path,
    input_xyz: Path,
    output_dir: Path,
    gdis: float,
    edis: float,
    temp_k: float,
    threads: int,
    energy_window: float
) -> IsostatResult:
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / "isostat.log"
    resolved_bin = _resolve_isostat_path(isostat_bin)
    if resolved_bin is None:
        message = (
            "isostat executable not found; skipping clustering. "
            "Set config.executables.isostat.path or ensure isostat on PATH."
        )
        logger.warning(message)
        log_file.write_text(message)
        return IsostatResult(
            cluster_xyz=input_xyz,
            log_file=log_file,
            n_within_window=None
        )
    cmd = [
        str(resolved_bin),
        str(input_xyz),
        "-Gdis",
        f"{gdis}",
        "-Edis",
        f"{edis}",
        "-nt",
        f"{threads}",
        "-T",
        f"{temp_k}"
    ]
    result = subprocess.run(
        cmd,
        cwd=output_dir,
        capture_output=True,
        text=True
    )
    log_file.write_text(result.stdout + result.stderr)
    cluster_xyz = output_dir / "cluster.xyz"
    if not cluster_xyz.exists():
        fallback = output_dir / "isomers.xyz"
        if fallback.exists():
            fallback.replace(cluster_xyz)
        else:
            raise RuntimeError(f"isostat 未生成 cluster.xyz: {log_file}")
    n_within = _count_within_window(log_file, energy_window)
    return IsostatResult(
        cluster_xyz=cluster_xyz,
        log_file=log_file,
        n_within_window=n_within
    )
