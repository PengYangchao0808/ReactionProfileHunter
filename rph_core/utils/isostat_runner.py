import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class IsostatResult:
    cluster_xyz: Path
    log_file: Path
    n_within_window: Optional[int]


def _count_within_window(log_file: Path, energy_window: float) -> Optional[int]:
    try:
        content = log_file.read_text(errors="ignore").splitlines()
    except Exception:
        return None
    count = 0
    for line in content:
        if "DE=" not in line:
            continue
        parts = line.replace("=", " = ").split()
        for i, token in enumerate(parts):
            if token == "DE" and i + 1 < len(parts):
                try:
                    val = float(parts[i + 1])
                    if val <= energy_window:
                        count += 1
                except ValueError:
                    continue
    return count if count > 0 else None


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
    cmd = [
        str(isostat_bin),
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
