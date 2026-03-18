from pathlib import Path

import numpy as np

from rph_core.steps.step2_retro.geometry_guard import check_scan_trajectory
from rph_core.utils.file_io import write_xyz


def test_check_scan_trajectory_reports_off_path_frames(tmp_path: Path) -> None:
    symbols = ["C", "C", "C", "C"]
    product_coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.45, 0.0, 0.0],
            [2.90, 0.0, 0.0],
            [4.35, 0.0, 0.0],
        ]
    )

    frame_ok = tmp_path / "frame_ok.xyz"
    frame_bad = tmp_path / "frame_bad.xyz"
    write_xyz(frame_ok, product_coords, symbols, title="ok")

    bad_coords = product_coords.copy()
    bad_coords[2] = np.array([1.05, 0.0, 0.0])
    write_xyz(frame_bad, bad_coords, symbols, title="bad")

    result = check_scan_trajectory(
        product_coords=product_coords,
        symbols=symbols,
        forming_bonds=((0, 1), (2, 3)),
        frame_paths=[frame_ok, frame_bad],
        graph_scale=1.25,
    )

    assert result["checked"] is True
    assert result["total_frames"] == 2
    assert result["off_path_count"] == 1
    assert result["off_path_indices"] == [1]
