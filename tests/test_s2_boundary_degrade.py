from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, cast

import numpy as np
import pytest

from rph_core.steps.step2_retro.retro_scanner import RetroScanner
from rph_core.utils.file_io import write_xyz


def _base_config() -> Dict[str, object]:
    return {
        "step2": {
            "scan": {
                "scan_start_distance": 3.2,
                "scan_end_distance": 1.8,
                "scan_steps": 8,
                "scan_mode": "concerted",
                "scan_force_constant": 0.5,
                "min_valid_points": 5,
                "reject_boundary_maximum": True,
                "boundary_retry_once": True,
                "boundary_retry_delta": 0.3,
                "boundary_retry_extra_steps": 4,
                "allow_boundary_degradation": True,
                "require_local_peak": False,
                "intermediate_validation_strict": False,
            },
            "xtb_settings": {"solvent": "acetone"},
            "intermediate_optimization": {"charge": 0, "multiplicity": 1},
        },
        "resources": {"nproc": 1},
        "theory": {"optimization": {"solvent": "acetone"}},
    }


def test_boundary_max_retry_then_degrade(monkeypatch, tmp_path: Path) -> None:
    scanner = RetroScanner(_base_config())

    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    symbols = ["C", "C", "C", "C"]

    product_xyz = tmp_path / "product.xyz"
    write_xyz(product_xyz, coords, symbols, title="product")

    intermediate_xyz = tmp_path / "intermediate.xyz"
    write_xyz(intermediate_xyz, coords, symbols, title="intermediate")

    ts_guess_primary = tmp_path / "ts_guess_primary.xyz"
    ts_guess_retry = tmp_path / "ts_guess_retry.xyz"
    write_xyz(ts_guess_primary, coords, symbols, title="ts1")
    write_xyz(ts_guess_retry, coords, symbols, title="ts2")

    monkeypatch.setattr(scanner, "_resolve_product_file", lambda p: Path(p))
    monkeypatch.setattr(
        scanner,
        "_optimize_intermediate",
        lambda seed, out_dir, forming_bonds, scan_start_distance: intermediate_xyz,
    )
    monkeypatch.setattr(scanner.bond_stretcher, "stretch_bonds", lambda c, _: c)

    from rph_core.utils.geometry_tools import LogParser

    monkeypatch.setattr(
        LogParser,
        "extract_last_converged_coords",
        staticmethod(lambda *_args, **_kwargs: (coords, symbols, None)),
    )

    scan_calls = {"count": 0}
    captured_constraints = []

    def _mock_execute_scan(self, start_xyz, output_dir, bonds, params, direction, charge=0, spin=1):
        scan_calls["count"] += 1
        # Return format: (result, energies, max_idx, boundary_max, local_peak)
        if scan_calls["count"] == 1:
            # First call: boundary maximum at index 4 (last position)
            result = SimpleNamespace(
                success=True,
                energies=[0.1, 0.2, 0.3, 0.35, 0.4],
                ts_guess_xyz=ts_guess_primary,
                geometries=[intermediate_xyz]
            )
            return result, [0.1, 0.2, 0.3, 0.35, 0.4], 4, True, False
        else:
            # Retry call: boundary maximum still persists at index 4
            result = SimpleNamespace(
                success=True,
                energies=[0.1, 0.2, 0.25, 0.29, 0.40],
                ts_guess_xyz=ts_guess_retry,
                geometries=[intermediate_xyz]
            )
            return result, [0.1, 0.2, 0.25, 0.29, 0.40], 4, True, False  # boundary max persists

    monkeypatch.setattr(RetroScanner, "_execute_scan", _mock_execute_scan)

    (
        ts_guess_xyz,
        reactant_xyz,
        dipolar_xyz,
        forming_bonds,
        scan_profile_json,
        status,
        confidence,
        degraded_reasons,
    ) = scanner.run(
        product_xyz=product_xyz,
        output_dir=tmp_path / "S2_Retro",
        forming_bonds=((0, 1), (2, 3)),
        scan_config=None,
    )

    assert scan_calls["count"] == 2
    assert ts_guess_xyz.exists()
    assert reactant_xyz.exists()
    assert dipolar_xyz.exists()
    assert forming_bonds == ((0, 1), (2, 3))
    assert scan_profile_json.exists()
    assert status == "DEGRADED"
    assert confidence == "low"
    assert "boundary_maximum_persisted_after_retry" in degraded_reasons


def test_intermediate_quality_warning_v51_behavior(monkeypatch, tmp_path: Path, caplog) -> None:
    """
    V5.1: S2 no longer fails on intermediate quality issues.
    Instead, it logs warnings and defers refinement to S3 (TSOptimizer).
    """
    cfg = _base_config()
    scanner = RetroScanner(cfg)

    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    symbols = ["C", "C", "C", "C"]

    product_xyz = tmp_path / "product.xyz"
    write_xyz(product_xyz, coords, symbols, title="product")

    # Create intermediate with unphysically short bonds to trigger warning
    intermediate_coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.2, 0.0, 0.0],  # Bond 0-1 = 0.2A (unphysically short)
            [0.0, 0.2, 0.0],  # Bond 0-2 = 0.2A (unphysically short)
            [0.0, 0.0, 0.2],
        ]
    )
    intermediate_xyz = tmp_path / "intermediate.xyz"
    write_xyz(intermediate_xyz, intermediate_coords, symbols, title="intermediate")

    monkeypatch.setattr(scanner, "_resolve_product_file", lambda p: Path(p))
    monkeypatch.setattr(
        scanner,
        "_optimize_intermediate",
        lambda seed, out_dir, forming_bonds, scan_start_distance: intermediate_xyz,
    )
    monkeypatch.setattr(scanner.bond_stretcher, "stretch_bonds", lambda c, _: c)

    from rph_core.utils.geometry_tools import LogParser

    monkeypatch.setattr(
        LogParser,
        "extract_last_converged_coords",
        staticmethod(lambda *_args, **_kwargs: (coords, symbols, None)),
    )

    # Mock the scan to avoid full execution
    scan_calls: Dict[str, Any] = {"count": 0}

    def mock_scan(*args, **kwargs):
        scan_calls["count"] += 1
        from types import SimpleNamespace

        class MockResult:
            success = True
            energies = [1.0, 2.0, 1.5]
            geometries = [intermediate_xyz]
            max_energy_index = 1
            ts_guess_xyz = intermediate_xyz
            scan_log = None

        return MockResult(), [1.0, 2.0, 1.5], 1, False, True

    monkeypatch.setattr(scanner, "_execute_scan", mock_scan)

    caplog.set_level("WARNING")

    # V5.1: Should NOT raise - just warn and continue
    (
        ts_guess_xyz,
        reactant_xyz,
        dipolar_xyz,
        forming_bonds,
        scan_profile_json,
        status,
        confidence,
        degraded_reasons,
    ) = scanner.run(
        product_xyz=product_xyz,
        output_dir=tmp_path / "S2_Retro",
        forming_bonds=((0, 1), (2, 3)),
        scan_config=None,
    )

    # Verify warning was logged
    assert any(
        "Intermediate geometry warning" in rec.message for rec in caplog.records
    ), "Expected warning about intermediate geometry"

    # Verify scan completed (did not raise)
    assert scan_calls["count"] == 1
    assert ts_guess_xyz.exists()


def test_forming_bond_preflight_warns_for_unusually_long_distance(monkeypatch, tmp_path: Path, caplog) -> None:
    scanner = RetroScanner(_base_config())

    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [4.2, 0.0, 0.0],
        ]
    )
    symbols = ["C", "C", "C", "C"]

    product_xyz = tmp_path / "product.xyz"
    write_xyz(product_xyz, coords, symbols, title="product")

    monkeypatch.setattr(scanner, "_resolve_product_file", lambda p: Path(p))
    monkeypatch.setattr(
        scanner,
        "_optimize_intermediate",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("stop_after_preflight")),
    )

    from rph_core.utils.geometry_tools import LogParser

    monkeypatch.setattr(
        LogParser,
        "extract_last_converged_coords",
        staticmethod(lambda *_args, **_kwargs: (coords, symbols, None)),
    )

    caplog.set_level("WARNING")
    with pytest.raises(RuntimeError, match="stop_after_preflight"):
        scanner.run(
            product_xyz=product_xyz,
            output_dir=tmp_path / "S2_Retro",
            forming_bonds=((0, 3), (1, 2)),
            scan_config=None,
        )

    assert any(
        "possible index mapping error" in rec.message and "Forming bond (0, 3)" in rec.message
        for rec in caplog.records
    )
