"""Degradation tests: Step4 should not crash on missing artifacts.

These tests are intentionally extract-only (V6.1 FeatureMiner) and should not
depend on external QC backends.
"""

from pathlib import Path
import tempfile


def _write_min_xyz(path: Path) -> None:
    path.write_text("1\n\nH 0.0 0.0 0.0\n")


def test_step4_generates_features_without_fchk() -> None:
    from rph_core.steps.step4_features.feature_miner import FeatureMiner

    with tempfile.TemporaryDirectory(prefix="test_degradation_") as tmp:
        tmpdir = Path(tmp)

        ts = tmpdir / "ts_final.xyz"
        reactant = tmpdir / "reactant.xyz"
        product = tmpdir / "product_min.xyz"
        _write_min_xyz(ts)
        _write_min_xyz(reactant)
        _write_min_xyz(product)

        out_dir = tmpdir / "S4_Data"
        miner = FeatureMiner(config={})
        features_raw_csv = miner.run(
            ts_final=ts,
            reactant=reactant,
            product=product,
            output_dir=out_dir,
            forming_bonds=None,
            fragment_indices=None,
            sp_matrix_report=None,
            ts_fchk=None,
            ts_orca_out=None,
            reactant_fchk=None,
            reactant_orca_out=None,
            product_fchk=None,
            product_orca_out=None,
        )

        assert features_raw_csv.exists()
        assert (out_dir / "features_mlr.csv").exists()
        assert (out_dir / "feature_meta.json").exists()


def test_step4_handles_missing_forming_bonds() -> None:
    from rph_core.steps.step4_features.feature_miner import FeatureMiner

    with tempfile.TemporaryDirectory(prefix="test_degradation_") as tmp:
        tmpdir = Path(tmp)

        ts = tmpdir / "ts_final.xyz"
        reactant = tmpdir / "reactant.xyz"
        product = tmpdir / "product_min.xyz"
        _write_min_xyz(ts)
        _write_min_xyz(reactant)
        _write_min_xyz(product)

        out_dir = tmpdir / "S4_Data"
        miner = FeatureMiner(config={})
        features_raw_csv = miner.run(
            ts_final=ts,
            reactant=reactant,
            product=product,
            output_dir=out_dir,
            forming_bonds=None,
            fragment_indices=None,
            sp_matrix_report=None,
        )

        assert features_raw_csv.exists()
