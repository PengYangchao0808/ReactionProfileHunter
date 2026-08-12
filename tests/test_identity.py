from __future__ import annotations

from pathlib import Path

from rph_core.utils.identity import classify_int, classify_int_v2, classify_minimum, classify_ts


TS_EXPECTED_KEYS = {
    "hessian_index",
    "curvature_class",
    "mode_identity",
    "stationary_point_class",
}


def test_classify_ts_strict_target() -> None:
    result = classify_ts(
        frequencies_cm1=[-100.0, 125.0],
        forming_bonds=[(0, 1)],
        normal_modes={
            "coordinates": [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            "displacements": [(-0.2, 0.0, 0.0), (0.2, 0.0, 0.0)],
        },
    )

    assert set(result.keys()) == TS_EXPECTED_KEYS
    assert result == {
        "hessian_index": 1,
        "curvature_class": "strict",
        "mode_identity": "target",
        "stationary_point_class": "valid_target_ts",
    }


def test_classify_ts_soft_target() -> None:
    result = classify_ts(
        frequencies_cm1=[-30.0, 150.0],
        forming_bonds=[(0, 1)],
        normal_modes={
            "coordinates": [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            "displacements": [(-0.2, 0.0, 0.0), (0.2, 0.0, 0.0)],
        },
    )

    assert result["hessian_index"] == 1
    assert result["curvature_class"] == "soft"
    assert result["mode_identity"] == "target"
    assert result["stationary_point_class"] == "soft_target_ts"


def test_classify_ts_wrong_mode() -> None:
    result = classify_ts(
        frequencies_cm1=[-120.0, 90.0],
        forming_bonds=[(0, 1)],
        normal_modes={
            "coordinates": [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            "displacements": [(0.0, -0.2, 0.0), (0.0, 0.2, 0.0)],
        },
    )

    assert result["mode_identity"] == "unrelated"
    assert result["stationary_point_class"] == "first_order_wrong_mode"


def test_classify_ts_higher_order_saddle() -> None:
    result = classify_ts(
        frequencies_cm1=[-120.0, -80.0, 95.0],
        forming_bonds=[(0, 1)],
        normal_modes=None,
    )

    assert result["hessian_index"] == 2
    assert result["curvature_class"] == "multi_imaginary"
    assert result["stationary_point_class"] == "higher_order_saddle"


def test_classify_ts_minimum_after_optts() -> None:
    result = classify_ts(
        frequencies_cm1=[25.0, 150.0],
        forming_bonds=[(0, 1)],
        normal_modes=None,
    )

    assert result["hessian_index"] == 0
    assert result["curvature_class"] == "none"
    assert result["stationary_point_class"] == "minimum_after_optts"


def test_classify_ts_no_normal_modes() -> None:
    result = classify_ts(
        frequencies_cm1=[-120.0, 50.0],
        forming_bonds=[(0, 1)],
        normal_modes=None,
    )

    assert result["mode_identity"] == "unavailable"


def test_classify_int_distinct_intermediate(tmp_path: Path) -> None:
    opt_xyz = _write_xyz(
        tmp_path / "opt.xyz",
        [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
    )
    precursor_ref = _write_xyz(
        tmp_path / "precursor.xyz",
        [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)],
    )
    product_ref = _write_xyz(
        tmp_path / "product.xyz",
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
    )

    result = classify_int(
        opt_xyz=opt_xyz,
        forming_bonds=[(0, 1)],
        precursor_ref=precursor_ref,
        product_ref=product_ref,
        atom_mapping=None,
        frequencies_cm1=[20.0, 40.0],
    )

    assert result["identity"] == "distinct_intermediate"
    assert result["avg_progress"] == 0.5
    assert result["rmsd_to_product"] is not None and result["rmsd_to_product"] > 0.3
    assert result["rmsd_to_precursor"] is not None and result["rmsd_to_precursor"] > 0.3


def test_classify_int_collapsed_to_product() -> None:
    result = classify_int(
        opt_xyz=None,
        forming_bonds=[(0, 1)],
        precursor_ref=None,
        product_ref=None,
        atom_mapping=None,
        frequencies_cm1=[10.0],
        coordinates=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
        precursor_coordinates=[(0.0, 0.0, 0.0), (2.5, 0.0, 0.0)],
        product_coordinates=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
    )

    assert result["identity"] == "collapsed_to_product"
    assert result["rmsd_to_product"] == 0.0


def test_classify_int_imaginary_frequency() -> None:
    result = classify_int(
        opt_xyz=None,
        forming_bonds=[(0, 1)],
        precursor_ref=None,
        product_ref=None,
        atom_mapping=None,
        frequencies_cm1=[-20.0, 100.0],
        coordinates=[(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
    )

    assert result["identity"] == "imaginary_frequency"
    assert result["imag_count"] == 1


def test_classify_int_no_hardcoded_1p7_threshold() -> None:
    result = classify_int(
        opt_xyz=None,
        forming_bonds=[(0, 1)],
        precursor_ref=None,
        product_ref=None,
        atom_mapping=None,
        frequencies_cm1=[25.0],
        coordinates=[(0.0, 0.0, 0.0), (1.6, 0.0, 0.0)],
        precursor_coordinates=[(0.0, 0.0, 0.0), (1.6, 0.0, 0.0)],
        product_coordinates=[(0.0, 0.0, 0.0), (0.6, 0.0, 0.0)],
    )

    assert result["per_bond_distance"][(0, 1)] == 1.6
    assert result["identity"] == "collapsed_to_precursor"


def test_classify_minimum_valid() -> None:
    result = classify_minimum(
        frequencies_cm1=[20.0, 100.0],
        expected_role="product",
        rmsd_to_expected=0.1,
    )

    assert result["identity"] == "valid_minimum"
    assert result["identity_status"] == "role_matched"
    assert result["imaginary_count"] == 0


def test_classify_minimum_imaginary() -> None:
    result = classify_minimum(
        frequencies_cm1=[-20.0, 50.0],
        expected_role="product",
        rmsd_to_expected=0.1,
    )

    assert result["identity"] == "imaginary_frequency"
    assert result["identity_status"] == "not_checked"
    assert result["imaginary_count"] == 1


def test_classify_minimum_identity_drift() -> None:
    result = classify_minimum(
        frequencies_cm1=[20.0, 50.0],
        expected_role="precursor",
        rmsd_to_expected=0.4,
    )

    assert result["identity"] == "identity_drift"
    assert result["identity_status"] == "mismatched"


def test_classify_int_v2_missing_refs_never_dipolar() -> None:
    result = classify_int_v2(
        forming_bonds=[(0, 1)],
        coordinates=[(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)],
        frequencies_cm1=[20.0, 40.0],
        energy_int_hartree=-1.0,
    )

    assert result.identity == "unclassified"
    assert result.usable_for_ml is False
    assert "product_ref" in result.missing_evidence
    assert result.as_dict()["classification_version"] == "int_identity_v2"


def test_classify_int_v2_collapsed_by_rmsd() -> None:
    result = classify_int_v2(
        forming_bonds=[(0, 1)],
        coordinates=[(0.0, 0.0, 0.0), (1.55, 0.0, 0.0)],
        product_coordinates=[(0.0, 0.0, 0.0), (1.53, 0.0, 0.0)],
        frequencies_cm1=[20.0, 40.0],
        energy_int_hartree=-1.0,
        energy_product_hartree=-1.0005,
    )

    assert result.identity == "collapsed_to_product"
    assert result.usable_for_ml is False


def test_classify_int_v2_collapsed_by_formed_bonds() -> None:
    result = classify_int_v2(
        forming_bonds=[(0, 1), (2, 3)],
        coordinates=[(0.0, 0.0, 0.0), (1.55, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 4.54, 0.0)],
        product_coordinates=[(0.0, 0.0, 0.0), (1.53, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 4.56, 0.0)],
        frequencies_cm1=[20.0, 40.0],
        energy_int_hartree=-1.0,
        energy_product_hartree=-1.0005,
    )

    assert result.identity == "collapsed_to_product"


def test_classify_int_v2_dipolar_intermediate() -> None:
    result = classify_int_v2(
        forming_bonds=[(0, 1)],
        coordinates=[(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)],
        product_coordinates=[(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)],
        ts_coordinates=[(0.0, 0.0, 0.0), (2.9, 0.0, 0.0)],
        frequencies_cm1=[20.0, 40.0],
        energy_int_hartree=-1.0,
        energy_ts_hartree=-0.998,
        energy_product_hartree=-1.1,
    )

    assert result.identity == "dipolar_intermediate"
    assert result.usable_for_ml is True
    assert result.metrics["well_depth_kcal"] > 0.0
    assert result.metrics["above_product_kcal"] > 60.0


def test_classify_int_v2_merged_with_ts() -> None:
    result = classify_int_v2(
        forming_bonds=[(0, 1)],
        coordinates=[(0.0, 0.0, 0.0), (2.9, 0.0, 0.0)],
        product_coordinates=[(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)],
        ts_coordinates=[(0.0, 0.0, 0.0), (2.91, 0.0, 0.0)],
        frequencies_cm1=[20.0, 40.0],
        energy_int_hartree=-1.0,
        energy_ts_hartree=-1.0001,
        energy_product_hartree=-1.1,
    )

    assert result.identity == "merged_with_ts"
    assert result.usable_for_ml is False


def test_classify_int_v2_above_ts_not_minimum() -> None:
    result = classify_int_v2(
        forming_bonds=[(0, 1)],
        coordinates=[(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)],
        product_coordinates=[(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)],
        ts_coordinates=[(0.0, 0.0, 0.0), (2.9, 0.0, 0.0)],
        frequencies_cm1=[20.0, 40.0],
        energy_int_hartree=-1.0,
        energy_ts_hartree=-1.002,
        energy_product_hartree=-1.1,
    )

    assert result.identity == "above_ts_not_minimum"
    assert result.usable_for_ml is False


def test_classify_int_v2_imaginary_frequency() -> None:
    result = classify_int_v2(
        forming_bonds=[(0, 1)],
        coordinates=[(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)],
        product_coordinates=[(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)],
        frequencies_cm1=[-20.0, 40.0],
        energy_int_hartree=-1.0,
        energy_product_hartree=-1.1,
    )

    assert result.identity == "imaginary_frequency"
    assert result.usable_for_ml is False


def _write_xyz(path: Path, coordinates: list[tuple[float, float, float]]) -> Path:
    lines = [str(len(coordinates)), "test"]
    for x, y, z in coordinates:
        lines.append(f"C {x:.6f} {y:.6f} {z:.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
