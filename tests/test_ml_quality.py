from rph_core.utils.ml_quality import (
    compute_legacy_usable_for_ml,
    compute_ml_usability,
)


EXPECTED_KEYS = {
    "geometry",
    "electronic_energy",
    "enthalpy",
    "gibbs_free_energy",
    "frequency_descriptors",
    "ts_descriptors",
    "intermediate_descriptors",
    "mechanism_label",
    "multifidelity_pair",
}


def test_compute_ml_usability_complete_ts() -> None:
    result = compute_ml_usability(
        requested_role="ts",
        requested_kind="ts",
        opt_converged=True,
        sp_complete=True,
        freq_available=True,
        identity_status="role_matched",
        ts_frequency_valid=True,
        geometry_consistent=True,
    )

    assert set(result.keys()) == EXPECTED_KEYS
    assert result["geometry"] is True
    assert result["electronic_energy"] is True
    assert result["enthalpy"] is True
    assert result["gibbs_free_energy"] is True
    assert result["frequency_descriptors"] is True
    assert result["ts_descriptors"] is True
    assert result["intermediate_descriptors"] is False


def test_compute_ml_usability_failed_opt() -> None:
    result = compute_ml_usability(
        requested_role="product",
        requested_kind="minimum",
        opt_converged=False,
        sp_complete=False,
        freq_available=False,
        identity_status="not_checked",
        ts_frequency_valid=None,
        geometry_consistent=None,
    )

    assert set(result.keys()) == EXPECTED_KEYS
    assert result["geometry"] is False
    assert result["electronic_energy"] is False
    assert result["enthalpy"] is False
    assert result["gibbs_free_energy"] is False
    assert result["frequency_descriptors"] is False
    assert result["ts_descriptors"] is False
    assert result["intermediate_descriptors"] is False


def test_compute_ml_usability_intermediate_role_matched() -> None:
    result = compute_ml_usability(
        requested_role="intermediate",
        requested_kind="minimum",
        opt_converged=True,
        sp_complete=True,
        freq_available=False,
        identity_status="role_matched",
        ts_frequency_valid=None,
        geometry_consistent=None,
    )

    assert result["intermediate_descriptors"] is True


def test_compute_ml_usability_intermediate_mismatched() -> None:
    result = compute_ml_usability(
        requested_role="intermediate",
        requested_kind="minimum",
        opt_converged=True,
        sp_complete=True,
        freq_available=False,
        identity_status="mismatched",
        ts_frequency_valid=None,
        geometry_consistent=None,
    )

    assert result["intermediate_descriptors"] is False


def test_compute_ml_usability_mechanism_label_always_true() -> None:
    result = compute_ml_usability(
        requested_role="intermediate",
        requested_kind="minimum",
        opt_converged=False,
        sp_complete=False,
        freq_available=False,
        identity_status="ambiguous",
        ts_frequency_valid=False,
        geometry_consistent=False,
    )

    assert result["mechanism_label"] is True


def test_compute_ml_usability_multifidelity_pair_always_false() -> None:
    result = compute_ml_usability(
        requested_role="ts",
        requested_kind="ts",
        opt_converged=True,
        sp_complete=True,
        freq_available=True,
        identity_status="role_matched",
        ts_frequency_valid=True,
        geometry_consistent=True,
    )

    assert result["multifidelity_pair"] is False


def test_compute_legacy_usable_for_ml_mapping() -> None:
    ml_usability = {
        "geometry": True,
        "electronic_energy": True,
        "enthalpy": False,
        "gibbs_free_energy": False,
        "frequency_descriptors": False,
        "ts_descriptors": False,
        "intermediate_descriptors": False,
        "mechanism_label": True,
        "multifidelity_pair": False,
    }

    assert compute_legacy_usable_for_ml(ml_usability) is True
