from __future__ import annotations

from pathlib import Path

import pytest

from rph_core.utils.config_loader import load_config
from rph_core.utils.qc_models import QCJobResult
from rph_core.utils.soft_mode_review import (
    SoftModeReviewResult,
    TS_FREQUENCY_CLASS_SOFT,
    TS_FREQUENCY_CLASS_STRICT,
    TS_FREQUENCY_CLASS_UNVERIFIED,
    classify_ts_frequencies,
    normalize_soft_mode_window,
    review_soft_mode,
)


def _write_ts_xyz(path: Path) -> Path:
    path.write_text(
        "2\nts\n"
        "C 0.000000 0.000000 0.000000\n"
        "O 1.200000 0.000000 0.000000\n",
        encoding="utf-8",
    )
    return path


class _FakeReviewCalculator:
    def __init__(
        self,
        *,
        opt_results: dict[str, QCJobResult],
        freq_results: dict[str, QCJobResult],
        max_attempts: int = 2,
    ) -> None:
        self._opt_results = opt_results
        self._freq_results = freq_results
        self._max_attempts = max_attempts

    def _soft_mode_review_max_attempts(self) -> int:
        return self._max_attempts

    def _soft_mode_review_strict_cutoff(self) -> float:
        return -50.0

    def _soft_mode_review_soft_window(self) -> tuple[float, float]:
        return (-50.0, -10.0)

    def _run_soft_mode_review_optimization(self, input_xyz: Path, *, direction: str):
        result = self._opt_results[direction]
        if result.output_xyz is None:
            result.output_xyz = Path(input_xyz)
        return {"result": result}

    def _run_soft_mode_review_frequency(self, input_xyz: Path, *, direction: str):
        _ = input_xyz
        return {"result": self._freq_results[direction]}


def test_frequency_classification_boundaries() -> None:
    assert classify_ts_frequencies((-50.0, 120.0))["ts_frequency_class"] == TS_FREQUENCY_CLASS_STRICT
    assert classify_ts_frequencies((-10.0, 120.0))["ts_frequency_class"] == TS_FREQUENCY_CLASS_SOFT
    assert classify_ts_frequencies((-9.9, 120.0))["ts_frequency_class"] == "invalid"
    assert classify_ts_frequencies((-55.0, -20.0, 120.0))["ts_frequency_class"] == "invalid"


def test_soft_mode_review_result_serialization() -> None:
    result = SoftModeReviewResult(
        reviewed=True,
        reviewed_at="2026-07-18T00:00:00+00:00",
        mode_alignment_score=0.42,
        perturbation_attempts=({"direction": "+", "freq_class_after": "soft"},),
        final_class=TS_FREQUENCY_CLASS_SOFT,
        eligible_for_s4_reoptimization=True,
        reason="aligned",
    )

    payload = result.to_dict()

    assert payload["reviewed"] is True
    assert payload["mode_alignment_score"] == pytest.approx(0.42)
    assert payload["perturbation_attempts"][0]["direction"] == "+"


def test_review_soft_mode_promotes_strict_when_one_perturbation_recovers_ts(tmp_path: Path) -> None:
    optimized_xyz = _write_ts_xyz(tmp_path / "opt.xyz")
    frequency_output = tmp_path / "freq.out"
    frequency_output.write_text("placeholder\n", encoding="utf-8")
    calculator = _FakeReviewCalculator(
        opt_results={
            "+": QCJobResult("complete", optimized_xyz, output_xyz=optimized_xyz),
            "-": QCJobResult("complete", optimized_xyz, output_xyz=optimized_xyz),
        },
        freq_results={
            "+": QCJobResult("complete", optimized_xyz, frequencies_cm1=(-120.0, 100.0)),
            "-": QCJobResult("complete", optimized_xyz, frequencies_cm1=(-20.0, 100.0)),
        },
    )

    result = review_soft_mode(
        calculator=calculator,
        structure={"id": "ts_soft"},
        structure_dir=tmp_path,
        forming_bonds=((0, 1),),
        imaginary_frequency_cm1=-25.0,
        frequency_output=frequency_output,
        optimized_xyz=optimized_xyz,
        mode_displacement_details={
            "mode_displacement_valid": True,
            "mode_alignment_score": 0.5,
            "atomic_displacements": [[0.2, 0.0, 0.0], [-0.2, 0.0, 0.0]],
        },
    )

    assert result.reviewed is True
    assert result.final_class == TS_FREQUENCY_CLASS_STRICT
    assert result.eligible_for_s4_reoptimization is True
    assert result.perturbation_attempts[0]["freq_class_after"] == TS_FREQUENCY_CLASS_STRICT


def test_review_soft_mode_keeps_aligned_soft_ts_s4_eligible(tmp_path: Path) -> None:
    optimized_xyz = _write_ts_xyz(tmp_path / "opt_soft.xyz")
    frequency_output = tmp_path / "freq_soft.out"
    frequency_output.write_text("placeholder\n", encoding="utf-8")
    calculator = _FakeReviewCalculator(
        opt_results={
            "+": QCJobResult("complete", optimized_xyz, output_xyz=optimized_xyz),
            "-": QCJobResult("complete", optimized_xyz, output_xyz=optimized_xyz),
        },
        freq_results={
            "+": QCJobResult("complete", optimized_xyz, frequencies_cm1=(-20.0, 100.0)),
            "-": QCJobResult("complete", optimized_xyz, frequencies_cm1=(-18.0, 120.0)),
        },
    )

    result = review_soft_mode(
        calculator=calculator,
        structure={"id": "ts_soft"},
        structure_dir=tmp_path,
        forming_bonds=((0, 1),),
        imaginary_frequency_cm1=-25.0,
        frequency_output=frequency_output,
        optimized_xyz=optimized_xyz,
        mode_displacement_details={
            "mode_displacement_valid": True,
            "mode_alignment_score": 0.35,
            "atomic_displacements": [[0.2, 0.0, 0.0], [-0.2, 0.0, 0.0]],
        },
    )

    assert result.reviewed is True
    assert result.final_class == TS_FREQUENCY_CLASS_SOFT
    assert result.eligible_for_s4_reoptimization is True


def test_review_soft_mode_marks_unverified_when_perturbations_do_not_help(tmp_path: Path) -> None:
    optimized_xyz = _write_ts_xyz(tmp_path / "opt_invalid.xyz")
    frequency_output = tmp_path / "freq_invalid.out"
    frequency_output.write_text("placeholder\n", encoding="utf-8")
    calculator = _FakeReviewCalculator(
        opt_results={
            "+": QCJobResult("complete", optimized_xyz, output_xyz=optimized_xyz),
            "-": QCJobResult("complete", optimized_xyz, output_xyz=optimized_xyz),
        },
        freq_results={
            "+": QCJobResult("complete", optimized_xyz, frequencies_cm1=(20.0, 100.0)),
            "-": QCJobResult("complete", optimized_xyz, frequencies_cm1=(-120.0, -20.0, 100.0)),
        },
    )

    result = review_soft_mode(
        calculator=calculator,
        structure={"id": "ts_soft"},
        structure_dir=tmp_path,
        forming_bonds=((0, 1),),
        imaginary_frequency_cm1=-25.0,
        frequency_output=frequency_output,
        optimized_xyz=optimized_xyz,
        mode_displacement_details={
            "mode_displacement_valid": True,
            "mode_alignment_score": 0.35,
            "atomic_displacements": [[0.2, 0.0, 0.0], [-0.2, 0.0, 0.0]],
        },
    )

    assert result.reviewed is True
    assert result.final_class == TS_FREQUENCY_CLASS_UNVERIFIED
    assert result.eligible_for_s4_reoptimization is False


def test_review_soft_mode_skips_when_mode_displacement_is_not_valid(tmp_path: Path) -> None:
    optimized_xyz = _write_ts_xyz(tmp_path / "opt_skip.xyz")
    frequency_output = tmp_path / "freq_skip.out"
    frequency_output.write_text("placeholder\n", encoding="utf-8")

    class FailIfCalledCalculator:
        def _soft_mode_review_max_attempts(self) -> int:
            return 2

        def _soft_mode_review_strict_cutoff(self) -> float:
            return -50.0

        def _soft_mode_review_soft_window(self) -> tuple[float, float]:
            return (-50.0, -10.0)

        def _run_soft_mode_review_optimization(self, input_xyz: Path, *, direction: str):
            _ = (input_xyz, direction)
            raise AssertionError("soft-mode review should not run")

        def _run_soft_mode_review_frequency(self, input_xyz: Path, *, direction: str):
            _ = (input_xyz, direction)
            raise AssertionError("soft-mode review should not run")

    result = review_soft_mode(
        calculator=FailIfCalledCalculator(),
        structure={"id": "ts_soft"},
        structure_dir=tmp_path,
        forming_bonds=((0, 1),),
        imaginary_frequency_cm1=-25.0,
        frequency_output=frequency_output,
        optimized_xyz=optimized_xyz,
        mode_displacement_details={"mode_displacement_valid": False},
    )

    assert result.reviewed is False
    assert result.eligible_for_s4_reoptimization is False


def test_config_soft_mode_window_parsing() -> None:
    config = load_config()

    s3_frequency = config["theory"]["s3_low_level"]["optimization"]["frequency"]
    s4_frequency = config["theory"]["s4_high_precision"]["optimization"]["frequency"]

    assert normalize_soft_mode_window(s3_frequency["soft_mode_window_cm1"], -50.0) == (-50.0, -10.0)
    assert s3_frequency["soft_mode_review"]["enabled"] is True
    assert s3_frequency["soft_mode_review"]["max_perturbation_attempts"] == 2
    assert normalize_soft_mode_window(s4_frequency["soft_mode_window_cm1"], -50.0) == (-50.0, -10.0)
