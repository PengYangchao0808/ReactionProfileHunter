"""Helpers for classifying and reviewing soft TS imaginary modes."""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Protocol, Sequence

import numpy as np

from rph_core.utils.file_io import read_xyz, write_xyz


logger = logging.getLogger(__name__)


class SoftModeReviewCalculator(Protocol):
    def _soft_mode_review_max_attempts(self) -> int: ...

    def _soft_mode_review_strict_cutoff(self) -> float: ...

    def _soft_mode_review_soft_window(self) -> tuple[float, float]: ...

    def _run_soft_mode_review_optimization(
        self,
        input_xyz: Path,
        *,
        direction: str,
    ) -> dict[str, Any]: ...

    def _run_soft_mode_review_frequency(
        self,
        input_xyz: Path,
        *,
        direction: str,
    ) -> dict[str, Any]: ...

TS_FREQUENCY_CLASS_STRICT = "strict"
TS_FREQUENCY_CLASS_SOFT = "soft"
TS_FREQUENCY_CLASS_INVALID = "invalid"
TS_FREQUENCY_CLASS_UNVERIFIED = "unverified"
TS_FREQUENCY_CLASS_NOT_APPLICABLE = "not_applicable"


def normalize_soft_mode_window(
    soft_mode_window_cm1: Sequence[float] | None,
    strict_cutoff_cm1: float,
) -> tuple[float, float]:
    """Normalize the soft-mode window bounds as (low, high)."""

    raw_window = soft_mode_window_cm1 or (strict_cutoff_cm1, -10.0)
    if len(raw_window) != 2:
        raise ValueError("soft_mode_window_cm1 must contain exactly two values")
    soft_low = float(raw_window[0])
    soft_high = float(raw_window[1])
    if soft_low > soft_high:
        soft_low, soft_high = soft_high, soft_low
    return soft_low, soft_high


def classify_ts_frequencies(
    frequencies_cm1: Sequence[float],
    *,
    strict_cutoff_cm1: float = -50.0,
    soft_mode_window_cm1: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Classify TS frequencies into strict / soft / invalid buckets."""

    soft_low, soft_high = normalize_soft_mode_window(
        soft_mode_window_cm1,
        strict_cutoff_cm1,
    )
    all_imaginaries = tuple(float(value) for value in frequencies_cm1 if float(value) < 0.0)
    strict_imaginaries = tuple(
        float(value) for value in frequencies_cm1 if float(value) <= strict_cutoff_cm1
    )
    soft_imaginaries = tuple(
        float(value)
        for value in frequencies_cm1
        if soft_low < float(value) <= soft_high
    )

    if len(all_imaginaries) == 1 and len(strict_imaginaries) == 1 and len(soft_imaginaries) == 0:
        ts_frequency_class = TS_FREQUENCY_CLASS_STRICT
        ts_frequency_valid = True
    elif len(all_imaginaries) == 1 and len(strict_imaginaries) == 0 and len(soft_imaginaries) == 1:
        ts_frequency_class = TS_FREQUENCY_CLASS_SOFT
        ts_frequency_valid = False
    else:
        ts_frequency_class = TS_FREQUENCY_CLASS_INVALID
        ts_frequency_valid = False

    return {
        "ts_frequency_class": ts_frequency_class,
        "ts_frequency_valid": ts_frequency_valid,
        "all_imaginaries_cm1": all_imaginaries,
        "strict_imaginaries_cm1": strict_imaginaries,
        "soft_imaginaries_cm1": soft_imaginaries,
        "soft_mode_window_cm1": (soft_low, soft_high),
    }


@dataclass(frozen=True)
class SoftModeReviewResult:
    reviewed: bool
    reviewed_at: str
    mode_alignment_score: float | None
    perturbation_attempts: tuple[dict[str, Any], ...]
    final_class: str
    eligible_for_s4_reoptimization: bool
    reason: str

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["perturbation_attempts"] = [dict(item) for item in self.perturbation_attempts]
        return payload


def review_soft_mode(
    *,
    calculator: SoftModeReviewCalculator,
    structure: Mapping[str, Any],
    structure_dir: Path,
    forming_bonds: Sequence[Sequence[int]],
    imaginary_frequency_cm1: float,
    frequency_output: Path,
    optimized_xyz: Path,
    mode_displacement_details: Mapping[str, Any] | None,
    projection_threshold: float = 0.3,
    perturbation_scale_bohr: float = 0.1,
) -> SoftModeReviewResult:
    """Perturb a soft TS guess along ±mode, then rerun OptTS/FREQ."""

    reviewed_at = _timestamp()
    try:
        if not Path(frequency_output).exists():
            return _not_reviewed(
                reviewed_at,
                _extract_alignment_score(mode_displacement_details or {}),
                f"Soft-mode review skipped because frequency output was missing: {frequency_output}",
            )
        details = dict(mode_displacement_details or {})
        if details.get("mode_displacement_valid") is not True:
            return _not_reviewed(
                reviewed_at,
                _extract_alignment_score(details),
                "Soft-mode review skipped because mode displacement alignment is unavailable or invalid.",
            )

        atomic_displacements = details.get("atomic_displacements")
        if not isinstance(atomic_displacements, Sequence) or not atomic_displacements:
            return _not_reviewed(
                reviewed_at,
                _extract_alignment_score(details),
                "Soft-mode review skipped because atomic displacements were not available.",
            )

        base_coordinates, symbols = read_xyz(Path(optimized_xyz))
        displacement_array = np.asarray(atomic_displacements, dtype=float)
        if displacement_array.shape != base_coordinates.shape:
            return _not_reviewed(
                reviewed_at,
                _extract_alignment_score(details),
                "Soft-mode review skipped because displacement and geometry shapes do not match.",
            )

        review_dir = Path(structure_dir) / "soft_mode_review"
        review_dir.mkdir(parents=True, exist_ok=True)
        max_attempts = max(0, int(calculator._soft_mode_review_max_attempts()))
        directions = ("+", "-")[:max_attempts]
        if not directions:
            return _not_reviewed(
                reviewed_at,
                _extract_alignment_score(details),
                "Soft-mode review skipped because max_perturbation_attempts is zero.",
            )

        attempt_payloads: list[dict[str, Any]] = []
        for direction in directions:
            sign = 1.0 if direction == "+" else -1.0
            perturbed_xyz = review_dir / f"perturb_{'plus' if direction == '+' else 'minus'}.xyz"
            perturbed_coordinates = base_coordinates + (sign * float(perturbation_scale_bohr) * displacement_array)
            write_xyz(
                perturbed_xyz,
                perturbed_coordinates,
                symbols,
                title=f"soft-mode perturbation {direction} from {structure.get('id', 'unknown')}",
            )

            opt_outcome = calculator._run_soft_mode_review_optimization(
                perturbed_xyz,
                direction=direction,
            )
            opt_result = opt_outcome["result"]
            freq_status = "not_run"
            freq_class_after = TS_FREQUENCY_CLASS_INVALID

            if opt_result.status == "complete":
                freq_input = Path(opt_result.output_xyz) if opt_result.output_xyz else perturbed_xyz
                freq_outcome = calculator._run_soft_mode_review_frequency(
                    freq_input,
                    direction=direction,
                )
                freq_result = freq_outcome["result"]
                freq_status = freq_result.status
                if freq_result.status == "complete" and freq_result.frequencies_cm1:
                    freq_class_after = str(
                        classify_ts_frequencies(
                            freq_result.frequencies_cm1,
                            strict_cutoff_cm1=calculator._soft_mode_review_strict_cutoff(),
                            soft_mode_window_cm1=calculator._soft_mode_review_soft_window(),
                        )["ts_frequency_class"]
                    )

            attempt_payloads.append(
                {
                    "direction": direction,
                    "opt_status": opt_result.status,
                    "freq_status": freq_status,
                    "freq_class_after": freq_class_after,
                }
            )

        alignment_score = _extract_alignment_score(details)
        if any(item["freq_class_after"] == TS_FREQUENCY_CLASS_STRICT for item in attempt_payloads):
            return SoftModeReviewResult(
                reviewed=True,
                reviewed_at=reviewed_at,
                mode_alignment_score=alignment_score,
                perturbation_attempts=tuple(attempt_payloads),
                final_class=TS_FREQUENCY_CLASS_STRICT,
                eligible_for_s4_reoptimization=True,
                reason=(
                    f"Soft TS review promoted {structure.get('id', 'structure')} to strict after perturbation."
                ),
            )

        if (
            attempt_payloads
            and all(item["freq_class_after"] == TS_FREQUENCY_CLASS_SOFT for item in attempt_payloads)
            and alignment_score is not None
            and alignment_score >= float(projection_threshold)
        ):
            return SoftModeReviewResult(
                reviewed=True,
                reviewed_at=reviewed_at,
                mode_alignment_score=alignment_score,
                perturbation_attempts=tuple(attempt_payloads),
                final_class=TS_FREQUENCY_CLASS_SOFT,
                eligible_for_s4_reoptimization=True,
                reason=(
                    "Soft TS remained soft after perturbation, but the imaginary mode stayed aligned with the forming bonds."
                ),
            )

        return SoftModeReviewResult(
            reviewed=True,
            reviewed_at=reviewed_at,
            mode_alignment_score=alignment_score,
            perturbation_attempts=tuple(attempt_payloads),
            final_class=TS_FREQUENCY_CLASS_UNVERIFIED,
            eligible_for_s4_reoptimization=False,
            reason=(
                "Soft TS perturbation review did not recover a strict TS or a consistently aligned soft TS."
            ),
        )
    except Exception as exc:  # pragma: no cover - defensive degradation path
        logger.warning(
            "Soft-mode review failed for %s with forming bonds %s and imaginary frequency %.3f cm-1: %s",
            structure.get("id", "structure"),
            list(forming_bonds or ()),
            float(imaginary_frequency_cm1),
            exc,
        )
        return _not_reviewed(
            reviewed_at,
            _extract_alignment_score(mode_displacement_details or {}),
            f"Soft-mode review failed and was ignored: {exc}",
        )


def _extract_alignment_score(details: Mapping[str, Any]) -> float | None:
    raw_value = details.get("mode_alignment_score")
    if raw_value is None:
        return None
    try:
        return float(raw_value)
    except (TypeError, ValueError):
        return None


def _not_reviewed(
    reviewed_at: str,
    alignment_score: float | None,
    reason: str,
) -> SoftModeReviewResult:
    return SoftModeReviewResult(
        reviewed=False,
        reviewed_at=reviewed_at,
        mode_alignment_score=alignment_score,
        perturbation_attempts=(),
        final_class=TS_FREQUENCY_CLASS_UNVERIFIED,
        eligible_for_s4_reoptimization=False,
        reason=reason,
    )


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()
