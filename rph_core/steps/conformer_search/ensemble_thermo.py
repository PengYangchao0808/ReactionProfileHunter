"""Canonical-ensemble thermodynamics for the V4 S1 conformer ensemble."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence, Tuple

from rph_core.utils.constants import HARTREE_TO_KCAL


GAS_CONSTANT_KCAL_MOL_K = 1.98720425864083e-3


@dataclass(frozen=True)
class EnsembleThermodynamics:
    """Partition-function result for one internally consistent ensemble."""

    temperature_k: float
    relative_free_energies_kcal: Tuple[float, ...]
    relative_partition_weights: Tuple[float, ...]
    boltzmann_populations: Tuple[float, ...]
    partition_function_relative: float
    conformational_free_energy_correction_kcal: float
    reference_free_energy_hartree: float
    ensemble_free_energy_hartree: float


def calculate_ensemble_thermodynamics(
    free_energies_hartree: Sequence[float],
    temperature_k: float,
    degeneracies: Sequence[int] | None = None,
) -> EnsembleThermodynamics:
    """Calculate ``Z_rel``, populations and ``G_conf_rel``.

    ``free_energies_hartree`` must contain comparable conformer free energies
    from one uniform method stack. Degeneracy defaults to one for every unique
    conformer; sampled duplicate frequency is deliberately not interpreted as
    physical degeneracy.
    """

    if not free_energies_hartree:
        raise ValueError("Cannot calculate ensemble thermodynamics for an empty ensemble")
    temperature = float(temperature_k)
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("temperature_k must be finite and greater than zero")

    energies = tuple(float(value) for value in free_energies_hartree)
    if not all(math.isfinite(value) for value in energies):
        raise ValueError("Conformer free energies must all be finite")

    if degeneracies is None:
        normalized_degeneracies = (1,) * len(energies)
    else:
        raw_degeneracies = tuple(degeneracies)
        if any(
            isinstance(value, bool)
            or not math.isfinite(float(value))
            or float(value) != int(value)
            for value in raw_degeneracies
        ):
            raise ValueError("Every conformer degeneracy must be a finite integer")
        normalized_degeneracies = tuple(int(value) for value in raw_degeneracies)
        if len(normalized_degeneracies) != len(energies):
            raise ValueError("degeneracies and free_energies_hartree must have equal length")
        if any(value < 1 for value in normalized_degeneracies):
            raise ValueError("Every conformer degeneracy must be at least one")

    reference = min(energies)
    relative = tuple((value - reference) * HARTREE_TO_KCAL for value in energies)
    rt = GAS_CONSTANT_KCAL_MOL_K * temperature
    log_weights = tuple(
        math.log(degeneracy) - delta_g / rt
        for degeneracy, delta_g in zip(normalized_degeneracies, relative)
    )
    max_log_weight = max(log_weights)
    scaled_weights = tuple(math.exp(value - max_log_weight) for value in log_weights)
    scaled_sum = math.fsum(scaled_weights)
    log_partition = max_log_weight + math.log(scaled_sum)
    partition = math.exp(log_partition)
    weights = tuple(math.exp(value) for value in log_weights)
    populations = tuple(value / scaled_sum for value in scaled_weights)
    g_conf_kcal = -rt * log_partition
    ensemble_free_energy = reference + g_conf_kcal / HARTREE_TO_KCAL

    return EnsembleThermodynamics(
        temperature_k=temperature,
        relative_free_energies_kcal=relative,
        relative_partition_weights=weights,
        boltzmann_populations=populations,
        partition_function_relative=partition,
        conformational_free_energy_correction_kcal=g_conf_kcal,
        reference_free_energy_hartree=reference,
        ensemble_free_energy_hartree=ensemble_free_energy,
    )
