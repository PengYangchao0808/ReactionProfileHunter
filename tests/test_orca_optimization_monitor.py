from __future__ import annotations

from rph_core.utils.orca_optimization_monitor import (
    OptimizationMonitorConfig,
    detect_oscillation,
    parse_orca_optimization_cycles,
)


def _cycle(cycle: int, energy_change: float, rms_gradient: float, max_step: float) -> str:
    return f"""
 * GEOMETRY OPTIMIZATION CYCLE {cycle} *
 Actually observed energy change .... {energy_change:.8f}
 ----------------------|Geometry convergence|-------------------------
 RMS gradient        {rms_gradient:.8f}            0.00010000      NO
 MAX gradient        0.00800000                    0.00030000      NO
 RMS step            0.01000000                    0.00200000      NO
 MAX step            {max_step:.8f}                0.00400000      NO
 """


def test_monitor_detects_confirmed_oscillation_only_after_full_window():
    content = "".join(
        _cycle(index, 2.0e-4 if index % 2 else -2.0e-4, 0.0010, 0.030)
        for index in range(1, 17)
    )
    records = parse_orca_optimization_cycles(content)

    trigger = detect_oscillation(records, OptimizationMonitorConfig())

    assert len(records) == 16
    assert trigger is not None
    assert trigger["reason"] == "oscillation"
    assert trigger["energy_reversals"] == 15


def test_monitor_does_not_stop_slow_but_progressing_optimization():
    records = [
        {
            "cycle": index,
            "energy_change": 2.0e-4,
            "rms_gradient": 0.0040 - index * 0.00015,
            "max_gradient": 0.008,
            "rms_step": 0.003,
            "max_step": 0.030,
            "rms_gradient_tolerance": 0.0001,
            "max_gradient_tolerance": 0.0003,
            "rms_step_tolerance": 0.002,
            "max_step_tolerance": 0.004,
        }
        for index in range(1, 39)
    ]

    assert detect_oscillation(records, OptimizationMonitorConfig()) is None


def test_parser_ignores_orca_diagnostic_ellipsis_lines():
    content = """
 * GEOMETRY OPTIMIZATION CYCLE 1 *
 RMS gradient                       ...    0.0064224164
 MAX gradient                       ...    0.0309534407
 RMS step                           ...    0.0109108945
 MAX step                           ...    0.0409846704
 Actually observed energy change         ....    -0.000356758
 ----------------------|Geometry convergence|-------------------------
 Energy change      -0.0003567585            0.0000050000      NO
 RMS gradient        0.0026297157            0.0001000000      NO
 MAX gradient        0.0285361003            0.0003000000      NO
 RMS step            0.0109108945            0.0020000000      NO
 MAX step            0.0422101151            0.0040000000      NO
 """

    records = parse_orca_optimization_cycles(content)

    assert len(records) == 1
    assert records[0]["energy_change"] == -0.0003567585
    assert records[0]["rms_gradient"] == 0.0026297157
    assert records[0]["max_gradient"] == 0.0285361003
    assert records[0]["rms_step"] == 0.0109108945
    assert records[0]["max_step"] == 0.0422101151
