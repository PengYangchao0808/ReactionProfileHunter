"""
Unified Thermochemistry Module
===============================

Provides a single source of truth for Shermo invocation, .sum parsing,
thermo.json generation, and cache self-healing.

This module supersedes the scattered Shermo parsing logic that previously
lived in handler.py, small_molecule_precompute.py, and inline orchestrator code.

Key classes:
    ThermoRecord   — canonical thermodynamic record (kcal/mol + Hartree)
    ShermoOptions  — Shermo CLI parameter bundle with config-based factory

Key functions:
    parse_shermo_sum()              — parse any Shermo .sum file
    run_shermo_record()             — invoke Shermo, return ThermoRecord
    write_thermo_json()             — write canonical thermo.json
    read_thermo_json()              — read thermo.json (v0 + v1 schema)
    derive_thermo_json_from_sum()   — .sum → thermo.json
    ensure_thermo_json_from_entry() — cache self-heal
    derive_shermo_summary()         — .sum → shermo_summary.json
"""

from rph_core.utils.thermo.schema import (
    ThermoRecord,
    write_thermo_json,
    read_thermo_json,
)
from rph_core.utils.thermo.options import ShermoOptions
from rph_core.utils.thermo.shermo import (
    parse_shermo_sum,
    run_shermo_record,
)
from rph_core.utils.thermo.artifacts import (
    derive_thermo_json_from_sum,
    ensure_thermo_json_from_entry,
    derive_shermo_summary,
)

__all__ = [
    "ThermoRecord",
    "ShermoOptions",
    "parse_shermo_sum",
    "run_shermo_record",
    "write_thermo_json",
    "read_thermo_json",
    "derive_thermo_json_from_sum",
    "ensure_thermo_json_from_entry",
    "derive_shermo_summary",
]
