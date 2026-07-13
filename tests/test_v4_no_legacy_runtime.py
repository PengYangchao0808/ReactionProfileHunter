"""Regression guard for the supported V4 S0--S4 runtime surface."""

from __future__ import annotations

import ast
from pathlib import Path

import yaml


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CORE_ROOT = REPOSITORY_ROOT / "rph_core"
FORBIDDEN_MODULES = {
    "rph_core.scheduling.v3_scheduler",
    "rph_core.utils.checkpoint_manager",
    "rph_core.utils.gau_xtb_interface",
    "rph_core.utils.geometry_preprocessor",
    "rph_core.utils.intra_reaction_scheduler",
    "rph_core.utils.oscillation_detector",
    "rph_core.utils.qc_task_runner",
    "rph_core.utils.v3_progress",
    "rph_core.utils.v3_stage_display",
    "rph_core.steps.anchor",
    "rph_core.steps.conformer_search.candidates",
    "rph_core.steps.conformer_search.engine",
    "rph_core.steps.conformer_search.funnel",
    "rph_core.steps.conformer_search.pipeline",
    "rph_core.steps.conformer_search.state_manager",
    "rph_core.steps.step2_retro.kinematic_stretcher",
    "rph_core.steps.step2_retro.retro_scanner",
    "rph_core.steps.step3_opt",
}
FORBIDDEN_CONFIG_KEYS = {"gau_xtb", "neutral_precursor", "path_search", "preoptimization"}


def _imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
    return names


def _config_keys(value: object) -> set[str]:
    if isinstance(value, dict):
        return set(value) | set().union(*(_config_keys(item) for item in value.values()))
    if isinstance(value, list):
        return set().union(*(_config_keys(item) for item in value)) if value else set()
    return set()


def test_v4_runtime_has_no_v3_module_imports() -> None:
    imports = set().union(*(_imports(path) for path in CORE_ROOT.rglob("*.py")))
    forbidden = sorted(
        imported
        for imported in imports
        if any(imported == module or imported.startswith(f"{module}.") for module in FORBIDDEN_MODULES)
    )
    assert not forbidden, f"V4 runtime imports retired V3 modules: {forbidden}"


def test_v4_defaults_exclude_retired_s2_and_preoptimization_options() -> None:
    defaults = yaml.safe_load((REPOSITORY_ROOT / "config" / "defaults.yaml").read_text(encoding="utf-8"))
    assert not (_config_keys(defaults) & FORBIDDEN_CONFIG_KEYS)
