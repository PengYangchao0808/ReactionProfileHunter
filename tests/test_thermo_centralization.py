from __future__ import annotations

"""
Verification tests for Shermo centralization architecture.
Validates that all Shermo calls are centralized in the condition layer.

Author: RPH Team
Date: 2026-05-19
"""

import ast
import importlib
import re
from pathlib import Path
from typing import Callable, cast


REPO_ROOT = Path(__file__).resolve().parent.parent


def _repo_file(*parts: str) -> Path:
    return REPO_ROOT.joinpath(*parts)


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _iter_run_shermo_call_sites() -> set[str]:
    call_sites: set[str] = set()
    for path in sorted((REPO_ROOT / "rph_core").rglob("*.py")):
        if "backup" in path.name.lower():
            continue
        source = _read_text(path)
        tree = ast.parse(source, filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if isinstance(func, ast.Name) and func.id == "run_shermo":
                call_sites.add(path.relative_to(REPO_ROOT).as_posix())
                break
            if isinstance(func, ast.Attribute) and func.attr == "run_shermo":
                call_sites.add(path.relative_to(REPO_ROOT).as_posix())
                break
    return call_sites


def _function_source(path: Path, function_name: str) -> str:
    source = _read_text(path)
    tree = ast.parse(source, filename=str(path))
    lines = source.splitlines()
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == function_name:
            start = node.lineno - 1
            end = node.end_lineno
            return "\n".join(lines[start:end])
    raise AssertionError(f"Function {function_name!r} not found in {path}")


def _method_return_constant(path: Path, class_name: str, method_name: str) -> str:
    tree = ast.parse(_read_text(path), filename=str(path))
    for node in tree.body:
        if not isinstance(node, ast.ClassDef) or node.name != class_name:
            continue
        for child in node.body:
            if not isinstance(child, ast.FunctionDef) or child.name != method_name:
                continue
            for stmt in child.body:
                if isinstance(stmt, ast.Return) and isinstance(stmt.value, ast.Constant):
                    constant_value = stmt.value.value
                    assert isinstance(constant_value, str)
                    return constant_value
    raise AssertionError(f"{class_name}.{method_name} not found in {path}")


def test_run_shermo_calls_are_centralized() -> None:
    allowed_modules = {
        "rph_core/steps/condition_thermo.py",
        "rph_core/utils/shermo_runner.py",
        "rph_core/utils/thermo/shermo.py",
        "rph_core/utils/thermo/artifacts.py",
    }

    call_sites = _iter_run_shermo_call_sites()
    unexpected = sorted(call_sites - allowed_modules)

    assert not unexpected, (
        "run_shermo() must stay centralized in the condition/thermo layer. "
        f"Unexpected call sites: {unexpected}"
    )
    assert "rph_core/steps/conformer_search/engine.py" not in call_sites
    assert "rph_core/steps/step3_opt/ts_optimizer.py" not in call_sites
    assert "rph_core/steps/condition_thermo.py" in call_sites


def test_s4_extractors_defer_gibbs_outputs() -> None:
    thermo_path = _repo_file("rph_core", "steps", "step4_features", "extractors", "thermo.py")
    step1_path = _repo_file("rph_core", "steps", "step4_features", "extractors", "step1_activation.py")

    thermo_source = _read_text(thermo_path)
    step1_source = _read_text(step1_path)

    assert _method_return_constant(thermo_path, "ThermoExtractor", "get_plugin_name") == "thermo"
    assert "deferred_to_condition" in thermo_source
    assert re.search(r'thermo\.dG_activation"\]\s*=\s*np\.nan', thermo_source)
    assert re.search(r'thermo\.dG_reaction"\]\s*=\s*np\.nan', thermo_source)
    assert re.search(r'energy_source_activation"\]\s*=\s*"deferred_to_condition"', thermo_source)
    assert re.search(r'energy_source_reaction"\]\s*=\s*"deferred_to_condition"', thermo_source)

    assert "W_S1_THERMO_DEFERRED_TO_CONDITION" in step1_source
    assert re.search(r"['\"]s1_dG_act['\"]\s*:\s*float\(['\"]nan['\"]\)", step1_source)


def test_parse_shermo_sum_temperature_handling(tmp_path: Path) -> None:
    from rph_core.utils.thermo.shermo import parse_shermo_sum

    explicit_sum = tmp_path / "explicit.sum"
    default_sum = tmp_path / "default.sum"

    _ = explicit_sum.write_text(
        """Some header info
Temperature:  228.150 K
Pressure:     1.00000 atm
Sum of electronic energy and thermal correction to U     -102.345678
Sum of electronic energy and thermal correction to H     -102.345678
Sum of electronic energy and thermal correction to G     -102.345678
Total S:      80.011 J/mol/K      19.123 cal/mol/K    -TS:   -10.379 kcal/mol
""",
        encoding="utf-8",
    )
    _ = default_sum.write_text(
        """Some header info
Sum of electronic energy and thermal correction to U     -102.345678
Sum of electronic energy and thermal correction to H     -102.345678
Sum of electronic energy and thermal correction to G     -102.345678
Total S:      80.011 J/mol/K      19.123 cal/mol/K    -TS:   -10.379 kcal/mol
""",
        encoding="utf-8",
    )

    explicit_record = parse_shermo_sum(explicit_sum)
    default_record = parse_shermo_sum(default_sum)

    assert explicit_record.temperature_K == 228.15
    assert default_record.temperature_K == 298.15


def test_ts_optimizer_compute_shermo_gibbs_is_noop() -> None:
    ts_optimizer_path = _repo_file("rph_core", "steps", "step3_opt", "ts_optimizer.py")
    method_source = _function_source(ts_optimizer_path, "_compute_shermo_gibbs")

    assert "deferred_to_condition" in method_source
    assert "sp_report.g_ts = None" in method_source
    assert "sp_report.g_reactant = None" in method_source
    assert "run_shermo" not in method_source


def test_dr_aggregator_handles_nan_barriers(tmp_path: Path) -> None:
    dr_aggregator = importlib.import_module("rph_core.steps.dr_aggregator")
    read_barrier = cast(
        Callable[[Path, str, list[str]], float | None],
        getattr(dr_aggregator, "_read_barrier"),
    )
    compute_predicted_dr = cast(
        Callable[[list[dict[str, object]], float, list[str]], dict[str, object]],
        getattr(dr_aggregator, "_compute_predicted_dr"),
    )

    csv_file = tmp_path / "features_raw.csv"
    _ = csv_file.write_text(
        "thermo.dG_activation,thermo.dE_activation\nnan,5.5\n",
        encoding="utf-8",
    )

    warnings: list[str] = []
    barrier = read_barrier(csv_file, "test", warnings)
    assert barrier == 5.5

    result = compute_predicted_dr(
        [
            {"branch_id": "BR1", "source": ".", "delta_g_act_kcal_mol": barrier},
            {"branch_id": "BR2", "source": ".", "delta_g_act_kcal_mol": float("nan")},
        ],
        298.15,
        warnings,
    )
    assert result["status"] == "insufficient_data"


def test_feature_context_temperature_defaults_to_none() -> None:
    from rph_features.context import FeatureContext

    context_path = _repo_file("rph_core", "steps", "step4_features", "context.py")
    context_source = _read_text(context_path)

    ctx = FeatureContext()
    assert ctx.temperature_K is None
    assert "temperature_K: Optional[float] = None" in context_source
