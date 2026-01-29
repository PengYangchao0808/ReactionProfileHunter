"""
Multiwfn Recipe Loader and Renderer (V4.2 Phase B/C)
================================================

Load Multiwfn recipe specifications from YAML and render menu sequences.
Phase A/B: Interface only (dry-run, no execution).
Phase C: Actual Multiwfn execution to be added.

Author: QC Descriptors Team
Date: 2026-01-18
"""

import yaml
from pathlib import Path
from typing import Dict, List, Any, Optional
from dataclasses import dataclass, field


@dataclass
class MultiwfnRecipe:
    """Multiwfn recipe specification."""
    description: str
    menu_sequence: List[str]
    outputs: List[str]
    version: str = "1.0"


def load_multiwfn_recipes(path: Path) -> Dict[str, MultiwfnRecipe]:
    """Load Multiwfn recipes from YAML file.

    Args:
        path: Path to multiwfn_recipes.yaml

    Returns:
        Dictionary mapping recipe name -> MultiwfnRecipe
    """
    if not path.exists():
        raise FileNotFoundError(f"Multiwfn recipes file not found: {path}")

    with open(path, "r") as f:
        data = yaml.safe_load(f)

    if not data:
        return {}

    recipes = {}
    for name, spec in data.items():
        if not isinstance(spec, dict):
            continue

        description = spec.get("description", "")
        menu_sequence = spec.get("menu_sequence", [])
        outputs = spec.get("outputs", [])
        version = spec.get("version", "1.0")

        if not menu_sequence or not outputs:
            continue

        recipes[name] = MultiwfnRecipe(
            description=description,
            menu_sequence=menu_sequence,
            outputs=outputs,
            version=version,
        )

    return recipes


def render_recipe(
    recipe: MultiwfnRecipe,
    *,
    input_fchk: Path,
    workdir: Path,
    output_prefix: str = "output",
) -> List[str]:
    """Render Multiwfn recipe into menu sequence (dry-run).

    Args:
        recipe: MultiwfnRecipe to render
        input_fchk: Path to Gaussian .fchk file
        workdir: Working directory
        output_prefix: Prefix for output files

    Returns:
        List of rendered menu choices (lines to send to Multiwfn)
    """
    input_str = str(input_fchk)
    workdir_str = str(workdir)

    rendered = []
    for item in recipe.menu_sequence:
        line = item.replace("{input_fchk}", input_str)
        line = line.replace("{workdir}", workdir_str)
        line = line.replace("{output_prefix}", output_prefix)
        rendered.append(line)

    return rendered
