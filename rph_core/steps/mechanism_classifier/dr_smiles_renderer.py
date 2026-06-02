r"""S0 DR SMILES rendering utilities for chiral-flip verification.

This module provides low-level RDKit rendering primitives used by
``MechanismVisualizer`` (visualizer.py) to build the DR-comparison row
of the composite mechanism_graph.png.  The only public API is
``render_dr_comparison_images`` (for ad-hoc scripting); the color
constants and ``_safe_parse`` are imported by the visualizer directly.

Typical direct usage (scripting only)::

    render_dr_comparison_images(dr_plan, output_dir)

The orchestrator no longer calls this function — the visualizer handles
DR rendering internally via ``graph.dr_completion``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, List, Optional

logger = logging.getLogger(__name__)

_CHIRAL_FLIPPED_COLOR = (0.9, 0.3, 0.3)
_CHIRAL_FIXED_COLOR = (0.3, 0.5, 0.9)
_BOND_FLIPPED_COLOR = (0.9, 0.5, 0.5)
_LEGEND_FONT_SIZE = 18
_ATOM_LABEL_FONT_SIZE = 14


def render_dr_comparison_images(
    dr_plan: Any,
    output_dir: Path,
    *,
    prefix: str = "smiles",
) -> bool:
    """Render DR comparison images from a DRBranchPlan.

    Attempts to render a side-by-side comparison PNG and individual per-branch
    PNGs. The function is designed to be called from the orchestrator's S0
    step, right after attach_dr_completion() populates the branch plan.

    The ``dr_plan`` argument must be a DRBranchPlan-like object with a
    ``branches`` attribute containing at least one DRProductBranch.

    Returns True if any PNG was written, False otherwise.
    """
    branches: List[Any] = getattr(dr_plan, "branches", []) or []
    if len(branches) < 2:
        return False

    output_dir.mkdir(parents=True, exist_ok=True)
    produced = False

    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        logger.warning("RDKit unavailable; SMILES comparison rendering skipped.")
        return False

    mols: List[Chem.Mol] = []
    legends: List[str] = []
    highlight_atoms: List[List[int]] = []
    highlight_colors: List[dict[int, tuple[float, float, float]]] = []

    for branch in branches:
        smi: Optional[str] = getattr(branch, "product_smiles", None)
        if not smi:
            continue

        mol = _safe_parse(smi)
        if mol is None:
            logger.warning(
                f"[S0] Cannot render {getattr(branch, 'branch_id', '?')}: "
                f"SMILES unparseable ({smi[:60]})"
            )
            continue

        rdDepictor.Compute2DCoords(mol)
        mols.append(mol)

        bid: str = getattr(branch, "branch_id", "?")
        role: str = getattr(branch, "role", "")
        legends.append(f"{bid} ({role})")

        flipped_set: set[int] = set(
            int(v) for v in (getattr(branch, "flipped_map_numbers", []) or [])
        )
        fixed_set: set[int] = set(
            int(v) for v in (getattr(branch, "fixed_stereocenters", []) or [])
        )

        atoms: List[int] = []
        colors: dict[int, tuple[float, float, float]] = {}
        for atom in mol.GetAtoms():
            map_num = int(atom.GetAtomMapNum())
            idx = int(atom.GetIdx())
            if map_num <= 0:
                continue
            if map_num in flipped_set:
                atoms.append(idx)
                colors[idx] = _CHIRAL_FLIPPED_COLOR
            elif map_num in fixed_set:
                atoms.append(idx)
                colors[idx] = _CHIRAL_FIXED_COLOR

        highlight_atoms.append(atoms)
        highlight_colors.append(colors)

    if not mols:
        return False

    for idx, mol in enumerate(mols):
        bid_str: str = legends[idx].split(" (")[0] if idx < len(legends) else "branch"
        try:
            img_path = output_dir / f"{prefix}_{bid_str}.png"
            _render_single_mol(
                mol, img_path,
                highlight_atoms[idx] if idx < len(highlight_atoms) else [],
                highlight_colors[idx] if idx < len(highlight_colors) else {},
                legend=legends[idx] if idx < len(legends) else "",
            )
            produced = True
        except Exception as exc:
            logger.warning(f"Failed to render {bid_str}: {exc}")

    if len(mols) >= 2:
        try:
            comparison_path = output_dir / f"{prefix}_comparison.png"
            _render_comparison(
                mols, legends,
                highlight_atoms, highlight_colors,
                comparison_path,
            )
            produced = True
        except Exception as exc:
            logger.warning(f"Failed to render comparison: {exc}")

    return produced


def _safe_parse(smiles: str):
    from rdkit import Chem

    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return None
    return mol


def _render_single_mol(
    mol,
    output_path: Path,
    highlight_atoms: List[int],
    atom_colors: dict[int, tuple[float, float, float]],
    *,
    legend: str = "",
    size: tuple[int, int] = (600, 500),
) -> None:
    from rdkit.Chem.Draw import rdMolDraw2D

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.addAtomIndices = False

    highlight_bonds: List[int] = []
    bond_colors: dict[int, tuple[float, float, float]] = {}
    if highlight_atoms:
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in highlight_atoms and bond.GetEndAtomIdx() in highlight_atoms:
                highlight_bonds.append(bond.GetIdx())
                bond_colors[bond.GetIdx()] = _BOND_FLIPPED_COLOR

    drawer.DrawMolecule(
        mol,
        legend=legend,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=atom_colors,
        highlightBonds=highlight_bonds,
        highlightBondColors=bond_colors,
    )
    drawer.FinishDrawing()

    with open(output_path, "wb") as f:
        f.write(drawer.GetDrawingText())

    logger.info(f"[S0] SMILES image saved: {output_path}")


def _render_comparison(
    mols,
    legends: List[str],
    highlight_atoms: List[List[int]],
    highlight_colors: List[dict[int, tuple[float, float, float]]],
    output_path: Path,
) -> None:
    from rdkit.Chem.Draw import rdMolDraw2D

    per_width = 500
    total_width = per_width * len(mols)
    height = 450

    drawer = rdMolDraw2D.MolDraw2DCairo(total_width, height, per_width, height)
    opts = drawer.drawOptions()

    highlight_bonds_lists: List[List[int]] = []
    bond_colors_lists: List[dict[int, tuple[float, float, float]]] = []
    for mol, atoms in zip(mols, highlight_atoms):
        bonds: List[int] = []
        bcolors: dict[int, tuple[float, float, float]] = {}
        if atoms:
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in atoms and bond.GetEndAtomIdx() in atoms:
                    bonds.append(bond.GetIdx())
                    bcolors[bond.GetIdx()] = _BOND_FLIPPED_COLOR
        highlight_bonds_lists.append(bonds)
        bond_colors_lists.append(bcolors)

    drawer.DrawMolecules(
        mols,
        legends=legends,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_colors,
        highlightBonds=highlight_bonds_lists,
        highlightBondColors=bond_colors_lists,
    )
    drawer.FinishDrawing()

    with open(output_path, "wb") as f:
        f.write(drawer.GetDrawingText())

    logger.info(f"[S0] SMILES comparison saved: {output_path}")
