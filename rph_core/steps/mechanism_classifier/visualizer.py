"""S0 mechanism visualization — single-row diastereomer reaction comparison.

Renders a single-row composite PNG with 2 or 3 molecules:
  No DR:  [Precursor (red=break)]  [Product (green=form)]
  With DR: [Precursor (red=break)]  [BR_MAJOR (green=form)]  [BR_DR_001 (green=form, red=flip)]

Bond changes and chiral-center highlights follow the same color scheme
as the original visualizer, with flipped chiral centres additionally
marked in red on the DR product.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional, Tuple

if TYPE_CHECKING:
    from rph_core.steps.mechanism_classifier.models import MechanismGraph

logger = logging.getLogger(__name__)

_BREAKING_COLOR = (1.0, 0.0, 0.0)
_FORMING_COLOR = (0.0, 0.8, 0.0)
_BREAKING_ATOM_COLOR = (1.0, 0.7, 0.7)
_FORMING_ATOM_COLOR = (0.7, 1.0, 0.7)


def _resolve_atom_idx(mol, atom_value, map_to_idx):
    atom_idx = int(atom_value)
    if 0 <= atom_idx < mol.GetNumAtoms():
        return atom_idx
    return map_to_idx.get(atom_idx)


def _annotate_map_numbers(mol):
    for atom in mol.GetAtoms():
        map_num = int(atom.GetAtomMapNum())
        if map_num > 0:
            atom.SetProp("atomNote", str(map_num))


def _build_map_to_idx(mol):
    result: dict[int, int] = {}
    for atom in mol.GetAtoms():
        map_num = int(atom.GetAtomMapNum())
        if map_num > 0:
            result[map_num] = atom.GetIdx()
    return result


def _collect_bond_highlights(
    mol, atom_pairs, map_to_idx, atom_color, bond_color
):
    atoms = set()
    bonds = []
    bcolors = {}
    acolors = {}
    for a_map, b_map in atom_pairs:
        a_idx = _resolve_atom_idx(mol, a_map, map_to_idx)
        b_idx = _resolve_atom_idx(mol, b_map, map_to_idx)
        if a_idx is None or b_idx is None:
            continue
        bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
        if bond is None:
            continue
        atoms.update((a_idx, b_idx))
        bonds.append(bond.GetIdx())
        bcolors[bond.GetIdx()] = bond_color
        acolors[a_idx] = atom_color
        acolors[b_idx] = atom_color
    return list(atoms), bonds, bcolors, acolors


class MechanismVisualizer:
    def render(
        self,
        graph: "MechanismGraph",
        output_path: Path,
        forming_bonds_map: Optional[List[Tuple[int, int]]] = None,
        breaking_bonds_map: Optional[List[Tuple[int, int]]] = None,
    ) -> Optional[Path]:
        try:
            from rdkit import Chem
            from rdkit.Chem import rdDepictor
        except ImportError:
            logger.warning("RDKit not available for mechanism visualization")
            return None

        try:
            precursor_smiles = None
            product_smiles = None
            for node in graph.nodes:
                if node.role == "reactants" and node.smiles:
                    precursor_smiles = node.smiles
                elif node.role == "product" and node.smiles:
                    product_smiles = node.smiles

            if not precursor_smiles or not product_smiles:
                logger.warning("Missing SMILES; visualization skipped")
                return None

            precursor_mol = Chem.MolFromSmiles(precursor_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)
            if precursor_mol is None or product_mol is None:
                logger.warning("Invalid SMILES; visualization skipped")
                return None

            precursor_mol = Chem.Mol(precursor_mol)
            product_mol = Chem.Mol(product_mol)
            rdDepictor.Compute2DCoords(precursor_mol)
            rdDepictor.Compute2DCoords(product_mol)
            _annotate_map_numbers(precursor_mol)
            _annotate_map_numbers(product_mol)

            precursor_map = _build_map_to_idx(precursor_mol)
            product_map = _build_map_to_idx(product_mol)

            if forming_bonds_map is None:
                forming_bonds_map = []
                for edge in graph.edges:
                    if getattr(edge, "pathway_id", "primary") != "primary":
                        continue
                    for pair in edge.forming_bonds:
                        forming_bonds_map.append((int(pair[0]), int(pair[1])))
            if breaking_bonds_map is None:
                breaking_bonds_map = []
                for edge in graph.edges:
                    if getattr(edge, "pathway_id", "primary") != "primary":
                        continue
                    for pair in edge.breaking_bonds:
                        breaking_bonds_map.append((int(pair[0]), int(pair[1])))

            forming_bonds_map = forming_bonds_map or []
            breaking_bonds_map = breaking_bonds_map or []

            precursor_atoms, precursor_bonds, precursor_bond_c, precursor_atom_c = (
                _collect_bond_highlights(
                    precursor_mol, breaking_bonds_map, precursor_map,
                    _BREAKING_ATOM_COLOR, _BREAKING_COLOR,
                )
            )
            product_atoms, product_bonds, product_bond_c, product_atom_c = (
                _collect_bond_highlights(
                    product_mol, forming_bonds_map, product_map,
                    _FORMING_ATOM_COLOR, _FORMING_COLOR,
                )
            )

            mols = [precursor_mol, product_mol]
            legends = ["Precursor (red=break)", "BR_MAJOR (green=form)"]
            hl_atoms = [precursor_atoms, product_atoms]
            hl_bonds = [precursor_bonds, product_bonds]
            hl_atom_colors = [precursor_atom_c, product_atom_c]
            hl_bond_colors = [precursor_bond_c, product_bond_c]

            dr_completion: dict = getattr(graph, "dr_completion", None) or {}
            dr_branches: list = dr_completion.get("branches", []) if isinstance(dr_completion, dict) else []

            if len(dr_branches) >= 2:
                try:
                    from rph_core.steps.mechanism_classifier.dr_smiles_renderer import (
                        _CHIRAL_FLIPPED_COLOR,
                        _CHIRAL_FIXED_COLOR,
                        _BOND_FLIPPED_COLOR,
                        _safe_parse,
                    )

                    dr_branch = dr_branches[1]
                    dr_smi = dr_branch.get("product_smiles")
                    if dr_smi:
                        dr_mol = _safe_parse(dr_smi)
                        if dr_mol is not None:
                            rdDepictor.Compute2DCoords(dr_mol)
                            _annotate_map_numbers(dr_mol)

                            flipped_set = set(dr_branch.get("flipped_map_numbers", []) or [])
                            fixed_set = set(dr_branch.get("fixed_stereocenters", []) or [])

                            dr_chiral_atoms = []
                            dr_chiral_colors = {}
                            for atom in dr_mol.GetAtoms():
                                map_num = int(atom.GetAtomMapNum())
                                idx = int(atom.GetIdx())
                                if map_num <= 0:
                                    continue
                                if map_num in flipped_set:
                                    dr_chiral_atoms.append(idx)
                                    dr_chiral_colors[idx] = _CHIRAL_FLIPPED_COLOR
                                elif map_num in fixed_set:
                                    dr_chiral_atoms.append(idx)
                                    dr_chiral_colors[idx] = _CHIRAL_FIXED_COLOR

                            dr_product_map = _build_map_to_idx(dr_mol)
                            dr_form_atoms, dr_form_bonds, dr_form_bc, dr_form_ac = (
                                _collect_bond_highlights(
                                    dr_mol, forming_bonds_map, dr_product_map,
                                    _FORMING_ATOM_COLOR, _FORMING_COLOR,
                                )
                            )

                            dr_product_atoms = list(dr_form_atoms)
                            dr_bond_ids = list(dr_form_bonds)
                            dr_bond_colors_all = dict(dr_form_bc)
                            dr_atom_colors_all = dict(dr_form_ac)

                            if dr_chiral_atoms:
                                for idx in dr_chiral_atoms:
                                    if idx not in dr_product_atoms:
                                        dr_product_atoms.append(idx)
                                dr_atom_colors_all.update(dr_chiral_colors)

                                for bond in dr_mol.GetBonds():
                                    if bond.GetBeginAtomIdx() in dr_chiral_atoms and \
                                       bond.GetEndAtomIdx() in dr_chiral_atoms:
                                        bid = bond.GetIdx()
                                        if bid not in dr_bond_ids:
                                            dr_bond_ids.append(bid)
                                        dr_bond_colors_all[bid] = _BOND_FLIPPED_COLOR

                            mols.append(dr_mol)
                            legends.append("BR_DR_001 (green=form, red=flip)")
                            hl_atoms.append(dr_product_atoms)
                            hl_bonds.append(dr_bond_ids)
                            hl_atom_colors.append(dr_atom_colors_all)
                            hl_bond_colors.append(dr_bond_colors_all)
                except Exception as exc:
                    logger.debug(f"DR molecule rendering skipped: {exc}")

            from rdkit.Chem.Draw import rdMolDraw2D

            n_mols = len(mols)
            cell_w = 600
            total_w = cell_w * n_mols
            height = 450

            drawer = rdMolDraw2D.MolDraw2DCairo(total_w, height, cell_w, height)
            drawer.DrawMolecules(
                mols, legends=legends,
                highlightAtoms=hl_atoms,
                highlightBonds=hl_bonds,
                highlightAtomColors=hl_atom_colors,
                highlightBondColors=hl_bond_colors,
            )
            drawer.FinishDrawing()

            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "wb") as f:
                f.write(drawer.GetDrawingText())

            logger.info(f"[S0] Mechanism visualization saved: {output_path}")
            return output_path

        except Exception as exc:
            logger.warning(f"[S0] Mechanism visualization failed (non-blocking): {exc}")
            return None
