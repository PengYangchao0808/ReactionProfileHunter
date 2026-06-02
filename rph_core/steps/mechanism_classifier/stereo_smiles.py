from __future__ import annotations

import logging
from typing import List, Optional

logger = logging.getLogger(__name__)


def generate_stereo_branch_smiles(
    mapped_product_smiles: str,
    flipped_map_numbers: List[int],
    fixed_map_numbers: List[int],
) -> Optional[str]:
    """Generate an explicit stereochemical branch SMILES by flipping chiral tags.

    Reads ``mapped_product_smiles``, identifies chiral centers by atom-map
    number, flips every tag for atoms in ``flipped_map_numbers``, leaves
    ``fixed_map_numbers`` untouched, and returns a canonical isomeric SMILES.

    Returns ``None`` when RDKit is unavailable, parsing fails, or no
    chiral center can be flipped. On partial failure the function degrades
    gracefully: it logs a warning and still returns the best-effort SMILES.
    """

    try:
        from rdkit import Chem
    except ImportError:
        logger.warning("RDKit unavailable; cannot generate stereo branch SMILES.")
        return None

    mol = _safe_parse(mapped_product_smiles)
    if mol is None:
        return None

    flipped_set = set(int(m) for m in flipped_map_numbers)
    fixed_set = set(int(m) for m in fixed_map_numbers)

    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    if not centers:
        return None

    map_to_idx: dict[int, int] = {}
    for atom in mol.GetAtoms():
        map_num = int(atom.GetAtomMapNum())
        if map_num > 0:
            map_to_idx[map_num] = int(atom.GetIdx())

    flipped_any = False
    for atom_idx, tag in centers:
        map_num = _map_num_of_atom(mol, atom_idx)
        if map_num is None:
            continue
        if map_num in fixed_set:
            continue
        if map_num in flipped_set:
            flipped_any = True
            _flip_chiral_tag(mol, atom_idx)

    if not flipped_any:
        logger.warning(
            "No chiral center matched flipped_map_numbers; returning original mapped SMILES."
        )

    try:
        smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception as exc:
        logger.warning(f"Failed to generate stereo branch SMILES: {exc}")
        return None

    if not smi:
        return None
    return smi


def _safe_parse(smiles: str):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning(f"Could not parse mapped product SMILES for stereo branch generation.")
        return None
    return mol


def _map_num_of_atom(mol, atom_idx: int) -> Optional[int]:
    atom = mol.GetAtomWithIdx(int(atom_idx))
    map_num = int(atom.GetAtomMapNum())
    if map_num <= 0:
        return None
    return map_num


def _flip_chiral_tag(mol, atom_idx: int) -> None:
    from rdkit import Chem
    from rdkit.Chem import ChiralType

    atom = mol.GetAtomWithIdx(int(atom_idx))
    current = atom.GetChiralTag()
    if current == ChiralType.CHI_TETRAHEDRAL_CW:
        atom.SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CCW)
    elif current == ChiralType.CHI_TETRAHEDRAL_CCW:
        atom.SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CW)
    else:
        atom.SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CW)
