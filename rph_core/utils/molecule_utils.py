"""
Molecule Utilities
==================

Utility functions for molecule processing, SMILES canonicalization, 
molecular formula generation, and small molecule detection.
"""

import logging
from typing import Optional
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

logger = logging.getLogger(__name__)

def canonicalize_smiles(smiles: str) -> Optional[str]:
    """
    Canonicalize a SMILES string using RDKit.

    Args:
        smiles: Input SMILES string.

    Returns:
        Canonical SMILES string or None if invalid.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
            return None
        mol = Chem.RemoveHs(mol)
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        logger.error(f"Error canonicalizing SMILES {smiles}: {e}")
        return None

def get_molecular_formula(smiles: str) -> Optional[str]:
    """
    Get molecular formula from a SMILES string.

    Args:
        smiles: Input SMILES string.

    Returns:
        Molecular formula string (e.g., 'C2H6O') or None if invalid.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        logger.error(f"Error getting formula for SMILES {smiles}: {e}")
        return None

def get_molecule_key(smiles: str) -> Optional[str]:
    """
    Generate a unique key for a molecule: {formula}_{canonical_smiles}.

    Args:
        smiles: Input SMILES string.

    Returns:
        Unique key string or None if invalid.
    """
    canonical = canonicalize_smiles(smiles)
    if canonical is None:
        return None
    formula = get_molecular_formula(canonical)
    if formula is None:
        return None
    return f"{formula}_{canonical}"

def is_small_molecule(smiles: str, threshold: int = 10) -> bool:
    """
    Determine if a molecule is a "small molecule" based on heavy atom count.

    Args:
        smiles: Input SMILES string.
        threshold: Heavy atom count threshold (exclusive).

    Returns:
        True if heavy atom count < threshold, False otherwise.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return mol.GetNumHeavyAtoms() < threshold
    except Exception as e:
        logger.error(f"Error checking if small molecule for SMILES {smiles}: {e}")
        return False
