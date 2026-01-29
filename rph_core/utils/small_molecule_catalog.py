"""
Small Molecule Catalog for Reference States

Provides mapping from leaving group labels (e.g., "AcOH", "TFE")
to SMILES and quantum chemistry parameters (charge, multiplicity).
"""

from dataclasses import dataclass
from typing import Optional, Dict, List


@dataclass
class SmallMolecule:
    """Small molecule definition for reference state calculation."""
    key: str
    smiles: str
    charge: int = 0
    multiplicity: int = 1


class UnknownSmallMoleculeError(Exception):
    """Raised when a small molecule key is not found in catalog."""
    pass


class SmallMoleculeCatalog:
    """
    Catalog of small molecules with configuration-driven mapping.

    Loaded from config['reference_states']['small_molecule_map'].
    """

    def __init__(self, config: Dict):
        """
        Initialize catalog from configuration.

        Args:
            config: Configuration dictionary containing
                     reference_states.small_molecule_map section
        """
        self._molecules: Dict[str, SmallMolecule] = {}

        ref_config = config.get('reference_states', {})
        map_config = ref_config.get('small_molecule_map', {})

        for key, mol_dict in map_config.items():
            smiles = mol_dict.get('smiles', '').strip()
            charge = int(mol_dict.get('charge', 0))
            multiplicity = int(mol_dict.get('multiplicity', 1))

            if not smiles:
                continue  # Skip entries without SMILES

            self._molecules[key] = SmallMolecule(
                key=key,
                smiles=smiles,
                charge=charge,
                multiplicity=multiplicity
            )

    def get(self, key: str) -> Optional[SmallMolecule]:
        """
        Get small molecule by key.

        Args:
            key: Molecule key (e.g., "AcOH", "TFE")

        Returns:
            SmallMolecule object or None if not found
        """
        return self._molecules.get(key)

    def require(self, key: str) -> SmallMolecule:
        """
        Get small molecule by key, raising error if not found.

        Args:
            key: Molecule key

        Returns:
            SmallMolecule object

        Raises:
            UnknownSmallMoleculeError: If key not in catalog
        """
        mol = self.get(key)
        if mol is None:
            raise UnknownSmallMoleculeError(f"Unknown small molecule key: {key}")
        return mol

    def validate_keys(self, keys: List[str]) -> List[str]:
        """
        Validate list of keys, returning unknown keys.

        Args:
            keys: List of molecule keys to validate

        Returns:
            List of unknown keys (empty if all valid)
        """
        unknown = [k for k in keys if k not in self._molecules]
        return unknown

    def list_keys(self) -> List[str]:
        """
        Get list of all available molecule keys.

        Returns:
            List of keys
        """
        return list(self._molecules.keys())

    def __len__(self) -> int:
        return len(self._molecules)
