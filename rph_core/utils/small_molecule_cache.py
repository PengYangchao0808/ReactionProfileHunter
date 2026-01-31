"""
Small Molecule Cache Manager
============================

Manages a global cache directory for small molecules to avoid redundant 
conformer searches and optimizations.
"""

import logging
from pathlib import Path
from typing import Optional

from rph_core.utils.molecule_utils import get_molecule_key

logger = logging.getLogger(__name__)

class SmallMoleculeCache:
    """
    Manages a global cache for small molecule conformers and geometries.
    
    Structure:
    cache_root/
        {molecule_key}/
            molecule_min.xyz
            ...
    """

    def __init__(self, cache_root: Path):
        """
        Initialize the cache manager.

        Args:
            cache_root: Path to the root directory of the cache.
        """
        self.cache_root = Path(cache_root).resolve()
        self.cache_root.mkdir(parents=True, exist_ok=True)
        logger.info(f"SmallMoleculeCache initialized at: {self.cache_root}")

    def get_path(self, smiles: str) -> Optional[Path]:
        """
        Get the path to the cache directory for a given SMILES.

        Args:
            smiles: SMILES string of the molecule.

        Returns:
            Path to the cache directory or None if SMILES is invalid.
        """
        key = get_molecule_key(smiles)
        if key is None:
            return None
        return self.cache_root / key

    def exists(self, smiles: str) -> bool:
        """
        Check if a valid cache entry exists for the given SMILES.
        A valid entry must contain 'molecule_min.xyz'.

        Args:
            smiles: SMILES string of the molecule.

        Returns:
            True if valid cache exists, False otherwise.
        """
        path = self.get_path(smiles)
        if path is None or not path.exists():
            return False
        
        min_xyz = path / "molecule_min.xyz"
        return min_xyz.exists()

    def get_or_create(self, smiles: str, name: str = "") -> Path:
        """
        Find or create a cache directory for the given SMILES.

        Args:
            smiles: SMILES string of the molecule.
            name: Optional name for logging.

        Returns:
            Path to the cache directory.

        Raises:
            ValueError: If SMILES is invalid and key cannot be generated.
        """
        path = self.get_path(smiles)
        if path is None:
            raise ValueError(f"Invalid SMILES provided to cache: {smiles}")
        
        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created new cache directory for {name or smiles}: {path}")
        else:
            logger.debug(f"Found existing cache directory for {name or smiles}: {path}")
        
        return path
