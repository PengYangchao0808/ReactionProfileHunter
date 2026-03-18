"""
Small Molecule Cache Manager
============================

Manages a global cache directory for small molecules to avoid redundant 
conformer searches and optimizations.
"""

import fcntl
import json
import logging
import os
import time
from pathlib import Path
from typing import Any, Dict, Optional

from rph_core.utils.molecule_utils import get_molecule_key

logger = logging.getLogger(__name__)


class SmallMoleculeCache:
    """
    Manages a global cache for small molecule conformers and geometries.
    
    Structure:
    cache_root/
        {molecule_key}/
            molecule_min.xyz
            dft/
            cache_meta.json
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

    def exists(self, smiles: str, theory_signature: Optional[Dict[str, Any]] = None) -> bool:
        """
        Check if a valid cache entry exists for the given SMILES.
        A valid entry must contain 'molecule_min.xyz'.
        If theory_signature is provided, also checks cache signature compatibility.

        Args:
            smiles: SMILES string of the molecule.
            theory_signature: Optional dictionary with theory parameters to validate cache compatibility.

        Returns:
            True if valid cache exists and is compatible, False otherwise.
        """
        path = self.get_path(smiles)
        if path is None or not path.exists():
            return False
        
        min_xyz = path / "molecule_min.xyz"
        if not min_xyz.exists():
            return False
        
        if theory_signature is not None:
            meta_file = path / "cache_meta.json"
            if meta_file.exists():
                try:
                    with open(meta_file, 'r') as f:
                        cache_meta = json.load(f)
                    cached_sig = cache_meta.get('theory_signature', {})
                    if not self._is_signature_compatible(theory_signature, cached_sig):
                        logger.warning(f"Cache signature mismatch for {smiles}, will recompute")
                        return False
                except Exception as e:
                    logger.warning(f"Failed to read cache meta for {smiles}: {e}")
                    return False
            else:
                logger.debug(f"No cache meta for {smiles}, skipping signature check")
        
        return True

    def _is_signature_compatible(self, requested: Dict[str, Any], cached: Dict[str, Any]) -> bool:
        """Check if requested theory signature is compatible with cached signature.
        
        Args:
            requested: Requested theory parameters.
            cached: Cached theory parameters.
            
        Returns:
            True if compatible (cached params match or are superset of requested).
        """
        critical_params = ['method', 'basis', 'solvent', 'engine']
        for param in critical_params:
            if param in requested:
                if param not in cached or cached[param] != requested[param]:
                    return False
        return True

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

    def acquire_compute_lock(self, smiles: str, timeout: float = 300.0) -> Optional[Path]:
        """Acquire a lock for computing a small molecule to prevent cache stampede.
        
        Uses file-based locking with fcntl (Unix) or a simple sentinel file.
        
        Args:
            smiles: SMILES string of the molecule.
            timeout: Maximum time to wait for lock in seconds.
            
        Returns:
            Path to the lock file if acquired, None if lock could not be acquired.
        """
        path = self.get_path(smiles)
        if path is None:
            return None
        
        lock_file = path / ".computing"
        start_time = time.time()
        
        while time.time() - start_time < timeout:
            try:
                fd = os.open(str(lock_file), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                with os.fdopen(fd, 'w') as f:
                    f.write(f"PID: {os.getpid()}\nTime: {time.time()}\n")
                return lock_file
            except FileExistsError:
                time.sleep(0.5)
                continue
        
        logger.warning(f"Could not acquire compute lock for {smiles} within {timeout}s")
        return None

    def release_compute_lock(self, lock_file: Path) -> None:
        """Release the compute lock.
        
        Args:
            lock_file: Path to the lock file returned by acquire_compute_lock.
        """
        try:
            if lock_file.exists():
                lock_file.unlink()
        except Exception as e:
            logger.warning(f"Failed to release lock {lock_file}: {e}")

    def write_cache_meta(self, smiles: str, theory_signature: Dict[str, Any], 
                        extra_meta: Optional[Dict[str, Any]] = None) -> None:
        """Write cache metadata including theory signature.
        
        Args:
            smiles: SMILES string of the molecule.
            theory_signature: Dictionary with theory parameters used for computation.
            extra_meta: Optional additional metadata to store.
        """
        path = self.get_path(smiles)
        if path is None:
            return
        
        meta_file = path / "cache_meta.json"
        meta_data = {
            'theory_signature': theory_signature,
            'timestamp': time.time(),
        }
        if extra_meta:
            meta_data.update(extra_meta)
        
        try:
            with open(meta_file, 'w') as f:
                json.dump(meta_data, f, indent=2)
        except Exception as e:
            logger.warning(f"Failed to write cache meta for {smiles}: {e}")

    def find_hoac_thermo(self) -> Optional[Path]:
        """Find HOAc thermo.json in the cache.
        
        Searches for AcOH/HOAc in the cache and returns the thermo.json path if found.
        
        Returns:
            Path to thermo.json if found, None otherwise.
        """
        hoac_smiles = "CC(=O)O"
        path = self.get_path(hoac_smiles)
        if path is None:
            return None
        
        thermo_file = path / "thermo.json"
        if thermo_file.exists():
            return thermo_file
        
        return None
