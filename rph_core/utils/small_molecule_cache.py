"""
Small Molecule Cache Manager
============================

Manages a global cache directory for small molecules to avoid redundant 
conformer searches and optimizations.
"""

import json
import logging
import os
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from rph_core.utils.molecule_utils import get_molecule_key

logger = logging.getLogger(__name__)

_TOXIC_CHARS_RE = re.compile(r"[ \[\](){}]")


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

    def __init__(self, cache_root: Path, lewis_acid_enabled: bool = False):
        """
        Initialize the cache manager.

        Args:
            cache_root: Path to the root directory of the cache.
            lewis_acid_enabled: If True, LA-specific geometry cache is used
                with a `_la` suffix to differentiate from non-LA cached geometries.
        """
        self.cache_root = Path(cache_root).resolve()
        self.lewis_acid_enabled = lewis_acid_enabled
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
        if self.lewis_acid_enabled:
            key = f"{key}_la"
        return self.cache_root / self._safe_key(key)

    @staticmethod
    def _safe_key(key: str) -> str:
        return _TOXIC_CHARS_RE.sub("_", key)

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

    def find_thermo(self, smiles: str) -> Optional[Path]:
        """Find thermo.json for a cached molecule by SMILES."""
        path = self.get_path(smiles)
        if path is None:
            return None

        thermo_file = path / "thermo.json"
        if thermo_file.exists():
            return thermo_file

        return None

    def find_thermo_by_key(self, key: str, catalog: Any) -> Optional[Path]:
        """Find thermo.json for a catalog key via its mapped SMILES."""
        mol = catalog.get(key)
        if mol is None:
            return None

        return self.find_thermo(mol.smiles)

    def _is_signature_compatible(self, requested: Dict[str, Any], cached: Dict[str, Any]) -> bool:
        """Check if requested theory signature is compatible with cached signature.
        
        Args:
            requested: Requested theory parameters.
            cached: Cached theory parameters.
            
        Returns:
            True if compatible (cached params match or are superset of requested).
        """
        requested_sig = self._normalize_theory_signature(requested)
        cached_sig = self._normalize_theory_signature(cached)
        critical_params = [
            'opt_method',
            'opt_basis',
            'opt_engine',
            'sp_method',
            'sp_basis',
            'sp_engine',
            'solvent',
        ]
        for param in critical_params:
            requested_value = requested_sig.get(param)
            if requested_value is not None:
                cached_value = cached_sig.get(param)
                if cached_value is None and param.startswith("sp_"):
                    # Older cache metadata did not record final-SP theory. Let
                    # those entries remain usable; new writes include sp_*.
                    continue
                if cached_value != requested_value:
                    return False
        return True

    @staticmethod
    def _normalize_theory_signature(signature: Dict[str, Any]) -> Dict[str, Any]:
        return {
            "opt_method": signature.get("opt_method", signature.get("method")),
            "opt_basis": signature.get("opt_basis", signature.get("basis")),
            "opt_engine": signature.get("opt_engine", signature.get("engine")),
            "sp_method": signature.get("sp_method"),
            "sp_basis": signature.get("sp_basis"),
            "sp_engine": signature.get("sp_engine"),
            "solvent": signature.get("solvent"),
        }

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

        # The per-molecule directory must exist before the atomic sentinel file
        # is opened. Fresh cache entries previously skipped get_or_create() and
        # failed here with FileNotFoundError under S1 molecule-level parallelism.
        path.mkdir(parents=True, exist_ok=True)
        
        lock_file = path / ".computing"
        start_time = time.time()
        
        while time.time() - start_time < timeout:
            try:
                fd = os.open(str(lock_file), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                with os.fdopen(fd, 'w') as f:
                    f.write(f"PID: {os.getpid()}\nTime: {time.time()}\n")
                return lock_file
            except FileExistsError:
                try:
                    with open(lock_file, 'r') as f:
                        content = f.read()
                    pid_line = next((line for line in content.split('\n') if line.startswith('PID: ')), None)
                    if pid_line:
                        pid = int(pid_line.split(': ')[1])
                        process_exists = False
                        try:
                            os.kill(pid, 0)
                            process_exists = True
                        except OSError:
                            pass
                        
                        if not process_exists:
                            logger.warning(f"Found stale lock for {smiles} (PID {pid} dead), removing it.")
                            lock_file.unlink(missing_ok=True)
                            continue
                except Exception:
                    pass
                
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
        return self.find_thermo("CC(=O)O")

    def is_complete(
        self,
        smiles: str,
        theory_signature: Optional[Dict[str, Any]] = None,
        require_thermo: bool = True,
    ) -> bool:
        """Check if a cache entry is fully complete.

        Unlike exists(), this requires all expected artifacts to be present.
        With require_thermo=True, both molecule_min.xyz AND thermo.json must exist.

        Args:
            smiles: SMILES string of the molecule.
            theory_signature: Optional theory parameters to validate cache compatibility.
            require_thermo: If True, also require thermo.json to exist.

        Returns:
            True if cache entry is complete and compatible, False otherwise.
        """
        path = self.get_path(smiles)
        if path is None or not path.exists():
            return False

        if not (path / "molecule_min.xyz").exists():
            return False

        if require_thermo and not (path / "thermo.json").exists():
            return False

        if theory_signature is not None:
            meta_file = path / "cache_meta.json"
            if meta_file.exists():
                try:
                    with open(meta_file, "r") as f:
                        cache_meta = json.load(f)
                    cached_sig = cache_meta.get('theory_signature', {})
                    if not self._is_signature_compatible(theory_signature, cached_sig):
                        return False
                except Exception:
                    return False
            else:
                return False

        return True

    def get_entry_manifest(self, smiles: str) -> Dict[str, Any]:
        """Read the cache_meta.json for a cached molecule.

        Args:
            smiles: SMILES string of the molecule.

        Returns:
            Dict of cache metadata, or empty dict if not found.
        """
        path = self.get_path(smiles)
        if path is None:
            return {}
        meta_file = path / "cache_meta.json"
        if not meta_file.exists():
            return {}
        try:
            with open(meta_file, "r") as f:
                data = json.load(f)
            return data if isinstance(data, dict) else {}
        except Exception:
            return {}


class CalibrationCache:
    """Cache manager for Lewis acid calibration artifacts."""

    def __init__(self, cache_root: Union[Path, str]) -> None:
        """Initialize with the Lewis acid calibration cache root."""
        self.cache_root = Path(cache_root).resolve()
        self.logger = logging.getLogger(__name__)

    @staticmethod
    def _safe_component(name: str) -> str:
        """Normalize path components for calibration cache directories."""
        return _TOXIC_CHARS_RE.sub("_", name)

    def get_free_molecule_dir(self, molecule_name: str) -> Path:
        """Get cache directory for a free molecule."""
        return self.cache_root / self._safe_component(molecule_name) / "free"

    def get_complex_dir(self, la_name: str) -> Path:
        """Get cache directory for an acetone-LA complex."""
        return self.cache_root / self._safe_component(la_name) / "acetone_complex"

    def free_molecule_exists(self, molecule_name: str) -> bool:
        """Check if free molecule cache exists."""
        return (self.get_free_molecule_dir(molecule_name) / "molecule_min.xyz").exists()

    def complex_exists(self, la_name: str) -> bool:
        """Check if acetone-LA complex cache exists."""
        return (self.get_complex_dir(la_name) / "best_complex.xyz").exists()

    def get_free_molecule_xyz(self, molecule_name: str) -> Optional[Path]:
        """Get cached free-molecule geometry path if present."""
        xyz_path = self.get_free_molecule_dir(molecule_name) / "molecule_min.xyz"
        return xyz_path if xyz_path.exists() else None

    def get_complex_xyz(self, la_name: str) -> Optional[Path]:
        """Get cached acetone-LA complex geometry path if present."""
        xyz_path = self.get_complex_dir(la_name) / "best_complex.xyz"
        return xyz_path if xyz_path.exists() else None

    def write_cache_meta(self, dir_path: Path, signature: Dict[str, Any]) -> None:
        """Write cache metadata with theory signature."""
        dir_path.mkdir(parents=True, exist_ok=True)
        meta_path = dir_path / "cache_meta.json"
        payload = {
            "theory_signature": signature,
            "created_at": datetime.now(timezone.utc).isoformat(),
        }
        meta_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    def get_cache_meta(self, dir_path: Path) -> Optional[Dict[str, Any]]:
        """Read cache_meta.json if present and valid."""
        meta_path = dir_path / "cache_meta.json"
        if not meta_path.exists():
            return None
        try:
            data = json.loads(meta_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as exc:
            self.logger.warning(f"Failed to read calibration cache meta {meta_path}: {exc}")
            return None
        return data if isinstance(data, dict) else None

    def compute_cache_signature(
        self,
        la_name: str,
        surrogate_smiles: str,
        mode: str,
        charge: int,
        multiplicity: int,
        opt_method: str,
        opt_basis: str,
        sp_method: str,
        sp_basis: str,
        solvent: str,
        charge_priority: List[str],
        calibration_schema_version: str = "la_calibration_v1",
        complex_builder_version: str = "v1",
    ) -> Dict[str, Any]:
        """Compute the cache signature for Lewis acid calibration artifacts."""
        return {
            "lewis_acid": la_name,
            "surrogate_smiles": surrogate_smiles,
            "mode": mode,
            "charge": charge,
            "multiplicity": multiplicity,
            "opt_method": opt_method,
            "opt_basis": opt_basis,
            "sp_method": sp_method,
            "sp_basis": sp_basis,
            "solvent": solvent,
            "charge_priority": list(charge_priority),
            "calibration_schema_version": calibration_schema_version,
            "complex_builder_version": complex_builder_version,
        }

    def is_signature_compatible(
        self, stored: Dict[str, Any], expected: Dict[str, Any]
    ) -> bool:
        """Check whether stored and expected cache signatures are compatible."""
        if not stored or not expected:
            return False
        for key in expected.keys():
            if stored.get(key) != expected.get(key):
                return False
        return True
