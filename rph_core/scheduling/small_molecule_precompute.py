from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path
from typing import Any

from rph_core.scheduling.models import ReactionJob
from rph_core.utils.small_molecule_cache import SmallMoleculeCache
from rph_core.utils.small_molecule_catalog import SmallMoleculeCatalog

logger = logging.getLogger(__name__)


class SmallMoleculePrecomputer:
    def __init__(self, config: dict[str, Any], cache_root: Path) -> None:
        self.config = config
        self.cache_root = Path(cache_root)
        self.cache = SmallMoleculeCache(self.cache_root)
        self.catalog = SmallMoleculeCatalog(config)

    def collect_keys(self, reaction_jobs: list[ReactionJob]) -> set[str]:
        keys: set[str] = set()
        for job in reaction_jobs:
            meta = job.representative.meta or {}
            sm_keys_raw = meta.get("small_molecular_keys")
            if isinstance(sm_keys_raw, (list, tuple)):
                for k in sm_keys_raw:
                    k_str = str(k).strip()
                    if k_str:
                        keys.add(k_str)

            leaving_key = meta.get("leaving_small_molecule_key")
            if isinstance(leaving_key, str) and leaving_key.strip():
                keys.add(leaving_key.strip())

        ref_terms = self.config.get("reaction_reference_terms", {}) or {}
        for _profile, profile_data in ref_terms.items():
            if not isinstance(profile_data, dict):
                continue
            for step_data in profile_data.values():
                if not isinstance(step_data, dict):
                    continue
                for species_dict in (step_data.get("reactants", {}), step_data.get("products", {})):
                    if not isinstance(species_dict, dict):
                        continue
                    for key in species_dict:
                        mol = self.catalog.get(str(key))
                        if mol is not None:
                            keys.add(str(key))

        return keys

    def ensure_all(self, keys: set[str]) -> dict[str, Path]:
        results: dict[str, Path] = {}
        for key in sorted(keys):
            try:
                path = self.ensure_one(key)
                results[key] = path
            except Exception as exc:
                logger.warning(f"Small molecule precompute failed for {key}: {exc}")
        return results

    def ensure_one(self, key: str) -> Path:
        mol = self.catalog.get(key)
        if mol is None:
            raise ValueError(f"Unknown small molecule key: {key}")

        theory_signature = self._build_theory_signature()
        if self.cache.is_complete(mol.smiles, theory_signature=theory_signature, require_thermo=True):
            path = self.cache.get_path(mol.smiles)
            if path is not None:
                logger.debug(f"Small molecule cache hit: {key} ({mol.smiles})")
                return path

        lock = self.cache.acquire_compute_lock(mol.smiles)
        if lock is None:
            import time
            for _ in range(60):
                if self.cache.is_complete(mol.smiles, theory_signature=theory_signature, require_thermo=True):
                    path = self.cache.get_path(mol.smiles)
                    if path is not None:
                        logger.info(f"Small molecule cache ready (waited): {key}")
                        return path
                time.sleep(0.5)
            raise RuntimeError(f"Timeout waiting for small molecule cache: {key}")

        try:
            cache_dir = self.cache.get_or_create(mol.smiles, name=key)
            cache_dir.mkdir(parents=True, exist_ok=True)

            if not self.cache.is_complete(
                mol.smiles,
                theory_signature=theory_signature,
                require_thermo=False,
            ):
                self._compute_molecule(key, mol.smiles, cache_dir)

            self.ensure_thermo_json(cache_dir)
            return cache_dir
        finally:
            self.cache.release_compute_lock(lock)

    def _compute_molecule(self, key: str, smiles: str, cache_dir: Path) -> None:
        from rph_core.steps.anchor.engine_factory import create_s1_engine

        protocol = str((self.config.get("step1", {}) or {}).get("protocol", "lite"))
        engine = create_s1_engine(
            protocol=protocol,
            config=self.config,
            work_dir=cache_dir.parent,
            molecule_name=key,
        )
        best_sp_out, _sp_energy = engine.run_optimization_only(smiles=smiles)
        shutil.copy(best_sp_out, cache_dir / "molecule_min.xyz")

        source_dft = cache_dir.parent / key / "finalDFT"
        if not source_dft.exists():
            source_dft = cache_dir.parent / key / "dft"
        if source_dft.exists():
            dest_dft = cache_dir / "finalDFT"
            if dest_dft.exists():
                shutil.rmtree(dest_dft)
            shutil.copytree(source_dft, dest_dft)

        work_dir = cache_dir.parent / key
        if work_dir.exists() and work_dir.resolve() != cache_dir.resolve():
            shutil.rmtree(work_dir, ignore_errors=True)

        self.cache.write_cache_meta(
            smiles,
            self._build_theory_signature(),
            extra_meta={"source": "v3_small_molecule_precompute", "key": key},
        )

    def _build_theory_signature(self) -> dict[str, Any]:
        theory = self.config.get("theory", {}) or {}
        opt = theory.get("optimization", {}) or {}
        return {
            "method": opt.get("method"),
            "basis": opt.get("basis"),
            "solvent": opt.get("solvent"),
            "engine": opt.get("engine"),
        }

    def ensure_thermo_json(self, cache_entry: Path) -> Path:
        thermo_path = cache_entry / "thermo.json"
        if thermo_path.exists():
            return thermo_path

        from rph_core.utils.thermo import ensure_thermo_json_from_entry

        result = ensure_thermo_json_from_entry(cache_entry)
        if result is not None:
            return result

        meta_path = cache_entry / "cache_meta.json"
        if meta_path.exists():
            try:
                with open(meta_path, "r") as f:
                    meta = json.load(f)
                warnings = meta.get("warnings", [])
                if "W_MISSING_SHERMO_SUM_FOR_THERMO_JSON" not in warnings:
                    warnings.append("W_MISSING_SHERMO_SUM_FOR_THERMO_JSON")
                    meta["warnings"] = warnings
                    with open(meta_path, "w") as f:
                        json.dump(meta, f, indent=2)
            except Exception:
                pass

        logger.warning(f"No Shermo .sum found for {cache_entry.name}; thermo.json not created")
        return thermo_path
