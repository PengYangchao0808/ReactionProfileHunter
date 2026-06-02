from __future__ import annotations
import logging
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol, TypedDict, cast

logger = logging.getLogger(__name__)


class TaskLike(Protocol):
    rx_id: str
    product_smiles: str
    meta: object


class TaskDict(TypedDict):
    rx_id: str
    product_smiles: str
    precursor_smiles: str | None
    small_molecular_keys: list[str]
    reaction_profile: str | None
    cleaner_data: dict[str, object] | None


class ReactionResultDict(TypedDict):
    success: bool
    rx_id: str
    error: str | None


ReactionGroupResults = dict[str, list[ReactionResultDict]]


@dataclass
class ResourcePlan:
    max_reaction_workers: int = 4
    qc_nproc: int = 8
    orca_maxcore: int = 3000
    scratch_root: Path | None = None


def _run_one_reaction_group(
    config_path: str,
    log_level: str,
    reaction_id: str,
    task_dicts: list[TaskDict],
    output_root: str,
    skip_steps: list[str],
    resume: bool,
) -> list[ReactionResultDict]:
    import sys
    from pathlib import Path as _Path

    if str(_Path(__file__).parent.parent) not in sys.path:
        sys.path.insert(0, str(_Path(__file__).parent.parent))
    from rph_core.orchestrator import ReactionProfileHunter
    from rph_core.utils.path_manager import get_reaction_root

    _ = resume
    hunter = ReactionProfileHunter(config_path=_Path(config_path) if config_path else None, log_level=log_level)
    reaction_root = get_reaction_root(_Path(output_root), reaction_id)
    results: list[ReactionResultDict] = []

    for td in task_dicts:
        try:
            result = hunter.run_pipeline(
                product_smiles=td["product_smiles"],
                work_dir=reaction_root,
                skip_steps=list(skip_steps) if skip_steps else [],
                precursor_smiles=td.get("precursor_smiles"),
                small_molecular_keys=td.get("small_molecular_keys", []),
                reaction_profile=td.get("reaction_profile"),
                cleaner_data=td.get("cleaner_data"),
            )
            results.append({"success": result.success, "rx_id": td["rx_id"], "error": result.error_message})
        except Exception as exc:
            logger.error(f"Reaction group {reaction_id} task {td['rx_id']} failed: {exc}")
            results.append({"success": False, "rx_id": td["rx_id"], "error": str(exc)})
    return results


@dataclass
class DatasetScheduler:
    config: dict[str, object]
    resource_plan: ResourcePlan = field(default_factory=ResourcePlan)

    def run_reaction_groups(
        self,
        reaction_groups: dict[str, list[TaskLike]],
        output_root: Path,
        skip_steps: list[str] | None = None,
        resume: bool = True,
    ) -> ReactionGroupResults:
        if self.resource_plan.max_reaction_workers <= 1:
            results: ReactionGroupResults = {}
            for reaction_id, tasks in reaction_groups.items():
                td_list = [self._task_to_dict(t) for t in tasks]
                results[reaction_id] = _run_one_reaction_group(
                    config_path=self._config_path(),
                    log_level=self._log_level(),
                    reaction_id=reaction_id,
                    task_dicts=td_list,
                    output_root=str(output_root),
                    skip_steps=list(skip_steps) if skip_steps else [],
                    resume=resume,
                )
            return results

        futures_map: dict[Future[list[ReactionResultDict]], str] = {}
        with ProcessPoolExecutor(max_workers=self.resource_plan.max_reaction_workers) as executor:
            for reaction_id, tasks in reaction_groups.items():
                td_list = [self._task_to_dict(t) for t in tasks]
                future = executor.submit(
                    _run_one_reaction_group,
                    config_path=self._config_path(),
                    log_level=self._log_level(),
                    reaction_id=reaction_id,
                    task_dicts=td_list,
                    output_root=str(output_root),
                    skip_steps=list(skip_steps) if skip_steps else [],
                    resume=resume,
                )
                futures_map[future] = reaction_id

            all_results: ReactionGroupResults = {}
            for future in as_completed(futures_map):
                reaction_id = futures_map[future]
                try:
                    all_results[reaction_id] = future.result()
                except Exception as exc:
                    logger.error(f"Reaction group {reaction_id} crashed: {exc}")
                    all_results[reaction_id] = [{"success": False, "rx_id": reaction_id, "error": str(exc)}]
            return all_results

    def prewarm_reference_cache(self, hunter: object) -> None:
        from rph_core.utils.small_molecule_catalog import SmallMoleculeCatalog
        from rph_core.utils.small_molecule_cache import SmallMoleculeCache

        del hunter
        ref_terms = self._config_dict("reaction_reference_terms")
        catalog = SmallMoleculeCatalog(self.config)
        global_cfg = self._config_dict("global")
        cache_dir_obj = global_cfg.get("small_molecule_cache_dir")
        cache_dir = cache_dir_obj if isinstance(cache_dir_obj, str) else None
        cache = SmallMoleculeCache(Path(cache_dir) if cache_dir else Path("./rph_output/small_molecules"))

        species_keys: set[str] = set()
        for profile_cfg in ref_terms.values():
            profile_cfg_dict = self._as_object_dict(profile_cfg)
            pti = self._as_object_dict(profile_cfg_dict.get("precursor_to_intermediate", {}))
            for side in ("reactants", "products"):
                side_dict = self._as_object_dict(pti.get(side, {}))
                species_keys.update(side_dict.keys())

        for key in species_keys:
            mol = catalog.get(key)
            if mol is None:
                logger.warning(f"Reference species '{key}' not in catalog, skip prewarm")
                continue
            existing = cache.find_thermo(mol.smiles)
            if existing:
                logger.info(f"Reference cache hit: {key} -> {existing}")
            else:
                logger.info(f"Reference cache miss: {key} ({mol.smiles}) — will compute on first use")

    @staticmethod
    def _task_to_dict(task: TaskLike) -> TaskDict:
        meta_dict = DatasetScheduler._as_object_dict(getattr(task, "meta", {}) or {})
        precursor_smiles_obj = meta_dict.get("precursor_smiles")
        small_molecular_keys_obj = meta_dict.get("small_molecular_keys")
        reaction_profile_obj = meta_dict.get("reaction_profile")
        cleaner_data_obj = meta_dict.get("cleaner_data")
        return {
            "rx_id": getattr(task, "rx_id", ""),
            "product_smiles": getattr(task, "product_smiles", ""),
            "precursor_smiles": precursor_smiles_obj if isinstance(precursor_smiles_obj, str) else None,
            "small_molecular_keys": small_molecular_keys_obj if isinstance(small_molecular_keys_obj, list) else [],
            "reaction_profile": reaction_profile_obj if isinstance(reaction_profile_obj, str) else None,
            "cleaner_data": cleaner_data_obj if isinstance(cleaner_data_obj, dict) else None,
        }

    def _config_dict(self, key: str) -> dict[str, object]:
        value = self.config.get(key)
        return self._as_object_dict(value)

    def _config_path(self) -> str:
        value = self.config.get("_config_path")
        return value if isinstance(value, str) else ""

    def _log_level(self) -> str:
        value = self._config_dict("global").get("log_level")
        return value if isinstance(value, str) else "INFO"

    @staticmethod
    def _as_object_dict(value: object) -> dict[str, object]:
        if not isinstance(value, dict):
            return {}
        typed_value = cast(dict[object, object], value)
        result: dict[str, object] = {}
        for key, item in typed_value.items():
            result[str(key)] = item
        return result
