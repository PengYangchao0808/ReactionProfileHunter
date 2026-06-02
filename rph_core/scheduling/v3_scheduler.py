from __future__ import annotations

import copy
import json
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Optional

from rph_core.scheduling.models import BranchJob, ConditionJob, ReactionJob
from rph_core.utils.json_io import read_json_dict, write_json
from rph_core.utils.path_compat import normalize_path
from rph_core.utils.path_manager import get_branch_root, get_reaction_root
from rph_core.utils.task_builder import TaskSpec, build_tasks_from_run_config, sanitize_rx_id
from rph_core.utils.v3_progress import (
    BranchProgress,
    FinalSummary,
    PlanningProgress,
    ReactionExecutionProgress,
    ReactionPlanRecord,
    SmallMoleculeCacheProgress,
    SmallMoleculeRecord,
    V3ItemState,
    V3RunProgress,
    V3StageState,
)
from rph_core.utils.v3_stage_display import V3StageDisplay

logger = logging.getLogger(__name__)
_SAFE_SEGMENT_RE = re.compile(r"^[A-Za-z0-9_-]+$")


def _ensure_dir(path: Path) -> None:
    if path.is_file():
        path.unlink()
    path.mkdir(parents=True, exist_ok=True)


def _safe_segment(value: Any, label: str) -> str:
    text = str(value or "").strip()
    if not text or not _SAFE_SEGMENT_RE.match(text):
        raise ValueError(f"Unsafe {label}: {value!r}")
    return text


def _load_dr_branch_plan(reaction_root: Path) -> Optional[dict[str, Any]]:
    plan_path = reaction_root / "S0_Mechanism" / "dr_branch_plan.json"
    if not plan_path.exists():
        return None
    try:
        with open(plan_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def _write_branch_manifest(
    branch_root: Path,
    branch_id: str,
    pathway_id: str,
    product_smiles: str,
    generation_policy: str,
    flipped_map_numbers: list[int],
    fixed_stereocenters: list[int],
    notes: list[str],
    dr_plan_path: str,
    parent_product_smiles: str,
) -> None:
    _ensure_dir(branch_root)
    manifest = {
        "version": "rph_v3.0",
        "branch_id": branch_id,
        "pathway_id": pathway_id,
        "product_smiles": product_smiles,
        "parent_product_smiles": parent_product_smiles,
        "generation_policy": generation_policy,
        "flipped_map_numbers": flipped_map_numbers,
        "fixed_stereocenters": fixed_stereocenters,
        "notes": notes,
        "dr_plan_path": dr_plan_path,
        "status": "PENDING",
    }
    write_json(branch_root / "branch_manifest.json", manifest)


def _to_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _theory_signature_string(config: dict[str, Any]) -> str:
    theory_opt = (config.get("theory", {}) or {}).get("optimization", {}) or {}
    theory_sp = (config.get("theory", {}) or {}).get("single_point", {}) or {}
    opt_method = str(theory_opt.get("method") or "").strip()
    opt_basis = str(theory_opt.get("basis") or "").strip()
    sp_method = str(theory_sp.get("method") or "").strip()
    sp_basis = str(theory_sp.get("basis") or "").strip()
    solvent_cfg = (config.get("solvent", {}) or {})
    solvent_name = str(solvent_cfg.get("name") or "").strip().lower()
    if not solvent_name:
        solvent_name = str(theory_opt.get("solvent") or theory_sp.get("solvent") or "").strip().lower()
    solvent_token = solvent_name or "dcm"
    tokens = [t for t in [opt_method, opt_basis, sp_method, sp_basis] if t]
    base = "_".join(tokens) if tokens else ""
    if base:
        return f"{base}_SMD_{solvent_token.upper()}"
    return f"SMD_{solvent_token.upper()}"


def _make_pipeline_failure_result(
    *,
    product_smiles: str,
    work_dir: Path,
    error_step: str,
    error_message: str,
) -> Any:
    orchestrator_module = sys.modules.get("rph_core.orchestrator")
    pipeline_result_cls = getattr(orchestrator_module, "PipelineResult", None)
    if pipeline_result_cls is None:
        from types import SimpleNamespace

        return SimpleNamespace(
            success=False,
            product_smiles=product_smiles,
            work_dir=work_dir,
            error_step=error_step,
            error_message=error_message,
        )
    return pipeline_result_cls(
        success=False,
        product_smiles=product_smiles,
        work_dir=work_dir,
        error_step=error_step,
        error_message=error_message,
    )


class V3Scheduler:
    def __init__(
        self,
        hunter: Any,
        run_cfg: dict[str, Any],
        display: Optional[V3StageDisplay] = None,
        progress: Optional[V3RunProgress] = None,
    ) -> None:
        self.hunter = hunter
        self.run_cfg = run_cfg
        self.config = hunter.config
        self.output_root = normalize_path(str(run_cfg.get("output_root", "./rph_output")))
        self.v3_cfg = (self.config.get("scheduler", {}) or {}).get("v3", {}) or {}
        self._theory_sig = _theory_signature_string(self.config)
        self._thermo_cfg = self.config.get("thermo", {}) or {}
        self._default_temperature_k = _to_float(self._thermo_cfg.get("temperature_k")) or 298.15
        self.display = display
        self.progress = progress if progress is not None else V3RunProgress()
        self._reaction_progress_history: list[ReactionExecutionProgress] = []

        global_cfg = self.config.setdefault("global", {})
        if not global_cfg.get("small_molecule_cache_dir"):
            global_cache_dir = self.output_root / "small_molecules"
            global_cfg["small_molecule_cache_dir"] = str(global_cache_dir)

    def plan(self) -> list[ReactionJob]:
        planning: PlanningProgress = self.progress.planning
        run_cfg_for_build = copy.deepcopy(self.run_cfg)
        if str(run_cfg_for_build.get("source", "dataset")) == "dataset":
            dataset_cfg_value = run_cfg_for_build.get("dataset") or {}
            dataset_cfg = dict(dataset_cfg_value) if isinstance(dataset_cfg_value, dict) else {}
            dataset_cfg["reaction_profiles"] = self.config.get("reaction_profiles", {}) or {}
            run_cfg_for_build["dataset"] = dataset_cfg

        tasks = build_tasks_from_run_config(run_cfg_for_build, theory_signature=self._theory_sig)
        planning.total_rows = len(tasks)
        planning.state = V3StageState.RUNNING
        if self.display:
            self.display.stage0_planning_header(self.progress)
        _ensure_dir(self.output_root)

        from collections import defaultdict
        grouped: dict[str, list[TaskSpec]] = defaultdict(list)
        for task in tasks:
            reaction_id = getattr(task, "reaction_id", None) or sanitize_rx_id(task.rx_id)
            grouped[_safe_segment(reaction_id, "reaction_id")].append(task)

        planning.total_reactions = len(grouped)
        planning.total_conditions = sum(len(v) for v in grouped.values())
        if self.display:
            self.display.stage0_planning_header(self.progress)

        jobs: list[ReactionJob] = []
        for reaction_id, condition_tasks in grouped.items():
            representative = condition_tasks[0]
            reaction_root = get_reaction_root(self.output_root, reaction_id)
            _ensure_dir(reaction_root)

            row_ids = [str(getattr(t, "row_id", t.rx_id)) for t in condition_tasks]
            condition_ids = [str(getattr(t, "condition_id", f"COND_{t.rx_id}")) for t in condition_tasks]

            reaction_manifest = {
                "version": "rph_v3.0",
                "reaction_id": reaction_id,
                "reaction_cache_key": str(getattr(representative, "reaction_cache_key", "")),
                "row_ids": row_ids,
                "condition_ids": condition_ids,
                "theory_signature": self._theory_sig,
                "status": "PENDING",
            }
            manifest_path = reaction_root / "reaction_manifest.json"
            if manifest_path.exists():
                try:
                    existing = read_json_dict(manifest_path)
                    if isinstance(existing, dict):
                        for k, v in reaction_manifest.items():
                            if k not in existing or existing.get(k) in (None, ""):
                                existing[k] = v
                        reaction_manifest = existing
                except Exception:
                    pass
            write_json(manifest_path, reaction_manifest)

            condition_jobs: list[ConditionJob] = []
            for ct in condition_tasks:
                cid = _safe_segment(getattr(ct, "condition_id", f"COND_{ct.rx_id}"), "condition_id")
                c_root = reaction_root / "conditions" / cid
                cleaner_data = (ct.meta or {}).get("cleaner_data")
                cleaner_data = cleaner_data if isinstance(cleaner_data, dict) else {}
                temperature_c = _to_float(cleaner_data.get("temperature_c") or cleaner_data.get("temp_celsius"))
                temperature_k = _to_float(cleaner_data.get("temperature_K") or cleaner_data.get("temperature_k"))
                if temperature_k is None and temperature_c is not None:
                    temperature_k = temperature_c + 273.15
                if temperature_k is None:
                    temperature_k = self._default_temperature_k
                condition_jobs.append(ConditionJob(
                    reaction_id=reaction_id,
                    condition_id=cid,
                    condition_root=c_root,
                    task=ct,
                    temperature_k=temperature_k,
                ))

            dr_plan_path = reaction_root / "S0_Mechanism" / "dr_branch_plan.json"
            if not dr_plan_path.exists() and not self.run_cfg.get("dry_run", False):
                self._run_s0_for_planning(reaction_root, representative)

            plan_rec = planning.reaction_plans.setdefault(
                reaction_id, ReactionPlanRecord(reaction_id=reaction_id)
            )
            if dr_plan_path.exists():
                plan_rec.s0_status = "S0_Mechanism ready"
                plan_rec.state = V3ItemState.COMPLETED
            else:
                plan_rec.s0_status = "pending"
                plan_rec.state = (
                    V3ItemState.RUNNING
                    if not self.run_cfg.get("dry_run")
                    else V3ItemState.PENDING
                )
            plan_rec.condition_count = len(condition_jobs)
            if self.display:
                self.display.stage0_planning_update(self.progress)

            job = ReactionJob(
                reaction_id=reaction_id,
                reaction_root=reaction_root,
                representative=representative,
                conditions=condition_jobs,
                dr_plan_path=dr_plan_path,
                branches=[],
            )

            s0_cfg = self.config.get("s0", {}) or {}
            dr_cfg = (s0_cfg.get("dr_completion", {}) or {})
            dr_enabled = bool(dr_cfg.get("enabled", True))
            loaded_branch_ids: set[str] = set()
            if dr_enabled and dr_plan_path.exists():
                dr_plan = _load_dr_branch_plan(reaction_root)
                if dr_plan and isinstance(dr_plan, dict):
                    branches_data = dr_plan.get("branches") or []
                    if isinstance(branches_data, list):
                        for branch_entry in branches_data:
                            raw_bid = branch_entry.get("branch_id", "")
                            if not raw_bid:
                                continue
                            bid = _safe_segment(raw_bid, "branch_id")
                            loaded_branch_ids.add(bid)
                            pathway_id = str(branch_entry.get("pathway_id", ""))
                            generation_policy = str(branch_entry.get("generation_policy", ""))
                            bp_smiles = (
                                str(branch_entry.get("product_smiles") or "").strip()
                                or representative.product_smiles
                            )
                            flipped = [int(m) for m in (branch_entry.get("flipped_map_numbers") or [])]
                            fixed = [int(m) for m in (branch_entry.get("fixed_stereocenters") or [])]
                            bnotes = [str(n) for n in (branch_entry.get("notes") or [])]
                            branch_root = get_branch_root(reaction_root, bid)
                            _write_branch_manifest(
                                branch_root=branch_root,
                                branch_id=bid,
                                pathway_id=pathway_id,
                                product_smiles=bp_smiles,
                                generation_policy=generation_policy,
                                flipped_map_numbers=flipped,
                                fixed_stereocenters=fixed,
                                notes=bnotes,
                                dr_plan_path=str(dr_plan_path),
                                parent_product_smiles=representative.product_smiles,
                            )
                            job.branches.append(BranchJob(
                                parent_reaction_id=reaction_id,
                                branch_id=bid,
                                pathway_id=pathway_id,
                                product_smiles=bp_smiles,
                                branch_root=branch_root,
                                generation_policy=generation_policy,
                                flipped_map_numbers=flipped,
                                fixed_stereocenters=fixed,
                                notes=bnotes,
                            ))

            if "BR_MAJOR" not in loaded_branch_ids:
                major_root = get_branch_root(reaction_root, "BR_MAJOR")
                _write_branch_manifest(
                    branch_root=major_root,
                    branch_id="BR_MAJOR",
                    pathway_id="primary",
                    product_smiles=representative.product_smiles,
                    generation_policy="primary",
                    flipped_map_numbers=[],
                    fixed_stereocenters=[],
                    notes=["default major branch"],
                    dr_plan_path=str(dr_plan_path),
                    parent_product_smiles=representative.product_smiles,
                )
                job.branches.insert(0, BranchJob(
                    parent_reaction_id=reaction_id,
                    branch_id="BR_MAJOR",
                    pathway_id="primary",
                    product_smiles=representative.product_smiles,
                    branch_root=major_root,
                    generation_policy="primary",
                    notes=["default major branch"],
                ))

            plan_rec.branch_count = len(job.branches)
            plan_rec.branch_ids = [b.branch_id for b in job.branches]
            if self.display:
                self.display.stage0_planning_update(self.progress)

            jobs.append(job)

        planning.state = V3StageState.COMPLETED
        planning.branch_sequence_formed = True
        if self.display:
            self.display.stage0_planning_complete(self.progress)

        return jobs

    def _run_s0_for_planning(self, reaction_root: Path, representative: TaskSpec) -> None:
        from rph_core.utils.checkpoint_manager import CheckpointManager

        cleaner_data = (representative.meta or {}).get("cleaner_data")
        cleaner_data = cleaner_data if isinstance(cleaner_data, dict) else None
        try:
            self.hunter._run_s0(
                work_dir=reaction_root,
                product_smiles=representative.product_smiles,
                cleaner_data=cleaner_data,
                checkpoint_mgr=CheckpointManager(reaction_root),
                resume_enabled=bool(self.run_cfg.get("resume", True)),
                pm=None,
            )
        except Exception as exc:
            logger.warning(f"V3: S0 planning failed for {reaction_root.name}: {exc}")

    def precompute_small_molecules(self, jobs: list[ReactionJob]) -> None:
        self.hunter.logger.info("DEBUG: Precomputing small molecules...")
        cache_progress: SmallMoleculeCacheProgress = self.progress.small_molecule_cache
        cache_progress.state = V3StageState.RUNNING
        cache_progress.cache_root = str(
            self.config.get("global", {}).get("small_molecule_cache_dir", "")
        )
        if not self.v3_cfg.get("precompute_small_molecules", True):
            self.hunter.logger.info("DEBUG: Precomputing disabled by config.")
            cache_progress.state = V3StageState.COMPLETED
            return
        from rph_core.scheduling.small_molecule_precompute import SmallMoleculePrecomputer
        precomputer = SmallMoleculePrecomputer(
            config=self.config,
            cache_root=Path(self.config.get("global", {}).get("small_molecule_cache_dir", "")),
        )
        keys = precomputer.collect_keys(jobs)
        sorted_keys = sorted(keys)
        cache_progress.summary_total = len(keys)
        cache_progress.summary_hit = 0
        cache_progress.summary_computed = 0
        cache_progress.summary_failed = 0
        cache_progress.summary_pending = len(keys)
        for key in sorted_keys:
            cache_progress.molecules.setdefault(
                key, SmallMoleculeRecord(key=key)
            )
        if self.display:
            self.display.stage1_sm_cache_header(self.progress)
            self.display.stage1_sm_cache_update(self.progress)
        self.hunter.logger.info(f"DEBUG: Precomputing keys: {keys}")

        theory_signature = precomputer._build_theory_signature()
        for key in sorted_keys:
            rec = cache_progress.molecules[key]
            mol = precomputer.catalog.get(key)
            was_cached = False
            if mol is not None:
                try:
                    was_cached = precomputer.cache.is_complete(
                        mol.smiles,
                        theory_signature=theory_signature,
                        require_thermo=True,
                    )
                except Exception:
                    was_cached = False
            try:
                rec.state = V3ItemState.RUNNING
                rec.detail = "checking cache / compute"
                if self.display:
                    self.display.stage1_sm_cache_update(self.progress)
                precomputer.ensure_one(key)
                if was_cached:
                    rec.state = V3ItemState.CACHED
                    rec.detail = "cache hit"
                    cache_progress.summary_hit += 1
                else:
                    rec.state = V3ItemState.COMPLETED
                    rec.detail = "thermo.json ready"
                    cache_progress.summary_computed += 1
                cache_progress.summary_pending -= 1
            except Exception as exc:
                rec.state = V3ItemState.FAILED
                rec.detail = str(exc)[:60]
                cache_progress.summary_failed += 1
                cache_progress.summary_pending -= 1
            if self.display:
                self.display.stage1_sm_cache_update(self.progress)

        cache_progress.state = V3StageState.COMPLETED
        if self.display:
            self.display.stage1_sm_cache_complete(self.progress)
        self.hunter.logger.info("DEBUG: Precomputing finished.")

    def run(self) -> list[Any]:
        self.hunter.logger.info("DEBUG: V3Scheduler.run() called")
        jobs = self.plan()
        self.hunter.logger.info("DEBUG: V3Scheduler.plan() finished")
        if self.run_cfg.get("dry_run", False):
            for job in jobs:
                logger.info(
                    f"[DRY RUN] Reaction {job.reaction_id}: "
                    f"{len(job.conditions)} conditions, {len(job.branches)} DR branches"
                )
            return []

        self.precompute_small_molecules(jobs)

        results: list[Any] = []
        self.progress.total_reactions = len(jobs)
        for idx, job in enumerate(jobs):
            self.progress.current_reaction_index = idx + 1
            self.progress.current_reaction_id = job.reaction_id
            self.hunter.logger.info(
                "V3: Processing reaction %s (%d conditions, %d DR branches)",
                job.reaction_id,
                len(job.conditions),
                len(job.branches),
            )
            try:
                reaction_results = self._run_reaction(job)
                results.extend(reaction_results)
            except Exception as exc:
                logger.error(f"V3: Reaction {job.reaction_id} failed: {exc}")
                plan_rec = self.progress.planning.reaction_plans.get(job.reaction_id)
                if plan_rec is not None:
                    plan_rec.state = V3ItemState.FAILED
                    plan_rec.detail = str(exc)[:80]
                if (
                    self.progress.reaction_execution is not None
                    and self.progress.reaction_execution.reaction_id == job.reaction_id
                ):
                    self.progress.reaction_execution.state = V3StageState.FAILED
                results.append(_make_pipeline_failure_result(
                    product_smiles=job.representative.product_smiles,
                    work_dir=job.reaction_root,
                    error_step="v3_reaction",
                    error_message=str(exc),
                ))

        self._build_final_summary(results)
        if self.display:
            self.display.stage3_summary(self.progress)
        return results

    def _run_reaction(self, job: ReactionJob) -> list[Any]:
        results: list[Any] = []

        precursor_smiles = (job.representative.meta or {}).get("precursor_smiles")
        small_molecular_keys = (job.representative.meta or {}).get("small_molecular_keys")

        rxn_progress = ReactionExecutionProgress(
            reaction_id=job.reaction_id,
            state=V3StageState.RUNNING,
            total_branches=len(job.branches),
        )
        setattr(rxn_progress, "small_molecule_refs_state", V3ItemState.PENDING)
        setattr(rxn_progress, "small_molecule_refs_detail", "waiting")
        setattr(rxn_progress, "condition_progress", {})
        setattr(rxn_progress, "condition_dr_progress", {})
        for branch in job.branches:
            rxn_progress.branches[branch.branch_id] = BranchProgress(
                branch_id=branch.branch_id,
                s0_skipped_by_plan=True,
            )
        self.progress.reaction_execution = rxn_progress

        if self.display:
            self.display.stage2_reaction_header(self.progress)
            self.display.stage2_reaction_setup(self.progress)

        try:
            if self.v3_cfg.get("precursor_s1_once_per_reaction", True) and precursor_smiles:
                rxn_progress.shared_precursor_s1_state = V3ItemState.RUNNING
                rxn_progress.shared_precursor_s1_detail = "conformer search / DFT opt"
                if self.display:
                    self.display.stage2_reaction_setup(self.progress)
                self._run_precursor_s1(job)
                rxn_progress.shared_precursor_s1_state = V3ItemState.COMPLETED
                setattr(rxn_progress, "sm_cache_refs_ready", True)
                setattr(rxn_progress, "small_molecule_refs_state", V3ItemState.COMPLETED)
                setattr(rxn_progress, "small_molecule_refs_detail", "ready")
                if self.display:
                    self.display.stage2_reaction_setup(self.progress)

            reaction_features_dir = job.reaction_root / "reaction_features"

            for branch in job.branches:
                branch_rec = rxn_progress.branches[branch.branch_id]
                try:
                    from rph_core.scheduling.artifact_refs import (
                        link_or_copy_s0_into_branch,
                        link_or_copy_precursor_into_branch,
                        materialize_small_molecule_refs,
                    )

                    link_or_copy_s0_into_branch(job.reaction_root, branch.branch_root)

                    if self.v3_cfg.get("precursor_s1_once_per_reaction", True) and precursor_smiles:
                        link_or_copy_precursor_into_branch(job.reaction_root, branch.branch_root)
                        branch_rec.precursor_linked = True

                        if small_molecular_keys:
                            cache_root = Path(self.config.get("global", {}).get("small_molecule_cache_dir", ""))
                            materialize_small_molecule_refs(
                                branch.branch_root,
                                cache_root,
                                small_molecular_keys,
                                config=self.config,
                            )
                            branch_rec.sm_refs_materialized = True

                    rxn_progress.current_branch_index = rxn_progress.current_branch_index + 1
                    branch_rec.state = V3ItemState.RUNNING
                    branch_rec.pipeline_step = "s1"
                    branch_rec.pipeline_step_state = V3ItemState.RUNNING
                    if self.display:
                        self.display.stage2_branch_start(self.progress, branch.branch_id)

                    branch_result = self._run_branch(job, branch)
                    results.append(branch_result)

                    branch_rec = rxn_progress.branches[branch.branch_id]
                    if branch_result.success:
                        branch_rec.state = V3ItemState.COMPLETED
                        branch_rec.pipeline_step = "s4"
                        branch_rec.pipeline_step_state = V3ItemState.COMPLETED
                        try:
                            sp_report = getattr(branch_result, "sp_matrix_report", None)
                            if sp_report:
                                dg = sp_report.get_activation_energy()
                                if dg is not None:
                                    branch_rec.dg_activation = dg
                        except Exception:
                            pass
                    else:
                        branch_rec.state = V3ItemState.FAILED
                        branch_rec.pipeline_step_state = V3ItemState.FAILED
                        error_step = getattr(branch_result, "error_step", "")
                        branch_rec.detail = f"{error_step} failed" if error_step else "failed"
                    if self.display:
                        self.display.stage2_branch_complete(self.progress, branch.branch_id)

                    bm_path = branch.branch_root / "branch_manifest.json"
                    if bm_path.exists():
                        with open(bm_path, "r", encoding="utf-8") as bf:
                            bm = json.load(bf)
                        bm["status"] = "COMPLETE" if branch_result.success else "FAILED"
                        with open(bm_path, "w", encoding="utf-8") as bf:
                            json.dump(bm, bf, indent=2, ensure_ascii=False)
                except Exception as exc:
                    branch_rec = rxn_progress.branches[branch.branch_id]
                    branch_rec.state = V3ItemState.FAILED
                    branch_rec.pipeline_step_state = V3ItemState.FAILED
                    branch_rec.detail = str(exc)[:80]
                    if self.display:
                        self.display.stage2_branch_complete(self.progress, branch.branch_id)
                    logger.warning(f"DR branch {branch.branch_id} failed for {job.reaction_id}: {exc}")

            if job.branches:
                try:
                    rxn_progress.dr_aggregation_state = V3ItemState.RUNNING
                    rxn_progress.dr_aggregation_detail = "writing dr_prediction_default.json"
                    if self.display:
                        self.display.stage2_dr_aggregation(self.progress)
                    from rph_core.steps.dr_aggregator import DRAggregator
                    DRAggregator().write(
                        reaction_root=job.reaction_root,
                        output_path=reaction_features_dir / "dr_prediction_default.json",
                        temperature_k=self._default_temperature_k,
                    )
                    rxn_progress.dr_aggregation_state = V3ItemState.COMPLETED
                    if self.display:
                        self.display.stage2_dr_aggregation(self.progress)
                except Exception as exc:
                    rxn_progress.dr_aggregation_state = V3ItemState.FAILED
                    rxn_progress.dr_aggregation_detail = str(exc)[:80]
                    if self.display:
                        self.display.stage2_dr_aggregation(self.progress)
                    logger.warning(f"DR aggregation failed for {job.reaction_id}: {exc}")

            rxn_progress.condition_post_processing = True
            rxn_progress.condition_total = len(job.conditions)
            rxn_progress.condition_branch_total = len(job.conditions) * len(job.branches)
            if self.display:
                self.display.stage2_condition_header(self.progress)

            self._run_conditions(job)

            try:
                manifest_path = job.reaction_root / "reaction_manifest.json"
                updated = read_json_dict(manifest_path)
                if isinstance(updated, dict):
                    updated["status"] = "COMPLETE"
                    write_json(manifest_path, updated)
            except Exception:
                pass

            rxn_progress.state = V3StageState.COMPLETED
            if self.display:
                self.display.stage2_reaction_complete(self.progress)

            return results
        except Exception:
            rxn_progress.state = V3StageState.FAILED
            raise
        finally:
            if not self._reaction_progress_history or self._reaction_progress_history[-1] is not rxn_progress:
                self._reaction_progress_history.append(rxn_progress)

    def _run_precursor_s1(self, job: ReactionJob) -> Path:
        from rph_core.utils.path_manager import get_precursor_s1_dir

        precursor_dir = get_precursor_s1_dir(job.reaction_root)
        precursor_smiles = (job.representative.meta or {}).get("precursor_smiles")
        if not precursor_smiles:
            return precursor_dir

        _ensure_dir(precursor_dir)
        min_xyz = precursor_dir / "precursor" / "precursor_min.xyz"
        global_min_xyz = precursor_dir / "precursor" / "precursor_global_min.xyz"
        if min_xyz.exists():
            return precursor_dir
        if global_min_xyz.exists():
            min_xyz.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(global_min_xyz, min_xyz)
            logger.info(f"V3: Reused existing precursor global minimum for {job.reaction_id}")
            return precursor_dir

        self.hunter.s1_engine.base_work_dir = precursor_dir
        anchor_result = self.hunter.s1_engine.run(
            molecules={"precursor": precursor_smiles},
        )
        if anchor_result.success:
            if global_min_xyz.exists() and not min_xyz.exists():
                shutil.copy2(global_min_xyz, min_xyz)
            logger.info(f"V3: Precursor S1 complete for reaction {job.reaction_id}")
        else:
            logger.warning(f"V3: Precursor S1 failed for reaction {job.reaction_id}: {anchor_result.error_message}")

        return precursor_dir

    def _run_branch(self, job: ReactionJob, branch: BranchJob) -> Any:
        configured_skip = [s for s in (self.run_cfg.get("skip_steps", [])) if str(s).strip()]
        skip_steps = ["s0"] + [s for s in configured_skip if str(s).lower() != "s4"]
        reaction_profile = (job.representative.meta or {}).get("reaction_profile") or self.run_cfg.get("reaction_profile")
        cleaner_data = (job.representative.meta or {}).get("cleaner_data")
        cleaner_data = cleaner_data if isinstance(cleaner_data, dict) else None

        branch_skip = skip_steps
        if self.v3_cfg.get("branch_product_only_s1", True):
            branch_skip = skip_steps

        logger.info(f"V3: Running branch {branch.branch_id} for reaction {job.reaction_id}")
        result = self.hunter.run_pipeline(
            product_smiles=branch.product_smiles,
            work_dir=branch.branch_root,
            skip_steps=branch_skip,
            precursor_smiles=None,
            leaving_group_key=None,
            small_molecular_keys=[],
            reaction_profile=reaction_profile,
            cleaner_data=cleaner_data,
            include_reference_small_molecules=False,
            reaction_id=job.reaction_id,
            branch_id=branch.branch_id,
            quiet_resume=True,
        )
        logger.info(f"V3: Branch {branch.branch_id} complete → {result.features_csv}")
        return result

    def _run_conditions(self, job: ReactionJob) -> None:
        rxn_progress = self.progress.reaction_execution
        condition_progress = getattr(rxn_progress, "condition_progress", None) if rxn_progress is not None else None
        condition_dr_progress = getattr(rxn_progress, "condition_dr_progress", None) if rxn_progress is not None else None
        if not self.v3_cfg.get("condition_branch_thermo", True):
            for cond in job.conditions:
                self._run_condition_for_root(job, cond, job.reaction_root)
            return

        all_branch_roots: list[tuple[str, Path]] = [
            (branch.branch_id, branch.branch_root) for branch in job.branches
        ]

        for cond in job.conditions:
            _ensure_dir(cond.condition_root)
            self._write_condition_manifest(cond, job)

            for branch_id, branch_root in all_branch_roots:
                condition_branch_root = cond.condition_root / "branches" / branch_id
                _ensure_dir(condition_branch_root)
                self._write_condition_manifest(cond, job, output_root=condition_branch_root)

                try:
                    from rph_core.steps.condition_thermo import ConditionThermoCalculator
                    ConditionThermoCalculator(
                        config=self.hunter.config,
                        reaction_root=job.reaction_root,
                        condition_root=condition_branch_root,
                        artifact_root=branch_root,
                    ).run()
                    if rxn_progress is not None:
                        cb_key = f"{cond.condition_id}/{branch_id}"
                        rxn_progress.condition_records[cb_key] = V3ItemState.COMPLETED
                        if isinstance(condition_progress, dict):
                            condition_progress[cb_key] = {
                                "condition_id": cond.condition_id,
                                "branch_id": branch_id,
                                "state": V3ItemState.COMPLETED,
                                "detail": "ready",
                            }
                        if self.display:
                            self.display.stage2_condition_update(
                                self.progress,
                                f"{cond.condition_id} × {branch_id}",
                            )
                except Exception as exc:
                    if rxn_progress is not None:
                        cb_key = f"{cond.condition_id}/{branch_id}"
                        rxn_progress.condition_records[cb_key] = V3ItemState.FAILED
                        if isinstance(condition_progress, dict):
                            condition_progress[cb_key] = {
                                "condition_id": cond.condition_id,
                                "branch_id": branch_id,
                                "state": V3ItemState.FAILED,
                                "detail": str(exc)[:80],
                            }
                        if self.display:
                            self.display.stage2_condition_update(
                                self.progress,
                                f"{cond.condition_id} × {branch_id}",
                            )
                    logger.warning(
                        f"Condition thermo failed for {cond.condition_id}/{branch_id}: {exc}"
                    )

                try:
                    from rph_core.steps.condition_feature_merger import ConditionFeatureMerger
                    ConditionFeatureMerger(
                        reaction_root=job.reaction_root,
                        condition_root=condition_branch_root,
                        feature_root=branch_root,
                        branch_id=branch_id,
                    ).run()
                except Exception as exc:
                    logger.warning(
                        f"Condition merge failed for {cond.condition_id}/{branch_id}: {exc}"
                    )

            try:
                from rph_core.steps.dr_aggregator import DRAggregator
                DRAggregator().run_for_condition(
                    reaction_root=job.reaction_root,
                    condition_root=cond.condition_root,
                    temperature_k=cond.temperature_k,
                )
                if rxn_progress is not None:
                    rxn_progress.condition_dr_state[cond.condition_id] = V3ItemState.COMPLETED
                    if isinstance(condition_dr_progress, dict):
                        condition_dr_progress[cond.condition_id] = {
                            "condition_id": cond.condition_id,
                            "state": V3ItemState.COMPLETED,
                            "detail": "ready",
                        }
            except Exception as exc:
                if rxn_progress is not None:
                    rxn_progress.condition_dr_state[cond.condition_id] = V3ItemState.FAILED
                    if isinstance(condition_dr_progress, dict):
                        condition_dr_progress[cond.condition_id] = {
                            "condition_id": cond.condition_id,
                            "state": V3ItemState.FAILED,
                            "detail": str(exc)[:80],
                        }
                logger.warning(
                    f"Condition DR aggregation failed for {cond.condition_id}: {exc}"
                )

    def _build_final_summary(self, results: list[Any]) -> None:
        summary: FinalSummary = self.progress.summary
        summary.total_reactions = self.progress.total_reactions
        summary.reaction_success = 0
        summary.reaction_failed = 0
        summary.reaction_records = {}
        summary.total_branches = 0
        summary.branch_success = 0
        summary.branch_failed = 0
        summary.branch_records = {}
        summary.total_conditions = 0
        summary.condition_complete = 0
        summary.condition_failed = 0
        summary.failed_planning = []
        summary.failed_cache = []
        summary.failed_branches = []
        summary.failed_conditions = []
        summary.output_root = str(self.output_root)

        for r in results:
            success = getattr(r, "success", False)
            if success:
                summary.branch_success += 1
            else:
                summary.branch_failed += 1
        summary.total_branches = summary.branch_success + summary.branch_failed

        for rxn_id, plan_rec in self.progress.planning.reaction_plans.items():
            if plan_rec.state == V3ItemState.FAILED:
                summary.reaction_failed += 1
                summary.failed_planning.append(rxn_id)
                summary.reaction_records[rxn_id] = "failed"
            else:
                summary.reaction_success += 1
                summary.reaction_records[rxn_id] = "completed"

        for key, rec in self.progress.small_molecule_cache.molecules.items():
            if rec.state == V3ItemState.FAILED:
                summary.failed_cache.append(key)

        for rxn in self._reaction_progress_history:
            for bid, bp in rxn.branches.items():
                branch_key = f"{rxn.reaction_id}/{bid}"
                summary.branch_records[branch_key] = bp.state.label
                if bp.state == V3ItemState.FAILED:
                    summary.failed_branches.append(branch_key)

            for cond_id, state in rxn.condition_dr_state.items():
                cond_key = f"{rxn.reaction_id}/{cond_id}"
                if state == V3ItemState.FAILED:
                    summary.condition_failed += 1
                    summary.failed_conditions.append(cond_key)
                elif state in {
                    V3ItemState.COMPLETED,
                    V3ItemState.CACHED,
                    V3ItemState.THERMO_DERIVED,
                }:
                    summary.condition_complete += 1

        if self._reaction_progress_history:
            summary.branch_success = 0
            summary.branch_failed = 0
            summary.total_branches = 0
            for rxn in self._reaction_progress_history:
                summary.total_branches += len(rxn.branches)
                for bid, bp in rxn.branches.items():
                    branch_key = f"{rxn.reaction_id}/{bid}"
                    if bp.state == V3ItemState.FAILED:
                        summary.branch_failed += 1
                        if branch_key not in summary.failed_branches:
                            summary.failed_branches.append(branch_key)
                    elif bp.state in {
                        V3ItemState.COMPLETED,
                        V3ItemState.CACHED,
                        V3ItemState.THERMO_DERIVED,
                    }:
                        summary.branch_success += 1
            summary.total_conditions = self.progress.planning.total_conditions

        summary.sm_cache_total = self.progress.small_molecule_cache.summary_total
        summary.sm_cache_hit = self.progress.small_molecule_cache.summary_hit
        summary.sm_cache_computed = self.progress.small_molecule_cache.summary_computed
        summary.sm_cache_failed = self.progress.small_molecule_cache.summary_failed

    def _run_condition_for_root(self, job: ReactionJob, cond: ConditionJob, artifact_root: Path) -> None:
        _ensure_dir(cond.condition_root)
        self._write_condition_manifest(cond, job)

        try:
            from rph_core.steps.condition_thermo import ConditionThermoCalculator
            ConditionThermoCalculator(
                config=self.hunter.config,
                reaction_root=job.reaction_root,
                condition_root=cond.condition_root,
                artifact_root=artifact_root,
            ).run()
        except Exception as exc:
            logger.warning(f"Condition thermo failed for {cond.condition_id}: {exc}")

        try:
            from rph_core.steps.condition_feature_merger import ConditionFeatureMerger
            ConditionFeatureMerger(
                reaction_root=job.reaction_root,
                condition_root=cond.condition_root,
                feature_root=artifact_root,
            ).run()
        except Exception as exc:
            logger.warning(f"Condition merge failed for {cond.condition_id}: {exc}")

    def _write_condition_manifest(
        self,
        cond: ConditionJob,
        job: ReactionJob,
        output_root: Optional[Path] = None,
    ) -> None:
        cleaner_data = (cond.task.meta or {}).get("cleaner_data")
        cleaner_data = cleaner_data if isinstance(cleaner_data, dict) else {}
        temperature_c = _to_float(cleaner_data.get("temperature_c") or cleaner_data.get("temp_celsius"))

        condition_manifest = {
            "version": "rph_v3.0",
            "condition_id": cond.condition_id,
            "reaction_id": job.reaction_id,
            "row_id": str(getattr(cond.task, "row_id", cond.task.rx_id)),
            "temperature_c": temperature_c,
            "temperature_K": cond.temperature_k,
            "solvent": cleaner_data.get("solvent") or "DCM",
            "catalyst": cleaner_data.get("catalyst"),
            "additive": cleaner_data.get("additive"),
            "oxidant": cleaner_data.get("oxidant"),
            "yield_pct": _to_float(cleaner_data.get("yield_pct") or cleaner_data.get("yield")),
            "dr_major": _to_float(cleaner_data.get("dr_major")),
            "dr_minor": _to_float(cleaner_data.get("dr_minor")),
            "ee_pct": _to_float(cleaner_data.get("ee_pct") or cleaner_data.get("ee")),
            "source_ref": cleaner_data.get("source_ref"),
        }
        target_root = output_root or cond.condition_root
        write_json(target_root / "condition_manifest.json", condition_manifest)
