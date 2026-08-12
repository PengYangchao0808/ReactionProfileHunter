"""Fixed V4 CENSO-LITE conformer workflow.

The workflow deliberately stops before geometry-level DFT refinement.  It uses
CREST/GFN2 geometries, B97-3c single-point ranking and optional xTB mRRHO
corrections.  High-level refinement is owned by S4.

Pipeline order:
    1.  RDKit embed → initial.xyz
    2.  CREST GFN2/ALPB → crest_conformers.xyz
    3.  Split ensemble → N conf_####.xyz
    4.  Parse xTB energies (fail-fast) + cross-validate with crest.energies
    5.  Torsion-aware dedup (geometric, cheap)
    6.  Optional prefilter (disabled by default)
    7.  Parallel ORCA B97-3c/CPCM SP on every unique member
    8.  Parallel xTB mRRHO on the same complete ensemble
    9.  Rank with G_i = E_B97-3c,i + G(RRHO)_xTB,i
    10. Calculate Z_rel, Boltzmann populations and G_conf_rel
    11. Select representative structures without truncating the ensemble
"""

from __future__ import annotations

import copy
import json
import logging
import re
import shutil
import sys
import threading
import time
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional, Sequence

from rdkit import Chem

from rph_core.steps.conformer_search.censo_lite_runtime import CensoLiteRuntime, CrestEnergyParseError
from rph_core.steps.conformer_search.deduplicator import DedupCandidate, TorsionAwareDeduplicator
from rph_core.steps.conformer_search.ensemble_thermo import calculate_ensemble_thermodynamics
from rph_core.steps.conformer_search.torsion_signature import signatures_equivalent
from rph_core.utils.atom_mapping import bind_atom_mapping_to_xyz
from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.utils.json_io import write_json_atomic
from rph_core.utils.provenance import build_provenance, sha256_of_file
from rph_core.utils.resource_utils import mem_to_mb
from rph_core.utils.run_id import RUN_ID_FIELD

logger = logging.getLogger(__name__)

S1_MANIFEST_SCHEMA_VERSION = "s1_censo_light_ranking_v4"


@dataclass(frozen=True)
class _WavePlan:
    """Operational layout for one homogeneous S1 batch."""

    bulk_count: int
    final_cores: tuple[int, ...]
    policy: str


class CensoLiteEngine:
    """Run the single supported S1 protocol and emit a stable manifest."""

    protocol_name = "censo_lite"
    manifest_schema_version = S1_MANIFEST_SCHEMA_VERSION

    def __init__(
        self,
        config: Dict[str, Any],
        work_dir: Path,
        molecule_name: str,
        event_callback: Callable[[str, Dict[str, Any]], None] | None = None,
        run_id: Optional[str] = None,
        s0_mechanism_json_path: Path | None = None,
        forming_bonds: Sequence[Sequence[int]] | None = None,
    ):
        self.config = copy.deepcopy(config)
        step1 = self.config.setdefault("step1", {})
        step1["protocol"] = self.protocol_name
        lite = dict(step1.get("censo_lite", {}) or {})
        ranking = dict(lite.get("ranking", {}) or {})
        profiles = dict(step1.get("fast_sp_profiles", {}) or {})
        profiles["censo_lite_gga"] = ranking
        step1["fast_sp_profiles"] = profiles
        self.config["step1"] = step1
        self.work_dir = Path(work_dir)
        self.molecule_name = molecule_name
        self.engine = CensoLiteRuntime(self.config, self.work_dir, molecule_name)
        self.deduplicator = TorsionAwareDeduplicator(lite.get("deduplication", {}))
        self.lite_config = lite
        self._event_callback = event_callback
        self.run_id = run_id
        self.s0_mechanism_json_path = (
            Path(s0_mechanism_json_path) if s0_mechanism_json_path is not None else None
        )
        self.forming_bonds = forming_bonds
        if event_callback is not None and hasattr(self.engine, "set_event_callback"):
            self.engine.set_event_callback(lambda event, payload: self._emit(event, **payload))

    def _emit(self, event: str, **fields: Any) -> None:
        if self._event_callback is None:
            return
        try:
            self._event_callback(event, {"variant": self.molecule_name, **fields})
        except Exception as exc:  # pragma: no cover - UI must not affect QC
            logger.warning("Ignoring S1 UI callback failure for %s: %s", event, exc)

    def set_event_callback(
        self, callback: Callable[[str, Dict[str, Any]], None] | None
    ) -> None:
        self._event_callback = callback
        if hasattr(self.engine, "set_event_callback"):
            self.engine.set_event_callback(
                (lambda event, payload: self._emit(event, **payload)) if callback is not None else None
            )

    @contextmanager
    def _step(
        self,
        step_id: str,
        index: int,
        label: str,
        *,
        output_path: Path | None = None,
        **detail: Any,
    ) -> Iterator[Dict[str, Any]]:
        started = time.time()
        result: Dict[str, Any] = {}
        self._emit(
            "step_started",
            step=step_id,
            index=index,
            total_steps=10,
            label=label,
            started_at=started,
            **detail,
        )
        stop = threading.Event()
        ui_cfg = dict(self.config.get("ui", {}) or {})
        interval_key = "heartbeat_seconds" if sys.stdout.isatty() else "non_tty_heartbeat_seconds"
        interval = max(1.0, float(ui_cfg.get(interval_key, ui_cfg.get("heartbeat_seconds", 30))))

        def _heartbeat() -> None:
            while not stop.wait(interval):
                activity = self._output_activity(output_path)
                self._emit(
                    "step_heartbeat",
                    step=step_id,
                    index=index,
                    total_steps=10,
                    label=label,
                    started_at=started,
                    elapsed_seconds=time.time() - started,
                    **activity,
                )

        thread: threading.Thread | None = None
        if self._event_callback is not None:
            thread = threading.Thread(target=_heartbeat, name=f"rph-ui-{step_id}", daemon=True)
            thread.start()
        try:
            yield result
        except Exception as exc:
            self._emit(
                "step_failed",
                step=step_id,
                index=index,
                total_steps=10,
                label=label,
                started_at=started,
                finished_at=time.time(),
                elapsed_seconds=time.time() - started,
                error=str(exc),
                **self._output_activity(output_path),
            )
            raise
        else:
            finished_payload = self._output_activity(output_path)
            finished_payload.update(result)
            self._emit(
                "step_finished",
                step=step_id,
                index=index,
                total_steps=10,
                label=label,
                started_at=started,
                finished_at=time.time(),
                elapsed_seconds=time.time() - started,
                status=str(finished_payload.pop("status", "complete")),
                **finished_payload,
            )
        finally:
            stop.set()
            if thread is not None:
                thread.join(timeout=0.2)

    @staticmethod
    def _output_activity(path: Path | None) -> Dict[str, Any]:
        if path is None:
            return {}
        target = Path(path)
        try:
            if target.is_file():
                stat = target.stat()
                return {
                    "output": str(target),
                    "output_size_bytes": stat.st_size,
                    "last_output_age_seconds": max(0.0, time.time() - stat.st_mtime),
                }
            if target.is_dir():
                files = [item for item in target.rglob("*") if item.is_file()]
                if not files:
                    return {"output": str(target), "output_size_bytes": 0}
                stats = [item.stat() for item in files]
                return {
                    "output": str(target),
                    "output_size_bytes": sum(stat.st_size for stat in stats),
                    "last_output_age_seconds": max(0.0, time.time() - max(stat.st_mtime for stat in stats)),
                }
        except OSError:
            pass
        return {"output": str(target)}

    def run(self, smiles: str) -> Dict[str, Any]:
        output_dir = self.engine.molecule_dir
        pipeline_log: List[Dict[str, Any]] = []

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES for CENSO-LITE: {smiles}")
        mol = Chem.AddHs(mol)

        with self._step(
            "rdkit_embed", 1, "RDKit ETKDG embedding",
            output_path=output_dir / "initial.xyz",
            purpose="generate initial 3D geometry",
            engine="rdkit", method="ETKDGv3",
        ) as step:
            initial_xyz = self.engine.embed(smiles)
            step.update(output=str(initial_xyz), atoms=mol.GetNumAtoms())

        crest_cfg = dict(self.lite_config.get("crest", {}) or {})
        with self._step(
            "crest_search", 2, "CREST/GFN2 conformer search",
            output_path=self.engine.crest_dir,
            purpose="sample the conformational ensemble",
            engine="crest", method=f"GFN{crest_cfg.get('gfn_level', 2)}-xTB",
            solvent=f"ALPB({crest_cfg.get('solvent', 'acetone')})",
            search_mode=crest_cfg.get("search_mode", "imtd_gc"),
            energy_window_kcal=crest_cfg.get("energy_window_kcal", 6.0),
            nprocs=(self.config.get("resources", {}) or {}).get("nproc"),
        ) as step:
            ensemble_xyz = self.engine.crest_search(initial_xyz)
            step.update(output=str(ensemble_xyz))

        with self._step(
            "split_ensemble", 3, "Split CREST ensemble",
            output_path=self.engine.raw_dir,
            purpose="write one XYZ file per CREST conformer",
        ) as step:
            candidate_paths = self.engine.split_ensemble(ensemble_xyz)
            if not candidate_paths:
                raise RuntimeError("CENSO-LITE produced no CREST candidates")
            step.update(input_conformers=len(candidate_paths), written=len(candidate_paths))
        pipeline_log.append({"step": "crest", "candidates": len(candidate_paths)})
        self._emit("funnel_progress", step="crest_generated", candidates=len(candidate_paths))

        # ---- Stage 4: parse xTB energies (fail-fast) + cross-validate ----
        with self._step(
            "energy_validation", 4, "Parse and validate CREST energies",
            output_path=self.engine.crest_dir,
            purpose="cross-check embedded xTB energies against crest.energies",
        ) as step:
            records = self._parse_and_annotate(mol, candidate_paths)
            step.update(parsed=len(records), invalid=len(candidate_paths) - len(records), validation="matched")
        pipeline_log.append({"step": "energy_parse", "candidates": len(records)})
        self._emit("funnel_progress", step="energy_valid", candidates=len(records))

        # ---- Stage 5: torsion-aware dedup (geometric, cheap) ----
        before_dedup = len(records)
        with self._step(
            "torsion_dedup", 5, "Torsion/RMSD deduplication",
            purpose="remove geometrically equivalent conformers",
        ) as step:
            extra_dedup = bool(self._prefilter_cfg().get("extra_deduplication", True))
            if extra_dedup:
                records = self.deduplicator.deduplicate(mol, records)
                pipeline_log.append({"step": "torsion_dedup", "candidates": len(records)})
            step.update(input=before_dedup, unique=len(records), removed=before_dedup - len(records))
        self._emit("funnel_progress", step="torsion_unique", candidates=len(records))

        with self._step(
            "prefilter", 6, "Optional conformer prefilter",
            purpose="apply configured xTB window and diversity cap",
        ) as step:
            prefilter_input = len(records)
            prefilter_cfg = self._prefilter_cfg()
            xtb_window_raw = prefilter_cfg.get("xtb_window_kcal")
            if xtb_window_raw is not None:
                xtb_window = float(xtb_window_raw)
                records = self._apply_energy_window(records, xtb_window)
                pipeline_log.append({"step": "xtb_window", "window_kcal": xtb_window, "candidates": len(records)})
            max_soft_raw = prefilter_cfg.get("max_candidates_soft")
            if max_soft_raw is not None and int(max_soft_raw) > 0:
                max_soft = int(max_soft_raw)
                records = self._apply_diversity_cap(records, max_soft)
                pipeline_log.append({"step": "diversity_cap_pre_sp", "max": max_soft, "candidates": len(records)})
            step.update(
                input=prefilter_input, output_candidates=len(records),
                window_kcal=xtb_window_raw, max_candidates=max_soft_raw,
                status="skipped" if xtb_window_raw is None and max_soft_raw is None else "complete",
            )
        self._emit("funnel_progress", step="prefiltered", candidates=len(records))

        # ---- Stage 8: parallel ORCA B97-3c SP ----
        with self._step(
            "b97_sp", 7, "ORCA B97-3c SP ranking",
            output_path=getattr(self.engine, "ranking_dir", output_dir / "ranking"),
            purpose="rank every unique CREST conformer with B97-3c",
            engine="orca", method="B97-3c",
        ) as step:
            records = self._run_b97_sp_parallel(mol, records)
            sp_failed = [r for r in records if r.metadata.get("sp_status") == "failed"]
            b97_failure_policy = str(
                self._b97_cfg().get("failure_policy", "strict")
            ).strip().lower()
            if b97_failure_policy != "strict":
                raise ValueError(
                    "step1.censo_lite.b97_3c.failure_policy must be 'strict'"
                )
            if sp_failed:
                raise RuntimeError(
                    f"B97-3c SP failed for {len(sp_failed)} ensemble candidates; "
                    "refusing to calculate a truncated partition function"
                )
            records = [r for r in records if r.metadata.get("sp_status") != "failed"]
            if not records:
                raise RuntimeError("All B97-3c SP calculations failed")
            step.update(valid=len(records), failed=len(sp_failed))
            pipeline_log.append({"step": "b97_sp", "candidates": len(records), "failed": len(sp_failed)})
        self._emit("funnel_progress", step="b97_valid", candidates=len(records))

        # ---- Stage 9: B97-3c energy window ----
        # No B97 electronic-energy filter is allowed before thermochemistry.

        # ---- Stage 10: diversity cap → thermo pool ----
        xtb_thermo_cfg = self.lite_config.get("xtb_thermo", {}) or {}
        # The complete B97-ranked ensemble enters mRRHO and the partition sum.

        # ---- Stage 11-12: mRRHO on entire thermo pool + final ranking ----
        thermo_enabled = bool(xtb_thermo_cfg.get("enabled", True))
        scoring_mode = "electronic_only"
        thermochemistry_status = "not_requested"
        with self._step(
            "mrrho", 8, "xTB mRRHO corrections",
            output_path=getattr(self.engine, "mrrho_dir", output_dir / "mrrho"),
            purpose="calculate per-conformer thermostatistical corrections",
            engine="xtb", method=f"GFN{xtb_thermo_cfg.get('gfn_level', 2)}-xTB SPH+mRRHO",
            temperature_k=xtb_thermo_cfg.get("temperature_k", 298.15),
        ) as step:
            if thermo_enabled:
                records, scoring_mode, thermochemistry_status = self._apply_mrrho_uniform(
                    records, xtb_thermo_cfg
                )
                pipeline_log.append(
                    {
                        "step": "mrrho",
                        "candidates": len(records),
                        "scoring_mode": scoring_mode,
                        "thermochemistry_status": thermochemistry_status,
                    }
                )
                step.update(
                    valid=len(records),
                    scoring_mode=scoring_mode,
                    thermochemistry_status=thermochemistry_status,
                )
            else:
                step.update(status="skipped", reason="xtb_thermo.enabled=false")
        self._emit("funnel_progress", step="mrrho_valid", candidates=len(records))

        with self._step(
            "ensemble_thermo", 9, "Ensemble thermodynamics",
            purpose="calculate partition weights, populations and Gconf",
            temperature_k=xtb_thermo_cfg.get("temperature_k", 298.15),
        ) as step:
            records.sort(key=lambda item: (item.score, str(item.path)))
            records, ensemble_thermo = self._annotate_ensemble_thermodynamics(
                records, thermochemistry_status
            )
            step.update(
                ensemble_members=len(records),
                partition_function_relative=ensemble_thermo.get("partition_function_relative"),
                conformational_free_energy_correction_kcal=ensemble_thermo.get(
                    "conformational_free_energy_correction_kcal"
                ),
            )

        # ---- Stage 13: final selection ----
        final_window = float(self._selection_cfg().get("final_window_kcal", 2.5))
        min_keep = int(self._selection_cfg().get("min_keep", 3))
        max_keep = int(self._selection_cfg().get("max_keep", 10))
        with self._step(
            "selection_manifest", 10, "Representative selection and manifest",
            output_path=output_dir,
            purpose="select downstream representatives and write the S1 contract",
            final_window_kcal=final_window, min_keep=min_keep, max_keep=max_keep,
        ) as step:
            representatives = self._apply_final_selection(records, final_window, min_keep, max_keep)
            pipeline_log.append(
                {
                    "step": "representative_selection",
                    "ensemble_candidates": len(records),
                    "representatives": len(representatives),
                }
            )
            result = self._write_manifest(
                records, representatives, output_dir, pipeline_log, scoring_mode, ensemble_thermo,
            )
            data = dict(result.get("data") or {})
            in_window = sum(
                float(
                    row.get("relative_free_energy_kcal")
                    if row.get("relative_free_energy_kcal") is not None
                    else row.get("relative_s1_score_kcal") or 0.0
                ) <= final_window
                for row in data.get("candidates", [])
            )
            selected = next(
                (row for row in data.get("candidates", []) if row.get("id") == data.get("selected")),
                {},
            )
            step.update(
                in_window=in_window, representatives=len(representatives),
                selected=data.get("selected"),
                selected_population=selected.get("boltzmann_population"),
                manifest=str(result.get("manifest")), geometry=str(result.get("selected_xyz")),
            )
        return result

    def _prefilter_cfg(self) -> Dict[str, Any]:
        return dict(self.lite_config.get("prefilter", {}) or {})

    def _b97_cfg(self) -> Dict[str, Any]:
        return dict(self.lite_config.get("b97_3c", {}) or {})

    def _selection_cfg(self) -> Dict[str, Any]:
        return dict(self.lite_config.get("selection", {}) or {})

    def _parse_and_annotate(self, mol: Chem.Mol, paths: List[Path]) -> List[DedupCandidate]:
        xtb_energies: List[float] = []
        for path in paths:
            xtb_energies.append(self.engine.extract_energy(path))

        if all(abs(e) < 1e-12 for e in xtb_energies):
            raise CrestEnergyParseError("All CREST energies parsed as zero")

        rel_energies = self.engine.read_crest_relative_energies(self.engine.crest_dir)
        self.engine.validate_energies(xtb_energies, rel_energies)

        records: List[DedupCandidate] = []
        for path, xtb_energy in zip(paths, xtb_energies):
            records.append(
                self.deduplicator.annotate(
                    mol,
                    path,
                    float(xtb_energy),
                    {
                        "crest_gfn2_electronic_alpb_hartree": float(xtb_energy),
                        # Deprecated v2 alias.
                        "xtb_energy_hartree": float(xtb_energy),
                        "b973c_sp_energy_hartree": None,
                        "g_rrho_correction_hartree": None,
                    },
                )
            )
        return records

    def _apply_energy_window(
        self, records: List[DedupCandidate], window_kcal: float
    ) -> List[DedupCandidate]:
        if not records:
            return records
        minimum = min(r.score for r in records)
        kept = [r for r in records if (r.score - minimum) * HARTREE_TO_KCAL <= window_kcal]
        return kept if kept else [min(records, key=lambda r: r.score)]

    def _apply_diversity_cap(
        self, records: List[DedupCandidate], max_candidates: int
    ) -> List[DedupCandidate]:
        if len(records) <= max_candidates:
            return list(records)

        # Do not group by integer torsion bins: values on opposite sides of a
        # bin edge can still be circularly equivalent. Reuse the same circular
        # distance rule as the geometric deduplicator.
        groups: List[List[DedupCandidate]] = []
        for record in sorted(records, key=lambda candidate: candidate.score):
            matching_group = next(
                (
                    group
                    for group in groups
                    if self.deduplicator.torsion_tolerance_deg >= 0
                    and signatures_equivalent(
                        record.signature,
                        group[0].signature,
                        self.deduplicator.torsion_tolerance_deg,
                    )
                ),
                None,
            )
            if matching_group is None:
                groups.append([record])
            else:
                matching_group.append(record)
        for group in groups:
            group.sort(key=lambda c: c.score)

        if len(groups) >= max_candidates:
            sorted_groups = sorted(groups, key=lambda g: g[0].score)
            return [g[0] for g in sorted_groups[:max_candidates]]

        selected: List[DedupCandidate] = [g[0] for g in groups]
        remaining_slots = max_candidates - len(selected)

        extras: List[DedupCandidate] = []
        for group in groups:
            extras.extend(group[1:])
        extras.sort(key=lambda c: c.score)
        selected.extend(extras[:remaining_slots])
        selected.sort(key=lambda c: c.score)
        return selected

    def _run_b97_sp_parallel(
        self, mol: Chem.Mol, records: List[DedupCandidate]
    ) -> List[DedupCandidate]:
        del mol
        b97_cfg = self._b97_cfg()
        parallel_jobs, cores_per_job = self._resolve_parallel_layout(b97_cfg, 1)
        total_cores = max(1, int((self.config.get("resources", {}) or {}).get("nproc", 1)))
        parallel_jobs = self._limit_orca_jobs_by_memory(
            b97_cfg,
            requested_jobs=parallel_jobs,
            cores_per_job=cores_per_job,
        )
        wave_plan = self._resolve_wave_plan(
            b97_cfg,
            total_cores=total_cores,
            parallel_jobs=parallel_jobs,
            base_cores_per_job=cores_per_job,
            candidate_count=len(records),
        )
        logger.debug(
            "B97-3c ensemble fan-out: jobs=%d cores_per_job=%d candidates=%d "
            "policy=%s final_wave=%s",
            parallel_jobs,
            cores_per_job,
            len(records),
            wave_plan.policy,
            list(wave_plan.final_cores),
        )
        batch_id = f"{self.molecule_name}:b97_sp"
        started = time.monotonic()
        self._emit(
            "batch_started",
            batch=batch_id,
            label="B97-3c SP",
            total=len(records),
            running=min(parallel_jobs, len(records)),
            cores_per_job=cores_per_job,
            total_cores=total_cores,
            scheduling_policy=wave_plan.policy,
            final_wave_cores=list(wave_plan.final_cores),
        )

        def _sp_worker(
            candidate: DedupCandidate,
            job_cores: int,
            job_maxcore: int | None,
            scheduler_layout: str,
            submitted_at: float,
        ) -> DedupCandidate:
            job_started = time.time()
            job_id = candidate.path.stem
            output = getattr(self.engine, "ranking_dir", self.engine.molecule_dir / "ranking") / job_id
            self._emit(
                "batch_job_started", batch=batch_id, job_id=job_id,
                engine="orca", method="B97-3c", nprocs=job_cores,
                maxcore_mb_per_rank=job_maxcore,
                scheduler_layout=scheduler_layout,
                submitted_at=submitted_at,
                queue_wait_seconds=max(0.0, job_started - submitted_at),
                started_at=job_started, input_xyz=str(candidate.path), output=str(output),
            )
            try:
                energy = self.engine.run_sp(
                    candidate.path,
                    nprocs=job_cores,
                    maxcore=job_maxcore,
                )
            except Exception as exc:
                self._emit(
                    "batch_job_failed", batch=batch_id, job_id=job_id,
                    elapsed_seconds=time.time() - job_started, error=str(exc), output=str(output),
                )
                raise
            metadata = dict(candidate.metadata)
            job_finished = time.time()
            metadata.update(
                {
                    "sp_scheduler_layout": scheduler_layout,
                    "sp_submitted_at_epoch": submitted_at,
                    "sp_started_at_epoch": job_started,
                    "sp_finished_at_epoch": job_finished,
                    "sp_queue_wait_seconds": max(0.0, job_started - submitted_at),
                    "sp_wall_seconds": max(0.0, job_finished - job_started),
                    "sp_nprocs": job_cores,
                    "sp_maxcore_mb_per_rank": job_maxcore,
                }
            )
            metadata.update(self._read_orca_task_metrics(output))
            if energy is not None:
                metadata["b973c_sp_energy_hartree"] = float(energy)
                metadata["b973c_electronic_cpcm_hartree"] = float(energy)
                ranking_cfg = dict(self.lite_config.get("ranking", {}) or {})
                ranking_solvent = ranking_cfg.get("solvent")
                metadata["b973c_electronic_model"] = (
                    f"B97-3c/CPCM({ranking_solvent})"
                    if ranking_solvent
                    else "B97-3c"
                )
                metadata["sp_status"] = "ok"
                self._emit(
                    "batch_job_finished", batch=batch_id, job_id=job_id,
                    elapsed_seconds=job_finished - job_started,
                    energy_hartree=float(energy), output=str(output), status="complete",
                )
                return DedupCandidate(candidate.path, float(energy), candidate.signature, metadata)
            metadata["sp_status"] = "failed"
            self._emit(
                "batch_job_failed", batch=batch_id, job_id=job_id,
                elapsed_seconds=time.time() - job_started,
                error="ORCA B97-3c SP did not converge", output=str(output),
            )
            return DedupCandidate(candidate.path, candidate.score, candidate.signature, metadata)

        def _progress(
            done: int,
            result: DedupCandidate,
            failed: int,
            running: int,
            active_cores: int,
        ) -> None:
            elapsed = max(time.monotonic() - started, 1e-9)
            rate = done / elapsed * 60.0
            self._emit(
                "batch_progress",
                batch=batch_id,
                label="B97-3c SP",
                total=len(records),
                done=done,
                failed=failed,
                running=running,
                active_cores=active_cores,
                total_cores=total_cores,
                current=result.path.stem,
                energy_hartree=result.metadata.get("b973c_sp_energy_hartree"),
                elapsed_seconds=elapsed,
                rate_per_minute=rate,
                eta_seconds=(len(records) - done) / max(rate / 60.0, 1e-9),
            )

        materialized, failures = self._execute_wave_plan(
            records,
            wave_plan=wave_plan,
            parallel_jobs=parallel_jobs,
            base_cores_per_job=cores_per_job,
            total_cores=total_cores,
            stage_cfg=b97_cfg,
            worker=_sp_worker,
            is_failure=lambda result: result.metadata.get("sp_status") == "failed",
            progress=_progress,
            memory_managed=True,
            stage_label="B97-3c",
        )
        self._emit(
            "batch_finished", batch=batch_id, label="B97-3c SP",
            total=len(records), done=len(materialized), failed=failures,
            status="complete" if failures == 0 else "degraded"
        )
        return materialized

    def _apply_mrrho_uniform(
        self, records: List[DedupCandidate], xtb_thermo_cfg: Dict[str, Any]
    ) -> tuple[List[DedupCandidate], str, str]:
        failure_policy = str(xtb_thermo_cfg.get("failure_policy", "defer")).strip().lower()
        if failure_policy == "ensemble_electronic_only":
            logger.warning(
                "Deprecated xTB mRRHO failure_policy=ensemble_electronic_only; "
                "using defer semantics without ensemble thermodynamics"
            )
            failure_policy = "defer"
        if failure_policy not in {"strict", "defer"}:
            raise ValueError(
                "step1.censo_lite.xtb_thermo.failure_policy must be 'strict' or 'defer'"
            )
        parallel_jobs, cores_per_job = self._resolve_parallel_layout(xtb_thermo_cfg, 1)
        total_cores = max(
            1, int((self.config.get("resources", {}) or {}).get("nproc", 1))
        )
        wave_plan = self._resolve_wave_plan(
            xtb_thermo_cfg,
            total_cores=total_cores,
            parallel_jobs=parallel_jobs,
            base_cores_per_job=cores_per_job,
            candidate_count=len(records),
        )
        logger.debug(
            "xTB mRRHO ensemble fan-out: jobs=%d cores_per_job=%d candidates=%d "
            "policy=%s final_wave=%s",
            parallel_jobs,
            cores_per_job,
            len(records),
            wave_plan.policy,
            list(wave_plan.final_cores),
        )
        batch_id = f"{self.molecule_name}:mrrho"
        started = time.monotonic()
        self._emit(
            "batch_started",
            batch=batch_id,
            label="xTB mRRHO",
            total=len(records),
            running=min(parallel_jobs, len(records)),
            cores_per_job=cores_per_job,
            total_cores=total_cores,
            scheduling_policy=wave_plan.policy,
            final_wave_cores=list(wave_plan.final_cores),
        )

        def _mrrho_worker(
            candidate: DedupCandidate,
            job_cores: int,
            _job_maxcore: int | None,
            scheduler_layout: str,
            submitted_at: float,
        ) -> DedupCandidate:
            job_started = time.time()
            job_id = candidate.path.stem
            output = getattr(self.engine, "mrrho_dir", self.engine.molecule_dir / "mrrho") / job_id
            self._emit(
                "batch_job_started", batch=batch_id, job_id=job_id,
                engine="xtb", method=f"GFN{xtb_thermo_cfg.get('gfn_level', 2)}-xTB mRRHO",
                nprocs=job_cores, started_at=job_started,
                scheduler_layout=scheduler_layout,
                submitted_at=submitted_at,
                queue_wait_seconds=max(0.0, job_started - submitted_at),
                input_xyz=str(candidate.path), output=str(output), attempt=1,
            )
            try:
                thermo_result = self.engine.run_mrrho(
                    candidate.path, nprocs=job_cores
                )
            except Exception as exc:
                self._emit(
                    "batch_job_failed", batch=batch_id, job_id=job_id,
                    elapsed_seconds=time.time() - job_started, error=str(exc), output=str(output),
                )
                raise
            metadata = dict(candidate.metadata)
            job_finished = time.time()
            metadata.update(
                {
                    "mrrho_scheduler_layout": scheduler_layout,
                    "mrrho_submitted_at_epoch": submitted_at,
                    "mrrho_started_at_epoch": job_started,
                    "mrrho_finished_at_epoch": job_finished,
                    "mrrho_queue_wait_seconds": max(0.0, job_started - submitted_at),
                    "mrrho_wall_seconds": max(0.0, job_finished - job_started),
                    "mrrho_nprocs": job_cores,
                }
            )
            # Runtime results carry the complete validated energy ledger. A
            # numeric value remains accepted only for lightweight test/runtime
            # adapters and is explicitly marked as unverified legacy input.
            if isinstance(thermo_result, (int, float)):
                correction = float(thermo_result)
                success = True
                ledger = {
                    "xtb_energy_ledger_status": "legacy_unverified",
                    "xtb_electronic_alpb_hartree": None,
                    "xtb_total_free_energy_alpb_hartree": None,
                    "xtb_energy_ledger_residual_hartree": None,
                    "xtb_solvation_hartree": None,
                    "xtb_mrrho_gfn_level": xtb_thermo_cfg.get("gfn_level", 2),
                    "xtb_mrrho_solvent_model": (
                        f"ALPB({xtb_thermo_cfg.get('solvent')})"
                        if xtb_thermo_cfg.get("solvent")
                        else None
                    ),
                    "explicit_xtb_solvation_correction_added": False,
                }
                error = None
            else:
                success = bool(
                    thermo_result is not None
                    and getattr(thermo_result, "success", False)
                    and getattr(thermo_result, "g_rrho_correction_hartree", None) is not None
                )
                thermo_correction = getattr(
                    thermo_result, "g_rrho_correction_hartree", None
                )
                correction = (
                    float(thermo_correction)
                    if success and thermo_correction is not None
                    else None
                )
                ledger = {
                    "xtb_energy_ledger_status": "validated" if success else "failed",
                    "xtb_electronic_alpb_hartree": getattr(
                        thermo_result, "xtb_electronic_energy_hartree", None
                    ),
                    "xtb_total_free_energy_alpb_hartree": getattr(
                        thermo_result, "xtb_total_free_energy_hartree", None
                    ),
                    "xtb_energy_ledger_residual_hartree": getattr(
                        thermo_result, "energy_ledger_residual_hartree", None
                    ),
                    # No independent ALPB solvation term is added to the CPCM
                    # electronic energy. The ALPB label describes the Hessian
                    # environment only.
                    "xtb_solvation_hartree": getattr(
                        thermo_result, "xtb_solvation_energy_hartree", None
                    ),
                    "xtb_mrrho_gfn_level": getattr(
                        thermo_result, "gfn_level", xtb_thermo_cfg.get("gfn_level", 2)
                    ),
                    "xtb_mrrho_solvent_model": getattr(
                        thermo_result, "solvent_model", None
                    ),
                    "explicit_xtb_solvation_correction_added": bool(
                        getattr(
                            thermo_result,
                            "explicit_solvation_correction_added",
                            False,
                        )
                    ),
                }
                error = getattr(thermo_result, "error", None) if thermo_result is not None else None
            metadata.update(ledger)
            if success and correction is not None:
                metadata["xtb_mrrho_thermal_correction_hartree"] = correction
                # Deprecated v2 alias.
                metadata["g_rrho_correction_hartree"] = correction
                metadata["mrrho_status"] = "ok"
                b97_energy = metadata.get("b973c_sp_energy_hartree")
                if b97_energy is None:
                    raise RuntimeError(f"Missing B97-3c energy for {candidate.path}")
                score = float(b97_energy) + correction
                self._emit(
                    "batch_job_finished", batch=batch_id, job_id=job_id,
                    elapsed_seconds=time.time() - job_started,
                    energy_hartree=float(correction), output=str(output), status="complete",
                )
            else:
                metadata["mrrho_status"] = "failed"
                metadata["mrrho_error"] = error or "xTB mRRHO failed after configured attempts"
                metadata["xtb_mrrho_thermal_correction_hartree"] = None
                metadata["g_rrho_correction_hartree"] = None
                b97_energy = metadata.get("b973c_sp_energy_hartree")
                if b97_energy is None:
                    b97_energy = candidate.score
                score = float(b97_energy)
                self._emit(
                    "batch_job_failed", batch=batch_id, job_id=job_id,
                    elapsed_seconds=time.time() - job_started,
                    error=metadata["mrrho_error"], output=str(output),
                )
            return DedupCandidate(candidate.path, score, candidate.signature, metadata)

        def _mrrho_progress(
            done: int,
            result: DedupCandidate,
            failed: int,
            running: int,
            active_cores: int,
        ) -> None:
            elapsed = max(time.monotonic() - started, 1e-9)
            rate = done / elapsed * 60.0
            self._emit(
                "batch_progress",
                batch=batch_id,
                label="xTB mRRHO",
                total=len(records),
                done=done,
                failed=failed,
                running=running,
                active_cores=active_cores,
                total_cores=total_cores,
                current=result.path.stem,
                energy_hartree=result.metadata.get("g_rrho_correction_hartree"),
                elapsed_seconds=elapsed,
                rate_per_minute=rate,
                eta_seconds=(len(records) - done) / max(rate / 60.0, 1e-9),
            )

        corrected, failures = self._execute_wave_plan(
            records,
            wave_plan=wave_plan,
            parallel_jobs=parallel_jobs,
            base_cores_per_job=cores_per_job,
            total_cores=total_cores,
            stage_cfg=xtb_thermo_cfg,
            worker=_mrrho_worker,
            is_failure=lambda result: result.metadata.get("mrrho_status") != "ok",
            progress=_mrrho_progress,
            memory_managed=False,
            stage_label="xTB mRRHO",
        )

        self._emit(
            "batch_finished", batch=batch_id, label="xTB mRRHO",
            total=len(records), done=len(corrected), failed=failures,
            status="complete" if failures == 0 else "degraded"
        )

        any_failure = any(
            candidate.metadata.get("mrrho_status") != "ok" for candidate in corrected
        )
        method_stacks = {
            (
                candidate.metadata.get("xtb_mrrho_gfn_level"),
                candidate.metadata.get("xtb_mrrho_solvent_model"),
            )
            for candidate in corrected
            if candidate.metadata.get("mrrho_status") == "ok"
        }
        method_mismatch = len(method_stacks) > 1
        any_failure = any_failure or method_mismatch
        scoring_mode = "b97_3c_plus_mrrho"
        thermochemistry_status = "complete"
        if any_failure:
            if failure_policy == "strict":
                raise RuntimeError(
                    "xTB mRRHO failed or produced an inconsistent method stack "
                    "and failure_policy is 'strict'"
                )
            logger.warning(
                "xTB mRRHO had partial failures; retaining B97-3c electronic "
                "ranking for geometry handoff while deferring ensemble thermodynamics"
            )
            deferred_records: List[DedupCandidate] = []
            for c in corrected:
                b97 = c.metadata.get("b973c_sp_energy_hartree")
                if b97 is not None:
                    metadata = dict(c.metadata)
                    metadata["thermochemistry_usable"] = False
                    metadata["thermochemistry_incomplete_reason"] = (
                        "inconsistent_mrrho_method_stack"
                        if method_mismatch
                        else "partial_mrrho_failure"
                    )
                    deferred_records.append(
                        DedupCandidate(c.path, float(b97), c.signature, metadata)
                    )
            corrected = deferred_records
            scoring_mode = "electronic_only_deferred"
            thermochemistry_status = "incomplete"
            self._emit(
                "decision",
                decision="thermochemistry_deferred",
                severity="warning",
                message=(
                    "inconsistent xTB mRRHO method stacks"
                    if method_mismatch
                    else (
                        f"{sum(c.metadata.get('mrrho_status') != 'ok' for c in corrected)} "
                        "mRRHO jobs failed"
                    )
                ),
                previous="E(B97-3c) + G(RRHO)xTB",
                active="provisional E(B97-3c) geometry ranking; no Zrel/population/Gconf",
                reason="partial thermochemistry cannot define a uniform partition function",
            )
        else:
            for candidate in corrected:
                candidate.metadata["thermochemistry_usable"] = True

        return corrected, scoring_mode, thermochemistry_status

    def _execute_wave_plan(
        self,
        records: Sequence[DedupCandidate],
        *,
        wave_plan: _WavePlan,
        parallel_jobs: int,
        base_cores_per_job: int,
        total_cores: int,
        stage_cfg: Dict[str, Any],
        worker: Callable[[DedupCandidate, int, int | None, str, float], DedupCandidate],
        is_failure: Callable[[DedupCandidate], bool],
        progress: Callable[[int, DedupCandidate, int, int, int], None],
        memory_managed: bool,
        stage_label: str,
    ) -> tuple[List[DedupCandidate], int]:
        """Execute a throughput bulk phase and an optional profiled final wave.

        A final wave is a real barrier: no widened task is launched while a
        bulk task is still running. This avoids the former static ``last N``
        queue split that mixed one- and two-rank ORCA jobs at batch start.
        """

        results: List[Optional[DedupCandidate]] = [None] * len(records)
        failures = 0
        done = 0

        phases: List[tuple[str, List[tuple[int, DedupCandidate, int]]]] = []
        if wave_plan.bulk_count:
            phases.append(
                (
                    "throughput_bulk",
                    [
                        (index, records[index], base_cores_per_job)
                        for index in range(wave_plan.bulk_count)
                    ],
                )
            )
        if wave_plan.final_cores:
            start = wave_plan.bulk_count
            phases.append(
                (
                    "profiled_final_wave",
                    [
                        (start + offset, records[start + offset], cores)
                        for offset, cores in enumerate(wave_plan.final_cores)
                    ],
                )
            )

        for layout_name, phase_jobs in phases:
            if not phase_jobs:
                continue
            phase_limit = min(parallel_jobs, len(phase_jobs))
            if layout_name == "profiled_final_wave":
                phase_limit = len(phase_jobs)
            phase_cores = [cores for _index, _candidate, cores in phase_jobs]
            phase_maxcore = (
                self._orca_wave_maxcore(stage_cfg, phase_cores[:phase_limit])
                if memory_managed
                else None
            )
            pending = list(phase_jobs)
            active_cores = 0
            futures: Dict[Future[DedupCandidate], tuple[int, int]] = {}

            with ThreadPoolExecutor(max_workers=max(1, phase_limit)) as pool:
                while pending or futures:
                    submitted = False
                    while pending and len(futures) < phase_limit:
                        index, candidate, job_cores = pending[0]
                        if active_cores + job_cores > total_cores:
                            break
                        pending.pop(0)
                        submitted_at = time.time()
                        future = pool.submit(
                            worker,
                            candidate,
                            job_cores,
                            phase_maxcore,
                            layout_name,
                            submitted_at,
                        )
                        futures[future] = (index, job_cores)
                        active_cores += job_cores
                        submitted = True

                    if not futures:
                        if pending:
                            raise RuntimeError(
                                f"{stage_label} wave scheduler stalled with pending jobs"
                            )
                        break

                    completed, _ = wait(tuple(futures), return_when=FIRST_COMPLETED)
                    for future in completed:
                        index, job_cores = futures.pop(future)
                        active_cores -= job_cores
                        result = future.result()
                        results[index] = result
                        done += 1
                        failures += int(is_failure(result))
                        progress(done, result, failures, len(futures), active_cores)

                    if not submitted and not completed and pending:
                        raise RuntimeError(
                            f"{stage_label} wave scheduler made no progress"
                        )

        return [candidate for candidate in results if candidate is not None], failures

    @staticmethod
    def _resolve_wave_plan(
        stage_cfg: Dict[str, Any],
        *,
        total_cores: int,
        parallel_jobs: int,
        base_cores_per_job: int,
        candidate_count: int,
    ) -> _WavePlan:
        """Resolve an exact, benchmark-provided final-wave layout.

        No heuristic widening is performed. Without an exact profile entry the
        complete batch stays in throughput mode, even if that leaves cores
        idle for a small ensemble.
        """

        if candidate_count <= 0:
            return _WavePlan(0, (), "throughput_first")

        scheduling = dict(stage_cfg.get("scheduling", {}) or {})
        final_cfg = dict(scheduling.get("final_wave", {}) or {})
        legacy = stage_cfg.get("queue_drain") or stage_cfg.get("adaptive_tail")
        if legacy and bool(dict(legacy).get("enabled", False)):
            logger.warning(
                "Ignoring deprecated queue_drain/adaptive_tail configuration; "
                "use scheduling.final_wave.profile_layouts after benchmarking"
            )
        if not bool(final_cfg.get("enabled", False)):
            return _WavePlan(candidate_count, (), "throughput_first")

        wave_capacity = max(1, parallel_jobs)
        remainder = candidate_count % wave_capacity
        final_count = remainder if remainder else min(candidate_count, wave_capacity)
        profiles = dict(final_cfg.get("profile_layouts", {}) or {})
        raw_layout = profiles.get(str(final_count), profiles.get(final_count))
        if raw_layout is None:
            return _WavePlan(candidate_count, (), "throughput_first_profile_fallback")

        layout = tuple(int(value) for value in raw_layout)
        max_cores = max(
            base_cores_per_job,
            int(final_cfg.get("max_cores_per_job", total_cores)),
        )
        valid = (
            len(layout) == final_count
            and all(base_cores_per_job <= value <= max_cores for value in layout)
            and sum(layout) <= total_cores
            and len(layout) <= parallel_jobs
        )
        if not valid:
            logger.warning(
                "Ignoring invalid final-wave profile for %d jobs: %s",
                final_count,
                list(layout),
            )
            return _WavePlan(candidate_count, (), "throughput_first_profile_invalid")

        if all(value == base_cores_per_job for value in layout):
            return _WavePlan(candidate_count, (), "throughput_first_profile_noop")
        return _WavePlan(
            candidate_count - final_count,
            layout,
            "throughput_bulk_with_profiled_final_wave",
        )

    def _limit_orca_jobs_by_memory(
        self,
        stage_cfg: Dict[str, Any],
        *,
        requested_jobs: int,
        cores_per_job: int,
    ) -> int:
        scheduling = dict(stage_cfg.get("scheduling", {}) or {})
        memory_cfg = dict(scheduling.get("memory", {}) or {})
        resources = dict(self.config.get("resources", {}) or {})
        total_mb = mem_to_mb(str(resources.get("mem", "32GB")))
        safety = float(resources.get("orca_maxcore_safety", 0.65))
        usable_mb = max(1, int(total_mb * safety))
        overhead_mb = max(0, int(memory_cfg.get("job_overhead_mb", 256)))
        minimum_mb = max(1, int(memory_cfg.get("min_maxcore_mb_per_rank", 1000)))

        for jobs in range(max(1, requested_jobs), 0, -1):
            required = jobs * overhead_mb + jobs * cores_per_job * minimum_mb
            if required <= usable_mb:
                if jobs != requested_jobs:
                    logger.warning(
                        "Reduced ORCA fan-out from %d to %d jobs to satisfy the "
                        "configured memory budget",
                        requested_jobs,
                        jobs,
                    )
                return jobs
        return 1

    def _orca_wave_maxcore(
        self, stage_cfg: Dict[str, Any], job_cores: Sequence[int]
    ) -> int:
        scheduling = dict(stage_cfg.get("scheduling", {}) or {})
        memory_cfg = dict(scheduling.get("memory", {}) or {})
        resources = dict(self.config.get("resources", {}) or {})
        total_mb = mem_to_mb(str(resources.get("mem", "32GB")))
        safety = float(resources.get("orca_maxcore_safety", 0.65))
        usable_mb = max(1, int(total_mb * safety))
        overhead_mb = max(0, int(memory_cfg.get("job_overhead_mb", 256)))
        minimum_mb = max(1, int(memory_cfg.get("min_maxcore_mb_per_rank", 1000)))
        ranks = max(1, sum(int(value) for value in job_cores))
        available = usable_mb - len(job_cores) * overhead_mb
        maxcore = available // ranks
        if maxcore < minimum_mb:
            raise RuntimeError(
                "ORCA wave does not fit the configured memory budget: "
                f"jobs={len(job_cores)} ranks={ranks} maxcore={maxcore} MB"
            )
        configured_cap = memory_cfg.get("maxcore_cap_mb_per_rank")
        if configured_cap is not None:
            maxcore = min(maxcore, max(1, int(configured_cap)))
        return int(maxcore)

    @staticmethod
    def _read_orca_task_metrics(output_dir: Path) -> Dict[str, Any]:
        """Read lightweight timing/convergence metrics from the latest ORCA output."""

        outputs = sorted(
            Path(output_dir).glob("*.out"),
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )
        if not outputs:
            return {}
        try:
            text = outputs[0].read_text(encoding="utf-8", errors="ignore")
        except OSError:
            return {}
        metrics: Dict[str, Any] = {"sp_output": str(outputs[0])}
        cycles = re.findall(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", text)
        if cycles:
            metrics["sp_scf_cycles"] = int(cycles[-1])
        runtimes = re.findall(
            r"TOTAL RUN TIME:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+"
            r"(\d+) seconds\s+(\d+) msec",
            text,
        )
        if runtimes:
            days, hours, minutes, seconds, milliseconds = map(int, runtimes[-1])
            metrics["sp_orca_total_run_seconds"] = (
                days * 86400
                + hours * 3600
                + minutes * 60
                + seconds
                + milliseconds / 1000.0
            )
        return metrics

    @staticmethod
    def _summarize_runtime_metrics(
        rows: Sequence[Dict[str, Any]], prefix: str
    ) -> Dict[str, Any]:
        durations = sorted(
            float(row[f"{prefix}_wall_seconds"])
            for row in rows
            if row.get(f"{prefix}_wall_seconds") is not None
        )
        layouts: Dict[str, int] = {}
        core_counts: Dict[str, int] = {}
        total_core_seconds = 0.0
        for row in rows:
            layout = row.get(f"{prefix}_scheduler_layout")
            if layout is not None:
                layouts[str(layout)] = layouts.get(str(layout), 0) + 1
            nprocs = row.get(f"{prefix}_nprocs")
            if nprocs is not None:
                core_counts[str(int(nprocs))] = core_counts.get(str(int(nprocs)), 0) + 1
                wall = row.get(f"{prefix}_wall_seconds")
                if wall is not None:
                    total_core_seconds += float(wall) * int(nprocs)
        if not durations:
            return {
                "jobs": 0,
                "layouts": layouts,
                "cores_per_job_counts": core_counts,
            }
        p95_index = min(len(durations) - 1, max(0, (95 * len(durations) + 99) // 100 - 1))
        return {
            "jobs": len(durations),
            "layouts": layouts,
            "cores_per_job_counts": core_counts,
            "p50_wall_seconds": durations[(len(durations) - 1) // 2],
            "p95_wall_seconds": durations[p95_index],
            "total_core_seconds": total_core_seconds,
        }

    def _resolve_parallel_layout(
        self, stage_cfg: Dict[str, Any], default_cores_per_job: int
    ) -> tuple[int, int]:
        """Fit process fan-out into the global S1 core budget."""

        resources = dict(self.config.get("resources", {}) or {})
        total_cores = max(1, int(resources.get("nproc", 1)))
        cores_per_job = max(
            1,
            min(total_cores, int(stage_cfg.get("cores_per_job") or default_cores_per_job)),
        )
        capacity = max(1, total_cores // cores_per_job)
        requested = stage_cfg.get("parallel_jobs", "auto")
        if requested is None or str(requested).strip().lower() == "auto":
            parallel_jobs = capacity
        else:
            parallel_jobs = max(1, min(int(requested), capacity))
        max_jobs = stage_cfg.get("max_parallel_jobs")
        if max_jobs is not None:
            parallel_jobs = min(parallel_jobs, max(1, int(max_jobs)))
        return parallel_jobs, cores_per_job

    def _annotate_ensemble_thermodynamics(
        self, records: List[DedupCandidate], thermochemistry_status: str
    ) -> tuple[List[DedupCandidate], Dict[str, Any]]:
        temperature = float(
            (self.lite_config.get("xtb_thermo", {}) or {}).get(
                "temperature_k",
                (self.config.get("thermo", {}) or {}).get("temperature_k", 298.15),
            )
        )
        failed_ids = [
            str(record.metadata.get("source_conformer_id") or record.path.stem)
            for record in records
            if record.metadata.get("mrrho_status") == "failed"
        ]
        incomplete_reasons = sorted(
            {
                str(record.metadata["thermochemistry_incomplete_reason"])
                for record in records
                if record.metadata.get("thermochemistry_incomplete_reason")
            }
        )
        if thermochemistry_status != "complete":
            minimum = min(record.score for record in records)
            annotated_records: List[DedupCandidate] = []
            for record in records:
                metadata = dict(record.metadata)
                metadata.update(
                    {
                        "degeneracy": int(metadata.get("degeneracy", 1)),
                        "relative_s1_score_kcal": (
                            record.score - minimum
                        ) * HARTREE_TO_KCAL,
                        "relative_free_energy_kcal": None,
                        "relative_partition_weight": None,
                        "boltzmann_population": None,
                    }
                )
                annotated_records.append(
                    DedupCandidate(record.path, record.score, record.signature, metadata)
                )
            return annotated_records, {
                "status": thermochemistry_status,
                "temperature_k": temperature,
                "ensemble_member_count": len(records),
                "partition_function_relative": None,
                "conformational_free_energy_correction_kcal": None,
                "reference_free_energy_hartree": None,
                "ensemble_free_energy_estimate_hartree": None,
                "population_sum": None,
                "failed_conformer_ids": failed_ids,
                "incomplete_reasons": incomplete_reasons,
                "required_downstream_action": (
                    "frequency_completion_or_review"
                    if thermochemistry_status == "incomplete"
                    else None
                ),
            }

        degeneracies = [int(record.metadata.get("degeneracy", 1)) for record in records]
        result = calculate_ensemble_thermodynamics(
            [record.score for record in records],
            temperature,
            degeneracies,
        )
        annotated: List[DedupCandidate] = []
        for index, record in enumerate(records):
            metadata = dict(record.metadata)
            metadata.update(
                {
                    "degeneracy": degeneracies[index],
                    "relative_s1_score_kcal": result.relative_free_energies_kcal[index],
                    "relative_free_energy_kcal": result.relative_free_energies_kcal[index],
                    "relative_partition_weight": result.relative_partition_weights[index],
                    "boltzmann_population": result.boltzmann_populations[index],
                }
            )
            annotated.append(DedupCandidate(record.path, record.score, record.signature, metadata))
        summary = {
            "status": "complete",
            "temperature_k": result.temperature_k,
            "ensemble_member_count": len(records),
            "partition_function_relative": result.partition_function_relative,
            "conformational_free_energy_correction_kcal": (
                result.conformational_free_energy_correction_kcal
            ),
            "reference_free_energy_hartree": result.reference_free_energy_hartree,
            "ensemble_free_energy_estimate_hartree": result.ensemble_free_energy_hartree,
            "population_sum": sum(result.boltzmann_populations),
            "failed_conformer_ids": [],
            "incomplete_reasons": [],
            "required_downstream_action": None,
        }
        return annotated, summary

    def _apply_final_selection(
        self,
        records: List[DedupCandidate],
        final_window_kcal: float,
        min_keep: int,
        max_keep: int,
    ) -> List[DedupCandidate]:
        if not records:
            return records
        minimum = min(r.score for r in records)
        in_window = [r for r in records if (r.score - minimum) * HARTREE_TO_KCAL <= final_window_kcal]
        if len(in_window) < min_keep:
            in_window = list(records)
        return in_window[:max_keep]

    def _write_manifest(
        self,
        ensemble: List[DedupCandidate],
        representatives: List[DedupCandidate],
        output_dir: Path,
        pipeline_log: List[Dict[str, Any]],
        scoring_mode: str,
        ensemble_thermo: Dict[str, Any],
    ) -> Dict[str, Any]:
        if not ensemble:
            raise RuntimeError("CENSO-LITE produced no valid candidates after filtering")
        ensemble.sort(key=lambda item: (item.score, str(item.path)))
        minimum = ensemble[0].score

        candidate_dir = output_dir / "candidates"
        candidate_dir.mkdir(parents=True, exist_ok=True)
        manifest_candidates: List[Dict[str, Any]] = []
        id_by_source: Dict[str, str] = {}
        for index, candidate in enumerate(ensemble, start=1):
            destination = candidate_dir / f"conf_{index:04d}.xyz"
            if candidate.path.resolve() != destination.resolve():
                destination.write_bytes(candidate.path.read_bytes())
            payload = dict(candidate.metadata)
            candidate_id = f"{self.molecule_name}_conf_{index:04d}"
            source_id = str(payload.get("source_conformer_id") or candidate.path.stem)
            merged_from = [str(value) for value in (payload.get("merged_from") or [source_id])]
            payload.update(
                {
                    "id": candidate_id,
                    "xyz": str(destination.relative_to(self.work_dir)),
                    # Backward-compatible alias; the value is now explicitly a
                    # relative conformer free energy when mRRHO is active.
                    "relative_energy_kcal": (candidate.score - minimum) * HARTREE_TO_KCAL,
                    "relative_s1_score_kcal": (candidate.score - minimum) * HARTREE_TO_KCAL,
                    "s1_score_hartree": candidate.score,
                    "source": "crest_gfn2_censo_light_ranking",
                    "source_conformer_id": source_id,
                    "merged_from": merged_from,
                    "merge_count": len(merged_from),
                    "degeneracy": int(payload.get("degeneracy", 1)),
                    "degeneracy_source": payload.get(
                        "degeneracy_source", "default_unique_minimum"
                    ),
                    "status": "valid",
                }
            )
            manifest_candidates.append(payload)
            id_by_source[str(candidate.path.resolve())] = candidate_id

        ranking_formula = (
            "E_B97-3c_SP + G(RRHO)_xTB"
            if scoring_mode == "b97_3c_plus_mrrho"
            else "E_B97-3c_SP (provisional geometry ranking only)"
        )
        representative_ids = [
            id_by_source[str(candidate.path.resolve())]
            for candidate in representatives
            if str(candidate.path.resolve()) in id_by_source
        ]
        thermo_payload = dict(ensemble_thermo)
        thermochemistry_status = str(thermo_payload.get("status") or "not_requested")
        thermochemistry_complete = thermochemistry_status == "complete"
        reference_g_rrho = (
            manifest_candidates[0].get("xtb_mrrho_thermal_correction_hartree")
            if thermochemistry_complete
            else None
        )
        ensemble_correction_hartree = None
        if reference_g_rrho is not None:
            ensemble_correction_hartree = float(reference_g_rrho) + (
                float(thermo_payload["conformational_free_energy_correction_kcal"])
                / HARTREE_TO_KCAL
            )
        thermo_payload.update(
            {
                "reference_conformer_id": manifest_candidates[0]["id"],
                "reference_g_rrho_correction_hartree": reference_g_rrho,
                "ensemble_thermochemistry_correction_hartree": ensemble_correction_hartree,
                "thermochemistry_status": thermochemistry_status,
                "thermochemistry_complete": thermochemistry_complete,
                "free_energy_formula": ranking_formula,
                "ensemble_formula": "G_ensemble = G_reference + G_conf_rel",
                "partition_formula": "Z_rel = sum_i d_i exp(-DeltaG_i/RT)",
                "conformational_correction_formula": "G_conf_rel = -RT ln(Z_rel)",
            }
        )
        selection_cfg = self._selection_cfg()
        selected_id = manifest_candidates[0]["id"]
        thermodynamic_rank1 = selected_id if thermochemistry_complete else None
        manifest = {
            "schema_version": self.manifest_schema_version,
            "stage": "S1",
            RUN_ID_FIELD: self.run_id,
            "protocol_id": "censo_light_inspired",
            "protocol": {
                "family": "CENSO-light-inspired",
                "stage_scope": "ensemble_ranking_and_representative_selection",
                "complete_censo_light": False,
                "selected_geometry_level": "GFN2-xTB/ALPB",
                "geometry_is_final_minimum": False,
                "requires_downstream_optimization": True,
                "downstream_protocol": "RPH_V4_S3_S4",
                "literature_exact_reproduction": False,
                "deferred_steps": [
                    "representative geometry optimization",
                    "higher-level representative single point",
                    "final electronic-plus-ensemble free-energy assembly",
                ],
            },
            "molecule": self.molecule_name,
            "ranking_formula": ranking_formula,
            "scoring_mode": scoring_mode,
            "thermochemistry_status": thermochemistry_status,
            "final_window_kcal": float(selection_cfg.get("final_window_kcal", 2.5)),
            "pipeline_log": pipeline_log,
            "scheduling_summary": {
                "b97_3c": self._summarize_runtime_metrics(manifest_candidates, "sp"),
                "xtb_mrrho": self._summarize_runtime_metrics(
                    manifest_candidates, "mrrho"
                ),
            },
            "ensemble_thermodynamics": thermo_payload,
            "candidates": manifest_candidates,
            "thermodynamic_rank1": thermodynamic_rank1,
            "reactivity_screening_candidates": representative_ids,
            "provisional_geometry_selection": (
                selected_id if not thermochemistry_complete else None
            ),
            "selection_basis": (
                "complete_s1_free_energy"
                if thermochemistry_complete
                else "b97_3c_electronic_energy_deferred_thermochemistry"
            ),
            # Deprecated v2 aliases retained for one migration cycle.
            "representative_candidates": representative_ids,
            "selected": selected_id,
            "selected_xyz": "selected.xyz",
            "atom_mapping_ref": "atom_mapping.json",
            "atom_mapping_status": "verified",
            "deprecated_fields": ["selected", "representative_candidates"],
        }
        selected_candidate = manifest_candidates[0]
        selected_src = self.work_dir / str(selected_candidate["xyz"])
        selected_dst = output_dir / "selected.xyz"
        _ = shutil.copy2(selected_src, selected_dst)
        mapping_path = bind_atom_mapping_to_xyz(
            output_dir / "initial_atom_mapping.json",
            selected_dst,
            output_dir / "atom_mapping.json",
        )
        mapping_payload = json.loads(mapping_path.read_text(encoding="utf-8"))
        manifest["provenance"] = build_provenance(
            run_id=self.run_id,
            schema_version=str(manifest["schema_version"]),
            protocol_version=None,
            parent_manifest_paths={"s0": self.s0_mechanism_json_path},
            atom_mapping_payload=mapping_payload,
            forming_bonds=self.forming_bonds,
            variant_manifest_path=None,
            extra={"selected_xyz_hash": sha256_of_file(selected_dst)},
        ).to_dict()
        manifest_path = output_dir / "manifest.json"
        write_json_atomic(manifest_path, manifest)
        top_population = [
            {
                "id": row.get("id"),
                "population": row.get("boltzmann_population"),
                "relative_free_energy_kcal": row.get("relative_free_energy_kcal"),
            }
            for row in sorted(
                manifest_candidates,
                key=lambda item: float(item.get("boltzmann_population") or 0.0),
                reverse=True,
            )[: max(1, int((self.config.get("ui", {}) or {}).get("top_population_count", 5)))]
        ]
        self._emit(
            "science_summary",
            funnel={row["step"]: row.get("candidates") for row in pipeline_log if row.get("candidates") is not None},
            ensemble_members=len(manifest_candidates),
            representatives=len(representative_ids),
            ranking_formula=ranking_formula,
            scoring_mode=scoring_mode,
            final_window_kcal=manifest["final_window_kcal"],
            selected=manifest["selected"],
            top_population=top_population,
            conformational_free_energy_correction_kcal=thermo_payload.get(
                "conformational_free_energy_correction_kcal"
            ),
            thermochemistry_complete=thermo_payload.get("thermochemistry_complete"),
        )
        return {
            "manifest": manifest_path,
            "selected_xyz": selected_dst,
            "atom_mapping": mapping_path,
            "data": manifest,
        }
