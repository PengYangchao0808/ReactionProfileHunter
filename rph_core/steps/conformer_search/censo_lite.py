"""Fixed V4 CENSO-LITE conformer workflow.

The workflow deliberately stops before geometry-level DFT refinement.  It uses
CREST/GFN2 geometries, B97-3c single-point ranking and optional xTB mRRHO
corrections.  High-level refinement is owned by S4.

Pipeline order (revised):
    1.  RDKit embed → initial.xyz
    2.  CREST GFN2/ALPB → crest_conformers.xyz
    3.  Split ensemble → N conf_####.xyz
    4.  Parse xTB energies (fail-fast) + cross-validate with crest.energies
    5.  Torsion-aware dedup (geometric, cheap)
    6.  xTB relative-energy window ≤ prefilter.xtb_window_kcal
    7.  Diversity-preserving soft cap ≤ prefilter.max_candidates_soft
    8.  Parallel ORCA B97-3c/CPCM SP
    9.  B97-3c energy window ≤ b97_3c.energy_window_kcal
    10. Diversity cap → thermo pool (≤ xtb_thermo.thermo_pool_max)
    11. mRRHO on entire thermo pool (uniform — no mixed scoring)
    12. Final ranking: G = E_B97-3c + G(RRHO)_xTB
    13. Final window ≤ selection.final_window_kcal, min_keep/max_keep
"""

from __future__ import annotations

import copy
import json
import logging
import shutil
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional

from rdkit import Chem

from rph_core.steps.conformer_search.censo_lite_runtime import CensoLiteRuntime, CrestEnergyParseError
from rph_core.steps.conformer_search.deduplicator import DedupCandidate, TorsionAwareDeduplicator
from rph_core.utils.constants import HARTREE_TO_KCAL

logger = logging.getLogger(__name__)


class CensoLiteEngine:
    """Run the single supported S1 protocol and emit a stable manifest."""

    protocol_name = "censo_lite"

    def __init__(self, config: Dict[str, Any], work_dir: Path, molecule_name: str):
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

    def run(self, smiles: str) -> Dict[str, Any]:
        output_dir = self.engine.molecule_dir
        pipeline_log: List[Dict[str, Any]] = []

        initial_xyz = self.engine.embed(smiles)
        ensemble_xyz = self.engine.crest_search(initial_xyz)
        candidate_paths = self.engine.split_ensemble(ensemble_xyz)
        if not candidate_paths:
            raise RuntimeError("CENSO-LITE produced no CREST candidates")

        pipeline_log.append({"step": "crest", "candidates": len(candidate_paths)})

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES for CENSO-LITE: {smiles}")
        mol = Chem.AddHs(mol)

        # ---- Stage 4: parse xTB energies (fail-fast) + cross-validate ----
        records = self._parse_and_annotate(mol, candidate_paths)
        pipeline_log.append({"step": "energy_parse", "candidates": len(records)})

        # ---- Stage 5: torsion-aware dedup (geometric, cheap) ----
        extra_dedup = bool(self._prefilter_cfg().get("extra_deduplication", True))
        if extra_dedup:
            records = self.deduplicator.deduplicate(mol, records)
            pipeline_log.append({"step": "torsion_dedup", "candidates": len(records)})

        # ---- Stage 6: xTB relative-energy window ----
        xtb_window = float(self._prefilter_cfg().get("xtb_window_kcal", 5.0))
        records = self._apply_energy_window(records, xtb_window)
        pipeline_log.append({"step": "xtb_window", "window_kcal": xtb_window, "candidates": len(records)})

        # ---- Stage 7: diversity-preserving soft cap ----
        max_soft = int(self._prefilter_cfg().get("max_candidates_soft", 60))
        records = self._apply_diversity_cap(records, max_soft)
        pipeline_log.append({"step": "diversity_cap_pre_sp", "max": max_soft, "candidates": len(records)})

        # ---- Stage 8: parallel ORCA B97-3c SP ----
        records = self._run_b97_sp_parallel(mol, records)
        sp_failed = [r for r in records if r.metadata.get("sp_status") == "failed"]
        records = [r for r in records if r.metadata.get("sp_status") != "failed"]
        if not records:
            raise RuntimeError("All B97-3c SP calculations failed")
        pipeline_log.append({"step": "b97_sp", "candidates": len(records), "failed": len(sp_failed)})

        # ---- Stage 9: B97-3c energy window ----
        b97_window = float(self._b97_cfg().get("energy_window_kcal", 4.0))
        records = self._apply_energy_window(records, b97_window)
        pipeline_log.append({"step": "b97_window", "window_kcal": b97_window, "candidates": len(records)})

        # ---- Stage 10: diversity cap → thermo pool ----
        xtb_thermo_cfg = self.lite_config.get("xtb_thermo", {}) or {}
        thermo_pool_max = int(xtb_thermo_cfg.get("thermo_pool_max", 30))
        records = self._apply_diversity_cap(records, thermo_pool_max)
        pipeline_log.append({"step": "thermo_pool", "max": thermo_pool_max, "candidates": len(records)})

        # ---- Stage 11-12: mRRHO on entire thermo pool + final ranking ----
        thermo_enabled = bool(xtb_thermo_cfg.get("enabled", True))
        if thermo_enabled:
            records = self._apply_mrrho_uniform(records, xtb_thermo_cfg)
            pipeline_log.append({"step": "mrrho", "candidates": len(records)})

        records.sort(key=lambda item: (item.score, str(item.path)))

        # ---- Stage 13: final selection ----
        final_window = float(self._selection_cfg().get("final_window_kcal", 2.5))
        min_keep = int(self._selection_cfg().get("min_keep", 3))
        max_keep = int(self._selection_cfg().get("max_keep", 10))
        kept = self._apply_final_selection(records, final_window, min_keep, max_keep)
        pipeline_log.append({"step": "final_selection", "candidates": len(kept)})

        return self._write_manifest(kept, output_dir, pipeline_log, thermo_enabled)

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

        groups: Dict[str, List[DedupCandidate]] = defaultdict(list)
        for r in records:
            groups[r.signature.key()].append(r)
        for group in groups.values():
            group.sort(key=lambda c: c.score)

        if len(groups) >= max_candidates:
            sorted_groups = sorted(groups.values(), key=lambda g: g[0].score)
            return [g[0] for g in sorted_groups[:max_candidates]]

        selected: List[DedupCandidate] = [g[0] for g in groups.values()]
        remaining_slots = max_candidates - len(selected)

        extras: List[DedupCandidate] = []
        for group in groups.values():
            extras.extend(group[1:])
        extras.sort(key=lambda c: c.score)
        selected.extend(extras[:remaining_slots])
        selected.sort(key=lambda c: c.score)
        return selected

    def _run_b97_sp_parallel(
        self, mol: Chem.Mol, records: List[DedupCandidate]
    ) -> List[DedupCandidate]:
        b97_cfg = self._b97_cfg()
        resources = dict(self.config.get("resources", {}) or {})
        total_cores = int(resources.get("nproc", 1))
        parallel_jobs = int(b97_cfg.get("parallel_jobs", 1))
        cores_per_job = int(b97_cfg.get("cores_per_job", total_cores))
        if parallel_jobs * cores_per_job > total_cores:
            logger.warning(
                "parallel_jobs(%d) * cores_per_job(%d) > nproc(%d); "
                "reducing parallel_jobs to %d",
                parallel_jobs,
                cores_per_job,
                total_cores,
                max(1, total_cores // cores_per_job),
            )
            parallel_jobs = max(1, total_cores // cores_per_job)

        def _sp_worker(candidate: DedupCandidate) -> DedupCandidate:
            energy = self.engine.run_sp(candidate.path, nprocs=cores_per_job)
            metadata = dict(candidate.metadata)
            if energy is not None:
                metadata["b973c_sp_energy_hartree"] = float(energy)
                metadata["sp_status"] = "ok"
                return DedupCandidate(candidate.path, float(energy), candidate.signature, metadata)
            metadata["sp_status"] = "failed"
            return DedupCandidate(candidate.path, candidate.score, candidate.signature, metadata)

        if parallel_jobs <= 1:
            return [_sp_worker(r) for r in records]

        results: List[DedupCandidate] = [None] * len(records)
        with ThreadPoolExecutor(max_workers=parallel_jobs) as pool:
            future_to_idx = {pool.submit(_sp_worker, r): i for i, r in enumerate(records)}
            for future in as_completed(future_to_idx):
                results[future_to_idx[future]] = future.result()
        return [r for r in results if r is not None]

    def _apply_mrrho_uniform(
        self, records: List[DedupCandidate], xtb_thermo_cfg: Dict[str, Any]
    ) -> List[DedupCandidate]:
        failure_policy = str(xtb_thermo_cfg.get("failure_policy", "ensemble_electronic_only"))

        corrected: List[DedupCandidate] = []
        any_failure = False
        for candidate in records:
            correction = self.engine.run_mrrho(candidate.path)
            metadata = dict(candidate.metadata)
            if correction is not None:
                metadata["g_rrho_correction_hartree"] = float(correction)
                b97_energy = metadata.get("b973c_sp_energy_hartree") or 0.0
                score = float(b97_energy) + float(correction)
            else:
                any_failure = True
                b97_energy = metadata.get("b973c_sp_energy_hartree")
                if b97_energy is None:
                    b97_energy = candidate.score
                score = float(b97_energy)
            corrected.append(DedupCandidate(candidate.path, score, candidate.signature, metadata))

        if any_failure:
            if failure_policy == "strict":
                raise RuntimeError(
                    "xTB mRRHO failed for one or more thermo-pool candidates "
                    "and failure_policy is 'strict'"
                )
            logger.warning(
                "xTB mRRHO had partial failures; falling back to B97-3c "
                "electronic-energy ranking for entire thermo pool"
            )
            for c in corrected:
                b97 = c.metadata.get("b973c_sp_energy_hartree")
                if b97 is not None:
                    c_with_b97_score = DedupCandidate(c.path, float(b97), c.signature, dict(c.metadata))
                    corrected[corrected.index(c)] = c_with_b97_score

        return corrected

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
        kept: List[DedupCandidate],
        output_dir: Path,
        pipeline_log: List[Dict[str, Any]],
        thermo_enabled: bool,
    ) -> Dict[str, Any]:
        if not kept:
            raise RuntimeError("CENSO-LITE produced no valid candidates after filtering")
        kept.sort(key=lambda item: (item.score, str(item.path)))
        minimum = kept[0].score

        candidate_dir = output_dir / "candidates"
        candidate_dir.mkdir(parents=True, exist_ok=True)
        manifest_candidates: List[Dict[str, Any]] = []
        for index, candidate in enumerate(kept, start=1):
            destination = candidate_dir / f"conf_{index:04d}.xyz"
            if candidate.path.resolve() != destination.resolve():
                destination.write_bytes(candidate.path.read_bytes())
            payload = dict(candidate.metadata)
            payload.update(
                {
                    "id": f"{self.molecule_name}_conf_{index:04d}",
                    "xyz": str(destination.relative_to(self.work_dir)),
                    "relative_energy_kcal": (candidate.score - minimum) * HARTREE_TO_KCAL,
                    "s1_score_hartree": candidate.score,
                    "source": "crest_gfn2_censo_lite",
                    "status": "valid",
                }
            )
            manifest_candidates.append(payload)

        ranking_formula = (
            "E_B97-3c_SP + G(RRHO)_xTB"
            if thermo_enabled
            else "E_B97-3c_SP"
        )
        selection_cfg = self._selection_cfg()
        manifest = {
            "schema_version": "s1_censo_lite_v1",
            "stage": "S1",
            "protocol": self.protocol_name,
            "molecule": self.molecule_name,
            "ranking_formula": ranking_formula,
            "final_window_kcal": float(selection_cfg.get("final_window_kcal", 2.5)),
            "pipeline_log": pipeline_log,
            "candidates": manifest_candidates,
            "selected": manifest_candidates[0]["id"],
            "selected_xyz": "selected.xyz",
        }
        manifest_path = output_dir / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        selected_candidate = manifest_candidates[0]
        selected_src = self.work_dir / str(selected_candidate["xyz"])
        selected_dst = output_dir / "selected.xyz"
        _ = shutil.copy2(selected_src, selected_dst)
        return {"manifest": manifest_path, "selected_xyz": selected_dst, "data": manifest}
