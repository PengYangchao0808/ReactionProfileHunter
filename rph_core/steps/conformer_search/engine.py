# pyright: reportImportCycles=false

"""
Unified Conformer Engine (UCE) - v3.0 (Refactored)
=====================================================
A rigorous, "Gold Standard" conformer search engine for ReactionProfileHunter v2.2.

This version implements:
1. Molecule-autonomous directory structure (flat architecture)
2. OPT-SP coupled execution loop (single-conformer atomic operation)
3. Robust Log Parser for extracting last converged coordinates
4. Naming-driven rescue and version control

Directory Structure:
    S1_ConfGeneration/[Molecule_Name]/
        ├── crest/
        ├── xtb/
        ├── cluster/
        ├── prescan/
        ├── fastsp/
        ├── finalDFT/
        └── [Molecule_Name]_global_min.xyz

Logic ported from: Original_Eddition/Conf_Search_20251222/config/confsearch.lib.sh
"""

# type: ignore

import logging
import os
import shutil
import subprocess
import json
import math
import threading
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdMolAlign

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import write_xyz, read_xyz
from rph_core.utils.qc_interface import (
    XTBInterface,
    CRESTInterface,
    GaussianInterface,
    QCResult,
    QCInterfaceFactory
)
from rph_core.utils.geometry_tools import GeometryUtils, LogParser
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.utils.qc_task_runner import QCTaskRunner
from rph_core.utils.qc_runner import run_with_timeout, QCTimeoutError
from rph_core.utils.optimization_config import build_gaussian_route_from_config
from rph_core.utils.isostat_runner import run_isostat
from rph_core.utils.ui import get_progress_manager
from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.steps.conformer_search.state_manager import ConformerStateManager
from rph_core.steps.conformer_search.candidates import CandidateSet, ConformerCandidate, candidate_set_from_paths
from rph_core.steps.conformer_search.pipeline.executor import PipelineExecutor
from rph_core.steps.conformer_search.protocols import ProtocolSpec, resolve_protocol_spec

from rph_core.utils.intra_reaction_scheduler import (
    IntraReactionScheduler,
    ParallelQCJob,
    ResourceLane,
    ParallelResult,
)


logger = logging.getLogger(__name__)
_STATE_MANAGER_LOCK = threading.Lock()

class ConformerEngine(LoggerMixin):
    """
    Unified Conformer Engine (UCE) - v3.1 (Two-Stage Architecture)

    This engine performs conformer search using a two-stage CREST workflow:
    - Stage 1: GFN0-xTB (fast, broad sampling) → ISOSTAT clustering
    - Stage 2: GFN2-xTB (refinement) → ISOSTAT clustering
    - Stage 3: DFT OPT-SP coupled optimization (single-conformer atomic operation)

    Supports backward-compatible single-stage mode (GFN2 only) via config toggle.

    Directory Structure (Two-Stage Mode):
        S1_ConfGeneration/[Molecule_Name]/
            ├── crest/
            │   ├── crest_conformers.xyz
            │   └── ensemble.xyz          # Final ensemble for downstream processing
            ├── xtb/
            │   ├── stage1_gfn0/          # GFN0 coarse search
            │   │   ├── crest_conformers.xyz
            │   │   └── cluster/
            │   │       ├── cluster.xyz   # Representatives for Stage 2
            │   │       └── isostat.log
            │   ├── stage2_gfn2/          # GFN2 refinement
            │   │   ├── crest_ensemble.xyz
            │   │   └── cluster/
            │   │       ├── cluster.xyz   # Final ensemble candidates
            │   │       └── isostat.log
            ├── cluster/                  # Single-stage clustering output
            ├── prescan/                  # FULL protocol fast SP pre-screening
            ├── fastsp/                   # FULL/LITE protocol fast SP screening
            └── finalDFT/                 # Final DFT OPT-SP jobs

    Note on enabled semantics:
    - ConformerEngine does not have an 'enabled' flag.
    - When instantiated, it always operates in real engine mode.
    - To disable conformer search: do not instantiate ConformerEngine.
    - Placeholder mode is not implemented; use conditional instantiation instead.
    """
    def __init__(
        self,
        config: Dict[str, Any],
        work_dir: Path,
        molecule_name: str,
        protocol_spec: Optional[ProtocolSpec] = None,
    ):
        super().__init__()
        self.config = config
        self.molecule_name = molecule_name

        # [v3.0 FIX] Resolve all paths to absolute paths immediately
        work_dir = Path(work_dir).resolve()

        # Create molecule-autonomous directory structure
        self.molecule_dir = (work_dir / molecule_name).resolve()
        self.molecule_dir.mkdir(parents=True, exist_ok=True)

        self.crest_dir = (self.molecule_dir / "crest").resolve()
        self.xtb_dir = (self.molecule_dir / "xtb").resolve()
        self.cluster_dir = (self.molecule_dir / "cluster").resolve()
        self.final_dft_dir = (self.molecule_dir / "finalDFT").resolve()
        self.prescan_dir = (self.molecule_dir / "prescan").resolve()
        self.fastsp_dir = (self.molecule_dir / "fastsp").resolve()
        self.crest_dir.mkdir(parents=True, exist_ok=True)
        self.cluster_dir.mkdir(parents=True, exist_ok=True)
        self.xtb_dir.mkdir(parents=True, exist_ok=True)
        self.final_dft_dir.mkdir(parents=True, exist_ok=True)
        self.prescan_dir.mkdir(parents=True, exist_ok=True)
        self.fastsp_dir.mkdir(parents=True, exist_ok=True)

        self.dft_dir = self.final_dft_dir
        self.state_manager = ConformerStateManager(self.molecule_dir, self.molecule_name)

        self.logger.info(f"📁 S1 workspace: {self.molecule_dir.name}/")

        self.step1_config = config.get('step1', {})
        self.protocol = str(self.step1_config.get('protocol', 'ext')).lower()
        self.protocol_spec = protocol_spec or resolve_protocol_spec(config, self.protocol)
        self.pipeline_executor = PipelineExecutor(self.protocol_spec)
        self.last_final_opt_sp_meta: Dict[str, Any] = {}
        self.crest_config = self.step1_config.get('crest', {})
        self.confsearch_config = self.step1_config.get('conformer_search', {})
        self.theory_opt = config.get('theory', {}).get('optimization', {})
        self.theory_sp = config.get('theory', {}).get('single_point', {})
        self.thermo_config = config.get('thermo', {})
        self.solvent_config = config.get('solvent', {})

        # ==============================================================
        # Two-Stage Configuration (v3.1)
        # ==============================================================
        self.two_stage_enabled = self.confsearch_config.get('two_stage_enabled', False)
        self.logger.info(f"🔀 Two-stage conformer search: {'enabled' if self.two_stage_enabled else 'disabled'}")

        # Stage 1: GFN0 configuration
        stage1_config = self.confsearch_config.get('stage1_gfn0', {})
        self.stage1_enabled = stage1_config.get('enabled', True)
        self.stage1_gfn_level = stage1_config.get('gfn_level', 0)
        self.stage1_energy_window = stage1_config.get('energy_window_kcal', 10.0)
        self.stage1_crest_flags = stage1_config.get('crest_flags', '')
        stage1_clustering = stage1_config.get('clustering', {})
        self.stage1_run_clustering = stage1_clustering.get('run_after', True)
        self.stage1_isostat_gdis = stage1_clustering.get('isostat_gdis', 0.125)
        self.stage1_isostat_edis = stage1_clustering.get('isostat_edis', 1.0)

        # Stage 2: GFN2 configuration
        stage2_config = self.confsearch_config.get('stage2_gfn2', {})
        self.stage2_enabled = stage2_config.get('enabled', True)
        self.stage2_gfn_level = stage2_config.get('gfn_level', 2)
        self.stage2_energy_window = stage2_config.get('energy_window_kcal', 3.0)
        self.stage2_crest_flags = stage2_config.get('crest_flags', '')
        stage2_clustering = stage2_config.get('clustering', {})
        self.stage2_run_clustering = stage2_clustering.get('run_after', True)
        self.stage2_isostat_gdis = stage2_clustering.get('isostat_gdis', 0.125)
        self.stage2_isostat_edis = stage2_clustering.get('isostat_edis', 1.0)

        # Common configuration
        common_config = self.confsearch_config.get('common', {})
        self.energy_window = self.confsearch_config.get('energy_window_kcal', common_config.get('energy_window_kcal', 3.0))
        self.solvent = self.solvent_config.get('name', common_config.get('solvent', self.crest_config.get('solvent', 'acetone')))
        self.nproc = common_config.get('threads', self.crest_config.get('threads', 8))
        self.max_conformers = self.confsearch_config.get('ngeom_max', common_config.get('ngeom_max', 20))
        self.ngeom_default = self.confsearch_config.get('ngeom_default', common_config.get('ngeom_default', 6))

        # ISOSTAT configuration (uses Stage 1 params as defaults for single-stage)
        self.isostat_gdis = self.confsearch_config.get('isostat_gdis', self.stage1_isostat_gdis)
        self.isostat_edis = self.confsearch_config.get('isostat_edis', self.stage1_isostat_edis)
        self.isostat_temp_k = self.thermo_config.get('temperature_k', 298.15)
        self.isostat_threads = self.confsearch_config.get('isostat_threads', self.nproc)
        executables = config.get('executables', {})
        self.isostat_bin = Path(executables.get('isostat', {}).get('path', 'isostat'))
        self.shermo_bin = Path(executables.get('shermo', {}).get('path', 'Shermo'))

        # Initialize CREST Interface (uses Stage 2 GFN level by default for single-stage)
        default_gfn = self.stage2_gfn_level if self.two_stage_enabled else self.crest_config.get('gfn_level', 2)
        self.crest_interface = CRESTInterface(
            gfn_level=default_gfn,
            solvent=self.solvent,
            nproc=self.nproc,
            config=config
        )

        # QC Interfaces for OPT and SP
        opt_engine_type = self.theory_opt.get('engine', 'gaussian').lower()
        sp_engine_type = self.theory_sp.get('engine', 'orca').lower()

        self.qc_opt_interface = QCInterfaceFactory.create_interface(
            opt_engine_type,
            nprocshared=self.theory_opt.get('nproc', 16),
            mem=self.theory_opt.get('mem', '16GB'),
            charge=0,
            multiplicity=1,
            config=self.config
        )

        # ORCA for high-precision SP
        self.orca_sp = ORCAInterface(
            method=self.theory_sp.get('method', 'M062X'),
            basis=self.theory_sp.get('basis', 'def2-TZVPP'),
            nprocs=self.theory_sp.get('nproc', 16),
            solvent=self.theory_sp.get('solvent', 'acetone'),
            config=config
        )

        # Intra-reaction parallel scheduler (P2)
        self._parallel_scheduler = IntraReactionScheduler(config)

    def _display_path(self, path: Path) -> str:
        path = Path(path)
        for root in (self.molecule_dir, Path.cwd()):
            try:
                return str(path.resolve().relative_to(root.resolve()))
            except Exception:
                continue
        parts = path.parts
        if len(parts) >= 3:
            return str(Path("...") / parts[-2] / parts[-1])
        return path.name or str(path)

    def _protocol_tier_label(self) -> str:
        mapping = {
            'ext': 'RPH baseline',
            'default': 'RPH baseline',
            'full': 'CENSO-like full funnel',
            'lite': 'CENSO-like lite funnel',
            'zero': 'CENSO-like zero funnel',
        }
        return mapping.get(self.protocol_spec.name, self.protocol_spec.name)

    def _protocol_funnel_summary(self) -> str:
        mapping = {
            'ext': 'GFN0→GFN2 two-stage ensemble → all candidates',
            'default': 'GFN0→GFN2 two-stage ensemble → all candidates',
            'full': 'GFN2 ensemble → prescreen fast-SP → screening fast-SP → survivor window',
            'lite': 'GFN2 ensemble → screening-level fast-SP → cutoff',
            'zero': 'GFN2 ensemble → narrow energy window',
        }
        return mapping.get(self.protocol_spec.name, 'protocol-defined funnel')

    def _protocol_handoff_summary(self) -> str:
        return 'shared DFT handoff → Gaussian OPT/FREQ → ORCA final SP'

    def _protocol_progress_label(self) -> str:
        mapping = {
            'ext': 'EXT funnel · two-stage baseline',
            'default': 'EXT funnel · two-stage baseline',
            'full': 'FULL funnel · prescreen + screening',
            'lite': 'LITE funnel · screening + cutoff',
            'zero': 'ZERO funnel · narrow window',
        }
        return mapping.get(self.protocol_spec.name, 'Protocol funnel')

    def run(self, smiles: str) -> Tuple[Path, float]:
        """
        Execute full conformer search workflow (OPT-SP coupled).

        Automatically dispatches to single-stage or two-stage workflow based on config.

        Args:
            smiles: SMILES string of the molecule

        Returns:
            Tuple[Path, float]: (Path to global min XYZ, Global min SP Energy in Hartree)
        """
        self.logger.info(f"[S1] 🚀 Starting V4.0 CENSO-like conformer search: {self.molecule_name}")
        pm = get_progress_manager()
        if pm:
            pm.update_step("s1", description=f"S1: [{self.molecule_name}] Starting conformer search")
            pm.log_event("S1", f"Starting conformer search for {self.molecule_name}")

        mode = "two-stage (GFN0→GFN2)" if self.two_stage_enabled else "single-stage (GFN2)"
        self.logger.info(f"[S1]    Tier         : {self.protocol_spec.name.upper()} · {self._protocol_tier_label()}")
        self.logger.info(f"[S1]    Funnel       : {self._protocol_funnel_summary()}")
        self.logger.info(f"[S1]    Handoff      : {self._protocol_handoff_summary()}")
        self.logger.info(f"[S1]    Search mode  : {mode}")
        self.logger.info(f"[S1]    Candidate cap: default={self.ngeom_default}, max={self.max_conformers}")
        self.state_manager.start_run(smiles=smiles, two_stage_enabled=self.two_stage_enabled)
        self.state_manager.set_protocol_signature(
            protocol=self.protocol_spec.name,
            funnel_signature={
                "search_mode": self.protocol_spec.funnel_policy.search_mode,
                "clustering_mode": self.protocol_spec.funnel_policy.clustering_mode,
                "prescreen_mode": self.protocol_spec.funnel_policy.prescreen_mode,
                "rerank_mode": self.protocol_spec.funnel_policy.rerank_mode,
            },
            handoff_signature={
                "mode": self.protocol_spec.handoff_policy.mode,
                "fallback_mode": self.protocol_spec.handoff_policy.fallback_mode,
            },
        )
        self._emit_s1_progress("run_started", {"mode": mode, "smiles": smiles})
        if pm:
            pm.log_event("S1", f"{self._protocol_tier_label()} · {mode}")

        try:
            global_min_path = self.molecule_dir / f"{self.molecule_name}_global_min.xyz"

            if self._is_valid_xyz(global_min_path):
                self.logger.info(f"[S1] ⏭️ Found {global_min_path.name}, skipping conformer search")
                if pm:
                    pm.log_event("S1", f"Found cache {global_min_path.name}, skipping conformer search")

                energy = self._extract_energy_from_xyz(global_min_path)
                self.state_manager.set_global_min(self.molecule_name, energy, global_min_path)
                self.state_manager.mark_run_complete()
                self._emit_s1_progress(
                    "run_reused",
                    {"global_min_xyz": str(global_min_path), "energy_hartree": energy},
                )
                self.logger.info(f"[S1]    Extracted energy: {energy:.6f} Hartree")
                if pm:
                    pm.log_event("S1", f"Extracted cached energy: {energy:.6f} Hartree")

                return global_min_path, energy

            if pm:
                pm.enter_phase("S1", "Stage 1: RDKit conformer generation")
                pm.update_step("s1", completed=5, description=f"S1: [{self.molecule_name}] RDKit Embedding")
            initial_xyz = self._step_rdkit_embed(smiles)

            if pm:
                pm.enter_phase("S1", "Stage 2: Funnel runtime")
                pm.update_step("s1", completed=15, description=f"S1: [{self.molecule_name}] {self._protocol_progress_label()}")
            candidate_set = self._run_protocol_funnel(initial_xyz)
            candidates = [candidate.source_xyz for candidate in candidate_set.candidates]
            self.logger.info(
                f"[S1] 🔍 Funnel summary: basis={candidate_set.ranking_basis}, "
                f"candidates={candidate_set.total_count}, survivors={len(candidates)}"
            )
            if pm:
                pm.log_event("S1", f"{self._protocol_progress_label()} → {len(candidates)} DFT candidates")

            self._write_crest_ensemble_energies_json()
            self._write_fastsp_ensemble_energies_json()

            if not candidates:
                raise RuntimeError("No valid conformers found after CREST processing.")

            if pm:
                pm.enter_phase("S1", "Stage 4: DFT optimization")
                pm.update_step("s1", completed=55, description=f"S1: [{self.molecule_name}] DFT optimization")
            best_log, min_energy, stage_meta = self._run_shared_dft_handoff(candidate_set)
            self.last_final_opt_sp_meta = dict(stage_meta)
            self.state_manager.set_funnel_summary(
                ranking_basis=str(candidate_set.ranking_basis),
                candidate_total=int(candidate_set.total_count),
                selected_candidate_count=int(stage_meta.get("selected_candidate_count", len(candidates))),
                fallback_triggered=bool(stage_meta.get("fallback_triggered", False)),
                window_count=int(candidate_set.metadata.get("window_count", len(candidate_set.candidates))),
                gap_rank1_rank2=candidate_set.metadata.get("gap_rank1_rank2"),
                mode_effective=str(stage_meta.get("mode_effective", stage_meta.get("handoff_mode", ""))),
                stages_executed=list(candidate_set.metadata.get("stages_executed", [])),
            )

            global_min_path = self.molecule_dir / f"{self.molecule_name}_global_min.xyz"

            coords, symbols, _ = LogParser.extract_last_converged_coords(
                best_log,
                engine_type='auto'
            )

            if coords is None:
                raise RuntimeError(f"Failed to extract coordinates from {best_log}")
            if symbols is None:
                _, symbols = read_xyz(best_log)

            if symbols is None:
                symbols = ["X"] * len(coords)

            write_xyz(global_min_path, coords, symbols, title=f"Global Min SP E={min_energy:.6f}")
            best_name = Path(best_log).stem.replace("_Res", "")
            self.state_manager.set_global_min(best_name, min_energy, global_min_path)
            self.state_manager.mark_run_complete()
            self._emit_s1_progress(
                "run_completed",
                {
                    "global_min_xyz": str(global_min_path),
                    "energy_hartree": min_energy,
                    "best_conformer": best_name,
                },
            )
            self.logger.info(f"[S1] 🏆 Global minimum: {self._display_path(global_min_path)}")
            self.logger.info(f"[S1]    SP Energy: {min_energy:.6f} Hartree")
            if pm:
                pm.update_step("s1", completed=98, description=f"S1: [{self.molecule_name}] Finalizing global minimum")
                pm.log_event("S1", f"Global minimum found: {global_min_path.name} (E={min_energy:.6f} Hartree)")

            return global_min_path, min_energy
        except Exception as exc:
            self.state_manager.mark_run_failed(str(exc))
            self._emit_s1_progress("run_failed", {"error": str(exc)})
            raise

    def _step_rdkit_embed(self, smiles: str, num_conf: int = 1) -> Path:
        """
        Generate initial 3D structure using RDKit.
        
        Args:
            smiles: Input SMILES string
            num_conf: Number of conformers to generate (default: 1)
            
        Returns:
            Path to generated XYZ file
        """
        self.logger.info(f"[S1]   [1/5] RDKit 3D embedding (num_conf={num_conf})...")
        pm = get_progress_manager()
        if pm:
            pm.log_event("S1", f"RDKit 3D embedding started (num_conf={num_conf})")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES for RDKit embedding: {smiles}")
        mol = Chem.AddHs(mol)

        params = rdDistGeom.ETKDG()
        embed_status = rdDistGeom.EmbedMolecule(mol, params)

        if embed_status != 0:
            self.logger.warning(
                f"[S1] ETKDG embedding failed for {self.molecule_name}; "
                "retrying with random coordinates"
            )
            params = rdDistGeom.ETKDGv3()
            params.useRandomCoords = True
            embed_status = rdDistGeom.EmbedMolecule(mol, params)

        if embed_status != 0:
            raise RuntimeError(
                f"RDKit embedding failed for {self.molecule_name} "
                f"(SMILES: {smiles}); stereochemistry may be geometrically unrealizable"
            )

        try:
            rdForceFieldHelpers.MMFFOptimizeMolecule(mol)
        except Exception as exc:
            self.logger.warning(
                f"[S1] MMFF optimization failed for {self.molecule_name}: {exc}; "
                "falling back to UFF"
            )
            try:
                rdForceFieldHelpers.UFFOptimizeMolecule(mol)
            except Exception as uff_exc:
                self.logger.warning(
                    f"[S1] UFF fallback also failed for {self.molecule_name}: {uff_exc}; "
                    "using unoptimized coordinates"
                )

        output_path = self.crest_dir / f"{self.molecule_name}_init.xyz"

        conf = mol.GetConformer()
        coords = conf.GetPositions()
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

        assert symbols is not None, "Symbols cannot be None after RDKit generation"

        write_xyz(output_path, coords, symbols, title="RDKit_Initial")
        if pm:
            pm.log_event("S1", f"RDKit embedding completed: {output_path.name}")
        return output_path

    def _step_crest_search(self, input_xyz: Path) -> Path:
        """Single-stage CREST search (GFN2 only) - backward compatible."""
        self.logger.info("[S1]   [2/5] GFN2 ensemble generation (single-stage)...")
        pm = get_progress_manager()
        if pm:
            pm.log_event("S1", "CREST global search started (single-stage GFN2)")

        crest_ensemble_file = self.crest_dir / "ensemble.xyz"

        if self._is_valid_xyz(crest_ensemble_file, min_size=500):
            self.logger.info("[S1] ⏭️ Found CREST cache, skipping search")
            self.state_manager.mark_crest_stage(
                stage_name="single_stage",
                status="cached",
                output_file=crest_ensemble_file,
            )
            if pm:
                pm.log_event("S1", "Found CREST cache, skipping single-stage search")
            self._emit_s1_progress("crest_single_stage_cached", {"ensemble": str(crest_ensemble_file)})
            return crest_ensemble_file

        best_xyz = self.crest_interface.run_conformer_search(input_xyz, self.crest_dir)

        ensemble_file = self.crest_dir / "crest_conformers.xyz"
        if ensemble_file.exists():
            shutil.copy(ensemble_file, crest_ensemble_file)
            self.state_manager.mark_crest_stage(
                stage_name="single_stage",
                status="completed",
                output_file=crest_ensemble_file,
            )
            if pm:
                pm.log_event("S1", f"CREST single-stage completed: {crest_ensemble_file.name}")
            self._emit_s1_progress("crest_single_stage_completed", {"ensemble": str(crest_ensemble_file)})
            return crest_ensemble_file

        self.logger.warning("crest_conformers.xyz not found, using best struct only.")
        shutil.copy(best_xyz, crest_ensemble_file)
        self.state_manager.mark_crest_stage(
            stage_name="single_stage",
            status="fallback",
            output_file=crest_ensemble_file,
        )
        if pm:
            pm.log_event("S1", f"CREST fallback to best structure: {crest_ensemble_file.name}")
        self._emit_s1_progress("crest_single_stage_fallback", {"ensemble": str(crest_ensemble_file)})
        return crest_ensemble_file

    def _step_two_stage_crest(self, input_xyz: Path) -> Path:
        """Execute two-stage CREST workflow: GFN0 → ISOSTAT → GFN2 → ISOSTAT.

        Stage 1: GFN0-xTB for fast, broad conformational space sampling
        Stage 2: GFN2-xTB for refinement of Stage 1 candidates

        Returns:
            Path to final ensemble file (after Stage 2 + clustering)
        """
        # ==============================================================
        # Stage 1: GFN0 Conformer Search
        # ==============================================================
        self.logger.info("[S1]   [2a/5] Stage 1: CREST GFN0 Conformer Search (Fast Sampling)...")
        pm = get_progress_manager()
        if pm:
            pm.enter_phase("S1", "Stage 2: CREST searching")
            pm.log_event("S1", "Stage 1 CREST GFN0 search started")

        stage1_dir = self.xtb_dir / "stage1_gfn0"
        stage1_dir.mkdir(parents=True, exist_ok=True)
        stage1_ensemble = stage1_dir / "crest_conformers.xyz"
        if self._is_valid_xyz(stage1_ensemble, min_size=500):
            gfn0_ensemble = stage1_ensemble
            self.logger.info("[S1]   [2a/5] Stage 1 ensemble cached, skipping GFN0 search")
            self.state_manager.mark_crest_stage("stage1_gfn0", "cached", output_file=stage1_ensemble)
        else:
            gfn0_ensemble = self._run_crest_stage(
                input_xyz=input_xyz,
                output_dir=stage1_dir,
                gfn_level=self.stage1_gfn_level,
                stage_name="GFN0",
                additional_flags=self.stage1_crest_flags
            )
            self.state_manager.mark_crest_stage("stage1_gfn0", "completed", output_file=gfn0_ensemble)

        # ==============================================================
        # Stage 1: ISOSTAT Clustering
        # ==============================================================
        if self.stage1_run_clustering:
            self.logger.info("[S1]   [2b/5] Stage 1: ISOSTAT Clustering...")
            stage1_cluster = stage1_dir / "cluster" / "cluster.xyz"
            if self._is_valid_xyz(stage1_cluster):
                gfn0_clustered = stage1_cluster
                self.logger.info("[S1]   [2b/5] Stage 1 cluster cached, skipping")
                self.state_manager.mark_crest_stage("stage1_cluster", "cached", output_file=stage1_cluster)
            else:
                gfn0_clustered = self._run_isostat_clustering(
                    input_ensemble=gfn0_ensemble,
                    output_dir=stage1_dir,
                    energy_window=self.stage1_energy_window,
                    gdis=self.stage1_isostat_gdis,
                    edis=self.stage1_isostat_edis
                )
                self.state_manager.mark_crest_stage("stage1_cluster", "completed", output_file=gfn0_clustered)
        else:
            gfn0_clustered = gfn0_ensemble

        # ==============================================================
        # Stage 2: GFN2 Batch Optimization
        # ==============================================================
        self.logger.info("[S1]   [2c/5] Stage 2: CREST GFN2 Batch Optimization (Refinement)...")
        if pm:
            pm.enter_phase("S1", "Stage 3: xTB optimization")
            pm.log_event("S1", "Stage 2 GFN2 xTB batch optimization started")

        stage2_dir = self.xtb_dir / "stage2_gfn2"
        stage2_dir.mkdir(parents=True, exist_ok=True)
        stage2_ensemble = stage2_dir / "crest_ensemble.xyz"
        if self._is_valid_xyz(stage2_ensemble, min_size=500):
            gfn2_ensemble = stage2_ensemble
            self.logger.info("[S1]   [2c/5] Stage 2 ensemble cached, skipping GFN2 optimization")
            self.state_manager.mark_crest_stage("stage2_gfn2", "cached", output_file=stage2_ensemble)
        else:
            gfn2_ensemble = self._run_crest_batch_optimization(
                input_xyz=gfn0_clustered,
                output_dir=stage2_dir,
                gfn_level=self.stage2_gfn_level,
                stage_name="GFN2",
                additional_flags=self.stage2_crest_flags
            )
            self.state_manager.mark_crest_stage("stage2_gfn2", "completed", output_file=gfn2_ensemble)

        # ==============================================================
        # Stage 2: ISOSTAT Clustering
        # ==============================================================
        if self.stage2_run_clustering:
            self.logger.info("[S1]   [2d/5] Stage 2: ISOSTAT Clustering...")
            stage2_cluster = stage2_dir / "cluster" / "cluster.xyz"
            if self._is_valid_xyz(stage2_cluster):
                gfn2_clustered = stage2_cluster
                self.logger.info("[S1]   [2d/5] Stage 2 cluster cached, skipping")
                self.state_manager.mark_crest_stage("stage2_cluster", "cached", output_file=stage2_cluster)
            else:
                gfn2_clustered = self._run_isostat_clustering(
                    input_ensemble=gfn2_ensemble,
                    output_dir=stage2_dir,
                    energy_window=self.stage2_energy_window,
                    gdis=self.stage2_isostat_gdis,
                    edis=self.stage2_isostat_edis
                )
                self.state_manager.mark_crest_stage("stage2_cluster", "completed", output_file=gfn2_clustered)
        else:
            gfn2_clustered = gfn2_ensemble

        # ==============================================================
        # Final: Copy to standard location for downstream processing
        # ==============================================================
        self.logger.info("[S1]   [2e/5] Finalizing ensemble for DFT processing...")

        final_ensemble = self.crest_dir / "ensemble.xyz"
        if not final_ensemble.exists() or final_ensemble.resolve() != gfn2_clustered.resolve():
            shutil.copy(gfn2_clustered, final_ensemble)
        self.state_manager.mark_crest_stage("final_ensemble", "completed", output_file=final_ensemble)
        self._emit_s1_progress("crest_two_stage_completed", {"ensemble": str(final_ensemble)})

        self.logger.info(f"[S1]   ✓ Two-stage funnel complete: {final_ensemble.name}")
        if pm:
            pm.log_event("S1", f"Two-stage CREST completed: {final_ensemble.name}")

        return final_ensemble

    def _run_crest_stage(
        self,
        input_xyz: Path,
        output_dir: Path,
        gfn_level: int,
        stage_name: str,
        additional_flags: Optional[str] = None
    ) -> Path:
        """Run CREST conformer search for a single stage.

        Args:
            input_xyz: Input structure file
            output_dir: Stage output directory
            gfn_level: GFN level (0, 1, or 2)
            stage_name: Display name for logging
            additional_flags: Additional CREST command-line flags

        Returns:
            Path to ensemble file
        """
        self.logger.info(f"    [{stage_name}] Running conformer search (GFN{gfn_level})...")

        try:
            result = self.crest_interface.run_conformer_search(
                xyz_file=input_xyz,
                output_dir=output_dir,
                gfn_override=gfn_level,
                additional_flags=additional_flags
            )

            ensemble_file = output_dir / "crest_conformers.xyz"
            if result != ensemble_file and ensemble_file.exists():
                shutil.copy(ensemble_file, result)

            self.logger.info(f"    [{stage_name}] ✓ Found {result.name}")
            return result

        except Exception as e:
            self.logger.warning(f"    [{stage_name}] CREST failed: {e}, using input structure")
            fallback = output_dir / input_xyz.name
            shutil.copy(input_xyz, fallback)
            return fallback

    def _run_crest_batch_optimization(
        self,
        input_xyz: Path,
        output_dir: Path,
        gfn_level: int,
        stage_name: str,
        additional_flags: Optional[str] = None
    ) -> Path:
        """Run CREST batch optimization on ensemble structures.

        Args:
            input_xyz: Input ensemble file (from Stage 1)
            output_dir: Stage output directory
            gfn_level: GFN level for optimization
            stage_name: Display name for logging
            additional_flags: Additional CREST command-line flags

        Returns:
            Path to optimized ensemble file
        """
        self.logger.info(f"    [{stage_name}] Running batch optimization on ensemble...")

        try:
            result = self.crest_interface.run_batch_optimization(
                ensemble_xyz=input_xyz,
                output_dir=output_dir,
                gfn_level=gfn_level,
                additional_flags=additional_flags
            )

            # Prefer crest_ensemble.xyz if it exists
            optimized = output_dir / "crest_ensemble.xyz"
            if optimized.exists():
                result = optimized

            self.logger.info(f"    [{stage_name}] ✓ Batch optimization complete: {result.name}")
            return result

        except Exception as e:
            self.logger.warning(f"    [{stage_name}] Batch optimization failed: {e}")
            return input_xyz

    def _run_isostat_clustering(
        self,
        input_ensemble: Path,
        output_dir: Path,
        energy_window: float,
        gdis: float,
        edis: float
    ) -> Path:
        """Run ISOSTAT clustering on ensemble structures.

        Args:
            input_ensemble: Input ensemble XYZ file
            output_dir: Directory for clustering output
            energy_window: Energy window for filtering (kcal/mol)
            gdis: RMSD distance threshold for clustering (Å)
            edis: Energy distance threshold for clustering (kcal/mol)

        Returns:
            Path to cluster.xyz (representative structures)
        """
        cluster_dir = output_dir / "cluster"
        cluster_dir.mkdir(parents=True, exist_ok=True)

        # Copy ensemble to cluster directory
        isomers_xyz = cluster_dir / "isomers.xyz"
        shutil.copy(input_ensemble, isomers_xyz)

        # Run ISOSTAT
        result = run_isostat(
            isostat_bin=self.isostat_bin,
            input_xyz=isomers_xyz,
            output_dir=cluster_dir,
            gdis=gdis,
            edis=edis,
            temp_k=self.isostat_temp_k,
            threads=self.isostat_threads,
            energy_window=energy_window
        )

        self.logger.info(f"    ✓ Clustering complete: {result.cluster_xyz.name}")
        return result.cluster_xyz

    def _step_process_ensemble(self, ensemble_file: Path) -> List[Path]:
        self.logger.info("[S1]   [3/5] Ensemble processing (clustering + filtering)...")

        isomers_xyz = self.cluster_dir / "isomers.xyz"
        shutil.copy(ensemble_file, isomers_xyz)

        result = run_isostat(
            isostat_bin=self.isostat_bin,
            input_xyz=isomers_xyz,
            output_dir=self.cluster_dir,
            gdis=self.isostat_gdis,
            edis=self.isostat_edis,
            temp_k=self.isostat_temp_k,
            threads=self.isostat_threads,
            energy_window=self.energy_window
        )

        if result.n_within_window is None:
            n_calc = min(self.ngeom_default, self.max_conformers)
        else:
            n_calc = min(result.n_within_window, self.max_conformers)

        conformers = self._split_xyz_ensemble(result.cluster_xyz, self.final_dft_dir, limit=n_calc)
        self.logger.info(f"    - Found {len(conformers)} conformers from isostat (Limit: {n_calc}).")
        self.state_manager.mark_crest_stage(
            stage_name="dft_candidates",
            status="completed",
            output_file=result.cluster_xyz,
            metadata={"selected": len(conformers), "limit": n_calc},
        )

        for idx, conf_path in enumerate(conformers):
            self.state_manager.upsert_conformer(
                conf_name=f"conf_{idx:03d}",
                source_xyz=conf_path,
                source_index=idx,
            )

        if not conformers:
            return []

        self.logger.info(f"    - Selected {len(conformers)} conformers for DFT.")
        return conformers

    def _split_xyz_ensemble(self, ensemble_file: Path, output_dir: Optional[Path] = None, limit: Optional[int] = None) -> List[Path]:

        split_dir = output_dir if output_dir is not None else self.crest_dir / "conf_split"
        split_dir.mkdir(parents=True, exist_ok=True)

        paths = []
        with open(ensemble_file, 'r') as f:
            lines = f.readlines()

        i = 0
        count = 0
        while i < len(lines):
            if limit is not None and count >= limit:
                break

            try:
                natoms = int(lines[i].strip())
                block = lines[i:i + natoms + 2]

                out_path = split_dir / f"conf_{count:03d}.xyz"
                with open(out_path, 'w') as out:
                    out.writelines(block)

                paths.append(out_path)

                i += natoms + 2
                count += 1
            except (ValueError, IndexError):
                break

        return paths

    def _rmsd_clustering(self, conf_list: List[Dict[str, Any]], threshold: float) -> List[Dict[str, Any]]:
        sorted_confs = sorted(conf_list, key=lambda x: x['energy'])

        selected = []
        mols_cache = {}

        def get_mol(path):
            if path in mols_cache:
                return mols_cache[path]
            try:
                mol = Chem.MolFromXYZFile(str(path))
                if mol is None:
                    return None
                mols_cache[path] = mol
                return mol
            except Exception:
                return None

        for candidate in sorted_confs:
            cand_mol = get_mol(candidate['path'])
            if cand_mol is None:
                continue

            is_unique = True
            for cluster_rep in selected:
                rep_mol = get_mol(cluster_rep['path'])
                if rep_mol is None:
                    continue

                try:
                    rmsd = rdMolAlign.GetBestRMS(cand_mol, rep_mol)
                    if rmsd < threshold:
                        is_unique = False
                        break
                except Exception:
                    pass

            if is_unique:
                selected.append(candidate)

        self.logger.info(f"    - Clustering reduced {len(conf_list)} -> {len(selected)} structures.")
        return selected

    def _step_dft_opt_sp_coupled(self, candidates: List[Path]) -> Tuple[Path, float, List[Dict[str, Any]]]:
        if self._parallel_scheduler.s1_conformer_parallel:
            return self._step_dft_opt_sp_coupled_parallel(candidates)

        self.logger.info(f"[S1]   [4/5] DFT OPT-SP Coupled Loop (Total {len(candidates)} conformers)...")

        pm = get_progress_manager()
        total_confs = len(candidates)
        if pm:
            pm.enter_phase("S1", "Stage 4: DFT optimization")
            pm.log_event("S1", f"DFT OPT-SP loop started for {total_confs} conformers")
        records: List[Dict[str, Any]] = []

        for idx, xyz_file in enumerate(candidates):
            conf_name = f"conf_{idx:03d}"
            self.state_manager.upsert_conformer(conf_name, Path(xyz_file).resolve(), idx)

            cached_record = self.state_manager.get_conformer_record(conf_name)
            if self.state_manager.is_conformer_completed(conf_name) and isinstance(cached_record, dict):
                cached_log_file = Path(str(cached_record.get("log_file", "")))
                if cached_log_file.exists():
                    self.logger.info(f"[S1]     ⏭️ [{idx+1}/{len(candidates)}] Reusing cached conformer: {conf_name}")
                    if pm:
                        pm.log_event("S1", f"Conformer cached {idx+1}/{total_confs}: {conf_name}")
                    self._emit_s1_progress("conformer_skipped", {"conformer": conf_name, "index": idx + 1, "total": total_confs})
                    records.append(self._restore_record_from_state(cached_record))
                    continue

            if pm:
                progress = 55 + int(((idx + 1) / max(total_confs, 1)) * 40)
                pm.update_step(
                    "s1",
                    completed=min(progress, 95),
                    description=f"S1: [{self.molecule_name}] Opt conformer {idx+1}/{total_confs}"
                )
                pm.set_subtask("S1", "DFT Conformers", idx + 1, total_confs)
                pm.log_event("S1", f"Optimizing conformer {idx+1}/{total_confs}: {conf_name}")

            self.logger.info(f"[S1]     🔄 [{idx+1}/{len(candidates)}] DFT handoff on conformer: {conf_name}")
            self.state_manager.mark_conformer_running(conf_name)
            self._emit_s1_progress("conformer_started", {"conformer": conf_name, "index": idx + 1, "total": total_confs})

            current_xyz_source = Path(xyz_file).resolve()
            converged_this_conf = False
            conf_failed_reason = ""

            for attempt in range(2):
                current_conf_name = f"{conf_name}_Res" if attempt > 0 else conf_name

                gjf_file = self.final_dft_dir / f"{current_conf_name}.gjf"
                log_file = self.final_dft_dir / f"{current_conf_name}.log"

                from rph_core.utils.qc_interface import GaussianInterface

                if attempt > 0:
                    rescue_route = self.theory_opt.get('rescue_route') or build_gaussian_route_from_config(
                        self.config,
                        rescue=True
                    )
                    gauss_int = GaussianInterface(
                        charge=self.theory_opt.get('charge', 0),
                        multiplicity=self.theory_opt.get('multiplicity', 1),
                        nprocshared=self.theory_opt.get('nproc', 16),
                        mem=self.theory_opt.get('mem', '16GB'),
                        config=self.config
                    )
                    gauss_int.write_input_file(
                        current_xyz_source,
                        gjf_file,
                        route=rescue_route,
                        title=f"{self.molecule_name}_{current_conf_name}_Rescue"
                    )
                else:
                    gauss_int = GaussianInterface(
                        charge=self.theory_opt.get('charge', 0),
                        multiplicity=self.theory_opt.get('multiplicity', 1),
                        nprocshared=self.theory_opt.get('nproc', 16),
                        mem=self.theory_opt.get('mem', '16GB'),
                        config=self.config
                    )
                    gauss_int.write_input_file(
                        current_xyz_source,
                        gjf_file,
                        route=self.theory_opt.get('route') or build_gaussian_route_from_config(
                            self.config,
                            rescue=False
                        ),
                        title=f"{self.molecule_name}_{current_conf_name}"
                    )

                opt_converged, next_xyz = self._run_gaussian_opt(
                    gjf_file,
                    log_file,
                    xyz_source=current_xyz_source,
                    attempt=attempt
                )

                if opt_converged:
                    self.state_manager.mark_opt_attempt(conf_name, attempt, "converged", log_file)
                    self._emit_s1_progress("opt_converged", {"conformer": conf_name, "attempt": attempt})
                else:
                    self.state_manager.mark_opt_attempt(
                        conf_name,
                        attempt,
                        "failed",
                        log_file,
                        note="gaussian_opt_not_converged",
                    )
                    self._emit_s1_progress("opt_failed", {"conformer": conf_name, "attempt": attempt})

                if not opt_converged:
                    if attempt == 0:
                        if next_xyz and next_xyz.exists():
                            current_xyz_source = next_xyz
                        conf_failed_reason = "opt_failed_attempt_0"
                        continue
                    conf_failed_reason = "opt_failed_after_rescue"
                    break

                final_coords, final_symbols, parse_error = LogParser.extract_last_converged_coords(
                    log_file,
                    engine_type='gaussian'
                )

                if final_coords is None:
                    conf_failed_reason = "opt_parse_failed"
                    if attempt == 0:
                        continue
                    break

                if final_symbols is None:
                    _, final_symbols = read_xyz(current_xyz_source)

                sp_in_file = self.final_dft_dir / f"{current_conf_name}_SP.inp"
                sp_out_file = self.final_dft_dir / f"{current_conf_name}_SP.out"

                sp_energy = self._run_orca_sp(
                    final_coords,
                    final_symbols,
                    sp_in_file,
                    sp_out_file
                )

                if sp_energy is None:
                    self.state_manager.mark_sp_result(
                        conf_name,
                        status="failed",
                        output_file=sp_out_file,
                        note="orca_sp_failed",
                    )
                    self._emit_s1_progress("sp_failed", {"conformer": conf_name, "output": str(sp_out_file)})
                    self.logger.warning(
                        "      ⚠️  ORCA SP failed after converged OPT; skipping rescue OPT retry and moving on"
                    )
                    conf_failed_reason = "sp_failed"
                    break

                self.state_manager.mark_sp_result(
                    conf_name,
                    status="completed",
                    output_file=sp_out_file,
                    sp_energy=sp_energy,
                )
                self._emit_s1_progress("sp_completed", {"conformer": conf_name, "energy_hartree": sp_energy})

                # Shermo removed — conformer ranking now uses electronic energy (sp_energy).
                # Gibbs/Boltzmann averaging at condition temperature is done in ConditionThermoEngine.
                record = {
                    "name": current_conf_name,
                    "g_used": sp_energy,  # Was Gibbs, now electronic energy for backward compat
                    "sp_energy": sp_energy,
                    "log_file": log_file,
                }
                records.append(record)
                self.state_manager.mark_conformer_completed(conf_name, self._serialize_record_for_state(record))
                self._emit_s1_progress("conformer_completed", {"conformer": conf_name, "energy_hartree": sp_energy})
                converged_this_conf = True
                break

            if not converged_this_conf:
                self.state_manager.mark_conformer_failed(conf_name, conf_failed_reason or "conformer_failed")
                self._emit_s1_progress(
                    "conformer_failed",
                    {"conformer": conf_name, "reason": conf_failed_reason or "conformer_failed"},
                )

        return self._finalize_conformer_results(records, total_confs, pm)

    def _step_dft_opt_sp_coupled_parallel(self, candidates: List[Path]) -> Tuple[Path, float, List[Dict[str, Any]]]:
        """Parallel variant of _step_dft_opt_sp_coupled for dual-lane conformer optimization.

        Processes conformers in pairs using IntraReactionScheduler.
        Falls back to sequential for odd-numbered last conformer.
        """
        self.logger.info(
            f"[S1]   [4/5] DFT OPT-SP Coupled Loop — PARALLEL MODE "
            f"({len(candidates)} conformers, {self._parallel_scheduler.s1_max_conformer_workers} workers)"
        )

        pm = get_progress_manager()
        total_confs = len(candidates)
        if pm:
            pm.enter_phase("S1", "Stage 4: DFT optimization (parallel)")
            pm.log_event("S1", f"DFT OPT-SP parallel loop started for {total_confs} conformers")

        records: List[Dict[str, Any]] = []
        max_workers = max(1, int(self._parallel_scheduler.s1_max_conformer_workers))
        lane: ResourceLane = self._parallel_scheduler.get_s1_conformer_lane()

        idx = 0
        while idx < total_confs:
            batch_end = min(idx + max_workers, total_confs)
            batch = candidates[idx:batch_end]

            to_run: List[Tuple[int, Path]] = []
            for j, xyz_file in enumerate(batch):
                conf_idx = idx + j
                conf_name = f"conf_{conf_idx:03d}"
                xyz_path = Path(xyz_file).resolve()

                with _STATE_MANAGER_LOCK:
                    self.state_manager.upsert_conformer(conf_name, xyz_path, conf_idx)
                    cached_record = self.state_manager.get_conformer_record(conf_name)
                    is_cached = self.state_manager.is_conformer_completed(conf_name)

                if is_cached and isinstance(cached_record, dict):
                    cached_log_file = Path(str(cached_record.get("log_file", "")))
                    if cached_log_file.exists():
                        self.logger.info(
                            f"[S1]     ⏭️ [{conf_idx+1}/{total_confs}] Reusing cached conformer: {conf_name}"
                        )
                        if pm:
                            pm.log_event("S1", f"Conformer cached {conf_idx+1}/{total_confs}: {conf_name}")
                        self._emit_s1_progress(
                            "conformer_skipped",
                            {"conformer": conf_name, "index": conf_idx + 1, "total": total_confs},
                        )
                        restored_record = self._restore_record_from_state(cached_record)
                        restored_record.setdefault("conf_name", conf_name)
                        restored_record.setdefault("g_energy", restored_record.get("g_used"))
                        restored_record["source_index"] = conf_idx
                        records.append(restored_record)
                        continue

                to_run.append((conf_idx, xyz_path))

            if not to_run:
                idx = batch_end
                continue

            if len(to_run) == 1:
                record = self._run_single_conformer_dft(*to_run[0], total_confs, pm)
                if not record.get("failed"):
                    records.append(record)
            else:
                jobs: List[ParallelQCJob] = []
                for conf_idx, xyz_path in to_run:
                    conf_name = f"conf_{conf_idx:03d}"
                    jobs.append(
                        ParallelQCJob(
                            job_id=conf_name,
                            lane=lane,
                            func=self._run_single_conformer_dft,
                            args=(conf_idx, xyz_path, total_confs, pm),
                            description=f"DFT OPT-SP for {conf_name}",
                        )
                    )

                parallel_results: List[ParallelResult] = self._parallel_scheduler.run_parallel(jobs)
                for pr in parallel_results:
                    if pr.success and isinstance(pr.result, dict):
                        if pr.result.get("failed"):
                            self.logger.warning(
                                f"[S1]     ❌ Parallel conformer job {pr.job_id} failed: "
                                f"{pr.result.get('error', 'conformer_failed')}"
                            )
                            continue
                        records.append(pr.result)
                        continue

                    error_message = pr.error or "parallel_conformer_job_failed"
                    self.logger.warning(f"[S1]     ❌ Parallel conformer job {pr.job_id} failed: {error_message}")
                    with _STATE_MANAGER_LOCK:
                        self.state_manager.mark_conformer_failed(pr.job_id, error_message)
                    self._emit_s1_progress(
                        "conformer_failed",
                        {"conformer": pr.job_id, "reason": error_message},
                    )

            idx = batch_end

        return self._finalize_conformer_results(records, total_confs, pm)

    def _run_single_conformer_dft(self, idx: int, xyz_file: Path, total_confs: int, pm: Any) -> Dict[str, Any]:
        """Run DFT OPT+SP for a single conformer. Thread-safe for parallel execution."""
        conf_name = f"conf_{idx:03d}"
        current_xyz_source = Path(xyz_file).resolve()
        s1_parallel_config = self._get_s1_parallel_config()
        opt_nproc = int(s1_parallel_config.get("opt_cores_per_job", self.theory_opt.get("nproc", 16)))
        opt_mem = str(s1_parallel_config.get("opt_mem_per_job", self.theory_opt.get("mem", "16GB")))
        sp_interface = self._build_parallel_orca_sp_interface()

        if pm:
            progress = 55 + int(((idx + 1) / max(total_confs, 1)) * 40)
            pm.update_step(
                "s1",
                completed=min(progress, 95),
                description=f"S1: [{self.molecule_name}] Opt conformer {idx+1}/{total_confs}",
            )
            pm.set_subtask("S1", "DFT Conformers", idx + 1, total_confs)
            pm.log_event("S1", f"Optimizing conformer {idx+1}/{total_confs}: {conf_name}")

        self.logger.info(f"[S1]     🔄 [{idx+1}/{total_confs}] DFT handoff on conformer: {conf_name}")
        with _STATE_MANAGER_LOCK:
            self.state_manager.mark_conformer_running(conf_name)
        self._emit_s1_progress("conformer_started", {"conformer": conf_name, "index": idx + 1, "total": total_confs})

        converged_this_conf = False
        conf_failed_reason = ""
        successful_record: Optional[Dict[str, Any]] = None

        for attempt in range(2):
            current_conf_name = f"{conf_name}_Res" if attempt > 0 else conf_name
            gjf_file = self.final_dft_dir / f"{current_conf_name}.gjf"
            log_file = self.final_dft_dir / f"{current_conf_name}.log"

            if attempt > 0:
                rescue_route = self.theory_opt.get('rescue_route') or build_gaussian_route_from_config(
                    self.config,
                    rescue=True,
                )
                gauss_int = GaussianInterface(
                    charge=self.theory_opt.get('charge', 0),
                    multiplicity=self.theory_opt.get('multiplicity', 1),
                    nprocshared=opt_nproc,
                    mem=opt_mem,
                    config=self.config,
                )
                gauss_int.write_input_file(
                    current_xyz_source,
                    gjf_file,
                    route=rescue_route,
                    title=f"{self.molecule_name}_{current_conf_name}_Rescue",
                )
            else:
                gauss_int = GaussianInterface(
                    charge=self.theory_opt.get('charge', 0),
                    multiplicity=self.theory_opt.get('multiplicity', 1),
                    nprocshared=opt_nproc,
                    mem=opt_mem,
                    config=self.config,
                )
                gauss_int.write_input_file(
                    current_xyz_source,
                    gjf_file,
                    route=self.theory_opt.get('route') or build_gaussian_route_from_config(
                        self.config,
                        rescue=False,
                    ),
                    title=f"{self.molecule_name}_{current_conf_name}",
                )

            opt_converged, next_xyz = self._run_gaussian_opt(
                gjf_file,
                log_file,
                xyz_source=current_xyz_source,
                attempt=attempt,
            )

            if opt_converged:
                with _STATE_MANAGER_LOCK:
                    self.state_manager.mark_opt_attempt(conf_name, attempt, "converged", log_file)
                self._emit_s1_progress("opt_converged", {"conformer": conf_name, "attempt": attempt})
            else:
                with _STATE_MANAGER_LOCK:
                    self.state_manager.mark_opt_attempt(
                        conf_name,
                        attempt,
                        "failed",
                        log_file,
                        note="gaussian_opt_not_converged",
                    )
                self._emit_s1_progress("opt_failed", {"conformer": conf_name, "attempt": attempt})

            if not opt_converged:
                if attempt == 0:
                    if next_xyz and next_xyz.exists():
                        current_xyz_source = next_xyz
                    conf_failed_reason = "opt_failed_attempt_0"
                    continue
                conf_failed_reason = "opt_failed_after_rescue"
                break

            final_coords, final_symbols, parse_error = LogParser.extract_last_converged_coords(
                log_file,
                engine_type='gaussian',
            )

            if final_coords is None:
                conf_failed_reason = "opt_parse_failed"
                self.logger.warning(f"      ⚠️  Failed to parse converged coordinates for {current_conf_name}: {parse_error}")
                if attempt == 0:
                    continue
                break

            if final_symbols is None:
                _, final_symbols = read_xyz(current_xyz_source)

            sp_in_file = self.final_dft_dir / f"{current_conf_name}_SP.inp"
            sp_out_file = self.final_dft_dir / f"{current_conf_name}_SP.out"
            sp_energy = self._run_orca_sp(
                final_coords,
                final_symbols,
                sp_in_file,
                sp_out_file,
                orca_interface=sp_interface,
            )

            if sp_energy is None:
                with _STATE_MANAGER_LOCK:
                    self.state_manager.mark_sp_result(
                        conf_name,
                        status="failed",
                        output_file=sp_out_file,
                        note="orca_sp_failed",
                    )
                self._emit_s1_progress("sp_failed", {"conformer": conf_name, "output": str(sp_out_file)})
                self.logger.warning(
                    "      ⚠️  ORCA SP failed after converged OPT; skipping rescue OPT retry and moving on"
                )
                conf_failed_reason = "sp_failed"
                break

            with _STATE_MANAGER_LOCK:
                self.state_manager.mark_sp_result(
                    conf_name,
                    status="completed",
                    output_file=sp_out_file,
                    sp_energy=sp_energy,
                )
            self._emit_s1_progress("sp_completed", {"conformer": conf_name, "energy_hartree": sp_energy})

            # Shermo removed — conformer ranking now uses electronic energy (sp_energy).
            # Gibbs/Boltzmann averaging at condition temperature is done in ConditionThermoEngine.
            successful_record = {
                "conf_name": conf_name,
                "name": current_conf_name,
                "source_index": idx,
                "g_energy": sp_energy,
                "g_used": sp_energy,  # Was Gibbs, now electronic energy for backward compat
                "sp_energy": sp_energy,
                "log_file": log_file,
            }
            with _STATE_MANAGER_LOCK:
                self.state_manager.mark_conformer_completed(
                    conf_name,
                    self._serialize_record_for_state(successful_record),
                )
            self._emit_s1_progress("conformer_completed", {"conformer": conf_name, "energy_hartree": sp_energy})
            converged_this_conf = True
            break

        if not converged_this_conf:
            failure_reason = conf_failed_reason or "conformer_failed"
            with _STATE_MANAGER_LOCK:
                self.state_manager.mark_conformer_failed(conf_name, failure_reason)
            self._emit_s1_progress(
                "conformer_failed",
                {"conformer": conf_name, "reason": failure_reason},
            )
            return {
                "conf_name": conf_name,
                "name": conf_name,
                "source_index": idx,
                "g_energy": float("nan"),
                "g_used": float("nan"),
                "sp_energy": float("nan"),
                "log_file": self.final_dft_dir / f"{conf_name}.log",
                "failed": True,
                "error": failure_reason,
            }

        if successful_record is None:
            raise RuntimeError(f"Conformer record missing for {conf_name}")

        return successful_record

    def _get_s1_parallel_config(self) -> Dict[str, Any]:
        intra_config = self.config.get("intra_reaction_parallel", {})
        if not isinstance(intra_config, dict):
            return {}
        s1_config = intra_config.get("s1", {})
        return s1_config if isinstance(s1_config, dict) else {}

    def _build_parallel_orca_sp_interface(self) -> ORCAInterface:
        s1_parallel_config = self._get_s1_parallel_config()
        return ORCAInterface(
            method=self.theory_sp.get("method", "M062X"),
            basis=self.theory_sp.get("basis", "def2-TZVPP"),
            aux_basis=self.theory_sp.get("aux_basis", getattr(self.orca_sp, "aux_basis", "def2/J")),
            nprocs=int(s1_parallel_config.get("sp_cores_per_job", getattr(self.orca_sp, "nprocs", 16))),
            maxcore=s1_parallel_config.get("sp_maxcore_per_job", getattr(self.orca_sp, "maxcore", None)),
            solvent=self.theory_sp.get("solvent", getattr(self.orca_sp, "solvent", "acetone")),
            route_extras=self.theory_sp.get("route_extras", ""),
            config=self.config,
        )

    def _finalize_conformer_results(
        self,
        records: List[Dict[str, Any]],
        total_confs: int,
        pm: Any,
    ) -> Tuple[Path, float, List[Dict[str, Any]]]:
        def _record_sort_key(item: Tuple[int, Dict[str, Any]]) -> int:
            source_index = item[1].get("source_index")
            if isinstance(source_index, int):
                return source_index
            if isinstance(source_index, str):
                try:
                    return int(source_index)
                except ValueError:
                    pass
            return item[0]

        successful_records: List[Dict[str, Any]] = []
        for record in records:
            g_used = record.get("g_used")
            sp_energy = record.get("sp_energy")
            try:
                if g_used is None or sp_energy is None:
                    continue
                if not math.isfinite(float(g_used)) or not math.isfinite(float(sp_energy)):
                    continue
            except (TypeError, ValueError):
                continue
            successful_records.append(record)

        if not successful_records:
            raise RuntimeError("所有 OPT-SP 循环均失败。")

        successful_records = [
            record
            for _, record in sorted(
                enumerate(successful_records),
                key=_record_sort_key,
            )
        ]

        best_weight = -1.0
        best_log: Optional[Path] = None
        best_sp_energy: Optional[float] = None

        # NOTE: After Shermo centralization, g_used is electronic energy (sp_energy).
        # Boltzmann weighting here is for conformer population estimation only.
        # Temperature-dependent Gibbs/Boltzmann is computed in ConditionThermoEngine.
        values = [float(record["g_used"]) for record in successful_records]
        min_value = min(values)
        rt = 0.0019872041 * float(self.thermo_config.get("temperature_k", 298.15))
        weights = [float(np.exp(-(value - min_value) / rt)) for value in values]
        total = float(np.sum(weights))
        if total > 0.0:
            weights = [float(weight / total) for weight in weights]
        else:
            weights = [0.0 for _ in weights]

        for record, weight in zip(successful_records, weights):
            record["weight"] = weight
            if weight > best_weight:
                best_weight = weight
                best_log = Path(str(record["log_file"]))
                best_sp_energy = float(record["sp_energy"])

        output_file = self.final_dft_dir / "conformer_thermo.csv"
        headers = [
            "name",
            "weight",
            "g_used",
            "g_sum",
            "g_conc",
            "h_sum",
            "u_sum",
            "s_total",
            "sp_energy",
            "log_file",
        ]
        lines = [",".join(headers)]
        for record in successful_records:
            lines.append(",".join([
                str(record.get("name", "")),
                f"{record.get('weight', 0.0):.6f}",
                f"{record.get('g_used', 0.0):.8f}",
                f"{record.get('g_sum', 0.0):.8f}",
                "" if record.get("g_conc") is None else f"{record.get('g_conc'):.8f}",
                f"{record.get('h_sum', 0.0):.8f}",
                "" if record.get("s_total") is None else f"{record.get('s_total'):.6f}",
                f"{record.get('sp_energy', 0.0):.8f}",
                str(record.get("log_file", "")),
            ]))
        output_file.write_text("\n".join(lines))

        self._write_conformer_energies_json(successful_records)

        if best_log is None or best_sp_energy is None:
            raise RuntimeError("无法确定最佳构象。")

        self.logger.info(f"  ✅ 最佳构象权重: {best_weight:.4f} ({best_log.name})")
        self._emit_s1_progress(
            "best_conformer_selected",
            {
                "conformer": best_log.stem,
                "weight": best_weight,
                "energy_hartree": best_sp_energy,
            },
        )
        if pm:
            pm.log_event(
                "S1",
                f"Best conformer selected after evaluating {len(successful_records)}/{total_confs} conformers: "
                f"{best_log.name} (weight={best_weight:.4f})",
            )

        return best_log, best_sp_energy, successful_records

    def _state_manager_call(self, method_name: str, *args: Any, **kwargs: Any) -> Any:
        with _STATE_MANAGER_LOCK:
            method = getattr(self.state_manager, method_name)
            return method(*args, **kwargs)

    def _emit_s1_progress(self, event: str, payload: Optional[Dict[str, Any]] = None) -> None:
        summary = self._state_manager_call("get_summary")
        message: Dict[str, Any] = {
            "schema": "s1_progress_v1",
            "event": event,
            "molecule": self.molecule_name,
            "state_file": str(self.state_manager.state_path),
            "summary": {
                "total_conformers": summary.get("total_conformers", 0),
                "completed": summary.get("completed", 0),
                "failed": summary.get("failed", 0),
                "running": summary.get("running", 0),
                "best_conformer": summary.get("best_conformer", ""),
                "global_min_energy": summary.get("global_min_energy"),
            },
        }
        if payload:
            message["payload"] = payload
        safe_message = self._json_safe(message)
        self.logger.info("S1_PROGRESS|" + json.dumps(safe_message, sort_keys=True, ensure_ascii=True, default=str, allow_nan=False))

    def _json_safe(self, value: Any) -> Any:
        if isinstance(value, float) and not math.isfinite(value):
            return None
        if isinstance(value, dict):
            return {str(k): self._json_safe(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [self._json_safe(item) for item in value]
        return value

    def _serialize_record_for_state(self, record: Dict[str, Any]) -> Dict[str, Any]:
        payload = dict(record)
        payload["log_file"] = str(record.get("log_file", ""))
        return payload

    def _restore_record_from_state(self, record: Dict[str, Any]) -> Dict[str, Any]:
        payload = dict(record)
        payload["log_file"] = Path(str(record.get("log_file", "")))
        return payload

    def _write_conformer_energies_json(self, records: List[Dict[str, Any]]) -> None:
        """Write conformer_energies.json for S4 Step1 activation extractor.

        Converts conformer thermo records to JSON format expected by
        step1_activation extractor (list of energies in kcal/mol).

        Args:
            records: List of conformer records with 'g_used' or 'sp_energy'
        """
        import json

        energies = []
        for record in records:
            # Prefer g_used (Gibbs free energy), fallback to sp_energy
            g_val = record.get("g_used")
            if g_val is not None:
                energies.append(float(g_val) * HARTREE_TO_KCAL)
            else:
                sp_val = record.get("sp_energy")
                if sp_val is not None:
                    energies.append(float(sp_val) * HARTREE_TO_KCAL)

        if energies:
            json_output = self.final_dft_dir / "conformer_energies.json"
            with open(json_output, 'w') as f:
                json.dump(energies, f, indent=2)
            self.logger.info(f"  ✅ Generated conformer_energies.json ({len(energies)} conformers)")

    def _write_crest_ensemble_energies_json(self) -> None:
        """Write crest/conformer_ensemble_energies.json from CREST GFN2 ensemble.

        Reads crest.energies (produced by CREST in the output directory),
        extracts relative GFN2-xTB energies for ALL conformers, and writes
        a structured JSON for the S4 conformer ensemble feature extractors.

        Search order:
        1. self.crest_dir / "crest.energies"           (single-stage)
        2. self.xtb_dir / "stage2_gfn2" / "crest.energies" (two-stage GFN2)
        3. self.xtb_dir / "stage1_gfn0" / "crest.energies" (two-stage GFN0 fallback)

        Output: crest/conformer_ensemble_energies.json with structure:
        {
            "source": "crest_gfn2",
            "n_conformers": N,
            "energies": [...],  // relative energies in kcal/mol
            "energy_unit": "kcal/mol",
            "temperature_K": 298.15
        }
        """
        import json

        search_dirs = [
            self.crest_dir,
            self.xtb_dir / "stage2_gfn2",
            self.xtb_dir / "stage1_gfn0",
        ]

        crest_energies_path = None
        for search_dir in search_dirs:
            candidate = search_dir / "crest.energies"
            if candidate.exists():
                crest_energies_path = candidate
                break

        if crest_energies_path is None:
            self.logger.debug("  crest.energies not found in any CREST output directory")
            return

        energies_kcal = []
        try:
            with open(crest_energies_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            # Format: <index> <energy_kcal>
                            energy = float(parts[1])
                            energies_kcal.append(energy)
                        except (ValueError, IndexError):
                            continue
                    elif len(parts) == 1:
                        try:
                            energy = float(parts[0])
                            energies_kcal.append(energy)
                        except ValueError:
                            continue
        except Exception as e:
            self.logger.warning(f"  Failed to parse crest.energies: {e}")
            return

        if not energies_kcal:
            self.logger.warning("  crest.energies found but no valid energies parsed")
            return

        # Convert to relative energies (subtract global minimum)
        min_energy = min(energies_kcal)
        relative_energies = [e - min_energy for e in energies_kcal]

        thermo_config = self.config.get('thermo', {}) or {}
        temperature_k = float(thermo_config.get('temperature_k', 298.15))

        payload = {
            "source": "crest_gfn2",
            "n_conformers": len(relative_energies),
            "energies": relative_energies,
            "energy_unit": "kcal/mol",
            "temperature_K": temperature_k,
            "source_file": str(crest_energies_path),
        }

        json_output = self.crest_dir / "conformer_ensemble_energies.json"
        with open(json_output, 'w') as f:
            json.dump(payload, f, indent=2)
        self.logger.info(
            f"  ✅ Generated crest/conformer_ensemble_energies.json "
            f"({len(relative_energies)} conformers, E_span={max(relative_energies):.3f} kcal/mol)"
        )

        csv_output = self.crest_dir / "conformer_thermo.csv"
        with open(csv_output, 'w', newline='') as f:
            f.write("index,energy_kcal_mol,rel_energy_kcal_mol\n")
            for i, (abs_energy, rel_energy) in enumerate(
                zip(energies_kcal, relative_energies), start=1
            ):
                f.write(f"{i},{abs_energy:.6f},{rel_energy:.6f}\n")
        self.logger.info(
            f"  ✅ Generated crest/conformer_thermo.csv "
            f"({len(relative_energies)} conformers, full ensemble)"
        )

    def _write_fastsp_ensemble_energies_json(self) -> None:
        """Write fastsp/conformer_ensemble_energies.json from fastSP screening.

        Scans all fastsp/conf_*/ directories for ORCA SP output files,
        extracts converged SP energies, and writes a structured JSON
        for the S4 conformer ensemble feature extractors.

        Output: fastsp/conformer_ensemble_energies.json with structure:
        {
            "source": "fastsp_dft",
            "n_conformers": N,
            "energies": [...],  // SP energies in kcal/mol (Hartree converted)
            "energy_unit": "kcal/mol",
            "temperature_K": 298.15
        }
        """
        import json
        import re

        if not self.fastsp_dir.exists():
            return

        # ORCA SP energy pattern: "FINAL SINGLE POINT ENERGY" followed by a float
        # Also try Gaussian SCF Done pattern as fallback
        energy_patterns = [
            re.compile(r"FINAL SINGLE POINT ENERGY\s+([-+]?\d+\.\d+)"),
            re.compile(r"SCF Done:\s+E\([^)]+\)\s*=\s*([-+]?\d+\.\d+)"),
        ]

        energies_hartree = []
        conf_dirs = sorted(
            [d for d in self.fastsp_dir.iterdir() if d.is_dir()],
            key=lambda d: d.name,
        )

        for conf_dir in conf_dirs:
            # Search for .log or .out files
            for log_path in sorted(conf_dir.glob("*.log")) + sorted(conf_dir.glob("*.out")):
                try:
                    content = log_path.read_text(encoding="utf-8", errors="ignore")
                except OSError:
                    continue

                for pattern in energy_patterns:
                    matches = pattern.findall(content)
                    if matches:
                        try:
                            # Use the last match (final energy)
                            energy = float(matches[-1])
                            energies_hartree.append(energy)
                            break
                        except (ValueError, IndexError):
                            continue
                else:
                    continue
                break  # Found energy for this conformer, move to next

        if not energies_hartree:
            self.logger.debug("  No fastSP energies extracted from fastsp/")
            return

        # Convert Hartree to kcal/mol
        relative_energies = [
            (e - min(energies_hartree)) * HARTREE_TO_KCAL
            for e in energies_hartree
        ]

        thermo_config = self.config.get('thermo', {}) or {}
        temperature_k = float(thermo_config.get('temperature_k', 298.15))

        payload = {
            "source": "fastsp_dft",
            "n_conformers": len(relative_energies),
            "energies": relative_energies,
            "energy_unit": "kcal/mol",
            "temperature_K": temperature_k,
            "screening_method": str(
                self.step1_config.get("fast_sp_profiles", {}).get(
                    "screening_sp", {}
                ).get("method", "r2SCAN-3c")
            ) if isinstance(self.step1_config, dict) else "r2SCAN-3c",
        }

        json_output = self.fastsp_dir / "conformer_ensemble_energies.json"
        with open(json_output, 'w') as f:
            json.dump(payload, f, indent=2)
        self.logger.info(
            f"  ✅ Generated fastsp/conformer_ensemble_energies.json "
            f"({len(relative_energies)} conformers from fastSP screening)"
        )

    def run_optimization_only(self, smiles: str) -> Tuple[Path, float]:
        """
        Run optimization-only workflow for rigid small molecules (Skip CREST).

        Flow:
        1. RDKit embed (generate 1 conformer)
        2. DFT Opt+Freq+SP (Standard/Rescue)
        3. Identify best result

        Args:
            smiles: SMILES string of the molecule

        Returns:
            Tuple of (path_to_global_min_xyz, energy_hartree)
        """
        import time
        start_time = time.time()
        self.logger.info(f"[S1] 🚀 启动小分子优化流程: {self.molecule_name}")
        self.logger.info(f"[S1]    SMILES: {smiles}")
        self.state_manager.start_run(smiles=smiles, two_stage_enabled=False)

        pm = get_progress_manager()
        if pm:
            pm.update_step("s1", description=f"S1: [{self.molecule_name}] Small molecule optimization")

        # Check for existing result (idempotency)
        global_min_path = self.molecule_dir / f"{self.molecule_name}_global_min.xyz"
        if self._is_valid_xyz(global_min_path):
            try:
                e_sp = self._extract_energy_from_xyz(global_min_path)
                self.state_manager.set_global_min(self.molecule_name, e_sp, global_min_path)
                self.state_manager.mark_run_complete()
                self.logger.info(f"[S1] ⏭️  Found existing global min: {global_min_path.name} (E={e_sp:.6f})")
                return global_min_path, e_sp
            except ValueError:
                self.logger.warning(f"[S1] Existing global min found but energy extract failed. Re-running.")

        # 1. RDKit Initialization (1 conformer)
        if pm:
            pm.update_step("s1", description=f"S1: [{self.molecule_name}] RDKit Initialization")
        self.logger.info("[S1]   [1/2] RDKit Initialization (Single Conformer)...")
        init_xyz = self._step_rdkit_embed(smiles, num_conf=1)

        # 2. DFT Optimization & SP
        # We treat the single RDKit structure as the "candidate list"
        self.logger.info("[S1]   [2/2] DFT Optimization & SP (Direct)...")
        candidates = [init_xyz]
        
        # Reuse the coupled DFT/SP logic
        # This handles Opt -> Freq -> SP and returns the best log/energy
        best_log, best_energy, _ = self._run_final_opt_sp_stage(candidates)

        # 3. Finalize
        duration = time.time() - start_time
        self.logger.info(f"[S1] ✅ Optimization-Only completed in {duration:.1f}s. Best E={best_energy:.6f}")


        # Create global min XYZ from the best log
        self._write_global_min_xyz(best_log, best_energy, global_min_path)
        self.state_manager.set_global_min(Path(best_log).stem.replace("_Res", ""), best_energy, global_min_path)
        self.state_manager.mark_run_complete()

        return global_min_path, best_energy

    def _run_final_opt_sp_stage(self, candidates: List[Path]) -> Tuple[Path, float, Dict[str, Any]]:
        selected = [Path(item).resolve() for item in candidates]
        candidate_set = candidate_set_from_paths(
            selected,
            source_protocol=self.protocol_spec.name,
            ranking_basis="input_order",
            ranking_unit="hartree",
            source_stage="post_ensemble",
        )
        stage_sequence = self.pipeline_executor.stage_sequence()
        self._emit_s1_progress(
            "final_opt_sp_started",
            {
                "protocol": self.protocol_spec.name,
                "selection_mode": self.protocol_spec.selection_mode,
                "pipeline_stages": stage_sequence,
                "selected_candidates": len(selected),
                "freq_requested": self.protocol_spec.freq_enabled,
            },
        )
        stage_result = self.pipeline_executor.execute_handoff(
            candidate_set=candidate_set,
            runner=self._run_dft_opt_sp_for_candidates,
        )
        self.last_final_opt_sp_meta = dict(stage_result.stage_meta)
        self._emit_s1_progress(
            "final_opt_sp_completed",
            {
                "protocol": self.protocol_spec.name,
                "energy_hartree": stage_result.best_energy,
                "stage_meta": self.last_final_opt_sp_meta,
            },
        )
        return stage_result.best_log, stage_result.best_energy, self.last_final_opt_sp_meta

    def _run_protocol_funnel(self, initial_xyz: Path) -> CandidateSet:
        from rph_core.steps.conformer_search.funnel import FunnelRunner

        runner = FunnelRunner(self, self.protocol_spec)
        candidate_set = runner.run(initial_xyz)
        self._write_funnel_candidate_summary(candidate_set)
        return candidate_set

    def _run_shared_dft_handoff(self, candidate_set: CandidateSet) -> Tuple[Path, float, Dict[str, Any]]:
        stage_result = self.pipeline_executor.execute_handoff(
            candidate_set=candidate_set,
            runner=self._run_dft_opt_sp_for_candidates,
        )
        stage_meta = dict(stage_result.stage_meta)
        stage_meta.setdefault("stages_executed", list(candidate_set.metadata.get("stages_executed", [])))
        stage_meta.setdefault("approx_thermo_applied", bool(candidate_set.metadata.get("approx_thermo_applied", False)))
        stage_meta.setdefault("boltzmann_cutoff", candidate_set.metadata.get("boltzmann_cutoff"))
        stage_meta.setdefault("window_count", int(candidate_set.metadata.get("window_count", len(candidate_set.candidates))))
        stage_meta.setdefault("gap_rank1_rank2", candidate_set.metadata.get("gap_rank1_rank2"))
        return stage_result.best_log, stage_result.best_energy, stage_meta

    def _write_funnel_candidate_summary(self, candidate_set: CandidateSet) -> Path:
        output = self.molecule_dir / "funnel_candidate_summary.json"
        payload = {
            "protocol": self.protocol_spec.name,
            "ranking_basis": candidate_set.ranking_basis,
            "ranking_unit": candidate_set.ranking_unit,
            "total_count": candidate_set.total_count,
            "selected_count_before_handoff": candidate_set.selected_count_before_handoff,
            "metadata": dict(candidate_set.metadata),
            "candidates": [
                {
                    "candidate_id": candidate.candidate_id,
                    "source_xyz": str(candidate.source_xyz),
                    "source_stage": candidate.source_stage,
                    "rank_initial": candidate.rank_initial,
                    "xtb_energy": candidate.xtb_energy,
                    "xtb_free_energy": candidate.xtb_free_energy,
                    "prescreen_energy": candidate.prescreen_energy,
                    "screening_energy": candidate.screening_energy,
                    "rerank_energy": candidate.rerank_energy,
                    "rerank_method": candidate.rerank_method,
                    "ranking_history": dict(candidate.ranking_history),
                    "included_by": list(candidate.included_by),
                    "within_window": candidate.within_window,
                    "cluster_id": candidate.cluster_id,
                }
                for candidate in candidate_set.candidates
            ],
        }
        output.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
        return output

    def _run_dft_opt_sp_for_candidates(self, candidates: List[ConformerCandidate]) -> Tuple[Path, float]:
        candidate_paths = [candidate.source_xyz for candidate in candidates]
        best_log, best_energy, _records = self._step_dft_opt_sp_coupled(candidate_paths)
        return best_log, best_energy



    def _run_gaussian_opt(
        self,
        gjf_file: Path,
        log_file: Path,
        xyz_source: Optional[Path] = None,
        attempt: int = 0
    ) -> Tuple[bool, Optional[Path]]:
        """
        [v3.0 FIX] Robust Gaussian OPT execution with absolute paths and input localization.

        Returns:
            Tuple[bool, Optional[Path]]:
                - bool: True if converged, False otherwise
                - Path: Coordinates for next attempt (None if success, or fallback path if failed)
        """
        import subprocess

        # [v3.0 FIX] Path Sanitization - Ensure absolute paths
        dft_dir_abs = self.final_dft_dir.resolve()
        gjf_path_abs = Path(gjf_file).resolve()
        log_path_abs = Path(log_file).resolve()
        xyz_source_abs = Path(xyz_source).resolve() if xyz_source else None

        # [v3.0 FIX] Input Localization - Copy files to execution directory
        local_gjf = dft_dir_abs / gjf_path_abs.name
        local_log = dft_dir_abs / log_path_abs.name

        # Copy input file if not already local
        if not local_gjf.exists() or local_gjf.resolve() != gjf_path_abs.resolve():
            try:
                shutil.copy(gjf_path_abs, local_gjf)
                self.logger.debug(f"  [Path-Fix] Copied {gjf_path_abs.name} to execution directory")
            except shutil.SameFileError:
                pass  # Already local

        exe_config = self.config.get('executables', {}).get('gaussian', {})
        use_wrapper = exe_config.get('use_wrapper', True)

        if use_wrapper:
            wrapper_path = exe_config.get('wrapper_path', './scripts/run_g16_worker.sh')
            if not Path(wrapper_path).is_absolute():
                wrapper_path = (Path.cwd() / wrapper_path).resolve()
            g16_cmd = str(wrapper_path)
            self.logger.info(f"      🔄 Running Gaussian OPT via wrapper: {local_gjf.name} (attempt {attempt})")
        else:
            g16_bin = exe_config.get('gaussian_bin', 'g16')
            g16_cmd = g16_bin
            self.logger.info(f"      🔄 Running Gaussian OPT: {local_gjf.name} (attempt {attempt})")

        self.logger.debug(f"  [Path-Fix] Executing in: {dft_dir_abs}")

        cmd = [g16_cmd, local_gjf.name, local_log.name]
        timeout_value = self.theory_opt.get('timeout')
        if timeout_value is None:
            timeout_cfg = self.config.get('optimization_control', {}).get('timeout', {})
            timeout_value = timeout_cfg.get('default_seconds', 21600)
        try:
            timeout_seconds = max(1, int(timeout_value))
        except (TypeError, ValueError):
            timeout_seconds = 21600

        try:
            run_with_timeout(cmd=cmd, timeout=timeout_seconds, cwd=dft_dir_abs)

            # Check if log file was created and is non-empty
            if not local_log.exists() or local_log.stat().st_size == 0:
                self.logger.warning(f"      ⚠️  OPT produced empty or missing log file")

                # [v3.0 FIX] Determine next coordinate source
                next_xyz = self._determine_rescue_xyz(xyz_source_abs, local_log, attempt)
                return False, next_xyz

            # Parse convergence status
            with open(local_log, 'r') as f:
                content = f.read()

            if "Normal termination" in content:
                # Check for imaginary frequencies
                import re
                freq_pattern = r"Frequencies --\s+([\d\-\.]+)"
                freqs = re.findall(freq_pattern, content)
                imag_count = sum(1 for f in freqs if float(f) < 0)

                if imag_count > 0:
                    self.logger.warning(f"      ⚠️  Found {imag_count} imaginary frequency(ies)")

                    # [v3.0 FIX] Try to rescue from last geometry even with imaginary frequencies
                    next_xyz = self._try_extract_rescue_coords(local_log, xyz_source_abs)
                    return False, next_xyz
                else:
                    self.logger.info(f"      ✓ OPT converged (no imaginary frequencies)")
                    return True, None
            else:
                self.logger.warning(f"      ⚠️  OPT did not terminate normally")

                # [v3.0 FIX] Try to rescue last coordinates from failed job
                next_xyz = self._try_extract_rescue_coords(local_log, xyz_source_abs)
                return False, next_xyz

        except QCTimeoutError:
            self.logger.error(f"      ❌ OPT timeout")

            # [v3.0 FIX] Try rescue from partial log
            next_xyz = self._try_extract_rescue_coords(local_log, xyz_source_abs)
            return False, next_xyz
        except subprocess.CalledProcessError as e:
            self.logger.error(f"      ❌ Gaussian execution failed (returncode {e.returncode})")

            # [v3.0 FIX] Determine next coordinate source on failure
            next_xyz = self._determine_rescue_xyz(xyz_source_abs, local_log, attempt)
            return False, next_xyz
        except Exception as e:
            self.logger.error(f"      ❌ OPT failed: {e}")

            # [v3.0 FIX] Return original source for retry
            return False, xyz_source_abs

    def _determine_rescue_xyz(
        self,
        original_xyz: Optional[Path],
        log_file: Path,
        attempt: int
    ) -> Optional[Path]:
        """
        [v3.0 FIX] Determine the coordinate source for rescue attempts.

        Strategy:
        - attempt 0: If job failed to start, use original CREST coordinates
        - attempt > 0: Try to rescue coordinates from failed log, else fallback to original
        """
        # Try to extract from log first (even if incomplete)
        rescue_xyz = self._try_extract_rescue_coords(log_file, original_xyz)

        if rescue_xyz is not None and rescue_xyz.exists():
            self.logger.info(f"      [Rescue] Using extracted coordinates from log: {rescue_xyz.name}")
            return rescue_xyz

        if original_xyz is not None:
            self.logger.info(f"      [Rescue] Falling back to original coordinates: {original_xyz.name}")

        return original_xyz

    def _try_extract_rescue_coords(
        self,
        log_file: Path,
        fallback_xyz: Optional[Path]
    ) -> Optional[Path]:
        if not log_file.exists():
            return fallback_xyz

        coords, symbols, error = LogParser.extract_last_converged_coords(log_file, 'gaussian')

        if coords is None or symbols is None:
            self.logger.debug(f"      [Rescue] Could not extract coordinates: {error}")
            return fallback_xyz

        # Create rescue XYZ file in dft directory
        rescue_xyz_path = self.final_dft_dir / f"{log_file.stem}_rescue.xyz"

        try:
            write_xyz(rescue_xyz_path, coords, symbols, title="Rescue Coordinates")
            self.logger.info(f"      [Rescue] Successfully extracted {len(coords)} atoms coordinates")
            return rescue_xyz_path
        except Exception as e:
            self.logger.warning(f"      [Rescue] Failed to save rescue coordinates: {e}")
            return fallback_xyz

    def _run_orca_sp(
        self,
        coords: np.ndarray,
        symbols: Optional[List[str]],
        inp_file: Path,
        out_file: Path,
        orca_interface: Optional[ORCAInterface] = None,
    ) -> Optional[float]:
        """
        [v3.0 FIX] Run ORCA single-point calculation with path localization.

        Args:
            coords: Coordinate array (N,3)
            symbols: Element symbols list (may be None from LogParser)
            inp_file: ORCA input file path
            out_file: ORCA output file path

        Returns:
            SP Energy in Hartree, or None if failed
        """
        # Handle None symbols
        if symbols is None or coords is None:
            self.logger.error("      ❌ Invalid coordinates or symbols")
            return None

        # Ensure symbols is List[str]
        symbols_list: List[str] = symbols
        sp_interface = orca_interface or self.orca_sp

        # [v3.0 FIX] Path localization
        dft_dir_abs = self.final_dft_dir.resolve()
        inp_path_abs = Path(inp_file).resolve()
        out_path_abs = Path(out_file).resolve()

        # Localize input file
        local_inp = dft_dir_abs / inp_path_abs.name
        local_out = dft_dir_abs / out_path_abs.name

        env = os.environ.copy()
        if hasattr(sp_interface, "_build_orca_runtime_env"):
            env = sp_interface._build_orca_runtime_env()

        nprocs_configured = int(getattr(sp_interface, "nprocs", 1) or 1)
        nprocs_for_sp = self._check_mpirun_compatibility(nprocs_configured, env)

        # Generate ORCA input directly (write to local path)
        self._generate_orca_sp_input(
            coords,
            symbols_list,
            local_inp,
            nprocs_for_sp,
            orca_interface=sp_interface,
        )

        if not local_inp.exists():
            self.logger.error(f"      ❌ ORCA input not created: {local_inp}")
            return None

        if sp_interface.orca_binary is None:
            self.logger.error("      ❌ ORCA binary not found")
            return None

        self.logger.info(f"      🔄 Running ORCA SP: {local_inp.name}")

        try:
            local_out = sp_interface._run_orca(
                local_inp,
                dft_dir_abs,
                timeout=self.theory_sp.get('timeout', 3600)
            )

            out_text = local_out.read_text(errors="ignore")
            if "ORCA TERMINATED NORMALLY" in out_text:
                energy_pattern = r"FINAL SINGLE POINT ENERGY\s+([\-\d\.]+)"
                match = __import__('re').search(energy_pattern, out_text)
                if match:
                    return float(match.group(1))
                energy_pattern = r"TOTAL SCF ENERGY\s+:\s+([\-\d\.]+)"
                match = __import__('re').search(energy_pattern, out_text)
                if match:
                    return float(match.group(1))
                self.logger.warning(f"      ⚠️  ORCA terminated normally but energy not found in output")
                return None

            self.logger.warning(f"      ⚠️  ORCA did not terminate normally")
            out_lines = out_text.split('\n')[-30:]
            err_summary = '\n'.join([f"         {line}" for line in out_lines if line.strip()])
            if err_summary:
                self.logger.debug(f"      ORCA output tail:\n{err_summary}")
            if "error" in out_text.lower() or "fatal" in out_text.lower():
                error_lines = [line for line in out_lines if any(kw in line.lower() for kw in ['error', 'fatal', 'abort'])]
                for line in error_lines[-3:]:
                    self.logger.error(f"         {line.strip()}")
            return None

        except subprocess.TimeoutExpired:
            self.logger.error(f"      ❌ ORCA SP timeout")
            return None
        except Exception as e:
            self.logger.error(f"      ❌ ORCA SP failed: {e}")
            return None

    def _run_fast_sp_profile(
        self,
        xyz_file: Path,
        *,
        profile_name: str,
        output_subdir: str,
    ) -> Optional[float]:
        profiles = self.step1_config.get("fast_sp_profiles", {}) if isinstance(self.step1_config, dict) else {}
        profile_cfg = profiles.get(profile_name, {}) if isinstance(profiles, dict) else {}
        profile_cfg = profile_cfg if isinstance(profile_cfg, dict) else {}
        defaults = {
            "prescreen_sp": {
                "method": "PBEh-3c",
                "basis": "",
                "aux_basis": "",
                "timeout": 1800,
            },
            "screening_sp": {
                "method": "r2SCAN-3c",
                "basis": "",
                "aux_basis": "",
                "timeout": 1800,
            },
        }
        merged_profile = {**defaults.get(profile_name, {}), **profile_cfg}
        if profile_name == "prescreen_sp":
            output_dir = self.prescan_dir / xyz_file.stem
        elif profile_name == "screening_sp":
            output_dir = self.fastsp_dir / xyz_file.stem
        else:
            output_dir = self.final_dft_dir / "fast_sp" / output_subdir / xyz_file.stem
        output_dir.mkdir(parents=True, exist_ok=True)
        fast_sp = ORCAInterface(
            method=str(merged_profile.get("method", self.theory_sp.get("method", "r2SCAN-3c"))),
            basis=str(merged_profile.get("basis", "")),
            aux_basis=str(merged_profile.get("aux_basis", "")),
            nprocs=int(merged_profile.get("nproc", self.theory_sp.get("nproc", 16))),
            solvent=str(merged_profile.get("solvent", self.theory_sp.get("solvent", self.solvent))),
            config=self.config,
        )
        result = fast_sp.single_point(
            xyz_file=Path(xyz_file).resolve(),
            output_dir=output_dir,
            timeout=int(merged_profile.get("timeout", 1800)),
            charge=int(self.theory_opt.get("charge", 0)),
            spin=int(self.theory_opt.get("multiplicity", 1)),
        )
        if result.converged and result.energy is not None:
            return float(result.energy)
        return None

    def _generate_orca_sp_input(
        self,
        coords: np.ndarray,
        symbols: List[str],
        inp_file: Path,
        nprocs: int,
        orca_interface: Optional[ORCAInterface] = None,
    ):
        """
        Generate ORCA SP input file from coordinates.

        Args:
            coords: Coordinate array (N, 3)
            symbols: Element symbols list
            inp_file: Output ORCA input file path
        """
        sp_interface = orca_interface or self.orca_sp
        route = sp_interface.render_simple_keywords(task_type="sp")
        rendered_blocks = sp_interface.render_blocks()
        block_text = ""
        if rendered_blocks:
            ordered_keys = ("basis", "method", "scf", "mdci")
            chunks = []
            for key in ordered_keys:
                value = rendered_blocks.get(key)
                if value:
                    chunks.append(value.rstrip())
            for key, value in rendered_blocks.items():
                if key in ordered_keys:
                    continue
                if value:
                    chunks.append(value.rstrip())
            if chunks:
                block_text = "\n" + "\n".join(chunks) + "\n"

        cpcm_block = ""
        if sp_interface.solvent and sp_interface.solvent.upper() != "NONE":
            cpcm_block = f"""
%cpcm
   smd true
   SMDsolvent "{sp_interface.solvent}"
end
"""

        # Coordinates as string
        coord_lines = []
        for symbol, coord in zip(symbols, coords):
            coord_lines.append(f"{symbol:2s} {coord[0]:15.10f} {coord[1]:15.10f} {coord[2]:15.10f}")

        coord_content = "\n".join(coord_lines)

        pal_block = f"%pal nprocs {nprocs} end\n" if nprocs > 1 else ""

        inp_content = f"""{route}
%maxcore {sp_interface.maxcore}
{pal_block}
        {cpcm_block}{block_text}
 * xyz 0 1
{coord_content}
 *
"""

        inp_file.write_text(inp_content)
        self.logger.debug(f"  ✓ ORCA SP input written: {inp_file}")

    def _check_mpirun_compatibility(self, nprocs: int, env: Optional[Dict[str, str]] = None) -> int:
        if nprocs <= 1:
            return 1

        search_path = (env or os.environ).get("PATH", "")
        mpirun_path = shutil.which("mpirun", path=search_path)
        if not mpirun_path:
            self.logger.warning("  mpirun not found, using single-core ORCA")
            return 1
        try:
            result = subprocess.run(
                [mpirun_path, "-np", "1", "/bin/true"],
                capture_output=True,
                text=True,
                timeout=5,
                env=env
            )
            if result.returncode != 0:
                stderr_text = (result.stderr or "").strip()
                stdout_text = (result.stdout or "").strip()
                probe_error = stderr_text or stdout_text or f"return code {result.returncode}"
                self.logger.warning(f"  mpirun probe failed ({probe_error}), using single-core ORCA")
                return 1
        except Exception:
            pass
        return nprocs

    def _is_valid_xyz(self, path: Path, min_size: int = 100) -> bool:
        """
        Validate XYZ file existence, size, and format.

        Args:
            path: Path to XYZ file
            min_size: Minimum file size in bytes (default 100 for single structures,
                    use 500+ for ensemble files like CREST outputs)

        Returns:
            True if file exists, has minimum size, and valid format
        """
        if not path.exists():
            return False

        if path.stat().st_size < min_size:
            return False

        try:
            with open(path, 'r') as f:
                first_line = f.readline().strip()
                return first_line.isdigit()
        except Exception:
            return False

    def _extract_energy_from_xyz(self, xyz_path: Path) -> float:
        """
        Extract energy from XYZ file with fail-fast validation.

        Strategy:
        1. Try to extract energy from XYZ comment line (2nd line)
           - Regex: r'[Ee](?:nergy)?[=:\\s]+(-?\\d+\\.?\\d*)'
        2. If not found, check companion log file (e.g., xyz_path.stem_SP.out)
           - Regex: r'FINAL SINGLE POINT ENERGY\\\\s+(-?\\d+\\.?\\d+)'

        Args:
            xyz_path: Path to XYZ file

        Returns:
            Energy in Hartree

        Raises:
            ValueError: If energy cannot be extracted (fail-fast - never returns 0.0)
        """
        import re

        if not self._is_valid_xyz(xyz_path):
            raise ValueError(f"CRITICAL: XYZ file invalid or empty: {xyz_path}")

        # Strategy A: Extract from XYZ comment line
        try:
            with open(xyz_path, 'r') as f:
                f.readline()  # Skip atom count
                comment_line = f.readline().strip()

                match = re.search(r'[Ee](?:nergy)?[=:\s]+(-?\d+\.?\d*)', comment_line)
                if match:
                    return float(match.group(1))
        except Exception as e:
            self.logger.debug(f"Energy extraction from XYZ comment failed: {e}")

        # Strategy B: Check companion log file
        log_patterns = [
            f"{xyz_path.stem}_SP.out",
            f"{xyz_path.stem}.log",
            f"{xyz_path.stem}.out"
        ]

        for log_pattern in log_patterns:
            log_path = xyz_path.parent / log_pattern
            if log_path.exists():
                try:
                    content = log_path.read_text()
                    # ORCA pattern
                    match = re.search(r'TOTAL SCF ENERGY\s+:\s+(-?\d+\.\d+)', content)
                    if match:
                        return float(match.group(1))

                    # Gaussian pattern
                    match = re.search(r'FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)', content)
                    if match:
                        return float(match.group(1))

                    # Generic pattern
                    match = re.search(r'SCF DONE:.*?(-?\d+\.\d+)', content)
                    if match:
                        return float(match.group(1))
                except Exception as e:
                    self.logger.debug(f"Energy extraction from log failed ({log_path}): {e}")

        # Fail Fast - do NOT return 0.0
        raise ValueError(
            f"CRITICAL: Cannot extract energy from {xyz_path}. "
            f"File may be corrupted or format unknown. "
            f"Tried: XYZ comment line and companion log files."
        )

    def _write_global_min_xyz(self, log_file: Path, energy: float, output_path: Path) -> None:
        from rph_core.utils.geometry_tools import LogParser
        
        coords, symbols, parse_error = LogParser.extract_last_converged_coords(
            log_file,
            engine_type='auto'
        )

        if coords is None:
            raise RuntimeError(f"Failed to extract coordinates from {log_file}")
        if symbols is None:
            _, symbols = read_xyz(log_file)

        if symbols is None:
            symbols = ["X"] * len(coords)

        write_xyz(output_path, coords, symbols, title=f"Global Min SP E={energy:.6f}")
        self.logger.info(f"[S1] ✓ Global minimum saved: {output_path.name}")
