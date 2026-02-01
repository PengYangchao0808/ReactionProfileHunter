# pyright: ignore

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
        ├── xtb2/
        ├── cluster/
        ├── dft/
        └── [Molecule_Name]_global_min.xyz

Logic ported from: Original_Eddition/Conf_Search_20251222/config/confsearch.lib.sh
"""

# type: ignore

import logging
import shutil
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
from rph_core.utils.resource_utils import setup_ld_library_path
from rph_core.utils.qc_task_runner import QCTaskRunner
from rph_core.utils.optimization_config import build_gaussian_route_from_config
from rph_core.utils.isostat_runner import run_isostat
from rph_core.utils.shermo_runner import run_shermo


logger = logging.getLogger(__name__)

HARTREE_TO_KCAL = 627.509


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
            ├── xtb2/
            │   ├── stage1_gfn0/          # GFN0 coarse search
            │   │   ├── crest_conformers.xyz
            │   │   └── cluster/
            │   │       ├── cluster.xyz   # Representatives for Stage 2
            │   │       └── isostat.log
            │   ├── stage2_gfn2/          # GFN2 refinement
            │   │   ├── crest_ensemble.xyz
            │   │   └── cluster/
            │   │       ├── cluster.xyz   # Final ensemble
            │   │       └── isostat.log
            │   └── ensemble.xyz          # Symlink to stage2/cluster/cluster.xyz
            ├── cluster/                  # Legacy (for single-stage)
            └── dft/                      # DFT OPT-SP jobs

    Note on enabled semantics:
    - ConformerEngine does not have an 'enabled' flag.
    - When instantiated, it always operates in real engine mode.
    - To disable conformer search: do not instantiate ConformerEngine.
    - Placeholder mode is not implemented; use conditional instantiation instead.
    """
    def __init__(self, config: Dict[str, Any], work_dir: Path, molecule_name: str):
        super().__init__()
        self.config = config
        self.molecule_name = molecule_name

        # [v3.0 FIX] Resolve all paths to absolute paths immediately
        work_dir = Path(work_dir).resolve()

        # Create molecule-autonomous directory structure
        self.molecule_dir = (work_dir / molecule_name).resolve()
        self.molecule_dir.mkdir(parents=True, exist_ok=True)

        self.crest_dir = (self.molecule_dir / "xtb2").resolve()
        self.cluster_dir = (self.molecule_dir / "cluster").resolve()
        self.dft_dir = (self.molecule_dir / "dft").resolve()
        self.crest_dir.mkdir(exist_ok=True)
        self.cluster_dir.mkdir(exist_ok=True)
        self.dft_dir.mkdir(exist_ok=True)

        self.logger.info(f"📁 Created molecule directory: {self.molecule_dir}")

        self.step1_config = config.get('step1', {})
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

    def run(self, smiles: str) -> Tuple[Path, float]:
        """
        Execute full conformer search workflow (OPT-SP coupled).

        Automatically dispatches to single-stage or two-stage workflow based on config.

        Args:
            smiles: SMILES string of the molecule

        Returns:
            Tuple[Path, float]: (Path to global min XYZ, Global min SP Energy in Hartree)
        """
        self.logger.info(f"🚀 Starting Unified Conformer Search for: {self.molecule_name}")
        mode = "two-stage (GFN0→GFN2)" if self.two_stage_enabled else "single-stage (GFN2)"
        self.logger.info(f"   Mode: {mode}")

        # Step 0: Final Output Check (Idempotency)
        global_min_path = self.molecule_dir / f"{self.molecule_name}_global_min.xyz"

        if self._is_valid_xyz(global_min_path):
            self.logger.info(f"⏭️ Found {global_min_path.name}, skipping conformer search")

            # Extract energy - let ValueError propagate if extraction fails
            energy = self._extract_energy_from_xyz(global_min_path)
            self.logger.info(f"   Extracted energy: {energy:.6f} Hartree")

            return global_min_path, energy

        # Step 1: RDKit Embedding
        initial_xyz = self._step_rdkit_embed(smiles)

        # Step 2: CREST Search (Single or Two-Stage)
        if self.two_stage_enabled and self.stage1_enabled and self.stage2_enabled:
            self.logger.info("🔀 Executing two-stage CREST workflow (GFN0→GFN2)...")
            crest_ensemble_file = self._step_two_stage_crest(initial_xyz)
        else:
            self.logger.info("📍 Executing single-stage CREST workflow (GFN2)...")
            crest_ensemble_file = self._step_crest_search(initial_xyz)

        # Step 3: Process Ensemble (Cluster & Filter)
        candidates = self._step_process_ensemble(crest_ensemble_file)
        self.logger.info(f"🔍 Selected {len(candidates)} conformers for DFT optimization.")

        if not candidates:
            raise RuntimeError("No valid conformers found after CREST processing.")

        # Step 4: DFT OPT-SP Coupled Loop
        best_log, min_energy = self._step_dft_opt_sp_coupled(candidates)

        # Step 5: Save Global Min
        global_min_path = self.molecule_dir / f"{self.molecule_name}_global_min.xyz"

        # Use final coordinates from the best conf
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
        self.logger.info(f"🏆 Global Minimum Found: {global_min_path}")
        self.logger.info(f"   SP Energy: {min_energy:.6f} Hartree")

        return global_min_path, min_energy

    def _step_rdkit_embed(self, smiles: str) -> Path:
        self.logger.info("  [1/5] RDKit 3D Embedding...")
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        params = rdDistGeom.ETKDG()
        rdDistGeom.EmbedMolecule(mol, params)
        rdForceFieldHelpers.MMFFOptimizeMolecule(mol)

        output_path = self.crest_dir / f"{self.molecule_name}_init.xyz"

        conf = mol.GetConformer()
        coords = conf.GetPositions()
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

        assert symbols is not None, "Symbols cannot be None after RDKit generation"

        write_xyz(output_path, coords, symbols, title="RDKit_Initial")
        return output_path

    def _step_crest_search(self, input_xyz: Path) -> Path:
        """Single-stage CREST search (GFN2 only) - backward compatible."""
        self.logger.info("  [2/5] CREST Global Search (Single-Stage)...")

        crest_ensemble_file = self.crest_dir / "ensemble.xyz"

        if self._is_valid_xyz(crest_ensemble_file, min_size=500):
            self.logger.info("⏭️ Found CREST cache, skipping search")
            return crest_ensemble_file

        best_xyz = self.crest_interface.run_conformer_search(input_xyz, self.crest_dir)

        ensemble_file = self.crest_dir / "crest_conformers.xyz"
        if ensemble_file.exists():
            shutil.copy(ensemble_file, crest_ensemble_file)
            return crest_ensemble_file

        self.logger.warning("crest_conformers.xyz not found, using best struct only.")
        shutil.copy(best_xyz, crest_ensemble_file)
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
        self.logger.info("  [2a/5] Stage 1: CREST GFN0 Conformer Search (Fast Sampling)...")

        stage1_dir = self.crest_dir / "stage1_gfn0"
        stage1_dir.mkdir(exist_ok=True)

        gfn0_ensemble = self._run_crest_stage(
            input_xyz=input_xyz,
            output_dir=stage1_dir,
            gfn_level=self.stage1_gfn_level,
            stage_name="GFN0",
            additional_flags=self.stage1_crest_flags
        )

        # ==============================================================
        # Stage 1: ISOSTAT Clustering
        # ==============================================================
        if self.stage1_run_clustering:
            self.logger.info("  [2b/5] Stage 1: ISOSTAT Clustering...")
            gfn0_clustered = self._run_isostat_clustering(
                input_ensemble=gfn0_ensemble,
                output_dir=stage1_dir,
                energy_window=self.stage1_energy_window,
                gdis=self.stage1_isostat_gdis,
                edis=self.stage1_isostat_edis
            )
        else:
            gfn0_clustered = gfn0_ensemble

        # ==============================================================
        # Stage 2: GFN2 Batch Optimization
        # ==============================================================
        self.logger.info("  [2c/5] Stage 2: CREST GFN2 Batch Optimization (Refinement)...")

        stage2_dir = self.crest_dir / "stage2_gfn2"
        stage2_dir.mkdir(exist_ok=True)

        gfn2_ensemble = self._run_crest_batch_optimization(
            input_xyz=gfn0_clustered,
            output_dir=stage2_dir,
            gfn_level=self.stage2_gfn_level,
            stage_name="GFN2",
            additional_flags=self.stage2_crest_flags
        )

        # ==============================================================
        # Stage 2: ISOSTAT Clustering
        # ==============================================================
        if self.stage2_run_clustering:
            self.logger.info("  [2d/5] Stage 2: ISOSTAT Clustering...")
            gfn2_clustered = self._run_isostat_clustering(
                input_ensemble=gfn2_ensemble,
                output_dir=stage2_dir,
                energy_window=self.stage2_energy_window,
                gdis=self.stage2_isostat_gdis,
                edis=self.stage2_isostat_edis
            )
        else:
            gfn2_clustered = gfn2_ensemble

        # ==============================================================
        # Final: Copy to standard location for downstream processing
        # ==============================================================
        self.logger.info("  [2e/5] Finalizing ensemble for DFT processing...")

        final_ensemble = self.crest_dir / "ensemble.xyz"
        shutil.copy(gfn2_clustered, final_ensemble)

        self.logger.info(f"  ✓ Two-stage CREST complete: {final_ensemble}")

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
        cluster_dir.mkdir(exist_ok=True)

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
        self.logger.info("  [3/5] Processing Ensemble (Clustering & Filtering)...")

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

        conformers = self._split_xyz_ensemble(result.cluster_xyz, self.dft_dir, limit=n_calc)
        self.logger.info(f"    - Found {len(conformers)} conformers from isostat (Limit: {n_calc}).")

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

    def _step_dft_opt_sp_coupled(self, candidates: List[Path]) -> Tuple[Path, float]:
        self.logger.info("  [4/5] DFT OPT-SP Coupled Loop (Single-Conformer Atomic Operation)...")

        best_weight = -1.0
        best_log = None
        best_sp_energy = None
        records: List[Dict[str, Any]] = []

        for idx, xyz_file in enumerate(candidates):
            conf_name = f"conf_{idx:03d}"
            self.logger.info(f"    [{idx+1}/{len(candidates)}] 处理 {conf_name}...")

            current_xyz_source = Path(xyz_file).resolve()

            for attempt in range(2):
                current_conf_name = f"{conf_name}_Res" if attempt > 0 else conf_name

                gjf_file = self.dft_dir / f"{current_conf_name}.gjf"
                log_file = self.dft_dir / f"{current_conf_name}.log"

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

                if not opt_converged:
                    if attempt == 0:
                        if next_xyz and next_xyz.exists():
                            current_xyz_source = next_xyz
                        continue
                    break

                final_coords, final_symbols, parse_error = LogParser.extract_last_converged_coords(
                    log_file,
                    engine_type='gaussian'
                )

                if final_coords is None:
                    if attempt == 0:
                        continue
                    break

                if final_symbols is None:
                    _, final_symbols = read_xyz(current_xyz_source)

                sp_in_file = self.dft_dir / f"{current_conf_name}_SP.inp"
                sp_out_file = self.dft_dir / f"{current_conf_name}_SP.out"

                sp_energy = self._run_orca_sp(
                    final_coords,
                    final_symbols,
                    sp_in_file,
                    sp_out_file
                )

                if sp_energy is None:
                    if attempt == 0:
                        continue
                    break

                shermo_out = self.dft_dir / f"{current_conf_name}_Shermo.sum"
                thermo = run_shermo(
                    shermo_bin=self.shermo_bin,
                    freq_output=log_file,
                    sp_energy=sp_energy,
                    output_file=shermo_out,
                    temperature_k=self.thermo_config.get('temperature_k', 298.15),
                    pressure_atm=self.thermo_config.get('pressure_atm', 1.0),
                    scl_zpe=self.thermo_config.get('scl_zpe', 0.9905),
                    ilowfreq=self.thermo_config.get('ilowfreq', 2),
                    imagreal=self.thermo_config.get('imagreal', 0),
                    conc=self.thermo_config.get('conc')
                )

                g_used = thermo.g_conc if thermo.g_conc is not None else thermo.g_sum

                records.append({
                    "name": current_conf_name,
                    "g_used": g_used,
                    "g_sum": thermo.g_sum,
                    "g_conc": thermo.g_conc,
                    "h_sum": thermo.h_sum,
                    "u_sum": thermo.u_sum,
                    "s_total": thermo.s_total,
                    "sp_energy": sp_energy,
                    "log_file": log_file
                })
                break

        if not records:
            raise RuntimeError("所有 OPT-SP 循环均失败。")

        values = [record["g_used"] for record in records]
        min_value = min(values)
        rt = 0.0019872041 * float(self.thermo_config.get("temperature_k", 298.15))
        weights = [float(np.exp(-(value - min_value) / rt)) for value in values]
        total = float(np.sum(weights))
        if total > 0.0:
            weights = [float(weight / total) for weight in weights]
        else:
            weights = [0.0 for _ in weights]

        for record, weight in zip(records, weights):
            record["weight"] = weight
            if weight > best_weight:
                best_weight = weight
                best_log = record["log_file"]
                best_sp_energy = record["sp_energy"]

        output_file = self.dft_dir / "conformer_thermo.csv"
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
            "log_file"
        ]
        lines = [",".join(headers)]
        for record in records:
            lines.append(",".join([
                str(record.get("name", "")),
                f"{record.get('weight', 0.0):.6f}",
                f"{record.get('g_used', 0.0):.8f}",
                f"{record.get('g_sum', 0.0):.8f}",
                "" if record.get("g_conc") is None else f"{record.get('g_conc'):.8f}",
                f"{record.get('h_sum', 0.0):.8f}",
                "" if record.get("s_total") is None else f"{record.get('s_total'):.6f}",
                f"{record.get('sp_energy', 0.0):.8f}",
                str(record.get("log_file", ""))
            ]))
        output_file.write_text("\n".join(lines))

        if best_log is None or best_sp_energy is None:
            raise RuntimeError("无法确定最佳构象。")

        self.logger.info(f"  ✅ 最佳构象权重: {best_weight:.4f} ({best_log.name})")

        return best_log, best_sp_energy


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
        import os

        # [v3.0 FIX] Path Sanitization - Ensure absolute paths
        dft_dir_abs = self.dft_dir.resolve()
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

        if use_wrapper:
            cmd = f'"{g16_cmd}" {local_gjf.name} {local_log.name}'
        else:
            cmd = f"{g16_cmd} {local_gjf.name} {local_log.name}"

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                cwd=dft_dir_abs,
                timeout=self.theory_opt.get('timeout', 3600),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # Check shell-level errors
            if result.returncode != 0:
                self.logger.error(f"      ❌ Gaussian execution failed (returncode {result.returncode})")
                self.logger.debug(f"      stderr: {result.stderr}")

                # [v3.0 FIX] Determine next coordinate source on failure
                next_xyz = self._determine_rescue_xyz(xyz_source_abs, local_log, attempt)
                return False, next_xyz

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

        except subprocess.TimeoutExpired:
            self.logger.error(f"      ❌ OPT timeout")

            # [v3.0 FIX] Try rescue from partial log
            next_xyz = self._try_extract_rescue_coords(local_log, xyz_source_abs)
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
        rescue_xyz_path = self.dft_dir / f"{log_file.stem}_rescue.xyz"

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
        out_file: Path
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
        import subprocess
        import os

        # Handle None symbols
        if symbols is None or coords is None:
            self.logger.error("      ❌ Invalid coordinates or symbols")
            return None

        # Ensure symbols is List[str]
        symbols_list: List[str] = symbols

        # [v3.0 FIX] Path localization
        dft_dir_abs = self.dft_dir.resolve()
        inp_path_abs = Path(inp_file).resolve()
        out_path_abs = Path(out_file).resolve()

        # Localize input file
        local_inp = dft_dir_abs / inp_path_abs.name
        local_out = dft_dir_abs / out_path_abs.name

        # Setup ORCA environment
        setup_ld_library_path([])

        # Generate ORCA input directly (write to local path)
        self._generate_orca_sp_input(coords, symbols_list, local_inp)

        if not local_inp.exists():
            self.logger.error(f"      ❌ ORCA input not created: {local_inp}")
            return None

        if self.orca_sp.orca_binary is None:
            self.logger.error("      ❌ ORCA binary not found")
            return None

        self.logger.info(f"      🔄 Running ORCA SP: {local_inp.name}")

        try:
            cmd = [str(self.orca_sp.orca_binary), local_inp.name]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=dft_dir_abs,
                timeout=self.theory_sp.get('timeout', 3600)
            )

            if not local_out.exists() or local_out.stat().st_size == 0:
                local_out.write_text(result.stdout + result.stderr)

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

            self.logger.warning(f"      ⚠️  ORCA did not terminate normally")
            return None

        except subprocess.TimeoutExpired:
            self.logger.error(f"      ❌ ORCA SP timeout")
            return None
        except Exception as e:
            self.logger.error(f"      ❌ ORCA SP failed: {e}")
            return None

    def _generate_orca_sp_input(
        self,
        coords: np.ndarray,
        symbols: List[str],
        inp_file: Path
    ):
        """
        Generate ORCA SP input file from coordinates.

        Args:
            coords: Coordinate array (N, 3)
            symbols: Element symbols list
            inp_file: Output ORCA input file path
        """
        route = f"! {self.orca_sp.method} {self.orca_sp.basis} {self.orca_sp.aux_basis} RIJCOSX tightSCF"
        route += " noautostart miniprint nopop"

        cpcm_block = ""
        if self.orca_sp.solvent and self.orca_sp.solvent.upper() != "NONE":
            cpcm_block = f"""
%cpcm
   smd true
   SMDsolvent "{self.orca_sp.solvent}"
end
"""

        # Coordinates as string
        coord_lines = []
        for symbol, coord in zip(symbols, coords):
            coord_lines.append(f"{symbol:2s} {coord[0]:15.10f} {coord[1]:15.10f} {coord[2]:15.10f}")

        coord_content = "\n".join(coord_lines)

        inp_content = f"""{route}
%maxcore {self.orca_sp.maxcore}
%pal nprocs {self.orca_sp.nprocs} end
{cpcm_block}
 * xyz 0 1
{coord_content}
 *
"""

        inp_file.write_text(inp_content)
        self.logger.debug(f"  ✓ ORCA SP input written: {inp_file}")

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
