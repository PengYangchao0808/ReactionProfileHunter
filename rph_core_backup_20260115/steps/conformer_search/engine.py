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
    S1_Product/[Molecule_Name]/
        ├── crest/          # CREST search results
        ├── dft/            # DFT OPT + SP calculations (no subdirectories)
        └── [Molecule_Name]_global_min.xyz

Logic ported from: Original_Eddition/Conf_Search_20251222/config/confsearch.lib.sh
"""

import logging
import shutil
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

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

logger = logging.getLogger(__name__)

HARTREE_TO_KCAL = 627.509


class ConformerEngine(LoggerMixin):
    """
    Unified Conformer Engine - v3.0 (Molecule-Autonomous Architecture)

    Key Features:
    - Flat directory structure per molecule
    - OPT-SP coupled execution loop
    - Robust log parsing for coordinate extraction
    - Naming-driven rescue workflow
    """

    def __init__(self, config: dict, work_dir: Path, molecule_name: str):
        """
        Initialize Conformer Engine with molecule-autonomous structure.

        [v3.0 FIX] All paths are resolved to absolute paths to eliminate ambiguity.

        Args:
            config: Full configuration dictionary
            work_dir: Base working directory (typically S1_Product/)
            molecule_name: Name of the molecule (creates subdirectory)
        """
        self.config = config
        self.molecule_name = molecule_name

        # [v3.0 FIX] Resolve all paths to absolute paths immediately
        work_dir = Path(work_dir).resolve()

        # Create molecule-autonomous directory structure
        self.molecule_dir = (work_dir / molecule_name).resolve()
        self.molecule_dir.mkdir(parents=True, exist_ok=True)

        self.crest_dir = (self.molecule_dir / "crest").resolve()
        self.dft_dir = (self.molecule_dir / "dft").resolve()
        self.crest_dir.mkdir(exist_ok=True)
        self.dft_dir.mkdir(exist_ok=True)

        self.logger.info(f"📁 Created molecule directory: {self.molecule_dir}")

        # Configuration Shortcuts
        self.step1_config = config.get('step1', {})
        self.crest_config = self.step1_config.get('crest', {})
        self.theory_opt = config.get('theory', {}).get('optimization', {})
        self.theory_sp = config.get('theory', {}).get('single_point', {})

        # Parameters
        self.energy_window = self.crest_config.get('energy_window', 6.0)  # kcal/mol
        self.dft_energy_window = self.step1_config.get('dft_energy_window', 3.0)  # kcal/mol
        self.rmsd_threshold = self.crest_config.get('rmsd_threshold', 0.125)  # Å
        self.solvent = self.crest_config.get('solvent', 'acetone')
        self.nproc = self.crest_config.get('threads', 8)
        self.max_conformers = self.crest_config.get('max_conformers', 20)

        # Initialize Interfaces
        self.xtb = XTBInterface(
            gfn_level=self.crest_config.get('gfn_level', 2),
            solvent=self.solvent,
            nproc=self.nproc,
            config=config
        )
        self.crest_interface = CRESTInterface(
            gfn_level=self.crest_config.get('gfn_level', 2),
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

        Args:
            smiles: SMILES string of the molecule

        Returns:
            Tuple[Path, float]: (Path to global min XYZ, Global min SP Energy in Hartree)
        """
        self.logger.info(f"🚀 Starting Unified Conformer Search for: {self.molecule_name}")

        # Step 1: RDKit Embedding
        initial_xyz = self._step_rdkit_embed(smiles)

        # Step 2: CREST Global Search
        crest_ensemble_file = self._step_crest_search(initial_xyz)

        # Step 3: Process Ensemble (Cluster & Filter)
        candidates = self._step_process_ensemble(crest_ensemble_file)
        self.logger.info(f"🔍 Selected {len(candidates)} conformers for DFT optimization.")

        if not candidates:
            raise RuntimeError("No valid conformers found after CREST processing.")

        # Step 4: DFT OPT-SP Coupled Loop
        best_xyz, min_energy = self._step_dft_opt_sp_coupled(candidates)

        # Step 5: Save Global Min
        global_min_path = self.molecule_dir / f"{self.molecule_name}_global_min.xyz"

        # Use final coordinates from the best conf
        coords, symbols, _ = LogParser.extract_last_converged_coords(
            best_xyz,
            engine_type='auto'
        )

        if coords is None:
            raise RuntimeError(f"Failed to extract coordinates from {best_xyz}")

        write_xyz(global_min_path, coords, symbols, title=f"Global Min SP E={min_energy:.6f}")
        self.logger.info(f"🏆 Global Minimum Found: {global_min_path}")
        self.logger.info(f"   SP Energy: {min_energy:.6f} Hartree")

        return global_min_path, min_energy

    def _step_rdkit_embed(self, smiles: str) -> Path:
        """Step 1: Generate initial 3D structure using RDKit."""
        self.logger.info("  [1/5] RDKit 3D Embedding...")
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        # Use ETKDGv3 for better initial geometries
        params = AllChem.ETKDGv3()
        params.useSmallRingTorsions = True
        AllChem.EmbedMolecule(mol, params)
        AllChem.MMFFOptimizeMolecule(mol)

        output_path = self.crest_dir / f"{self.molecule_name}_init.xyz"

        # Get coords and write
        conf = mol.GetConformer()
        coords = conf.GetPositions()
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

        # Ensure symbols is List[str]
        assert symbols is not None, "Symbols cannot be None after RDKit generation"

        write_xyz(output_path, coords, symbols, title="RDKit_Initial")
        return output_path

    def _step_crest_search(self, input_xyz: Path) -> Path:
        """Step 2: Run CREST global search."""
        self.logger.info("  [2/5] CREST Global Search...")
        best_xyz = self.crest_interface.run_conformer_search(input_xyz, self.crest_dir)

        # Look for ensemble file
        ensemble_file = self.crest_dir / "crest_conformers.xyz"
        if not ensemble_file.exists():
            self.logger.warning("crest_conformers.xyz not found, using best struct only.")
            # Create dummy ensemble with single structure
            shutil.copy(best_xyz, ensemble_file)

        return ensemble_file

    def _step_process_ensemble(self, ensemble_file: Path) -> List[Path]:
        """
        Step 3: Parse, Cluster, and Filter ensemble.

        Returns:
            List of Path objects to individual XYZ files of selected conformers.
        """
        self.logger.info("  [3/5] Processing Ensemble (Clustering & Filtering)...")

        # Split ensemble file
        conformers = self._split_xyz_ensemble(ensemble_file)
        self.logger.info(f"    - Found {len(conformers)} raw conformers from CREST.")

        if not conformers:
            return []

        # Filter by Energy (relative to min)
        parsed_confs = []
        min_e = float('inf')

        for conf_path in conformers:
            try:
                with open(conf_path, 'r') as f:
                    f.readline()  # atoms
                    title_line = f.readline().strip()

                e_match = __import__('re').search(r"(-?\d+\.\d+)", title_line)
                if e_match:
                    e = float(e_match.group(1))
                    if e < min_e:
                        min_e = e
                    parsed_confs.append({'path': conf_path, 'energy': e})
                else:
                    parsed_confs.append({'path': conf_path, 'energy': 0.0})
            except Exception:
                continue

        # Filter by energy window
        filtered_confs = []
        if min_e != float('inf'):
            for c in parsed_confs:
                rel_e = (c['energy'] - min_e) * HARTREE_TO_KCAL
                if rel_e <= self.energy_window:
                    filtered_confs.append(c)

            # Limit to max_conformers
            if len(filtered_confs) > self.max_conformers:
                filtered_confs = filtered_confs[:self.max_conformers]
        else:
            filtered_confs = parsed_confs

        self.logger.info(f"    - {len(filtered_confs)} conformers within {self.energy_window} kcal/mol window.")

        # Cluster (RMSD Pruning)
        final_selection = self._rmsd_clustering(filtered_confs, self.rmsd_threshold)

        return [item['path'] for item in final_selection]

    def _split_xyz_ensemble(self, ensemble_file: Path) -> List[Path]:
        """Splits a multi-structure XYZ file into individual files."""
        split_dir = self.crest_dir / "conf_split"
        split_dir.mkdir(exist_ok=True)

        paths = []
        with open(ensemble_file, 'r') as f:
            lines = f.readlines()

        i = 0
        count = 0
        while i < len(lines):
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

    def _rmsd_clustering(self, conf_list: List[dict], threshold: float) -> List[dict]:
        """
        Greedy RMSD clustering.

        Args:
            conf_list: [{'path': Path, 'energy': float}, ...] sorted by energy
            threshold: RMSD threshold in Angstroms

        Returns:
            Filtered list with unique conformers
        """
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
        """
        Step 4: DFT OPT-SP Coupled Execution Loop (Single-Conformer Atomic Operation).

        For each conformer:
        1. Run Gaussian OPT (with potential rescue)
        2. Extract last converged coordinates from .log
        3. Generate ORCA SP input (immediately)
        4. Run ORCA SP
        5. Record SP energy

        [NEW v3.0] 添加救援逻辑：使用 _Res 后缀重试失败的 OPT/SP

        Returns:
            Tuple[Path, float]: (Path to best SP output file, Best SP Energy)
        """
        self.logger.info("  [4/5] DFT OPT-SP Coupled Loop (Single-Conformer Atomic Operation)...")

        best_sp_energy = float('inf')
        best_sp_out = None

        for idx, xyz_file in enumerate(candidates):
            conf_name = f"conf_{idx:03d}"
            self.logger.info(f"    [{idx+1}/{len(candidates)}] 处理 {conf_name}...")

            # [v3.0 FIX] Initialize coordinate source for this conformer
            # Start with original CREST coordinates
            current_xyz_source = Path(xyz_file).resolve()

            # [v3.0 FIX] 尝试主流程和救援流程
            for attempt in range(2):  # 最多尝试2次（主流程 + 救援）
                # 确定文件名
                current_conf_name = f"{conf_name}_Res" if attempt > 0 else conf_name

                # Step A: Run Gaussian OPT
                gjf_file = self.dft_dir / f"{current_conf_name}.gjf"
                log_file = self.dft_dir / f"{current_conf_name}.log"

                # [v3.1 FIX] Generate GJF using GaussianInterface (unified logic)
                from rph_core.utils.qc_interface import GaussianInterface

                if attempt > 0:
                    self.logger.warning(f"      🔧 触发救援策略 (_Res)...")
                    rescue_route = self.theory_opt.get('rescue_route', '# B3LYP/def2-SVP Opt=CalcFC Freq NoEigenTest')
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
                        route=self.theory_opt.get('route', '# B3LYP/def2-SVP Opt=CalcFC Freq'),
                        title=f"{self.molecule_name}_{current_conf_name}"
                    )

                # [v3.0 FIX] Run OPT with coordinate source tracking
                opt_converged, next_xyz = self._run_gaussian_opt(
                    gjf_file,
                    log_file,
                    xyz_source=current_xyz_source,
                    attempt=attempt
                )

                if not opt_converged:
                    if attempt == 0:
                        self.logger.warning(f"      ⚠️  {conf_name} OPT 失败，尝试救援 (_Res)...")

                        # [v3.0 FIX] Update coordinate source for rescue attempt
                        if next_xyz and next_xyz.exists():
                            current_xyz_source = next_xyz
                            self.logger.debug(f"      [Rescue] Using extracted coordinates: {next_xyz.name}")
                        else:
                            # Fallback: keep using original coordinates
                            self.logger.debug(f"      [Rescue] Keeping original coordinates: {current_xyz_source.name}")

                        continue  # 尝试救援
                    else:
                        self.logger.error(f"      ❌ {conf_name} OPT 失败（救援也失败），跳过 SP...")
                        break  # 放弃这个构象

                # Step B: Extract last converged coordinates from .log
                final_coords, final_symbols, parse_error = LogParser.extract_last_converged_coords(
                    log_file,
                    engine_type='gaussian'
                )

                if final_coords is None:
                    if attempt == 0:
                        self.logger.warning(f"      ⚠️  无法提取坐标: {parse_error}，尝试救援 (_Res)...")
                        continue  # 尝试救援
                    else:
                        self.logger.error(f"      ❌ {conf_name} 坐标提取失败（救援也失败），跳过 SP...")
                        break  # 放弃这个构象

                # Ensure symbols is valid
                if final_symbols is None:
                    self.logger.warning(f"      ⚠️  未提取到符号，从 XYZ 文件读取")
                    coords, final_symbols = read_xyz(current_xyz_source)

                self.logger.info(f"      ✓ 提取到 {len(final_coords)} 个原子坐标（来自 {current_conf_name}.log）")

                # Step C: Run ORCA SP (immediately)
                sp_in_file = self.dft_dir / f"{current_conf_name}_SP.inp"
                sp_out_file = self.dft_dir / f"{current_conf_name}_SP.out"

                # [NEW] 救援时可能使用不同的 SP 参数
                if attempt > 0:
                    self.logger.info(f"      🔧 使用救援 SP 参数 (_Res)...")
                    # 可以在这里设置更宽松的 SCF 收敛标准等

                sp_energy = self._run_orca_sp(
                    final_coords,
                    final_symbols,
                    sp_in_file,
                    sp_out_file
                )

                if sp_energy is not None:
                    self.logger.info(f"      ✓ SP 能量: {sp_energy:.6f} Hartree ({current_conf_name})")

                    if sp_energy < best_sp_energy:
                        best_sp_energy = sp_energy
                        best_sp_out = sp_out_file
                        self.logger.info(f"      🏆 更新最佳能量: {best_sp_energy:.6f} Hartree")
                    break  # SP 成功，无需救援
                else:
                    if attempt == 0:
                        self.logger.warning(f"      ⚠️  {conf_name} SP 失败，尝试救援 (_Res)...")
                        # 使用更宽松的 SP 参数重试
                        continue  # 尝试救援
                    else:
                        self.logger.error(f"      ❌ {conf_name} SP 失败（救援也失败）")
                        break  # 放弃这个构象

        if best_sp_out is None:
            raise RuntimeError("所有 OPT-SP 循环均失败。")

        self.logger.info(f"  ✅ 最佳 SP 能量: {best_sp_energy:.6f} Hartree ({best_sp_out.name})")

        return best_sp_out, best_sp_energy

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
        original_xyz: Path,
        log_file: Path,
        attempt: int
    ) -> Path:
        """
        [v3.0 FIX] Determine the coordinate source for rescue attempts.

        Strategy:
        - attempt 0: If job failed to start, use original CREST coordinates
        - attempt > 0: Try to rescue coordinates from failed log, else fallback to original
        """
        # Try to extract from log first (even if incomplete)
        rescue_xyz = self._try_extract_rescue_coords(log_file, original_xyz)

        if rescue_xyz.exists():
            self.logger.info(f"      [Rescue] Using extracted coordinates from log: {rescue_xyz.name}")
            return rescue_xyz

        # Fallback to original
        self.logger.info(f"      [Rescue] Falling back to original coordinates: {original_xyz.name}")
        return original_xyz

    def _try_extract_rescue_coords(
        self,
        log_file: Path,
        fallback_xyz: Path
    ) -> Path:
        """
        [v3.0 FIX] Attempt to extract last coordinates from log file.

        If extraction fails, returns fallback_xyz path without creating file.
        """
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

        # Copy input if not already local
        if not local_inp.exists() or local_inp.resolve() != inp_path_abs.resolve():
            try:
                shutil.copy(inp_path_abs, local_inp)
                self.logger.debug(f"  [Path-Fix] Copied ORCA input to execution directory")
            except shutil.SameFileError:
                pass  # Already local

        # Setup ORCA environment
        setup_ld_library_path([])

        # Generate ORCA input directly (write to local path)
        self._generate_orca_sp_input(coords, symbols_list, local_inp)

        self.logger.info(f"      🔄 Running ORCA SP: {local_inp.name}")

        # Run ORCA in dft directory
        try:
            cmd = [self.orca_sp.orca_binary, local_inp.name]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=dft_dir_abs,
                timeout=self.theory_sp.get('timeout', 3600)
            )

            # Write output to local path
            with open(local_out, 'w') as f:
                f.write(result.stdout + result.stderr)

            # Parse energy
            if "ORCA TERMINATED NORMALLY" in result.stdout:
                energy_pattern = r"TOTAL SCF ENERGY\s+:\s+([\-\d\.]+)"
                match = __import__('re').search(energy_pattern, result.stdout)
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
