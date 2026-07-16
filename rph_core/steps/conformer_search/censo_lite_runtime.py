"""Minimal runtime primitives for the fixed V4 CENSO-LITE protocol."""

from __future__ import annotations

import logging
import re
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom

from rph_core.utils.constants import HARTREE_TO_KCAL
from rph_core.utils.file_io import write_xyz
from rph_core.utils.qc_interface import CRESTInterface, XTBInterface
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.steps.conformer_search.xtb_thermo import XTBThermoResult

logger = logging.getLogger(__name__)


class CrestEnergyParseError(RuntimeError):
    """Raised when xTB/CREST energies cannot be parsed from output files."""


class CensoLiteRuntime:
    """CENSO-LITE execution primitives without the V3 ConformerEngine."""

    def __init__(self, config: Dict[str, Any], work_dir: Path, molecule_name: str):
        self.config = config
        self.work_dir = Path(work_dir)
        self.molecule_name = molecule_name
        self._event_callback = None
        self.molecule_dir = self.work_dir / molecule_name
        self.crest_dir = self.molecule_dir / "crest"
        self.raw_dir = self.molecule_dir / "candidates_raw"
        self.ranking_dir = self.molecule_dir / "ranking"
        self.mrrho_dir = self.molecule_dir / "mrrho"
        for directory in (self.crest_dir, self.raw_dir, self.ranking_dir, self.mrrho_dir):
            directory.mkdir(parents=True, exist_ok=True)

    def set_event_callback(self, callback) -> None:
        self._event_callback = callback

    def _emit(self, event: str, **fields: Any) -> None:
        if self._event_callback is None:
            return
        try:
            self._event_callback(event, fields)
        except Exception as exc:  # pragma: no cover - UI isolation
            logger.warning("Ignoring CENSO-LITE runtime UI callback failure: %s", exc)

    def embed(self, smiles: str) -> Path:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES for CENSO-LITE: {smiles}")
        mol = Chem.AddHs(mol)
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = 20260711
        params.useRandomCoords = True
        conformer_id = rdDistGeom.EmbedMolecule(mol, params)
        if conformer_id < 0:
            raise RuntimeError("RDKit failed to generate an initial 3D conformer")
        coordinates = np.asarray(mol.GetConformer(conformer_id).GetPositions(), dtype=float)
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        output = self.molecule_dir / "initial.xyz"
        write_xyz(output, coordinates, symbols, title="CENSO-LITE initial geometry")
        return output

    def crest_search(self, input_xyz: Path) -> Path:
        cfg = dict(self.config.get("step1", {}).get("censo_lite", {}).get("crest", {}) or {})
        resources = dict(self.config.get("resources", {}) or {})
        total_cores = max(1, int(resources.get("nproc", 1)))
        requested_nproc = cfg.get("nproc", "auto")
        if requested_nproc is None or str(requested_nproc).strip().lower() == "auto":
            nproc = total_cores
        else:
            requested = max(1, int(requested_nproc))
            nproc = min(requested, total_cores)
            if requested > total_cores:
                logger.warning(
                    "Clamping CREST nproc=%d to global S1 core budget=%d",
                    requested,
                    total_cores,
                )
        search_mode = cfg.get("search_mode")
        ewin = cfg.get("energy_window_kcal")
        interface = CRESTInterface(
            gfn_level=int(cfg.get("gfn_level", 2)),
            solvent=cfg.get("solvent"),
            nproc=nproc,
            config=self.config,
            additional_flags=str(cfg.get("additional_flags", "")),
            search_mode=search_mode,
            energy_window_kcal=float(ewin) if ewin is not None else None,
        )
        logger.debug(
            "CREST: GFN%d/ALPB(%s) search_mode=%s ewin=%s nproc=%d",
            int(cfg.get("gfn_level", 2)),
            cfg.get("solvent"),
            search_mode,
            ewin,
            nproc,
        )
        interface.crest_timeout_seconds = int(cfg.get("timeout", 21600))
        ensemble = self.crest_dir / "crest_conformers.xyz"
        if ensemble.exists() and ensemble.stat().st_size > 0:
            return ensemble
        result = interface.run_conformer_search(
            Path(input_xyz).resolve(),
            self.crest_dir,
            gfn_override=int(cfg.get("gfn_level", 2)),
            additional_flags=str(cfg.get("additional_flags", "")),
            search_mode=search_mode,
            energy_window_kcal=float(ewin) if ewin is not None else None,
        )
        if ensemble.exists():
            return ensemble
        result = Path(result)
        if not result.exists():
            raise RuntimeError("CREST produced no ensemble or best conformer")
        shutil.copy2(result, ensemble)
        return ensemble

    def split_ensemble(self, ensemble: Path) -> List[Path]:
        lines = Path(ensemble).read_text(encoding="utf-8", errors="ignore").splitlines()
        paths: List[Path] = []
        index = 0
        while index < len(lines):
            if not lines[index].strip():
                index += 1
                continue
            try:
                atom_count = int(lines[index].strip())
            except ValueError:
                index += 1
                continue
            end = index + 2 + atom_count
            if end > len(lines):
                break
            block = lines[index:end]
            output = self.raw_dir / f"conf_{len(paths) + 1:04d}.xyz"
            output.write_text("\n".join(block) + "\n", encoding="utf-8")
            paths.append(output)
            index = end
        return paths

    @staticmethod
    def extract_energy(xyz_path: Path) -> float:
        """Parse the xTB electronic energy (Hartree) from an XYZ comment line.

        Raises ``CrestEnergyParseError`` if parsing fails — never returns 0.0
        as a silent fallback.
        """
        try:
            comment = Path(xyz_path).read_text(encoding="utf-8", errors="ignore").splitlines()[1]
            number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
            match = re.search(
                rf"(?:energy|E)\s*[:=]?\s*({number})",
                comment,
                re.IGNORECASE,
            )
            if match is None:
                match = re.fullmatch(rf"\s*({number})\s*", comment)
            if match is None:
                raise ValueError("energy token not found")
            energy = float(match.group(1).replace("D", "E").replace("d", "e"))
            if not np.isfinite(energy):
                raise ValueError(f"non-finite energy: {energy}")
            return energy
        except (AttributeError, OSError, IndexError, ValueError) as exc:
            raise CrestEnergyParseError(
                f"Could not parse xTB energy from {xyz_path}: {exc}"
            ) from exc

    @staticmethod
    def read_crest_relative_energies(crest_dir: Path) -> Optional[List[float]]:
        """Read ``crest.energies`` (relative energies in kcal/mol).

        Returns ``None`` if the file does not exist.  The values are
        sorted by energy ascending and are relative to the lowest conformer.
        """
        energies_file = Path(crest_dir) / "crest.energies"
        if not energies_file.exists():
            return None
        try:
            text = energies_file.read_text(encoding="utf-8", errors="ignore")
            rel_energies: List[float] = []
            for line in text.splitlines():
                parts = line.split()
                if not parts:
                    continue
                try:
                    rel_energies.append(float(parts[-1]))
                except ValueError:
                    continue
            return rel_energies if rel_energies else None
        except OSError:
            return None

    @staticmethod
    def validate_energies(
        xyz_energies: List[float],
        rel_energies_kcal: Optional[List[float]],
        tolerance_kcal: float = 1.0,
    ) -> None:
        """Cross-validate XYZ absolute energies against ``crest.energies``.

        Derives relative energies from absolute Hartree values and checks
        consistency with CREST's reported relative energies.  Logs a warning
        on mismatch but does not raise — the XYZ absolute energies are the
        primary source of truth.
        """
        if rel_energies_kcal is None or not xyz_energies:
            return
        if len(rel_energies_kcal) != len(xyz_energies):
            logger.warning(
                "crest.energies count (%d) != candidate count (%d); "
                "skipping cross-validation",
                len(rel_energies_kcal),
                len(xyz_energies),
            )
            return
        min_hartree = min(xyz_energies)
        derived_rel = [(e - min_hartree) * HARTREE_TO_KCAL for e in xyz_energies]
        max_diff = max(abs(d - r) for d, r in zip(derived_rel, rel_energies_kcal))
        if max_diff > tolerance_kcal:
            logger.warning(
                "crest.energies cross-validation mismatch (max diff %.3f kcal/mol "
                "> tolerance %.3f); using XYZ absolute energies as primary source",
                max_diff,
                tolerance_kcal,
            )

    def run_sp(
        self,
        xyz_path: Path,
        nprocs: Optional[int] = None,
        maxcore: Optional[int] = None,
    ) -> Optional[float]:
        """Run an ORCA B97-3c single-point calculation.

        Args:
            xyz_path: Path to the input XYZ file.
            nprocs: Override the number of ORCA cores (for parallel fan-out).
                If ``None``, uses the config default.
            maxcore: Per-rank ORCA memory ceiling resolved for the complete
                concurrent wave. If ``None``, uses the ranking/config default.
        """
        cfg = dict(self.config.get("step1", {}).get("censo_lite", {}).get("ranking", {}) or {})
        resources = dict(self.config.get("resources", {}) or {})
        effective_nprocs = int(nprocs or cfg.get("nproc") or resources.get("nproc", 1))
        logger.debug("B97-3c SP: %s/CPCM(%s) nprocs=%d", cfg.get("method", "B97-3c"), cfg.get("solvent", "acetone"), effective_nprocs)
        interface = ORCAInterface(
            method=str(cfg.get("method", "B97-3c")),
            basis=str(cfg.get("basis", "")),
            aux_basis=str(cfg.get("aux_basis", "")),
            nprocs=effective_nprocs,
            maxcore=maxcore if maxcore is not None else cfg.get("maxcore"),
            solvent=str(cfg.get("solvent", "acetone")),
            route_extras=str(cfg.get("route_extras", "")),
            config=self.config,
        )
        sp_output_dir = self.ranking_dir / Path(xyz_path).stem
        sp_output_dir.mkdir(parents=True, exist_ok=True)
        result = interface.single_point(
            Path(xyz_path).resolve(),
            sp_output_dir,
            timeout=cfg.get("timeout"),
            charge=int(cfg.get("charge", 0)),
            spin=int(cfg.get("multiplicity", 1)),
        )
        return float(result.energy) if result.converged and result.energy is not None else None

    def run_mrrho(
        self, xyz_path: Path, nprocs: Optional[int] = None
    ) -> Optional[XTBThermoResult]:
        """Run xTB SPH+mRRHO and return its validated energy ledger.

        Returns the thermostatistical free-energy correction ``G(T)`` which
        should be added to a DFT electronic energy::

            G_ranked = E_B97-3c + result.g_rrho_correction_hartree

        Returns ``None`` only when mRRHO is disabled. Failed attempts return the
        final ``XTBThermoResult(success=False)`` so diagnostics are not erased.
        """
        cfg = dict(self.config.get("step1", {}).get("censo_lite", {}).get("xtb_thermo", {}) or {})
        if not bool(cfg.get("enabled", True)):
            return None
        resources = dict(self.config.get("resources", {}) or {})
        solvent = str(cfg.get("solvent", self.config.get("step1", {}).get("censo_lite", {}).get("crest", {}).get("solvent", "acetone")))
        primary_gfn = int(cfg.get("gfn_level", 1))
        primary_nproc = int(nprocs or cfg.get("cores_per_job") or cfg.get("nproc") or resources.get("nproc", 1))
        retry_nproc = int(cfg.get("retry_nproc") or 1)
        primary_scc_iterations = cfg.get("scc_max_iterations")
        retry_scc_iterations = int(cfg.get("scc_retry_max_iterations") or 1000)
        fallback_enabled = bool(cfg.get("fallback_enabled", False))
        fallback_gfn = int(cfg.get("fallback_gfn_level", primary_gfn))
        fallback_nproc = int(cfg.get("fallback_nproc") or 1)

        attempts: List[Tuple[int, int, Optional[int]]] = [
            (primary_gfn, primary_nproc, primary_scc_iterations)
        ]
        retry_attempt = (primary_gfn, retry_nproc, retry_scc_iterations)
        if retry_attempt not in attempts:
            attempts.append(retry_attempt)
        fallback_attempt = (fallback_gfn, fallback_nproc, retry_scc_iterations)
        if fallback_enabled and fallback_attempt not in attempts:
            attempts.append(fallback_attempt)

        base_dir = self.mrrho_dir / Path(xyz_path).stem
        last_result: Optional[XTBThermoResult] = None
        for attempt_index, (gfn_level, nproc, max_scc_iterations) in enumerate(attempts):
            attempt_dir = base_dir if attempt_index == 0 else base_dir / (
                f"retry_{attempt_index}_gfn{gfn_level}_nproc{nproc}_scc{max_scc_iterations}"
            )
            interface = XTBInterface(
                gfn_level=gfn_level,
                solvent=solvent,
                nproc=nproc,
                config=self.config,
            )
            if attempt_index:
                self._emit(
                    "batch_job_retry",
                    batch=f"{self.molecule_name}:mrrho",
                    job_id=Path(xyz_path).stem,
                    attempt=attempt_index + 1,
                    gfn_level=gfn_level,
                    nprocs=nproc,
                    max_scc_iterations=max_scc_iterations,
                    output=str(attempt_dir),
                    reason="previous xTB mRRHO attempt failed",
                )
                logger.debug(
                    "Retrying xTB mRRHO for %s with GFN%d, nproc=%d, SCC iterations=%s",
                    Path(xyz_path).stem,
                    gfn_level,
                    nproc,
                    max_scc_iterations,
                )
            result = interface.enso_thermo(
                xyz_file=Path(xyz_path).resolve(),
                output_dir=attempt_dir,
                charge=int(cfg.get("charge", 0)),
                spin=int(cfg.get("multiplicity", 1)),
                temperature_k=float(cfg.get("temperature_k", 298.15)),
                sthr=float(cfg.get("sthr", 50.0)),
                imagthr=cfg.get("imagthr", -100.0),
                solvent=solvent,
                timeout=cfg.get("timeout"),
                max_scc_iterations=max_scc_iterations,
            )
            last_result = result
            if result.success and result.g_rrho_correction_hartree is not None:
                return result

            error_message = result.error or "xTB mRRHO failed without diagnostics"
            if attempt_index < len(attempts) - 1 and self._is_mrrho_retryable_failure(error_message):
                continue
            logger.debug(
                "xTB mRRHO failed for %s without retry: %s",
                Path(xyz_path).stem,
                error_message[:300],
            )
            break
        return last_result

    @staticmethod
    def _is_mrrho_retryable_failure(error_message: str) -> bool:
        """Return whether an xTB crash-like mRRHO failure is retryable."""

        message = (error_message or "").lower()
        if any(token in message for token in ("timed out", "timeout", "binary not found")):
            return False
        return any(
            token in message
            for token in (
                "sigsegv",
                "segmentation fault",
                "forrtl:",
                "rc=174",
                "rc=-11",
                "rc=128",
                ".sccnotconverged",
                "self consistent charge iterator did not converge",
                "longjmp causes uninitialized stack frame",
                "abnormal termination",
                "access violation",
            )
        )
