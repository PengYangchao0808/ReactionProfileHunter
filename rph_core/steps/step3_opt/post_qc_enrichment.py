"""
Post-QC Enrichment
===================

DIAS/ASM enrichment with tether-cut fragmenter, executed after TS success.

Author: QCcalc Team
Date: 2026-01-27
"""

import logging
import hashlib
import json
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, TYPE_CHECKING
from dataclasses import asdict

import numpy as np

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import write_xyz, read_xyz
from rph_core.utils.orca_interface import ORCAInterface
from rph_core.steps.step3_opt.intramolecular_fragmenter import IntramolecularFragmenter, FragmenterResult
from rph_core.utils.fragment_manipulation import get_fragment_charges, get_fragment_multiplicities

logger = logging.getLogger(__name__)

# Import SPMatrixReport from ts_optimizer to avoid circular import
try:
    from rph_core.steps.step3_opt.ts_optimizer import SPMatrixReport
except ImportError:
    SPMatrixReport = None
    logger.warning("SPMatrixReport not available - enrichment will be limited")


class PostQCEnrichment(LoggerMixin):
    """
    Post-QC enrichment handler for DIAS/ASM features.

    Executes after TS optimization success:
    1. Fragments system using tether-cut strategy
    2. Runs SP calculations on capped fragments (R and TS geometries)
    3. Writes enrichment contract JSON
    """

    def __init__(self, config: dict = None):
        """
        Initialize enrichment handler.

        Args:
            config: Configuration dict with enrichment settings
        """
        self.config = config or {}
        self.fragmenter = IntramolecularFragmenter()

    def run(
        self,
        s3_dir: Path,
        reactant_complex_xyz: Path,
        ts_final_xyz: Path,
        sp_report: "SPMatrixReport",
        forming_bonds: Tuple[Tuple[int, int], Tuple[int, int]]
    ) -> None:
        """
        Run enrichment workflow.

        Args:
            s3_dir: Step3 output directory (e.g., work_dir/S3_TransitionAnalysis/)
            reactant_complex_xyz: Path to reactant complex XYZ
            ts_final_xyz: Path to TS final XYZ
            sp_report: SPMatrixReport with electronic energies
            forming_bonds: ((u1, v1), (u2, v2)) forming bond atom pairs

        Writes:
            <s3_dir>/S3_PostQCEnrichment/enrichment.json
            <s3_dir>/S3_PostQCEnrichment/enrichment_status.json
        """
        if not self.config.get('enabled', False):
            self.logger.info("Enrichment disabled by config, skipping")
            self._write_status(s3_dir, "disabled", "Enrichment disabled by config")
            return

        enrichment_dir = s3_dir / self.config.get('write_dirname', 'S3_PostQCEnrichment')
        enrichment_dir.mkdir(parents=True, exist_ok=True)

        try:
            self.logger.info(f"Starting post-QC enrichment in {enrichment_dir}")

            reactant_coords, reactant_symbols = read_xyz(reactant_complex_xyz)
            ts_coords, ts_symbols = read_xyz(ts_final_xyz)

            total_charge = self._extract_system_charge(sp_report)
            total_multiplicity = getattr(sp_report, "total_multiplicity", 1) if sp_report else 1

            fragment_result = self.fragmenter.fragment(
                reactant_coords=reactant_coords,
                reactant_symbols=reactant_symbols,
                ts_coords=ts_coords,
                ts_symbols=ts_symbols,
                forming_bonds=forming_bonds,
                config=self.config.get('fragmenter', {})
            )

            if fragment_result.status != "ok":
                self.logger.warning(f"Fragmenter failed: {fragment_result.reason}")
                self._write_status(
                    s3_dir,
                    "fragmenter_failed",
                    f"Fragmentation failed: {fragment_result.reason}"
                )
                return

            self._write_fragment_geometries(
                enrichment_dir,
                fragment_result
            )

            enrichment_data = self._build_enrichment_data(
                fragment_result,
                total_charge,
                sp_report,
                forming_bonds
            )

            self._write_enrichment_json(enrichment_dir, enrichment_data)

            self._write_status(s3_dir, "ok", "Enrichment completed successfully")

            self.logger.info("Post-QC enrichment completed successfully")

        except Exception as e:
            self.logger.error(f"Post-QC enrichment failed: {e}", exc_info=True)
            self._write_status(
                s3_dir,
                "failed",
                f"Exception: {str(e)}"
            )

    def _extract_system_charge(self, sp_report: "SPMatrixReport") -> int:
        """
        Extract total system charge from SPMatrixReport.

        Args:
            sp_report: SPMatrixReport instance

        Returns:
            Total charge (defaults to 0)
        """
        return getattr(sp_report, 'total_charge', 0)

    def _extract_system_multiplicity(self, sp_report: "SPMatrixReport") -> int:
        """
        Extract total system multiplicity from SPMatrixReport.

        Args:
            sp_report: SPMatrixReport instance

        Returns:
            Total multiplicity (defaults to 1)
        """
        return getattr(sp_report, 'total_multiplicity', 1)

    def _build_orca_template_id(self) -> str:
        """
        Generate canonical ORCA template identifier.

        Format: "orca:{method}/{basis}:{aux_basis}:{solvent}:{nprocs}:{maxcore}"

        Returns:
            Canonical template identifier string
        """
        method = self.config.get('orca_method', 'M062X')
        basis = self.config.get('orca_basis', 'def2-TZVPP')
        aux_basis = self.config.get('orca_aux_basis', 'def2/J')
        solvent = self.config.get('orca_solvent', 'acetone')
        nprocs = self.config.get('orca_nprocs', 16)
        maxcore = self.config.get('orca_maxcore', 8000)

        return f"orca:{method}/{basis}:{aux_basis}:{solvent}:{nprocs}:{maxcore}"

    def _write_fragment_geometries(
        self,
        enrichment_dir: Path,
        fragment_result: FragmenterResult
    ) -> None:
        """
        Write capped fragment XYZ files for inspection/debugging.

        Args:
            enrichment_dir: Enrichment output directory
            fragment_result: FragmenterResult with capped geometries
        """
        for name, coords, symbols in [
            ('fragA_R', fragment_result.fragA_coords_R, fragment_result.fragA_symbols_R),
            ('fragB_R', fragment_result.fragB_coords_R, fragment_result.fragB_symbols_R),
            ('fragA_TS', fragment_result.fragA_coords_TS, fragment_result.fragA_symbols_TS),
            ('fragB_TS', fragment_result.fragB_coords_TS, fragment_result.fragB_symbols_TS),
        ]:
            if len(coords) > 0 and len(symbols) > 0:
                xyz_path = enrichment_dir / f"{name}.xyz"
                write_xyz(xyz_path, coords, symbols, energy=0.0)

    def _build_enrichment_data(
        self,
        fragment_result: FragmenterResult,
        total_charge: int,
        sp_report: SPMatrixReport,
        forming_bonds: Tuple[Tuple[int, int], Tuple[int, int]] = None,
    ) -> Dict[str, Any]:
        """
        Build enrichment contract data structure.

        Args:
            fragment_result: FragmenterResult
            total_charge: Total system charge
            sp_report: SPMatrixReport

        Returns:
            Dictionary with enrichment contract data
        """
        n_fragA = len(fragment_result.fragA_indices)
        n_fragB = len(fragment_result.fragB_indices)

        charge_fragA, charge_fragB = get_fragment_charges(
            total_charge, n_fragA, n_fragB, dipole_in_fragA=True
        )

        total_multiplicity = getattr(sp_report, 'total_multiplicity', 1)
        mult_fragA, mult_fragB = get_fragment_multiplicities(
            total_multiplicity, n_fragA, n_fragB, dipole_in_fragA=True
        )

        fragmenter_config = self.config.get('fragmenter', {})
        orca_template_id = self._build_orca_template_id()

        return {
            'fragmenter': {
                'type': fragmenter_config.get('type', 'intramolecular_tether_cut'),
                'status': fragment_result.status,
                'reason': fragment_result.reason if fragment_result.status != 'ok' else None,
                'cut_bond': fragment_result.cut_bond_indices,
                'dipole_core': fragment_result.dipole_core_indices,
                'dipole_ends': fragment_result.dipole_end_indices,
                'alkene_ends': fragment_result.alkene_end_indices,
                'fragment_sizes': {
                    'fragA': n_fragA,
                    'fragB': n_fragB
                }
            },
            'charges': {
                'total': total_charge,
                'fragA': charge_fragA,
                'fragB': charge_fragB
            },
            'multiplicities': {
                'total': total_multiplicity,
                'fragA': mult_fragA,
                'fragB': mult_fragB
            },
            'fragments': {
                'fragA_R': {
                    'atoms': fragment_result.fragA_indices,
                    'charge': charge_fragA,
                    'multiplicity': mult_fragA,
                    'file': 'fragA_R.xyz'
                },
                'fragB_R': {
                    'atoms': fragment_result.fragB_indices,
                    'charge': charge_fragB,
                    'multiplicity': mult_fragB,
                    'file': 'fragB_R.xyz'
                },
                'fragA_TS': {
                    'atoms': fragment_result.fragA_indices,
                    'charge': charge_fragA,
                    'multiplicity': mult_fragA,
                    'file': 'fragA_TS.xyz'
                },
                'fragB_TS': {
                    'atoms': fragment_result.fragB_indices,
                    'charge': charge_fragB,
                    'multiplicity': mult_fragB,
                    'file': 'fragB_TS.xyz'
                }
            },
            'sp_results': {
                'fragA_relaxed': sp_report.e_frag_a_relaxed,
                'fragB_relaxed': sp_report.e_frag_b_relaxed,
                'fragA_ts': sp_report.e_frag_a_ts,
                'fragB_ts': sp_report.e_frag_b_ts
            },
            'orca_config': {
                'template_id': orca_template_id,
                'template_hash': None,
                'method': self.config.get('orca_method', 'M062X'),
                'basis': self.config.get('orca_basis', 'def2-TZVPP'),
                'aux_basis': self.config.get('orca_aux_basis', 'def2/J'),
                'solvent': self.config.get('orca_solvent', 'acetone'),
                'nprocs': self.config.get('orca_nprocs', 16),
                'maxcore': self.config.get('orca_maxcore', 8000)
            },
            'inputs': {
                'forming_bonds': forming_bonds_to_list(forming_bonds)
            },
            'metadata': {
                'version': '6.2',
                'timestamp': None,
                'sp_report_hash': self._compute_sp_report_hash(sp_report)
            }
        }

    def _compute_sp_report_hash(self, sp_report: SPMatrixReport) -> str:
        """
        Compute hash of SPMatrixReport for provenance.

        Args:
            sp_report: SPMatrixReport instance

        Returns:
            SHA256 hash hex string
        """
        report_dict = asdict(sp_report)
        report_str = json.dumps(report_dict, sort_keys=True)
        return hashlib.sha256(report_str.encode()).hexdigest()

    def _write_enrichment_json(
        self,
        enrichment_dir: Path,
        data: Dict[str, Any]
    ) -> None:
        """
        Write enrichment contract JSON file.

        Args:
            enrichment_dir: Enrichment output directory
            data: Enrichment contract data
        """
        enrichment_path = enrichment_dir / 'enrichment.json'
        with open(enrichment_path, 'w') as f:
            json.dump(data, f, indent=2)

        self.logger.info(f"Wrote enrichment contract to {enrichment_path}")

    def _write_status(
        self,
        s3_dir: Path,
        status: str,
        message: str
    ) -> None:
        """
        Write enrichment status JSON file.

        Args:
            s3_dir: Step3 output directory
            status: Status string ("ok", "failed", "disabled", etc.)
            message: Status message
        """
        write_dir = s3_dir / self.config.get('write_dirname', 'S3_PostQCEnrichment')
        write_dir.mkdir(parents=True, exist_ok=True)

        status_path = write_dir / 'enrichment_status.json'
        status_data = {
            'status': status,
            'message': message
        }

        with open(status_path, 'w') as f:
            json.dump(status_data, f, indent=2)

        self.logger.info(f"Wrote enrichment status to {status_path}: {status}")


def forming_bonds_to_list(forming_bonds: Tuple[Tuple[int, int], Tuple[int, int]]) -> list:
    """Convert forming_bonds tuple to list format."""
    return [list(bond) for bond in forming_bonds]
