"""
Reference States Runner

Manages reference state calculations for precursors and leaving small molecules.
Creates directory structure and index files (without QC in V5.1).
"""

from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict
import json

from rph_core.utils.tsv_dataset import ReactionRecord, collect_leaving_small_molecule_keys
from rph_core.utils.small_molecule_catalog import SmallMoleculeCatalog, UnknownSmallMoleculeError


@dataclass
class ReferenceStateEntry:
    """Reference state for a single species (precursor or small molecule)."""
    smiles: str
    charge: int = 0
    multiplicity: int = 1
    global_min_xyz: str = ""
    sp_energy_hartree: Optional[float] = None
    g_used_hartree: Optional[float] = None

    thermo: Optional[Dict[str, Any]] = None
    artifacts: Optional[Dict[str, str]] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ReactionReferenceState:
    """Reference state for a single reaction record."""
    rx_id: str
    precursor: ReferenceStateEntry
    raw_meta: Dict[str, str]
    leaving_small_molecule_key: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        if self.precursor:
            d['precursor'] = self.precursor.to_dict()
        return d


class ReferenceStatesRunner:
    """
    Runner for reference state calculations.

    V5.1: Creates directory structure and writes index files
    without running actual QC calculations.
    """

    def __init__(self, config: Dict, work_dir: Path):
        """
        Initialize reference states runner.

        Args:
            config: Configuration dictionary
            work_dir: Base working directory
        """
        self.config = config
        self.work_dir = Path(work_dir).resolve()

        ref_config = config.get('reference_states', {})
        self.enabled = ref_config.get('enabled', True)
        self.base_dirname = ref_config.get('base_dirname', 'reference_states')

        # Initialize small molecule catalog
        self.small_mol_catalog = SmallMoleculeCatalog(config)

        # Results storage
        self.reaction_states: Dict[str, ReactionReferenceState] = {}
        self.small_molecule_states: Dict[str, ReferenceStateEntry] = {}
        self.errors: Dict[str, List[str]] = {}

    def run(
        self,
        records: List[ReactionRecord],
        require_thermo: Optional[bool] = None
    ) -> Dict[str, Any]:
        """
        Run reference states generation for given reaction records.

        Args:
            records: List of ReactionRecord objects
            require_thermo: Whether thermo is required (from config if None)

        Returns:
            Summary dict with stats and paths
        """
        if require_thermo is None:
            ref_config = self.config.get('reference_states', {})
            require_thermo_val = ref_config.get('require_thermo', True)
        else:
            require_thermo_val = require_thermo

        # Create base directory
        ref_base = self.work_dir / self.base_dirname
        ref_base.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        rx_dir = ref_base / 'reactions'
        rx_dir.mkdir(exist_ok=True)
        smol_dir = ref_base / 'small_molecules'
        smol_dir.mkdir(exist_ok=True)

        # Collect leaving small molecule keys
        leaving_keys = collect_leaving_small_molecule_keys(records)

        # Validate small molecule keys
        unknown_keys = self.small_mol_catalog.validate_keys(list(leaving_keys))
        if unknown_keys:
            policy = self._get_unknown_small_mol_policy()
            if policy == 'error':
                raise UnknownSmallMoleculeError(
                    f"Unknown small molecules: {unknown_keys}"
                )
            elif policy == 'warn':
                self.errors['unknown_small_molecules'] = unknown_keys
            else:  # 'skip'
                pass

        # Process each reaction record
        for record in records:
            self._process_reaction(record, rx_dir)

        # Process each small molecule
        for key in leaving_keys:
            if key not in unknown_keys:
                self._process_small_molecule(key, smol_dir)

        # Write index.json
        index_path = ref_base / 'index.json'
        self._write_index(index_path, require_thermo_val)

        return {
            'success': True,
            'total_reactions': len(self.reaction_states),
            'total_small_molecules': len(self.small_molecule_states),
            'errors': self.errors,
            'index_path': index_path
        }

    def _process_reaction(self, record: ReactionRecord, rx_dir: Path):
        """Process a single reaction record (V5.1: write structure only)."""
        rx_subdir = rx_dir / f'rx_{record.rx_id}'
        rx_subdir.mkdir(parents=True, exist_ok=True)

        precursor_dir = rx_subdir / 'precursor'
        precursor_dir.mkdir(exist_ok=True)

        # Write placeholder energy.json (V5.1: no QC yet)
        placeholder_entry = ReferenceStateEntry(
            smiles=record.precursor_smiles,
            global_min_xyz='precursor_global_min.xyz'
        )
        energy_path = precursor_dir / 'energy.json'
        energy_path.write_text(json.dumps(placeholder_entry.to_dict(), indent=2))

        # Write meta.json
        meta_path = rx_subdir / 'meta.json'
        meta_path.write_text(json.dumps(record.raw, indent=2))

        # Store in memory
        self.reaction_states[record.rx_id] = ReactionReferenceState(
            rx_id=record.rx_id,
            precursor=placeholder_entry,
            leaving_small_molecule_key=record.get_leaving_small_molecule_key(),
            raw_meta=record.raw
        )

    def _process_small_molecule(self, key: str, smol_dir: Path):
        """Process a single small molecule (V5.1: write structure only)."""
        mol = self.small_mol_catalog.require(key)

        mol_dir = smol_dir / key
        mol_dir.mkdir(parents=True, exist_ok=True)

        # Write placeholder energy.json
        placeholder_entry = ReferenceStateEntry(
            smiles=mol.smiles,
            charge=mol.charge,
            multiplicity=mol.multiplicity,
            global_min_xyz=f'{key}_global_min.xyz'
        )
        energy_path = mol_dir / 'energy.json'
        energy_path.write_text(json.dumps(placeholder_entry.to_dict(), indent=2))

        # Store in memory
        self.small_molecule_states[key] = placeholder_entry

    def _write_index(self, index_path: Path, require_thermo: bool):
        """Write index.json with all reference state data."""
        index_data = {
            'version': '1',
            'csv_schema': {
                'id_col': 'rx_id',
                'precursor_smiles_col': 'precursor_smiles',
                'ylide_leaving_group_col': 'ylide_leaving_group',
                'leaving_group_col_fallback': 'leaving_group'
            },
            'config_snapshot': {
                'theory_opt': self.config.get('theory', {}).get('optimization', {}),
                'theory_sp': self.config.get('theory', {}).get('single_point', {})
            },
            'reactions': {
                rx_id: state.to_dict()
                for rx_id, state in self.reaction_states.items()
            },
            'small_molecules': {
                key: entry.to_dict()
                for key, entry in self.small_molecule_states.items()
            },
            'errors': self.errors
        }

        if require_thermo:
            # Validate that thermo fields are present (V5.2: machine-detectable validation)
            missing_thermo_rx = [
                rx_id for rx_id, state in self.reaction_states.items()
                if state.precursor.thermo is None
            ]
            missing_thermo_sm = [
                key for key, entry in self.small_molecule_states.items()
                if entry.thermo is None
            ]
            if missing_thermo_rx or missing_thermo_sm:
                self.errors['missing_thermo_reactions'] = missing_thermo_rx
                self.errors['missing_thermo_small_molecules'] = missing_thermo_sm

                # Mark validation failed for machine detection
                index_data['thermo_validation_passed'] = False
                index_data['thermo_validation_errors'] = {
                    'missing_reactions': len(missing_thermo_rx),
                    'missing_small_molecules': len(missing_thermo_sm)
                }
            else:
                index_data['thermo_validation_passed'] = True
        else:
            index_data['thermo_validation_passed'] = None

        index_path.write_text(json.dumps(index_data, indent=2))

        index_path.write_text(json.dumps(index_data, indent=2))

    def _get_unknown_small_mol_policy(self) -> str:
        """Get policy for unknown small molecules from config."""
        ref_config = self.config.get('reference_states', {})
        return ref_config.get('unknown_small_molecule_policy', 'error')
