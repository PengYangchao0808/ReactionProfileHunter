"""Artifact-first validation and no-QC migration for V4 resume planning."""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

from rdkit import Chem

from rph_core.utils.atom_mapping import (
    PebResolutionError,
    load_s1_smiles_to_xyz_map,
    write_smiles_to_xyz_mapping,
)

logger = logging.getLogger(__name__)


S1_CONTRACT_FIELDS = frozenset(
    {"signature_schema", "manifest_schema", "protocol_version", "schema_version"}
)


@dataclass(frozen=True)
class ArtifactAssessment:
    node: str
    status: str
    reason: str
    manifest: Optional[Path] = None
    selected_xyz: Optional[Path] = None
    atom_mapping: Optional[Path] = None
    evidence: Dict[str, Any] = field(default_factory=dict)

    @property
    def reusable(self) -> bool:
        return self.status in {"reusable_exact", "reusable_validated", "reusable_after_migration"}


def canonical_molecule_identity(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f"Invalid SMILES in artifact provenance: {smiles}")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


class ArtifactReconciler:
    """Validate existing artifacts and materialize current contracts without QC."""

    def __init__(self, work_dir: Path, current_s1_schema: str):
        self.work_dir = Path(work_dir).resolve()
        self.current_s1_schema = current_s1_schema

    @staticmethod
    def _read_json(path: Path) -> Dict[str, Any]:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError(f"Manifest must contain an object: {path}")
        return payload

    @staticmethod
    def _selected_xyz(manifest_path: Path, payload: Mapping[str, Any]) -> Path:
        direct = payload.get("selected_xyz")
        if isinstance(direct, str) and direct.strip():
            return manifest_path.parent / direct
        selected = payload.get("selected")
        for candidate in payload.get("candidates", []) or []:
            if isinstance(candidate, dict) and candidate.get("id") == selected:
                raw = Path(str(candidate["xyz"]))
                local = manifest_path.parent / raw
                if local.is_file():
                    return local
                return manifest_path.parent.parent / raw
        raise ValueError(f"S1 manifest has no resolvable selected XYZ: {manifest_path}")

    @staticmethod
    def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(
            json.dumps(dict(payload), indent=2, ensure_ascii=False, default=str),
            encoding="utf-8",
        )
        temporary.replace(path)

    def assess_s1(
        self,
        manifest_path: Path,
        *,
        expected_smiles: str,
        recorded_smiles: Optional[str] = None,
        changed_paths: Sequence[str] = (),
    ) -> ArtifactAssessment:
        manifest_path = Path(manifest_path)
        if not manifest_path.is_file():
            return ArtifactAssessment("s1", "missing", "manifest_missing", manifest=manifest_path)
        try:
            payload = self._read_json(manifest_path)
            selected_xyz = self._selected_xyz(manifest_path, payload)
        except (OSError, ValueError, json.JSONDecodeError, KeyError) as exc:
            return ArtifactAssessment("s1", "corrupt", str(exc), manifest=manifest_path)
        if not selected_xyz.is_file() or selected_xyz.stat().st_size == 0:
            return ArtifactAssessment(
                "s1", "corrupt", "selected_xyz_missing", manifest_path, selected_xyz
            )

        provenance_smiles = str(recorded_smiles or expected_smiles)
        try:
            same_molecule = canonical_molecule_identity(provenance_smiles) == canonical_molecule_identity(expected_smiles)
        except ValueError as exc:
            return ArtifactAssessment("s1", "unsupported_legacy", str(exc), manifest_path, selected_xyz)
        if not same_molecule:
            return ArtifactAssessment(
                "s1", "dependency_changed", "molecular_identity_changed", manifest_path, selected_xyz
            )

        leaf_paths = {str(path).rsplit(".", 1)[-1] for path in changed_paths}
        incompatible = leaf_paths - S1_CONTRACT_FIELDS - {"smiles"}
        if incompatible:
            return ArtifactAssessment(
                "s1",
                "dependency_changed",
                f"scientific_signature_changed:{','.join(sorted(incompatible))}",
                manifest_path,
                selected_xyz,
            )

        mapping_ref = str(payload.get("atom_mapping_ref") or "atom_mapping.json")
        mapping_path = manifest_path.parent / mapping_ref
        if mapping_path.is_file() and payload.get("schema_version") == self.current_s1_schema:
            try:
                load_s1_smiles_to_xyz_map(mapping_path, selected_xyz)
            except PebResolutionError as exc:
                return ArtifactAssessment("s1", "invalid_mapping", str(exc), manifest_path, selected_xyz, mapping_path)
            return ArtifactAssessment(
                "s1",
                "reusable_exact",
                "current_contract_verified",
                manifest_path,
                selected_xyz,
                mapping_path,
            )

        return ArtifactAssessment(
            "s1",
            "reusable_after_migration",
            "qc_complete_contract_upgrade_required",
            manifest_path,
            selected_xyz,
            mapping_path,
            evidence={"provenance_smiles": provenance_smiles},
        )

    def migrate_s1(self, assessment: ArtifactAssessment) -> ArtifactAssessment:
        if assessment.status != "reusable_after_migration":
            return assessment
        assert assessment.manifest is not None
        assert assessment.selected_xyz is not None
        provenance_smiles = str(assessment.evidence["provenance_smiles"])
        payload = self._read_json(assessment.manifest)
        old_bytes = assessment.manifest.read_bytes()
        old_hash = hashlib.sha256(old_bytes).hexdigest()
        migration_dir = assessment.manifest.parent / "migrations"
        migration_dir.mkdir(parents=True, exist_ok=True)
        backup = migration_dir / f"manifest_{payload.get('schema_version', 'legacy')}_{old_hash[:12]}.json"
        if not backup.exists():
            backup.write_bytes(old_bytes)

        mapping_path = write_smiles_to_xyz_mapping(
            provenance_smiles,
            assessment.selected_xyz,
            assessment.manifest.parent / "atom_mapping.json",
        )
        payload.update(
            {
                "schema_version": self.current_s1_schema,
                "selected_xyz": str(assessment.selected_xyz.relative_to(assessment.manifest.parent)),
                "atom_mapping_ref": mapping_path.name,
                "atom_mapping_status": "verified",
                "artifact_reconciliation": {
                    "status": "reusable_after_migration",
                    "qc_recomputed": False,
                    "source_manifest": str(backup),
                    "source_manifest_sha256": old_hash,
                    "geometry_reference_smiles": provenance_smiles,
                },
            }
        )
        self._atomic_write_json(assessment.manifest, payload)
        load_s1_smiles_to_xyz_map(mapping_path, assessment.selected_xyz)
        logger.info("[V4] Migrated S1 artifact without QC: %s", assessment.manifest)
        return ArtifactAssessment(
            "s1",
            "reusable_validated",
            "manifest_and_mapping_migrated_without_qc",
            assessment.manifest,
            assessment.selected_xyz,
            mapping_path,
            evidence={"source_manifest_sha256": old_hash, "qc_recomputed": False},
        )
