from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

from rph_core.utils.bond_pairs import canonicalize_bond_pairs


@dataclass(frozen=True)
class Provenance:
    run_id: str | None
    schema_version: str | None
    protocol_version: str | None
    parent_stage_manifest_hashes: dict[str, str]
    atom_mapping_sha256: str | None
    forming_bonds_hash: str | None
    variant_manifest_hash: str | None
    extra: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def sha256_of_file(path: Path | None) -> str:
    """SHA-256 hexdigest of a file's bytes, or empty string if path is None/missing."""

    if path is None:
        return ""
    candidate = Path(path)
    if not candidate.is_file():
        return ""
    try:
        return hashlib.sha256(candidate.read_bytes()).hexdigest()
    except OSError:
        return ""


def sha256_of_json(payload: Any) -> str:
    """SHA-256 hexdigest of json.dumps(payload, sort_keys=True, default=str)."""

    serialized = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def canonicalize_forming_bonds_for_hash(bonds: Sequence[Sequence[int]] | None) -> str:
    """Return canonicalized JSON string for unordered forming-bond pairs."""

    canonical_pairs = canonicalize_bond_pairs(bonds or [])
    return json.dumps([list(pair) for pair in canonical_pairs])


def _normalize_index_value(value: Any) -> int | str:
    if isinstance(value, bool):
        return str(value)
    try:
        return int(value)
    except (TypeError, ValueError):
        return str(value)


def _sort_key(value: int | str) -> tuple[int, int | str]:
    return (0, value) if isinstance(value, int) else (1, value)


def _sorted_mapping_items(mapping: Mapping[Any, Any]) -> list[list[int | str]]:
    items = [
        [_normalize_index_value(key), _normalize_index_value(value)]
        for key, value in mapping.items()
    ]
    items.sort(key=lambda row: (_sort_key(row[0]), _sort_key(row[1])))
    return items


def _coerce_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _pair_mapping(left_pairs: Any, right_pairs: Any) -> dict[int, int]:
    if not isinstance(left_pairs, list) or not isinstance(right_pairs, list):
        return {}
    if len(left_pairs) != len(right_pairs):
        return {}
    mapping: dict[int, int] = {}
    for left_pair, right_pair in zip(left_pairs, right_pairs):
        if not isinstance(left_pair, (list, tuple)) or not isinstance(right_pair, (list, tuple)):
            continue
        if len(left_pair) != 2 or len(right_pair) != 2:
            continue
        for left_value, right_value in zip(left_pair, right_pair):
            left_index = _coerce_int(left_value)
            right_index = _coerce_int(right_value)
            if left_index is None or right_index is None:
                continue
            mapping[left_index] = right_index
    return mapping


def _extract_smiles_to_xyz(payload: Mapping[str, Any]) -> dict[int, int]:
    explicit = payload.get("smiles_to_xyz")
    if isinstance(explicit, Mapping):
        mapping = {
            key_int: value_int
            for raw_key, raw_value in explicit.items()
            if (key_int := _coerce_int(raw_key)) is not None
            and (value_int := _coerce_int(raw_value)) is not None
        }
        if mapping:
            return mapping

    atoms = payload.get("atoms")
    if isinstance(atoms, list):
        mapping: dict[int, int] = {}
        for row in atoms:
            if not isinstance(row, Mapping):
                continue
            smiles_idx = _coerce_int(row.get("smiles_idx"))
            xyz_idx = _coerce_int(row.get("xyz_idx"))
            if smiles_idx is None or xyz_idx is None:
                continue
            mapping[smiles_idx] = xyz_idx
        if mapping:
            return mapping

    forming_bonds = payload.get("forming_bonds")
    if isinstance(forming_bonds, Mapping):
        mapping = _pair_mapping(
            forming_bonds.get("product_smiles_idx"),
            forming_bonds.get("product_xyz_0based"),
        )
        if mapping:
            return mapping

    map_to_product_smiles = payload.get("map_to_product_smiles")
    if (
        isinstance(map_to_product_smiles, Mapping)
        and payload.get("product_smiles_idx_space") == "geometry_product_smiles_idx"
    ):
        identity_mapping = {
            smiles_idx: smiles_idx
            for raw_value in map_to_product_smiles.values()
            if (smiles_idx := _coerce_int(raw_value)) is not None
        }
        if identity_mapping:
            return identity_mapping

    return {}


def _extract_map_to_xyz(
    payload: Mapping[str, Any],
    smiles_to_xyz: Mapping[int, int],
) -> dict[int, int]:
    explicit = payload.get("map_to_xyz")
    if isinstance(explicit, Mapping):
        mapping = {
            key_int: value_int
            for raw_key, raw_value in explicit.items()
            if (key_int := _coerce_int(raw_key)) is not None
            and (value_int := _coerce_int(raw_value)) is not None
        }
        if mapping:
            return mapping

    atoms = payload.get("atoms")
    if isinstance(atoms, list):
        atom_mapping: dict[int, int] = {}
        for row in atoms:
            if not isinstance(row, Mapping):
                continue
            map_number = _coerce_int(row.get("map_number"))
            xyz_idx = _coerce_int(row.get("xyz_idx"))
            if map_number is None or xyz_idx is None:
                continue
            atom_mapping[map_number] = xyz_idx
        if atom_mapping:
            return atom_mapping

    forming_bonds = payload.get("forming_bonds")
    if isinstance(forming_bonds, Mapping):
        mapping = _pair_mapping(
            forming_bonds.get("map_space"),
            forming_bonds.get("product_xyz_0based"),
        )
        if mapping:
            return mapping

    map_to_product_smiles = payload.get("map_to_product_smiles")
    if isinstance(map_to_product_smiles, Mapping):
        mapping: dict[int, int] = {}
        for raw_map_idx, raw_smiles_idx in map_to_product_smiles.items():
            map_idx = _coerce_int(raw_map_idx)
            smiles_idx = _coerce_int(raw_smiles_idx)
            if map_idx is None or smiles_idx is None:
                continue
            mapping[map_idx] = smiles_to_xyz.get(smiles_idx, smiles_idx)
        if mapping:
            return mapping

    return {}


def hash_atom_mapping(payload: Mapping[str, Any]) -> str:
    """SHA-256 over a canonical, index-space-stable subset of an atom mapping."""

    smiles_to_xyz = _extract_smiles_to_xyz(payload)
    map_to_xyz = _extract_map_to_xyz(payload, smiles_to_xyz)
    canonical_payload: dict[str, Any] = {
        "mapped_smiles": payload.get("mapped_smiles")
        or payload.get("mapped_product_smiles")
        or payload.get("reference_smiles"),
        "smiles_to_xyz": _sorted_mapping_items(smiles_to_xyz),
    }
    if map_to_xyz:
        canonical_payload["map_to_xyz"] = _sorted_mapping_items(map_to_xyz)
    return sha256_of_json(canonical_payload)


def build_provenance(
    *,
    run_id: str | None,
    schema_version: str | None,
    protocol_version: str | None,
    parent_manifest_paths: Mapping[str, Path | None],
    atom_mapping_payload: Mapping[str, Any] | None = None,
    forming_bonds: Sequence[Sequence[int]] | None = None,
    variant_manifest_path: Path | None = None,
    extra: Mapping[str, str] | None = None,
) -> Provenance:
    """Build a Provenance record from on-disk parent manifests and canonical inputs."""

    parent_hashes = {
        str(stage): sha256_of_file(path)
        for stage, path in parent_manifest_paths.items()
    }
    atom_mapping_sha256 = hash_atom_mapping(atom_mapping_payload) if atom_mapping_payload is not None else None
    forming_bonds_hash = None
    if forming_bonds is not None:
        forming_bonds_hash = hashlib.sha256(
            canonicalize_forming_bonds_for_hash(forming_bonds).encode("utf-8")
        ).hexdigest()
    normalized_extra = {
        str(key): "" if value is None else str(value)
        for key, value in (extra or {}).items()
    }
    return Provenance(
        run_id=run_id,
        schema_version=schema_version,
        protocol_version=protocol_version,
        parent_stage_manifest_hashes=parent_hashes,
        atom_mapping_sha256=atom_mapping_sha256,
        forming_bonds_hash=forming_bonds_hash,
        variant_manifest_hash=(
            None if variant_manifest_path is None else sha256_of_file(variant_manifest_path)
        ),
        extra=normalized_extra,
    )


def verify_provenance_chain(
    manifest_path: Path,
    *,
    expected_parent_paths: Mapping[str, Path],
) -> tuple[bool, list[str]]:
    """Re-hash expected parent manifests and compare them to recorded provenance."""

    try:
        payload = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return False, [f"failed to read manifest {manifest_path}: {exc}"]

    provenance = payload.get("provenance")
    if not isinstance(provenance, Mapping):
        return False, [f"{manifest_path}: missing provenance payload"]

    recorded_hashes = provenance.get("parent_stage_manifest_hashes")
    if not isinstance(recorded_hashes, Mapping):
        return False, [f"{manifest_path}: missing provenance.parent_stage_manifest_hashes"]

    mismatches: list[str] = []
    for stage, expected_path in expected_parent_paths.items():
        actual_hash = sha256_of_file(expected_path)
        recorded_hash = recorded_hashes.get(stage)
        if recorded_hash != actual_hash:
            mismatches.append(
                f"{stage}: recorded={recorded_hash!r} actual={actual_hash!r} path={expected_path}"
            )
    return not mismatches, mismatches
