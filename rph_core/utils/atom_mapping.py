"""
PEB atom-mapping resolver.

Composes S0 atom-map sidecars with S1 SMILES→XYZ sidecars to resolve
forming bonds for Product Edge Breaking (PEB) without post-hoc remapping.
"""

import hashlib
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem

from rph_core.utils.file_io import read_xyz
from rph_core.utils.path_compat import normalize_path

logger = logging.getLogger(__name__)

EXPECTED_PRODUCT_SMILES_IDX_SPACE = "geometry_product_smiles_idx"

REJECTED_S1_MAPPING_SOURCES = {
    "assumed_heavy_order",
    "substructure_match",
    "kabsch_hungarian",
    "legacy",
}

ADDITIVE_OR_H_ELEMENTS = {"H", "Mg", "Cl", "Al", "Zn", "Li", "Na", "K", "Ca", "Br", "I"}
ORGANIC_HEAVY_ELEMENTS = {"C", "N", "O", "S", "P", "F", "B"}

__all__ = [
    "PebMappingResult",
    "PebResolutionError",
    "load_s0_atom_map_smiles",
    "load_s1_smiles_to_xyz_map",
    "build_product_smiles_to_xyz",
    "resolve_peb_forming_bonds",
    "write_peb_mapping",
    "write_smiles_to_xyz_mapping",
    "bind_atom_mapping_to_xyz",
    "verify_atom_mapping_table",
]


@dataclass(frozen=True)
class PebMappingResult:
    """Resolved PEB forming-bond mapping across map/SMILES/XYZ index spaces."""

    forming_bonds_map_space: Tuple[Tuple[int, int], ...]
    forming_bonds_product_smiles: Tuple[Tuple[int, int], ...]
    forming_bonds_product_xyz_0based: Tuple[Tuple[int, int], ...]
    forming_bonds_product_xyz_1based: Tuple[Tuple[int, int], ...]
    confidence: str
    diagnostics: Dict[str, object] = field()


class PebResolutionError(Exception):
    """Fail-fast error for PEB sidecar loading, resolution, or validation."""


def _load_json_object(path: Path, *, missing_message: str, label: str) -> Dict[str, object]:
    normalized_path = normalize_path(path)
    if not normalized_path.exists():
        raise PebResolutionError(missing_message)

    try:
        with normalized_path.open("r", encoding="utf-8") as handle:
            payload: Any = json.load(handle)
    except json.JSONDecodeError as exc:
        raise PebResolutionError(f"{label} is not valid JSON: {normalized_path}: {exc}") from exc
    except OSError as exc:
        raise PebResolutionError(f"Failed to read {label}: {normalized_path}: {exc}") from exc

    if not isinstance(payload, dict):
        raise PebResolutionError(f"{label} must contain a JSON object: {normalized_path}")

    return payload


def _compute_file_hash(file_path: Path) -> Optional[str]:
    normalized_path = normalize_path(file_path)
    if not normalized_path.exists():
        return None

    try:
        with normalized_path.open("rb") as handle:
            return hashlib.sha256(handle.read()).hexdigest()[:16]
    except OSError as exc:
        logger.debug(f"Failed to compute hash for {normalized_path}: {exc}")
        return None


def _normalize_element_symbol(value: object) -> str:
    element = str(value).strip()
    if not element:
        return ""
    if len(element) == 1:
        return element.upper()
    return element[0].upper() + element[1:].lower()


def _xyz_elements(path: Path) -> List[str]:
    try:
        _coordinates, elements = read_xyz(normalize_path(path))
    except (FileNotFoundError, OSError, ValueError, IndexError) as exc:
        raise PebResolutionError(f"Failed to read XYZ atom order from {path}: {exc}") from exc
    return [_normalize_element_symbol(element) for element in elements]


def verify_atom_mapping_table(payload: Dict[str, object], xyz_path: Path) -> None:
    """Fail fast unless a mapping table exactly matches an XYZ atom sequence."""

    atoms = payload.get("atoms")
    if not isinstance(atoms, list) or not atoms:
        raise PebResolutionError("Atom mapping sidecar has no atoms table")
    xyz_elements = _xyz_elements(xyz_path)
    if len(atoms) != len(xyz_elements):
        raise PebResolutionError(
            f"Atom mapping/XYZ atom-count mismatch: mapping={len(atoms)} xyz={len(xyz_elements)}"
        )
    seen: set[int] = set()
    for row_index, row in enumerate(atoms):
        if not isinstance(row, dict):
            raise PebResolutionError(f"atoms[{row_index}] must be an object")
        try:
            xyz_idx = int(row.get("xyz_idx"))
        except (TypeError, ValueError) as exc:
            raise PebResolutionError(f"atoms[{row_index}] has invalid xyz_idx") from exc
        if xyz_idx in seen or xyz_idx != row_index:
            raise PebResolutionError("Atom mapping XYZ indices must be unique and contiguous")
        seen.add(xyz_idx)
        expected = _normalize_element_symbol(row.get("element"))
        if expected != xyz_elements[xyz_idx]:
            raise PebResolutionError(
                f"Atom order drift at XYZ idx {xyz_idx}: mapping={expected} xyz={xyz_elements[xyz_idx]}"
            )


def write_smiles_to_xyz_mapping(smiles: str, xyz_path: Path, output_path: Path) -> Path:
    """Write the generation-time SMILES/QC atom table used by all later stages."""

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise PebResolutionError(f"Cannot build atom mapping for invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    hydrogen_ordinals: Dict[int, int] = {}
    atoms: List[Dict[str, object]] = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        map_number = atom.GetAtomMapNum() or None
        if atom.GetAtomicNum() == 1:
            heavy_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() != 1]
            parent = heavy_neighbors[0].GetIdx() if heavy_neighbors else -1
            ordinal = hydrogen_ordinals.get(parent, 0) + 1
            hydrogen_ordinals[parent] = ordinal
            stable_id = f"H:parent:{parent}:{ordinal}"
            atom_type = "hydrogen"
        else:
            stable_id = f"map:{map_number}" if map_number else f"product_local:{idx}"
            atom_type = "organic_heavy"
        atoms.append(
            {
                "stable_id": stable_id,
                "map_number": map_number,
                "smiles_idx": idx,
                "xyz_idx": idx,
                "element": atom.GetSymbol(),
                "type": atom_type,
            }
        )
    payload: Dict[str, object] = {
        "schema_version": "s1_atom_mapping_v1",
        "mapping_source": "rdkit_generation_sidecar",
        "mapping_status": "verified",
        "confidence": "high",
        "reference_smiles": smiles,
        "product_smiles_idx_space": EXPECTED_PRODUCT_SMILES_IDX_SPACE,
        "xyz_ref": str(normalize_path(xyz_path)),
        "xyz_sha256": _compute_file_hash(xyz_path),
        "atoms": atoms,
    }
    verify_atom_mapping_table(payload, xyz_path)
    normalized_output = normalize_path(output_path)
    normalized_output.parent.mkdir(parents=True, exist_ok=True)
    normalized_output.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return normalized_output


def bind_atom_mapping_to_xyz(source_path: Path, xyz_path: Path, output_path: Path) -> Path:
    """Verify unchanged atom order and bind an existing table to a new geometry."""

    payload = _load_json_object(
        normalize_path(source_path),
        missing_message=f"Atom mapping sidecar not found: {source_path}",
        label="atom mapping sidecar",
    )
    verify_atom_mapping_table(payload, xyz_path)
    payload["xyz_ref"] = str(normalize_path(xyz_path))
    payload["xyz_sha256"] = _compute_file_hash(xyz_path)
    payload["mapping_status"] = "verified"
    normalized_output = normalize_path(output_path)
    normalized_output.parent.mkdir(parents=True, exist_ok=True)
    normalized_output.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return normalized_output


def _coerce_int_pair(raw_pair: object, *, context: str) -> Tuple[int, int]:
    if not isinstance(raw_pair, (list, tuple)) or len(raw_pair) != 2:
        raise PebResolutionError(f"{context} must be a 2-item pair, got: {raw_pair!r}")

    try:
        return int(raw_pair[0]), int(raw_pair[1])
    except (TypeError, ValueError) as exc:
        raise PebResolutionError(f"{context} must contain integer-compatible values: {raw_pair!r}") from exc


def _coerce_int_pair_list(raw_pairs: object, *, context: str) -> List[Tuple[int, int]]:
    if not isinstance(raw_pairs, list):
        raise PebResolutionError(f"{context} must be a list of 2-item pairs")

    pairs: List[Tuple[int, int]] = []
    for index, raw_pair in enumerate(raw_pairs):
        pairs.append(_coerce_int_pair(raw_pair, context=f"{context}[{index}]"))
    return pairs


def _extract_map_to_product_smiles(s0_payload: Dict[str, object]) -> Dict[int, int]:
    raw_mapping: Optional[Any] = None

    smiles_atom_mapping = s0_payload.get("smiles_atom_mapping")
    if isinstance(smiles_atom_mapping, dict):
        raw_mapping = smiles_atom_mapping.get("map_to_product_smiles")

    if raw_mapping is None:
        raw_mapping = s0_payload.get("map_to_product_smiles")

    if raw_mapping is None:
        raise PebResolutionError("S0 atom_map_smiles.json missing map_to_product_smiles")
    if not isinstance(raw_mapping, dict):
        raise PebResolutionError("S0 map_to_product_smiles must be a JSON object")

    mapping: Dict[int, int] = {}
    for raw_map_idx, raw_smiles_idx in raw_mapping.items():
        try:
            map_idx = int(raw_map_idx)
            smiles_idx = int(raw_smiles_idx)
        except (TypeError, ValueError) as exc:
            raise PebResolutionError(
                f"Invalid map_to_product_smiles entry: {raw_map_idx!r} -> {raw_smiles_idx!r}"
            ) from exc
        mapping[map_idx] = smiles_idx

    return mapping


def _extract_forming_bonds_map_space(s0_payload: Dict[str, object]) -> List[Tuple[int, int]]:
    raw_map_space = s0_payload.get("forming_bonds_map_space")
    if raw_map_space is not None:
        pairs = _coerce_int_pair_list(raw_map_space, context="forming_bonds_map_space")
        if pairs:
            return pairs

    raw_annotated = s0_payload.get("forming_bonds_annotated")
    if raw_annotated is not None:
        if not isinstance(raw_annotated, list):
            raise PebResolutionError("forming_bonds_annotated must be a list")

        annotated_pairs: List[Tuple[int, int]] = []
        for index, raw_entry in enumerate(raw_annotated):
            if not isinstance(raw_entry, dict):
                raise PebResolutionError(f"forming_bonds_annotated[{index}] must be an object")
            if "map_space" not in raw_entry:
                raise PebResolutionError(f"forming_bonds_annotated[{index}] missing map_space")
            annotated_pairs.append(
                _coerce_int_pair(raw_entry["map_space"], context=f"forming_bonds_annotated[{index}].map_space")
            )
        if annotated_pairs:
            return annotated_pairs

    if "forming_bonds" in s0_payload:
        if s0_payload.get("forming_bonds_index_space") != "map_space":
            raise PebResolutionError(
                "Legacy S0 forming_bonds can only be used when forming_bonds_index_space == 'map_space'"
            )
        return _coerce_int_pair_list(s0_payload.get("forming_bonds"), context="forming_bonds")

    raise PebResolutionError("S0 atom_map_smiles.json missing forming bonds in map space")


def _collect_lewis_acid_additive_indices(smiles_to_xyz_payload: Dict[str, object]) -> List[int]:
    lewis_acid = smiles_to_xyz_payload.get("lewis_acid")
    if lewis_acid is None:
        return []
    if not isinstance(lewis_acid, dict):
        raise PebResolutionError("S1 smiles_to_xyz_map lewis_acid payload must be an object")

    raw_indices = lewis_acid.get("additive_atom_indices")
    if raw_indices is None:
        return []
    if not isinstance(raw_indices, list):
        raise PebResolutionError("lewis_acid.additive_atom_indices must be a list")

    additive_indices: List[int] = []
    for index, raw_value in enumerate(raw_indices):
        try:
            additive_indices.append(int(raw_value))
        except (TypeError, ValueError) as exc:
            raise PebResolutionError(
                f"lewis_acid.additive_atom_indices[{index}] must be integer-compatible"
            ) from exc

    return additive_indices


def _validate_xyz_index(xyz_idx: int, atom_count: int, product_xyz_path: Path) -> None:
    if xyz_idx < 0 or xyz_idx >= atom_count:
        raise PebResolutionError(
            f"XYZ idx {xyz_idx} out of range for {product_xyz_path} (atom_count={atom_count})"
        )


def load_s0_atom_map_smiles(path: Path) -> Dict[str, object]:
    """Load S0 atom_map_smiles.json.

    Args:
        path: Path to S0 atom_map_smiles.json.

    Returns:
        Parsed JSON payload as a dictionary.

    Raises:
        PebResolutionError: If the file is missing, unreadable, invalid JSON,
            or does not contain a top-level object.
    """

    normalized_path = normalize_path(path)
    return _load_json_object(
        normalized_path,
        missing_message=f"S0 atom_map_smiles.json not found: {normalized_path}",
        label="S0 atom_map_smiles.json",
    )


def load_s1_smiles_to_xyz_map(path: Path, product_xyz: Path) -> Dict[str, object]:
    """Load and gate the S1 smiles_to_xyz_map.json sidecar.

    Args:
        path: Path to S1 smiles_to_xyz_map.json.
        product_xyz: Path to the S1 product XYZ file used for optional hash checks.

    Returns:
        Parsed JSON payload as a dictionary.

    Raises:
        PebResolutionError: If the file is missing, unreadable, invalid, has
            non-high confidence, or uses a rejected mapping source.
    """

    normalized_path = normalize_path(path)
    normalized_product_xyz = normalize_path(product_xyz)
    payload = _load_json_object(
        normalized_path,
        missing_message=f"S1 smiles_to_xyz_map.json not found: {normalized_path}",
        label="S1 smiles_to_xyz_map.json",
    )

    confidence = payload.get("confidence")
    if confidence != "high":
        raise PebResolutionError(f"S1 smiles_to_xyz_map confidence must be 'high', got {confidence!r}")

    mapping_source: Any = payload.get("mapping_source")
    if mapping_source is None:
        raise PebResolutionError("S1 smiles_to_xyz_map.json missing mapping_source")
    if not isinstance(mapping_source, str):
        raise PebResolutionError(
            f"S1 smiles_to_xyz_map mapping_source must be a string, got {type(mapping_source).__name__}"
        )
    if mapping_source in REJECTED_S1_MAPPING_SOURCES:
        raise PebResolutionError(
            f"S1 smiles_to_xyz_map mapping_source {mapping_source!r} is rejected for PEB resolution"
        )

    expected_hash = payload.get("xyz_sha256")
    if expected_hash is not None:
        actual_hash = _compute_file_hash(normalized_product_xyz)
        if actual_hash is not None and str(expected_hash) != actual_hash:
            raise PebResolutionError(
                "S1 smiles_to_xyz_map xyz_sha256 mismatch for "
                f"{normalized_product_xyz}: sidecar={expected_hash} actual={actual_hash}"
            )

    # V4 generation sidecars contain a complete atom table and therefore gate
    # the whole XYZ sequence.  Older validation fixtures may intentionally
    # contain only mapped heavy atoms; their resolved indices are still gated
    # below by resolve_peb_forming_bonds.
    if payload.get("schema_version") == "s1_atom_mapping_v1":
        verify_atom_mapping_table(payload, normalized_product_xyz)

    return payload


def build_product_smiles_to_xyz(smiles_to_xyz_payload: Dict[str, object]) -> Dict[int, int]:
    """Build a product-SMILES-index to XYZ-index map from the S1 sidecar.

    Args:
        smiles_to_xyz_payload: Parsed S1 smiles_to_xyz_map.json payload.

    Returns:
        Dictionary mapping 0-based product SMILES atom indices to 0-based XYZ indices.

    Raises:
        PebResolutionError: If the payload lacks a valid atoms list or contains
            malformed atom entries.
    """

    raw_atoms = smiles_to_xyz_payload.get("atoms")
    if raw_atoms is None:
        raise PebResolutionError("S1 smiles_to_xyz_map.json missing atoms")
    if not isinstance(raw_atoms, list):
        raise PebResolutionError("S1 smiles_to_xyz_map.json atoms must be a list")

    product_smiles_to_xyz: Dict[int, int] = {}
    for index, raw_entry in enumerate(raw_atoms):
        if not isinstance(raw_entry, dict):
            raise PebResolutionError(f"atoms[{index}] must be an object")

        smiles_idx = raw_entry.get("smiles_idx")
        if smiles_idx is None:
            continue

        raw_xyz_idx = raw_entry.get("xyz_idx")
        if raw_xyz_idx is None:
            raise PebResolutionError(f"atoms[{index}] missing xyz_idx for smiles_idx {smiles_idx!r}")

        try:
            smiles_idx_int = int(smiles_idx)
            xyz_idx_int = int(raw_xyz_idx)
        except (TypeError, ValueError) as exc:
            raise PebResolutionError(
                f"atoms[{index}] has non-integer smiles_idx/xyz_idx: {smiles_idx!r}, {raw_xyz_idx!r}"
            ) from exc

        existing_xyz_idx = product_smiles_to_xyz.get(smiles_idx_int)
        if existing_xyz_idx is not None and existing_xyz_idx != xyz_idx_int:
            raise PebResolutionError(
                f"Duplicate smiles_idx {smiles_idx_int} maps to multiple xyz_idx values: {existing_xyz_idx} and {xyz_idx_int}"
            )

        product_smiles_to_xyz[smiles_idx_int] = xyz_idx_int

    return product_smiles_to_xyz


def resolve_peb_forming_bonds(
    *,
    atom_map_smiles_path: Path,
    smiles_to_xyz_map_path: Path,
    product_xyz_path: Path,
    allow_legacy_index_space: bool = False,
) -> PebMappingResult:
    """Resolve PEB forming bonds from S0/S1 sidecars into product XYZ space.

    Args:
        atom_map_smiles_path: Path to S0 atom_map_smiles.json.
        smiles_to_xyz_map_path: Path to S1 smiles_to_xyz_map.json.
        product_xyz_path: Path to S1 product_min.xyz.
        allow_legacy_index_space: When True, bypass fail-fast index-space
            validation for legacy migration/validation workflows.

    Returns:
        Frozen PEB mapping result containing map-space, SMILES-space, and XYZ-space bonds.

    Raises:
        PebResolutionError: If any required sidecar is missing, malformed,
            low-confidence, unresolved, or fails PEB validation.
    """

    normalized_s0_path = normalize_path(atom_map_smiles_path)
    normalized_s1_path = normalize_path(smiles_to_xyz_map_path)
    normalized_product_xyz_path = normalize_path(product_xyz_path)

    s0_payload = load_s0_atom_map_smiles(normalized_s0_path)
    s1_payload = load_s1_smiles_to_xyz_map(normalized_s1_path, normalized_product_xyz_path)

    s0_idx_space = s0_payload.get("product_smiles_idx_space")
    s1_idx_space = s1_payload.get("product_smiles_idx_space")

    if not allow_legacy_index_space:
        if s0_idx_space is None:
            raise PebResolutionError(
                "S0 atom_map_smiles.json missing 'product_smiles_idx_space'; cannot verify index space consistency. "
                + "Set allow_legacy_index_space=True for validation/migration only."
            )
        if s1_idx_space is None:
            raise PebResolutionError(
                "S1 smiles_to_xyz_map.json missing 'product_smiles_idx_space'; cannot verify index space consistency."
            )
        if s0_idx_space != s1_idx_space:
            raise PebResolutionError(
                f"S0/S1 product_smiles_idx_space mismatch: S0='{s0_idx_space}', S1='{s1_idx_space}'"
            )
        if s0_idx_space != EXPECTED_PRODUCT_SMILES_IDX_SPACE:
            raise PebResolutionError(
                f"S0/S1 product_smiles_idx_space must both be '{EXPECTED_PRODUCT_SMILES_IDX_SPACE}'; "
                + f"got S0='{s0_idx_space}', S1='{s1_idx_space}'"
            )
    elif (
        s0_idx_space != s1_idx_space
        or s0_idx_space != EXPECTED_PRODUCT_SMILES_IDX_SPACE
        or s1_idx_space != EXPECTED_PRODUCT_SMILES_IDX_SPACE
    ):
        logger.warning(
            "Index space validation bypassed (allow_legacy_index_space=True): expected '%s'; S0='%s', S1='%s'",
            EXPECTED_PRODUCT_SMILES_IDX_SPACE,
            s0_idx_space,
            s1_idx_space,
        )

    map_to_product_smiles = _extract_map_to_product_smiles(s0_payload)
    forming_bonds_map_space = _extract_forming_bonds_map_space(s0_payload)
    smiles_to_xyz = build_product_smiles_to_xyz(s1_payload)

    try:
        coordinates, elements = read_xyz(normalized_product_xyz_path)
    except (FileNotFoundError, OSError, ValueError, IndexError) as exc:
        raise PebResolutionError(
            f"Failed to read product XYZ for PEB resolution: {normalized_product_xyz_path}: {exc}"
        ) from exc

    atom_count = len(coordinates)

    resolved_product_smiles: List[Tuple[int, int]] = []
    resolved_product_xyz_0based: List[Tuple[int, int]] = []
    resolved_product_xyz_1based: List[Tuple[int, int]] = []

    for left_map, right_map in forming_bonds_map_space:
        left_product_smiles = map_to_product_smiles.get(left_map)
        if left_product_smiles is None:
            raise PebResolutionError(f"Map# {left_map} has no product SMILES idx")

        right_product_smiles = map_to_product_smiles.get(right_map)
        if right_product_smiles is None:
            raise PebResolutionError(f"Map# {right_map} has no product SMILES idx")

        left_xyz = smiles_to_xyz.get(left_product_smiles)
        if left_xyz is None:
            raise PebResolutionError(f"product SMILES idx {left_product_smiles} has no XYZ idx")

        right_xyz = smiles_to_xyz.get(right_product_smiles)
        if right_xyz is None:
            raise PebResolutionError(f"product SMILES idx {right_product_smiles} has no XYZ idx")

        resolved_product_smiles.append((left_product_smiles, right_product_smiles))
        resolved_product_xyz_0based.append((left_xyz, right_xyz))
        resolved_product_xyz_1based.append((left_xyz + 1, right_xyz + 1))

    if len(forming_bonds_map_space) != 2:
        raise PebResolutionError(
            f"Expected exactly 2 forming bonds for PEB resolution, got {len(forming_bonds_map_space)}"
        )

    additive_indices = set(_collect_lewis_acid_additive_indices(s1_payload))

    all_xyz_indices: List[int] = []
    for left_xyz, right_xyz in resolved_product_xyz_0based:
        all_xyz_indices.extend([left_xyz, right_xyz])

    unique_xyz_indices = set(all_xyz_indices)
    if len(unique_xyz_indices) != 4:
        raise PebResolutionError(
            f"PEB forming bonds must span 4 unique XYZ atoms, got {len(unique_xyz_indices)}"
        )

    for xyz_idx in all_xyz_indices:
        _validate_xyz_index(xyz_idx, atom_count, normalized_product_xyz_path)

        element = _normalize_element_symbol(elements[xyz_idx])
        if xyz_idx in additive_indices:
            raise PebResolutionError(f"XYZ idx {xyz_idx} is marked as a Lewis acid additive atom")
        if element in ADDITIVE_OR_H_ELEMENTS:
            raise PebResolutionError(
                f"XYZ idx {xyz_idx} resolves to disallowed additive/hydrogen element {element}"
            )
        if element not in ORGANIC_HEAVY_ELEMENTS:
            raise PebResolutionError(
                f"XYZ idx {xyz_idx} resolves to unexpected non-organic heavy element {element}"
            )

    result = PebMappingResult(
        forming_bonds_map_space=tuple(forming_bonds_map_space),
        forming_bonds_product_smiles=tuple(resolved_product_smiles),
        forming_bonds_product_xyz_0based=tuple(resolved_product_xyz_0based),
        forming_bonds_product_xyz_1based=tuple(resolved_product_xyz_1based),
        confidence=str(s1_payload.get("confidence", "unknown")),
        diagnostics={
            "s0_path": str(normalized_s0_path),
            "s1_path": str(normalized_s1_path),
            "xyz_path": str(normalized_product_xyz_path),
            "atom_count": atom_count,
            "s0_schema_version": s0_payload.get("schema_version"),
            "s0_product_smiles_idx_space": s0_idx_space,
            "s1_mapping_source": s1_payload.get("mapping_source"),
            "s1_confidence": s1_payload.get("confidence"),
            "s1_product_smiles_idx_space": s1_idx_space,
            "allow_legacy_index_space": allow_legacy_index_space,
            "map_to_product_smiles_count": len(map_to_product_smiles),
            "smiles_to_xyz_count": len(smiles_to_xyz),
            "validation_passed": True,
        },
    )
    return result


def write_peb_mapping(result: PebMappingResult, output_path: Path, *, product_xyz: Path) -> Path:
    """Write a resolved PEB mapping JSON payload.

    Args:
        result: Resolved PEB mapping result.
        output_path: Target JSON file path.
        product_xyz: Product XYZ path to record in the output payload.

    Returns:
        The output path that was written.

    Raises:
        PebResolutionError: If the payload cannot be written.
    """

    normalized_output_path = normalize_path(output_path)
    normalized_product_xyz = normalize_path(product_xyz)

    unique_atom_count = len(
        {xyz_idx for pair in result.forming_bonds_product_xyz_0based for xyz_idx in pair}
    )

    payload = {
        "schema_version": "3.0",
        "mapping_source": "s0_smiles_map_plus_s1_smiles_xyz_sidecar",
        "confidence": result.confidence,
        "product_xyz": str(normalized_product_xyz),
        "forming_bonds": {
            "map_space": [list(pair) for pair in result.forming_bonds_map_space],
            "product_smiles_idx": [list(pair) for pair in result.forming_bonds_product_smiles],
            "product_xyz_0based": [list(pair) for pair in result.forming_bonds_product_xyz_0based],
            "product_xyz_1based": [list(pair) for pair in result.forming_bonds_product_xyz_1based],
        },
        "validation": {
            "all_map_numbers_resolved": True,
            "all_smiles_indices_resolved": True,
            "all_xyz_indices_in_range": True,
            "no_additive_atoms_in_forming_bonds": True,
            "forming_bond_atom_count": unique_atom_count,
        },
    }

    normalized_output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with normalized_output_path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
    except OSError as exc:
        raise PebResolutionError(f"Failed to write PEB mapping JSON: {normalized_output_path}: {exc}") from exc

    return normalized_output_path
