"""Comprehensive unit tests for the PEB atom-mapping resolver."""

# pyright: reportDeprecated=false, reportMissingImports=false, reportUnknownVariableType=false, reportUnknownMemberType=false, reportUnknownArgumentType=false, reportUnknownParameterType=false, reportExplicitAny=false, reportAny=false

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pytest

from rph_core.utils.atom_mapping import (
    PebMappingResult,
    PebResolutionError,
    build_product_smiles_to_xyz,
    load_s0_atom_map_smiles,
    load_s1_smiles_to_xyz_map,
    resolve_peb_forming_bonds,
    write_peb_mapping,
)


def _write_s0_atom_map_smiles(
    path: Path,
    *,
    map_to_product_smiles: Optional[Dict[str, int]] = None,
    forming_bonds_map_space: Optional[List[List[int]]] = None,
    forming_bonds_annotated: Optional[List[Dict[str, List[int]]]] = None,
    forming_bonds: Optional[List[List[int]]] = None,
    forming_bonds_index_space: Optional[str] = None,
    product_smiles_idx_space: Optional[str] = "geometry_product_smiles_idx",
    schema_version: str = "3.0",
) -> Path:
    """Write a synthetic S0 atom_map_smiles.json."""
    payload: Dict[str, object] = {
        "schema_version": schema_version,
        "index_spaces": {
            "map_numbers": "1-based atom map numbers from mapped SMILES",
            "geometry_product_smiles_idx": "0-based RDKit atom index in S1 geometry product SMILES",
        },
        "smiles_atom_mapping": {
            "map_to_product_smiles": map_to_product_smiles or {},
            "product_smiles_to_map": {},
        },
    }
    if product_smiles_idx_space is not None:
        payload["product_smiles_idx_space"] = product_smiles_idx_space
    if forming_bonds_map_space is not None:
        payload["forming_bonds_map_space"] = forming_bonds_map_space
    if forming_bonds_annotated is not None:
        payload["forming_bonds_annotated"] = forming_bonds_annotated
    if forming_bonds is not None:
        payload["forming_bonds"] = forming_bonds
    if forming_bonds_index_space is not None:
        payload["forming_bonds_index_space"] = forming_bonds_index_space
    path.parent.mkdir(parents=True, exist_ok=True)
    _ = path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def _write_s1_smiles_to_xyz_map(
    path: Path,
    *,
    atoms: Optional[List[Dict[str, object]]] = None,
    confidence: str = "high",
    mapping_source: str = "rdkit_generation_sidecar",
    lewis_acid: Optional[Dict[str, List[int]]] = None,
    product_smiles_idx_space: Optional[str] = "geometry_product_smiles_idx",
) -> Path:
    """Write a synthetic S1 smiles_to_xyz_map.json."""
    payload: Dict[str, object] = {
        "schema_version": "3.0",
        "mapping_source": mapping_source,
        "confidence": confidence,
        "atoms": atoms or [],
    }
    if product_smiles_idx_space is not None:
        payload["product_smiles_idx_space"] = product_smiles_idx_space
    if lewis_acid is not None:
        payload["lewis_acid"] = lewis_acid
    path.parent.mkdir(parents=True, exist_ok=True)
    _ = path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def _write_product_xyz(path: Path, elements: List[str]) -> Path:
    """Write a minimal valid XYZ file with the given element symbols."""
    lines = [str(len(elements)), "test product"]
    for element in elements:
        lines.append(f"{element}  0.000  0.000  0.000")
    path.parent.mkdir(parents=True, exist_ok=True)
    _ = path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_case(
    tmp_path: Path,
    *,
    s0_kwargs: Dict[str, Any],
    s1_kwargs: Dict[str, Any],
    elements: List[str],
) -> Tuple[Path, Path, Path]:
    s0_path = _write_s0_atom_map_smiles(tmp_path / "S0" / "atom_map_smiles.json", **s0_kwargs)
    s1_path = _write_s1_smiles_to_xyz_map(tmp_path / "S1" / "smiles_to_xyz_map.json", **s1_kwargs)
    xyz_path = _write_product_xyz(tmp_path / "S1" / "product_min.xyz", elements)
    return s0_path, s1_path, xyz_path


def _resolve(
    s0_path: Path,
    s1_path: Path,
    xyz_path: Path,
    *,
    allow_legacy_index_space: bool = False,
) -> PebMappingResult:
    return resolve_peb_forming_bonds(
        atom_map_smiles_path=s0_path,
        smiles_to_xyz_map_path=s1_path,
        product_xyz_path=xyz_path,
        allow_legacy_index_space=allow_legacy_index_space,
    )


def _basic_four_atom_atoms() -> List[Dict[str, object]]:
    return [
        {"smiles_idx": 0, "xyz_idx": 0, "element": "C", "type": "organic_heavy"},
        {"smiles_idx": 1, "xyz_idx": 1, "element": "C", "type": "organic_heavy"},
        {"smiles_idx": 2, "xyz_idx": 2, "element": "C", "type": "organic_heavy"},
        {"smiles_idx": 3, "xyz_idx": 3, "element": "C", "type": "organic_heavy"},
    ]


class TestPebResolver:
    def test_peb_resolves_map_to_xyz_correctly(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 2, "11": 3, "3": 7, "8": 6},
                "forming_bonds_map_space": [[1, 11], [3, 8]],
            },
            s1_kwargs={
                "atoms": [
                    {"smiles_idx": 2, "xyz_idx": 2, "element": "C", "type": "organic_heavy"},
                    {"smiles_idx": 3, "xyz_idx": 3, "element": "C", "type": "organic_heavy"},
                    {"smiles_idx": 7, "xyz_idx": 6, "element": "C", "type": "organic_heavy"},
                    {"smiles_idx": 6, "xyz_idx": 21, "element": "C", "type": "organic_heavy"},
                ],
            },
            elements=["C"] * 22,
        )

        assert load_s0_atom_map_smiles(s0_path)["schema_version"] == "3.0"
        assert load_s1_smiles_to_xyz_map(s1_path, xyz_path)["confidence"] == "high"

        result = _resolve(s0_path, s1_path, xyz_path)

        assert isinstance(result, PebMappingResult)
        assert result.forming_bonds_map_space == ((1, 11), (3, 8))
        assert result.forming_bonds_product_smiles == ((2, 3), (7, 6))
        assert result.forming_bonds_product_xyz_0based == ((2, 3), (6, 21))
        assert result.forming_bonds_product_xyz_1based == ((3, 4), (7, 22))
        assert result.confidence == "high"

    def test_peb_fails_when_map_number_unresolved(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0},
                "forming_bonds_map_space": [[1, 99]],
            },
            s1_kwargs={"atoms": [{"smiles_idx": 0, "xyz_idx": 0, "element": "C", "type": "organic_heavy"}]},
            elements=["C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "Map# 99" in str(excinfo.value)

    def test_peb_fails_when_smiles_idx_unresolved(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 5, "2": 6},
                "forming_bonds_map_space": [[1, 2]],
            },
            s1_kwargs={"atoms": [{"smiles_idx": 6, "xyz_idx": 1, "element": "C", "type": "organic_heavy"}]},
            elements=["C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "product SMILES idx 5 has no XYZ idx" in str(excinfo.value)

    def test_peb_fails_when_xyz_out_of_range(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
            },
            s1_kwargs={
                "atoms": [
                    {"smiles_idx": 0, "xyz_idx": 0, "element": "C", "type": "organic_heavy"},
                    {"smiles_idx": 1, "xyz_idx": 1, "element": "C", "type": "organic_heavy"},
                    {"smiles_idx": 2, "xyz_idx": 2, "element": "C", "type": "organic_heavy"},
                    {"smiles_idx": 3, "xyz_idx": 100, "element": "C", "type": "organic_heavy"},
                ],
            },
            elements=["C"] * 10,
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "XYZ idx 100 out of range" in str(excinfo.value)

    def test_peb_fails_when_additive_atom_in_bond(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
            },
            s1_kwargs={"atoms": _basic_four_atom_atoms()},
            elements=["C", "C", "C", "Mg"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "disallowed additive/hydrogen element Mg" in str(excinfo.value)

    def test_peb_fails_when_hydrogen_in_bond(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
            },
            s1_kwargs={"atoms": _basic_four_atom_atoms()},
            elements=["C", "C", "H", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "disallowed additive/hydrogen element H" in str(excinfo.value)

    def test_peb_rejects_low_confidence_sidecar(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
            },
            s1_kwargs={
                "atoms": _basic_four_atom_atoms(),
                "confidence": "low",
                "mapping_source": "assumed_heavy_order",
            },
            elements=["C", "C", "C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "confidence must be 'high'" in str(excinfo.value)

    def test_peb_rejects_missing_s0_index_space(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
                "product_smiles_idx_space": None,
            },
            s1_kwargs={"atoms": _basic_four_atom_atoms()},
            elements=["C", "C", "C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "S0 atom_map_smiles.json missing 'product_smiles_idx_space'" in str(excinfo.value)

    def test_peb_rejects_missing_s1_index_space(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
            },
            s1_kwargs={
                "atoms": _basic_four_atom_atoms(),
                "product_smiles_idx_space": None,
            },
            elements=["C", "C", "C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "S1 smiles_to_xyz_map.json missing 'product_smiles_idx_space'" in str(excinfo.value)

    def test_peb_rejects_mismatched_index_space(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
            },
            s1_kwargs={
                "atoms": _basic_four_atom_atoms(),
                "product_smiles_idx_space": "mapped_product_smiles_idx",
            },
            elements=["C", "C", "C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "S0/S1 product_smiles_idx_space mismatch" in str(excinfo.value)

    def test_peb_rejects_noncanonical_matching_index_space(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
                "product_smiles_idx_space": "mapped_product_smiles_idx",
            },
            s1_kwargs={
                "atoms": _basic_four_atom_atoms(),
                "product_smiles_idx_space": "mapped_product_smiles_idx",
            },
            elements=["C", "C", "C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "must both be 'geometry_product_smiles_idx'" in str(excinfo.value)

    def test_peb_legacy_mode_allows_index_space_mismatch(self, caplog: pytest.LogCaptureFixture, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
                "product_smiles_idx_space": None,
            },
            s1_kwargs={"atoms": _basic_four_atom_atoms()},
            elements=["C", "C", "C", "C"],
        )

        with caplog.at_level("WARNING"):
            result = _resolve(
                s0_path,
                s1_path,
                xyz_path,
                allow_legacy_index_space=True,
            )

        assert result.forming_bonds_product_xyz_0based == ((0, 1), (2, 3))
        assert "allow_legacy_index_space=True" in caplog.text

    def test_peb_rejects_assumed_heavy_order_source(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
            },
            s1_kwargs={
                "atoms": _basic_four_atom_atoms(),
                "mapping_source": "assumed_heavy_order",
            },
            elements=["C", "C", "C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "mapping_source 'assumed_heavy_order' is rejected" in str(excinfo.value)

    def test_peb_rejects_kabsch_hungarian_source(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
            },
            s1_kwargs={
                "atoms": _basic_four_atom_atoms(),
                "mapping_source": "kabsch_hungarian",
            },
            elements=["C", "C", "C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "mapping_source 'kabsch_hungarian' is rejected" in str(excinfo.value)

    def test_peb_rejects_substructure_match_source(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
            },
            s1_kwargs={
                "atoms": _basic_four_atom_atoms(),
                "mapping_source": "substructure_match",
            },
            elements=["C", "C", "C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "mapping_source 'substructure_match' is rejected" in str(excinfo.value)

    def test_peb_rejects_legacy_forming_bonds_without_index_space(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds": [[1, 2], [3, 4]],
            },
            s1_kwargs={"atoms": _basic_four_atom_atoms()},
            elements=["C", "C", "C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "forming_bonds_index_space == 'map_space'" in str(excinfo.value)

    def test_peb_accepts_legacy_forming_bonds_with_map_space(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds": [[1, 2], [3, 4]],
                "forming_bonds_index_space": "map_space",
            },
            s1_kwargs={"atoms": _basic_four_atom_atoms()},
            elements=["C", "C", "C", "C"],
        )

        result = _resolve(s0_path, s1_path, xyz_path)

        assert result.forming_bonds_map_space == ((1, 2), (3, 4))
        assert result.forming_bonds_product_xyz_0based == ((0, 1), (2, 3))

    def test_peb_rejects_wrong_bond_count(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1},
                "forming_bonds_map_space": [[1, 2]],
            },
            s1_kwargs={
                "atoms": [
                    {"smiles_idx": 0, "xyz_idx": 0, "element": "C", "type": "organic_heavy"},
                    {"smiles_idx": 1, "xyz_idx": 1, "element": "C", "type": "organic_heavy"},
                ],
            },
            elements=["C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "Expected exactly 2 forming bonds" in str(excinfo.value)

    def test_peb_rejects_duplicate_atoms_in_bonds(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2},
                "forming_bonds_map_space": [[1, 2], [2, 3]],
            },
            s1_kwargs={
                "atoms": [
                    {"smiles_idx": 0, "xyz_idx": 0, "element": "C", "type": "organic_heavy"},
                    {"smiles_idx": 1, "xyz_idx": 1, "element": "C", "type": "organic_heavy"},
                    {"smiles_idx": 2, "xyz_idx": 2, "element": "C", "type": "organic_heavy"},
                ],
            },
            elements=["C", "C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "must span 4 unique XYZ atoms" in str(excinfo.value)

    def test_peb_rejects_la_additive_in_forming_bonds(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
            },
            s1_kwargs={
                "atoms": _basic_four_atom_atoms(),
                "lewis_acid": {"additive_atom_indices": [3]},
            },
            elements=["C", "C", "C", "C"],
        )

        with pytest.raises(PebResolutionError) as excinfo:
            _resolve(s0_path, s1_path, xyz_path)

        assert "Lewis acid additive atom" in str(excinfo.value)

    def test_peb_missing_s0_file(self, tmp_path: Path):
        missing_s0_path = tmp_path / "missing_atom_map_smiles.json"

        with pytest.raises(PebResolutionError) as excinfo:
            load_s0_atom_map_smiles(missing_s0_path)

        assert "S0 atom_map_smiles.json not found" in str(excinfo.value)

    def test_peb_missing_s1_file(self, tmp_path: Path):
        xyz_path = _write_product_xyz(tmp_path / "product_min.xyz", ["C"])
        missing_s1_path = tmp_path / "missing_smiles_to_xyz_map.json"

        with pytest.raises(PebResolutionError) as excinfo:
            load_s1_smiles_to_xyz_map(missing_s1_path, xyz_path)

        assert "S1 smiles_to_xyz_map.json not found" in str(excinfo.value)

    def test_write_peb_mapping_produces_valid_json(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_map_space": [[1, 2], [3, 4]],
            },
            s1_kwargs={"atoms": _basic_four_atom_atoms()},
            elements=["C", "C", "C", "C"],
        )
        result = _resolve(s0_path, s1_path, xyz_path)
        output_path = tmp_path / "S2_Retro" / "atom_mapping_peb.json"

        written_path = write_peb_mapping(result, output_path, product_xyz=xyz_path)
        payload = json.loads(written_path.read_text(encoding="utf-8"))

        assert written_path == output_path
        assert payload["schema_version"] == "3.0"
        assert payload["mapping_source"] == "s0_smiles_map_plus_s1_smiles_xyz_sidecar"
        assert payload["forming_bonds"]["product_xyz_0based"] == [[0, 1], [2, 3]]
        assert payload["forming_bonds"]["product_xyz_1based"] == [[1, 2], [3, 4]]
        assert payload["validation"]["forming_bond_atom_count"] == 4

    def test_build_product_smiles_to_xyz_skips_additive(self):
        payload = {
            "atoms": [
                {"smiles_idx": 0, "xyz_idx": 0, "element": "C", "type": "organic_heavy"},
                {"smiles_idx": None, "xyz_idx": 4, "element": "Mg", "type": "additive"},
                {"smiles_idx": 3, "xyz_idx": 5, "element": "O", "type": "organic_heavy"},
                {"smiles_idx": None, "xyz_idx": 6, "element": "H", "type": "hydrogen"},
            ]
        }

        result = build_product_smiles_to_xyz(payload)

        assert result == {0: 0, 3: 5}

    def test_forming_bonds_annotated_fallback(self, tmp_path: Path):
        s0_path, s1_path, xyz_path = _write_case(
            tmp_path,
            s0_kwargs={
                "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                "forming_bonds_annotated": [
                    {"map_space": [1, 2]},
                    {"map_space": [3, 4]},
                ],
            },
            s1_kwargs={"atoms": _basic_four_atom_atoms()},
            elements=["C", "C", "C", "C"],
        )

        result = _resolve(s0_path, s1_path, xyz_path)

        assert result.forming_bonds_map_space == ((1, 2), (3, 4))
        assert result.forming_bonds_product_smiles == ((0, 1), (2, 3))
