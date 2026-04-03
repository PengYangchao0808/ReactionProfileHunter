from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
import csv
import json
import re

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdDetermineBonds

from rph_core.utils.file_io import read_xyz
from rph_core.utils.tsv_dataset import ReactionRecord


class CleanerAdapterError(Exception):
    pass


def load_cleaner_records(
    path: Path,
    reaction_profiles: Optional[Dict[str, Dict[str, Any]]] = None,
    delimiter: str = "\t",
) -> List[ReactionRecord]:
    if not path.exists():
        raise CleanerAdapterError(f"Cleaner file not found: {path}")

    rows = list(_iter_cleaner_rows(path, delimiter=delimiter))
    records: List[ReactionRecord] = []
    profiles = reaction_profiles or {}

    for idx, row in enumerate(rows, start=1):
        record = convert_cleaner_row_to_record(
            row=row,
            row_index=idx,
            reaction_profiles=profiles,
        )
        if record is not None:
            records.append(record)

    if not records:
        raise CleanerAdapterError(f"No valid cleaner records found in {path}")

    return records


def convert_cleaner_row_to_record(
    row: Dict[str, Any],
    row_index: int,
    reaction_profiles: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Optional[ReactionRecord]:
    rx_id = _first_nonempty(row, "rx_id", "reaction_id", "record_id", "id") or f"cleaner_{row_index:06d}"

    precursor_smiles = _first_nonempty(
        row,
        "precursor_smiles",
        "reactant_smiles",
        "substrate_smiles",
    )
    if not precursor_smiles:
        return None

    map_status = (_first_nonempty(row, "map_status") or "").strip().upper()
    map_confidence = _to_float(_first_nonempty(row, "map_confidence", "map_score", "mapping_score"))
    use_mapped = map_status == "OK" and map_confidence is not None and map_confidence >= 0.8

    mapped_product_smiles = _first_nonempty(
        row,
        "mapped_product_smiles",
        "product_smiles_mapped",
        "mapped_product",
    )
    mapped_precursor_smiles = _first_nonempty(
        row,
        "mapped_precursor_smiles",
        "precursor_smiles_mapped",
        "mapped_reactant_smiles",
        "reactant_smiles_mapped",
    )
    mapped_smiles = mapped_product_smiles or mapped_precursor_smiles or precursor_smiles

    mapped_reaction = _first_nonempty(row, "rxn_smiles_mapped", "reaction_smiles_mapped")
    if not mapped_product_smiles:
        mapped_product_smiles = extract_product_smiles(mapped_reaction)
    if not mapped_smiles:
        mapped_smiles = extract_reactant_smiles(mapped_reaction)

    index_source_smiles = mapped_product_smiles or mapped_smiles

    formed_map_pairs: List[Tuple[int, int]] = []
    formed_index_pairs: List[Tuple[int, int]] = []
    bond_source = "none"

    if use_mapped:
        formed_map_pairs = parse_formed_pairs_from_core_bond_changes(
            _first_nonempty(row, "core_bond_changes")
        )
        formed_index_pairs = map_pairs_to_internal_indices(index_source_smiles, formed_map_pairs)
        bond_source = "core_bond_changes"
    else:
        smarts = _first_nonempty(row, "reaction_smarts", "smarts", "fallback_smarts")
        formed_map_pairs = extract_formed_pairs_from_reaction_smarts(smarts)
        formed_index_pairs = map_pairs_to_internal_indices(index_source_smiles, formed_map_pairs)
        bond_source = "smarts_fallback"

    reaction_type = _first_nonempty(row, "reaction_type", "rxn_type", "reaction_family")
    profile_key = match_reaction_profile_key(reaction_type, reaction_profiles or {})

    raw = {k: str(v) for k, v in row.items() if v is not None and str(v).strip()}
    raw["map_status"] = map_status
    if map_confidence is not None:
        raw["map_confidence"] = f"{map_confidence:.6g}"
    raw["bond_change_source"] = bond_source
    if mapped_product_smiles:
        raw["mapped_product_smiles"] = mapped_product_smiles
    raw["formed_bond_map_pairs"] = _pairs_to_str(formed_map_pairs)
    raw["formed_bond_index_pairs"] = _pairs_to_str(formed_index_pairs)
    raw["forming_bonds"] = raw["formed_bond_index_pairs"]
    if formed_index_pairs:
        raw["forming_bonds_index_base"] = "0"
        raw["index_base"] = "0"
    if profile_key:
        raw["reaction_profile"] = profile_key

    return ReactionRecord(
        rx_id=str(rx_id),
        precursor_smiles=str(precursor_smiles),
        product_smiles_main=_first_nonempty(row, "product_smiles_main", "product_smiles"),
        solvent=_first_nonempty(row, "solvent"),
        base=_first_nonempty(row, "base"),
        temp_celsius=_first_nonempty(row, "temp_celsius", "temperature"),
        yield_=_first_nonempty(row, "yield", "yield_percent"),
        ee=_first_nonempty(row, "ee"),
        dr_major=_first_nonempty(row, "dr_major"),
        dr_minor=_first_nonempty(row, "dr_minor"),
        raw=raw,
    )


def parse_formed_pairs_from_core_bond_changes(core_bond_changes: Optional[str]) -> List[Tuple[int, int]]:
    if not core_bond_changes:
        return []

    pairs: List[Tuple[int, int]] = []
    for item in core_bond_changes.split(";"):
        piece = item.strip()
        if not piece or ":" not in piece:
            continue
        bond_part, change_type = piece.split(":", 1)
        if change_type.strip().lower() != "formed":
            continue
        match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", bond_part)
        if not match:
            continue
        a, b = int(match.group(1)), int(match.group(2))
        if a != b:
            pairs.append((a, b))
    return _deduplicate_pairs(pairs)


def map_pairs_to_internal_indices(
    mapped_smiles: Optional[str],
    map_pairs: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    if not mapped_smiles or not map_pairs:
        return []

    mol = Chem.MolFromSmiles(mapped_smiles)
    if mol is None:
        return []

    map_to_idx: Dict[int, int] = {}
    for atom in mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0 and map_num not in map_to_idx:
            map_to_idx[map_num] = atom.GetIdx()

    idx_pairs: List[Tuple[int, int]] = []
    for a_map, b_map in map_pairs:
        if a_map in map_to_idx and b_map in map_to_idx:
            idx_pairs.append((map_to_idx[a_map], map_to_idx[b_map]))
    return _deduplicate_pairs(idx_pairs)


def map_pairs_to_xyz_indices(
    mapped_smiles: Optional[str],
    map_pairs: Sequence[Tuple[int, int]],
    xyz_file: Path,
) -> List[Tuple[int, int]]:
    if not mapped_smiles or not map_pairs:
        return []

    mol = Chem.MolFromSmiles(mapped_smiles)
    if mol is None:
        return []

    map_to_query_idx: Dict[int, int] = {}
    for atom in mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0 and map_num not in map_to_query_idx:
            map_to_query_idx[map_num] = atom.GetIdx()

    if not map_to_query_idx:
        return []

    coords, symbols = read_xyz(Path(xyz_file))
    xyz_lines = [str(len(symbols)), "xyz"]
    for symbol, coord in zip(symbols, coords):
        xyz_lines.append(f"{symbol} {float(coord[0]):.10f} {float(coord[1]):.10f} {float(coord[2]):.10f}")
    xyz_block = "\n".join(xyz_lines) + "\n"

    xyz_mol = Chem.MolFromXYZBlock(xyz_block)
    if xyz_mol is None:
        return []

    try:
        rdDetermineBonds.DetermineBonds(xyz_mol)
    except Exception:
        try:
            rdDetermineBonds.DetermineConnectivity(xyz_mol)
        except Exception:
            return []

    query_mol = Chem.Mol(mol)
    for atom in query_mol.GetAtoms():
        atom.SetAtomMapNum(0)

    matches = xyz_mol.GetSubstructMatches(query_mol, uniquify=False, useChirality=False)
    if not matches:
        if query_mol.GetNumAtoms() != xyz_mol.GetNumAtoms():
            return []
        mapped_symbols = [atom.GetSymbol() for atom in query_mol.GetAtoms()]
        if mapped_symbols != symbols:
            return []
        matches = [tuple(range(query_mol.GetNumAtoms()))]

    # V5.1 Symmetry Disambiguation
    best_match = matches[0]
    if len(matches) > 1:
        from rph_core.utils.geometry_tools import GeometryUtils
        best_score = float('inf')
        for match in matches:
            query_to_xyz_cand = {q_idx: int(x_idx) for q_idx, x_idx in enumerate(match)}
            penalty = 0.0
            valid_pairs = 0
            for a_map, b_map in map_pairs:
                a_query = map_to_query_idx.get(int(a_map))
                b_query = map_to_query_idx.get(int(b_map))
                if a_query is None or b_query is None:
                    continue
                a_xyz = query_to_xyz_cand.get(a_query)
                b_xyz = query_to_xyz_cand.get(b_query)
                if a_xyz is None or b_xyz is None:
                    continue
                # For product forming bonds, atoms should be close. We penalize larger distances.
                dist = GeometryUtils.calculate_distance(coords, a_xyz, b_xyz)
                penalty += dist
                valid_pairs += 1
            if valid_pairs > 0 and penalty < best_score:
                best_score = penalty
                best_match = match

    query_to_xyz: Dict[int, int] = {q_idx: int(x_idx) for q_idx, x_idx in enumerate(best_match)}

    xyz_pairs: List[Tuple[int, int]] = []
    for a_map, b_map in map_pairs:
        a_query = map_to_query_idx.get(int(a_map))
        b_query = map_to_query_idx.get(int(b_map))
        if a_query is None or b_query is None:
            continue
        a_xyz = query_to_xyz.get(a_query)
        b_xyz = query_to_xyz.get(b_query)
        if a_xyz is None or b_xyz is None:
            continue
        xyz_pairs.append((a_xyz, b_xyz))

    return _deduplicate_pairs(xyz_pairs)

def get_map_to_xyz_dict(
    mapped_smiles: Optional[str],
    xyz_file: Path,
) -> Dict[int, int]:
    """Get mapping from MapId to MolIdx for V5.1 contracts"""
    if not mapped_smiles:
        return {}

    mol = Chem.MolFromSmiles(mapped_smiles)
    if mol is None:
        return {}

    map_to_query_idx: Dict[int, int] = {}
    for atom in mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0 and map_num not in map_to_query_idx:
            map_to_query_idx[map_num] = atom.GetIdx()

    if not map_to_query_idx:
        return {}

    coords, symbols = read_xyz(Path(xyz_file))
    xyz_lines = [str(len(symbols)), "xyz"]
    for symbol, coord in zip(symbols, coords):
        xyz_lines.append(f"{symbol} {float(coord[0]):.10f} {float(coord[1]):.10f} {float(coord[2]):.10f}")
    xyz_block = "\n".join(xyz_lines) + "\n"

    xyz_mol = Chem.MolFromXYZBlock(xyz_block)
    if xyz_mol is None:
        return {}

    try:
        rdDetermineBonds.DetermineBonds(xyz_mol)
    except Exception:
        try:
            rdDetermineBonds.DetermineConnectivity(xyz_mol)
        except Exception:
            return {}

    query_mol = Chem.Mol(mol)
    for atom in query_mol.GetAtoms():
        atom.SetAtomMapNum(0)

    matches = xyz_mol.GetSubstructMatches(query_mol, uniquify=False, useChirality=False)
    if not matches:
        if query_mol.GetNumAtoms() != xyz_mol.GetNumAtoms():
            return {}
        mapped_symbols = [atom.GetSymbol() for atom in query_mol.GetAtoms()]
        if mapped_symbols != symbols:
            return {}
        matches = [tuple(range(query_mol.GetNumAtoms()))]

    query_to_xyz: Dict[int, int] = {q_idx: int(x_idx) for q_idx, x_idx in enumerate(matches[0])}
    
    map_to_xyz: Dict[int, int] = {}
    for map_num, query_idx in map_to_query_idx.items():
        if query_idx in query_to_xyz:
            map_to_xyz[map_num] = query_to_xyz[query_idx]
            
    return map_to_xyz

def parse_pairs_text(pairs_text: Optional[str]) -> List[Tuple[int, int]]:
    if not pairs_text:
        return []
    pairs: List[Tuple[int, int]] = []
    for chunk in re.split(r"[;,]", str(pairs_text)):
        piece = chunk.strip()
        if not piece:
            continue
        match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", piece)
        if not match:
            continue
        a, b = int(match.group(1)), int(match.group(2))
        if a != b:
            pairs.append((a, b))
    return _deduplicate_pairs(pairs)


def extract_formed_pairs_from_reaction_smarts(reaction_smarts: Optional[str]) -> List[Tuple[int, int]]:
    if not reaction_smarts:
        return []

    try:
        rxn = rdChemReactions.ReactionFromSmarts(reaction_smarts)
    except Exception:
        return []
    if rxn is None:
        return []

    reactant_bonds = _collect_mapped_bonds(rxn.GetReactants())
    product_bonds = _collect_mapped_bonds(rxn.GetProducts())
    formed = product_bonds.difference(reactant_bonds)
    return sorted(formed)


def extract_reactant_smiles(reaction_smiles: Optional[str]) -> Optional[str]:
    if not reaction_smiles:
        return None

    if ">>" not in reaction_smiles:
        return reaction_smiles

    left = reaction_smiles.split(">>", 1)[0].strip()
    return left or None


def extract_product_smiles(reaction_smiles: Optional[str]) -> Optional[str]:
    if not reaction_smiles:
        return None

    if ">>" not in reaction_smiles:
        return reaction_smiles

    right = reaction_smiles.split(">>", 1)[1].strip()
    return right or None


def match_reaction_profile_key(
    reaction_type: Optional[str],
    reaction_profiles: Dict[str, Dict[str, Any]],
) -> Optional[str]:
    if not reaction_profiles:
        return None

    if not reaction_type:
        return "_universal" if "_universal" in reaction_profiles else None

    rt = reaction_type.strip()
    candidates: List[str] = []

    def _add(candidate: str) -> None:
        if candidate and candidate not in candidates:
            candidates.append(candidate)

    normalized = rt.replace(" ", "")
    _add(rt)
    _add(f"{rt}_default")
    _add(normalized)
    _add(f"{normalized}_default")

    if normalized.startswith("[") and normalized.endswith("]") and len(normalized) > 2:
        inner = normalized[1:-1]
        _add(inner)
        _add(f"{inner}_default")
    else:
        bracketed = f"[{normalized}]"
        _add(bracketed)
        _add(f"{bracketed}_default")

    bracket_match = re.match(r"\s*(\[[^\]]+\])", rt)
    if bracket_match:
        bracket_part = bracket_match.group(1)
        bracket_normalized = bracket_part.replace(" ", "")
        _add(bracket_normalized)
        _add(f"{bracket_normalized}_default")
        inner = bracket_normalized[1:-1] if bracket_normalized.startswith("[") else bracket_normalized
        _add(f"[{inner}]_default")

    for key in candidates:
        if key in reaction_profiles:
            return key

    return "_universal" if "_universal" in reaction_profiles else None


def _iter_cleaner_rows(path: Path, delimiter: str) -> Iterable[Dict[str, Any]]:
    suffix = path.suffix.lower()
    if suffix == ".json":
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, list):
            for item in data:
                if isinstance(item, dict):
                    yield item
            return
        if isinstance(data, dict):
            if isinstance(data.get("records"), list):
                for item in data["records"]:
                    if isinstance(item, dict):
                        yield item
                return
            raise CleanerAdapterError("JSON cleaner data must be a list or contain a 'records' list")
        raise CleanerAdapterError("Unsupported JSON cleaner format")

    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            yield row


def _collect_mapped_bonds(mols: Sequence[Chem.Mol]) -> set[Tuple[int, int]]:
    bonds: set[Tuple[int, int]] = set()
    for mol in mols:
        for bond in mol.GetBonds():
            a = bond.GetBeginAtom().GetAtomMapNum()
            b = bond.GetEndAtom().GetAtomMapNum()
            if a > 0 and b > 0 and a != b:
                bonds.add((a, b) if a < b else (b, a))
    return bonds


def _deduplicate_pairs(pairs: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    seen: set[Tuple[int, int]] = set()
    out: List[Tuple[int, int]] = []
    for a, b in pairs:
        pair = (a, b) if a < b else (b, a)
        if pair not in seen:
            seen.add(pair)
            out.append(pair)
    return out


def _first_nonempty(row: Dict[str, Any], *keys: str) -> Optional[str]:
    for key in keys:
        value = row.get(key)
        if value is None:
            continue
        text = str(value).strip()
        if text:
            return text
    return None


def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _pairs_to_str(pairs: Sequence[Tuple[int, int]]) -> str:
    if not pairs:
        return ""
    return ";".join(f"{a}-{b}" for a, b in pairs)
