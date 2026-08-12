"""Resolve trusted reaction records into the V4 atom-mapping contract.

The dataset stores mechanistic bond changes in atom-map space.  V4 keeps the
mapped product SMILES as the geometry reference so RDKit, S1 XYZ files and S2
all share one atom order.  Canonical/unmapped SMILES are display and identity
fields only; their atom indices must never drive a QC constraint.
"""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Sequence, Tuple

from rdkit import Chem


BondPair = Tuple[int, int]


@dataclass(frozen=True)
class S0ReactionRecord:
    """Trusted, product-side mechanism data ready for V4 S0."""

    rx_id: str
    product_smiles: str
    reaction_type: str
    mapped_product_smiles: str
    mapped_forming_bonds: Tuple[BondPair, ...]
    forming_bonds: Tuple[BondPair, ...]
    mapping_confidence: Optional[float]
    mapping_trusted: bool
    source_csv: Path
    canonical_precursor_smiles: str = ""
    mapped_precursor_smiles: str = ""
    topology: str = "INTER"
    cyclo_mode: str = "UNKNOWN"
    precursor_type: str = ""
    source_row_hash: str = ""
    raw_row: Dict[str, str] = field(default_factory=dict)

    def signature_payload(self) -> Dict[str, object]:
        """Return the fields that define the S0 result for checkpointing."""

        return {
            "rx_id": self.rx_id,
            "product_smiles": self.product_smiles,
            "reaction_type": self.reaction_type,
            "canonical_precursor_smiles": self.canonical_precursor_smiles,
            "mapped_product_smiles": self.mapped_product_smiles,
            "mapped_precursor_smiles": self.mapped_precursor_smiles,
            "mapped_forming_bonds": self.mapped_forming_bonds,
            "forming_bonds": self.forming_bonds,
            "topology": self.topology,
            "cyclo_mode": self.cyclo_mode,
            "precursor_type": self.precursor_type,
            "source_row_hash": self.source_row_hash,
            "mapping_confidence": self.mapping_confidence,
            "mapping_trusted": self.mapping_trusted,
        }


def load_s0_reaction_record(csv_path: Path, rx_id: str) -> S0ReactionRecord:
    """Load one trusted dataset reaction and convert bonds to product indices.

    ``mechanistic_forming_bonds`` is authoritative for S0.  It is recorded
    in atom-map-number space and must never be treated as XYZ indices.
    """

    path = Path(csv_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Reaction CSV does not exist: {path}")

    row = _find_row(path, rx_id)
    if row is None:
        raise ValueError(f"rx_id={rx_id!r} was not found in {path}")

    product_smiles = _required(row, "product_smiles_main")
    mapped_reaction = _required(row, "rxn_smiles_mapped")
    mapped_precursor_smiles, mapped_product_smiles = _split_mapped_reaction(mapped_reaction)
    reaction_type = _required(row, "reaction_type")
    mapping_trusted = _as_bool(row.get("mapping_trusted_for_mechanistic"))
    if not mapping_trusted:
        raise ValueError(
            f"rx_id={rx_id!r} has untrusted mechanistic mapping; "
            "refuse to use it as the S0 forming-bond source"
        )

    mapped_forming_bonds = _bond_pairs(row.get("mechanistic_forming_bonds"), rx_id)
    if len(mapped_forming_bonds) != 2:
        raise ValueError(
            f"rx_id={rx_id!r} must provide exactly two mechanistic_forming_bonds; "
            f"got {len(mapped_forming_bonds)}"
        )
    _validate_mapped_reaction_bonds(
        mapped_precursor_smiles,
        mapped_product_smiles,
        mapped_forming_bonds,
    )
    forming_bonds = _map_bonds_to_product_indices(
        product_smiles,
        mapped_product_smiles,
        mapped_forming_bonds,
    )
    source_row_hash = hashlib.sha256(
        json.dumps(row, sort_keys=True).encode("utf-8")
    ).hexdigest()
    return S0ReactionRecord(
        rx_id=str(row["rx_id"]),
        product_smiles=product_smiles,
        reaction_type=reaction_type,
        mapped_product_smiles=mapped_product_smiles,
        mapped_forming_bonds=mapped_forming_bonds,
        forming_bonds=forming_bonds,
        mapping_confidence=_optional_float(row.get("mapping_confidence")),
        mapping_trusted=mapping_trusted,
        source_csv=path,
        canonical_precursor_smiles=_optional_str(row, "precursor_smiles"),
        mapped_precursor_smiles=mapped_precursor_smiles,
        topology=_optional_str(row, "topology", default="INTER"),
        cyclo_mode=_optional_str(row, "cyclo_mode", default="UNKNOWN"),
        precursor_type=_optional_str(row, "precursor_type"),
        source_row_hash=source_row_hash,
        raw_row=dict(row),
    )


def _find_row(path: Path, rx_id: str) -> Optional[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if str(row.get("rx_id", "")).strip() == str(rx_id).strip():
                return row
    return None


def _required(row: Dict[str, str], field: str) -> str:
    value = str(row.get(field, "")).strip()
    if not value:
        raise ValueError(f"Reaction record has no {field}")
    return value


def _optional_str(row: Dict[str, str], field: str, default: str = "") -> str:
    value = str(row.get(field, "")).strip()
    return value or default


def _mapped_product(mapped_reaction: str) -> str:
    return _split_mapped_reaction(mapped_reaction)[1]


def _split_mapped_reaction(mapped_reaction: str) -> Tuple[str, str]:
    try:
        reactants, product = mapped_reaction.split(">>", 1)
    except ValueError as exc:
        raise ValueError("rxn_smiles_mapped must contain a reactants>>products separator") from exc
    if not product:
        raise ValueError("rxn_smiles_mapped has no mapped product")
    return reactants.strip(), product.strip()


def _bond_pairs(value: Optional[str], rx_id: str) -> Tuple[BondPair, ...]:
    try:
        raw_pairs: Any = json.loads(value or "[]")
    except json.JSONDecodeError as exc:
        raise ValueError(f"rx_id={rx_id!r} has invalid mechanistic_forming_bonds JSON") from exc
    if not isinstance(raw_pairs, list):
        raise ValueError(f"rx_id={rx_id!r} has invalid mechanistic_forming_bonds JSON")
    pairs = []
    for pair in raw_pairs:
        if not isinstance(pair, Sequence) or len(pair) != 2:
            raise ValueError(f"rx_id={rx_id!r} has invalid forming-bond pair: {pair!r}")
        raw_atom_a: Any = pair[0]
        raw_atom_b: Any = pair[1]
        atom_a, atom_b = int(raw_atom_a), int(raw_atom_b)
        if atom_a == atom_b:
            raise ValueError(f"rx_id={rx_id!r} has a self-bond: {pair!r}")
        pairs.append((atom_a, atom_b))
    return tuple(pairs)


def _map_bonds_to_product_indices(
    product_smiles: str,
    mapped_product_smiles: str,
    mapped_pairs: Iterable[BondPair],
) -> Tuple[BondPair, ...]:
    map_to_product = _build_map_to_product_smiles(product_smiles, mapped_product_smiles)
    product = Chem.MolFromSmiles(mapped_product_smiles)
    if product is None:
        raise ValueError("Could not parse product SMILES while resolving S0 atom indices")

    product_pairs = []
    for map_a, map_b in mapped_pairs:
        try:
            product_a = map_to_product[map_a]
            product_b = map_to_product[map_b]
        except KeyError as exc:
            raise ValueError(f"Forming-bond map number {exc.args[0]} is absent from the mapped product") from exc
        if product.GetBondBetweenAtoms(product_a, product_b) is None:
            raise ValueError(
                f"Mapped forming bond ({map_a}, {map_b}) is absent from mapped product"
            )
        product_pairs.append((product_a, product_b))
    return tuple(product_pairs)


def _build_map_to_product_smiles(
    product_smiles: str,
    mapped_product_smiles: str,
) -> Dict[int, int]:
    product = Chem.MolFromSmiles(product_smiles)
    mapped_product = Chem.MolFromSmiles(mapped_product_smiles)
    if product is None or mapped_product is None:
        raise ValueError("Could not parse product SMILES while resolving S0 atom indices")

    atom_match = product.GetSubstructMatch(mapped_product, useChirality=True)
    if not atom_match:
        atom_match = product.GetSubstructMatch(mapped_product, useChirality=False)
    if len(atom_match) != mapped_product.GetNumAtoms():
        raise ValueError("Mapped product does not map one-to-one onto product_smiles_main")

    # The substructure match above is validation only.  Executable indices are
    # the indices of the exact mapped SMILES handed to RDKit in S1.  Returning
    # atom_match[...] here was the V4 regression that changed rx_id=2 C7-C8
    # into C7-C22 after canonical SMILES reordering.
    map_to_product: Dict[int, int] = {}
    for mapped_atom in mapped_product.GetAtoms():
        atom_map = mapped_atom.GetAtomMapNum()
        # Mappers occasionally leave spectator atoms (for example, a carbonyl
        # oxygen) unnumbered.  They are irrelevant unless a declared forming
        # bond references them; that is checked explicitly below.
        if atom_map > 0:
            if atom_map in map_to_product:
                raise ValueError(f"Duplicate atom-map number {atom_map} in mapped product")
            map_to_product[atom_map] = mapped_atom.GetIdx()
    return map_to_product


def _map_number_to_atom(mol: Chem.Mol) -> Dict[int, int]:
    result: Dict[int, int] = {}
    for atom in mol.GetAtoms():
        map_number = atom.GetAtomMapNum()
        if map_number <= 0:
            continue
        if map_number in result:
            raise ValueError(f"Duplicate atom-map number {map_number} in mapped reaction")
        result[map_number] = atom.GetIdx()
    return result


def _validate_mapped_reaction_bonds(
    mapped_precursor_smiles: str,
    mapped_product_smiles: str,
    mapped_pairs: Iterable[BondPair],
) -> None:
    """Require every declared PEB bond to resolve and exist in the product."""

    precursor = Chem.MolFromSmiles(mapped_precursor_smiles)
    product = Chem.MolFromSmiles(mapped_product_smiles)
    if precursor is None or product is None:
        raise ValueError("Could not parse mapped reaction while validating forming bonds")
    _map_number_to_atom(precursor)
    product_map = _map_number_to_atom(product)
    for map_a, map_b in mapped_pairs:
        if map_a not in product_map or map_b not in product_map:
            raise ValueError(f"Forming bond ({map_a}, {map_b}) is unresolved in mapped product")
        if product.GetBondBetweenAtoms(product_map[map_a], product_map[map_b]) is None:
            raise ValueError(f"Forming bond ({map_a}, {map_b}) is absent from mapped product")


def _as_bool(value: Optional[str]) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y", "on"}


def _optional_float(value: Optional[str]) -> Optional[float]:
    text = str(value or "").strip()
    return float(text) if text else None
