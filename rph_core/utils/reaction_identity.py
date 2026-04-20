from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from collections.abc import Mapping

from rph_core.utils.molecule_utils import canonicalize_smiles


_ID_SANITIZE_RE = re.compile(r"[^A-Za-z0-9_-]+")


def sanitize_row_id(row_id: str) -> str:
    cleaned = _ID_SANITIZE_RE.sub("_", str(row_id).strip())
    return cleaned or "unknown"


def generate_condition_id(row_id: str) -> str:
    return f"COND_{sanitize_row_id(row_id)}"


def generate_reaction_id(
    *,
    reactant_smiles_canon: str,
    product_smiles_canon: str,
    reaction_type: str,
    cyclo_mode: str,
) -> str:
    payload = (
        f"{reactant_smiles_canon}>>{product_smiles_canon}"
        f"|{reaction_type.strip()}"
        f"|{cyclo_mode.strip()}"
    )
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:8]
    return f"RXN_{digest}"


def _first_str(data: Mapping[str, object], *keys: str) -> str | None:
    for key in keys:
        value = data.get(key)
        if value is None:
            continue
        text = str(value).strip()
        if text:
            return text
    return None


def _canonical_or_none(smiles: str | None) -> str | None:
    if smiles is None:
        return None
    text = str(smiles).strip()
    if not text:
        return None
    return canonicalize_smiles(text)


@dataclass(frozen=True)
class RecordIdentity:
    row_id: str
    reaction_id: str
    condition_id: str
    reactant_smiles_canon: str
    product_smiles_canon: str
    reaction_type: str
    cyclo_mode: str


def normalize_record_identity(record: Mapping[str, object]) -> RecordIdentity:
    row_id = _first_str(record, "row_id", "rx_id", "id")
    if not row_id:
        raise ValueError("Missing row_id/rx_id in dataset record")

    product_smiles = _first_str(
        record,
        "product_smiles_canon",
        "product_smiles_main",
        "product_smiles",
        "mapped_product_smiles",
    )
    product_canon = _canonical_or_none(product_smiles)
    if not product_canon:
        raise ValueError(f"Invalid or missing product SMILES for row_id={row_id}")

    reactant_or_precursor_smiles = _first_str(
        record,
        "reactant_smiles_canon",
        "reactant_smiles",
        "precursor_smiles",
        "mapped_reactant_smiles",
    )
    reactant_canon = _canonical_or_none(reactant_or_precursor_smiles)
    if not reactant_canon:
        raise ValueError(f"Invalid or missing reactant/precursor SMILES for row_id={row_id}")

    reaction_type = _first_str(record, "reaction_type", "rxn_type", "reaction_family") or ""
    cyclo_mode = _first_str(record, "cyclo_mode", "mode") or "concerted"

    reaction_id = generate_reaction_id(
        reactant_smiles_canon=reactant_canon,
        product_smiles_canon=product_canon,
        reaction_type=reaction_type,
        cyclo_mode=cyclo_mode,
    )
    condition_id = generate_condition_id(row_id)

    return RecordIdentity(
        row_id=str(row_id),
        reaction_id=reaction_id,
        condition_id=condition_id,
        reactant_smiles_canon=reactant_canon,
        product_smiles_canon=product_canon,
        reaction_type=reaction_type,
        cyclo_mode=cyclo_mode,
    )
