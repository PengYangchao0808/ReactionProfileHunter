import re

import pytest

from rph_core.utils.reaction_identity import (
    generate_condition_id,
    generate_reaction_id,
    normalize_record_identity,
    sanitize_row_id,
)


def test_sanitize_row_id_basic() -> None:
    assert sanitize_row_id("row 001") == "row_001"
    assert sanitize_row_id("  ") == "unknown"


def test_generate_condition_id() -> None:
    assert generate_condition_id("row 001") == "COND_row_001"


def test_generate_reaction_id_format_and_stability() -> None:
    rid1 = generate_reaction_id(
        reactant_smiles_canon="C=C",
        product_smiles_canon="C=CC",
        reaction_type="[4+3]",
        cyclo_mode="concerted",
    )
    rid2 = generate_reaction_id(
        reactant_smiles_canon="C=C",
        product_smiles_canon="C=CC",
        reaction_type="[4+3]",
        cyclo_mode="concerted",
    )
    assert rid1 == rid2
    assert re.fullmatch(r"RXN_[0-9a-f]{8}", rid1)


def test_generate_reaction_id_changes_with_mode() -> None:
    rid1 = generate_reaction_id(
        reactant_smiles_canon="C=C",
        product_smiles_canon="C=CC",
        reaction_type="[4+3]",
        cyclo_mode="concerted",
    )
    rid2 = generate_reaction_id(
        reactant_smiles_canon="C=C",
        product_smiles_canon="C=CC",
        reaction_type="[4+3]",
        cyclo_mode="stepwise",
    )
    assert rid1 != rid2


def test_normalize_record_identity_smiles_canonicalization() -> None:
    ident = normalize_record_identity(
        {
            "rx_id": "row_0012",
            "precursor_smiles": "C=C",
            "product_smiles_main": "C=CC",
            "reaction_type": "[4+3]",
            "cyclo_mode": "concerted",
        }
    )
    assert ident.row_id == "row_0012"
    assert ident.condition_id == "COND_row_0012"
    assert re.fullmatch(r"RXN_[0-9a-f]{8}", ident.reaction_id)


def test_normalize_record_identity_missing_smiles_raises() -> None:
    with pytest.raises(ValueError):
        normalize_record_identity({"rx_id": "row_1", "product_smiles_main": "C"})
