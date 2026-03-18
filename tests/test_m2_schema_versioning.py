import json
import re

import pytest


@pytest.fixture
def mech_index_payload() -> dict[str, object]:
    return {
        "version": "1.0.0",
        "schema_version": "mech_index_v1",
        "generated_at": "2026-01-21T00:00:00.000Z",
        "mechanism_status": "COMPLETE",
        "quality_flags": {
            "atom_count_ok": True,
            "forming_bond_window_ok": None,
            "suspect_optimized_to_product": "unknown",
        },
    }


def test_schema_version_present(mech_index_payload: dict[str, object]) -> None:
    assert mech_index_payload["schema_version"] == "mech_index_v1"


def test_generated_at_is_iso8601_utc(mech_index_payload: dict[str, object]) -> None:
    generated_at = mech_index_payload["generated_at"]
    assert isinstance(generated_at, str)
    assert re.match(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z$", generated_at)


def test_quality_flags_three_state_serializable(mech_index_payload: dict[str, object], tmp_path) -> None:
    out = tmp_path / "mech_index.json"
    out.write_text(json.dumps(mech_index_payload, indent=2), encoding="utf-8")
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["quality_flags"]["atom_count_ok"] is True
    assert loaded["quality_flags"]["forming_bond_window_ok"] is None
