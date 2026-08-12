from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

import rph_core.steps.refinement.engine as refinement_engine_module
from rph_core.steps.refinement.manifest_io import REFINEMENT_MANIFEST_V1
from rph_core.steps.step3_lowlevel import LowLevelEngine
from rph_core.utils.config_loader import load_config
from rph_core.utils.provenance import (
    Provenance,
    build_provenance,
    canonicalize_forming_bonds_for_hash,
    hash_atom_mapping,
    sha256_of_file,
    sha256_of_json,
    verify_provenance_chain,
)
from rph_core.utils.qc_models import QCJobResult


def _write_json(path: Path, payload: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def test_sha256_of_file_existing_and_missing(tmp_path: Path) -> None:
    existing = tmp_path / "sample.txt"
    existing.write_text("abc", encoding="utf-8")

    assert sha256_of_file(existing) == hashlib.sha256(b"abc").hexdigest()
    assert sha256_of_file(tmp_path / "missing.txt") == ""


def test_sha256_of_json_is_deterministic() -> None:
    left = {"b": 2, "a": {"d": 4, "c": 3}}
    right = {"a": {"c": 3, "d": 4}, "b": 2}

    assert sha256_of_json(left) == sha256_of_json(right)


def test_canonicalize_forming_bonds_for_hash() -> None:
    canonical = canonicalize_forming_bonds_for_hash([[3, 2], [1, 4]])

    assert canonical == json.dumps([[1, 4], [2, 3]])
    assert hashlib.sha256(canonical.encode("utf-8")).hexdigest() == hashlib.sha256(
        json.dumps([[1, 4], [2, 3]]).encode("utf-8")
    ).hexdigest()


def test_provenance_to_dict_roundtrip() -> None:
    provenance = Provenance(
        run_id="run-123",
        schema_version=REFINEMENT_MANIFEST_V1,
        protocol_version=None,
        parent_stage_manifest_hashes={"s2": "abc"},
        atom_mapping_sha256="def",
        forming_bonds_hash="ghi",
        variant_manifest_hash="jkl",
        extra={"source": "test"},
    )

    assert Provenance(**provenance.to_dict()) == provenance


def test_build_provenance_with_full_inputs(tmp_path: Path) -> None:
    s0_manifest = _write_json(tmp_path / "S0_Mechanism" / "mechanism.json", {"stage": "S0"})
    s1_manifest = _write_json(tmp_path / "S1_ConfSearch" / "product" / "manifest.json", {"stage": "S1"})
    variant_manifest = _write_json(tmp_path / "variant.json", {"variant": "product_major"})
    atom_mapping_payload = {
        "mapped_product_smiles": "[CH2:1][CH2:2]",
        "map_to_product_smiles": {"1": 0, "2": 1},
        "product_smiles_idx_space": "geometry_product_smiles_idx",
    }

    provenance = build_provenance(
        run_id="11111111-1111-4111-8111-111111111111",
schema_version="s2_peb_manifest_v9",
        protocol_version=None,
        parent_manifest_paths={"s0": s0_manifest, "s1": s1_manifest},
        atom_mapping_payload=atom_mapping_payload,
        forming_bonds=[[3, 2], [1, 4]],
        variant_manifest_path=variant_manifest,
        extra={"selected_xyz_hash": "xyz"},
    )

    assert provenance.parent_stage_manifest_hashes == {
        "s0": sha256_of_file(s0_manifest),
        "s1": sha256_of_file(s1_manifest),
    }
    assert provenance.atom_mapping_sha256 == hash_atom_mapping(atom_mapping_payload)
    assert provenance.forming_bonds_hash == hashlib.sha256(
        canonicalize_forming_bonds_for_hash([[3, 2], [1, 4]]).encode("utf-8")
    ).hexdigest()
    assert provenance.variant_manifest_hash == sha256_of_file(variant_manifest)
    assert provenance.extra == {"selected_xyz_hash": "xyz"}


def test_build_provenance_with_none_inputs() -> None:
    provenance = build_provenance(
        run_id=None,
        schema_version="rph_run_manifest_v1",
        protocol_version=None,
        parent_manifest_paths={},
        atom_mapping_payload=None,
        forming_bonds=None,
        variant_manifest_path=None,
        extra=None,
    )

    assert provenance.atom_mapping_sha256 is None
    assert provenance.forming_bonds_hash is None
    assert provenance.variant_manifest_hash is None
    assert provenance.parent_stage_manifest_hashes == {}


def test_build_provenance_tolerates_missing_parent_file(tmp_path: Path) -> None:
    provenance = build_provenance(
        run_id=None,
        schema_version="pipeline_result",
        protocol_version=None,
        parent_manifest_paths={"s0": tmp_path / "missing.json"},
    )

    assert provenance.parent_stage_manifest_hashes["s0"] == ""


def test_s0_manifest_includes_provenance_field(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    from rph_core.steps.mechanism_classifier.s0_artifact_writer import write_s0_artifacts
    from rph_core.steps.mechanism_classifier.s0_record import S0ReactionRecord

    source_csv = tmp_path / "dataset.csv"
    source_csv.write_text("rx_id\nRXN_PROVENANCE\n", encoding="utf-8")
    record = S0ReactionRecord(
        rx_id="RXN_PROVENANCE",
        product_smiles="C1CCC1",
        reaction_type="4+3",
        mapped_product_smiles="[CH2:1]1[CH2:2][CH2:3][CH2:4]1",
        mapped_forming_bonds=((1, 2), (3, 4)),
        forming_bonds=((0, 1), (2, 3)),
        mapping_confidence=0.99,
        mapping_trusted=True,
        source_csv=source_csv,
        canonical_precursor_smiles="CC",
        mapped_precursor_smiles="[CH3:1][CH3:2]",
        topology="INTER",
        cyclo_mode="UNKNOWN",
        source_row_hash="row-hash",
    )

    _context, _variants, mechanism_path = write_s0_artifacts(
        tmp_path / "run",
        record,
        run_id="22222222-2222-4222-8222-222222222222",
    )
    mechanism = json.loads(mechanism_path.read_text(encoding="utf-8"))

    assert mechanism["provenance"]["schema_version"] == "s0_mechanism_v3"
    assert mechanism["provenance"]["parent_stage_manifest_hashes"] == {}
    assert mechanism["provenance"]["atom_mapping_sha256"]
    assert mechanism["provenance"]["forming_bonds_hash"]
    assert mechanism["provenance"]["extra"]["source_row_hash"] == "row-hash"


def test_s3_manifest_includes_provenance_field_with_parent_s2_hash(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    xyz = tmp_path / "seed.xyz"
    xyz.write_text("1\nseed\nH 0.0 0.0 0.0\n", encoding="utf-8")
    s2_manifest = _write_json(tmp_path / "S2_PEB" / "manifest.json", {"stage": "S2"})

    def fake_opt(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_xyz=out_xyz,
            output_file=output_dir / "opt.out",
            energy_hartree=-1.0,
        )

    def fake_frequency(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / "freq.out"
        out_file.write_text("mock freq\n", encoding="utf-8")
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=out_file,
            energy_hartree=-0.9,
            frequencies_cm1=(25.0, 125.0, 325.0),
        )

    monkeypatch.setattr(refinement_engine_module, "run_optimization", fake_opt)
    monkeypatch.setattr(refinement_engine_module, "run_frequency", fake_frequency)

    config = load_config()
    manifest_path = LowLevelEngine(
        config,
        run_id="33333333-3333-4333-8333-333333333333",
        parent_manifest_paths={"s2": s2_manifest},
    ).run(
        [{"id": "product", "role": "product", "kind": "minimum", "input_xyz": str(xyz)}],
        tmp_path / "S3_LowLevel",
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert manifest["provenance"]["parent_stage_manifest_hashes"]["s2"] == sha256_of_file(s2_manifest)
    assert manifest["provenance"]["variant_manifest_hash"] == sha256_of_file(s2_manifest)


def test_verify_provenance_chain_detects_tampering(tmp_path: Path) -> None:
    parent_manifest = _write_json(tmp_path / "parent.json", {"value": 1})
    child_manifest = tmp_path / "child.json"
    provenance = build_provenance(
        run_id="44444444-4444-4444-8444-444444444444",
        schema_version="child_v1",
        protocol_version=None,
        parent_manifest_paths={"s0": parent_manifest},
    )
    _write_json(
        child_manifest,
        {"schema_version": "child_v1", "provenance": provenance.to_dict()},
    )

    assert verify_provenance_chain(child_manifest, expected_parent_paths={"s0": parent_manifest}) == (True, [])

    parent_manifest.write_text(json.dumps({"value": 2}, indent=2), encoding="utf-8")
    ok, mismatches = verify_provenance_chain(
        child_manifest,
        expected_parent_paths={"s0": parent_manifest},
    )

    assert not ok
    assert mismatches
    assert "s0" in mismatches[0]
