from __future__ import annotations

import json
from pathlib import Path

from rph_core.steps.fidelity_profile import FidelityProfile
from rph_core.steps.refinement import RefinementEngine
from rph_core.steps.refinement.models import StructureRequest
from rph_core.utils.config_loader import load_config


def _engine(stage: str = "S3") -> RefinementEngine:
    config = load_config()
    return RefinementEngine(config, FidelityProfile.from_config(config, stage))


def _write_xyz(path: Path, distance: float = 1.5) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"2\nmock\nH 0.0 0.0 0.0\nH {distance:.6f} 0.0 0.0\n",
        encoding="utf-8",
    )
    return path


def test_pass0_valid_structure_succeeds(tmp_path: Path):
    engine = _engine()
    input_xyz = _write_xyz(tmp_path / "inputs" / "valid.xyz")
    request = StructureRequest(
        id="valid",
        role="precursor",
        kind="minimum",
        input_xyz=input_xyz,
        forming_bonds=[(0, 1)],
    )

    outcome = engine._run_pass0_preflight([request], tmp_path / "stage")[0]

    assert outcome.status == "ok"
    assert outcome.input_xyz == input_xyz


def test_pass0_missing_xyz_fails(tmp_path: Path):
    engine = _engine()
    request = StructureRequest(
        id="missing",
        role="precursor",
        kind="minimum",
        input_xyz=tmp_path / "inputs" / "missing.xyz",
    )

    outcome = engine._run_pass0_preflight([request], tmp_path / "stage")[0]

    assert outcome.status == "failed_preflight"
    assert "input_xyz not found" in str(outcome.error)


def test_pass0_fallback_xyz_used_when_primary_missing(tmp_path: Path):
    engine = _engine()
    fallback_xyz = _write_xyz(tmp_path / "inputs" / "fallback.xyz")
    request = StructureRequest(
        id="fallback",
        role="product",
        kind="minimum",
        input_xyz=tmp_path / "inputs" / "missing.xyz",
        fallback_xyz=fallback_xyz,
    )

    outcome = engine._run_pass0_preflight([request], tmp_path / "stage")[0]

    assert outcome.status == "ok"
    assert outcome.input_xyz == fallback_xyz


def test_pass0_forming_bond_out_of_range_fails(tmp_path: Path):
    engine = _engine()
    input_xyz = _write_xyz(tmp_path / "inputs" / "bad_bond.xyz")
    request = StructureRequest(
        id="bad-bond",
        role="intermediate",
        kind="minimum",
        input_xyz=input_xyz,
        forming_bonds=[(0, 2)],
    )

    outcome = engine._run_pass0_preflight([request], tmp_path / "stage")[0]

    assert outcome.status == "failed_preflight"
    assert "out of range" in str(outcome.error)


def test_pass0_writes_provenance_per_structure(tmp_path: Path):
    engine = _engine()
    input_xyz = _write_xyz(tmp_path / "inputs" / "provenance.xyz")
    request = StructureRequest(
        id="with-provenance",
        role="ts",
        kind="ts",
        input_xyz=input_xyz,
        forming_bonds=[(0, 1)],
    )

    engine._run_pass0_preflight([request], tmp_path / "stage")
    provenance_path = tmp_path / "stage" / "with-provenance" / "provenance.json"

    assert provenance_path.is_file()
    payload = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert payload["structure_id"] == "with-provenance"
    assert payload["input_xyz"] == str(input_xyz)


def test_pass0_creates_attempt_recorder_directory(tmp_path: Path):
    engine = _engine()
    input_xyz = _write_xyz(tmp_path / "inputs" / "attempts.xyz")
    request = StructureRequest(
        id="attempted",
        role="intermediate",
        kind="minimum",
        input_xyz=input_xyz,
        forming_bonds=[(0, 1)],
    )

    engine._run_pass0_preflight([request], tmp_path / "stage")

    assert (tmp_path / "stage" / "attempted" / "attempts").is_dir()
