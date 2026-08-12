from __future__ import annotations

import json
from pathlib import Path

from rph_core.steps.fidelity_profile import FidelityProfile
from rph_core.steps.refinement import RefinementEngine
from rph_core.steps.refinement.engine import RefinementEngine as ModuleRefinementEngine
from rph_core.steps.refinement.manifest_io import (
    LEGACY_S3_SCHEMA,
    LEGACY_S4_SCHEMA,
    REFINEMENT_MANIFEST_V1,
    read_refinement_manifest,
    write_refinement_manifest,
)
from rph_core.utils.config_loader import load_config
from rph_core.utils.provenance import sha256_of_file
from rph_core.utils.stage_progress import StageProgressReporter


def _load_config() -> dict:
    return load_config()


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_refinement_engine_writes_s3_manifest(tmp_path: Path):
    config = _load_config()
    engine = RefinementEngine(config, FidelityProfile.from_config(config, "S3"))

    manifest_path = engine.run([], tmp_path)
    manifest = _read_json(manifest_path)

    assert manifest_path == tmp_path / "manifest.json"
    assert manifest["schema_version"] == REFINEMENT_MANIFEST_V1
    assert manifest["stage"] == "S3"
    assert manifest["fidelity"] == "low"
    assert manifest["profile_id"] == "b97_3c_r2scan_3c_v1"
    assert manifest["structures"] == []


def test_refinement_engine_writes_s4_manifest(tmp_path: Path):
    config = _load_config()
    engine = RefinementEngine(config, FidelityProfile.from_config(config, "S4"))

    manifest_path = engine.run([], tmp_path)
    manifest = _read_json(manifest_path)

    assert manifest_path == tmp_path / "manifest.json"
    assert manifest["schema_version"] == REFINEMENT_MANIFEST_V1
    assert manifest["stage"] == "S4"
    assert manifest["fidelity"] == "high"
    assert manifest["profile_id"] == "m062x_wb97mv_v1"
    assert manifest["structures"] == []


def test_refinement_engine_manifest_has_provenance(tmp_path: Path):
    config = _load_config()
    engine = RefinementEngine(config, FidelityProfile.from_config(config, "S3"))

    manifest = _read_json(engine.run([], tmp_path))
    provenance = manifest["provenance"]

    assert isinstance(provenance, dict)
    assert provenance["schema_version"] == REFINEMENT_MANIFEST_V1
    assert "run_id" in provenance
    assert isinstance(provenance["parent_stage_manifest_hashes"], dict)


def test_refinement_engine_accepts_parent_manifest_paths(tmp_path: Path):
    config = _load_config()
    parent_manifest = tmp_path / "s2.json"
    parent_manifest.write_text('{"schema_version": "s2_v1"}', encoding="utf-8")
    engine = RefinementEngine(
        config,
        FidelityProfile.from_config(config, "S3"),
        parent_manifest_paths={"s2": parent_manifest},
    )

    manifest = _read_json(engine.run([], tmp_path / "s3"))
    hashes = manifest["provenance"]["parent_stage_manifest_hashes"]

    assert hashes["s2"] == sha256_of_file(parent_manifest)


def test_refinement_engine_archive_called(tmp_path: Path):
    config = _load_config()
    output_dir = tmp_path / "stage"
    output_dir.mkdir()
    stale_file = output_dir / "stale.txt"
    stale_file.write_text("stale", encoding="utf-8")
    engine = RefinementEngine(config, FidelityProfile.from_config(config, "S3"))

    manifest = _read_json(engine.run([], output_dir))

    archived_to = manifest["stale_outputs_archived_to"]
    assert archived_to is not None
    assert Path(archived_to).is_dir()
    assert not stale_file.exists()
    assert (Path(archived_to) / "stale.txt").is_file()


def test_refinement_engine_restores_live_status_after_archiving(tmp_path: Path):
    config = _load_config()
    output_dir = tmp_path / "stage"
    output_dir.mkdir()
    (output_dir / "stale.txt").write_text("stale", encoding="utf-8")
    reporter = StageProgressReporter(output_dir, "S3")
    engine = RefinementEngine(config, FidelityProfile.from_config(config, "S3"))
    engine.set_progress_reporter(reporter)

    engine.run([], output_dir)

    assert (output_dir / "status.json").is_file()
    assert "stage_progress_reinitialized" in (output_dir / "events.jsonl").read_text(encoding="utf-8")


def test_refinement_engine_accepts_metadata_paths_in_mapping_request(tmp_path: Path):
    config = _load_config()
    engine = RefinementEngine(config, FidelityProfile.from_config(config, "S3"))
    input_xyz = tmp_path / "input.xyz"
    input_xyz.write_text("1\nH\nH 0.0 0.0 0.0\n", encoding="utf-8")
    s1_manifest = tmp_path / "s1" / "manifest.json"
    mapping_audit = tmp_path / "s2" / "atom_mapping.json"

    request = engine._coerce_to_structure_request(
        {
            "id": "product_major",
            "role": "product",
            "kind": "minimum",
            "input_xyz": str(input_xyz),
            "s1_manifest": str(s1_manifest),
            "mapping_audit": str(mapping_audit),
        }
    )

    assert request.s1_manifest == str(s1_manifest)
    assert request.mapping_audit == str(mapping_audit)


def test_read_refinement_manifest_v1(tmp_path: Path):
    manifest_path = write_refinement_manifest(
        tmp_path / "manifest.json",
        stage="S3",
        fidelity="low",
        profile_id="b97_3c_r2scan_3c_v1",
        structures=[],
        run_id="run-123",
        extra={"note": "ok"},
    )

    manifest = read_refinement_manifest(manifest_path)

    assert manifest["schema_version"] == REFINEMENT_MANIFEST_V1
    assert manifest["stage"] == "S3"
    assert manifest["fidelity"] == "low"
    assert manifest["profile_id"] == "b97_3c_r2scan_3c_v1"
    assert manifest["structures"] == []
    assert manifest["run_id"] == "run-123"
    assert manifest["note"] == "ok"


def test_read_refinement_manifest_legacy_s3(tmp_path: Path):
    manifest_path = tmp_path / "legacy_s3.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema_version": LEGACY_S3_SCHEMA,
                "structures": {
                    "a": {"id": "a", "status": "complete"},
                    "b": {"id": "b", "status": "failed"},
                },
            }
        ),
        encoding="utf-8",
    )

    manifest = read_refinement_manifest(manifest_path)

    assert manifest["schema_version"] == LEGACY_S3_SCHEMA
    assert manifest["stage"] == "S3"
    assert manifest["fidelity"] == "low"
    assert manifest["profile_id"] == "b97_3c_r2scan_3c_v1_legacy"
    assert manifest["structures"] == [
        {"id": "a", "status": "complete"},
        {"id": "b", "status": "failed"},
    ]


def test_read_refinement_manifest_legacy_s4(tmp_path: Path):
    manifest_path = tmp_path / "legacy_s4.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema_version": LEGACY_S4_SCHEMA,
                "structures": [{"id": "ts", "status": "complete"}],
            }
        ),
        encoding="utf-8",
    )

    manifest = read_refinement_manifest(manifest_path)

    assert manifest["schema_version"] == LEGACY_S4_SCHEMA
    assert manifest["stage"] == "S4"
    assert manifest["fidelity"] == "high"
    assert manifest["profile_id"] == "m062x_wb97mv_v1_legacy"
    assert manifest["structures"] == [{"id": "ts", "status": "complete"}]


def test_read_refinement_manifest_unknown_schema_returns_raw(tmp_path: Path):
    manifest_path = tmp_path / "future.json"
    payload = {
        "schema_version": "future_v99",
        "stage": "S9",
        "structures": ["raw"],
    }
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")

    manifest = read_refinement_manifest(manifest_path)

    assert manifest == payload


def test_refinement_config_section_present():
    config = _load_config()

    assert ModuleRefinementEngine is RefinementEngine
    assert config["refinement"]["common"]["workflow"]["warmup"]["enabled_roles"] == [
        "intermediate",
        "ts",
    ]
    assert config["refinement"]["s3"]["profile_id"] == "b97_3c_r2scan_3c_v1"
    assert config["refinement"]["s4"]["profile_id"] == "m062x_wb97mv_v1"


def test_refinement_config_does_not_break_theory_section():
    config = _load_config()

    assert config["theory"]["s3_low_level"]["optimization"]["method"] == "B97-3c"
