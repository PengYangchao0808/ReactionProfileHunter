import json
from pathlib import Path

from rph_core.steps.conformer_search.censo_lite import S1_MANIFEST_SCHEMA_VERSION
from rph_core.steps.conformer_search.censo_lite_runtime import CensoLiteRuntime
from rph_core.utils.artifact_reconciler import ArtifactReconciler


def test_legacy_s1_is_migrated_without_qc(tmp_path: Path):
    smiles = "[CH2:1]1[CH2:2][CH2:3][CH2:4]1"
    runtime = CensoLiteRuntime({}, tmp_path / "S1_ConfSearch", "product_major")
    selected = runtime.embed(smiles)
    selected.rename(runtime.molecule_dir / "selected.xyz")
    selected = runtime.molecule_dir / "selected.xyz"
    (runtime.molecule_dir / "initial_atom_mapping.json").unlink()
    manifest = runtime.molecule_dir / "manifest.json"
    manifest.write_text(
        json.dumps({
            "schema_version": "s1_censo_light_ranking_v3",
            "selected": "product_major_conf_0001",
            "selected_xyz": "selected.xyz",
            "candidates": [{"id": "product_major_conf_0001", "xyz": "selected.xyz"}],
        }),
        encoding="utf-8",
    )
    reconciler = ArtifactReconciler(tmp_path, S1_MANIFEST_SCHEMA_VERSION)

    assessment = reconciler.assess_s1(
        manifest,
        expected_smiles=smiles,
        recorded_smiles=smiles,
        changed_paths=("manifest_schema", "protocol_version"),
    )
    migrated = reconciler.migrate_s1(assessment)

    assert assessment.status == "reusable_after_migration"
    assert migrated.status == "reusable_validated"
    assert migrated.evidence["qc_recomputed"] is False
    assert (runtime.molecule_dir / "atom_mapping.json").is_file()
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload["schema_version"] == S1_MANIFEST_SCHEMA_VERSION
    assert payload["artifact_reconciliation"]["qc_recomputed"] is False
    assert list((runtime.molecule_dir / "migrations").glob("manifest_*.json"))


def test_scientific_s1_change_is_not_migrated(tmp_path: Path):
    manifest = tmp_path / "manifest.json"
    selected = tmp_path / "selected.xyz"
    selected.write_text("1\ntest\nH 0 0 0\n", encoding="utf-8")
    manifest.write_text(
        json.dumps({"schema_version": "legacy", "selected_xyz": "selected.xyz"}),
        encoding="utf-8",
    )
    assessment = ArtifactReconciler(tmp_path, S1_MANIFEST_SCHEMA_VERSION).assess_s1(
        manifest,
        expected_smiles="[H][H]",
        recorded_smiles="[H][H]",
        changed_paths=("config.ranking.method",),
    )

    assert assessment.status == "dependency_changed"
