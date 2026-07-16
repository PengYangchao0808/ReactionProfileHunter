import json
from pathlib import Path

from rph_core.v4_orchestrator import V4Orchestrator
from rph_core.utils.v4_checkpoint import V4Checkpoint


def test_s1_scientific_signature_ignores_queue_scheduling():
    baseline = {
        "crest": {"gfn_level": 2, "nproc": "auto"},
        "b97_3c": {
            "method": "B97-3c",
            "parallel_jobs": 16,
            "cores_per_job": 1,
            "queue_drain": {
                "enabled": True,
                "remaining_jobs": 8,
                "min_cores_per_job": 2,
                "max_cores_per_job": 16,
            },
            "scheduling": {
                "policy": "throughput_first",
                "final_wave": {"enabled": False, "profile_layouts": {}},
            },
        },
    }
    alternate_layout = {
        "crest": {"gfn_level": 2, "nproc": 8},
        "b97_3c": {
            "method": "B97-3c",
            "parallel_jobs": 8,
            "cores_per_job": 2,
            "queue_drain": {
                "enabled": False,
                "remaining_jobs": 4,
                "min_cores_per_job": 1,
                "max_cores_per_job": 4,
            },
            "scheduling": {
                "policy": "profiled_final_wave",
                "final_wave": {
                    "enabled": True,
                    "profile_layouts": {"4": [4, 4, 4, 4]},
                },
            },
        },
    }

    assert V4Orchestrator._scientific_s1_config(
        baseline
    ) == V4Orchestrator._scientific_s1_config(alternate_layout)


def test_v4_checkpoint_roundtrip(tmp_path: Path):
    checkpoint = V4Checkpoint(tmp_path)
    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}", encoding="utf-8")
    signature = checkpoint.signature({"stage": "s3"})
    checkpoint.mark("s3", signature, manifest)
    assert checkpoint.reusable("s3", signature, manifest)
    assert not checkpoint.reusable("s3", "different", manifest)


def test_v4_checkpoint_invalidates_downstream_stages(tmp_path: Path):
    checkpoint = V4Checkpoint(tmp_path)
    manifests = {}
    for stage in ("s0", "s1", "s2", "s3", "s4"):
        manifest = tmp_path / f"{stage}.json"
        manifest.write_text("{}", encoding="utf-8")
        manifests[stage] = manifest
        checkpoint.mark(stage, checkpoint.signature({"stage": stage}), manifest)

    replacement = tmp_path / "s1_replacement.json"
    replacement.write_text("{}", encoding="utf-8")
    signature = checkpoint.signature({"stage": "s1", "replacement": True})
    checkpoint.mark("s1", signature, replacement)

    assert checkpoint.reusable("s1", signature, replacement)
    assert not checkpoint.reusable("s2", checkpoint.signature({"stage": "s2"}), manifests["s2"])
    assert not checkpoint.reusable("s4", checkpoint.signature({"stage": "s4"}), manifests["s4"])


def test_record_s0_update_preserves_independent_s1_and_invalidates_s2(tmp_path: Path):
    checkpoint = V4Checkpoint(tmp_path)
    manifests = {}
    signatures = {}
    for stage in ("s0", "s1", "s2"):
        manifest = tmp_path / f"{stage}.json"
        manifest.write_text("{}", encoding="utf-8")
        manifests[stage] = manifest
        signatures[stage] = checkpoint.signature({"stage": stage})
        checkpoint.mark(stage, signatures[stage], manifest)

    updated_s0 = tmp_path / "s0_updated.json"
    updated_s0.write_text("{}", encoding="utf-8")
    updated_signature = checkpoint.signature({"stage": "s0", "record": "updated"})
    checkpoint.mark_s0(updated_signature, updated_s0)

    assert checkpoint.reusable("s0", updated_signature, updated_s0)
    assert checkpoint.reusable("s1", signatures["s1"], manifests["s1"])
    assert not checkpoint.reusable("s2", signatures["s2"], manifests["s2"])


def test_scoped_checkpoint_explains_signature_difference(tmp_path: Path):
    checkpoint = V4Checkpoint(tmp_path)
    manifest = tmp_path / "S1_ConfSearch" / "product_major" / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("{}", encoding="utf-8")
    recorded_payload = {
        "stage": "s1",
        "config": {"xtb_thermo": {"gfn_level": 1}},
    }
    current_payload = {
        "stage": "s1",
        "config": {"xtb_thermo": {"gfn_level": 2}},
    }
    checkpoint.mark_scope(
        "product_major",
        "s1",
        checkpoint.signature(recorded_payload),
        manifest,
        signature_payload=recorded_payload,
    )

    decision = checkpoint.check_reuse(
        "s1",
        checkpoint.signature(current_payload),
        manifest,
        scope="product_major",
        signature_payload=current_payload,
    )

    assert decision.reusable is False
    assert decision.reason == "signature_mismatch"
    assert decision.differences == ({
        "path": "config.xtb_thermo.gfn_level",
        "recorded": 1,
        "current": 2,
    },)


def test_start_from_s3_accepts_signature_mismatch_for_all_required_upstream_stages(
    tmp_path: Path,
):
    orchestrator = object.__new__(V4Orchestrator)
    orchestrator.resume_policy = "use-existing-upstream"
    orchestrator.start_from = "s3"
    orchestrator._reused_signature_overrides = {}

    for stage, scope in (("s0", None), ("s1", "product_major"), ("s2", "product_major")):
        checkpoint = V4Checkpoint(tmp_path / stage)
        manifest = tmp_path / stage / "manifest.json"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        manifest.write_text("{}", encoding="utf-8")
        recorded_payload = {"stage": stage, "config_revision": "recorded"}
        current_payload = {"stage": stage, "config_revision": "current"}
        recorded_signature = checkpoint.signature(recorded_payload)
        current_signature = checkpoint.signature(current_payload)
        if scope is None:
            checkpoint.mark(
                stage,
                recorded_signature,
                manifest,
                signature_payload=recorded_payload,
            )
        else:
            checkpoint.mark_scope(
                scope,
                stage,
                recorded_signature,
                manifest,
                signature_payload=recorded_payload,
            )

        assert orchestrator._checkpoint_reusable(
            checkpoint,
            stage,
            current_signature,
            manifest,
            scope=scope,
            signature_payload=current_payload,
        )
        assert orchestrator._reused_signature_overrides[(stage, scope)] == recorded_signature


def test_imports_legacy_variant_checkpoint_into_pipeline_state(tmp_path: Path):
    manifest = tmp_path / "S1_ConfSearch" / "product_major" / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("{}", encoding="utf-8")
    legacy_path = tmp_path / ".rph" / "checkpoint.json"
    legacy_path.parent.mkdir(parents=True)
    legacy_path.write_text(
        json.dumps({
            "schema_version": "rph_variant_checkpoint_v1",
            "reaction_id": "RXN_1",
            "variants": {
                "product_major": {
                    "s1": {
                        "signature": "legacy-s1",
                        "manifest": str(manifest.resolve()),
                        "status": "complete",
                    }
                }
            },
        }),
        encoding="utf-8",
    )

    checkpoint = V4Checkpoint(tmp_path)
    state = checkpoint.load()

    assert checkpoint.reusable_scope(
        "product_major", "s1", "legacy-s1", manifest
    )
    assert state["schema_version"] == "rph_v4_checkpoint_v2"
    assert state["reaction_id"] == "RXN_1"
    assert state["stages"]["s1"]["scopes"]["product_major"]["status"] == "complete"
    assert legacy_path.exists()


def test_run_config_snapshot_is_immutable(tmp_path: Path):
    checkpoint = V4Checkpoint(tmp_path)
    path = checkpoint.write_config_snapshot({"step1": {"protocol": "censo_lite"}})
    checkpoint.write_config_snapshot({"step1": {"protocol": "changed"}})

    saved = json.loads(path.read_text(encoding="utf-8"))
    assert saved["config"]["step1"]["protocol"] == "censo_lite"


def test_incomplete_s4_checkpoint_is_not_reused(tmp_path: Path):
    checkpoint = V4Checkpoint(tmp_path)
    manifest = tmp_path / "S4_HighLevel" / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text('{"status": "incomplete"}', encoding="utf-8")
    payload = {"signature_schema": "s4_signature_v3", "stage": "s4"}
    signature = checkpoint.signature(payload)
    checkpoint.mark(
        "s4",
        signature,
        manifest,
        status="incomplete",
        signature_payload=payload,
    )

    orchestrator = object.__new__(V4Orchestrator)
    orchestrator.resume_policy = "strict"
    orchestrator.start_from = None
    orchestrator._reused_signature_overrides = {}

    assert not orchestrator._checkpoint_reusable(
        checkpoint,
        "s4",
        signature,
        manifest,
        signature_payload=payload,
    )
