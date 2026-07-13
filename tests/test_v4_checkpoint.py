from pathlib import Path

from rph_core.utils.v4_checkpoint import V4Checkpoint


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
