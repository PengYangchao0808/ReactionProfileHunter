import re
from pathlib import Path

from rph_core.utils.superseded_archive import archive_stale_stage_outputs


def test_archive_stale_stage_outputs_returns_none_for_empty_stage_dir(tmp_path: Path):
    stage_dir = tmp_path / "S3_LowLevel"
    stage_dir.mkdir()

    archived_to = archive_stale_stage_outputs(stage_dir, stage_name="S3")

    assert archived_to is None
    assert not (stage_dir / "superseded").exists()


def test_archive_stale_stage_outputs_moves_structure_directory(tmp_path: Path):
    stage_dir = tmp_path / "S3_LowLevel"
    stale_output = stage_dir / "product_ts" / "opt"
    stale_output.mkdir(parents=True)
    (stale_output / "orca.out").write_text("stale", encoding="utf-8")

    archived_to = archive_stale_stage_outputs(stage_dir, stage_name="S3")

    assert archived_to is not None
    assert re.fullmatch(r"\d{8}T\d{12}Z", archived_to.name)
    assert not (stage_dir / "product_ts").exists()
    assert (archived_to / "product_ts" / "opt" / "orca.out").read_text(encoding="utf-8") == "stale"


def test_archive_stale_stage_outputs_skips_superseded_directory_and_is_idempotent(tmp_path: Path):
    stage_dir = tmp_path / "S4_HighLevel"
    existing_archive = stage_dir / "superseded" / "20240101T000000000000Z"
    existing_archive.mkdir(parents=True)
    (existing_archive / "manifest.json").write_text("old", encoding="utf-8")
    stale_output = stage_dir / "product"
    stale_output.mkdir(parents=True)
    (stale_output / "status.txt").write_text("stale", encoding="utf-8")

    archived_to = archive_stale_stage_outputs(stage_dir, stage_name="S4")

    assert archived_to is not None
    assert existing_archive.exists()
    assert not (archived_to / "superseded").exists()
    assert archive_stale_stage_outputs(stage_dir, stage_name="S4") is None
    archive_names = sorted(path.name for path in (stage_dir / "superseded").iterdir())
    assert archive_names == ["20240101T000000000000Z", archived_to.name]


def test_archive_stale_stage_outputs_preserves_requested_entries(tmp_path: Path):
    stage_dir = tmp_path / "S3_LowLevel"
    stage_dir.mkdir()
    (stage_dir / "status.json").write_text("keep", encoding="utf-8")
    stale_output = stage_dir / "product_int"
    stale_output.mkdir()

    archived_to = archive_stale_stage_outputs(
        stage_dir,
        stage_name="S3",
        preserve=["status.json"],
    )

    assert archived_to is not None
    assert (stage_dir / "status.json").read_text(encoding="utf-8") == "keep"
    assert not (stage_dir / "product_int").exists()
    assert (archived_to / "product_int").is_dir()


def test_archive_stale_stage_outputs_includes_run_id_in_archive_directory_name(tmp_path: Path):
    stage_dir = tmp_path / "S4_HighLevel"
    stage_dir.mkdir()
    (stage_dir / "manifest.json").write_text("stale", encoding="utf-8")

    archived_to = archive_stale_stage_outputs(
        stage_dir,
        stage_name="S4",
        run_id="rxn-001",
    )

    assert archived_to is not None
    assert re.fullmatch(r"\d{8}T\d{12}Z__run-rxn-001", archived_to.name)
    assert (archived_to / "manifest.json").read_text(encoding="utf-8") == "stale"
