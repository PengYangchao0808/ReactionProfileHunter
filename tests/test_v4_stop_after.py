import json
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

import rph_core.v4_orchestrator as v4_module


REPO_ROOT = Path(__file__).resolve().parents[1]
DATASET = REPO_ROOT / "data" / "reaxys_cleaned.csv"


def test_cli_forwards_stop_after(monkeypatch, tmp_path: Path):
    captured = {}

    class FakeOrchestrator:
        def __init__(
            self,
            config_path=None,
            color=False,
            ui_mode=None,
            resume_policy=None,
            start_from=None,
            recompute_from=None,
        ):
            captured["config"] = config_path
            captured["color"] = color
            captured["ui_mode"] = ui_mode
            captured["resume_policy"] = resume_policy
            captured["start_from"] = start_from
            captured["recompute_from"] = recompute_from

        def run(self, s0_record, output, stop_after="s4"):
            captured.update(
                s0_record=s0_record,
                output=output,
                stop_after=stop_after,
            )
            return {}

    monkeypatch.setattr(v4_module, "V4Orchestrator", FakeOrchestrator)
    assert v4_module.main(
        [
                "--csv",
                str(DATASET),
                "--rx-id",
                "1",
                "--output",
                str(tmp_path),
            "--stop-after",
            "s3",
        ]
    ) == 0
    assert captured["stop_after"] == "s3"
    assert captured["s0_record"].reaction_type == "4+3"


def test_cli_forwards_explicit_resume_policy(monkeypatch, tmp_path: Path):
    captured = {}

    class FakeOrchestrator:
        def __init__(self, _config=None, **kwargs):
            captured.update(kwargs)

        def run(self, _record, _output, stop_after="s4"):
            captured["stop_after"] = stop_after
            return {}

    monkeypatch.setattr(v4_module, "V4Orchestrator", FakeOrchestrator)
    assert v4_module.main([
        "--csv", str(DATASET),
        "--rx-id", "1",
        "--output", str(tmp_path),
        "--stop-after", "s2",
        "--resume-policy", "use-existing-upstream",
        "--start-from", "s2",
    ]) == 0

    assert captured["resume_policy"] == "use-existing-upstream"
    assert captured["start_from"] == "s2"
    assert captured["stop_after"] == "s2"


def test_cli_forwards_dashboard_ui_mode(monkeypatch, tmp_path: Path):
    captured = {}

    class FakeOrchestrator:
        def __init__(self, _config=None, **kwargs):
            captured.update(kwargs)

        def run(self, _record, _output, stop_after="s4"):
            return {}

    monkeypatch.setattr(v4_module, "V4Orchestrator", FakeOrchestrator)
    assert v4_module.main([
        "--csv", str(DATASET),
        "--rx-id", "1",
        "--output", str(tmp_path),
        "--stop-after", "s1",
        "--ui", "dashboard",
    ]) == 0
    assert captured["ui_mode"] == "dashboard"


def test_cli_requires_explicit_policy_for_start_from(tmp_path: Path):
    with pytest.raises(SystemExit):
        v4_module.main([
            "--csv", str(DATASET),
            "--rx-id", "1",
            "--output", str(tmp_path),
            "--stop-after", "s2",
            "--start-from", "s2",
        ])


def test_cli_rejects_product_smiles_mode(tmp_path: Path):
    with pytest.raises(SystemExit):
        v4_module.main([
            "--smiles", "C=C",
            "--output", str(tmp_path),
        ])


def test_partial_pipeline_result_contains_only_completed_manifests(tmp_path: Path):
    manifest = tmp_path / "S3_LowLevel" / "manifest.json"
    manifest.parent.mkdir()
    manifest.write_text("{}", encoding="utf-8")

    result = v4_module.V4Orchestrator._write_pipeline_result(
        tmp_path,
        "s3",
        {"s1": tmp_path / "s1.json", "s3": manifest},
    )
    saved = json.loads((tmp_path / "pipeline.result.json").read_text(encoding="utf-8"))

    assert result["completed_through"] == "s3"
    assert "s3_manifest" in saved
    assert "s4_manifest" not in saved
