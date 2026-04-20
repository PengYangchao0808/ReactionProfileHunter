import argparse
import json
from typing import cast
from pathlib import Path

import pytest
import yaml

from benchmark import phasea_driver


XYZ_TEXT = """1
generated
H 0.0 0.0 0.0
"""


def _write_required_phasea_artifacts(run_dir: Path) -> Path:
    s1_dir = run_dir / "S1_ConfGeneration"
    (s1_dir / "product" / "crest").mkdir(parents=True, exist_ok=True)
    _ = (s1_dir / "provenance.json").write_text('{"ok": true}', encoding="utf-8")
    ensemble = s1_dir / "product" / "crest" / "ensemble.xyz"
    _ = ensemble.write_text(XYZ_TEXT, encoding="utf-8")
    return s1_dir


def test_validate_artifacts_accepts_canonical_flat_product_path(tmp_path: Path):
    run_dir = tmp_path / "run"
    s1_dir = _write_required_phasea_artifacts(run_dir)
    flat_product = s1_dir / "product_min.xyz"
    _ = flat_product.write_text(XYZ_TEXT, encoding="utf-8")

    artifacts = phasea_driver.artifact_candidates(run_dir)

    assert artifacts["product_xyz"] == flat_product
    assert phasea_driver.validate_artifacts(run_dir) == {
        "provenance_json": True,
        "product_min_xyz": True,
        "ensemble_found": True,
    }


def test_validate_artifacts_falls_back_to_legacy_nested_product_path(tmp_path: Path):
    run_dir = tmp_path / "run"
    s1_dir = _write_required_phasea_artifacts(run_dir)
    legacy_product = s1_dir / "product" / "product_min.xyz"
    _ = legacy_product.write_text(XYZ_TEXT, encoding="utf-8")

    artifacts = phasea_driver.artifact_candidates(run_dir)

    assert artifacts["product_xyz"] == legacy_product
    assert phasea_driver.artifacts_complete(run_dir) is True


def _make_phasea_plan(tmp_path: Path, rx_id: str = "42") -> phasea_driver.Paths:
    paths = phasea_driver.build_paths(tmp_path, "phasea_test")
    phasea_driver.ensure_plan_dirs(paths)
    phasea_driver.write_csv(
        phasea_driver.panel_path(paths),
        [{"rx_id": rx_id, "enabled": "1", "tier": "primary", "notes": "", "rxn_key": f"key-{rx_id}"}],
        ["rx_id", "enabled", "tier", "notes", "rxn_key"],
    )
    phasea_driver.write_csv(
        phasea_driver.task_manifest_path(paths),
        [],
        [
            "order_index",
            "task_id",
            "rx_id",
            "protocol",
            "replicate",
            "tier",
            "config_relpath",
            "run_relpath",
            "log_relpath",
            "cache_relpath",
        ],
    )
    phasea_driver.write_json(
        phasea_driver.rep_config_path(paths),
        {"targets": {}, "updated_at": phasea_driver.now_iso()},
    )
    return paths


def _write_task_manifest_rows(
    paths: phasea_driver.Paths,
    rx_id: str,
    replicate: int,
    protocols: list[str],
) -> None:
    phasea_driver.write_csv(
        phasea_driver.task_manifest_path(paths),
        [
            {
                "order_index": str(index),
                "task_id": phasea_driver.normalize_task_id(rx_id, protocol, replicate),
                "rx_id": rx_id,
                "protocol": protocol,
                "replicate": str(replicate),
                "tier": "primary",
                "config_relpath": f"configs/{phasea_driver.normalize_task_id(rx_id, protocol, replicate)}.yaml",
                "run_relpath": f"runs/{phasea_driver.normalize_task_id(rx_id, protocol, replicate)}",
                "log_relpath": f"logs/{phasea_driver.normalize_task_id(rx_id, protocol, replicate)}.log",
                "cache_relpath": f"cache/{phasea_driver.normalize_task_id(rx_id, protocol, replicate)}",
            }
            for index, protocol in enumerate(protocols, start=1)
        ],
        [
            "order_index",
            "task_id",
            "rx_id",
            "protocol",
            "replicate",
            "tier",
            "config_relpath",
            "run_relpath",
            "log_relpath",
            "cache_relpath",
        ],
    )


def test_run_rx_uses_default_protocol_order_and_creates_suite_report(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    paths = _make_phasea_plan(tmp_path)
    _write_task_manifest_rows(paths, "42", 1, phasea_driver.DEFAULT_PROTOCOLS)
    seen_task_ids: list[str] = []
    statuses: dict[str, str] = {}

    def fake_run_task(run_paths: phasea_driver.Paths, task_id: str, preclaimed: bool = False) -> int:
        assert run_paths == paths
        assert preclaimed is False
        seen_task_ids.append(task_id)
        cache_dir = phasea_driver.task_cache_dir(run_paths, task_id)
        cache_dir.mkdir(parents=True, exist_ok=True)
        phasea_driver.write_json(cache_dir / "done.json", {"run_seconds": 1.25})
        phasea_driver.write_json(
            cache_dir / "metrics.json",
            {
                "task_id": task_id,
                "raw_nconf": 4,
                "selected_nconf": 2,
                "selected_candidate_count": 2,
                "mean_native_energy": -10.5,
                "two_stage_enabled": True,
                "selection_mode": "rank1",
                "final_opt_sp_stage_status": "complete",
            },
        )
        statuses[task_id] = phasea_driver.STATUS_COLLECTED
        return 0

    monkeypatch.setattr(phasea_driver, "run_task", fake_run_task)

    def fake_task_status(_paths: phasea_driver.Paths, task_id: str) -> str:
        return statuses[task_id]

    monkeypatch.setattr(phasea_driver, "task_status", fake_task_status)

    rc = phasea_driver.cmd_run_rx(paths, argparse.Namespace(rx_id="42", replicate=1))

    expected_task_ids = [phasea_driver.normalize_task_id("42", protocol, 1) for protocol in phasea_driver.DEFAULT_PROTOCOLS]
    assert rc == 0
    assert seen_task_ids == expected_task_ids

    report = cast(dict[str, object], json.loads((paths.reports_dir / "rx42__rep01_suite_report.json").read_text(encoding="utf-8")))
    protocol_results = cast(list[dict[str, object]], report["protocol_results"])
    first_metrics = cast(dict[str, object], protocol_results[0]["metrics"])
    first_comparison = cast(dict[str, object], protocol_results[0]["comparison_to_ext"])
    assert report["suite_status"] == "DONE"
    assert report["protocol_order"] == phasea_driver.DEFAULT_PROTOCOLS
    assert [item["protocol"] for item in protocol_results] == phasea_driver.DEFAULT_PROTOCOLS
    assert [item["task_id"] for item in protocol_results] == expected_task_ids
    assert first_metrics["raw_nconf"] == 4
    assert first_comparison["delta_mean_native_energy"] == 0.0


def test_run_rx_returns_nonzero_on_partial_failure_but_still_writes_suite_report(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    paths = _make_phasea_plan(tmp_path, rx_id="7")
    _write_task_manifest_rows(paths, "7", 2, phasea_driver.DEFAULT_PROTOCOLS)
    seen_protocols: list[str] = []
    statuses: dict[str, str] = {}

    def fake_run_task(run_paths: phasea_driver.Paths, task_id: str, preclaimed: bool = False) -> int:
        assert preclaimed is False
        _, protocol, _ = phasea_driver.parse_task_id(task_id)
        seen_protocols.append(protocol)
        cache_dir = phasea_driver.task_cache_dir(run_paths, task_id)
        cache_dir.mkdir(parents=True, exist_ok=True)
        if protocol == "lite":
            phasea_driver.write_json(
                cache_dir / "failed.json",
                {"error_message": "lite failed", "run_seconds": 2.5},
            )
            statuses[task_id] = phasea_driver.STATUS_FAILED
            return 1
        phasea_driver.write_json(cache_dir / "done.json", {"run_seconds": 1.0})
        phasea_driver.write_json(cache_dir / "metrics.json", {"task_id": task_id, "mean_native_energy": -9.0})
        statuses[task_id] = phasea_driver.STATUS_COLLECTED
        return 0

    monkeypatch.setattr(phasea_driver, "run_task", fake_run_task)

    def fake_task_status(_paths: phasea_driver.Paths, task_id: str) -> str:
        return statuses[task_id]

    monkeypatch.setattr(phasea_driver, "task_status", fake_task_status)

    rc = phasea_driver.cmd_run_rx(paths, argparse.Namespace(rx_id="7", replicate=2))

    report = cast(dict[str, object], json.loads((paths.reports_dir / "rx7__rep02_suite_report.json").read_text(encoding="utf-8")))
    protocol_results = cast(list[dict[str, object]], report["protocol_results"])
    lite_result = next(item for item in protocol_results if item["protocol"] == "lite")
    zero_result = next(item for item in protocol_results if item["protocol"] == "zero")

    assert rc == 1
    assert seen_protocols == phasea_driver.DEFAULT_PROTOCOLS
    assert report["suite_status"] == "FAILED"
    assert report["completed_protocols"] == 3
    assert report["failed_protocols"] == 1
    assert lite_result["exit_code"] == 1
    assert lite_result["error_message"] == "lite failed"
    assert lite_result["metrics_path"] is None
    assert zero_result["exit_code"] == 0


def test_run_rx_rejects_missing_manifest_rows(tmp_path: Path):
    paths = _make_phasea_plan(tmp_path, rx_id="11")

    with pytest.raises(RuntimeError, match="missing task_manifest rows"):
        _ = phasea_driver.cmd_run_rx(paths, argparse.Namespace(rx_id="11", replicate=1))


def test_collect_task_accepts_nested_reaction_root_layout(tmp_path: Path):
    paths = _make_phasea_plan(tmp_path)
    task_id = phasea_driver.normalize_task_id("42", "ext", 1)
    _write_task_manifest_rows(paths, "42", 1, ["ext"])

    cache_dir = phasea_driver.task_cache_dir(paths, task_id)
    cache_dir.mkdir(parents=True, exist_ok=True)
    phasea_driver.write_json(cache_dir / "done.json", {"run_seconds": 0.25})

    reaction_root = paths.runs_dir / task_id / "RXN_deadbeef"
    s1_dir = _write_required_phasea_artifacts(reaction_root)
    product_min_xyz = s1_dir / "product_min.xyz"
    _ = product_min_xyz.write_text(XYZ_TEXT, encoding="utf-8")

    metrics = phasea_driver.collect_task(paths, task_id)

    assert metrics["run_dir"] == str(paths.runs_dir / task_id)
    assert metrics["product_min_xyz"] == str(product_min_xyz)
    assert metrics["ensemble_path"] == str(s1_dir / "product" / "crest" / "ensemble.xyz")
    assert metrics["artifact_check"] == {
        "provenance_json": True,
        "product_min_xyz": True,
        "ensemble_found": True,
    }


def test_resolve_task_roots_ignores_none_reaction_id_and_picks_nested_root(tmp_path: Path):
    paths = _make_phasea_plan(tmp_path)
    task_id = phasea_driver.normalize_task_id("42", "ext", 1)
    cache_dir = phasea_driver.task_cache_dir(paths, task_id)
    cache_dir.mkdir(parents=True, exist_ok=True)
    phasea_driver.write_json(
        cache_dir / "request.json",
        {
            "task_id": task_id,
            "reaction_id": None,
            "reaction_root": str(paths.runs_dir / task_id),
        },
    )

    reaction_root = paths.runs_dir / task_id / "RXN_deadbeef"
    _write_required_phasea_artifacts(reaction_root)

    resolved = phasea_driver.resolve_task_roots(paths, task_id)

    assert resolved.reaction_root == reaction_root
    assert resolved.reaction_id == "RXN_deadbeef"


def test_resolve_task_roots_rejects_reaction_id_path_traversal(tmp_path: Path):
    paths = _make_phasea_plan(tmp_path)
    task_id = phasea_driver.normalize_task_id("42", "ext", 1)
    cache_dir = phasea_driver.task_cache_dir(paths, task_id)
    cache_dir.mkdir(parents=True, exist_ok=True)
    task_root = paths.runs_dir / task_id
    task_root.mkdir(parents=True, exist_ok=True)
    phasea_driver.write_json(
        cache_dir / "request.json",
        {
            "task_id": task_id,
            "reaction_id": "../../escape",
            "reaction_root": str(task_root),
        },
    )

    resolved = phasea_driver.resolve_task_roots(paths, task_id)

    assert resolved.reaction_root == task_root
    assert resolved.reaction_id is None


def test_candidate_reaction_roots_rejects_symlink_escape(tmp_path: Path):
    base_root = tmp_path / "run"
    base_root.mkdir(parents=True, exist_ok=True)
    external_root = tmp_path / "external" / "RXN_escape"
    external_root.mkdir(parents=True, exist_ok=True)
    symlink_path = base_root / "RXN_escape"
    symlink_path.symlink_to(external_root, target_is_directory=True)

    assert phasea_driver._candidate_reaction_roots(base_root) == []


def test_render_task_config_uses_task_scoped_output_root_per_protocol(tmp_path: Path):
    paths = _make_phasea_plan(tmp_path)
    config_dir = tmp_path / "config"
    config_dir.mkdir(parents=True, exist_ok=True)
    defaults_path = config_dir / "defaults.yaml"
    _ = defaults_path.write_text(
        yaml.safe_dump(
            {
                "run": {"workdir_naming": "rx_{rx_id}"},
                "step1": {},
            },
            allow_unicode=True,
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    defaults = cast(dict[str, object], phasea_driver.load_defaults(paths))
    run_defaults = cast(dict[str, object], defaults["run"])
    default_workdir_naming = cast(str, run_defaults["workdir_naming"])
    ext_task_id = phasea_driver.normalize_task_id("42", "ext", 1)
    lite_task_id = phasea_driver.normalize_task_id("42", "lite", 1)

    ext_cfg = cast(
        dict[str, object],
        yaml.safe_load(phasea_driver.render_task_config(paths, ext_task_id, "42", "ext").read_text(encoding="utf-8")),
    )
    lite_cfg = cast(
        dict[str, object],
        yaml.safe_load(phasea_driver.render_task_config(paths, lite_task_id, "42", "lite").read_text(encoding="utf-8")),
    )
    ext_run_cfg = cast(dict[str, object], ext_cfg["run"])
    lite_run_cfg = cast(dict[str, object], lite_cfg["run"])

    assert ext_run_cfg["output_root"] == str(paths.runs_dir / ext_task_id)
    assert lite_run_cfg["output_root"] == str(paths.runs_dir / lite_task_id)
    assert ext_run_cfg["output_root"] != lite_run_cfg["output_root"]
    assert ext_run_cfg["workdir_naming"] == default_workdir_naming
    assert lite_run_cfg["workdir_naming"] == default_workdir_naming


def test_build_rx_suite_report_writes_expected_files_and_basic_markdown_content(tmp_path: Path):
    paths = _make_phasea_plan(tmp_path, rx_id="9")

    report = phasea_driver.build_rx_suite_report(
        paths,
        "9",
        3,
        [
            {
                "protocol": "ext",
                "task_id": "rx9__ext__rep03",
                "exit_code": 0,
                "task_status": phasea_driver.STATUS_COLLECTED,
                "run_seconds": 1.1,
                "metrics_path": "metrics/ext.json",
                "log_file": "logs/ext.log",
                "error_message": None,
                "metrics": {"raw_nconf": 12, "mean_native_energy": -20.0},
            },
            {
                "protocol": "full",
                "task_id": "rx9__full__rep03",
                "exit_code": 1,
                "task_status": phasea_driver.STATUS_FAILED,
                "run_seconds": 2.2,
                "metrics_path": None,
                "log_file": "logs/full.log",
                "error_message": "full failed",
                "metrics": {},
            },
        ],
    )

    json_path = paths.reports_dir / "rx9__rep03_suite_report.json"
    md_path = paths.reports_dir / "rx9__rep03_suite_report.md"
    markdown = md_path.read_text(encoding="utf-8")

    assert report["report_paths"] == {"json": str(json_path), "markdown": str(md_path)}
    assert json_path.exists()
    assert md_path.exists()
    assert "# Phase A Per-Rx Suite Report: rx9 rep03" in markdown
    assert "Suite status: **FAILED**" in markdown
    assert "### full" in markdown
    assert "- error: full failed" in markdown
    assert "- raw_nconf: 12" in markdown
    assert "- delta_mean_native_energy_vs_ext: 0.0" in markdown
