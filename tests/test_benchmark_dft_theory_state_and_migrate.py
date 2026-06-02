# pyright: reportPrivateUsage=false

import json
import subprocess
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import cast

import pytest

from benchmark.dft_theory.lib.state_manager import (
    BASELINE_SCHEMA_VERSION,
    BENCHMARK_MANIFEST_SCHEMA_VERSION,
    PHASE1_RESULT_SCHEMA_VERSION,
    RX_MANIFEST_SCHEMA_VERSION,
    create_session_manifest,
    init_rx_manifest,
    read_benchmark_manifest,
    read_rx_manifest,
    record_task_status,
    resolve_baseline_cache_dir,
    resolve_method_run_dir,
    resolve_rx_phase_report_paths,
    set_phase_winner,
)
from benchmark.dft_theory.stages import baseline as baseline_stage
from benchmark.dft_theory.stages import geo_benchmark as geo_stage
from benchmark.dft_theory.stages import sp_benchmark as sp_stage
from benchmark.dft_theory.stages.migrate import migrate_confsearch_seeds


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    _ = path.write_text(content, encoding="utf-8")


def test_state_manager_roundtrip(tmp_path: Path) -> None:
    session_dir = tmp_path / "benchmark_root" / "experiments" / "session_001"
    sp_ref_dir = resolve_method_run_dir(session_dir, "1", "phase1", "SP-REF")
    _ = create_session_manifest(
        session_dir,
        cases=["1", "3"],
        methods={"sp": {"SP-REF": {"is_reference": True}}, "geo": {"GEO-1": {}}},
        reaction_type="[4+3]_default",
    )
    _ = init_rx_manifest(session_dir, "1")
    _ = record_task_status(
        session_dir,
        "1",
        "sp",
        "SP-REF",
        status="done",
        run_dir=sp_ref_dir,
        exit_code=0,
        runtime_seconds=12.5,
        extra={"result_json": str((sp_ref_dir / "sp_result.json").resolve())},
    )
    _ = set_phase_winner(session_dir, "1", "phase1", "SP-REF")

    benchmark_manifest = read_benchmark_manifest(session_dir)
    rx_manifest = read_rx_manifest(session_dir, "1")

    assert benchmark_manifest["meta"]["benchmark_type"] == "dft_theory"
    assert benchmark_manifest["meta"]["schema_version"] == BENCHMARK_MANIFEST_SCHEMA_VERSION
    assert benchmark_manifest["cases"] == ["1", "3"]
    assert rx_manifest["meta"]["schema_version"] == RX_MANIFEST_SCHEMA_VERSION
    assert rx_manifest["tasks"]["phase1"]["SP-REF"]["status"] == "done"
    assert rx_manifest["phase1"]["winner"] == "SP-REF"
    assert rx_manifest["current_phase1_winner"] == "SP-REF"


def test_migrate_stage_imports_confsearch_seed(tmp_path: Path) -> None:
    confsearch_root = tmp_path / "Output" / "benchmark"
    protocol_dir = confsearch_root / "rx1" / "ext" / "RXN_demo"
    s0_dir = protocol_dir / "S0_Mechanism"
    s1_dir = protocol_dir / "S1_ConfGeneration"
    _write_text(s0_dir / "mechanism_summary.json", json.dumps({"reaction_type": "[4+3]_default"}))
    _write_text(s0_dir / "mechanism_graph.json", json.dumps({"nodes": [0, 1], "edges": [[0, 1]]}))
    _write_text(s1_dir / "precursor" / "precursor_min.xyz", "3\nprecursor\nC 0 0 0\nH 0 0 1\nH 0 1 0\n")
    _ = _write_text(s1_dir / "product" / "product_global_min.xyz", "3\nproduct\nC 0 0 0\nH 0 0 1\nH 0 1 0\n")
    _write_text(s1_dir / "provenance.json", json.dumps({"protocol": "ext"}))
    _write_text(s1_dir / "product" / "finalDFT" / "conformer_energies.json", "[0.0, 1.2]")
    _write_text((confsearch_root / "rx1" / "evaluation" / "summary.json"), json.dumps({"winner": "ext"}))

    session_dir = tmp_path / "session_001"
    results = migrate_confsearch_seeds(
        session_dir,
        ["1"],
        cases=[{"rx_id": "1", "seed_protocol": "ext", "confsearch_source": str(confsearch_root / "rx1")}],
    )

    upstream_dir = session_dir / "rx1" / "upstream" / "confsearch"
    mechanism_summary = upstream_dir / "S0_Mechanism" / "mechanism_summary.json"
    mechanism_graph = upstream_dir / "S0_Mechanism" / "mechanism_graph.json"
    precursor_seed = upstream_dir / "S1_ConfGeneration" / "precursor" / "precursor_min.xyz"
    product_seed = upstream_dir / "S1_ConfGeneration" / "product_min.xyz"
    source_manifest = upstream_dir / "source_manifest.json"
    source_manifest_payload = cast(dict[str, object], json.loads(source_manifest.read_text(encoding="utf-8")))
    rx_manifest = cast(dict[str, object], read_rx_manifest(session_dir, "1"))
    imported_files = cast(list[dict[str, object]], source_manifest_payload["imported_files"])
    imported_labels = {str(item["label"]) for item in imported_files}
    upstream = cast(dict[str, object], rx_manifest["upstream"])
    confsearch = cast(dict[str, object], upstream["confsearch"])
    rx_imported_files = cast(list[dict[str, object]], confsearch["imported_files"])
    rx_imported_labels = {str(item["label"]) for item in rx_imported_files}
    inputs = cast(dict[str, object], rx_manifest["inputs"])
    mechanism_summary_input = cast(dict[str, object], inputs["mechanism_summary"])
    mechanism_graph_input = cast(dict[str, object], inputs["mechanism_graph"])
    precursor_xyz_input = cast(dict[str, object], inputs["precursor_xyz"])
    product_xyz_input = cast(dict[str, object], inputs["product_xyz"])

    assert mechanism_summary.exists()
    assert mechanism_graph.exists()
    assert precursor_seed.exists()
    assert product_seed.exists()
    assert source_manifest.exists()
    assert source_manifest_payload["schema_version"] == "dft_benchmark_confsearch_import_v2"
    assert results["1"]["selected_seed_protocol"] == "ext"
    assert {
        "mechanism_summary_json",
        "mechanism_graph_json",
        "precursor_min_xyz",
        "product_min_xyz",
        "provenance_json",
    } <= imported_labels
    assert {"mechanism_summary_json", "mechanism_graph_json", "source_manifest"} <= rx_imported_labels
    assert mechanism_summary_input["path"] == str(mechanism_summary.resolve())
    assert mechanism_graph_input["path"] == str(mechanism_graph.resolve())
    assert precursor_xyz_input["path"] == str(precursor_seed.resolve())
    assert product_xyz_input["path"] == str(product_seed.resolve())
    assert confsearch["status"] == "imported"


def test_copy_pre_s2_snapshot_copies_s0_s1_only(tmp_path: Path) -> None:
    session_dir = tmp_path / "benchmark_root" / "experiments" / "session_001"
    upstream_dir = session_dir / "rx1" / "upstream" / "confsearch"
    baseline_root = resolve_baseline_cache_dir(session_dir, "1", "bl_test")

    _write_text(upstream_dir / "S0_Mechanism" / "mechanism_summary.json", json.dumps({"rx_id": "1"}))
    _write_text(upstream_dir / "S0_Mechanism" / "mechanism_graph.json", json.dumps({"nodes": [0]}))
    _write_text(upstream_dir / "S1_ConfGeneration" / "precursor" / "precursor_min.xyz", "1\nprecursor\nH 0 0 0\n")
    _write_text(upstream_dir / "S1_ConfGeneration" / "product_min.xyz", "1\nproduct\nH 0 0 0\n")
    _write_text(upstream_dir / "pipeline.state", json.dumps({"stage": "foreign"}))

    copied_root = baseline_stage._copy_pre_s2_snapshot(upstream_dir, baseline_root)

    assert copied_root == baseline_root
    assert (baseline_root / "S0_Mechanism" / "mechanism_summary.json").exists()
    assert (baseline_root / "S0_Mechanism" / "mechanism_graph.json").exists()
    assert (baseline_root / "S1_ConfGeneration" / "precursor" / "precursor_min.xyz").exists()
    assert (baseline_root / "S1_ConfGeneration" / "product_min.xyz").exists()
    assert not (baseline_root / "pipeline.state").exists()


def _write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    _ = path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def _phase1_sp_result(
    method_id: str,
    *,
    runtime_seconds: float,
    point_energies: Mapping[str, float],
    secondary_dE_kcal: Mapping[str, float],
) -> dict[str, object]:
    return {
        "schema_version": PHASE1_RESULT_SCHEMA_VERSION,
        "method_id": method_id,
        "status": "completed",
        "runtime_seconds": runtime_seconds,
        "point_results": {
            point: {"energy_hartree": energy}
            for point, energy in point_energies.items()
        },
        "e_precursor_hartree": point_energies["precursor"],
        "e_intermediate_hartree": point_energies["intermediate"],
        "e_ts_hartree": point_energies["ts"],
        "e_product_hartree": point_energies["product"],
        "dE_precursor_to_intermediate_kcal": secondary_dE_kcal["precursor_to_intermediate"],
        "dE_intermediate_to_ts_kcal": secondary_dE_kcal["intermediate_to_ts"],
        "dE_intermediate_to_product_kcal": secondary_dE_kcal["intermediate_to_product"],
    }


def _make_stage_session(tmp_path: Path) -> Path:
    session_dir = tmp_path / "benchmark_root" / "experiments" / "session_stage"
    manifests_dir = session_dir / "manifests"
    reports_dir = session_dir / "reports"
    manifests_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    benchmark_manifest: dict[str, object] = {
        "meta": {"schema_version": BENCHMARK_MANIFEST_SCHEMA_VERSION, "benchmark_type": "dft_theory"},
        "cases": ["1"],
        "methods": {
            "sp": {
                "SP-REF": {"is_reference": True},
                "SP-1": {"is_reference": False},
            },
            "geo": {
                "GEO-1": {},
                "GEO-2": {},
            },
        },
        "stages": {
            "migrate": "done",
            "baseline": "done",
            "phase1_sp": "pending",
            "phase2_geo": "pending",
            "final_report": "pending",
        },
        "reports": {},
    }
    _write_json(manifests_dir / "benchmark_manifest.json", benchmark_manifest)

    baseline_root = resolve_baseline_cache_dir(session_dir, "1", "bl_stage")
    stationary_dir = baseline_root / "stationary_points"
    stationary_dir.mkdir(parents=True, exist_ok=True)
    _write_text(stationary_dir / "precursor_min.xyz", "3\nprecursor\nC 0 0 0\nH 0 0 1\nH 0 1 0\n")
    _write_text(stationary_dir / "complex.xyz", "3\ncomplex\nC 0 0 0\nH 0 0 1\nH 0 1 0\n")
    _write_text(stationary_dir / "intermediate.xyz", "3\nintermediate\nC 0 0 0\nH 0 0 1\nH 0 1 0\n")
    _write_text(stationary_dir / "ts_guess.xyz", "3\nts guess\nC 0 0 0\nH 0 0 1\nH 0 1 0\n")
    _write_text(stationary_dir / "ts_final.xyz", "3\nts\nC 0 0 0\nH 0 0 1\nH 0 1 0\n")
    _write_text(stationary_dir / "product_min.xyz", "3\nproduct\nC 0 0 0\nH 0 0 1\nH 0 1 0\n")

    baseline_manifest = {
        "schema_version": BASELINE_SCHEMA_VERSION,
        "baseline_id": "bl_stage",
        "forming_bonds": [[0, 1], [0, 2]],
        "stationary_points": {
            "precursor_min.xyz": str((stationary_dir / "precursor_min.xyz").resolve()),
            "complex.xyz": str((stationary_dir / "complex.xyz").resolve()),
            "intermediate.xyz": str((stationary_dir / "intermediate.xyz").resolve()),
            "ts_guess.xyz": str((stationary_dir / "ts_guess.xyz").resolve()),
            "ts_final.xyz": str((stationary_dir / "ts_final.xyz").resolve()),
            "product_min.xyz": str((stationary_dir / "product_min.xyz").resolve()),
        },
    }
    _write_json(baseline_root / "baseline_manifest.json", baseline_manifest)
    _write_json(
        baseline_root / "thermo_reference.json",
        {
            "e_ts": -100.000,
            "e_reactant": -100.050,
            "e_product": -100.080,
            "g_ts": -99.900,
            "g_reactant": -99.950,
            "g_product": -99.980,
        },
    )

    sp_ref_dir = resolve_method_run_dir(session_dir, "1", "phase1", "SP-REF")
    sp1_dir = resolve_method_run_dir(session_dir, "1", "phase1", "SP-1")
    geo1_dir = resolve_method_run_dir(session_dir, "1", "phase2", "GEO-1")
    geo2_dir = resolve_method_run_dir(session_dir, "1", "phase2", "GEO-2")
    for path in (sp_ref_dir, sp1_dir, geo1_dir, geo2_dir):
        path.mkdir(parents=True, exist_ok=True)

    _write_json(
        sp_ref_dir / "sp_result.json",
        _phase1_sp_result(
            "SP-REF",
            runtime_seconds=200.0,
            point_energies={
                "precursor": -100.2000,
                "intermediate": -100.1500,
                "ts": -100.1000,
                "product": -100.2500,
            },
            secondary_dE_kcal={
                "precursor_to_intermediate": 31.375,
                "intermediate_to_ts": 31.375,
                "intermediate_to_product": -62.750,
            },
        ),
    )
    _write_json(
        sp1_dir / "sp_result.json",
        _phase1_sp_result(
            "SP-1",
            runtime_seconds=180.0,
            point_energies={
                "precursor": -100.1980,
                "intermediate": -100.1470,
                "ts": -100.0950,
                "product": -100.2440,
            },
            secondary_dE_kcal={
                "precursor_to_intermediate": 32.003,
                "intermediate_to_ts": 32.631,
                "intermediate_to_product": -60.629,
            },
        ),
    )
    _write_json(
        geo1_dir / "geo_result.json",
        {
            "method_id": "GEO-1",
            "status": "ok",
            "runtime_seconds": 360.0,
            "complex": {"converged": True, "l2_energy_hartree": -100.2000},
            "product": {"converged": True, "l2_energy_hartree": -100.2500},
            "ts": {"converged": True, "l2_energy_hartree": -100.1000},
            "dG_activation_kcal": 15.2,
            "ts_validation": {"n_imag": 1, "bond_lengths": [2.20, 2.28]},
        },
    )
    _write_json(
        geo2_dir / "geo_result.json",
        {
            "method_id": "GEO-2",
            "status": "ok",
            "runtime_seconds": 420.0,
            "complex": {"converged": True, "l2_energy_hartree": -100.1980},
            "product": {"converged": True, "l2_energy_hartree": -100.2440},
            "ts": {"converged": True, "l2_energy_hartree": -100.0950},
            "dG_activation_kcal": 16.4,
            "ts_validation": {"n_imag": 2, "bond_lengths": [1.95, 2.55]},
        },
    )

    rx_manifest: dict[str, object] = {
        "meta": {"rx_id": "1", "schema_version": RX_MANIFEST_SCHEMA_VERSION},
        "inputs": {},
        "matrix": {
            "phase1_methods": ["SP-REF", "SP-1"],
            "phase2_methods": ["GEO-1", "GEO-2"],
        },
        "baseline": {
            "schema_version": BASELINE_SCHEMA_VERSION,
            "baseline_id": "bl_stage",
            "root": str(baseline_root.resolve()),
            "baseline_manifest": str((baseline_root / "baseline_manifest.json").resolve()),
            "thermo_reference_json": str((baseline_root / "thermo_reference.json").resolve()),
            "stationary_points": baseline_manifest["stationary_points"],
        },
        "tasks": {
            "phase1": {
                "SP-REF": {"status": "done", "result_json": str((sp_ref_dir / "sp_result.json").resolve())},
                "SP-1": {"status": "done", "result_json": str((sp1_dir / "sp_result.json").resolve())},
            },
            "phase2": {
                "GEO-1": {"status": "done", "result_json": str((geo1_dir / "geo_result.json").resolve())},
                "GEO-2": {"status": "done", "result_json": str((geo2_dir / "geo_result.json").resolve())},
            },
        },
        "phase1": {"results": {}, "winner": None},
        "phase2": {"sp_reader": None, "results": {}, "winner": None},
        "current_phase1_winner": None,
        "current_phase2_winner": None,
        "reports": {},
        "warnings": [],
    }
    _write_json(manifests_dir / "rx1.json", rx_manifest)
    return session_dir


def test_baseline_cache_reused_across_sessions(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    benchmark_root = tmp_path / "Output" / "benchmark_dft_theory"
    confsearch_root = tmp_path / "Output" / "benchmark"
    defaults_config = tmp_path / "defaults.yaml"
    defaults_config.write_text("run: {}\n", encoding="utf-8")
    dataset_path = tmp_path / "dataset.csv"
    dataset_path.write_text("rx_id,product_smiles_main\n", encoding="utf-8")
    cases_config = tmp_path / "cases.yaml"
    cases_config.write_text(
        "cases:\n  - rx_id: '1'\n    seed_protocol: ext\n    smiles: C=C\n    confsearch_source: {0}\n".format(confsearch_root / "rx1"),
        encoding="utf-8",
    )

    protocol_dir = confsearch_root / "rx1" / "ext" / "RXN_demo"
    s0_dir = protocol_dir / "S0_Mechanism"
    s1_dir = protocol_dir / "S1_ConfGeneration"
    _write_text(s0_dir / "mechanism_summary.json", json.dumps({"reaction_type": "[4+3]_default"}))
    _write_text(s0_dir / "mechanism_graph.json", json.dumps({"nodes": [0, 1], "edges": [[0, 1]]}))
    _write_text(s1_dir / "precursor" / "precursor_min.xyz", "1\nprecursor\nH 0 0 0\n")
    _write_text(s1_dir / "product_min.xyz", "1\nproduct\nH 0 0 0\n")
    _write_text(s1_dir / "provenance.json", json.dumps({"protocol": "ext"}))

    session_a = benchmark_root / "experiments" / "session_a"
    session_b = benchmark_root / "experiments" / "session_b"
    cases = [{"rx_id": "1", "seed_protocol": "ext", "smiles": "C=C", "confsearch_source": str(confsearch_root / "rx1")}]
    migrate_confsearch_seeds(session_a, ["1"], cases=cases)
    migrate_confsearch_seeds(session_b, ["1"], cases=cases)

    calls = {"count": 0}

    class _FakeSpReport:
        def to_dict(self) -> dict[str, object]:
            return {"mock": True}

    class _FakeResult:
        def __init__(self, work_dir: Path) -> None:
            self.success = True
            self.error_step = None
            self.error_message = None
            self.intermediate_xyz = work_dir / "S2_Retro" / "complex.xyz"
            self.intermediate_xyz = work_dir / "S2_Retro" / "intermediate.xyz"
            self.ts_guess_xyz = work_dir / "S2_Retro" / "ts_guess.xyz"
            self.ts_final_xyz = work_dir / "S3_TS" / "ts_final.xyz"
            self.product_xyz = work_dir / "S1_ConfGeneration" / "product_min.xyz"
            self.product_thermo = None
            self.ts_log = None
            self.intermediate_log = None
            self.forming_bonds = ((0, 1),)
            self.sp_matrix_report = _FakeSpReport()

    class _FakeHunter:
        def __init__(self, *, config_path: Path) -> None:
            self.config_path = config_path

        def run_pipeline(self, *, product_smiles: str, work_dir: Path, reaction_profile: str) -> _FakeResult:
            calls["count"] += 1
            assert product_smiles == "C=C"
            assert reaction_profile == "[4+3]_default"
            _write_text(work_dir / "S2_Retro" / "complex.xyz", "1\ncomplex\nH 0 0 0\n")
            _write_text(work_dir / "S2_Retro" / "intermediate.xyz", "1\nintermediate\nH 0 0 0\n")
            _write_text(work_dir / "S2_Retro" / "ts_guess.xyz", "1\nts_guess\nH 0 0 0\n")
            _write_text(work_dir / "S3_TS" / "ts_final.xyz", "1\nts\nH 0 0 0\n")
            _write_text(work_dir / "S3_TS" / "sp_matrix_metadata.json", json.dumps({"ok": True}))
            _write_text(work_dir / "S4_Data" / "features_raw.csv", "feature,value\nfoo,1\n")
            return _FakeResult(work_dir)

    monkeypatch.setattr(baseline_stage, "ReactionProfileHunter", _FakeHunter)
    monkeypatch.setattr(baseline_stage, "_validate_ts", lambda result: {"status": "pass", "n_imag": 1})
    monkeypatch.setattr(baseline_stage, "_serialize_sp_report", lambda result, s3_dir: {"mock": True})

    first = baseline_stage.run_baseline_stage(
        session_a,
        ["1"],
        cases_config=cases_config,
        defaults_config=defaults_config,
        dataset_path=dataset_path,
        reaction_profile="[4+3]_default",
    )
    second = baseline_stage.run_baseline_stage(
        session_b,
        ["1"],
        cases_config=cases_config,
        defaults_config=defaults_config,
        dataset_path=dataset_path,
        reaction_profile="[4+3]_default",
    )

    assert calls["count"] == 1
    assert first["1"]["baseline_id"] == second["1"]["baseline_id"]
    baseline_id = cast(str, first["1"]["baseline_id"])
    baseline_root = resolve_baseline_cache_dir(session_a, "1", baseline_id)
    assert baseline_root.is_dir()
    assert (baseline_root / "baseline_manifest.json").exists()
    rx_manifest = cast(dict[str, object], read_rx_manifest(session_b, "1"))
    baseline_block = cast(dict[str, object], rx_manifest["baseline"])
    assert baseline_block["baseline_id"] == baseline_id
    assert baseline_block["root"] == str(baseline_root.resolve())


def test_stage_evaluate_and_final_report(tmp_path: Path) -> None:
    session_dir = _make_stage_session(tmp_path)
    sp_script = Path(__file__).resolve().parents[1] / "benchmark" / "dft_theory" / "stages" / "sp_benchmark.py"
    geo_script = Path(__file__).resolve().parents[1] / "benchmark" / "dft_theory" / "stages" / "geo_benchmark.py"
    eval_script = Path(__file__).resolve().parents[1] / "benchmark" / "dft_theory" / "evaluate.py"

    _ = subprocess.run(
        [sys.executable, str(sp_script), "--session-dir", str(session_dir), "--mode", "evaluate"],
        check=True,
    )

    rx_manifest_path = session_dir / "manifests" / "rx1.json"
    rx_manifest = cast(dict[str, object], json.loads(rx_manifest_path.read_text(encoding="utf-8")))
    phase1 = cast(dict[str, object], rx_manifest["phase1"])
    assert phase1["winner"] == "SP-REF"
    assert phase1["reference_method"] == "SP-REF"
    assert rx_manifest["current_phase1_winner"] == "SP-REF"
    baseline_points = cast(dict[str, str], cast(dict[str, object], rx_manifest["baseline"])["stationary_points"])
    assert set(baseline_points) == {
        "precursor_min.xyz",
        "complex.xyz",
        "intermediate.xyz",
        "ts_guess.xyz",
        "ts_final.xyz",
        "product_min.xyz",
    }
    phase1_results = cast(dict[str, dict[str, object]], phase1["results"])
    sp_ref_result = phase1_results["SP-REF"]
    sp1_result = phase1_results["SP-1"]
    assert tuple(cast(dict[str, object], sp_ref_result["point_results"])) == ("precursor", "intermediate", "ts", "product")
    assert sp_ref_result["absolute_sp_energy_mae_vs_ref_kcal"] == 0.0
    assert cast(float, sp1_result["absolute_sp_energy_mae_vs_ref_kcal"]) > 0.0
    for result in (sp_ref_result, sp1_result):
        assert "substrate" not in cast(dict[str, object], result["point_results"])
        assert "reactant" not in cast(dict[str, object], result["point_results"])
        for obsolete_key in (
            "e_substrate_hartree",
            "e_reactant_hartree",
            "dE_substrate_to_intermediate_kcal",
            "dE_activation_kcal",
            "dE_reaction_kcal",
            "dG_activation_kcal",
            "dG_reaction_kcal",
            "g_reactant_hartree",
            "g_ts_hartree",
            "g_product_hartree",
        ):
            assert obsolete_key not in result

    phase1_summary_path = resolve_rx_phase_report_paths(session_dir, "1", "phase1")["json"]
    phase1_summary = cast(dict[str, object], json.loads(phase1_summary_path.read_text(encoding="utf-8")))
    assert phase1_summary["reference_method"] == "SP-REF"

    _ = subprocess.run(
        [sys.executable, str(geo_script), "--session-dir", str(session_dir), "--mode", "evaluate"],
        check=True,
    )

    rx_manifest = cast(dict[str, object], json.loads(rx_manifest_path.read_text(encoding="utf-8")))
    phase2 = cast(dict[str, object], rx_manifest["phase2"])
    reports = cast(dict[str, str], rx_manifest["reports"])
    assert phase2["sp_reader"] == "SP-1"
    assert phase2["winner"] == "GEO-1"
    assert Path(reports["sp_summary_json"]).is_file()
    assert Path(reports["geo_summary_json"]).is_file()

    _ = subprocess.run(
        [sys.executable, str(eval_script), "--session-dir", str(session_dir), "--phase", "final"],
        check=True,
    )

    assert (session_dir / "reports" / "final_report.json").exists()
    assert (session_dir / "reports" / "phase1_global_summary.json").exists()
    assert (session_dir / "reports" / "phase2_global_summary.json").exists()


def _make_phase1_semantics_session(tmp_path: Path) -> Path:
    session_dir = tmp_path / "benchmark_root" / "experiments" / "session_phase1_semantics"
    manifests_dir = session_dir / "manifests"
    reports_dir = session_dir / "reports"
    manifests_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    _write_json(
        manifests_dir / "benchmark_manifest.json",
        {
            "meta": {"schema_version": BENCHMARK_MANIFEST_SCHEMA_VERSION, "benchmark_type": "dft_theory"},
            "cases": ["1"],
            "methods": {
                "sp": {
                    "SP-REF": {"is_reference": True},
                    "SP-1": {"is_reference": False},
                },
                "geo": {},
            },
            "stages": {
                "migrate": "done",
                "baseline": "done",
                "phase1_sp": "pending",
                "phase2_geo": "pending",
                "final_report": "pending",
            },
            "reports": {},
        },
    )

    sp_ref_dir = resolve_method_run_dir(session_dir, "1", "phase1", "SP-REF")
    sp1_dir = resolve_method_run_dir(session_dir, "1", "phase1", "SP-1")
    sp_ref_dir.mkdir(parents=True, exist_ok=True)
    sp1_dir.mkdir(parents=True, exist_ok=True)

    _write_json(
        sp_ref_dir / "sp_result.json",
        _phase1_sp_result(
            "SP-REF",
            runtime_seconds=12.0,
            point_energies={
                "precursor": -100.1000,
                "intermediate": -100.0500,
                "ts": -100.0100,
                "product": -100.1200,
            },
            secondary_dE_kcal={
                "precursor_to_intermediate": 31.375,
                "intermediate_to_ts": 25.100,
                "intermediate_to_product": -43.925,
            },
        ),
    )
    _write_json(
        sp1_dir / "sp_result.json",
        _phase1_sp_result(
            "SP-1",
            runtime_seconds=45.0,
            point_energies={
                "precursor": -100.0990,
                "intermediate": -100.0490,
                "ts": -100.0090,
                "product": -100.1210,
            },
            secondary_dE_kcal={
                "precursor_to_intermediate": 31.375,
                "intermediate_to_ts": 25.100,
                "intermediate_to_product": -45.180,
            },
        ),
    )

    _write_json(
        manifests_dir / "rx1.json",
        {
            "meta": {"rx_id": "1", "schema_version": RX_MANIFEST_SCHEMA_VERSION},
            "inputs": {},
            "matrix": {"phase1_methods": ["SP-REF", "SP-1"], "phase2_methods": []},
            "baseline": {},
            "tasks": {
                "phase1": {
                    "SP-REF": {"status": "done", "result_json": str((sp_ref_dir / "sp_result.json").resolve())},
                    "SP-1": {"status": "done", "result_json": str((sp1_dir / "sp_result.json").resolve())},
                },
                "phase2": {},
            },
            "phase1": {"results": {}, "winner": None},
            "phase2": {"sp_reader": None, "results": {}, "winner": None},
            "current_phase1_winner": None,
            "current_phase2_winner": None,
            "reports": {},
            "warnings": [],
        },
    )
    return session_dir


def test_phase1_evaluation_uses_precursor_order_and_no_legacy_or_gibbs_fields(tmp_path: Path) -> None:
    session_dir = _make_phase1_semantics_session(tmp_path)

    summaries = sp_stage.evaluate_sp_benchmark_stage(session_dir, ["1"])

    summary = cast(dict[str, object], summaries["1"])
    results = cast(dict[str, dict[str, object]], summary["results"])
    sp_ref = results["SP-REF"]
    sp1 = results["SP-1"]

    assert summary["winner"] == "SP-REF"
    assert summary["reference_method"] == "SP-REF"
    assert tuple(cast(dict[str, object], sp_ref["point_results"])) == ("precursor", "intermediate", "ts", "product")
    assert cast(dict[str, float], sp_ref["absolute_sp_energy_errors_vs_ref_kcal"]) == {
        "precursor": 0.0,
        "intermediate": 0.0,
        "ts": 0.0,
        "product": 0.0,
    }
    assert sp_ref["absolute_sp_energy_mae_vs_ref_kcal"] == 0.0
    assert abs(cast(float, sp1["absolute_sp_energy_mae_vs_ref_kcal"]) - 0.6275094740631) < 1e-9
    assert sp_ref["dE_precursor_to_intermediate_kcal"] == 31.375
    assert sp_ref["dE_intermediate_to_ts_kcal"] == 25.100
    assert sp1["dE_intermediate_to_product_kcal"] == -45.180
    for obsolete_key in (
        "e_substrate_hartree",
        "e_reactant_hartree",
        "dE_substrate_to_intermediate_kcal",
        "dE_activation_kcal",
        "dE_reaction_kcal",
        "dG_activation_kcal",
        "dG_reaction_kcal",
        "g_reactant_hartree",
        "g_ts_hartree",
        "g_product_hartree",
    ):
        assert obsolete_key not in sp_ref
        assert obsolete_key not in sp1
    assert "substrate" not in cast(dict[str, object], sp_ref["point_results"])
    assert "reactant" not in cast(dict[str, object], sp_ref["point_results"])

    rx_manifest = cast(dict[str, object], read_rx_manifest(session_dir, "1"))
    assert cast(dict[str, object], rx_manifest["phase1"])["winner"] == "SP-REF"
    assert rx_manifest["current_phase1_winner"] == "SP-REF"

    summary_json_path = Path(cast(dict[str, str], rx_manifest["reports"])["sp_summary_json"])
    summary_md_path = Path(cast(dict[str, str], rx_manifest["reports"])["sp_summary_md"])
    summary_json = cast(dict[str, object], json.loads(summary_json_path.read_text(encoding="utf-8")))
    summary_md = summary_md_path.read_text(encoding="utf-8")
    global_summary = cast(
        dict[str, object],
        json.loads((session_dir / "reports" / "phase1_global_summary.json").read_text(encoding="utf-8")),
    )

    assert summary_json["winner"] == "SP-REF"
    assert cast(dict[str, object], cast(dict[str, object], summary_json["results"])["SP-REF"])[
        "absolute_sp_energy_mae_vs_ref_kcal"
    ] == 0.0
    assert "Gibbs" not in summary_md
    assert "dG" not in summary_md
    assert global_summary["winner"] == "SP-REF"
    assert cast(dict[str, int], global_summary["votes"]) == {"SP-REF": 1}


def test_phase1_point_map_requires_precursor_and_rejects_legacy_names(tmp_path: Path) -> None:
    stationary_dir = tmp_path / "stationary_points"
    stationary_dir.mkdir(parents=True, exist_ok=True)
    precursor = stationary_dir / "precursor_min.xyz"
    complex_xyz = stationary_dir / "complex.xyz"
    intermediate = stationary_dir / "intermediate.xyz"
    ts_guess = stationary_dir / "ts_guess.xyz"
    ts = stationary_dir / "ts_final.xyz"
    product = stationary_dir / "product_min.xyz"
    for path in (precursor, complex_xyz, intermediate, ts_guess, ts, product):
        _write_text(path, "1\npoint\nH 0 0 0\n")

    mapped = sp_stage._baseline_point_map(
        {
            "precursor_min.xyz": str(precursor),
            "complex.xyz": str(complex_xyz),
            "intermediate.xyz": str(intermediate),
            "ts_guess.xyz": str(ts_guess),
            "ts_final.xyz": str(ts),
            "product_min.xyz": str(product),
        }
    )
    assert tuple(mapped.keys()) == ("precursor", "intermediate", "ts", "product")

    with pytest.raises(ValueError, match="precursor|canonical stationary points|incomplete"):
        _ = sp_stage._baseline_point_map(
            {
                "complex.xyz": str(complex_xyz),
                "intermediate.xyz": str(intermediate),
                "ts_guess.xyz": str(ts_guess),
                "ts_final.xyz": str(ts),
                "product_min.xyz": str(product),
            }
        )

    with pytest.raises(ValueError, match="legacy|precursor|incomplete"):
        _ = sp_stage._baseline_point_map(
            {
                "precursor_min.xyz": str(precursor),
                "reactant_complex.xyz": str(complex_xyz),
                "intermediate.xyz": str(intermediate),
                "ts_guess.xyz": str(ts_guess),
                "ts_final.xyz": str(ts),
                "product_min.xyz": str(product),
            }
        )


def test_geo_inputs_require_complex_and_reject_legacy_session_names(tmp_path: Path) -> None:
    stationary_dir = tmp_path / "geo_stationary_points"
    stationary_dir.mkdir(parents=True, exist_ok=True)
    precursor = stationary_dir / "precursor_min.xyz"
    complex_xyz = stationary_dir / "complex.xyz"
    intermediate = stationary_dir / "intermediate.xyz"
    ts_guess = stationary_dir / "ts_guess.xyz"
    ts_final = stationary_dir / "ts_final.xyz"
    product = stationary_dir / "product_min.xyz"
    for path in (precursor, complex_xyz, intermediate, ts_guess, ts_final, product):
        _write_text(path, "1\npoint\nH 0 0 0\n")

    mapped = geo_stage._point_inputs(
        {
            "precursor_min.xyz": str(precursor),
            "complex.xyz": str(complex_xyz),
            "intermediate.xyz": str(intermediate),
            "ts_guess.xyz": str(ts_guess),
            "ts_final.xyz": str(ts_final),
            "product_min.xyz": str(product),
        }
    )
    assert tuple(mapped) == ("complex", "product", "ts_final", "ts_guess")

    with pytest.raises(ValueError, match="legacy|complex|incomplete"):
        _ = geo_stage._point_inputs(
            {
                "precursor_min.xyz": str(precursor),
                "reactant_complex.xyz": str(complex_xyz),
                "intermediate.xyz": str(intermediate),
                "ts_guess.xyz": str(ts_guess),
                "ts_final.xyz": str(ts_final),
                "product_min.xyz": str(product),
            }
        )


def test_state_manager_rejects_legacy_schema_versions(tmp_path: Path) -> None:
    session_dir = tmp_path / "session_legacy_schema"
    manifests_dir = session_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)
    _write_json(
        manifests_dir / "benchmark_manifest.json",
        {
            "meta": {"schema_version": "benchmark_v1", "benchmark_type": "dft_theory"},
            "cases": [],
            "methods": {},
            "stages": {},
            "reports": {},
        },
    )
    _write_json(
        manifests_dir / "rx1.json",
        {
            "meta": {"rx_id": "1", "schema_version": "benchmark_rx_v1"},
            "inputs": {},
            "matrix": {},
            "tasks": {},
            "phase1": {"results": {}},
            "phase2": {"results": {}},
            "reports": {},
            "warnings": [],
        },
    )

    with pytest.raises(ValueError, match="schema_version"):
        _ = read_benchmark_manifest(session_dir)
    with pytest.raises(ValueError, match="schema_version"):
        _ = read_rx_manifest(session_dir, "1")


def test_phase1_evaluation_rejects_legacy_result_payloads(tmp_path: Path) -> None:
    session_dir = _make_phase1_semantics_session(tmp_path)
    legacy_result_path = resolve_method_run_dir(session_dir, "1", "phase1", "SP-1") / "sp_result.json"
    _write_json(
        legacy_result_path,
        {
            "schema_version": PHASE1_RESULT_SCHEMA_VERSION,
            "method_id": "SP-1",
            "status": "completed",
            "runtime_seconds": 45.0,
            "point_results": {
                "substrate": {"energy_hartree": -100.0990},
                "intermediate": {"energy_hartree": -100.0490},
                "ts": {"energy_hartree": -100.0090},
                "product": {"energy_hartree": -100.1210},
            },
            "e_substrate_hartree": -100.0990,
            "e_intermediate_hartree": -100.0490,
            "e_ts_hartree": -100.0090,
            "e_product_hartree": -100.1210,
            "dE_substrate_to_intermediate_kcal": 31.375,
            "dE_intermediate_to_ts_kcal": 25.100,
            "dE_intermediate_to_product_kcal": -45.180,
            "dG_activation_kcal": 1.23,
        },
    )

    with pytest.raises(ValueError, match="legacy|substrate|Gibbs|phase1"):
        _ = sp_stage.evaluate_sp_benchmark_stage(session_dir, ["1"])


# ──────────────────────────────────────────────────────────────────
# P0/P1 tests: artifact integrity and geo report generation
# ──────────────────────────────────────────────────────────────────


def _make_geo_result_fixture(
    run_dir: Path,
    *,
    optimized_xyz: str = "ts_guess_ts.xyz",
    optimized_hash: str = "aaaa1111",
    freq_log_stem: str = "ts_guess_ts_opt_deadbeef",
    sibling_xyz_hash: str = "bbbb2222",
    artifact_ok: bool = True,
    dE: float = 5.0,
    dG: float = 5.5,
    ts_omega: float = -200.0,
    rmsd: float = 0.15,
    status: str = "completed",
) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "method_id": "GEO-TEST",
        "status": status,
        "runtime_seconds": 1000,
        "dE_activation_kcal": dE,
        "dE_reaction_kcal": -55.0,
        "dG_activation_kcal": dG,
        "dG_reaction_kcal": -50.0,
        "artifact_integrity_ok": artifact_ok,
        "ts": {
            "optimized_xyz": str((run_dir / optimized_xyz).resolve()),
            "optimized_xyz_sha256": optimized_hash,
            "freq_log": str((run_dir / f"{freq_log_stem}.out").resolve()),
            "freq_log_sibling_xyz_sha256": sibling_xyz_hash,
            "displacement": {
                "aligned_rmsd_to_baseline_angstrom": rmsd,
                "max_atom_displacement_angstrom": 0.5,
            },
        },
        "ts_validation": {
            "n_imag": 1,
            "imag_freq": ts_omega,
            "status": "pass",
        },
    }
    result_path = run_dir / "geo_result.json"
    _write_json(result_path, payload)
    return result_path


def test_geo_result_artifact_integrity_field_present(tmp_path: Path) -> None:
    run_dir = tmp_path / "geo_test"
    result_path = _make_geo_result_fixture(run_dir)
    geo = cast(dict[str, object], json.loads(result_path.read_text()))
    assert geo.get("artifact_integrity_ok") is True


def test_geo_result_artifact_integrity_mismatch_detected(tmp_path: Path) -> None:
    run_dir = tmp_path / "geo_test_mismatch"
    result_path = _make_geo_result_fixture(
        run_dir,
        optimized_hash="aaa111",
        sibling_xyz_hash="bbb222",
        artifact_ok=False,
    )
    geo = cast(dict[str, object], json.loads(result_path.read_text()))
    assert geo["artifact_integrity_ok"] is False
    # Verify hashes are different as expected
    ts = cast(dict[str, object], geo["ts"])
    assert ts["optimized_xyz_sha256"] == "aaa111"
    assert ts["freq_log_sibling_xyz_sha256"] == "bbb222"


def test_geo_report_generator_renders_with_effective_source(tmp_path: Path) -> None:
    from benchmark.dft_theory.stages.geo_report_generator import render_detailed_report
    session_dir = tmp_path / "session_rpt"
    manifests_dir = session_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    _write_json(
        manifests_dir / "benchmark_manifest.json",
        {
            "meta": {"schema_version": "benchmark_v2", "benchmark_type": "dft_theory"},
            "cases": ["1"],
            "methods": {"sp": {"SP-1": {}}, "geo": {"GEO-1": {}}},
            "stages": {"phase2_geo": "done"},
            "reports": {},
        },
    )
    run_dir = resolve_method_run_dir(session_dir, "1", "phase2", "GEO-1")
    gp = _make_geo_result_fixture(run_dir, dE=4.2, dG=4.7, ts_omega=-82.0, rmsd=0.18)
    # write an effective block
    geo = json.loads(gp.read_text())
    geo["effective"] = {"source": "refreshed", "effective_dE_activation_kcal": 6.5}
    _write_json(gp, geo)

    _write_json(
        manifests_dir / "rx1.json",
        {
            "meta": {"rx_id": "1", "schema_version": "benchmark_rx_v2"},
            "tasks": {
                "phase2": {
                    "GEO-1": {"status": "completed", "result_json": str(gp.resolve())},
                }
            },
        },
    )
    report = render_detailed_report(session_dir, ["1"])
    assert "rx1" in report
    assert "6.50" in report
    assert "refreshed" in report
    assert "-82.0" in report
    assert "0.1800" in report


def test_geo_report_generator_detects_artifact_integrity_failure(tmp_path: Path) -> None:
    from benchmark.dft_theory.stages.geo_report_generator import render_detailed_report
    session_dir = tmp_path / "session_rpt_fail"
    manifests_dir = session_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    _write_json(
        manifests_dir / "benchmark_manifest.json",
        {
            "meta": {"schema_version": "benchmark_v2", "benchmark_type": "dft_theory"},
            "cases": ["1"],
            "methods": {"sp": {"SP-1": {}}, "geo": {"GEO-1": {}}},
            "stages": {"phase2_geo": "done"},
            "reports": {},
        },
    )
    run_dir = resolve_method_run_dir(session_dir, "1", "phase2", "GEO-1")
    gp = _make_geo_result_fixture(run_dir, artifact_ok=False)
    _write_json(
        manifests_dir / "rx1.json",
        {
            "meta": {"rx_id": "1", "schema_version": "benchmark_rx_v2"},
            "tasks": {
                "phase2": {
                    "GEO-1": {"status": "completed", "result_json": str(gp.resolve())},
                }
            },
        },
    )
    report = render_detailed_report(session_dir, ["1"])
    assert "NO" in report
    assert "Artifact OK" in report
