from __future__ import annotations

import json
from pathlib import Path

import pytest

import rph_core.steps.refinement.engine as engine_module
from rph_core.steps.fidelity_profile import FidelityProfile
from rph_core.steps.refinement import RefinementEngine
from rph_core.steps.refinement.models import StructureRequest
from rph_core.utils.config_loader import load_config
from rph_core.utils.qc_models import QCJobResult
from rph_core.utils.stage_progress import StageProgressReporter


def _engine(stage: str = "S3") -> RefinementEngine:
    config = load_config()
    return RefinementEngine(config, FidelityProfile.from_config(config, stage))


def _write_xyz(path: Path, distance: float) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"2\nmock\nH 0.0 0.0 0.0\nH {distance:.6f} 0.0 0.0\n",
        encoding="utf-8",
    )
    return path


def _request(
    inputs_dir: Path,
    *,
    structure_id: str,
    role: str,
    kind: str,
    input_distance: float = 1.8,
    seed_distance: float | None = None,
    source_stage: str = "S2",
    missing: bool = False,
) -> StructureRequest:
    input_xyz = inputs_dir / f"{structure_id}.xyz"
    if not missing:
        _write_xyz(input_xyz, input_distance)
    seed_xyz = None
    if seed_distance is not None:
        seed_xyz = _write_xyz(inputs_dir / f"{structure_id}_seed.xyz", seed_distance)
    return StructureRequest(
        id=structure_id,
        role=role,
        kind=kind,
        input_xyz=input_xyz,
        original_seed_xyz=seed_xyz,
        source_stage=source_stage,
        forming_bonds=[(0, 1)] if role in {"intermediate", "ts"} else [],
    )


def _manifest_for(
    engine: RefinementEngine,
    requests: list[StructureRequest],
    output_dir: Path,
) -> dict:
    manifest_path = engine.run(requests, output_dir)
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def test_partial_s3_rerun_preserves_complete_structure_and_recalculates_failed(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    engine = _engine()
    requests = [
        _request(tmp_path / "inputs", structure_id="complete", role="product", kind="minimum"),
        _request(tmp_path / "inputs", structure_id="failed", role="product", kind="minimum"),
    ]
    initial_opt_calls, _ = _install_qc_mocks(monkeypatch)
    output_dir = tmp_path / "s3"
    initial_manifest = _manifest_for(engine, requests, output_dir)
    preserved_dir = output_dir / "complete" / "preserved"
    preserved_dir.mkdir(parents=True)
    canonical_xyz = _write_xyz(preserved_dir / "canonical.xyz", 1.0)
    opt_output = preserved_dir / "output.out"
    opt_output.write_text("preserved\n", encoding="utf-8")
    initial_complete = {
        "id": "complete",
        "status": "complete",
        "canonical_xyz": str(canonical_xyz),
        "opt_output": str(opt_output),
        "canonical_frequency_status": "not_run",
        "sp_status": "not_run",
        "marker": "preserved",
    }
    initial_manifest["structures"] = [
        initial_complete if row["id"] == "complete" else row
        for row in initial_manifest["structures"]
    ]
    initial_opt_calls.clear()

    failed_row = next(row for row in initial_manifest["structures"] if row["id"] == "failed")
    failed_row["status"] = "failed"
    failed_row["canonical_xyz"] = None
    (output_dir / "manifest.json").write_text(
        json.dumps(initial_manifest), encoding="utf-8"
    )

    partial_manifest = json.loads(
        engine.run(requests, output_dir, resume_incomplete=True).read_text(encoding="utf-8")
    )

    entries = {row["id"]: row for row in partial_manifest["structures"]}
    assert len(initial_opt_calls) == 1
    assert "failed" in str(initial_opt_calls[0]["output_dir"])
    assert entries["complete"] == initial_complete
    assert entries["failed"]["status"] in {"complete", "degraded"}
    assert partial_manifest["partial_rerun"]["rerun_structure_ids"] == ["failed"]
    assert (output_dir / "superseded").exists() is False


def _install_qc_mocks(
    monkeypatch: pytest.MonkeyPatch,
    *,
    fail_primary_for: set[str] | None = None,
):
    opt_calls: list[dict] = []
    freq_calls: list[dict] = []
    fail_primary_for = fail_primary_for or set()

    def fake_run_optimization(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
        out_file = output_dir / "opt.out"
        out_file.write_text("mock opt\n", encoding="utf-8")
        call = {
            "spec": spec,
            "input_xyz": Path(input_xyz),
            "output_dir": output_dir,
            "output_xyz": out_xyz,
            "config": config,
        }
        opt_calls.append(call)
        is_primary = spec.route_extras == "" and spec.task in {"opt", "opt_ts"}
        should_fail = is_primary and any(marker in str(output_dir) for marker in fail_primary_for)
        if should_fail:
            return QCJobResult(
                "failed",
                Path(input_xyz),
                output_xyz=out_xyz,
                output_file=out_file,
                energy_hartree=-1.0,
                error="primary opt failed",
            )
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_xyz=out_xyz,
            output_file=out_file,
            energy_hartree=-1.0,
        )

    def fake_run_frequency(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / "freq.out"
        out_file.write_text("mock freq\n", encoding="utf-8")
        frequencies = (-321.4, 120.0, 311.2) if "ts" in str(input_xyz) else (25.0, 125.0, 325.0)
        call = {
            "spec": spec,
            "input_xyz": Path(input_xyz),
            "output_dir": output_dir,
            "config": config,
        }
        freq_calls.append(call)
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_file=out_file,
            energy_hartree=-0.9,
            frequencies_cm1=frequencies,
        )

    def fake_run_single_point(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / "sp.out"
        out_file.write_text("mock sp\n", encoding="utf-8")
        return QCJobResult(
            "complete",
            Path(input_xyz),
            output_file=out_file,
            energy_hartree=-1.1,
        )

    monkeypatch.setattr(engine_module, "run_optimization", fake_run_optimization)
    monkeypatch.setattr(engine_module, "run_frequency", fake_run_frequency)
    monkeypatch.setattr(engine_module, "run_single_point", fake_run_single_point)
    return opt_calls, freq_calls


def test_pass1_precursor_no_warmup_model_hessian(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls = _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    request = _request(
        tmp_path / "inputs",
        structure_id="precursor",
        role="precursor",
        kind="minimum",
    )

    manifest = _manifest_for(engine, [request], tmp_path / "stage")
    entry = manifest["structures"][0]

    assert entry["warmup_status"] == "not_requested"
    assert opt_calls[0]["spec"].route == "Opt"
    assert opt_calls[0]["spec"].initial_hessian == "model"


def test_pass1_product_no_warmup_model_hessian(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls = _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    request = _request(
        tmp_path / "inputs",
        structure_id="product",
        role="product",
        kind="minimum",
    )

    manifest = _manifest_for(engine, [request], tmp_path / "stage")
    entry = manifest["structures"][0]

    assert entry["warmup_status"] == "not_requested"
    assert opt_calls[0]["spec"].route == "Opt"
    assert opt_calls[0]["spec"].initial_hessian == "model"


def test_refinement_engine_emits_live_stage_progress(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    output_dir = tmp_path / "stage"
    reporter = StageProgressReporter(output_dir, "S3")
    engine.set_progress_reporter(reporter)
    request = _request(
        tmp_path / "inputs",
        structure_id="product",
        role="product",
        kind="minimum",
    )
    request.s1_thermochemistry_status = "complete"

    _manifest_for(engine, [request], output_dir)

    status = json.loads((output_dir / "status.json").read_text(encoding="utf-8"))
    structure = status["structures"]["product"]
    events = (output_dir / "events.jsonl").read_text(encoding="utf-8")
    assert status["total_structures"] == 1
    assert structure["status"] == "complete"
    assert structure["tasks"]["optimization"]["status"] == "complete"
    assert structure["tasks"]["frequency"]["status"] == "complete"
    assert structure["tasks"]["single_point"]["status"] == "complete"
    assert "structure_queued" in events
    assert "structure_started" in events
    assert "optimization_started" in events


def test_pass1_intermediate_warmup_calc_hess(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls = _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    request = _request(
        tmp_path / "inputs",
        structure_id="intermediate",
        role="intermediate",
        kind="minimum",
        seed_distance=2.4,
    )

    manifest = _manifest_for(engine, [request], tmp_path / "stage")
    entry = manifest["structures"][0]

    assert entry["warmup_status"] == "complete"
    assert opt_calls[0]["spec"].route_extras == "LooseOpt"
    assert opt_calls[1]["spec"].route == "Opt"
    assert opt_calls[1]["spec"].initial_hessian == "calculate"
    assert opt_calls[1]["spec"].max_cycles == 60


def test_pass1_ts_warmup_optts_calc_hess(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls = _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    request = _request(
        tmp_path / "inputs",
        structure_id="ts",
        role="ts",
        kind="ts",
        seed_distance=2.5,
    )

    manifest = _manifest_for(engine, [request], tmp_path / "stage")
    entry = manifest["structures"][0]

    assert entry["warmup_status"] == "complete"
    assert opt_calls[0]["spec"].route_extras == "LooseOpt"
    assert opt_calls[0]["spec"].trust is None
    assert opt_calls[1]["spec"].route == "OptTS"
    assert opt_calls[1]["spec"].initial_hessian == "calculate"
    assert opt_calls[1]["spec"].trust == pytest.approx(0.15)
    assert opt_calls[1]["spec"].max_cycles == 60


def test_pass1_s3_warmup_uses_s2_seed_distance(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls = _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    request = _request(
        tmp_path / "inputs",
        structure_id="intermediate",
        role="intermediate",
        kind="minimum",
        input_distance=1.8,
        seed_distance=2.4,
        source_stage="S2",
    )

    _manifest_for(engine, [request], tmp_path / "stage")

    assert opt_calls[0]["spec"].bond_constraints[0][2] == pytest.approx(2.4)


def test_pass1_s4_warmup_uses_s3_input_distance(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, _freq_calls = _install_qc_mocks(monkeypatch)
    engine = _engine("S4")
    request = _request(
        tmp_path / "inputs",
        structure_id="intermediate",
        role="intermediate",
        kind="minimum",
        input_distance=1.8,
        seed_distance=2.4,
        source_stage="S3",
    )

    _manifest_for(engine, [request], tmp_path / "stage")

    assert opt_calls[0]["spec"].bond_constraints[0][2] == pytest.approx(1.8)


def test_pass1_frequency_independent_of_opt(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _opt_calls, freq_calls = _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    request = _request(
        tmp_path / "inputs",
        structure_id="precursor",
        role="precursor",
        kind="minimum",
    )

    manifest = _manifest_for(engine, [request], tmp_path / "stage")

    assert freq_calls[0]["input_xyz"] == Path(manifest["structures"][0]["opt_xyz"])
    assert freq_calls[0]["spec"].initial_hessian == "model"


def test_pass1_classification_uses_identity_module(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    inputs_dir = tmp_path / "inputs"
    requests = [
        _request(inputs_dir, structure_id="precursor", role="precursor", kind="minimum"),
        _request(inputs_dir, structure_id="product", role="product", kind="minimum"),
        _request(
            inputs_dir,
            structure_id="intermediate",
            role="intermediate",
            kind="minimum",
            seed_distance=2.4,
        ),
        _request(inputs_dir, structure_id="ts", role="ts", kind="ts", seed_distance=2.5),
    ]
    called = {"ts": 0, "int": 0, "minimum": 0}

    def fake_classify_ts(**kwargs):
        called["ts"] += 1
        return {"kind": "ts", "frequencies": list(kwargs["frequencies_cm1"])}

    def fake_classify_int(**kwargs):
        called["int"] += 1
        return {"kind": "int", "opt_xyz": str(kwargs["opt_xyz"]) if kwargs["opt_xyz"] else None}

    def fake_classify_minimum(**kwargs):
        called["minimum"] += 1
        return {"kind": "minimum", "role": kwargs.get("expected_role")}

    monkeypatch.setattr(engine_module.identity_module, "classify_ts", fake_classify_ts)
    monkeypatch.setattr(engine_module.identity_module, "classify_int", fake_classify_int)
    monkeypatch.setattr(engine_module.identity_module, "classify_minimum", fake_classify_minimum)

    _manifest_for(engine, requests, tmp_path / "stage")

    assert called == {"ts": 1, "int": 1, "minimum": 4}


def test_pass1_opt_failure_does_not_block_frequency_attempt(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    opt_calls, freq_calls = _install_qc_mocks(monkeypatch, fail_primary_for={"precursor"})
    engine = _engine("S3")
    request = _request(
        tmp_path / "inputs",
        structure_id="precursor",
        role="precursor",
        kind="minimum",
    )

    manifest = _manifest_for(engine, [request], tmp_path / "stage")
    entry = manifest["structures"][0]

    assert len(opt_calls) == 1
    assert freq_calls == []
    assert entry["frequency_status"] == "skipped"


def test_pass1_opt_failure_status_degraded(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch, fail_primary_for={"precursor"})
    engine = _engine("S3")
    request = _request(
        tmp_path / "inputs",
        structure_id="precursor",
        role="precursor",
        kind="minimum",
    )

    manifest = _manifest_for(engine, [request], tmp_path / "stage")
    entry = manifest["structures"][0]

    assert entry["status"] in {"degraded", "failed"}


def test_pass1_attempt_history_recorded(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    request = _request(
        tmp_path / "inputs",
        structure_id="ts",
        role="ts",
        kind="ts",
        seed_distance=2.5,
    )

    manifest = _manifest_for(engine, [request], tmp_path / "stage")
    entry = manifest["structures"][0]

    assert len(entry["attempt_history"]) == 3
    assert all(item["attempt_id"] for item in entry["attempt_history"])
    assert entry["current_attempt_id"] == entry["attempt_history"][-1]["attempt_id"]


def test_pass1_manifest_contains_all_pass1_fields(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    inputs_dir = tmp_path / "inputs"
    manifest = _manifest_for(
        engine,
        [
            _request(inputs_dir, structure_id="precursor", role="precursor", kind="minimum"),
            _request(inputs_dir, structure_id="product", role="product", kind="minimum"),
            _request(
                inputs_dir,
                structure_id="intermediate",
                role="intermediate",
                kind="minimum",
                seed_distance=2.4,
            ),
            _request(inputs_dir, structure_id="ts", role="ts", kind="ts", seed_distance=2.5),
        ],
        tmp_path / "stage",
    )

    assert manifest["schema_version"] == "refinement_manifest_v1"
    assert len(manifest["structures"]) == 4
    for entry in manifest["structures"]:
        assert "warmup_status" in entry
        assert "opt_status" in entry
        assert "frequency_status" in entry
        assert "ts_classification" in entry
        assert "int_classification" in entry
        assert "minimum_classification" in entry


def test_pass1_pass2_pass3_stubs_return_placeholders(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    request = _request(
        tmp_path / "inputs",
        structure_id="product",
        role="product",
        kind="minimum",
    )

    manifest = _manifest_for(engine, [request], tmp_path / "stage")
    entry = manifest["structures"][0]

    assert entry["pass2_rescue_attempts"] == []
    assert entry["canonical_attempt_id"] is None
    assert entry["canonical_xyz"] is None
    assert entry["sp_status"] == "not_run"
    assert entry["sp_energy_hartree"] is None
    assert entry["sp_output"] is None
    assert entry["thermochemistry"] is None
    assert entry["ml_usability"] is None
    assert entry["irc_status"] == "not_run"
    assert entry["irc_endpoints"] is None


def test_pass1_s3_and_s4_profiles_differ_in_warmup_cycles(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_qc_mocks(monkeypatch)
    s3_engine = _engine("S3")
    s4_engine = _engine("S4")
    s3_request = _request(
        tmp_path / "inputs_s3",
        structure_id="ts_s3",
        role="ts",
        kind="ts",
        seed_distance=2.5,
    )
    s4_request = _request(
        tmp_path / "inputs_s4",
        structure_id="ts_s4",
        role="ts",
        kind="ts",
        seed_distance=2.5,
        source_stage="S3",
    )

    s3_manifest = _manifest_for(s3_engine, [s3_request], tmp_path / "stage_s3")
    s4_manifest = _manifest_for(s4_engine, [s4_request], tmp_path / "stage_s4")

    assert s3_manifest["structures"][0]["warmup_max_cycles"] == 50
    assert s4_manifest["structures"][0]["warmup_max_cycles"] == 6


def test_pass1_failed_preflight_structure_in_manifest_with_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    _install_qc_mocks(monkeypatch)
    engine = _engine("S3")
    inputs_dir = tmp_path / "inputs"
    requests = [
        _request(inputs_dir, structure_id="precursor", role="precursor", kind="minimum"),
        _request(
            inputs_dir,
            structure_id="broken",
            role="product",
            kind="minimum",
            missing=True,
        ),
    ]

    manifest = _manifest_for(engine, requests, tmp_path / "stage")
    structures = {entry["id"]: entry for entry in manifest["structures"]}

    assert structures["broken"]["preflight_status"] == "failed_preflight"
    assert structures["broken"]["status"] == "failed"
    assert structures["broken"]["error"]
    assert "opt_status" not in structures["broken"]
