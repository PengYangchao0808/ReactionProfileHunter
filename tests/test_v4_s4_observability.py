import json
from pathlib import Path
from types import SimpleNamespace

import rph_core.steps.refinement.engine as refinement_engine_module
from rph_core.steps.step4_highlevel import HighLevelEngine
from rph_core.utils import qc_jobs
from rph_core.utils.config_loader import load_config
from rph_core.utils.qc_models import QCJobResult, QCJobSpec
from rph_core.utils.s4_progress import S4ProgressReporter
from rph_core.v4_watch import load_status, render_status


def _install_successful_s4_qc_mocks(monkeypatch):
    def fake_opt(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_xyz = output_dir / "opt.xyz"
        out_xyz.write_text(Path(input_xyz).read_text(encoding="utf-8"), encoding="utf-8")
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_xyz=out_xyz,
            output_file=output_dir / "opt.out",
            energy_hartree=-1.0,
        )

    def fake_frequency(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=output_dir / "freq.out",
            energy_hartree=-1.0,
            frequencies_cm1=(25.0, 125.0, 325.0),
            enthalpy_hartree=-0.98,
            gibbs_free_energy_hartree=-1.04,
            gibbs_correction_hartree=-0.04,
        )

    def fake_sp(spec, input_xyz, output_dir, config, subprocess_callback=None):
        del spec, config, subprocess_callback
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        return QCJobResult(
            status="complete",
            input_xyz=Path(input_xyz),
            output_file=output_dir / "sp.out",
            energy_hartree=-1.23,
        )

    monkeypatch.setattr(refinement_engine_module, "run_optimization", fake_opt)
    monkeypatch.setattr(refinement_engine_module, "run_frequency", fake_frequency)
    monkeypatch.setattr(refinement_engine_module, "run_single_point", fake_sp)
def test_s4_writes_live_status_events_and_manifest(monkeypatch, tmp_path: Path):
    xyz = tmp_path / "seed.xyz"
    xyz.write_text("1\nseed\nH 0.0 0.0 0.0\n", encoding="utf-8")
    _install_successful_s4_qc_mocks(monkeypatch)

    config = load_config()
    stage_dir = tmp_path / "S4_HighLevel"
    structures = [
        {
            "id": "product_conf_0001",
            "role": "product",
            "kind": "minimum",
            "input_xyz": str(xyz),
            "s1_thermochemistry_status": "complete",
        }
    ]
    manifest = HighLevelEngine(config).run(structures, stage_dir)
    saved_manifest = json.loads(manifest.read_text(encoding="utf-8"))

    reporter = S4ProgressReporter(stage_dir, structures, config=config)
    reporter.start_structure("product_conf_0001")
    reporter.calculator_event(
        "optimization_started",
        {
            "structure_id": "product_conf_0001",
            "engine": "orca",
            "method": "M062X",
            "solvent": "acetone",
            "solvent_model": "CPCM",
        },
    )
    reporter.calculator_event(
        "optimization_finished",
        {
            "structure_id": "product_conf_0001",
            "engine": "orca",
            "method": "M062X",
            "solvent": "acetone",
            "solvent_model": "CPCM",
            "status": "complete",
        },
    )
    reporter.calculator_event(
        "single_point_started",
        {
            "structure_id": "product_conf_0001",
            "engine": "orca",
            "method": "wB97M-V",
        },
    )
    reporter.calculator_event(
        "single_point_finished",
        {
            "structure_id": "product_conf_0001",
            "engine": "orca",
            "method": "wB97M-V",
            "status": "complete",
            "energy_hartree": -1.23,
        },
    )
    reporter.finish_structure(saved_manifest["structures"][0])
    reporter.finish_stage(manifest, saved_manifest["structures"])

    status = json.loads((stage_dir / "status.json").read_text(encoding="utf-8"))
    events = (stage_dir / "events.jsonl").read_text(encoding="utf-8").splitlines()
    assert status["status"] == "completed"
    assert status["summary"]["usable_for_ml"] == 1
    assert status["structures"][0]["tasks"]["optimization"]["status"] == "complete"
    assert status["structures"][0]["tasks"]["single_point"]["status"] == "complete"
    assert any('"event": "optimization_started"' in line for line in events)
    assert any('"event": "stage_finished"' in line for line in events)
    assert saved_manifest["schema_version"] == "refinement_manifest_v1"
    assert saved_manifest["summary"]["complete"] == 1
    assert (stage_dir / "s4.log").exists()

    viewer_status = load_status(tmp_path)
    rendered = render_status(viewer_status, tmp_path, colour=False, event_count=2)
    assert "RPH V4 S4: COMPLETED" in rendered
    assert "product_conf_0001" in rendered


def test_gaussian_frequency_route_uses_explicit_cpcm_and_parses_frequencies(monkeypatch, tmp_path: Path):
    xyz = tmp_path / "ts.xyz"
    xyz.write_text("1\nts\nH 0.0 0.0 0.0\n", encoding="utf-8")
    captured = {}

    class FakeGaussian:
        def __init__(self, **kwargs):
            captured["init"] = kwargs

        def optimize(self, input_xyz, output_dir, route, timeout, charge, spin):
            captured["route"] = route
            return SimpleNamespace(
                converged=True,
                output_file=Path(output_dir) / "freq.log",
                energy=-10.0,
                frequencies=(-345.6, 120.0, 345.0),
                error_message=None,
            )

    monkeypatch.setattr(qc_jobs, "GaussianInterface", FakeGaussian)
    spec = QCJobSpec(
        engine="gaussian",
        task="freq",
        method="M062X",
        basis="def2-SVP",
        solvent="acetone",
        solvent_model="CPCM",
        grid="UltraFine",
        scf="XQC",
    )
    result = qc_jobs.run_frequency(spec, xyz, tmp_path / "freq", {"resources": {"nproc": 4, "mem": "4GB"}})

    assert result.status == "complete"
    assert result.frequencies_cm1 == (-345.6, 120.0, 345.0)
    assert "Freq" in captured["route"]
    assert "SCRF=(CPCM,Solvent=acetone)" in captured["route"]
    assert "M062X/def2SVP" in captured["route"]
    assert "Int=UltraFine" in captured["route"]
    assert "SCF=XQC" in captured["route"]


def test_s4_marks_failed_optimization_incomplete_and_unusable(monkeypatch, tmp_path: Path):
    xyz = tmp_path / "seed.xyz"
    xyz.write_text("1\nseed\nH 0.0 0.0 0.0\n", encoding="utf-8")

    monkeypatch.setattr(
        refinement_engine_module,
        "run_optimization",
        lambda spec, input_xyz, output_dir, config, subprocess_callback=None: QCJobResult(
            status="failed",
            input_xyz=Path(input_xyz),
            output_file=Path(output_dir) / "opt.out",
            error="optimization failed",
        ),
    )
    monkeypatch.setattr(
        refinement_engine_module,
        "run_single_point",
        lambda *args, **kwargs: QCJobResult(status="complete", input_xyz=Path(args[1]), energy_hartree=-1.23),
    )
    manifest = HighLevelEngine(load_config()).run(
        [
            {
                "id": "product",
                "role": "product",
                "kind": "minimum",
                "input_xyz": str(xyz),
                "s1_thermochemistry_status": "complete",
            }
        ],
        tmp_path / "S4_HighLevel",
    )

    saved = json.loads(manifest.read_text(encoding="utf-8"))
    assert saved["summary"].get("complete", 0) == 0
    assert saved["summary"].get("failed", 0) == 1
    assert saved["structures"][0]["usable_for_ml"] is False


def test_highlevel_engine_uses_s4_profile_defaults():
    engine = HighLevelEngine(load_config())

    assert engine.profile.stage == "S4"
    assert engine.profile.fidelity == "high"
    assert engine.profile.profile_id == "m062x_wb97mv_v1"
