import json
from pathlib import Path
from types import SimpleNamespace

from rph_core.steps.step4_highlevel import engine as highlevel_module
from rph_core.steps.step4_highlevel.engine import HighLevelEngine
from rph_core.utils import qc_jobs
from rph_core.utils.qc_models import QCJobSpec
from rph_core.v4_watch import load_status, render_status


def test_s4_writes_live_status_events_and_manifest(monkeypatch, tmp_path: Path):
    xyz = tmp_path / "seed.xyz"
    xyz.write_text("1\nseed\nH 0.0 0.0 0.0\n", encoding="utf-8")

    class FakeCalculator:
        def __init__(self, config, theory, event_callback=None):
            self.event_callback = event_callback

        def run_structure(self, structure, output_dir):
            payload = {
                "structure_id": structure["id"],
                "engine": "gaussian",
                "method": "M062X",
                "solvent": "acetone",
                "solvent_model": "CPCM",
            }
            self.event_callback("optimization_started", payload)
            self.event_callback("optimization_finished", {**payload, "status": "complete"})
            self.event_callback("single_point_started", {**payload, "engine": "orca", "method": "wB97M-V"})
            self.event_callback(
                "single_point_finished",
                {**payload, "engine": "orca", "method": "wB97M-V", "status": "complete", "energy_hartree": -1.23},
            )
            return {
                "id": structure["id"],
                "kind": structure["kind"],
                "input_xyz": structure["input_xyz"],
                "opt_status": "complete",
                "sp_status": "complete",
                "sp_energy_hartree": -1.23,
                "frequency_status": "not_requested",
                "ts_frequency_valid": None,
                "status": "complete",
            }

    monkeypatch.setattr(highlevel_module, "StageCalculator", FakeCalculator)
    stage_dir = tmp_path / "S4_HighLevel"
    manifest = HighLevelEngine({"theory": {"s4_high_precision": {}}}).run(
        [{"id": "product_conf_0001", "kind": "minimum", "input_xyz": str(xyz)}],
        stage_dir,
    )

    status = json.loads((stage_dir / "status.json").read_text(encoding="utf-8"))
    events = (stage_dir / "events.jsonl").read_text(encoding="utf-8").splitlines()
    saved_manifest = json.loads(manifest.read_text(encoding="utf-8"))
    assert status["status"] == "completed"
    assert status["summary"]["usable_for_ml"] == 1
    assert status["structures"][0]["tasks"]["optimization"]["status"] == "complete"
    assert status["structures"][0]["tasks"]["single_point"]["status"] == "complete"
    assert any('"event": "optimization_started"' in line for line in events)
    assert any('"event": "stage_finished"' in line for line in events)
    assert saved_manifest["schema_version"] == "s4_high_level_v2"
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
    assert "Integral=UltraFine" in captured["route"]
    assert "SCF=XQC" in captured["route"]
