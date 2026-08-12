from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

import rph_core.v4_orchestrator as v4_orchestrator
import rph_core.steps.refinement.engine as refinement_engine_module
from rph_core.steps.conformer_search import censo_lite as censo_lite_module
from rph_core.steps.conformer_search.censo_lite import CensoLiteEngine
from rph_core.steps.conformer_search.deduplicator import DedupCandidate
from rph_core.steps.conformer_search.torsion_signature import TorsionSignature
from rph_core.steps.conformer_search.xtb_thermo import XTBThermoResult
from rph_core.steps.mechanism_classifier.s0_record import S0ReactionRecord
from rph_core.steps.step3_lowlevel import LowLevelEngine
from rph_core.steps.step4_highlevel import HighLevelEngine
from rph_core.utils.config_loader import load_config
from rph_core.utils.qc_models import QCJobResult
from rph_core.utils.run_id import RUN_ID_FIELD, is_valid_run_id, new_run_id
from rph_core.utils.s4_progress import S4ProgressReporter
from rph_core.utils.stage_progress import StageProgressReporter
from rph_core.utils.v4_checkpoint import STATE_SCHEMA_V2, V4Checkpoint


def _record(tmp_path: Path) -> S0ReactionRecord:
    return S0ReactionRecord(
        rx_id="RXN_RUN_ID",
        product_smiles="C1CCC1",
        reaction_type="4+3",
        mapped_product_smiles="[CH2:1]1[CH2:2][CH2:3][CH2:4]1",
        mapped_forming_bonds=((1, 2), (3, 4)),
        forming_bonds=((0, 1), (2, 3)),
        mapping_confidence=0.99,
        mapping_trusted=True,
        source_csv=tmp_path / "dataset.csv",
        canonical_precursor_smiles="CC",
        mapped_precursor_smiles="[CH3:1][CH3:2]",
        topology="INTER",
        cyclo_mode="UNKNOWN",
        source_row_hash="row-hash",
    )


def test_new_run_id_returns_unique_uuids() -> None:
    first = new_run_id()
    second = new_run_id()

    assert first != second
    assert is_valid_run_id(first)
    assert is_valid_run_id(second)


def test_is_valid_run_id_accepts_valid_uuid_and_rejects_garbage() -> None:
    assert is_valid_run_id("11111111-1111-4111-8111-111111111111")
    assert not is_valid_run_id("not-a-uuid")
    assert not is_valid_run_id(None)
    assert not is_valid_run_id(123)


def test_orchestrator_generates_run_id_once_per_run_and_persists_outputs(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    run_ids = [
        "11111111-1111-4111-8111-111111111111",
        "22222222-2222-4222-8222-222222222222",
    ]
    calls: list[str] = []

    def fake_new_run_id() -> str:
        value = run_ids[len(calls)]
        calls.append(value)
        return value

    monkeypatch.setattr(v4_orchestrator, "new_run_id", fake_new_run_id)

    orchestrator = v4_orchestrator.V4Orchestrator(color=False)
    result_one = orchestrator.run(_record(tmp_path), tmp_path / "run_one", stop_after="s0")
    result_two = orchestrator.run(_record(tmp_path), tmp_path / "run_two", stop_after="s0")

    assert calls == run_ids
    assert result_one[RUN_ID_FIELD] == run_ids[0]
    assert result_two[RUN_ID_FIELD] == run_ids[1]

    for run_dir, expected in (
        (tmp_path / "run_one", run_ids[0]),
        (tmp_path / "run_two", run_ids[1]),
    ):
        for relative_path in (
            Path("pipeline.result.json"),
            Path("run.manifest.json"),
            Path("pipeline.state"),
            Path(".rph_schema.json"),
            Path("S0_Mechanism/mechanism.json"),
            Path("S0_Mechanism/status.json"),
        ):
            payload = json.loads((run_dir / relative_path).read_text(encoding="utf-8"))
            assert payload[RUN_ID_FIELD] == expected


def test_s1_manifest_includes_run_id(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    class FakeRuntime:
        def __init__(self, _config: object, work_dir: Path, molecule_name: str):
            self.molecule_dir = Path(work_dir) / molecule_name
            self.raw_dir = self.molecule_dir / "candidates_raw"
            self.crest_dir = self.molecule_dir / "crest"
            self.molecule_dir.mkdir(parents=True, exist_ok=True)
            self.raw_dir.mkdir(parents=True, exist_ok=True)
            self.crest_dir.mkdir(parents=True, exist_ok=True)

        def embed(self, _smiles: str) -> Path:
            path = self.molecule_dir / "initial.xyz"
            path.write_text("1\ninitial\nH 0.0 0.0 0.0\n", encoding="utf-8")
            (self.molecule_dir / "initial_atom_mapping.json").write_text(
                json.dumps(
                    {
                        "schema_version": "s1_atom_mapping_v1",
                        "atoms": [
                            {"smiles_idx": 0, "xyz_idx": 0, "element": "H", "type": "hydrogen"}
                        ],
                    }
                ),
                encoding="utf-8",
            )
            return path

        def crest_search(self, input_xyz: Path) -> Path:
            return input_xyz

        def split_ensemble(self, _ensemble: Path) -> list[Path]:
            paths = []
            for index in (1, 2):
                path = self.raw_dir / f"conf_{index:04d}.xyz"
                path.write_text(
                    f"1\nenergy = -10.{index}\nH {index}.0 0.0 0.0\n",
                    encoding="utf-8",
                )
                paths.append(path)
            return paths

        def extract_energy(self, xyz_path: Path) -> float:
            return -10.0 if xyz_path.stem.endswith("0001") else -9.0

        def read_crest_relative_energies(self, _crest_dir: Path):
            return None

        def validate_energies(self, xyz_energies, rel_energies, tolerance_kcal=1.0):
            del xyz_energies, rel_energies, tolerance_kcal

        def run_sp(self, xyz_path: Path, nprocs=None, maxcore=None) -> float:
            del nprocs, maxcore
            return -20.0 if xyz_path.stem.endswith("0001") else -19.0

        def run_mrrho(self, _xyz_path: Path, nprocs=None) -> XTBThermoResult:
            del nprocs
            return XTBThermoResult(
                g_rrho_correction_hartree=0.1,
                xtb_electronic_energy_hartree=-10.0,
                xtb_total_free_energy_hartree=-9.9,
                energy_ledger_residual_hartree=0.0,
                temperature_k=298.15,
                gfn_level=2,
                solvent_model="ALPB(acetone)",
            )

    class FakeDeduplicator:
        def annotate(self, _mol, path: Path, score: float, metadata: dict[str, object]) -> DedupCandidate:
            return DedupCandidate(
                path=Path(path),
                score=float(score),
                signature=TorsionSignature((), (), ()),
                metadata=dict(metadata),
            )

        def deduplicate(self, _mol, candidates: list[DedupCandidate]) -> list[DedupCandidate]:
            return list(candidates)

    run_id = "33333333-3333-4333-8333-333333333333"
    monkeypatch.setattr(censo_lite_module, "CensoLiteRuntime", FakeRuntime)

    engine = CensoLiteEngine(
        {"step1": {"protocol": "censo_lite", "censo_lite": {"retention": {"energy_window_kcal": 10.0}}}},
        tmp_path / "S1_ConfSearch",
        "product",
        run_id=run_id,
    )
    monkeypatch.setattr(engine, "deduplicator", FakeDeduplicator(), raising=False)

    result = engine.run("CCCC")
    manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))

    assert manifest[RUN_ID_FIELD] == run_id


def test_s2_s3_s4_manifests_and_status_include_run_id(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    run_id = "44444444-4444-4444-8444-444444444444"

    ts_guess = tmp_path / "ts_guess.xyz"
    intermediate = tmp_path / "intermediate.xyz"
    scan_profile = tmp_path / "scan_profile.json"
    for path in (ts_guess, intermediate):
        path.write_text("1\nH\nH 0.0 0.0 0.0\n", encoding="utf-8")
    scan_profile.write_text(json.dumps({"selections": {}}), encoding="utf-8")

    s2_manifest = v4_orchestrator.V4Orchestrator._write_s2(
        tmp_path / "S2_PEB",
        ts_guess,
        intermediate,
        [(0, 1), (2, 3)],
        scan_profile,
        "COMPLETE",
        "high",
        [],
        run_id=run_id,
    )
    assert json.loads(s2_manifest.read_text(encoding="utf-8"))[RUN_ID_FIELD] == run_id

    xyz = tmp_path / "seed.xyz"
    xyz.write_text("1\nseed\nH 0.0 0.0 0.0\n", encoding="utf-8")
    structures = [{"id": "product", "kind": "minimum", "input_xyz": str(xyz)}]

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
            energy_hartree=-0.9,
            frequencies_cm1=(25.0, 125.0, 325.0),
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
    config = load_config()
    s3_dir = tmp_path / "S3_LowLevel"
    s3_engine = LowLevelEngine(
        config,
        run_id=run_id,
    )
    s3_manifest = s3_engine.run(
        [{**structures[0], "role": "product", "s1_thermochemistry_status": "complete"}],
        s3_dir,
    )
    StageProgressReporter(s3_dir, "S3", default_fields={"reaction_id": "RXN"}, run_id=run_id)
    assert json.loads(s3_manifest.read_text(encoding="utf-8"))[RUN_ID_FIELD] == run_id
    assert json.loads((s3_dir / "status.json").read_text(encoding="utf-8"))[RUN_ID_FIELD] == run_id

    s4_dir = tmp_path / "S4_HighLevel"
    s4_engine = HighLevelEngine(config, run_id=run_id)
    s4_manifest = s4_engine.run(
        [{**structures[0], "role": "product", "s1_thermochemistry_status": "complete"}],
        s4_dir,
    )
    S4ProgressReporter(
        s4_dir,
        [{**structures[0], "role": "product", "s1_thermochemistry_status": "complete"}],
        run_id=run_id,
        config=config,
    )
    assert json.loads(s4_manifest.read_text(encoding="utf-8"))[RUN_ID_FIELD] == run_id
    assert json.loads((s4_dir / "status.json").read_text(encoding="utf-8"))[RUN_ID_FIELD] == run_id


def test_cross_stage_gate_rejects_stale_s2_manifest(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    class FakeCensoLiteEngine:
        def __init__(self, _config: object, stage_dir: Path, molecule_name: str):
            self.stage_dir = Path(stage_dir)
            self.molecule_name = molecule_name

        def run(self, smiles: str) -> dict[str, object]:
            molecule_dir = self.stage_dir / self.molecule_name
            candidate_dir = molecule_dir / "candidates"
            candidate_dir.mkdir(parents=True, exist_ok=True)
            xyz = candidate_dir / "conf_0001.xyz"
            xyz.write_text(
                f"4\n{smiles}\nC 0 0 0\nC 1 0 0\nC 1 1 0\nC 0 1 0\n",
                encoding="utf-8",
            )
            (molecule_dir / "atom_mapping.json").write_text(
                json.dumps(
                    {
                        "schema_version": "s1_atom_mapping_v1",
                        "mapping_source": "rdkit_generation_sidecar",
                        "confidence": "high",
                        "product_smiles_idx_space": "geometry_product_smiles_idx",
                        "atoms": [
                            {"smiles_idx": i, "xyz_idx": i, "element": "C", "type": "organic_heavy"}
                            for i in range(4)
                        ],
                    }
                ),
                encoding="utf-8",
            )
            manifest = molecule_dir / "manifest.json"
            payload = {
                "selected": f"{self.molecule_name}_conf_0001",
                "selected_xyz": "selected.xyz",
                "atom_mapping_ref": "atom_mapping.json",
                "candidates": [
                    {
                        "id": f"{self.molecule_name}_conf_0001",
                        "xyz": f"{self.molecule_name}/candidates/conf_0001.xyz",
                    }
                ],
            }
            manifest.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            selected_xyz = molecule_dir / "selected.xyz"
            selected_xyz.write_text(xyz.read_text(encoding="utf-8"), encoding="utf-8")
            return {"manifest": manifest, "selected_xyz": selected_xyz, "data": payload}

    class FakePEBScanner:
        def __init__(self, _config: object, molecule_name: str | None = None):
            self.molecule_name = molecule_name or "unknown"
            self.event_callback = None

        def run(self, product_xyz: Path, output_dir: Path, forming_bonds, scan_config=None):
            del scan_config
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            ts_seed = output_dir / "ts_guess.xyz"
            intermediate_seed = output_dir / "intermediate_seed.xyz"
            profile = output_dir / "scan_profile.json"
            contents = Path(product_xyz).read_text(encoding="utf-8")
            ts_seed.write_text(contents, encoding="utf-8")
            intermediate_seed.write_text(contents, encoding="utf-8")
            profile.write_text("{}", encoding="utf-8")
            return ts_seed, intermediate_seed, intermediate_seed, tuple(forming_bonds), profile, "COMPLETE", "high", ()

    class FailIfCalledRefinementEngine:
        def __init__(self, *args, **kwargs):
            pass

        def set_progress_reporter(self, reporter):
            del reporter

        def set_event_callback(self, callback):
            del callback

        def run(self, structures, output_dir: Path) -> Path:
            del structures, output_dir
            raise AssertionError("S3 must not start after stale S2 detection")

    original_write_s2 = v4_orchestrator.V4Orchestrator._write_s2

    def stale_write_s2(*args, **kwargs):
        kwargs[RUN_ID_FIELD] = "55555555-5555-4555-8555-555555555555"
        return original_write_s2(*args, **kwargs)

    monkeypatch.setattr(v4_orchestrator, "CensoLiteEngine", FakeCensoLiteEngine)
    monkeypatch.setattr(v4_orchestrator, "PEBScanner", FakePEBScanner)
    monkeypatch.setattr(v4_orchestrator, "RefinementEngine", FailIfCalledRefinementEngine)
    monkeypatch.setattr(v4_orchestrator.V4Orchestrator, "_write_s2", staticmethod(stale_write_s2))

    orchestrator = v4_orchestrator.V4Orchestrator(color=False)
    # This test exercises the stale-run-id gate via the legacy PEB path; pin the
    # S2 method so the V4.1 default (neb_assisted) does not reroute the fakes.
    orchestrator.config.setdefault("step2", {})["method"] = "legacy_peb"
    with pytest.raises(RuntimeError, match=r"Stale S2 output detected\. Refusing to start S3\."):
        orchestrator.run(_record(tmp_path), tmp_path / "stale_s2", stop_after="s3")


def test_cross_stage_gate_allows_checkpoint_validated_reused_manifest(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "S2_PEB" / "manifest.json"
    manifest.parent.mkdir(parents=True)
    old_run_id = "77777777-7777-4777-8777-777777777777"
    current_run_id = "88888888-8888-4888-8888-888888888888"
    manifest.write_text(json.dumps({RUN_ID_FIELD: old_run_id}), encoding="utf-8")
    orchestrator = v4_orchestrator.V4Orchestrator(color=False)
    orchestrator.run_id = current_run_id

    orchestrator._assert_stage_results_run_id(
        {
            "product_major": {
                "status": "complete",
                "manifest": manifest,
                "reused": True,
            }
        },
        "S2",
        "S3",
    )

    with pytest.raises(RuntimeError, match=r"Stale S2 output detected\. Refusing to start S3\."):
        orchestrator._assert_stage_results_run_id(
            {
                "product_major": {
                    "status": "complete",
                    "manifest": manifest,
                    "reused": False,
                }
            },
            "S2",
            "S3",
        )


def test_cross_stage_gate_allows_reused_s3_manifest(tmp_path: Path) -> None:
    manifest = tmp_path / "S3_LowLevel" / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(
        json.dumps({RUN_ID_FIELD: "99999999-9999-4999-8999-999999999999"}),
        encoding="utf-8",
    )
    orchestrator = v4_orchestrator.V4Orchestrator(color=False)
    orchestrator.run_id = "aaaaaaaa-aaaa-4aaa-8aaa-aaaaaaaaaaaa"

    orchestrator._assert_manifest_run_id(
        manifest,
        "S3",
        "S4",
        allow_reused=True,
    )


def test_legacy_state_without_run_id_is_tolerated(
    caplog: pytest.LogCaptureFixture,
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "S2_PEB" / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("{}", encoding="utf-8")
    state_path = tmp_path / "pipeline.state"
    state_path.write_text(
        json.dumps(
            {
                "schema_version": STATE_SCHEMA_V2,
                "reaction_id": "RXN_LEGACY",
                "legacy_variant_checkpoint_imported": True,
                "stages": {
                    "s2": {
                        "signature": "legacy-s2",
                        "manifest": str(manifest.resolve()),
                        "status": "complete",
                    }
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    checkpoint = V4Checkpoint(tmp_path)
    state = checkpoint.load()

    assert state[RUN_ID_FIELD] is None
    assert checkpoint.run_id is None
    assert checkpoint.reusable("s2", "legacy-s2", manifest)

    new_run_id_value = "66666666-6666-4666-8666-666666666666"
    with caplog.at_level(logging.WARNING):
        checkpoint.set_run_id(new_run_id_value)

    saved = json.loads(state_path.read_text(encoding="utf-8"))
    assert saved[RUN_ID_FIELD] == new_run_id_value
    assert f"Resuming across runs (old run_id=None, new run_id={new_run_id_value})" in caplog.text
