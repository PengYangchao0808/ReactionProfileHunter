from __future__ import annotations

import json
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("rdkit")

from rdkit import Chem
from rdkit.Chem import rdDistGeom

import rph_core.v4_orchestrator as v4_orchestrator
from rph_core.steps.conformer_search import censo_lite as censo_lite_module
from rph_core.steps.conformer_search.censo_lite import CensoLiteEngine
from rph_core.steps.conformer_search.censo_lite_runtime import CensoLiteRuntime, CrestEnergyParseError
from rph_core.steps.conformer_search.deduplicator import (
    DedupCandidate,
    TorsionAwareDeduplicator,
)
from rph_core.steps.conformer_search.ensemble_thermo import calculate_ensemble_thermodynamics
from rph_core.steps.mechanism_classifier.s0_record import S0ReactionRecord
from rph_core.steps.conformer_search.torsion_signature import TorsionSignature, build_signature, signatures_equivalent
from rph_core.steps.conformer_search.xtb_thermo import XTBThermoResult, _xyz_to_coord
from rph_core.utils.constants import BOHR_TO_ANGSTROM
from rph_core.utils.orca_interface import ORCAInterface


def _embed_coordinates(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    mol = Chem.AddHs(mol)
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = 42
    assert rdDistGeom.EmbedMolecule(mol, params) == 0
    return mol, mol.GetConformer().GetPositions()


def test_censo_lite_engine_writes_selected_xyz(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
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
            return path

        def crest_search(self, input_xyz: Path) -> Path:
            return input_xyz

        def split_ensemble(self, _ensemble: Path) -> list[Path]:
            paths = []
            for index in (1, 2):
                path = self.raw_dir / f"conf_{index:04d}.xyz"
                path.write_text(f"1\nenergy = -10.{index}\nH {index}.0 0.0 0.0\n", encoding="utf-8")
                paths.append(path)
            return paths

        def extract_energy(self, xyz_path: Path) -> float:
            return -10.0 if xyz_path.stem.endswith("0001") else -9.0

        def read_crest_relative_energies(self, crest_dir: Path):
            return None

        def validate_energies(self, xyz_energies, rel_energies, tolerance_kcal=1.0):
            pass

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
        def annotate(self, _mol: Chem.Mol, path: Path, score: float, metadata: dict[str, object]) -> DedupCandidate:
            return DedupCandidate(
                path=Path(path),
                score=float(score),
                signature=TorsionSignature((), (), ()),
                metadata=dict(metadata),
            )

        def deduplicate(self, _mol: Chem.Mol, candidates: list[DedupCandidate]) -> list[DedupCandidate]:
            return list(candidates)

    monkeypatch.setattr(censo_lite_module, "CensoLiteRuntime", FakeRuntime)

    engine = CensoLiteEngine(
        {"step1": {"protocol": "censo_lite", "censo_lite": {"retention": {"energy_window_kcal": 10.0}}}},
        tmp_path / "S1_ConfSearch",
        "product",
    )
    ui_events: list[tuple[str, dict[str, object]]] = []
    engine.set_event_callback(lambda event, payload: ui_events.append((event, dict(payload))))
    monkeypatch.setattr(engine, "deduplicator", FakeDeduplicator(), raising=False)

    result = engine.run("CCCC")

    manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
    selected_xyz = tmp_path / "S1_ConfSearch" / "product" / "selected.xyz"
    candidate_xyz = tmp_path / "S1_ConfSearch" / "product" / "candidates" / "conf_0001.xyz"
    assert manifest["selected_xyz"] == "selected.xyz"
    assert Path(result["selected_xyz"]) == selected_xyz
    assert selected_xyz.exists()
    assert selected_xyz.read_text(encoding="utf-8") == candidate_xyz.read_text(encoding="utf-8")
    started_steps = [payload["step"] for event, payload in ui_events if event == "step_started"]
    finished_steps = [payload["step"] for event, payload in ui_events if event == "step_finished"]
    assert started_steps == [
        "rdkit_embed",
        "crest_search",
        "split_ensemble",
        "energy_validation",
        "torsion_dedup",
        "prefilter",
        "b97_sp",
        "mrrho",
        "ensemble_thermo",
        "selection_manifest",
    ]
    assert finished_steps == started_steps
    assert any(event == "batch_job_started" for event, _payload in ui_events)
    assert any(event == "batch_job_finished" for event, _payload in ui_events)
    assert manifest["schema_version"] == "s1_censo_light_ranking_v3"
    assert manifest["protocol"]["complete_censo_light"] is False
    assert manifest["protocol"]["requires_downstream_optimization"] is True
    assert manifest["thermodynamic_rank1"] == manifest["selected"]
    assert manifest["reactivity_screening_candidates"] == manifest["representative_candidates"]
    assert len(manifest["candidates"]) == 2
    assert manifest["ensemble_thermodynamics"]["ensemble_member_count"] == 2
    assert manifest["scheduling_summary"]["b97_3c"]["jobs"] == 2
    assert manifest["scheduling_summary"]["b97_3c"]["cores_per_job_counts"] == {
        "1": 2
    }
    assert manifest["scheduling_summary"]["xtb_mrrho"]["jobs"] == 2
    assert sum(row["boltzmann_population"] for row in manifest["candidates"]) == pytest.approx(1.0)
    first = manifest["candidates"][0]
    assert first["s1_score_hartree"] == pytest.approx(
        first["b973c_electronic_cpcm_hartree"]
        + first["xtb_mrrho_thermal_correction_hartree"]
    )
    assert first["xtb_energy_ledger_status"] == "validated"
    assert first["explicit_xtb_solvation_correction_added"] is False
    assert first["merge_count"] == 1
    assert first["degeneracy"] == 1


def test_v4_orchestrator_runs_precursor_s1_before_variants(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    record = S0ReactionRecord(
        rx_id="RXN_P3",
        product_smiles="C=C",
        reaction_type="4+3",
        mapped_product_smiles="[CH2:1]=[CH2:2]",
        mapped_forming_bonds=((1, 2), (1, 2)),
        forming_bonds=((0, 1), (0, 1)),
        mapping_confidence=None,
        mapping_trusted=True,
        source_csv=tmp_path / "dataset.csv",
        canonical_precursor_smiles="CC",
        mapped_precursor_smiles="[CH3:1][CH3:2]",
    )

    s1_calls: list[str] = []
    s3_structures: list[dict[str, object]] = []

    def fake_write_s0(work_dir: Path, record: S0ReactionRecord) -> Path:
        stage_dir = Path(work_dir) / "S0_Mechanism"
        stage_dir.mkdir(parents=True, exist_ok=True)
        manifest = stage_dir / "mechanism.json"
        manifest.write_text(
            json.dumps({
                "schema_version": "s0_mechanism_v2",
                "stage": "S0",
                "source": "trusted_reaction_record",
                "rx_id": record.rx_id,
                "forming_bonds": [list(pair) for pair in record.forming_bonds],
            }, indent=2),
            encoding="utf-8",
        )
        return manifest

    class FakeCensoLiteEngine:
        def __init__(self, _config: object, stage_dir: Path, molecule_name: str):
            self.stage_dir = Path(stage_dir)
            self.molecule_name = molecule_name

        def run(self, smiles: str) -> dict[str, object]:
            s1_calls.append(self.molecule_name)
            molecule_dir = self.stage_dir / self.molecule_name
            candidate_dir = molecule_dir / "candidates"
            candidate_dir.mkdir(parents=True, exist_ok=True)
            xyz = candidate_dir / "conf_0001.xyz"
            xyz.write_text(f"1\n{smiles}\nH 0.0 0.0 0.0\n", encoding="utf-8")
            manifest = molecule_dir / "manifest.json"
            payload = {
                "selected": f"{self.molecule_name}_conf_0001",
                "candidates": [
                    {
                        "id": f"{self.molecule_name}_conf_0001",
                        "xyz": f"{self.molecule_name}/candidates/conf_0001.xyz",
                    }
                ],
            }
            manifest.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            return {"manifest": manifest, "selected_xyz": xyz, "data": payload}

    class FakePEBScanner:
        def __init__(self, _config: object, molecule_name: str | None = None):
            self.molecule_name = molecule_name or "unknown"

        def run(self, product_xyz: Path, output_dir: Path, forming_bonds, scan_config=None):
            del forming_bonds, scan_config
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            ts_seed = output_dir / "ts_guess.xyz"
            intermediate_seed = output_dir / "intermediate_seed.xyz"
            profile = output_dir / "scan_profile.json"
            ts_seed.write_text(f"1\n{product_xyz}\nH 0.0 0.0 0.0\n", encoding="utf-8")
            intermediate_seed.write_text(f"1\n{product_xyz}\nH 1.0 0.0 0.0\n", encoding="utf-8")
            profile.write_text("{}", encoding="utf-8")
            return ts_seed, intermediate_seed, intermediate_seed, ((0, 1), (0, 1)), profile, "COMPLETE", "high", ()

    class FakeLowLevelEngine:
        def __init__(self, _config: object):
            pass

        def run(self, structures, output_dir: Path) -> Path:
            rows = list(structures)
            s3_structures.extend(rows)
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            manifest = output_dir / "manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "schema_version": "s3_low_level_v1",
                        "stage": "S3",
                        "structures": rows,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            return manifest

    monkeypatch.setattr(v4_orchestrator.V4Orchestrator, "_write_s0_from_record", staticmethod(fake_write_s0))
    monkeypatch.setattr(v4_orchestrator, "CensoLiteEngine", FakeCensoLiteEngine)
    monkeypatch.setattr(v4_orchestrator, "PEBScanner", FakePEBScanner)
    monkeypatch.setattr(v4_orchestrator, "LowLevelEngine", FakeLowLevelEngine)

    orchestrator = object.__new__(v4_orchestrator.V4Orchestrator)
    orchestrator.config = {
        "step1": {"protocol": "censo_lite"},
        "step2": {},
        "theory": {"s3_low_level": {}, "s4_high_precision": {}},
    }

    result = orchestrator.run(record, tmp_path / "run", stop_after="s3")

    assert result["completed_through"] == "s3"
    assert s1_calls == ["precursor", "product_major"]
    assert (tmp_path / "run" / "S1_ConfSearch" / "precursor" / "manifest.json").exists()
    assert [row["id"] for row in s3_structures] == ["precursor", "product_major", "product_major_int", "product_major_ts"]


def test_atom_map_torsion_signatures_use_map_numbers_in_keys():
    mol, coordinates = _embed_coordinates("[CH3:11][CH2:22][CH2:33][CH3:44]")

    signature = build_signature(mol, coordinates)

    assert signature.bond_space == "atom_map"
    assert "11-22" in signature.key()
    assert "22-33" in signature.key()
    assert "33-44" in signature.key()
    assert "0-1" not in signature.key()


def test_torsion_signatures_with_different_atom_maps_are_not_equivalent():
    mol_a, coordinates_a = _embed_coordinates("[CH3:11][CH2:22][CH2:33][CH3:44]")
    mol_b, coordinates_b = _embed_coordinates("[CH3:101][CH2:202][CH2:303][CH3:404]")

    signature_a = build_signature(mol_a, coordinates_a)
    signature_b = build_signature(mol_b, coordinates_b)

    assert signature_a.bond_space == signature_b.bond_space == "atom_map"
    assert signature_a.key() != signature_b.key()
    assert signatures_equivalent(signature_a, signature_b) is False


def test_extract_energy_raises_on_parse_failure(tmp_path: Path):
    xyz = tmp_path / "broken.xyz"
    xyz.write_text("1\nnot_an_energy_comment\nH 0.0 0.0 0.0\n", encoding="utf-8")

    with pytest.raises(CrestEnergyParseError, match="Could not parse xTB energy"):
        CensoLiteRuntime.extract_energy(xyz)


def test_extract_energy_accepts_crest_bare_numeric_comment(tmp_path: Path):
    xyz = tmp_path / "crest.xyz"
    xyz.write_text("1\n-62.10299138\nH 0.0 0.0 0.0\n", encoding="utf-8")

    assert CensoLiteRuntime.extract_energy(xyz) == pytest.approx(-62.10299138)


def test_xyz_to_xtb_coord_converts_angstrom_to_bohr(tmp_path: Path):
    xyz = tmp_path / "one_angstrom.xyz"
    xyz.write_text("1\ngeometry in angstrom\nH 1.0 0.0 0.0\n", encoding="utf-8")
    coord = tmp_path / "coord"

    _xyz_to_coord(xyz, coord)

    atom_fields = coord.read_text(encoding="utf-8").splitlines()[1].split()
    assert float(atom_fields[0]) == pytest.approx(1.0 / BOHR_TO_ANGSTROM)


def test_mrrho_retries_crash_with_configured_gfn_and_nproc(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    calls = []

    class FakeXTBInterface:
        def __init__(self, *, gfn_level, solvent, nproc, config):
            calls.append((gfn_level, nproc))

        def enso_thermo(self, **kwargs):
            if calls[-1] == (1, 4):
                return SimpleNamespace(
                    success=False,
                    g_rrho_correction_hartree=0.0,
                    error="xTB SPH+MRRHO failed (rc=174): SIGSEGV",
                )
            return SimpleNamespace(success=True, g_rrho_correction_hartree=0.25, error=None)

    monkeypatch.setattr(
        "rph_core.steps.conformer_search.censo_lite_runtime.XTBInterface",
        FakeXTBInterface,
    )
    xyz = tmp_path / "conf.xyz"
    xyz.write_text("1\n-1.0\nH 0.0 0.0 0.0\n", encoding="utf-8")
    config = {
        "resources": {"nproc": 16},
        "step1": {
            "censo_lite": {
                "xtb_thermo": {
                    "gfn_level": 1,
                    "nproc": 4,
                    "fallback_enabled": True,
                    "fallback_gfn_level": 0,
                    "fallback_nproc": 1,
                }
            }
        },
    }

    runtime = CensoLiteRuntime(config, tmp_path / "work", "molecule")
    result = runtime.run_mrrho(xyz)
    assert result is not None
    assert result.g_rrho_correction_hartree == pytest.approx(0.25)
    assert calls == [(1, 4), (1, 1)]


def test_ensemble_partition_function_and_boltzmann_populations():
    # Two equally populated conformers give Z_rel=2 and G_conf=-RT ln(2).
    result = calculate_ensemble_thermodynamics([-100.0, -100.0], 298.15)

    assert result.partition_function_relative == pytest.approx(2.0)
    assert result.boltzmann_populations == pytest.approx((0.5, 0.5))
    assert result.conformational_free_energy_correction_kcal == pytest.approx(
        -0.00198720425864083 * 298.15 * 0.6931471805599453
    )
    assert result.ensemble_free_energy_hartree < result.reference_free_energy_hartree


def test_ensemble_thermodynamics_physical_invariants():
    single = calculate_ensemble_thermodynamics([-100.0], 298.15)
    assert single.partition_function_relative == pytest.approx(1.0)
    assert single.conformational_free_energy_correction_kcal == pytest.approx(0.0)
    assert single.relative_free_energies_kcal == pytest.approx((0.0,))
    assert sum(single.boltzmann_populations) == pytest.approx(1.0)

    base = calculate_ensemble_thermodynamics([-100.0, -99.999], 298.15)
    shifted = calculate_ensemble_thermodynamics([-50.0, -49.999], 298.15)
    assert shifted.relative_free_energies_kcal == pytest.approx(
        base.relative_free_energies_kcal
    )
    assert shifted.boltzmann_populations == pytest.approx(base.boltzmann_populations)
    assert shifted.conformational_free_energy_correction_kcal == pytest.approx(
        base.conformational_free_energy_correction_kcal
    )

    degenerate = calculate_ensemble_thermodynamics([-100.0], 298.15, [2])
    assert degenerate.partition_function_relative == pytest.approx(2.0)
    assert degenerate.conformational_free_energy_correction_kcal == pytest.approx(
        -0.00198720425864083 * 298.15 * 0.6931471805599453
    )
    with pytest.raises(ValueError, match="finite integer"):
        calculate_ensemble_thermodynamics([-100.0], 298.15, [1.5])


def test_dedup_tracks_merge_provenance_without_inventing_degeneracy(tmp_path: Path):
    mol, coordinates = _embed_coordinates("CC")
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    left = tmp_path / "crest_conf_001.xyz"
    right = tmp_path / "crest_conf_012.xyz"
    from rph_core.utils.file_io import write_xyz

    write_xyz(left, coordinates, symbols, title="duplicate 1")
    write_xyz(right, coordinates, symbols, title="duplicate 2")
    deduplicator = TorsionAwareDeduplicator(
        {"heavy_atom_rmsd_prefilter_A": 0.25, "torsion_rmsd_deg": 25.0}
    )
    records = [
        deduplicator.annotate(mol, left, -10.0),
        deduplicator.annotate(mol, right, -9.9),
    ]

    unique = deduplicator.deduplicate(mol, records)

    assert len(unique) == 1
    assert unique[0].metadata["merge_count"] == 2
    assert unique[0].metadata["merged_from"] == [
        "crest_conf_001",
        "crest_conf_012",
    ]
    assert unique[0].metadata["degeneracy"] == 1
    assert unique[0].metadata["degeneracy_source"] == "default_unique_minimum"


def test_mrrho_defer_preserves_candidates_but_suppresses_partition_function(
    tmp_path: Path,
):
    engine = CensoLiteEngine(
        {
            "resources": {"nproc": 1},
            "step1": {"censo_lite": {"xtb_thermo": {"failure_policy": "defer"}}},
        },
        tmp_path,
        "product",
    )
    signature = TorsionSignature((), (), ())
    records = []
    for index, energy in enumerate((-100.0, -99.9), start=1):
        path = tmp_path / f"conf_{index:04d}.xyz"
        path.write_text("1\ntest\nH 0.0 0.0 0.0\n", encoding="utf-8")
        records.append(
            DedupCandidate(
                path,
                energy,
                signature,
                {
                    "source_conformer_id": path.stem,
                    "b973c_sp_energy_hartree": energy,
                    "b973c_electronic_cpcm_hartree": energy,
                },
            )
        )

    def fake_mrrho(path: Path, nprocs=None):
        del nprocs
        if path.stem.endswith("0002"):
            return XTBThermoResult(
                g_rrho_correction_hartree=None,
                success=False,
                error="SCC did not converge",
            )
        return XTBThermoResult(
            g_rrho_correction_hartree=0.2,
            xtb_electronic_energy_hartree=-10.0,
            xtb_total_free_energy_hartree=-9.8,
            energy_ledger_residual_hartree=0.0,
            gfn_level=2,
            solvent_model="ALPB(acetone)",
        )

    engine.engine.run_mrrho = fake_mrrho
    deferred, scoring_mode, status = engine._apply_mrrho_uniform(
        records, {"failure_policy": "defer", "gfn_level": 2}
    )
    annotated, thermo = engine._annotate_ensemble_thermodynamics(deferred, status)

    assert scoring_mode == "electronic_only_deferred"
    assert status == "incomplete"
    assert [row.score for row in deferred] == pytest.approx([-100.0, -99.9])
    assert thermo["partition_function_relative"] is None
    assert thermo["conformational_free_energy_correction_kcal"] is None
    assert thermo["failed_conformer_ids"] == ["conf_0002"]
    assert all(row.metadata["boltzmann_population"] is None for row in annotated)
    assert all(row.metadata["relative_free_energy_kcal"] is None for row in annotated)

    result = engine._write_manifest(
        annotated,
        annotated[:1],
        engine.engine.molecule_dir,
        [],
        scoring_mode,
        thermo,
    )
    manifest = result["data"]
    assert manifest["thermodynamic_rank1"] is None
    assert manifest["provisional_geometry_selection"] == manifest["selected"]
    assert manifest["ensemble_thermodynamics"]["thermochemistry_complete"] is False
    assert manifest["ensemble_thermodynamics"][
        "ensemble_thermochemistry_correction_hartree"
    ] is None


def test_mrrho_strict_rejects_partial_thermochemistry(tmp_path: Path):
    engine = CensoLiteEngine(
        {
            "resources": {"nproc": 1},
            "step1": {"censo_lite": {"xtb_thermo": {"failure_policy": "strict"}}},
        },
        tmp_path,
        "product",
    )
    path = tmp_path / "conf_0001.xyz"
    path.write_text("1\ntest\nH 0.0 0.0 0.0\n", encoding="utf-8")
    record = DedupCandidate(
        path,
        -100.0,
        TorsionSignature((), (), ()),
        {"b973c_sp_energy_hartree": -100.0},
    )
    engine.engine.run_mrrho = lambda *_args, **_kwargs: XTBThermoResult(
        g_rrho_correction_hartree=None,
        success=False,
        error="SCC did not converge",
    )

    with pytest.raises(RuntimeError, match="failure_policy is 'strict'"):
        engine._apply_mrrho_uniform([record], {"failure_policy": "strict"})


def test_parallel_layout_respects_global_core_budget(tmp_path: Path):
    engine = CensoLiteEngine(
        {
            "resources": {"nproc": 16},
            "step1": {"censo_lite": {}},
        },
        tmp_path,
        "product",
    )

    assert engine._resolve_parallel_layout(
        {"parallel_jobs": "auto", "cores_per_job": 4}, 1
    ) == (4, 4)
    assert engine._resolve_parallel_layout(
        {"parallel_jobs": 99, "cores_per_job": 3}, 1
    ) == (5, 3)


def test_orca_runtime_environment_prevents_hidden_thread_oversubscription(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setenv("OMP_NUM_THREADS", "64")
    monkeypatch.setenv("MKL_NUM_THREADS", "64")
    interface = ORCAInterface.__new__(ORCAInterface)
    interface.config = {"resources": {"orca_math_threads_per_rank": 1}}
    interface.orca_binary = None

    env = interface._build_orca_runtime_env()

    assert env["OMP_NUM_THREADS"] == "1"
    assert env["MKL_NUM_THREADS"] == "1"
    assert env["OPENBLAS_NUM_THREADS"] == "1"
    assert env["OMP_THREAD_LIMIT"] == "1"
    assert env["OMP_MAX_ACTIVE_LEVELS"] == "1"


def test_b97_legacy_tail_is_ignored_in_throughput_mode(tmp_path: Path):
    engine = CensoLiteEngine(
        {
            "resources": {"nproc": 16},
            "step1": {
                "censo_lite": {
                    "b97_3c": {
                        "parallel_jobs": 16,
                        "cores_per_job": 1,
                        "max_parallel_jobs": 16,
                        "adaptive_tail": {
                            "enabled": True,
                            "jobs": 8,
                            "min_cores_per_job": 2,
                            "max_cores_per_job": 16,
                        },
                    }
                }
            },
        },
        tmp_path,
        "product",
    )
    signature = TorsionSignature((), (), ())
    records = []
    for index in range(20):
        path = tmp_path / f"conf_{index:04d}.xyz"
        path.write_text("1\n-1.0\nH 0.0 0.0 0.0\n", encoding="utf-8")
        records.append(DedupCandidate(path, -1.0, signature, {}))

    lock = threading.Lock()
    active_cores = 0
    peak_cores = 0
    assigned: dict[str, int] = {}

    def fake_run_sp(path: Path, nprocs=None, maxcore=None):
        nonlocal active_cores, peak_cores
        del maxcore
        cores = int(nprocs or 1)
        with lock:
            active_cores += cores
            peak_cores = max(peak_cores, active_cores)
            assigned[path.stem] = cores
        time.sleep(0.05)
        with lock:
            active_cores -= cores
        return -100.0

    engine.engine.run_sp = fake_run_sp
    completed = engine._run_b97_sp_parallel(Chem.MolFromSmiles("C"), records)

    assert len(completed) == 20
    assert peak_cores == 16
    assert peak_cores <= 16
    assert all(value == 1 for value in assigned.values())


def test_b97_profiled_final_wave_uses_exact_benchmarked_layout(tmp_path: Path):
    engine = CensoLiteEngine(
        {
            "resources": {"nproc": 16},
            "step1": {
                "censo_lite": {
                    "b97_3c": {
                        "parallel_jobs": 16,
                        "cores_per_job": 1,
                        "scheduling": {
                            "final_wave": {
                            "enabled": True,
                                "max_cores_per_job": 4,
                                "profile_layouts": {
                                    "6": [3, 3, 3, 3, 2, 2],
                                },
                            },
                        },
                    }
                }
            },
        },
        tmp_path,
        "product",
    )
    signature = TorsionSignature((), (), ())
    records = []
    for index in range(6):
        path = tmp_path / f"small_{index:04d}.xyz"
        path.write_text("1\n-1.0\nH 0.0 0.0 0.0\n", encoding="utf-8")
        records.append(DedupCandidate(path, -1.0, signature, {}))

    lock = threading.Lock()
    active_cores = 0
    peak_cores = 0
    assigned: dict[str, int] = {}

    def fake_run_sp(path: Path, nprocs=None, maxcore=None):
        nonlocal active_cores, peak_cores
        del maxcore
        cores = int(nprocs or 1)
        with lock:
            active_cores += cores
            peak_cores = max(peak_cores, active_cores)
            assigned[path.stem] = cores
        time.sleep(0.05)
        with lock:
            active_cores -= cores
        return -100.0

    engine.engine.run_sp = fake_run_sp
    engine._run_b97_sp_parallel(Chem.MolFromSmiles("C"), records)

    assert peak_cores == 16
    assert sum(assigned.values()) == 16
    assert sorted(assigned.values()) == [2, 2, 3, 3, 3, 3]


def test_b97_small_batch_without_profile_keeps_one_core_per_job(tmp_path: Path):
    engine = CensoLiteEngine(
        {
            "resources": {"nproc": 16, "mem": "32GB"},
            "step1": {
                "censo_lite": {
                    "b97_3c": {
                        "parallel_jobs": 16,
                        "cores_per_job": 1,
                        "scheduling": {
                            "final_wave": {
                                "enabled": True,
                                "profile_layouts": {},
                            }
                        },
                    }
                }
            },
        },
        tmp_path,
        "product",
    )
    signature = TorsionSignature((), (), ())
    records = []
    for index in range(12):
        path = tmp_path / f"safe_{index:04d}.xyz"
        path.write_text("1\n-1.0\nH 0.0 0.0 0.0\n", encoding="utf-8")
        records.append(DedupCandidate(path, -1.0, signature, {}))

    lock = threading.Lock()
    active_cores = 0
    peak_cores = 0
    assigned: dict[str, tuple[int, int | None]] = {}

    def fake_run_sp(path: Path, nprocs=None, maxcore=None):
        nonlocal active_cores, peak_cores
        cores = int(nprocs or 1)
        with lock:
            active_cores += cores
            peak_cores = max(peak_cores, active_cores)
            assigned[path.stem] = (cores, maxcore)
        time.sleep(0.05)
        with lock:
            active_cores -= cores
        return -100.0

    engine.engine.run_sp = fake_run_sp
    engine._run_b97_sp_parallel(Chem.MolFromSmiles("C"), records)

    assert peak_cores == 12
    assert all(cores == 1 for cores, _maxcore in assigned.values())
    assert all(maxcore == 1518 for _cores, maxcore in assigned.values())


def test_b97_final_wave_waits_for_bulk_barrier(tmp_path: Path):
    engine = CensoLiteEngine(
        {
            "resources": {"nproc": 16, "mem": "32GB"},
            "step1": {
                "censo_lite": {
                    "b97_3c": {
                        "parallel_jobs": 16,
                        "cores_per_job": 1,
                        "scheduling": {
                            "final_wave": {
                                "enabled": True,
                                "max_cores_per_job": 4,
                                "profile_layouts": {"4": [4, 4, 4, 4]},
                            }
                        },
                    }
                }
            },
        },
        tmp_path,
        "product",
    )
    signature = TorsionSignature((), (), ())
    records = []
    for index in range(20):
        path = tmp_path / f"barrier_{index:04d}.xyz"
        path.write_text("1\n-1.0\nH 0.0 0.0 0.0\n", encoding="utf-8")
        records.append(DedupCandidate(path, -1.0, signature, {}))

    lock = threading.Lock()
    active_bulk = 0
    final_started_during_bulk = False
    assigned: dict[str, int] = {}

    def fake_run_sp(path: Path, nprocs=None, maxcore=None):
        nonlocal active_bulk, final_started_during_bulk
        del maxcore
        index = int(path.stem.rsplit("_", 1)[-1])
        cores = int(nprocs or 1)
        with lock:
            if index < 16:
                active_bulk += 1
            elif active_bulk:
                final_started_during_bulk = True
            assigned[path.stem] = cores
        time.sleep(0.03)
        with lock:
            if index < 16:
                active_bulk -= 1
        return -100.0

    engine.engine.run_sp = fake_run_sp
    engine._run_b97_sp_parallel(Chem.MolFromSmiles("C"), records)

    assert final_started_during_bulk is False
    assert all(assigned[f"barrier_{index:04d}"] == 1 for index in range(16))
    assert all(assigned[f"barrier_{index:04d}"] == 4 for index in range(16, 20))


def test_b97_memory_budget_reduces_parallel_job_count(tmp_path: Path):
    engine = CensoLiteEngine(
        {
            "resources": {
                "nproc": 16,
                "mem": "8GB",
                "orca_maxcore_safety": 0.65,
            },
            "step1": {"censo_lite": {}},
        },
        tmp_path,
        "product",
    )

    jobs = engine._limit_orca_jobs_by_memory(
        {
            "scheduling": {
                "memory": {
                    "job_overhead_mb": 256,
                    "min_maxcore_mb_per_rank": 1000,
                }
            }
        },
        requested_jobs=16,
        cores_per_job=1,
    )

    assert jobs == 4


def test_orca_task_metrics_are_parsed_for_manifest(tmp_path: Path):
    output_dir = tmp_path / "ranking" / "conf_0001"
    output_dir.mkdir(parents=True)
    output = output_dir / "conf_0001.out"
    output.write_text(
        "* SCF CONVERGED AFTER  15 CYCLES *\n"
        "TOTAL RUN TIME: 0 days 0 hours 2 minutes 6 seconds 327 msec\n",
        encoding="utf-8",
    )

    metrics = CensoLiteEngine._read_orca_task_metrics(output_dir)

    assert metrics["sp_scf_cycles"] == 15
    assert metrics["sp_orca_total_run_seconds"] == pytest.approx(126.327)
    assert metrics["sp_output"] == str(output)


def test_mrrho_profiled_final_wave_uses_core_budget_without_oversubscription(tmp_path: Path):
    thermo_cfg = {
        "failure_policy": "strict",
        "gfn_level": 2,
        "parallel_jobs": 16,
        "cores_per_job": 1,
        "max_parallel_jobs": 16,
        "scheduling": {
            "final_wave": {
                "enabled": True,
                "max_cores_per_job": 4,
                "profile_layouts": {"6": [3, 3, 3, 3, 2, 2]},
            }
        },
    }
    engine = CensoLiteEngine(
        {
            "resources": {"nproc": 16},
            "step1": {"censo_lite": {"xtb_thermo": thermo_cfg}},
        },
        tmp_path,
        "product",
    )
    signature = TorsionSignature((), (), ())
    records = []
    for index in range(6):
        path = tmp_path / f"thermo_{index:04d}.xyz"
        path.write_text("1\n-1.0\nH 0.0 0.0 0.0\n", encoding="utf-8")
        records.append(
            DedupCandidate(
                path,
                -100.0 + index * 0.001,
                signature,
                {"b973c_sp_energy_hartree": -100.0 + index * 0.001},
            )
        )

    lock = threading.Lock()
    active_cores = 0
    peak_cores = 0
    assigned: dict[str, int] = {}

    def fake_run_mrrho(path: Path, nprocs=None):
        nonlocal active_cores, peak_cores
        cores = int(nprocs or 1)
        with lock:
            active_cores += cores
            peak_cores = max(peak_cores, active_cores)
            assigned[path.stem] = cores
        time.sleep(0.05)
        with lock:
            active_cores -= cores
        return 0.1

    engine.engine.run_mrrho = fake_run_mrrho
    corrected, scoring_mode, status = engine._apply_mrrho_uniform(records, thermo_cfg)

    assert len(corrected) == 6
    assert scoring_mode == "b97_3c_plus_mrrho"
    assert status == "complete"
    assert peak_cores == 16
    assert sum(assigned.values()) == 16
    assert sorted(assigned.values()) == [2, 2, 3, 3, 3, 3]


def test_crest_nproc_is_clamped_to_global_budget(monkeypatch, tmp_path: Path):
    captured = {}

    class FakeCRESTInterface:
        def __init__(self, *, nproc, **kwargs):
            del kwargs
            captured["nproc"] = nproc
            self.crest_timeout_seconds = None

        def run_conformer_search(self, _input, output_dir, **kwargs):
            del kwargs
            output = Path(output_dir) / "crest_conformers.xyz"
            output.write_text("1\n-1.0\nH 0.0 0.0 0.0\n", encoding="utf-8")
            return output

    monkeypatch.setattr(
        "rph_core.steps.conformer_search.censo_lite_runtime.CRESTInterface",
        FakeCRESTInterface,
    )
    runtime = CensoLiteRuntime(
        {
            "resources": {"nproc": 16},
            "step1": {"censo_lite": {"crest": {"nproc": 64}}},
        },
        tmp_path,
        "product",
    )
    initial = runtime.molecule_dir / "initial.xyz"
    initial.write_text("1\ntest\nH 0.0 0.0 0.0\n", encoding="utf-8")

    runtime.crest_search(initial)

    assert captured["nproc"] == 16
