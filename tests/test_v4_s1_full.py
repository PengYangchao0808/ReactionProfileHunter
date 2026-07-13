from __future__ import annotations

import json
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
from rph_core.steps.conformer_search.deduplicator import DedupCandidate
from rph_core.steps.mechanism_classifier.s0_record import S0ReactionRecord
from rph_core.steps.conformer_search.torsion_signature import TorsionSignature, build_signature, signatures_equivalent


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

        def run_sp(self, xyz_path: Path, nprocs=None) -> float:
            return -20.0 if xyz_path.stem.endswith("0001") else -19.0

        def run_mrrho(self, _xyz_path: Path) -> float | None:
            return None

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
    monkeypatch.setattr(engine, "deduplicator", FakeDeduplicator(), raising=False)

    result = engine.run("CCCC")

    manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
    selected_xyz = tmp_path / "S1_ConfSearch" / "product" / "selected.xyz"
    candidate_xyz = tmp_path / "S1_ConfSearch" / "product" / "candidates" / "conf_0001.xyz"
    assert manifest["selected_xyz"] == "selected.xyz"
    assert Path(result["selected_xyz"]) == selected_xyz
    assert selected_xyz.exists()
    assert selected_xyz.read_text(encoding="utf-8") == candidate_xyz.read_text(encoding="utf-8")


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
    assert runtime.run_mrrho(xyz) == pytest.approx(0.25)
    assert calls == [(1, 4), (0, 1)]
