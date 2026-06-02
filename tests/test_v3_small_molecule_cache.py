from __future__ import annotations

import json
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

from rph_core.utils.small_molecule_cache import SmallMoleculeCache
from rph_core.scheduling.small_molecule_precompute import SmallMoleculePrecomputer


class TestSmallMoleculeCacheIsComplete:
    def test_incomplete_when_no_entry(self, tmp_path: Path):
        cache = SmallMoleculeCache(tmp_path / "cache")
        assert not cache.is_complete("CC(=O)O")

    def test_incomplete_when_no_molecule_xyz(self, tmp_path: Path):
        cache = SmallMoleculeCache(tmp_path / "cache")
        entry = cache.get_or_create("CC(=O)O", name="AcOH")
        (entry / "thermo.json").write_text("{}")
        assert not cache.is_complete("CC(=O)O", require_thermo=False)
        assert not cache.is_complete("CC(=O)O", require_thermo=True)

    def test_complete_with_xyz_only(self, tmp_path: Path):
        cache = SmallMoleculeCache(tmp_path / "cache")
        entry = cache.get_or_create("CC(=O)O", name="AcOH")
        (entry / "molecule_min.xyz").write_text("2\n\nC 0 0 0\nO 1 0 0\n")
        assert cache.is_complete("CC(=O)O", require_thermo=False)
        assert not cache.is_complete("CC(=O)O", require_thermo=True)

    def test_complete_with_xyz_and_thermo(self, tmp_path: Path):
        cache = SmallMoleculeCache(tmp_path / "cache")
        entry = cache.get_or_create("CC(=O)O", name="AcOH")
        (entry / "molecule_min.xyz").write_text("2\n\nC 0 0 0\nO 1 0 0\n")
        (entry / "thermo.json").write_text('{"g_kcal": -100.0}')
        assert cache.is_complete("CC(=O)O", require_thermo=True)

    def test_theory_signature_mismatch(self, tmp_path: Path):
        cache = SmallMoleculeCache(tmp_path / "cache")
        entry = cache.get_or_create("CC(=O)O", name="AcOH")
        (entry / "molecule_min.xyz").write_text("2\n\nC 0 0 0\nO 1 0 0\n")
        (entry / "thermo.json").write_text("{}")
        cache.write_cache_meta("CC(=O)O", {"method": "B3LYP", "basis": "def2-SVP"})
        assert not cache.is_complete("CC(=O)O", theory_signature={"method": "M062X", "basis": "def2-SVP"})


class TestSmallMoleculeCacheGetEntryManifest:
    def test_missing_entry(self, tmp_path: Path):
        cache = SmallMoleculeCache(tmp_path / "cache")
        assert cache.get_entry_manifest("nonexistent") == {}

    def test_existing_manifest(self, tmp_path: Path):
        cache = SmallMoleculeCache(tmp_path / "cache")
        cache.get_or_create("CC(=O)O", name="AcOH")
        cache.write_cache_meta("CC(=O)O", {"method": "B3LYP"}, extra_meta={"custom": "value"})
        manifest = cache.get_entry_manifest("CC(=O)O")
        assert manifest["theory_signature"]["method"] == "B3LYP"
        assert manifest["custom"] == "value"


class TestSmallMoleculePrecomputerSignature:
    def test_recomputes_when_theory_signature_mismatches(self, tmp_path: Path, monkeypatch):
        config = {
            "reference_states": {"small_molecule_map": {"AcOH": {"smiles": "CC(=O)O"}}},
            "theory": {"optimization": {"method": "M062X", "basis": "def2-SVP", "solvent": "acetone", "engine": "gaussian"}},
            "step1": {"protocol": "lite"},
        }
        precomputer = SmallMoleculePrecomputer(config=config, cache_root=tmp_path / "cache")
        cache = precomputer.cache
        entry = cache.get_or_create("CC(=O)O", name="AcOH")
        (entry / "molecule_min.xyz").write_text("2\n\nC 0 0 0\nO 1 0 0\n")
        (entry / "thermo.json").write_text("{}")
        cache.write_cache_meta("CC(=O)O", {"method": "B3LYP", "basis": "def2-SVP", "solvent": "acetone", "engine": "gaussian"})

        called = {"value": False}

        def fake_compute(key: str, smiles: str, cache_dir: Path) -> None:
            called["value"] = True
            (cache_dir / "molecule_min.xyz").write_text("2\n\nC 0 0 0\nO 1 0 0\n")
            cache.write_cache_meta(smiles, precomputer._build_theory_signature())

        monkeypatch.setattr(precomputer, "_compute_molecule", fake_compute)

        precomputer.ensure_one("AcOH")

        assert called["value"] is True
