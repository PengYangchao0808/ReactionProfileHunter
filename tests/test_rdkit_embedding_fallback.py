"""Tests for RDKit embedding robustness in ConformerEngine.

Validates that _step_rdkit_embed handles embedding failures gracefully
with proper fallback to random coordinates and UFF optimization.
"""
import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


class TestEmbeddingFallback:
    """Test that engine._step_rdkit_embed recovers from EmbedMolecule failures."""

    @staticmethod
    def _make_mock_mol():
        mol = MagicMock()
        mol.GetConformer.return_value.GetPositions.return_value = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ]
        atom = MagicMock()
        atom.GetSymbol.return_value = "C"
        mol.GetAtoms.return_value = [atom, atom]
        return mol

    def test_normal_embedding_succeeds(self):
        from rdkit import Chem
        from rph_core.steps.conformer_search.engine import ConformerEngine

        config = {
            "step1": {"conformer_search": {"two_stage_enabled": False}},
            "theory": {"preoptimization": {"gfn_level": 2}},
            "executables": {},
        }
        engine = ConformerEngine.__new__(ConformerEngine)
        engine.config = config
        engine.molecule_name = "test_mol"
        engine.molecule_dir = Path("/tmp/test_rdkit_embed")
        engine.crest_dir = engine.molecule_dir / "crest"
        engine.crest_dir.mkdir(parents=True, exist_ok=True)

        try:
            result = engine._step_rdkit_embed("C")
            assert result.exists()
            assert result.name == "test_mol_init.xyz"
        finally:
            import shutil
            shutil.rmtree(engine.molecule_dir, ignore_errors=True)

    @patch("rph_core.steps.conformer_search.engine.rdForceFieldHelpers")
    @patch("rph_core.steps.conformer_search.engine.rdDistGeom")
    @patch("rph_core.steps.conformer_search.engine.Chem")
    def test_etkdg_failure_triggers_random_coords_fallback(
        self, mock_chem, mock_dist_geom, mock_ff
    ):
        from rph_core.steps.conformer_search.engine import ConformerEngine
        import logging

        mock_mol = self._make_mock_mol()
        mock_chem.MolFromSmiles.return_value = mock_mol
        mock_chem.AddHs.return_value = mock_mol

        call_count = {"n": 0}

        def fake_embed(mol, params=None, **kwargs):
            call_count["n"] += 1
            if call_count["n"] == 1:
                return -1
            return 0

        mock_dist_geom.EmbedMolecule.side_effect = fake_embed
        mock_dist_geom.ETKDG.return_value = MagicMock()
        mock_dist_geom.ETKDGv3.return_value = MagicMock()
        mock_ff.MMFFOptimizeMolecule.return_value = 0

        config = {
            "step1": {"conformer_search": {"two_stage_enabled": False}},
            "theory": {"preoptimization": {"gfn_level": 2}},
            "executables": {},
        }
        engine = ConformerEngine.__new__(ConformerEngine)
        engine.config = config
        engine.molecule_name = "test_fallback"
        engine.molecule_dir = Path("/tmp/test_rdkit_fallback")
        engine.crest_dir = engine.molecule_dir / "crest"
        engine.crest_dir.mkdir(parents=True, exist_ok=True)
        try:
            engine._step_rdkit_embed("C")
            assert call_count["n"] == 2
            mock_dist_geom.ETKDGv3.assert_called_once()
        finally:
            import shutil
            shutil.rmtree(engine.molecule_dir, ignore_errors=True)

    @patch("rph_core.steps.conformer_search.engine.rdForceFieldHelpers")
    @patch("rph_core.steps.conformer_search.engine.rdDistGeom")
    @patch("rph_core.steps.conformer_search.engine.Chem")
    def test_all_embedding_fail_raises(self, mock_chem, mock_dist_geom, mock_ff):
        from rph_core.steps.conformer_search.engine import ConformerEngine
        import logging

        mock_mol = self._make_mock_mol()
        mock_chem.MolFromSmiles.return_value = mock_mol
        mock_chem.AddHs.return_value = mock_mol

        mock_dist_geom.EmbedMolecule.return_value = -1
        mock_dist_geom.ETKDG.return_value = MagicMock()
        mock_dist_geom.ETKDGv3.return_value = MagicMock()

        config = {
            "step1": {"conformer_search": {"two_stage_enabled": False}},
            "theory": {"preoptimization": {"gfn_level": 2}},
            "executables": {},
        }
        engine = ConformerEngine.__new__(ConformerEngine)
        engine.config = config
        engine.molecule_name = "test_total_fail"
        engine.molecule_dir = Path("/tmp/test_rdkit_total_fail")
        engine.crest_dir = engine.molecule_dir / "crest"
        engine.crest_dir.mkdir(parents=True, exist_ok=True)
        try:
            with pytest.raises(RuntimeError, match="embedding failed"):
                engine._step_rdkit_embed("C")
        finally:
            import shutil
            shutil.rmtree(engine.molecule_dir, ignore_errors=True)

    @patch("rph_core.steps.conformer_search.engine.rdForceFieldHelpers")
    @patch("rph_core.steps.conformer_search.engine.rdDistGeom")
    @patch("rph_core.steps.conformer_search.engine.Chem")
    def test_mmff_failure_falls_back_to_uff(self, mock_chem, mock_dist_geom, mock_ff):
        from rph_core.steps.conformer_search.engine import ConformerEngine
        import logging

        mock_mol = self._make_mock_mol()
        mock_chem.MolFromSmiles.return_value = mock_mol
        mock_chem.AddHs.return_value = mock_mol

        mock_dist_geom.EmbedMolecule.return_value = 0
        mock_dist_geom.ETKDG.return_value = MagicMock()
        mock_ff.MMFFOptimizeMolecule.side_effect = ValueError("MMFF failed")
        mock_ff.UFFOptimizeMolecule.return_value = 0

        config = {
            "step1": {"conformer_search": {"two_stage_enabled": False}},
            "theory": {"preoptimization": {"gfn_level": 2}},
            "executables": {},
        }
        engine = ConformerEngine.__new__(ConformerEngine)
        engine.config = config
        engine.molecule_name = "test_uff_fallback"
        engine.molecule_dir = Path("/tmp/test_rdkit_uff")
        engine.crest_dir = engine.molecule_dir / "crest"
        engine.crest_dir.mkdir(parents=True, exist_ok=True)
        try:
            engine._step_rdkit_embed("C")
            mock_ff.UFFOptimizeMolecule.assert_called_once()
        finally:
            import shutil
            shutil.rmtree(engine.molecule_dir, ignore_errors=True)

    @patch("rph_core.steps.conformer_search.engine.Chem")
    def test_invalid_smiles_raises(self, mock_chem):
        from rph_core.steps.conformer_search.engine import ConformerEngine
        import logging

        mock_chem.MolFromSmiles.return_value = None

        config = {
            "step1": {"conformer_search": {"two_stage_enabled": False}},
            "theory": {"preoptimization": {"gfn_level": 2}},
            "executables": {},
        }
        engine = ConformerEngine.__new__(ConformerEngine)
        engine.config = config
        engine.molecule_name = "test_invalid"
        engine.molecule_dir = Path("/tmp/test_rdkit_invalid")
        engine.crest_dir = engine.molecule_dir / "crest"
        engine.crest_dir.mkdir(parents=True, exist_ok=True)
        try:
            with pytest.raises(ValueError, match="Invalid SMILES"):
                engine._step_rdkit_embed("INVALID_SMILES_STRING")
        finally:
            import shutil
            shutil.rmtree(engine.molecule_dir, ignore_errors=True)
