import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
import sys

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from rph_core.orchestrator import ReactionProfileHunter, PipelineResult
from rph_core.steps.anchor.handler import AnchorPhaseResult
from rph_core.utils.task_builder import TaskSpec

class TestOrchestratorMultiMolecule:
    @pytest.fixture
    def config(self):
        return {
            "global": {"log_level": "DEBUG"},
            "run": {"resume": False, "output_root": "./rph_output", "workdir_naming": "rx_{rx_id}"},
            "reference_states": {
                "small_molecule_map": {
                    "AcOH": {"smiles": "CC(=O)O", "charge": 0, "multiplicity": 1},
                    "MeOH": {"smiles": "CO", "charge": 0, "multiplicity": 1}
                }
            }
        }

    @patch("rph_core.orchestrator.load_config")
    @patch("rph_core.orchestrator.normalize_qc_config")
    @patch("rph_core.orchestrator.ui.print_pipeline_header")
    def test_run_pipeline_multi_molecule(self, mock_ui, mock_normalize, mock_load, config, tmp_path):
        mock_load.return_value = config
        mock_normalize.return_value = (config, [])
        
        hunter = ReactionProfileHunter()
        
        # Mock S1 engine
        mock_s1 = MagicMock()
        mock_s1.run.return_value = AnchorPhaseResult(
            success=True,
            anchored_molecules={
                "product": {"xyz": tmp_path / "S1_ConfGeneration" / "prod.log", "e_sp": -100.0}
            }
        )
        
        # Use property mock
        with patch.object(ReactionProfileHunter, 's1_engine', mock_s1):
            # Create a fake log file to satisfy LogParser
            prod_log = tmp_path / "S1_ConfGeneration" / "prod.log"
            prod_log.parent.mkdir(parents=True, exist_ok=True)
            prod_log.write_text("Dummy Gaussian log")
            
            # Mock LogParser to return dummy coords
            with patch("rph_core.utils.geometry_tools.LogParser.extract_last_converged_coords") as mock_extract:
                mock_extract.return_value = ([[0, 0, 0]], ["C"], None)
                
                # Mock write_xyz to do nothing
                with patch("rph_core.utils.file_io.write_xyz"):
                    result = hunter.run_pipeline(
                        product_smiles="C1CCCCC1",
                        work_dir=tmp_path,
                        skip_steps=['s2', 's3', 's4'],
                        precursor_smiles="C=C",
                        leaving_group_key="AcOH"
                    )
        
        assert result.success
        # Verify molecules passed to S1
        expected_molecules = {
            "product": "C1CCCCC1",
            "precursor": "C=C",
            "leaving_group": "CC(=O)O"
        }
        mock_s1.run.assert_called_once_with(molecules=expected_molecules)

    @patch("rph_core.orchestrator.load_config")
    @patch("rph_core.orchestrator.normalize_qc_config")
    @patch("rph_core.orchestrator.ui.print_pipeline_header")
    def test_run_pipeline_missing_lg(self, mock_ui, mock_normalize, mock_load, config, tmp_path):
        mock_load.return_value = config
        mock_normalize.return_value = (config, [])
        
        hunter = ReactionProfileHunter()
        
        mock_s1 = MagicMock()
        mock_s1.run.return_value = AnchorPhaseResult(
            success=True,
            anchored_molecules={
                "product": {"xyz": tmp_path / "S1_ConfGeneration" / "prod.log", "e_sp": -100.0}
            }
        )
        
        with patch.object(ReactionProfileHunter, 's1_engine', mock_s1):
            prod_log = tmp_path / "S1_ConfGeneration" / "prod.log"
            prod_log.parent.mkdir(parents=True, exist_ok=True)
            prod_log.write_text("Dummy Gaussian log")
            
            with patch("rph_core.utils.geometry_tools.LogParser.extract_last_converged_coords") as mock_extract:
                mock_extract.return_value = ([[0, 0, 0]], ["C"], None)
                with patch("rph_core.utils.file_io.write_xyz"):
                    result = hunter.run_pipeline(
                        product_smiles="C1CCCCC1",
                        work_dir=tmp_path,
                        skip_steps=['s2', 's3', 's4'],
                        leaving_group_key="UNKNOWN"
                    )
        
        assert result.success
        # Should only have product because UNKNOWN was skipped
        expected_molecules = {
            "product": "C1CCCCC1"
        }
        mock_s1.run.assert_called_once_with(molecules=expected_molecules)

    @patch("rph_core.orchestrator.build_tasks_from_run_config")
    def test_run_tasks_extracts_meta(self, mock_build_tasks, config, tmp_path):
        # We need a hunter instance
        with patch("rph_core.orchestrator.load_config", return_value=config):
            with patch("rph_core.orchestrator.normalize_qc_config", return_value=(config, [])):
                with patch("rph_core.orchestrator.ui.print_pipeline_header"):
                    from rph_core.orchestrator import ReactionProfileHunter, _run_tasks
                    hunter = ReactionProfileHunter()
        
        # Mock task with meta
        task = TaskSpec(
            rx_id="test_rx",
            product_smiles="C1CCCCC1",
            meta={
                "precursor_smiles": "C=C",
                "leaving_small_molecule_key": "AcOH"
            }
        )
        mock_build_tasks.return_value = [task]
        
        run_cfg = config["run"]
        run_cfg["output_root"] = str(tmp_path)
        
        with patch.object(hunter, 'run_pipeline') as mock_run:
            mock_run.return_value = PipelineResult(success=True)
            _run_tasks(hunter, run_cfg)
            
            mock_run.assert_called_once_with(
                product_smiles="C1CCCCC1",
                work_dir=tmp_path / "rx_test_rx",
                skip_steps=[],
                precursor_smiles="C=C",
                leaving_group_key="AcOH"
            )
