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
from rph_core.utils.checkpoint_manager import CheckpointManager
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
            
            mock_run.assert_called_once()
            call_kwargs = mock_run.call_args[1]  # keyword arguments
            assert call_kwargs["product_smiles"] == "C1CCCCC1"
            assert call_kwargs["precursor_smiles"] == "C=C"
            assert call_kwargs["leaving_group_key"] == "AcOH"
            assert call_kwargs["skip_steps"] == []

    @patch("rph_core.orchestrator.load_config")
    @patch("rph_core.orchestrator.normalize_qc_config")
    @patch("rph_core.orchestrator.ui.print_pipeline_header")
    @patch("rph_core.orchestrator.get_progress_manager")
    @patch("rph_core.orchestrator.notify.notify_completion")
    def test_run_pipeline_marks_s0_skipped_reason(
        self,
        _mock_notify,
        mock_get_pm,
        _mock_ui,
        mock_normalize,
        mock_load,
        config,
        tmp_path,
    ):
        cfg = dict(config)
        cfg["s0"] = {"enabled": True}
        cfg["run"] = {"resume": False}
        mock_load.return_value = cfg
        mock_normalize.return_value = (cfg, [])
        pm = MagicMock()
        mock_get_pm.return_value = pm

        hunter = ReactionProfileHunter()
        result = hunter.run_pipeline(
            product_smiles="C=C",
            work_dir=tmp_path,
            skip_steps=["s1", "s2", "s3", "s4"],
            cleaner_data=None,
        )

        assert result.success is True
        descriptions = [
            call.kwargs.get("description", "")
            for call in pm.update_step.call_args_list
            if call.args and call.args[0] == "s0"
        ]
        assert any("Step 0: Mechanism [SKIPPED: cleaner_data_unavailable]" in d for d in descriptions)

    @patch("rph_core.orchestrator.load_config")
    @patch("rph_core.orchestrator.normalize_qc_config")
    @patch("rph_core.orchestrator.ui.print_pipeline_header")
    def test_run_s0_fallback_adds_rxn_key_hash_and_completes(
        self,
        _mock_ui,
        mock_normalize,
        mock_load,
        config,
        tmp_path,
    ):
        cfg = dict(config)
        cfg["s0"] = {"enabled": True}
        cfg["run"] = {"resume": False}
        mock_load.return_value = cfg
        mock_normalize.return_value = (cfg, [])

        hunter = ReactionProfileHunter()
        checkpoint_mgr = CheckpointManager(tmp_path)
        mock_s0 = MagicMock()
        mock_s0.classify_from_dict.return_value = None

        cleaner_data = {
            "reaction_id": "rx_manual",
            "reaction_type": "cycloaddition",
            "precursor_smiles": "C=C",
            "forming_bonds": [[0, 1]],
        }

        with patch.object(ReactionProfileHunter, "s0_engine", mock_s0):
            summary = hunter._run_s0(
                work_dir=tmp_path,
                product_smiles="CC",
                cleaner_data=cleaner_data,
                checkpoint_mgr=checkpoint_mgr,
                resume_enabled=False,
            )

        assert summary is not None
        assert summary.get("status") == "complete"
        assert summary.get("forming_bonds") == [[0, 1]]

        called_row = mock_s0.classify_from_dict.call_args[0][0]
        assert called_row.get("rxn_key_hash") == "rx_manual"

        s0_dir = tmp_path / "S0_Mechanism"
        assert (s0_dir / "mechanism_graph.json").exists()
        assert (s0_dir / "mechanism_summary.json").exists()

    @patch("rph_core.orchestrator.load_config")
    @patch("rph_core.orchestrator.normalize_qc_config")
    @patch("rph_core.orchestrator.ui.print_pipeline_header")
    def test_s0_status_json_written_on_every_exit_path(
        self,
        _mock_ui,
        mock_normalize,
        mock_load,
        config,
        tmp_path,
    ):
        cfg = dict(config)
        cfg["s0"] = {"enabled": True}
        cfg["run"] = {"resume": False}
        mock_load.return_value = cfg
        mock_normalize.return_value = (cfg, [])
        hunter = ReactionProfileHunter()
        checkpoint_mgr = CheckpointManager(tmp_path)

        result = hunter._run_s0(
            work_dir=tmp_path,
            product_smiles="CC",
            cleaner_data=None,
            checkpoint_mgr=checkpoint_mgr,
            resume_enabled=False,
        )
        assert result is not None
        assert result["status"] == "skipped"
        s0_dir = tmp_path / "S0_Mechanism"
        status_file = s0_dir / "s0_status.json"
        assert status_file.exists()
        import json
        artifact = json.loads(status_file.read_text())
        assert artifact["status"] == "skipped"
        assert artifact["reason"] == "cleaner_data_unavailable"
        assert "timestamp" in artifact

        cfg2 = dict(config)
        cfg2["s0"] = {"enabled": False}
        cfg2["run"] = {"resume": False}
        mock_load.return_value = cfg2
        mock_normalize.return_value = (cfg2, [])
        hunter2 = ReactionProfileHunter()
        result2 = hunter2._run_s0(
            work_dir=tmp_path,
            product_smiles="CC",
            cleaner_data={"reaction_id": "x"},
            checkpoint_mgr=checkpoint_mgr,
            resume_enabled=False,
        )
        assert result2 is not None
        assert result2["status"] == "skipped"
        assert result2["reason"] == "disabled_by_config"
        artifact2 = json.loads(status_file.read_text())
        assert artifact2["status"] == "skipped"
        assert artifact2["reason"] == "disabled_by_config"

    @patch("rph_core.orchestrator.load_config")
    @patch("rph_core.orchestrator.normalize_qc_config")
    @patch("rph_core.orchestrator.ui.print_pipeline_header")
    def test_normalize_cleaner_data_for_pipeline_maps_dataset_json_fields(
        self,
        _mock_ui,
        mock_normalize,
        mock_load,
        config,
    ):
        cfg = dict(config)
        cfg["run"] = {"resume": False}
        mock_load.return_value = cfg
        mock_normalize.return_value = (cfg, [])

        hunter = ReactionProfileHunter()
        cleaner_data = {
            "record_id": "X2003-Sch2-1A",
            "substrate_smiles": "C=C",
            "product_smiles": "C1CCO1",
            "reaction_family": "[4+3] cycloaddition",
        }

        normalized = hunter._normalize_cleaner_data_for_pipeline(
            cleaner_data,
            product_smiles="C1CCO1",
            precursor_smiles=None,
            reaction_profile="[4+3]_default",
        )

        assert normalized is not None
        assert normalized.get("rxn_key_hash") == "X2003-Sch2-1A"
        assert normalized.get("precursor_smiles") == "C=C"
        assert normalized.get("product_smiles_main") == "C1CCO1"
        assert normalized.get("reaction_type") == "[4+3]"
        assert normalized.get("reaction_profile") == "[4+3]_default"
        assert isinstance(normalized.get("raw"), dict)
        assert normalized["raw"].get("reaction_type") == "[4+3]"
