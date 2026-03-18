from rph_core.steps.step4_features.feature_miner import FeatureMiner


def test_feature_miner_defaults_to_disallow_job_runs() -> None:
    miner = FeatureMiner(config={"step4": {"job_run_policy": "disallow"}})
    assert miner.config["step4"]["job_run_policy"] == "disallow"
