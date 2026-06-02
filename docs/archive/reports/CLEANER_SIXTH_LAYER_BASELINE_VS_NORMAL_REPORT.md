# CLEANER SIXTH LAYER — Baseline vs Normal Run Deep-Dive (rx1)

## Scope / invariant

This layer assumes the **forming-bond gold standard** is defined purely by the input XYZ (1-based indexing):

- **Gold standard forming bonds (XYZ 1-based): 12–19 and 15–16**

Everything below is a **baseline vs normal-run differential** intended to explain why the **benchmark baseline** for rx1 can look dramatically worse (S2 degradation, topology drift, negative barrier) even though both pipelines are still propagating the *wrong* forming-bond pairs in their metadata.

## Artifacts used (primary evidence)

### Baseline (benchmark_dft_theory)

Baseline cache directory:

- `Output/benchmark_dft_theory/baselines/rx1/bl_3f31ea3a42fe5184/`

Key files:

- `baseline_manifest.json`
- `pipeline.state`
- `S2_Retro/scan_profile.json`
- `S3_TS/mechanism_meta.json`
- `S3_TS/sp_matrix_metadata.json`

### Normal run (rph_output)

- `rph_output/rx_1/pipeline.state`
- `rph_output/rx_1/S2_Retro/scan_profile.json`
- `rph_output/rx_1/S3_TS/mechanism_meta.json`
- `rph_output/rx_1/S3_TS/sp_matrix_metadata.json`

## Finding A — S2.2 xTB path_search is the most visible pipeline divergence

### Baseline: no S2.2 path_search

Baseline `pipeline.state` indicates:

- `steps.step_s2.metadata.step2_signature.path_search.enabled = false`

Baseline `S2_Retro/scan_profile.json` indicates:

- `generation_method = "retro_scan_from_product"`
- `scan_quality.status = "DEGRADED"`
- `scan_quality.boundary_maximum = true`
- `scan_quality.trajectory_check.off_path_count = 9` (topology drift)

### Normal run: S2.2 path_search executed

Normal `pipeline.state` indicates:

- `steps.step_s2.metadata.s2_generation_method = "xtb_path_search"`
- `steps.step_s2.metadata.step2_signature.path_search.enabled = true`
- `config_snapshot.step2.path_search.enabled = true`

Normal `S2_Retro/scan_profile.json` indicates:

- `generation_method = "xtb_path_search"`
- `ts_quality.status = "COMPLETE"` and `ts_guess_confidence = "high"`
- It does **not** carry the baseline-style topology-drift diagnostics.

### Implication

Even with the forming-bond metadata still wrong, **running S2.2 path_search tends to produce a smoother TS-guessing pipeline** (at least for rx1), while the baseline’s relaxed-scan route is much more prone to **topology drift + boundary peak** failure modes.

## Finding B — Baseline TS validation is too weak and can “pass” with empty bond-length checks

Baseline `baseline_manifest.json` contains:

- `ts_validation.status = "pass"`
- `ts_validation.bond_lengths = []`

So baseline acceptance is currently dominated by “imaginary frequency exists” while **geometric verification can be skipped entirely**.

## Finding C — The baseline produces an unphysical negative activation energy; the normal run does not

From `S3_TS/sp_matrix_metadata.json`:

### Baseline

- `activation_energy_kcal = -44.5287625399651`
- `reaction_energy_kcal = -52.09099515737395`

### Normal run

- `activation_energy_kcal = 26.6102295081364`
- `reaction_energy_kcal = -25.707389210652167`

### Interpretation

The baseline’s TS/reactant ordering is inconsistent with a conventional barrier definition, and that inconsistency is **not caught** by baseline validation (Finding B).

## Finding D — Index-base / forming-bond representation is inconsistent across layers (still)

Normal `S3_TS/mechanism_meta.json` shows:

- Top-level `index_base = 0`
- Top-level `forming_bonds = [[11,12],[14,15]]` (0-based)
- But `source.step2_signature.forming_bonds = [[12,13],[15,16]]` (appears 1-based)

Baseline `S3_TS/mechanism_meta.json` shows:

- `index_base = 0` but `forming_bonds = [[12,13],[15,16]]` (appears 1-based)

This reinforces the earlier layers’ conclusion: **index-base corruption and/or mixed conventions are still present in the metadata chain**, and for this layer we treat that as an upstream invariant to hold constant while explaining baseline-vs-normal severity.

## Critical open question (needs code-provenance confirmation)

Repo `config/defaults.yaml` currently has `step2.path_search.enabled: false`, yet the normal run `pipeline.state` snapshot records `step2.path_search.enabled: true`.

### Closed: this is a defaults.yaml *revision* mismatch, not an in-memory override

There is no code path in current `rph_core` that flips `step2.path_search.enabled` to true at runtime; the value is read from the loaded YAML and snapshotted directly into `pipeline.state`:

- `rph_core/steps/runners.py` reads it via `path_search_enabled = bool(step2_cfg.get("path_search", {}).get("enabled", False))`.
- `rph_core/orchestrator.py` snapshots config via `checkpoint_mgr.initialize_state(..., config=self.config)`.
- `rph_core/utils/checkpoint_manager.py:compute_step2_signature()` embeds `path_search.enabled` from `config["step2"]["path_search"]["enabled"]`.

The refactor plan explicitly documents a default flip:

- `docs/v211_refactor_plan.md`: “`step2.path_search.enabled` 默认值从 `true` 改为 `false`”.

And we can directly confirm that an earlier `defaults.yaml` revision had it enabled:

- `git show 921b8d6e:config/defaults.yaml` contains `step2.path_search.enabled: true`.

This aligns with timestamps:

- Normal run `rph_output/rx_1/pipeline.state` timestamp is 2026-04-06 (pre-refactor flip), and records `enabled=true`.
- Baseline run is 2026-04-21 and records `enabled=false`.

Therefore the baseline-vs-normal divergence in S2.2 is explained by **config revision provenance**, not hidden CLI flags or resume logic.

## Recommended next checks / experiments (no-code-change)

1. **Provenance check:** locate in code how `step2.path_search.enabled` can become `true` at runtime even if `defaults.yaml` has it `false`.
2. **Compare TS validation rules:** baseline manifest passes with empty `bond_lengths`; normal-run packaging does not produce such a manifest, so decide whether baseline should require at least one geometric check.
3. **Reproduce minimal delta:** re-run baseline once with only `step2.path_search.enabled=true` (holding everything else constant) to see if baseline’s S2 drift and negative barrier disappear.

## Minimal mitigations suggested by evidence

- In benchmark baseline validation, treat `bond_lengths=[]` as **not validated** (i.e., downgrade to degraded/fail) rather than “pass”.
- If `scan_quality.status=DEGRADED` due to topology drift, baseline should not accept downstream TS as “pass” without additional structural verification.
