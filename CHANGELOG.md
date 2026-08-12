# Changelog

## [4.0.1] - 2026-08-12

### Added

- Unified S3/S4 `RefinementEngine` (3-pass DAG: preflight → primary → rescue → canonical) under `rph_core/steps/refinement/`.
- `FidelityProfile` dataclass as the single source of S3 vs S4 stage differences (methods, bases, routes, grid/SCF, rescue flags, resources, thermo).
- P0–P2 pipeline hardening modules: `stale_recovery`, `superseded_archive`, `attempt_recorder`, `artifact_reconciler`, `identity`, `run_id`, `provenance`, `soft_mode_review`, `orca_failure_classifier`, `orca_optimization_monitor`, `irc_trajectory`, `ml_quality`.
- S2 PEB refinements: `path_selector` (TS/INT seed selection), `path_profile`, `relaxed_scan_rescue`.
- Hierarchical AGENTS.md refresh (root, `rph_core/`, `rph_core/utils/`, new `tests/AGENTS.md`) with stale V3 references corrected.
- Design docs: V4 atom-mapping contract, V4 ML contract, S3 optimization/rescue normalization plan, S2 ORCA rescue unified selection plan.

### Changed

- Standardized the supported S4 runtime on ORCA: M062X/def2-SVP OPT/OptTS/FREQ followed by wB97M-V/def2-TZVPP SP.
- Upgraded the S4 manifest and checkpoint signature schemas to v3.
- Simplified S2 to one root-level `scan_profile.png` and retired legacy publication/diagnostic plot variants.
- Selected S2 TS guesses from refined B97-3c PATH extrema; xTB PATH estimates are diagnostics only.
- Reframed S2 INT output as an S3 search seed: a pre-TS basin when resolved, otherwise a pre-TS arc-length midpoint.
- CI gate `py_compile` target updated from removed `stage_calculator.py` to `refinement/engine.py`.
- `rph_core/__init__.py` package docstring corrected from stale "v5" to "v4 - S0–S4 DFT Pipeline".

### Fixed

- Propagated ORCA grid and SCF controls into optimization and frequency routes.
- Prevented failed S4 optimizations from being marked usable for ML or checkpointed as complete.

### Removed

- Retired legacy Kneedle, xTB-TS-priority, post-TS plateau, corridor, and multi-panel plotting logic from the active S2 flow.
- Deleted V3 remnants: `stage_calculator.py`, `step3_lowlevel/engine.py`, `step4_highlevel/engine.py` (replaced by `refinement/engine.py` + alias `__init__.py` subclasses).
- Removed stale `pytest.ini` `norecursedirs` entries for nonexistent `tests/deprecated{,_v3}`.
- Excluded `scripts/audit/` ML training scripts from the release scope (out-of-bounds per AGENTS.md; belong in `RPH_Postprocess`).

## [4.0.0] - 2026-07-13

### Added

- S0 mechanism validation and manifest generation.
- S1 CENSO-LITE conformer search with torsion-signature deduplication.
- S2 PEB scanning.
- S3 low-level ORCA B97-3c optimization and r2SCAN-3c single-point calculation.
- S4 high-level Gaussian M062X optimization and ORCA wB97M-V single-point calculation.
- Versioned manifests, checkpoint/resume state, structured UI events, and V4 contract tests.

### Changed

- V4 is the supported runtime: S0 mechanism → S1 CENSO-LITE → S2 PEB → S3 low-level → S4 high-level.
- Feature extraction and ML training remain external responsibilities of `RPH_Postprocess`.

### Removed

- V3 orchestration, `QCTaskRunner`, Berny/QST2/IRC rescue, and S2 pre-optimization.

### Validation

- Automated tests cover V4 protocol contracts, checkpointing, stage calculation, and UI adapters.
- Real WSL/QC validation remains governed by `docs/WSL_TEST_PLAN_V4.md`.
