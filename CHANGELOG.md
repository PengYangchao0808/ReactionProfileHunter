# Changelog

## [Unreleased]

### Changed

- Standardized the supported S4 runtime on ORCA: M062X/def2-SVP OPT/OptTS/FREQ followed by wB97M-V/def2-TZVPP SP.
- Upgraded the S4 manifest and checkpoint signature schemas to v3.

### Fixed

- Propagated ORCA grid and SCF controls into optimization and frequency routes.
- Prevented failed S4 optimizations from being marked usable for ML or checkpointed as complete.

## [4.0.1] - 2026-07-20

### Changed

- Simplified S2 to one root-level `scan_profile.png` and retired legacy publication/diagnostic plot variants.
- Selected S2 TS guesses from refined B97-3c PATH extrema; xTB PATH estimates are diagnostics only.
- Reframed S2 INT output as an S3 search seed: a pre-TS basin when resolved, otherwise a pre-TS arc-length midpoint.

### Removed

- Retired legacy Kneedle, xTB-TS-priority, post-TS plateau, corridor, and multi-panel plotting logic from the active S2 flow.

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
