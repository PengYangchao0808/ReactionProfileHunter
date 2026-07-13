# Changelog

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
