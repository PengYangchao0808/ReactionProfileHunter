# ARCHIVE_V3.md — V3 Behaviour Archive (Frozen)

**Status:** Frozen reference document. Do not edit.
**Purpose:** Preserve the historical V3 behaviour so it is available for
reference but does not influence the V4 runtime.
**Applies to:** ReactionProfileHunter v3.0.0–v3.1.0

---

## 1. Purpose & Status

This document is a **frozen archive** of the V3 pipeline behaviour. V4 has
replaced every subsystem described here. Per AGENTS.md, V3 orchestration,
QCTaskRunner, Berny/QST2/IRC rescue, S2 pre-optimization and all other V3
subsystems catalogued below are removed and **must not be reintroduced**.

The V4 supported runtime is:

```
S0 mechanism → S1 CENSO-LITE → S2 PEB → S3 low-level → S4 high-level
```

Nothing in this document describes current V4 behaviour. For V4
implementation details, see `docs/RPH_V4_COMPLETION_MASTER_PLAN.md` and
`AGENTS.md`.

---

## 2. V3 Input Model

V3 accepted pipeline input through three mutually exclusive modes,
configured via `run.source` in `config/defaults.yaml`:

| Mode | Config key | How it worked |
|------|-----------|---------------|
| Single | `run.source: single` | Product SMILES string passed via `--smiles` CLI flag or `run.single.product_smiles` in config. Optional reaction type via `--reaction-type` (e.g. `[4+3]_default`). |
| Batch | `run.source: batch` | TSV file with `rx_id` and `product_smiles` columns. Reactions processed sequentially. |
| Dataset | `run.source: dataset` | CSV with configurable column names for product SMILES, optional precursor SMILES, and row ID. |

**S0 was optional.** When S0 failed to classify a mechanism, S2 fell back to
SMARTS auto-detection of forming bonds directly from the SMILES string.

**V4 replaces all three modes** with a single dataset-only entry:
`--csv <file> --rx-id <id>`. SMILES-only runs are no longer accepted.

---

## 3. V3 Pipeline Stages

### S0 — Mechanism Classification (optional)

- Built a reaction graph from SMILES input.
- Identified cycloaddition mode ([4+3], [5+2], [4+2], [3+2]).
- Detected topology and forming/breaking bonds.
- Produced `mechanism_summary.json` with forming bond indices and reaction type.
- Could plan diastereomeric ratio (DR) branches.
- When S0 failed or was skipped, S2 used SMARTS auto-detection instead.

### S1 — Anchor / Conformer Search

- UCE v3.1 two-stage conformer engine (GFN0 → GFN2 → DFT).
- Multi-protocol funnel: `ext` / `default` / `full` / `lite` / `zero`.
- Protocol-aware selection of conformer search depth.
- DFT optimization and single-point ranking on selected conformers
  (`finalDFT` directory).
- Produced `product_min.xyz` and `precursor_min.xyz` (optional).
- Each molecule had its own working directory under
  `S1_ConfGeneration/[Molecule]/`.

### S2 — Retro Scan

- SMARTS-based forming bond detection.
- xTB constrained scan along forming bonds.
- Gau_XTB TS pre-optimization.
- S2 pre-optimization driver (integrated Gaussian/xTB pre-opt).
- Produced `ts_guess.xyz` and `intermediate.xyz` (legacy alias
  `reactant_complex.xyz`).

### S3 — TS Optimization

- **Berny** as primary optimizer (Gaussian).
- **QST2 rescue** on Berny failure.
- **Optional IRC** validation.
- Optimized reactant geometry in `S3_intermediate_opt/`.
- Dual-level single point: B3LYP/def2-SVP (low) → wB97X-D3BJ/def2-TZVPP (high).
- NBO analysis (optional, controlled by `step3.reactant_opt.enable_nbo`).
- Outputs in `S3_TransitionAnalysis/`:
  - `ts_final.xyz`, `reactant_sp.xyz`
  - `ts_opt/berny/`, `ts_opt/qst2_rescue/`
  - `L2_SP/`, `irc/`

### S4 — Feature Extraction (External in V3, Different Semantics)

In V3, S4 was an **external post-processing step** run by the companion
`RPH_Postprocess` package:

```
rph-features extract --rph-run ./Output/rx_001 --output ./Output/rx_001/S4_Data
```

It consumed S3 artifacts and produced:
- `features_raw.csv` — all extracted features
- `features_mlr.csv` — ML-ready feature matrix
- `feature_meta.json` — provenance and config metadata
- `qc_nbo.37` (optional NBO file)

V3 S4 semantics: **feature extraction and descriptor computation for ML
training**.

V4 S4 semantics have changed: V4 S4 is a **high-precision QC stage**
(Gaussian M062X OPT/OptTS → ORCA wB97M-V SP). Feature extraction remains
the responsibility of `RPH_Postprocess` as a separate step.

---

## 4. V3 Subsystems Removed in V4

| Subsystem | V3 behaviour | V4 replacement |
|-----------|-------------|----------------|
| SMILES-only entry | Pipeline launched from bare SMILES string; S0 optional | Dataset-only `--csv --rx-id`. SMILES-only removed. |
| Multi-protocol conformer search | ext/default/full/lite/zero protocols with DFT ranking | CENSO-LITE only. No protocol variants. |
| S1 DFT OPT/SP (finalDFT) | DFT optimization and single-point ranking after xTB conformer search | CENSO-LITE produces xTB energy + mRRHO ranking only. No DFT OPT/SP in S1. |
| Berny/QST2/IRC rescue chain | Berny → QST2 → IRC with fallback chain in S3 | Intentionally removed. V4 S3 uses ORCA B97-3c OptTS with independent Freq. |
| QCTaskRunner orchestration | Centralized task submission, monitoring, and sandbox management | V4 uses `v4_orchestrator.py` with stage calculators. QCTaskRunner is quarantined. |
| S2 pre-optimization | Gaussian/xTB integrated pre-optimization driver in S2 | Removed. V4 S2 is PEB only. |
| V3 linear CheckpointManager | Hash-validated linear step persistence with partial rehydration | V4 uses `v4_checkpoint.py` with manifest-based checkpoint. |
| Lewis acid calibration | `LewisAcidCalibrationRunner` computing acetone_delta_q_c via DFT | Intentionally removed as V4 runtime feature. Lewis acid models remain on disk due to transitive dependency in `retro_scanner.py`. |
| Lewis acid acidity/detector/modes | Three-mode system (none/strong_single_site/weak_multi_site), structured detection | Same as above: archived but not removed due to imports. |
| DR aggregator | `dr_aggregator.py` — diastereomeric ratio aggregation logic | V4 will reimplement DR as variant-aware (`product_major`/`product_minor`) in P1–P2. |
| condition_feature_merger / condition_thermo | Merged condition descriptors into feature space; thermochemistry under different conditions | Externalized to RPH_Postprocess. Not part of V4 pipeline. |
| V3 feature extraction as S4 | External feature extraction (geometric, electronic, NBO descriptors) | V4 S4 is high-precision QC. Feature extraction stays in RPH_Postprocess. |
| `reaction_profiles` config section | YAML block defining scan parameters per reaction type | Replaced by step2 PEB configuration. |
| V3 `orchestrator.py` | `ReactionProfileHunter` class with full pipeline wiring, CLI entry | Deleted. V4 uses `v4_orchestrator.py`. |

---

## 5. V3 Source Modules Still on Disk (Quarantined)

These modules exist on disk but are **not called by the V4 orchestrator**.
They remain due to transitive import dependencies or because they have not
yet been scheduled for removal under P-phases. They will be removed in later
P-phases once their last import sites are refactored.

### Not imported by V4 code (safe to delete once dependencies resolved)

| Module | Lines | Notes |
|--------|-------|-------|
| `rph_core/utils/qc_task_runner.py` | 1617 | Centralized QC task orchestration. Replaced by QC interfaces + v4_orchestrator. |
| `rph_core/utils/oscillation_detector.py` | 423 | Detected geometry oscillations during optimization. Not used in V4. |
| `rph_core/utils/checkpoint_manager.py` | 908 | V3 linear checkpoint/resume. Replaced by `v4_checkpoint.py`. |
| `rph_core/utils/task_progress.py` | 445 | V3 task progress table. Replaced by events.jsonl + status.json. |
| `rph_core/utils/v3_progress.py` | 390 | V3 progress reporting. Not used in V4. |
| `rph_core/utils/v3_stage_display.py` | 1458 | V3 stage display formatting. Not used in V4. |
| `rph_core/utils/small_molecule_cache.py` | 479 | V3 small molecule caching. Not used in V4. |
| `rph_core/steps/anchor/` | Multiple | S1 anchor phase: conformer routing. Replaced by CENSO-LITE. |
| `rph_core/steps/conformer_search/engine.py` | 3195 | The V3 UCE v3.1 two-stage conformer engine (138KB). NOT `censo_lite.py`. |
| `rph_core/steps/conformer_search/funnel.py` | 352 | V3 protocol-aware conformer funnel. |
| `rph_core/steps/conformer_search/pipeline/` | Multiple | V3 conformer search pipeline stages. |
| `rph_core/steps/conformer_search/protocols.py` | 229 | V3 multi-protocol definitions (ext/default/full/lite/zero). |
| `rph_core/steps/step3_opt/` | Multiple | V3 TS optimization with Berny/QST2/IRC. Replaced by stage_calculator. |
| `rph_core/steps/dr_aggregator.py` | 296 | V3 DR aggregation. To be replaced by variant-aware DR plan in P1. |
| `rph_core/steps/condition_thermo.py` | 1227 | V3 condition thermochemistry. |
| `rph_core/steps/condition_feature_merger.py` | 132 | V3 condition feature merging. |
| `rph_core/steps/contracts.py` | 86 | V3 pipeline contracts. |
| `rph_core/scheduling/v3_scheduler.py` | 1004 | V3 batch scheduling and task orchestration. |

### Imported transitively by V4 code (cannot delete until dependency refactored)

| Module | Lines | Dependency |
|--------|-------|------------|
| `rph_core/lewis_acid/models.py` | Part of `rph_core/lewis_acid/` | Imported by `rph_core/steps/step2_retro/retro_scanner.py` for `LewisAcidQualityFlags`, `CenterCoordinationTrace`, `HALOGEN_SYMBOLS`. |
| `rph_core/lewis_acid/__init__.py` | Part of `rph_core/lewis_acid/` | Re-exports from `models.py`. |
| `rph_core/lewis_acid/mode_resolver.py` | Part of `rph_core/lewis_acid/` | Imports `LewisAcidModeSpec` from models. |
| `rph_core/lewis_acid/calibration.py` | Part of `rph_core/lewis_acid/` | Imports from models. |
| `rph_core/lewis_acid/acidity.py` | Part of `rph_core/lewis_acid/` | Lewis acid acidity functions. |

The entire `rph_core/lewis_acid/` package remains on disk because
`retro_scanner.py` imports `LewisAcidQualityFlags`,
`CenterCoordinationTrace`, and `HALOGEN_SYMBOLS` from
`rph_core.lewis_acid.models`. Until that import is refactored, the lewis_acid
package must stay.

### Scheduler support modules

| Module | Notes |
|--------|-------|
| `rph_core/scheduling/__init__.py` | Package init. |
| `rph_core/scheduling/artifact_refs.py` | Artifact reference models. |
| `rph_core/scheduling/models.py` | Scheduling models. |
| `rph_core/scheduling/small_molecule_precompute.py` | Small molecule precomputation. |

---

## 6. Archived V3 Tests

All V3 tests have been moved to `tests/deprecated_v3/`. They are excluded
from pytest collection via `pytest.ini` (`norecursedirs`) and are not
maintained.

**30 files** in `tests/deprecated_v3/`:

| Test file | V3 subsystem tested |
|-----------|---------------------|
| `test_atom_mapping_index_space.py` | Atom mapping and index space conventions |
| `test_checkpoint_partial_resume.py` | V3 CheckpointManager partial rehydration |
| `test_e2e_precursor_leaving_group.py` | End-to-end precursor with leaving group |
| `test_fast_config.yaml` | Fast test configuration |
| `test_gaussian_env.py` | Gaussian environment detection |
| `test_la_quality_flags.py` | Lewis acid quality flags |
| `test_lewis_acid_acidity.py` | Lewis acid acidity subsystem |
| `test_lewis_acid_calibration.py` | Lewis acid calibration (acetone probe) |
| `test_lewis_acid_detector.py` | Lewis acid structured detection |
| `test_lewis_acid_modes.py` | Lewis acid three-mode system |
| `test_lewis_acid_phase1.py` | Lewis acid phase 1 integration |
| `test_mock_integration.py` | Mock QC integration tests |
| `test_mock_qc_e2e.py` | Mock QC end-to-end |
| `test_molecule_utils.py` | Molecule utility functions |
| `test_mrrho.py` | mRRHO thermochemistry |
| `test_orca_interface.py` | ORCA interface (V3 variant) |
| `test_preopt_driver.py` | S2 pre-optimization driver |
| `test_qc_interface_v52.py` | QC interface v5.2 compatibility |
| `test_retro_scanner_v52.py` | V3 retro scanner v5.2 |
| `test_s2_la_trace.py` | S2 Lewis acid coordination trace |
| `test_s3_checkpoint.py` | S3 checkpoint/resume |
| `test_s3_la_validator.py` | S3 Lewis acid validator |
| `test_s3_ts_rescue_policy.py` | S3 Berny/QST2/IRC rescue policy |
| `test_sandbox_toxic_paths.py` | Sandbox path security |
| `test_small_molecule_cache.py` | Small molecule cache |
| `test_step1_protocol_contract.py` | S1 protocol contract (multi-protocol) |
| `test_step2_path_compat.py` | S2 path compatibility |
| `test_two_stage_conformer.py` | UCE v3.1 two-stage conformer engine |
| `test_xtb_scan_input.py` | xTB scan input generation |

See `tests/deprecated_v3/README.md` for the archived test policy.

---

## 7. V3 to V4 Theory Mapping

| Stage | V3 method | V4 method | Change |
|-------|-----------|-----------|--------|
| S1 conformer ranking | GFN0 → GFN2 → **B3LYP/def2-SVP DFT OPT + SP** (finalDFT) | CREST/GFN2 → **B97-3c SP ranking + xTB mRRHO** | DFT ranking removed. All S1 is semi-empirical + minimal DFT SP only. |
| S2 scan | xTB constrained scan + **Gau_XTB TS pre-optimization** | xTB **PEB** scan | Pre-optimization driver removed. |
| S3 optimization | **Gaussian B3LYP/def2-SVP** OPT/Berny → QST2 rescue | **ORCA B97-3c** OPT/OptTS | Gaussian → ORCA. Rescue chain removed. Single protocol. |
| S3 single point | **Gaussian wB97X-D3BJ/def2-TZVPP** | **ORCA r2SCAN-3c** | Functional changed from wB97X-D3BJ to r2SCAN-3c. Engine changed from Gaussian to ORCA. |
| S4 | **External feature extraction** (geometric, electronic, NBO) via RPH_Postprocess | **Gaussian M062X/def2-SVP OPT/OptTS** then **ORCA wB97M-V/def2-TZVPP SP** | Semantic change: S4 is now high-precision QC, not feature extraction. |

Solvent model: V3 assumed implicit solvation on a per-config basis. V4
mandates `CPCM(acetone)` in all QC routes.

---

## 8. Must-Not-Reintroduce Clause

The following V3 behaviours are **intentionally removed** and **must not** be
reintroduced as default behaviour (from `AGENTS.md` and
`docs/RPH_V4_COMPLETION_MASTER_PLAN.md` section 2.1):

1. SMILES-only pipeline entry.
2. High-precision geometry optimization and high-accuracy single-point
   energy in S1.
3. Implicit, unlimited Berny/QST2/IRC rescue chains.
4. Guessing mechanism or DR variants from temporary directories, stdout, or
   file names.
5. Writing Lewis acid/additive back into canonical SMILES.
6. Writing real Lewis acid descriptors as empirical energy barrier corrections.
7. V3 orchestration via `QCTaskRunner` and linear `CheckpointManager`.
8. S2 pre-optimization drivers.
9. Multi-protocol conformer search (ext/default/full/lite/zero) with DFT
   ranking in S1.
10. Any `subprocess.run` call in stage modules (all QC must route through
    `qc_jobs.py` and QC interfaces).

Violations of this clause are blocking issues for V4 acceptance.

---

## Appendix: Key document references

| Document | Location |
|----------|----------|
| V4 coding guide (AGENTS.md) | `/AGENTS.md` |
| V4 completion master plan | `/docs/RPH_V4_COMPLETION_MASTER_PLAN.md` |
| V3 README (historical reference) | `/README.md` |
| Deprecated V3 tests README | `/tests/deprecated_v3/README.md` |
| S4 completion & observability | `/docs/RPH_V4_S4_COMPLETION_AND_OBSERVABILITY.md` |
| WSL test plan | `/docs/WSL_TEST_PLAN_V4.md` |
