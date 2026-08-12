# rph_core/AGENTS.md

## OVERVIEW
Core package: `v4_orchestrator.py` wires S0→S4; `steps/` holds per-step business logic; `utils/` provides the QC/IO/logging/checkpoint infrastructure shared by all steps. The only supported runtime is `S0 mechanism → S1 CENSO-LITE → S2 PEB → S3 low-level → S4 high-level`.

## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| Pipeline orchestration | `v4_orchestrator.py` | `V4Orchestrator` class; `run()` + CLI `main()` |
| CLI entry | `v4_orchestrator.py` | `main()` — argparse; flags: `--csv`, `--rx-id`, `--output`, `--config`, `--stop-after` |
| Module run | `__main__.py` | `python -m rph_core` → calls `v4_orchestrator.main()` |
| S0 mechanism record | `steps/mechanism_classifier/s0_record.py` | `S0ReactionRecord` loaded from trusted dataset CSV |
| S0 mechanism graph/models | `steps/mechanism_classifier/models.py`, `graph_builder.py`, `dr_completion.py` | `MechanismGraph`, `DRBranchPlan` (V3-era, not yet wired into V4 orchestrator — see P1 in master plan) |
| S1 CENSO-LITE | `steps/conformer_search/censo_lite.py`, `censo_lite_runtime.py` | CREST/GFN2 + B97-3c SP + xTB mRRHO; no DFT OPT/FREQ |
| S1 torsion dedup | `steps/conformer_search/torsion_signature.py`, `deduplicator.py` | Torsion-aware conformer deduplication |
| S2 PEB | `steps/step2_retro/peb_scanner.py` → `peb_engine.py` | `PEBScanner` is a 10-LOC facade over `PEBScanEngine` (4,341 LOC); xTB PEB backward scan from S1 selected product |
| S3/S4 unified engine | `steps/refinement/engine.py` | `RefinementEngine` (4,389 LOC) — 3-pass DAG (preflight → primary → rescue → canonical); no `if stage == "S3"/"S4"` branches |
| S3/S4 stage policy | `steps/fidelity_profile.py` | `FidelityProfile` dataclass — single source of all S3 vs S4 differences; built via `FidelityProfile.from_config(config, stage)` |
| S3/S4 backward-compat aliases | `steps/step3_lowlevel/__init__.py`, `steps/step4_highlevel/__init__.py` | 58-LOC empty subclasses (`LowLevelEngine`/`HighLevelEngine`); orchestrator imports them only into `_BACKWARD_COMPAT_STAGE_ENGINES` for downstream consumers |
| QC job routing | `utils/qc_jobs.py`, `qc_models.py` | `run_optimization()`, `run_frequency()`, `run_single_point()` |
| QC interfaces | `utils/orca_interface.py`, `qc_interface.py` | ORCA + Gaussian + xTB + CREST runners |
| Checkpoint/resume | `utils/v4_checkpoint.py` | `V4Checkpoint` — hash-validated stage resume via `pipeline.state` |
| S4 observability | `utils/s4_progress.py` | `S4ProgressReporter` — writes `status.json` + `events.jsonl` |
| WSL live viewer | `v4_watch.py` | `bin/rph_watch --output <run> --watch` |
| PEB atom mapping | `utils/atom_mapping.py` | Resolves forming bonds from map-space → SMILES-space → XYZ-space |
| Config loading | `utils/config_loader.py` | Loads `config/defaults.yaml` |

## V4 CONFIG CONTRACT
All methods, paths, resources and timeouts belong in `config/defaults.yaml`. Key sections:
- `theory.s3_low_level` / `theory.s4_high_precision` — legacy per-stage theory keys (still read by `FidelityProfile` dual-write during Phase 1)
- `refinement.common` / `refinement.s3` / `refinement.s4` — additive V4 unified-engine profiles (the forward path)
- `step1.protocol: censo_lite` (only allowed protocol)
- `step1.censo_lite` — CREST, ranking, xTB thermo, deduplication parameters
- `step2.scan` — PEB backward scan parameters

## V4 OUTPUT CONTRACT
```
S0_Mechanism/mechanism.json
S1_ConfSearch/product/manifest.json
S2_PEB/manifest.json
S3_LowLevel/manifest.json
S4_HighLevel/manifest.json
pipeline.state
pipeline.result.json
```

## ARCHITECTURE NOTES
- `V4Orchestrator` owns one `V4Checkpoint`; stages are hashed and persisted after each step.
- The pipeline is dataset-only: `--csv <trusted.csv> --rx-id <id>`. No SMILES-only entry.
- `--stop-after {s0,s1,s2,s3,s4}` allows partial pipeline runs.
- Forming bonds are authoritative in the S0 manifest and use 0-based XYZ indices.
- Each S3/S4 structure record retains input, optimized geometry, output files, energies, status and `usable_for_ml`.
- Feature extraction and ML training are external — handled by `RPH_Postprocess`.

## CONVENTIONS
- Absolute `rph_core...` imports only (no multi-dot relative imports).
- `pathlib.Path` everywhere; no string path concatenation in core code.
- `logging.getLogger(__name__)`; no `print()` in library code.
- Route every external QC job through `qc_jobs.py` and the existing QC interfaces.
- Never add `subprocess.run` to a stage module.
- Preserve XYZ, output logs and failure diagnostics when a job degrades.
- A failed S4 job must never delete or replace usable S3 outputs.
- Every stage writes a versioned manifest and is resumable through `pipeline.state`.

## ANTI-PATTERNS
- Direct `subprocess.run` for QC binaries inside steps — bypasses `qc_jobs.py` routing.
- Hardcoding theory methods, paths, or resources in stage code — use `config/defaults.yaml`.
- Reintroducing V3 orchestration (`orchestrator.py`), `QCTaskRunner`, Berny/QST2/IRC rescue, or S2 pre-optimisation.
- Reintroducing SMILES-only entry, multi-protocol S1, or S1 DFT OPT/SP.
- Treating S4 as feature extraction — it is high-precision QC only.
- Creating `COND_xxx/` condition directories inside S1–S4.

## V3 REMOVAL STATUS

V3 modules are **deleted from disk**, not merely deprecated. The list below is historical context only — none of these paths exist in the current tree.

**Authoritative enforcement** (live, runs in CI):
- `tests/test_v4_no_legacy_runtime.py` — `FORBIDDEN_MODULES` (25 paths) + `FORBIDDEN_CONFIG_KEYS` (`gau_xtb`, `neutral_precursor`, `path_search`, `preoptimization`). Importing any forbidden module or adding any forbidden config key fails CI.
- `scripts/ci/check_imports.py` — bans `from ...utils` / `from ....utils` multi-dot relative imports.

**Frozen history** (do not edit, do not trust for current state):
- `docs/ARCHIVE_V3.md` — V3 archive snapshot; §5 "quarantined modules on disk" overstates what remains.
- `docs/V4_CLEANUP_LOG.md` — concise list of actually-removed surfaces; reflects current disk state.

**Removed V3 paths** (all gone from disk): `utils/qc_task_runner.py`, `utils/checkpoint_manager.py`, `utils/oscillation_detector.py`, `utils/gau_xtb_interface.py`, `steps/stage_calculator.py`, `steps/anchor/`, `steps/step3_opt/`, `steps/conformer_search/{engine,funnel,candidates,state_manager}.py`, `steps/conformer_search/pipeline/`, `steps/step2_retro/{kinematic_stretcher,retro_scanner}.py`, `scheduling/`, `lewis_acid/`, `tests/deprecated_v3/`.

**Repurposed survivors**: `steps/conformer_search/protocols.py` is now a V4-only contract (`SUPPORTED_PROTOCOLS = {"censo_lite"}`), not V3 multi-protocol machinery.

**IRC scope rule**: The V3 implicit Berny→QST2→IRC fallback chain is banned. A scoped `RescueMethod.IRC_MIDPOINT_RECOVERY` survives in `steps/refinement/` for recovering a dipolar intermediate, gated by `defaults.yaml` keys `irc_midpoint_recovery` / `irc_max_iter` / `irc_direction`, and runs only after `valid_target_ts` validation — never as INT rescue, never with `UseHess`.
