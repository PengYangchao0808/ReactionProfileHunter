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
| S2 PEB | `steps/step2_retro/peb_scanner.py` → `retro_scanner.py` | xTB PEB backward scan from S1 selected product |
| S3 low-level | `steps/step3_lowlevel/engine.py` | ORCA B97-3c OPT/OptTS → r2SCAN-3c SP |
| S4 high-level | `steps/step4_highlevel/engine.py` | ORCA M062X OPT/OptTS/FREQ → ORCA wB97M-V SP |
| Shared stage engine | `steps/stage_calculator.py` | OPT/Freq/SP dispatch for both S3 and S4 |
| QC job routing | `utils/qc_jobs.py`, `qc_models.py` | `run_optimization()`, `run_frequency()`, `run_single_point()` |
| QC interfaces | `utils/orca_interface.py`, `qc_interface.py` | ORCA + Gaussian + xTB + CREST runners |
| Checkpoint/resume | `utils/v4_checkpoint.py` | `V4Checkpoint` — hash-validated stage resume via `pipeline.state` |
| S4 observability | `utils/s4_progress.py` | `S4ProgressReporter` — writes `status.json` + `events.jsonl` |
| WSL live viewer | `v4_watch.py` | `bin/rph_watch --output <run> --watch` |
| PEB atom mapping | `utils/atom_mapping.py` | Resolves forming bonds from map-space → SMILES-space → XYZ-space |
| Config loading | `utils/config_loader.py` | Loads `config/defaults.yaml` |

## V4 CONFIG CONTRACT
All methods, paths, resources and timeouts belong in `config/defaults.yaml`. Key sections:
- `theory.s3_low_level` — ORCA B97-3c OPT/Freq + r2SCAN-3c SP (with CPCM acetone)
- `theory.s4_high_precision` — ORCA M062X OPT/Freq + ORCA wB97M-V SP (with CPCM acetone)
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

## DEPRECATED V3 CODE (see docs/ARCHIVE_V3.md)
The following V3 modules remain on disk but are **not called by the V4 orchestrator**. They are preserved transitively (some V4 modules import from them) and will be removed in later P-phases:
- `utils/qc_task_runner.py`, `utils/oscillation_detector.py`, `utils/checkpoint_manager.py`
- `steps/anchor/`, `steps/conformer_search/engine.py` (V3 138KB engine, not `censo_lite.py`)
- `steps/step3_opt/`, `steps/dr_aggregator.py`, `steps/condition_thermo.py`
- `lewis_acid/`, `scheduling/v3_scheduler.py`

Do not extend or depend on these modules for V4 work.
