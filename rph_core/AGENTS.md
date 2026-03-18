# rph_core/AGENTS.md

## OVERVIEW
Core package: `orchestrator.py` wires S1→S4; `steps/` holds per-step business logic; `utils/` (41 files) provides the QC/IO/logging/checkpoint infrastructure shared by all steps.

## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| Pipeline orchestration | `orchestrator.py:79` | `ReactionProfileHunter` class; `run_pipeline()` + `run_batch()` |
| CLI entry | `orchestrator.py:1287` | `main()` — argparse; flags: --smiles, --output, --config, --log-level, --reaction-type |
| Module run | `__main__.py` | `python -m rph_core` → calls `orchestrator.main()` |
| S1 anchor/conformer | `steps/anchor/`, `steps/conformer_search/` | Two-stage UCE v3.1 engine |
| S2 retro scan | `steps/step2_retro/` | SMARTS matching + bond stretching; supports forward_scan (xTB $scan) |
| S3 TS optimization | `steps/step3_opt/` | Berny + QST2 rescue + IRC + validation |
| S4 feature extraction | `steps/step4_features/` | Plugin-based; 14 extractors |
| QC facade | `utils/qc_interface.py` | ALL subprocess QC calls route here |
| XTB scan (NEW) | `utils/qc_interface.py:991` | `XTBInterface.scan()` — forward scan facade |
| Checkpoint/resume | `utils/checkpoint_manager.py` | Hash-validated step resume (603 lines) |
| Forming bonds S3→S4 | `utils/forming_bonds_resolver.py` | Resolves forming bond indices post-S3 |

## KEY NEW FEATURES
- `reaction_profiles` config drives forward_scan parameters: `scan_start_distance`, `scan_end_distance`, `scan_steps`, `scan_mode`, `scan_force_constant`
- `--reaction-type` CLI arg selects reaction profile (e.g., `[4+3]_default`)
- Forming bonds from S2 preserved through S3→S4 without recomputation

## ARCHITECTURE NOTES
- `ReactionProfileHunter` owns one `CheckpointManager`; step results are hashed and persisted after each step.
- `run_batch()` uses `ProcessPoolExecutor` — each reaction gets its own `work_dir`.
- Forming bonds resolution runs *between* S3 and S4 inside `run_pipeline()` — not inside either step.
- `_resolve_s1_artifacts()` inside orchestrator handles v2.1 vs v3.0/v6.1 directory layout differences.

## CONVENTIONS
- `pathlib.Path` everywhere; no string path concatenation in core code.
- `LoggerMixin` or `logging.getLogger(__name__)`; no `print()` in library code.
- Output directories are idempotent: skip recomputation if prior output exists and checkpoint is valid.
- All QC subprocess calls through `utils/qc_interface.py` (sandbox + toxic-path enforcement).

## ANTI-PATTERNS
- Direct `subprocess.run` for QC binaries inside steps — bypasses sandbox + logging.
- Hardcoding layout assumptions outside `orchestrator._resolve_s1_artifacts()` or `path_compat.py`.
- Treating `rph_output/` or `test_tmpdir/` as source code.
- Forward scan: hardcoded scan params instead of `reaction_profiles` config
- Forward scan: overwriting S2-derived `forming_bonds`
