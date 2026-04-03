# rph_core/steps/AGENTS.md

## OVERVIEW
Step implementations: S1 anchor/conformer → S2 retro scan → S3 TS optimization/rescue → S4 feature extraction. Each step is a self-contained directory with its own `AGENTS.md`.

## WHERE TO LOOK
| Step | Directory | Entry point |
|------|-----------|-------------|
| S1 anchor | `anchor/` | `handler.py` → `AnchorPhase.run()` |
| S1 conformer engine | `conformer_search/` | `engine.py` → `ConformerEngine` |
| S2 retro scan | `step2_retro/` | `retro_scanner.py` → `RetroScanner.run()` |
| S3 TS optimization | `step3_opt/` | `ts_optimizer.py` → `TSOptimizer.run()` |
| S4 feature mining | `step4_features/` | `feature_miner.py:31` → `FeatureMiner.run()` |

## NAMING CONVENTION (v6.3+)
Unified naming for all reaction types ([4+3], [3+2], [4+2], etc.):

| Concept | New Name | Legacy Alias | Notes |
|---------|----------|--------------|-------|
| S2 output | `intermediate.xyz` | `reactant_complex.xyz` | Backward compatible |
| S3 subdirectory | `S3_intermediate_opt/` | `S3_Intermediate/` | Intermediate DFT optimization |
| S4 source label | `S2_intermediate` | `s2_reactant_complex` | Resolved via `naming_compat.py` |

Use `rph_core.utils.naming_compat` for path resolution.

## OUTPUT CONTRACT (strict — orchestrator asserts these)
| Step | Required outputs | Notes |
|------|-----------------|-------|
| S1 | `S1_ConfGeneration/product/product_min.xyz` | `precursor_min.xyz` optional |
| S2 | `S2_Retro/ts_guess.xyz` + `intermediate.xyz` | Both required; forming bonds indices; `reactant_complex.xyz` alias created |
| S3 | `S3_TransitionAnalysis/ts_final.xyz` | Must preserve fchk/log for S4 |
| S4 | `S4_Data/features_raw.csv`, `features_mlr.csv`, `feature_meta.json` | |

## INTER-STEP HANDOFFS
- S1 → S2: `product_min.xyz` (Path), `e_sp` (float)
- S2 → S3: `ts_guess.xyz`, `intermediate.xyz` (or legacy `reactant_complex.xyz`)
- S3 → S4 (via orchestrator): `ts_final.xyz`, `sp_report` (SPMatrixReport), `ts/reactant fchk/log`, `forming_bonds` metadata
- S4 also receives S1 artifacts: `product_thermo.csv`, precursor xyz, shermo summary

## ANTI-PATTERNS
- Missing `intermediate.xyz` and silently continuing into S3/S4 — must fail-fast or rescue.
- Changing output directory layout without updating orchestrator `_resolve_s1_artifacts()` and `path_compat.py`.
- Step implementations calling QC tools directly — all QC goes through `rph_core/utils/qc_interface.py`.
