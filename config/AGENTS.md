# config/AGENTS.md

## OVERVIEW
Single source of truth for all runtime configuration: `defaults.yaml` (V4 canonical config — executable paths, theory levels, step parameters, resources). All QC input files (Gaussian `.gjf`, ORCA `.inp`) are generated programmatically from this config.

## WHERE TO LOOK
| File | Role |
|------|------|
| `defaults.yaml` | All V4 config keys; never duplicate or fork this file |
| `mechanism_templates.yaml` | S0 reaction-type template parameters ([4+3], [5+2], etc.) |

## V4 CONFIG SECTIONS
| Section | Purpose |
|---------|---------|
| `schema_version` | `rph_v4_config_v1` |
| `executables.*` | Paths for Gaussian, ORCA, xTB, CREST |
| `executables.discovery.*` | Auto-discovery control |
| `resources.*` | `mem`, `nproc`, `orca_maxcore_safety` |
| `theory.s3_low_level.optimization` | ORCA B97-3c OPT/OptTS method/basis/solvent (CPCM acetone) |
| `theory.s3_low_level.optimization.frequency` | `enabled_for_ts`, `imaginary_cutoff_cm1` (default -50.0) |
| `theory.s3_low_level.single_point` | ORCA r2SCAN-3c SP method/basis/solvent |
| `theory.s4_high_precision.optimization` | Gaussian M062X OPT/OptTS method/basis/solvent (CPCM acetone) |
| `theory.s4_high_precision.optimization.frequency` | TS frequency validation |
| `theory.s4_high_precision.single_point` | ORCA wB97M-V SP method/basis/solvent |
| `step1.protocol` | Must be `censo_lite` (only allowed protocol) |
| `step1.allowed_protocols` | `[censo_lite]` — V4 gate |
| `step1.censo_lite.*` | CREST search, B97-3c SP ranking, xTB mRRHO, deduplication, energy window, retention |
| `step2.scan.*` | PEB backward scan: start/end distance, steps, constraints |
| `run.*` | Resume, output root |

## V3 SECTIONS REMOVED (do not reintroduce)
- `theory.preoptimization` — xTB pre-optimisation (removed in V4)
- `theory.optimization` / `theory.single_point` — flat theory sections (replaced by `s3_low_level` / `s4_high_precision`)
- `reaction_profiles.*` — per-reaction-type scan parameters (removed)
- `step1.conformer_search` — V3 two-stage UCE engine config (replaced by `step1.censo_lite`)
- `step4.enabled_plugins` / `step4.mlr.columns` — V3 feature extraction (now external in `RPH_Postprocess`)
- `optimization_control` — V3 optimisation tuning (removed)

## CONVENTIONS
- `defaults.yaml` is the **only** config file — do not create per-variant, per-test, or per-machine forks. Conditions are a run-root identity (see master plan §3.3), not a config fork.
- All code reads config through `rph_core/utils/config_loader.py` (passed as dict at runtime).
- All QC route lines (including `CPCM(acetone)` solvent) must be generated from config — never hardcoded in stage code.
- `pathlib.Path` for all path values.

## ANTI-PATTERNS
- Maintaining multiple `defaults.yaml` forks with diverged behaviour.
- Reading config files directly from steps — config dict is passed from orchestrator.
- Hardcoding theory methods, basis sets, or solvent models in stage modules.
- Reintroducing V3 config sections (`reaction_profiles`, flat `theory.*`, `step1.conformer_search`).
- Bypassing `resolve_executable_config()` / `declare_resolved_executables()` for executable path lookups.
