# config/AGENTS.md

## OVERVIEW
Single source of truth for all runtime configuration: `defaults.yaml` (executable paths, theory levels, step parameters, resources) + Gaussian job templates (`templates/*.gjf`/`.com`).

## WHERE TO LOOK
| File | Role |
|------|------|
| `defaults.yaml` | All config keys; never duplicate or fork this file |
| `templates/gaussian_*.gjf` | Gaussian input templates; read at runtime via `config_loader.py` |
| `templates/gaussian_*.com` | Alternative Gaussian input templates |

## KEY CONFIG SECTIONS
| Section | Purpose |
|---------|---------|
| `executables.*` | Paths for Gaussian, ORCA, xTB, CREST, Multiwfn, ISOSTAT, Shermo, wrapper_path |
| `resources.*` | `mem`, `nproc`, `orca_maxcore_safety` |
| `theory.preoptimization` | xTB GFN level, overlap threshold |
| `theory.optimization` | DFT method/basis/dispersion/engine |
| `theory.single_point` | SP method/basis/engine |
| `step1.conformer_search` | `two_stage_enabled`, stage configs, energy windows |
| `reaction_profiles.*.s2_strategy` | S2 strategy routing (`forward_scan` or `retro_scan`) |
| `step2.xtb_settings` | xTB scan controls (`gfn_level`, `solvent`, optional `etemp`) |
| `step3.reactant_opt.enable_nbo` | NBO disabled by default — enable explicitly |
| `run.*` | `source` (single/batch/dataset), `resume`, `output_root` |
| `step4.enabled_plugins` | Filter extractor plugins; if absent, all registered run |
| `step4.mlr.columns` | Override `DEFAULT_MLR_COLUMNS` for `features_mlr.csv` |

## CONVENTIONS
- `defaults.yaml` is the **only** config file — do not create forks. `.bak` files exist but must not diverge in behavior.
- All code reads config through `rph_core/utils/config_loader.py` (passed as dict at runtime).
- Templates are read at runtime; do not hardcode template strings in Python source.
- `oldchk`/`Guess=Read` continuations: if enabling, check template traceability.

## ANTI-PATTERNS
- Copying template strings inline in Python code.
- Maintaining multiple `defaults.yaml` forks with diverged behavior.
- Reading config files directly from steps — config dict is passed from orchestrator.
