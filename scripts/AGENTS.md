# scripts/AGENTS.md

## OVERVIEW
Utility scripts for the RPH pipeline. Three categories:
- **QC wrappers** (`run_g16_worker.sh`): Gaussian g16 wrapper for per-job scratch isolation, low-disk protection, environment detection, and auto-cleanup.
- **CI tools** (`ci/check_imports.py`): grep-based import style checker, exit 1 on violation.
- **ML utilities** (`rebuild_all.sh`, `train_yield_ml.py`): dataset construction pipeline and Stage 1 branch-yield training.

## WHERE TO LOOK
| File | Role |
|------|------|
| `run_g16_worker.sh` | Default Gaussian wrapper (`executables.gaussian.wrapper_path` in `config/defaults.yaml`) |
| `ci/check_imports.py` | Scans `rph_core/` for `from ...utils`-style imports; exit 0 = pass, exit 1 = fail |
| `rebuild_all.sh` | End-to-end rebuild: Phase 1 (S4 features) → Phase 2 (condition thermo) → Phase 3 (`rph-features build-dataset`) |
| `build_tables.py` | Upstream raw table generator from dataset/config (retained for P2 review) |
| `backfill_condition_thermo.py` | Condition thermo backfill utility (retained) |
| `rescue_only_rerun.py` | S3 rescue rerun utility (retained) |

## INTEGRATION POINTS
`run_g16_worker.sh` is invoked by:
- `rph_core/utils/qc_interface.py` (`GaussianRunner`)
- `rph_core/steps/step3_opt/qst2_rescue.py`, `irc_driver.py`
- `rph_core/steps/conformer_search/engine.py`

`ci/check_imports.py` is invoked by:
- CI pipeline: `python scripts/ci/check_imports.py rph_core`
- Pre-commit hook: see `docs/CI_INTEGRATION_GUIDE.md`

`rebuild_all.sh` Phase 3 invokes:
- `python -m rph_features.cli build-dataset --rph-root rph_output_backup --output data/processed_rph_v2.2`
- Module: `rph_features.dataset_builder` (single source of truth for ALL ML tables)

`train_yield_ml.py` is invoked by:
- Manual: `python scripts/train_yield_ml.py --data-dir data/processed_rph_v2.2 --output data/ml_yield_results`
- Reads: `stage1_branch_yield_dataset.csv` (produced by `rph-features build-dataset`)
- Outputs: `structured_training_summary.json`, `structured_predictions.csv`

## WHAT ci/check_imports.py CHECKS
| Pattern | Why forbidden | Fix |
|---------|--------------|-----|
| `from ...utils` | Resolves to `rph_core.steps.utils` (doesn't exist) | `from rph_core.utils import ...` |
| `from ....utils` | Ambiguous depth, breaks on restructure | Same fix |

Expected pass output: `✅ No forbidden import patterns found`
Expected fail output: lists file + line + suggested fix

## DATA FLOW: rph-features build-dataset → train_yield_ml.py
```
rph_output_backup/
  └── RXN_*/           → discover_samples()
                         ↓
              rph-features build-dataset
                ├── sample_index.csv
                ├── targets.csv
                ├── features_pre_s3.csv
                ├── condition_table.csv
                ├── condition_features.csv
                ├── retrospective_dataset.csv
                ├── deployable_pre_s3_dataset.csv
                ├── stage1_branch_yield_dataset.csv
                ├── stage1_branch_yield_quality_report.json
                ├── stage1_branch_yield_feature_registry.json
                └── feature_registry.json / quality_report.json / leakage_report.json
                         ↓
              train_yield_ml.py
                ├── structured_training_summary.json
                └── structured_predictions.csv
```

## CONVENTIONS
- `run_g16_worker.sh`: wrapper handles per-job `GAUSS_SCRDIR` isolation — Python side must NOT reuse or clean up the same scratch directory.
- `ci/check_imports.py`: multi-dot relative imports are banned repo-wide — use absolute `from rph_core.utils...`.
- CI smoke tests also live in `tests/test_imports_*.py` (faster than grep-based check for import chain errors).
- Keep machine-specific paths out of `run_g16_worker.sh` — use config/PATH-based discovery.
- `rph-features build-dataset`: single source of truth for ALL ML tables — training scripts must NOT re-merge raw tables.
- `train_yield_ml.py`: only consumes `stage1_branch_yield_dataset.csv` — no raw table merging.

## ANTI-PATTERNS
- Python code manually cleaning or reusing `GAUSS_SCRDIR` set by the wrapper — wrapper owns scratch lifecycle.
- Hardcoding machine-specific paths inside `run_g16_worker.sh`.
- Patching `sys.path` to "fix" an import — change the import statement itself.
- Modifying `rph_core_backup_20260115/` imports then wholesale replacing the source dir.
- Training scripts re-merging condition_table + targets + features — use `rph-features build-dataset` output directly.
- Adding new scripts when functionality belongs in `rph_features.dataset_builder` or `train_yield_ml.py`.
