# scripts/AGENTS.md

## OVERVIEW
Utility scripts for the RPH pipeline. Two categories:
- **QC wrappers** (`run_g16_worker.sh`): Gaussian g16 wrapper for per-job scratch isolation, low-disk protection, environment detection, and auto-cleanup.
- **CI tools** (`ci/check_imports.py`): grep-based import style checker, exit 1 on violation.
- **ML utilities** (`train_mlr_loocv.py`): post-pipeline ML scripts, not part of the QC pipeline.

## WHERE TO LOOK
| File | Role |
|------|------|
| `run_g16_worker.sh` | Default Gaussian wrapper (`executables.gaussian.wrapper_path` in `config/defaults.yaml`) |
| `ci/check_imports.py` | Scans `rph_core/` for `from ...utils`-style imports; exit 0 = pass, exit 1 = fail |
| `train_mlr_loocv.py` | ML regression training script (LOOCV); not part of the QC pipeline |

## INTEGRATION POINTS
`run_g16_worker.sh` is invoked by:
- `rph_core/utils/qc_interface.py` (`GaussianRunner`)
- `rph_core/steps/step3_opt/qst2_rescue.py`, `irc_driver.py`
- `rph_core/steps/conformer_search/engine.py`

`ci/check_imports.py` is invoked by:
- CI pipeline: `python scripts/ci/check_imports.py rph_core`
- Pre-commit hook: see `docs/CI_INTEGRATION_GUIDE.md`

## WHAT ci/check_imports.py CHECKS
| Pattern | Why forbidden | Fix |
|---------|--------------|-----|
| `from ...utils` | Resolves to `rph_core.steps.utils` (doesn't exist) | `from rph_core.utils import ...` |
| `from ....utils` | Ambiguous depth, breaks on restructure | Same fix |

Expected pass output: `✅ No forbidden import patterns found`
Expected fail output: lists file + line + suggested fix

## CONVENTIONS
- `run_g16_worker.sh`: wrapper handles per-job `GAUSS_SCRDIR` isolation — Python side must NOT reuse or clean up the same scratch directory.
- `ci/check_imports.py`: multi-dot relative imports are banned repo-wide — use absolute `from rph_core.utils...`.
- CI smoke tests also live in `tests/test_imports_*.py` (faster than grep-based check for import chain errors).
- Keep machine-specific paths out of `run_g16_worker.sh` — use config/PATH-based discovery.

## ANTI-PATTERNS
- Python code manually cleaning or reusing `GAUSS_SCRDIR` set by the wrapper — wrapper owns scratch lifecycle.
- Hardcoding machine-specific paths inside `run_g16_worker.sh`.
- Patching `sys.path` to "fix" an import — change the import statement itself.
- Modifying `rph_core_backup_20260115/` imports then wholesale replacing the source dir.
