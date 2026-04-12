# AGENTS.md — Agentic Coding Guide

**Project**: ReactionProfileHunter  
**Language**: Python 3.8+  
**Lines**: ~54k Python (212 .py files), 723 total files  
**Purpose**: Product-driven reaction mechanism pipeline (S1→S4)

---

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| Add a new pipeline step | `rph_core/orchestrator.py` + `steps/runners.py` | Wire in `run_pipeline()`, add to checkpoint |
| Add a new feature extractor | `rph_core/steps/step4_features/extractors/` | Subclass `BaseExtractor`, call `register_extractor()` |
| Change QC theory level | `config/defaults.yaml` | Never hardcode — single source of truth |
| Fix a path handling bug | `rph_core/utils/path_compat.py` | `normalize_path()` + `is_toxic_path()` |
| Add a new Gaussian template | `config/templates/*.gjf` | Read at runtime, never hardcode |
| Debug checkpoint resume | `rph_core/utils/checkpoint_manager.py` | Hash-validated step state in `pipeline.state` |
| Change conformer search | `rph_core/steps/conformer_search/engine.py` | UCE v3.1: GFN0→GFN2→DFT |
| Add SMARTS for new reaction | `rph_core/steps/step2_retro/smarts_matcher.py` | Registry: [5+2], [4+3], [4+2], [3+2] |
| Fix import violations | `scripts/ci/check_imports.py` | Blocks `from ...utils` |
| Run a specific reaction | `bin/rph_run --smiles "C=C(C)C(=O)O"` | Or `python -m rph_core` |

## PIPELINE FLOW

```
SMILES → S0 (classify) → S1 (anchor) → S2 (scan) → S3 (TS opt) → S4 (features)
         mechanism_       product_      ts_guess     ts_final      features_raw.csv
         summary.json     min.xyz       intermediate  sp_report     features_mlr.csv
                                        .xyz                        feature_meta.json
```

**Inter-step handoffs** (orchestrator wires these):
- S0→S1: `mechanism_summary.json` (forming_bonds, reaction_type)
- S1→S2: `product_min.xyz` (Path), `e_sp` (float), fchk/log/thermo
- S2→S3: `ts_guess.xyz`, `intermediate.xyz`, forming_bonds indices
- S3→S4 (via orchestrator): `ts_final.xyz`, SPMatrixReport, fchk/log, forming_bonds metadata
- Forming bonds resolution runs **between S3 and S4** in orchestrator — not inside either step

## CODE MAP (Top Symbols)

| Symbol | Type | Location | Role |
|--------|------|----------|------|
| `ReactionProfileHunter` | Class | `rph_core/orchestrator.py:86` | Pipeline orchestrator; owns S0-S4 engines + CheckpointManager |
| `PipelineResult` | Dataclass | `rph_core/orchestrator.py:39` | Carries all stage artifacts through pipeline |
| `ConformerEngine` | Class (LoggerMixin) | `rph_core/steps/conformer_search/engine.py:64` | UCE v3.1 two-stage conformer search |
| `TSOptimizer` | Class (LoggerMixin) | `rph_core/steps/step3_opt/ts_optimizer.py:220` | Berny + QST2 rescue + IRC + validation |
| `BaseExtractor` | ABC | `rph_core/steps/step4_features/extractors/base.py:40` | Plugin interface; 14 subclasses |
| `CheckpointManager` | Class | `rph_core/utils/checkpoint_manager.py` | Hash-validated step resume |
| `LinuxSandbox` | Class | `rph_core/utils/qc_interface.py:671` | Context manager for safe QC execution |
| `TaskKind` | Enum | `rph_core/utils/qc_interface.py:161` | OPTIMIZATION/SP/FREQ/TS/IRC/NBO/SCAN |
| `MechanismGraph` | Dataclass | `rph_core/steps/mechanism_classifier/models.py:137` | Graph-based reaction classification |
| `FeatureMiner` | Class | `rph_core/steps/step4_features/feature_miner.py:31` | Discovers plugins, builds context, runs extractors |

---

## Build / Test / Lint Commands

### Run Tests
```bash
# All tests
pytest -v tests/

# Single test file
pytest tests/test_s4_no_qc_execution.py -v

# Single test function
pytest tests/test_s4_no_qc_execution.py::test_extractor_degrades_gracefully -v

# Fast CI gate (import smoke + no-QC tests)
pytest tests/test_imports_step4_features.py tests/test_s4_no_qc_execution.py -v

# S4 contract tests only
pytest tests/test_s4_*.py tests/test_m2_*.py tests/test_m4_*.py -v

# With coverage
pytest --cov=rph_core --cov-report=html
```

### Import Style Check (MANDATORY CI Gate)
```bash
# Block multi-dot relative imports (exit 1 on violation)
python scripts/ci/check_imports.py rph_core

# What it catches:
# ❌ from ...utils.file_io import X  → resolves to wrong path
# ✅ from rph_core.utils.file_io import X  → absolute only
```

### Run Pipeline (Manual Testing)
```bash
# CLI entrypoint
bin/rph_run --smiles "C=C(C)C(=O)O" --output ./Output/rx_001

# Python module
python -m rph_core --smiles "C=C(C)C(=O)O" --output ./Output/rx_001
```

---

## Code Style Guidelines

### Imports (CRITICAL — CI Enforced)
- **ABSOLUTE ONLY**: `from rph_core.utils.file_io import read_xyz`
- **NEVER multi-dot relative**: `from ...utils.file_io` → CI will reject
- **NEVER patch sys.path** to fix imports — fix the import statement
- **NEVER modify rph_core_backup_20260115/** — reference only

### Path Handling
- **ALWAYS use `pathlib.Path`** — no string paths in core code
- **ALWAYS use `normalize_path()`** from `path_compat` before passing to QC tools
- **NEVER use paths with spaces/brackets** — `is_toxic_path()` guards this

### Logging (Library Code)
- **ALWAYS use `logging.getLogger(__name__)`** or `LoggerMixin`
- **NEVER use `print()`** in library code — only in CLI entrypoints
- Log at appropriate levels: DEBUG (details), INFO (milestones), WARNING (degradation), ERROR (failures)

### Type Hints
- **USE typing annotations** for function signatures and dataclasses
- **PREFER `Optional[X]`** over `Union[X, None]`
- **USE `Path`** type for all file/directory parameters

### Naming Conventions
- **Modules**: `snake_case.py`
- **Classes**: `PascalCase`
- **Functions/Methods**: `snake_case()`
- **Constants**: `UPPER_SNAKE_CASE`
- **Private**: `_leading_underscore` for internal use

### Error Handling
- **Fail-fast for config errors** — raise early with clear messages
- **Degrade gracefully for missing artifacts** — log warning + record NaN in `feature_meta.json`
- **NEVER suppress exceptions silently** — always log or raise
- **NEVER use bare `except:`** — catch specific exceptions

### QC Execution (CRITICAL)
- **ALL QC calls MUST go through `utils/qc_interface.py`** — no direct subprocess in steps
- **ALWAYS use sandbox isolation** — each QC task runs in isolated environment
- **ALWAYS check `is_toxic_path()`** before QC execution — paths with `[](){}` break tools

### Configuration
- **NEVER hardcode paths/theory levels** — use `config/defaults.yaml`
- **NEVER fork `defaults.yaml`** — single source of truth (`.bak` files are frozen)
- **Templates**: read at runtime from `config/templates/` — never hardcode in Python

### S4 Extractor Plugins
- **Subclass `BaseExtractor`** and implement required methods
- **Call `register_extractor()`** at module level for discovery
- **Return `FeatureResultStatus`** (COMPLETE / DEGRADED / FAILED)
- **Missing data → degrade with NaN + warning** — never skip silently

---

## Project Structure

```
rph_core/              # Core source — absolute imports only
├── orchestrator.py    # Main pipeline + CLI
├── steps/             # S1-S4 implementations
│   └── step4_features/
│       ├── extractors/     # Plugin subclasses
│       └── schema.py       # FIXED_COLUMNS contract
└── utils/             # 41-file QC/IO/checkpoint infra
    ├── qc_interface.py     # ALL QC calls go here
    ├── config_loader.py    # defaults.yaml reader
    └── checkpoint_manager.py

tests/                 # 53 pytest files (refactored from 75)
├── conftest.py        # Adds repo root to sys.path
├── tmp_v2_2_test/     # Committed S1 fixtures (DON'T DELETE)
├── deprecated/        # 16 archived redundant tests
└── test_imports_*.py  # Fast CI gate

config/                # Runtime configuration
├── defaults.yaml      # Single source of truth
└── templates/         # Gaussian .gjf/.com templates

scripts/               # Utility scripts
├── ci/
│   └── check_imports.py   # Import style gate (exit 1 on violation)
├── run_g16_worker.sh  # Gaussian wrapper
└── train_mlr_loocv.py # ML training script

docs/                  # Supplementary documentation
```

---

## Quick Reference

```python
# Good: Absolute import
from rph_core.utils.file_io import read_xyz

# Good: Path handling
from pathlib import Path
from rph_core.utils.path_compat import normalize_path
path = normalize_path(Path(work_dir) / "file.xyz")

# Good: Logging
import logging
logger = logging.getLogger(__name__)
logger.info("Step completed")

# Good: Error handling with degradation
if not fchk_path.exists():
    logger.warning(f"Missing fchk: {fchk_path}")
    return {key: float('nan') for key in required_keys}

# Bad: Relative import (CI will reject)
from ...utils.file_io import read_xyz

# Bad: String paths
with open(str(path), 'r') as f: ...

# Bad: Print in library code
print("Step completed")

# Bad: Direct subprocess in step
subprocess.run(["g16", "input.gjf"])  # USE qc_interface instead
```

---

## NOTES / GOTCHAS

- **Project path contains brackets** (`[4+3]` in parent dir) — `is_toxic_path()` guards against this; sandbox execution handles it
- **No `pyproject.toml`/`setup.py`** — project uses `sys.path` manipulation via `conftest.py` and `bin/` wrappers; no `pip install` required
- **`rph_core_backup_20260115/`** referenced in docs but **does not exist** in current workspace — legacy reference
- **Backup files in source**: `engine.py.backup_20260117`, `ts_optimizer.py.backup_20260117`, `utils/backup_20260122/` — frozen, do not modify
- **`__init__.py` exports are empty** in `rph_core/utils/` — all imports use absolute paths to specific files
- **`is_path_toxic()` exists in two places**: `qc_interface.py:1415` and `path_compat.py:38` — keep in sync
- **`tests/tmp_v2_2_test/`** contains committed Diels-Alder S1 fixtures — **DO NOT DELETE**
- **Batch mode** uses `ProcessPoolExecutor` — each reaction gets its own `work_dir` and process
- **S0 mechanism classifier** is optional — when it fails, S2 falls back to SMARTS auto-detection
