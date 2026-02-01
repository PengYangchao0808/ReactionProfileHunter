# AGENTS.md (ReactionProfileHunter)

## OVERVIEW
ReactionProfileHunter is a product-first reaction mechanistic pipeline:
S1(anchor/conformer) → S2(retro scan) → S3(TS optimize/rescue) → S4(feature extraction)

**Version**: v6.1.0  
**Language**: Python 3.9+  
**Architecture**: Pipeline orchestrator + step modules + utils infrastructure

## ENTRYPOINTS
- **CLI wrapper**: `bin/rph_run`
- **Python CLI main**: `rph_core/orchestrator.py` (contains `main()` / `__main__`)
- **Integration runner**: `run_auto_test.py` (run by step / resume)

## DIRECTORY STRUCTURE
```
rph_core/                   Core library (orchestrator + steps + utils)
├── orchestrator.py         Pipeline wiring & data flow
├── steps/                  Step implementations (S1-S4)
│   ├── anchor/             S1: Product anchoring (v3.0 molecular autonomy)
│   ├── conformer_search/   Conformer search engine
│   ├── step2_retro/        S2: Retro scanner (TS guess + reactant complex)
│   ├── step3_opt/          S3: TS optimization & rescue (Berny/QST2)
│   └── step4_features/     S4: Feature extraction & packaging
└── utils/                  QC execution, IO, logging, checkpoint
    ├── qc_interface.py     QC execution + sandbox/toxic path handling
    ├── qc_task_runner.py   OPT/Freq driver & route assembly
    ├── qc_runner.py        Retry/failure typing
    ├── orca_interface.py   ORCA runner
    ├── xtb_runner.py       xTB runner
    ├── checkpoint_manager.py  Resume/checkpoint logic
    └── log_manager.py      Logging infrastructure

config/                     Default configs + Gaussian templates
tests/                      pytest (toxic path / degradation / integration)
bin/                        CLI wrapper (adds repo root to sys.path)
scripts/                    External program wrappers (e.g., Gaussian g16)
ci/                         CI scripts & import style checker
test_results/               Run artifacts/examples (NOT source code)
rph_core_backup_20260115/   Historical frozen backup (reference only)
```

## BUILD / INSTALL / TEST COMMANDS

### Install (Development Mode)
```bash
pip install -e ".[dev]"
```

### Run Tests
```bash
# Run all tests
pytest -v tests/

# Run a single test file
pytest -v tests/test_sandbox_toxic_paths.py

# Run a specific test function
pytest -v tests/test_integration.py::TestReactionProfileHunter::test_pipeline_complete

# Run with coverage
pytest --cov=rph_core --cov-report=html tests/

# Run only fast tests (skip those marked with @pytest.mark.slow)
pytest -v -m "not slow" tests/
```

### Lint / Style Check
```bash
# Import style checker (blocks multi-dot relative imports)
python ci/check_imports.py

# Run on specific directory
python ci/check_imports.py rph_core/steps/
```

### Run Pipeline
```bash
# Single reaction (override config)
bin/rph_run --smiles "C=C(C)C(=O)O" --output ./Output/rx_manual

# Config-driven run (edit config/defaults.yaml first)
bin/rph_run

# With custom config
bin/rph_run --config my_config.yaml

# With logging level
bin/rph_run --smiles "CCO" --log-level DEBUG
```

### Integration Test Runner
```bash
# Run all steps
python run_auto_test.py --step ALL

# Resume from checkpoint
python run_auto_test.py --step ALL --resume
```

## CODE STYLE GUIDELINES

### Imports
**CRITICAL**: NEVER use multi-dot relative imports (e.g., `from ...utils`)

✅ **Correct**:
```python
from rph_core.utils.qc_interface import run_gaussian_optimization
from rph_core.utils.log_manager import setup_logger
from rph_core.steps.step2_retro import RetroScanner
```

❌ **Wrong**:
```python
from ...utils.qc_interface import run_gaussian_optimization  # Will fail!
from ..utils import qc_interface  # Ambiguous, avoid
```

**Import Order** (follow PEP 8):
1. Standard library imports
2. Third-party imports (numpy, scipy, rdkit, etc.)
3. Local imports (`rph_core.*`)

### Paths
- **ALWAYS** use `pathlib.Path` (never string paths)
- Use `Path.resolve()` for absolute paths
- Use `normalize_path()` from `rph_core.utils.path_compat` for config paths

### Logging
- Use `logging.getLogger(__name__)` (NEVER `print()` in core library)
- Logging levels: DEBUG (verbose), INFO (progress), WARNING (recoverable), ERROR (failure)
- For orchestrator/steps: inherit from `LoggerMixin` if available

### Type Hints
- Use type hints for all function signatures
- Use `Optional[T]` for nullable types
- Use `Tuple`, `List`, `Dict` from `typing` (Python 3.9+ native also OK)
- Example:
```python
def extract_features(
    ts_xyz: Path,
    reactant_xyz: Optional[Path] = None,
    forming_bonds: Optional[Tuple[Tuple[int, int], ...]] = None
) -> Dict[str, Any]:
    ...
```

### Naming Conventions
- **Functions/methods**: `snake_case` (e.g., `extract_last_converged_coords`)
- **Classes**: `PascalCase` (e.g., `ReactionProfileHunter`, `TSOptimizer`)
- **Constants**: `UPPER_SNAKE_CASE` (e.g., `MAX_RETRIES`, `DEFAULT_TIMEOUT`)
- **Private methods**: `_leading_underscore` (e.g., `_validate_input`)
- **Module-level "private"**: `_leading_underscore` (e.g., `_format_gaussian_route_block`)

### Error Handling
- Define custom exceptions for domain errors:
  ```python
  class GaussianError(Exception): pass
  class TSValidationError(Exception): pass
  ```
- Catch specific exceptions (avoid bare `except:`)
- Log exceptions with `logger.error(..., exc_info=True)` for stack traces
- Re-raise with context when needed:
  ```python
  try:
      result = run_qc_task()
  except QCRunnerError as e:
      logger.error(f"QC task failed: {e}", exc_info=True)
      raise RuntimeError(f"Step3 failed: {e}") from e
  ```

### Docstrings
- Use triple-quoted docstrings for all public functions/classes
- Include Args, Returns, Raises sections (Google/NumPy style)
- Example:
```python
def run_gaussian_optimization(
    input_xyz: Path,
    output_dir: Path,
    route: str
) -> Tuple[Path, float]:
    """Run Gaussian optimization in a safe sandbox.

    Args:
        input_xyz: Input geometry file
        output_dir: Output directory (will be created)
        route: Gaussian route card (e.g., "opt b3lyp/6-31g(d)")

    Returns:
        Tuple of (output_xyz_path, final_energy_hartree)

    Raises:
        GaussianError: If optimization fails or output invalid
    """
```

## PROJECT GOTCHAS (HIGH SIGNAL)

### Toxic Paths
**CRITICAL**: Spaces and `[](){}` in paths break QC programs on WSL/Linux
- Related: `rph_core/utils/qc_interface.py` + `tests/test_sandbox_toxic_paths.py`
- **Always** use `is_toxic_path()` / `LinuxSandbox` for QC execution
- **Never** write custom `subprocess.run(cwd=...)` logic in steps

### Step Output Contracts
**S2 MUST output BOTH** `ts_guess.xyz` AND `reactant_complex.xyz`
- `reactant_complex.xyz` is required by S3 (QST2 rescue) and S4 (distortion/fragment)
- Missing either file should fail-fast, not silently proceed

**S4 Fixed Outputs** (v6.1):
- `S4_Data/features_raw.csv`
- `S4_Data/features_mlr.csv`
- `S4_Data/feature_meta.json`

### Key Directories
- `test_results/` is **run artifacts/examples**, NOT source code
  - Fix bugs in `rph_core/`, not in `test_results/`
- `rph_core_backup_20260115/` is frozen backup for reference
  - **Never** modify; use only for comparison

### Test Discovery
- `pytest.ini` restricts `testpaths=tests`
- Root-level `test_*.py` are manual verification scripts, NOT pytest tests

## ANTI-PATTERNS

❌ **Don't**: Implement QC runners inside step modules  
✅ **Do**: Reuse `rph_core/utils/qc_interface.py` or `qc_task_runner.py`

❌ **Don't**: Fix bugs with large refactors (repo has backups + large artifact dirs)  
✅ **Do**: Make minimal, targeted fixes

❌ **Don't**: Use `sys.path` hacks to fix imports  
✅ **Do**: Fix import statements to use absolute imports

❌ **Don't**: Hardcode file paths in extractors  
✅ **Do**: Get paths from `FeatureContext` / `PathAccessor`

❌ **Don't**: Silently skip missing artifacts  
✅ **Do**: Degrade gracefully with logging + status tracking

❌ **Don't**: Write to stdout/print in library code  
✅ **Do**: Use logger (stdout is for CLI/user-facing messages only)

## WHERE TO LOOK (Quick Reference)

| Task | Location |
|------|----------|
| Pipeline wiring / data flow | `rph_core/orchestrator.py` |
| Step2 Retro scan | `rph_core/steps/step2_retro/retro_scanner.py` |
| Step3 TS optimize/rescue | `rph_core/steps/step3_opt/ts_optimizer.py` |
| Step4 Feature extraction | `rph_core/steps/step4_features/feature_miner.py` |
| QC execution + sandbox/toxic path | `rph_core/utils/qc_interface.py` |
| ORCA / xTB runners | `rph_core/utils/orca_interface.py`, `rph_core/utils/xtb_runner.py` |
| Checkpoint/resume | `rph_core/utils/checkpoint_manager.py` |
| Import style rules | `ci/check_imports.py` |

## ADDITIONAL RESOURCES
- **README.md**: High-level architecture & I/O contract
- **Sub-AGENTS.md files**: Domain-specific details in each subdirectory
- **IMPORT_GUIDELINES.md** (if exists): Repository-wide import conventions
