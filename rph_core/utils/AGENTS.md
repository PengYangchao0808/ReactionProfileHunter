# rph_core/utils/AGENTS.md

## OVERVIEW
41-file cross-step infrastructure: QC execution facade, sandbox/toxic-path enforcement, ORCA/xTB/CREST/Multiwfn runners, geometry parsing, checkpoint/resume, logging, and IO helpers.

## WHERE TO LOOK
| File | Role |
|------|------|
| `qc_interface.py` (1614 lines) | **ALL QC calls route here** — Gaussian/xTB/CREST interfaces, `LinuxSandbox`, `is_path_toxic()`, `try_formchk()`, `TaskKind` enum (includes SCAN), `QCInterfaceFactory` |
| `qc_runner.py` | Retry logic, failure type classification (316 lines) |
| `qc_task_runner.py` | Orchestration helper: OPT/SP coupling loop (989 lines) |
| `orca_interface.py` | ORCA input generation, parsing, execution (1124 lines) |
| `xtb_runner.py` (588 lines) | xTB subprocess wrapper; includes `.run_scan()` (NEW), `_write_scan_input()`, `_parse_scan_log()` |
| `shermo_runner.py` | Shermo thermochemistry runner (312 lines) |
| `isostat_runner.py` | ISOSTAT clustering runner |
| `multiwfn_runner.py` | Multiwfn non-interactive batch runner (661 lines) |
| `checkpoint_manager.py` | Step-level resume with artifact hash validation (603 lines) |
| `log_manager.py` | Logging setup; `LoggerMixin` base class |
| `geometry_tools.py` | XYZ parsing, geometry manipulation (837 lines) |
| `fchk_reader.py` | Gaussian .fchk parser (468 lines) |
| `forming_bonds_resolver.py` | S3→S4 forming bonds resolution |
| `fragment_cut.py` | Fragment cutting utilities (497 lines) |
| `optimization_config.py` | Optimization parameter building (369 lines) |
| `oscillation_detector.py` | Geometry oscillation detection → rescue escalation (362 lines) |
| `semantic_slicer.py` | Semantic log slicing (548 lines) |
| `path_compat.py` | Legacy/new directory layout compatibility |
| `small_molecule_cache.py` | Global cache to avoid re-running S1 on common small molecules |
| `data_types.py` | Shared dataclasses (`QCResult`, `ScanResult` (NEW)) |

## KEY QC INTERFACE API (`qc_interface.py`)
```python
# Enums
TaskKind: OPTIMIZATION, SINGLE_POINT, FREQUENCY, TS_OPTIMIZATION, IRC, NBO, SCAN (NEW)

# Utilities
is_path_toxic(path: Path) -> bool          # detects spaces / [](){} in path
try_formchk(chk_path: Path) -> Optional[Path]  # chk → fchk conversion
harvest_nbo_files(output_dir, jobname, sub_dir) -> Dict[str, Path]

# High-level Gaussian entrypoints
run_gaussian_task(task_kind, xyz_file, output_dir, config, ...) -> QCTaskResult
run_gaussian_optimization(route, atoms, charge, mult, output_dir, config) -> dict

# Interface classes
GaussianInterface   # write_input_file(), optimize(), constrained_optimize()
XTBInterface       # optimize(), scan() (NEW) → delegates to XTBRunner
CRESTInterface     # run_conformer_search(), run_batch_optimization()
QCInterfaceFactory  # create_interface(engine_type, **kwargs)

# Sandbox
LinuxSandbox        # context manager: disk check, isolation, cleanup
GaussianRunner       # .run(sandbox_path, input_content, timeout)
ResultHarvester     # .harvest(sandbox_path, destination_dir) -> Dict[str, Path]
```

## XTB SCAN API (NEW)
```python
# XTBRunner.run_scan() signature
runner.run_scan(
    input_xyz: Path,
    constraints: Dict[str, float],      # e.g., {"0 1": 2.2}
    scan_range: Tuple[float, float],   # (start, end) distances
    scan_steps: int,
    scan_mode: str = "concerted",      # or "sequential"
    scan_force_constant: float = 1.0,  # NEW - constraint force constant
    solvent: Optional[str] = None,
    charge: int = 0,
    uhf: int = 0
) -> ScanResult

# XTBInterface.scan() signature
xtb.scan(
    xyz_file: Path,
    output_dir: Path,
    constraints: Dict[str, float],
    scan_range: Tuple[float, float],
    scan_steps: int,
    scan_mode: str = "concerted",
    scan_force_constant: float = 1.0,   # NEW
    charge: int = 0,
    spin: int = 1
) -> ScanResult
```

## PROJECT-SPECIFIC RULES
- **Toxic path**: any path with spaces or `[](){}` must use sandbox execution. `is_path_toxic()` is the gate — check before any QC subprocess.
- **Executable lookup priority**: `config['executables']` → PATH → fallback. xTB and CREST implement this via `resolve_executable_config`. Never hardcode.
- **NBO whitelist**: `NBO_WHITELIST` in `qc_interface.py` controls accepted NBO extensions (`*.37`, `*.nbo`, `*.nbo7`). Do not extend without updating whitelist.
- **Checkpoint hashing**: `CheckpointManager` validates artifact hashes on resume. When adding new step outputs, register them with the checkpoint manager or resume may reject/ignore them.
- **LoggerMixin**: every class in core code should inherit `LoggerMixin` — provides `self.logger` without boilerplate.

## ANTI-PATTERNS
- Separate executable-lookup / retry / log-format implementations in multiple modules — centralize in utils.
- Swallowing `stderr`/`stdout` on exception paths — always preserve for diagnostics.
- Direct `subprocess.run` for QC in steps — use `qc_interface.py` wrappers.
- Hardcoding xTB scan force constant — use `scan_force_constant` parameter from config
