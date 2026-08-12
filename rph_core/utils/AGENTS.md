# rph_core/utils/AGENTS.md

## OVERVIEW
78-file cross-step infrastructure: QC execution facade, sandbox/toxic-path enforcement, ORCA/xTB/CREST/Multiwfn/Shermo/ISOSTAT runners, geometry/atom-mapping/fragment utilities, V4 checkpoint/resume, provenance, attempt recording, logging, Rich UI, and IO helpers. Plus nested `thermo/` subpackage (541 LOC, Shermo thermochemistry).

## WHERE TO LOOK
| File | Role |
|------|------|
| `qc_interface.py` (2,350 lines) | **ALL QC calls route here** — Gaussian/xTB/CREST interfaces, `LinuxSandbox`, `is_path_toxic()`, `try_formchk()`, `TaskKind` enum (OPTIMIZATION, SINGLE_POINT, FREQUENCY, TS_OPTIMIZATION, IRC, NBO, SCAN), `QCInterfaceFactory` |
| `qc_jobs.py` (645 lines) | QC job spec mapping — `run_optimization()`, `run_frequency()`, `run_single_point()` |
| `orca_interface.py` (2,877 lines) | ORCA input generation, parsing, execution, optimization monitor |
| `xtb_runner.py` (1,018 lines) | xTB subprocess wrapper; `.run_scan()`, `_write_scan_input()`, `_parse_scan_log()`, `enso_thermo()` |
| `shermo_runner.py` (271 lines) | Shermo thermochemistry runner |
| `isostat_runner.py` (162 lines) | ISOSTAT clustering runner |
| `multiwfn_runner.py` (669 lines) | Multiwfn non-interactive batch runner |
| `orca_failure_classifier.py` (576 lines) | ORCA failure parsing/classification for manifest error reporting |
| `v4_checkpoint.py` (535 lines) | `V4Checkpoint` — hash-validated stage resume via `pipeline.state` |
| `stale_recovery.py` (672 lines) | Detect interrupted/stale runs from heartbeat artifacts before resume |
| `superseded_archive.py` (87 lines) | Archive stale stage outputs before replacement |
| `attempt_recorder.py` (86 lines) | Structured attempt/audit records for QC retries and degradations |
| `artifact_reconciler.py` (204 lines) | Reconcile artifacts after interruptions |
| `run_id.py` (22 lines) | Per-run identity propagated into manifests, status snapshots, UI filtering |
| `provenance.py` (293 lines) | Parent-manifest, atom-mapping and run-chain provenance utilities |
| `stage_scheduler.py` (252 lines) | Stage scheduling and task orchestration |
| `stage_progress.py` (499 lines), `s4_progress.py` (375 lines) | Progress reporters — structured events via `event_callback` |
| `log_manager.py` | Logging setup; `setup_v4_logging()`; `LoggerMixin` base class |
| `geometry_tools.py` (874 lines) | XYZ parsing, geometry manipulation |
| `atom_mapping.py` (667 lines) | Resolves forming bonds from map-space → SMILES-space → XYZ-space |
| `identity.py` (714 lines) | Structure identity / canonical comparison |
| `fchk_reader.py` (468 lines) | Gaussian .fchk parser |
| `forming_bonds_resolver.py` (340 lines) | S3→S4 forming bonds resolution |
| `fragment_cut.py` (515 lines) | Fragment cutting utilities |
| `optimization_config.py` (468 lines) | Optimization parameter building |
| `semantic_slicer.py` (548 lines) | Semantic log slicing |
| `resource_utils.py` (766 lines) | Resource resolution (mem, nproc, maxcore) from config |
| `path_compat.py` | Legacy/new directory layout compatibility |
| `small_molecule_cache.py` (479 lines) | Global cache to avoid re-running S1 on common small molecules |
| `data_types.py` (131 lines) | Shared dataclasses (`QCResult`, `ScanResult`) |
| `ui.py` (1,047 lines), `ui_reporter.py` (819 lines), `ui_adapter.py` (347 lines), `ui_state.py` (56 lines) | UI state/adapters/reporters — stage status payloads normalized through `ui_adapter` before rendering |
| `shared_console.py` (59 lines) | `get_console()` + shared `RPH_THEME` |
| `thermo/` (541 lines, 5 files) | Shermo thermochemistry subpackage — `ThermoRecord`, `ShermoOptions`, `parse_shermo_sum`; see `thermo/__init__.py` docstring |

## KEY QC INTERFACE API (`qc_interface.py`)
```python
# Enums
TaskKind: OPTIMIZATION, SINGLE_POINT, FREQUENCY, TS_OPTIMIZATION, IRC, NBO, SCAN

# Utilities
is_path_toxic(path: Path) -> bool          # detects spaces / [](){} in path
try_formchk(chk_path: Path) -> Optional[Path]  # chk → fchk conversion
harvest_nbo_files(output_dir, jobname, sub_dir) -> Dict[str, Path]

# High-level Gaussian entrypoints
run_gaussian_task(task_kind, xyz_file, output_dir, config, ...) -> QCTaskResult
run_gaussian_optimization(route, atoms, charge, mult, output_dir, config) -> dict

# Interface classes
GaussianInterface   # write_input_file(), optimize(), constrained_optimize()
XTBInterface       # optimize(), scan() (NEW), enso_thermo() (NEW) → delegates to XTBRunner
CRESTInterface     # run_conformer_search(), run_batch_optimization()
QCInterfaceFactory  # create_interface(engine_type, **kwargs)

# Sandbox
LinuxSandbox        # context manager: disk check, isolation, cleanup
GaussianRunner       # .run(sandbox_path, input_content, timeout)
ResultHarvester     # .harvest(sandbox_path, destination_dir) -> Dict[str, Path]
```

## XTB SCAN API
```python
# XTBRunner.run_scan() signature
runner.run_scan(
    input_xyz: Path,
    constraints: Dict[str, float],      # e.g., {"0 1": 2.2}
    scan_range: Tuple[float, float],   # (start, end) distances
    scan_steps: int,
    scan_mode: str = "concerted",      # or "sequential"
    scan_force_constant: float = 1.0,  # constraint force constant
    solvent: Optional[str] = None,
    charge: int = 0,
    uhf: int = 0,
    fixed_constraints: Optional[Dict[str, float]] = None,  # constrained but not scanned
) -> ScanResult

# XTBInterface.scan() signature
xtb.scan(
    xyz_file: Path,
    output_dir: Path,
    constraints: Dict[str, float],
    scan_range: Tuple[float, float],
    scan_steps: int,
    scan_mode: str = "concerted",
    scan_force_constant: float = 1.0,
    charge: int = 0,
    spin: int = 1,
    fixed_constraints: Optional[Dict[str, float]] = None,
) -> ScanResult
```

## XTB ENSO THERMO API
```python
# XTBInterface.enso_thermo() signature
xtb.enso_thermo(
    xyz_file: Path,
    output_dir: Path,
    *,
    charge: int = 0,
    spin: int = 1,
    temperature_k: float = 298.15,
    sthr: float = 50.0,
    imagthr: float = -100.0,
    solvent: Optional[str] = None,
    timeout: Optional[int] = None,
) -> XTBThermoResult
# CENSO-style xTB SPH+mRRHO thermochemistry. Wraps run_xtb_enso() from
# rph_core.steps.conformer_search.xtb_thermo. Delegates to xtb --bhess --enso.
# Returns G(T), ZPVE, H(T) on success, or success=False with error on failure.
```

## PROJECT-SPECIFIC RULES
- **Toxic path**: any path with spaces or `[](){}` must use sandbox execution. `is_path_toxic()` is the gate — check before any QC subprocess.
- **Executable lookup priority**: `config['executables']` → PATH → fallback. xTB and CREST implement this via `resolve_executable_config`. Never hardcode.
- **NBO whitelist**: `NBO_WHITELIST` in `qc_interface.py` controls accepted NBO extensions (`*.37`, `*.nbo`, `*.nbo7`). Do not extend without updating whitelist.
- **Checkpoint hashing**: `V4Checkpoint` (`v4_checkpoint.py`) validates stage signatures on resume via `pipeline.state`. When adding new stage outputs, ensure they are covered by the stage signature or resume may reject/ignore them.
- **LoggerMixin**: every class in core code should inherit `LoggerMixin` — provides `self.logger` without boilerplate.
- **UI failures are swallowed**: progress reporters emit events via `event_callback`; UI failures must never abort a QC job.
- **Plain-text log**: `rph_v4.log` must never contain ANSI escape sequences or Rich markup (Rich is TTY-only).

## ANTI-PATTERNS
- Separate executable-lookup / log-format implementations in multiple modules — centralize in utils.
- Swallowing `stderr`/`stdout` on exception paths — always preserve for diagnostics.
- Direct `subprocess.run` for QC in steps — use `qc_interface.py` wrappers.
- Hardcoding xTB scan force constant — use `scan_force_constant` parameter from config
