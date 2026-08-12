# tests/AGENTS.md

## OVERVIEW
53 flat-layout pytest files + `conftest.py`; no `__init__.py`; all QC calls mocked — no real ORCA/xTB/CREST binaries required.

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| V3 reintroduction guard | `test_v4_no_legacy_runtime.py` | `FORBIDDEN_MODULES` (25 paths) + `FORBIDDEN_CONFIG_KEYS` (4 keys); the live enforcement of root AGENTS.md anti-patterns |
| V4 protocol contract | `test_v4_protocol_contract.py` | Legacy protocol rejection, output schema validation |
| Checkpoint/resume | `test_v4_checkpoint.py` | Hash-validated stage resume via `pipeline.state` |
| Stage engine | `test_v4_stage_calculator.py` | S3/S4 dispatch contract (legacy test name; tests the refinement engine surface) |
| UI rendering | `test_v4_ui.py` | Rich console, dashboard, ui_adapter normalization |
| Refinement passes | `test_refinement_pass{0,1,2,3}_*.py`, `test_refinement_engine_phase1.py`, `test_refinement_rescue_only.py` | 6 files covering preflight → primary → rescue → canonical |
| S2 PEB | `test_s2_*.py` (7 files) | PEB contract, intermediate selection, profile figures, ORCA GFN2 workflow |
| S0 / S1 | `test_s0_record_restore.py`, `test_v4_s0_context.py`, `test_v4_s1_full.py` | Mechanism record loading, S1 CENSO-LITE full flow |
| Run identity / provenance | `test_run_id_propagation.py`, `test_provenance.py` | Per-run identity in manifests and status snapshots |
| Stale/interrupted recovery | `test_stale_recovery.py`, `test_status_pid_heartbeat.py`, `test_superseded_archive.py`, `test_v4_artifact_reconciler.py` | Heartbeat-based stale detection, archive-before-replace |
| ORCA specifics | `test_orca6_normal_modes.py`, `test_orca_irc.py`, `test_orca_optimization_monitor.py`, `test_v4_orca_convergence.py`, `test_v4_frequency_job.py`, `test_qc_hessian_rendering.py`, `test_orca_failure_classifier` (via resource resolution) | ORCA interface parsing, convergence, frequency validation |

## CONVENTIONS

- **sys.path hack**: `conftest.py` inserts repo root into `sys.path` — tests run without `pip install -e .` (there is no packaging metadata).
- **Flat layout**: no `tests/__init__.py`, no subpackages. All 53 `test_*.py` files sit at top level.
- **Mocked QC**: integration tests mock `qc_interface.py` / `orca_interface.py` / `xtb_runner.py`; no real binaries needed. Do not add tests that shell out to real ORCA/xTB.
- **Naming**: `test_v4_*.py` = V4 contract/architecture tests (CI-relevant). `test_<feature>_*.py` = unit tests. `test_refinement_pass{N}_*.py` = per-pass refinement DAG tests.
- **Fixtures**: `tests/fixtures/` exists but is empty (untracked). `tests/data/` likewise empty. Tests synthesize fixtures in-code or skip if data missing.

## CI GATE vs FULL SUITE

CI (`.github/workflows/ci.yml`) runs **only 4 files**:
```
pytest -q tests/test_v4_protocol_contract.py tests/test_v4_checkpoint.py tests/test_v4_stage_calculator.py tests/test_v4_ui.py
```
The other 49 test files never run in CI. Run the full suite locally:
```
pytest -q tests/
```

## ANTI-PATTERNS

- **Do not** add tests that require real ORCA/xTB/CREST binaries — mock the QC interfaces.
- **Do not** import V3 modules listed in `test_v4_no_legacy_runtime.py::FORBIDDEN_MODULES` — CI fails on any such import.
- **Do not** create `tests/deprecated_v3/` or `tests/deprecated/` — both are referenced in `pytest.ini` `norecursedirs` but neither exists; the dirs were removed with V3.
- **Do not** trust `scripts/AGENTS.md` references to `tests/test_imports_*.py` smoke tests — those files do not exist.
- **Do not** add `__init__.py` to `tests/` — flat layout is intentional; adding one changes pytest collection semantics.

## NOTES

- `pytest.ini`: `testpaths=tests`, `norecursedirs=tests/deprecated tests/deprecated_v3` (stale — neither dir exists, currently a no-op).
- Python matrix: CI pins 3.11 + 3.12. Tests rely on f-strings, `pathlib`, `dataclasses` — 3.8+ compatible in practice.
- `tests/__pycache__/` shows compiled artifacts across pytest 7.4.4, 8.3.4, 9.1.1 — multiple pytest versions work.
- The legacy test name `test_v4_stage_calculator.py` predates the refinement-engine refactor; it now tests the S3/S4 dispatch surface (alias engines + FidelityProfile), not a `stage_calculator.py` file (which is deleted).
