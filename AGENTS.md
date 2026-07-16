# AGENTS.md — ReactionProfileHunter V4 Coding Guide

## Scope

This repository owns the V4 S0–S4 DFT pipeline. Feature extraction and ML training are external responsibilities of `RPH_Postprocess`.

The only supported runtime is:

```text
S0 mechanism → S1 CENSO-LITE → S2 PEB → S3 low-level → S4 high-level
```

V3 orchestration, `QCTaskRunner`, Berny/QST2/IRC rescue and S2 pre-optimization are removed and must not be reintroduced.

## Authoritative locations

| Concern | Location |
|---|---|
| Pipeline orchestration and resume | `rph_core/v4_orchestrator.py`, `rph_core/utils/v4_checkpoint.py` |
| Single source of configuration | `config/defaults.yaml` |
| S1 runtime | `rph_core/steps/conformer_search/censo_lite.py`, `censo_lite_runtime.py` |
| S1 torsion deduplication | `rph_core/steps/conformer_search/torsion_signature.py`, `deduplicator.py` |
| S0 SMARTS mechanism validation | `rph_core/steps/step2_retro/smarts_matcher.py` |
| S2 PEB | `rph_core/steps/step2_retro/peb_scanner.py` |
| S3/S4 stage policy | `rph_core/steps/step3_lowlevel/`, `step4_highlevel/`, `stage_calculator.py` |
| QC job mapping | `rph_core/utils/qc_models.py`, `qc_jobs.py` |
| Low-level QC interfaces | `rph_core/utils/orca_interface.py`, `qc_interface.py` |
| UI state / status adapters | `rph_core/utils/ui_state.py`, `rph_core/utils/ui_adapter.py` |
| Terminal UI reporter | `rph_core/utils/ui_reporter.py` |
| Shared Rich console / theme | `rph_core/utils/shared_console.py` |
| Logging setup | `rph_core/utils/log_manager.py` |
| Live terminal watcher | `rph_core/v4_watch.py` |

## Theory contract

```text
S1: CREST/GFN2 + B97-3c SP ranking + xTB mRRHO; no DFT OPT/FREQ
S3: ORCA B97-3c OPT/OptTS → ORCA r2SCAN-3c SP
S4: ORCA M062X OPT/OptTS/FREQ → ORCA wB97M-V SP
```

All methods, paths, resources and timeouts belong in `config/defaults.yaml`. Do not hardcode them in stage code.

## Coding rules

- Use absolute `rph_core...` imports only.
- Use `pathlib.Path` for paths.
- Route every external QC job through `qc_jobs.py` and the existing QC interfaces.
- Never add `subprocess.run` to a stage module.
- Preserve XYZ, output logs and failure diagnostics when a job degrades.
- A failed S4 job must never delete or replace usable S3 outputs.
- Every stage writes a versioned manifest and is resumable through `pipeline.state`.
- **UI / logging**
  - Use `logging.getLogger(__name__)` or `LoggerMixin`; no `print()` in library code.
  - V4 logging is configured via `rph_core.utils.log_manager.setup_v4_logging()`. It configures the `root` logger, preserves host handlers, and avoids duplicate `FileHandler` registrations.
  - Console output is Rich-enhanced in TTY mode; it falls back to plain text for `--no-color`, non-TTY, or when `rich` is unavailable.
  - `rph_v4.log` is always plain text and must not contain ANSI escape sequences or Rich markup.
  - Terminal UI uses `rph_core.utils.shared_console.get_console()` and the shared `RPH_THEME`.
  - Progress reporters emit structured events via `event_callback`; UI failures are swallowed and must never abort a QC job.
  - Stage status payloads are normalized through `rph_core.utils.ui_adapter` before rendering; do not assume S3/S4 structure layouts are identical.
  - `rph_core/v4_watch.py` is the only CLI module that may render a full-screen dashboard; stage modules must not use `rich.live.Live`.

## Output contract

```text
S0_Mechanism/mechanism.json
S1_ConfSearch/product/manifest.json
S2_PEB/manifest.json
S3_LowLevel/manifest.json
S4_HighLevel/manifest.json
pipeline.state
pipeline.result.json
```

`forming_bonds` are authoritative in the S0 manifest and use 0-based XYZ indices. Each S3/S4 structure record must retain input, optimized geometry when available, output files, energies, status and `usable_for_ml`.

## Verification commands

```bash
python -m py_compile rph_core/v4_orchestrator.py rph_core/steps/stage_calculator.py
python scripts/ci/check_imports.py rph_core
pytest -q tests/test_v4_protocol_contract.py tests/test_v4_checkpoint.py tests/test_v4_stage_calculator.py tests/test_v4_ui.py
```

Real WSL execution follows `docs/WSL_TEST_PLAN_V4.md`. Do not run a large benchmark until the T0–T6 gates in that document pass.
