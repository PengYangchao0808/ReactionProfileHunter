# rph_benchmark/AGENTS.md — Benchmark Infrastructure

## OVERVIEW
Conformer search protocol benchmarks + DFT theory level benchmarks. Not part of the core pipeline — stand-alone evaluation suite for protocol/method selection decisions.

## WHERE TO LOOK

| Module | Focus | Entry Point |
|--------|-------|-------------|
| `confsearch/` | Conformer search protocol comparison (ext/full/lite/zero x multiple molecules) | `python -m rph_benchmark.confsearch.evaluate` |
| `dft_theory/` | DFT method benchmark (geo optimization + SP methods x multiple reactions) | `python -m rph_benchmark.dft_theory.evaluate` |
| `dft_theory/stages/` | 9 stage scripts: migrate → baseline → geo_benchmark → sp_benchmark → shermo → etc. | `bash rph_benchmark/dft_theory/run.sh --stage all` |
| `dft_theory/lib/` | Shared state management (`state_manager.py`) | Internal lib |
| `phasea_driver.py` | Phase A orchestration driver | `python -m rph_benchmark.phasea_driver` |

## KEY DESIGN
- **Confsearch benchmarks** compare CREST/xTB protocol variants (two-stage vs single-stage, energy windows).
- **DFT theory benchmarks** are multi-stage pipelines: migrate confsearch seeds → run baseline → benchmark geo + SP → compute corrections → report.
- Results are written to `rph_benchmark/logs/` and `rph_benchmark/dft_theory/logs/`.
- `dft_theory/stages/` scripts share state via `state_manager.py` (immutable incremental fields pattern).

## CONVENTIONS
- Each stage script is independently runnable (`main()` + `parse_args()`) — no CLI registry.
- Stage scripts MUST NOT overwrite previous stage results — always add `refreshed_*` fields.
- Bash dispatchers (`run.sh`) are the intended invocation path, not direct Python calls.

## ANTI-PATTERNS
- Adding benchmark stages that depend on unreleased core pipeline features.
- Benchmark results stored in `logs/` — do NOT commit large result files to git.
- Stage scripts calling QC tools directly — use `rph_core.utils.qc_interface` when possible.
