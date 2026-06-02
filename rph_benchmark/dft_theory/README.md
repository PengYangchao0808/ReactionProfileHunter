# DFT Theory Benchmark

## Overview

This module benchmarks DFT method combinations for `[4+3]` cycloaddition workflows in RPH.

It now uses a staged architecture:

- `run.sh`: thin stage dispatcher
- `lib/state_manager.py`: manifest/state I/O
- `stages/migrate.py`: import confsearch S1 seeds
- `stages/baseline.py`: one-time baseline RPH run per reaction
- `stages/sp_benchmark.py`: fixed-geometry SP benchmark
- `stages/geo_benchmark.py`: fixed-SP geometry benchmark
- `evaluate.py`: final/global report aggregation

## Design Boundaries

- Reuse RPH QC stack (`rph_core/utils/qc_interface.py`, `rph_core/utils/orca_interface.py`)
- Reuse existing templates in `config/templates/`
- Reuse existing parsers/extractors (no local parser reimplementation)
- No `day0_freeze/`, no local `templates/`, no `winner.txt`
- Single source of benchmark state: `manifests/*.json`

## Config Files

- `config/methods_sp.yaml`: Phase-1 single-point candidates
- `config/methods_geo.yaml`: Phase-2 geometry candidates
- `config/benchmark_cases.yaml`: benchmark case list

## Quick Start

```bash
# import seeds, run baseline, SP, GEO, then final report
bash benchmark/dft_theory/run.sh --all --stage all

# import confsearch seeds only
bash benchmark/dft_theory/run.sh --rx-id 1,3,8 --stage migrate

# run one-time baseline RPH stage
bash benchmark/dft_theory/run.sh --rx-id 1,3,8 --stage baseline

# run only SP phase for selected methods
bash benchmark/dft_theory/run.sh --rx-id 1,3,8 --stage sp --methods SP-REF,SP-1

# run GEO phase using phase1 winner as SP reader
bash benchmark/dft_theory/run.sh --session-dir /tmp/rph_dft_bench/session_20260420_120000 --stage geo --geo-methods GEO-1,GEO-3
```

## Main Arguments

- `--rx-id ID[,ID...]`: run selected reaction ids
- `--all`: run all rx_ids from `benchmark_cases.yaml`
- `--stage migrate|baseline|sp|geo|report|all`: run one stage or the full flow
- `--output DIR`: experiment output root (default: `Output/benchmark_dft_theory/experiments`)
- `--session-dir DIR`: reuse existing session directory
- `--confsearch-root DIR`: override imported confsearch benchmark root
- `--methods LIST`: restrict SP methods
- `--geo-methods LIST`: restrict GEO methods
- `--force-clean`: clean stage output before rerun

## Output Layout

```text
Output/benchmark_dft_theory/
├── baselines/
│   └── rx{ID}/
│       └── bl_<signature>/
│           ├── baseline_manifest.json
│           ├── thermo_reference.json
│           ├── stationary_points/
│           ├── S2_Retro/
│           ├── S3_TS/
│           └── S4_Data/
└── experiments/
    └── session_<timestamp>/
        ├── manifests/
        │   ├── benchmark_manifest.json
        │   └── rx{ID}.json
        ├── rx{ID}/
        │   ├── upstream/confsearch/
        │   ├── sp/SP-*/
        │   ├── geo/GEO-*/
        │   └── reports/
        │       ├── sp_summary.json
        │       ├── sp_summary.md
        │       ├── geo_summary.json
        │       └── geo_summary.md
        └── reports/
            ├── phase1_global_summary.json
            ├── phase1_global_summary.md
            ├── phase2_global_summary.json
            ├── phase2_global_summary.md
            ├── final_report.json
            └── final_report.md
```

## Manifest Contract (Brief)

- `benchmark_manifest.json`: global benchmark metadata + stage/report status
- `rx{ID}.json`: per-reaction inputs, matrix, tasks, phase results, winners
- imported confsearch seeds live under `rx{ID}/upstream/confsearch/`
- reusable baseline artifacts live under `baselines/rx{ID}/bl_<signature>/`
- winner propagation:
  - `phase1.winner`
  - `current_phase1_winner`

Phase-2 always reads `current_phase1_winner` from the same `rx{ID}.json`.
