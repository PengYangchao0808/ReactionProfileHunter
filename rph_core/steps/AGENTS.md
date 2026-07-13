# rph_core/steps/AGENTS.md

## OVERVIEW
V4 step implementations: S0 mechanism → S1 CENSO-LITE conformer search → S2 PEB backward scan → S3 low-level QC → S4 high-level QC. Each step is a self-contained directory with versioned manifest output.

## WHERE TO LOOK
| Step | Directory | Entry point | Theory |
|------|-----------|-------------|--------|
| S0 mechanism | `mechanism_classifier/` | `s0_record.py` → `load_s0_reaction_record()` | Trusted dataset CSV input |
| S1 CENSO-LITE | `conformer_search/` | `censo_lite.py` → `CensoLiteEngine.run()` | CREST/GFN2 + B97-3c SP + xTB mRRHO |
| S2 PEB | `step2_retro/` | `peb_scanner.py` → `RetroScanner.run()` | xTB PEB backward scan |
| S3 low-level | `step3_lowlevel/` | `engine.py` → `LowLevelEngine.run()` | ORCA B97-3c OPT/OptTS → r2SCAN-3c SP |
| S4 high-level | `step4_highlevel/` | `engine.py` → `HighLevelEngine.run()` | Gaussian M062X OPT/OptTS → ORCA wB97M-V SP |
| Shared engine | — | `stage_calculator.py` → `StageCalculator` | OPT/Freq/SP dispatch for S3 + S4 |

## V4 OUTPUT CONTRACT
```
S0_Mechanism/mechanism.json
S1_ConfSearch/product/manifest.json
S2_PEB/manifest.json
S3_LowLevel/manifest.json
S4_HighLevel/manifest.json
```

Each S3/S4 structure uses flat directory naming:
- `<variant>` — product (e.g. `product_major`)
- `<variant>_int` — PEB intermediate
- `<variant>_ts` — PEB TS candidate
- `precursor` — precursor minimum

No nested `intermediate/` or `ts/` subdirectories. No `conf_0001/` directories in S2–S4.

## INTER-STEP HANDOFFS
- S0 → S1: `mechanism.json` (forming bonds, mapped SMILES, reaction type)
- S1 → S2: `manifest.json` with `selected` candidate → resolved `selected.xyz` path
- S2 → S3: `ts_guess.xyz`, `intermediate.xyz`, forming bonds, scan profile
- S3 → S4: per-structure `{opt,freq,sp}/` outputs with energies and status
- S4 → external: high-precision energies + geometries for `RPH_Postprocess`

## ANTI-PATTERNS
- Reintroducing Berny/QST2/IRC rescue, `QCTaskRunner`, or V3 `step3_opt/` logic.
- Treating S4 as feature extraction — it is high-precision QC (M062X/wB97M-V).
- Creating `conf_0001/` directories downstream of S1 — S1 exports `selected.xyz`, S2–S4 consume only that.
- Step implementations calling QC tools directly — all QC goes through `rph_core/utils/qc_jobs.py`.
- Nested variant subdirectories (`product_major/intermediate/`) — use flat names (`product_major_int/`).

## DEPRECATED V3 DIRECTORIES (see docs/ARCHIVE_V3.md)
- `anchor/` — V3 S1 anchor phase (not used by V4)
- `conformer_search/engine.py` — V3 two-stage UCE engine (V4 uses `censo_lite.py`)
- `step3_opt/` — V3 Berny/QST2/IRC (V4 uses `step3_lowlevel/`)
- `dr_aggregator.py`, `condition_thermo.py`, `condition_feature_merger.py`, `contracts.py` — V3-era
