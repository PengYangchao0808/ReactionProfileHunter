# rph_core/steps/conformer_search/AGENTS.md

## OVERVIEW
V4 S1 conformer search uses **CENSO-LITE** protocol only: CREST/GFN2 conformer
sampling → B97-3c single-point ranking → xTB mRRHO correction. No DFT OPT or
DFT FREQ is performed in S1.

> **Note:** This directory also contains V3-era code (`engine.py`, `funnel.py`,
> `pipeline/`, `protocols.py`, `candidates.py`, `state_manager.py`) that is
> **not called by V4**. See `docs/ARCHIVE_V3.md`. Do not extend V3 modules.

## V4 FILES (active)

| File | Role |
|------|------|
| `censo_lite.py` | `CensoLiteEngine` — V4 S1 orchestrator: embed → CREST → SP → mRRHO → dedup → manifest |
| `censo_lite_runtime.py` | Runtime primitives: `embed()`, `crest_search()`, `split_ensemble()`, `extract_energy()`, `run_sp()`, `run_mrrho()` |
| `torsion_signature.py` | `TorsionSignature` — rotatable-bond dihedral binning for dedup keys |
| `deduplicator.py` | `TorsionAwareDeduplicator` — torsion signature + heavy-atom RMSD dedup |
| `xtb_thermo.py` | `run_xtb_enso()` — xTB `--bhess --enso` mRRHO thermochemistry |

## V4 PROTOCOL (CENSO-LITE)

```
RDKit ETKDGv3 embed → initial.xyz
    ↓
CREST GFN2 conformer search → crest_conformers.xyz
    ↓
Split ensemble → candidates_raw/conf_####.xyz
    ↓
Per candidate: ORCA B97-3c SP + xTB mRRHO → s1_score = E_SP + G_mRRHO
    ↓
Torsion-aware dedup (signature equivalence + RMSD ≤ 0.25 Å)
    ↓
Energy window filter (default 6.0 kcal/mol)
    ↓
Sort by s1_score → manifest.json with `selected` = lowest-energy candidate
```

**No DFT OPT. No DFT FREQ.** S1 only does conformer search + cheap SP ranking.

## CONFIG

| Key | Default | Notes |
|-----|---------|-------|
| `step1.protocol` | `censo_lite` | Only allowed protocol |
| `step1.allowed_protocols` | `[censo_lite]` | V4 gate |
| `step1.censo_lite.energy_window_kcal` | `6.0` | Energy window for valid candidates |
| `step1.censo_lite.crest.*` | — | CREST search parameters |
| `step1.censo_lite.ranking.*` | — | B97-3c SP ranking parameters |
| `step1.censo_lite.xtb_thermo.*` | — | xTB mRRHO parameters |
| `step1.censo_lite.deduplication.*` | — | Torsion + RMSD thresholds |
| `step1.censo_lite.retention.*` | — | Candidate retention policy |

## OUTPUT

```
S1_ConfSearch/<molecule>/
├── manifest.json          # schema_version: s1_censo_lite_v1
├── initial.xyz            # RDKit embed
├── crest/                 # CREST outputs
├── candidates_raw/        # split conformers
└── candidates/            # ranked + deduplicated conf_####.xyz
```

The `manifest.json` `selected` field names the lowest-energy candidate. The V4
orchestrator resolves this to an XYZ path for S2 PEB consumption.

## V3 FILES (deprecated — see docs/ARCHIVE_V3.md)

| File | Status |
|------|--------|
| `engine.py` (138 KB) | V3 `ConformerEngine` — two-stage UCE, NOT used by V4 |
| `funnel.py` | V3 protocol funnel |
| `pipeline/` | V3 pipeline stages (incl. `final_opt_sp.py` — V3 DFT) |
| `protocols.py` | V3 multi-protocol definitions (ext/default/full/lite/zero) |
| `candidates.py` | V3 candidate management |
| `state_manager.py` | V3 state manager |

## ANTI-PATTERNS
- Reintroducing multi-protocol S1 (`ext`/`default`/`full`/`lite`/`zero`) — V4 is CENSO-LITE only.
- Adding DFT OPT or DFT FREQ to S1 — forbidden by V4 theory contract.
- Importing from `engine.py`, `funnel.py`, or `protocols.py` in V4 code.
- Hardcoding binary paths — use `config['executables']`.
