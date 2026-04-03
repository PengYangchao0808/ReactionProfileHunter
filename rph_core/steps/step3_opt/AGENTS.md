# step3_opt/AGENTS.md

## OVERVIEW
TS optimization: Berny (primary) → QST2 rescue (on Berny failure) → optional IRC. Validator enforces convergence + imaginary frequency consistency. Post-QC enrichment and intramolecular fragmentation also live here.

## WHERE TO LOOK
| File | Role |
|------|------|
| `ts_optimizer.py` | Coordinator (1000+ lines); owns retry/rescue logic |
| `intermediate_driver.py` | Intermediate DFT optimization driver |
| `berny_driver.py` | Berny TS optimization driver |
| `qst2_rescue.py` | QST2 rescue strategy |
| `irc_driver.py` | IRC computation (430 lines); optional |
| `validator.py` | Convergence + imaginary-freq validation |
| `post_qc_enrichment.py` | Post-optimization QC enrichment (407 lines); ASM/DIAS contract |
| `intramolecular_fragmenter.py` | Fragment splitting for intramolecular reactions (565 lines) |
| `rph_core/utils/oscillation_detector.py` | Detects geometry oscillation; triggers rescue level escalation |

## CONTRACT
**Inputs:**
- `ts_guess.xyz` (from S2)
- `intermediate.xyz` (from S2) — used for intermediate optimization and QST2; legacy: `reactant_complex.xyz`
- `product_xyz` (from S1) — QST2 endpoint
- `e_product_l2` (float) — S1 SP energy
- `product_thermo` (optional Path) — S1 thermo file
- `old_checkpoint` (optional) — reuse S1 artifacts from checkpoint

**Outputs (required):**
- `ts_final.xyz` — validated TS structure
- `sp_report` (SPMatrixReport) — ΔG‡ and ΔG_rxn energies
- `ts_fchk`, `ts_log`, `intermediate_fchk`, `intermediate_log` (Paths, for S4)

**Output directory layout:**
```
S3_TransitionAnalysis/
├── S3_intermediate_opt/    # NEW: Intermediate DFT optimization (was S3_Intermediate)
│   └── intermediate_opt.xyz
├── ts_opt/
│   ├── berny/             # primary Berny attempt
│   └── L2_SP/              # L2 single-point on TS
├── ASM_SP_Mat/             # activation strain matrix
└── (legacy directories preserved for backward compatibility)
```

## RESCUE ESCALATION
Berny fails → `OscillationDetector` checks for geometry oscillation → QST2 rescue → if QST2 fails → mark step as failed (no silent continuation).

## ANTI-PATTERNS
- Silently continuing after Berny failure without triggering rescue or explicit failure.
- Not preserving `log`/`chk`/`fchk` files — S4 requires them for QC-derived features.
- Skipping imaginary frequency validation after convergence.
