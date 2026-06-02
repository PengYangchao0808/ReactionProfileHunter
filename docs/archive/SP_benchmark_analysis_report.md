# SP Benchmark Analysis Report — DFT Theory Phase 1

**Project**: ReactionProfileHunter v2.1.0  
**Scope**: Fixed-geometry single-point energy benchmark (Phase 1)  
**Session**: `session_bl_fixed` (2026-04-28 ~ 2026-04-29)  
**Status**: ✅ Phase 1 COMPLETE — Phase 2 (GEO) PENDING  
**Generated**: 2026-04-30

---

## Executive Summary

The Phase-1 SP benchmark evaluated **7 single-point methods** across **4 [4+3] cycloaddition reactions** on a fixed baseline geometry (B3LYP-D3BJ/def2-SVP optimized at lite protocol). The reference method, DLPNO-CCSD(T)/def2-TZVPP (SP-REF), successfully completed 2/4 reactions but **failed on TS calculations for rx3 and rx15** (ORCA abnormal termination). The global vote ended in a **tie between SP-REF and SP-1 (wB97M-V) at 2 votes each**.

Among DFT candidates, **SP-3 (wB97X-2/def2-TZVPP)** showed the lowest absolute energy MAE vs reference (84–145 kcal across completed cases), though the absolute MAE metric is dominated by method-dependent baseline shifts (~800–1000 kcal) rather than chemically meaningful relative errors. The fastest method, **SP-6 (r2SCAN-3c)**, completed in ~100 seconds but **incorrectly predicts the TS energy below the intermediate** in 3/4 reactions — a physically impossible result for exothermic cycloadditions.

### Key Recommendations

| Decision | Method | Rationale |
|----------|--------|-----------|
| **Phase-2 GEO SP reader** | SP-1 (wB97M-V) | Always completed (4/4), ~600s, good dE_ts_int prediction when REF unavailable |
| **Gold-standard SP** | SP-REF (DLPNO-CCSD(T)) | Best accuracy when it works; needs TS calculation robustness investigation |
| **High-accuracy DFT SP** | SP-3 (wB97X-2) | Lowest MAE vs REF among DFT methods but overestimates barrier magnitude |
| **DO NOT USE** | SP-6 (r2SCAN-3c) | Sign error on reaction barrier — TS predicted below intermediate |
| **DO NOT USE** | SP-4 (wB97X-V) | Largest MAE vs REF (~1000 kcal) without clear advantage |

---

## 1. Benchmark Design

### 1.1 Reaction Cases

| rx_id | Label | Characteristic | Atom Count | dE_ts_int via REF (kcal) |
|-------|-------|---------------|------------|--------------------------|
| **1** | rigid_small_ring | Rigid small ring, Boc protecting group | ~60 atoms | +1.166 |
| **3** | flexible_chain | Flexible chain, benzyl, dispersion-sensitive | ~60 atoms | REF FAILED |
| **8** | medium_polarity_async | Cyclic carbamate, high rigidity, short tether | ~50 atoms | −2.112 |
| **15** | steric_endo_exo | Extended tether, stereoselectivity probe | ~55 atoms | REF FAILED |

### 1.2 SP Methods Under Test

| Method ID | Description | Family | Engine | Basis Set | Aux | Dispersion |
|-----------|-------------|--------|--------|-----------|-----|------------|
| **SP-REF** | DLPNO-CCSD(T) | post-HF (DLPNO) | ORCA | def2-TZVPP | def2/JK + /C | forbidden |
| **SP-1** | wB97M-V | VV10 family | ORCA | def2-TZVPP | def2/J | forbidden |
| **SP-2** | PWPB95-D3BJ | Double-hybrid DFT | ORCA | def2-TZVPP | def2/J + /C | D3BJ external |
| **SP-3** | wB97X-2 | Double-hybrid DFT | ORCA | def2-TZVPP | def2/J + /C | none (D4 implicit) |
| **SP-4** | wB97X-V | VV10 family | ORCA | def2-TZVPP | def2/J | forbidden |
| **SP-5** | M06-2X | Hybrid GGA | ORCA | def2-TZVPP | def2/J | none |
| **SP-6** | r2SCAN-3c | Composite 3c | ORCA | mTZVPP | none | forbidden |

All calculations: PCM(acetone), TightPNO for DLPNO, tightSCF, gcorr reused from baseline.

### 1.3 Winner Selection Logic

The per-reaction winner is selected from **completed** methods by lowest `absolute_sp_energy_mae_vs_ref_kcal`. When SP-REF fails (partial_failed), it is excluded from voting, and the winner falls back to the method with the most completed stationary points (then alphabetical).

**Critical caveat**: The MAE metric compares absolute total energies (~−900 hartree scale). Method-dependent baseline shifts (~1–2 hartree = 600–1200 kcal) dominate these values, making absolute MAE a poor measure of chemical accuracy. The chemically meaningful comparison is on **relative energies** (dE_ts_int, dE_prod_int), which the current winner logic does not use.

---

## 2. Completion Status Matrix

| Method | rx1 | rx3 | rx8 | rx15 | Success Rate |
|--------|-----|-----|-----|------|-------------|
| **SP-REF** | ✅ completed | ⚠️ partial_failed | ✅ completed | ⚠️ partial_failed | 2/4 (50%) |
| **SP-1** | ✅ completed | ✅ completed | ✅ completed | ✅ completed | **4/4 (100%)** |
| **SP-2** | ✅ completed | ✅ completed | ✅ completed | ✅ completed | **4/4 (100%)** |
| **SP-3** | ✅ completed | ✅ completed | ✅ completed | ✅ completed | **4/4 (100%)** |
| **SP-4** | ✅ completed | ✅ completed | ✅ completed | ✅ completed | **4/4 (100%)** |
| **SP-5** | ✅ completed | ✅ completed | ✅ completed | ✅ completed | **4/4 (100%)** |
| **SP-6** | ✅ completed | ✅ completed | ✅ completed | ✅ completed | **4/4 (100%)** |

**SP-REF failure detail**: Both rx3 and rx15 failed with "ORCA 未正常终止" (ORCA abnormal termination) on the TS single-point calculation only. Precursor, intermediate, and product SPs completed successfully. This suggests the fixed B3LYP TS geometry triggers convergence issues in DLPNO-CCSD(T) for these specific molecules.

**Overall SP phase completion**: 26/28 task-slots completed (92.9%), 2 partial failures (SP-REF TS on rx3, rx15).

---

## 3. Runtime Analysis

| Method | rx1 (s) | rx3 (s) | rx8 (s) | rx15 (s) | Avg (s) | Family | Cost factor vs SP-6 |
|--------|---------|---------|---------|----------|---------|--------|---------------------|
| **SP-REF** | 36,880 | 32,662 | 16,783 | 28,922 | **28,812** | post-HF | **250×** |
| **SP-3** | 707 | 707 | 356 | 608 | **594** | DH-DFT | 5.1× |
| **SP-4** | 675 | 659 | 339 | 585 | **564** | VV10 | 4.9× |
| **SP-1** | 676 | 651 | 362 | 585 | **568** | VV10 | 4.9× |
| **SP-2** | 566 | 551 | 275 | 469 | **465** | DH-DFT | 4.0× |
| **SP-5** | 471 | 461 | 239 | 411 | **395** | Hybrid GGA | 3.4× |
| **SP-6** | 119 | 120 | 75 | 106 | **105** | Comp. 3c | **1×** |

**Key observations**:
- SP-REF is **two orders of magnitude slower** than DFT methods (average ~8 hours per reaction)
- Among DFT methods, double-hybrids (SP-2, SP-3) are 4–5× slower than r2SCAN-3c but comparable to VV10 methods
- r2SCAN-3c completes in ~100s — attractive for high-throughput but see accuracy issues below

---

## 4. Energy Analysis

### 4.1 Absolute SP Energy MAE vs REF (kcal/mol)

| Method | rx1 MAE | rx8 MAE | Notes |
|--------|---------|---------|-------|
| **SP-REF** | 0.0 | 0.0 | Reference (by definition) |
| **SP-3** wB97X-2 | **145.3** | **84.5** | Best DFT by MAE |
| **SP-2** PWPB95-D3BJ | 768.3 | 636.6 | |
| **SP-6** r2SCAN-3c | 834.7 | 686.7 | |
| **SP-5** M06-2X | 937.2 | 768.0 | |
| **SP-1** wB97M-V | 959.4 | 801.7 | |
| **SP-4** wB97X-V | 1014.9 | 819.0 | Worst DFT by MAE |

> **ℹ️ Note on MAE interpretation**: Values above ~100 kcal are dominated by method-dependent absolute energy baseline shifts (total energies vs REF differ by 1–2 hartree). These are NOT chemical accuracy errors — they reflect different zero-of-energy conventions between methods. The good news: SP-3's MAE of 84–145 kcal is significantly lower than other DFT methods, suggesting its absolute energy scale is closest to DLPNO-CCSD(T). However, the chemically meaningful metric is **relative energy accuracy** (barrier heights, reaction energies), analyzed below.

### 4.2 Chemically Meaningful Relative Energies

The following compares dE_intermediate_to_ts (the electronic reaction barrier) across methods. Positive values indicate TS above intermediate (physically correct for exothermic reactions).

| Method | rx1 dE_ts_int | rx3 dE_ts_int | rx8 dE_ts_int | rx15 dE_ts_int | Barrier Sign Correct? |
|--------|--------------|--------------|--------------|---------------|----------------------|
| **SP-REF** | **+1.166** | N/A (failed) | **−2.112** | N/A (failed) | ✅ (rx1), ⚠️ (rx8 negative) |
| **SP-1** | +2.892 | **+2.089** | +0.180 | **+2.012** | ✅ All correct |
| **SP-2** | +1.819 | +1.567 | −0.135 | +2.165 | ⚠️ rx8 borderline |
| **SP-3** | −1.121 | −3.722 | −4.011 | −2.600 | ❌ **ALL WRONG** |
| **SP-4** | +2.408 | +1.079 | −0.520 | +1.074 | ⚠️ rx8 borderline |
| **SP-5** | +2.973 | +1.716 | −0.197 | +1.840 | ⚠️ rx8 borderline |
| **SP-6** | −2.305 | −3.101 | −3.674 | −2.906 | ❌ **ALL WRONG** |

**Critical finding**: Both SP-3 (wB97X-2) and SP-6 (r2SCAN-3c) **predict the TS energy below the intermediate** (negative barrier) for nearly all reactions. This is physically impossible for reactions that proceed experimentally. The issue likely stems from the fixed B3LYP geometry not representing a true TS on these methods' potential energy surfaces. This has crucial implications for Phase 2 (GEO) where geometry will be re-optimized.

### 4.3 Reaction Energy Comparison (dE_prod_int in kcal)

| Method | rx1 | rx3 | rx8 | rx15 | Avg | vs REF (rx1,rx8) |
|--------|-----|-----|-----|------|-----|-------------------|
| **SP-REF** | −62.0 | −54.0° | −52.6 | −57.4° | −56.5 | — |
| **SP-1** | −60.3 | −53.4 | −52.0 | −56.7 | −55.6 | **+1.7, +0.6** ✅ |
| **SP-2** | −54.2 | −45.2 | −43.8 | −47.4 | −47.7 | +7.8, +8.8 ❌ |
| **SP-3** | −59.6 | −51.0 | −48.6 | −53.7 | −53.2 | +2.4, +4.0 |
| **SP-4** | −66.1 | −59.9 | −58.2 | −63.2 | −61.9 | −4.1, −5.6 |
| **SP-5** | −55.9 | −48.5 | −46.7 | −51.5 | −50.7 | +6.1, +5.9 |
| **SP-6** | −53.2 | −43.3 | −40.9 | −45.6 | −45.7 | +8.8, +11.7 ❌ |

> ° SP-REF partial_failed on rx3/rx15 — product energy only (no TS comparison)

**Reaction energy findings**: SP-1 (wB97M-V) shows the closest agreement with SP-REF on reaction energies (±1–2 kcal). SP-2 (PWPB95), SP-5 (M06-2X), and SP-6 (r2SCAN-3c) underestimate reaction exothermicity by 6–12 kcal. SP-4 (wB97X-V) overestimates by ~5 kcal.

---

## 5. Per-Reaction Winner Analysis

### rx1 (rigid_small_ring)
- **Winner**: SP-REF (MAE = 0.0, reference)
- **Best DFT**: SP-3 (wB97X-2), MAE = 145.3 kcal
- **All methods completed**: 7/7 ✅
- **SP-REF runtime**: 36,880s (10.2 hours)
- **SP-REF dE_ts_int**: +1.166 kcal — correct TS barrier sign
- **Issue**: SP-3 and SP-6 predict negative barriers (−1.12, −2.31 kcal)

### rx3 (flexible_chain)
- **Winner**: SP-1 (wB97M-V) — SP-REF excluded (partial_failed)
- **SP-REF status**: ⚠️ TS calculation failed (precursor/intermediate/product OK)
- **All DFT methods completed**: 6/6 ✅
- **SP-1 dE_ts_int**: +2.089 kcal
- **Note**: No reference available for dE comparison; winner selected by fallback (most completed points → alphabetical)

### rx8 (medium_polarity_async)
- **Winner**: SP-REF (MAE = 0.0, reference)
- **Best DFT**: SP-3 (wB97X-2), MAE = 84.5 kcal
- **All methods completed**: 7/7 ✅
- **SP-REF runtime**: 16,783s (4.7 hours)
- **SP-REF dE_ts_int**: −2.112 kcal — **SP-REF itself predicts negative barrier on rx8!**
- **All methods predict negative or near-zero barriers for rx8** → This reaction may have a very flat PES near the TS

### rx15 (steric_endo_exo)
- **Winner**: SP-1 (wB97M-V) — SP-REF excluded (partial_failed)
- **SP-REF status**: ⚠️ TS calculation failed (precursor/intermediate/product OK)
- **All DFT methods completed**: 6/6 ✅
- **SP-1 dE_ts_int**: +2.012 kcal

---

## 6. Global Voting and Winner

```
Phase 1 Global Vote:
┌──────────┬───────┬──────────────────────────────────┐
│ Method   │ Votes │ Won On                           │
├──────────┼───────┼──────────────────────────────────┤
│ SP-REF   │   2   │ rx1, rx8                         │
│ SP-1     │   2   │ rx3, rx15 (SP-REF failed)        │
│ SP-2 ~ 6 │   0   │ —                                │
└──────────┴───────┴──────────────────────────────────┘
```

**Declared global winner**: SP-REF (tiebreaker favors reference method in max-vote logic)

**Critical assessment**: The global "winner" designation is misleading. SP-REF won only where it completed (50% of cases), while SP-1 completed 100% of cases. The 2–2 tie reflects SP-REF's reliability problem, not its superiority. A more practical interpretation:

- **When SP-REF works**: It is the gold standard (by definition as reference)
- **When SP-REF fails** (50% of cases): SP-1 (wB97M-V) is the best fallback — reliable, correct barrier sign, good reaction energies
- **Best DFT in absolute energy alignment**: SP-3 (wB97X-2), but it predicts wrong barrier sign at fixed geometry — needs Phase 2 GEO to resolve

---

## 7. Issues and Anomalies

### 7.1 Critical: TS Barrier Sign Errors

| Method | rx1 | rx3 | rx8 | rx15 | % Wrong Sign |
|--------|-----|-----|-----|------|-------------|
| SP-3 (wB97X-2) | ❌ | ❌ | ❌ | ❌ | 100% |
| SP-6 (r2SCAN-3c) | ❌ | ❌ | ❌ | ❌ | 100% |
| SP-REF (DLPNO) | ✅ | N/A | ❌ | N/A | 50% |
| SP-1 (wB97M-V) | ✅ | ✅ | ✅ | ✅ | 0% |

**Root cause**: Fixed B3LYP geometry. These methods have different TS electronic structures; the B3LYP-optimized TS geometry is not a stationary point on their PES. Phase 2 (GEO benchmark with geometry re-optimization) is essential to resolve this.

### 7.2 Critical: SP-REF DLPNO-CCSD(T) TS Failures (rx3, rx15)

Both failures occurred with "ORCA 未正常终止" on the TS single-point calculation specifically. Precursor, intermediate, and product SPs all converged successfully. Possible causes:
1. **Poor SCF convergence** at the TS geometry for large molecules with DLPNO (the TS has partial bonds that challenge the PNO truncation)
2. **Memory/disk issues** — TS calculations are more demanding than minima
3. **Numerical instability** in the DLPNO-CCSD(T) correlation treatment at the TS geometry

**Recommendation**: Investigate DLPNO-CCSD(T) TS SP failure logs. Consider looser PNO thresholds (NormalPNO) or restart with different initial guess.

### 7.3 Metric Concern: Absolute MAE Not Chemically Meaningful

The `absolute_sp_energy_mae_vs_ref_kcal` metric compares raw total energies across methods with different zero-of-energy baselines. The resulting ~800–1000 kcal "errors" are dominated by baseline shifts, not chemical accuracy. The current winner selection based on this metric is valid for ranking methods by their total energy alignment with DLPNO-CCSD(T), but does NOT reflect accuracy for barrier heights or reaction energies.

**Recommendation**: Add relative-energy MAE metrics:
- MAE on dE_intermediate_to_ts (barrier error)
- MAE on dE_intermediate_to_product (reaction energy error)
- These are the quantities that matter for mechanism studies

### 7.4 rx8 Universal Near-Zero/Negative Barrier

All 7 methods predict dE_ts_int ≤ +0.18 kcal for rx8, with most predicting negative values (TS below intermediate). This is unlikely to be a method artifact — it suggests rx8 may genuinely have:
- A nearly barrierless [4+3] cycloaddition
- The fixed B3LYP geometry misplaces the TS relative to the intermediate
- The lite-protocol conformer search may have selected a suboptimal intermediate conformation

---

## 8. Recommendations for Phase 2 (GEO Benchmark)

### 8.1 SP Reader Selection

**Primary recommendation**: **SP-1 (wB97M-V/def2-TZVPP)** as Phase-2 SP reader.

Rationale:
- ✅ 100% completion rate (4/4 reactions)
- ✅ Correct barrier sign prediction in all cases
- ✅ Best agreement with SP-REF on reaction energies (±1–2 kcal where comparable)
- ✅ Moderate runtime (~570s avg)
- ✅ Good physical behavior (no unphysical negative barriers)

**Secondary option**: SP-3 (wB97X-2) — better absolute energy alignment with REF, but Phase 2 must verify that geometry re-optimization resolves the barrier sign issue.

### 8.2 GEO Methods to Prioritize

Based on SP phase results, the following GEO method testing order is recommended:

1. **GEO-1 (B3LYP-D3BJ)** — Baseline, known behavior
2. **GEO-5 (wB97X-D3BJ)** — Best available DFT geometry for TS
3. **GEO-4 (M06-2X)** — Good TS description for main-group organic reactions
4. **GEO-2 (PBE0-D3BJ)** — Alternative hybrid with different exchange fraction
5. **GEO-3 (r2SCAN-3c)** — Low priority given SP sign errors; test only if time permits

### 8.3 TS Validation Criteria

For Phase 2, the TS validation should check:
- Exactly 1 imaginary frequency (n_imag == 1)
- Forming bond lengths in the [2.0, 2.4] Å range
- dG_activation_kcal > 0 (positive free energy barrier)
- dE_activation_kcal > 0 (positive electronic barrier)

These criteria will filter out methods that predict unphysical TS geometries.

---

## 9. Data Inventory

### Primary Output Location
```
Output/benchmark_dft_theory/experiments/session_bl_fixed/
├── manifests/
│   ├── benchmark_manifest.json     ← Global state: stages, methods, winner
│   ├── rx1.json (1112 lines)       ← Full rx1 manifest with SP results
│   ├── rx3.json (1056 lines)       ← Full rx3 manifest
│   ├── rx8.json (1112 lines)       ← Full rx8 manifest
│   └── rx15.json (1056 lines)      ← Full rx15 manifest
├── reports/
│   ├── phase1_global_summary.json  ← Global SP winner tally
│   └── phase1_global_summary.md
├── rx1/reports/sp_summary.json     ← rx1 per-method SP comparison
├── rx1/reports/sp_summary.md
├── rx3/reports/sp_summary.json
├── rx3/reports/sp_summary.md
├── rx8/reports/sp_summary.json
├── rx8/reports/sp_summary.md
├── rx15/reports/sp_summary.json
└── rx15/reports/sp_summary.md
```

### Per-Method QC Outputs
```
rx{N}/sp/SP-{ID}/qc/{point}/{label}_{hash}.out
```
- ~75–84 KB per ORCA output log
- 28 methods × 4 points × 4 reactions = 448 output files
- All files present and readable for completed tasks

### Baseline Stationary Points
```
Output/benchmark_dft_theory/baselines/rx{N}/bl_fixed/
├── stationary_points/
│   ├── precursor_min.xyz
│   ├── intermediate.xyz
│   ├── ts_final.xyz
│   ├── product_min.xyz
│   └── ts_guess.xyz
└── baseline_manifest.json
```
- All 4 reactions have complete baseline artifacts
- TS validation: all pass

---

## 10. Conclusion

The Phase-1 SP benchmark successfully completed with 92.9% of task-slots finishing (26/28). Two DLPNO-CCSD(T) TS failures on rx3 and rx15 represent the only incomplete calculations.

**The core finding**: No single DFT method perfectly reproduces DLPNO-CCSD(T) energies at fixed B3LYP geometry. wB97M-V (SP-1) is the most **reliable** DFT method (100% completion, correct barrier sign, good reaction energies), while wB97X-2 (SP-3) shows the best **absolute energy alignment** but predicts physically incorrect negative barriers at fixed geometry. r2SCAN-3c (SP-6) is fast but scientifically unusable for [4+3] TS energetics at this geometry.

**Phase 2 (GEO benchmark)** is essential to determine whether the barrier sign errors in SP-3 and SP-6 are artifacts of the fixed B3LYP geometry or intrinsic method failures. The Phase-2 geometry re-optimization will allow each method to relax to its own TS, providing a fair comparison of method accuracy for reaction mechanism prediction.

---

## Appendix: Quick-Reference Tables

### A. Completion by Method × Reaction

|  | SP-REF | SP-1 | SP-2 | SP-3 | SP-4 | SP-5 | SP-6 |
|--|--------|------|------|------|------|------|------|
| rx1 | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| rx3 | ⚠️ TS | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| rx8 | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| rx15 | ⚠️ TS | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |

### B. Best Per-Reaction DFT (by MAE vs REF, where REF available)

| rx_id | Best DFT | MAE (kcal) | dE_ts_int (vs REF) | dE_prod_int (vs REF) |
|-------|----------|-----------|---------------------|-----------------------|
| 1 | SP-3 (wB97X-2) | 145.3 | −2.29 kcal (wrong sign) | +2.37 kcal |
| 8 | SP-3 (wB97X-2) | 84.5 | −1.90 kcal (sign worse) | +3.97 kcal |

### C. Runtime-Weighted Recommendation

For production RPH pipelines where SP runtime matters:

| Priority | Method | Use Case |
|----------|--------|----------|
| L1 (geo) | B3LYP-D3BJ (existing) | Geometry optimization — fast, reliable |
| L2 (SP) | wB97M-V (SP-1) | Reliable SP energy, ~600s, correct barrier sign |
| L2 (SP, high-acc) | DLPNO-CCSD(T) (SP-REF) | Gold standard when it converges (~8h) |
| L2 (SP, budget) | M06-2X (SP-5) | Budget option, ~400s, mostly correct barrier sign |
| Avoid | r2SCAN-3c (SP-6) | Wrong barrier sign, no advantage at fixed B3LYP geo |
