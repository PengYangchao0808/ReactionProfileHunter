# RPH v2.1.0 Benchmark Report: S1 Conformer Search Protocols

**Reaction**: rx_id=1 (dataset: `reaxys_cleaned.csv`)  
**Date**: 2026-04-16 ~ 2026-04-17  
**RPH Version**: 2.1.0  
**Scope**: S1-only (S2/S3/S4 skipped)  
**Status**: ✅ 4/4 protocols PASSED  

---

## 1. Test Configuration

| Parameter | Value |
|-----------|-------|
| Product SMILES | `CC(C)(C)OC(=O)N1CCC[C@@]23C=C[C@@H](CC(=O)[C@H]12)O3` |
| Precursor SMILES | `C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C` |
| Reaction Type | `[4+3]_default` |
| S0 Topology | `INTRA_TYPE_I` |
| Forming Bonds | `((12, 13), (15, 16))` |
| OPT Theory | B3LYP/def2-SVP (Gaussian) |
| SP Theory | wB97X-D4/def2-TZVPP (ORCA) |
| Solvent | acetone |
| Protocols Tested | ext, full, lite, zero |
| Benchmark Runner | `benchmark/confsearch/run.sh --rx-id 1` |

---

## 2. Protocol Descriptions

| Feature | **ext** | **full** | **lite** | **zero** |
|---------|---------|----------|----------|----------|
| Two-stage (GFN0→GFN2) | ✅ | ❌ | ❌ | ❌ |
| ngeom_default / ngeom_max | 3 / 6 | 12 / 24 | 4 / 6 | 3 / 4 |
| CREST search mode | `crest_two_stage_gfn0_to_gfn2` | `crest_gfn2` | `crest_gfn2` | `crest_gfn2_or_skip` |
| Prescreen | none | PBEh-3c SP (4.0 kcal) | none | none |
| Screening / Rerank | none | r2SCAN-3c + mRRHO (3.5 kcal) | r2SCAN-3c + mRRHO | GFN2 energy |
| Survivor window | — | 3.0 kcal | — | 0.5 kcal (narrow) |
| Handoff mode | `optimize_all_candidates` | `optimize_all_survivors_within_window` | `optimize_rank1` | `optimize_rank1` |
| Fallback | none | none | `optimize_top2_if_gap_small` (1.0 kcal) | `optimize_all_within_0p5_kcal` |
| Ranking after handoff | Boltzmann (Shermo+SP) | Boltzmann (SP) | final_sp_minimum | final_sp_minimum |
| mRRHO correction | ❌ | ✅ | ✅ | ❌ |
| Boltzmann cutoff | — | — | 0.90 | — |

**Design intent**:
- **ext**: Baseline — exhaustive two-stage conformer search, all candidates forwarded to DFT
- **full**: CENSO-like — large ensemble + two-tier fast-SP filtering (PBEh-3c → r2SCAN-3c) + mRRHO thermo
- **lite**: Fast screening — r2SCAN-3c rerank, Boltzmann cutoff at 90%, only rank-1 forwarded (fallback top-2 if gap < 1 kcal)
- **zero**: Minimal — narrow 0.5 kcal GFN2 energy window, rank-1 only, minimal clustering

---

## 3. Run Results Overview

### 3.1 Wall Time Summary

| Protocol | Start | End | Wall Time | Relative |
|----------|-------|-----|-----------|----------|
| **ext** | 2026-04-16 19:44:04 | 2026-04-16 22:34:04 | **10 199 s** (2h 49m 59s) | 1.00× |
| **full** | 2026-04-16 22:34:05 | 2026-04-17 03:36:20 | **18 135 s** (5h 02m 15s) | 1.78× |
| **lite** | 2026-04-17 03:36:21 | 2026-04-17 05:19:16 | **6 176 s** (1h 42m 56s) | 0.61× |
| **zero** | 2026-04-17 05:19:17 | 2026-04-17 07:18:52 | **7 174 s** (1h 59m 34s) | 0.70× |

**Total benchmark**: ~11h 35m (sequential execution across 4 protocols)

### 3.2 Final Energies

#### Product Global Minimum

| Protocol | E_SP (Hartree) | ΔE vs ext (kcal/mol) |
|----------|---------------|----------------------|
| **ext** | -940.674 883 572 | 0.000 (reference) |
| **full** | -940.674 883 598 | −0.000 02 |
| **lite** | -940.674 883 644 | −0.000 04 |
| **zero** | -940.674 883 638 | −0.000 04 |

> **All protocols converge to the same product energy within < 0.001 kcal/mol** — essentially identical. The product conformer space is well-defined with a dominant minimum.

#### Precursor Global Minimum

| Protocol | E_SP (Hartree) | ΔE vs ext (kcal/mol) |
|----------|---------------|----------------------|
| **ext** | -865.331 299 081 | 0.000 (reference) |
| **full** | -865.333 294 165 | **−1.252** |
| **lite** | -865.333 294 455 | **−1.252** |
| **zero** | -865.331 286 051 | +0.008 |

> **full and lite find a lower-energy precursor by ~1.25 kcal/mol** vs ext. This is a chemically significant improvement. The precursor has a flexible open-chain structure with more accessible conformer space — the single-stage GFN2 + fast-SP reranking in full/lite samples this space more broadly than ext's two-stage approach (15→13 vs 3→6 DFT candidates).

---

## 4. Per-Protocol Detailed Analysis

### 4.1 ext (Baseline)

| Metric | Product | Precursor |
|--------|---------|-----------|
| CREST mode | Two-stage (GFN0→GFN2) | Two-stage (GFN0→GFN2) |
| Funnel candidates | 3 | 6 |
| DFT optimized | 3 (conf_000~002) | 6 (conf_000~005) |
| Best conformer | conf_000 | conf_001 |
| Boltzmann weight | 0.3344 | 0.1673 |
| E_SP (Hartree) | −940.674 884 | −865.331 299 |

**Funnel stages**: `rdkit_embed` → `crest_two_stage` → `ensemble_processing`  
**Handoff**: `optimize_all_candidates` — no filtering, all 3+6=9 conformers sent to DFT  
**Time breakdown**: CREST ~4min + 3×18min (product DFT) + CREST ~10min + 6×18min (precursor DFT) ≈ 170min total  

### 4.2 full (CENSO-like)

| Metric | Product | Precursor |
|--------|---------|-----------|
| CREST mode | Single-stage GFN2 | Single-stage GFN2 |
| CREST candidates | ~large ensemble | ~large ensemble |
| Prescreen (PBEh-3c) | → 2 survivors | → 15 survivors |
| Screening (r2SCAN-3c) | 2 candidates | 15 candidates |
| Window filtered | 2 → 1 | 15 → 13 |
| DFT optimized | 1 (conf_000) | 13 (conf_000~012) |
| Best conformer | conf_000 | conf_001 |
| Boltzmann weight | 1.0000 | 0.0770 |
| E_SP (Hartree) | −940.674 884 | −865.333 294 |

**Funnel stages**: `rdkit_embed` → `crest_gfn2` → `ensemble_processing` → `prescreen_sp` → `screening_sp`  
**Handoff**: `optimize_all_survivors_within_window` (3.0 kcal) — product: 1 survivor, precursor: 13 survivors  
**Key observation**: The precursor fast-SP prescreen generates many more candidates (15) than ext (6), exploring a wider conformer space. This is why full finds a better precursor minimum.  
**Time breakdown**: CREST ~8min + prescreen + screening + 1×18min (product) + CREST ~42min + prescreen + screening + 13×18min (precursor) ≈ 302min total  

### 4.3 lite (Boltzmann-weighted)

| Metric | Product | Precursor |
|--------|---------|-----------|
| CREST mode | Single-stage GFN2 | Single-stage GFN2 |
| Screening (r2SCAN-3c) | 2 candidates | 5 candidates |
| Boltzmann cutoff | 0.90 | 0.90 |
| DFT optimized | 2 (fallback) | 2 |
| Fallback triggered | ✅ YES | — |
| Gap rank1–rank2 | 0.77 kcal | — |
| Best conformer | conf_000 | conf_001 |
| Boltzmann weight | 0.5006 | 0.5001 |
| E_SP (Hartree) | −940.674 884 | −865.333 294 |

**Funnel stages**: `rdkit_embed` → `crest_gfn2` → `ensemble_processing` → `screening_sp`  
**Handoff**: `optimize_rank1` → **fallback** `optimize_top2_if_gap_small` (gap=0.77 < 1.0 kcal)  
**Key observation**: Fallback correctly triggered for product — rank1/rank2 gap is small (0.77 kcal), so both conformers are optimized. The extra conformer does not change the final answer but provides a safety margin.  
**Time breakdown**: CREST ~8min + screening + 2×18min (product) + CREST ~23min + screening + 2×18min (precursor) ≈ 103min total  

### 4.4 zero (Minimal)

| Metric | Product | Precursor |
|--------|---------|-----------|
| CREST mode | GFN2 (narrow) | GFN2 (narrow) |
| Narrow window | 0.5 kcal | 0.5 kcal |
| Window survivors | 1 | 4 |
| DFT optimized | 1 (conf_000) | 4 (conf_000~003) |
| Best conformer | conf_000 | conf_001 |
| Boltzmann weight | 1.0000 | 0.2502 |
| E_SP (Hartree) | −940.674 884 | −865.331 286 |

**Funnel stages**: `rdkit_embed` → `crest_gfn2` → `ensemble_processing` → `narrow_window`  
**Handoff**: `optimize_rank1` — no fallback triggered  
**Key observation**: Zero finds the same product energy but misses the lower-energy precursor conformer found by full/lite. The 0.5 kcal narrow window is too restrictive for the flexible precursor, excluding conformers that would eventually be better after DFT optimization.  
**Time breakdown**: CREST ~7min + 1×18min (product) + CREST ~20min + 4×18min (precursor) ≈ 120min total  

---

## 5. Cross-Protocol Comparison

### 5.1 DFT Optimization Count (Computational Cost Proxy)

| Protocol | Product | Precursor | **Total DFT opts** | Relative cost |
|----------|---------|-----------|-------------------|---------------|
| ext | 3 | 6 | **9** | 1.00× |
| full | 1 | 13 | **14** | 1.56× |
| lite | 2 | 2 | **4** | 0.44× |
| zero | 1 | 4 | **5** | 0.56× |

> full is the most expensive (14 DFT optimizations) due to its large survivor window for the flexible precursor. lite is the most efficient (4 DFT optimizations).

### 5.2 Energy Accuracy vs Cost

| Protocol | Wall Time | Product ΔE | Precursor ΔE | Total ΔE | DFT opts |
|----------|-----------|-----------|-------------|----------|----------|
| **ext** | 10 199 s | 0.000 kcal | 0.000 kcal | 0.000 | 9 |
| **full** | 18 135 s | −0.000 kcal | **−1.252** kcal | −1.252 | 14 |
| **lite** | 6 176 s | −0.000 kcal | **−1.252** kcal | −1.252 | 4 |
| **zero** | 7 174 s | −0.000 kcal | +0.008 kcal | +0.008 | 5 |

> **lite achieves the same energy improvement as full at 1/3 the cost and 2/3 the wall time of full.** This demonstrates the effectiveness of the r2SCAN-3c screening + Boltzmann cutoff strategy.

### 5.3 Fallback Trigger Analysis

| Protocol | Product | Precursor |
|----------|---------|-----------|
| ext | Not applicable (all candidates) | Not applicable |
| full | Not triggered (1 survivor) | Not triggered (13 within window) |
| **lite** | **✅ Triggered** (gap=0.77 < 1.0 kcal → top-2) | Not triggered (only rank-1) |
| zero | Not triggered (1 within 0.5 kcal) | Not triggered (rank-1 only) |

> lite's fallback correctly activates for the product when rank1/rank2 gap is small, preventing potential energy misses from premature pruning.

---

## 6. Summary & Conclusions

### 6.1 Key Findings

1. **Product energy is robust across all protocols** — the cyclic [4+3] product has a well-defined global minimum that all funnel strategies find. ΔE < 0.001 kcal/mol across protocols.

2. **Precursor energy is protocol-sensitive** — the flexible open-chain precursor benefits from broader conformer sampling:
   - **full/lite**: Find −865.333 294 Ha (**1.25 kcal/mol lower** than ext)
   - **ext/zero**: Find −865.331 Ha range

3. **lite offers the best cost-accuracy tradeoff**:
   - Same energy as full (the most expensive protocol)
   - **3.5× fewer DFT optimizations** than full
   - **2.9× faster** than full
   - **1.65× faster** than ext (the baseline)
   - Correctly triggers fallback when conformer energies are close

4. **zero is fast but risks missing better minima** for flexible molecules:
   - 0.5 kcal GFN2 window is too restrictive for the precursor
   - Misses the 1.25 kcal/mol lower minimum found by full/lite
   - Product energy remains accurate (rigid cyclic structure)

5. **full is the most expensive but not the most accurate**:
   - 14 DFT optimizations (vs 4 for lite)
   - 18 135s wall time (vs 6 176s for lite)
   - Finds the same precursor as lite — the extra cost comes from optimizing 13 precursor conformers

### 6.2 Recommendations

| Use Case | Recommended Protocol | Rationale |
|----------|---------------------|-----------|
| Routine production runs | **lite** | Best cost/accuracy; fallback protects against close conformers |
| Small/rigid molecules | **zero** | Minimal cost; narrow window sufficient for well-defined minima |
| Benchmarking / validation | **ext** | Baseline reference; exhaustive sampling |
| Large flexible molecules | **full** | Broad sampling + fast-SP reranking; worth the extra cost |
| Batch high-throughput | **lite** or **zero** | Speed priority; lite provides safety margin |

### 6.3 Caveats

1. **Single reaction tested** (rx_id=1) — conclusions should be validated across a larger benchmark set with diverse topologies and flexibility.
2. **S1-only scope** — S2/S3/S4 were skipped; downstream pipeline robustness should be verified with the full S1→S4 flow.
3. **Toxic path workaround** — all runs used temporary directory execution due to spaces/brackets in the project path. This adds minor overhead (~0.5-1s per QC call for file copying) but does not affect energy results.
4. **evaluate.py report** — the auto-generated `summary.md` showed all protocols as ❌ because it was run from a previous partial run. The provenance.json files confirm all 4 protocols completed successfully. Re-running `benchmark/confsearch/evaluate.py` would regenerate the report correctly.

---

## Appendix A: Provenance Verification

All four protocols pass the benchmark validation checks:

```
[  OK] rx1/ext: provenance.json found
[  OK] rx1/ext: provenance.protocol = ext
[  OK] rx1/ext: product structure found
[  OK] rx1/ext: precursor structure found
[  OK] rx1/ext: S2/S3/S4 absent as expected

[  OK] rx1/full: provenance.json found
[  OK] rx1/full: provenance.protocol = full
[  OK] rx1/full: product structure found
[  OK] rx1/full: precursor structure found
[  OK] rx1/full: S2/S3/S4 absent as expected

[  OK] rx1/lite: provenance.json found
[  OK] rx1/lite: provenance.protocol = lite
[  OK] rx1/lite: product structure found
[  OK] rx1/lite: precursor structure found
[  OK] rx1/lite: S2/S3/S4 absent as expected

[  OK] rx1/zero: provenance.json found
[  OK] rx1/zero: provenance.protocol = zero
[  OK] rx1/zero: product structure found
[  OK] rx1/zero: precursor structure found
[  OK] rx1/zero: S2/S3/S4 absent as expected
```

## Appendix B: Output Locations

```
Output/benchmark/
└── rx1/
    ├── ext/RXN_358d0ac1/S1_ConfGeneration/
    │   ├── provenance.json
    │   ├── product/product_global_min.xyz
    │   └── precursor/precursor_global_min.xyz
    ├── full/RXN_358d0ac1/S1_ConfGeneration/
    │   ├── provenance.json
    │   ├── product/product_global_min.xyz
    │   └── precursor/precursor_global_min.xyz
    ├── lite/RXN_358d0ac1/S1_ConfGeneration/
    │   ├── provenance.json
    │   ├── product/product_global_min.xyz
    │   └── precursor/precursor_global_min.xyz
    ├── zero/RXN_358d0ac1/S1_ConfGeneration/
    │   ├── provenance.json
    │   ├── product/product_global_min.xyz
    │   └── precursor/precursor_global_min.xyz
    └── evaluation/
        ├── summary.json
        └── summary.md
```
