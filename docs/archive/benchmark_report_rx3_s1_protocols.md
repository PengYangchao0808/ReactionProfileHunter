# RPH v2.1.0 Benchmark Report: rx3 S1 Conformer Search Protocols

**Reaction**: rx_id=3 (dataset: `reaxys_cleaned.csv`)  
**Date**: 2026-04-17 ~ 2026-04-18  
**RPH Version**: 2.1.0  
**Scope**: S1-only (S2/S3/S4 skipped)  
**Status**: ✅ 4/4 protocols PASSED  

---

## 1. Test Configuration

| Parameter | Value |
|-----------|-------|
| Product SMILES | `O=C1C[C@@H]2C=C[C@@]3(CCN(C(=O)OCc4ccccc4)[C@@H]13)O2` |
| Precursor SMILES | `C=C=CN(CCc1ccfo1)C(=O)OCc1ccccc1` |
| Reaction Type | `[4+3]_default` |
| S0 Classification | `INTRA_TYPE_I`, Forming Bonds `((12, 13), (15, 16))` |
| Map Status | `LOW_CONFIDENCE` (0.785, sanity check passed) |
| OPT Theory | B3LYP/def2-SVP (Gaussian) |
| SP Theory | wB97X-D4/def2-TZVPP (ORCA) |
| Solvent | acetone |
| Protocols Tested | ext, full, lite, zero |

**与 rx1 的结构差异**：rx3 使用 benzyl carbamate 保护基（—OCc₆H₅）替代 rx1 的 Boc（—OC(C)(₃)C），同时 precursor tether 少一个 —CH₂—（—CCc₆H₅ vs —CCCc₆H₅）。产物骨架相同但取代基不同，构象空间特征有显著差异。

---

## 2. Wall Time Summary

| Protocol | Start | End | Wall Time | Relative (vs ext) |
|----------|-------|-----|-----------|-------------------|
| **ext** | 2026-04-17 13:05 | 2026-04-17 15:06 | **7 294 s** (2h 02m) | 1.00× |
| **full** | 2026-04-18 01:06 | 2026-04-18 11:27 | **37 316 s** (10h 22m) | 5.12× |
| **lite** | 2026-04-18 12:35 | 2026-04-18 13:12 | **2 248 s** (0h 37m) | 0.31× |
| **zero** | 2026-04-18 14:55 | 2026-04-18 16:05 | **4 197 s** (1h 10m) | 0.58× |

> rx3 的 full 协议比 rx1 慢得多（10h22m vs 5h02m），因为 CREST 为 benzyl 取代基生成了更多构象（product 16 + precursor 20 = 36 个 DFT 优化）。

---

## 3. Final Energies

### 3.1 Product Global Minimum

| Protocol | E_SP (Hartree) | ΔE vs ext (kcal/mol) | DFT conformers |
|----------|---------------|----------------------|----------------|
| **ext** | −1014.502 037 830 | **0.000** (reference) | 4 |
| **full** | −1014.501 490 937 | **+0.343** | 16 |
| **lite** | −1014.503 440 281 | **−0.880** ✨ | 2 |
| **zero** | −1014.503 438 957 | **−0.879** ✨ | 4 |

> ⚠️ **重大发现：lite 和 zero 找到了比 ext 基准更低 0.88 kcal/mol 的产物全局最低点！**
>
> 这与 rx1（所有协议收敛到相同能量，ΔE < 0.001 kcal/mol）形成鲜明对比。原因：rx3 的 benzyl 取代基引入额外自由度，xTB 能量排序与 DFT 排序不一致——ext 的 GFN0→GFN2 两阶段采样 + xTB 排名将真正的全局最低点遗漏在外，而 lite/zero 的 r2SCAN-3c 或 GFN2 energy 重排能更好地识别低能构象。
>
> full 虽然 DFT 优化了 16 个产物构象，但最终找到的能量比 ext 还高 0.34 kcal/mol，说明 full 的 screening 阶段虽然保留了 17 个候选者，但仍未包含 lite 发现的那个最优构象。

### 3.2 Precursor Global Minimum

| Protocol | E_SP (Hartree) | ΔE vs ext (kcal/mol) | DFT conformers |
|----------|---------------|----------------------|----------------|
| **ext** | −939.163 942 867 | **0.000** (reference) | 6 |
| **full** | −939.165 231 673 | **−0.809** | 20 |
| **lite** | −939.165 230 273 | **−0.808** | 2 |
| **zero** | −939.163 350 620 | **+0.372** | 4 |

> 与 rx1 一致：full 和 lite 在 precursor 上比 ext 低约 0.81 kcal/mol。zero 仍然偏高（+0.37 kcal/mol），0.5 kcal/mol 窄窗对柔性 precursor 仍然不够宽。

---

## 4. Per-Protocol Detailed Analysis

### 4.1 ext (Baseline)

| Metric | Product | Precursor |
|--------|---------|-----------|
| CREST mode | Two-stage (GFN0→GFN2) | Two-stage (GFN0→GFN2) |
| Funnel candidates | 4 | 6 |
| DFT optimized | 4 (conf_000~003) | 6 (conf_000~005) |
| Best conformer | conf_000 (weight=0.2504) | — |
| Gap rank1–rank2 (xTB) | **0.025 kcal/mol** | — |
| E_SP (Hartree) | −1014.502 038 | −939.163 943 |

**关键特征**：产物 xTB 能隙极小（0.025 kcal/mol），4 个候选者几乎简并。这意味着 xTB 的排序置信度很低——真正的 DFT 全局最低点可能不在这 4 个之中。

### 4.2 full (CENSO-like)

| Metric | Product | Precursor |
|--------|---------|-----------|
| CREST mode | Single-stage GFN2 | Single-stage GFN2 |
| Prescreen (PBEh-3c) | → large ensemble | → large ensemble |
| Screening (r2SCAN-3c) | 17 candidates | 20 candidates |
| Window (3.0 kcal) | 17 survivors | 20 survivors |
| DFT optimized | **16** | **20** |
| Best conformer | — | — |
| Gap rank1–rank2 (screening) | 0.656 kcal/mol | — |
| E_SP (Hartree) | −1014.501 491 | −939.165 232 |

**关键特征**：full 的 3.0 kcal survivor 窗口将所有 17 个候选者全部送入 DFT，没有起到筛选作用。尽管 DFT 优化了 16 个产物构象，但仍未找到 lite 发现的最优构象。这说明 CREST 的 GFN2 初始采样空间就遗漏了该构象（而 lite 的 screening 阶段恰好捕获了它）。

### 4.3 lite (Boltzmann-weighted)

| Metric | Product | Precursor |
|--------|---------|-----------|
| CREST mode | Single-stage GFN2 | Single-stage GFN2 |
| Screening (r2SCAN-3c) | 5 candidates | 2 candidates |
| Boltzmann cutoff | 0.90 | — |
| Fallback triggered | ✅ **YES** (gap=0.794 < 1.0 kcal) | — |
| DFT optimized | **2** | **2** |
| Best conformer | cand_002 | — |
| Gap rank1–rank2 (screening) | 0.794 kcal/mol | — |
| E_SP (Hartree) | **−1014.503 440** ✨ | **−939.165 230** |

**关键特征**：
1. Fallback 正确触发（rank1/rank2 gap=0.794 < 1.0 kcal），将 rank-1 和 rank-2 都送入 DFT。
2. **仅用 2 个 DFT 优化就找到了所有协议中的最优产物能量**。r2SCAN-3c screening 的重排能力在此发挥了关键作用。
3. Precursor 也以仅 2 个 DFT 优化达到了与 full（20 个 DFT）几乎相同的能量（Δ = 0.001 kcal/mol）。

### 4.4 zero (Minimal)

| Metric | Product | Precursor |
|--------|---------|-----------|
| CREST mode | GFN2 (narrow) | GFN2 (narrow) |
| Narrow window | 0.5 kcal | 0.5 kcal |
| Window survivors | **4** (fallback!) | **4** (fallback!) |
| Fallback triggered | ✅ **YES** (all within 0.5 kcal) | ✅ **YES** (all within 0.5 kcal) |
| DFT optimized | 4 | 4 |
| Gap rank1–rank2 (GFN2) | **0.017 kcal/mol** | — |
| E_SP (Hartree) | −1014.503 439 | −939.163 351 |

**关键特征**：
1. 产物 GFN2 能隙极小（0.017 kcal/mol），4 个候选者几乎完全简并 → fallback 将全部 4 个送入 DFT。
2. 产物能量与 lite 几乎相同（ΔE = 0.001 kcal/mol），说明 GFN2 排序在这个特定分子上恰好是正确的。
3. 但 precursor 的 zero 仍然比 full/lite 差 1.18 kcal/mol，说明窄窗对柔性前体构象采样不足。

---

## 5. Cross-Protocol Comparison

### 5.1 DFT Optimization Count (Computational Cost)

| Protocol | Product | Precursor | **Total DFT** | Relative cost |
|----------|---------|-----------|--------------|---------------|
| ext | 4 | 6 | **10** | 1.00× |
| full | 16 | 20 | **36** | 3.60× |
| lite | 2 | 2 | **4** | 0.40× |
| zero | 4 | 4 | **8** | 0.80× |

### 5.2 Energy Accuracy vs Cost

| Protocol | Wall Time | Product ΔE | Precursor ΔE | Total ΔE | DFT opts |
|----------|-----------|-----------|-------------|----------|----------|
| **ext** | 7 294 s | 0.000 kcal | 0.000 kcal | 0.000 | 10 |
| **full** | 37 316 s | +0.343 kcal | **−0.809** kcal | −0.466 | 36 |
| **lite** | 2 248 s | **−0.880** kcal | **−0.808** kcal | **−1.688** | 4 |
| **zero** | 4 197 s | **−0.879** kcal | +0.372 kcal | −0.507 | 8 |

> **lite 在所有指标上全面胜出**：最低的产物能量（−0.880 kcal）、几乎最低的 precursor 能量（−0.808 kcal）、最少的 DFT 优化数（4 次）、最快的壁钟时间（37 分钟）。

### 5.3 Fallback Trigger Summary

| Protocol | Product | Precursor |
|----------|---------|-----------|
| ext | N/A (all candidates) | N/A |
| full | Not triggered (all within 3 kcal) | Not triggered |
| **lite** | **✅ Triggered** (gap=0.794 < 1.0 → top-2) | Not triggered |
| **zero** | **✅ Triggered** (gap=0.017, all within 0.5 kcal) | **✅ Triggered** (all within 0.5 kcal) |

---

## 6. rx1 vs rx3 对比分析

| 维度 | rx1 (Boc) | rx3 (Bn) |
|------|-----------|----------|
| 产物 ΔE (ext vs lite) | < 0.001 kcal/mol | **−0.880 kcal/mol** |
| 前体 ΔE (ext vs lite) | −1.252 kcal/mol | −0.808 kcal/mol |
| ext 产物 DFT 数 | 3 | 4 |
| full 产物 DFT 数 | 1 | **16** |
| lite 产物 DFT 数 | 2 | 2 |
| lite wall time | 1h 43m | **37 min** |
| full wall time | 5h 02m | **10h 22m** |
| ext product xTB gap | 2.17 kcal/mol | **0.025 kcal/mol** |
| lite product fallback | ✅ (gap=0.77) | ✅ (gap=0.79) |

**核心差异**：
1. **rx1 产物构象空间简单**（Boc 刚性大），所有协议收敛到同一能量；**rx3 产物构象空间复杂**（Bn 引入额外旋转自由度），协议间出现显著能量差异。
2. rx3 的 full 协议成本暴增（36 个 DFT vs rx1 的 14 个），因为 CREST 对 benzyl 取代基生成了大量低能构象。
3. **rx3 中 xTB 排序质量明显下降**：ext 的 xTB gap 仅 0.025 kcal/mol，说明 xTB 无法可靠区分 rx3 产物构象。这解释了为什么 lite 的 r2SCAN-3c 重排能找到更好的构象——DFT 级别的重排对 benzyl 取代基的构象排序至关重要。

---

## 7. Summary & Conclusions

### 7.1 Key Findings

1. **lite 全面优于 ext 基准**：产物能量低 0.88 kcal/mol，前体能量低 0.81 kcal/mol，DFT 优化数仅 4 次（ext 的 40%），壁钟时间仅 37 分钟（ext 的 31%）。

2. **ext 基准在 rx3 上存在盲区**：xTB 两阶段采样 + xTB 排名遗漏了产物最优构象。对于含 benzyl 等较柔性取代基的分子，xTB 排序的可靠性下降。

3. **full 协议性价比最差**：36 个 DFT 优化、10 小时壁钟时间，产物能量反而比 ext 高 0.34 kcal/mol。3.0 kcal survivor 窗口在此案例中未起到有效筛选作用。

4. **zero 对产物表现意外地好**（−0.879 kcal/mol），但对前体仍然不够（+0.37 kcal/mol）。产物 GFN2 能隙极小导致 fallback 将全部 4 个候选送入 DFT，无意中覆盖了最优构象。

5. **Fallback 机制在 rx3 中频繁触发**：lite 和 zero 都触发了产物 fallback，说明 rx3 的构象空间存在多个近简并构象。

### 7.2 rx1+rx3 综合建议

| Use Case | Recommended Protocol | Rationale |
|----------|---------------------|-----------|
| **一般生产运行** | **lite** | 2/2 反应中均为最优或并列最优；成本最低 |
| 小分子 / 刚性骨架 | lite 或 zero | 构象空间简单时 zero 足够 |
| 大柔性分子验证 | lite + ext 对照 | ext 可作为备用参考，但不应作为唯一依赖 |
| 高通量筛选 | **lite** | 壁钟时间最短（37 min），能量最优 |
| ❌ 不推荐 | full | 成本最高，能量不是最优 |

### 7.3 Caveats

1. **2 个反应的结论需更多数据验证**：rx3 展现了与 rx1 不同的行为模式（协议间能量差异显著），但仅凭 2 个样本不足以得出普适结论。rx8（刚性环状）和 rx15（高柔性长链）将是关键的验证案例。

2. **S1-only scope**：S2/S3/S4 被跳过；下游管线鲁棒性需要全流程验证。

3. **Path workaround**：所有 ORCA 运行使用临时目录执行（"检测到路径包含空格/特殊字符"），每次 QC 调用有 ~0.5-1s 的文件复制开销，不影响能量结果。

---

## Appendix A: Provenance Verification

```
rx3/ext:  ✅ provenance.json found, protocol = ext
rx3/full: ✅ provenance.json found, protocol = full
rx3/lite: ✅ provenance.json found, protocol = lite
rx3/zero: ✅ provenance.json found, protocol = zero
```

All 4 protocols: product structure found, precursor structure found, S2/S3/S4 absent as expected.

## Appendix B: Output Locations

```
Output/benchmark/rx3/
├── ext/RXN_8c098878/S1_ConfGeneration/
│   ├── provenance.json
│   ├── product/  (4 DFT conformers)
│   └── precursor/  (6 DFT conformers)
├── full/RXN_8c098878/S1_ConfGeneration/
│   ├── provenance.json
│   ├── product/  (16 DFT conformers)
│   └── precursor/  (20 DFT conformers)
├── lite/RXN_8c098878/S1_ConfGeneration/
│   ├── provenance.json
│   ├── product/  (2 DFT conformers)
│   └── precursor/  (2 DFT conformers)
└── zero/RXN_8c098878/S1_ConfGeneration/
    ├── provenance.json
    ├── product/  (4 DFT conformers)
    └── precursor/  (4 DFT conformers)
```
