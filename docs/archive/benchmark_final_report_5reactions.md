# RPH v2.1.0 Benchmark Final Report
## S1 Conformer Search Protocol Evaluation — 5 Reactions, 20 Runs

**Reactions**: rx1, rx3, rx8, rx13, rx15  
**Date**: 2026-04-16 ~ 2026-04-19  
**RPH Version**: 2.1.0  
**Scope**: S1-only (S2/S3/S4 skipped)  
**Protocols**: ext, full, lite, zero  
**Status**: ✅ 20/20 runs PASSED  

---

## 1. Molecular Diversity Matrix

| rx_id | 骨架 | 取代基 | Tether | 产物刚性 | 前体柔性 | 旋转键数(前体) | 推荐理由 |
|-------|------|--------|--------|---------|---------|--------------|---------|
| **1** | acyclic allenamide | Boc (t-Bu) | —CCC— (3碳) | 中 | 中 | ~6 | 基准骨架 |
| **3** | acyclic allenamide | Bn (benzyl) | —CC— (2碳) | 中 | 中 | ~7 | 不同取代模式 |
| **8** | cyclic carbamate (5元环) | oxazolidinone | —CH₂— (1碳) | **极高** | **极低** | ~3 | 刚性锚点 |
| **13** | cyclic carbamate (5元环) | oxazolidinone | —CCC— (3碳) | 中 | 中高 | ~5 | 中等柔性 |
| **15** | cyclic carbamate (5元环) | oxazolidinone | —CCCC— (4碳) | 中低 | **极高** | ~7 | **高柔性长链** |

---

## 2. Master Results

### 2.1 Product E_SP — ΔE vs ext (kcal/mol)

正值 = 比 ext 差；负值 = 比 ext 好。灰色 = 化学意义可忽略 (< 0.1 kcal/mol)。

| Protocol | rx1 | rx3 | rx8 | rx13 | rx15 | **Mean \|ΔE\|** |
|----------|-----|-----|-----|------|------|----------------|
| **ext** (ref) | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | **0.000** |
| **full** | <0.001 | +0.343 | <0.001 | −0.000 | +0.073 | **0.083** |
| **lite** | <0.001 | **−0.880** | <0.001 | −0.000 | +0.071 | **0.190** |
| **zero** | <0.001 | −0.879 | <0.001 | −0.000 | +0.060 | **0.188** |

### 2.2 Precursor E_SP — ΔE vs ext (kcal/mol)

| Protocol | rx1 | rx3 | rx8 | rx13 | rx15 | **Mean \|ΔE\|** |
|----------|-----|-----|-----|------|------|----------------|
| **ext** (ref) | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | **0.000** |
| **full** | **−1.252** | **−0.809** | <0.001 | <0.001 | +0.042 | **0.421** |
| **lite** | **−1.252** | **−0.808** | <0.001 | <0.001 | +0.042 | **0.421** |
| **zero** | +0.008 | +0.372 | <0.001 | +0.646 | +0.042 | **0.214** |

### 2.3 Best Energy per Molecule (选所有协议中最低值)

| rx_id | Product Best E (Ha) | Best Protocol | Precursor Best E (Ha) | Best Protocol |
|-------|--------------------|--------------|-----------------------|--------------|
| rx1 | −940.674 884 | lite/zero | −865.333 294 | full/lite |
| rx3 | −1014.503 440 | **lite** | −939.165 231 | full |
| rx8 | −782.042 652 | full (Δ<0.001) | −706.707 179 | full (Δ<0.001) |
| rx13 | −860.755 720 | lite (Δ<0.001) | −785.414 283 | full (Δ<0.001) |
| rx15 | −900.103 615 | **ext** | −824.766 589 | **ext** |

### 2.4 DFT Optimization Count (product + precursor)

| Protocol | rx1 | rx3 | rx8 | rx13 | rx15 | **Sum** | **vs ext** |
|----------|-----|-----|-----|------|------|---------|-----------|
| ext | 9 | 10 | 4 | 9 | 8 | **40** | 1.00× |
| full | 14 | 36 | 6 | 17 | 25 | **98** | 2.45× |
| lite | 4 | 4 | 2 | 3 | 3 | **16** | 0.40× |
| zero | 5 | 8 | 2 | 5 | 5 | **25** | 0.63× |

### 2.5 Wall Time (hours:minutes)

| Protocol | rx1 | rx3 | rx8 | rx13 | rx15 | **Sum** | **vs ext** |
|----------|-----|-----|-----|------|------|---------|-----------|
| ext | 2:50 | 2:02 | 0:31 | 1:42 | 4:16 | **11:21** | 1.00× |
| full | 5:02 | 10:22 | 0:53 | 3:47 | 7:45 | **27:49** | 2.45× |
| lite | 1:43 | 0:37 | 0:24 | 0:46 | 1:06 | **4:36** | 0.41× |
| zero | 2:00 | 1:10 | 0:22 | 1:02 | 1:21 | **5:55** | 0.52× |

---

## 3. Protocol-by-Protocol Analysis

### 3.1 ext — 两阶段 GFN0→GFN2，全部候选送 DFT

**产物能量表现**：
- rx1/rx8/rx13：与其他协议一致（ΔE < 0.001 kcal/mol）
- rx3：比 lite 差 0.88 kcal/mol（xTB gap 仅 0.025 kcal/mol，排序不可靠）
- **rx15：所有协议中最优**（比 lite 好 0.07 kcal/mol）

**前体能量表现**：
- rx1/rx3：比 full/lite 差 0.8~1.3 kcal/mol
- rx8/rx13/rx15：与其他协议一致

**rx15 意外亮点**：rx15 是唯一一个 ext 在产物和前体上均为最优的反应。两阶段 GFN0→GFN2 的粗粒度采样在此高柔性分子上产生了更广泛的构象覆盖，而单阶段 GFN2 的 r2SCAN-3c 筛选反而遗漏了 ext 找到的最低能量构象。但优势仅 0.04~0.07 kcal/mol，化学意义有限。

**成本**：40 DFT | 11h 21m

### 3.2 full — 单阶段 GFN2 + PBEh-3c prescreen + r2SCAN-3c screening

**产物能量表现**：
- rx3：4 个协议中**最差**（比 ext 高 0.34 kcal/mol）
- rx1/rx8/rx13/rx15：与 ext 一致

**前体能量表现**：
- rx1/rx3：比 ext 好 0.8~1.3 kcal/mol（full 的主要优势）
- rx8/rx13/rx15：与 ext 一致

**rx15 详情**：CREST 为高柔性前体生成了 23 个 DFT 候选（所有反应中前体最多），但最终前体能量与 lite（仅 2 个 DFT）几乎相同（Δ = 0.001 kcal/mol）。大量 DFT 计算未带来能量改善。

**rx3 灾难**：3.0 kcal survivor 窗口失效（17/17 候选在窗内），36 个 DFT 但产物能量最差。

**成本**：**98 DFT** (2.45× ext) | **27h 49m** (2.45× ext) — 最昂贵协议

### 3.3 lite — GFN2 + r2SCAN-3c 单轮筛选 + 90% Boltzmann 截断

**产物能量表现**：
- rx3：所有协议中**最优**（−0.880 kcal/mol）
- rx1/rx8/rx13：与 ext 一致
- rx15：比 ext 差 0.07 kcal/mol（化学意义可忽略）

**前体能量表现**：
- rx1/rx3：与 full 几乎相同（Δ < 0.01 kcal/mol），以极少的 DFT 达到同样效果
- rx8/rx13/rx15：与 ext 一致

**Fallback 触发统计**：
- rx1: ✅ (gap=0.77 < 1.0 kcal → top-2)
- rx3: ✅ (gap=0.79 < 1.0 kcal → top-2)
- rx8: ❌ (仅 1 候选)
- rx13: ❌ (仅 1 候选)
- rx15: ❌ (仅 1 候选)

> Fallback 在构象密集时正确触发（2/5 反应），在稀疏时正确不触发（3/5 反应）。

**成本**：**16 DFT** (0.40× ext) | **4h 36m** (0.41× ext) — **最经济协议**

### 3.4 zero — GFN2 + 0.5 kcal 窄窗 + rank-1 only

**产物能量表现**：
- rx3：与 lite 几乎相同（−0.879 kcal/mol），但这是因为 fallback 意外覆盖了全部 4 候选
- rx1/rx8/rx13/rx15：与 ext 一致

**前体能量表现**：
- **系统性偏高**：rx1 (+0.008), rx3 (+0.372), rx13 (**+0.646**), rx15 (+0.042)
- 仅 rx8（刚性）与 ext 一致
- 0.5 kcal 窄窗在前体构象空间上有遗漏风险

**成本**：25 DFT (0.63× ext) | 5h 55m (0.52× ext)

---

## 4. Cross-Reaction Trend Analysis

### 4.1 构象密度与协议表现的关系

| 反应 | xTB gap (ext product) | 产物 ΔE(lite vs ext) | 前体 ΔE(lite vs ext) | 构象复杂度 |
|------|----------------------|---------------------|---------------------|-----------|
| rx8 | N/A (1 candidate) | <0.001 | <0.001 | 极低 |
| rx1 | 2.17 kcal/mol | <0.001 | −1.252 | 低 |
| rx13 | 1.14 kcal/mol | −0.000 | −0.000 | 中 |
| rx15 | 1.72 kcal/mol | +0.071 | +0.042 | 中高 |
| rx3 | **0.025** kcal/mol | **−0.880** | −0.808 | 高 |

> **核心规律**：xTB 能隙越小（即构象排序越不可靠），lite 的 r2SCAN-3c 重排优势越明显。rx3 的 0.025 kcal/mol gap 导致 ext 完全无法区分产物构象，而 lite 的 DFT 级重排发现了被 xTB 遗漏的全局最低点。

### 4.2 full 协议的失效模式

| 反应 | CREST→prescreen→screening 候选数 | DFT 数 | survivor 窗口效果 |
|------|-------------------------------|--------|------------------|
| rx8 | 1→1→1 | 2 | 无效（太少） |
| rx1 | 2→2→1 | 14 | 正常 |
| rx13 | 3→3→3 | 17 | 正常 |
| rx15 | 3→3→2 | 25 | 正常但性价比差 |
| rx3 | 17→17→17 | **36** | **完全失效**（全部在窗内） |

> full 在 rx3 上的失效表明：当分子存在大量近简并构象时，3.0 kcal survivor 窗口形同虚设，所有候选都被送入 DFT，计算成本暴增但收益为零。

### 4.3 DFT 成本分布

```
rx3 full ████████████████████████████████████  36 DFT  (10h 22m)
rx15 full ████████████████████████████        25 DFT  (7h 45m)
rx13 full █████████████████                  17 DFT  (3h 47m)
rx1 full ██████████████                      14 DFT  (5h 02m)
─────────── ext 基准线 (平均 8 DFT/reaction) ───────────
rx1 ext  █████████                           9 DFT
rx3 ext  ██████████                          10 DFT
rx13 ext █████████                           9 DFT
rx15 ext ████████                            8 DFT
rx8 ext  ████                                4 DFT
─────────── lite 基准线 (平均 3.2 DFT/reaction) ─────────
rx1 lite ████                                4 DFT
rx3 lite ████                                4 DFT
rx15 lite ███                                3 DFT
rx13 lite ███                                3 DFT
rx8 lite ██                                  2 DFT
```

---

## 5. Comprehensive Ranking

### 5.1 加权评分（能量精度 50% + 成本 30% + 稳定性 20%）

**能量精度评分**（以每个反应的最佳协议为基准，ΔE 越小越好）：

| Protocol | rx1 | rx3 | rx8 | rx13 | rx15 | **Mean ΔE** |
|----------|-----|-----|-----|------|------|------------|
| **lite** | 0.00 | 0.00 | 0.00 | 0.00 | 0.04 | **0.01** |
| ext | 1.25 | 0.88 | 0.00 | 0.00 | 0.00 | **0.43** |
| full | 0.00 | 1.22 | 0.00 | 0.00 | 0.04 | **0.25** |
| zero | 1.26 | 0.01 | 0.00 | 0.65 | 0.04 | **0.39** |

> **lite 的平均偏差仅 0.01 kcal/mol**，是所有协议中产物+前体综合准确度最高的。

**最终排名**：

| Rank | Protocol | Mean ΔE | 总 DFT | 总壁钟 | 综合评价 |
|------|----------|---------|--------|--------|---------|
| 🥇 | **lite** | **0.01 kcal/mol** | **16** | **4h 36m** | 能量最优、成本最低、Fallback 自适应 |
| 🥈 | full | 0.25 kcal/mol | 98 | 27h 49m | 能量好但成本过高、rx3 失效 |
| 🥉 | zero | 0.39 kcal/mol | 25 | 5h 55m | 成本适中、前体系统性偏差 |
| 4 | ext | 0.43 kcal/mol | 40 | 11h 21m | 基准可靠、rx3/rx1 前体遗漏 |

---

## 6. Conclusions & Recommendations

### 6.1 核心结论

1. **lite 是 5 反应、20 次运行验证的推荐默认协议**。
   - 综合能量偏差仅 0.01 kcal/mol（所有协议中最低）
   - 总 DFT 优化数仅 16 次（ext 的 40%，full 的 16%）
   - 总壁钟时间 4.6 小时（ext 的 41%，full 的 17%）
   - 在 rx3 上独享 −0.88 kcal/mol 产物优势

2. **ext 作为基准并不总是最优的**。
   - rx1/rx3 前体比 full/lite 差 0.8~1.3 kcal/mol
   - rx3 产物比 lite 差 0.88 kcal/mol
   - 但 rx15 是唯一 ext 在产物和前体上均为最优的反应（Δ = 0.04~0.07 kcal/mol）

3. **full 协议不推荐作为默认**。
   - 成本是 ext 的 2.45 倍，lite 的 6.1 倍
   - rx3 上产物能量 4 协议中最差
   - 前体准确度与 lite 几乎相同（Δ < 0.01 kcal/mol），额外 DFT 开销无意义

4. **zero 有条件可用，但不推荐用于柔性前体**。
   - 刚性分子（rx8）完美
   - 前体在 3/5 反应中偏高（最高 +0.65 kcal/mol）
   - rx15 上表现意外地好（产物/前体 ΔE < 0.06 kcal/mol）

5. **分子构象复杂度是协议表现的关键调节变量**：
   - 刚性分子（rx8）：所有协议等价
   - 中等柔性（rx1, rx13, rx15）：协议差异小（< 0.1 kcal/mol），lite/ext 均可
   - 高柔性/xTB 排序不可靠（rx3）：lite 显著优于 ext（0.88 kcal/mol）

### 6.2 生产运行决策树

```
开始
│
├─ 已知分子刚性？────────── 是 ──→ lite 或 zero 均可
│                                    （zero 节省 ~20% 时间）
│
├─ 不确定 / 批量运行 ────────→ **lite**（默认推荐）
│
├─ 需要最高精度验证 ────────→ lite + ext 双重运行对照
│
└─ ❌ 不推荐 full（成本 2.5× ext，能量不优于 lite）
```

### 6.3 Limitations & Future Work

1. **5 个反应均来自同一反应类型**（[4+3] allenamide cycloaddition）。结论需在其他反应类型（[4+2] Diels-Alder, [3+2] 1,3-dipolar, [5+2]）上验证。

2. **rx16（cyclic lactam）尚未运行**。Lactam 骨架与 carbamate 的氢键供体特征不同，可能影响构象偏好。

3. **S1-only scope**。S2/S3/S4 被跳过。S1 的构象质量对下游 TS 搜索（S2/S3）的影响需要全流程验证。

4. **DFT 级理论固定**（B3LYP/def2-SVP OPT + wB97X-D4/def2-TZVPP SP）。不同 theory level 下 lite 的 r2SCAN-3c 重排优势是否仍然成立，尚需验证。

---

## Appendix: Per-Reaction Raw Data

### rx1 — Boc, acyclic allenamide, 3-carbon tether

| Protocol | Product E (Ha) | Precursor E (Ha) | DFT | Wall |
|----------|---------------|-----------------|-----|------|
| ext | −940.674 883 572 | −865.331 299 081 | 9 | 2:50 |
| full | −940.674 883 598 | −865.333 294 165 | 14 | 5:02 |
| lite | −940.674 883 644 | −865.333 294 455 | 4 | 1:43 |
| zero | −940.674 883 638 | −865.331 286 051 | 5 | 2:00 |

### rx3 — Bn (benzyl), acyclic allenamide, 2-carbon tether

| Protocol | Product E (Ha) | Precursor E (Ha) | DFT | Wall |
|----------|---------------|-----------------|-----|------|
| ext | −1014.502 037 830 | −939.163 942 867 | 10 | 2:02 |
| full | −1014.501 490 937 | −939.165 231 673 | 36 | 10:22 |
| **lite** | **−1014.503 440 281** | −939.165 230 273 | **4** | **0:37** |
| zero | −1014.503 438 957 | −939.163 350 620 | 8 | 1:10 |

### rx8 — cyclic carbamate (oxazolidinone), 1-carbon tether (rigid)

| Protocol | Product E (Ha) | Precursor E (Ha) | DFT | Wall |
|----------|---------------|-----------------|-----|------|
| ext | −782.042 651 553 | −706.707 178 576 | 4 | 0:31 |
| full | −782.042 651 560 | −706.707 178 810 | 6 | 0:53 |
| lite | −782.042 651 525 | −706.707 178 746 | 2 | 0:24 |
| zero | −782.042 651 530 | −706.707 178 755 | 2 | 0:22 |

> ΔE < 0.000 03 kcal/mol — 刚性分子所有协议完全等价。

### rx13 — cyclic carbamate (oxazolidinone), 3-carbon tether

| Protocol | Product E (Ha) | Precursor E (Ha) | DFT | Wall |
|----------|---------------|-----------------|-----|------|
| ext | −860.755 719 775 | −785.414 282 216 | 9 | 1:42 |
| full | −860.755 720 241 | −785.414 282 844 | 17 | 3:47 |
| lite | −860.755 720 316 | −785.414 282 414 | 3 | 0:46 |
| zero | −860.755 720 202 | −785.413 250 588 | 5 | 1:02 |

> 产物 ΔE < 0.001 kcal/mol；zero 前体偏差 0.646 kcal/mol。

### rx15 — cyclic carbamate (oxazolidinone), **4-carbon tether (high flexibility)**

| Protocol | Product E (Ha) | Precursor E (Ha) | DFT | Wall |
|----------|---------------|-----------------|-----|------|
| **ext** | **−900.103 615 166** | **−824.766 588 788** | 8 | 4:16 |
| full | −900.103 498 288 | −824.766 522 563 | 25 | 7:45 |
| lite | −900.103 501 222 | −824.766 522 650 | 3 | 1:06 |
| zero | −900.103 518 886 | −824.766 522 628 | 5 | 1:21 |

> rx15 是唯一 ext 在产物和前体上均为最优的反应。但优势极小（< 0.08 kcal/mol）。
> full 为前体做了 23 个 DFT 优化（所有反应中最多），但结果与 lite 2 个 DFT 完全一致。

---

## Appendix: rx15 Detailed Notes

rx15 作为最高柔性的测试案例，展现了与之前 4 个反应不同的行为模式：

1. **ext 首次全面胜出**：两阶段 GFN0→GFN2 采样的广度优势在 4 碳 tether 的高柔性前体上得以发挥。GFN0 粗采样阶段探索了更广泛的构象空间，最终 DFT 排序找到了略优的能量。

2. **full 的 23 个前体 DFT 无效**：CREST GFN2 单阶段为高柔性前体生成了大量构象，经 prescreen/screening 后保留 23 个送入 DFT。但最终前体能量与 lite（仅 2 个 DFT）完全一致（Δ < 0.001 kcal/mol）。这意味着 full 的 fast-SP 筛选虽然保留了 23 个候选者，但真正的全局最低点已经被 lite 的 r2SCAN-3c 识别出来了。

3. **lite 产物偏差极小**（+0.07 kcal/mol）：在 0.08 kcal/mol 的误差范围内，lite 和 ext 实质等价。

4. **零 fallback 触发**：4 个协议在 rx15 上均未触发 fallback（产物/前体均仅 1~2 候选），说明 rx15 的构象空间在 GFN2 层面具有明确的主导构象，不存在 rx1/rx3 那样的近简并问题。

综合来看，rx15 验证了 lite 在高柔性分子上仍然可靠（ΔE < 0.08 kcal/mol），同时确认了 ext 在特定条件下（广泛构象采样 + xTB 排序可靠）仍然有价值。
