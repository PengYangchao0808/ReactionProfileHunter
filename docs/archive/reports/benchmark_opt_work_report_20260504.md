# Benchmark OPT 工作情况报告（修正版）

**报告日期**：2026-05-04  
**修正说明**：上一版报告依据较早的 `docs/weekly_report_20260503.md` 判断 OPT 仍处于 `7/20` 进行中状态；经直接核对仓库落地目录 `Output/benchmark_dft_theory/experiments/session_opt_full_orca_sp1`，当前 Phase 2 GEO/OPT 主矩阵已经完整完成。本文以该会话目录中的 `manifest`、`geo_result.json`、`phase2_global_summary`、`GEO_BENCHMARK_FINAL_REPORT` 和 `final_report` 为准。  
**核心会话**：`Output/benchmark_dft_theory/experiments/session_opt_full_orca_sp1`  
**Benchmark 类型**：固定 `SP-1` reader 的几何优化方法比较（`fixed_sp1_opt_comparison`）

---

## 1. 总体结论

Benchmark OPT 的 **Phase 2 GEO/OPT 主矩阵已经完成**：4 个 `[4+3]` 反应（`rx1/rx3/rx8/rx15`）× 5 个几何优化方法（`GEO-1`~`GEO-5`）共 **20/20 个任务全部完成，0 failed**。

所有 TS 优化均收敛到有效的一阶鞍点，`GEO_BENCHMARK_FINAL_REPORT.md` 明确记录成功率为 **100%**。当前已经不再是“计算未完成/err125 阻塞”的状态；旧报告中关于 `rx3/rx8` ORCA `err125`、`rx15` 未运行、`7/20 完成` 的描述已过期，应全部替换为“主矩阵已落地完成”。

当前结论可以分为两个层次：

1. **落地状态**：Phase 1 SP、Shermo、Phase 2 GEO、final report 均已完成；
2. **方法推荐**：详细 GEO 报告推荐 `GEO-2` 作为 routine TS optimization 首选，`GEO-3` 用于 high-throughput screening，`GEO-5` 用于更高精度；而 `final_report.md/json` 的 legacy 投票聚合给出 `recommended_l1_geo_method = GEO-1`。这两个推荐口径需要在后续报告生成逻辑中统一。

---

## 2. 证据与落地文件

### 2.1 全局 manifest 状态

`manifests/benchmark_manifest.json` 显示：

| 阶段 | 状态 |
|---|---|
| migrate | done |
| baseline | done |
| phase1_sp | done |
| shermo_correction | done |
| phase2_geo | done |
| final_report | done |

会话更新时间为 `2026-05-04T12:07:11+00:00`，覆盖反应为 `1, 3, 8, 15`。

### 2.2 主矩阵结果文件

会话目录下已落地 **20 个** `geo_result.json`：

```text
rx1/geo/GEO-1..5/geo_result.json
rx3/geo/GEO-1..5/geo_result.json
rx8/geo/GEO-1..5/geo_result.json
rx15/geo/GEO-1..5/geo_result.json
```

全局报告文件已经生成：

- `reports/phase2_global_summary.md/json`
- `reports/GEO_BENCHMARK_FINAL_REPORT.md`
- `reports/final_report.md/json`
- `reports/geo_sp_refresh_audit.md/json`
- `reports/geo_sp_refresh_summary.json`

### 2.3 Phase 2 全局完成矩阵

`reports/phase2_global_summary.md` 记录如下：

| rx_id | SP reader | methods completed |
|---|---|---:|
| rx1 | SP-1 | 5 / 5 |
| rx3 | SP-1 | 5 / 5 |
| rx8 | SP-1 | 5 / 5 |
| rx15 | SP-1 | 5 / 5 |

因此当前主线判断应为：**OPT 主矩阵已完成，不存在未补齐的 rx/GEO 组合。**

---

## 3. 方法矩阵与运行条件

### 3.1 固定 SP reader

OPT 阶段统一使用：

- `SP-1 = wB97M-V/def2-TZVPP`
- ORCA / RIJCOSX
- CPCM acetone
- 10 cores

这一点由 `GEO_BENCHMARK_FINAL_REPORT.md` 和各 `geo_summary.md` 共同确认。

### 3.2 GEO 方法矩阵

| 方法 | 理论级别 | 定位 |
|---|---|---|
| GEO-1 | B3LYP-D3BJ/def2-SVP | industry/baseline standard |
| GEO-2 | PBE0-D3BJ/def2-SVP | routine all-rounder |
| GEO-3 | r2SCAN-3c/mTZVPP | high-throughput / low-cost screening |
| GEO-4 | M06-2X/def2-SVP | Minnesota high-HF-exchange method |
| GEO-5 | wB97X-D3BJ/def2-SVP | higher-accuracy, slower option |

---

## 4. 完整 OPT 结果概览

### 4.1 每个反应的完成状态

| 反应 | GEO-1 | GEO-2 | GEO-3 | GEO-4 | GEO-5 | 状态 |
|---|---:|---:|---:|---:|---:|---|
| rx1 | ✅ | ✅ | ✅ | ✅ | ✅ | 5/5 completed |
| rx3 | ✅ | ✅ | ✅ | ✅ | ✅ | 5/5 completed |
| rx8 | ✅ | ✅ | ✅ | ✅ | ✅ | 5/5 completed |
| rx15 | ✅ | ✅ | ✅ | ✅ | ✅ | 5/5 completed |

合计：**20/20 completed, 0 failed**。

### 4.2 完整结果表（来自 `GEO_BENCHMARK_FINAL_REPORT.md`）

| rx | 方法 | dE‡ | dG‡ | dE_rxn | dG_rxn | TS ω₁ | RMSD | Runtime |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rx1 | GEO-1 | 4.28 | 4.67 | -59.2 | -54.6 | -82 | 0.1846 | 2.2h |
| rx1 | GEO-2 | 4.28 | 4.43 | -59.2 | -54.5 | -88 | 0.1846 | 2.1h |
| rx1 | GEO-3 | 4.28 | 6.15 | -59.2 | -53.0 | -110 | 0.1846 | 1.7h |
| rx1 | GEO-4 | 4.28 | 5.59 | -59.2 | -53.6 | -226 | 0.1846 | 1.7h |
| rx1 | GEO-5 | 4.28 | 3.39 | -59.2 | -54.5 | -13 | 0.1846 | 3.2h |
| rx3 | GEO-1 | 5.29 | 5.95 | -51.2 | -46.5 | -12 | 0.6526 | 1.5h |
| rx3 | GEO-2 | 5.29 | 6.37 | -51.2 | -46.9 | -205 | 0.6526 | 1.6h |
| rx3 | GEO-3 | 5.29 | 6.82 | -51.2 | -45.9 | -178 | 0.6526 | 1.6h |
| rx3 | GEO-4 | 5.29 | 6.23 | -51.2 | -47.1 | -255 | 0.6526 | 2.0h |
| rx3 | GEO-5 | 5.29 | 5.30 | -51.2 | -46.4 | -20 | 0.6526 | 2.6h |
| rx8 | GEO-1 | 2.66 | 4.38 | -49.9 | -45.8 | -253 | 0.0922 | 0.6h |
| rx8 | GEO-2 | 2.66 | 4.07 | -49.9 | -45.8 | -191 | 0.0922 | 0.6h |
| rx8 | GEO-3 | 2.66 | 3.74 | -49.9 | -45.9 | -143 | 0.0922 | 0.6h |
| rx8 | GEO-4 | 2.66 | 4.04 | -49.9 | -46.0 | -248 | 0.0922 | 0.7h |
| rx8 | GEO-5 | 2.66 | 4.37 | -49.9 | -45.5 | -291 | 0.0922 | 0.7h |
| rx15 | GEO-1 | 2.68 | 3.97 | -54.8 | -50.0 | -257 | 0.1144 | 1.1h |
| rx15 | GEO-2 | 2.68 | 3.91 | -54.8 | -49.8 | -208 | 0.1144 | 1.1h |
| rx15 | GEO-3 | 2.68 | 3.87 | -54.8 | -49.9 | -195 | 0.1144 | 1.4h |
| rx15 | GEO-4 | 2.68 | 4.38 | -54.8 | -50.0 | -286 | 0.1144 | 1.3h |
| rx15 | GEO-5 | 2.68 | 4.73 | -54.8 | -49.0 | -328 | 0.1144 | 1.3h |

---

## 5. 关键分析结论

### 5.1 成功率与 TS 合法性

所有 20 个 TS 优化均成功，且均为有效一阶鞍点。`GEO_BENCHMARK_FINAL_REPORT.md` 明确记录：

- `20/20 completed, 0 failed`
- all TS optimizations converged to valid first-order saddle points
- 总 wall time 约 **29.4 h**

这说明此前 `err125` 属于中间过程问题，当前落地结果中已不再构成阻塞。

### 5.2 dE‡ 在同一反应内几乎方法不敏感

详细报告指出，固定 `SP-1` reader 后，同一反应内各 GEO 方法的 dE‡ 基本一致：

- rx1：4.28 kcal/mol
- rx3：5.29 kcal/mol
- rx8：2.66 kcal/mol
- rx15：2.68 kcal/mol

因此 OPT 方法差异主要体现在：

- 频率/热力学校正导致的 dG‡ 差异；
- TS 虚频大小；
- runtime；
- 几何偏移与稳定性。

### 5.3 dG‡ 排名

按平均 dG‡（越低越好）：

| 排名 | 方法 | Mean dG‡ |
|---:|---|---:|
| 1 | GEO-5 | 4.45 |
| 2 | GEO-2 | 4.69 |
| 3 | GEO-1 | 4.74 |
| 4 | GEO-4 | 5.06 |
| 5 | GEO-3 | 5.14 |

`GEO-5` 给出最低平均 dG‡，尤其在 `rx1/rx3` 上最低，但总 runtime 最高。

### 5.4 Runtime 排名

按总 runtime（越低越好）：

| 排名 | 方法 | Total runtime |
|---:|---|---:|
| 1 | GEO-2 | 5.2h |
| 2 | GEO-3 | 5.3h |
| 3 | GEO-1 | 5.4h |
| 4 | GEO-4 | 5.7h |
| 5 | GEO-5 | 7.8h |

`GEO-2` 是整体最快方案，`GEO-5` 约慢 50%。

### 5.5 几何一致性

各方法在同一反应内收敛到非常接近的 TS 几何；报告按反应给出的 baseline-aligned RMSD 为：

| Reaction | RMSD | 评价 |
|---|---:|---|
| rx1 | 0.1846 Å | good agreement |
| rx3 | 0.6526 Å | significant deviation |
| rx8 | 0.0922 Å | excellent agreement |
| rx15 | 0.1144 Å | excellent agreement |

其中 `rx3` 的 TS 相对 baseline 偏移最大，说明它仍是几何敏感性最高的案例。

### 5.6 SP refresh 审计已完成

`geo_sp_refresh_audit.md` 显示：

- refreshable: 20
- failed_matches: 0

这说明 20 个 GEO 结果均可匹配到刷新后的 SP 计算，且匹配审计没有失败。审计表中每个 rx/GEO/point 都有独立 matched hash，说明“指标相近/部分 dE 相同”不应简单解读为计算缺失。

需要注意的是：`GEO_BENCHMARK_FINAL_REPORT.md` 和 `geo_sp_refresh_summary.json` 是两个不同层级的产物。前者是当前人读的最终详细报告；后者记录了 refresh 后的 SP 能量与 hash。后续如果要以 refresh 后能量作为最终口径，应重新聚合一次 detailed GEO report，避免报告口径不一致。

---

## 6. 方法推荐口径

### 6.1 详细 GEO 报告的推荐

`GEO_BENCHMARK_FINAL_REPORT.md` 给出的实际使用建议是：

| 场景 | 推荐方法 | 理由 |
|---|---|---|
| Routine TS optimization | GEO-2 | 最快、dG‡ 稳定、SCF 表现好 |
| High-throughput screening | GEO-3 | 3c composite，适合快速筛选 |
| High accuracy | GEO-5 | rx1/rx3 最低 dG‡，但更慢 |
| Industry standard | GEO-1 | 最广泛验证、行为可预期 |

因此从实际工作流角度，当前最合理的主推荐是：

> **常规生产使用 GEO-2；大规模筛选用 GEO-3；需要更高精度时用 GEO-5；保守对照保留 GEO-1。**

### 6.2 `final_report` 的推荐

`reports/final_report.md/json` 记录：

- `recommended_l2_sp_method = SP-REF`
- `recommended_l1_geo_method = GEO-1`
- Phase2 votes: `GEO-1 = 2`, `GEO-4 = 1`, `GEO-2 = 1`

这里存在一个需要清理的口径问题：`phase2_global_summary.json` 明确写了 `ranking_policy = comparison_only_no_winner`，但 `final_report` 仍然通过 legacy vote 聚合生成了 `recommended_l1_geo_method = GEO-1`。因此不应只看 `final_report` 的 `GEO-1`，而应把它视为旧聚合逻辑的结果；对实际方法选择，应优先采用详细 GEO 报告中的场景化推荐。

---

## 7. 与上一版报告相比的修正点

| 上一版说法 | 修正后说法 |
|---|---|
| OPT 仅 7/20 完成 | OPT 主矩阵 20/20 完成 |
| rx3 只有 GEO-1 成功 | rx3 GEO-1~5 全部完成 |
| rx8 全部 err125 阻塞 | rx8 GEO-1~5 全部完成 |
| rx15 尚未运行 | rx15 GEO-1~5 全部完成 |
| rx1 结果一致性是 P0 阻塞 | 当前已有完整详细报告与 SP refresh audit；指标相近需解释，但不是“计算未完成” |
| 下一步主要是补齐计算 | 下一步主要是统一报告推荐口径、决定采用 GEO-2/GEO-3/GEO-5 的生产策略 |

---

## 8. 后续建议

### P0：统一 final report 与 detailed GEO report 的推荐逻辑

当前 `final_report` 给出 `GEO-1`，而详细 GEO 报告推荐 routine 使用 `GEO-2`。建议修改最终聚合逻辑：在 `ranking_policy = comparison_only_no_winner` 时，不应输出单一 `recommended_l1_geo_method`，或应输出场景化推荐字段。

### P1：决定生产默认几何方法

基于现有结果，建议默认策略为：

- `GEO-2`：常规 TS 优化默认；
- `GEO-3`：高通量初筛；
- `GEO-5`：重点样本或复核；
- `GEO-1`：baseline 对照和保守兼容。

### P1：重新聚合 refresh 后 SP 结果（如采用 refresh 口径）

`geo_sp_refresh_summary.json` 已落地且审计通过。如果后续要使用 refreshed L2 SP 能量作为最终判据，应生成对应的 refreshed GEO final report，避免 `geo_result.json`、refresh summary、detailed report 三者出现解释负担。

### P2：明确 Phase 2b 是否纳入当前 benchmark 范围

当前会话未发现 `p2b_result.json`。这不影响 20/20 GEO 主矩阵完成结论；但如果要把 precursor/intermediate 也按每个 GEO 方法补成完整 5-point profile，需要另行执行并报告 `p2b_complement.py` 的结果。

---

## 9. 当前一句话状态

`session_opt_full_orca_sp1` 中 Benchmark OPT 主矩阵已经完整落地：**4 reactions × 5 GEO methods = 20/20 completed, 0 failed**；常规使用建议倾向 `GEO-2`，高通量倾向 `GEO-3`，高精度复核倾向 `GEO-5`，后续重点不再是补计算，而是统一 final report 的推荐逻辑与是否采用 refreshed SP 结果作为最终口径。
