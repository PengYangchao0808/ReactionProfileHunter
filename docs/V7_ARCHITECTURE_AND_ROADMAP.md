# V7 完整数据集架构与后续设计路线图

**版本**: V7.1 | **日期**: 2026-05-26 | **状态**: Phase 1 完成，Phase 2-5 设计中

---

## 1. 总体架构概览

V7 的核心设计原则：

- **分层构建**：单次运行特征 (`rph_features`) → 条件层合并 (`condition_thermo`/`condition_feature_merger`) → Stage1 分支-产率数据集 (`build_ml_dataset`)
- **特征族隔离**：每个特征族有独立的命名空间前缀（`s1_*`, `s2_*`, `int_*`, `ts_*`, `ts_freq.*`, `ts_opt.*`, `int_mw_*`），避免跨族命名冲突
- **退化优先**：所有提取器 `get_required_inputs() → []`，缺失文件 → NaN + warning code，不抛异常
- **QC 隔离**：特征提取阶段 extract-only（`job_run_policy: disallow`），不运行新 QC 计算

### 1.1 三层数据架构

```
Layer 1: 单次运行特征 (rph_features extract)
  ┌─────────────────────────────────────────────┐
  │ 每个 (reaction, condition, branch) 三元组   │
  │ 产生一份 features_raw.csv                   │
  │ ~250 列 (固定 + 计算 + QA)                  │
  └──────────────┬──────────────────────────────┘
                 │
                 ▼
Layer 2: 条件层合并 (condition_feature_merger)
  ┌─────────────────────────────────────────────┐
  │ 同一条件的两个 branch 合并                   │
  │ 产生 merged_features.csv                    │
  │ + 条件-level 热力学 Gibbs 自由能             │
  │ + major_*/minor_*/delta_* 前缀化            │
  └──────────────┬──────────────────────────────┘
                 │
                 ▼
Layer 3: Stage1 分支-产率数据集 (build_ml_dataset)
  ┌─────────────────────────────────────────────┐
  │ 跨 reaction × condition 聚合                │
  │ 合并实验数据 (产率/DR)                       │
  │ 输出: stage1_branch_yield_dataset.csv       │
  │ ~217 列 → ML 训练                           │
  └─────────────────────────────────────────────┘
```

---

## 2. 特征族完整地图

### 2.1 全量特征族目录 (12 族)

| # | 特征族 | 来源阶段 | 提取器 | 列数 | 完整度现状 |
|---|--------|---------|--------|:---:|-----------|
| 1 | Condition | 实验文献 | `reaxys_cleaned.csv` | 14 | ✅ 100% |
| 2 | S1_conformational_energy | S1 DFT | `conformational_population` | 32 | ⚠ 系综太小 |
| 3 | S1_conformational_population | S1 DFT | `conformational_population` | 24 | ⚠ 系综太小 |
| 4 | S1_preorganization | S1 + product | `conformational_preorganization` / `precursor_preorganization` | 26 | ⚠ 命名漂移 |
| 5 | S1_activation | S1 DFT + Shermo | `step1_activation` | ⚠ | ⚠ 缺失 |
| 6 | S1_tether_flexibility | S1 conf ensemble | `torsion_flexibility` | ⚠ | ⚠ 缺失 |
| 7 | S1_topology | product SMILES | `geometry` | 3 | ✅ |
| 8 | S2_electronic | S2/S3 fchk | `fmo_cdft_dipolar` + `step2_cyclization` | 26 | ✅ |
| 9 | **S3_intermediate_geometry** | **S3 opt** | **`intermediate_features`** | **9** | **🆕 V7.1** |
| 10 | **S3_intermediate_electronic** | **S3 fchk** | **`intermediate_features`** | **6** | **🆕 V7.1** |
| 11 | **S3_intermediate_wavefunction** | **S3 fchk + Multiwfn** | **`intermediate_multiwfn`** | **12** | **🆕 V7.1** |
| 12 | **S3_ts_quality** | **S3 opt + freq** | **`ts_quality`** + **`ts_optimization_metadata`** | **13** | **🆕 V7.1** |

> 🆕 = V7.1 新增，当前实现完成

### 2.2 S3 特征族列明细

#### S3_intermediate_geometry (9 列)

| 列 | 含义 | 退化策略 |
|----|------|---------|
| `int_geom.rg` | 回转半径 (Å) | 无 XYZ → NaN |
| `int_geom.min_nonbonded` | 最小非键接触 (Å) | 同上 |
| `int_geom.close_contacts` | close contact 计数 | 同上 |
| `int_geom.close_contacts_density` | close_contacts / natoms | 同上 |
| `int_geom.forming_distance_1` | forming bond 1 距离 (Å) | 无 forming_bonds → NaN |
| `int_geom.forming_distance_2` | forming bond 2 距离 (Å) | 同上 |
| `int_geom.asynch` | |d1 − d2| (Å) | 同上 |
| `int_geom.rmsd_to_s2_guess` | S3 vs S2 guess RMSD (Å) | 缺任一方 → NaN |
| `int_geom.source` | `s3_intermediate_xyz` / `s2_intermediate_fallback` / `missing` | 总是存在 |

#### S3_intermediate_electronic (6 列)

| 列 | 含义 | 来源 |
|----|------|------|
| `int_cdft.eps_homo` | HOMO (eV) | `intermediate_fchk` via `read_fchk_cdft_indices()` |
| `int_cdft.eps_lumo` | LUMO (eV) | 同上 |
| `int_cdft.mu` | 化学势 (eV) | 同上 |
| `int_cdft.eta` | 化学硬度 (eV) | 同上 |
| `int_cdft.omega` | 亲电指数 (eV) | 同上 |
| `int_cdft.gap` | HOMO-LUMO gap (eV) | 同上 |

#### S3_intermediate_wavefunction (12 列)

| 列 | 含义 | 来源 |
|----|------|------|
| `int_mw_fukui_fplus_atomA/B` | forming atom 亲核 Fukui | Multiwfn (per-atom) |
| `int_mw_fukui_fminus_atomA/B` | forming atom 亲电 Fukui | 同上 |
| `int_mw_fukui_f0_atomA/B` | forming atom 自由基 Fukui | 同上 |
| `int_mw_dual_descriptor_atomA/B` | forming atom dual descriptor | 同上 |
| `int_mw_rho_bcp_forming1` | BCP 电子密度 | QTAIM (bond A-B) |
| `int_mw_laplacian_bcp_forming1` | BCP Laplacian | 同上 |
| `int_mw_status` | `ok` / `degraded` / `disabled` / `failed` | QA |
| `int_mw_missing_reason` | 退化原因 | QA |
| `int_mw_warnings_count` | 警告数 | QA |

> BCP 分析针对 forming bond (atomA, atomB)。Intermediate 中 forming contacts 不一定形成真正的 BCP，此时 `degraded` 而非 `failed`。

#### S3_ts_quality (13 列)

| 列 | 含义 | 来源 |
|----|------|------|
| `ts_freq.n_imag` | 虚频数 | TS freq log |
| `ts_freq.imag_cm1` | 最负虚频 (cm⁻¹，负值) | 同上 |
| `ts_freq.imag_abs_cm1` | 虚频绝对值 (cm⁻¹) | 同上 |
| `ts_freq.mode_overlap_forming` | 虚频模式与 forming bond 方向重叠 [0,1] | freq log + TS XYZ |
| `ts_freq.has_single_imag` | 恰好 1 个虚频 → 1 | QA |
| `ts_opt.method_used` | Berny / QST2 / Intermediate_Reuse | s3_resume.json |
| `ts_opt.rescue_attempted` | 是否触发 rescue | 同上 |
| `ts_opt.rescue_route` | qst2 / oscillation / none | 同上 |
| `ts_opt.failure_kind` | oscillation / convergence / none | 同上 |
| `ts_opt.n_cycles` | 优化迭代次数 | 同上 |
| `ts_opt.converged` | 最终收敛标志 | .rph_step_status.json |
| `ts_opt.has_l2_sp` | L2 SP 产物存在 | 文件检查 |
| `ts_opt.has_freq` | freq 计算产物存在 | 文件检查 |

> `ts_freq.*` 为解释性特征，不替代现有 `ts.*`（`ts.n_imag`, `ts.imag1_cm1_abs`, `ts.dipole_debye`）。两者并存。

### 2.3 Schema.py 中的 V7 Named Feature Sets

```python
FEATURE_SETS_V7 = {
    "F0_condition":      ["temperature_K", "LA_intensity", ...],         # 4 col
    "F1_core":           ["major_s1_E_span", "major_s2_eta", ...],       # 12 col
    "F2_yield":          ["major_s1_dG_act", "major_s2_gedt_value", ...],# 15 col
    "F3_selectivity":    ["delta_s1_d_reactive_q90", "delta_s2_eta", ...],# 22 col
    "F4_pca_candidates": [],                                              # dynamic
}
```

---

## 3. 提取器注册表

### 3.1 当前活动提取器 (20 个)

| # | 提取器 | 插件名 | 输出前缀 | 状态 |
|---|--------|--------|---------|:---:|
| 1 | `FmoCdftDipolarParser` | `fmo_cdft_dipolar` | `s2_*` | active |
| 2 | `TSQualityExtractor` | `ts_quality` | `ts.*`, `ts_freq.*` | active |
| 3 | `QCChecksExtractor` | `qc_checks` | `qc.*` | active |
| 4 | `GeometryExtractor` | `geometry` | `geom.*` | active |
| 5 | `ThermoExtractor` | `thermo` | `thermo.*` | active |
| 6 | `NBOE2Extractor` | `nbo_e2` | `nbo.*` | active |
| 7 | `NICSExtractor` | `nics` | `nics.*` | active |
| 8 | `InteractionAnalysisExtractor` | `interaction_analysis` | `ia_*` | active |
| 9 | `MultiwfnFeaturesExtractor` | `multiwfn_features` | `mw_*` | active |
| 10 | `ConformationalPopulationExtractor` | `conformational_population` | `s1_Hpop`, ... | active |
| 11 | `ConformationalPreorganizationExtractor` | `conformational_preorganization` | `s1_d_reactive_*` | active |
| 12 | `PrecursorPreorganizationExtractor` | `precursor_preorganization` | `s1_preorg_*` | active |
| 13 | `TorsionFlexibilityExtractor` | `torsion_flexibility` | `s1_tether_*` | active |
| 14 | `Step1ActivationExtractor` | `step1_activation` | `s1_dG_act`, ... | active |
| 15 | `Step2CyclizationExtractor` | `step2_cyclization` | `s2_*` | active |
| 16 | `TSOptimizationMetadataExtractor` | `ts_optimization_metadata` | `ts_opt.*` | 🆕 active |
| 17 | `IntermediateFeaturesExtractor` | `intermediate_features` | `int_geom.*`, `int_energy.*`, `int_cdft.*`, `int_charge.*` | 🆕 active |
| 18 | `IntermediateMultiwfnExtractor` | `intermediate_multiwfn` | `int_mw_*` | 🆕 active |
| 19 | `ASMEnrichmentExtractor` | `asm_enrichment` | `asm_*` | active |
| 20 | `PrecursorGeometryExtractor` | `precursor_geometry` | `s1_precursor_*` | active |

> 🆕 = V7.1 新增

### 3.2 提取器与特征族映射

```
提取器                         → 特征族
─────────────────────────────────────────────────────
conformational_population      → S1_conformational_energy + S1_conformational_population
conformational_preorganization → S1_preorganization
precursor_preorganization      → S1_preorganization (precursor 侧)
torsion_flexibility            → S1_tether_flexibility
step1_activation               → S1_activation
geometry                       → S1_topology
fmo_cdft_dipolar               → S2_electronic (dipolar + GEDT)
step2_cyclization              → S2_electronic (CDFT)
multiwfn_features              → S2_electronic (波函数补充)
intermediate_features           → S3_intermediate_geometry + S3_intermediate_electronic 🆕
intermediate_multiwfn           → S3_intermediate_wavefunction 🆕
ts_quality                      → TS_baseline + S3_ts_quality 🆕
ts_optimization_metadata        → S3_ts_quality 🆕
thermo                          → thermo.* (热力学基线)
qc_checks                       → qc.* (QA)
```

---

## 4. 数据流全链路

```
SMILES ──→ S1 (anchor/conformer) ──→ S2 (retro scan) ──→ S3 (TS opt) ──→ S4 (features)
                                                                              │
              ┌───────────────────────────────────────────────────────────────┘
              ▼
     rph_features extract
     ├── 构造 FeatureContext (20+ optional path handles)
     ├── 运行 20 个提取器 (parallel-safe, extract-only)
     ├── 输出 features_raw.csv (~250 列)
     ├── 输出 features_mlr.csv (ML-ready 子集)
     └── 输出 feature_meta.json (provenance, warnings, status)

              │  × (reaction × condition × branch)
              ▼
     condition_feature_merger
     ├── 合并 major/minor branch 特征
     ├── 计算 delta_* = major − minor
     ├── 注入条件热力学 Gibbs (Shermo at actual T)
     └── 输出 merged_features.csv

              │  × (reaction × condition)
              ▼
     rph_features build-dataset (或 rph_ml build-dataset)
     ├── 合并实验数据 (产率, DR)
     ├── 生成 support_count (每个 reaction 的条件数)
     ├── 标记 structured_trainable / failure_aware_trainable
     └── 输出 stage1_branch_yield_dataset.csv

              │
              ▼
     rph_ml train-yield / select-features
     ├── LORO CV (Leave-One-Reaction-Out)
     ├── Two-Head Ridge: yield_head + selectivity_head
     ├── 特征选择: F0→F1→F2→F3 增量比较
     └── 输出 structured_predictions.csv + summary.json
```

---

## 5. ML 训练架构

### 5.1 双头模型

```
Input: X ∈ R^(N×D)  (N samples, D features)
       ├── condition 特征 (temperature, LA, ...)
       ├── major_* 特征 (主分支计算特征)
       ├── minor_* 特征 (副分支计算特征)
       └── delta_* 特征 (major − minor)

Model: Ridge(alpha) with two heads

  Head 1 — yield_head:
    Target: yield_fraction ∈ [0, 1]  (scale-transformed)
    Metric: R², MAE

  Head 2 — selectivity_head:
    Target: p_major ∈ [0, 1]  (major branch probability)
    Metric: R², MAE

CV: Leave-One-Reaction-Out (LORO)
    ─ 确保同一 reaction 的多个 condition 不同时出现在 train/val
```

### 5.2 特征集演进路径

```
F0 (condition-only)       → 4 特征, baseline
F1 (core low-collinearity) → 12 特征, 已验证
F2 (yield-biased)          → 15 特征, yield head 优化
F3 (selectivity-biased)    → 22 特征, DR head 优化
F4 (PCA candidates)        → 动态, 降维探索
F5 (interaction)           → 条件×计算交互项
```

### 5.3 后续 F6 候选 (S3 intermediate 特征加入 ML)

```python
# F6_yield_with_intermediate (16 特征)
FEATURE_SET_F6_YIELD_WITH_INTERMEDIATE = FEATURE_SET_F2_YIELD + [
    "major_int_geom.forming_distance_1",   # 中间体 forming bond 距离 → 预反应接近程度
    "major_int_geom.asynch",               # 异步性 → 反应模式判别
    "major_int_energy.rel_to_ts_kcal",     # 中间体→TS 能垒 (中间体侧视角)
    "delta_int_cdft.eta",                  # branch 间中间体硬度差异
]

# F7_selectivity_with_intermediate (24 特征)
FEATURE_SET_F7_SELECTIVITY_WITH_INTERMEDIATE = FEATURE_SET_F3_SELECTIVITY + [
    "major_int_geom.rmsd_to_s2_guess",     # S3 vs S2 几何变化量 → branch 稳定性
    "major_int_cdft.omega",                # 中间体亲电性
    "delta_int_cdft.omega",                # branch 间亲电性差异
    "major_int_mw_dual_descriptor_atomA",  # forming atom 波函数特征
    "delta_int_mw_dual_descriptor_atomA",  # branch 间 forming atom 差异
    "major_ts_freq.mode_overlap_forming",  # TS 模式与反应坐标对齐度
]
```

---

## 6. 后续设计路线图

### Phase 2: 条件层热力学富集 (预计 1 周)

**目标**: 把温度依赖的 Gibbs 自由能从条件层注入 Stage1 数据集，替代当前恒温 298K 假设。

```
当前状态: thermo.dG_activation 基于 298K Shermo
          → 低温反应 (-78°C) 的 Gibbs 被高估

修复方案:
  1. condition_thermo.py 已实现 per-condition Shermo 计算
  2. backfill_condition_thermo.py 已实现批量回填
  3. condition_feature_merger 注入条件温度下的 G_ts, G_int, G_product
  4. Stage1 数据集使用 "thermo.dG_activation_at_T" 替代固定 298K 值

关键产出:
  - condition_features_mlr.csv 新增列:
    * thermo.G_ts_at_T       (kcal/mol, at reaction T)
    * thermo.G_int_at_T      (kcal/mol)
    * thermo.G_product_at_T  (kcal/mol)
    * thermo.dG_activation_at_T  (= G_ts − G_reactants)
    * thermo.dG_reaction_at_T    (= G_product − G_reactants)
```

### Phase 3: NBO/NICS 高级电子特征 (低优先级, 仅在需求明确时启动)

**当前状态**: NBO/NICS extractor 已注册但默认 disabled。Multiwfn extractor 仅做 Fukui/Dual/QTAIM。

```
候选扩展:
  1. NBO E2 供体-受体相互作用能
     - 需要 S3 中 enable_nbo=true
     - 可产出: nbo_e2_max (最强供体→受体稳定化能), nbo_e2_sum
  2. NICS(0) / NICS(1) 芳香性指标
     - 需要 Multiwfn NICS 模块
     - 可产出: nics_0_forming_ring, nics_1_forming_ring
  3. 电荷分解分析 (CDA)
     - fragment → fragment 电荷转移分解

决策门: 先验证 S3 intermediate + TS freq 特征对 yield/DR 的增量预测力。
        如果 F6/F7 已显著提升, NBO/NICS 优先级降低。
```

### Phase 4: 多反应类型泛化 (中期, 需要数据积累)

**目标**: RPH 当前只支持 [4+3] 环加成。扩展到 [4+2], [3+2], [5+2] 需要：

```
1. 反应类型特征化
   - reaction_type: [4+3] / [4+2] / [3+2] / [5+2]
   - forming_bond_count: 2 (concerted) / 1 (stepwise)
   - 环大小: product_ring_size

2. 跨反应类型训练
   - 问题: 不同反应类型的产率/DR 分布不同
   - 方案: multi-task learning (shared encoder + type-specific heads)
     ＊ 或: reaction_type 作为 categorical feature + interaction terms

3. 数据集需求
   - 每个反应类型 ≥ 20 个 reaction (当前 [4+3] 有 10 个)
   - 需要扩大实验数据采集

4. 化学空间多样性
   - 当前所有反应共享相同骨架 (oxallyl cation 前体)
   - 需要多样化前体 (不同取代基, 不同离去基团)
```

### Phase 5: 端到端主动学习管线 (远期愿景)

```
SMILES 库 ──→ RPH 批量运行 ──→ 特征提取 ──→ ML 预测 ──→ 不确定性估计
                                                              │
                                              ┌───────────────┘
                                              ▼
                                    高不确定性样本 ──→ 实验验证 ──→ 反馈到训练集
```

核心组件：
1. **批量调度器**: 按优先级和资源约束调度 RPH 运行
2. **不确定性量化**: ensemble variance / conformal prediction
3. **采集函数**: expected improvement, UCB, Thompson sampling
4. **实验反馈闭环**: 实验室产率/DR 自动回写到训练集

---

## 7. 文件与模块地图

| 模块 | 文件 | 职责 |
|------|------|------|
| **特征提取** | `rph_features/rph_features/feature_miner.py` | extract-only 管线入口 |
| | `rph_features/rph_features/context.py` | FeatureContext 数据容器 |
| | `rph_features/rph_features/extractors/*.py` | 20 个提取器 |
| | `rph_features/rph_features/schema.py` | 列定义 / 特征族 / V7 命名集 |
| **条件层** | `rph_core/steps/condition_thermo.py` | 条件温度下 Shermo 热力学 |
| | `rph_core/steps/condition_feature_merger.py` | branch 合并 + delta 计算 |
| **上游数据** | `rph_core/steps/step3_opt/ts_optimizer.py` | S3 中间体数据产出 |
| | `rph_core/steps/runners.py` | step→step 数据传递 |
| | `rph_core/orchestrator.py` | 主流程编排 |
| **波函数** | `rph_core/utils/multiwfn_runner.py` | Multiwfn 封装 (Fukui/Dual/QTAIM) |
| **ML 训练** | `rph_ml/rph_ml/train_yield.py` | 双头 Ridge + LORO CV |
| | `rph_ml/rph_ml/feature_selection.py` | 增量特征选择 |
| | `rph_ml/rph_ml/diagnostics.py` | 特征质量诊断 |
| **数据构建** | `rph_features/rph_features/dataset_builder.py` | Stage1 数据集构建 |

---

## 8. 验收标准

### Phase 1 (当前 — 已完成 ✅)

- [x] `features_raw.csv` 包含 `int_geom.*`, `int_energy.*`, `int_cdft.*`, `int_charge.*`, `int_mw_*`, `ts_freq.*`, `ts_opt.*` 列
- [x] 所有新提取器在文件缺失时输出 NaN 而非崩溃
- [x] `feature_meta.json` 包含 `artifact_presence.s3_intermediate_xyz`
- [x] LSP error = 0
- [x] 快速 CI gate 通过 (import + no-QC)
- [x] 9 个专项 degrade 测试通过

### Phase 2 (待实施)

- [ ] `condition_features_mlr.csv` 包含 `thermo.dG_activation_at_T`
- [ ] 低温条件下的 Gibbs 值与 298K 基准值可比较
- [ ] Stage1 数据集中产率与 `thermo.dG_activation_at_T` 的 Spearman ρ > 与 `thermo.dG_activation` 的 ρ

### Phase 3-5 (未来)

- [ ] F6/F7 特征集在 LORO CV 上的增量 R² > 0.05 vs F2/F3
- [ ] 跨反应类型模型 RMSE < 0.15 (yield) 和 < 0.20 (DR)
- [ ] 主动学习管线每轮实验后模型 R² 单调递增
