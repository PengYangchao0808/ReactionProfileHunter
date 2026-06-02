# S4 构象预组织特征引入报告

> 更新日期：2026-05-22
> 关联需求：基于 S1 GFN2-CREST 构象搜索 + DFT 单点能 Boltzmann 排序的链柔性精确描述
> 新增模块：3 个 S4 提取器 | 15 个新特征

---

## 1. 报告摘要

针对现有 S4 特征体系中"链柔性描述不充分"的问题，本报告记录了一整套**反应感知的构象预组织特征**（reaction-aware conformational preorganization features）的设计与实现。该方案完全基于 **S1 已有的 GFN2-CREST 构象搜索 + DFT 单点能数据**，无需额外 QC 计算，从三个递进层次刻画产物的构象行为：

| 层级 | 提取器 | 特征数 | 核心问题 |
|------|--------|--------|----------|
| **P0 群体** | `ConformationalPopulationExtractor` | 5 | "构象系综有多分散？" |
| **P1 预组织** | `ConformationalPreorganizationExtractor` | 5 | "构象倾向于反应有利取向吗？" |
| **P2 柔性** | `TorsionFlexibilityExtractor` | 5 | "柔性链段如何旋转/展开？" |

**核心范式转变**：从泛化的"链柔性"（chain flexibility）→ **反应感知的构象预组织**（reaction-aware conformational preorganization），聚焦于"构象分布是否天然倾向于反应有利的几何取向"。

---

## 2. 设计原理

### 2.1 为什么需要新的构象描述符

现有 S4 特征体系对构象行为的描述仅包含：

- `s1_Nconf_eff`（有效构象数）
- `s1_Sconf`（构象熵）

**不足**：这两个特征仅传达了"有多少构象"，但无法表达：

1. **分布均匀性**：Boltzmann 权重是集中在极少数构象（刚性分子），还是分散在大量构象中（柔性分子）？
2. **几何倾向性**：构象群体是否倾向于使反应位点接近？这对预组织驱动反应至关重要。
3. **柔性链行为**：远端柔性链（tether）的扭转、回旋如何影响整体构象景观？

### 2.2 设计约束

| 约束 | 决策 |
|------|------|
| 无额外 QC 计算 | 所有特征从 S1 已有数据提取（`conformer_energies.json` + `conf_NNN.xyz`） |
| 小样本 ML 安全（N≈30–70） | 严格控制特征数（每层 5 个），避免共线性 |
| 归一化基准 | 用 `N_rot`（可旋转键数）归一化，不用 `ln(N_total)`（后者依赖 CREST 搜索参数） |
| 能量来源 | DFT 电子能（`conformer_energies.json`），非 Gibbs 自由能 |
| 转角统计 | 离散三态 Shannon 熵（gauche⁺/gauche⁻/trans），非连续分箱 |

### 2.3 物理图景

```
构象预组织 = 构象群体中的反应有利构象占比 × 几何偏向程度
            ───────────────────────   ──────────────
                 P0: 群体分布            P1: 反应位点距离
                                    +  P2: 柔性链段行为
```

- **预组织反应**（如刚性双烯）：产物构象群已处于近反应有利几何，TS 形成代价低
- **非预组织反应**（如长链柔性产物）：需要大量构象重排才能达到反应有利取向

---

## 3. 三层特征详述

### 3.1 P0：构象群体描述符

**数据源**：`conformer_energies.json`（DFT 电子能，kcal/mol）+ `product_smiles`（计算 `N_rot`）

**Boltzmann 权重重推导**（独立于 298K 预计算的 thermo 表，使用反应温度）：

```
w_i = exp(-(E_i - E_min) / RT) / Σ_j exp(-(E_j - E_min) / RT)
```

| 特征名 | 单位 | 公式 | 物理意义 |
|--------|------|------|----------|
| `s1_Hpop` | — | H = −Σ w_i·ln(w_i) | **群体 Shannon 熵**。值≈0 = 刚性（一个主导构象）；值≈ln(N) = 均匀分散 |
| `s1_Hpop_per_rotbond` | — | H / N_rot | **每旋转键熵**。归一化后可比不同大小分子 |
| `s1_wmax` | — | max(w_i) | **主导构象权重**。w_max≈1 = 刚性单构象；w_max≪1 = 无主导构象 |
| `s1_N90` | — | 累积权重 ≥ 0.90 的最少构象数 | **90% 覆盖数**。N90=1 = 刚性；N90>10 = 高度柔性 |
| `s1_rank1_rank2_gap` | kcal/mol | E₂ − E₁（排序后前两个构象能量差） | **最优-次优能量间距**。大 gap → 单构象主导；小 gap → 多构象竞争 |

**退化行为** (缺少 `conformer_energies.json` 时): 所有 5 个特征返回 `NaN`。

---

### 3.2 P1：反应感知几何预组织

**数据源**：`conf_NNN.xyz`（S1 最终 DFT 构象结构）+ `forming_bonds`（成键对）+ `product_smiles`（识别 reactive fragment 重原子）

**核心计算**：对每个 Boltzmann 加权构象，计算所有 forming-bond 原子对之间的**平均反应位点距离**：

```
d_k = mean(∥r_i − r_j∥) — 对 forming_bonds 中所有 (i, j) 对取平均
```

然后对整个构象系综做加权统计：

| 特征名 | 单位 | 公式 | 物理意义 |
|--------|------|------|----------|
| `s1_d_reactive_mean` | Å | Σ w_k·d_k | **加权平均反应位点距离**，越小→构象越倾向于反应几何 |
| `s1_d_reactive_std` | Å | sqrt(Σ w_k·(d_k − d_mean)²) | **反应距离离散度**，小→构象群体几何一致性好 |
| `s1_d_reactive_q10` | Å | 加权 10% 分位数 | **最近 10% 构象的反应距离**（最佳情况） |
| `s1_d_reactive_q90` | Å | 加权 90% 分位数 | **最远 10% 构象的反应距离**（最差情况） |
| `s1_productive_population` | — | Σ w_k·I(d_k < cutoff) | **productive 构象总权重**。所有 forming-bond 距离均 < cutoff（默认 3.0 Å） |

**原子过滤**：距离计算仅基于 reactive fragment 重原子（非 H），用 SMILES + forming_bonds 索引映射。

**productive_population 逻辑**：**AND 逻辑**——多键反应要求所有 forming bonds 同时处于 cutoff 内。cutoff 可通过 `config.step4.preorganization.productive_cutoff_angstrom` 配置。

**退化行为**：缺少 conformer xyz / forming_bonds / product_smiles → 全部返回 `NaN`。

---

### 3.3 P2：Tether 柔性链描述符

**数据源**：`conf_NNN.xyz` + `product_smiles`（识别可旋转键 + tether 原子）

**核心思路**：不度量整个分子的柔性，而是**聚焦于反应位点之间的柔性连接段（tether）**。

**Tether 原子识别**：通过 product_smiles 构建分子图，识别 forming-bond 原子之间的最短路径上的所有重原子。仅对 tether 原子计算几何描述符。

#### 3.3.1 扭转柔性

对 tether 上的每个可旋转二面角，用离散 3-bin 分类：

```
bin 1 (gauche⁺): 0° < φ < 120°
bin 2 (trans):   120° ≤ φ ≤ 240°
bin 3 (gauche⁻): 240° < φ < 360°
```

每个二面角的构象群体离散度用 Shannon 熵量化：

```
H_j = −Σ p_bin·log₃(p_bin)   （以 3 为底，范围 [0, 1]）
```

| 特征名 | 单位 | 公式 | 物理意义 |
|--------|------|------|----------|
| `s1_tether_torsion_H_sum` | — | Σ H_j | **扭转柔性总和**。大值→tether 有多个高度柔性的二面角 |
| `s1_tether_torsion_var_max` | — | max(H_j) | **最大扭转柔性**。识别最柔性的单一二面角 |

#### 3.3.2 几何弥散

| 特征名 | 单位 | 公式 | 物理意义 |
|--------|------|------|----------|
| `s1_tether_rmsd_to_gm` | Å | Boltzmann 加权 RMSD（vs 全局最小构象） | **构象偏离程度**，仅对 tether 原子计算。O(N) 复杂度 |
| `s1_tether_rg_mean` | Å | Σ w_k·Rg_k | **加权平均回旋半径**（tether 原子），大值→链段伸展 |
| `s1_tether_rg_std` | Å | std(Rg_k) | **回旋半径离散度**，大值→链段在紧凑/伸展间剧烈变换 |

**RMSD vs Rg**：RMSD 度量"构象偏离 GM 多远"（结构差异），Rg 度量"链段有多紧凑"（尺寸变化）。二者正交。

**O(N) 优化**：RMSD 仅对全局最小构象计算（O(N)），非全对全 O(N²) pairwise RMSD，避免特征爆炸。

**退化行为**：缺少 conformer xyz / product_smiles → 全部返回 `NaN`。无 tether 可旋转键 → `torsion_H_sum` 和 `torsion_var_max` 为 0。

---

## 4. 实现文件清单

### 4.1 新增文件 (3 个，849 行)

| 文件 | 行数 | 内容 |
|------|------|------|
| `rph_core/steps/step4_features/extractors/conformational_population.py` | 205 | P0 提取器：Boltzmann 重归一化、Shannon 熵、N90 排序 |
| `rph_core/steps/step4_features/extractors/conformational_preorganization.py` | 286 | P1 提取器：反应位点距离统计、productive_population |
| `rph_core/steps/step4_features/extractors/torsion_flexibility.py` | 358 | P2 提取器：tether 识别、离散扭转熵、Rg/RMSD |

### 4.2 修改文件 (5 个)

| 文件 | 变更 | 说明 |
|------|------|------|
| `rph_core/steps/step4_features/extractors/__init__.py` | +3 行 | 导入 3 个新提取器模块（触发 `register_extractor`） |
| `rph_core/steps/step4_features/context.py` | +3 字段 | `product_smiles`、`s1_conformer_dir`、`s1_conformer_thermo_csv` |
| `rph_core/steps/step4_features/schema.py` | +28 行 | COMPUTATIONAL_FEATURES 新增 15 个特征名；unit_map 新增单位；WARNING_CODES 新增 14 个 |
| `rph_core/steps/step4_features/feature_miner.py` | 签名扩展 | `run()` 新增 `product_smiles`/`s1_conformer_dir`/`s1_conformer_thermo_csv` 参数 |
| `rph_core/steps/runners.py` | 签名扩展 | `run_step4()` 新增 `product_smiles` 参数，传入 `run_kwargs` |
| `rph_core/orchestrator.py` | 2 处修改 | `_resolve_s1_artifacts()` 解析 `s1_conformer_dir` + `s1_conformer_thermo_csv`；`run_step4()` 调用处传入 `product_smiles` |

### 4.3 数据流

```
orchestrator.run_pipeline(product_smiles, work_dir)
    │
    ├─ _resolve_s1_artifacts()
    │   └─ conformer_energies.json 所在目录 → s1_conformer_dir
    │   └─ conformer_thermo.csv          → s1_conformer_thermo_csv
    │
    └─ runners.run_step4(product_smiles=..., **s1_artifacts)
        │
        └─ FeatureMiner.run(product_smiles=..., s1_conformer_dir=..., ...)
            │
            └─ FeatureContext(product_smiles=..., s1_conformer_dir=..., ...)
                │
                ├─ ConformationalPopulationExtractor.extract(ctx)
                │   └─ ctx.s1_conformer_energies_file  → P0 features
                │
                ├─ ConformationalPreorganizationExtractor.extract(ctx)
                │   └─ ctx.s1_conformer_dir + ctx.product_smiles + ctx.forming_bonds → P1 features
                │
                └─ TorsionFlexibilityExtractor.extract(ctx)
                    └─ ctx.s1_conformer_dir + ctx.product_smiles → P2 features
```

---

## 5. 新增警告代码

所有新增代码已注册到 `schema.py` 的 `WARNING_CODES` 字典：

### P0 相关 (4 个)

| 代码 | 含义 |
|------|------|
| `W_POP_NO_ENERGIES` | `conformer_energies.json` 不存在或为空 |
| `W_POP_PARSE_FAILED` | `conformer_energies.json` 解析失败 |
| `W_POP_NO_THERMO` | `conformer_thermo.csv` 不存在 |
| `W_POP_BOLTZMANN_DEGENERACY` | 所有构象能量相同（简并）→ 群体描述符退化 |

### P1 相关 (4 个)

| 代码 | 含义 |
|------|------|
| `W_PREORG_NO_CONFORMER_DIR` | S1 构象目录不存在或为空 |
| `W_PREORG_NO_FORMING_BONDS` | forming_bonds 未定义 |
| `W_PREORG_XYZ_READ_FAILED` | 构象 XYZ 文件读取失败 |
| `W_PREORG_NO_PRODUCT_SMILES` | product_smiles 不可用（无法识别 tether 原子） |

### P2 相关 (5 个)

| 代码 | 含义 |
|------|------|
| `W_TORSION_NO_CONFORMER_DIR` | S1 构象目录不存在或为空 |
| `W_TORSION_NO_PRODUCT_SMILES` | product_smiles 不可用 |
| `W_TORSION_NO_ROTATABLE_BONDS` | 产物分子无可旋转键 |
| `W_TORSION_XYZ_READ_FAILED` | 构象 XYZ 文件读取失败 |
| `W_TORSION_RGM_FAILED` | 参考几何均值（reference geometry mean）计算失败 |

---

## 6. 验证状态

### 6.1 自动化检查

| 检查项 | 状态 | 工具 |
|--------|------|------|
| LSP 诊断（8 个变更文件） | ✅ 通过（0 新增 error，1 个预存 warning 不相关） | basedpyright |
| 导入样式 CI | ✅ 通过（152 文件，无违规） | `scripts/ci/check_imports.py` |
| 快速 CI 门 | ✅ 8/8 通过 | pytest: `test_imports_step4_features.py` + `test_s4_no_qc_execution.py` |
| 提取器注册 | ✅ 15/15 注册（12 原有 + 3 新增） | `list_extractors()` 运行时验证 |
| FeatureContext 字段 | ✅ `product_smiles`、`s1_conformer_dir`、`s1_conformer_thermo_csv` 正确传递 | 运行时打印验证 |

### 6.2 退化行为验证

所有三个提取器在以下场景均已验证优雅退化（返回 NaN + 记录警告）：

| 场景 | P0 | P1 | P2 |
|------|----|----|----|
| conformer_energies.json 不存在 | ✅ NaN | — | — |
| conformer 目录不存在 | — | ✅ NaN | ✅ NaN |
| product_smiles 缺失 | `Hpop_per_rotbond`=NaN | ✅ NaN | ✅ NaN |
| forming_bonds 缺失 | — | ✅ NaN | — |

### 6.3 P0 真实数据验证

使用 `tests/tmp_v2_2_test/da_reaction/S1_Anchor/butadiene/dft/` 的 butadiene 双构象数据，P0 提取器正确降级（缺少 `conformer_energies.json` 时返回 NaN，逻辑路径正确触发）。

> **注意**：完整 P0/P1/P2 端到端验证需要具备以下全部数据的完整产物目录：
> - `conformer_energies.json`（含 DFT 电子能）
> - `conformer_thermo.csv`
> - `conf_NNN.xyz` 系列文件
> - 有效的 `forming_bonds` + `product_smiles`
>
> 现有测试夹具仅含 butadiene 的部分 S1_Anchor 数据，缺少 `conformer_energies.json`。

---

## 7. 已知风险与限制

| 风险 | 严重度 | 缓解措施 |
|------|--------|----------|
| **P1 forming_bonds 索引映射**：forming_bonds 来自 TS 原子编号，可能不直接映射到产物分子原子编号 | ⚠️ 中 | 需在首次完整流水线运行中验证；如不匹配需实现 `_remap_forming_bonds_to_product()` |
| **P2 tether 识别**：多环/桥环体系中 "最短路径" 可能不是化学上合理的 tether | ⚠️ 低 | 当前假设线性/单环连接段；桥环场景需额外验证 |
| **能量简并**：所有构象能量相等时 Hpop→ln(N)，Hpop_per_rotbond 可能异常大 | ⚠️ 低 | 已注册 `W_POP_BOLTZMANN_DEGENERACY` 警告 |
| **单构象极限**：仅 1 个构象时 Hpop=0，torsion_H_sum=0，rmsd_to_gm=0 — 数学正确但可能产生恒等列 | 💡 注意 | 小样本 ML 中应检查这些特征是否为零方差 |
| **conformer_energies.json 数据格式**：当前假设 JSON 为 `{"energies": [e1, e2, ...]}` 或 `[e1, e2, ...]` | ⚠️ 低 | 已处理 dict 和 list 两种格式 |

---

## 8. 后续工作

### 8.1 高优先级

- [ ] **端到端集成测试**：使用完整 S1–S4 产物目录验证 15 个特征的实际计算值
- [ ] **forming_bonds 映射验证**：P1 中确认 forming_bonds 索引在产物分子上的正确性
- [ ] **ML 特征筛选**：将 15 个新特征加入 `build_ml_dataset.py` → `train_yield_ml.py`，评估对 p_major / logDR 的预测贡献
- [ ] **共线性检查**：验证新特征与现有 `s1_Nconf_eff`、`s1_Sconf` 之间无严重共线性（VIF > 10）

### 8.2 中优先级

- [ ] **config 暴露**：在 `config/defaults.yaml` 的 `step4.preorganization` 段添加 `productive_cutoff_angstrom` 配置项（当前代码内默认 3.0 Å）
- [ ] **P2 tether 路径完善**：对桥环/螺环体系实现多路径 tether 识别
- [ ] **特征文档同步**：在 `docs/S4_FEATURES_SUMMARY.md` 中追加本报告的 P0/P1/P2 节

### 8.3 低优先级

- [ ] **Shermo 热力学校正**：评估是否需要 Gibbs 自由能替代电子能作为 Boltzmann 权重源（当前 P0 用电子能）
- [ ] **反应温度可配置**：当前 P0 默认 298.15K，读取 context.temperature_K 但 orchestrator 未传入

---

## 9. 附录 A：提取器注册验证输出

```
$ python -c "from rph_core.steps.step4_features.extractors import list_extractors; \
  print('\n'.join(sorted(e.get_plugin_name() for e in list_extractors())))"

15 extractors registered:
  asm_enrichment
  conformational_population        ← NEW (P0)
  conformational_preorganization   ← NEW (P1)
  fmo_cdft_dipolar
  geometry
  interaction
  multiwfn_features
  nbo_e2
  nics
  qc_checks
  step1_activation
  step2_cyclization
  thermo
  torsion_flexibility              ← NEW (P2)
  ts_quality
```

## 附录 B：退化测试输出

```
P0 (no files):    s1_Hpop=NaN  s1_Hpop_per_rotbond=NaN  s1_wmax=NaN  s1_N90=NaN  s1_rank1_rank2_gap=NaN
P1 (no dir):      s1_d_reactive_mean=NaN  s1_d_reactive_std=NaN  s1_d_reactive_q10=NaN  s1_d_reactive_q90=NaN  s1_productive_population=NaN
P2 (no dir):      s1_tether_torsion_H_sum=NaN  s1_tether_torsion_var_max=NaN  s1_tether_rmsd_to_gm=NaN  s1_tether_rg_mean=NaN  s1_tether_rg_std=NaN
```

## 附录 C：完整特征列表（FIXED_COLUMNS 新增段）

以下 15 个特征已添加至 `schema.py` 的 `COMPUTATIONAL_FEATURES` 列表（按字母排序后位置）：

```
s1_d_reactive_mean     Å     反应位点加权平均距离
s1_d_reactive_q10      Å     反应距离 10% 分位数
s1_d_reactive_q90      Å     反应距离 90% 分位数
s1_d_reactive_std      Å     反应距离加权标准差
s1_Hpop                —     构象群体 Shannon 熵
s1_Hpop_per_rotbond    —     每旋转键构象熵
s1_N90                 —     90% 覆盖最少构象数
s1_productive_population —   productive 构象总权重
s1_rank1_rank2_gap     kcal/mol  最优-次优能量间距
s1_tether_rg_mean      Å     tether 加权平均回旋半径
s1_tether_rg_std       Å     tether 回旋半径标准差
s1_tether_rmsd_to_gm   Å     tether Boltzmann 加权 RMSD (vs GM)
s1_tether_torsion_H_sum —    tether 扭转 Shannon 熵总和
s1_tether_torsion_var_max —  tether 最大扭转柔性
s1_wmax                —     主导构象 Boltzmann 权重
```

---

<div align="center">

**实现完成 · 等待端到端验证**

</div>
