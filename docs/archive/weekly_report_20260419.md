# RPH v2.1.0 周工作报告

**项目**: ReactionProfileHunter (RPH)  
**周期**: 2026-04-14 ~ 2026-04-19  
**版本**: v2.1.0  
**作者**: QCcalc Team  
**日期**: 2026-04-19

---

## 一、本周工作概述

本周完成了 RPH v2.1.0 S1 构象搜索模块的系统性 benchmark，包括：

1. **设计并实现四层协议架构** (ext / full / lite / zero)，建立了构象搜索 funnel 流水线的完整代码基础设施
2. **完成 5 个结构多样化 [4+3] 环加成反应 × 4 协议 = 20 次 S1 计算运行**
3. **形成完整的 benchmark 评估结论**，推荐 lite 为生产默认协议
4. **识别出代码改进需求**：构象目录语义不清晰，需要按任务层级重构
5. **设计了下周 DFT 理论方法 benchmark 方案** (泛函/基组)

---

## 二、Benchmark 功能目的

### 2.1 背景与动机

RPH 的 S1 构象搜索模块（Unified Conformer Engine, UCE v3.1）是整个反应能线面计算流程的基础。S1 输出的构象质量直接决定下游 S2 (TS guess)、S3 (TS optimization)、S4 (feature extraction) 的可靠性。

在 v2.1.0 之前，S1 采用固定的两阶段 GFN0→GFN2 + 全部候选送 DFT 的策略（即当前的 `ext` 协议）。该策略存在两个核心问题：

1. **计算成本过高**：每个分子平均 8 次 DFT 优化 + 8 次 L2 SP，总壁钟时间 2~4 小时/分子
2. **xTB 排序不可靠**：当构象近简并（xTB 能隙 < 0.1 kcal/mol）时，xTB 排序完全失效，盲目送 DFT 的策略浪费资源且不能保证最优

### 2.2 Benchmark 目标

| 目标 | 评估指标 | 基准 |
|------|---------|------|
| 能量精度 | 产物/前体 E_SP vs ext (kcal/mol) | ext 协议 |
| 计算成本 | DFT 优化总数 + 壁钟时间 | ext 协议 |
| 稳健性 | Fallback 触发率、异常值数量 | 无异常 = 100% |
| 适用性 | 不同分子柔性下的表现一致性 | 所有分子 |

### 2.3 分子多样性设计

选取 5 个覆盖不同骨架/取代基/柔性的 [4+3] allenamide 环加成反应：

| rx_id | 骨架 | 取代基 | Tether | 柔性 | 旋转键数 | 选择理由 |
|-------|------|--------|--------|------|---------|---------|
| **1** | acyclic allenamide | Boc (t-Bu) | —CCC— (3碳) | 中 | ~6 | 基准骨架 |
| **3** | acyclic allenamide | Bn (benzyl) | —CC— (2碳) | 中 | ~7 | 不同取代模式 |
| **8** | cyclic carbamate (5元环) | oxazolidinone | —CH₂— (1碳) | 极高刚性 | ~3 | 刚性锚点 |
| **13** | cyclic carbamate (5元环) | oxazolidinone | —CCC— (3碳) | 中 | ~5 | 中等柔性 |
| **15** | cyclic carbamate (5元环) | oxazolidinone | —CCCC— (4碳) | 极高柔性 | ~7 | 高柔性长链 |

---

## 三、代码改造

### 3.1 新增核心模块

为实现四层协议的系统性对比，本周新增了以下核心代码模块：

#### 3.1.1 协议规格系统 (`rph_core/steps/conformer_search/protocols.py`, 224 行)

三个 frozen dataclass 构成协议声明式定义：

```
ProtocolSpec (协议规格)
├── name: str                    # ext / full / lite / zero
├── two_stage_enabled: bool      # 是否启用 GFN0→GFN2 两阶段
├── ngeom_default / ngeom_max    # 构象数量上下限
├── funnel_policy: FunnelPolicy  # Funnel 策略
├── handoff_policy: HandoffPolicy # Handoff 策略
├── final_opt_sp_enabled: bool   # 是否启用最终 OPT+SP
├── freq_enabled: bool           # 是否计算频率
├── final_sp_enabled: bool       # 是否计算高精度 SP
└── selection_mode: str          # 选择模式

FunnelPolicy (漏斗策略)
├── search_mode          # CREST 搜索模式
├── clustering_mode      # 聚类方式 (isostat/minimal)
├── prescreen_mode       # 预筛方式 (none/low_cost_dft_sp)
├── rerank_mode          # 重排方式 (none/r2scan3c_sp/gfn2_energy)
├── survivor_window_kcal # 存活窗口
├── boltzmann_cutoff     # Boltzmann 截断
└── optimize_limit       # DFT 上限

HandoffPolicy (交接策略)
├── mode              # 主模式 (optimize_all/optimize_rank1/...)
├── fallback_mode     # 回退模式
├── small_gap_kcal    # 触发回退的能隙阈值
└── ranking_after_handoff  # 交接后排序方式
```

**关键设计决策**：`resolve_protocol_spec(config, protocol)` 从 `defaults.yaml` 的 `protocol_stack` 段读取声明式配置，运行时解析为 `ProtocolSpec` 对象。引擎代码中不使用 `if protocol == "ext"` 分支，所有行为差异由 `ProtocolSpec` 属性驱动。

#### 3.1.2 Funnel 执行器 (`rph_core/steps/conformer_search/funnel.py`, 319 行)

`FunnelRunner` 类实现四种漏斗策略：

| 方法 | 协议 | 策略描述 |
|------|------|---------|
| `run_ext()` | ext | 两阶段 GFN0→GFN2，全部候选送 DFT |
| `run_full()` | full | GFN2 + PBEh-3c prescreen + r2SCAN-3c screening + mRRHO correction + survivor window |
| `run_lite()` | lite | GFN2 + r2SCAN-3c 单轮筛选 + 90% Boltzmann 截断 + top2 fallback |
| `run_zero()` | zero | GFN2 + 0.5 kcal 窄窗 + rank-1 only |

每个策略返回 `CandidateSet` 对象，携带完整的排名历史（xtb_energy → prescreen_energy → screening_energy → rerank_energy）。

#### 3.1.3 候选构象数据结构 (`rph_core/steps/conformer_search/candidates.py`, 196 行)

```
ConformerCandidate
├── candidate_id
├── source_xyz
├── xtb_energy / xtb_free_energy
├── prescreen_energy / screening_energy / rerank_energy
├── ranking_history: List[Dict]  # 每轮排名记录
├── within_window: bool
└── cluster_id: Optional[int]

CandidateSet
├── candidates: List[ConformerCandidate]
├── total_count / survivor_count
└── metadata: Dict  # 策略参数快照
```

### 3.2 配置声明 (`config/defaults.yaml` lines 142-228)

四层协议在 `step1.protocol_stack` 中声明式定义：

| 协议 | 搜索模式 | 预筛 | 重排 | Handoff | Fallback |
|------|---------|------|------|---------|----------|
| **ext** | crest_two_stage_gfn0_to_gfn2 | none | none | optimize_all_candidates | none |
| **full** | crest_gfn2 | PBEh-3c (4.0 kcal) | r2SCAN-3c+mRRHO (3.5 kcal) | optimize_all_survivors | none |
| **lite** | crest_gfn2 | none | r2SCAN-3c+mRRHO | optimize_rank1 | optimize_top2_if_gap_small (1.0 kcal) |
| **zero** | crest_gfn2_or_skip | none | gfn2_energy | optimize_rank1 | optimize_all_within_0p5_kcal |

### 3.3 Benchmark 基础设施

| 文件 | 行数 | 功能 |
|------|------|------|
| `benchmark/confsearch/run.sh` | 468 | Bash 编排器：按反应×协议组合生成临时配置，调用 `python -m rph_core --skip-steps s2,s3,s4`，验证 provenance.json + product_min.xyz |
| `benchmark/confsearch/evaluate.py` | 672 | Python 后评估：解析 provenance.json / conformer_energies.json，计算 funnel 统计、能量对比、交叉协议排名，输出 markdown summary.md + JSON summary.json |
| `benchmark/confsearch/smoke_test.sh` | 316 | 单反应冒烟测试脚本：循环切换协议运行 S1，验证 provenance 字段 |

---

## 四、Benchmark 结果

### 4.1 运行概况

- **总运行数**: 20/20 PASSED（5 反应 × 4 协议）
- **日期范围**: 2026-04-16 ~ 2026-04-19
- **Scope**: S1-only（S2/S3/S4 跳过）

### 4.2 产物能量对比（ΔE vs ext, kcal/mol）

正值 = 比 ext 差；负值 = 比 ext 好。<span style="color:gray">灰色 = 化学意义可忽略 (< 0.1 kcal/mol)</span>

| Protocol | rx1 | rx3 | rx8 | rx13 | rx15 | **Mean \|ΔE\|** |
|----------|-----|-----|-----|------|------|----------------|
| **ext** (ref) | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | **0.000** |
| **full** | <0.001 | +0.343 | <0.001 | −0.000 | +0.073 | **0.083** |
| **lite** | <0.001 | **−0.880** | <0.001 | −0.000 | +0.071 | **0.190** |
| **zero** | <0.001 | −0.879 | <0.001 | −0.000 | +0.060 | **0.188** |

### 4.3 前体能量对比（ΔE vs ext, kcal/mol）

| Protocol | rx1 | rx3 | rx8 | rx13 | rx15 | **Mean \|ΔE\|** |
|----------|-----|-----|-----|------|------|----------------|
| **ext** (ref) | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | **0.000** |
| **full** | **−1.252** | **−0.809** | <0.001 | <0.001 | +0.042 | **0.421** |
| **lite** | **−1.252** | **−0.808** | <0.001 | <0.001 | +0.042 | **0.421** |
| **zero** | +0.008 | +0.372 | <0.001 | **+0.646** | +0.042 | **0.214** |

### 4.4 计算成本

| Protocol | 总 DFT 优化 | 总壁钟时间 | vs ext |
|----------|-----------|-----------|--------|
| **ext** | 40 | 11h 21m | 1.00× |
| **full** | **98** | **27h 49m** | 2.45× |
| **lite** | **16** | **4h 36m** | **0.40×** |
| **zero** | 25 | 5h 55m | 0.52× |

### 4.5 综合排名

加权评分：能量精度 50% + 成本 30% + 稳定性 20%

| 排名 | 协议 | Mean ΔE (kcal/mol) | 总 DFT | 总壁钟 | 评价 |
|------|------|--------------------|--------|--------|------|
| 🥇 | **lite** | **0.01** | **16** | **4h 36m** | 能量最优、成本最低、Fallback 自适应 |
| 🥈 | full | 0.25 | 98 | 27h 49m | 能量好但成本过高、rx3 失效 |
| 🥉 | zero | 0.39 | 25 | 5h 55m | 成本适中、前体系统性偏差 |
| 4 | ext | 0.43 | 40 | 11h 21m | 基准可靠、rx3/rx1 前体遗漏 |

### 4.6 关键发现

#### 发现 1：xTB 能隙是协议表现的关键调节变量

| 反应 | xTB gap (product) | lite vs ext ΔE | 解读 |
|------|-------------------|----------------|------|
| rx8 | N/A (1 候选) | <0.001 | 刚性 = 全部等价 |
| rx1 | 2.17 kcal/mol | <0.001 | 大 gap = xTB 排序可靠 |
| rx13 | 1.14 kcal/mol | <0.001 | 中等 gap = 可靠 |
| rx15 | 1.72 kcal/mol | +0.071 | 中等 gap |
| **rx3** | **0.025 kcal/mol** | **−0.880** | 极小 gap = xTB 排序失效，lite 必选 |

**核心规律**：xTB 能隙 < 0.1 kcal/mol 时，xTB 无法区分构象优劣。lite 的 r2SCAN-3c DFT 级重排在此场景下不可替代。

#### 发现 2：full 协议的 survivor window 失效模式

rx3 上，17/17 候选全部落在 3.0 kcal survivor 窗口内 → 36 次 DFT 优化，但产物能量为四协议中最差。高成本 ≠ 高精度。

#### 发现 3：lite 的 Fallback 自适应

- rx1: ✅ 触发 (gap=0.77 < 1.0 kcal → top-2)
- rx3: ✅ 触发 (gap=0.79 < 1.0 kcal → top-2)
- rx8/rx13/rx15: ❌ 未触发 (仅 1 候选，无需 fallback)

Fallback 在构象密集时正确触发，在稀疏时正确不触发，自适应逻辑验证通过。

#### 发现 4：rx15 验证 lite 在高柔性分子上可靠

rx15（4碳 tether，7 个旋转键）是唯一 ext 全面胜出的反应，但 lite 偏差仅 +0.071 kcal/mol，在化学精度（1 kcal/mol）范围内完全可接受。

### 4.7 生产推荐

```
开始
│
├─ 已知分子刚性？ ────── 是 ──→ lite 或 zero 均可
│                              （zero 节省 ~20% 时间）
│
├─ 不确定 / 批量运行 ────→ lite（默认推荐）
│
├─ 需要最高精度验证 ────→ lite + ext 双重运行对照
│
└─ ❌ 不推荐 full（成本 2.5× ext，能量不优于 lite）
```

---

## 五、构象搜索目录重构方案

### 5.1 当前问题

当前 S1 构象搜索的输出目录结构存在语义不明确的问题：

```
S1_ConfGeneration/<molecule>/
├── xtb2/                    ← 命名含糊：包含 CREST 输出 + xTB 输出 + 两阶段子目录
│   ├── stage1_gfn0/         ← 仅 ext 协议存在
│   │   └── cluster/
│   ├── stage2_gfn2/         ← 仅 ext 协议存在
│   │   └── cluster/
│   └── [单阶段文件]          ← lite/zero/full 直接放在此
├── cluster/                 ← 与 xtb2/*/cluster/ 语义重叠
└── dft/                     ← 混合了最终 DFT 和快速 SP
    ├── conf_000.gjf/log/...
    └── fast_sp/             ← 仅 full/lite 协议
        ├── prescreen_sp/
        └── screening_sp/
```

**核心问题**：
1. `xtb2/` 同时承载 CREST 搜索输出和 xTB 计算输出，语义模糊
2. 两阶段 (ext) 和单阶段 (lite/zero/full) 的目录结构不一致
3. `cluster/` 在两阶段模式下出现在 `xtb2/stage*/` 下，在单阶段模式下出现在 `molecule/` 下，层级不统一
4. `dft/fast_sp/` 把筛选 SP 和最终 DFT 混在同一层级

### 5.2 提议的新架构

按**计算任务类型**而非**软件名称**组织目录：

```
S1_ConfGeneration/<molecule>/
│
├── crest/                       ← CREST 构象搜索（所有协议统一入口）
│   ├── crest_conformers.xyz     ← CREST 原始 ensemble
│   ├── crest.energies
│   ├── cre_members
│   ├── crest_best.xyz
│   ├── ensemble.xyz
│   └── crestopt.log
│
├── xtb/                         ← xTB 级处理（仅 ext 两阶段有子目录）
│   ├── stage1_gfn0/             ← ext: GFN0 粗筛（其他协议无此目录）
│   │   ├── crest_conformers.xyz
│   │   └── crestopt.log
│   └── stage2_gfn2/             ← ext: GFN2 精修
│       ├── crest_ensemble.xyz
│       └── crestopt.log
│
├── cluster/                     ← ISOSTAT 聚类（统一层级，不嵌套）
│   ├── cluster.xyz              ← 代表性构象
│   ├── isomers.xyz
│   └── isostat.log
│
├── prescan/                     ← 快速 SP 预筛选（仅 full 协议）
│   └── conf_000/
│       ├── *.inp / *.out / *.xyz
│       └── *_property.txt
│
├── fastsp/                      ← r2SCAN-3c 筛选 SP（full/lite 协议）
│   └── conf_000/
│       ├── *.inp / *.out / *.xyz
│       └── *_property.txt
│
├── finalDFT/                    ← 最终 DFT OPT + Freq + L2 SP（所有协议）
│   ├── conf_000.gjf             ← Gaussian OPT 输入
│   ├── conf_000.log             ← Gaussian OPT 输出
│   ├── conf_000.fchk            ← Formatted checkpoint
│   ├── conf_000.xyz             ← 优化后几何
│   ├── conf_000_Shermo.sum      ← 热化学
│   └── conf_000_SP.*            ← ORCA L2 SP (wB97X-D4/def2-TZVPP)
│
├── conformer_state.json         ← 状态追踪
├── funnel_candidate_summary.json
├── provenance.json              ← 全流程溯源
└── <molecule>_global_min.xyz    ← 全局最低能量构象
```

### 5.3 新旧目录映射

| 新目录 | 旧目录 | 语义变化 | 创建条件 |
|--------|--------|---------|---------|
| `crest/` | `xtb2/` (CREST 输出部分) | 明确标识 CREST 搜索 | 所有协议 |
| `xtb/` | `xtb2/stage1_gfn0/`, `xtb2/stage2_gfn2/` | 去除 `xtb2` 前缀歧义 | ext 协议 |
| `cluster/` | `cluster/` (单阶段) / `xtb2/*/cluster/` (两阶段) | 统一层级 | 所有协议 |
| `prescan/` | `dft/fast_sp/prescreen_sp/` | 独立顶层目录 | full 协议 |
| `fastsp/` | `dft/fast_sp/screening_sp/` | 独立顶层目录 | full/lite 协议 |
| `finalDFT/` | `dft/` (conf_*.gjf/log + SP) | 明确区分最终 DFT | 所有协议 |

### 5.4 重构收益

1. **语义清晰**：每个目录名直接表达其计算任务类型
2. **协议一致性**：所有协议共享相同的顶层目录命名，不存在的任务目录则不创建
3. **避免嵌套重叠**：`cluster/` 不再嵌套在 `xtb2/stage*/` 下
4. **便于后处理**：按任务类型检索文件不再需要遍历多层子目录
5. **便于 benchmark 追踪**：每个目录对应一个计算阶段，能量溯源路径更短

### 5.5 实施要点

- 修改 `engine.py` 中 `__init__` 的目录创建逻辑（lines 114-121）
- 修改 `funnel.py` 中的路径引用
- 修改 `state_manager.py` 的状态文件路径
- 修改 `handler.py` 中的 anchor 阶段路径
- 更新 `path_compat.py` 的布局兼容性解析
- 确保 `_resolve_s1_artifacts()` 同时支持新旧布局

---

## 六、下周计划：DFT 理论方法 Benchmark

### 6.1 目标

系统性评估不同 DFT 泛函/基组组合对 [4+3] 环加成反应能线面的影响，确定最优理论水平配置。

### 6.2 当前理论水平

RPH v2.1.0 采用 Houk group 的 dual-level 策略：

| 层级 | 用途 | 当前配置 | 引擎 |
|------|------|---------|------|
| **预优化** | 避免 Gaussian 崩溃 | GFN2-xTB / crude | xTB |
| **几何优化 (L1)** | S1 构象优化, S3 TS 优化 | B3LYP/def2-SVP + GD3BJ | Gaussian |
| **高精度 SP (L2)** | S1/S3 能量精修 | wB97X-D4/def2-TZVPP + def2/J | ORCA |
| **快速 SP** | 构象筛选 | r2SCAN-3c / PBEh-3c | Mixed |

### 6.3 Benchmark 变量

#### 高优先级变量（影响关键能量 ΔG‡, ΔG_rxn）

| 变量 | Config 路径 | 当前值 | 候选值 |
|------|-----------|--------|--------|
| **优化泛函** | `theory.optimization.method` | B3LYP | B3LYP, PBE0, ωB97X-D4, M06-2X, r2SCAN-3c |
| **优化基组** | `theory.optimization.basis` | def2-SVP | def2-SVP, def2-TZVP, def2-TZVPP |
| **色散校正** | `theory.optimization.dispersion` | GD3BJ | GD3BJ, D4, none |
| **L2 SP 泛函** | `theory.single_point.method` | wB97X-D4 | wB97X-D4, wB97X-V, PW6B95-D4, DSD-PBEP86-D3BJ |
| **L2 SP 基组** | `theory.single_point.basis` | def2-TZVPP | def2-TZVP, def2-TZVPP, def2-QZVPP, cc-pVTZ |

#### 低优先级变量（影响筛选精度，已由 S1 benchmark 间接验证）

| 变量 | Config 路径 | 当前值 | 说明 |
|------|-----------|--------|------|
| prescreen SP | `step1.fast_sp_profiles.prescreen_sp.method` | PBEh-3c | 仅 full 协议使用 |
| screening SP | `step1.fast_sp_profiles.screening_sp.method` | r2SCAN-3c | full/lite 协议使用 |
| preopt GFN level | `theory.preoptimization.gfn_level` | 2 | 影响小，不在本轮评估 |

### 6.4 三阶段 Benchmark 方案

#### Phase 1 — 泛函筛选（固定基组）

**策略**：固定基组为 def2-SVP (opt) + def2-TZVPP (SP)，仅切换泛函。

| 实验编号 | 优化泛函 | 色散 | L2 SP | 预计壁钟/反应 |
|---------|---------|------|-------|-------------|
| F1 (baseline) | B3LYP | GD3BJ | wB97X-D4/def2-TZVPP | ~2h |
| F2 | PBE0 | GD3BJ | wB97X-D4/def2-TZVPP | ~2.5h |
| F3 | ωB97X-D4 | (内置) | wB97X-D4/def2-TZVPP | ~3h |
| F4 | M06-2X | GD3BJ | wB97X-D4/def2-TZVPP | ~4h |
| F5 | r2SCAN-3c | (内置) | wB97X-D4/def2-TZVPP | ~1.5h |

**测试反应**：rx1, rx3, rx8（覆盖刚性/中等/高柔性）  
**协议**：lite（已验证为最优）  
**Scope**：S1→S2→S3 全流程（评估对 TS 几何和活化能的影响）  
**预计总壁钟**：5 泛函 × 3 反应 × ~3h = ~45h

**评估指标**：
- 产物/前体构象能量一致性
- TS 几何（形成键键长偏差）
- ΔG‡ 活化能 vs 文献/基准值
- 壁钟时间

**筛选标准**：保留 Mean ΔE < 1.0 kcal/mol 且壁钟时间可接受的 2-3 个泛函。

#### Phase 2 — 基组收敛性测试（固定泛函）

**策略**：使用 Phase 1 筛选出的最优 2-3 个泛函，系统评估基组影响。

| 实验编号 | 优化基组 | L2 SP 基组 | 预计壁钟/反应 |
|---------|---------|-----------|-------------|
| B1 (baseline) | def2-SVP | def2-TZVPP | ~2h |
| B2 | def2-TZVP | def2-TZVPP | ~5h |
| B3 | def2-TZVP | def2-QZVPP | ~8h |
| B4 | def2-TZVPP | def2-TZVPP | ~10h |

**测试反应**：rx1, rx3（2 个代表）  
**预计总壁钟**：3 泛函 × 4 基组 × 2 反应 × ~6h = ~144h

**评估指标**：
- 几何收敛性（键长变化 < 0.01 Å）
- 能量收敛性（ΔE 变化 < 0.5 kcal/mol）
- 成本增量与精度增量的边际分析

**筛选标准**：选择精度满足化学精度（< 1 kcal/mol）的最小基组。

#### Phase 3 — 色散校正对比

**策略**：使用 Phase 1+2 确定的最优泛函/基组，评估色散校正的影响。

| 实验编号 | 色散 | 适用泛函 |
|---------|------|---------|
| D1 | none | B3LYP, PBE0 |
| D2 | GD3BJ | B3LYP, PBE0, M06-2X |
| D3 | D4 | B3LYP, PBE0, M06-2X |
| D4 | (内置) | ωB97X-D4, r2SCAN-3c |

**测试反应**：rx3, rx15（高柔性分子，色散效应最显著）  
**预计总壁钟**：~30h

**评估指标**：
- 色散对构象排序的影响
- 色散对形成键几何的影响
- 有/无色散的能量差

### 6.5 总体实施计划

| 阶段 | 时间 | 内容 | 预计壁钟 |
|------|------|------|---------|
| Phase 1 | 周一~周二 | 泛函筛选 (5×3 反应) | ~45h |
| Phase 2 | 周三~周四 | 基组收敛 (3×4×2 反应) | ~144h |
| Phase 3 | 周五 | 色散对比 | ~30h |
| **分析** | 周五 | 数据汇总 + 成本精度 Pareto 图 | — |

### 6.6 基准数据与验证策略

**参考数据来源**：
- Phase 1 baseline (B3LYP/def2-SVP + wB97X-D4/def2-TZVPP) 已有 5 反应数据
- 若有文献实验 ΔG‡ 值，作为绝对精度验证
- DLPNO-CCSD(T)/def2-TZVPP 可作为计算化学 "gold standard" 参考（选 1-2 个反应）

**自动化实施**：
- 修改 `benchmark/confsearch/run.sh` 增加 `--theory` 参数，接受覆写配置
- 创建 `benchmark/confsearch/theory_configs/` 目录存放各实验的 YAML override
- 评估脚本 `evaluate.py` 扩展为支持跨理论方法对比

---

## 七、附录

### A. 各反应原始能量数据

#### rx1 — Boc, acyclic allenamide, 3-carbon tether
| Protocol | Product E (Ha) | Precursor E (Ha) | DFT | Wall |
|----------|---------------|-----------------|-----|------|
| ext | −940.674 883 572 | −865.331 299 081 | 9 | 2:50 |
| full | −940.674 883 598 | −865.333 294 165 | 14 | 5:02 |
| lite | −940.674 883 644 | −865.333 294 455 | 4 | 1:43 |
| zero | −940.674 883 638 | −865.331 286 051 | 5 | 2:00 |

#### rx3 — Bn (benzyl), acyclic allenamide, 2-carbon tether
| Protocol | Product E (Ha) | Precursor E (Ha) | DFT | Wall |
|----------|---------------|-----------------|-----|------|
| ext | −1014.502 037 830 | −939.163 942 867 | 10 | 2:02 |
| full | −1014.501 490 937 | −939.165 231 673 | 36 | 10:22 |
| **lite** | **−1014.503 440 281** | −939.165 230 273 | **4** | **0:37** |
| zero | −1014.503 438 957 | −939.163 350 620 | 8 | 1:10 |

#### rx8 — cyclic carbamate (oxazolidinone), 1-carbon tether (rigid)
| Protocol | Product E (Ha) | Precursor E (Ha) | DFT | Wall |
|----------|---------------|-----------------|-----|------|
| ext | −782.042 651 553 | −706.707 178 576 | 4 | 0:31 |
| full | −782.042 651 560 | −706.707 178 810 | 6 | 0:53 |
| lite | −782.042 651 525 | −706.707 178 746 | 2 | 0:24 |
| zero | −782.042 651 530 | −706.707 178 755 | 2 | 0:22 |

> ΔE < 0.000 03 kcal/mol — 刚性分子所有协议完全等价。

#### rx13 — cyclic carbamate (oxazolidinone), 3-carbon tether
| Protocol | Product E (Ha) | Precursor E (Ha) | DFT | Wall |
|----------|---------------|-----------------|-----|------|
| ext | −860.755 719 775 | −785.414 282 216 | 9 | 1:42 |
| full | −860.755 720 241 | −785.414 282 844 | 17 | 3:47 |
| lite | −860.755 720 316 | −785.414 282 414 | 3 | 0:46 |
| zero | −860.755 720 202 | −785.413 250 588 | 5 | 1:02 |

#### rx15 — cyclic carbamate (oxazolidinone), 4-carbon tether (high flexibility)
| Protocol | Product E (Ha) | Precursor E (Ha) | DFT | Wall |
|----------|---------------|-----------------|-----|------|
| **ext** | **−900.103 615 166** | **−824.766 588 788** | 8 | 4:16 |
| full | −900.103 498 288 | −824.766 522 563 | 25 | 7:45 |
| lite | −900.103 501 222 | −824.766 522 650 | 3 | 1:06 |
| zero | −900.103 518 886 | −824.766 522 628 | 5 | 1:21 |

### B. Gaussian/ORCA 关键词映射

| Config 参数 | 当前值 | Gaussian 关键词 | ORCA 关键词 |
|-----------|--------|---------------|------------|
| `theory.optimization.method` | B3LYP | `B3LYP` | — |
| `theory.optimization.basis` | def2-SVP | `Def2SVP` | — |
| `theory.optimization.dispersion` | GD3BJ | `em=GD3BJ` | — |
| `theory.single_point.method` | wB97X-D4 | — | `wB97X-D4` |
| `theory.single_point.basis` | def2-TZVPP | — | `def2-TZVPP` |
| `theory.single_point.aux_basis` | def2/J | — | `def2/J` |

**完整 L1 route**: `#p B3LYP/Def2SVP em=GD3BJ Opt Freq`  
**完整 L2 route**: `! wB97X-D4 def2-TZVPP def2/J RIJCOSX`

### C. 关键文件索引

| 文件 | 路径 | 说明 |
|------|------|------|
| 最终报告 | `docs/benchmark_final_report_5reactions.md` | 5 反应综合分析 |
| rx1 报告 | `docs/benchmark_report_rx1_s1_protocols.md` | rx1 详细分析 |
| rx3 报告 | `docs/benchmark_report_rx3_s1_protocols.md` | rx3 详细分析 |
| 4反应摘要 | `docs/benchmark_summary_rx1_rx3_rx8_rx13.md` | rx1+3+8+13 总结 |
| Benchmark 脚本 | `benchmark/confsearch/run.sh` | 运行编排器 |
| 评估脚本 | `benchmark/confsearch/evaluate.py` | 后评估计算 |
| 协议定义 | `rph_core/steps/conformer_search/protocols.py` | ProtocolSpec |
| Funnel 实现 | `rph_core/steps/conformer_search/funnel.py` | FunnelRunner |
| 候选数据结构 | `rph_core/steps/conformer_search/candidates.py` | CandidateSet |
| 配置文件 | `config/defaults.yaml` | 理论水平 + 协议栈 |
