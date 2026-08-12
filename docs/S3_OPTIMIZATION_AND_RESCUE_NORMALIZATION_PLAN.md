# S3 流程优化 + 救援方法族矩阵方案（修订 v2）

> 状态：**方案稿 v2（Plan Only，未实施）** — 按 2026-08-05 评审意见修订
> 关联：`RPH_V4_S3_S4_UNIFIED_REFINEMENT_PLAN.md`（本方案在其之上做 S3 流程优化与命名规范化）
> 证据基准：`RPH_Test_Results/s2_gfn2_benchmark_alpb_v2/rx_id_1`（S3 首轮全量运行，7 结构 5 成 2 败）

---

## 目录

1. [评审意见与修订对照](#1-评审意见与修订对照)
2. [运行数据分析结论（方案依据）](#2-运行数据分析结论方案依据)
3. [目标](#3-目标)
4. [S3.x 规范化步骤体系（v2）](#4-s3x-规范化步骤体系v2)
5. [救援方法族与救援矩阵（v2 核心变更）](#5-救援方法族与救援矩阵v2-核心变更)
6. [S3 计算流程优化（v2）](#6-s3-计算流程优化v2)
7. [Manifest / 日志 / UI / 事件集成](#7-manifest--日志--ui--事件集成)
8. [向后兼容](#8-向后兼容)
9. [分阶段实施计划](#9-分阶段实施计划)
10. [验证门](#10-验证门)
11. [关键决策点（待确认）](#11-关键决策点待确认)
12. [附录：相关代码位置](#12-附录相关代码位置)

---

## 1. 评审意见与修订对照

| # | 评审意见 | 修订措施 | 所在节 |
|---|---|---|---|
| 1 | Warmup cycle 过小，应给长到 LooseOpt 收敛限 | `warmup_max_cycles` INT 8→**40** / TS 12→**50**；ORCA 原生提前退出保护易收敛点 | §4.2 S3.1, §6 W1 |
| 2 | 稳定点不需这么多轮；建议全局统一 ~60 轮收敛限 | **已验证**：成功结构实际收敛 3-17 轮（precursor 11 / product_major 8 / product_major_int 17 / product_minor_001 11 / product_major_ts 3）→ `max_cycles_*` 全局统一 **60** | §2.4, §6 W2 |
| 3 | S3.3 需要更多详细子分类层级 | S3.4 救援细分为 **S3.4.0 失败诊断 / S3.4.1 重启族 / S3.4.2 模式族 / S3.4.3 精确 Hessian 族** | §4.3 |
| 4 | S3.4 意义不明 | 原"Canonical"重命名 **S3.5 Canonical Selection**，明确为"选择"非"计算"（纯内存排序+落地，非 QC） | §4.2 S3.5 |
| 5 | FREQ 为何不是独立层级 | **S3.3 Frequency & Classification** 独立成层（主几何 Freq + 驻点分类）；canonical Freq 作为 S3.5 子记录 | §4.2 S3.3 |
| 6 | 救援顺序太复杂、状态机过多；R1/R2/R3 应定义为不同救援方法而非 S3 子流程，与不同报错构成救援矩阵 | **R1/R2/R3 重定义为救援方法族**（Restart / Mode / Exact），废除逐角色嵌套决策树，改为 **错误类型 × 结构类型 二维救援矩阵**，单元格 = 有序方法列表 | §5 |

**第二轮评审修订（2026-08-05）**：

| # | 评审意见 | 修订措施 | 所在节 |
|---|---|---|---|
| 7 | R1 不应读取旧 Hessian，应**重算一次精确 Hessian** 再继续（更新曲率、避免震荡环境） | R1 由 `read_hessian_restart`（InHess Read）改为 **`fresh_hessian_restart`**（`Calc_Hess true` 重算精确 Hessian，默认步长）；**O2 Hessian 复用优化废除** | §5.1, §5.3, §6 O2 |
| 8 | 起始结构无虚频时，应看后续优化能否找到虚频并确认方向，回退 warmup 无意义 | F3 单元格由 `warmup_reseed` 改为 **R1 继续下降 + 周期 `Recalc_Hess` 监测虚频涌现** → 模式确认方向后 R2 定向 | §5.2, §5.3 |
| 9 | 回退原始 S2 种子无意义 | `original_seed_fallback` **删除**（checkpoint 已在下降轨迹上，回退种子必然重演失败） | §5.1 |
| 10 | R2 mode 基本合理 | 保留（模式选择 + Trust 0.03 保守脊线跟随） | §5.1 |
| 11 | R3 calcall 不应加 Trust 0.05，默认步长即可，否则耗时爆炸 | R3 **移除 Trust 覆盖**（精确 Hessian 已保证步长质量，限步长徒增轮数） | §5.1, §5.3 |

**第三轮评审修订（2026-08-05，终版）**：

| # | 评审意见 | 修订措施 | 所在节 |
|---|---|---|---|
| 12 | R2 trust 0.03 太小，应到 0.15，否则计算时间太长 | `ts_mode_directed` trust **0.03 → 0.15**（与主 OptTS 一致；模式定向 + 大步长沿脊线快速收敛） | §5.1, §5.3 |
| 13 | recalc = 5（约 6 次精确 Hessian 评估仍不收敛必有问题） | R1 系列 `recalc_hessian_interval` **统一 5**；方法级预算 = MaxIter 30 / Recalc 5 ≈ **6 次精确 Hessian 评估**，仍不收敛即判败升级下一方法 | §5.1 |
| 14 | F5/F6 直接报错退出，具体问题具体分析 | F5/F6 单元格由"R0 级重试/重跑"改为 **直接失败退出（S3.7）**：不重试、不进入矩阵，记录错误待人工分析 | §5.3 |

---

## 2. 运行数据分析结论（方案依据）

> 数据来源：`rx_id_1/S3_LowLevel/product_minor_001_{int,ts}/diagnostics/attempt_002_*` ORCA 输出 + 5 个成功结构 opt 输出。

### 2.1 收敛轨迹：无能量震荡，梯度未收敛

**TS 30 轮能量（FINAL SINGLE POINT ENERGY，Eh）**：`-939.4124 → -939.4173`，**每步单调下降**，无反弹；末段 Δ≈1e-4 Eh/步。

**INT 30 轮能量**：`-939.4133 → -939.4272`，单调下降（仅 4 处 ~1e-4 微小上翘，属 BFGS 正常过冲，非震荡）。

**梯度轨迹（决定性证据）**：

| 结构 | RMS 梯度范围 (Eh/bohr) | 收敛阈值 | 状态 |
|---|---|---|---|
| TS | 0.66e-3 ~ 1.37e-3（7~14×） | 1e-4 | **钉死不降**，`RMS step 0.00724` 反复重复 |
| INT | 1.9e-3 → 0.47e-3（稳定下降） | 1e-4 | 收敛中，仅步数不够 |

### 2.2 首步 Hessian：TS 是二阶鞍点（2 个负本征值）

首步 Hessian 最低本征值（Hartree/bohr²）：

```
-0.024849365   -0.003067565   0.002320652   0.003237920   0.004097098
 └─ mode 0 ─┘   └─ mode 1 ─┘
 （反应坐标）   （浅虚频 ~30 cm⁻¹ 量级）
```

- 每次 Hessian 重算点（cycle 1/8/16/24）均输出 `Hessian has 2 negative eigenvalues`，**mode 1 从未变正**。
- **P-RFO 算法限制**：OptTS 沿 mode 0 最大化、沿所有其他模式（含负曲率 mode 1）"最小化"；负曲率方向的"最小化"在调和模型中无驻点，仅靠 trust radius 钳制步长 → 爬行。**标准 OptTS 在真正的二阶鞍点上结构性无法收敛，cycle 数再多也没用。**
- 附加异常：第 1 轮 P-RFO 的 λ 二分搜索失败（`UNABLE TO DETERMINE LAMBDA / Bisect method failed`）。

### 2.3 Warmup 实测：LooseOpt 收敛快于 cap（评审点 1 依据）

`product_minor_001_int` warmup（LooseOpt, TolRMSG=5e-4）：**第 8 轮 RMSG 0.00032 < 5e-4，已收敛**。说明 LooseOpt 对多数种子在 ~10 轮内收敛，cap 增大（40/50）不拖慢易收敛点（ORCA 原生提前退出），但给差种子留足空间。

### 2.4 稳定点收敛轮数实测（评审点 2 依据）

| 结构 | 类型 | 实际几何收敛轮数 | 原 cap |
|---|---|---|---|
| precursor | minimum | **11** | 100 |
| product_major | minimum | **8** | 100 |
| product_major_int | intermediate | **17** | 30 |
| product_minor_001 | minimum | **11** | 100 |
| product_major_ts | TS | **3** | 30 |
| product_minor_001_int | intermediate | 未收敛（30 轮梯度 0.00047） | 30 |
| product_minor_001_ts | TS | 未收敛（30 轮梯度钉死） | 30 |

**结论**：收敛点普遍 ≤17 轮（TS 可低至 3 轮）——**cap 只是失败检测上界，不是实际耗轮**。全局 60 cap 的代价仅在"该失败但没失败"的结构上多烧时间；而慢 INT（需 ~50-60 轮）获得足够空间。配合 O1 二阶鞍点提前中止（~cycle 16-24 即断），TS 的 60 cap 实际被提前终止覆盖。

### 2.5 对方案的影响

- **全局 60 cap 成立**（评审点 2）；TS 二阶鞍点靠 O1 提前中止 + 救援矩阵处理。
- **Warmup cap 提升成立**（评审点 1）。
- **R1 mode-directed OptTS 对本案例只能部分解决**：mode 0 会被选中，mode 1 残留需 R2（mode 位移/破鞍点）或 R3（精确 Hessian）接力——这正是救援矩阵设计的动因。

---

## 3. 目标

1. **S3 子步骤规范化**：S3.0~S3.7 标准步骤体系（FREQ 独立、救援子层级、Canonical 语义明确），统一日志/UI/manifest/事件/恢复。
2. **救援方法族矩阵**：R1/R2/R3 重定义为方法族；错误类型 × 结构类型的二维矩阵取代逐角色嵌套状态机。
3. **计算流程优化**：全局 60 cap、Warmup 至 LooseOpt 收敛、二阶鞍点提前中止、救援 Hessian 复用、角色并发、SP 缓存持久化。
4. **全量回归**：5 个成功结构零重算，2 个失败结构经救援矩阵收敛，`steps[]` 契约完整。

---

## 4. S3.x 规范化步骤体系（v2）

### 4.1 单一事实源：`StepCode` 枚举

新建 `rph_core/steps/refinement/step_codes.py`，全 S3 步骤身份的**唯一来源**：

```python
class StepCode(str, Enum):
    S3_0_PREFLIGHT = "S3.0"   # 预检：XYZ/电荷/自旋/forming_bonds/目录/provenance
    S3_1_WARMUP    = "S3.1"   # 条件化热身（仅 INT/TS）：LooseOpt + 键约束 → LooseOpt 收敛
    S3_2_PRIMARY   = "S3.2"   # 主 Opt/OptTS（全局 60 cap，二阶鞍点可提前中止）
    S3_3_FREQ      = "S3.3"   # 独立 Freq + 驻点分类（minimum/TS/INT + alignment）
    S3_4_RESCUE    = "S3.4"   # 救援矩阵（子层级 S3.4.0~S3.4.3，见 §4.3）
    S3_5_CANONICAL = "S3.5"   # Canonical Selection：多候选→唯一 canonical.xyz（纯选择，非 QC）
    S3_6_PROPERTIES = "S3.6"  # r2SCAN-3c SP + 复合热化学 + ml_usability 9-key
    S3_7_FINALIZE  = "S3.7"   # manifest 写盘 + 陈旧归档 + summary + fail-fast 信号
```

> S4 建议同构（`S4.0~S4.7`）——见 §11 决策点 3。

### 4.2 每步的规范化契约

| code | name | entry_gate | exit_gate（分岔） | outputs |
|---|---|---|---|---|
| S3.0 | preflight | 结构请求就绪 | 校验通过 → S3.1 | 校验记录、provenance.json |
| S3.1 | warmup | INT/TS 角色（precursor/product 直通 S3.2） | warmup_xyz 存在（complete 或 partial）→ S3.2 | warmup/opt.xyz、约束明细 |
| S3.2 | primary | warmup_xyz（或原种子） | converged → S3.3；maxiter → S3.4；二阶鞍点提前中止 → S3.4；其他失败 → S3.7 | opt.xyz/out、stop_reason |
| S3.3 | freq | 主优化收敛几何 | 分类完成 → S3.5（无救援需求）或 S3.4（分类异常） | freq/out、hess、分类结果 |
| S3.4 | rescue | S3.2/S3.3 失败或分类异常 | 任一方法族产出有效候选 → S3.5；耗尽 → S3.7 | 方法执行记录（矩阵单元格） |
| S3.5 | canonical | 有效候选存在（主优化或救援） | canonical.xyz 落地 + canonical Freq（复用或重跑）→ S3.6 | canonical.xyz、选择依据 |
| S3.6 | properties | canonical 完成 | ml_usability 判定 → S3.7 | sp_energy、thermo、ml 9-key |
| S3.7 | finalize | 所有前置 | manifest 写盘 → 批次 fail-fast | manifest、summary、健康报告 |

**FREQ 独立理由（评审点 5）**：Freq 是独立 QC 作业类型（B97-3c，~15min/次），在主几何与 canonical 几何各跑一次；独立成层使分类逻辑、Freq 失败处理、canonical Freq 复用/重跑都有明确归属，日志/UI 可直接呈现"正在做频率分析"。

**Canonical 语义（评审点 4）**：S3.5 是**选择而非计算**——多候选（主优化 + 救援方法产物）按 `stationary_rank → mode_rank → hessian_index → converged → alignment_score → gradient_norm` 排序取唯一，落地 canonical.xyz，并决定 canonical Freq 是复用主 Freq（几何不变）还是重跑。全程无新 QC 计算（Freq 复用判定除外），毫秒级。

### 4.3 S3.4 救援子分类层级（评审点 3）

```
S3.4 Rescue
 ├─ S3.4.0  失败诊断     : failure_type 分类 (F1~F6) + trigger 记录 → 查救援矩阵
 ├─ S3.4.1  R1 重启族    : 重算精确 Hessian 后继续优化（Calc_Hess true，周期刷新可配）
 ├─ S3.4.2  R2 模式族    : TS_Mode 定向 OptTS / 虚频 mode ± 位移 / 鞍点破除位移
 └─ S3.4.3  R3 精确族    : Recalc_Hess 1 (CalcAll) 兜底
```

矩阵决定 S3.4.x 的触发顺序与参数（§5.3），每次方法执行记录为标准 rescue record。

### 4.4 步骤流转图（v2）

```
S3.0 Preflight ──► S3.1 Warmup(INT/TS) ──► S3.2 Primary(60轮cap)
                                              │
                          ┌───────────────────┼───────────────────┐
                          ▼                   ▼                   ▼
                     converged            maxiter 耗竭/       其他失败
                          │               二阶鞍点提前
                          ▼                   │                   │
                     S3.3 Freq&分类            ▼                   │
                          │               S3.4 Rescue             │
                     ┌────┴─────┐        ├─S3.4.0 诊断            │
                     ▼         ▼        ├─S3.4.1 R1 重启族        │
                 分类正常   分类异常     ├─S3.4.2 R2 模式族        │
                     │         │        └─S3.4.3 R3 精确族        │
                     └────┬────┘              │                   │
                          ▼                   ▼                   │
                     S3.5 Canonical Selection ◄─ 有效候选           │
                          │                   │ (耗尽→failed)      │
                          ▼                   ▼                   ▼
                     S3.6 Properties      S3.7 Finalize(failed)
                          │
                          ▼
                     S3.7 Finalize
```

---

## 5. 救援方法族与救援矩阵（v2 核心变更）

### 5.1 R1/R2/R3 重定义：救援方法族（非子流程序号）

| 方法族 | 名称 | 原子方法 | 固定参数 |
|---|---|---|---|
| **R1 Restart 族** | 重算精确 Hessian 后继续优化 | `fresh_hessian_restart`（`Calc_Hess true` 在 checkpoint 重算精确 Hessian，**不用**主运行 BFGS 更新的旧 Hessian；`Recalc_Hess 5` 周期刷新） | Calc_Hess true, Recalc_Hess 5, max_cycles 30, 默认步长；**方法级预算 ≈6 次精确 Hessian 评估**（30/5），仍不收敛即判败升级 |
| **R2 Mode 族** | 沿虚频方向 | `ts_mode_directed`（TS_Mode {M n} 定向 OptTS，模式来自 checkpoint 处**新算**的局部 Hessian）、`mode_displacement`（虚频 mode ± 位移）、`saddle_break`（mode-1 ± 位移破除二阶鞍点） | **trust 0.15**（与主 OptTS 一致）, max_cycles 12, overlap ≥0.35, margin ≥0.08 |
| **R3 Exact 族** | 每步精确 Hessian 兜底 | `calcall_opt`（Recalc_Hess 1，**默认步长**，不设 Trust 覆盖——精确 Hessian 已保证步长质量，限步长徒增轮数） | max_cycles 30, recalc 1, 默认 trust |

> 一个方法族 = 一类**方法论**，不是流程阶段。同一族内多个方法按序尝试；不同族的组合由矩阵单元格指定。

### 5.2 失败类型分类（任务报错）

| 代码 | failure_type / 证据 | 来源 |
|---|---|---|
| F1 | `geometry_optimization_not_converged` / `native_maxiter_checkpoint` | `_is_geometry_nonconvergence` |
| F2 | `Hessian ≥2 负本征值`（连续 2 个 Hessian 重算点确认） | ORCA 输出解析（新增） |
| F3 | TS 无虚频（checkpoint Hessian 无负本征值） | `_select_target_ts_mode` 结果（**策略**：不换种子——继续下降，周期 Hessian 重算中监测虚频涌现，模式出现且方向确认后转 R2 定向） |
| F4 | 极小点带虚频（precursor/product 出现虚频） | `_classify_structure` |
| F5 | `scf_not_converged` | failure_classifier |
| F6 | 崩溃 / 超时 / 其他 | returncode / timed_out |

### 5.3 救援矩阵（核心交付：错误类型 × 结构类型 → 有序方法列表）

```
                              │   TS           │   INT          │   minimum (precursor/product)
──────────────────────────────┼────────────────┼────────────────┼──────────────────────────────
F1 未收敛 (maxiter)           │ R1: fresh_hessian (OptTS)     │ R1: fresh_hessian (Opt)       │ R1: fresh_hessian (Opt)
                              │  → R2: ts_mode_directed       │  → R3: calcall_opt            │  → R3: calcall_opt
                              │  → R3: calcall_opt            │                              │
──────────────────────────────┼────────────────┼────────────────┼──────────────────────────────
F2 高阶鞍点 (≥2 虚频)         │ R2: saddle_break (mode-1 ±)   │ —（INT 不应出现）             │ —（minimum 不应出现）
                              │  → R2: ts_mode_directed       │                              │
                              │  → R3: calcall_opt            │                              │
──────────────────────────────┼────────────────┼────────────────┼──────────────────────────────
F3 TS 无虚频                  │ R1: fresh_hessian (OptTS,     │ —                             │ —
                              │     Recalc_Hess 5 监测虚频涌现)│                              │
                              │  → R2: ts_mode_directed       │                              │
                              │  → R3: calcall_opt            │                              │
──────────────────────────────┼────────────────┼────────────────┼──────────────────────────────
F4 极小点带虚频               │ —                              │ R2: mode_displacement ±      │ R2: mode_displacement ±
──────────────────────────────┼────────────────┼────────────────┼──────────────────────────────
F5 SCF 不收敛                 │ 直接报错退出（S3.7 failed）：不重试、不进入矩阵，记录错误待人工分析         │
──────────────────────────────┼────────────────┼────────────────┼──────────────────────────────
F6 崩溃/超时                  │ 直接报错退出（S3.7 failed）：不重试、不进入矩阵，记录错误待人工分析         │
```

**执行语义**：单元格 = 有序方法列表，顺序执行，任一方法产出有效候选（`_candidate_is_valid_for_request`）即终止并进入 S3.5；耗尽则 S3.7(failed)。方法之间无嵌套分支、无角色专属状态机——**单表驱动，维护 = 改表单元格**。

**对比旧设计（评审点 6）**：

| 维度 | v1（废除） | v2（采用） |
|---|---|---|
| 结构 | 逐角色嵌套决策树（TS→前置判断→R1→R2→兜底…） | 二维矩阵查表（failure × kind → 方法列表） |
| R1/R2/R3 含义 | S3 内部轮次序号（level 字段） | **方法族**（Restart/Mode/Exact），与流程解耦 |
| 新增救援 | 修改分支逻辑 | 追加矩阵单元格 / 新增方法族 |
| 状态机数量 | 多（每个角色一条链） | **1 个**（矩阵查表 + 列表迭代） |
| 可测试性 | 每链一条测试 | 每单元格一条测试，可穷举 |

### 5.4 统一 rescue record（manifest 内）

```json
{
  "method_family": "R2",
  "method": "ts_mode_directed",
  "round_ordinal": 1,
  "trigger": {"failure_type": "F1", "reason": "native_maxiter_checkpoint", "primary_cycle_limit": 60},
  "status": "complete | failed | skipped",
  "input_geometry": "...", "output_geometry": "...",
  "energy_hartree": -939.417, "frequencies_cm1": [...],
  "mode_analysis": {"selected_mode": 0, "selected_overlap": 0.82, "mode_scores": [...]},
  "attempt_id": "...", "duration_seconds": 1731
}
```

---

## 6. S3 计算流程优化（v2）

### W1. Warmup 延长至 LooseOpt 收敛（评审点 1）

- `warmup_max_cycles`: INT 8→**40**，TS 12→**50**。
- ORCA 原生提前退出（收敛即停），实测 INT warmup 第 8 轮已收敛（RMSG 0.00032 < 5e-4）——易收敛点零损失，差种子获得足够空间。
- `allow_unconverged_geometry=True` 保留（partial 仍可作主优化起点，但不作为首选路径）。

### W2. 全局统一 60 轮 cap（评审点 2）

- `max_cycles_minimum` 100→**60**，`max_cycles_intermediate` 30→**60**，`max_cycles_ts` 30→**60**。
- 依据 §2.4：收敛点实际 3-17 轮，cap 仅是失败检测上界。统一语义 = 单数值心智模型，易维护。
- TS 二阶鞍点由 O1 提前中止（~cycle 16-24）覆盖，60 cap 不产生额外空耗。

### O1. 二阶鞍点提前中止 + `saddle_break`

- S3.2 首个 `Recalc_Hess` 点解析 `Hessian has N negative eigenvalues`；N≥2 且第 2 负本征值 < -1e-3，**连续 2 个重算点确认** → 立即中止主优化 → S3.4 矩阵查 F2 行。
- `saddle_break` = R2 族内方法：沿 mode-1 ± 位移（复用 `_build_mode_displaced_xyz`）打破高阶鞍点 → OptTS。

### O2. ~~救援 Hessian 复用~~（已废除，评审点 7）

- **废除**：R1 一律重算精确 Hessian（`Calc_Hess true`），不复用主运行 BFGS 更新的旧 `.hess`——旧 Hessian 在震荡/近平坦区质量退化，正是 R1 要解决的问题。
- R2 模式定向仍需 checkpoint 处**新算**的局部 Hessian（`_run_local_hessian_analysis` 的 FREQ 产物），属模式选择的必要输入，不是"复用"。

### O3. 角色感知并发

- S3.2 按角色分池：precursor/product（minima）并行；INT/TS 保守串行。`stage_scheduler` 的 `max_concurrent_by_role` 实际生效，`max_workers` ≥2。

### O4. SP 缓存持久化

- `_sp_cache` 磁盘持久化（geometry hash + method/basis/solvent 键控，run 根目录），跨 rx_id 复用共用片段。

### O5. 失败诊断健康报告

- S3.7 输出结构化报告：救援矩阵命中分布、各方法族成功率、平均方法数、cycle 分布、二阶鞍点清单 → 驱动批次决策。

---

## 7. Manifest / 日志 / UI / 事件集成

### 7.1 manifest 新增 `steps[]` 数组（每结构）

```json
"steps": [
  {"code": "S3.0", "name": "preflight", "status": "complete", "duration_s": 0.4},
  {"code": "S3.1", "name": "warmup", "status": "complete", "attempt_id": "attempt_001_warmup_opt", "duration_s": 742},
  {"code": "S3.2", "name": "primary_opt", "status": "failed", "attempt_id": "attempt_002_opt_ts", "stop_reason": "native_maxiter_checkpoint", "duration_s": 3521},
  {"code": "S3.3", "name": "freq", "status": "skipped", "reason": "primary_opt_failed", "duration_s": 0},
  {"code": "S3.4", "name": "rescue", "status": "complete", "sub_steps": [
      {"code": "S3.4.0", "failure_type": "F2", "trigger": "higher_order_saddle"},
      {"code": "S3.4.2", "method_family": "R2", "method": "saddle_break", "status": "complete", "attempt_id": "attempt_004_r2_saddle_break", "duration_s": 1200},
      {"code": "S3.4.2", "method_family": "R2", "method": "ts_mode_directed", "status": "complete", "attempt_id": "attempt_005_r2_mode_directed", "duration_s": 630}
  ], "duration_s": 1830},
  {"code": "S3.5", "name": "canonical", "status": "complete", "selection": {"rank": 1, "source": "r2_mode_directed"}, "freq_reused": false, "duration_s": 900},
  {"code": "S3.6", "name": "properties", "status": "complete", "sp_energy_hartree": -939.6438, "duration_s": 900},
  {"code": "S3.7", "name": "finalize", "status": "complete", "duration_s": 3}
]
```

### 7.2 日志与事件

- 日志前缀 `[S3.2 primary_opt]`；`rph_v4.log` 保持纯文本。
- 事件：`s3.4.r2_mode_directed.start|complete|failed`（step code + method + attempt_id）；废弃 `phase: wave2_pass2_pass3`。

### 7.3 UI 状态

- `ui_adapter.adapt_s3_structures` 增量：每结构步骤进度（S3.0→S3.7）+ 救援方法族徽标。现适配器命名中立，改动为纯增量。

### 7.4 恢复粒度

- `pipeline.state` 记录到 S3.x 步骤 + S3.4 方法级：resume 从失败步骤/方法继续。

---

## 8. 向后兼容

| 变更 | 兼容策略 |
|---|---|
| config 键改名（`ts_nonconvergence_rescue` → 救援矩阵 `matrix`/`methods`） | `config_loader` 读旧键 → 新键 + DeprecationWarning（一个版本窗口） |
| manifest `level` → `round_ordinal` + `method_family` | `manifest_io` 读适配器映射 |
| 旧 strategy id（`ts_nonconvergence_l1_mode_directed`） | 读适配器映射到 `{method_family, method}` 对 |
| `_run_pass2_rescue_one` / 逐角色嵌套逻辑 | 删除，统一由矩阵驱动器执行（已全库验证无外部引用） |

---

## 9. 分阶段实施计划

| 阶段 | 内容 | 交付物 | 验证 |
|---|---|---|---|
| **P1 命名底座** | `step_codes.py`（S3.0~S3.7）+ `StepSpec`；manifest `steps[]` writer | 新枚举、schema v2（读兼容 v1） | `test_v4_protocol_contract` 扩展 |
| **P2 矩阵引擎** | 失败分类器（F1~F6）+ 救援矩阵表 + 方法族执行器（R1/R2/R3 原子方法）；删除旧嵌套逻辑 | `rescue_matrix` 模块 + 执行器 | 新测试：矩阵逐单元格穷举；`test_refinement_pass2_rescue` 改造 |
| **P3 流程参数** | W1 warmup 40/50；W2 全局 60；O1 二阶鞍点提前中止 + saddle_break；O2 Hessian 复用 | 配置 + 检测逻辑 | 新增测试：提前中止、Hessian 复用 |
| **P4 运维** | O3 角色并发；O4 SP 缓存持久化；O5 健康报告；UI/事件步骤化 | 报告生成器 + 缓存层 | rx_id=1 回归：2 失败结构经矩阵收敛，`steps[]` 完整 |

---

## 10. 验证门

```bash
python -m py_compile rph_core/steps/refinement/engine.py rph_core/steps/refinement/step_codes.py
pytest -q tests/test_v4_protocol_contract.py tests/test_v4_checkpoint.py tests/test_v4_stage_calculator.py tests/test_v4_ui.py tests/test_refinement_pass2_rescue.py
python scripts/ci/check_imports.py rph_core
bash scripts/run_s3_serial.sh data/reaxys_cleaned.csv RPH_Test_Results/s2_gfn2_benchmark_alpb_v2   # 回归：5 结构复用 + 2 失败结构经矩阵收敛
```

---

## 11. 关键决策点（待确认）

1. **全局 60 cap 的 TS 例外**：TS 二阶鞍点靠 O1 提前中止（~16-24 轮）兜底，60 cap 仅在一阶鞍点慢收敛时实际消耗——确认接受"统一 60，不加例外"。
2. **救援矩阵的数据驱动程度**：矩阵表放 `config/defaults.yaml`（科学评审可见、免改代码）还是代码常量（类型安全、测试友好）？建议：**矩阵结构在代码（dataclass + 测试穷举），参数在 config**。
3. **S3.x 是否延伸到 S4**：同一 `RefinementEngine` 驱动 S4，建议 `S4.0~S4.7` 同构（profile 参数不同）——确认纳入。
4. **~~R1 周期刷新默认值~~（已定，评审 13）**：R1 统一 `Recalc_Hess 5`（≈6 次精确 Hessian 评估为方法级预算）；如个别震荡案例需更保守，可按需调小间隔。
5. **~~F5/F6 归属~~（已定，评审 14）**：直接报错退出（S3.7 failed），不重试、不进入矩阵，具体问题具体分析。

---

## 12. 附录：相关代码位置

| 关注点 | 位置 |
|---|---|
| 生产 pass2 入口（改造为矩阵驱动器） | `rph_core/steps/refinement/engine.py` L253, L1426 `_run_pass2_rescue` |
| 死代码（删除/并入） | `engine.py` L2031 `_run_pass2_rescue_one`；L2062 `_run_ts_rescue`；L2132 `_run_int_rescue`；L2228 `_run_minimum_rescue`；L2300 `_maybe_run_ts_endpoint_discovery` |
| 现有 ladder（并入方法族） | `engine.py` L1781 `_run_ts_nonconvergence_rescue`；L1791 `_run_int_nonconvergence_rescue`；L1801 `_run_stationary_point_nonconvergence_rescue` |
| 模式选择（R2 用） | `engine.py` L1309 `_select_target_ts_mode` |
| 局部 Hessian（O2 复用点） | `engine.py` L1363 `_run_local_hessian_analysis` |
| 位移工具（saddle_break 复用） | `engine.py` L3053 `_build_mode_displaced_xyz` |
| 触发判定（F1 分类源） | `engine.py` L1998 `_is_geometry_nonconvergence`；L2014 `_native_checkpoint_trigger` |
| Warmup 构建（W1 改参数处） | `engine.py` L829 `_maybe_warmup`；L861 warmup_spec |
| 主优化构建（W2 改参数处） | `engine.py` L914 `_run_primary_opt`；profile max_cycles |
| 配置解析 | `rph_core/steps/fidelity_profile.py` L421-425（rescue cfg 路径）、L558-681（字段映射）、L29-78（默认值） |
| 配置默认值 | `config/defaults.yaml` L74-93（warmup/max_cycles）、L191-234（rescue 区）、L273-293（S3 profile） |
| ORCA 渲染（TS_Mode/InHess/Recalc_Hess） | `rph_core/utils/qc_jobs.py` L69-108 `_orca_geom_block` |
| manifest 读写（兼容层） | `rph_core/steps/refinement/manifest_io.py` |
| UI 适配（增量扩展） | `rph_core/utils/ui_adapter.py` L211 `adapt_s3_structures`；`ui_state.py` |
| 事件残留（需清理） | `engine.py` L3426 `extra={"phase": "wave2_pass2_pass3"}` |
| 批次调度（O3 并发） | `rph_core/utils/stage_scheduler.py`（`max_concurrent_by_role` 钳制警告） |
