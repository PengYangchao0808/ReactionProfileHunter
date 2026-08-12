# RPH V4 S3/S4 统一双精度优化方案（v1.2 — S4 Warmup + INT Calc_Hess 修正版）

> **状态**: 审核修订完成，待确认后实施
> **日期**: 2026-07-20
> **修订历史**:
> - v1.0: 原始方案（已废弃）
> - v1.1: 19 项 P0 阻断性问题修正
> - **v1.2: S4 保留 Warmup + INT 正式优化使用 Calc_Hess true**

---

## 一、设计原则

### 1.1 核心思想

S3 和 S4 不是"两个共享逻辑的 stage"，而是**同一个 `RefinementEngine` 在同一套统一状态机下，用两套不同的 `FidelityProfile` 参数运行两次**。

```
S3 = RefinementEngine(profile=low_fidelity).run(requests, output_dir)
S4 = RefinementEngine(profile=high_fidelity).run(requests, output_dir)
```

### 1.2 硬性约束（审核修订后）

1. **唯一执行类**：`RefinementEngine` 是唯一实现；`LowLevelEngine` 和 `HighLevelEngine` 退化为空子类别名，不覆写任何方法
2. **统一 manifest schema**：`refinement_manifest_v1`，S3/S4 仅 `stage` 和 `fidelity` 字段不同
3. **核心引擎不出现 `if stage == "S3"/"S4"`**：所有差异来自 `FidelityProfile` 参数
4. **3-pass DAG 调度**：不依赖串行顺序和 `_ts_results` 字典
5. **不新建目录结构**：保持 `S3_LowLevel/` 和 `S4_HighLevel/` 输出目录，内部命名区分路径

### 1.3 理论合同（不变）

```
S3: ORCA B97-3c OPT/OptTS + independent Freq → ORCA r2SCAN-3c SP
S4: ORCA M062X/def2-SVP OPT/OptTS + independent Freq → ORCA wB97M-V/def2-TZVPP SP
Solvent: CPCM(acetone)
```

---

## 二、代码结构

### 2.1 文件清单

```
rph_core/
├── steps/
│   ├── refinement/
│   │   ├── __init__.py           # 导出 RefinementEngine, FidelityProfile
│   │   └── engine.py             # 🔴 唯一实现类 RefinementEngine
│   │
│   ├── step3_lowlevel/
│   │   ├── __init__.py           # LowLevelEngine = RefinementEngine (别名)
│   │   └── engine.py             # 保留或删除（空壳）
│   ├── step4_highlevel/
│   │   ├── __init__.py           # HighLevelEngine = RefinementEngine (别名)
│   │   └── engine.py             # 保留或删除（空壳）
│   │
│   ├── stage_calculator.py       # 🔴 核心修改：重构为 per-pass actions
│   └── fidelity_profile.py       # 🆕 FidelityProfile 数据类
│
├── utils/
│   ├── qc_jobs.py                # 🔴 新增 run_irc(), 修改 run_optimization() Hessian 参数
│   ├── qc_models.py              # 🔴 新增 IRCJobSpec, OptimizationSpec 类型化字段
│   ├── orca_interface.py         # 🔴 新增 IRC 接口, 修改 Hessian 渲染
│   └── ml_quality.py             # 🆕 属性级 ml_usability 计算
│
├── v4_orchestrator.py            # 🟡 S3/S4 统一为 RefinementEngine 调度
│
config/
└── defaults.yaml                 # 🟡 新增 refinement profile 块 + rescue 参数
```

### 2.2 `LowLevelEngine` / `HighLevelEngine` 生命周期

```python
# rph_core/steps/step3_lowlevel/__init__.py
from rph_core.steps.refinement.engine import RefinementEngine

class LowLevelEngine(RefinementEngine):
    """向后兼容别名 — 不覆写任何方法"""
    pass

# rph_core/steps/step4_highlevel/__init__.py
from rph_core.steps.refinement.engine import RefinementEngine

class HighLevelEngine(RefinementEngine):
    """向后兼容别名 — 不覆写任何方法"""
    pass
```

旧测试可继续使用 `LowLevelEngine` / `HighLevelEngine`，它们实际实例化 `RefinementEngine`。

---

## 三、统一 Manifest Schema

```
S3_LowLevel/manifest.json:
{
  "schema_version": "refinement_manifest_v1",
  "stage": "S3",
  "fidelity": "low",
  "profile_id": "b97_3c_r2scan_3c_v1",
  "structures": [...]
}

S4_HighLevel/manifest.json:
{
  "schema_version": "refinement_manifest_v1",
  "stage": "S4",
  "fidelity": "high",
  "profile_id": "m062x_wb97mv_v1",
  "structures": [...]
}
```

如需兼容旧代码读取旧 manifest，增加适配器（不修改旧文件）：

```python
def read_refinement_manifest(path: Path):
    data = json.loads(path.read_text())
    if data.get("schema_version") == "refinement_manifest_v1":
        return data
    if data.get("schema_version") == "s3_low_level_v3":
        return adapt_s3_to_refinement(data)
    if data.get("schema_version") == "s4_high_level_v4":
        return adapt_s4_to_refinement(data)
```

---

## 四、FidelityProfile 数据类

```python
@dataclass(frozen=True)
class FidelityProfile:
    """S3/S4 的全部差异参数"""

    stage: str                          # "S3" | "S4"
    fidelity: str                       # "low" | "high"
    profile_id: str                     # 标识符，进入 checkpoint 签名

    # 几何优化
    geometry_method: str                # "B97-3c" | "M062X"
    geometry_basis: str                 # "" | "def2-SVP"
    geometry_aux_basis: str             # "" | "def2/J"
    geometry_grid: str | None           # None | "DefGrid3"
    geometry_scf: str | None            # None | "TightSCF"
    route_minimum: str                  # "Opt"
    route_ts: str                       # "OptTS"
    max_cycles_minimum: int             # 100 | 200
    max_cycles_ts: int                  # 150 | 200

    # 初始 Hessian（按角色）
    initial_hessian_precursor: str      # "model"
    initial_hessian_product: str        # "model"
    initial_hessian_intermediate: str   # "calculate"  ← INT 也用 Calc_Hess
    initial_hessian_ts: str             # "calculate"

    # Warmup
    warmup_max_cycles_int: int          # S3=8, S4=4
    warmup_max_cycles_ts: int           # S3=12, S4=6

    # 频率
    frequency_method: str
    frequency_basis: str

    # 单点能
    sp_method: str                      # "r2SCAN-3c" | "wB97M-V"
    sp_basis: str                       # "" | "def2-TZVPP"
    sp_aux_basis: str                   # "" | "def2/J"

    # 溶剂
    solvent: str                        # "acetone"
    solvent_model: str                  # "CPCM"

    # 救援策略
    ts_rescue_enabled: bool
    ts_rescue_recalc_hessian_interval: int  # 5
    ts_rescue_trust_radius: float | None    # 0.3

    int_rescue_enabled: bool
    int_rescue_mode_displacement: bool

    irc_enabled: bool
    irc_max_iter: int                   # 50
    irc_direction: str                  # "both"

    # 资源
    cores_per_worker: int
    memory_gb_per_worker: int
    max_workers: int
    timeout_seconds: int

    # 热化学
    temperature_k: float                # 247.55
    standard_state: str                 # "1M"
    qrrho: bool                         # True
    ensemble_correction_source: str     # "s1"
```

渲染 ORCA 输入时，按角色选择初始 Hessian 策略：

```python
def render_orca_geom_block(spec: OptimizationSpec, role: str, profile: FidelityProfile) -> str:
    lines = []
    # 按角色确定初始 Hessian
    hessian_mode = profile.initial_hessian.get(role, "model")
    if hessian_mode == "calculate":
        lines.append("  Calc_Hess true")
    elif hessian_mode == "read":
        lines.append("  InHess Read")
        lines.append(f'  InHessName "{spec.hessian_filename}"')
    if spec.ts_mode is not None:
        lines.append(f"  TS_Mode {{M {spec.ts_mode}}}")
    if spec.recalc_hessian is not None:
        lines.append(f"  Recalc_Hess {spec.recalc_hessian}")
    if spec.trust is not None:
        lines.append(f"  Trust {spec.trust}")
    return "\n".join(lines)
```

**注意：Calc_Hess 不能阻止 INT 滑入 product。** 它能改善初始曲率描述和数值稳定性，但不能在不存在独立 INT 的势能面上创造 INT。如果该理论水平下 INT 本就不是独立最低点，即使 Warmup + Calc_Hess + TightOpt，最终仍可能收敛到 product。此时结果应标记为 `converged_to_product_minimum`（有效化学结论），而非 `optimization_failed`。

---

## 五、S3 / S4 启动差异

### 5.1 S3 启动

```
S1 selected.xyz (precursor, products)
S2 ts_guess.xyz + intermediate_xyz
       │
       ▼
    build StructureRequest[]
       │
       │ seed_state:
       │ · precursor / product  → stable_minimum_seed
       │ · intermediate / TS    → path_seed (带 forming_bonds)
       │
       ▼
    RefinementEngine(profile=low_fidelity)
       │
       ▼
    S3_LowLevel/manifest.json
```

### 5.2 S4 启动

```
S3_LowLevel/manifest.json
       │
       │ 对每个 S3 structure:
       │
       │ ┌─ S3 opt_status=complete AND identity_role_match:
       │ │    input   = S3 opt_xyz
       │ │    seed    = optimized_parent
       │ │    warmup  = true           ← S4 也执行 warmup
       │ │    constraints = S3 几何中 forming_bonds 的当前距离
       │ │    max_cycles = 4 (INT) / 6 (TS) — 比 S3 短
       │ │
       │ │    目的: ∇E_S3(x_S3) ≈ 0，但 ∇E_S4(x_S3) ≠ 0
       │ │          在保持 S3 已获得反应进度的前提下，
       │ │          让其他坐标适应 S4 势能面
       │ │
       │ ├─ S3 opt_status=complete AND identity_role_mismatch:
       │ │    input   = S2 original_seed_xyz
       │ │    seed    = parent_identity_mismatch
       │ │    warmup  = true (S2 seed 距离约束)
       │ │
       │ └─ S3 opt_status!=complete:
       │      input   = S2 original_seed_xyz
       │      seed    = parent_opt_failed
       │      warmup  = true (S2 seed 距离约束)
       │
       ▼
    RefinementEngine(profile=high_fidelity)
       │
       ▼
    S4_HighLevel/manifest.json
```

**S4 为何保留 Warmup**：S3 优化后结构仅在 S3 势能面上接近驻点（∇E_S3 ≈ 0），但在 S4 势能面上不一定（∇E_S4 ≠ 0）。尤其对于浅势能面、高度异步 INT、Lewis acid 配位结构、S3/S4 对色散描述不同的体系，直接无约束高层级优化可能滑入相邻盆地。ORCA 约束优化将梯度投影到允许自由度子空间，先放松其余坐标同时暂时保持关键反应坐标。

---

## 六、3-Pass DAG 调度

**不再依赖串行 `_ts_results` 字典。** 统一为三遍调度，S3/S4 完全相同：

```
═══════════════════════════════════════════════════════════
  PASS 0: PREFLIGHT (独立，所有结构并行)
═══════════════════════════════════════════════════════════

  每个结构独立执行:
  · validate XYZ / atom_count / atom_mapping
  · resolve charge / multiplicity
  · write provenance.json
  · create AttemptRecorder

  失败: hard failure → 标记 failed_preflight，不进入后续

═══════════════════════════════════════════════════════════
  PASS 1: PRIMARY CALCULATION (独立，所有结构并行)
═══════════════════════════════════════════════════════════

  ┌─ WARMUP (仅 intermediate / TS，S3/S4 都执行)
  │    · 约束距离 S3: S2 seed 距离 / S4: S3 几何中的当前距离
  │    · LooseOpt 收敛
  │    · S3 max_cycles: INT=8, TS=12
  │    · S4 max_cycles: INT=4, TS=6
  │    · 失败不阻断
  │
  ├─ PRIMARY OPT
  │    · precursor/product  → Opt + model Hessian
  │    · intermediate        → Opt + Calc_Hess true
  │    · TS                  → OptTS + Calc_Hess true
  │    · OPT 成功 → 写入 opt.xyz
  │    · OPT 失败 → 记录 last_xyz，不阻断
  │
  │    ⚠ 正确顺序: WARMUP → 在 warmup 输出几何上 Calc_Hess → OPT
  │      不能在 warmup 之前的几何上算 Hessian
  │
  ├─ FINAL-GEOMETRY FREQ (OPT 成功才执行)
  │    · 独立 Freq/AnFreq/NumFreq 在 opt.xyz 上
  │    · 不依赖初始 Hessian（独立计算）
  │
  └─ INITIAL CLASSIFICATION
       · TS: hessian_index / curvature_class / mode_identity
       · INT: frequency check / forming_bond distances
       · ALL ROLES (含 product/precursor): identity 验证

═══════════════════════════════════════════════════════════
  GLOBAL BARRIER — 等待所有 Pass 1 结果
═══════════════════════════════════════════════════════════

═══════════════════════════════════════════════════════════
  PASS 2: RESCUE DAG (依赖型，仅需要时执行)
═══════════════════════════════════════════════════════════

  DAG 依赖:
    variant_TS_validated ──→ variant_INT_endpoint_discovery

  ┌─ TS RESCUE (stationary_point_class ∈ {soft, wrong_mode,
  │            higher_order_saddle})
  │    Level 1: 读取 Prelim Freq Hessian + TS_Mode {M} 指定模式
  │             → OptTS + InHess Read + TS_Mode
  │             → 独立 Freq
  │    Level 2: Recalc_Hess N (默认 N=5)
  │             → OptTS + Calc_Hess true + Recalc_Hess N
  │             → 独立 Freq
  │
  ├─ INT RESCUE (有虚频 — 优先自身 mode displacement)
  │    Level 1: TightOpt + 更严格网格 + 重新 Freq
  │    Level 2: 沿 INT 自身虚模 ± displacement → Opt → Freq
  │
  ├─ MINIMUM RESCUE (任何 minimum 有显著虚频)
  │    Level 1: TightOpt + 更严格网格 + 重新 Freq
  │    Level 2: 沿最低虚模 ± displacement → Opt → Freq
  │
  └─ OPTIONAL: TS ENDPOINT DISCOVERY (通过 IRC)
       · 仅当 TS 已确认正确 (validated_ts)
       · IRC Direction both → 产生 endpoint_A 和 endpoint_B
       · 对两个 endpoint 分别 OPT+Freq
       · 分类: precursor / product / intermediate / other
       · 不预设"一定能找到 INT"

═══════════════════════════════════════════════════════════
  CANONICAL ATTEMPT SELECTION
═══════════════════════════════════════════════════════════

  从所有 attempt（primary + rescue）中选 canonical geometry:
  · 候选排序: (mode_identity, hessian_index, converged,
                alignment_score, -gradient_norm)
  · 确定 resolved_kind / resolved_identity / geometry_hash

═══════════════════════════════════════════════════════════
  PASS 3: CANONICAL RESULT
═══════════════════════════════════════════════════════════

  ┌─ FINAL FREQ (仅当 canonical geometry ≠ Pass 1 Freq 几何)
  │    · 独立 Freq 在 canonical.xyz 上
  │
  ├─ FINAL SP (在 canonical.xyz 上)
  │    · OPT 成功 → stationary_point_energy = true
  │    · OPT 失败 → diagnostic_sp, stationary_point_energy = false
  │
  ├─ THERMOCHEMISTRY
  │    · G_composite = E_SP + (G_freq - E_freq) + ΔG_ensemble
  │    · H_composite = E_SP + (H_freq - E_freq) + ΔH_ensemble
  │    · 记录 geometry_level / frequency_level / sp_level
  │
  └─ PROPERTY-LEVEL ML USABILITY
       · 拆分为属性级 dict，不单一 bool
```

---

## 七、核心计算流程（per-pass 详情）

### 7.1 阶段 A：WARMUP

**S3/S4 使用完全相同的 Warmup 逻辑，根据化学角色判断，不根据 seed_state。**

| 结构 | S3 | S4 |
|------|-----|-----|
| precursor | 不 warmup | 不 warmup |
| product | 不 warmup | 不 warmup |
| intermediate | warmup | warmup |
| TS | warmup | warmup |

```python
def warmup_applies(structure: StructureRequest) -> bool:
    """Warmup 基于化学角色，S3/S4 完全相同"""
    return (
        structure.role in {"intermediate", "ts"}
        and bool(structure.forming_bonds)
    )

def build_warmup_constraints(request, input_xyz, profile):
    """
    S3: 约束距离 = S2 seed 几何中的 forming_bond 距离
    S4: 约束距离 = S3 输入几何中的当前距离 (保持 S3 已获得的反应进度)
    """
    constraints = []
    for (i, j) in request.forming_bonds:
        if request.source_stage == "S3":
            # S4: 使用 S3 几何中的当前 forming_bond 距离
            dist = measure_distance(input_xyz, i, j)
        else:
            # S3: 使用 S2 seed 距离
            dist = measure_distance(request.original_seed_xyz, i, j)
        constraints.append(BondConstraint(i, j, dist))
    return constraints

def run_warmup(request, input_xyz, output_dir, profile):
    constraints = build_warmup_constraints(request, input_xyz, profile)
    result = run_optimization(
        spec=OptimizationSpec(
            task="minimum",
            route="Opt",
            route_extras="LooseOpt",
            bond_constraints=constraints,
            allow_unconverged_geometry=True,
            max_cycles=profile.warmup_max_cycles.get(
                request.role,
                12 if request.role == "ts" else 8,
            ),
        ),
        input_xyz=input_xyz,
        output_dir=output_dir / "warmup",
    )
    return result.output_xyz if result.status in ("complete", "partial") else input_xyz
```

**关键差异（纯参数，非代码逻辑）**：

| 参数 | S3 | S4 | 原因 |
|------|-----|-----|------|
| 约束距离 | S2 seed 距离 | S3 opt 几何中的当前距离 | S4 已在 S3 势能面接近驻点，不需要拉回 S2 坐标 |
| max_cycles (INT) | 8 | 4 | S4 从更接近驻点的几何出发 |
| max_cycles (TS) | 12 | 6 | 同上 |

### 7.2 阶段 B：PRIMARY OPT + FREQ

**Warmup 与 Calc_Hess 解决不同问题：**

- **Warmup**：在约束条件下让结构靠近目标势能面，减少后续无约束优化中的大位移和盆地滑落风险
- **Calc_Hess**：为优化器提供准确的初始曲率信息，改善步长方向和收敛稳定性

**正确顺序**：WARMUP → 在 warmup 输出几何上 Calc_Hess → 正式 OPT

不能在 warmup 之前在原始几何上算 Hessian，因为 warmup 已改变了几何。

#### 7.2.1 各角色优化策略

| 角色 | 优化类型 | 初始 Hessian | ORCA 输入 |
|------|---------|-------------|-----------|
| precursor | Opt | model (默认) | `! Opt B97-3c` |
| product | Opt | model (默认) | `! Opt B97-3c` |
| intermediate | Opt | **Calc_Hess true** | `! Opt B97-3c` + `%geom Calc_Hess true` |
| TS | OptTS | **Calc_Hess true** | `! OptTS B97-3c` + `%geom Calc_Hess true` |

**为何 INT 需要 Calc_Hess**：INT 不是普通稳定 minimum——

- 来自反应路径谷点，可能靠近 TS
- 曲率很浅，形成键高度异步
- 容易滑入 product
- 可能包含较软的 Lewis acid 配位自由度
- 模型 Hessian 可能无法准确描述初始局部曲率

#### 7.2.2 INT 优化（正确 ORCA 输入）

```
S3:
! B97-3c Opt TightSCF

%geom
  Calc_Hess true
  MaxIter 100
end

S4:
! M062X def2-SVP Opt TightSCF

%geom
  Calc_Hess true
  MaxIter 200
end
```

```python
int_opt_spec = OptimizationSpec(
    task="minimum",
    route="Opt",
    method=profile.geometry_method,
    basis=profile.geometry_basis,
    initial_hessian="calculate",      # Calc_Hess true
    max_cycles=profile.max_cycles_minimum,
)
```

#### 7.2.3 TS 优化（正确 ORCA 输入）

```
S3:
! B97-3c OptTS TightSCF

%geom
  Calc_Hess true
  MaxIter 150
end

S4:
! M062X def2-SVP OptTS TightSCF

%geom
  Calc_Hess true
  MaxIter 200
end
```

```python
# TS 优化
ts_opt_spec = OptimizationSpec(
    task="ts",
    route="OptTS",
    method=profile.geometry_method,
    basis=profile.geometry_basis,
    initial_hessian="calculate",      # Calc_Hess true
    max_cycles=profile.max_cycles_ts,
)
```

#### 7.2.4 频率（独立运行，不继承初始 Hessian）

```python
# ❌ 错误: 不存在的 UseHess
# ! Freq UseHess B97-3c

# ✅ 正确: 独立 Freq
# ! Freq B97-3c
# 或数值: ! NumFreq B97-3c

freq_result = run_frequency(
    spec=QCJobSpec(
        task="freq",            # 或 "numfreq"
        method=profile.frequency_method,
        basis=profile.frequency_basis,
    ),
    input_xyz=opt_result.output_xyz,    # 独立在优化后几何上运行
    output_dir=output_dir / "freq",
)
```

**关键规则**：初始 `Calc_Hess true` 计算的是优化初期初猜结构的 Hessian。优化完成后几何已变，不能将初始 Hessian 当作最终频率。优化过程中即便发生过 `Recalc_Hess`，最后一次 Hessian 也不保证恰好对应最终收敛几何。**频率必须独立在最终几何上运行。**

### 7.3 阶段 C：RESCUE

#### 7.3.1 TS 救援梯度

```
Level 1: 读取 Hessian + 指定目标模式

  OptTS
  %geom
    InHess Read
    InHessName "prelim_freq.hess"
    TS_Mode {M 0}
  end

Level 2: 周期性重算 Hessian

  OptTS
  %geom
    Calc_Hess true
    Recalc_Hess 5
    Trust 0.3
  end
```

配置对应：

```yaml
ts_rescue:
  enabled: true
  strategies:
    - read_hessian_target_mode    # Level 1
    - recalc_hessian              # Level 2
  recalc_hessian_interval: 5
  trust_radius_update: true
```

#### 7.3.2 INT 救援梯度

```
情况 A: OPT 成功，有小虚频
  → TightOpt + 更严格积分网格 + 重新 SCF + 重新 Freq
  → 若仍未解决: Recalc_Hess 10

情况 B: OPT 成功，有一个显著虚频
  → 沿 INT 自身虚模 ± displacement → Opt + Calc_Hess true → Freq
  → 选择 zero-imag 候选
  → 若仍未解决: Recalc_Hess 5

情况 C: OPT 滑入已知 product/precursor
  → 非数值失败，是 identity mismatch
  → 结果标记: converged_to_product_minimum（有效化学结论）
  → 可选: 通过已验证的 TS 执行 IRC endpoint discovery

情况 D: OPT 未收敛
  → 继续优化 / 收紧步长 / 更换坐标系
  → 无有效 TS 时 IRC 无法执行
```

**INT Primary 使用 Calc_Hess true，Rescue 使用 Recalc_Hess**：

| 阶段 | Hessian 策略 | 适用场景 |
|------|-------------|---------|
| Primary | `Calc_Hess true` | 所有 INT 的初始优化 |
| Rescue L1 | `Recalc_Hess 10` | 收敛慢、振荡或 Hessian 更新失效 |
| Rescue L2 | `Recalc_Hess 5` | 特别不稳定结构 |

不建议对所有 INT 默认 `Recalc_Hess 1`（每步重算 Hessian），成本极高，且如果势能面本不存在独立 INT，频繁重算也无济于事。

#### 7.3.3 TS Endpoint Discovery（IRC，不再是"INT rescue"）

```python
def run_ts_endpoint_discovery(validated_ts, profile, output_dir):
    """IRC 的作用: 判定 TS 连接哪些最低点"""
    irc_result = run_irc(
        spec=IRCJobSpec(
            direction="both",               # 同时跑 forward 和 backward
            max_iter=profile.irc_max_iter,  # 50
            init_hessian="read",            # 读取已验证 TS 的 Hessian
            hessian_filename=validated_ts.hess_file,
            hessian_mode=0,
            init_displacement_mode="energy",
            initial_delta_energy_mEh=2.0,
            scale_initial_displacement=0.1,
            scale_steepest_descent=0.15,
            adaptive_step=True,
        ),
        input_xyz=validated_ts.opt_xyz,
        output_dir=output_dir / "irc_endpoints",
    )

    if irc_result.status != "complete":
        return {"irc_status": "failed"}

    # 获取两侧终点
    endpoint_a = irc_result.endpoint_a
    endpoint_b = irc_result.endpoint_b

    # 分别执行 OPT+Freq
    opt_a = run_optimization(OptSpec(task="minimum"), endpoint_a.xyz, ...)
    opt_b = run_optimization(OptSpec(task="minimum"), endpoint_b.xyz, ...)

    # 分类
    identity_a = classify_identity(opt_a, ...)
    identity_b = classify_identity(opt_b, ...)

    return {
        "irc_status": "complete",
        "endpoint_a_identity": identity_a,   # "precursor" | "product" | "intermediate" | "other"
        "endpoint_b_identity": identity_b,
        "distinct_intermediate_found": "intermediate" in {identity_a, identity_b},
        "mechanism_implication": (
            "concerted_at_this_fidelity" if not distinct_intermediate_found
            else "stepwise_ts_int_product"
        ),
    }
```

**关键语义变化**：不预设"从 TS 做 IRC 一定能找回 INT"。IRC 可能证明该 TS 直接连接 precursor 和 product，不存在独立 INT。

正确的 ORCA IRC 输入：

```
! IRC B97-3c

%irc
  MaxIter 50
  Direction both

  InitHess read
  Hess_Filename "validated_ts.hess"
  HessMode 0

  Init_Displ_DE
  DE_Init_Displ 2.0
  Scale_Init_Displ 0.1

  Scale_Displ_SD 0.15
  Adapt_Scale_Displ true
end
```

注意参数名：`MaxIter`（不是 MaxPoints），`Direction both`（不是 forward/reverse）。

### 7.4 阶段 D：FINAL SP + 热化学 + ML 可用性

```python
def run_final_sp_and_thermo(canonical, profile, output_dir):
    # 1. FINAL SP
    sp_result = run_single_point(
        spec=QCJobSpec(task="sp", method=profile.sp_method, basis=profile.sp_basis),
        input_xyz=canonical.xyz,
        output_dir=output_dir / "sp",
    )

    # 2. THERMOCHEMISTRY
    # G_composite = E_SP(canonical) + [G_freq(canonical) - E_freq(canonical)] + ΔG_ensemble
    if canonical.freq_available and sp_result.complete:
        composite = {
            "gibbs_free_energy_hartree":
                sp_result.energy_hartree
                + canonical.G_correction_hartree
                + canonical.ensemble_correction_hartree,
            "enthalpy_hartree":
                sp_result.energy_hartree
                + canonical.H_correction_hartree
                + canonical.ensemble_enthalpy_correction_hartree,
            "geometry_level": profile.geometry_method,
            "frequency_level": profile.frequency_method,
            "single_point_level": profile.sp_method,
            "temperature_K": profile.temperature_k,
            "standard_state": profile.standard_state,
            "qrrho": profile.qrrho,
            "ensemble_correction_id": canonical.ensemble_correction_id,
            "geometry_consistent": canonical.freq_geometry_hash == canonical.sp_geometry_hash,
        }

    # 3. ML USABILITY (属性级，非单一 bool)
    ml_usability = {
        "geometry": canonical.opt_converged,
        "electronic_energy": sp_result.complete,
        "enthalpy": sp_result.complete and canonical.freq_available,
        "gibbs_free_energy": sp_result.complete and canonical.freq_available,
        "frequency_descriptors": canonical.freq_available,
        "ts_descriptors": (
            canonical.requested_kind == "ts"
            and canonical.ts_frequency_valid
        ),
        "intermediate_descriptors": (
            canonical.requested_role == "intermediate"
            and canonical.identity_status == "role_matched"
        ),
        "mechanism_label": True,           # 即使 INT 未找到，分类仍有价值
        "multifidelity_pair": False,       # 由 S3/S4 pairing 阶段填写
    }
```

---

## 八、结构检测规则

### 8.1 TS 分类

```python
def classify_ts(freq_result, forming_bonds):
    imag_freqs = [f for f in freq_result.frequencies_cm1 if f < -10.0]
    imag_count = len(imag_freqs)

    # 1. Hessian 阶数
    hessian_index = imag_count

    # 2. 曲率分类
    if imag_count == 1:
        if imag_freqs[0] < -50.0:
            curvature_class = "strict"
        else:
            curvature_class = "soft"
    elif imag_count == 0:
        curvature_class = "none"
    else:
        curvature_class = "multi_imaginary"

    # 3. 模式身份
    mode_identity = check_mode_alignment(
        freq_result.normal_modes,
        imag_freqs,
        forming_bonds,
    ) if imag_count == 1 else "unavailable"

    # 4. 驻点分类
    if hessian_index == 1 and curvature_class == "strict" and mode_identity == "target":
        sp_class = "valid_target_ts"
    elif hessian_index == 1 and curvature_class == "soft" and mode_identity == "target":
        sp_class = "soft_target_ts"
    elif hessian_index == 1 and mode_identity == "unrelated":
        sp_class = "first_order_wrong_mode"
    elif hessian_index >= 2:
        sp_class = "higher_order_saddle"
    elif hessian_index == 0:
        sp_class = "minimum_after_optts"
    else:
        sp_class = "unclassifiable"

    return {
        "hessian_index": hessian_index,
        "curvature_class": curvature_class,
        "mode_identity": mode_identity,
        "stationary_point_class": sp_class,
    }
```

**保留 soft-mode review**：`Calc_Hess true` 帮助优化的初猜找到负曲率方向，不能替代最终驻点的 soft-mode 判定。两者回答不同问题。

### 8.2 INT 检测（废除固定 1.7 Å）

```python
def classify_int(opt_result, freq_result, forming_bonds,
                 precursor_ref, product_ref, atom_mapping):
    """联合多判据，不依赖单一固定距离阈值"""

    if opt_result.status != "complete":
        return {"identity": "opt_failed"}

    # 1. 频率检查
    sig_imag = [f for f in freq_result.frequencies_cm1 if f < -10.0]
    if sig_imag:
        return {"identity": "imaginary_frequency", "imag_count": len(sig_imag)}

    # 2. 归一化反应进度
    progress = {}
    for bond in forming_bonds:
        d_current = measure_bond_distance(opt_result.opt_xyz, bond)
        d_precursor = measure_bond_distance(precursor_ref, bond)
        d_product = measure_bond_distance(product_ref, bond)
        if d_precursor != d_product:
            progress[bond] = (d_precursor - d_current) / (d_precursor - d_product)
    avg_progress = sum(progress.values()) / len(progress) if progress else 0.5

    # 3. 分子图拓扑
    mol_graph = build_molecular_graph(opt_result.opt_xyz, forming_bonds)

    # 4. Mapped RMSD
    rmsd_to_product = compute_mapped_rmsd(opt_result.opt_xyz, product_ref, atom_mapping)
    rmsd_to_precursor = compute_mapped_rmsd(opt_result.opt_xyz, precursor_ref, atom_mapping)

    # 5. 综合判定
    if mol_graph.matches(product_ref) and rmsd_to_product < 0.3:
        return {"identity": "collapsed_to_product"}
    if mol_graph.matches(precursor_ref) and rmsd_to_precursor < 0.3:
        return {"identity": "collapsed_to_precursor"}
    if 0.2 < avg_progress < 0.8:
        return {"identity": "distinct_intermediate", "avg_progress": avg_progress}

    return {"identity": "topology_ambiguous", "avg_progress": avg_progress}
```

注意：**所有 role（含 product/precursor）都必须经过 identity 验证**。product 和 precursor 也可能解离 Lewis acid、改变配位位点、发生异构化。

### 8.3 ORCA 关键词规范（重要）

| 意图 | ❌ 错误写法 | ✅ 正确写法 |
|------|-----------|-----------|
| TS 优化初算 Hessian | `CalcFC` | `%geom Calc_Hess true` |
| 周期性重算 Hessian | `Calcall` | `%geom Recalc_Hess N` |
| 复用 Hessian 做频率 | `UseHess` | 独立运行 `! Freq` |
| IRC 最大步数 | `MaxPoints` | `MaxIter` |
| IRC 方向 | `forward` / `reverse` | `forward` / `backward` / `both` |
| IRC 步长 | `StepSize` | `Scale_Init_Displ` + `Scale_Displ_SD` |

---

## 九、完整配置文件

```yaml
# config/defaults.yaml

refinement:
  engine_schema: refinement_engine_v1

  # ===== 通用策略（S3/S4 共享，不可被覆盖） =====
  common:
    workflow:
      warmup:
        enabled_roles:                  # 按化学角色，S3/S4 完全相同
          - intermediate
          - ts
        constraint_mode: freeze_forming_bonds_at_input_distance
        fallback_to_input_geometry: true
        accept_partial_geometry: true

      initial_hessian:                  # 按角色，S3/S4 完全相同
        precursor: model
        product: model
        intermediate: calculate         # ← INT 也用 Calc_Hess true
        ts: calculate

      frequency_always_independent: true
      continue_on_structure_failure: true

    rescue:
      ts_rescue:
        enabled: true
        strategies:
          - read_hessian_target_mode      # Level 1
          - recalc_hessian                 # Level 2
        recalc_hessian_interval: 5
        trust_radius_update: true

      int_rescue:
        enabled: true
        mode_displacement_first: true      # 有虚频优先自身 mode displacement

      minimum_rescue:
        enabled: true

    irc:
      direction: both
      max_iter: 50
      init_hessian: read
      init_displacement_mode: energy
      initial_delta_energy_mEh: 2.0
      scale_initial_displacement: 0.1
      scale_steepest_descent: 0.15
      adaptive_step: true

    identity:
      use_topology: true
      use_mapped_rmsd: true
      use_reaction_progress: true
      validate_all_roles: true           # product/precursor 也验证

    soft_mode:
      enabled: true                       # 保留，不与 Calc_Hess 冲突
      imaginary_cutoff_cm1: -50.0
      soft_mode_window_cm1: [-50.0, -10.0]
      mode_alignment_threshold: 0.30

    thermochemistry:
      temperature_K: 247.55
      standard_state: "1M"
      qrrho: true
      ensemble_correction_source: s1
      composite_formula: "E_SP + (G_freq - E_freq) + ΔG_ensemble"

    ml_quality:
      granular: true                     # 属性级，非单一 bool
      geometry_consistent_check: true

    output:
      manifest_schema: refinement_manifest_v1

  # ===== S3 Profile (low fidelity) =====
  s3:
    fidelity: low
    profile_id: b97_3c_r2scan_3c_v1

    warmup:
      max_cycles:
        intermediate: 8                # S2 seed 较粗糙，需更多循环
        ts: 12

    geometry:
      method: B97-3c
      basis: ""
      route_minimum: Opt
      route_ts: OptTS
      max_cycles_minimum: 100
      max_cycles_ts: 150
      timeout: 864000

    frequency:
      method: B97-3c
      basis: ""

    single_point:
      method: r2SCAN-3c
      basis: ""
      timeout: 864000

    solvent:
      solvent: acetone
      solvent_model: CPCM

    resources:
      cores_per_worker: 16
      memory_gb_per_worker: 32
      max_workers: 1

  # ===== S4 Profile (high fidelity) =====
  s4:
    fidelity: high
    profile_id: m062x_wb97mv_v1

    warmup:
      max_cycles:
        intermediate: 4                # S3 几何更接近驻点，可更短
        ts: 6

    geometry:
      method: M062X
      basis: def2-SVP
      aux_basis: def2/J
      grid: DefGrid3
      scf: TightSCF
      route_minimum: Opt
      route_ts: OptTS
      max_cycles_minimum: 200
      max_cycles_ts: 200
      timeout: 864000

    frequency:
      method: M062X
      basis: def2-SVP

    single_point:
      method: wB97M-V
      basis: def2-TZVPP
      aux_basis: def2/J
      timeout: 864000

    solvent:
      solvent: acetone
      solvent_model: CPCM

    resources:
      cores_per_worker: 32
      memory_gb_per_worker: 64
      max_workers: 1

# ===== 向后兼容 =====
theory:
  s3_low_level:    # 别名，实际引用 refinement.s3
    optimization: { engine: orca }
    single_point:  { engine: orca }
  s4_high_precision:
    optimization: { engine: orca }
    single_point:  { engine: orca }
```

---

## 十、关键代码修改点（修订后）

### 10.1 `rph_core/steps/refinement/engine.py`（🆕 新建）

```python
class RefinementEngine:
    """S3/S4 唯一实现类。LowLevelEngine 和 HighLevelEngine 是空子类别名。"""

    def __init__(self, config: dict, profile: FidelityProfile,
                 parent_manifest_paths: dict = None):
        self.config = config
        self.profile = profile
        self.parent_manifest_paths = parent_manifest_paths or {}
        self._progress_reporter = None

    def run(self, requests: List[StructureRequest],
            output_dir: Path) -> Path:
        """
        三遍调度:
        Pass 0: Preflight (并行)
        Pass 1: Primary Calc (并行) → GLOBAL BARRIER
        Pass 2: Rescue DAG (依赖型)
        Pass 3: Canonical Result (独立)
        """
        # ... 实现 ...
```

### 10.2 `rph_core/utils/qc_models.py`（类型化字段）

```python
@dataclass(frozen=True)
class OptimizationSpec:
    task: str                           # "opt" | "opt_ts"
    method: str
    basis: str = ""
    route: str = "Opt"
    route_extras: str = ""
    initial_hessian: str = "model"      # "model" | "calculate" | "read" | "hybrid"
    hessian_filename: str | None = None
    ts_mode: int | None = None          # 用于 TS_Mode {M N}
    recalc_hessian: int | None = None   # Recalc_Hess N
    trust: float | None = None          # Trust radius
    # ... 其余字段 ...

@dataclass(frozen=True)
class IRCJobSpec:
    direction: str = "both"             # "both" | "forward" | "backward"
    max_iter: int = 50
    init_hessian: str = "read"          # "read" | "calc_anfreq" | "calc_numfreq"
    hessian_filename: str | None = None
    hessian_mode: int = 0
    init_displacement_mode: str = "energy"  # "energy" | "length"
    initial_delta_energy_mEh: float = 2.0
    scale_initial_displacement: float = 0.1
    scale_steepest_descent: float = 0.15
    adaptive_step: bool = True
    method: str = ""
    basis: str = ""
```

### 10.3 `rph_core/utils/orca_interface.py`（Hessian 渲染 + IRC）

```python
class ORCAInterface:

    def _render_geom_block(self, spec: OptimizationSpec) -> str:
        lines = []
        lines.append("%geom")
        if spec.initial_hessian == "calculate":
            lines.append("  Calc_Hess true")
        elif spec.initial_hessian == "read":
            lines.append("  InHess Read")
            lines.append(f'  InHessName "{spec.hessian_filename}"')
        if spec.ts_mode is not None:
            lines.append(f"  TS_Mode {{M {spec.ts_mode}}}")
        if spec.recalc_hessian is not None:
            lines.append(f"  Recalc_Hess {spec.recalc_hessian}")
        if spec.trust is not None:
            lines.append(f"  Trust {spec.trust}")
        lines.append("end")
        return "\n".join(lines)

    def run_irc(self, spec: IRCJobSpec, input_xyz: Path,
                output_dir: Path) -> dict:
        """执行 ORCA IRC 计算"""
        input_lines = [
            f"! IRC {spec.method} {spec.basis}",
            "",
            "%irc",
            f"  MaxIter {spec.max_iter}",
            f"  Direction {spec.direction}",
            f"  InitHess {spec.init_hessian}",
        ]
        if spec.hessian_filename:
            input_lines.append(f'  Hess_Filename "{spec.hessian_filename}"')
        input_lines.append(f"  HessMode {spec.hessian_mode}")
        input_lines.append("")
        if spec.init_displacement_mode == "energy":
            input_lines.append("  Init_Displ_DE")
            input_lines.append(f"  DE_Init_Displ {spec.initial_delta_energy_mEh}")
        input_lines.append(f"  Scale_Init_Displ {spec.scale_initial_displacement}")
        input_lines.append(f"  Scale_Displ_SD {spec.scale_steepest_descent}")
        if spec.adaptive_step:
            input_lines.append("  Adapt_Scale_Displ true")
        input_lines.append("end")
        # ... 执行和解析 ...
```

### 10.4 `rph_core/utils/qc_jobs.py`

```python
def run_irc(spec: IRCJobSpec, input_xyz: Path, output_dir: Path,
            config: dict, subprocess_callback=None) -> QCJobResult:
    """执行 IRC 计算"""
    # ...

def run_optimization(spec: OptimizationSpec, ...):
    # 移除: CalcFC / Calcall 字符串处理
    # 改为: 读取 spec.initial_hessian / spec.recalc_hessian 渲染 ORCA 输入
    # ...
```

### 10.5 `rph_core/v4_orchestrator.py`

```python
# 统一入口
from rph_core.steps.refinement import RefinementEngine
from rph_core.steps.fidelity_profile import FidelityProfile

def _run_refinement_stage(self, stage: str, structures, output_dir):
    profile = FidelityProfile.from_config(self.config, stage)
    engine = RefinementEngine(self.config, profile,
                              parent_manifest_paths=self._stage_manifest_paths())
    return engine.run(structures, output_dir)

# S3
s3_manifest = self._run_refinement_stage("S3", s3_structures,
                                          self.work_dir / "S3_LowLevel")

# S4
s4_manifest = self._run_refinement_stage("S4", s4_structures,
                                          self.work_dir / "S4_HighLevel")
```

---

## 十一、实施计划

### 11.1 P0 — 阻断性修复（必须先完成）

| # | 修改 | 文件 |
|---|------|------|
| 1 | 创建 `RefinementEngine` 唯一实现类 | `rph_core/steps/refinement/engine.py` |
| 2 | LowLevelEngine/HighLevelEngine → 空子类别名 | `step3_lowlevel/`, `step4_highlevel/` |
| 3 | 统一 manifest schema → `refinement_manifest_v1` | `engine.py` |
| 4 | `Calc_Hess true` 替代 `CalcFC` | `qc_models.py`, `orca_interface.py` |
| 5 | `Recalc_Hess N` 替代 `Calcall` | `qc_models.py`, `orca_interface.py`, `stage_calculator.py` |
| 6 | 删除 `UseHess` 逻辑，Freq 独立运行 | `stage_calculator.py`, `qc_jobs.py` |
| 7 | IRC `MaxIter` / `Direction both` / 正式位移参数 | `qc_models.py`, `orca_interface.py` |
| 8 | INT rescue → TS endpoint discovery | `stage_calculator.py` |
| 9 | TS rescue 梯度: read Hessian + TS_Mode → Recalc_Hess | `stage_calculator.py` |
| 10 | INT 救援梯度: mode displacement → IRC discovery | `stage_calculator.py` |
| 11 | 保留 soft-mode review | `stage_calculator.py` |
| 12 | 废除固定 1.7 Å，改用联合判据 | `stage_calculator.py`, `identity.py` |
| 13 | 串行 `_ts_results` → 3-pass DAG 调度 | `stage_calculator.py`, `engine.py` |
| 14 | ~~删除 S4 warmup~~ → S4 也保留 Warmup（role-based，非 seed_state-based） | `stage_calculator.py` |
| 15 | 增加完整 thermochemistry 阶段 | `stage_calculator.py` |
| 16 | `usable_for_ml: bool` → `ml_usability: dict` | `stage_calculator.py`, `ml_quality.py` |
| 17 | ALL roles identity validation | `identity.py` |
| 18 | 配置文件重构（含 warmup role-based + initial_hessian per-role） | `defaults.yaml` |
| 19 | 删除未经基准验证的固定耗时百分比 | `docs/` |
| **20** | **INT 正式优化使用 Calc_Hess true** | `stage_calculator.py` |
| **21** | **S4 Warmup 约束距离 = S3 输入几何中的 forming_bond 当前距离** | `stage_calculator.py` |
| **22** | **Warmup → Calc_Hess → OPT 正确顺序（Hessian 在 warmup 输出几何上算）** | `stage_calculator.py` |
| **23** | **INT Rescue 增加 Recalc_Hess N 分级** | `stage_calculator.py` |

### 11.2 P1 — 后续增强

- S3/S4 pairing 逻辑（垂直/松弛/总差值）
- Checkpoint 签名分离（calculation vs scheduling）
- 读适配器兼容旧 manifest schema
- 测试全面更新

### 11.3 不在本方案中的功能

以下功能明确不属于 S3/S4 统一方案：
- IRC 验证 TS 连接关系（属于 S5 Connectivity Validation 后处理阶段）
- Δ-ML 训练（属于 RPH_Postprocess）
- Lewis acid 校准（V3 遗留，V4 不恢复）

---

## 十二、与审核前的差异总结

| 项目 | 审核前 (v1.0) | 审核后 (v1.1) | v1.2 修正 |
|------|--------------|--------------|-----------|
| **引擎** | LowLevelEngine + HighLevelEngine 分别修改 | 唯一 RefinementEngine，旧名称是空子类 | — |
| **Manifest** | s3_low_level_v3 + s4_high_level_v4 | 统一 refinement_manifest_v1 | — |
| **TS Hessian** | `CalcFC` (错误关键词) | `%geom Calc_Hess true` | — |
| **INT Hessian** | model (默认) | model (默认) | **Calc_Hess true** |
| **TS 救援** | `Calcall` (不存在) | read Hessian + TS_Mode / Recalc_Hess N | — |
| **频率** | `UseHess` 复用初始 Hessian | 独立 Freq 在 canonical 几何上 | — |
| **IRC** | `MaxPoints` / `StepSize` / `reverse` | `MaxIter` / 正式位移参数 / `Direction both` | — |
| **INT rescue** | IRC 强制找回 INT | TS endpoint discovery | + mode displacement 优先 + Recalc_Hess 分级 |
| **INT 检测** | 固定 1.7 Å | 归一化反应进度 + molecular graph + RMSD | — |
| **Soft-mode** | 删除 | 保留（与 Calc_Hess 正交） | — |
| **调度** | 串行 `_ts_results` | 3-pass DAG | — |
| **S4 warmup** | 删除 | 删除 | **保留（全部 INT/TS 都 warmup，role-based）** |
| **S4 warmup 约束距离** | — | — | **S3 opt 几何中的当前 forming_bond 距离** |
| **S4 warmup 循环数** | — | — | **INT=4, TS=6（比 S3 的 8/12 短）** |
| **Warmup 顺序** | 未明确 | 未明确 | **WARMUP → Calc_Hess → OPT（Hessian 在 warmup 输出上算）** |
| **调度** | 串行 `_ts_results` 字典 | 3-pass DAG 调度 |
| **S4 warmup** | 一步 sanity optimization | optimized_parent → 无 warmup |
| **热化学** | 未描述 | 完整的 G_composite / H_composite |
| **ML 质量** | `usable_for_ml: bool` | `ml_usability: {property: bool}` |
| **Identity** | product/precursor 跳过 | ALL roles 验证 |
| **百分比** | 固定 10-30% | 删除，待基准测试 |

---

*文档版本: v1.2 — S4 Warmup + INT Calc_Hess 修正版 · 最后更新: 2026-07-20*
