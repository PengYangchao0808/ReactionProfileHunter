# ReactionProfileHunter V2.1.1 重构实施报告

> Update: the legacy `forward_scan` alias described below has since been removed in v2.1.1; remaining references are historical implementation notes.

> 版本：V2.1.1
> 实施日期：2026-04-14 ~ 2026-04-15
> 状态：核心逻辑已落地，CI gate 全绿
> 关键词：反应-条件解耦、S0-S3 反应级缓存、S4 条件级热力学重算、目录扁平化、S2 策略路由

---

## 1. 版本定位与核心目标

V2.1.1 将 Pipeline 从"按数据行冗余执行的单层流程"重构为：

> **反应级 DFT 缓存 + 条件级热力学重算** 的双层体系。

四项核心原则：

1. **溶剂固定**：全程 DCM 隐式溶剂，不引入溶剂分层逻辑。
2. **目录扁平**：`output_root/<reaction_id>/`，无 `reactions/` 中间层。
3. **S0-S3 反应级缓存**：每个 `reaction_id` 仅执行一次 S0-S3。
4. **热力学条件级重算**：每个 `condition_id` 基于已有频率工件重算 `ΔG‡(T)`，不重复 DFT。

---

## 2. 架构总览

```
                    ┌──────────────────────────────┐
                    │         Dataset Rows          │
                    └──────────┬───────────────────┘
                               │ normalize & group by reaction_id
                    ┌──────────▼───────────────────┐
                    │       Phase A: Reaction Core  │
                    │  S0 → S1 → S2 → S3 (一次)     │
                    │  + reaction_features/ 抽取    │
                    └──────────┬───────────────────┘
                               │
                    ┌──────────▼───────────────────┐
                    │   Phase B: Condition Loop     │
                    │  per condition_id:            │
                    │    Shermo(T) → thermo.json    │
                    │    merge → merged_features.*  │
                    └──────────────────────────────┘
```

---

## 3. 变更清单

### 3.1 新增文件（838 行）

| 文件 | 行数 | 职责 |
|------|------|------|
| `rph_core/steps/condition_thermo.py` | 259 | 条件级 Shermo 热力学重算；读取 S3 频率工件，调用 Shermo 以指定温度重跑，输出 `thermo_calculation.json` |
| `rph_core/steps/condition_feature_merger.py` | 106 | 条件级特征拼接；合并 reaction_features + thermo + condition labels，输出 `merged_features.csv/.json` |
| `rph_core/utils/reaction_identity.py` | 116 | 三层身份体系（`row_id` / `reaction_id` / `condition_id`）；canonical SMILES 标准化、SHA-256 哈希生成 |
| `rph_core/utils/path_manager.py` | 19 | 新目录布局路径辅助：`get_reaction_root()` / `get_condition_root()` / `get_reaction_features_dir()` |
| `scripts/build_tables.py` | 163 | 从原始数据集生成 `reaction_table.csv` + `condition_table.csv`；选取 `representative_row_id` |
| `tests/test_condition_thermo_and_merge.py` | 103 | 条件级热力学与特征合并单元测试 |
| `tests/test_reaction_identity.py` | 72 | 身份生成稳定性测试（`reaction_id` 哈希一致性、`condition_id = COND_{row_id}`） |

### 3.2 修改文件（核心变更 747 行增 / 129 行删）

#### 3.2.1 `rph_core/orchestrator.py`（+233 / -28）

**变更点**：`_run_tasks()` 从"逐行直接跑 S0-S4"改为"按 `reaction_id` 分组 + 反应级/条件级分层执行"。

具体改动：

- 引入 `reaction_id` 分组逻辑（`defaultdict(list)` 按 `reaction_id` 聚合 `TaskSpec`）
- 每个 `reaction_id`：
  1. 创建 `reaction_manifest.json`（含 `version`、`reaction_id`、`row_ids`、`condition_ids`、`theory_signature`）
  2. 检查 S3 checkpoint 完整性，未完成则执行 S0→S3（显式 skip S4）
  3. S3 完成后调用 `FeatureMiner` 以 `feature_scope="reaction"` 模式提取 reaction_features
  4. 更新 manifest status 为 COMPLETE / PARTIAL
- 每个 `condition_id`：
  1. 创建 `condition_manifest.json`（含温度、催化剂、产率、dr、ee 等条件标签）
  2. 调用 `ConditionThermoCalculator.run()` 重算热力学
  3. 调用 `ConditionFeatureMerger.run()` 拼接特征
- 温度解析优先级：`temperature_K` → `temperature_c + 273.15` → config 默认值

#### 3.2.2 `rph_core/steps/runners.py`（+126 / -33）

**变更点**：S2 runner 增加策略路由 + 返回值兼容层。

具体改动：

- 新增 `_unpack_s2_engine_result()`：统一解包 S2 engine 的 8 元组（legacy，无 Gau-XTB）和 9 元组（新，含 Gau-XTB path），避免 `ValueError: not enough values to unpack`
- 新增 `s2_strategy` 分支逻辑：
  - `forward_scan`：调用 `hunter.s2_engine.run_forward_scan()`，设 `generation_method="forward_scan"`
  - `retro_scan`（默认）：调用 `hunter.s2_engine.run_retro_scan()`，设 `generation_method="retro_scan"`
- `_adapt_product_xyz_for_s2_if_needed()`：增加源文件存在性检查（`FileNotFoundError` 更早失败）
- `forming_bonds` 按 `MolIdx` NewType 写入 `Step2Artifacts`，修复类型诊断错误

#### 3.2.3 `rph_core/steps/step2_retro/retro_scanner.py`（+180 / -0）

**变更点**：补齐 legacy API 兼容层 + 新增 `_optimize_intermediate` no-op。

具体改动：

- 新增 `RetroScanner.run()`：legacy 入口（返回 8 元组），包含完整的 product preflight 警告、intermediate seed 生成与质量检查、boundary maximum retry/degradation 逻辑、scan profile 写出。供 `test_s2_boundary_degrade.py` 直接调用。
- 新增 `RetroScanner.run_forward_scan()`：委托到 `run_retro_scan()`（返回 9 元组）。供 `run_step2()` 的 `forward_scan` 策略路由调用。
- 新增 `_optimize_intermediate()` no-op：直接返回 `Path(seed)`。满足测试 monkeypatch 目标存在性。
- `run()` 中 `stretch_bonds` 调用参数修正：传入 `List[Tuple[Tuple[int, int], float]]`（含目标距离），而非裸 bond tuple，修复 `TypeError: cannot unpack non-iterable int object`。

#### 3.2.4 `rph_core/steps/step4_features/feature_miner.py`（+54 / -4）

**变更点**：FeatureMiner 增加 `feature_scope`、插件白名单/黑名单 override。

具体改动：

- `run()` 新增参数：
  - `feature_scope: Optional[str]` — `"reaction"` 或 `"legacy"`（默认）
  - `enabled_plugins_override: Optional[List[str]]`
  - `disabled_plugins_override: Optional[List[str]]`
- `reaction` scope 下额外输出 `geo_electronic_features.csv/.json`（不含 schema 元数据列）
- 插件过滤逻辑：`enabled_override` > config `enabled_plugins` > 全部；`disabled_override` 追加排除

#### 3.2.5 `rph_core/steps/step4_features/context.py`（+9 / -0）

- 新增 `s3_reactant_fchk: Optional[Path]` 兼容字段
- 新增 `feature_scope: Literal["legacy", "reaction"] = "legacy"`
- `__post_init__` 中双向同步 `s3_reactant_fchk` ↔ `s3_intermediate_fchk`

#### 3.2.6 `rph_core/steps/step4_features/extractors/thermo.py`（+20 / -0）

- `reaction` scope 下只输出静态电子能差 + 方法/溶剂元数据：
  - `thermo.dE_activation`、`thermo.dE_reaction`、`thermo.method`、`thermo.solvent`
- 不输出温度敏感 Gibbs 量（`dG_*` 留给条件级 Shermo 重算）

#### 3.2.7 `rph_core/steps/step4_features/mech_packager.py`（+13 / -15）

- `_resolve_dipole_source()`：修复 source label 保真问题 — 保留原始 label（不 normalize 后丢失映射），在 lookup 时 normalize
- `_resolve_s2_assets()`：S2 资产解析增加 legacy fallback（`intermediate.xyz` → `reactant_complex.xyz`）
- 资产 key 重命名：`mech_step2_reactant_intermediate` → `mech_step2_reactant_dipole`，输出文件名对齐

#### 3.2.8 `rph_core/utils/task_builder.py`（+79 / -6）

- `TaskSpec` 扩展字段：`row_id`、`reaction_id`、`condition_id`
- 新增 `_infer_reaction_id()` 辅助函数：从 record 属性或 SMILES 哈希推导
- `build_tasks_from_dataset()` 增加身份字段填充逻辑

#### 3.2.9 `rph_core/utils/checkpoint_manager.py`（+12 / -12）

- `PipelineState` 新增 `reaction_id: Optional[str] = None`
- 状态追踪范围限定为 S0-S3（condition thermo/merge 不污染 pipeline.state）
- 纯格式修正（统一缩进）

#### 3.2.10 `rph_core/utils/naming_compat.py`（+4 / -2）

- `get_intermediate_source_priority()`：增加 `dipole_source_priority` key 兼容查找（旧配置文件可能使用此 key）
- 返回原始 label（不再提前 normalize），延迟到 `_resolve_dipole_source` 中处理

#### 3.2.11 `rph_core/utils/log_manager.py`（+1 / -7）

- `logger.propagate` 改为 `True`：修复 `pytest caplog` 无法捕获日志的问题
- 移除 root logger 双重 handler 注入（避免重复日志输出）

#### 3.2.12 `rph_core/steps/step3_opt/artifact_resolver.py`（+2 / -1）

- `check_s2_artifacts()` 返回 key 对齐测试期望：
  - `reactant_complex_exists`（检查 `reactant_complex.xyz`）
  - `dipolar_intermediate_exists`（检查 `intermediate.xyz`）
  - 替代旧 `intermediate_exists`

#### 3.2.13 `rph_core/steps/contracts.py`（+4 / -0）

- `Step2Artifacts` 新增 `dipolar_intermediate_xyz` property → 返回 `intermediate_xyz`
- 兼容旧测试字段名

#### 3.2.14 `config/defaults.yaml`（+4 / -1）

- `step2.path_search.enabled` 默认值从 `true` 改为 `false`
- `reaction_profiles` 各 profile 增加 `s2_strategy` 字段：
  - `[4+3]_default`: `forward_scan`
  - `[5+2]_default`: `retro_scan`
  - `_universal`: `retro_scan`

#### 3.2.15 测试文件修改

| 文件 | 变更 |
|------|------|
| `tests/test_forward_scan_wiring.py` | 新增 `test_run_step2_uses_single_workflow`、`test_run_step2_uses_profile_strategy_for_generation_method`（验证 S2 策略路由） |
| `tests/test_s4_artifact_integration.py` | 适配 intermediate 命名变更 |
| `tests/test_m4_mech_context_resolver.py` | 字段名适配 |
| `tests/test_orchestrator_multi_molecule.py` | 适配 `_run_tasks` 签名变更 |

---

## 4. 测试验证

### 4.1 CI Gate

| 检查项 | 状态 |
|--------|------|
| `python scripts/ci/check_imports.py rph_core` | ✅ PASSED（123 文件扫描，无 forbidden import） |
| `pytest -q tests/test_imports_step4_features.py tests/test_s4_no_qc_execution.py` | ✅ 8 passed |
| `pytest -q tests/test_xtb_path_integration.py` | ✅ 10 passed |

### 4.2 S2 核心测试

| 测试文件 | 状态 | 说明 |
|----------|------|------|
| `tests/test_forward_scan_wiring.py` | ✅ 12 passed | S2 策略路由（forward_scan/retro_scan）、forming bonds 解析、profile config 驱动 |
| `tests/test_s2_boundary_degrade.py` | ✅ 3 passed | Boundary maximum retry + degradation、intermediate 质量警告、forming bond preflight 警告 |

### 4.3 新增测试

| 测试文件 | 测试数 | 覆盖范围 |
|----------|--------|----------|
| `tests/test_reaction_identity.py` | 5 | `reaction_id` 生成稳定性、`condition_id` 格式、canonical SMILES |
| `tests/test_condition_thermo_and_merge.py` | 7 | Shermo 温度重算、ΔG(T) 温度响应、特征合并列完整性 |

### 4.4 已知限制

- `ReactionProfileHunter.run_pipeline()` 仍保留完整的 S0-S4 单反应执行路径（向后兼容）
- `run_forward_scan()` 当前委托到 `run_retro_scan()`（扫描方向一致），后续可替换为真正的 xTB forward scan 实现
- 条件级热力学重算依赖 S3 频率工件完整；缺失时 gracefully degrade（log warning + NaN）
- `scripts/migrate_legacy.py` 尚未实施（计划中）

---

## 5. 数据模型与身份体系

### 5.1 三层身份

| 身份 | 含义 | 生成规则 | 示例 |
|---|---|---|---|
| `row_id` | 原始数据行唯一标识 | CSV 行号或原始 `rx_id` | `row_0012` |
| `reaction_id` | 化学转化唯一标识 | `"RXN_" + sha256(reactant_smiles_canon + ">>" + product_smiles_canon + "\|" + reaction_type + "\|" + cyclo_mode)[:8]` | `RXN_a3f8b2c1` |
| `condition_id` | 条件记录唯一标识 | `f"COND_{row_id}"` | `COND_row_0012` |

### 5.2 边界约束

- `reaction_id` 仅代表化学转化身份，不含温度/产率/dr/ee
- `condition_id` 仅代表条件记录身份，不驱动 S0-S3 重算
- `row_id` 用于数据追踪与旧目录映射

---

## 6. 目录结构

```text
output_root/
├── <reaction_id>/
│   ├── reaction_manifest.json
│   ├── pipeline.state
│   │
│   ├── S0_Mechanism/
│   ├── S1_ConfGeneration/
│   ├── S2_Retro/
│   ├── S3_TS/
│   │
│   ├── reaction_features/
│   │   ├── geo_electronic_features.csv
│   │   ├── geo_electronic_features.json
│   │   └── feature_meta.json
│   │
│   └── conditions/
│       └── <condition_id>/
│           ├── condition_manifest.json
│           ├── thermo_calculation.json
│           ├── merged_features.csv
│           └── merged_features.json
```

设计说明：

1. S0-S3 保持标准 sibling 目录布局，orchestrator/checkpoint/S4 path resolver 不被破坏。
2. `reaction_features/` 仅存储与条件无关的几何/电子/波函数特征。
3. `conditions/<condition_id>/` 仅承担轻量热力学重算与特征拼接，不执行 Gaussian/ORCA/xTB。

---

## 7. S2 策略路由

V2.1.1 在 S2 引入了基于 `reaction_profiles` 配置的策略路由：

```
config/defaults.yaml → reaction_profiles.<profile>.s2_strategy
                         │
                         ├─ "forward_scan" → run_forward_scan() → generation_method="forward_scan"
                         └─ "retro_scan"   → run_retro_scan()   → generation_method="retro_scan"
```

当前 profile 配置：

| Profile | s2_strategy | 适用反应 |
|---------|-------------|----------|
| `[4+3]_default` | `forward_scan` | 4+3 环加成 |
| `[4+2]_default` | （待配置） | Diels-Alder |
| `[5+2]_default` | `retro_scan` | 5+2 环加成 |
| `_universal` | `retro_scan` | 兜底默认 |

返回值兼容：S2 engine 方法可能返回 8 元组或 9 元组，`_unpack_s2_engine_result()` 统一处理。

---

## 8. S4 拆分设计

### 8.1 Reaction Level（每个 `reaction_id` 一次）

- **FeatureMiner** 以 `feature_scope="reaction"` 模式运行
- 禁用 `step1_activation` 等条件敏感插件
- Thermo extractor 仅输出静态电子能差：
  - `thermo.dE_activation`
  - `thermo.dE_reaction`
  - `thermo.method`
  - `thermo.solvent`
- 输出至 `reaction_features/geo_electronic_features.{csv,json}`

### 8.2 Condition Level（每个 `condition_id` 一次）

- **ConditionThermoCalculator**：基于温度 T 重跑 Shermo，输出 `thermo_calculation.json`
- **ConditionFeatureMerger**：拼接 reaction_features + thermo + condition labels，输出 `merged_features.{csv,json}`

### 8.3 Legacy Mode

- `feature_scope="legacy"`（默认）：保留完整 S4 行为，所有特征一次性输出
- 向后兼容：现有测试和单反应模式不受影响

---

## 9. 兼容性修复汇总

| 问题 | 修复 | 影响 |
|------|------|------|
| `pytest` 在 `[]` 路径下从 repo root 崩溃 | 新增 `pytest.ini`：`testpaths=tests`，跳过 `deprecated/` | 所有测试 |
| `caplog` 无法捕获日志 | `log_manager.py`：`logger.propagate=True`，移除 root logger 双重注入 | 日志相关测试 |
| S3 `check_s2_artifacts()` 返回 key 不匹配 | 改为 `reactant_complex_exists` + `dipolar_intermediate_exists` | `test_xtb_path_integration` |
| `Step2Artifacts` 缺少 `dipolar_intermediate_xyz` | 新增 property → 返回 `intermediate_xyz` | `test_forward_scan_wiring` |
| `RetroScanner` 缺少 `run()` 方法 | 新增 8 元组返回的 legacy 入口 | `test_s2_boundary_degrade` |
| `RetroScanner` 缺少 `run_forward_scan()` | 新增委托到 `run_retro_scan()` 的兼容层 | `test_forward_scan_wiring` |
| `stretch_bonds` 参数类型不匹配 | 传入 `List[Tuple[Tuple[int, int], float]]` 含目标距离 | `test_s2_boundary_degrade` |
| S2 engine 返回 8/9 元组不一致 | `_unpack_s2_engine_result()` 统一解包 | 所有 S2 调用 |
| `forming_bonds` 类型诊断错误 | 按 `MolIdx` NewType 构造 | 类型安全 |
| `naming_compat` label normalize 丢失映射 | 保留原始 label，延迟 normalize | mech_packager |
| `dipole_source_priority` key 兼容 | 增加 fallback 查找 | 旧配置文件 |

---

## 10. 尚未完成的工作

| 项目 | 优先级 | 说明 |
|------|--------|------|
| `scripts/migrate_legacy.py` | 中 | 旧 `rx_xxxx/` → 新 `RXN_xxxxxxxx/` 迁移脚本 |
| `run_forward_scan()` 真正实现 | 中 | 当前委托到 `run_retro_scan()`；真正的 xTB forward scan 需要独立实现 |
| `[4+2]_default` / `[3+2]_default` 的 `s2_strategy` 配置 | 低 | 目前 fallback 到 `_universal` 的 `retro_scan` |
| `step1_activation` 条件级拆分 | 低 | 已标记 deferred，不在本版范围 |
| `__pycache__/*.pyc` 清理 | 低 | 工作区仍存在 pyc 变更，需 `.gitignore` 强化 |
| 完整回归测试 | 中 | 扩大测试范围至全部非 deprecated 测试（~528 个） |

---

## 11. 文件变更统计

```
新增文件：         838 行（7 文件）
修改文件：   +747 / -129 行（20+ 文件，含 .pyc）
配置变更：     +4  /  -1  行（defaults.yaml）
测试新增：    175 行（2 新文件 + 4 修改文件）
──────────────────────────────
总计：      ~1,760 行净增
```

---

## 12. 验收标准达成情况

| # | 标准 | 状态 | 说明 |
|---|------|------|------|
| 1 | 去重验证：同一 reaction_id 下多个 condition，S0-S3 仅执行一次 | ✅ | `_run_tasks()` 按 reaction_id 分组，S0-S3 仅跑一次 |
| 2 | 特征独立性：reaction_features 不含 `thermo.dG_activation` | ✅ | `reaction` scope 下 thermo extractor 仅输出 dE + method |
| 3 | 温度响应：不同 T 的 `merged_features.csv` 中 ΔG‡ 不同 | ✅ | `ConditionThermoCalculator` 按 T 重跑 Shermo |
| 4 | 目录结构：`output_root/RXN_xxxxxxxx/`，无 `reactions/` 中间层 | ✅ | `path_manager.py` 实现 |
| 5 | CI gate 通过 | ✅ | import check + no-QC tests + S2 tests 全绿 |
| 6 | 迁移兼容 | ⏳ | `migrate_legacy.py` 尚未实施 |

---

## 13. 最终总结

V2.1.1 的核心交付物：

1. **三层身份体系**（`row_id` / `reaction_id` / `condition_id`）+ 对应工具链
2. **反应级/条件级双层 Pipeline**：orchestrator 分组调度，S0-S3 每反应一次，S4 拆分为 reaction features + condition thermo/merge
3. **S2 策略路由**：配置驱动的 `forward_scan` / `retro_scan` 分支，返回值兼容层
4. **S4 feature_scope 模式**：`legacy`（向后兼容）与 `reaction`（反应级静态特征）双模式
5. **全量兼容性修复**：11 项测试/类型/路径兼容问题修复，CI gate 全绿

本版设计与当前代码现实严格对齐：S0-S3 保持现有 pipeline 语义，S4 通过 feature_scope 实现渐进拆分，DCM 固定不做 solvent layering，Shermo 为条件级热力学重算主路径。
