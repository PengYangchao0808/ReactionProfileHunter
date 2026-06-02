# RPH V3.0 特征提取模块优化方案

> 日期: 2026-05-23  
> 范围: `rph_features` 模块、`external_data/rph_output_backup` 历史输出、与 `rph_ml` 的数据集构建衔接

---

## 1. 审核结论

上一版“从 `rph_output_backup` 批量重新提取特征”的方案需要修正。当前数据并不是大面积缺少 S4，而是：

| 项目 | 实际状态 |
|------|----------|
| RXN 数 | 12 |
| 分支数 | 24 (`BR_MAJOR` + `BR_DR_001`) |
| 条件数 | 20 |
| 已有分支级 S4 特征 | 22/24 |
| 已有条件级特征 | 40 个 `condition_features_mlr.csv` |
| 缺失 S4 的分支 | 2 个，且均缺少可用 S3 `ts_final.xyz` |

因此，当前最优策略不是“大量重跑 S4”，而是将 `rph_features` 从单次 `extract` 工具升级为：

```text
inspect  →  normalize  →  repair  →  validate  →  handoff to rph_ml
审计       标准化         必要时补提     质量门控      交给 ML 构建器
```

即：

1. **优先复用已有 S4 和 condition 特征**，避免破坏已完成结果。
2. **只对“有完整 S3 但缺 S4”的分支补提**；当前数据中没有这种分支。
3. **对“缺 S3”的分支只记录降级状态**，不要尝试伪造特征。
4. **数据集构建继续交给 `rph-features build-dataset`**，不要在 `rph_features` 中重复实现 ML 聚合逻辑。

---

## 2. 当前数据结构与完成度

### 2.1 数据结构

`external_data/rph_output_backup` 是 V3.0 风格输出：

```text
RXN_xxx/
├── reaction_manifest.json
├── S0_Mechanism/
├── branches/
│   ├── BR_MAJOR/
│   │   ├── S1_ConfGeneration/
│   │   ├── S2_Retro/
│   │   ├── S3_TS/
│   │   └── S4_Data/
│   │       ├── features_raw.csv
│   │       ├── features_mlr.csv
│   │       ├── feature_meta.json
│   │       ├── deployable_features.csv
│   │       └── qa_metadata.csv
│   └── BR_DR_001/
│       ├── S1_ConfGeneration/
│       ├── S2_Retro/
│       ├── S3_TS/
│       └── S4_Data/       # 10/12 存在
└── conditions/
    └── COND_x/
        ├── condition_manifest.json
        ├── dr_prediction.json
        └── branches/
            ├── BR_MAJOR/
            │   ├── condition_features_mlr.csv
            │   ├── merged_features.csv
            │   ├── merged_features.json
            │   └── thermo_calculation.json
            └── BR_DR_001/
                ├── condition_features_mlr.csv
                ├── merged_features.csv
                ├── merged_features.json
                └── thermo_calculation.json
```

### 2.2 完成度

| 层级 | 完成状态 |
|------|----------|
| `BR_MAJOR` | 12/12 有 S3 + S4 |
| `BR_DR_001` | 10/12 有 S3 + S4 |
| 缺失分支 | `RXN_0b71b9e9/BR_DR_001`, `RXN_1be0da47/BR_DR_001` |
| 条件级分支特征 | 40/40 存在 |

两个缺失 S4 的 DR 分支都缺少有效 `S3_TS/ts_final.xyz`，所以不能通过 `rph_features` 单独补救；必须回到 S3 或标记为 `INCOMPLETE_S3`。

---

## 3. 当前 `rph_features` 的问题

### 3.1 CLI 过于单次化

当前 `rph_features/rph_features/cli.py` 只有：

```bash
rph-features extract --rph-run <single-run-dir> --output <output-dir>
```

问题：

- 不会自动发现 `RXN_*/branches/BR_*`。
- 不会识别 `conditions/COND_*`。
- 没有 `inspect`、`validate`、`normalize`、`repair`。
- 没有 `--skip-existing`，容易覆盖已有 S4。

### 3.2 `FeatureMiner` 是单样本执行器

`FeatureMiner.run()` 一次只写一行 `features_raw.csv` / `features_mlr.csv`，适合“一个分支的一次提取”，不适合批量管控。

应保留其定位：**低层单样本 extractor engine**。批处理应放在新建的 batch/orchestrator 层。

### 3.3 路径解析仍偏 V2/S4_Data 假设

`mech_packager.resolve_mechanism_context()` 通过 `s4_dir.parent` 寻找 `S1_ConfGeneration/S2_Retro/S3_TS`。这对 `branches/BR_MAJOR/S4_Data` 有效，但对 `conditions/COND_x/branches/BR_MAJOR` 无效，因为条件目录只有 post-processing 文件，没有 S1/S2/S3。

### 3.4 S3 中间体命名兼容不足

当前 `_resolve_s3_assets()` 偏好 `reactant_sp.xyz`，但 V3 输出常见位置是：

```text
S3_TS/S3_intermediate_opt/standard/intermediate_opt.xyz
S3_TS/S3_intermediate_opt/standard/intermediate.fchk
S3_TS/S3_intermediate_opt/standard/intermediate.log
```

需要加入候选路径，否则中间体相关特征会降级。

### 3.5 condition metadata 没有进入 `rph_features`

`condition_manifest.json` 中已有：

- `temperature_K`
- `solvent`
- `yield_pct`
- `dr_major`
- `dr_minor`

但 `rph_features` 当前没有显式的 condition ingestion 层。这个设计不一定错误，因为条件级聚合已经由 `rph_features.dataset_builder` 处理；关键是边界要明确。

---

## 4. 优化后的模块边界

### 4.1 职责划分

| 模块 | 应负责 | 不应负责 |
|------|--------|----------|
| `rph_core` | S0–S3 机理计算；输出分支级 raw artifacts | 不再默认跑 S4；不生成 ML 数据集 |
| `rph_features` | 审计、标准化、补提、验证分支级 S4 artifacts | 不训练模型；不重新实现 ML dataset builder |
| `rph_ml` | 读取 S4 + condition features，构建 ML 表并训练 | 不读取 `.log/.fchk/.xyz`，不做 QC 特征提取 |
| `rph_schemas` | 列定义、artifact contract、manifest schema | 不做 IO 和计算 |

### 4.2 推荐数据流

```text
external_data/rph_output_backup
        │
        ▼
rph-features inspect
        │  生成 inventory + quality report
        ▼
rph-features normalize
        │  复用已有 S4_Data，标准化到 rph_features_output
        ▼
rph-features repair --only-missing
        │  只补提有完整 S3 但缺 S4 的分支
        ▼
rph-features validate
        │  验证 schema、manifest、缺失列、NaN 比例
        ▼
rph-features build-dataset
        │  构建 stage1_branch_yield_dataset 等 ML 表
        ▼
rph-ml train-yield / diagnose / select-features / visualize
```

---

## 5. `rph_features` 应新增的核心能力

### 5.1 新增 CLI

建议将 `rph-features` 扩展为以下命令：

```bash
# 1. 只读审计，不写特征
rph-features inspect \
  --rph-root external_data/rph_output_backup \
  --output rph_features_output/audit

# 2. 标准化已有 S4，不重跑提取
rph-features normalize \
  --rph-root external_data/rph_output_backup \
  --output rph_features_output \
  --mode symlink \
  --skip-existing

# 3. 仅补提缺失但可恢复的分支
rph-features repair \
  --rph-root external_data/rph_output_backup \
  --output rph_features_output \
  --only-missing \
  --skip-incomplete-s3

# 4. 单分支提取，保留现有能力但增强参数
rph-features extract \
  --rph-run external_data/rph_output_backup/RXN_0d5f57a5/branches/BR_DR_001 \
  --output rph_features_output/RXN_0d5f57a5/branches/BR_DR_001 \
  --sample-id RXN_0d5f57a5__BR_DR_001 \
  --skip-existing

# 5. 验证标准化后的输出
rph-features validate \
  --features-root rph_features_output \
  --schema rph_schemas.feature_schema
```

### 5.2 新增内部模块

建议新增文件：

```text
rph_features/rph_features/
├── layout.py          # 自动识别 flat / branch / condition layout
├── inventory.py       # 生成 extraction_inventory.csv/json
├── normalizer.py      # 复用已有 S4_Data，标准化输出结构
├── repair.py          # 对可恢复分支调用 FeatureMiner.run()
├── validator.py       # 检查 schema、文件完整性、NaN、状态
├── manifest.py        # raw_manifest/computed_manifest/provenance
└── batch.py           # 批处理调度与错误聚合
```

### 5.3 核心数据模型

建议用 dataclass 表示扫描结果：

```python
@dataclass
class BranchFeatureRecord:
    reaction_id: str
    branch_id: str
    branch_root: Path
    s1_dir: Optional[Path]
    s2_dir: Optional[Path]
    s3_dir: Optional[Path]
    s4_dir: Optional[Path]
    has_ts_final: bool
    has_features_raw: bool
    has_features_mlr: bool
    has_feature_meta: bool
    has_qa_metadata: bool
    status: str  # COMPLETE | EXISTING_S4 | INCOMPLETE_S3 | MISSING_S4_REPAIRABLE | DEGRADED | FAILED
```

条件层单独建模：

```python
@dataclass
class ConditionFeatureRecord:
    reaction_id: str
    condition_id: str
    branch_id: str
    condition_root: Path
    condition_manifest: Path
    condition_features_mlr: Optional[Path]
    merged_features_json: Optional[Path]
    thermo_calculation_json: Optional[Path]
```

---

## 6. 针对当前备份数据的推荐操作

当前数据中 22/24 分支已有 S4，40 个条件分支特征也已存在。因此建议按以下方式处理：

### Step 1: 审计，不重跑

目标输出：

```text
rph_features_output/audit/
├── branch_inventory.csv
├── condition_inventory.csv
├── missing_artifacts.json
├── quality_summary.json
└── extraction_plan.json
```

预期结论：

```text
branch_count: 24
existing_s4: 22
missing_s4: 2
repairable_missing_s4: 0
incomplete_s3: 2
condition_feature_rows: 40
```

### Step 2: 标准化已有 S4

不要覆盖原始 `external_data`，只在 `rph_features_output` 建立标准化视图：

```text
rph_features_output/
├── RXN_0d5f57a5/
│   └── branches/
│       ├── BR_MAJOR/
│       │   ├── features/features_raw.csv
│       │   ├── features/features_mlr.csv
│       │   ├── features/feature_meta.json
│       │   └── provenance/source_manifest.json
│       └── BR_DR_001/
│           └── features/...
└── _manifests/
    ├── extraction_inventory.csv
    ├── normalized_manifest.json
    └── checksum_report.json
```

实现方式优先 `symlink`，避免复制大量文件。

### Step 3: 对两个缺失 DR 分支降级记录

这两个分支没有完整 S3，不能通过 `rph_features` 补提：

```text
RXN_0b71b9e9 / BR_DR_001 → INCOMPLETE_S3
RXN_1be0da47 / BR_DR_001 → INCOMPLETE_S3
```

应生成占位 manifest，而不是生成伪 CSV：

```json
{
  "reaction_id": "RXN_0b71b9e9",
  "branch_id": "BR_DR_001",
  "status": "INCOMPLETE_S3",
  "reason": "missing S3_TS/ts_final.xyz",
  "repairable": false
}
```

### Step 4: 交给 `rph_ml` 构建数据集

当前 `rph_features/rph_features/dataset_builder.py` 已经支持 `rph_output_backup` 格式：

```bash
rph-features build-dataset \
  --rph-root external_data/rph_output_backup \
  --output rph_ml_output/dataset \
  --barrier-threshold 10.0
```

可选只看主分支：

```bash
rph-features build-dataset \
  --rph-root external_data/rph_output_backup \
  --output rph_ml_output/dataset_major \
  --branch BR_MAJOR
```

---

## 7. 必须修复的代码点

### 7.1 修复 S3 中间体路径解析

位置：`rph_features/rph_features/mech_packager.py::_resolve_s3_assets()`

增加候选路径：

```text
S3_TS/S3_intermediate_opt/standard/intermediate_opt.xyz
S3_TS/S3_intermediate_opt/standard/intermediate.xyz
S3_TS/S3_intermediate_opt/standard/intermediate.log
S3_TS/S3_intermediate_opt/standard/intermediate.fchk
```

原因：V3 输出中的中间体不是 `reactant_sp.xyz`，而是 `intermediate_opt.xyz`。

### 7.2 `FeatureMiner.run()` 增加 `sample_id`

当前 `sample_id = output_dir.parent.name`，对嵌套结构不可靠。建议增加显式参数：

```python
def run(..., sample_id: Optional[str] = None, source_manifest: Optional[Path] = None):
    resolved_sample_id = sample_id or output_dir.parent.name
```

推荐格式：

```text
RXN_0d5f57a5__BR_MAJOR
RXN_0d5f57a5__COND_15__BR_MAJOR
```

### 7.3 `extract` 增加安全选项

```bash
--skip-existing
--force
--sample-id
--feature-scope branch|condition
--dry-run
```

默认行为必须是 **不覆盖已有 S4**。

### 7.4 明确 condition 与 branch 的边界

不建议让 `FeatureMiner` 直接处理 `conditions/COND_*` 目录，因为那里没有 raw QC artifacts。条件层应由 `rph_features.dataset_builder` 或未来的 `rph_features.condition_normalizer` 读取：

- `condition_manifest.json`
- `condition_features_mlr.csv`
- `merged_features.json`
- `thermo_calculation.json`

---

## 8. 质量门控

`rph-features validate` 应至少检查：

| 检查 | 规则 |
|------|------|
| 文件完整性 | `features_raw.csv`, `features_mlr.csv`, `feature_meta.json`, `qa_metadata.csv` |
| Schema | 与 `rph_schemas.feature_schema` 对齐 |
| 行数 | 每个 feature CSV 应为 1 行 |
| 状态 | `feature_status` 只能是 COMPLETE/DEGRADED/FAILED |
| NaN 比例 | 输出每列缺失率；关键列缺失报警 |
| S3-S4 对应 | 有 S4 的分支必须能追溯到 S3 `ts_final.xyz` 或记录 legacy source |
| 条件层 | `condition_manifest.json` 与 `condition_features_mlr.csv` 数量匹配 |
| 不可恢复分支 | 缺 S3 的分支只允许 manifest 占位，不允许伪造 CSV |

---

## 9. 目标命令流

最终希望用户操作变成：

```bash
# 1. 审计历史输出
rph-features inspect \
  --rph-root external_data/rph_output_backup \
  --output rph_features_output/audit

# 2. 标准化已有特征，不重算
rph-features normalize \
  --rph-root external_data/rph_output_backup \
  --output rph_features_output \
  --mode symlink \
  --skip-existing

# 3. 验证标准化输出
rph-features validate \
  --features-root rph_features_output

# 4. 直接从原始 RPH 输出构建 ML 表
rph-features build-dataset \
  --rph-root external_data/rph_output_backup \
  --output rph_ml_output/dataset

# 5. 模型诊断与训练
rph-features diagnose --output rph_ml_output/diagnostics
rph-ml train-yield \
  --dataset rph_ml_output/dataset/stage1_branch_yield_dataset.csv \
  --output rph_ml_output/models
```

说明：第 4 步当前已经可以用 `rph_features dataset_builder` 的 existing logic 支持；第 1–3 步是需要补强的 `rph_features` 能力。

---

## 10. 实施优先级

### P0: 立即修正方案认知

- 当前备份数据已有 22/24 分支特征，不要默认重跑。
- 2 个缺失分支不可由 S4 补救，需回到 S3 或降级。
- ML 聚合优先使用 `rph-features build-dataset --rph-root external_data/rph_output_backup`。

### P1: 最小可用优化

1. 新增 `rph-features inspect`。
2. 新增 `rph-features normalize`。
3. 修复 `_resolve_s3_assets()` 对 `S3_intermediate_opt/standard` 的支持。
4. `extract` 增加 `--skip-existing`。
5. 增加 `sample_id` 参数。

### P2: 质量和可复现性

1. 新增 `computed_manifest.json`。
2. 新增 checksum report。
3. 新增 schema validation。
4. 新增 dry-run 计划输出。

### P3: 长期改造

1. 将 condition normalizer 独立为 `rph_features.condition`。
2. 将 shared path/artifact schemas 下沉到 `rph_schemas`。
3. 未来可将 `rph_features` 完全去除对 `rph_core.utils` 的依赖。

---

## 11. 推荐最终定位

`rph_features` 不应只是“运行一次 FeatureMiner 的 CLI”，而应成为 **RPH 后处理资产管理层**：

```text
raw RPH output  →  audited feature assets  →  validated feature tables  →  rph_ml datasets
```

它的核心价值不是重复计算，而是：

- 看清哪些分支可用；
- 安全复用已有 S4；
- 对缺失分支做明确降级；
- 统一 manifests 和 schema；
- 给 `rph_ml` 一个可追溯、可验证、可复现的输入基础。
