# RPH V3.0 历史错误 S4 清理与专用特征提取器重算方案

> 日期: 2026-05-23  
> 数据源: `external_data/rph_output_backup/`  
> 新前提: 历史 `S4_Data` 与 condition-level feature 输出为错误派生结果，不能复用，只能作为需要隔离的 contaminated artifacts。

---

## 1. 核心纠偏

上一版方案默认“已有 S4 可复用”，这与当前目标不符。新的判断如下：

1. **历史 S4 不可信**：`branches/*/S4_Data/*` 全部视为错误派生数据。
2. **历史 condition features 不可信**：`conditions/*/branches/*/condition_features_mlr.csv`、`merged_features.*`、`thermo_calculation.json` 均视为由旧 S4/旧热力学流程派生，必须重算或隔离。
3. **可信输入只来自 S0–S3 raw artifacts 与实验条件标签**。
4. **新 `rph_features` 必须从 S1/S2/S3/QC 原始产物重新提取**，并生成新的 manifest、checksum、schema signature。
5. **不要直接删除历史 S4**：先 quarantine，生成审计清单，确认新特征通过验证后再决定是否永久删除。

---

## 2. Artifact 信任边界

### 2.1 可信输入层

这些文件可以作为新特征提取输入：

| 层级 | 路径/文件 | 用途 |
|------|-----------|------|
| RXN metadata | `reaction_manifest.json` | reaction_id, condition_ids, theory_signature |
| S0 | `S0_Mechanism/*` | reaction_type, topology, branch plan, atom map |
| Branch metadata | `branches/*/branch_manifest.json` | branch_id, branch status |
| S1 raw | `branches/*/S1_ConfGeneration/` | product geometry, conformer energies, S1 activation features |
| S2 raw | `branches/*/S2_Retro/` | ts_guess, intermediate, scan_profile, forming_bonds source |
| S3 raw | `branches/*/S3_TS/` | ts_final, ts/intermediate logs/fchk/out, sp_matrix_metadata |
| Condition labels | `conditions/*/condition_manifest.json` | temperature, solvent, yield, dr labels |
| Small molecule cache | `small_molecules/*` | small molecule Gibbs/thermo reference，需 checksum 记录 |
| Shermo raw cache | `conditions/*/branches/*/*_Shermo.sum` | 可作为热力学原始缓存，但必须重新解析并记录来源 |

### 2.2 不可信派生层

这些文件不应被新提取器读取为输入：

| 类型 | 路径 |
|------|------|
| 历史分支 S4 | `branches/*/S4_Data/` |
| 历史 condition features | `conditions/*/branches/*/condition_features_mlr.csv` |
| 历史 merged features | `conditions/*/branches/*/merged_features.csv`, `merged_features.json` |
| 历史 condition thermo result | `conditions/*/branches/*/thermo_calculation.json` |
| 历史 DR prediction | `conditions/*/dr_prediction.json`, `reaction_features/dr_prediction_default.json` |
| stale state | RXN-level `pipeline.state`，只可做参考，不可作为完成度判据 |

原则：新提取流程中，任何读取 `S4_Data/features_*.csv` 或 `merged_features.*` 的行为都应触发测试失败。

---

## 3. 清理策略：quarantine-first，不直接删除

### 3.1 新 CLI

新增命令：

```bash
rph-features clean-derived \
  --rph-root external_data/rph_output_backup \
  --quarantine-root external_data/rph_output_backup_invalid_s4_20260523 \
  --dry-run

rph-features clean-derived \
  --rph-root external_data/rph_output_backup \
  --quarantine-root external_data/rph_output_backup_invalid_s4_20260523 \
  --apply
```

### 3.2 迁移对象

`clean-derived` 应移动而不是删除：

```text
RXN_*/branches/BR_*/S4_Data/
RXN_*/conditions/COND_*/branches/BR_*/condition_features_mlr.csv
RXN_*/conditions/COND_*/branches/BR_*/merged_features.csv
RXN_*/conditions/COND_*/branches/BR_*/merged_features.json
RXN_*/conditions/COND_*/branches/BR_*/thermo_calculation.json
RXN_*/conditions/COND_*/dr_prediction.json
RXN_*/reaction_features/dr_prediction_default.json
```

可保留但标记为 raw-cache 的文件：

```text
RXN_*/conditions/COND_*/branches/BR_*/*_Shermo.sum
```

这些 `.sum` 文件是否复用由新 `thermo` 子模块决定，但必须重新解析并写入新 provenance。

### 3.3 清理 manifest

清理后生成：

```text
external_data/rph_output_backup_invalid_s4_20260523/
├── quarantine_manifest.json
├── quarantine_manifest.csv
├── checksum_before_move.json
└── restore.sh
```

每条记录至少包含：

```json
{
  "source_path": "RXN_xxx/branches/BR_MAJOR/S4_Data/features_raw.csv",
  "quarantine_path": "RXN_xxx/branches/BR_MAJOR/S4_Data/features_raw.csv",
  "artifact_type": "legacy_s4_features",
  "sha256": "...",
  "size_bytes": 12345,
  "reason": "legacy_wrong_s4_output",
  "moved_at": "2026-05-23T..."
}
```

### 3.4 清理验证门

清理后必须满足：

```bash
# 旧 S4 不再存在于 active root
find external_data/rph_output_backup -path '*/S4_Data/*' -type f | wc -l

# 历史 merged/condition feature 不再存在于 active root
find external_data/rph_output_backup -name 'merged_features.*' -o -name 'condition_features_mlr.csv' | wc -l

# S0-S3 raw artifacts 仍存在
find external_data/rph_output_backup -path '*/S3_TS/ts_final.xyz' | wc -l
```

预期：

```text
active legacy S4 files = 0
active merged/condition feature files = 0
ts_final.xyz ≈ 22
```

---

## 4. 新专用特征提取器设计

### 4.1 模块结构

建议将 `rph_features` 从“迁移过来的旧 S4”升级为专用后处理包：

```text
rph_features/rph_features/
├── cli.py
├── layout.py              # 识别 V3 branch/condition/root layout
├── inventory.py           # 扫描 S0-S3 raw artifacts，生成提取计划
├── cleaner.py             # clean-derived / quarantine
├── resolver.py            # 从 branch_root 解析 S1/S2/S3/QC artifacts
├── context_builder.py     # 构造 BranchFeatureContext / ConditionFeatureContext
├── engine.py              # 调度插件，不直接关心文件布局
├── extractors/
│   ├── geometry.py
│   ├── ts_quality.py
│   ├── s1_activation.py
│   ├── s2_electronic.py
│   ├── thermo_condition.py
│   └── qc_provenance.py
├── thermo/
│   ├── shermo_parser.py
│   ├── condition_thermo.py
│   └── reference_state.py
├── writer.py              # features_raw/mlr/meta/qa/provenance
├── validator.py           # schema + source purity validation
└── manifest.py            # raw_manifest / computed_manifest / checksum
```

### 4.2 输入解析层

新 resolver 只接受 `branch_root`，例如：

```text
external_data/rph_output_backup/RXN_0d5f57a5/branches/BR_MAJOR
```

解析输出：

```python
BranchArtifacts(
    reaction_id="RXN_0d5f57a5",
    branch_id="BR_MAJOR",
    s1_product_xyz=...,
    s1_conformer_energies=...,
    s2_ts_guess=...,
    s2_intermediate=...,
    s3_ts_final=...,
    s3_ts_log=...,
    s3_ts_fchk=...,
    s3_intermediate_xyz=...,
    s3_intermediate_log=...,
    s3_intermediate_fchk=...,
    mechanism_meta=...,
    sp_matrix_metadata=...,
)
```

必须支持 V3 路径：

```text
S3_TS/ts_final.xyz
S3_TS/ts_opt/berny/ts_guess.log
S3_TS/ts_opt/berny/ts_guess.fchk
S3_TS/S3_intermediate_opt/standard/intermediate_opt.xyz
S3_TS/S3_intermediate_opt/standard/intermediate.log
S3_TS/S3_intermediate_opt/standard/intermediate.fchk
```

### 4.3 禁止读取旧 S4

新提取器必须有硬性 guard：

```python
FORBIDDEN_INPUT_PATTERNS = [
    "*/S4_Data/features_raw.csv",
    "*/S4_Data/features_mlr.csv",
    "*/S4_Data/feature_meta.json",
    "*/merged_features.csv",
    "*/merged_features.json",
    "*/condition_features_mlr.csv",
]
```

任何 extractor 若读取这些路径，应直接报错：

```text
LegacyDerivedArtifactError: old S4/merged features are forbidden as extraction inputs
```

---

## 5. 新提取流程

### 5.1 Branch-level extraction

对每个 `RXN_*/branches/BR_*`：

1. 读取 S0/S1/S2/S3 raw artifacts。
2. 如果缺少 `S3_TS/ts_final.xyz`：
   - 输出 `branch_status.json`，状态为 `INCOMPLETE_S3`。
   - 不生成 fake `features_raw.csv`。
3. 如果 S3 完整：
   - 调用专用 extractor engine。
   - 生成新 S4：

```text
rph_features_output_clean/RXN_xxx/branches/BR_MAJOR/
├── raw_manifest.json
├── computed_manifest.json
├── features/
│   ├── features_raw.csv
│   ├── features_mlr.csv
│   ├── deployable_features.csv
│   ├── qa_metadata.csv
│   └── feature_meta.json
└── provenance/
    ├── source_checksums.json
    ├── extractor_config.json
    └── schema_signature.json
```

### 5.2 Condition-level extraction

对每个 `RXN_*/conditions/COND_*/branches/BR_*`：

1. 读取 `condition_manifest.json`。
2. 读取或重算 Shermo/thermal correction：
   - 优先解析 `*_Shermo.sum` 作为 raw cache；
   - 如果 `.sum` 缺失，则从对应 branch S3/S1 logs/fchk 重算。
3. 与新 branch-level features 合并。
4. 输出：

```text
rph_features_output_clean/RXN_xxx/conditions/COND_x/branches/BR_MAJOR/
├── condition_features_mlr.csv
├── condition_features_raw.csv
├── thermo_results.json
├── merged_features.csv
├── merged_features.json
└── provenance/condition_source_manifest.json
```

### 5.3 预期输出数量

基于当前数据：

| 输出 | 预期数量 |
|------|----------|
| branch-level 新特征 | 22 条 COMPLETE |
| branch-level incomplete manifest | 2 条 INCOMPLETE_S3 |
| condition branch rows | 35 条可完整计算 + 5 条 failed/upstream-missing 记录 |

说明：2 个缺 S3 的 DR 分支对应的条件数约为 5，因此 condition-level 中应明确标记为 upstream failed，而不是填入历史错误特征。

---

## 6. 新 CLI 设计

### 6.1 审计

```bash
rph-features inspect-raw \
  --rph-root external_data/rph_output_backup \
  --output rph_features_output_clean/audit
```

输出：

```text
audit/
├── raw_branch_inventory.csv
├── condition_inventory.csv
├── extractability_report.json
├── forbidden_legacy_artifacts.json
└── plan.json
```

### 6.2 清理历史错误 S4

```bash
rph-features clean-derived \
  --rph-root external_data/rph_output_backup \
  --quarantine-root external_data/rph_output_backup_invalid_s4_20260523 \
  --dry-run

rph-features clean-derived \
  --rph-root external_data/rph_output_backup \
  --quarantine-root external_data/rph_output_backup_invalid_s4_20260523 \
  --apply
```

### 6.3 重算 branch-level features

```bash
rph-features extract-batch \
  --rph-root external_data/rph_output_backup \
  --output rph_features_output_clean \
  --source raw-s1s2s3 \
  --ignore-existing-s4 \
  --write-incomplete-manifest
```

### 6.4 重算 condition-level features

```bash
rph-features extract-conditions \
  --rph-root external_data/rph_output_backup \
  --branch-features-root rph_features_output_clean \
  --output rph_features_output_clean \
  --thermo-source shermo-sum-or-recompute
```

### 6.5 验证新输出

```bash
rph-features validate \
  --features-root rph_features_output_clean \
  --fail-on-legacy-source \
  --schema rph_schemas.feature_schema
```

### 6.6 导出给 rph_ml

两种方式：

**方案 A：输出兼容 RPH root layout**

```bash
rph-features export-rph-layout \
  --features-root rph_features_output_clean \
  --output rph_features_output_clean_rph_layout

rph-features build-dataset \
  --rph-root rph_features_output_clean_rph_layout \
  --output rph_ml_output/dataset_clean
```

**方案 B：扩展 rph_ml 支持 `--features-root`**

```bash
rph-features build-dataset \
  --features-root rph_features_output_clean \
  --output rph_ml_output/dataset_clean
```

短期推荐方案 A，因为现有 `rph_features.dataset_builder` 已经稳定支持 RPH root layout。

---

## 7. 输出 manifest 规范

### 7.1 `raw_manifest.json`

记录所有原始输入：

```json
{
  "reaction_id": "RXN_0d5f57a5",
  "branch_id": "BR_MAJOR",
  "source_root": "external_data/rph_output_backup/RXN_0d5f57a5/branches/BR_MAJOR",
  "trusted_inputs": {
    "s1_product_xyz": {"path": "...", "sha256": "..."},
    "s2_intermediate_xyz": {"path": "...", "sha256": "..."},
    "s3_ts_final_xyz": {"path": "...", "sha256": "..."},
    "s3_ts_log": {"path": "...", "sha256": "..."}
  },
  "forbidden_inputs_detected": [],
  "legacy_s4_ignored": true
}
```

### 7.2 `computed_manifest.json`

记录新提取器版本：

```json
{
  "extractor_package": "rph_features",
  "extractor_version": "0.2.0",
  "schema_signature": "...",
  "git_commit": "...",
  "status": "COMPLETE",
  "warnings": [],
  "outputs": {
    "features_raw_csv": "features/features_raw.csv",
    "features_mlr_csv": "features/features_mlr.csv",
    "feature_meta_json": "features/feature_meta.json"
  }
}
```

### 7.3 `feature_meta.json` 必须包含

```json
{
  "source_policy": "raw_s1s2s3_only",
  "legacy_s4_read": false,
  "condition_thermo_source": "shermo_sum_reparsed",
  "feature_status": "COMPLETE|DEGRADED|FAILED_UPSTREAM",
  "degradation_reasons": []
}
```

---

## 8. 测试计划

### 8.1 清理测试

- `test_clean_derived_dry_run_no_mutation`
- `test_clean_derived_quarantine_preserves_checksums`
- `test_clean_derived_restore_script_roundtrip`

### 8.2 输入纯度测试

- `test_extractor_never_reads_legacy_s4`
- `test_extractor_fails_on_forbidden_merged_features_input`
- `test_feature_meta_marks_legacy_s4_ignored`

### 8.3 V3 layout 解析测试

- `test_resolve_v3_branch_major_artifacts`
- `test_resolve_s3_intermediate_opt_standard`
- `test_missing_ts_final_yields_incomplete_s3_status`

### 8.4 批处理测试

- `test_extract_batch_counts_complete_and_failed`
- `test_extract_conditions_skips_upstream_failed_branch`
- `test_export_rph_layout_compatible_with_rph_ml_builder`

### 8.5 数值/Schema 测试

- `test_features_raw_schema_matches_rph_schemas`
- `test_features_mlr_no_target_leakage`
- `test_condition_features_include_temperature_solvent_yield_dr`

---

## 9. 推荐实施顺序

### Phase 1: 只读审计

实现 `inspect-raw`，确认：

```text
24 branches
22 extractable branches
2 incomplete S3 branches
20 conditions
40 condition-branch pairs
```

### Phase 2: 清理历史派生物

实现 `clean-derived --dry-run` 与 `--apply`，将错误 S4 quarantine。

### Phase 3: 专用 branch extractor

实现 V3 branch artifact resolver + branch-level extract-batch。

### Phase 4: condition thermo extractor

重新解析 Shermo 或重算热力学，生成 condition-level features。

### Phase 5: validation + ML handoff

验证新特征无历史 S4 来源，导出兼容 layout，交给 `rph-features build-dataset`。

---

## 10. 最终目标

完成后应形成：

```text
external_data/rph_output_backup/                 # 只保留可信 S0-S3/raw artifacts
external_data/rph_output_backup_invalid_s4_*/    # 历史错误 S4 quarantine
rph_features_output_clean/                       # 新专用提取器产物
rph_features_output_clean_rph_layout/            # 给 rph_ml 的兼容视图
rph_ml_output/dataset_clean/                     # 新 ML 数据集
```

质量判据：

1. 新特征的 provenance 显示 `legacy_s4_read=false`。
2. 所有新 features 均可追溯到 S1/S2/S3 raw artifact checksum。
3. 2 个缺 S3 分支明确标记为 upstream failed，不生成伪特征。
4. condition-level thermodynamics 由 `.sum` 重新解析或从 raw logs 重新计算。
5. `rph-features build-dataset` 只读取新输出，不读取历史错误 S4。
