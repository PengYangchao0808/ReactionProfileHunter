# RPH V3.0 错误 S4 清理与新特征提取：后续实施建议

> 日期: 2026-05-23  
> 依据: `docs/RPH_FEATURE_REEXTRACTION_CLEAN_PLAN.md` 与当前代码库/数据实测

---

## 1. 先统一统计口径

调研报告中出现了 `~603` 个非 MAJOR 分支、`80` 个 conditions 等数字。这类数字大概率来自递归扫描，把 xTB/scan 子目录、condition branch 副本或重复 manifest 混入了统计。

后续实现必须使用**严格口径**：

```text
RXN root:        external_data/rph_output_backup/RXN_*
真实 branch:    RXN_*/branches/BR_*       # 仅 immediate children
真实 condition: RXN_*/conditions/COND_*   # 仅 immediate children
condition pair: RXN_*/conditions/COND_*/branches/BR_*
```

当前严格口径实测：

| 指标 | 数量 |
|------|------|
| RXN | 12 |
| 真实 branch | 24 |
| branch names | `BR_MAJOR`, `BR_DR_001` |
| 含 `S3_TS/ts_final.xyz` | 22 |
| 含历史 `S4_Data/features_raw.csv` | 22 |
| 缺 S3 的真实分支 | 2 |
| root condition dirs | 20 |
| condition-branch pairs | 40 |
| recursive `condition_manifest.json` | 60，含 branch 内副本，不应当作 60 个条件 |

结论：`clean-derived` 和 `extract-batch` 必须只遍历 immediate branch，不允许递归 `**/BR*`。

---

## 2. 方案方向修正

你的新前提是：**历史 S4 是错误特征输出，必须清理并用新专用提取器重算**。

因此后续建议为：

1. **保留现有 `rph_features` 包，不新建包。** 现有包已有 `FeatureMiner`、`FeatureContext`、`BaseExtractor`、schema、extractors、thermo/nbo/multiwfn 子模块。
2. **新增模块，而不是推倒重写。** 新增 `cleaner.py`、`inventory.py`、`resolver.py`、`validator.py`、`batch.py` 即可。
3. **历史 S4 先 quarantine，后删除。** 不直接删除，避免不可逆损失。
4. **新提取器必须禁止读取旧派生物。** 任何读取旧 `S4_Data/features_*`、`merged_features.*`、`condition_features_mlr.csv` 的行为都应报错。
5. **缺 S3 的 2 个分支不补提。** 它们只能写 `INCOMPLETE_S3` manifest，不能伪造 CSV。

---

## 3. 推荐实施顺序

### Phase 1: 只读审计 `inspect-raw`

目标：把真实可提取范围确定下来，不改任何文件。

新增：

```text
rph_features/rph_features/inventory.py
```

CLI：

```bash
rph-features inspect-raw \
  --rph-root external_data/rph_output_backup \
  --output rph_features_output_clean/audit
```

输出：

```text
rph_features_output_clean/audit/
├── branch_inventory.csv
├── condition_inventory.csv
├── forbidden_legacy_artifacts.csv
├── extractability_report.json
└── extraction_plan.json
```

关键规则：

- 只扫描 `RXN_*/branches/BR_*` immediate child。
- 只扫描 `RXN_*/conditions/COND_*` immediate child。
- 检测但不读取历史 S4 内容。
- 标记每个分支状态：
  - `READY_RAW_S3`
  - `INCOMPLETE_S3`
  - `HAS_LEGACY_S4_CONTAMINATED`

### Phase 2: 隔离历史派生物 `clean-derived`

新增：

```text
rph_features/rph_features/cleaner.py
```

CLI：

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

移动对象：

```text
RXN_*/branches/BR_*/S4_Data/
RXN_*/conditions/COND_*/branches/BR_*/condition_features_mlr.csv
RXN_*/conditions/COND_*/branches/BR_*/merged_features.csv
RXN_*/conditions/COND_*/branches/BR_*/merged_features.json
RXN_*/conditions/COND_*/branches/BR_*/thermo_calculation.json
RXN_*/conditions/COND_*/dr_prediction.json
RXN_*/reaction_features/dr_prediction_default.json
```

不要移动：

```text
S0_Mechanism/
S1_ConfGeneration/
S2_Retro/
S3_TS/
condition_manifest.json
reaction_manifest.json
branch_manifest.json
small_molecules/
*_Shermo.sum   # 可作为 raw cache，但后续必须重新解析并记录 provenance
```

必须生成：

```text
quarantine_manifest.csv
quarantine_manifest.json
checksum_before_move.json
restore.sh
```

### Phase 3: V3 raw artifact resolver

新增：

```text
rph_features/rph_features/resolver.py
```

职责：从 `branch_root` 解析 S1/S2/S3 raw artifacts，而不是从旧 `S4_Data` 反推。

必须支持：

```text
S1_ConfGeneration/product_min.xyz
S1_ConfGeneration/product/product_global_min.xyz
S1_ConfGeneration/product/finalDFT/conformer_energies.json
S2_Retro/ts_guess.xyz
S2_Retro/intermediate.xyz
S2_Retro/scan_profile.json
S3_TS/ts_final.xyz
S3_TS/ts_opt/berny/ts_guess.log
S3_TS/ts_opt/berny/ts_guess.fchk
S3_TS/S3_intermediate_opt/standard/intermediate_opt.xyz
S3_TS/S3_intermediate_opt/standard/intermediate.log
S3_TS/S3_intermediate_opt/standard/intermediate.fchk
S3_TS/sp_matrix_metadata.json
S3_TS/mechanism_meta.json
```

同步修复：`mech_packager.py` 中对 `S3_intermediate_opt/standard/intermediate_opt.xyz` 的兼容。

### Phase 4: 新专用提取 batch

新增：

```text
rph_features/rph_features/batch.py
```

CLI：

```bash
rph-features extract-batch \
  --rph-root external_data/rph_output_backup \
  --output rph_features_output_clean \
  --source raw-s1s2s3 \
  --ignore-existing-s4 \
  --write-incomplete-manifest
```

原则：

- 复用现有 `FeatureMiner` 和 `BaseExtractor`。
- 新 extractor 继续继承现有 `BaseExtractor`。
- batch 层只负责发现样本、构造 context、写 manifest、聚合错误。
- 单样本特征仍由 `FeatureMiner.run()` 或新 `FeatureEngine` 完成。

建议先给 `FeatureMiner.run()` 添加两个参数：

```python
sample_id: Optional[str] = None
source_manifest: Optional[Path] = None
```

避免当前 `sample_id = output_dir.parent.name` 在嵌套输出中失真。

### Phase 5: condition-level 重算

新增或扩展：

```text
rph_features/rph_features/thermo/condition_thermo.py
rph_features/rph_features/condition.py
```

CLI：

```bash
rph-features extract-conditions \
  --rph-root external_data/rph_output_backup \
  --branch-features-root rph_features_output_clean \
  --output rph_features_output_clean \
  --thermo-source shermo-sum-or-recompute
```

注意：

- `condition_manifest.json` 是可信标签源。
- 历史 `condition_features_mlr.csv`、`merged_features.*` 不可信，不能读取。
- `*_Shermo.sum` 可作为 raw thermo cache，但要重新解析并记录 checksum。

### Phase 6: validate + ML handoff

新增：

```text
rph_features/rph_features/validator.py
```

CLI：

```bash
rph-features validate \
  --features-root rph_features_output_clean \
  --fail-on-legacy-source
```

最终给 `rph_ml` 的输入建议采用兼容 layout：

```bash
rph-features export-rph-layout \
  --features-root rph_features_output_clean \
  --output rph_features_output_clean_rph_layout

rph-features build-dataset \
  --rph-root rph_features_output_clean_rph_layout \
  --output rph_ml_output/dataset_clean
```

---

## 4. 与现有基础设施的整合点

| 现有组件 | 后续用法 |
|----------|----------|
| `FeatureMiner` | 保留为单样本提取 engine，增加 `sample_id/source_manifest` |
| `FeatureContext` | 继续作为插件上下文，新增 condition/raw manifest 字段可选 |
| `BaseExtractor` | 新提取器继续继承它 |
| `schema.py` | 继续写 `features_raw/mlr/deployable/qa`，逐步迁移列定义到 `rph_schemas` |
| `thermo/` | 承接 condition-level thermodynamics 重算 |
| `rph_features.dataset_builder` | 最终 ML 构建器，短期通过 `export-rph-layout` 适配 |
| `ConditionThermoEngine` / `ConditionFeatureMerger` | 旧实现可参考，但不能直接复用历史输出；逻辑需迁入/重写到 `rph_features` |

---

## 5. 必须避免的坑

1. **不要递归扫所有 `BR*` 目录。** 必须限定 `RXN_*/branches/BR_*` immediate child。
2. **不要把 `condition_manifest.json` 的副本数当作 condition 数。** root condition 是 20，branch 内副本会重复。
3. **不要读取旧 S4 CSV 作为输入。** 这是本轮重算的核心约束。
4. **不要为缺 S3 的分支生成伪特征。** 只能写 upstream failed manifest。
5. **不要直接删除历史 S4。** 先 quarantine，保留 restore 能力。
6. **不要在 `rph_features` 重造完整 ML builder。** 数据集构建仍交给 `rph_ml`。

---

## 6. 建议立即执行的最小闭环

最小可交付目标：

```text
inspect-raw → clean-derived --dry-run → clean-derived --apply → extract-batch → validate
```

第一轮不做 condition-level 重算，只做 branch-level clean extraction：

1. 确认 24 个真实 branch。
2. quarantine 22 个历史 S4。
3. 对 22 个有 `ts_final.xyz` 的 branch 重算 branch-level features。
4. 对 2 个缺 S3 的 branch 写 `INCOMPLETE_S3` manifest。
5. validate 确认新输出没有读取旧 S4。

然后第二轮再做 condition-level thermo/features 重算。
