# V5 验收标准文档

> 本文档定义 ReactionProfileHunter V5 阶段的验收标准（Acceptance Criteria），对应 `AGENTS.md` 中 V5.1/V5.2 + M1–M4 里程碑。

## 1. 验收门禁（Gate）

| Gate | 命令 | 验收标准 | 备注 |
|------|------|----------|------|
| **Gate-0** | `python -m compileall rph_core` | 编译无错误 | 确保代码语法正确 |
| **Gate-1** | `python -c "import rph_core.utils.qc_interface; import rph_core.steps.step4_features.mech_packager"` | 最小导入链成功 | 确保核心模块可导入 |
| **Gate-2** | `pytest --collect-only -q` | pytest **可收集**（无 ERROR） | **必须在 safe path 下运行**（见下方说明） |
| **Gate-3** | `pytest -q tests/test_s4_mech_packager.py tests/test_m4_mech_index_contract.py tests/test_m4_qc_artifacts_mech_index.py` | 纯 mock 测试全部通过 | 禁止调用外部 QC（二进制） |
| **Gate-4** | 手动/脚本运行最小 pipeline + 断言 | 行为合同验证（见 §3） | 不跑真实 Gaussian/ORCA/xTB |

### Gate-2 Safe Path 说明

**问题背景**：pytest 核心对路径中的 `[ ]` 字符有硬限制（会把 `[]` 当作 nodeid 参数化语法的一部分），导致在包含 `[]` 的路径下无法收集测试。

**验收策略**：

- **Gate-2 必须在 safe path（不含 `[]` 的绝对路径）下执行**。  
  示例：在同一磁盘/WSL2 内建立软链接：
  ```
  实际路径：/mnt/e/Calculations/[5+2] Mechain learning/Scripts/ReactionProfileHunter/ReactionProfileHunter_20260121
  Safe 链接：/mnt/e/Calculations/RPH_safe -> 上述实际路径
  ```

- **Toxic Path Safety 验证**：与 pytest 运行无关。  
  “toxic path safety”指的是 QC runner（Gaussian/ORCA/xTB）对路径字符（空格、`[](){}` 等）的敏感性，这些在 `tests/test_sandbox_toxic_paths.py` 中验证。

- **验收证据**：Gate-2 命令输出（截图或日志）显示 `ERROR` 为 0。

---

## 2. 里程碑验收（按 V5.1/V5.2 + M1–M4）

### V5.1（数据集与参考态“输入侧”能力）

| 目标 | 验收要点 |
|------|----------|
| TSV/CSV 反应数据加载 | `ReactionRecord` / `TSVLoader` 支持列名映射、过滤、必填字段校验 |
| 离去小分子目录 | `SmallMoleculeCatalog` + unknown policy（error/warn/skip） |
| 参考态目录/索引 | `ReferenceStatesRunner` 能创建目录结构与 `index.json`（mock 模式） |

### V5.2（S2/QC 框架扩展的“接口侧”能力）

| 目标 | 验收要点 |
|------|----------|
| RetroScanner 扩展 | 支持“neutral precursor（可选）”的契约输出（不破坏旧 run） |
| QC 框架扩展 | 新增 NBO/NMR/Hirshfeld 白名单收集/解析函数（文件级逻辑，不调用外部二进制） |
| 模板占位符一致性 | `gaussian_ts.gjf` 等模板不硬编码 charge/mult |

### M1（S4 机理资产打包基础能力）

| 目标 | 验收要点 |
|------|----------|
| S4 根目录固定命名资产 | `features_raw.csv`、`mech_index.json`、`mech_step2_ts2.xyz`、`mech_step2_reactant_dipole.xyz`、`mech_step2_product.xyz`、`mech_step1_precursor.xyz`（可缺失） |
| 缺失输入不阻断 | `mechanism_status=INCOMPLETE`，`missing_inputs`/`degradation_reasons` 有解释 |
| Toxic path safety | 轻量质量检查（atom count、forming bond window） |

### M2（Resume/Skip 语义 + Schema v1 + 前体回退链）

| 目标 | 验收要点 |
|------|----------|
| `is_step4_complete` | 必须同时满足 `features.csv` + `mech_index.json` + `schema_version` 匹配 |
| 前体来源优先级 | `S1_precursor` → `S2_neutral_precursor` → `S2_reactant_complex` |
| mech_index v1 字段 | `schema_version`、`generated_at`（UTC ISO8601）、`mechanism_status`、`missing_inputs`、`degradation_reasons`、`quality_flags`（三态） |

### M3（Backfill hook + Schema migration + 模板修复）

| 目标 | 验收要点 |
|------|----------|
| feature_miner backfill hook | `features_raw.csv` 存在但 `mech_index` 过期 → 只运行 packager |
| `migrate_mech_index` | 支持旧字段到 v1 映射，保留 deprecated alias |
| QC mock 测试 | 文件系统 + regex 级测试，不调用外部 QC |

### M4（QC Artifacts 集成 + 质量标志扩展 + 更新原因常量化）

| 目标 | 验收要点 |
|------|----------|
| `qc_artifacts` 字段 | NMR/Hirshfeld/NBO 输出文件的打包结果（相对文件名 + 可选 source_paths） |
| 质量标志扩展 | `ts_imag_freq_ok`、`asset_hash_ok`（三态） |
| `is_mech_index_up_to_date` | 返回 `UpdateReason` 常量（无副作用，不自动写回） |
| 纯 mock 测试 | 所有测试为 mock/文件级，不触发 subprocess |

---

## 3. Gate-4 行为合同验证（最小 pipeline 树）


## 4. V6.1: S4 Output Contract Update
### V6.1 S4 根目录固定命名资产
| S4 根目录固定命名资产 | `features_raw.csv`、`features_mlr.csv`、`feature_meta.json`、`mech_index.json`、`mech_step2_ts2.xyz`、`mech_step2_reactant_dipole.xyz`、`mech_step2_product.xyz`、`mech_step1_precursor.xyz`（可缺失） |
| 缺失输入不阻断 | `mechanism_status=INCOMPLETE`，`missing_inputs`/`degradation_reasons` 有解释 |
| Toxic path safety | 潤量质量检查（atom count、forming bond window） |

### V6.1 S4 QC Guardrails
**"S4 Never Runs QC" Enforcement:**
- `can_submit_jobs` attribute added to `BaseExtractor` (default: False)
- Extract-only plugins (thermo, geometry, qc_checks, ts_quality) have `can_submit_jobs=False`
- Phase C/B plugins (interaction_analysis, nics, nbo_e2) may have `can_submit_jobs=True` to generate job_specs
- In job_run_policy=disallow mode, plugins may only generate job_specs, never execute them
- Hard safety test added to detect S4 QC execution attempts (tests/test_s4_no_qc_execution.py)

### S4 输出文件变更
| 文件名 | 旧文件名 | 新文件名 | 说明 |
|---------|---------|---------|---------|
| features_raw.csv | features.csv | V6.1 主输出（full features + QC columns） |
| features_mlr.csv | N/A | V6.1 MLR就绪输出（≤10列） |
| feature_meta.json | feature_meta.json | V6.1 扩展（units/sources/missing_stats/mlr_columns） |

### 验收要点（V6.1）
- `python -m pytest -q` passes
- `features_raw.csv`, `features_mlr.csv`, `feature_meta.json` 存在且非空
- `features.csv` 不再被写入（已被 `features_raw.csv` 替代）
- `DEFAULT_ENABLED_PLUGINS = ["thermo", "geometry", "qc_checks", "ts_quality"]`
- 新增 `ts_quality.py` extractor提供 `ts.n_imag`, `ts.imag1_cm1_abs`, `ts.dipole_debye`



在临时目录构造以下树（使用 mock 文件），调用 `pack_mechanism_assets()` 后断言：

```
tmp_path/
├── pipeline_root/
│   ├── S1_Product/
│   │   └── product_min.xyz
│   ├── S2_Retro/
│   │   ├── ts_guess.xyz
│   │   └── reactant_complex.xyz
│   ├── S3_TS/
│   │   ├── ts_final.xyz
│   │   ├── reactant_sp.xyz
│   │   └── hirshfeld/
│   │       └── job_hirshfeld.dat
│   └── S4_Data/
│       └── (packager 输出)
```

### 断言点

1. **固定命名资产**：
   - `S4_Data/mech_step2_ts2.xyz`
   - `S4_Data/mech_step2_reactant_dipole.xyz`
   - `S4_Data/mech_step2_product.xyz`
   - `S4_Data/mech_step1_meta.json`
   - `S4_Data/mech_step2_meta.json`
   - `S4_Data/mech_index.json`

2. **mech_index.json 结构**：
   - 顶层字段：`schema_version`、`generated_at`、`mechanism_status`
   - `assets` 字典：键为固定命名，值为 `{"filename": "...", "source_label": "..."}` 或 `null`
   - `quality_flags`：`atom_count_ok`、`forming_bond_window_ok`、`suspect_optimized_to_product`（三态）
   - `missing_inputs` / `degradation_reasons`：解释缺失/降级原因
   - `qc_artifacts`：`hirshfeld_outputs` 等键，`filename` 为相对路径，`meta` 含 `picked` / `candidates` / `reason`

3. **缺失输入降级**：
   - 删除 `S1_Product/product_min.xyz` → `mechanism_status=INCOMPLETE`，`missing_inputs` 包含 `mech_step2_product.xyz`

4. **qc_artifacts 收集**：
   - `S4_Data/qc_hirshfeld.dat` 存在
   - `mech_index.json['qc_artifacts']['hirshfeld_outputs']['filename'] == "qc_hirshfeld.dat"`
   - `meta.reason == "picked_by_mtime"`

5. **Subprocess 调用禁止**：
   - 测试中 monkeypatch `subprocess.run/Popen` 为直接 raise，确保 packager 路径不触发

---

## 4. 测试矩阵

| 里程碑 | 必须通过的测试类型 | 禁止项 |
|--------|-------------------|--------|
| V5.1 | TSVLoader / Catalog / ReferenceStatesRunner mock | 不读外部大数据、不跑 QC |
| V5.2 | 模板纯渲染 + whitelist 收集/解析 | 禁止调用 gaussian/orca/xtb |
| M1 | packager 合同测试（固定命名、降级） | 不依赖真实 S1/S2/S3 计算 |
| M2 | resume/skip 语义测试 + schema v1 字段测试 | 禁止真实 pipeline 运行 |
| M3 | backfill hook + migrate 测试 | 禁止真实 QC |
| M4 | qc_artifacts/UpdateReason/flags 扩展测试（mock） | 禁止 subprocess、禁止真实 QC |

---

## 5. 证据输出（PR/记录应提交）

每次里程碑验收应附：

1. **Gate-0/1/2 命令输出**（Gate-2 必须是 safe path 截图/日志）
2. **关键测试集通过摘要**（`pytest -q` 输出）
3. **最小样例 `mech_index.json`**（脱敏路径，仅展示关键字段结构）
4. **Schema 变化迁移说明**（如果涉及）

---

## 6. “通过/不通过”判定口径

V5 视为“达成”的必要且充分条件：

1. **仓库健康**：Gate-0/1/2 全部通过（safe path）
2. **S4 机理打包合同达成**：M1+M2+M3 的合同测试全部通过；`mech_index_v1` 稳定输出
3. **可恢复性达成**：任意入口（resume/backfill/feature_miner）能补齐 `mech_index`，而不重复跑昂贵步骤
4. **QC artifacts 可选增强达成（M4）**：
   - `qc_artifacts` 字段存在且可为“空/部分填充”
   - `quality_flags` 扩展字段在 `mech_index` 中存在（值允许 null）
   - 所有相关测试为 mock，不触发真实 QC
5. **与机理假设一致**：第一步不做 TS，不引入重复计算；偶极子态来源以 S3 reactant 为主、S2 为兜底，并在 meta 中记录来源与降级原因

---

## 7. 常见问题（FAQ）

### Q1: 为什么 Gate-2 必须用 safe path？

A: pytest 核心对路径中的 `[ ]` 字符有硬限制，会把 `[]` 当作 nodeid 参数化语法的一部分，导致 `ERROR: path cannot contain [] parametrization`。这不是代码 bug，而是 pytest 工具的限制。为避免 pytest 限制污染“toxic path safety”的业务目标，我们把 Gate-2 验收限定在 safe path。

### Q2: 如何创建 safe path？

A: 在同一磁盘/WSL2 内创建软链接：
```bash
# 假设实际路径包含 []
ln -s "/mnt/e/Calculations/[5+2] Mechain learning/Scripts/ReactionProfileHunter/ReactionProfileHunter_20260121" "/mnt/e/Calculations/RPH_safe"
```
然后在 `/mnt/e/Calculations/RPH_safe` 下运行 Gate-2/3/4。

### Q3: `is_mech_index_up_to_date` 是否会写回文件？

A: **不会**。`is_mech_index_up_to_date()` 是纯检查函数（无副作用），返回 `(is_up_to_date: bool, reason: UpdateReason)`。若需要迁移并写回，使用新增的 `ensure_mech_index_schema()` 函数。

### Q4: qc_artifacts 的路径格式是什么？

A: `filename` 为 **S4 根目录下的相对文件名**（如 `"qc_hirshfeld.dat"`），绝对路径只在 `meta.source_paths` 中作为追溯信息出现，以保证 S4 包的**可移植性**。

### Q5: 如何验证测试没有触发真实 QC？

A: 测试中应使用 `monkeypatch.setattr("subprocess.run", lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("Subprocess not allowed in mock test")))` 等方式显式禁止 subprocess 调用。
