# RPH v3.0.0 → v3.0.1 重算方案

**日期**: 2026-06-29  
**版本**: RPH v3.0.0 → v3.0.1  
**范围**: `rph_output_backup/` 中 12 个 RXN 的选择性重算 + LA 体系全新计算

---

## 0. 变更摘要

| 变更项 | 影响范围 | 严重度 |
|--------|---------|--------|
| mRRHO GFN2→GFN1 回退机制 | S1 构象排序（`G_total = E_sp + G_mRRHO`） | High — 可能使全局最小构象选择错误 |
| CENSO 模块 SCF 收敛检测 | S1 funnel 阶段个别构象 mRRHO 失败 | Medium — 部分构象缺失校正 |
| 路易斯酸原子映射修复 | S1–S3 全链路原子索引 | Critical — LA 体系完全不可用（已修复） |
| Checkpoint 签名纳入 LA metadata | S2/S3 断点续算 | Medium — LA 配置变更可正确触发重算 |

---

## 1. 现有数据评估

### 1.1 备份目录结构

```
rph_output_backup/
├── RXN_0b71b9e9/          ← 12 个反应，全部无 LA (additive: null)
├── RXN_0d5f57a5/
├── ...
├── RXN_f5d5c7b9/
├── rph_output_backup/     ← 嵌套备份（二次安全副本）
└── small_molecules/       ← 小分子缓存（C3H6O2, C3H6O 等）
```

### 1.2 每个反应的内部结构（v3 调度器）

```
RXN_xxxx/
├── pipeline.state                ← 顶层状态（S0=✓, S1–S4 跟踪在 condition/branch 级）
├── reaction_manifest.json        ← 理论级别、状态、条件 ID 列表
├── S0_Mechanism/                 ← 机制分类（可复用）
├── branches/
│   ├── BR_DR_001/                ← 非对映异构体分支 1
│   │   ├── S1_ConfGeneration/product/finalDFT/conformer_thermo.csv
│   │   ├── S1_ConfGeneration/product_min.xyz
│   │   ├── S2_Retro/ts_guess.xyz
│   │   └── S3_TS/...
│   └── BR_MAJOR/                 ← 主产物分支
├── conditions/
│   ├── COND_8/                   ← 实验条件 1 (温度/溶剂/收率/DR)
│   │   ├── condition_manifest.json   ← catalyst: null, additive: null
│   │   └── branches/BR_DR_001/, BR_MAJOR/
│   ├── COND_9/
│   └── COND_10/
└── precursor/S1_ConfGeneration/  ← 前驱体构象搜索
```

### 1.3 关键发现

| 检查项 | 结果 |
|--------|------|
| Lewis acid 数据 | **全部 12 个 RXN 均无 LA**（`additive: null`, `catalyst: null`） |
| 理论级别 | M062X/def2-SVP + wB97M-V/def2-TZVPP + SMD/ACETONE |
| mRRHO 文件 | `conformer_thermo.csv` 存在，含 `g_used` 列（Gibbs 校正值） |
| 条件数 | 每反应 3 个条件（COND_8/9/10），不同温度/溶剂/DR |
| 分支数 | 每条件 2 个分支（BR_DR_001 + BR_MAJOR） |
| 小分子缓存 | 2 个（C3H6O2, C3H6O）—— 需检查是否受 mRRHO 影响 |

---

## 2. 重算分类

### Category A：无 LA 体系 — mRRHO-only 选择性重算

**目标**: 修复 mRRHO/CENSO bug 对热力学校正和模型训练特征的影响，同时最大限度复用既有 S1/S2/S3 结构结果。

**核心判断**: 默认不删除 S1 checkpoint，也不重跑 GFN0/GFN2/DFT 构象搜索；先基于现有构象和既有 DFT/SP 结果补算 mRRHO，并生成旁路结果文件。只有当新 mRRHO 自由能排序明确改变当前全局最小构象时，才对对应 branch 做结构级升级重算。
- **排序不变或差异很小** → `product_min.xyz` 保持 byte-level hash 不变，S2/S3 直接复用
- **排序显著改变** → 只升级受影响 branch，优先复用已有 DFT OPT/SP 构象，再按需级联 S2/S3

### Category B：有 LA 体系 — 全新计算

**目标**: 对需要路易斯酸模拟的反应进行完整管线计算，同时验证 LA 管线兼容性。

---

## 3. Category A：mRRHO-only 选择性重算方案

### 3.1 核心修正

全量重算 S1 是过度保守的默认方案。mRRHO GFN2→GFN1 回退主要影响热力学校正项和 funnel 自由能排序；GFN0/GFN2 初筛、既有 DFT OPT/SP 结构、S2 retro scan 和 S3 TS 优化本身并不因为 mRRHO 补算而必然失效。

Category A 默认改为：

```
保留现有 S1/S2/S3 canonical artifacts
        ↓
基于现有构象和既有 DFT/SP 能量补算 mRRHO
        ↓
生成 sidecar thermochemistry 结果和排序对比报告
        ↓
仅当新 mRRHO 排序显著改变 S1 全局最小构象时，升级少数 branch
```

**默认禁止操作**:
- 不删除 `S1_ConfGeneration/`
- 不修改 `pipeline.state`
- 不重写 `product_min.xyz`
- 不运行完整 `bin/rph_run` 触发 resume 级联

这样可以保持 `product_min.xyz` 的 byte-level hash 不变，S2/S3 checkpoint 复用语义不会被误触发。mRRHO 修复结果作为训练特征侧的旁路 thermochemistry 数据进入后处理流程。

### 3.2 mRRHO-only 重算范围

对每个无 LA RXN 的每个 S1 branch，读取既有产物构象集合和 DFT/SP 能量：

| 输入 | 用途 |
|------|------|
| `S1_ConfGeneration/product/finalDFT/conformer_thermo.csv` | 旧排序、旧 `g_used`、旧 conformer ID |
| `S1_ConfGeneration/product/finalDFT/*` | 现有 DFT OPT/SP 构象与能量来源 |
| `S1_ConfGeneration/product_min.xyz` | 当前 canonical product，仅读取 hash 与 conformer identity |
| `pipeline.state` / provenance | 记录当前 S1/S2/S3 完成状态，不修改 |

补算输出只写 sidecar 文件：

| 输出 | 说明 |
|------|------|
| `S1_ConfGeneration/product/finalDFT/conformer_mrrho_v301.csv` | 新 mRRHO 结果，含 GFN2/GFN1 attempt、`G_mRRHO`、`G_total_v301` |
| `S1_ConfGeneration/product/finalDFT/mrrho_recompute_manifest.json` | 输入 hash、代码版本、xTB 版本、fallback 记录、失败原因 |
| `S1_ConfGeneration/product/finalDFT/mrrho_rank_compare_v301.json` | 新旧排序、top conformer 是否变化、能量差阈值判定 |
| `reports/category_a_mrrho_only_summary.csv` | 全局汇总，供人工验收和 RPH_Postprocess 消费 |

### 3.3 计算规则

对每个已有 conformer：
1. 优先使用已完成 DFT OPT 后的几何作为 mRRHO 几何输入。
2. 使用既有 DFT SP electronic energy 作为 `E_sp`。
3. 运行 xTB mRRHO：先 GFN2；若出现 SCF convergence failure，则回退 GFN1。
4. 计算 `G_total_v301 = E_sp + G_mRRHO_v301`。
5. 不改变原始 `conformer_thermo.csv`，除非后续明确进入“受控覆盖”阶段。

推荐 sidecar schema：

```json
{
  "schema_version": "mrrho_recompute_v1",
  "rph_version": "3.0.1",
  "branch_path": "RXN_xxxx/branches/BR_MAJOR",
  "product_xyz_hash_before": "sha256[:16]",
  "canonical_artifacts_modified": false,
  "conformers": [
    {
      "conformer_id": "conf_001",
      "geometry_source": "finalDFT",
      "sp_energy_hartree": -123.456,
      "mrrho_gfn_attempted": [2, 1],
      "mrrho_gfn_used": 1,
      "g_mrrho_hartree": 0.0123,
      "g_total_v301_hartree": -123.4437,
      "status": "ok"
    }
  ]
}
```

### 3.4 排序门控

使用新旧排序对比决定是否升级：

| 结果 | 判定 | 后续动作 |
|------|------|----------|
| 新 top conformer 与当前 `product_min.xyz` 对应 conformer 一致 | Safe reuse | 不重算 S1/S2/S3；训练特征读取 sidecar mRRHO |
| 新 top 改变，但与旧 top 的差值小于 0.3 kcal/mol | Ambiguous | 不自动级联；标记人工复核或使用 ensemble/低置信度标签 |
| 新 top 改变，且新 top 低于旧 top 超过 0.3-0.5 kcal/mol | Escalate | 仅对该 branch 做结构级升级 |
| mRRHO 对多个关键 conformer 失败 | Incomplete | 补跑失败 conformer 的 mRRHO；失败仍存在则标记不可自动验收 |

阈值建议先用 0.3 kcal/mol 作为敏感门控，0.5 kcal/mol 作为强升级门控。最终阈值应在 1-2 个 RXN 试算后根据排序稳定性调整。

### 3.5 结构级升级路径

只有 `Escalate` branch 才进入结构级重算，且优先级如下：

1. **Promotion 优先**：如果新 top conformer 已有 DFT OPT/SP 输出，则复制该优化后几何到新的受控重算目录，生成新的 `product_min.xyz`，再只对该 branch 重跑 S2/S3。
2. **Targeted DFT**：如果新 top 缺少 DFT OPT/SP，则仅对该 conformer 补跑 DFT OPT/SP，不重启 GFN0/GFN2 全构象搜索。
3. **Full S1 fallback**：只有当 conformer identity 丢失、DFT artifacts 不完整、或 promotion 无法建立可审计链路时，才删除 S1 checkpoint 并完整重跑 S1。

结构级升级必须放在新的工作目录，例如 `rph_output_recompute_escalated/`，原 `rph_output_backup/` 和 mRRHO-only sidecar 结果保持只读可追溯。

### 3.6 Hash 语义说明

S2 signature 使用 `product_min.xyz` 的 byte-level SHA256[:16]，不是 RMSD：

- 标题行、能量精度、浮点末位、原子顺序任一变化都会改变 hash。
- mRRHO-only 方案不写 `product_min.xyz`，所以不会触发 S2/S3 checkpoint 级联。
- RMSD 只用于结构级升级后的事后解释，不参与 resume 自动复用判定。

### 3.7 小分子缓存处理

`small_molecules/` 中的 C3H6O2、C3H6O 不应默认删除缓存。先做同样的 mRRHO-only sidecar 检查：

- 若只有热力学校正值变化，保留结构缓存并更新训练特征侧 mRRHO 数据。
- 若缓存本身缺少必要几何或能量 provenance，再对小分子做 targeted mRRHO/DFT 补算。

### 3.8 预期工作量估算

**v3 调度器结构确认**: S1 构象搜索仅在顶层 branches/ 下，conditions/ 下主要是不同温度/条件的后处理。每个 RXN 的核心 S1 对象约为 2 branches + 1 precursor。

| 项目 | 估算 | 说明 |
|------|------|------|
| mRRHO-only sidecar | 12 RXN × 约 3 S1 × 若干 conformer × GFN mRRHO | GFN 级别，通常小时级而非天级 |
| 排序对比与报告 | < 1 h | 纯文件解析与能量排序 |
| Escalate branch | 预计少数 | 仅排序显著改变时触发 |
| S2/S3 级联 | 预计少数 | 仅结构级升级后触发 |
| **默认总计** | **~2-8 h** | 取决于 conformer 数量和 xTB 并行度 |

> 资源规划应按“mRRHO-only 先快速完成，结构级重算作为少数例外”执行，而不是按“多数 branch 级联 S2/S3”执行。

---

## 4. Category B：路易斯酸体系全新计算方案

### 4.1 适用场景

对实验条件中明确使用了路易斯酸（如 LiCl、MgCl₂、AlCl₃、BF₃ 等）的反应，需要引入 LA 模型进行完整 DFT 计算。

### 4.2 LA 配置

```yaml
# config/defaults.yaml 中启用 LA
lewis_acid:
  enabled: true
  surrogate: "LiCl"                    # 默认 LiCl
  
  surrogate_registry:
    LiCl:
      name: "LiCl"
      smiles: "[Li]Cl"
      charge: 0
      multiplicity: 1
    MgCl2:
      name: "MgCl2"
      smiles: "Cl[Mg]Cl"
      charge: 0
      multiplicity: 1
    AlCl3:
      name: "AlCl3"
      smiles: "Cl[Al](Cl)Cl"
      charge: 0
      multiplicity: 1
  
  s1_sampling:
    use_nci_mode: true
    use_wall_potential: true
    max_li_cl_distance: 3.5
    max_li_substrate_distance: 4.0
  
  s2_quality_gate:
    li_cl_max: 3.5
    li_o_max: 2.6
    cl_reaction_center_min: 3.2
    coordination_switch_allowed: false
  
  s3_quality_gate:
    reaction_mode_projection_min: 0.35
    additive_mode_projection_max: 0.25
    require_irc_for_flagged_ts: true
  
  append_to_precursor: true
```

### 4.3 执行策略

#### Phase 1: 小规模验证（1-3 个反应）

选择 1-3 个实验条件明确使用 LA 的反应，使用 LiCl surrogate 进行完整管线测试：

```bash
# 单反应测试
bin/rph_run \
    --smiles "PRODUCT_SMILES" \
    --output ./Output/la_test_001 \
    --config config/defaults_la_test.yaml
```

**验证检查清单**：
- [ ] S1: `product_min.xyz` 包含 LA 原子（Li + Cl 在末尾）
- [ ] S1: `smiles_to_xyz_map.json` 有 `lewis_acid` 段（v2.0 schema）
- [ ] S1: LA 几何校验通过（metal-Cl < 3.5 Å, metal-O < 4.0 Å）
- [ ] S2: forming bonds 不含 LA 原子索引
- [ ] S2: LA coordination trace 正常记录
- [ ] S2: quality gate 通过（无 dissociation/coordination_switch）
- [ ] S3: TS 虚频模式投影在反应键上（非 LA 运动）
- [ ] S3: `ts_final.xyz` 包含 LA 原子
- [ ] Checkpoint: `pipeline.state` 含 LA metadata

#### Phase 2: 批量计算

验证通过后，对全部 LA 反应进行批量计算：

```bash
# 批量模式
# config 中设置:
#   run.source: dataset
#   run.dataset.path: "data/la_reactions.tsv"
#   run.resume: false  ← 全新计算，不复用旧 checkpoint

bin/rph_run --config config/defaults_la_test_dataset.yaml
```

#### Phase 3: 多 surrogate 对比（可选）

对关键反应，使用不同 surrogate（LiCl vs MgCl₂ vs AlCl₃）进行对比。

> **注意**: 当前 CLI 不支持 `--override` 参数。需为每种 surrogate 生成独立 config 文件：

```bash
# 生成 3 个 config 文件（手动修改 surrogate 字段）
for surrogate in LiCl MgCl2 AlCl3; do
    sed "s/surrogate: \"LiCl\"/surrogate: \"${surrogate}\"/" \
        config/defaults_la_test.yaml > config/defaults_la_test_${surrogate}.yaml
done

# 分别运行
for surrogate in LiCl MgCl2 AlCl3; do
    bin/rph_run \
        --smiles "PRODUCT_SMILES" \
        --output ./Output/la_compare_${surrogate} \
        --config config/defaults_la_test_${surrogate}.yaml
done
```

### 4.4 LA 质量门控验证

每个 LA 反应完成后，检查质量标志：

```python
# 检查 PipelineResult.lewis_acid_quality_flags
flags = result.lewis_acid_quality_flags
if not flags.s2_passed:
    print(f"S2 质量门控失败: {flags.s2_violations}")
    # 常见问题: LA 解离、配位切换、Cl 靠近反应中心
if not flags.s3_passed:
    print(f"S3 质量门控失败: {flags.s3_reason}")
    # 常见问题: 虚频模式在 LA 上而非反应键上
```

---

## 5. 执行时间线

```
Week 1:
  Day 1:   Category A — 盘点既有 S1 conformer/DFT/SP artifacts，生成输入清单
  Day 2:   Category A — 执行 mRRHO-only sidecar 试算（1-2 个 RXN）
  Day 3:   Category A — 批量执行 mRRHO-only sidecar 补算
  Day 4:   Category A — 生成排序对比报告，标记 Safe/Ambiguous/Escalate
  Day 5:   Category A — 仅对 Escalate branch 制定 targeted promotion/DFT 方案

Week 2:
  Day 1-3: Category A — 仅对 Escalate branch 执行结构级升级与必要 S2/S3 级联
  Day 4-5: Category B Phase 1 — LA 小规模验证（1-3 反应）
  Day 6-7: Category B Phase 1 — 质量门控验证

Week 3:
  Day 1-3: Category B Phase 2 — LA 批量计算
  Day 4-5: Category B Phase 3 — 多 surrogate 对比（可选）
  Day 6-7: 数据汇总、特征提取（RPH_Postprocess）
```

---

## 6. 数据完整性保障

### 6.1 备份策略

```bash
# 原始备份永不修改
rph_output_backup/          ← 只读，永久保留

# 重算工作目录
rph_output_recompute/       ← Category A mRRHO-only sidecar 与必要升级结果
rph_output_la/              ← Category B LA 全新计算结果
```

### 6.2 版本标记

```bash
# 重算前，记录版本信息
echo "RPH v3.0.1 (mRRHO GFN fallback + LA atom mapping fix)" > rph_output_recompute/VERSION
echo "Recompute date: $(date -I)" >> rph_output_recompute/VERSION

# 同时更新 version.py
# __version__ = "3.0.1"
```

### 6.3 对比验证

对 Category A，建立旧/新结果对比报告：

| 指标 | 对比方法 | 容差 |
|------|---------|------|
| `product_min.xyz` hash | SHA256[:16] | mRRHO-only 阶段必须完全不变 |
| canonical artifacts modified | manifest 布尔值 | mRRHO-only 阶段必须为 `false` |
| mRRHO completion | conformer 级状态统计 | 关键 conformer 必须成功或有明确失败原因 |
| top conformer identity | 新旧 `G_total` 排序对比 | 一致 → Safe reuse；改变 → 进入门控 |
| ΔΔG(new_top - old_top) | kcal/mol | < 0.3 ambiguous；> 0.3-0.5 才考虑 Escalate |
| S2/S3 artifacts | hash/provenance 抽查 | Safe/Ambiguous branch 不应被改写 |

---

## 7. 风险评估

| 风险 | 概率 | 影响 | 缓解措施 |
|------|------|------|---------|
| sidecar mRRHO 未被 RPH_Postprocess 正确消费 | Medium | 训练特征仍使用旧热力学值 | 在 RPH_Postprocess 增加显式输入优先级或 overlay manifest 校验 |
| conformer identity 无法从旧 artifacts 稳定追踪 | Medium | 无法可靠判断 top 是否改变 | 用 geometry hash、文件名、旧 CSV 行号三重匹配；失败 branch 标记人工复核 |
| mRRHO 修复导致少数 top conformer 显著改变 | Low-Medium | 需要 targeted 结构升级 | 只升级 Escalate branch，优先 promotion/targeted DFT，最后才 full S1 |
| LA S1 CREST 解离 LiCl | Medium | S1 失败 | 已有 NCI mode + wall potential + 几何校验 |
| LA S3 虚频在 LiCl 上 | Medium | TS 无效 | 已有 mode projection 质量门控 |
| Checkpoint 签名不匹配导致意外重算 | Low | 浪费计算 | 已在 v3.0.1 修复签名包含 LA metadata |
| 小分子缓存使用了错误 mRRHO | Low | 间接影响 | 先做 mRRHO-only sidecar 检查，必要时 targeted 补算 |

---

## 8. 验收标准

### Category A 验收
- [ ] 全部 12 个 RXN 的 mRRHO-only sidecar 已用 v3.0.1 代码生成
- [ ] `mrrho_recompute_manifest.json` 记录输入 hash、GFN fallback、失败原因和 `canonical_artifacts_modified=false`
- [ ] `mrrho_rank_compare_v301.json` 已标注每个 branch 的 Safe/Ambiguous/Escalate 状态
- [ ] Safe/Ambiguous branch 的 `product_min.xyz`、S2、S3 artifacts hash 未改变
- [ ] Escalate branch 已完成 targeted promotion/DFT 或明确进入 full S1 fallback
- [ ] RPH_Postprocess 已能显式读取 v3.0.1 mRRHO sidecar，避免继续使用旧热力学值

### Category B 验收
- [ ] LA 管线端到端测试通过（S0→S1→S2→S3 全步骤成功）
- [ ] 全部 LA 质量门控标志已检查（s2_passed, s3_passed）
- [ ] `smiles_to_xyz_map.json` v2.0 schema 正确（含 LA 段）
- [ ] Checkpoint metadata 含完整 LA 信息（metal_atom_index 等）
- [ ] 至少 1 个反应完成了多 surrogate 对比

---

## 附录 A: Checkpoint 操作速查

| 目标 | 操作 |
|------|------|
| 强制全量重算 | `run.resume: false` 或删除整个工作目录 |
| mRRHO-only 默认路径 | 不改 checkpoint；只生成 `conformer_mrrho_v301.csv` 和 manifest sidecar |
| 仅重算 S1 | 删除 `S1_ConfGeneration/` + pipeline.state 中 `step_s1.completed = false`；仅作为 Full S1 fallback |
| 仅重算 S2 | 删除 `S2_Retro/` + pipeline.state 中 `step_s2` 的 `step2_signature` |
| 仅重算 S3 | 删除 `S3_TS/` — S3 签名不匹配会自动触发 |
| 保留 S1 重算 S2+S3 | `--skip-steps s1`（确保 `product_min.xyz` 存在） |

### Resume 机制关键点

1. **S2 复用条件**: `product_xyz_hash`（SHA256[:16]）+ `forming_bonds` + scan params + LA config 全部匹配。**注意: hash 是 byte-level 的——即使坐标数值相同，标题行、能量精度、浮点格式差异都会使 hash 变化，导致 S2 保守重算。**
2. **S3 复用条件**: `input_hashes`（ts_guess/intermediate/product 的 SHA256）+ `step3_signature` + `upstream_step2_signature` 全部匹配
3. **Provenance 重建**: 删除 `pipeline.state` 后，系统**优先**从 `step2_provenance.json` + `step3_provenance.json` 重建状态（精确恢复）。但如果 provenance 文件也缺失，fallback 到 `best_effort` 模式：此时 `reaction_profile=None`、`scan_config=None`、无 `la_additive` → 重算的 signature 不完整 → **不应将 best-effort 恢复视为精确签名等价**。建议操作时保留 provenance 文件。
4. **S1 产物变化 → S2 自动重算**: `product_xyz_hash` 不匹配 → S2 signature 不匹配 → 重算
5. **LA metadata 已纳入签名**: `compute_step2_signature` 现包含 `lewis_acid` 段（enabled/surrogate/elements/metal_index/indices）。LA 配置变更会正确触发 S2 重算。

> Category A 默认 mRRHO-only 路径不应触发上述 resume 机制；只有 Escalate branch 的结构级升级才使用 checkpoint 级联。

---

## 附录 B: 辅助脚本清单

需创建的辅助脚本：

| 脚本 | 用途 |
|------|------|
| `scripts/recompute/inventory_s1_conformers.py` | 盘点现有 branch、conformer、DFT/SP artifacts 和 `product_min.xyz` hash |
| `scripts/recompute/recompute_mrrho_sidecar.py` | 对既有 conformer 执行 GFN2→GFN1 mRRHO-only 补算并写 sidecar |
| `scripts/recompute/compare_mrrho_rankings.py` | 对比旧排序与 `G_total_v301` 排序，输出 Safe/Ambiguous/Escalate |
| `scripts/recompute/generate_recompute_report.py` | 汇总 mRRHO-only 状态、失败 conformer、Escalate branch |
| `scripts/recompute/promote_escalated_conformer.py` | 仅对 Escalate branch 生成受控结构升级目录 |
| `scripts/recompute/clear_s1_checkpoint.py` | Full S1 fallback 专用；禁止作为 Category A 默认入口 |

---

*本方案基于 RPH v3.0.1 代码库分析，结合 checkpoint resume 机制设计和 mRRHO 修复影响评估制定。*
