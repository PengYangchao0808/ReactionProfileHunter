# 20 条数据集完整运行方案

**最后更新**: 2026-05-07  
**RPH 版本**: v2.1.1 + S0-DR branch  
**数据集**: `data/reaxys_cleaned.csv` (20 行, 12 个独立反应)

---

## 一、数据集概况

| 属性 | 值 |
|------|-----|
| 总行数 | 20 |
| 独立反应 (unique rxn_key_hash) | 12 |
| 反应类型 | 全部 [4+3] |
| 拓扑类型 | 全部 INTRA_TYPE_I (分子内) |
| 有报告 DR 数据 | 16/20 (rx_id=5–20) |
| 无报告 DR 数据 | 4/20 (rx_id=1–4) |
| 前体类型 | allenamide (全部) |
| 产物已有手性标记 | 是（`[C@]`, `[C@@H]`, `[C@H]`） |
| Mapped SMILES 可用 | 是（`rxn_smiles_mapped` 列完整） |

### DR 分布

```
rx_id=1–4:  无报告 DR（仍应生成 BR_DR_001）
rx_id=5:    major=58,  minor=42
rx_id=6:    major=60,  minor=40
rx_id=7:    major=53,  minor=47
rx_id=8–11: major=96,  minor=4
rx_id=12:   major=95,  minor=5
rx_id=13:   major=87,  minor=13
rx_id=14:   major=83,  minor=17
rx_id=15:   major=70,  minor=30
rx_id=16:   major=87,  minor=13
rx_id=17:   major=93,  minor=7
rx_id=18:   major=52,  minor=48
rx_id=19:   major=60,  minor=40
rx_id=20:   major=62,  minor=38
```

---

## 二、运行前配置

### 2.1 确认 config/defaults.yaml

```yaml
run:
  source: dataset
  resume: true
  resume_rehydrate: true
  resume_rehydrate_policy: best_effort
  output_root: ./rph_output
  workdir_naming: rx_{rx_id}
  dataset:
    format: csv
    path: data/reaxys_cleaned.csv
    product_smiles_col: product_smiles_main
    id_col: rx_id
    delimiter: ","
    small_molecular_primary_col: ylide_leaving_group
    small_molecular_fallback_col: leaving_group

s0:
  enabled: true
  persist_graph: true
  fail_policy: degrade
  generate_intermediates: false
  enable_dr_paths: false
  dr_completion:
    enabled: true                     # ★ 必须为 true
    default_branch_when_unreported: true  # ★ 无报告 DR 时仍生成分支

checkpoint:
  gaussian_oldchk_reuse:
    enabled: false
    allow_cross_theory: false
```

### 2.2 验证导入和测试

```bash
# 1. 导入风格检查
python scripts/ci/check_imports.py rph_core

# 2. 运行全部 DR 相关测试
pytest tests/test_s0_dr_*.py tests/test_dr_aggregator.py \
       tests/test_branch_task_expansion.py \
       tests/test_mechanism_classifier_graph_builder.py \
       tests/test_checkpoint_partial_resume.py -v
```

预期: 全部通过。

---

## 三、预期输出结构

每个独立反应 (rxn_key_hash) 输出如下:

```
rph_output/RXN_xxxxxxxx/
│
├── S0_Mechanism/
│   ├── mechanism_graph.json          # 含 dr_completion 字段
│   ├── mechanism_summary.json        # 含 dr_completion 摘要
│   ├── dr_branch_plan.json           # ★ 唯一 DR 语义源
│   ├── atom_map_smiles.json          # SMILES↔Map 映射
│   ├── s0_status.json
│   └── mechanism_graph.png
│
├── S1_ConfGeneration/
│   └── product/product_min.xyz
│
├── S2_Retro/
│   ├── ts_guess.xyz
│   └── intermediate.xyz
│
├── S3_TS/
│   ├── ts_final.xyz
│   └── sp_matrix_metadata.json
│
├── S4_Data/
│   ├── features_raw.csv              # BR_MAJOR 的特征
│   ├── features_mlr.csv
│   └── feature_meta.json
│
├── branches/                         # ★ DR 分支目录
│   └── BR_DR_001/
│       ├── branch_manifest.json
│       ├── S1_ConfGeneration/
│       ├── S2_Retro/
│       ├── S3_TS/
│       └── S4_Data/
│           ├── features_raw.csv      # BR_DR_001 的特征
│           ├── features_mlr.csv
│           └── feature_meta.json
│
├── reaction_features/
│   ├── geo_electronic_features.csv
│   ├── geo_electronic_features.json
│   └── dr_prediction.json           # ★ ΔΔG‡ → DR 预测
│
├── reaction_manifest.json
│
└── conditions/
    └── COND_xxx/
        ├── condition_manifest.json
        ├── thermo_calculation.json
        └── merged_features.csv
```

---

## 四、执行步骤

### Step 1: 试跑单条反应（rx_id=1 验证 DR 逻辑）

```bash
python -m rph_core --config config/defaults.yaml --rx-id 1
```

**验证点**:

1. `S0_Mechanism/dr_branch_plan.json` 存在且 status 为 `"complete"`:
   ```bash
   cat rph_output/RXN_*/S0_Mechanism/dr_branch_plan.json | python -m json.tool | head -30
   ```
   预期: `status: "complete"`, branches 含 `BR_MAJOR` + `BR_DR_001`

2. `branches/BR_DR_001/branch_manifest.json` 存在且 product_smiles 非空:
   ```bash
   cat rph_output/RXN_*/branches/BR_DR_001/branch_manifest.json | python -m json.tool
   ```

3. `reaction_features/dr_prediction.json` 存在:
   ```bash
   cat rph_output/RXN_*/reaction_features/dr_prediction.json | python -m json.tool
   ```

### Step 2: 试跑有报告 DR 的反应（rx_id=5, major=58/minor=42）

```bash
python -m rph_core --config config/defaults.yaml --rx-id 5
```

**验证点**:

1. `dr_branch_plan.json["reported_dr"]` 应包含:
   ```json
   {"dr_major": "58", "dr_minor": "42"}
   ```

2. `dr_prediction.json` 应计算预测 DR（从 ΔG‡ 差值）。

### Step 3: 全量 20 条批量运行

```bash
# 确保 run 配置中没有 filter_ids，或改为注释
python -m rph_core --config config/defaults.yaml
```

或设置批量参数:

```yaml
run:
  source: dataset
  max_tasks: null          # 不限制
  # filter_ids: 若需要分批，设置具体 rx_id
```

---

## 五、验证清单

### 5.1 S0 DR Plan 验证

对每个反应目录:

```bash
for rxn in rph_output/RXN_*/; do
    plan="$rxn/S0_Mechanism/dr_branch_plan.json"
    if [ -f "$plan" ]; then
        echo "=== $(basename $rxn) ==="
        python3 -c "
import json
with open('$plan') as f:
    p = json.load(f)
branches = p.get('branches', [])
print(f'  status: {p.get(\"status\")}')
print(f'  branches: {len(branches)} ({[b[\"branch_id\"] for b in branches]})')
print(f'  new_stereocenters: {p.get(\"stereocenters\",{}).get(\"new_map_numbers\",[])}')
print(f'  reported_dr: {p.get(\"reported_dr\",{})}')
for b in branches:
    smi = b.get('product_smiles')
    print(f'  {b[\"branch_id\"]}: policy={b[\"generation_policy\"]}, has_smiles={smi is not None}')
"
    else
        echo "=== $(basename $rxn): MISSING dr_branch_plan.json ==="
    fi
done
```

### 5.2 Branch 执行验证

```bash
# 检查所有 BR_DR_001 分支是否完成
for rxn in rph_output/RXN_*/branches/BR_DR_001/; do
    if [ -f "$rxn/S4_Data/features_raw.csv" ]; then
        echo "$(basename $(dirname $(dirname $rxn))): BR_DR_001 OK"
    else
        echo "$(basename $(dirname $(dirname $rxn))): BR_DR_001 MISSING features"
    fi
done
```

### 5.3 DR Prediction 验证

```bash
for rxn in rph_output/RXN_*/reaction_features/dr_prediction.json; do
    rxn_id=$(basename $(dirname $(dirname $rxn)))
    python3 -c "
import json
with open('$rxn') as f:
    p = json.load(f)
pred = p.get('predicted_dr', {})
print(f'{rxn_id}: status={pred.get(\"status\")}, ddg={pred.get(\"ddg_kcal_mol\")}, ratio={pred.get(\"ratio\")}')
" 2>/dev/null || echo "$rxn_id: MISSING or BROKEN"
done
```

---

## 六、预期行为矩阵

### 6.1 rx_id=1–4（无报告 DR）

| 项目 | 预期 |
|------|------|
| `dr_branch_plan.json` status | `"complete"` |
| branches 数量 | 2 (BR_MAJOR + BR_DR_001) |
| `BR_DR_001.generation_policy` | `"concerted_bridgehead_pair_flip"` |
| `BR_DR_001.product_smiles` | 显式手性 SMILES（产物已有 @ 标记，stereo_smiles 可翻转）或降级 |
| `stereocenters.inferred_from_forming_bonds` | `true`（因为 mapped SMILES 中的手性中心可被 RDKit 识别）或 `false` |
| `reported_dr` | `{}` (空) |
| `branches/BR_DR_001/` | 存在且 S4 features 完整 |
| `dr_prediction.json` status | `"complete"` (若 BR_MAJOR 和 BR_DR_001 都有 barrier) |

### 6.2 rx_id=5–20（有报告 DR）

| 项目 | 预期 |
|------|------|
| `reported_dr` | 含 `{"dr_major": "...", "dr_minor": "..."}` |
| branches 数量 | 2 |
| `dr_prediction.json.predicted_dr` | 含 `ddg_kcal_mol` 和 `ratio` |
| 报告 DR vs 预测 DR | 可作为对比验证（见第八节） |

---

## 七、DR Plan 关键字段说明

### `dr_branch_plan.json`

```json
{
  "version": "s0-dr-branch-plan-v1",
  "status": "complete",           // complete | not_applicable | disabled
  "stereocenters": {
    "fixed_map_numbers": [],       // 反应前后不变的手性中心 map 编号
    "new_map_numbers": [1,2,4],    // 新生成手性中心 map 编号
    "inferred_from_forming_bonds": false  // 是否从成键推断
  },
  "branches": [
    {
      "branch_id": "BR_MAJOR",
      "product_smiles": "...",     // 原产物 SMILES
      "generation_policy": "reference_product"
    },
    {
      "branch_id": "BR_DR_001",
      "product_smiles": "...",     // 翻转后 SMILES（或 null 降级）
      "flipped_map_numbers": [1,2],
      "generation_policy": "concerted_bridgehead_pair_flip"
    }
  ],
  "reported_dr": {},              // 数据集原始 DR 值
  "warnings": []
}
```

### `dr_prediction.json`

```json
{
  "version": "rph-dr-prediction-v1",
  "temperature_k": 298.15,
  "branches": [
    {"branch_id": "BR_MAJOR", "delta_g_act_kcal_mol": 12.3},
    {"branch_id": "BR_DR_001", "delta_g_act_kcal_mol": 13.7}
  ],
  "predicted_dr": {
    "major_branch_id": "BR_MAJOR",
    "minor_branch_id": "BR_DR_001",
    "ddg_kcal_mol": 1.4,
    "ratio": 10.6,                // major/minor
    "status": "complete"
  }
}
```

---

## 八、DR 预测 vs 实验对比参考

| rx_id | 实验 DR (major:minor) | 预测 DR | 期望趋势 |
|-------|----------------------|---------|----------|
| 5 | 58:42 (≈1.4:1) | 由 ΔΔG‡ 计算 | ddg ≈ 0–0.5 kcal/mol |
| 6 | 60:40 (1.5:1) | 同上 | ddg ≈ 0.3–0.5 |
| 7 | 53:47 (≈1.1:1) | 同上 | ddg ≈ 0–0.2 (接近平手) |
| 8–11 | 96:4 (24:1) | 同上 | ddg ≈ 1.5–2.0+ kcal/mol |
| 12 | 95:5 (19:1) | 同上 | ddg > 1.5 kcal/mol |
| 18 | 52:48 (≈1.1:1) | 同上 | ddg ≈ 0–0.1 (几乎平手) |

**DR 对比说明**: 当前阶段 ΔG‡ 来自 S3 TS optimization + S4 thermo 提取。若不同分支的 TS barrier 差值太小（< 1 kcal/mol），预测 DR 的不确定性会很高。这是 DFT 精度限制，不是代码问题。建议后续 benchmark 时记录 ΔΔG‡ 误差范围。

---

## 九、故障恢复指南

### 9.1 部分反应失败

```bash
# 只重新运行失败的反应
python -m rph_core --config config/defaults.yaml --rx-id <失败的rx_id>
```

### 9.2 Checkpoint 清除

若需要完全重新计算某反应:

```bash
rm -rf rph_output/RXN_xxxxxxxx/S1_ConfGeneration/
rm -rf rph_output/RXN_xxxxxxxx/S2_Retro/
rm -rf rph_output/RXN_xxxxxxxx/S3_TS/
rm -rf rph_output/RXN_xxxxxxxx/S4_Data/
rm -rf rph_output/RXN_xxxxxxxx/branches/
rm rph_output/RXN_xxxxxxxx/.pipeline.state
rm rph_output/RXN_xxxxxxxx/.pipeline.state.lock
rm rph_output/RXN_xxxxxxxx/.rph_run.lock
```

然后重新运行即可。

### 9.3 DR branch 单独重跑

若主反应已完成但 branch 失败:

```bash
rm -rf rph_output/RXN_xxxxxxxx/branches/BR_DR_001/
# 删除 S0 checkpoint 以触发重新调度
rm rph_output/RXN_xxxxxxxx/.pipeline.state
# 重新运行
python -m rph_core --config config/defaults.yaml --rx-id xxx
```

### 9.4 DR plan 缺失

如果 `dr_branch_plan.json` 不存在:
- 检查 `config/defaults.yaml` 中 `s0.dr_completion.enabled` 是否为 `true`
- 检查 S0 是否执行成功（查看 `s0_status.json`）
- 重启: 删除 `.pipeline.state` 后重新运行

---

## 十、预期运行时间

| 阶段 | 每个反应（估算） | 备注 |
|------|-----------------|------|
| S0 | < 1 sec | 纯 Python/RDKit |
| S1 | 1–4 h | CREST + DFT opt（取决于构象复杂度） |
| S2 | 0.5–2 h | xTB scan |
| S3 | 2–6 h | TS opt + rescue + IRC + SP |
| S4 | < 5 min | 特征提取 |
| BR_DR_001 | ≈ S1+S2+S3+S4 | 与主分支同量级 |
| **总计/反应** | **10–30 h** | 含 DR 分支 |

建议分批运行（如每次 3–4 个独立反应），或使用多进程 `ProcessPoolExecutor`。

---

## 十一、最终验收标准

完成以下全部检查即为成功:

- [ ] 12 个独立反应全部完成 S0（`dr_branch_plan.json` 存在）
- [ ] 12 个独立反应中至少 8 个生成 `branches/BR_DR_001/` 且 S4 features 完整
- [ ] 每个反应的 `reaction_features/dr_prediction.json` 存在
- [ ] `dr_prediction.json["predicted_dr"]["status"]` 为 `"complete"` 或 `"insufficient_data"`（非 error）
- [ ] 所有 `conditions/` 挂载完成（condition_manifest + thermo + merged_features）
- [ ] `python scripts/ci/check_imports.py rph_core` 通过
- [ ] 相关 pytest 全部通过（47 items）
