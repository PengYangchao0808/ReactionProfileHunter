# S3 INT 数据质量问题与后续修改方向报告

> 生成日期：2026-08-08
> 关联修复：ORCA6 NORMAL MODES 解析缺陷（`parse_orca6_normal_mode_vectors`）与
> `parent_structure` 参考几何注入（`_reference_geometry_bundle`）——已完成并验证。
> 本报告聚焦**中间体（INT）极小值正确性**的剩余问题与修改方向。

---

## 1. 背景与结论速览

上一轮 S3 全量重算后（140/140 `usable_for_ml`，20 反应 × 7 结构），独立核查发现：

| 检查项 | 结论 |
|---|---|
| TS 虚频方向 vs 成键方向 | **修复完成**。35/36 primary TS 方向一致（独立重算），仅 rx4 major TS 例外；4 个被救 TS 经质量解权重后方向一致（引擎此前误标 `wrong_mode`） |
| **INT 极小值正确性** | **问题未解决（本次仅修复了检测通路，14 个异常 INT 的处置待定）** |
| INT 分类标签可信度 | `distinct_intermediate` 标签此前全部不可信（`parent_structure` 未注入）；**注入修复后需重算分类验证** |

独立核查（forming-bond 距离 + 产物几何 Kabsch-RMSD + 虚频）发现 **40 个 INT 中 14 个异常**：

| 类型 | 数量 | 判定标准 | 结构清单 |
|---|---|---|---|
| **坍缩到产物**（COLLAPSED） | 5 | 全部 forming-bond 距离 < 1.7 Å（= 产物成键距离）；与产物几何 RMSD < 0.35 Å | rx3/minor、rx7/major、rx7/minor、rx11/minor、rx12/minor |
| **解离/拓扑畸变**（DISSOCIATED） | 9 | 任一 forming-bond 距离 > 4.0 Å（超过解离距离，中间体"散开"） | rx1/major、rx4/minor、rx5/major、rx6/major、rx13/minor、rx14/minor、rx18/minor、rx19/minor、rx20/minor |
| 正常中间体 | 26 | forming-bond 2.6-3.4 Å，无显著虚频 | — |

---

## 2. 根因分析（已完成修复部分）

### 2.1 坍缩/解离检测在生产中从未真正执行（已修复）

**代码缺陷**：`RefinementEngine._classify_structure`（`rph_core/steps/refinement/engine.py` L1131-1139）
调用 `classify_int` 时传入 `precursor_ref` / `product_ref`，二者来自
`_resolve_reference_path(request, ...)`，而后者只读取 `request.parent_structure`。

但 **orchestrator 构建 S3 请求时从未设置 `parent_structure`**（`v4_orchestrator.py`
L997-1061 无此字段），S4 也只是转发 S3 的 `None`。

→ `classify_int` 的 RMSD 坍缩分支被跳过，`avg_progress` 回退默认值 0.5，
**所有 INT 无条件输出 `distinct_intermediate`**。

**修复（已完成）**：
- `v4_orchestrator.py` 新增 `_reference_geometry_bundle()`，为 INT/TS 请求注入
  `parent_structure = {"product_ref": <S1 selected.xyz>, "precursor_ref": <S1 precursor selected.xyz>}`。
- 参考几何选 S1 selected conformer：与 PEB 中间体种子共享原子排序，保证
  `forming_bonds` 索引可直接用于 `_measure_distance` 与 `_compute_mapped_rmsd`。

### 2.2 遗留：修复后需要重算分类才能确认标签

注入修复只打通了**检测通路**，已生成的 40 个 INT manifest 仍是旧的
`distinct_intermediate`。需要：
1. 重跑 S3 classify（纯分类、不重算 QC），或
2. 离线用注入后的 `classify_int` 重算 40 个 INT 的标签。

---

## 3. 14 个异常 INT 的化学含义

### 3.1 坍缩组（5 个）：rx3/minor、rx7/major、rx7/minor、rx11/minor、rx12/minor

- 特征：forming-bond 距离 1.53-1.57 Å（与产物一致），与产物几何 RMSD 0.14-0.35 Å。
- 化学含义：S2 PEB 扫描给出的中间体种子在 B97-3c 优化下**直接滑入产物极小值**——
  该反应路径上不存在稳定中间体，或 PEB 中间体种子与真实势能面极小不匹配。
- rx7（此前 TS 救援过）的 major/minor INT 双双坍缩，与 rx7 的 F2 高阶鞍点问题
  同源——该反应路径整体质量堪忧。

### 3.2 解离组（9 个）：rx1/major、rx4/minor、rx5/major、rx6/major、rx13-14/minor、rx18/minor、rx19-20/minor

- 特征：至少一个 forming-bond 距离 4.2-5.2 Å（如 rx19 minor 达 5.24 Å）。
- 化学含义：中间体"散开"，forming-bond 原子对分离到超出化学键/氢键范围——
  可能是：(a) 该反应无真实中间体，PEB 的 intermediate_seed 是扫描曲线上的
  非驻点；(b) B97-3c 优化逃逸到非化学合理的宽泛区域；(c) 原子排序/索引错位
  （需用 atom_mapping 复核后确认）。
- rx13≈rx14、rx19≈rx20 是重复反应对，解离现象同步出现，去重后独立位点为 7 个。

### 3.3 数据影响

若这 14 个异常 INT 直接进入 ML 特征提取：
- 坍缩组：INT 特征（per-bond distance/progress）与产物几乎相同 → 训练集中出现
  "产物级"中间体样本，污染反应进度标签。
- 解离组：forming-bond 距离 4-5 Å 超出合理反应坐标范围 → 异常特征值，
  且可能与 TS 特征矛盾（TS 的 forming-bond 距离反而更短）。

---

## 4. 后续修改方向（按优先级）

### 方向 A：重算 INT 分类标签（必做，低成本）

用已注入 `parent_structure` 的 `classify_int` 重算 40 个 INT 的
`int_classification`，确认 14 个异常 INT 被正确标记为
`collapsed_to_product` / `collapsed_to_precursor` / `topology_ambiguous`。

- **做法**：重跑 S3（分类阶段）或写离线脚本读 manifest + S1 selected.xyz 重算。
- **注意**：`classify_int` 的 RMSD 依赖 `atom_mapping`（INT→product 映射）。
  需验证 S2 的 `intermediate_atom_mapping.json` 与 INT canonical.xyz 的原子排序一致；
  若不一致需先补映射。

### 方向 B：坍缩组 INT 处置（5 个）

- **B1（推荐）**：接受"无中间体"结论——将 5 个坍缩 INT 从 ML 中间体集剔除，
  在 manifest 中标记 `usable_for_ml=false` + `identity=collapsed_to_product`，
  保留结构供分析。
- **B2**：用约束优化（固定 forming-bond 距离在 PEB 中间体值）重算，确认是否
  存在被绕过的高鞍中间体。成本高、收益不确定，仅在机理研究需要时做。

### 方向 C：解离组 INT 处置（9 个）

- **C1（必做）**：先用 atom_mapping 复核原子排序，排除索引错位假阳性。
- **C2（推荐）**：对复核后仍解离的 INT，检查 PEB `intermediate_seed` 与
  B97-3c 优化轨迹（`S3_LowLevel/<int>/diagnostics/attempt_*_opt/*_trj.xyz`）——
  判断是"种子即散开"还是"优化逃逸"。前者标记无中间体，后者可考虑从
  `intermediate_seed` 用更紧的优化参数（max_cycles、初始 Hessian）重算。
- **C3**：若确认无中间体，标记 `usable_for_ml=false`，保证 ML 数据一致性。

### 方向 D：管线级预防（中长期）

- **D1**：`classify_int` 的坍缩判定已可用（修复后），但建议增加
  **forming-bond 距离硬边界**（如 < 1.7 Å → collapsed、> 4.0 Å → dissociated），
  作为 RMSD/progress 之外的显式守卫——目前 `topology_ambiguous` 语义太宽。
- **D2**：S2 PEB 生成 intermediate_seed 时，若扫描曲线两翼单调（无局域极小），
  直接标记 `has_independent_int=false`（该字段已存在于 S2 manifest），
  避免向 S3 提交注定坍缩/解离的种子。
- **D3**：TS 与 INT 联动质量检查——坍缩组 INT 对应的 TS 若也为 F2/wrong_mode
  （如 rx7），整条路径标记"低置信"，供 ML 数据筛选。

### 方向 E：TS 侧遗留（与 INT 联动）

- rx4 major TS 是唯一方向不一致 TS（alignment 0.002）——建议复核其 S2 种子
  或标记剔除。
- rx6/7/16/17 被救 TS 的 `canonical_frequency_output=None`（calcall 复用 hessian），
  manifest 未记录 freq 输出路径——建议在 `_build_structure_payload` 中补充
  `canonical_frequency_output` 指向 hessian 派生的 freq 结果，保证下游可读。

---

## 5. 建议的验证门（进入 ML 前）

```bash
# 1. 回归（确认修复无副作用）
python -m pytest -q tests/test_orca6_normal_modes.py \
  tests/test_identity.py tests/test_refinement_pass1_primary.py \
  tests/test_refinement_pass2_rescue.py tests/test_refinement_rescue_only.py

# 2. 重算 INT 分类（方向 A 落地后）——验证 14 个异常 INT 标签正确
#    预期：5 个 collapsed_to_product，9 个 topology_ambiguous/dissociated

# 3. 独立复核抽查（方向 C1）——确认解离组非原子索引错位
#    对比 S2 intermediate.xyz 与 canonical.xyz 的 forming-bond 距离

# 4. 最终统计（去重 rx9/10、13/14、16/17、19/20 后）：
#    正常 INT >= 26，坍缩/解离标记 <= 14，全部 usable_for_ml 与 identity 一致
```

---

## 6. 附录：本次已完成的代码修复

| 文件 | 变更 |
|---|---|
| `rph_core/utils/orca_interface.py` | 新增 `parse_orca6_normal_mode_vectors()`：解析 ORCA6 `NORMAL MODES` 块（多 6 列分块 + 单列块），mass-weighted → 按 `sqrt(m)` 解权重并归一化，返回与 ORCA5 `_parse_orca_displacement_vectors` 语义一致的笛卡尔位移向量 |
| `rph_core/steps/refinement/engine.py` | `_build_ts_normal_mode_payload`：freq 路径优先 ORCA5 解析、失败回退 ORCA6 解析；hessian 路径新增 `_demass_weight_mode_vectors()` 解权重；import `_ELEMENT_MASS`、`parse_orca6_normal_mode_vectors` |
| `rph_core/v4_orchestrator.py` | 新增 `_reference_geometry_bundle()`；INT 请求注入 `parent_structure`（product_ref + precursor_ref = S1 selected.xyz） |
| `tests/test_orca6_normal_modes.py` | 新增 6 个测试：解析形状/归一化、sqrt(mass) 解权重数学、行数不匹配拒绝、解权重归一化、未知元素回退、形状不匹配返回 None |

**验证结果**：py_compile ✓；82 passed / 5 xfailed / 2 xpassed（xfail/xpass 为既有标记）✓；
check_imports PASSED ✓；端到端重算 36 个 primary TS → 23 `valid_target_ts` /
13 `first_order_wrong_mode`（alignment 0.19-0.29 边界）/ 1 例外（rx4 major）。
