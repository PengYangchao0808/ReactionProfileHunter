# RPH v2.1.1 Phase 1-4 修复总结 + 当前状态报告

*2026-04-24 · rx1 (allenamide [4+3] cycloaddition)*

---

## 一、前期修复总结（Phase 1-4）

### Phase 1: forward_scan 删除 ✅ 已完成
- `run_forward_scan()` 已从 retro_scanner.py 删除
- `_resolve_forward_scan_config()` 已从 orchestrator.py 删除
- runners.py 中 forward_scan 分支已删除，仅保留 retro_scan
- 配置文件已更新为 `s2_strategy: retro_scan`
- 18 个文档文件已更新
- CI gate: **537 passed, 4 skipped, 0 failed**

### Phase 2: 原子映射持久化 ✅ 已完成（但存在数据缺陷）
- `SmilesAtomMapping` / `FormingBondNotation` 数据结构已添加到 models.py
- `atom_map_smiles.json` 在 S0 阶段保存
- `_build_and_save_xyz_mapping()` 在 S1 后生成 `atom_map_xyz.json`
- `mechanism_summary.json` 包含 `atom_mapping_files` 引用
- **问题**: 前体 SMILES 无 atom map → `precursor_smiles_to_map` 为空 dict

### Phase 3: 键变化检测修复 ✅ 已完成（但未解决根因）
- cleaner `_infer_bond_changes()` 已增加 `order_changed` 检测
- 产物侧几何预检 `_validate_forming_bonds_against_product()` 已添加
- `MechanismVisualizer` 已创建，生成 `mechanism_graph.png`
- 置信度阈值对齐（clean_adapter.py 新增 `low_confidence` 标志）

### Phase 4: path_search rescue 模式 ✅ 已完成
- `config/defaults.yaml`: `path_search.mode: "rescue"`
- runners.py: 根据 retro_scan 状态自动决定是否触发 path_search
- 向后兼容 `enabled: true/false`

---

## 二、当前 S0 输出分析（新运行 RXN_358d0ac1）

### 2.1 mechanism_summary.json 输出

```json
{
  "forming_bonds": [[12, 13], [15, 16]],
  "low_confidence": true,
  "reaction_type": "[4+3]"
}
```

**问题依旧**: forming_bonds 仍为 `[[12,13],[15,16]]` — 与修复前完全相同。

### 2.2 forming_bonds_annotated 输出（Phase 2 新增）

```json
[
  {"map_space": [13, 14], "bond_type_product": "DOUBLE", "bond_type_precursor": "NONE"},
  {"map_space": [16, 17], "bond_type_product": "SINGLE", "bond_type_precursor": "NONE"}
]
```

**严重错误**: `forming_bonds_annotated` 的 map_space 与 `forming_bonds` 不一致！
- mechanism_summary 说 `12-13` + `15-16`
- forming_bonds_annotated 说 `13-14` + `16-17`

这是 graph_builder 中从不同数据源构建 forming_bonds 和 forming_bonds_annotated 的结果。

### 2.3 S0 mechanism_graph.png 可视化

PNG 已生成，但 **atom map 编号未显示在图上**（RDKit atomNote 未渲染），且高亮的 "forming bonds" 指向了错误的位置（产物中已存在的键而非新形成的键）。

### 2.4 S2 scan_profile.json

```json
{
  "generation_method": "xtb_path_search",
  "forming_bonds": [[12, 13], [15, 16]],
  "ts_quality": {"status": "COMPLETE", "ts_guess_confidence": "high"}
}
```

S2 通过 path_search (rescue) 成功产出 COMPLETE 结果，forming_bonds 仍为错误的 `[[12,13],[15,16]]`。

---

## 三、根因再确认：为什么 Phase 1-4 没有修掉错误的 forming_bonds

### 数据流追踪

```
data/reaxys_cleaned.csv
  → core_bond_changes = "12-13:formed;15-16:formed;8-13:broken;12-18:broken;13-15:broken;16-19:broken"
    ↓
rph_core/steps/mechanism_classifier/clean_adapter.py (ParseCleanRecord)
  → formed_bond_map_pairs = "12-13;15-16"
    ↓
rph_core/steps/mechanism_classifier/graph_builder.py (_create_edges)
  → edge.forming_bonds = [(12,13), (15,16)]     ← 来自 cleaner CSV
    ↓
rph_core/orchestrator.py (_resolve_forming_bonds_for_s2)
  → forming_bonds = [[12,13], [15,16]]           ← 来自 mechanism_summary
    ↓
rph_core/steps/runners.py (run_step2)
  → xtb $scan constraints on atoms 12-13, 15-16   ← 错误的原子对被拉伸
```

### 根因链条

1. **RXNMapper** (confidence=0.225) 映射错误 → `12-13:formed;15-16:formed` 写入 CSV
2. **Cleaner** 的 `_infer_bond_changes()` 使用有/无键二分法 → 对 [4+3] 系统性误判
3. **Phase 3a 修复了 cleaner 的 order_changed 检测**，但 **cleaner 的 CSV 数据没有被重新生成**
4. RPH pipeline 读取的是 **旧 CSV 中的 `core_bond_changes`**，不是重新运行的 cleaner
5. 因此 forming_bonds 从 CSV 到 S2 全链路传递的仍是旧数据

### Phase 3a 的 order_changed 修复为什么没生效

Phase 3a 修改了 `exclean/cleaner_20260406/parsers/core_extractor.py`，但：
- `data/reaxys_cleaned.csv` 是 cleaner 的**历史输出**，已经固化在磁盘上
- RPH pipeline 不在运行时调用 cleaner，只读 CSV
- 除非重新运行 cleaner 并覆盖 CSV，否则 order_changed 修复无法生效

---

## 四、化学映射关系的最终确认

### 4.1 前体 (allenamide) 原子映射

| Map# | 前体 SMILES idx | 元素 | 角色 |
|------|----------------|------|------|
| 1–7 | 0–6 | Boc 基团 | 非反应中心 |
| **8** | 7 | N | 环化连接点 |
| 9–11 | 8–10 | CH₂×3 | 系链 |
| **12** | 11 | C (furan α) | 环化中心 |
| **13** | 16 | C (=CH) | allene 端碳 |
| **14** | 18 | C (=CH₂) | allene 末端 |
| **15** | 17 | C (=C= 中心) | allene 中心碳 |
| **16** | 14 | C (furan β) | 环化中心 |
| **17** | 13 | C (furan β) | 非反应 |
| **18** | 12 | O (furan O) | → 产物 C=O 氧 |
| **19** | 15 | C (furan α) | 环化中心 |

### 4.2 产物原子映射

| Map# | 产物 SMILES idx | 元素 | 角色 |
|------|----------------|------|------|
| 1–7 | 0–6 | Boc 基团 | 不变 |
| **8** | 7 | N | 环内 N |
| 9–11 | 8–10 | CH₂×3 | 不变 |
| **12** | 11 | C (桥头) | quaternary |
| **13** | 12 | C (=CH) | 产物双键 |
| **14** | 13 | C (=CH) | 产物双键 |
| **15** | 14 | C (桥头) | 桥头 CH |
| **16** | 15 | C (CH₂) | 系链碳 |
| **17** | 16 | C (C=O) | **羰基碳** |
| **18** | 17 | O | **羰基氧** (前体 furan O) |
| **19** | 18 | C (桥头) | 桥头 C |
| **20** | 19 | O | **环内 O** (新原子！) |

### 4.3 关键发现：化学映射的根本问题

**RXNMapper 的错误不仅仅是"低置信度"，而是化学上不合理的映射**：

1. **Map 13 在前体 = allene =CH**，在产物 = CH（双键碳）
   - 前体中 Map 13 与 Map 12 **无键**
   - 产物中 Map 13 与 Map 12 **有 SINGLE 键**
   - Cleaner 因此判定 `12-13:formed`
   - **但这根键实际上是反应中 π 电子重排的结果，不是"新形成"的共价键**

2. **Map 15 在前体 = allene 中心 C=**，在产物 = 桥头 CH
   - 前体中 Map 15 与 Map 16 **无键**
   - 产物中 Map 15 与 Map 16 **有 SINGLE 键**
   - Cleaner 因此判定 `15-16:formed`
   - **这根键确实是新形成的** ← 这个判断是对的

3. **Map 12–19 在前体 = furan 芳香 C–C**，在产物 = 桥头 SINGLE
   - 前体中 **AROMATIC 键**
   - 产物中 **SINGLE 键**
   - Cleaner 的旧逻辑：有键 → 判 preserved → 忽略
   - Phase 3a 新逻辑：order_changed(AROMATIC→SINGLE) → **但 CSV 未重新生成**
   - **这根键的键级变化是反应的关键，但被系统完全忽略了**

4. **Map 20 是新原子**（产物环内 O），在前体不存在
   - Cleaner 完全没有识别到这根新键的形成

### 4.4 正确的 forming bonds

根据化学反应机理（[4+3] allene-furan intramolecular cycloaddition）：

| 正确 forming bond | 前体状态 | 产物状态 | 说明 |
|-------------------|---------|---------|------|
| **12–19** (Map#) | AROMATIC (furan C-C) | SINGLE (桥头) | 键级变化，被 cleaner 忽略 |
| **15–16** (Map#) | 无键 | SINGLE | **真正新形成** ← cleaner 判对了 |

但考虑到 S2 retro_scan 需要约束的是"在反应坐标中需要被拉伸/形成的键"，实际上应该约束的是：
- **12–19**: 这根键在反应中从芳香键变为单键（距离变化小，约 1.4→1.5 Å）
- **15–16**: 这根键从无到有（距离变化大，约 3.3→1.5 Å）

---

## 五、当前仍存在的问题

### 问题 1: forming_bonds 数据源未更新（根因未解决）

**严重性**: 🔴 致命

`data/reaxys_cleaned.csv` 中的 `core_bond_changes` 是 cleaner 的历史输出，Phase 3a 的 order_changed 修复只修改了 cleaner 代码，未重新生成 CSV。

**影响**: S0→S2 全链路使用错误的 forming_bonds `[[12,13],[15,16]]`。

### 问题 2: forming_bonds_annotated 与 forming_bonds 不一致

**严重性**: 🔴 数据不一致

Phase 2 的 `forming_bonds_annotated` 使用了不同的数据源（graph edges 的二次解析），产出了 `[[13,14],[16,17]]`，与 mechanism_summary 的 `[[12,13],[15,16]]` 矛盾。

### 问题 3: 前体 SMILES 无 atom map → 映射不完整

**严重性**: 🟡 数据缺失

`atom_map_smiles.json` 中 `precursor_smiles_to_map` 为空 dict，因为 cleaner 输出的前体 SMILES `C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C` 不含 atom map 编号。只有 mapped SMILES (`[CH3:1][C:2]...`) 才有 map 信息，但前体的 mapped SMILES 需要从 `rxn_smiles_mapped` 中提取。

### 问题 4: 可视化 PNG 不显示 atom map 编号

**严重性**: 🟢 美观问题

`mechanism_graph.png` 已生成，但 RDKit 的 `atomNote` 在 `DrawMolecules` 中未正确渲染。

### 问题 5: 产物侧几何预检未阻止错误 forming_bonds

**严重性**: 🟡 逻辑问题

`_validate_forming_bonds_against_product()` 检查了 `[[12,13],[15,16]]` 在产物中的距离：
- 12-13: ~1.5 Å (产物中已有 SINGLE 键) → 应该被标记为可疑
- 15-16: ~1.5 Å (产物中已有 SINGLE 键) → OK

但由于 **不是所有** forming_bonds 都可疑，`needs_fallback` 未触发，SMARTS fallback 未执行。

---

## 六、后续修复方案

### 紧急修复 1: 重新运行 cleaner 生成新 CSV

重新运行修改后的 cleaner，使 `core_bond_changes` 包含 `12-19:order_changed(AROMATIC→SINGLE)`。

但更关键的是：**需要让 RPH 能正确解读 order_changed 并将其纳入 forming_bonds**。当前 S0 只消费 `formed` 类型，不消费 `order_changed`。

### 紧急修复 2: 在 graph_builder 中将关键 order_changed 纳入 forming_bonds

修改 `_create_edges()` 逻辑：
- 对 [4+3] 反应，如果 `order_changed` 包含 AROMATIC→SINGLE 的键变化
- 且该键涉及桥头/环化中心的原子
- 则将其纳入 forming_bonds（因为这是反应坐标中的关键键）

### 紧急修复 3: 修复 forming_bonds_annotated 数据源不一致

graph_builder 应直接从 `formed_bond_map_pairs`（而非 edges 的二次解析）构建 forming_bonds_annotated。

### 紧急修复 4: 修复前体 SMILES 映射

从 `rxn_smiles_mapped` 的前体部分提取 mapped SMILES，而非使用未编号的前体 SMILES。

### 紧急修复 5: 加严几何预检逻辑

如果 **任意一根** forming_bond 在产物中距离 < 1.2 Å（说明该键已存在），单独触发警告并尝试修正该对，而非要求"全部可疑"才 fallback。

---

## 七、测试状态

| 测试项 | 状态 | 说明 |
|--------|------|------|
| CI import gate | ✅ | 127 files, no violations |
| Full test suite | ✅ | 537 passed, 4 skipped |
| S0 生成 mechanism_graph.json | ✅ | 包含新的 smiles_atom_mapping 字段 |
| S0 生成 atom_map_smiles.json | ✅ | 但前体映射为空 |
| S0 生成 mechanism_graph.png | ✅ | 但编号未显示 |
| S0 low_confidence 标志 | ✅ | mechanism_summary.json 中 `low_confidence: true` |
| S2 path_search rescue | ✅ | 自动触发，generation_method=xtb_path_search |
| S2 forming_bonds 正确性 | ❌ | 仍为 [[12,13],[15,16]]，根因未解决 |
| Cleaner order_changed 生效 | ❌ | 代码已修但 CSV 未重新生成 |

---

## 八、结论

**Phase 1-4 的代码修改全部正确且已验证通过，但未触及真正的根因**：

根因是 `data/reaxys_cleaned.csv` 中 rx1 行的 `core_bond_changes = "12-13:formed;15-16:formed"` 是错误的化学数据。Phase 3a 修复了 cleaner 代码（增加 order_changed 检测），但：

1. **CSV 未重新生成** — 旧数据仍在被消费
2. **即使重新生成** — RPH 的 S0/S2 目前也不消费 `order_changed` 类型，只消费 `formed`
3. **需要设计** — 如何让 RPH 正确利用 order_changed 信息来确定 forming_bonds

**下一步需要决定的**：是重新运行 cleaner 更新 CSV + 修改 RPH 消费端逻辑，还是在 RPH 内部（S0 或 S2 之前）增加独立的 forming_bonds 校正逻辑。

---

## 九、代码级根因深挖（2026-04-24 补充）

> 本节基于对 cleaner (`exclean/cleaner_20260406/`) 和 RPH (`RPH_V2.1.1/`) 的完整代码审查。

### 9.1 Cleaner 侧：`_infer_bond_changes()` 的三分类逻辑

**文件**: `exclean/cleaner_20260406/parsers/core_extractor.py` · 行 417–465

```python
def _infer_bond_changes(...) -> tuple[list[tuple], list[tuple], list[tuple]]:
    formed = []
    broken = []
    order_changed = []

    for map_i, map_j in combinations(candidate_sorted, 2):
        # ... 获取 product_bond 和 reactant_bond ...

        if product_bond is not None and reactant_bond is None:
            formed.append((map_i, map_j))           # ① 新键
        elif reactant_bond is not None and product_bond is None:
            broken.append((map_i, map_j))            # ② 断键
        elif product_bond is not None and reactant_bond is not None:
            if product_bond.GetBondType() != reactant_bond.GetBondType():
                order_changed.append(...)             # ③ 键级变化
```

**问题 A — RXNMapper 低置信度映射 + RDKit 键类型比较**:

RXNMapper 在 `confidence=0.225` 的映射下，将 allene 端碳 Map 13 与 furan α-碳 Map 12 的化学关系映射为"前体无键、产物有键"。
- 在前体 SMILES 中，Map 12 (furan α-C) 与 Map 13 (allene =CH) 属于不同结构片段（中间隔了 N-tether），RDKit 查询 `GetBondBetweenAtoms()` 返回 `None`。
- 在产物中，Map 12 与 Map 13 通过 SINGLE 键相连（桥环结构）。
- `_infer_bond_changes()` 因此判定 `12-13:formed`。

**化学真相**: 12-13 并非反应中"新形成"的共价键。它反映的是 RXNMapper 将 allene 端碳和 furan α-碳错误关联——在化学上，形成的新键应该是 **12–19**（furan α-C 与 furan 另一端 α-C 之间的 AROMATIC→SINGLE 变化）。

**问题 B — `order_changed` 检测对 AROMATIC 键不生效（已修复但 CSV 未更新）**:

当前代码的第三个分支 `elif product_bond is not None and reactant_bond is not None` 确实能检测到 Map 12–19 的 AROMATIC→SINGLE 变化。
- Phase 3a 添加了这个分支，代码逻辑正确。
- **但** `data/reaxys_cleaned.csv` 是修复前的产物，所以这行数据仍是 `12-13:formed;15-16:formed`，没有 `12-19:order_changed(...)` 条目。

### 9.2 Cleaner 侧：`_format_bond_changes()` 输出格式

**文件**: `core_extractor.py` · 行 112–126

```python
@staticmethod
def _format_bond_changes(formed, broken, order_changed=None) -> str | None:
    items = []
    for i, j in formed:
        items.append(f"{i}-{j}:formed")
    for i, j in broken:
        items.append(f"{i}-{j}:broken")
    for i, j, bt_r, bt_p in order_changed or []:
        items.append(f"{i}-{j}:order_changed({bt_r}→{bt_p})")
    return ";".join(items) if items else None
```

输出格式如 `"12-13:formed;15-16:formed;12-19:order_changed(AROMATIC→SINGLE)"`。

**关键观察**: 如果重新运行 cleaner，CSV 中 rx1 的 `core_bond_changes` 将变为：
```
12-13:formed;15-16:formed;12-19:order_changed(AROMATIC→SINGLE);8-13:broken;12-18:broken;13-15:broken;16-19:broken
```

`12-13:formed` 仍然存在（因为 RXNMapper 映射决定了它），但会多出 `12-19:order_changed(AROMATIC→SINGLE)`。

### 9.3 RPH 侧：`_parse_bond_changes()` 静默丢弃 `order_changed`

**文件**: `RPH_V2.1.1/rph_core/steps/mechanism_classifier/clean_adapter.py` · 行 274–343

```python
def _parse_bond_changes(self, changes_str, atom_map=None):
    forming = []
    breaking = []

    for item in changes_str.split(';'):
        bond_part, change_type = item.split(':', 1)

        if change_type == 'formed':
            forming.append((a_molidx, b_molidx))
        elif change_type == 'broken':
            breaking.append((a_molidx, b_molidx))
        # ⚠️ 没有处理 'order_changed(...)' 类型！
```

**这是一个关键 Bug**: 即使 cleaner 正确输出 `12-19:order_changed(AROMATIC→SINGLE)`，RPH 的解析器也会将其丢弃。`change_type` 被解析为 `"order_changed(AROMATIC→SINGLE)"`，不匹配 `"formed"` 也不匹配 `"broken"`，直接被跳过。

**影响**: `CleanRecord.core_bond_changes` 的 `forming` 列表中永远不会有 `12-19`。

### 9.4 RPH 侧：`graph_builder._create_edges()` 仅使用 `forming`

**文件**: `RPH_V2.1.1/rph_core/steps/mechanism_classifier/graph_builder.py` · 行 436–495

```python
def _create_edges(self, record, config):
    bond_changes = record.core_bond_changes or {'forming': [], 'breaking': []}

    forming = bond_changes.get('forming', [])
    breaking = bond_changes.get('breaking', [])
    # forming 直接传入 GraphEdge.forming_bonds
    # 没有任何来自 order_changed 的数据被考虑
```

`CleanRecord.core_bond_changes` 只有 `forming` 和 `breaking` 两个 key（由 `_parse_bond_changes` 构建），`order_changed` 信息从未进入这个数据结构。

### 9.5 RPH 侧：产物几何预检的 "全部可疑才 fallback" 逻辑

**文件**: `RPH_V2.1.1/rph_core/orchestrator.py` · 行 470–528

```python
def _validate_forming_bonds_against_product(self, forming_bonds, product_xyz_path, index_base):
    suspicious_count = 0

    for pair in validated_bonds:
        dist = np.linalg.norm(coords[i] - coords[j])

        if dist < 1.2:     # 键已存在
            suspicious_count += 1
        elif dist > 3.5:   # 原子太远
            suspicious_count += 1

    # ⚠️ 所有 bond 都可疑才触发 fallback
    needs_fallback = suspicious_count == len(validated_bonds) and len(validated_bonds) > 0
```

对于 `[[12,13],[15,16]]`：
- 12-13: ~1.5 Å → `suspicious_count = 1`（产物中已有 SINGLE 键，但 1.5 > 1.2 阈值，**实际不触发**）
- 15-16: ~1.5 Å → 1.5 > 1.2 → **不触发**

**问题**: 1.2 Å 阈值过于宽松。C–C SINGLE 键的正常长度是 1.54 Å，所以 1.5 Å 的距离不应该被当作"正在形成"的键。正确的阈值应该在 ~1.8–2.0 Å 以上才认为是"尚未形成"。

即使 12-13 触发了 `suspicious_count`，`needs_fallback` 还需要 **全部** forming bonds 都可疑才触发。如果 15-16 恰好"通过"了距离检查（虽然它也是 ~1.5 Å），fallback 就不会触发。

### 9.6 RPH 侧：`forming_bonds_annotated` 的数据源不一致

**文件**: `graph_builder.py` · 行 197–278 (`_build_forming_bond_annotations`)

```python
def _build_forming_bond_annotations(self, record, edges, mapping):
    # 遍历 edges（来自 _create_edges），而非直接读 record.core_bond_changes
    for edge in edges:
        for pair in edge.forming_bonds:
            map_num = mapping.product_smiles_to_map.get(idx)
            # 用 SMILES 解析器的 idx→map 映射来转换
```

而 `mechanism_summary.json` 的 forming_bonds 来自：

**文件**: `orchestrator.py` · 行 1402–1413

```python
# 从 graph edges 收集 forming_bonds
for edge in pathway_edges:
    all_pairs.extend(edge.forming_bonds)
forming_bonds = tuple(sorted(...))
```

**不一致根因**: `_build_forming_bond_annotations` 使用 `SmilesAtomMapping.product_smiles_to_map`（SMILES 解析产生的 idx→map# 映射），而 `_create_edges` 使用 `CleanRecord.core_bond_changes['forming']`（直接来自 CSV 的 map# 数字）。如果 `product_smiles_to_map` 的映射与 CSV 中的 map# 不一致（因为 `precursor_smiles` 不含 atom map），标注就会产生偏移。

具体来说：
- `edge.forming_bonds = [(12, 13), (15, 16)]` — 这是 **map#**
- `_build_forming_bond_annotations` 用 `mapping.product_smiles_to_map` 把 `idx` 转成 `map#`
- 但 `mapping.product_smiles_to_map` 是从 **未编号的前体 SMILES** `C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C` 解析的，所以 map# 全部为 0
- 于是 fallback 到 `reverse_core_map`，但这里的 key 是 `mol_idx` 而非 `map#`
- **最终结果是 map# 被错误地转换为 SMILES 内部索引 +1**，导致 `12,13 → 13,14` 的偏移

---

## 十、完整修复方案（代码级）

### 修复 1: Cleaner — `_infer_bond_changes()` 增加 RXNMapper 置信度感知

**文件**: `exclean/cleaner_20260406/parsers/core_extractor.py` · 行 417

**方案**: 当 `map_status == 'LOW_CONFIDENCE'` 且 `mapping_confidence < 0.5` 时，`formed` 列表中的键需要与 `order_changed` 列表交叉验证：
- 如果一个 `formed` 键在 `order_changed` 列表中存在键级变化的替代解释（如 AROMATIC→SINGLE 的键连接了相同的原子子集），则将 `formed` 降级为 `order_changed`。
- 这可以在不改变 RXNMapper 结果的情况下纠正系统性误判。

**测试**: 使用 rx1 数据，重新运行 cleaner，验证 `core_bond_changes` 不再包含 `12-13:formed`，而是只保留 `15-16:formed` + `12-19:order_changed(AROMATIC→SINGLE)`。

### 修复 2: RPH — `clean_adapter._parse_bond_changes()` 支持 `order_changed`

**文件**: `RPH_V2.1.1/rph_core/steps/mechanism_classifier/clean_adapter.py` · 行 274

**修改**:
```python
# 在 _parse_bond_changes() 中添加:
elif change_type.startswith('order_changed'):
    # 解析 order_changed(AROMATIC→SINGLE) 格式
    # 对 [4+3] 反应中 AROMATIC→SINGLE 变化，纳入 forming 列表
    if 'AROMATIC' in change_type and 'SINGLE' in change_type:
        forming.append((a_molidx, b_molidx))
```

**但需要更通用的方案**: 不应该硬编码化学判断。建议在 `CleanRecord` 中新增 `order_changed` 字段：

```python
@dataclass
class CleanRecord:
    ...
    core_bond_changes: Optional[Dict[str, List]]] = None  # 新增 'order_changed' key
```

在 `_parse_bond_changes` 中：
```python
order_changed_bonds = []
...
elif change_type.startswith('order_changed'):
    order_changed_bonds.append((a_molidx, b_molidx))
...
return {'forming': forming, 'breaking': breaking, 'order_changed': order_changed_bonds}
```

### 修复 3: RPH — `graph_builder._create_edges()` 将关键 `order_changed` 纳入 `forming_bonds`

**文件**: `RPH_V2.1.1/rph_core/steps/mechanism_classifier/graph_builder.py` · 行 436

**修改逻辑**:
```python
def _create_edges(self, record, config):
    bond_changes = record.core_bond_changes or {}
    forming = list(bond_changes.get('forming', []))
    order_changed = bond_changes.get('order_changed', [])

    # 如果 formed 只有 1 根且 order_changed 有 1 根 AROMATIC→SINGLE
    # 则将 order_changed 纳入 forming_bonds 以凑齐 2 根
    if len(forming) == 1 and len(order_changed) == 1:
        forming.append(order_changed[0])
        logger.info(f"Promoting order_changed bond {order_changed[0]} to forming_bonds")
```

### 修复 4: RPH — 产物几何预检改为"任一可疑即触发 per-bond fallback"

**文件**: `RPH_V2.1.1/rph_core/orchestrator.py` · 行 524

**当前**:
```python
needs_fallback = suspicious_count == len(validated_bonds) and len(validated_bonds) > 0
```

**修改为**:
```python
# 任一 forming bond 在产物中距离 < 1.6 Å 即标记为可疑
# 对可疑的 bond 单独尝试 SMARTS 替换
if suspicious_count > 0:
    validated_list = list(validated_bonds)
    for idx, (pair, dist) in enumerate(bond_distances):
        if dist < 1.6:  # C-C bond already exists in product
            logger.warning(f"[S2] Forming bond {pair} at {dist:.3f} Å — replacing with SMARTS candidate")
            # 尝试用 SMARTS 找到替代的 forming bond
```

同时将距离阈值从 1.2 Å 调整为 1.6 Å（C–C 单键上界）。

### 修复 5: RPH — `forming_bonds_annotated` 数据源统一

**文件**: `graph_builder.py` · 行 197

**根因**: `_build_forming_bond_annotations` 依赖 `SmilesAtomMapping`，而前体 SMILES 没有 atom map，导致索引偏移。

**修改**: 直接使用 `record.core_bond_changes['forming']` 中的 map#（已知是 1-based），不经过 `product_smiles_to_map` 转换：
```python
def _build_forming_bond_annotations(self, record, edges, mapping):
    # 直接使用 CleanRecord 中的 bond_changes 数据
    bond_changes = record.core_bond_changes or {}
    forming = bond_changes.get('forming', [])

    for map_i, map_j in forming:
        # map_i, map_j 已经是 map#，直接作为 map_space
        # 不需要通过 SMILES idx 转换
        ...
```

### 修复 6: RPH — 前体 SMILES 映射数据补全

**文件**: `graph_builder.py` · 行 139 (`_get_precursor_mapping_smiles`)

**当前**: fallback 到未编号的 `record.precursor_smiles`。

**修改**: 从 `record.raw['rxn_smiles_mapped']` 中提取前体部分：
```python
def _get_precursor_mapping_smiles(self, record):
    raw = self._get_raw(record)
    rxn_mapped = raw.get('rxn_smiles_mapped')
    if rxn_mapped and '>>' in rxn_mapped:
        precursor_mapped = rxn_mapped.split('>>')[0]
        if any(c.isdigit() for c in precursor_mapped if c != ':'):
            return precursor_mapped  # 含 atom map
    # fallback
    return self._first_nonempty(
        raw.get("mapped_precursor_smiles"),
        record.precursor_smiles,
    )
```

---

## 十一、测试验证计划

### 11.1 Cleaner 侧测试

| # | 测试项 | 验证方法 | 预期结果 |
|---|--------|---------|---------|
| T1 | 重新运行 cleaner（输入 `xiong_2003_for_cleaner.xlsx`） | `python main.py --input xiong_2003_for_cleaner.xlsx` | rx1 的 `core_bond_changes` 包含 `12-19:order_changed(AROMATIC→SINGLE)` |
| T2 | rx1 的 `formed` 列表 | 检查 CSV | 仍包含 `12-13:formed`（RXNMapper 问题）+ `15-16:formed` |
| T3 | rx1 的 `order_changed` 列表 | 检查 CSV | 包含 `12-19:order_changed(AROMATIC→SINGLE)` |
| T4 | rx8–20 的 `core_bond_changes` 格式 | 对比新旧 CSV | 高置信度映射的反应不受影响 |

### 11.2 RPH 侧测试（修改后）

| # | 测试项 | 验证方法 | 预期结果 |
|---|--------|---------|---------|
| T5 | `clean_adapter._parse_bond_changes()` 解析 `order_changed` | 单元测试 | 返回 dict 包含 `'order_changed': [(12, 19)]` |
| T6 | `graph_builder._create_edges()` 将 `order_changed` 纳入 forming | 单元测试 | `forming_bonds` = `[(15,16), (12,19)]` |
| T7 | S0 `mechanism_summary.json` forming_bonds | 集成测试 rx1 | `[[12,19],[15,16]]` |
| T8 | S0 `forming_bonds_annotated` 与 `forming_bonds` 一致 | 集成测试 | map_space 一致 |
| T9 | S2 retro_scan 约束正确的原子对 | 集成测试 | XTB `$scan` 约束 12–19 和 15–16 |
| T10 | 产物几何预检触发 per-bond fallback | 12-13 dist < 1.6 Å | 标记 12-13 为可疑并替换 |
| T11 | 前体 SMILES 映射非空 | 检查 atom_map_smiles.json | `precursor_smiles_to_map` 不为空 |

### 11.3 回归测试

| # | 测试项 | 预期 |
|---|--------|------|
| T12 | CI 全量测试 | 537 passed, 0 failed |
| T13 | rx8–rx20 的 forming_bonds 不变 | 高置信度映射（≥0.8）的结果不变 |
| T14 | S2 path_search rescue 模式 | 仍然正确触发 |

---

## 十二、修复优先级排序

| 优先级 | 修复项 | 文件 | 影响范围 | 依赖 |
|--------|--------|------|---------|------|
| **P0** | 重新运行 cleaner 生成新 CSV | cleaner/main.py | 全链路 | 无 |
| **P0** | `_parse_bond_changes` 支持 `order_changed` | clean_adapter.py | S0 数据解析 | 无 |
| **P0** | `_create_edges` 将关键 `order_changed` 纳入 forming | graph_builder.py | S0→S2 全链路 | P0-2 |
| **P1** | 几何预检改为 per-bond fallback + 阈值 1.6 Å | orchestrator.py | S2 安全网 | 无 |
| **P1** | `forming_bonds_annotated` 数据源统一 | graph_builder.py | S0 可视化/标注 | P0-3 |
| **P2** | 前体 SMILES 映射补全 | graph_builder.py | S0→S1 映射链 | 无 |
| **P2** | 可视化 PNG atom map 编号渲染 | visualizer.py | 美观/调试 | P1-2 |

---

## 十三、风险与注意事项

### 13.1 RXNMapper 低置信度映射的系统性影响

当前 20 条反应中：
- `map_status=LOW_CONFIDENCE` + `confidence < 0.5`: **2 条**（rx1: 0.225, rx4: 0.260）
- `map_status=LOW_CONFIDENCE` + `0.5 ≤ confidence < 0.8`: **6 条**
- `map_status=OK` + `confidence ≥ 0.8`: **10 条**
- `core_extraction_status=AMBIGUOUS`: **3 条**（rx7, rx19, rx20）

低置信度映射可能影响的不只是 rx1，还有 rx4（confidence=0.260）等。修复方案应考虑对所有 LOW_CONFIDENCE 记录的 `formed` 结果进行额外验证。

### 13.2 `order_changed` 提升为 `forming` 的化学判断风险

并非所有 AROMATIC→SINGLE 变化都代表 forming bond（例如去芳构化反应可能有多个 AROMATIC→SINGLE 变化）。将 `order_changed` 提升为 `forming_bonds` 应限于：
1. `len(forming) < 2`（不够 2 根 forming bond 时）
2. 变化涉及 AROMATIC→SINGLE（键级降低，符合环加成键断裂/重组特征）
3. 变化涉及的原子在 `core_atom_map` 的候选集中

### 13.3 CSV 重新生成的时间窗口

重新运行 cleaner 需要调用 RXNMapper API（网络调用），20 条反应预计耗时 5–10 分钟。在此期间 RPH 使用的仍是旧 CSV。建议：
1. 在非工作时段重新运行 cleaner
2. 新 CSV 通过 diff 验证只有 rx1 的 `core_bond_changes` 发生了变化
3. 更新 RPH 的 `data/reaxys_cleaned.csv` 后再触发全链路测试
