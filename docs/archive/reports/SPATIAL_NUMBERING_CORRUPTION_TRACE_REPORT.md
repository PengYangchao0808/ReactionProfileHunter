# 空间编号错误传播追踪报告

## 1. 任务定义

本报告建立在上一份 `SPATIAL_NUMBERING_ISSUE_REPORT.md` 的编号规范之上。

本次固定前提如下：

- **唯一金标准编号空间**：`product_min.xyz` / `ts_final.xyz` 的 **XYZ 1-based** 编号
- **唯一金标准 forming bonds**：
  - **12-19**
  - **15-16**

本报告只回答一个问题：

> 在 **cleaner → S0 → S2** 这条链路里，**哪一步第一次把金标准 `12-19` 写坏**？

---

## 2. 先给结论

### 最终结论

**第一次把金标准 `12-19` 写坏的，不是 S0，也不是 S2，而是 cleaner 输出本身。**

更准确地说：

1. **cleaner 输入文件 `data/reaxys_cleaned.csv` 的 `core_bond_changes` 已经写成了 `12-13:formed;15-16:formed`**
2. repo 内部的 loader / adapter 只是把这个错误继续正规化、投影并传播
3. **S0 不是第一次写坏，而是第一次把这个错误提升为“权威机制输出”**
4. **S2 不是第一次写坏，而是第一次把这个错误真正用于扫描执行，并且按 0-based 内部索引消费，造成更深一层语义漂移**

---

## 3. 金标准与系统当前值对比

### 3.1 金标准

| 编号空间 | 正确 forming bonds |
|---|---|
| XYZ 1-based / map | `12-19`, `15-16` |
| 内部 0-based | `11-18`, `14-15` |

### 3.2 当前系统实际传播值

| 编号空间 | 当前传播值 |
|---|---|
| map / XYZ 1-based 表观数字 | `12-13`, `15-16` |
| 内部 0-based | `11-12`, `14-15` |

即：

- 正确键 `15-16` 被保留下来了
- 正确键 `12-19` 被替换成了错误键 `12-13`

---

## 4. cleaner → S0 → S2 逐级审查

## 4.1 Stage A：cleaner 输入文件（第一次写坏发生在这里）

### 证据文件

- `data/reaxys_cleaned.csv:2`

### 关键字段

该行已经明确写出：

```text
rxn_smiles_mapped = [CH3:1]...[C@@:12]23...[C@H:19]12...[O:20]3
core_bond_changes = 12-13:formed;15-16:formed;8-13:broken;12-18:broken;13-15:broken;16-19:broken
map_status = LOW_CONFIDENCE
mapping_confidence = 0.225214...
```

### 关键观察

1. 当前 reaction row 中，`map 19` 并没有丢失——它仍然清楚地存在于 `rxn_smiles_mapped` 中
2. 但 `core_bond_changes` 已经把 forming bonds 写成：
   - `12-13`
   - `15-16`
3. 金标准要求的是：
   - `12-19`
   - `15-16`

### 本阶段判定

> **第一次把 `12-19` 写坏的阶段 = cleaner 输出 `core_bond_changes` 的阶段。**

这是最早、最原始、最上游的错误落点。

### 为什么这一步可判定为“首次写坏”

因为在 cleaner 输入文件中：

- `12` 和 `19` 仍然作为 atom-map 号存在
- 但“forming bond”字段已经不再包含 `12-19`
- 它已经被替换成 `12-13`

也就是说：

> **错误并不是在 S0 才生成，而是在 cleaner 已经生成完成。**

---

## 4.2 Stage B：repo loader / cleaner adapter（第一次 repo 内正规化错误）

这一层虽然不属于你要求的三大阶段名称之一，但它是 cleaner 和 S0 之间真正的数据桥，因此必须单独列出。

### 关键代码 1：`rph_core/utils/dataset_loader.py:190-232`

核心逻辑：

```python
formed_map_pairs = parse_formed_pairs_from_core_bond_changes(core_changes)
formed_index_pairs = map_pairs_to_internal_indices(index_source_smiles, formed_map_pairs)
raw["formed_bond_map_pairs"] = _pairs_to_str(formed_map_pairs)
raw["formed_bond_index_pairs"] = _pairs_to_str(formed_index_pairs)
raw["forming_bonds"] = raw["formed_bond_index_pairs"]
raw["forming_bonds_index_base"] = "0"
```

### 关键代码 2：`rph_core/utils/cleaner_adapter.py:95-122`

同样逻辑也在这里出现：

```python
formed_map_pairs = parse_formed_pairs_from_core_bond_changes(...)
formed_index_pairs = map_pairs_to_internal_indices(index_source_smiles, formed_map_pairs)
raw["formed_bond_map_pairs"] = ...
raw["formed_bond_index_pairs"] = ...
raw["forming_bonds"] = raw["formed_bond_index_pairs"]
raw["forming_bonds_index_base"] = "0"
```

### 实际效果

loader / adapter 并没有重新判断化学正确性，它只是把 cleaner 的错误正规化：

- cleaner 给 `12-13;15-16`
- loader 继续生成：
  - `formed_bond_map_pairs = 12-13;15-16`
  - `formed_bond_index_pairs = 11-12;14-15`
  - `forming_bonds = 11-12;14-15`
  - `index_base = 0`

### 本阶段判定

> 这一层**不是第一次写坏**，但它是**第一次在 repo 内把错误键对正规化成内部字段**。

也就是说：

- cleaner 首次写坏
- loader 首次把坏值转成 repo 内的正式输入

---

## 4.3 Stage C：S0（不是首次写坏，而是首次“官方确认错误”）

### 证据文件

- `Output/.../S0_Mechanism/mechanism_graph.json`
- `Output/.../S0_Mechanism/mechanism_summary.json`

### 证据 1：`mechanism_graph.json`

可见：

```json
"core_bond_changes": "12-13:formed;15-16:formed;..."
"formed_bond_map_pairs": "12-13;15-16"
"formed_bond_index_pairs": "11-12;14-15"
"forming_bonds": "11-12;14-15"
```

同时 edge 中又被写成：

```json
"forming_bonds": [[12,13]]
"forming_bonds": [[15,16]]
```

### 对应代码路径

#### `rph_core/steps/mechanism_classifier/clean_adapter.py:169-204`

S0 CleanAdapter 直接读取 cleaner 的 `core_bond_changes`：

```python
bond_changes = self._parse_bond_changes(row.get('core_bond_changes', ''), atom_map=core_atom_map)
```

#### `rph_core/steps/mechanism_classifier/graph_builder.py:260-303`

GraphBuilder 直接使用 `record.core_bond_changes`：

```python
forming = bond_changes.get('forming', [])
...
step1_forming = [forming[0]]
step2_forming = forming[1:]
```

### `mechanism_summary.json`

当前 S0 summary 最终写成：

```json
"forming_bonds": [[12,13],[15,16]]
```

其生成代码在：

- `rph_core/orchestrator.py:890-909`

```python
for edge in pathway_edges:
    all_pairs.extend(edge.forming_bonds)
summary["forming_bonds"] = [list(pair) for pair in forming_bonds]
```

### 本阶段判定

> **S0 不是第一次把 `12-19` 写坏。**

S0 做的事情是：

1. 直接接收 cleaner 已经写坏的键对
2. 将其写入 `mechanism_graph.json`
3. 再把它提升为 `mechanism_summary.json` 中的“权威 forming_bonds”

因此：

> **S0 是第一次“官方确认并发布这个错误”的阶段。**

---

## 4.4 Stage D：S2（不是首次写坏，而是首次把错误拿去真正执行扫描）

### 关键代码路径

#### `rph_core/orchestrator.py:488-506`

S2 优先读取 S0 summary：

```python
s0_bonds_raw = s0_data.get("forming_bonds")
s0_bonds = self._normalize_forming_bonds(
    s0_bonds_raw,
    atom_count=atom_count,
    index_base=0,
    require_exact_two=True,
)
return s0_bonds
```

这一步非常关键：

> **S2 明确把 S0 summary 中的数字解释为 0-based internal indices。**

也就是说：

- S0 summary 写 `12-13`
- S2 会按 0-based 理解成内部 `(12,13)`
- 这在 XYZ 1-based 空间中实际对应的是 **13-14**

同理：

- `15-16` 会被按 0-based 理解成内部 `(15,16)`
- 在 XYZ 1-based 空间中对应 **16-17**

### 证据文件：`S2_Retro/scan_profile.json`

当前内容：

```json
"forming_bonds": [[12,13],[15,16]]
```

其写入位置在：

- `rph_core/steps/step2_retro/retro_scanner.py:341-348`

```python
"forming_bonds": [list(pair) for pair in bonds]
```

而 `bonds` 在 Step2 约定中是内部索引。

### 本阶段判定

> **S2 不是第一次写坏 `12-19`。**

但 S2 是：

1. 第一次把 S0 的错误键对真正用于扫描执行
2. 第一次把这些数字按 0-based internal indices 严格消费
3. 第一次把“化学错误 + 编号空间错误”叠加在一起

因此：

> **S2 是第一次让这个错误真正变成计算行为错误的阶段。**

---

## 5. 分阶段责任判定

| 阶段 | 是否第一次写坏 `12-19` | 本阶段作用 |
|---|---|---|
| cleaner | **是** | 首次把正确键 `12-19` 写成 `12-13` |
| repo loader / cleaner adapter | 否 | 把 cleaner 的错误正规化为 repo 内部字段 |
| S0 | 否 | 把错误提升为“权威机制输出” |
| S2 | 否 | 把错误真正用于扫描执行，并按 0-based 消费 |

---

## 6. 必须区分的两类“写坏”

为了避免后续讨论混乱，这里必须把两类错误分开：

### 6.1 第一类：化学键身份写坏

即：

- 正确应为 `12-19`
- 却被写成 `12-13`

这一类错误的首次出现点是：

> **cleaner 的 `core_bond_changes`**

### 6.2 第二类：编号空间语义写坏

即：

- 原本数字长得像 map / XYZ1
- 却被下游当成 0-based internal 使用

这一类错误最关键的发生点是：

> **S2 在 `_resolve_forming_bonds_for_s2()` 中以 `index_base=0` 消费 S0 summary**

因此：

- cleaner 首先犯“键身份错误”
- S2 再叠加“编号语义错误”

---

## 7. 最终裁决

### Brutal Verdict

如果金标准固定为：

- `12-19`
- `15-16`

那么 cleaner → S0 → S2 这条链路的真实情况是：

1. **cleaner 先把 bond identity 写错了**
2. **S0 没有纠错，反而把它提升成权威输出**
3. **S2 再把这组错误数字当成 0-based internal indices 执行，进一步加重偏差**

所以本问题不能再模糊描述为“某处编号有点乱”，而应明确写成：

> **第一次把 `12-19` 写坏的是 cleaner。**  
> **第一次把这个坏值制度化的是 S0。**  
> **第一次把这个坏值真正拿去驱动错误扫描的是 S2。**

---

## 8. 后续建议

在后续修复报告中，应明确拆成两个独立问题：

### 问题 A：bond identity corruption

- cleaner 为什么把 `12-19` 推导成 `12-13`

### 问题 B：index-space corruption

- 为什么 S0 / S2 输出没有始终显式标明这组数字到底属于：
  - map
  - XYZ 1-based
  - XYZ 0-based / MolIdx

只有把这两个问题拆开，后续修复方案才不会继续互相污染。
