# 空间编号问题报告

## 1. 任务目标

本报告仅聚焦当前链路中的**空间编号混乱问题**，不讨论化学机理正确性本身。

本次排查采用用户指定的**金标准**：

- 以 **产物 / TS 输入文件**（`product_min.xyz` / `ts_final.xyz`）的原子编号为最高标准
- 在该编号体系下，真实 forming bonds 为：
  - **15-16**
  - **12-19**

本报告目标是逐步确认当前系统中各类编号空间的映射关系，并定位编号混乱发生的位置。

---

## 2. 金标准：产物 / TS 输入文件编号

### 2.1 权威输入文件

- `Output/benchmark_dft_theory/baselines/rx1/bl_3f31ea3a42fe5184/S1_ConfGeneration/product_min.xyz`
- `Output/benchmark_dft_theory/baselines/rx1/bl_3f31ea3a42fe5184/stationary_points/ts_final.xyz`

### 2.2 XYZ 文件编号规则

XYZ 文件对用户可见的编号是 **1-based**：

- 第 3 行原子 = atom 1
- 第 4 行原子 = atom 2
- ...

### 2.3 金标准成键距离验证

直接对 `product_min.xyz` 和 `ts_final.xyz` 计算得到：

| Bond (XYZ 1-based) | product | TS |
|---|---:|---:|
| 12-19 | 1.571 Å | 1.571 Å |
| 15-16 | 1.545 Å | 1.545 Å |

因此本报告后续所有“正确编号”均以：

- **12-19**
- **15-16**

作为金标准表达。

---

## 3. 当前系统中实际存在的编号空间

当前链路中至少存在 5 套编号空间：

| 编号空间 | 含义 | 基数 |
|---|---|---|
| XYZ 显示编号 | `.xyz` 文件中的人类可见编号 | 1-based |
| XYZ 内部编号 | `read_xyz()` 读入后的 `coords[i]` 下标 | 0-based |
| mapped SMILES atom-map | `[C:12]` 里的 `:12` | 通常 1-based |
| RDKit MolIdx | `atom.GetIdx()` | 0-based |
| S0/S2/S3 JSON 中 forming_bonds | 代码声称为内部索引，但实际来源混杂 | 混乱 |

---

## 4. 金标准编号与内部编号的精确对应

### 4.1 XYZ 1-based → XYZ 0-based

内部数组编号永远比文件显示编号小 1：

- XYZ **12-19** → 内部 **11-18**
- XYZ **15-16** → 内部 **14-15**

即：

| 金标准（XYZ 1-based） | 内部 XYZ / MolIdx（0-based） |
|---|---|
| 12-19 | 11-18 |
| 15-16 | 14-15 |

---

## 5. mapped product 与 product XYZ 的映射关系

### 5.1 使用的数据源

来自 `data/reaxys_cleaned.csv` 当前目标反应行的 `rxn_smiles_mapped`：

```text
[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][CH2:11][C@@:12]23[CH:13]=[CH:14][C@@H:15]([CH2:16][C:17](=[O:18])[C@H:19]12)[O:20]3
```

### 5.2 map number → RDKit MolIdx

程序实际解析结果：

- map 1 → molIdx 0
- map 2 → molIdx 1
- ...
- map 20 → molIdx 19

即当前这条反应里：

> **map number = MolIdx + 1**

### 5.3 map number → product XYZ

通过 `rph_core.utils.cleaner_adapter.get_map_to_xyz_dict()` 对 `mapped_product_smiles` 与 `product_min.xyz` 做子结构匹配，得到：

- map 1 → xyz0 0 → xyz1 1
- map 2 → xyz0 1 → xyz1 2
- ...
- map 20 → xyz0 19 → xyz1 20

即当前这条反应里：

> **map number = product XYZ 1-based heavy-atom 编号**

这是本次编号调查最重要的中间结论。

### 5.4 对金标准成键的直接映射

据此可得：

| 表达空间 | 金标准 bond 1 | 金标准 bond 2 |
|---|---|---|
| XYZ 1-based | 12-19 | 15-16 |
| XYZ 0-based / MolIdx | 11-18 | 14-15 |
| mapped atom-map | 12-19 | 15-16 |

因此，**在当前这条反应里，map 编号与 product/TS 文件的 1-based heavy-atom 编号是同一套编号。**

---

## 6. 当前错误结果如何在系统中出现

### 6.1 当前 cleaner 产物中的键对

当前 cleaner / repo 行里写的是：

- `core_bond_changes = 12-13:formed;15-16:formed;...`
- `formed_bond_map_pairs = 12-13;15-16`
- `formed_bond_index_pairs = 11-12;14-15`
- `forming_bonds = 11-12;14-15`

这说明当前链路已经把“错误的 map 键对”继续投影成了内部 0-based 键对：

| 当前系统使用值 | 实际对应的 XYZ 1-based |
|---|---|
| 11-12 | 12-13 |
| 14-15 | 15-16 |

也就是说：

- 当前系统保留了正确的 **15-16**
- 但把金标准 **12-19** 错写成了 **12-13**

---

## 7. 各模块中的编号处理链路

## 7.1 `rph_core/utils/cleaner_adapter.py`

### `map_pairs_to_internal_indices()`

位置：`rph_core/utils/cleaner_adapter.py:162-183`

功能：

- 输入：map pairs（如 `12-19`）
- 输出：MolIdx / 0-based pairs（如 `11-18`）

关键逻辑：

```python
map_to_idx[map_num] = atom.GetIdx()
```

### `map_pairs_to_xyz_indices()`

位置：`rph_core/utils/cleaner_adapter.py:186-278`

功能：

- 输入：map pairs + mapped SMILES + product XYZ
- 输出：product XYZ 0-based pairs

关键事实：

- `12-19` 会被映射为 `11-18`
- `15-16` 会被映射为 `14-15`
- 当前 cleaner 的 `12-13` 会被映射为 `11-12`

这一步本身没有错，**它只是忠实地把上游给出的 map pair 投影到 XYZ 空间。**

---

## 7.2 `rph_core/utils/dataset_loader.py`

位置：`rph_core/utils/dataset_loader.py:190-232`

功能：

- 从 `core_bond_changes` 解析形成键（map 编号）
- 再通过 `map_pairs_to_internal_indices()` 转成 0-based
- 最终写入：
  - `formed_bond_index_pairs`
  - `forming_bonds`
  - `forming_bonds_index_base = 0`

关键风险：

> 一旦上游 map 键对错了，这里会把错误稳定地转写成 0-based 正式输入。

---

## 7.3 `rph_core/steps/mechanism_classifier/clean_adapter.py`

位置：`rph_core/steps/mechanism_classifier/clean_adapter.py:166-295`

功能：

- 解析 `core_atom_map`
- 解析 `core_bond_changes`
- 将 cleaner 里的 bond changes 直接变成 S0 的 `forming` / `breaking`

关键问题：

```python
a_molidx = atom_map.get(a, a)
b_molidx = atom_map.get(b, b)
```

这里假定：

- `core_bond_changes` 中的数值是 MapId
- `core_atom_map` 能正确完成 MapId → MolIdx

但当前这条反应的 `core_atom_map` 本质上是恒等映射：

```json
{"11": 11, "12": 12, ..., "19": 19, "8": 8}
```

因此这一步没有真正消除编号空间，只是把同样的数字原样带入了 S0。

---

## 7.4 `rph_core/orchestrator.py`

### `_normalize_forming_bonds()`

位置：`rph_core/orchestrator.py:350-469`

规则：

- 若 `index_base == 1`，则减 1
- 若 `index_base == 0`，则原样使用

```python
if index_base_flag == 1:
    normalized = [(i - 1, j - 1) for i, j in raw_pairs]
else:
    normalized = list(raw_pairs)
```

### S0 summary 被强行当成 0-based

位置：`rph_core/orchestrator.py:488-506`

```python
s0_bonds = self._normalize_forming_bonds(
    s0_bonds_raw,
    atom_count=atom_count,
    index_base=0,
    require_exact_two=True,
)
```

这意味着：

> **只要 S0 summary 里写了 `[12,13]`，S2 就会无条件把它当作内部 0-based `(12,13)` 使用。**

而内部 0-based `(12,13)` 实际对应的是：

- XYZ 1-based **13-14**

这与用户看到的文件编号 **12-13** 已经发生了一次静默漂移。

### S0 summary 又被写回 effective cleaner data

位置：`rph_core/orchestrator.py:1042-1047`

```python
effective_cleaner_data["forming_bonds"] = _fb
effective_cleaner_data["forming_bonds_index_base"] = 0
```

即：

> S0 summary 的值会再次被贴上“0-based”标签，继续下游传播。

---

## 7.5 `rph_core/steps/step2_retro/retro_scanner.py`

### `_map_bonds()`

位置：`rph_core/steps/step2_retro/retro_scanner.py:230-240`

逻辑：

- 如果给了 `atom_map`，尝试把 bond 映射
- 如果没有，则默认输入就是 MolIdx / 0-based

```python
if not atom_map:
    return self._validate_forming_bonds(forming_bonds)
```

因此只要上游把 map-style 数字贴成 0-based，S2 就会直接拿去扫描。

### `scan_profile.json` 写入的是 `bonds`

位置：`retro_scanner.py:341-348`

```python
"forming_bonds": [list(pair) for pair in bonds]
```

而根据 Step2 约定，这里的 `bonds` 是内部索引（0-based）。

所以当前 `scan_profile.json` 里的：

```json
"forming_bonds": [[12,13],[15,16]]
```

若按代码约定解释，应当表示：

- XYZ 1-based **13-14**
- XYZ 1-based **16-17**

但文件本身没有再次标注 index_base，极易被误读为文件编号 12-13 / 15-16。

---

## 7.6 `S3_TS/mechanism_meta.json`

当前位置内容：

```json
{
  "index_base": 0,
  "forming_bonds": [[12,13],[15,16]]
}
```

而 `orchestrator.py:1528-1539` 确实明确写成：

```python
"index_base": 0,
"forming_bonds": [list(b) for b in result.forming_bonds]
```

所以这里的解释也应为：

- 0-based `(12,13)` → XYZ 1-based **13-14**
- 0-based `(15,16)` → XYZ 1-based **16-17**

这与用户当前金标准 **12-19 / 15-16** 完全不一致。

---

## 8. 当前编号混乱的核心结论

### 结论 1：这条反应里，map 编号与 product/TS 的 1-based heavy-atom 编号完全一致

这是当前最稳定的对照基准：

> **map 12 = xyz1 12**  
> **map 19 = xyz1 19**  
> **map 15 = xyz1 15**  
> **map 16 = xyz1 16**

因此，若金标准是 **12-19 / 15-16**，那么对应的 map pair 也必须是：

- **12-19**
- **15-16**

### 结论 2：当前 cleaner / S0 / S2 / S3 链路没有保持这个关系

当前链路实际使用的是：

- map pair: `12-13;15-16`
- 0-based pair: `11-12;14-15`

这只保留了一个正确键（15-16），把另一个正确键（12-19）替换成了错误键（12-13）。

### 结论 3：系统存在“同一组数字，被不同模块当成不同编号空间”的问题

最典型的例子：

- S0 summary 写 `12-13`
- orchestrator 把它认作 0-based
- S2 / S3 继续把它写入 JSON
- 但用户看到这些数字时，会自然按 product/TS 文件编号理解成 1-based

这就造成了：

> **同一个 `12-13`，在文件阅读者眼中是 XYZ atom 12-13，**  
> **在代码里却可能被当成内部 0-based 的 12-13（即 XYZ atom 13-14）。**

这正是当前空间编号混乱的核心。

---

## 9. 需要明确区分的三条规则

后续分析中必须固定使用以下规则，否则结论会反复漂移：

### 规则 A：用户讨论化学结构时

统一使用：

- **product / TS 输入文件的 XYZ 1-based 编号**

### 规则 B：代码内部形成键计算时

统一使用：

- **0-based internal XYZ / MolIdx**

即：

- 12-19（文件） ↔ 11-18（内部）
- 15-16（文件） ↔ 14-15（内部）

### 规则 C：只要看到 mapped atom-map 数字

对本条反应可直接等同于：

- **product XYZ 1-based heavy-atom 编号**

即：

- map 12 = xyz1 12
- map 19 = xyz1 19

---

## 10. 本阶段结论

本次“空间编号问题”排查已经确认：

1. **金标准编号体系**应固定为 product / TS 输入文件的 XYZ 1-based 编号。
2. 在当前目标反应中，**mapped product atom-map 编号与 XYZ 1-based heavy-atom 编号完全一致**。
3. 因而金标准 forming bonds 的跨空间表达应为：
   - map / XYZ1: **12-19, 15-16**
   - internal 0-based: **11-18, 14-15**
4. 当前 cleaner / S0 / S2 / S3 链路实际在传播：
   - map / XYZ1: **12-13, 15-16**
   - internal 0-based: **11-12, 14-15**
5. 当前系统存在明确的**编号空间混用**：
   - 同样的数值对在不同模块中被当作不同编号空间解释
   - 且 JSON 产物并未始终显式标明 index space

因此，**在继续排查化学成键逻辑之前，必须先把以下对应关系固定下来作为后续唯一编号规范：**

### 后续唯一编号规范

- 对外报告、人工分析、截图讨论：**XYZ 1-based**
- 对内代码、扫描、S4 特征：**0-based internal**
- 本条反应的桥接关系：
  - `12-19 (XYZ1 / map)` ↔ `11-18 (0-based)`
  - `15-16 (XYZ1 / map)` ↔ `14-15 (0-based)`

---

## 11. 下一步建议

在“空间编号”这一层已经厘清后，下一步应继续排查：

1. **当前 cleaner 为什么把 map pair `12-19` 错写成 `12-13`**
2. **S0 summary 为什么输出未声明空间的 `12-13 / 15-16`，却又在 S2 被按 0-based 消费**
3. **scan_profile.json / mechanism_meta.json 是否应该强制同时写出：**
   - `forming_bonds_xyz1`
   - `forming_bonds_xyz0`
   - `forming_bonds_map`

以彻底消除当前歧义。
