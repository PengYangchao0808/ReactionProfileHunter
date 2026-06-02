# Cleaner 内部根因报告

## 1. 任务范围

本报告只分析外部 cleaner：

- 路径：`exclean/cleaner_20260406`

目标：

> 解释为什么在用户已固定的金标准下，cleaner 内部会把正确 forming bond `12-19` 写成 `12-13`。

本报告继续沿用唯一编号规范：

- **金标准编号空间**：product / TS 输入文件的 **XYZ 1-based** 编号
- **金标准 forming bonds**：
  - `12-19`
  - `15-16`

---

## 2. 先给结论

### Cleaner 内部最可能的真实根因

cleaner 内部不是简单的“写错了一列 CSV”，而是有两层更早的问题：

1. **RXNMapper 低置信度映射被接受并继续进入主流程**
2. **CoreExtractor 完全依赖 `rxn_smiles_mapped` 的图连通关系来推断 formed/broken bonds**，并不参考 product/TS 几何金标准

因此：

> 如果 `rxn_smiles_mapped` 本身已经把反应图写成“12-13 在产物中成键、12-19 不属于 formed bond”，那么 CoreExtractor 就会稳定地产出 `12-13:formed`。

也就是说：

> **在 cleaner 内部，`12-19 → 12-13` 的首次失真点，高概率已经发生在 atom mapping 结果本身；CoreExtractor 只是把这个错误忠实地形式化为 `core_bond_changes`。**

---

## 3. cleaner 内部数据传播链

## 3.1 Step A：RXNMapper 产出 `rxn_smiles_mapped`

文件：

- `parsers/atom_mapper.py:323-368`

关键逻辑：

```python
outputs = self._mapper.get_attention_guided_atom_maps([normalized_reaction_str])
item = outputs[0]
mapped_rxn = item["mapped_rxn"]
confidence = float(item["confidence"])
```

对于低置信度 mapping：

```python
if confidence >= self.confidence_threshold:
    status = "OK"
else:
    sanity_result = self._check_mapping_sanity(mapped_rxn)
    if sanity_result['usable']:
        status = "LOW_CONFIDENCE"
        final_mapped = mapped_rxn
```

### 结论

当前 cleaner 的策略不是：

- 低置信度就拒绝

而是：

- **低置信度但“sanity 通过”就继续使用**

而本反应的 cleaner 输出明确显示：

- `map_status = LOW_CONFIDENCE`
- `mapping_confidence = 0.225214...`

这一步是 cleaner 内部的**第一处结构性风险源**。

---

## 3.2 Step B：main_v2 允许 LOW_CONFIDENCE mapping 进入 CoreExtractor

文件：

- `main_v2.py:340-386`

关键逻辑：

```python
if map_status in usable_statuses and rxn_smiles_mapped:
    sanity_result = compute_map_sanity(rxn_smiles_mapped)
    ...

extract_result = self.core_extractor.extract(
    rxn_smiles_mapped=rxn_smiles_mapped,
    reactants_all=row.get('reactants_all'),
    product_main_smiles=row.get('product_main_smiles'),
    map_status=map_status,
    precursor_type=precursor_type,
)
```

而 `CoreExtractor.USABLE_STATUSES` 定义为：

- `{"OK", "LOW_CONFIDENCE"}`

文件：

- `parsers/core_extractor.py:30-32`

### 结论

> cleaner 主流程明确允许 **LOW_CONFIDENCE** 的 atom mapping 直接进入核心成键推断。

因此，一旦低置信度 mapping 本身已经把反应图写歪，后续 `core_bond_changes` 就会在坏图上继续推理。

---

## 3.3 Step C：CoreExtractor 以 `rxn_smiles_mapped` 的 product-side 图为唯一真相

文件：

- `parsers/core_extractor.py:562-583`

关键逻辑：

```python
reactant_side, product_side = rxn_smiles_mapped.split(">>", 1)
...
product_mol = self._pick_main_product_mol(product_side)
...
product_map_idx = self._build_map_idx(product_mol)
reactant_map_idx_list = [self._build_map_idx(m) for m in reactant_mols]
```

### 关键事实

CoreExtractor 并没有用：

- product XYZ
- TS XYZ
- 任何几何文件

它只使用：

- `rxn_smiles_mapped` 的 reactant/product graph

所以：

> cleaner 内部对 bond change 的判断，是**纯图论**、纯 SMILES 层面的，不受几何金标准约束。

---

## 3.4 Step D：候选原子集来自 SMARTS 模式匹配，不来自最终几何

文件：

- `parsers/core_extractor.py:585-634`

关键逻辑：

```python
strategy_patterns = self._get_active_patterns(precursor_type)
...
matches = reactant_mol.GetSubstructMatches(patt)
...
atom_maps_for_match.add(map_num)
...
candidate_set_options.add(frozenset(merged))
```

### 结论

候选核心原子是从：

- reactant side 的 SMARTS 命中

拼出来的。

这意味着 cleaner 的核心问题不是“看了 product/TS 结果后选错”，而是：

> 它在**映射后的反应图**上先定义一组候选核心原子，再在这组原子上推 formed/broken bonds。

如果 atom mapping 本身有问题，或者候选集把原子 13 纳入并把 19 的角色定义错，那么后面推导就会系统性偏向 `12-13`。

---

## 3.5 Step E：真正生成 `core_bond_changes` 的位置

文件：

- `parsers/core_extractor.py:411-445`

这是最关键的函数：

```python
for map_i, map_j in combinations(candidate_sorted, 2):
    product_has_bond = False
    if map_i in product_map_idx and map_j in product_map_idx:
        idx_i = product_map_idx[map_i]
        idx_j = product_map_idx[map_j]
        product_has_bond = product_mol.GetBondBetweenAtoms(idx_i, idx_j) is not None

    reactant_has_bond = False
    for reactant_mol, reactant_map_idx in zip(reactant_mols, reactant_map_idx_list):
        ...
        if reactant_mol.GetBondBetweenAtoms(idx_i, idx_j) is not None:
            reactant_has_bond = True

    if product_has_bond and not reactant_has_bond:
        formed.append((map_i, map_j))
```

### 这一段意味着什么

Cleaner 对 formed bond 的定义非常机械：

- 在 `product_mol` 里有键
- 在所有 `reactant_mol` 里都没有键
- 就标成 `formed`

所以，如果当前 `rxn_smiles_mapped` 对应的 product-side graph 中：

- `12-13` 被表示为新出现的键
- `12-19` 没有满足“产物有 / 反应物无”的判定

那么 cleaner 必然输出：

- `12-13:formed`
- 而不是 `12-19:formed`

### 结论

> `core_bond_changes` 本身不是一个独立推理结果，它只是对 `mapped reaction graph` 的直接图差分。

---

## 3.6 Step F：CSV writer 只是原样写出，不制造新错误

文件：

- `outputs/csv_writer.py:63-95`

关键逻辑：

```python
df_output.to_csv(filepath, index=False, encoding='utf-8', quoting=1)
```

### 结论

CSV writer 只负责写盘，不参与 bond change 计算。

所以：

> **writer 不是错误源。**

---

## 4. cleaner 内部最早的真正 defect 在哪里？

## 4.1 直接触发点

**直接生成错误 `core_bond_changes` 的函数是：**

- `parsers/core_extractor.py:_infer_bond_changes()`

因为 `12-13:formed` 就是在这里被写入 `formed` 列表的。

但它并不是最根的根因。

---

## 4.2 真正更早的内部根因

Cleaner 内部最早的高概率根因是：

### 根因 A：错误或不可靠的 atom mapping 被接受

证据：

- `atom_mapper.py:339-345`
- `map_status = LOW_CONFIDENCE`
- `mapping_confidence = 0.225214`

即：

> cleaner 明知 mapping 极低置信度，却只要 `sanity_result['usable']` 就继续放行。

这使得一个可能已经错误编码了反应键连通关系的 `mapped_rxn` 进入下游。

### 根因 B：CoreExtractor 只看 graph difference，不看化学/几何金标准

证据：

- `core_extractor.py:562-583`
- `core_extractor.py:411-445`

即：

> CoreExtractor 不会问“真实 forming bond 是不是 12-19”，它只会问“mapped product graph 里哪条键存在而 reactant graph 里不存在”。

所以如果 mapping 图错了，bond change 必错。

---

## 5. 为什么 `12-19` 会被 cleaner 看成“不是 formed bond”？

当前 cleaner 代码最值得警惕的语义缺陷是：

### 5.1 它只做“有键 / 无键”差分，不处理更深层键语义

Cleaner 判断 formed/broken 的核心条件是：

- `GetBondBetweenAtoms(...) is not None`

它并不深入区分：

- aromatic → single
- 环闭合中的旧边角色变化
- 同一对 map atom 在不同反应图语义中的“新成键意义”

如果你的金标准认为：

> `12-19` 在化学机理上应被视为新的 cycloaddition forming bond

那么 cleaner 当前实现并不会做这种“机理语义判定”。

它只会做：

> “mapped reactant 图里有没有这条边？”
> “mapped product 图里有没有这条边？”

若图差分不满足，它就不会把 `12-19` 记为 formed。

### 5.2 这意味着什么

Cleaner 当前 `core_bond_changes` 的语义其实不是：

- “真实机理成键”

而是：

- “当前映射图上的拓扑差分”

这两者在你这条反应上已经发生了分离。

---

## 6. 内部根因排序

### 排名 1：低置信度 RXNMapper 结果被接受（最可能根因）

位置：

- `parsers/atom_mapper.py:339-345`

原因：

- `rxn_smiles_mapped` 是后续一切的唯一反应图基础
- 一旦这个图把 `12-19` 的角色写错，后面所有步骤都会一致地错

### 排名 2：CoreExtractor 的 bond-change 语义过于“图差分化”

位置：

- `parsers/core_extractor.py:411-445`

原因：

- 它只识别 `product has bond && reactant no bond`
- 不识别“机理意义上的 forming bond”
- 对复杂环化、芳香体系、桥连体系尤其脆弱

### 排名 3：候选核心原子依赖 SMARTS 命中组合，容易把错误 mapping 放大

位置：

- `parsers/core_extractor.py:585-634`

原因：

- 候选集不是从 product/TS 几何回推
- 而是从 mapped reactant SMARTS 命中拼出来
- 所以上游 mapping 偏差会直接污染候选集和后续 formed/broken 判断

### 排名 4：writer 导出污染（基本可排除）

位置：

- `outputs/csv_writer.py`

原因：

- 这里只写盘
- 没有改写 bond pair

---

## 7. 最终裁决

### Brutal Verdict

Cleaner 的问题不是 CSV 写错，也不是简单编号偏移。

Cleaner 当前真正的问题是：

> **它把一个低置信度 atom-mapped reaction 当成可信反应图，再用纯图差分逻辑去定义“forming bonds”。**

所以当你的金标准与这个 graph-diff 语义不一致时，cleaner 会稳定地产出看起来格式正确、实际上化学意义错误的 `core_bond_changes`。

对这条反应而言，最核心的一句话是：

> **Cleaner 内部第一次真正出问题的地方，不在 CSV，不在 writer，不在 repo；而是在“LOW_CONFIDENCE mapping 被接受”之后，`_infer_bond_changes()` 对这个坏图做了无条件图差分。**

---

## 8. 后续建议

下一步如果继续追 cleaner，应继续拆成两个子问题：

### 子问题 A：mapping 问题

- `rxn_smiles_mapped` 本身为什么会把这条反应编码成当前这种图连通关系？

### 子问题 B：语义问题

- cleaner 是否应该继续把 `core_bond_changes` 定义成纯 graph-diff
- 还是要引入“机理意义上的 forming bond”层

如果不拆开，后续修复会继续把 mapping 问题和 bond semantics 问题混成一团。
