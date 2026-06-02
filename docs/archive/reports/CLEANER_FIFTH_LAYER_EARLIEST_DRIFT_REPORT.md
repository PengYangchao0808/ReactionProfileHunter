# Cleaner 第五层最早漂移点报告

## 1. 任务定义

本报告回答第五层问题：

> 从 `raw reaction SMILES -> normalized reaction -> RXNMapper output -> accepted mapped reaction` 这条链里，atom 19 的角色最早是在哪里第一次漂移的？

本报告继续采用固定前提：

- **唯一金标准编号空间**：product / TS 输入文件 XYZ 1-based
- **唯一金标准 forming bonds**：
  - `12-19`
  - `15-16`

---

## 2. 先给结论

### 最终裁决

第五层最早漂移点已经可以明确判定：

> **atom 19 的角色第一次漂移，不是在 cleaner 的 pre-mapper normalization 阶段，而是在 RXNMapper 产出 `mapped_rxn` 的那一刻。**

更准确地说：

1. `raw reaction SMILES` 与传给 mapper 的 reaction string 在本条记录上基本一致
2. pre-mapper 阶段没有证据表明 atom 19 的身份在这里就已经被改写
3. **第一次出现“12-19 被编码为保留键、12-13 被编码为新键”的，是 RXNMapper 输出的 `mapped_rxn` 本身**
4. 后续 cleaner 只是通过一个过于宽松的 sanity gate 把这个低置信度错误映射放行

---

## 3. 直接断点证据

## 3.1 cleaner 当前记录

来自：

- `exclean/cleaner_20260406/cleaned_output/reaxys_cleaned.csv:2`

关键字段：

```text
rxn_key = C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C>>CC(C)(C)OC(=O)N1CCC[C@@]23C=C[C@@H](CC(=O)[C@H]12)O3
map_status = LOW_CONFIDENCE
mapping_confidence = 0.225214...
rxn_smiles_mapped = [CH3:1]...[c:12]1...[cH:19]1 >> [CH3:1]...[C@H:19]12...
```

---

## 3.2 重新计算 mapper 输入

直接调用 cleaner 内部代码复现：

- `AtomMapper._build_rxn_key(rxn_smiles_raw, product_main_smiles)`

得到：

```text
rxn_smiles_raw = C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C>>CC(C)(C)OC(=O)N1CCC[C@@]23C=C[C@@H](CC(=O)[C@H]12)O3
built_rxn_key = C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C>>CC(C)(C)OC(=O)N1CCC[C@@]23C=C[C@@H](CC(=O)[C@H]12)O3
stored_rxn_key = C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C>>CC(C)(C)OC(=O)N1CCC[C@@]23C=C[C@@H](CC(=O)[C@H]12)O3
```

### 这说明什么

在本条记录上：

> **传给 RXNMapper 的 reaction string 与 cleaner 当前存储的 rxn_key 是一致的。**

因此当前没有证据支持：

- pre-mapper normalization 先把 atom 19 的角色写歪

---

## 3.3 sanity gate 结果

直接复现两个 sanity gate：

```text
compute_map_sanity = {
  sanity_pass: True,
  overlap_count: 19,
  overlap_ratio: 1.0,
  precursor_map_count: 19,
  product_map_count: 20,
  reason: OK
}

_check_mapping_sanity = {
  usable: True,
  intersection_size: 19,
  n_min: 6,
  reason: OK: 19 >= 6
}
```

### 这说明什么

这两个 gate 都只验证：

- atom-map overlap 是否足够多

它们都不验证：

- atom 19 在反应中心的角色是否合理
- 关键 forming bond 是否符合机理预期

因此：

> **sanity gate 不是最早 defect，但它是让 RXNMapper 错图继续流入 cleaner 主流程的关键放行器。**

---

## 4. 第五层按阶段裁决

## 4.1 Pre-mapper normalization 阶段

### 涉及代码

- `parsers/smiles_surgeon.py:134-148` `_normalize_reaction_smiles`
- `parsers/smiles_surgeon.py:241-253` `product_main_smiles / rxn_key` 生成
- `parsers/atom_mapper.py:130-150` `_build_rxn_key`
- `parsers/atom_mapper.py:152-166` `_apply_length_guard`

### 本条记录上的实际判定

#### (1) arrow normalization

只是规范 `>` / `>>`，不涉及 atom-role。

#### (2) product selection

虽然代码支持多产品选择，但本条记录实测：

- `built_rxn_key == stored_rxn_key`

没有观察到产品选择歧义导致的 role drift。

#### (3) length guard

本条反应不在超长 SMILES 场景，未触发会改变结构语义的压缩逻辑。

### 结论

> **在本条记录上，没有证据表明 atom 19 的角色在 pre-mapper normalization 阶段已经漂移。**

也就是说：

- pre-mapper 可能是一般性风险点
- 但不是本案例的最早 concrete defect

---

## 4.2 RXNMapper 输出阶段

### 涉及代码

- `parsers/atom_mapper.py:323-345`

关键逻辑：

```python
outputs = self._mapper.get_attention_guided_atom_maps([normalized_reaction_str])
item = outputs[0]
mapped_rxn = item["mapped_rxn"]
confidence = float(item["confidence"])
```

### 本条记录上的实际判定

Cleaner 最终接受的 `rxn_smiles_mapped` 就是这个 `mapped_rxn`。

而第四层已确认：

- 在这个 accepted mapped reaction 中
  - `12-19 = react 有键, prod 有键`
  - `12-13 = react 无键, prod 有键`

因此：

> **atom 19 的角色第一次被写成“仍与 12 相连”的，就是 RXNMapper 输出本身。**

这是第五层最核心结论。

---

## 4.3 Post-mapper acceptance 阶段

### 涉及代码

- `parsers/atom_mapper.py:339-345`
- `parsers/atom_mapper.py:182-256`
- `main_v2.py:340-386`

关键逻辑：

```python
if confidence < threshold:
    sanity_result = self._check_mapping_sanity(mapped_rxn)
    if sanity_result['usable']:
        status = 'LOW_CONFIDENCE'
        final_mapped = mapped_rxn
```

### 本条记录上的实际判定

这一步并没有改写 atom 19 的角色，而是：

- **把已经漂移过的 mapped_rxn 放行了**

因此它的定位应为：

- **trigger / gate failure**
- 不是 earliest defect

---

## 5. 第五层 症状 / 触发点 / 根因分离

## 5.1 表面症状

最终 cleaner 输出：

```text
12-13:formed;15-16:formed
```

而不是金标准：

```text
12-19:formed;15-16:formed
```

## 5.2 直接触发点

LOW_CONFIDENCE 的 `mapped_rxn` 被 sanity gate 放行。

## 5.3 真正最早 defect

> **RXNMapper 输出本身第一次把 atom 19 的角色写歪。**

不是：

- csv writer
- CoreExtractor
- S0 / S2
- 本案例中的 pre-mapper normalization

---

## 6. 为什么第五层不是 pre-mapper 问题

虽然代码中存在几个理论上可能造成漂移的 pre-mapper 风险点：

- canonicalization
- product_main 选择
- reactant 排序
- length guard

但在本案例上：

1. `built_rxn_key == stored_rxn_key`
2. mapper 输入与当前 cleaner 存储的 reaction key 一致
3. 没有观察到“在 mapper 前就已经把 atom 19 改成别的角色”的证据

所以必须区分：

- **一般性风险点**
- **本案例最早 concrete defect**

本案例的最早 concrete defect 明显属于后者：

> **RXNMapper output stage**

---

## 7. 结合外部证据的解释

外部资料显示，RXNMapper / 类似 mapper 在以下体系上容易出现 atom-role drift：

- pericyclic cycloaddition
- intramolecular ring closure
- aromatic / furan systems
- bridge-forming products

这些失真模式的共同点是：

- atom overlap 仍然很高
- 但反应中心个别原子的角色归属会错

这与本案例完全一致：

- overlap = 19 / 19
- ratio = 1.0
- 但 atom 19 的反应中心角色被编码错了

因此：

> 本案例非常符合“mapper 输出局部角色漂移，但全局 overlap 看起来正常”的典型失败模式。

---

## 8. 根因排序

### 排名 1：RXNMapper 输出本身（最早 defect）

这是第五层结论。

### 排名 2：post-mapper overlap-only sanity（放行机制）

它不制造漂移，但让漂移结果继续存活。

### 排名 3：pre-mapper normalization / product selection（本案例弱证据）

理论上存在风险，但本案例没有直接证据支持它是最早 defect。

---

## 9. 最终裁决

### Brutal Verdict

第五层已经足够明确：

> **本案例的 atom 19 最早不是在 cleaner 前处理里被改坏的，而是在 RXNMapper 返回 `mapped_rxn` 时第一次被改坏。**

随后 cleaner 又犯了第二个错误：

> **它用一个只检查 overlap 的 sanity gate，把这个低置信度错误映射当成“可用”放行了。**

所以第五层的最短结论就是：

> **earliest drift = RXNMapper output**  
> **survival mechanism = overlap-only LOW_CONFIDENCE acceptance**

---

## 10. 下一步建议

如果继续第六层，最有价值的问题将变成：

> **能否直接在这条 raw reaction 上重放 / 对比多个 mapper 输出，证明 atom 19 的角色漂移是 RXNMapper 的局部错误而非 cleaner 前处理副作用？**

也就是：

1. 同一 raw reaction
2. 不同 mapper / 不同参数 / 不同 canonicalization 入口
3. 比较 atom 19 在 reaction center 的角色是否稳定

这会把问题进一步压到“模型输出稳定性”层面。
