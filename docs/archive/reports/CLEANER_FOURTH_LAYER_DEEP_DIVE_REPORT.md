# Cleaner 第四层深挖报告

## 1. 任务定义

本报告只回答第四层问题：

> 在 cleaner 已经接受的这条 `rxn_smiles_mapped` 内部，`12 / 13 / 19` 到底是怎样被编码的？为什么 `12-19` 无法在 cleaner 内部被识别为 formed，而 `12-13` 会被稳定识别为 formed？

本报告仍采用固定金标准：

- **唯一金标准编号空间**：product / TS 输入文件 XYZ 1-based
- **唯一金标准 forming bonds**：
  - `12-19`
  - `15-16`

但本报告的分析对象不是几何文件，而是：

- cleaner 已接受的 `rxn_smiles_mapped`
- cleaner 内部 `CoreExtractor` 对这条 map 的图语义解释

---

## 2. 先给最终结论

### 最终裁决

第四层最深结论已经明确：

> **在 cleaner 已接受的 `rxn_smiles_mapped` 中，`12-19` 被编码成“reactant 已有键，product 仍有键（只是 bond type 变化）”；而 `12-13` 被编码成“reactant 无键，product 有键”。**

因此对 cleaner 当前实现来说：

- `12-19` 不可能被识别为 formed
- `12-13` 必然被识别为 formed

换句话说：

> **第四层最深 defect 不是 candidate 评分，也不是 S2，而是 cleaner 接受的 mapped reaction 自身已经把 atom-role 语义写成了“12-19 是保留键 / 键级变化，12-13 是新键”。**

---

## 3. 直接图审计证据

### 3.1 审计对象

Cleaner 当前接受的 `rxn_smiles_mapped`：

```text
[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]([CH2:9][CH2:10][CH2:11][c:12]1[o:18][cH:17][cH:16][cH:19]1)[CH:13]=[C:15]=[CH2:14]
>>
[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][CH2:11][C@@:12]23[CH:13]=[CH:14][C@@H:15]([CH2:16][C:17](=[O:18])[C@H:19]12)[O:20]3
```

### 3.2 关键原子对的图差分结果

直接对 accepted mapped reaction 做 RDKit 图审计得到：

| Bond | Reactant | Product | cleaner 图语义 |
|---|---|---|---|
| 12-19 | AROMATIC | SINGLE | **preserved / bond-order changed** |
| 12-13 | None | SINGLE | **formed** |
| 15-16 | None | SINGLE | **formed** |
| 8-13 | SINGLE | None | broken |
| 12-18 | AROMATIC | None | broken |
| 13-15 | DOUBLE | None | broken |
| 16-19 | AROMATIC | None | broken |

这张表已经足够说明 cleaner 为什么会输出：

```text
12-13:formed;15-16:formed;8-13:broken;12-18:broken;13-15:broken;16-19:broken
```

因为对于 cleaner 而言，这就是该 mapped graph 的直接拓扑差分。

---

## 4. 为什么 `12-19` 在 cleaner 内永远不可能被判成 formed

## 4.1 `_infer_bond_changes()` 的定义

文件：

- `exclean/cleaner_20260406/parsers/core_extractor.py:411-445`

核心逻辑：

```python
if product_has_bond and not reactant_has_bond:
    formed.append((map_i, map_j))
elif reactant_has_bond and not product_has_bond:
    broken.append((map_i, map_j))
```

### 对 `12-19` 代入

根据图审计：

- reactant `12-19` = 有键（AROMATIC）
- product `12-19` = 有键（SINGLE）

所以：

- `product_has_bond = True`
- `reactant_has_bond = True`

最终既不满足：

- `product_has_bond and not reactant_has_bond`

也不满足：

- `reactant_has_bond and not product_has_bond`

因此：

> **`12-19` 在 cleaner 当前语义中只能被当作 preserved / bond-order-changed，绝不会进入 formed。**

---

## 4.2 对 `12-13` 代入

根据图审计：

- reactant `12-13` = 无键
- product `12-13` = SINGLE

所以：

- `product_has_bond = True`
- `reactant_has_bond = False`

满足：

```python
if product_has_bond and not reactant_has_bond:
    formed.append((12, 13))
```

即：

> **只要当前 accepted mapped reaction 不变，`12-13` 就一定会被 cleaner 判成 formed。**

---

## 5. 第四层对 candidate generation 的裁决

## 5.1 candidate set 并没有把 19 丢掉

直接复现 `CoreExtractor` 的候选集生成后，得到：

### active patterns

来自：

- `chemistry_templates.py:64-69`
- `core_extractor.py:43-50, 585`

对于 `allenamide`，启用：

- `furan_diene`
- `furan_substituted`
- `allenamide_3c`
- `alkene_generic`

### 实际 match sets

复现结果：

- `furan_diene` → `[12,16,17,18,19]`
- `furan_substituted` → `[11,12,16,17,18,19]`
- `allenamide_3c` → `[8,13,14,15]`
- `alkene_generic` → `[13,15]`, `[14,15]`

### candidate_set_options

最终只有 **1 个 candidate set**：

```text
[8,11,12,13,14,15,16,17,18,19]
```

### 结论

这一步非常关键：

> **19 没有在 candidate generation 阶段丢失。**

因此可以明确排除：

- “12-19 因为候选原子集没包含 19，所以没被识别”

这是错误假设。

---

## 5.2 candidate competition 也不存在

由于 `candidate_set_options count = 1`，所以：

- 不存在多个候选集互相竞争
- 不存在 `12-19` 候选集在 ranking 中输给 `12-13` 候选集

因此也可以明确排除：

- “12-19 在多候选评分里输给 12-13”

这同样不是根因。

---

## 6. 第四层对 ring-size / pruning 的裁决

文件：

- `core_extractor.py:447-528`

当前这条反应复现结果：

- `candidate_sorted = [8,11,12,13,14,15,16,17,18,19]`
- `_evaluate_candidate_option(...)` 返回：
  - `formed = [(12,13),(15,16)]`
  - `broken = [(8,13),(12,18),(13,15),(16,19)]`
  - `ring_size = 7`

### 这意味着什么

1. `formed` 数量刚好是 2
2. 因此不会进入 `len(formed) > 2` 的 pair-pruning 分支
3. 不存在“12-19 原本在 formed 里，后来被剪掉”的过程

### 结论

> **ring-size pruning 不是 12-19 丢失的原因。**

因为 `12-19` 在 `_infer_bond_changes()` 阶段就没有进入 formed 列表。

---

## 7. 第四层对 mapping sanity 的裁决

文件：

- `parsers/atom_mapper.py:323-348`
- `parsers/atom_mapper.py:182-256`

### 已确认事实

1. 该 reaction 的 `mapping_confidence = 0.225214`
2. cleaner 仍然接受该 mapping，因为 `_check_mapping_sanity(mapped_rxn)` 返回 usable
3. `_check_mapping_sanity()` 只检查：
   - reactant/product atom-map overlap 数量
   - 最低阈值是否满足

### 它没有检查的内容

它不检查：

- 关键反应中心 atom-role 是否正确
- 关键 forming bond 是否符合机理预期
- 芳香体系中的 bond-role 是否发生错误归属

### 结论

> **mapping sanity 太弱，足以放行一个“atom label overlap 足够，但 bond-role 语义已经错了”的 mapped reaction。**

这正是当前第四层根因的入口。

---

## 8. 第四层症状 / 触发点 / 根因分离

## 8.1 表面症状

Cleaner 输出：

```text
12-13:formed;15-16:formed
```

而不是金标准：

```text
12-19:formed;15-16:formed
```

## 8.2 直接触发点

`CoreExtractor._infer_bond_changes()`

文件：

- `core_extractor.py:411-445`

这里根据 accepted mapped graph 的“有键 / 无键”直接写出 formed/broken。

## 8.3 真正根因

**真正根因不是 `_infer_bond_changes()` 本身。**

真正根因是：

> **cleaner 接受了一个已经把 `12-19` 编码成“保留键”的 mapped reaction。**

从这一刻起，`12-19` 在 cleaner 当前语义里就已经不可能再被恢复成 formed。

---

## 9. 第四层根因排名

### 排名 1：accepted mapped reaction 自身的 atom-role 语义错误

这是最早、最深的 defect。

证据：

- `12-19` 在 accepted mapped graph 中是 `AROMATIC -> SINGLE`
- 而不是 `None -> SINGLE`

### 排名 2：mapping sanity 只检查 overlap，不检查关键 bond-role

这是放行坏 mapping 的门禁缺陷。

### 排名 3：CoreExtractor 的 graph-diff semantics

这是错误的直接执行机制，但不是最早 defect。

### 排名 4：candidate scoring / ring-size pruning

当前反应中基本可排除，因为：

- candidate 只有 1 个
- formed 只有 2 个
- 不存在竞争 / 剪枝淘汰 `12-19`

---

## 10. 最终裁决

### Brutal Verdict

第四层已经可以做出非常明确的最终裁决：

> **Cleaner 不是在 bond extraction 阶段“误选了 12-13”，而是在更早的 accepted mapped reaction 里，就已经把 12-19 的角色编码成了“旧键仍存在”。**

于是 cleaner 后续所有逻辑都会一致地得出：

- `12-13` 是 formed
- `15-16` 是 formed
- `12-19` 不是 formed

所以第四层最重要的一句话是：

> **问题不在“12-19 为什么输给 12-13”，而在“12-19 从一开始就没被 cleaner 的 accepted map 编码成 formed bond”。**

---

## 11. 下一步建议

如果继续第五层，应该直接追问：

1. **为什么 RXNMapper 给出的 accepted mapped reaction 会把 `12-19` 编码成 reactant 已有键？**
2. **这是不是 atom identity / aromatic ring correspondence / bridgehead role assignment 的具体错误？**

也就是说，下一层该追的是：

> **从 raw reaction SMILES 到 accepted mapped reaction，这条反应的 atom-role correspondence 到底在哪一步具体歪掉了。**
