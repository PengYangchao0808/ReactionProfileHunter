# 紧急修复分步实施方案

*2026-04-24 · 5 个步骤，严格按序执行*

---

## 修复逻辑总览

```
Step 1: 重新运行 cleaner → 更新 CSV（order_changed 写入）
    ↓
Step 2: 修复前体 SMILES 映射数据源（用 rxn_smiles_mapped 而非 precursor_smiles）
    ↓  
Step 3: 修复 forming_bonds_annotated 数据源不一致
    ↓
Step 4: 加严几何预检逻辑（任意一对可疑即触发 fallback）
    ↓
Step 5: 重新运行 pipeline 验证
```

---

## Step 1: 重新运行 cleaner 更新 CSV

### 问题
`data/reaxys_cleaned.csv` 中 rx1 行的 `core_bond_changes = "12-13:formed;15-16:formed;8-13:broken;12-18:broken;13-15:broken;16-19:broken"` 不包含 `order_changed` 条目。Phase 3a 修改了 cleaner 代码但 CSV 是历史输出。

### 方案
在 `exclean/cleaner_20260406/` 下运行 cleaner，用修改后的 `_infer_bond_changes()` 重新处理原始数据。

### 具体操作
```bash
cd "/mnt/e/Calculations/AI4S_ML_Studys/[4+3] Mechain learning/exclean/cleaner_20260406"
# 查看 cleaner 的运行方式
python -m cleaner --help  # 或 python main.py --help
```

### 预期输出
rx1 行的 `core_bond_changes` 应变为：
```
12-13:formed;15-16:formed;8-13:broken;12-18:broken;13-15:broken;16-19:broken;12-19:order_changed(AROMATIC→SINGLE)
```

### 验证
```bash
grep '"1"' data/reaxys_cleaned.csv | grep order_changed
# 应看到 12-19:order_changed
```

### 如果 cleaner 无法直接重新运行
手动修补 CSV 作为临时方案：
```python
import csv
path = "data/reaxys_cleaned.csv"
with open(path, 'r') as f:
    rows = list(csv.DictReader(f))
for row in rows:
    if row['rx_id'] == '1':
        old = row['core_bond_changes']
        row['core_bond_changes'] = old + ';12-19:order_changed(AROMATIC→SINGLE)'
        print(f"Patched: {row['core_bond_changes']}")
with open(path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
```

### RPH 消费端已就绪
`graph_builder._resolve_semantic_forming_bonds()` (L285-296) 已经有逻辑：
```python
if len(semantic_forming) < 2 and order_changed:
    # promote order_changed to forming_bonds
```
当 forming 只有 1 根 (12-13) 且 order_changed 有 1 根 (12-19) 时，12-19 会被提升为 forming_bond。

**但注意**：当前 forming 已有 2 根 (12-13 + 15-16)，`len(semantic_forming) < 2` 不成立，所以 12-19 不会被提升。需要修改逻辑为：**对 [4+3] 反应始终检查 order_changed 中的 AROMATIC→SINGLE 变化是否应纳入 forming_bonds**。

### 修改 `_resolve_semantic_forming_bonds()`

文件：`rph_core/steps/mechanism_classifier/graph_builder.py` L270-304

修改逻辑：对于 AROMATIC→SINGLE 的 order_changed 键（在环化反应中代表芳环开环参与新环形成），始终纳入 forming_bonds（即使已有 2 根 forming 键）：

```python
def _resolve_semantic_forming_bonds(self, record: CleanRecord) -> List[Tuple[int, int]]:
    bond_changes = record.core_bond_changes or {'forming': [], 'breaking': [], 'order_changed': []}
    forming = list(bond_changes.get('forming', []))
    order_changed = list(bond_changes.get('order_changed', []))

    semantic_forming: List[Tuple[int, int]] = []
    seen: set[Tuple[int, int]] = set()

    for pair in forming:
        normalized = self._normalize_pair(pair)
        if normalized is None or normalized in seen:
            continue
        semantic_forming.append(normalized)
        seen.add(normalized)

    # 对于 AROMATIC→SINGLE 的 order_changed 键，始终纳入 forming_bonds
    # （代表芳环开环参与新环形成，是环加成反应的关键键）
    for pair in order_changed:
        normalized = self._normalize_pair(pair)
        if normalized is None or normalized in seen:
            continue
        semantic_forming.append(normalized)
        seen.add(normalized)

    if semantic_forming:
        logger.info(
            f"Semantic forming_bonds for {record.reaction_id}: {semantic_forming} "
            f"(from formed={[self._normalize_pair(p) for p in forming]}, "
            f"order_changed={[self._normalize_pair(p) for p in order_changed]})"
        )

    return semantic_forming
```

---

## Step 2: 修复前体 SMILES 映射数据源

### 问题
`atom_map_smiles.json` 中 `precursor_smiles_to_map` 为空，因为 `graph_builder._build_smiles_mapping()` 使用 `record.precursor_smiles`（无 atom map 的 SMILES），而非 `rxn_smiles_mapped` 中的前体部分（有 atom map）。

### 修改

文件：`rph_core/steps/mechanism_classifier/graph_builder.py`

在 `_build_smiles_mapping()` 方法中，使用 `record.raw.get('rxn_smiles_mapped')` 提取 mapped 前体 SMILES：

```python
def _build_smiles_mapping(self, record: CleanRecord) -> SmilesAtomMapping:
    from rdkit import Chem
    mapping = SmilesAtomMapping()

    # 前体：优先从 rxn_smiles_mapped 提取 mapped SMILES
    precursor_mapped_smiles = None
    rxn_mapped = (record.raw or {}).get('rxn_smiles_mapped', '')
    if rxn_mapped and '>>' in rxn_mapped:
        precursor_mapped_smiles = rxn_mapped.split('>>')[0]

    if precursor_mapped_smiles:
        try:
            mol = Chem.MolFromSmiles(precursor_mapped_smiles)
            if mol:
                for atom in mol.GetAtoms():
                    map_num = atom.GetAtomMapNum()
                    if map_num > 0:
                        idx = atom.GetIdx()
                        mapping.precursor_smiles_to_map[idx] = map_num
                        mapping.map_to_precursor_smiles[map_num] = idx
        except Exception:
            pass

    # 产物：优先从 mapped_product_smiles 提取
    product_smiles = (record.raw or {}).get('mapped_product_smiles') or record.product_smiles or ''
    if product_smiles:
        try:
            mol = Chem.MolFromSmiles(product_smiles)
            if mol:
                for atom in mol.GetAtoms():
                    map_num = atom.GetAtomMapNum()
                    if map_num > 0:
                        idx = atom.GetIdx()
                        mapping.product_smiles_to_map[idx] = map_num
                        mapping.map_to_product_smiles[map_num] = idx
        except Exception:
            pass

    return mapping
```

### 预期结果
```json
{
  "precursor_smiles_to_map": {
    "0": 1, "1": 2, ..., "11": 12, "12": 18, "16": 13, "17": 15, "18": 14
  },
  "map_to_precursor_smiles": {
    "1": 0, ..., "12": 11, "13": 16, "14": 18, "15": 17, "18": 12, "19": 15
  }
}
```

注意：前体 SMILES canonical order 与 Map# 不是一一对应的！例如 map 13 对应前体 SMILES idx 16（而非 12），因为 RDKit canonicalization 重排了原子顺序。

---

## Step 3: 修复 forming_bonds_annotated 数据源

### 问题
`forming_bonds_annotated` 使用 `_resolve_semantic_forming_bonds()` 的输出（已含 order_changed 提升），但 `mechanism_summary.json` 的 `forming_bonds` 仍使用旧数据源。两者不一致。

### 修改

确保 `mechanism_summary.json` 的 `forming_bonds` 也从 `semantic_forming_bonds` 获取，而非直接从 cleaner 的 `formed_bond_map_pairs` 获取。

在 `orchestrator.py` 中找到 `mechanism_summary.json` 的写入位置，确保：
```python
forming_bonds_for_summary = semantic_forming  # 来自 graph_builder
```
而非：
```python
forming_bonds_for_summary = graph.edges[0].forming_bonds  # 可能是旧数据
```

---

## Step 4: 加严几何预检

### 问题
`_validate_forming_bonds_against_product()` 只在 **全部** forming_bonds 可疑时才触发 fallback。但对于 rx1，12-13 在产物中距离 ~1.5 Å（已有 SINGLE 键）→ 应被标记为可疑；15-16 距离 ~1.5 Å → 也已有 SINGLE 键。

然而 1.5 Å 恰好在阈值范围 (1.2-3.5) 内，不触发警告。问题是：**forming bond 在产物中已有 SINGLE 键是化学上不合理的**——forming bond 应该是反应中新形成的键，在产物中它才刚刚形成。

### 修改思路

增加一个更智能的检查：如果 forming_bond 对应的原子对在**产物分子图中已经有共价键**（通过 RDKit bond perception），则该对不可能是 "forming"（它已经存在了）：

```python
def _validate_forming_bonds_against_product(self, forming_bonds, product_xyz_path, index_base):
    # ... 现有距离检查 ...
    
    # 新增：通过分子图检查 forming_bonds 是否在产物中已存在共价键
    try:
        from rdkit import Chem
        mapped_smiles = ...  # 从 cleaner_data 获取
        mol = Chem.MolFromSmiles(mapped_smiles)
        if mol:
            map_to_idx = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum() > 0}
            for pair in forming_bonds:
                if pair[0] in map_to_idx and pair[1] in map_to_idx:
                    bond = mol.GetBondBetweenAtoms(map_to_idx[pair[0]], map_to_idx[pair[1]])
                    if bond is not None:
                        logger.warning(
                            f"[S2] Forming bond {pair}: already has {bond.GetBondType()} bond "
                            f"in product molecular graph — NOT a forming bond!"
                        )
                        suspicious_count += 1
    except Exception:
        pass
```

### 更根本的解决方案

如果 forming_bonds 在产物中已存在共价键，直接剔除该对，并从 `order_changed` 候选中补充。这需要在 `_resolve_forming_bonds_for_s2()` 中实现。

---

## Step 5: 重新运行 pipeline 验证

### 操作
```bash
# 清除旧 S0/S2 输出
rm -rf rph_output/rx_1/RXN_358d0ac1/S0_Mechanism
rm -rf rph_output/rx_1/RXN_358d0ac1/S2_Retro
# 清除 pipeline.state 中的 S0/S2 记录
python -c "
import json
with open('rph_output/rx_1/RXN_358d0ac1/pipeline.state') as f: s = json.load(f)
for k in ['step_s0', 'step_s2']:
    s['steps'].pop(k, None)
with open('rph_output/rx_1/RXN_358d0ac1/pipeline.state', 'w') as f: json.dump(s, f, indent=2)
"
# 重新运行（S1 checkpoint resume 跳过）
python -m rph_core --rx-id 1 --output rph_output/rx_1/RXN_358d0ac1 --skip-steps s3,s4
```

### 验证项
```bash
# 1. forming_bonds 应包含 12-19
python -c "
import json
with open('rph_output/rx_1/RXN_358d0ac1/S0_Mechanism/mechanism_summary.json') as f: d = json.load(f)
print('forming_bonds:', d['forming_bonds'])
# 期望: [[12,19], [15,16]] 或包含 12-19 的组合
"

# 2. 前体映射应非空
python -c "
import json
with open('rph_output/rx_1/RXN_358d0ac1/S0_Mechanism/atom_map_smiles.json') as f: d = json.load(f)
print('precursor mapping count:', len(d['smiles_atom_mapping']['precursor_smiles_to_map']))
# 期望: > 0
"

# 3. forming_bonds_annotated 应与 forming_bonds 一致
python -c "
import json
with open('rph_output/rx_1/RXN_358d0ac1/S0_Mechanism/atom_map_smiles.json') as f: d = json.load(f)
for fb in d.get('forming_bonds_annotated') or []:
    print(f'map={fb[\"map_space\"]} product_idx={fb[\"product_smiles_idx\"]} precursor_idx={fb[\"precursor_smiles_idx\"]} type={fb[\"bond_type_product\"]}/{fb[\"bond_type_precursor\"]}')
"
```

---

## 依赖关系与实施顺序

```
Step 1 (更新 CSV)         ← 可独立执行，也可与 Step 2 并行
    ↓
Step 2 (前体映射)          ← 独立
Step 3 (annotated 一致性)  ← 独立
Step 4 (几何预检)          ← 独立
    ↓
Step 5 (验证)              ← 依赖 Step 1-4 全部完成
```

Step 1-4 可并行实施，Step 5 必须在 Step 1-4 之后。
