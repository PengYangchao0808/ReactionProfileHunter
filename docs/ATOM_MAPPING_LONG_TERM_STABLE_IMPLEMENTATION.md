# LA/PEB 原子映射长期稳定实现方案

## 1. 当前最新结论

PEB 重写后，S2 已经不再需要 Kabsch + Hungarian、SubstructMatch、SMARTS 拓扑识别来重建完整原子对应。新的主链路是：

```text
S0 atom_map_smiles.json
  Map# <-> product SMILES idx

S1 smiles_to_xyz_map.json
  product SMILES idx -> product XYZ idx

S2 PEB
  Map# -> product SMILES idx -> product XYZ idx
```

验证脚本 `scripts/validate_peb_on_existing.py` 证明 PEB resolver 自身可以正确组合 sidecar，但也暴露了一个更底层的问题：

```text
S0 atom_map_smiles.json 使用 mapped product SMILES 的 RDKit atom idx
S1 smiles_to_xyz_map.json 使用 unmapped product SMILES 的 RDKit atom idx
```

复杂环系中 mapped SMILES 与 unmapped SMILES 的 RDKit 原子序号不等价，所以直接组合会得到错误结果：

```text
RXN_8ebf416f
PEB 初始错误结果: ((2, 3), (6, 7))   # C7-C8
期望正确结果:     ((2, 3), (6, 21))  # C7-C22
```

RDKit 2D 同构翻译证明：

```text
Map#3: mapped_idx=7 -> unmapped_idx=21 -> XYZ_idx=21
Map#8: mapped_idx=6 -> unmapped_idx=6  -> XYZ_idx=6
```

因此长期稳定实现的核心不是再换一个后验匹配方法，而是统一并显式声明所有索引空间。

## 2. 长期设计原则

1. **唯一 canonical product atom index space**

   PEB 使用的 `product_smiles_idx` 必须定义为 S1 几何生成所使用的 product molecule 的 RDKit atom idx。本文称为：

   ```text
   geometry_product_smiles_idx
   ```

   S0 写出的 `map_to_product_smiles` 必须已经翻译到这个空间。

2. **S0 不再暴露 mapped SMILES idx 作为默认 product_smiles_idx**

   mapped SMILES 的 RDKit idx 只能作为诊断字段：

   ```text
   mapped_product_smiles_idx
   ```

   不能被 PEB resolver 直接消费。

3. **S1 只消费 geometry product SMILES**

   S1 的 `smiles_to_xyz_map.json` 必须声明：

   ```json
   "product_smiles_idx_space": "geometry_product_smiles_idx"
   ```

4. **PEB resolver 必须检查 S0/S1 index_space 一致**

   若 S0 是 mapped idx、S1 是 unmapped idx，必须 fail-fast，不能静默组合。

5. **新计算从源头传递 atom identity**

   长期目标不是每次用 RDKit substructure match 翻译 mapped/unmapped，而是在 S0 生成 geometry product molecule 时就写出完整 atom identity table，并让 S1 直接使用该 molecule 或其确定性 SMILES。

## 3. 推荐分阶段实现

## Phase 1: 立即修复 S0 索引空间

### 3.1 修改位置

```text
rph_core/steps/mechanism_classifier/graph_builder.py
  GraphBuilder._build_smiles_mapping()
```

当前逻辑：

```python
product_smiles = self._get_product_mapping_smiles(record)
mol = Chem.MolFromSmiles(product_smiles)
for atom in mol.GetAtoms():
    map_num = atom.GetAtomMapNum()
    idx = atom.GetIdx()
    mapping.map_to_product_smiles[map_num] = idx
```

问题是 `product_smiles` 优先来自 mapped product SMILES，因此 `idx` 是 mapped product SMILES idx。

### 3.2 新逻辑

同时解析：

```text
mapped_product_smiles     # 带 Map#
geometry_product_smiles   # S1 实际用于几何生成的无 Map# product SMILES
```

建立翻译：

```text
mapped_product_idx -> geometry_product_smiles_idx
```

然后输出：

```text
map_to_product_smiles[Map#] = geometry_product_smiles_idx
product_smiles_to_map[geometry_product_smiles_idx] = Map#
```

保留诊断字段：

```json
"mapped_product_smiles_to_map": {"7": 3},
"map_to_mapped_product_smiles": {"3": 7},
"mapped_to_geometry_product_smiles_idx": {"7": 21},
"geometry_to_mapped_product_smiles_idx": {"21": 7}
```

### 3.3 翻译算法

建议新增内部 helper：

```python
def _build_mapped_to_geometry_product_idx(
    mapped_product_smiles: str,
    geometry_product_smiles: str,
) -> Dict[int, int]:
    ...
```

算法：

1. `mapped_mol = Chem.MolFromSmiles(mapped_product_smiles)`
2. `geometry_mol = Chem.MolFromSmiles(geometry_product_smiles)`
3. 复制 `mapped_mol`，清除 atom map number，得到 `query_mol`
4. `geometry_mol.GetSubstructMatches(query_mol, uniquify=False, useChirality=True)`
5. 若无匹配，用 `useChirality=False` 重试，并记录 warning
6. 若唯一匹配，直接使用
7. 若多匹配，使用以下消歧顺序：
   - 形成键 Map# 对应的 atoms 必须映射到 geometry product 中有对应 product bond 的 atom pair
   - 原子元素、degree、formal charge、aromatic、ring membership 必须一致
   - 手性中心 CIP/chiral tag 尽量一致
8. 若仍多解，不猜，fail-fast

注意：这里的 RDKit matching 只在同一 2D product SMILES 的 mapped/unmapped 之间建立索引翻译，不用于最终 XYZ 结构重建。风险远低于在 LA 复合物 XYZ 上做后验匹配，但仍必须唯一或可证明消歧。

### 3.4 S0 输出 schema

`atom_map_smiles.json` 升级为：

```json
{
  "schema_version": "3.1",
  "index_spaces": {
    "map_numbers": "1-based atom map numbers",
    "mapped_product_smiles_idx": "0-based RDKit atom index in mapped product SMILES",
    "geometry_product_smiles_idx": "0-based RDKit atom index in S1 geometry product SMILES"
  },
  "product_smiles_idx_space": "geometry_product_smiles_idx",
  "smiles_atom_mapping": {
    "product_smiles_to_map": {"21": 3, "6": 8},
    "map_to_product_smiles": {"3": 21, "8": 6},
    "mapped_product_smiles_to_map": {"7": 3, "6": 8},
    "map_to_mapped_product_smiles": {"3": 7, "8": 6},
    "mapped_to_geometry_product_smiles_idx": {"7": 21, "6": 6},
    "geometry_to_mapped_product_smiles_idx": {"21": 7, "6": 6}
  },
  "forming_bonds_map_space": [[1, 11], [3, 8]],
  "forming_bonds_index_space": "map_space"
}
```

## Phase 2: 强化 PEB resolver 的 index-space gate

### 4.1 修改位置

```text
rph_core/utils/atom_mapping.py
```

当前 `_extract_map_to_product_smiles()` 只读取 `map_to_product_smiles`，不验证它到底是哪种 SMILES idx。

### 4.2 新校验

`resolve_peb_forming_bonds()` 必须检查：

```text
S0 product_smiles_idx_space == S1 product_smiles_idx_space
```

推荐字段：

```json
S0:
"product_smiles_idx_space": "geometry_product_smiles_idx"

S1:
"product_smiles_idx_space": "geometry_product_smiles_idx"
```

如果缺字段：

- 新 schema 下直接 fail-fast
- 兼容模式仅允许验证脚本显式传 `allow_legacy_index_space=True`

### 4.3 禁止误用旧字段

若 S0 payload 中只有旧 `map_to_product_smiles`，且没有：

```json
"product_smiles_idx_space": "geometry_product_smiles_idx"
```

PEB resolver 必须拒绝：

```text
S0 map_to_product_smiles index space is undeclared; refusing to compose with S1 smiles_to_xyz_map
```

## Phase 3: S1 sidecar 真实化

### 5.1 当前问题

`write_smiles_to_xyz_map()` 目前会在验证通过时写：

```json
"mapping_source": "rdkit_generation_sidecar",
"confidence": "high"
```

但在验证脚本中，它是在已有 `product_min.xyz` 上后验重建 sidecar。即使元素序列匹配，也不等于真正追踪了坐标生成时的 atom identity。

### 5.2 长期改法

区分两种来源：

```text
rdkit_generation_sidecar
  坐标生成时同步写出，真实高置信

existing_xyz_reconstructed_sidecar
  对已有 XYZ 后验重建，只能用于迁移验证，不能默认 high
```

建议 `write_smiles_to_xyz_map()` 增加参数：

```python
def write_smiles_to_xyz_map(
    smiles: str,
    xyz_path: Path,
    output_dir: Path,
    *,
    atom_identity_source: str = "existing_xyz_reconstructed",
    product_smiles_idx_space: str = "geometry_product_smiles_idx",
    ...
) -> Path:
```

只有当 `atom_identity_source == "generation"` 时才允许：

```json
"mapping_source": "rdkit_generation_sidecar",
"confidence": "high"
```

验证脚本可以显式开启：

```text
allow_existing_xyz_reconstruction_for_validation = true
```

## Phase 4: 让 S1 使用 S0 生成的 geometry product SMILES

### 6.1 当前风险

S0 用一个 product SMILES 建映射，S1 可能用另一个 product SMILES 生成构象。即使两者表示同一分子，RDKit atom idx 也可能不同。

### 6.2 长期稳定路径

S0 输出：

```json
"geometry_product_smiles": "...",
"geometry_product_smiles_source": "mapped_product_smiles_with_maps_removed_and_atom_order_recorded",
"geometry_product_smiles_idx_space": "geometry_product_smiles_idx"
```

S1 必须读取并使用这个 exact `geometry_product_smiles`，不能重新从 cleaner raw product SMILES 选择另一个字符串。

如果 S1 因 LA/surrogate 需要扩展 SMILES，也必须保留 organic geometry product molecule 的 atom order，并把 LA additive 作为附加原子区间写入 sidecar。

## Phase 5: Atom identity table 作为最终形态

长期最稳的结构是统一写一个 atom identity table，贯穿 S0-S3：

```json
{
  "schema_version": "4.0",
  "atom_identity_space": "rph_product_atom_uid",
  "atoms": [
    {
      "atom_uid": "prod:00021",
      "map_num": 3,
      "mapped_product_smiles_idx": 7,
      "geometry_product_smiles_idx": 21,
      "s1_product_xyz_idx": 21,
      "element": "C",
      "role": "organic_heavy"
    }
  ]
}
```

S2 PEB 不再组合多个 dict，而是查询：

```text
Map# -> atom_uid -> s1_product_xyz_idx
```

每个阶段只做一件事：

- S0: 建立 `Map# -> atom_uid`
- S1: 建立 `atom_uid -> product_xyz_idx`
- S2: 建立 `forming_bonds_map_space -> forming_bonds_xyz_idx`
- S3: 继承并验证 atom order

## 7. Existing results 迁移方案

对已有 `Output/la_zncl2_as_mgcl2` 结果，不触发 DFT/QC。

### 7.1 验证脚本继续保留

```bash
python scripts/validate_peb_on_existing.py
```

但脚本应明确标记：

```json
"mapping_source": "existing_xyz_reconstructed_sidecar",
"validation_mode": true
```

### 7.2 回填 S0

对已有 `atom_map_smiles.json`：

1. 读取 mapped product SMILES
2. 读取 S1 sidecar 的 `source_smiles` 作为 geometry product SMILES
3. 建立 `mapped_to_geometry_product_smiles_idx`
4. 重写 `map_to_product_smiles` 到 geometry idx 空间
5. 保留旧 mapped idx 诊断字段

### 7.3 验证目标

`RXN_8ebf416f`：

```text
Map# [3,8] -> product_smiles_idx [21,6]
Map# [3,8] -> product_xyz_0based [21,6]
canonicalized forming bond pair -> (6,21)
```

`RXN_f5d5c7b9`：

```text
PEB result must be independent from kabsch_hungarian
```

## 8. Checkpoint 与重算策略

S2 checkpoint signature 必须包含：

- `S0_Mechanism/atom_map_smiles.json` hash
- `S1_ConfGeneration/smiles_to_xyz_map.json` hash
- `S1_ConfGeneration/product_min.xyz` hash
- `S2_Retro/atom_mapping_peb.json` hash
- `product_smiles_idx_space`
- `atom_mapping_schema_version`

修复后：

- `RXN_8ebf416f`: 如果无有效 S2，可直接跑 S2
- `RXN_f5d5c7b9`: 如果旧 S2 forming bonds 来自旧映射，checkpoint 应失效并重跑 S2
- 不需要重跑 S1 DFT，只需要重跑依赖 forming bonds 的 xTB scan/TS guess

## 9. 测试计划

### 9.1 GraphBuilder 单元测试

新增测试：

```text
tests/test_atom_mapping_index_space.py
```

覆盖：

- mapped product SMILES idx 与 geometry product SMILES idx 不同时，S0 输出 geometry idx
- `mapped_to_geometry_product_smiles_idx` 正确写出
- `product_smiles_idx_space == "geometry_product_smiles_idx"`
- `forming_bonds_index_space == "map_space"`

### 9.2 PEB resolver 单元测试

覆盖：

- S0/S1 index space 一致时可解析
- S0/S1 index space 不一致时 fail-fast
- S0 缺 `product_smiles_idx_space` 时 fail-fast
- 兼容模式必须显式开启

### 9.3 Existing results 回归

固定断言：

```text
RXN_8ebf416f:
  expected forming_bonds_product_xyz_0based == ((2,3), (6,21))
  forbidden == ((2,3), (6,7))

RXN_f5d5c7b9:
  PEB resolver 不读取 atom_map_xyz 的 kabsch_hungarian 结果
```

### 9.4 End-to-end dry run

新增无 QC dry-run：

```bash
python scripts/validate_peb_on_existing.py
pytest tests/test_atom_mapping_index_space.py -v
pytest tests/test_atom_mapping_peb.py -v
python scripts/ci/check_imports.py rph_core
```

## 10. 实施顺序

1. 在 S0 schema 中新增 index space 字段。
2. 修改 `GraphBuilder._build_smiles_mapping()`，将 Map# 映射翻译到 geometry product SMILES idx。
3. PEB resolver 增加 S0/S1 index space 一致性检查。
4. 更新 `write_smiles_to_xyz_map()`，区分 generation sidecar 与 existing XYZ reconstruction。
5. 更新 `validate_peb_on_existing.py`，回填时写出 mapped/unmapped 翻译表。
6. 更新 S2 checkpoint signature。
7. 加入 `RXN_8ebf416f` 和 `RXN_f5d5c7b9` 回归测试。
8. 后续再引入 schema 4.0 atom identity table，逐步替代分散 sidecar。

## 11. 成功标准

- S0/S1 sidecar 明确使用同一个 `geometry_product_smiles_idx` 空间。
- PEB resolver 不再能组合两个未声明或不一致的 index space。
- `RXN_8ebf416f` 稳定得到 `((2,3), (6,21))`。
- `RXN_f5d5c7b9` 不依赖 Kabsch/Hungarian 也能得到正确 PEB forming bonds。
- Existing results 验证不触发任何 DFT/QC。
- 新计算从 S0 到 S2 不需要后验完整原子映射猜测。

