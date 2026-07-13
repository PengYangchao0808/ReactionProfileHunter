# 原子映射模块重写与 S2 PEB 方案

## 1. 背景与目标

当前 S2 的形成键解析仍依赖后验猜测链路：

```text
P1 atom_map_xyz full notation
P2 S0/S1 atom-map annotations
P3 mechanism_graph edges + atom_map_xyz
P4 S0 mechanism_summary legacy fallback
P5 SMARTS topology fallback
```

在引入 Lewis acid 后，Kabsch + Hungarian、SubstructMatch、SMARTS 拓扑识别都无法稳定处理复杂体系中的完整原子对应。典型失败是 `RXN_8ebf416f` 中正确断裂边应为 `C7-C22`，但后验匹配把它解析成 `C7-C8`。根因不是单个阈值，而是这些方法都在最终 LA 复合物 XYZ 上重新猜测 Map# 与 XYZ 行号的关系。

重写目标：

- 移除 S2 中 P1/P2/P3 这三条依赖 `atom_map_xyz.json` 的旧解析路径。
- 废弃 Kabsch + Hungarian、SubstructMatch、`_smilesAtomOutputOrder`、SMARTS 作为主映射来源。
- 直接组合现有完整链路：
  - S0 `atom_map_smiles.json`: `Map Number <-> product SMILES idx`
  - S1 `smiles_to_xyz_map.json`: `product SMILES idx -> product XYZ idx`
- 在 S2 阶段执行 PEB（Product Edge Breaking）时，唯一主流程为：

```text
forming_bonds_map_space
  -> map_to_product_smiles
  -> smiles_to_xyz_map
  -> forming_bonds_product_xyz_0based
  -> RetroScanner/xTB scan constraints
```

核心原则：S2 不再重建原子映射，只消费 S0/S1 明确写出的 atom identity sidecar。

## 2. 新数据契约

### 2.1 S0: `atom_map_smiles.json`

S0 必须明确声明 `forming_bonds` 的 index space 是 Map#，不要再让调用方猜。

推荐 schema：

```json
{
  "schema_version": "3.0",
  "index_spaces": {
    "map_numbers": "1-based atom map numbers from mapped SMILES",
    "product_smiles_idx": "0-based RDKit atom index in mapped product SMILES"
  },
  "smiles_atom_mapping": {
    "product_smiles_to_map": {"7": 3, "6": 8},
    "map_to_product_smiles": {"3": 7, "8": 6}
  },
  "forming_bonds_map_space": [[1, 11], [3, 8]],
  "forming_bonds_annotated": [
    {
      "map_space": [3, 8],
      "product_smiles_idx": [7, 6],
      "bond_type_product": "SINGLE",
      "bond_type_precursor": "NONE"
    }
  ]
}
```

兼容期可以继续写旧字段 `forming_bonds`，但必须同时写：

```json
"forming_bonds_index_space": "map_space"
```

### 2.2 S1: `smiles_to_xyz_map.json`

S1 必须作为主映射 sidecar，不只是诊断文件。

推荐 schema：

```json
{
  "schema_version": "3.0",
  "mapping_source": "rdkit_generation_sidecar",
  "confidence": "high",
  "source_smiles": "...",
  "xyz_file": "product_min.xyz",
  "xyz_sha256": "...",
  "index_spaces": {
    "product_smiles_idx": "0-based RDKit atom index",
    "product_xyz_idx": "0-based full XYZ row index"
  },
  "atoms": [
    {
      "smiles_idx": 7,
      "xyz_idx": 21,
      "element": "C",
      "type": "organic_heavy"
    }
  ],
  "lewis_acid": {
    "enabled": true,
    "organic_atom_count": 23,
    "additive_atom_indices": [42, 43, 44],
    "additive_elements": ["Mg", "Cl", "Cl"]
  },
  "validation": {
    "atom_count_matches_xyz": true,
    "element_sequence_matches": true,
    "no_smiles_idx_on_additive": true
  }
}
```

关键要求：

- `smiles_idx -> xyz_idx` 必须来自坐标生成时的原子身份传递。
- 如果当前实现仍只是按重原子顺序补写，必须先标记为 `mapping_source="assumed_heavy_order"`、`confidence="low"`，不得进入 S2 PEB 主流程。
- LA 原子必须 `smiles_idx = null`，`type = "additive"`，永不参与 Map#。

### 2.3 S2: `atom_mapping_peb.json`

S2 解析完成后写一个新的 PEB 映射产物，作为 S2/S3/外部后处理的唯一来源。

```json
{
  "schema_version": "3.0",
  "mapping_source": "s0_smiles_map_plus_s1_smiles_xyz_sidecar",
  "confidence": "high",
  "product_xyz": "S1_ConfGeneration/product_min.xyz",
  "forming_bonds": {
    "map_space": [[1, 11], [3, 8]],
    "product_smiles_idx": [[2, 3], [7, 6]],
    "product_xyz_0based": [[2, 3], [6, 21]],
    "product_xyz_1based": [[3, 4], [7, 22]]
  },
  "validation": {
    "all_map_numbers_resolved": true,
    "all_smiles_indices_resolved": true,
    "all_xyz_indices_in_range": true,
    "no_additive_atoms_in_forming_bonds": true,
    "forming_bond_atom_count": 4
  }
}
```

## 3. 新模块设计

新增文件：

```text
rph_core/utils/atom_mapping.py
```

职责：

- 读取并校验 `atom_map_smiles.json`
- 读取并校验 `smiles_to_xyz_map.json`
- 组合 Map# -> SMILES idx -> XYZ idx
- 为 S2 PEB 输出 forming bonds
- 写入 `atom_mapping_peb.json`

建议 API：

```python
from pathlib import Path
from typing import Dict, List, Tuple
from dataclasses import dataclass

@dataclass(frozen=True)
class PebMappingResult:
    forming_bonds_map_space: Tuple[Tuple[int, int], ...]
    forming_bonds_product_smiles: Tuple[Tuple[int, int], ...]
    forming_bonds_product_xyz_0based: Tuple[Tuple[int, int], ...]
    forming_bonds_product_xyz_1based: Tuple[Tuple[int, int], ...]
    confidence: str
    diagnostics: Dict[str, object]

def load_s0_atom_map_smiles(path: Path) -> Dict[str, object]:
    ...

def load_s1_smiles_to_xyz_map(path: Path, product_xyz: Path) -> Dict[str, object]:
    ...

def build_product_smiles_to_xyz(smiles_to_xyz_payload: Dict[str, object]) -> Dict[int, int]:
    ...

def resolve_peb_forming_bonds(
    *,
    atom_map_smiles_path: Path,
    smiles_to_xyz_map_path: Path,
    product_xyz_path: Path,
) -> PebMappingResult:
    ...

def write_peb_mapping(result: PebMappingResult, output_path: Path) -> Path:
    ...
```

## 4. PEB 解析算法

### 4.1 输入

必需输入：

- `S0_Mechanism/atom_map_smiles.json`
- `S1_ConfGeneration/smiles_to_xyz_map.json`
- `S1_ConfGeneration/product_min.xyz`

可选输入：

- `S0_Mechanism/mechanism_summary.json`
- `S0_Mechanism/mechanism_graph.json`

可选输入只能用于补齐 `forming_bonds_map_space`，不能直接作为 XYZ 索引来源。

### 4.2 步骤

1. 读取 S0 `atom_map_smiles.json`。
2. 提取 `map_to_product_smiles`。
3. 提取形成键，优先级：
   - `forming_bonds_map_space`
   - `forming_bonds_annotated[].map_space`
   - 兼容期旧字段 `forming_bonds`，但必须要求 `forming_bonds_index_space == "map_space"` 或来自 S0 文件。
4. 读取 S1 `smiles_to_xyz_map.json`。
5. 从 `atoms[]` 构建 `product_smiles_idx -> product_xyz_idx`。
6. 对每条 Map# 形成键执行：

```text
left_map -> left_product_smiles_idx -> left_product_xyz_idx
right_map -> right_product_smiles_idx -> right_product_xyz_idx
```

7. 生成 0-based XYZ forming bonds，给 S2 RetroScanner 使用。
8. 生成 1-based XYZ forming bonds，只用于 JSON 可读性和人工检查。
9. 写入 `S2_Retro/atom_mapping_peb.json`。

### 4.3 失败策略

PEB 是主流程，不做猜测式 fallback。

必须 fail-fast 的情况：

- 缺少 S0 `atom_map_smiles.json`
- 缺少 S1 `smiles_to_xyz_map.json`
- `smiles_to_xyz_map.json` 的 `confidence != "high"`
- `mapping_source` 是 `assumed_heavy_order`、`substructure_match`、`kabsch_hungarian`、`legacy`
- 任一 Map# 找不到 product SMILES idx
- 任一 product SMILES idx 找不到 XYZ idx
- XYZ idx 越界
- forming bond 指向 additive/H 原子
- 对需要两条 forming bonds 的流程，解析结果不是两条键

可以 warning 但继续的情况：

- 旧字段仍存在但新字段完整
- 1-based 辅助字段缺失
- `mechanism_summary.json` 仍写旧字段 `forming_bonds`

## 5. Orchestrator 修改方案

### 5.1 移除 P1/P2/P3

修改 `ReactionProfileHunter._resolve_forming_bonds_for_s2()`。

删除或停用以下解析器：

- `_resolve_forming_bonds_from_xyz_mapping_artifact`
- `_resolve_forming_bonds_from_mapping_annotations`
- `_resolve_forming_bonds_from_graph_edges_and_xyz_mapping`

保留这些函数一段时间可以用于兼容测试，但不再被 S2 主流程调用。

新流程：

```text
_resolve_forming_bonds_for_s2
  -> resolve_peb_forming_bonds(...)
  -> _finalize_forming_bonds_for_s2(...)
  -> return product_xyz_0based forming_bonds
```

### 5.2 禁止 S0 legacy XYZ fallback

当前旧逻辑会把 `mechanism_summary.json` 的 `forming_bonds` 当作 0-based XYZ 索引。这必须删除。

新规则：

- S0 `forming_bonds` 只允许解释为 Map#。
- S0 不能直接返回 S2 xTB scan constraints。
- 没有 S1 `smiles_to_xyz_map.json` 时直接报错。

### 5.3 `_build_and_save_xyz_mapping()` 降级

`atom_map_xyz.json` 可以暂时保留给旧工具读取，但它必须由 PEB 主链路生成，而不是 Kabsch/SubstructMatch 生成。

修改方向：

```text
_build_and_save_xyz_mapping
  -> load atom_map_smiles.json
  -> load smiles_to_xyz_map.json
  -> compose map_to_product_xyz_1based
  -> write atom_map_xyz.json with source="s0_s1_sidecar_composition"
```

不再调用：

- `_build_mol_to_xyz_index_mapping`
- `get_map_to_xyz_dict`
- Kabsch + Hungarian
- SubstructMatch

### 5.4 S2 checkpoint signature

S2 checkpoint signature 必须加入：

- `S0_Mechanism/atom_map_smiles.json` hash
- `S1_ConfGeneration/smiles_to_xyz_map.json` hash
- `S1_ConfGeneration/product_min.xyz` hash
- `S2_Retro/atom_mapping_peb.json` hash
- `atom_mapping_schema_version`

否则修复后旧 checkpoint 可能继续复用错误 forming bonds。

## 6. Conformer/S1 修改方案

### 6.1 强化 `write_smiles_to_xyz_map()`

当前实现的风险在于按 XYZ 重原子顺序分配 SMILES idx。重写时必须区分两种来源：

1. `rdkit_generation_sidecar`: 坐标生成时直接携带 atom identity，可 high confidence。
2. `assumed_heavy_order`: 后验按元素/顺序推断，只能 low confidence。

S2 PEB 只接受第一种。

### 6.2 LA 追加规则

LA 原子必须在 sidecar 中清楚标记：

```json
{
  "smiles_idx": null,
  "xyz_idx": 42,
  "element": "Mg",
  "type": "additive"
}
```

验证规则：

- `additive_atom_indices` 必须等于 XYZ 中 LA 原子实际行号。
- `organic_atom_count + hydrogen_count + additive_count == xyz_atom_count`
- additive 原子不得出现在 `smiles_idx -> xyz_idx` 主映射中。

## 7. S2 PEB 输出与 RetroScanner 接口

S2 传给 `RetroScanner` 的 forming bonds 必须是：

```text
product_xyz_0based
```

示例：

```python
forming_bonds = ((2, 3), (6, 21))
```

对应 JSON：

```json
"product_xyz_1based": [[3, 4], [7, 22]]
```

`RetroScanner` 不再接收 Map#、SMILES idx 或 legacy organic-heavy idx。

## 8. 验证规则

PEB 解析后必须做以下校验：

- `len(forming_bonds) == 2`
- 两条键共 4 个唯一原子
- 每个 XYZ idx 在 `product_min.xyz` 范围内
- 每个 XYZ idx 指向 organic heavy atom
- 不涉及 H、Mg、Cl、Al、Zn、Li 等 additive 或氢
- 元素与 S0 `product_smiles_idx` 对应元素一致
- 如果 `bond_type_product == "SINGLE"`，product XYZ 中距离应位于合理 C-C/C-X 成键区间
- 对 LA case，forming bonds 不得包含 `la_additive.additive_atom_indices`

注意：距离校验只是辅助，不能用它重新选择原子。

## 9. 回归测试计划

### 9.1 单元测试

新增：

```text
tests/test_atom_mapping_peb.py
```

测试点：

1. `map_to_product_smiles + smiles_to_xyz_map` 可以正确组合 Map# -> XYZ。
2. 缺少 Map# 时 fail-fast。
3. 缺少 SMILES idx 时 fail-fast。
4. 指向 additive/H 时 fail-fast。
5. `mapping_source="assumed_heavy_order"` 时 S2 拒绝使用。
6. 旧 `forming_bonds` 字段没有 `forming_bonds_index_space` 时不允许当 XYZ idx 使用。

### 9.2 案例回归

`RXN_8ebf416f`：

```text
Map# [3,8] -> XYZ [7,22] 1-based
Map# [3,8] -> XYZ [6,21] 0-based
禁止 XYZ [7,8] / 0-based [6,7]
```

`RXN_f5d5c7b9`：

```text
Map# [3,12] -> 原当前正确 XYZ pair
不能依赖 kabsch_hungarian 才成功
```

### 9.3 Checkpoint 回归

修改 `smiles_to_xyz_map.json` 后，S2 checkpoint 必须失效并重算 forming bonds。

## 10. 迁移步骤

### Phase 0: 文档与测试先行

- 新增本方案文档。
- 新增 `tests/test_atom_mapping_peb.py`，先写 red tests。
- 固定 `RXN_8ebf416f` 的 C7-C22 断言。

### Phase 1: 新建 PEB 解析模块

- 新增 `rph_core/utils/atom_mapping.py`。
- 实现 loader、validator、resolver、writer。
- 不改 orchestrator 主流程，只跑单元测试。

### Phase 2: 接入 S2

- 修改 `_resolve_forming_bonds_for_s2()`，优先且唯一使用 PEB resolver。
- 删除 P1/P2/P3 调用。
- 禁止 S0 legacy fallback 返回 XYZ idx。
- S2 生成 `S2_Retro/atom_mapping_peb.json`。

### Phase 3: 重写 `atom_map_xyz.json` 生成

- `_build_and_save_xyz_mapping()` 不再调用 Kabsch/SubstructMatch。
- 改为从 S0/S1 sidecar 组合生成。
- provenance 写 `mapping_source="s0_s1_sidecar_composition"`。

### Phase 4: S1 sidecar 强化

- `write_smiles_to_xyz_map()` 增加 `mapping_source/confidence/validation`。
- 只有坐标生成时传递 atom identity 的 sidecar 可以 high confidence。
- 后验顺序推断全部 low confidence。

### Phase 5: 清理旧代码

- `_build_mol_to_xyz_index_mapping()` 标记 deprecated 或删除。
- `get_map_to_xyz_dict()` 不再被 S2 映射链路调用。
- SMARTS fallback 只保留为诊断工具，不作为主 PEB 映射来源。

## 11. 成功标准

完成后必须满足：

- S2 不再调用 P1/P2/P3。
- S2 forming bonds 只来自 S0/S1 sidecar 组合。
- `RXN_8ebf416f` 输出 C7-C22，不再输出 C7-C8。
- `RXN_f5d5c7b9` 不依赖 Kabsch 也能得到正确 forming bonds。
- LA 原子不可能进入 forming bonds。
- 修改 S0/S1 映射文件会触发 S2 checkpoint 失效。
- `scripts/ci/check_imports.py rph_core` 通过。
- `pytest tests/test_atom_mapping_peb.py -v` 通过。

## 12. 非目标

本次重写不做：

- 不恢复 S4 功能，本仓库仍只负责 S0-S3。
- 不用 ML/SMARTS 猜形成键作为主流程。
- 不在最终 LA 复合物 XYZ 上重新推断完整 Map#。
- 不放宽 Kabsch 阈值作为修复手段。

