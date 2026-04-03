# S2-S3 问题报告：形成键识别错误 & 中间体优化失败

**日期**: 2026-03-14 (更新)  
**分析范围**: `cleaner_adapter.py`, `dataset_loader.py`, `retro_scanner.py`, `qc_task_runner.py`, `defaults.yaml`

---

## 核心问题

1. **🔴 P0-新发现：形成键原子索引识别错误** — 拉伸/扫描的键并非 [4+3] 实际形成键
2. **🔴 P0：中间体几何优化直接坍缩到产物结构**

---

## 1. 形成键识别错误（P0 — 最高优先级）

### 1.1 实际数据验证

以 rx_9425282 为例：

- **产物 SMILES**: `O=C1C[C@@H]2C=C[C@@]3(CCN(C(=O)OCc4ccccc4)[C@@H]13)O2`
- **CSV `core_bond_changes`**: `3-4:formed;7-21:formed` ← 这是**原子映射编号 (atom map number)**
- **pipeline.state 中记录的 forming_bonds**: `[5, 10]` 和 `[8, 12]` ← 这是**0-based 原子索引**

### 1.2 根因分析：原子索引坐标系错配

数据流如下：

```
CSV core_bond_changes: "3-4:formed; 7-21:formed"
         │  (atom map numbers)
         ▼
cleaner_adapter.py: parse_formed_pairs_from_core_bond_changes()
         │  → [(3, 4), (7, 21)]  (atom map pairs)
         ▼
cleaner_adapter.py: map_pairs_to_internal_indices(mapped_smiles, ...)
         │  → 用 mapped_precursor_smiles 的 RDKit 原子顺序转换
         │  → 得到 (5, 10) 和 (8, 12)  ← ⚠️ 这是前体 SMILES 的原子序号！
         ▼
dataset_loader.py: raw["formed_bond_index_pairs"] = "5-10;8-12"
         │  raw["forming_bonds_index_base"] = "0"
         ▼
orchestrator._resolve_forming_bonds_for_s2()
         │  → 直接使用 (5, 10) 和 (8, 12) 作为产物 XYZ 的原子索引
         ▼
retro_scanner.run(): stretch_bonds(coords, [((5,10), 3.2), ((8,12), 3.2)])
         ← ⚠️ 错误！拉伸的是产物 XYZ 中第5-10、8-12号原子之间的键
            但这些索引来自前体 SMILES，不是产物 XYZ 的原子编号
```

**问题本质**: `map_pairs_to_internal_indices()` 将 atom map 编号转换为**前体 (precursor) SMILES** 的 RDKit 原子索引。但 S2 使用这些索引去操作**产物 XYZ 文件**，而产物 XYZ 的原子顺序与前体 SMILES 的 RDKit 原子顺序**完全不同**。

### 1.3 代码位置

| 文件 | 位置 | 问题 |
|------|------|------|
| `cleaner_adapter.py` L143-164 | `map_pairs_to_internal_indices()` | 用**前体** SMILES 做 map→index 转换 |
| `dataset_loader.py` L92-105 | `_enrich_cleaner_metadata()` | `mapped_smiles` 取自 precursor，非 product |
| `dataset_loader.py` L119 | `formed_index_pairs = map_pairs_to_internal_indices(mapped_smiles, ...)` | 转换结果是前体的原子序号 |
| `orchestrator.py` L449-496 | `_resolve_forming_bonds_for_s2()` | 直接使用前体索引去操作产物 XYZ |

### 1.4 正确做法

应使用**产物 SMILES (product_smiles_main)** 的 atom map → RDKit index 映射，或直接使用**产物 XYZ 的 atom map 对应表**来转换。具体需要：

1. 从 CSV 的 `rxn_smiles_mapped` 字段中提取产物端的 mapped SMILES
2. 用产物 mapped SMILES 构建 `atom_map → RDKit index` 映射
3. 将 RDKit index 与产物 XYZ 的原子顺序对齐（这需要 SMILES→XYZ 的原子对应关系）

> **注意**: 即使使用产物 SMILES 的 RDKit 索引，也不一定与 S1 输出的 XYZ 原子顺序一致。S1 经过 CREST 构象搜索和 DFT 优化后，原子顺序可能被重排。最可靠的方法是通过 **SMILES→3D 原子匹配** 来建立映射。

---

## 2. 中间体优化缺少约束冻结（P0）

### 现状

`retro_scanner.py` `_optimize_intermediate()` 执行**无约束** B3LYP/def2-SVP 优化，形成键原子完全自由 → 坍缩到产物。

### 应有逻辑

1. xTB 约束优化（冻结形成键 ~3.5Å，松弛其他自由度）
2. DFT 约束优化（modredundant 冻结形成键）
3. SP 确认中间体能量

---

## 3. 其他问题

### P1: 键拉伸距离 3.2Å 偏短 🟡

`defaults.yaml` `scan_start_distance = 3.2`，但 `[5+2]_default` 使用 `3.5`。建议统一为 ~3.5Å。

### P1: BondStretcher 默认参数混乱 🟡

`StretchingParams.target_length_A = 2.2` 从未使用，`RetroScanner` 使用 `scan_start_distance`。

### P2: S3 对 ts_guess 质量无预检 🟡

无检查形成键距离范围、产物 RMSD、能量梯度方向。

---

## 4. 根因链（更新）

```
[P0-根因] 形成键索引来自前体 SMILES → 操作的是错误的原子对
    ↓
拉伸了错误的键 → seed_reactant_like.xyz 化学意义错误
    ↓
[P0-根因] _optimize_intermediate 无约束 → 坍缩到产物
    ↓
xTB 扫描起点错误 → ts_guess.xyz 质量差
    ↓
S3 Berny/QST2 优化失败
```

---

## 5. 涉及文件（更新）

| 文件 | 修改内容 |
|------|----------|
| **`cleaner_adapter.py`** | `map_pairs_to_internal_indices()` 应使用**产物** mapped SMILES + SMILES↔XYZ 原子匹配 |
| **`dataset_loader.py`** | `_enrich_cleaner_metadata()` 中 `mapped_smiles` 应取产物端，非前体端 |
| `retro_scanner.py` | `_optimize_intermediate()` 增加 xTB 约束预优化 + DFT 约束冻结 |
| `qc_interface.py` | `XTBInterface.optimize()` 支持 distance constraint |
| `defaults.yaml` | `scan_start_distance → 3.5`, 增加 `freeze_forming_bonds: true` |

---

## 6. 修复优先级（更新）

| 级别 | 问题 | 影响 |
|------|------|------|
| **P0** | 形成键索引映射：前体→产物坐标系修正 | **所有反应的形成键均可能错误** |
| **P0** | `_optimize_intermediate` 增加约束冻结 | 阻止中间体坍缩 |
| **P0** | 增加 xTB 预松弛步骤 | 提供合理 DFT 起点 |
| **P1** | 拉伸距离改为 3.5Å | 远离产物势能阱 |
| **P1** | 中间体结构验证 | 检测坍缩并报错 |
| **P2** | 统一 StretchingParams 与 scan config | 消除混乱 |
| **P2** | S3 增加 ts_guess 质量预检 | 节省计算资源 |
