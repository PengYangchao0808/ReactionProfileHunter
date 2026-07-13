# Lewis Acid 原子映射问题排查报告

**项目**: ReactionProfileHunter v3.0.1  
**日期**: 2026-06-29  
**状态**: 排查完成，待修复  
**排查范围**: 路易斯酸（LiCl/MgCl₂/AlCl₃）引入后的原子索引映射全链路（S0→S1→S2→S3）

---

## 0. 执行摘要

经过对原子映射全链路的逐函数排查，共发现 **9 个确认问题**，其中 **4 个为 CRITICAL（会导致计算结果错误或静默失败）**，**2 个为 HIGH**，**3 个为 MEDIUM/LOW**。

**核心结论**：当前实现**仅对 LiCl（2 原子） surrogate 有效**；一旦切换到配置中已声明的 MgCl₂（3 原子）或 AlCl₃（4 原子），多处硬编码假设会崩溃。此外，S1 阶段的 LA 几何校验因索引时序问题**完全失效（空跑）**。

| 编号 | 严重度 | 位置 | 一句话描述 |
|------|--------|------|-----------|
| #1 | 🔴 CRITICAL | `engine.py:877-907` | S1 阶段 LA 几何校验空跑（`additive_atom_indices` 为空） |
| #2 | 🔴 CRITICAL | `engine.py:524,559` | 哑键创建/删除硬编码 `la_atoms[0]` 和 `N-2`，MgCl₂/AlCl₃ 下键连错误原子 |
| #3 | 🔴 CRITICAL | `engine.py:2860-2866` | `_strip_la_suffix` 正则无法处理分支 SMILES，AlCl₃ 后缀不剥离 |
| #4 | 🔴 CRITICAL | `engine.py:2922-2938` | `write_smiles_to_xyz_map` additive 计数硬编码检测 "Li"，MgCl₂/AlCl₃ 检测失败 |
| #5 | 🟠 HIGH | 全局 | 无显式 `organic_atom_count` 字段，有机/LA 边界到处重复推导 |
| #6 | 🟠 HIGH | `smarts_matcher.py:166,205` | SMARTS Mol（仅重原子）接收完整 coords（重+H+LA），索引空间隐式耦合 |
| #7 | 🟡 MEDIUM | `orchestrator.py:2599-2610` | additive_atom_indices 位置推导，不匹配时仅 warning 不 fail-fast |
| #8 | 🟡 MEDIUM | `orchestrator.py:245-246` | charge/multiplicity 存储但未传入 QC 输入生成 |
| #9 | 🟢 LOW | `lewis_acid/models.py` | `LiCoordinationTrace` 字段名保留 `li_` 前缀，与通用金属语义不符 |

---

## 1. 原子映射全链路追踪

下图展示了引入 LA 后，原子索引从 SMILES 到最终 TS 结构的完整流转路径。**绿色**表示索引空间一致；**红色**表示存在断裂风险的位置。

```
SMILES (canonical, 不含 LA)
  │
  │  orchestrator._prepare_lewis_acid_system()  [orchestrator.py:193-255]
  │  拼接: geometry_smiles = f"{product_smiles}.{la_smiles}"
  │  创建: LewisAcidAdditive(additive_atom_indices=())  ← 初始为空!
  ▼
geometry_smiles (含 LA 后缀, 如 "C=CCO.[Li]Cl")
  │
  │  ════════════ S1: ConformerEngine ════════════
  │
  │  engine._step_rdkit_embed()  [engine.py:460-592]
  │  ├─ [474-480] _strip_la_suffix → organic_smiles          ← 问题 #3
  │  ├─ [482-485] MolFromSmiles(organic) + AddHs
  │  ├─ [511-528] AddAtom(LA) + AddBond(O, la_atoms[0])      ← 问题 #2
  │  ├─ [537-547] ETKDG 3D embedding
  │  ├─ [556-564] RemoveBond(O, N-2)  ← 硬编码 N-2            ← 问题 #2
  │  ├─ [567]    MMFF optimize
  │  └─ [583-589] write_xyz → 原子序: [重原子(SMILES序), H, LA]
  ▼
_init.xyz  → xTB/CREST → product_min.xyz  (原子序保持)
  │
  │  engine._validate_la_geometry()  [engine.py:877-907]      ← 问题 #1
  │  additive_atom_indices 仍为 () → 校验直接 return True (空跑!)
  │
  │  ════════════ S1 结束 → orchestrator 推导索引 ════════════
  │
  │  orchestrator [2590-2612]
  │  读取 product_min.xyz → additive_atom_indices = (N-2, N-1)  ← 问题 #7
  │  假设: LA 原子总在 XYZ 末尾，仅校验末尾 N 个元素
  │
  │  engine.write_smiles_to_xyz_map()  [engine.py:2869-3003]   ← 问题 #4
  │  写出 smiles_to_xyz_map.json
  ▼
product_min.xyz (含 LA) + additive_atom_indices (已填充)
  │
  │  ════════════ S2: RetroScanner + SMARTSMatcher ════════════
  │
  │  smarts_matcher.find_reactive_bonds()  [smarts_matcher.py:131-254]
  │  ├─ [161] MolFromSmiles(canonical) → 仅有机重原子 Mol
  │  ├─ [166] _embed_coords_into_mol(mol, coords)             ← 问题 #6
  │  │         coords 含 全部原子(重+H+LA), Mol 仅含重原子
  │  │         按 i in range(N_heavy) 取 coords[i] → 恰好正确
  │  ├─ [205] identify(mol, coords)  ← coords 仍是完整数组
  │  └─ [232-251] LA 过滤: 形成键含 additive → 拒绝
  ▼
forming_bonds (0-based, 有机重原子索引空间)
  │
  │  retro_scanner.run_retro_scan()
  │  ├─ forming_bonds 直接用作 XYZ 索引（含 LA 的完整 XYZ）
  │  │  因有机重原子在 XYZ 前部 → 索引恰好一致            ← 隐式耦合
  │  ├─ xTB scan 用 forming_bonds 做距离约束
  │  └─ 输出 ts_guess.xyz, intermediate.xyz (含全部原子)
  ▼
ts_guess.xyz + intermediate.xyz + forming_bonds
  │
  │  ════════════ S3: TSOptimizer ════════════
  │
  │  ts_optimizer.run()
  │  ├─ 读取 ts_guess/intermediate (含 LA)
  │  ├─ forming_bonds 传入优化器
  │  └─ validator.validate_la_mode_projection()  [validator.py:109-208]
  │     使用 additive_atom_indices 排除 LA 振动模式
  ▼
ts_final.xyz (含 LA)
```

**索引不变量（全链路隐式依赖）**：
```
XYZ 原子排列: [有机重原子 (SMILES 顺序)] [有机 H 原子] [LA 原子]
索引区间:      0 ... H₀-1               H₀ ... O₀-1    O₀ ... N-1
```
- LA 原子**始终**在 XYZ 末尾
- `additive_atom_indices = (N - L, N - L + 1, ..., N - 1)`，其中 L = LA 原子数
- forming_bonds 索引**始终**指向有机重原子区间 `[0, H₀)`
- **无显式 `organic_atom_count` 字段**——每个消费方独立重推边界

---

## 2. 确认问题详解

### 🔴 问题 #1：S1 阶段 LA 几何校验空跑

**位置**: `rph_core/steps/conformer_search/engine.py:877-907` + 调用点 `638, 701`

**现象**: `_validate_la_geometry()` 在 S1 conformer 搜索过程中被调用两次（best_xyz 和 gfn0_ensemble 阶段），但此时 `la_additive.additive_atom_indices` 仍为空元组 `()`（该字段直到 S1 完成后在 `orchestrator.py:2600` 才填充）。

**代码证据**:
```python
# engine.py:882-884
additive_indices = la_additive.additive_atom_indices   # ← 此时 = ()
if coords is None or not additive_indices or len(coords) <= max(additive_indices):
    return True   # ← not () == True → 直接返回，校验被跳过!
```

**影响**: 设计文档（§7.3）要求的 S1 后几何校验（Li–Cl ≤ 3.5 Å, Li–O ≤ 4.0 Å）**从未执行**。LiCl 在 CREST 过程中解离也无法被检测，导致下游 S2/S3 在错误几何上运行。

**设计意图 vs 实际**:
| 阶段 | 设计要求 | 实际行为 |
|------|---------|---------|
| S1 conformer 搜索中 | 校验每个候选构型的 LA 配位 | ❌ 空跑（indices 为空） |
| S1 完成后 (orchestrator) | 推导 indices | ✅ 正确推导 |
| S2 quality gate | 使用 indices 校验 scan 轨迹 | ✅ 可用（indices 已填充） |

**根因**: 时序设计缺陷——indices 在 S1 **完成后**推导，但校验需要在 S1 **进行中**执行。

---

### 🔴 问题 #2：哑键创建/删除硬编码，多原子 surrogate 键连错误

**位置**: `rph_core/steps/conformer_search/engine.py:524, 559`

**现象**: LA 原子添加时，创建 dummy bond 连接羰基 O 和 `la_atoms[0]`；embedding 后删除该 bond 时硬编码 `GetNumAtoms() - 2`。这两处假设 **surrogate 第一个原子是金属** 且 **surrogate 恰好 2 个原子**。

**代码证据**:
```python
# engine.py:524 — dummy bond 创建
rw_mol.AddBond(carbonyl_o_idx, la_atoms[0], Chem.BondType.SINGLE)
#                                     ^^^^^^^^^^
# la_atoms[0] = surrogate SMILES 的第一个原子

# engine.py:559 — dummy bond 删除
li_idx = rw_mol.GetNumAtoms() - 2  # ← 硬编码：假设金属在倒数第二
rw_mol.RemoveBond(carbonyl_o_idx, li_idx)
```

**各 surrogate 的表现**:

| Surrogate | SMILES | 原子序 | `la_atoms[0]` | dummy bond 连接 | 删除目标 (N-2) | 结果 |
|-----------|--------|--------|---------------|----------------|---------------|------|
| LiCl | `[Li]Cl` | Li(0), Cl(1) | **Li** ✓ | O–Li ✓ | Li (N-2) ✓ | ✅ 正确 |
| MgCl₂ | `Cl[Mg]Cl` | Cl(0), Mg(1), Cl(2) | **Cl** ✗ | O–**Cl** ✗ | Mg (N-2) ✗ | ❌ 错误 |
| AlCl₃ | `Cl[Al](Cl)Cl` | Cl(0), Al(1), Cl(2), Cl(3) | **Cl** ✗ | O–**Cl** ✗ | Cl (N-2) ✗ | ❌ 错误 |

**MgCl₂ 具体失败链**:
1. `la_atoms = [N, N+1, N+2]` 对应 `[Cl, Mg, Cl]`
2. `AddBond(O, la_atoms[0])` → 创建 **O–Cl** 键（应为 O–Mg）
3. embedding 后 `RemoveBond(O, N-2)` → 尝试删除 **O–Mg** 键（不存在！）
4. O–Cl dummy 键**残留在结构中** → CREST/DFT 看到错误的 connectivity
5. `Chem.SanitizeMol` 可能报错或静默产生错误拓扑

**影响**: MgCl₂/AlCl₃ surrogate 完全不可用——RDKit embedding 产生错误拓扑，后续所有几何都是垃圾。

---

### 🔴 问题 #3：`_strip_la_suffix` 正则无法处理分支 SMILES

**位置**: `rph_core/steps/conformer_search/engine.py:2860-2866`

**代码**:
```python
def _strip_la_suffix(smiles: str) -> str:
    match = re.search(r'\.(?:\[[A-Z][a-z]?\]|[A-Z][a-z]?)+$', smiles)
    if match:
        return smiles[:match.start()]
    return smiles
```

**正则分析**: `\.(?:AtomToken)+$` 要求点号后跟随连续的原子 token 直到字符串末尾。

| Surrogate SMILES | 后缀 | 正则匹配 | 剥离结果 |
|-----------------|------|---------|---------|
| `.[Li]Cl` | `[Li]` + `Cl` | ✅ 全部匹配 | ✅ 正确剥离 |
| `.Cl[Mg]Cl` | `Cl` + `[Mg]` + `Cl` | ✅ 全部匹配 | ✅ 正确剥离 |
| `.Cl[Al](Cl)Cl` | `Cl` + `[Al]` + **`(`** ← 断裂 | ❌ 遇 `(` 无法匹配 | ❌ **不剥离** |

**AlCl₃ 失败链**:
1. `_strip_la_suffix("C=CCO.Cl[Al](Cl)Cl")` → 正则在 `(` 处中断，无法匹配到 `$`
2. 返回原始 SMILES `"C=CCO.Cl[Al](Cl)Cl"`（**含 LA 后缀**）
3. `MolFromSmiles("C=CCO.Cl[Al](Cl)Cl")` → Mol 包含 AlCl₃ 原子
4. 重原子计数错误 → SMILES→XYZ 映射错误 → forming bonds 索引错位

**影响**: AlCl₃ surrogate 下，organic_smiles 仍含 LA，整个 S1 embedding 和 SMARTS 匹配索引空间混乱。

---

### 🔴 问题 #4：`write_smiles_to_xyz_map` additive 计数硬编码 LiCl

**位置**: `rph_core/steps/conformer_search/engine.py:2921-2938`

**代码**:
```python
additive_count = 0
if la_enabled and total_count >= 2:
    last_two = [xyz_symbols[-2], xyz_symbols[-1]]
    # ...
    has_li_in_smiles = "Li" in smoles_elements
    if last_two[0] == "Li" and not has_li_in_smiles:   # ← 硬编码 "Li"
        additive_count = total_count - xyz_symbols.index("Li") if "Li" in xyz_symbols else 0
        # ... 反向扫描逻辑 ...
```

**问题**:
1. **硬编码 `"Li"`**：仅检测 Li 基 additive。MgCl₂ 末尾两原子是 `[Cl, Cl]`，`last_two[0] == "Li"` 为 False → `additive_count = 0`
2. **硬编码 `[-2]`**：假设恰好 2 个 additive 原子。AlCl₃ 有 4 个 additive 原子
3. **`xyz_symbols.index("Li")`**：找第一个 Li——若有机分子本身含 Li（有机锂试剂），误判

**各 surrogate 的表现**:

| Surrogate | XYZ 末尾 | `last_two[0]` | additive_count | 映射文件正确？ |
|-----------|---------|---------------|----------------|--------------|
| LiCl | `[..., Li, Cl]` | `"Li"` | 2 | ✅ |
| MgCl₂ | `[..., Cl, Mg, Cl]` | `"Cl"` | **0** | ❌ LA 原子被标为 heavy |
| AlCl₃ | `[..., Cl, Al, Cl, Cl]` | `"Cl"` | **0** | ❌ LA 原子被标为 heavy |

**影响**: MgCl₂/AlCl₃ 下，`smiles_to_xyz_map.json` 将 LA 原子误标为 `"type": "heavy"` 并分配 `smiles_idx`，导致下游所有依赖该映射的索引解析错误。

---

### 🟠 问题 #5：无显式 `organic_atom_count` 字段

**位置**: 全局架构问题

**现象**: `LewisAcidAdditive` 数据类有 `additive_atom_indices` 和 `additive_elements`，但**没有 `organic_atom_count` 字段**。每个消费方独立重推有机/LA 边界：

| 消费方 | 推导方式 | 代码位置 |
|--------|---------|---------|
| orchestrator | `len(symbols) - len(additive_elements)` | `orchestrator.py:2600` |
| retro_scanner | `metal_idx = additive_atom_indices[0]` | `retro_scanner.py:294` |
| validator | 直接用 `additive_atom_indices` | `validator.py:176` |
| forming_bonds_resolver | `set(additive_indices)` 过滤 | `forming_bonds_resolver.py:145` |
| write_smiles_to_xyz_map | 反向扫描（问题 #4） | `engine.py:2929-2936` |

**影响**: 边界推导逻辑散落各处，任何一处推导错误都会导致索引偏移，且难以追踪。新增第二种 additive 类型（如溶剂分子）时会全面崩溃。

---

### 🟠 问题 #6：SMARTS 匹配索引空间隐式耦合

**位置**: `rph_core/steps/step2_retro/smarts_matcher.py:159-166, 205`

**代码**:
```python
# line 161: 从 canonical SMILES 构建 Mol（仅重原子，无 H，无 LA）
mol = Chem.MolFromSmiles(product_smiles)

# line 166: 将 XYZ 坐标嵌入 Mol
mol = self._embed_coords_into_mol(mol, coords)
# coords 来源: read_xyz(product_xyz) → 含 [重原子, H, LA] 全部坐标

# _embed_coords_into_mol 内部:
for i in range(mol.GetNumAtoms()):        # mol 仅有 N_heavy 个原子
    conf.SetAtomPosition(i, coords[i])    # coords[i] 取前 N_heavy 个
# 因 XYZ 中重原子在前 → coords[0..N_heavy-1] 恰好是重原子 → 正确

# line 205: identify(mol, coords) ← coords 仍是完整数组（重+H+LA）
topo_candidate = identify(mol, coords)
```

**问题**: `identify` 函数接收的 `coords` 是完整数组（含 H 和 LA 坐标），但 Mol 仅有重原子。如果 `identify` 内部基于 `len(coords)` 做索引或访问 `coords[N_heavy:]`，会取到 H 或 LA 的坐标。

**当前是否出错**: 经追踪，`_topological_core_identification` 等函数基于 Mol 的 atom index 索引 coords，而 Mol atom index ∈ [0, N_heavy)，恰好落在重原子区间 → **当前恰好正确**。

**风险**: 这是**巧合性正确**，依赖 XYZ 重原子在前的不变量。任何打破该不变量的改动（如 H 原子穿插、LA 原子不在末尾）都会静默出错。

---

### 🟡 问题 #7：additive_atom_indices 位置推导，不匹配仅 warning

**位置**: `rph_core/orchestrator.py:2596-2610`

**代码**:
```python
num_additive = len(la_additive.additive_elements)
if num_additive > 0 and len(_la_symbols) >= num_additive:
    last_n = list(_la_symbols[-num_additive:])
    expected = list(la_additive.additive_elements)
    if last_n == expected:
        la_additive.additive_atom_indices = tuple(...)
    else:
        self.logger.warning(
            f"[LA] Expected {expected} at end of XYZ, got {last_n}"
        )
        # ← additive_atom_indices 保持为 () → 下游全部失效但仅 warning
```

**问题**: 当 XYZ 末尾元素不匹配预期时，仅打印 warning，`additive_atom_indices` 保持为空。下游 S2/S3 的所有 LA 相关逻辑（coordination trace, mode projection）因空 indices 而静默跳过。

**影响**: 如果 LA 原子在 DFT 优化中位置变动（如被排序算法重排），indices 无法推导，整个 LA 质量门控静默失效，但 pipeline 继续运行并产出**不可信结果**。

---

### 🟡 问题 #8：charge/multiplicity 未传入 QC 输入

**位置**: `rph_core/orchestrator.py:245-246`（存储）vs QC 输入生成（未消费）

**现象**: `LewisAcidAdditive` 存储 `charge` 和 `multiplicity` 字段，从 `surrogate_registry` 读取。但 QC 输入生成（`qc_interface.py` 的 task 创建路径）**未读取这些字段**来修改系统电荷/自旋。

**影响**:
- LiCl（charge=0, mult=1）：恰好与默认值一致 → 无影响
- 未来若引入带电 surrogate（如裸 `Li+`）：QC 输入会用错误的 charge/multiplicity
- 当前为 **Phase 2+ 占位符**，但文档未明确标注

---

### 🟢 问题 #9：`LiCoordinationTrace` 字段名保留 `li_` 前缀

**位置**: `rph_core/lewis_acid/models.py:27-45`

**现象**: 数据类字段名如 `li_cl_distances`, `li_o_distances`, `li_coordination_sites` 保留 `li_` 前缀，但 docstring 和实际用途是通用金属。

**影响**: 代码可读性和可维护性降低；配置键已通用化（`metal_cl_max`/`metal_o_max`），但数据类未同步。

---

## 3. 各 Surrogate 可用性矩阵

基于以上确认的问题，各 surrogate 的实际可用性：

| Surrogate | SMILES | 原子数 | S1 embedding | SMILES 剥离 | XYZ 映射 | S2 SMARTS | 综合可用性 |
|-----------|--------|--------|-------------|------------|---------|-----------|-----------|
| **LiCl** | `[Li]Cl` | 2 | ✅ | ✅ | ✅ | ✅ | ✅ **可用** |
| **MgCl₂** | `Cl[Mg]Cl` | 3 | ❌ 问题#2 | ✅ | ❌ 问题#4 | ⚠️ | ❌ **不可用** |
| **AlCl₃** | `Cl[Al](Cl)Cl` | 4 | ❌ 问题#2 | ❌ 问题#3 | ❌ 问题#4 | ❌ | ❌ **不可用** |

> **结论**：虽然 `config/defaults.yaml:768` 声明 `surrogate` 可选 `LiCl | MgCl2 | AlCl3`，但当前实现**仅 LiCl 可正确运行**。

---

## 4. S1 校验门控失效分析

设计文档（§7.3, §11 Gate A）要求 S1 后校验 LA 几何完整性。实际执行链路：

```
S1 conformer search
  ├─ _step_rdkit_embed() → _init.xyz
  ├─ CREST/xTB stages → best.xyz
  │   └─ _validate_la_geometry(best_xyz)     ← 问题#1: indices=() → return True
  ├─ gfn0 ensemble
  │   └─ _validate_la_geometry(ensemble)     ← 问题#1: indices=() → return True
  └─ DFT optimization → product_min.xyz

orchestrator post-S1
  └─ 推导 additive_atom_indices               ← 此时才填充

S2 quality gate
  └─ _check_s2_quality_gate()                 ← indices 可用 → 生效
```

**结论**：S1 阶段的 LA 几何校验**完全空跑**（两个调用点都因空 indices 而提前返回 True）。仅 S2 quality gate 实际生效。

---

## 5. 修复建议（按优先级）

### P0 — 立即修复（阻断 MgCl₂/AlCl₃ 使用）

#### Fix #2a: 通用化 LA 哑键创建
```python
# engine.py:511-528 — 替换 la_atoms[0] 硬编码
# 找到 surrogate 中的金属原子（最高原子序数），而非第一个原子
metal_idx_in_la = max(range(len(la_atoms)), 
                      key=lambda i: la_mol.GetAtomWithIdx(i).GetAtomicNum())
rw_mol.AddBond(carbonyl_o_idx, la_atoms[metal_idx_in_la], Chem.BondType.SINGLE)
# 记录用于后续删除
self._la_metal_graft_idx = la_atoms[metal_idx_in_la]
```

#### Fix #2b: 通用化哑键删除
```python
# engine.py:556-564 — 替换 N-2 硬编码
if self.la_enabled and carbonyl_o_idx is not None and hasattr(self, '_la_metal_graft_idx'):
    rw_mol = Chem.RWMol(mol)
    rw_mol.RemoveBond(carbonyl_o_idx, self._la_metal_graft_idx)
    mol = rw_mol.GetMol()
    Chem.FastFindRings(mol)
```

#### Fix #3: 修复 `_strip_la_suffix` 正则
```python
def _strip_la_suffix(smiles: str) -> str:
    """剥离 LA 后缀——基于最后一个 '.' 分割，而非正则匹配原子 token。"""
    # LA 后缀格式: "{organic}.{la_smiles}"
    # 最安全的做法：从配置获取 la_smiles 后精确匹配
    idx = smiles.rfind('.')
    if idx > 0:
        suffix = smiles[idx:]
        # 校验后缀是合法的 LA SMILES（不含有机骨架特征）
        # 或直接从 la_additive.surrogate_smiles 精确比对
        return smiles[:idx]
    return smiles
```
> 更优方案：不在 `_strip_la_suffix` 中猜正则，而是直接传入 `la_additive.surrogate_smiles` 做精确字符串匹配剥离。

#### Fix #4: 通用化 additive 计数
```python
# engine.py:2921-2938 — 从 la_additive.additive_elements 精确推导
additive_count = 0
if la_enabled:
    expected_tail = list(getattr(la_additive, 'additive_elements', []))
    n = len(expected_tail)
    if n > 0 and len(xyz_symbols) >= n and xyz_symbols[-n:] == expected_tail:
        additive_count = n
```
> 需要将 `la_additive` 对象传入 `write_smiles_to_xyz_map`（当前仅传 `la_enabled: bool`）。

### P1 — 尽快修复（提升健壮性）

#### Fix #1: S1 校验提前推导 indices
```python
# engine.py:_validate_la_geometry — 不依赖预填充的 indices，就地推导
def _validate_la_geometry(self, xyz_path, la_additive):
    coords, atoms = read_xyz(xyz_path)
    # 就地推导：末尾 N 个原子 = additive
    n_add = len(la_additive.additive_elements)
    expected = list(la_additive.additive_elements)
    if atoms[-n_add:] != expected:
        return True  # 无法确定，跳过
    additive_indices = tuple(range(len(atoms) - n_add, len(atoms)))
    # ... 后续校验逻辑 ...
```

#### Fix #5: 添加 `organic_atom_count` 字段
```python
@dataclass
class LewisAcidAdditive:
    # ... existing fields ...
    organic_atom_count: Optional[int] = None  # S1 完成后填充
    # additive_atom_indices 的补集：[0, organic_atom_count) 为有机原子
```

#### Fix #7: 不匹配时 fail-fast
```python
# orchestrator.py:2607-2610 — 改 warning 为 raise
else:
    raise RuntimeError(
        f"[LA] XYZ 末尾元素 {last_n} 与预期 {expected} 不匹配，"
        f"无法确定 additive 原子索引。中止以避免下游索引错误。"
    )
```

### P2 — 后续优化

- Fix #6: SMARTS matcher 的 `identify` 函数应仅接收 `coords[:N_heavy]`，并添加注释说明索引空间
- Fix #8: 在 `qc_interface.py` task 创建路径中消费 `la_additive.charge/multiplicity`
- Fix #9: 将 `LiCoordinationTrace` 字段名通用化（`metal_cl_distances` 等），保留 `li_` 别名做向后兼容

---

## 6. 验证检查清单

修复后应逐项验证：

- [ ] **LiCl 单元测试**: 端到端运行（S1→S2→S3），确认 indices 全链路一致
- [ ] **MgCl₂ 单元测试**: 确认 Fix #2 后 O–Mg 哑键正确创建和删除
- [ ] **AlCl₃ 单元测试**: 确认 Fix #3 后 SMILES 正确剥离
- [ ] **S1 校验生效**: `_validate_la_geometry` 在 indices 为空时也能就地推导并执行
- [ ] **additive 计数**: `smiles_to_xyz_map.json` 对三种 surrogate 都正确标记 additive 原子
- [ ] **形成键过滤**: forming bonds 不含任何 additive 索引（三种 surrogate）
- [ ] **fail-fast**: XYZ 末尾不匹配时 pipeline 中止而非静默继续
- [ ] **断点续算**: checkpoint identity hash 包含完整 LA 配置

---

## 7. 附录：关键文件索引

| 文件 | 关键位置 | 角色 |
|------|---------|------|
| `rph_core/lewis_acid/models.py` | L6-67 | 3 个数据类定义 |
| `rph_core/orchestrator.py` | L193-255 | `_prepare_lewis_acid_system()` |
| `rph_core/orchestrator.py` | L2590-2612 | additive_atom_indices 推导 |
| `rph_core/steps/conformer_search/engine.py` | L460-592 | `_step_rdkit_embed()` LA 原子添加 |
| `rph_core/steps/conformer_search/engine.py` | L877-907 | `_validate_la_geometry()` |
| `rph_core/steps/conformer_search/engine.py` | L2860-2866 | `_strip_la_suffix()` |
| `rph_core/steps/conformer_search/engine.py` | L2869-3003 | `write_smiles_to_xyz_map()` |
| `rph_core/steps/step2_retro/smarts_matcher.py` | L131-254 | `find_reactive_bonds()` LA 模式 |
| `rph_core/steps/step2_retro/smarts_matcher.py` | L456-462 | `_embed_coords_into_mol()` |
| `rph_core/steps/step2_retro/retro_scanner.py` | L270-355 | `_compute_la_coordination_trace()` |
| `rph_core/steps/step3_opt/validator.py` | L109-208 | `validate_la_mode_projection()` |
| `rph_core/utils/forming_bonds_resolver.py` | L129-159 | `_filter_additive_atoms()` |
| `config/defaults.yaml` | L765-807 | `lewis_acid:` 配置段 |

---

## 8. 二次复核补充（2026-06-29）

本轮复核将报告结论与当前代码重新对齐，确认原报告的大方向成立：当前实现仍主要只对 LiCl 可靠，MgCl₂/AlCl₃ 会在原子顺序、金属定位、映射持久化和 checkpoint 复用上形成连锁风险。但有 1 处需要修正，另有 3 条问题链路需要补入修复范围。

### 8.1 需要修正的原报告判断

**问题 #3 的作用范围应收窄**：当前 `ConformerEngine._step_rdkit_embed()` 已在 `engine.py:473-480` 使用 `la_additive.surrogate_smiles` 做精确后缀匹配：

```python
additive_smi = getattr(self.la_additive, 'surrogate_smiles', '')
suffix = f".{additive_smi}"
if smiles.endswith(suffix):
    organic_smiles = smiles[:-len(suffix)]
```

因此 AlCl₃ 的 S1 embedding 不再走 `_strip_la_suffix()` 正则路径；原报告中“S1 embedding 无法剥离 AlCl₃”的判断对当前代码不成立。  
但 `_strip_la_suffix()` 仍被 `write_smiles_to_xyz_map()` 使用，且仍无法处理 `Cl[Al](Cl)Cl` 这类带分支 SMILES，所以问题 #3 应改为：**映射文件生成阶段的 LA 后缀剥离仍有 AlCl₃ 风险**。

### 8.2 新增问题 #10：金属原子被假定为 additive 的第一个原子

**严重度**：CRITICAL  
**位置**：
- `engine.py:885`：`metal_idx = additive_indices[0]`
- `retro_scanner.py:294`：`metal_idx = la_additive.additive_atom_indices[0]`
- `retro_scanner.py:309`：校验 `symbols[metal_idx] == additive_elements[0]`

`additive_elements` 当前来自 RDKit 原子顺序：

| Surrogate | SMILES | RDKit additive_elements | 真正金属位置 |
|-----------|--------|--------------------------|--------------|
| LiCl | `[Li]Cl` | `("Li", "Cl")` | 0 |
| MgCl₂ | `Cl[Mg]Cl` | `("Cl", "Mg", "Cl")` | 1 |
| AlCl₃ | `Cl[Al](Cl)Cl` | `("Cl", "Al", "Cl", "Cl")` | 1 |

这意味着即使 `additive_atom_indices` 成功推导，S1 几何校验和 S2 coordination trace 仍会把 Cl 当作 metal，导致：

- metal-Cl 距离变成 Cl-Mg/Cl-Al 或 Cl-Cl 距离，物理含义错误；
- metal-O 距离变成 Cl-O，coordination gate 误判；
- MgCl₂/AlCl₃ 的质量门控可能产生假阳性或假阴性。

**修复方向**：`LewisAcidAdditive` 必须显式保存 `metal_symbol`、`metal_idx_in_additive`、`metal_atom_index`、`halide_atom_indices`，所有消费方禁止再用 `additive_atom_indices[0]` 推断金属。

### 8.3 新增问题 #11：S2 SMARTS fallback 的 LA 参数没有贯穿

**严重度**：HIGH  
**位置**：`orchestrator.py:716-722`, `orchestrator.py:884-889`

`_finalize_forming_bonds_for_s2()` 和 `_resolve_forming_bonds_for_s2()` 的部分 SMARTS fallback 调用没有传入 `la_additive` 和 canonical `product_smiles`：

```python
smarts_result = matcher.find_reactive_bonds(
    product_xyz=Path(product_xyz_file),
    cleaner_data=self._build_smarts_fallback_context(cleaner_data),
)
```

这会让 `SMARTSMatcher` 走非 LA 模式，用完整 XYZ 推断 Mol 连接性。若 LA 原子与底物距离较近，几何建图有机会把 LA 边引入拓扑判断，形成“主路径 LA-aware，fallback LA-unaware”的行为分叉。

**修复方向**：所有 SMARTS fallback 调用必须传入：

- `la_additive=la_additive`
- `product_smiles=<canonical product smiles>`
- `coords_for_matching=coords[:organic_heavy_atom_count]` 或由 atom mapping 显式切片

### 8.4 新增问题 #12：checkpoint / provenance 未纳入 atom mapping 身份

**严重度**：HIGH  
**位置**：
- `checkpoint_manager.py:393-396`：LA signature 只含 enabled/default_surrogate，未含完整 registry entry、additive_elements、metal position
- `checkpoint_manager.py:428-447`：S2 signature 不含 `atom_mapping_hash`
- `orchestrator.py:2723-2740`, `orchestrator.py:2890-2917`：S2/S3 checkpoint metadata 未持久化 atom mapping 文件路径和 hash

风险是：修复映射逻辑后，旧的 S2/S3 checkpoint 仍可能被认为有效并复用，从而继续使用历史错误的 `forming_bonds`、scan profile 或 LA quality flags。

**修复方向**：S1/S2/S3 的 checkpoint metadata 和 provenance 必须包含：

- `atom_mapping_json`
- `atom_mapping_hash`
- `lewis_acid_mapping_hash`
- `mapping_schema_version`

并将这些字段纳入 S2/S3 signature。只要 LA 配置、surrogate 原子顺序、metal_idx 或 S1 atom order 改变，就必须触发 S2/S3 重算。

---

## 9. 完整修复方案：将 atom mapping 提升为一等 artifact

用户目标是“原子映射链路每一步都清晰可知，并有专门文件储存相应原子映射关系”。建议不要继续在各模块内重复猜测 `[organic heavy][H][LA]` 边界，而是新增统一的 mapping 层。

### 9.1 新增统一模块

新增文件：`rph_core/utils/atom_mapping.py`

职责：

1. 从 canonical product SMILES、geometry SMILES、S1 XYZ 和 `LewisAcidAdditive` 构建标准映射。
2. 校验 XYZ 尾部 additive 元素、金属原子、卤素原子、organic/heavy/H/additive 分区。
3. 写入和读取每一步的 `atom_mapping.json`。
4. 提供唯一 API 给 S1/S2/S3 使用，替代散落的边界推导。

建议核心 API：

```python
def build_s1_atom_mapping(
    canonical_smiles: str,
    geometry_smiles: str,
    product_xyz: Path,
    la_additive: Optional[LewisAcidAdditive],
) -> AtomMapping

def write_atom_mapping(mapping: AtomMapping, path: Path) -> Path

def load_atom_mapping(path: Path) -> AtomMapping

def resolve_additive_indices(symbols: Sequence[str], additive: LewisAcidAdditive) -> LewisAcidAdditive

def organic_heavy_coords(coords: np.ndarray, mapping: AtomMapping) -> np.ndarray
```

### 9.2 扩展 `LewisAcidAdditive`

在 `rph_core/lewis_acid/models.py` 中添加显式字段：

```python
organic_atom_count: Optional[int] = None
organic_heavy_atom_count: Optional[int] = None
metal_symbol: Optional[str] = None
metal_idx_in_additive: Optional[int] = None
metal_atom_index: Optional[int] = None
halide_atom_indices: Tuple[int, ...] = ()
atom_mapping_json: Optional[str] = None
atom_mapping_hash: Optional[str] = None
```

构建 additive 时，从 surrogate RDKit Mol 推导 `metal_idx_in_additive`，规则为“候选金属元素优先，其次最高原子序数”。候选金属至少包含 `Li, Na, K, Mg, Ca, Zn, Al, B`；MgCl₂/AlCl₃ 应得到 `metal_idx_in_additive=1`。

### 9.3 标准 `atom_mapping.json` schema

每一步写独立文件，但 schema 保持一致。建议版本 `2.0`：

```json
{
  "schema_version": "2.0",
  "stage": "S1",
  "index_base": 0,
  "canonical_product_smiles": "...",
  "geometry_product_smiles": "...",
  "xyz_path": "S1_ConfGeneration/product/product_min.xyz",
  "xyz_sha256": "...",
  "atom_order_contract": "[organic_heavy][organic_hydrogen][lewis_acid_additive]",
  "counts": {
    "organic_heavy": 12,
    "organic_hydrogen": 14,
    "organic_total": 26,
    "additive": 3,
    "total": 29
  },
  "ranges": {
    "organic_heavy": [0, 12],
    "organic_hydrogen": [12, 26],
    "additive": [26, 29]
  },
  "lewis_acid": {
    "enabled": true,
    "surrogate_name": "MgCl2",
    "surrogate_smiles": "Cl[Mg]Cl",
    "additive_elements": ["Cl", "Mg", "Cl"],
    "additive_atom_indices": [26, 27, 28],
    "metal_symbol": "Mg",
    "metal_atom_index": 27,
    "halide_atom_indices": [26, 28],
    "charge": 0,
    "multiplicity": 1
  },
  "atoms": [
    {"xyz_idx": 0, "smiles_idx": 0, "element": "C", "role": "organic_heavy"},
    {"xyz_idx": 12, "smiles_idx": null, "element": "H", "role": "organic_hydrogen"},
    {"xyz_idx": 26, "smiles_idx": null, "element": "Cl", "role": "lewis_acid_halide"},
    {"xyz_idx": 27, "smiles_idx": null, "element": "Mg", "role": "lewis_acid_metal"}
  ],
  "forming_bonds": {
    "index_space": "xyz_0based",
    "bonds": [[1, 5], [3, 9]],
    "source": "S2_Retro"
  },
  "parent_mapping": null,
  "mapping_sha256": "..."
}
```

### 9.4 每一步应写出的专门文件

| 阶段 | 文件 | 内容 |
|------|------|------|
| S0 | `S0_Classification/atom_mapping.json` | canonical SMILES atom order、S0 forming_bonds 的 SMILES/organic-heavy index 空间、反应中心注释 |
| S1 | `S1_ConfGeneration/product/atom_mapping.json` | SMILES→XYZ 映射、organic/H/LA ranges、metal/halide indices、product_min.xyz hash |
| S2 | `S2_Retro/atom_mapping.json` | 继承 S1 mapping，增加 `forming_bonds`、scan input/ts_guess/intermediate atom count/hash 校验 |
| S3 | `S3_TransitionAnalysis/atom_mapping.json` | 继承 S2 mapping，增加 ts_final/intermediate_opt/reactant_sp atom order 校验、LA mode projection 输入索引 |
| root | `atom_mapping_manifest.json` | 汇总各阶段 mapping 文件路径、hash、schema version，供外部 S4/RPH_Postprocess 消费 |
| checkpoint | `pipeline.state` metadata | 每一步保存 mapping path/hash；S2/S3 signature 纳入上游 mapping hash |

旧文件 `smiles_to_xyz_map.json` 可保留为兼容别名，但应由新 `atom_mapping.json` 生成，不能再作为权威来源。

### 9.5 按优先级实施

**P0：阻断错误结果**

1. 为 `LewisAcidAdditive` 增加 metal/halide 显式字段。
2. 修复 `_step_rdkit_embed()`：dummy bond 连接到 `metal_idx_in_additive` 对应的新增原子，并记录本地 `graft_atom_idx`；删除时使用同一个 idx，失败时 warning，不再静默 `pass`。
3. 修复 `_validate_la_geometry()` 和 `_compute_la_coordination_trace()`：使用 `metal_atom_index` / `halide_atom_indices`，禁止 `additive_atom_indices[0]`。
4. 修复 `write_smiles_to_xyz_map()`：改为接收 `la_additive` 或迁移到 `atom_mapping.py`；用 `additive_elements` 精确尾部匹配；AlCl₃ 后缀用 surrogate 精确匹配，不再用正则猜。
5. orchestrator 推导 additive indices 失败时 fail-fast；LA enabled 但 indices 为空时禁止进入 S2/S3。

**P1：建立清晰映射链路**

1. 写出 S1/S2/S3 `atom_mapping.json` 和 root `atom_mapping_manifest.json`。
2. S2 forming_bonds resolution 全部读取 S1 mapping，并在 metadata 中声明 `index_space="xyz_0based"`。
3. 所有 SMARTS fallback 传入 canonical SMILES、`la_additive` 和 mapping；LA 模式只使用 organic heavy 坐标。
4. S3 mode projection 从 S3 mapping 读取 additive/forming bond indices，并先校验 TS displacement 原子数等于 mapping total。
5. `mechanism_meta.json` 增加 `atom_mapping_hash` 和 `forming_bonds.index_space`。

**P2：checkpoint 与兼容**

1. S2 signature 加入 `atom_mapping_hash`、`lewis_acid_mapping_hash`、`surrogate_smiles`、`metal_idx_in_additive`。
2. S3 signature 加入 upstream S2 mapping hash 和 ts/intermediate/product input hashes。
3. pipeline.state 每步 metadata 保存 mapping path/hash；resume 时 hash mismatch 必须重算。
4. 保留 `LiCoordinationTrace.li_*` 字段别名，同时新增 `metal_*` 字段，后续版本再去 Li 命名。
5. QC charge/multiplicity 保持从 config 单一入口流入；如果未来 surrogate 非中性，再将 LA charge/multiplicity 合并进系统 charge policy。

### 9.6 测试方案

新增或扩展测试：

- `tests/test_lewis_acid_atom_mapping.py`
  - LiCl/MgCl₂/AlCl₃ 三种 surrogate 的 `metal_idx_in_additive`、`metal_atom_index`、`halide_atom_indices`。
  - S1 mapping 的 ranges/counts/roles 正确。
  - AlCl₃ `Cl[Al](Cl)Cl` 不污染 canonical organic SMILES。
- `tests/test_lewis_acid_embedding.py`
  - dummy bond 对 Mg/Al 创建并删除；删除失败会记录 warning 或抛出受控异常。
- `tests/test_lewis_acid_s2_mapping.py`
  - SMARTS fallback 在 LA 模式下只消费 organic heavy 坐标。
  - forming_bonds 不包含 additive indices，且 index_space 写入 `xyz_0based`。
- `tests/test_lewis_acid_checkpoint_mapping.py`
  - 改变 surrogate 或 mapping hash 后 S2/S3 checkpoint 不可复用。

推荐验证命令：

```bash
pytest -v tests/test_lewis_acid_atom_mapping.py tests/test_lewis_acid_s2_mapping.py
python scripts/ci/check_imports.py rph_core
```

---

*报告生成方式：2 个 explore agent 并行扫描 + 关键代码段逐行验证 + 配置/SMILES 交叉比对*
