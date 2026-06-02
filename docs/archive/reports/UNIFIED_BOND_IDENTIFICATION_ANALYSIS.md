# 统一分析报告：键识别问题全链路诊断

> 注：文中涉及的 `forward_scan` 为历史别名；该别名已在 v2.1.1 中删除。

## 0. 金标准

以 **产物 `product_min.xyz`** 的 XYZ 1-based 原子编号为唯一权威编号空间：

| 金标准 forming bond | product 距离 | TS 距离 |
|---|---:|---:|
| **12–19** | 1.571 Å | 1.571 Å |
| **15–16** | 1.545 Å | 1.545 Å |

---

## 1. 化学反应全景图

### 1.1 反应类型

分子内 [4+3] 环加成（intramolecular allenamide–furan [4+3] cycloaddition）

```
前体（allenamide）：C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C
      ↓  [4+3] 环加成
产物（环化产物）：CC(C)(C)OC(=O)N1CCC[C@@]23C=C[C@@H](CC(=O)[C@H]12)O3
```

### 1.2 反应中心原子对照表

mapped SMILES 提供了前体 → 产物中各原子的对应关系（`rxn_smiles_mapped`）：

| Map# | 前体角色 | 产物角色 | 参与反应？ |
|---|---|---|---|
| 1–4 | Boc 甲基 C | Boc 甲基 C（不变） | 否 |
| 5 | Boc 酯 O | Boc 酯 O（不变） | 否 |
| 6 | 羰基 C | 羰基 C（不变） | 否 |
| 7 | 羰基 O | 羰基 O（不变） | 否 |
| 8 | N | N（从 allenamide N 变成环内 N） | 断键 |
| **9** | **CH₂** | **CH₂（不变）** | 桥连链 |
| **10** | **CH₂** | **CH₂（不变）** | 桥连链 |
| **11** | **CH₂** | **CH₂（不变）** | 桥连链 |
| **12** | **furan C（α位）** | **桥头 C（quaternary）** | ✅ 成键中心 |
| **13** | **allene =CH** | **C=（双键）** | ✅ 成键中心 |
| 14 | allene =CH₂ | CH（双键旁） | 是 |
| **15** | **allene =C=（中心）** | **CH（桥头）** | ✅ 成键中心 |
| **16** | **furan CH** | **CH₂** | ✅ 成键中心 |
| 17 | furan CH | C=O | 断键/重排 |
| 18 | furan O | C=O 中的 O | 断键/重排 |
| **19** | **furan CH（α位）** | **桥头 C** | ✅ 成键中心 |
| 20 | （不存在） | 环内 O | 新原子 |

**核心观察**：Map 12 和 Map 19 在前体中是 furan 环的 α-碳，通过 furan 环上的 O–C 键直接相邻（AROMATIC 键）。在产物中，12 和 19 是桥头碳，通过新的 SINGLE 键连接。

### 1.3 金标准成键的化学意义

| 金标准 bond | 化学意义 |
|---|---|
| **12–19** | furan α-C → furan α-C：在前体中通过芳香键相连（furan 环的一部分），反应后变成桥头 C–C 单键——**这是环化反应中"旧键保留但角色改变"的键** |
| **15–16** | allene 中心 C → furan CH：在前体中不相邻，反应后形成新的 C–C 单键——**这是典型的"新形成键"** |

---

## 2. 三套编号系统及其映射关系

### 2.1 前体 XYZ 编号（precursor_min.xyz，40 个原子，19 个重原子）

前体由 SMILES `C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C` 构建，RDKit canonicalization 产生的原子顺序与 mapped SMILES 不同：

```
XYZ#  元素  结构图中的位置          对应 Map#
  1    C    allene 终端 =CH₂       map 13
  2    C    allene 中心 =C=         map 15
  3    C    allene 终端 =CH₂       map 14
  4    N    氮（连接 allene+链+Boc）map 8
  5    C    链 CH₂                  map 9
  6    C    链 CH₂                  map 10
  7    C    链 CH₂                  map 11
  8    C    furan α-C               map 12
  9    C    furan CH                map 19
 10    C    furan CH                map 16
 11    C    furan CH                map 17
 12    O    furan O                 map 18
 13    C    羰基 C                   map 6
 14    O    羰基 O                   map 7
 15    O    酯桥 O                   map 5
 16    C    tBu 中心 C              map 2
 17    C    tBu CH₃                 map 1
 18    C    tBu CH₃                 map 3
 19    C    tBu CH₃                 map 4
```

**关键发现**：前体中 map 12（furan α-C）= XYZ **8**，map 19（furan α-C）= XYZ **9**。

因此：
- **前体 forming bond 1**：map 13–12（allene→furan）= 前体 XYZ **1–8**（距离 5.012 Å，未成键）
- **前体 forming bond 2**：map 15–16（allene 中心→furan）= 前体 XYZ **2–10**（距离 ?）

但用户报告前体构象搜索给出的键对是 **[3,8]** 和 **[1,11]**。对应 map：
- 前体 XYZ 3 = map 14，XYZ 8 = map 12 → **map 14–12**
- 前体 XYZ 1 = map 13，XYZ 11 = map 17 → **map 13–17**

### 2.2 产物 XYZ 编号（product_min.xyz，41 个原子，20 个重原子）

```
XYZ#  元素  对应 Map#
  1    C    map 1  (Boc CH₃)
  2    C    map 2  (Boc 中心 C)
  3    C    map 3  (Boc CH₃)
  4    C    map 4  (Boc CH₃)
  5    O    map 5  (酯桥 O)
  6    C    map 6  (羰基 C)
  7    O    map 7  (羰基 O)
  8    N    map 8  (氮)
  9    C    map 9  (链 CH₂)
 10    C    map 10 (链 CH₂)
 11    C    map 11 (链 CH₂)
 12    C    map 12 (桥头 C)       ← 成键中心
 13    C    map 13 (C= 双键)      ← 成键中心
 14    C    map 14 (CH)
 15    C    map 15 (桥头 CH)      ← 成键中心
 16    C    map 16 (CH₂)          ← 成键中心
 17    C    map 17 (C=O)
 18    O    map 18 (C=O 的 O)
 19    C    map 19 (桥头 C)       ← 成键中心
 20    O    map 20 (环内 O)
```

**关键发现**：产物 XYZ 编号 = map 编号 = 1-based。本反应中 **产物 XYZ 1-based = map number = map number**（三者合一）。

### 2.3 中间体编号（intermediate.xyz）

中间体是 S2.1 retro_scan 产出的拉伸构型。其原子顺序与产物 XYZ 相同（41 个原子），因为它是从产物拉伸出来的。

中间体中的键长：

| 键对（XYZ 1-based） | 中间体距离 | 产物距离 | 状态 |
|---|---:|---:|---|
| 12–19（金标准） | 3.314 Å | 1.571 Å | 已拉伸（接近断裂） |
| 15–16（金标准） | 1.525 Å | 1.545 Å | 仍接近产物值 |
| 12–13（错误键） | 1.290 Å | 1.521 Å | 反而更短了 |
| 13–14 | 3.207 Å | 1.335 Å | 已拉伸断裂 |
| 16–17 | 3.421 Å | 1.524 Å | 已拉伸断裂 |

**关键观察**：中间体中 12–13 距离（1.290 Å）比产物中还短——这说明 retro_scan 把 12–13 当作 forming bond 去拉伸时，实际上在压缩这对原子而非拉伸。

---

## 3. 为什么前体编号 [3,8]/[1,11] 与产物编号 [12,19]/[15,16] 不吻合

### 根本原因：前体和产物是不同的分子，有不同数目的原子

| | 前体 | 产物 |
|---|---|---|
| 重原子数 | 19 | 20 |
| 原子顺序 | RDKit 从前体 SMILES 构建 | RDKit 从产物 SMILES 构建 |
| 编号规则 | 各自独立的 1-based | 各自独立的 1-based |

**同一对化学键在不同分子中拥有完全不同的 XYZ 索引。** 例如：

| 化学键 | 前体 XYZ (1-based) | 产物 XYZ (1-based) | map 编号 |
|---|---|---|---|
| map 12 ↔ map 19 | 8–9（furan 环内键） | 12–19（桥头键） | 12–19 |
| map 15 ↔ map 16 | 2–10（allene–furan 跨距） | 15–16（环内键） | 15–16 |

这不是"编号混乱"，而是**前体和产物各自独立编号的必然结果**。问题出在 pipeline 没有在两者之间建立正确的原子映射。

---

## 4. Cleaner 为什么输出 12–13 / 15–16（而非 12–19 / 15–16）

### 4.1 Cleaner 的判断逻辑

Cleaner 的 `CoreExtractor._infer_bond_changes()` 对 mapped reaction graph 做纯拓扑差分：

```
对每对核心原子 (i, j):
    product_has_bond = (mapped_product_graph 中 i-j 有键?)
    reactant_has_bond = (mapped_reactant_graph 中 i-j 有键?)
    
    if product_has_bond and NOT reactant_has_bond → formed
    if reactant_has_bond and NOT product_has_bond → broken
```

### 4.2 各原子对在 mapped graph 中的状态

| 原子对 | 前体 graph | 产物 graph | Cleaner 判定 | 实际化学意义 |
|---|---|---|---|---|
| **12–19** | **有键**（furan 芳香键） | **有键**（桥头 SINGLE） | preserved（≠ formed） | ❌ 环化反应中角色转变的键 |
| **12–13** | **无键** | **有键**（SINGLE） | **formed** ✅ | 新形成的 C–C 键（但非金标准定义的 forming bond） |
| **15–16** | **无键** | **有键**（SINGLE） | **formed** ✅ | 正确的新形成键 |
| 8–13 | 有键 | 无键 | broken | 正确 |
| 12–18 | 有键 | 无键 | broken | 正确 |
| 13–15 | 有键 | 无键 | broken | 正确 |
| 16–19 | 有键 | 无键 | broken | 正确 |

### 4.3 问题的本质

Cleaner 的"有键/无键"二分法无法处理以下化学语义：

> **12–19 在前体中是 furan 芳香键，在产物中是桥头 C–C 单键。从环化反应机理角度看，这应当被归类为"反应中心键的角色转变"而非简单的"保留"。**

但对于 Cleaner 的纯图差分来说，12–19 两端都有键，所以它**永远不可能**被判为 formed。与此同时，12–13 确实在产物中才出现，所以它**必然**被判为 formed。

这不是 Cleaner 的 bug——这是其 bond-change 语义定义的**系统性盲区**。

### 4.4 更深层：RXNMapper 的映射错误

Cleaner 接受的 mapped reaction 本身来自 RXNMapper（confidence=0.225，LOW_CONFIDENCE），该映射将前体的 furan 环原子对应关系编码为：
- 前体 `c:12` → 产物 `C@@:12`（furan C → 桥头 C）
- 前体 `cH:19` → 产物 `C@H:19`（furan CH → 桥头 C）

在这个映射中，12–19 的 furan 环键被保留为"两端都有键"，从而使其不可能被识别为 forming bond。

---

## 5. S2 三种算法的拉伸机制详解

### 5.1 算法概览

| 算法 | 代码入口 | 起始几何 | 核心机制 | forming_bonds 的角色 |
|---|---|---|---|---|
| **retro_scan** | `retro_scanner.py:242` | 产物 `product_min.xyz` | XTB `$scan`：逐点约束优化 | **直接决定**约束哪对原子 |
| **forward_scan** | `retro_scanner.py:551` | 产物 `product_min.xyz` | 与 retro_scan **完全相同**（别名） | **直接决定**约束哪对原子 |
| **xtb_path_search** | `retro_scanner.py:655` | 中间体 `intermediate.xyz` → 产物 | XTB `--path`：链态路径优化 | **不参与约束**，仅用于事后距离计算 |

### 5.2 retro_scan / forward_scan 的详细拉伸过程

#### 步骤 1：准备

从 `product_min.xyz` 开始。forming_bonds（例如 `[12, 13]` 和 `[15, 16]`，当前系统传播的错误值）被传入。

#### 步骤 2：约束生成（`scan_policies.py`）

```python
# scan_policies.py:64–69
for bond_idx in config.constrained_bonds:
    bond = bonds[bond_idx]   # bonds = forming_bonds = [[12,13], [15,16]]
    constraints[f"{bond[0]} {bond[1]}"] = target_distance  # target = 1.8 Å
```

这生成 XTB 输入文件的约束块：
```
$constrain
  force constant=0.5
  distance: 12, 13, 1.800    ← 把原子 12-13 约束在 1.8 Å
  distance: 15, 16, 1.800    ← 把原子 15-16 约束在 1.8 Å
$scan
  mode=concerted
  1: 1.800, 3.500, 20        ← 20 步从 1.8 Å 扫描到 3.5 Å
```

#### 步骤 3：扫描执行（`_execute_scan`，direction="outward"）

```
retro_scanner.py:185–192:
    direction="outward":
        start_dist = scan_end_distance = 1.8 Å   ← 从产物的短键距离开始
        end_dist   = scan_start_distance = 3.5 Å  ← 向长键距离拉伸
```

物理意义：**从产物（键已形成）出发，逐步拉伸键直到断裂，模拟逆反应路径**。

#### 步骤 4：每步做什么

在每一步 k（k=0..19）：
- 当前约束距离 = 1.8 + k × (3.5−1.8)/20
- XTB 在该约束下做 constrained optimization
- 所有其他原子自由弛豫
- 记录能量 → 构成 scan profile

#### 步骤 5：TS guess 提取

从 scan profile 中找 TS guess：
- 找能量最大值（peak）或 knee point（拐点）
- 该点的几何即为 TS 初猜

### 5.3 xtb_path_search 的详细路径搜索过程

#### 步骤 1：准备

两个端点几何：
- **start_xyz** = `intermediate.xyz`（S2.1 retro_scan 产出的拉伸构型）
- **end_xyz** = `product_min.xyz`（S1 的产物构型）

#### 步骤 2：路径生成

调用 XTB 的 `--path` 算法：
```
xtb intermediate.xyz --path product.xyz --input path.inp
```

XTB `--path` 的核心参数：
```
$path
  nrun=1        ← 优化轮数
  npoint=25     ← 路径上的点数
  anopt=10      ← 每点优化步数
  kpush=0.003   ← 推力（排斥势）
  kpull=-0.015  ← 拉力（吸引势）
$end
```

**关键区别**：`--path` 不使用 `$constrain` 或 `$scan`。它使用 penalty function 在两个端点之间寻找**最小能量路径（MEP）**。

#### 步骤 3：forming_bonds 不参与路径计算

在 `run_path_search()` 中，forming_bonds 参数被传入但**仅用于事后分析**：

```python
# retro_scanner.py:769–772 — 仅为计算 TS 点处的键长
if returned_forming_bonds:
    dist = GeometryUtils.calculate_distance(coords, atom_i, atom_j)
```

**路径搜索的 TS guess 完全由几何端点决定，不受 forming_bonds 指定错误的影响。**

#### 步骤 4：TS guess 提取

从路径中取能量最高点作为 TS guess（`xtbpath_ts.xyz`）。

### 5.4 三种算法对错误 forming_bonds 的敏感度

| 算法 | forming_bonds 错误时的后果 | 敏感度 |
|---|---|---|
| **retro_scan** | 拉伸错误的原子对 → 能量曲线偏离真实反应坐标 → TS guess 质量差 → 可能 DEGRADED/topology drift | **高** |
| **forward_scan** | 与 retro_scan 相同（它是 retro_scan 的别名） | **高** |
| **xtb_path_search** | forming_bonds 不影响路径计算 → TS guess 由端点几何驱动 → 对错误 forming_bonds **几乎免疫** | **低** |

---

## 6. 为什么正常 run 即使 forming_bonds 错误也能找到正确 TS

### 6.1 正常 run 的 S2 执行路径

```
S2.1 retro_scan(product_min.xyz, forming_bonds=[[12,13],[15,16]])
  → 产出 intermediate.xyz（拉伸构型，可能不理想）
  → scan_profile.json: status=可能 DEGRADED

S2.2 xtb_path_search(start=intermediate.xyz, end=product_min.xyz)
  → forming_bonds 不影响路径计算
  → 路径由 intermediate → product 的几何端点驱动
  → 产出更好的 TS guess
  → 覆盖 scan_profile.json
```

### 6.2 核心机制：S2.2 的鲁棒性

xtb_path_search 不依赖 forming_bonds 做约束——它只看两个端点几何。只要 intermediate.xyz 和 product_min.xyz 在物理上合理（中间体接近解离态、产物接近环化产物），路径搜索就能找到合理的 TS guess。

**这就是为什么正常 run（使用旧版 defaults.yaml，`path_search.enabled=true`）即使 forming_bonds 是 [12,13]/[15,16]（错误），仍能得到 activation_energy = +26.61 kcal/mol（合理正势垒）。**

### 6.3 baseline 为什么失败

```
S2.1 retro_scan(product_min.xyz, forming_bonds=[[12,13],[15,16]])
  → 拉伸 12-13 和 15-16
  → 但 12-13 在产物中是已形成的 C-C 单键（1.521 Å），拉伸它会破坏已有环结构
  → 能量单调上升，9/20 帧 topology drift
  → scan_quality: DEGRADED

S2.2 未执行（path_search.enabled=false）
  → TS guess 来自 S2.1 的 knee point（质量极差）

S3 TS optimization
  → 从差质量的 TS guess 出发优化
  → 可能"挽救"但结果不可靠
  → activation_energy = -44.53 kcal/mol（负势垒，说明 TS/reactant 能量关系反转）
```

---

## 7. 全链路根因总结

### 7.1 问题链条

```
① RXNMapper 低置信度映射（confidence=0.225）
   ↓ 把 12-19 编码为"保留键"而非"forming bond"
② Cleaner 纯图差分逻辑
   ↓ 输出 core_bond_changes = "12-13:formed;15-16:formed"
③ Repo loader / adapter 正规化
   ↓ forming_bonds = [[11,12],[14,15]]（0-based）但标注 index_base=0
④ S0 机制分类器直接采用
   ↓ mechanism_summary.json: forming_bonds = [[12,13],[15,16]]
⑤ Orchestrator 以 index_base=0 消费
   ↓ S2 把 [12,13] 当作 0-based（实际对应 XYZ 1-based 的 13-14）
⑥ S2 retro_scan 拉伸错误原子对
   ↓ 能量单调上升 + topology drift + DEGRADED
⑦ baseline 未启用 S2.2 path_search
   ↓ 差质量 TS guess 进入 S3
⑧ S3 "挽救"但产出负势垒 TS
   ↓ baseline validation 过弱（bond_lengths=[] 仍 pass）
⑨ 最终结果：负势垒被当作合法 baseline 缓存
```

### 7.2 两类独立问题

| 问题 | 描述 | 首次出现点 | 影响范围 |
|---|---|---|---|
| **A. 化学键身份错误** | 12–19 被替换为 12–13 | Cleaner `core_bond_changes`（根源：RXNMapper 低置信度映射） | 所有下游步骤 |
| **B. 编号空间语义错误** | 同一组数字在不同模块被当作不同编号空间 | S0→S2 交接（S0 写 map-space 数字但标注 index_base=0） | scan_profile.json, mechanism_meta.json |
| **C. S2 算法选择错误** | forward_scan 是 retro_scan 的别名；baseline 未启用 path_search | `retro_scanner.py:551` + `defaults.yaml` v2.1.1 | baseline 特有 |

### 7.3 为什么正常 run 能工作但 baseline 不能

| 因素 | 正常 run | baseline |
|---|---|---|
| forming_bonds 错误 | 有（12-13, 15-16） | 有（12-13, 15-16） |
| S2.2 path_search | ✅ 已执行（旧版默认 enabled=true） | ❌ 未执行（新版默认 enabled=false） |
| TS guess 来源 | xtb_path_search（对错误 forming_bonds 免疫） | retro_scan knee point（受错误 forming_bonds 严重影响） |
| 势垒 | +26.61 kcal/mol（合理） | -44.53 kcal/mol（不合理） |
| TS validation | 不适用（无 manifest 机制） | bond_lengths=[] 仍 pass |

---

## 8. 修复方向建议

### P0（最高优先级）

1. **在 S2 执行前增加几何预检**：检查 forming_bonds 指定的原子对在 product_min.xyz 中的实际距离。如果距离 < 1.8 Å（说明键已形成），应发出警告或触发 fallback。

2. **默认开启 `step2.path_search.enabled`**：证据表明 xtb_path_search 对 forming_bonds 错误有显著鲁棒性优势。当前的 `enabled: false` 使得所有新 run 都暴露在 retro_scan 的脆弱性下。

### P1

3. **加固 baseline TS validation**：`bond_lengths=[]` 不应 pass；`activation_energy < 0` 应标记为异常。

4. **修正 cleaner 的 bond-change 语义**：对环化反应中"aromatic→single"的键对，应引入"角色转变"分类而非简单二分法。

5. **统一编号空间标注**：所有 JSON 产物强制写出 `forming_bonds_xyz1` 和 `forming_bonds_xyz0`。

### P2

6. **实现真正的 forward_scan**：当前 `run_forward_scan()` 是 `run_retro_scan()` 的别名。

7. **收紧 RXNMapper 置信度阈值**：当前 confidence=0.225 仍被接受。

---

*报告生成时间：2026-04-23*  
*基于版本：RPH v2.1.1*  
*分析样本：benchmark_dft_theory/baselines/rx1/bl_3f31ea3a42fe5184 及 rph_output/rx_1*
