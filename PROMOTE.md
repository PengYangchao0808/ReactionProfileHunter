# 🎯 Reaction Profile Hunter v2.1

## 架构重设计文档 (Serial Architecture)

**文档日期**: 2026-01-09  
**核心理念**: "以终为始" (End-to-Begin) - 从产物逆推过渡态  
**目标反应**: [5+2] 1,3-偶极环加成  
**架构模式**: 瀑布式串行依赖 (Waterfall Sequential)

---

## 📋 目录

1. [设计理念](#1-设计理念)
2. [串行四步架构](#2-串行四步架构)
3. [详细实施细节](#3-详细实施细节)
4. [现有代码映射](#4-现有代码映射)
5. [需要新增的功能](#5-需要新增的功能)
6. [测试策略](#6-测试策略)
7. [开发路线图](#7-开发路线图)
8. [设计复审：效率与合理性优化](#七设计复审效率与合理性优化2026-01-10)
9. [检查测试流程设计](#八检查测试流程设计2026-01-10)
10. [体系化补全迁移方案](#十体系化补全迁移方案)

---

## 1. 设计理念

### 1.1 v2.0 并行架构的致命缺陷

| 问题 | 化学解释 | 后果 |
|:---:|---|---|
| **盲目 TS 搜索** | 不知道产物精确构象时寻找 TS | 可能找到错误的非对映异构体 |
| **Hammond 假说违背** | [5+2] 是吸热反应，TS 是"晚期过渡态" | TS 几何更接近产物，需要产物信息 |
| **原子映射混乱** | 底物+试剂拼凑导致原子索引不一致 | 畸变能计算错误 |

### 1.2 v2.1 核心变更

```
┌─────────────────────────────────────────────────────────────────┐
│                    v2.0 → v2.1 关键变更                        │
├─────────────────────────────────────────────────────────────────┤
│ ❌ 取消并行分叉                                                │
│    TS 生成模块必须等待产物构象搜索完成                         │
│                                                                 │
│ ✅ 产物导向策略 (Product-First)                                │
│    TS 初猜由产物通过"逆向拉伸 (Retro-Stretching)"生成         │
│                                                                 │
│ ✅ 双向锚定                                                    │
│    确定产物和底物构象后，才能构建 NEB/QST3 路径                │
│                                                                 │
│ ✅ 原子映射零成本                                              │
│    底物、TS、产物共享相同原子索引                              │
└─────────────────────────────────────────────────────────────────┘
```

### 1.3 化学优势

1. **立体选择性保真**: [5+2] 反应产生特定手性中心。从产物逆推，保证 TS 对应正确的立体异构体
2. **原子映射零成本**: 所有结构由产物坐标直接变形得到，原子索引始终一致
3. **成功率更高**: "Unzipping"(解拉链) 总是比 "Assembling"(组装) 容易收敛

### 1.4 核心公式

$$\Delta G^\ddagger = G_{TS} - G_{Reactant}$$

$$E_{dist} = E(R_{distorted}) - E(R_{relaxed})$$

$$E_{int} = \Delta E^\ddagger - E_{dist}$$

---

## 2. 串行四步架构

### 2.1 架构概览 (Sequential Workflow)

```
┌─────────────────────────────────────────────────────────────────────────┐
│              Reaction Profile Hunter v2.1 (Serial Architecture)        │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│    输入: SMILES_Product (产物 SMILES)                                  │
│                           │                                             │
│                           ▼                                             │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │           Step 1: Product_Anchor (产物锚定)                      │  │
│  │  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐               │  │
│  │  │ RDKit 3D   │→ │ CREST 搜索 │→ │ DFT 优化   │               │  │
│  │  └─────────────┘  └─────────────┘  └─────────────┘               │  │
│  │                           │                                      │  │
│  │                           ▼                                      │  │
│  │               Product_Global_Min.xyz                             │  │
│  └──────────────────────────┬───────────────────────────────────────┘  │
│                             │ (串行阻塞)                               │
│                             ▼                                          │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │           Step 2: Retro_Scanner (逆向扫描)                       │  │
│  │  ┌─────────────────────────────────────────────────────────┐     │  │
│  │  │              基于 Product_Min 分叉                      │     │  │
│  │  └─────────────────────────┬───────────────────────────────┘     │  │
│  │              ┌─────────────┴─────────────┐                       │  │
│  │              ▼                           ▼                       │  │
│  │     ┌────────────────┐          ┌────────────────┐               │  │
│  │     │ 路径A: 拉伸   │          │ 路径B: 断键   │               │  │
│  │     │ 键→2.2Å      │          │ 键→3.5Å→松弛 │               │  │
│  │     │ (受限优化)    │          │ (无限制优化)  │               │  │
│  │     └───────┬────────┘          └───────┬────────┘               │  │
│  │             ▼                           ▼                        │  │
│  │     TS_Guess.xyz               Reactant_Complex.xyz              │  │
│  └──────────────────────────┬───────────────────────────────────────┘  │
│                             │ (串行阻塞)                               │
│                             ▼                                          │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │           Step 3: TS_Optimizer (TS 精准优化)                     │  │
│  │  ┌─────────────┐     ┌─────────────┐     ┌─────────────┐         │  │
│  │  │ Berny TS   │ ──→ │ 虚频检查   │ ──→ │ IRC 验证   │         │  │
│  │  │ Opt=TS     │     │ (1个虚频?) │     │ (路径确认) │         │  │
│  │  └─────────────┘     └──────┬──────┘     └─────────────┘         │  │
│  │                             │                                    │  │
│  │              ┌──────────────┴──────────────┐                     │  │
│  │              │ 失败? → QST2 救援           │                     │  │
│  │              │ (有 R 和 P 两个端点)        │                     │  │
│  │              └─────────────────────────────┘                     │  │
│  │                             │                                    │  │
│  │                             ▼                                    │  │
│  │                     TS_Final.xyz                                 │  │
│  └──────────────────────────┬───────────────────────────────────────┘  │
│                             │ (串行阻塞)                               │
│                             ▼                                          │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │           Step 3.5: High_Precision_SP (高精度单点矩阵)           │  │
│  │  ┌─────────────┐     ┌─────────────┐     ┌─────────────┐         │  │
│  │  │ 覆盖矩阵生成 │ ──→ │ 批量 SP 提交 │ ──→ │ Checkpoint复用│         │  │
│  │  └─────────────┘     └──────┬──────┘     └─────────────┘         │  │
│  │                             │                                    │  │
│  │                             ▼                                    │  │
│  │                   [High_Level_SP_Reports]                        │  │
│  └──────────────────────────┬───────────────────────────────────────┘  │
│                             │ (串行阻塞)                               │
│                             ▼                                          │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │           Step 4: Feature_Miner (特征提取 - 四模块架构)          │  │
│  │  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐               │  │
│  │  │ Electronic  │  │ FMO indices │  │ Steric/Geo  │               │  │
│  │  └─────────────┘  └─────────────┘  └─────────────┘               │  │
│  │                           │                                      │  │
│  │                           ▼                                      │  │
│  │                    features.csv                                  │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

### 2.2 步骤定义

| Step | 名称 | 英文全称 | 依赖 | 核心职责 |
|:---:|:---:|---|:---:|---|
| **S1** | Product_Anchor | Product Conformer Anchor | 输入 | 确定产物全局最低构象 (L1 Opt + L2 SP) |
| **S2** | Retro_Scanner | Retro Transition State Scanner | S1 | 逆向拉伸生成 TS 初猜 + 底物 |
| **S3** | TS_Optimizer | Transition State Optimizer | S2 | Berny 优化 + IRC 验证 |
| **S3.5** | **SP_Matrix** | **High-Precision SP Matrix** | S1,S2,S3 | 强制高精度单点校验 (**ORCA RI 加速**) |
| **S4** | Feature_Miner | Distortion/Interaction Feature Miner | S1-S3.5 | 四模块解耦特征提取 |

### 2.3 数据流图

```
                    SMILES_Product
                          │
                          ▼
                ┌─────────────────┐
                │  Step 1: S1     │
                │ Product_Anchor  │ ──┐
                └────────┬────────┘   │
                         │            │
                         ▼            │
              Product_Global_Min.xyz ─┼───────────────────┐
                         │            │                   │
                         ▼            │                   │
                ┌─────────────────┐   │                   │
                │  Step 2: S2     │   │                   │
                │ Retro_Scanner   │   │                   │
                └────────┬────────┘   │                   │
                         │            │                   │
           ┌─────────────┴─────────────┐                  │
           ▼                           ▼                  │
    TS_Guess.xyz              Reactant_Complex.xyz        │
           │                           │                  │
           └─────────────┬─────────────┘                  │
                         ▼                                │
                ┌─────────────────┐                       │
                │  Step 3: S3     │                       │
                │ TS_Optimizer    │                       │
                └────────┬────────┘                       │
                         │                                │
                         ▼                                │
                  TS_Final.xyz                            │
                         │                                │
           ┌─────────────┴────────────────────────────────┘
           ▼
  ┌──────────────────┐
  │  Step 3.5: S3.5  │ (High-Precision SP Matrix - ORCA)
  │    SP_Matrix     │
  └────────┬─────────┘
           │
           ▼
  [SP_Report_Matrix] ──────────────────────────┐
           │                                   │
           ▼                                   ▼
  ┌─────────────────┐                 ┌─────────────────┐
  │  Step 4: S4     │                 │ 机器学习特征集  │
  │ Feature_Miner   │ ──────────────→ │ features.csv    │
  └─────────────────┘                 └─────────────────┘
```

---

## 3. 精度策略：双层级协议 (Dual-Level Protocol)

为了确保物理特征的准确性，v2.1 引入了 **Dual-Level (DL)** 协议，并引入 **ORCA** 作为高精度计算后端以利用其 **RI/RIJCOSX** 加速特性。

### 3.1 精度配置

| 级别 | 泛函/基组 | 软件后端 | 适用阶段 | 目的 |
|:---:|---|:---:|---|---|
| **L1 (Low)** | B3LYP/def2-SVP | Gaussian | S1, S2, S3 | 几何优化、路径预探索 |
| **L2 (High)** | M06-2X/def2-TZVPP | **ORCA** | S1, S3.5 | **能量基准、特征提取核心** |

> **为什么选择 ORCA?**: 对于 L2 层级（如 def2-TZVPP），传统的 Gaussian 计算开销较大。ORCA 的 **RIJCOSX** 技术能显著加速双电子积分计算，结合辅助基组（/J, /C）在保证精度的前提下提升 5-10 倍计算效率。

### 3.2 SP 覆盖矩阵 (SP Coverage Matrix)

在 **S1** 和 **Step 3.5** 中，系统会自动生成并执行 L2 级别的单点能任务：

| 结构节点 | 几何来源 | 任务类型 | 物理用途 |
|---|---|---|---|
| **Product** | **S1 (DFT-OPT)** | **ORCA SP (L2)** | **产物绝对能量基准 $G_{prod}$** |
| **Reactant** | S2 (无限制OPT) | ORCA SP (L2) | 活化能/反应能基准 $G_{rea}$ |
| **Cation/Anion**| 独立构建 | ORCA SP (L2) | 计算 HBDE (异裂解离能) |
| **TS_Guess** | S2 (受限OPT) | ORCA SP (L2) | 畸变能计算的碎片基准 (Fragment A/B) |
| **TS_Final** | S3 (Berny TS) | ORCA SP (L2) | 活化能 $\Delta G^\ddagger$ (L2级别) |

**关键规则**: Step 4 特征读取时，**强制优先读取 S3.5 的高精度能量结果**，而非 S1-S3 的优化级别能量。

---

## 4. 详细实施细节

### 3.1 Step 1: Product_Anchor (产物锚定)

> **"一切始于终点。"**

#### 设计原理

[5+2] 产物通常是刚性的桥环体系，但仍存在构象异构（如侧链取向、环的皱褶）。
必须先找到能量最低的产物构象，因为这是反应坐标的**终点**，也是逆推 TS 的**起点**。

#### 输入/输出

```
INPUT:  SMILES_Product (产物 SMILES)
OUTPUT: 
  - Product_Global_Min.xyz   (产物全局最低构象)
  - E_product                (产物能量)
```

#### 处理流程

```python
def step1_product_anchor(product_smiles: str) -> Tuple[Path, float]:
    """
    Step 1: 产物锚定工作流
    
    Steps:
    1. RDKit 3D 嵌入 (ETKDGv3)
    2. CREST 全局搜索 (GFN2-xTB)
    3. isostat 聚类 (能量窗口 6 kcal/mol)
    4. DFT 优化 + 频率 (验证无虚频)
    5. 返回全局最低构象
    """
    # Step 1: 3D 嵌入
    mol = Chem.MolFromSmiles(product_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    product_3d = save_to_xyz(mol, "product_initial.xyz")
    
    # Step 2-3: CREST + 聚类
    conformers = run_crest_confsearch(
        product_3d, 
        method="gfn2", 
        solvent="acetone",
        energy_window=6.0  # kcal/mol
    )
    
    # Step 4: DFT 优化 (仅最低能构象)
    product_min = run_dft_opt_freq(
        conformers[0],  # 最低能构象
        method="B3LYP-D3BJ/def2-SVP",
        verify_no_imaginary=True  # 确认是极小点
    )
    
    # Step 5: [NEW] L2 高精度单点能 (ORCA)
    e_product_l2 = run_orca_sp(
        product_min,
        method="M062X/def2-TZVPP",
        accel="RIJCOSX"
    )
    
    return product_min, e_product_l2

# 串行阻塞：Step 2 必须等待 Step 1 完成
product_global_min, e_product_l2 = step1_product_anchor(smiles)
```

#### 关键参数

```yaml
product_anchor:
  rdkit:
    embedding: "ETKDGv3"
    force_field: "MMFF"
  crest:
    method: "gfn2"
    solvent: "acetone"
    energy_window: 6.0  # kcal/mol
    threads: 16
  dft:
    method: "B3LYP"
    basis: "def2-SVP"
    dispersion: "D3BJ"
    freq: true  # 验证无虚频
```

#### 现有代码映射

| 功能 | 现有文件 | 函数/类 | 状态 |
|---|---|---|:---:|
| RDKit 3D | `retrots_main.py` | `_phase_a_prepare_product()` | ✅ 可用 |
| CREST 调用 | `confsearch.sh` | `run_crest_phase()` | ✅ 可用 |
| isostat 聚类 | `confsearch.sh` | `run_clustering()` | ✅ 可用 |
| DFT 优化 | `opt.sh` | 包装 OPT-Freq | ✅ 可用 |

---

### 3.2 Step 2: Retro_Scanner (逆向扫描与生成)

> **"从山顶（产物）沿山脊倒推找到鞍点（TS）。"**

#### 设计原理

这是脚本逻辑改变**最大**的地方。我们不再试图把两个分离的分子（底物+试剂）凑在一起，
而是把产物"拆"开：

```
传统方法 (❌ 易失败):           逆向方法 (✅ 推荐):
底物 + 试剂 → 拼凑 → TS        产物 → 拉伸 → TS → 松弛 → 底物
```

#### 逆向拉伸 (Retro-Stretching) 原理

```
        产物 (P_min)                    TS 初猜                    底物复合物
        ┌───O───┐                     ┌───O───┐                   ┌───O───┐
       /   ‖     \                   /   ‖     \                 /   ‖     \
      C─────────C                   C ← ─ ─ → C                 C           C
      │  1.54Å  │                   │  2.2Å   │                 │  >3.5Å    │
      R─────────R                   R─────────R                 R───────────R
      
      形成键: ~1.5Å          →     形成键: ~2.2Å          →     形成键: 断裂
      (已成键)                      (TS 典型距离)                 (未成键)
```

#### 输入/输出

```
INPUT:  Product_Global_Min.xyz (来自 Step 1)
OUTPUT: 
  - TS_Guess.xyz           (TS 初猜，保留产物立体化学)
  - Reactant_Complex.xyz   (底物复合物，用于能量参考)
```

#### 处理流程

```python
def step2_retro_scanner(product_min: Path) -> Tuple[Path, Path]:
    """
    Step 2: 逆向扫描工作流
    
    Steps:
    1. SMARTS 识别 [5+2] 形成的两个 σ 键
    2. 路径 A: 拉伸至 TS 距离 (2.2 Å) + 受限优化
    3. 路径 B: 拉伸至断键 (3.5 Å) + 无限制优化
    """
    # 读取产物结构
    coords, symbols = read_xyz(product_min)
    mol = xyz_to_rdkit_mol(product_min)
    
    # Step 1: SMARTS 识别形成键
    match_result = smarts_engine.find_reactive_bonds(mol, coords)
    bond_1 = (match_result.bond_1.atom_idx_1, match_result.bond_1.atom_idx_2)
    bond_2 = (match_result.bond_2.atom_idx_1, match_result.bond_2.atom_idx_2)
    
    print(f"识别到形成键: {bond_1} (当前: {match_result.bond_1.current_length:.2f} Å)")
    print(f"识别到形成键: {bond_2} (当前: {match_result.bond_2.current_length:.2f} Å)")
    
    # ==========================================
    # 路径 A: 生成 TS 初猜 (拉伸至 2.2 Å)
    # ==========================================
    ts_guess = geometry_ops.stretch_and_constrain(
        coords, symbols,
        bonds=[bond_1, bond_2],
        target_distance=2.2,  # TS 典型距离
        constraints=True      # 保持距离约束
    )
    
    # 受限优化 (xTB)
    ts_guess_opt = xtb_constrained_opt(
        ts_guess,
        constraints=[
            f"$constrain",
            f"  distance: {bond_1[0]+1}, {bond_1[1]+1}, 2.2",
            f"  distance: {bond_2[0]+1}, {bond_2[1]+1}, 2.2",
            f"$end"
        ]
    )
    
    # ==========================================
    # 路径 B: 生成底物 (拉伸至 3.5 Å + 松弛)
    # ==========================================
    stretched = geometry_ops.stretch_bonds(
        coords, symbols,
        bonds=[bond_1, bond_2],
        target_distance=3.5  # 断裂距离
    )
    
    # 无限制优化
    reactant_complex = xtb_free_opt(stretched)
    
    return (ts_guess_opt, reactant_complex)

# 串行阻塞：Step 3 必须等待 Step 2 完成
ts_guess, reactant_complex = step2_retro_scanner(product_global_min)
```

#### 关键优势

| 优势 | 说明 |
|---|---|
| **立体化学保真** | TS 初猜保留了产物的手性中心，不会产生差向异构 |
| **原子映射零成本** | 所有结构共享相同的原子索引，无需重新映射 |
| **高成功率** | "解拉链"比"拼凑"更容易收敛 |

#### 现有代码映射

| 功能 | 现有文件 | 函数/类 | 状态 |
|---|---|---|:---:|
| SMARTS 匹配 | `smarts_engine.py` | `OxabicycloSMARTSEngine` | ✅ 可用 |
| 键拉伸 | `geometry_ops.py` | `BondStretcher.stretch_two_bonds()` | ✅ 可用 |
| 受限优化 | `constrained_opt.py` | `ConstrainedOptimizer` | ✅ 可用 |
| 无限制优化 | `qc_interface.py` | `XTBInterface.optimize()` | ✅ 可用 |

---

### 3.3 Step 3: TS_Optimizer (过渡态精准优化)

> **"在确定的路径上精修鞍点。"**

#### 设计原理

由于 `TS_Guess` 是从产物逆推而来的，它已经**非常接近真实的 TS**。
如果 Berny 优化失败，我们可以使用 **QST2 策略**，因为我们已经同时获得了：
- `Product_Min` (Step 1)
- `Reactant_Complex` (Step 2)

这两个确定的端点是 QST2 算法的完美输入。

#### 输入/输出

```
INPUT:  
  - TS_Guess.xyz           (来自 Step 2)
  - Reactant_Complex.xyz   (来自 Step 2, 备用)
  - Product_Global_Min.xyz (来自 Step 1, 备用)
OUTPUT: 
  - TS_Final.xyz           (验证后的 TS)
  - imaginary_frequency    (虚频值)
  - IRC_path.xyz           (可选: IRC 路径)
```

#### 处理流程

```python
def step3_ts_optimizer(
    ts_guess: Path,
    reactant: Path,
    product: Path
) -> Tuple[Path, float]:
    """
    Step 3: TS 精准优化工作流
    
    Steps:
    1. Berny TS 优化 (Opt=TS, CalcFC, NoEigenTest)
    2. 虚频检验 (必须恰好 1 个虚频)
    3. IRC 验证 (可选但推荐)
    4. 失败救援: QST2 策略
    """
    
    # ==========================================
    # 主策略: Berny TS 优化
    # ==========================================
    try:
        ts_result = gaussian_ts_opt(
            ts_guess,
            keywords="Opt=(TS, CalcFC, NoEigenTest) Freq",
            method="B3LYP-D3BJ/def2-SVP"
        )
        
        # 虚频检验
        imaginary_freqs = [f for f in ts_result.frequencies if f < 0]
        
        if len(imaginary_freqs) != 1:
            raise TSValidationError(
                f"期望 1 个虚频，实际 {len(imaginary_freqs)} 个"
            )
        
        print(f"✓ TS 优化成功，虚频 = {imaginary_freqs[0]:.1f} cm⁻¹")
        
    except (OptimizationError, TSValidationError) as e:
        print(f"Berny 优化失败: {e}")
        print(">>> 启用 QST2 救援策略...")
        
        # ==========================================
        # 救援策略: QST2
        # ==========================================
        ts_result = gaussian_qst2(
            reactant=reactant,
            product=product,
            keywords="Opt=(QST2, CalcFC) Freq",
            method="B3LYP-D3BJ/def2-SVP"
        )
        
        imaginary_freqs = [f for f in ts_result.frequencies if f < 0]
        
        if len(imaginary_freqs) != 1:
            raise FatalError("QST2 也失败了，需要人工干预")
    
    # ==========================================
    # 可选: IRC 验证
    # ==========================================
    if VERIFY_IRC:
        irc_result = gaussian_irc(
            ts_result.xyz,
            keywords="IRC=(CalcFC, MaxPoints=50, StepSize=10)",
            direction="both"
        )
        
        # 验证 IRC 连接正确的端点
        validate_irc_endpoints(irc_result, reactant, product)
    
    return (ts_result.xyz, imaginary_freqs[0])

# 串行阻塞：Step 3.5 必须等待 Step 3 完成
ts_final, imag_freq = step3_ts_optimizer(ts_guess, reactant_complex, product_global_min)

---

### 3.3.5 Step 3.5: High_Precision_SP (高精度单点矩阵 - ORCA 驱动)

> **"几何决定定性，能量决定定量。"**

#### 设计原理
在最终特征提取前，对 S1-S3 产生的关键几何点进行 M06-2X/def2-TZVPP 级别的 SP 计算。为了应对大规模计算，采用 **ORCA** 作为主要后端，利用其 **RIJCOSX** 算法加快计算速度。

#### 处理流程
1. **收集节点**: 获取 `Product_Min`, `Reactant_Complex`, `TS_Final`, 以及由 TS 拆解得到的 `Frag_A_at_TS`, `Frag_B_at_TS`。
2. **生成 ORCA 输入**: 使用 `! M062X def2-TZVPP def2/J RIJCOSX` 路由，并根据泛函类型自动匹配辅助基组。
3. **RIJCOSX 加速**: 启用 `! RIJCOSX` 与 `GridX` 系列关键字，确保双电子积分在大基组下的效率。
4. **解析结果**: 通过 ORCA 输出解析能量值，统一整合为 `SP_Matrix_Report`。

---

### 3.4 Step 4: Feature_Miner (特征提取 - 四模块重构)

> **"解构复杂性，职责分离。"**

#### 4.1 目录结构
为了应对 25+ 个特征的复杂性，Step 4 弃用单一文件，采用以下解耦架构：

```text
step4_features/
├── __init__.py
├── feature_miner.py      # [主控] 协调各模块，汇总 CSV
├── electronic.py         # [电子] HBDE, NBO, 电荷分布
├── fmo_reactivity.py     # [前线轨道] HOMO/LUMO, Parrish指数
├── steric_geometry.py    # [立体/几何] 畸变能, Sterimol, Vbur
└── entropy.py            # [熵/热力学] 振动熵, 构象布居
```

#### 4.2 核心算法伪代码

**3.4.1 HBDE 计算 (electronic.py)**
```python
def calculate_hbde(e_precursor: float, e_cation: float, e_anion: float) -> float:
    """
    异裂键解离能 (HBDE) = E(Cation) + E(Anion) - E(Precursor)
    必须使用相同层级的 L2 SP 能量。
    """
    return (e_cation + e_anion - e_precursor) * 627.5095
```

**3.4.2 双片段畸变能 (steric_geometry.py)**
```python
def calculate_dual_fragment_distortion(e_ts_sp: float, e_frag_a_ts: float, e_frag_b_ts: float, e_rea_relaxed: float) -> float:
    """
    E_dist_total = (E_frag_A_at_TS - E_frag_A_relaxed) + (E_frag_B_at_TS - E_frag_B_relaxed)
    简化版: E_dist = E_frag_A_at_TS + E_frag_B_at_TS - E_reactant_relaxed
    """
    return (e_frag_a_ts + e_frag_b_ts - e_rea_relaxed) * 627.5095
```

**3.4.3 NICS Ghost 原子添加 (geometry_tools.py)**
```python
def add_ghost_atom_at_centroid(coords: np.ndarray, ring_atom_indices: list, distance: float = 1.0):
    """
    SVD 拟合环平面，并在法线 1.0Å 处添加 Ghost 原子用于 NICS(1) 计算。
    """
    ring_coords = coords[ring_atom_indices]
    centroid = np.mean(ring_coords, axis=0)
    _, _, vh = np.linalg.svd(ring_coords - centroid)
    normal = vh[2] # 最小主轴即法线
    return centroid + normal * distance
```

#### 4.3 外部工具集成
- **morfeus**: 用于计算 Sterimol ($B_1, B_5, L$) 参数及 Buried Volume ($V_{bur}$)。
- **cclib**: 统一不同量化软件的 log 解析接口。

---

## 4. 现有代码映射与 GAP 分析 (Updated)

### 4.1 开发状态概览

| 功能模块 | 归属步骤 | 状态 | 关键缺失 |
|---|---|---|---|
| Product Optimization | S1 | ✅ | 无需修改 |
| TS Guess Generation | S2 | ✅ | 需确保产物锚定逻辑 |
| TS Opt (Berny) | S3 | ⚠️ | 需增加 QST2 自动救援 |
| High-Precision SP | S3.5 | ❌ | **NEW** 需实现单点矩阵 |
| Electronic Features | S4 | ⚠️ | HBDE 计算需对齐 L2 能量 |
| FMO Features | S4 | ⚠️ | Parrish 指数提取 |
| Steric Features | S4 | ⚠️ | Morfeus Sterimol 集成 |

---

## 5. 任务清单 (Phase 2.1)

### 5.1 Step 3.5: SP_Matrix 核心开发 [P0]
1. [ ] 创建 `rph_core/steps/step3_5_sp/`
2. [ ] 实现 Gaussian SP GJF 生成器 (支持 `%oldchk`)
3. [ ] 实现 SP 任务批处理器

### 5.2 Step 4: 四模块重构 [P0]
1. [ ] 完成 `electronic.py`: HBDE, NICS, Charges
2. [ ] 完成 `steric_geometry.py`: Distortion, Sterimol
3. [ ] 完成 `fmo_reactivity.py`: HOMO/LUMO, Orbital Analysis
4. [ ] 完成 `entropy.py`: Vibration/Entropy extraction

### 5.3 工程端 [P1]
1. [ ] 支持 Batch CSV 输入
2. [ ] 实现 JSON Checkpoint 机制 (S1-S4 意外中断恢复)

> **v2.0 vs v2.1 对比**: 
> - v2.0 需要 ~60h (含复杂的模板对齐和并行协调)
> - v2.1 仅需 ~42h (复用大量现有代码)
| **总计** | **28h** | **24h** | **8h** | **60h** |

---

## 6. 测试策略

### 6.1 单元测试

每个模块独立测试：

```python
# test_gs_anchor.py
def test_crest_confsearch():
    """测试 CREST 构象搜索"""
    result = run_crest_confsearch("c1ccccc1", method="gfn0")  # 苯，快速测试
    assert len(result.conformers) >= 1
    assert result.lowest_energy is not None

def test_boltzmann_weighting():
    """测试 Boltzmann 加权"""
    energies = [0.0, 1.0, 2.0]  # kcal/mol
    weights = boltzmann_weighting(energies, T=298.15)
    assert abs(sum(weights) - 1.0) < 1e-6
```

```python
# test_ts_generator.py
def test_smarts_matching():
    """测试 SMARTS 匹配"""
    smiles = "O=C1CC2CCCCC2O1"  # 氧杂桥环产物
    result = smarts_engine.find_reactive_bonds(smiles)
    assert result.matched
    assert result.bond_1 is not None
    assert result.bond_2 is not None

def test_retro_stretching():
    """测试逆向拉伸 (v2.1 核心功能)"""
    product = create_5plus2_product()
    ts_guess, reactant = retro_scanner_workflow(product)
    
    # TS 初猜的形成键距离应在 2.0-2.4 Å
    r1, r2 = get_forming_bond_lengths(ts_guess)
    assert 2.0 < r1 < 2.4
    assert 2.0 < r2 < 2.4
    
    # 底物的"形成键"应已断裂 (>3.0 Å)
    r1_r, r2_r = get_forming_bond_lengths(reactant)
    assert r1_r > 3.0
    assert r2_r > 3.0
```

```python
# test_ts_optimizer.py
def test_ts_validation():
    """测试虚频验证"""
    freqs = [-450.2, 50.1, 120.3, 200.5]  # 一个虚频
    assert validate_imaginary_freq(freqs) == True
    
    freqs_bad = [-450.2, -30.5, 120.3, 200.5]  # 两个虚频
    assert validate_imaginary_freq(freqs_bad) == False

def test_qst2_fallback():
    """测试 QST2 救援策略"""
    # Berny 失败后应自动切换到 QST2
    reactant = Path("test_data/reactant.xyz")
    product = Path("test_data/product.xyz")
    ts_guess = Path("test_data/bad_ts_guess.xyz")  # 故意用差的初猜
    
    ts_final = ts_optimizer_workflow(ts_guess, reactant, product)
    assert ts_final.method_used in ["Berny", "QST2"]
```

```python
# test_feature_miner.py
def test_distortion_energy():
    """测试畸变能计算"""
    E_relaxed = -100.0  # Hartree
    E_distorted = -99.95  # Hartree
    E_dist = calculate_distortion_energy(E_relaxed, E_distorted)
    assert abs(E_dist - 31.4) < 1.0  # ~31 kcal/mol

def test_asynchronicity():
    """测试非同步性指数"""
    r1, r2 = 2.1, 2.3  # Å
    async_idx = calculate_asynchronicity(r1, r2)
    assert 0 < async_idx < 1
```

### 6.2 集成测试

端到端测试完整串行流程：

```python
# test_integration.py
def test_full_pipeline_serial():
    """v2.1 串行流程全测试"""
    product_smiles = "O=C1CC2CCCCC2O1"  # [5+2] 氧杂桥环产物
    
    # ========================================
    # Step 1: Product Anchor (阻塞)
    # ========================================
    product_min = step1_product_anchor(product_smiles)
    assert product_min.exists()
    
    # ========================================
    # Step 2: Retro-Scanner (依赖 Step 1)
    # ========================================
    ts_guess, reactant = step2_retro_scanner(product_min)
    assert ts_guess.exists()
    assert reactant.exists()
    
    # 验证原子数量一致 (v2.1 关键优势!)
    n_atoms_product = count_atoms(product_min)
    n_atoms_ts = count_atoms(ts_guess)
    n_atoms_reactant = count_atoms(reactant)
    assert n_atoms_product == n_atoms_ts == n_atoms_reactant
    
    # ========================================
    # Step 3: TS Optimizer (依赖 Step 2)
    # ========================================
    ts_final, imag_freq = step3_ts_optimizer(
        ts_guess, reactant, product_min
    )
    assert ts_final.exists()
    assert imag_freq < -100  # 合理的虚频值
    
    # ========================================
    # Step 4: Feature Miner (依赖 Step 1-3)
    # ========================================
    features = step4_feature_miner(product_min, reactant, ts_final)
    
    # 验证关键特征
    assert "dG_activation" in features.columns
    assert features["dG_activation"].iloc[0] > 0  # 活化能应为正
    assert 0 < features["asynchronicity"].iloc[0] < 1
```

### 6.3 测试数据

```
test_data/
├── products/                      # v2.1: 从产物开始
│   ├── oxabicyclo_1.xyz           # 氧杂桥环产物示例 1
│   ├── oxabicyclo_2.xyz           # 氧杂桥环产物示例 2
│   └── oxabicyclo_fused.xyz       # 稠环产物示例
├── reference/
│   ├── ts_reference_1.xyz         # 参考 TS 结构
│   ├── reactant_reference_1.xyz   # 参考底物结构
│   └── features_reference_1.csv   # 参考特征值
└── configs/
    └── test_config.yaml           # 测试配置
```

---

## 7. 开发路线图 (v2.1 串行架构)

### 7.1 Phase 1: Product Anchor (3天)

```
目标: 完成 Step 1 产物锚定模块
任务:
[x] 1.1 CREST 调用封装 (已有 confsearch.sh)
[x] 1.2 isostat 聚类 (已有)
[x] 1.3 DFT 优化 (已有 opt.sh)
[ ] 1.4 verify_no_imaginary() 虚频检查
[ ] 1.5 编写 test_product_anchor.py
```

### 7.2 Phase 2: Retro-Scanner (3天)

```
目标: 完成 Step 2 逆向扫描模块
任务:
[x] 2.1 SMARTS 匹配 (已有 smarts_engine.py)
[x] 2.2 键拉伸 (已有 geometry_ops.py)
[x] 2.3 受限优化 (已有 constrained_opt.py)
[ ] 2.4 generate_dual_outputs() 同时生成 TS_Guess + Reactant
[ ] 2.5 编写 test_retro_scanner.py
```

### 7.3 Phase 3: TS Optimizer (1周)

```
目标: 完成 Step 3 过渡态优化模块 (核心开发)
任务:
[ ] 3.1 gaussian_ts_opt() Berny TS 接口
[ ] 3.2 validate_imaginary_freq() 虚频验证
[ ] 3.3 gaussian_qst2() QST2 救援策略
[ ] 3.4 gaussian_irc() IRC 验证 (可选)
[ ] 3.5 编写 test_ts_optimizer.py
```

### 7.4 Phase 4: Feature Miner (4天)

```
目标: 完成 Step 4 特征提取模块
任务:
[ ] 4.1 extract_fragment() 从 TS 提取分子片段
[ ] 4.2 calculate_distortion_energy() 畸变能计算
[ ] 4.3 calculate_asynchronicity() 非同步性指数
[ ] 4.4 run_nbo_analysis() NBO 接口 (可选)
[ ] 4.5 CSV 输出格式化
[ ] 4.6 编写 test_feature_miner.py
```

### 7.5 Phase 5: 集成验证 (3天)

```
目标: 端到端集成与真实数据验证
任务:
[ ] 5.1 ReactionProfileHunter 主类实现
[ ] 5.2 编写 test_integration.py
[ ] 5.3 使用真实 [5+2] 反应数据验证
[ ] 5.4 性能基准测试
[ ] 5.5 文档更新
```

### 7.6 v2.1 优势总结

| 对比项 | v2.0 并行架构 | v2.1 串行架构 |
|---|:---:|:---:|
| 化学正确性 | ❌ 可能产生错误构象 | ✅ Hammond 原理保证 |
| 原子映射 | 需要复杂对齐 | 零成本 (共享索引) |
| 新增代码量 | ~60h | ~42h |
| 现有代码复用 | 60% | **80%** |
| 开发周期 | 6周 | **3周** |

---

## 8. 附录

### 8.1 Gaussian 输入模板

#### TS 优化

```
%chk=ts.chk
%mem=32GB
%nproc=16
# B3LYP/def2-SVP EmpiricalDispersion=GD3BJ 
# Opt=(TS, CalcFC, NoEigenTest) Freq

TS Optimization for [5+2] Cycloaddition

0 1
[geometry]
```

#### IRC 计算

```
%chk=ts.chk
%mem=32GB
%nproc=16
# B3LYP/def2-SVP EmpiricalDispersion=GD3BJ
# IRC=(CalcFC, MaxPoints=50, StepSize=10) Geom=Check Guess=Read

IRC from TS

0 1

```

### 8.2 ORCA 输入模板

#### TS 优化

```
! B3LYP D3BJ def2-SVP TightSCF OptTS Freq

%geom
  Calc_Hess true
  Recalc_Hess 5
end

%pal nprocs 16 end
%maxcore 2000

* xyzfile 0 1 ts_guess.xyz
```

---

## 9. 🏗️ Python 骨架重构方案 (Modularization Surgery)

> **设计理念**: 从 "按工具分类" 转变为 "按步骤分类" (One Folder, One Step)

### 9.1 问题诊断：当前架构的 "Growing Pains"

#### 当前物理文件结构 vs 逻辑架构的错位

```
┌───────────────────────────────────────────────────────────────────────────┐
│                     ❌ 当前结构 (按工具分类)                              │
├───────────────────────────────────────────────────────────────────────────┤
│                                                                           │
│  ConfsearchIII/                                                           │
│  ├── retrots_main.py         # 混合了 S1+S2 的逻辑                       │
│  ├── qc_interface.py         # 底层工具 (无业务逻辑)                      │
│  ├── descriptor_extractor.py # S4 的部分逻辑                              │
│  └── QCDescriptors/                                                       │
│      └── core/                                                            │
│          ├── smarts_engine.py      # S2 的一部分                         │
│          ├── geometry_ops.py       # S2 的一部分                         │
│          ├── constrained_opt.py    # S2 的一部分                         │
│          ├── validators.py         # S3 的一部分                         │
│          └── qc_interface.py       # 底层 (重复!)                        │
│                                                                           │
│  问题:                                                                    │
│  1. S2 的逻辑散落在 3+ 个文件中                                          │
│  2. retrots_main.py 承担了太多职责 (God Object)                          │
│  3. S3 (TS 优化) 的代码几乎不存在                                        │
│  4. 调试某个 Step 时，需要在多个目录间跳转                               │
│                                                                           │
└───────────────────────────────────────────────────────────────────────────┘
```

#### 目标架构：物理隔离 + 职责明确

```
┌───────────────────────────────────────────────────────────────────────────┐
│                     ✅ 目标结构 (按步骤分类)                              │
├───────────────────────────────────────────────────────────────────────────┤
│                                                                           │
│  ReactionProfileHunter/                                                   │
│  ├── bin/rph_run                 # 单一入口                              │
│  ├── rph_core/                                                            │
│  │   ├── orchestrator.py         # 总指挥 (无业务逻辑)                   │
│  │   ├── steps/                                                          │
│  │   │   ├── step1_anchor/       # S1 所有代码在这                       │
│  │   │   ├── step2_retro/        # S2 所有代码在这                       │
│  │   │   ├── step3_opt/          # S3 所有代码在这                       │
│  │   │   └── step4_features/     # S4 所有代码在这                       │
│  │   └── utils/                  # 共享工具 (无业务逻辑)                 │
│  └── tests/                      # 按步骤组织的测试                      │
│                                                                           │
│  优势:                                                                    │
│  1. 调试 S2? 只看 steps/step2_retro/                                     │
│  2. 新增 S3 功能? 只改 steps/step3_opt/                                  │
│  3. 职责隔离: "技工" (utils) vs "化学家" (steps)                         │
│                                                                           │
└───────────────────────────────────────────────────────────────────────────┘
```

---

### 9.2 新项目结构蓝图

```
ReactionProfileHunter/
├── 📂 bin/                        # 命令行入口
│   └── rph_run                    # 启动命令 (entry point)
│
├── 📂 config/                     # 配置文件
│   ├── defaults.yaml              # 默认参数 (B3LYP, def2-SVP...)
│   └── templates/                 # 输入文件模板 (Gaussian/ORCA header)
│       ├── gaussian_ts.gjf
│       ├── gaussian_qst2.gjf
│       ├── gaussian_irc.gjf
│       └── orca_ts.inp
│
├── 📂 rph_core/                   # 核心源代码包 (Python Package)
│   ├── __init__.py
│   │
│   ├── 📜 orchestrator.py         # [总指挥] 串行工作流控制器
│   │
│   ├── 📂 steps/                  # [核心] 四大步骤模块
│   │   ├── __init__.py
│   │   │
│   │   ├── 📂 step1_anchor/       # Step 1: 产物锚定
│   │   │   ├── __init__.py
│   │   │   ├── anchor_manager.py  # 主逻辑: ProductAnchor 类
│   │   │   └── conf_searcher.py   # 调用 CREST/isostat
│   │   │
│   │   ├── 📂 step2_retro/        # Step 2: 逆向扫描
│   │   │   ├── __init__.py
│   │   │   ├── retro_scanner.py   # 主逻辑: RetroScanner 类
│   │   │   ├── smarts_matcher.py  # SMARTS 匹配 (从 smarts_engine.py)
│   │   │   └── bond_stretcher.py  # 键拉伸 (从 geometry_ops.py)
│   │   │
│   │   ├── 📂 step3_opt/          # Step 3: 精准优化
│   │   │   ├── __init__.py
│   │   │   ├── ts_optimizer.py    # 主逻辑: TSOptimizer 类
│   │   │   ├── berny_driver.py    # Berny TS 优化
│   │   │   ├── qst2_rescue.py     # QST2 救援策略
│   │   │   └── validator.py       # 虚频与 IRC 验证
│   │   │
│   │   └── 📂 step4_features/     # Step 4: 特征提取
│   │       ├── __init__.py
│   │       ├── feature_miner.py   # 主逻辑: FeatureMiner 类
│   │       ├── distortion_calc.py # 畸变能计算
│   │       └── nbo_analyzer.py    # NBO 分析 (可选)
│   │
│   └── 📂 utils/                  # [工具] 通用底层库
│       ├── __init__.py
│       ├── qc_interface.py        # Gaussian/ORCA/xTB 接口 (核心底层)
│       ├── geometry_tools.py      # 通用几何操作 (非业务逻辑)
│       ├── file_io.py             # 文件读写辅助 (XYZ, GJF, Log 解析)
│       └── log_manager.py         # 日志系统
│
├── 📂 tests/                      # 单元测试
│   ├── test_s1_anchor.py
│   ├── test_s2_retro.py
│   ├── test_s3_optimizer.py
│   ├── test_s4_features.py
│   └── test_integration.py
│
├── 📄 requirements.txt
├── 📄 setup.py
├── 📄 pyproject.toml
└── 📘 README.md
```

---

### 9.3 现有代码 → 新架构映射表

#### 9.3.1 文件迁移清单

| 现有文件 | 目标位置 | 迁移策略 |
|---|---|---|
| `ConfsearchIII/retrots_main.py` | 拆分 → `step1_anchor/` + `step2_retro/` | **拆解** |
| `ConfsearchIII/qc_interface.py` | `rph_core/utils/qc_interface.py` | **移动** |
| `ConfsearchIII/descriptor_extractor.py` | `rph_core/steps/step4_features/feature_miner.py` | **移动 + 重构** |
| `QCDescriptors/core/smarts_engine.py` | `rph_core/steps/step2_retro/smarts_matcher.py` | **移动 + 简化** |
| `QCDescriptors/core/geometry_ops.py` | 拆分 → `utils/geometry_tools.py` + `step2_retro/bond_stretcher.py` | **拆解** |
| `QCDescriptors/core/constrained_opt.py` | `step2_retro/retro_scanner.py` (合并) | **合并** |
| `QCDescriptors/core/validators.py` | `step3_opt/validator.py` | **移动** |
| `QCDescriptors/data/reaxys_loader.py` | `rph_core/utils/reaxys_loader.py` | **移动** |
| `QCDescriptors/core/parallel_processor.py` | `rph_core/utils/parallel_processor.py` | **移动** |
| `QCDescriptors/integration/bash_bridge.py` | `rph_core/utils/bash_bridge.py` | **移动** |

#### 9.3.2 详细拆解规划

##### `retrots_main.py` 拆解 (616 行 → 4 个模块)

```python
# ============== 当前 retrots_main.py 中的逻辑分布 ==============

class RetroTSPipeline:
    def process_single_molecule(self, ...):
        # Phase A: 产物准备 ─────────────────────────────> step1_anchor/
        product_xyz = self._phase_a_prepare_product(...)  # → anchor_manager.py
        
        # Phase B: 键位点识别 ───────────────────────────> step2_retro/
        match_result = self.smarts_engine.find_reactive_bonds(...)  # → smarts_matcher.py
        
        # Phase C: 逆向拉伸 + 受限优化 ──────────────────> step2_retro/
        prets_xyz = self._phase_c_stretch_and_optimize(...)  # → retro_scanner.py
        
        # Phase D: 结构验证 ─────────────────────────────> step3_opt/
        validation_result = self.validator.validate_structure(...)  # → validator.py
```

**拆解行动计划:**

| 行号范围 | 内容 | 迁移目标 |
|:---:|---|---|
| 1-55 | 导入 + 日志配置 | `rph_core/__init__.py` |
| 56-100 | `RetroTSPipeline.__init__` | `orchestrator.py` |
| 101-250 | `_phase_a_prepare_product()` | `step1_anchor/anchor_manager.py` |
| 251-320 | Phase B: SMARTS 匹配 | `step2_retro/smarts_matcher.py` |
| 321-450 | `_phase_c_stretch_and_optimize()` | `step2_retro/retro_scanner.py` |
| 451-500 | Phase D: 验证 | `step3_opt/validator.py` |
| 501-616 | 批量处理 + CLI | `bin/rph_run` |

##### `geometry_ops.py` 拆解 (597 行 → 2 个模块)

```python
# ============== 当前 geometry_ops.py 中的逻辑分布 ==============

# 业务逻辑 (与 Retro-TS 强相关) ────────────────> step2_retro/bond_stretcher.py
class BondStretcher:
    def stretch_bond(self, ...)        # TS 键拉伸
    def stretch_two_bonds(self, ...)   # [5+2] 双键拉伸
    def generate_scan_structures(...)  # 扫描结构生成

# 通用工具 (无业务逻辑) ─────────────────────────> utils/geometry_tools.py
class GeometryUtils:
    @staticmethod
    def calculate_distance(...)
    @staticmethod
    def calculate_angle(...)
    @staticmethod
    def calculate_dihedral(...)
    @staticmethod
    def get_com(...)  # 质心
```

---

## 10. 体系化补全迁移方案（2026-01-10）

本节系统性总结目前 ReactionProfileHunter v2.1 与 PROMOTE.md 设计的差异，并给出基于 Original_Eddition/Auto_Calc/ConfsearchIII 的体系化补全迁移方案，确保新架构既具备高通量批量处理能力，也能实现物理化学意义上的畸变/碎片能量分析。

### 一、差异与缺失点总结

1. **批量/并行处理能力缺失**
   - 现有 orchestrator 仅支持单分子串行处理，未能复用原有 Auto_Calc/ConfsearchIII 的批量与并行机制。
   - 原版 `parallel_processor.py`、`reaxys_loader.py` 提供了高效的多进程与数据集加载能力。

2. **畸变能/碎片能量分析未闭环**
   - Step 4 仅有能量差与几何特征，未实现“TS碎片提取→SP能量→碎片优化→畸变/相互作用能”全流程。
   - 原版 `descriptor_extractor.py`/Auto_Calc 已有成熟的碎片提取与能量分析逻辑。

3. **产物锚定精度与DFT精修**
   - 目前 CREST已接入，但DFT优化部分仍为占位，未实现高精度产物锚定。

4. **NBO/轨道相互作用分析缺失**
   - Gaussian/NBO相关的输入与解析尚未集成，无法提取轨道相互作用等高级电子特征。

### 二、体系化补全迁移方案

#### 1. 并行与批量处理能力迁移
- 迁移内容：
  - 将 `Original_Eddition/ConfsearchIII/parallel_processor.py` 迁移至 `rph_core/utils/parallel_processor.py`。
  - 将 `Original_Eddition/ConfsearchIII/reaxys_loader.py` 迁移至 `rph_core/utils/reaxys_loader.py`。
- 集成方式：
  - 在 orchestrator.py 增加 `run_batch()` 方法，支持批量 SMILES/CSV 输入，自动分批并行调度。
  - CLI 增加 `--batch` 参数，支持批量任务与单分子任务切换。

#### 2. 畸变能/碎片能量分析闭环
- 迁移内容：
  - 在 `rph_core/steps/step4_features/distortion_calculator.py` 增加 `extract_fragments()` 方法，实现基于 forming_bonds 的 TS 结构碎片化。
  - 在 `feature_miner.py` 中，调用 `GaussianInterface` 对碎片分别进行单点（SP）与优化（OPT）计算，获得 $E_{frag}^{TS}$ 与 $E_{frag}^{opt}$。
  - 复用原有 `calculate_distortion_interaction()` 公式，输出 Edist/ Eint/ Efrag 等物理特征。
- 实现思路：
  - forming_bonds 由 Step 2 传递，自动切分 TS 结构为两个碎片。
  - 先对碎片做 SP（冻结构型），再做 OPT（全松弛），与 Outcome/Auto_Calc 体系一致。

#### 3. 产物锚定与DFT精修
- 补全内容：
  - 在 `ProductAnchor` 中补全 DFT 优化分支，调用 `GaussianInterface.optimize()` 对 CREST 最优构象做高精度 DFT 优化。
  - 产物全局最低构象输出为 DFT 优化结果，提升后续反应路径精度。

#### 4. NBO/轨道相互作用分析
- 补全内容：
  - 在 `GaussianInterface` 增加 NBO 相关 route 支持与 log 解析。
  - 在 `feature_miner.py` 增加 NBO 特征提取分支，输出 donor-acceptor 相互作用等高级电子结构特征。

### 三、迁移与补全实施路线
1. 基础设施迁移
   - [x] 迁移 `parallel_processor.py`、`reaxys_loader.py` 至 utils。
   - [x] 在 orchestrator/CLI 增加批量模式支持。
2. 畸变能/碎片能量分析闭环
   - [x] 在 distortion_calculator/feature_miner 实现碎片提取与能量分析全流程。
   - [x] 复用 forming_bonds 自动切分与能量计算。
3. 产物锚定与DFT精修
   - [ ] 在 ProductAnchor 实现 DFT 优化分支，输出高精度产物构象。
4. NBO/轨道相互作用分析
   - [ ] 在 GaussianInterface/feature_miner 实现 NBO route 与特征提取。

### 四、参考原有体系的设计要点
- 目录分层与数据流：严格复用 Auto_Calc 的 OPT/SP/Outcome 目录与文件命名规范，便于后续数据追溯与批量管理。
- 批量与并行：所有批量任务均通过 ParallelProcessor/ChunkedParallelProcessor 调度，充分利用多核资源。
- 碎片能量分析：所有 TS/产物/底物均可自动切分为碎片，支持高通量畸变/相互作用能分析。
- 可扩展性：所有新功能均以 utils/steps/分层实现，便于后续维护与团队协作。

**结论**：本方案将最大程度复用原有体系的工程化与物理化学分析优势，补全 PROMOTE.md 设计与当前实现的所有差距，确保 ReactionProfileHunter v2.1 具备高通量、体系化、物理可信的自动化反应特征挖掘能力。


---

### 五、补充细节：已识别的潜在缺失点与解决方案

#### 5.1 碎片提取算法细节

**问题**：`extract_fragments()` 如何基于 `forming_bonds` 切分 TS 结构？

**解决方案**：
```python
def extract_fragments(coords, symbols, forming_bonds):
        """
        基于形成键切分 TS 结构为两个碎片
    
        算法:
        1. 构建分子邻接图 (基于共价半径判断键连)
        2. 移除 forming_bonds 中的两条边
        3. BFS/DFS 遍历获取连通分量
        4. 返回两个碎片的坐标与符号
    
        电荷/多重度确定:
        - 碎片A (oxidopyrylium): 电荷=0, 多重度=1
        - 碎片B (dienophile): 电荷=0, 多重度=1
        - 若切分后出现自由基, 需特殊处理
        """
        from scipy.sparse.csgraph import connected_components
        # ... 实现细节
```

**边界情况处理**：
| 情况 | 处理策略 |
|---|---|
| 切出 3+ 碎片 | 警告并合并最小碎片到最近邻 |
| 自由基碎片 | 多重度设为 2 (doublet) |
| 电荷不守恒 | 基于电负性分配电荷 |

#### 5.2 IRC 验证闭环

**问题**：IRC 端点验证失败后如何处理？

**解决方案**：
```python
def step3_ts_optimizer_with_irc(...):
        ts_final = berny_or_qst2_optimize(...)
    
        if VERIFY_IRC:
                irc_result = run_irc(ts_final)
        
                if not validate_irc_endpoints(irc_result, reactant, product):
                        # 策略1: 尝试从IRC端点重新优化
                        ts_from_irc = reoptimize_from_irc_saddle(irc_result)
            
                        # 策略2: 标记为"IRC未验证"但继续
                        if ts_from_irc is None:
                                logger.warning("IRC 验证失败，TS 可能不连接正确端点")
                                result.irc_verified = False
```

#### 5.3 错误处理与断点续传

**问题**：批量处理时某个分子失败，如何恢复？

**解决方案**：
```python
class BatchOrchestrator:
        def run_batch_with_checkpoint(self, smiles_list, checkpoint_file):
                """
                支持断点续传的批量处理
        
                机制:
                1. 每完成一个分子, 写入 checkpoint.json
                2. 重启时读取 checkpoint, 跳过已完成的
                3. 失败的分子记录到 failed.csv, 可后续重试
                """
                completed = self._load_checkpoint(checkpoint_file)
        
                for smiles in smiles_list:
                        if smiles in completed:
                                continue
            
                        try:
                                result = self.run_pipeline(smiles)
                                self._save_checkpoint(smiles, result, checkpoint_file)
                        except Exception as e:
                                self._log_failure(smiles, e, "failed.csv")
                                continue  # 不中断整个批次
```

#### 5.4 配置管理完善

**问题**：`defaults.yaml` 的具体内容未定义。

**建议内容**：
```yaml
# config/defaults.yaml
global:
    log_level: INFO
    scratch_dir: /tmp/rph_scratch
    keep_intermediates: false

step1:
    crest:
        gfn_level: 2
        solvent: acetone
        energy_window: 6.0  # kcal/mol
        threads: 16
    dft:
        method: B3LYP
        basis: def2-SVP
        dispersion: GD3BJ
        verify_no_imaginary: true

step2:
    ts_distance: 2.2  # Å
    reactant_distance: 3.5  # Å
    xtb_gfn_level: 2

step3:
    method: B3LYP
    basis: def2-SVP
    dispersion: GD3BJ
    nprocshared: 16
    mem: 32GB
    verify_irc: false
    irc_max_points: 50
    irc_step_size: 10

step4:
    enable_nbo: false
    enable_distortion: true
    fragment_sp_method: B3LYP/def2-SVP
```

#### 5.5 架构细节统一

**问题**：文档与代码不一致。

| 项目 | 文档描述 | 实际代码 | 统一方案 |
|---|---|---|---|
| `BondStretcher` 位置 | `utils/` | `step2_retro/` | **保留在 `step2_retro/`**（业务相关） |
| `forming_bonds` 传递 | 未说明 | S2返回→Orchestrator存储→S4使用 | **在此文档化** |
| `CRESTInterface` | 未提及 | 已实现于 `qc_interface.py` | **补充到文档** |

**数据传递机制图**：
```
S2.run() → (ts_guess, reactant, forming_bonds)
                            │                    │
                            ▼                    ▼
                Orchestrator         PipelineResult.forming_bonds
                            │                    │
                            ▼                    ▼
                S4.run(..., forming_bonds=result.forming_bonds)
                            │
                            ▼
                FeatureMiner._extract_features(..., forming_bonds)
                            │
                            ▼
                features["r1_ts"], features["r2_ts"], features["asynchronicity"]
```

#### 5.6 NBO 分析详细规划

**问题**：NBO route 和解析逻辑未详细说明。

**解决方案**：
```python
# Gaussian NBO 关键词
NBO_ROUTE = "# B3LYP/def2-SVP Pop=NBO"

# 需要提取的 NBO 特征
NBO_FEATURES = {
        "nbo_charge_c1": "形成键 C1 的 NBO 电荷",
        "nbo_charge_c2": "形成键 C2 的 NBO 电荷",
        "wiberg_bond_order": "形成键的 Wiberg 键级",
        "donor_acceptor_e2": "主要 donor→acceptor E(2) 能量 (kcal/mol)",
        "homo_nbo": "基于 NBO 的 HOMO 能量 (eV)",
        "lumo_nbo": "基于 NBO 的 LUMO 能量 (eV)"
}

# 解析逻辑 (正则模式)
NBO_PATTERNS = {
        "charge": r"^\s*(\d+)\.\s+(\w+)\s+([\-\d\.]+)\s*$",
        "wiberg": r"Wiberg bond index.*?(\d+)\s+(\d+)\s+([\d\.]+)",
        "e2": r"E\(2\)=\s*([\d\.]+)\s+kcal/mol"
}
```

---

### 六、更新后的实施路线（修订版）

| 阶段 | 任务 | 优先级 | 状态 |
|:---:|---|:---:|:---:|
| **基础设施** | 迁移 parallel_processor/reaxys_loader | P0 | ✅ 已完成 |
| **基础设施** | 实现断点续传机制 | P1 | ⬜ 待开发 |
| **S1 精修** | ProductAnchor DFT 优化分支 | P1 | ⬜ 待开发 |
| **S2 增强** | forming_bonds 传递机制文档化 | P0 | ✅ 已完成 |
| **S3 闭环** | IRC 验证失败回退策略 | P2 | ⬜ 待开发 |
| **S4 核心** | extract_fragments() 图切分算法 | P0 | ⬜ 待开发 |
| **S4 核心** | 碎片 SP/OPT 能量计算 | P0 | ⬜ 待开发 |
| **S4 扩展** | NBO 分析集成 | P2 | ⬜ 待开发 |
| **配置** | defaults.yaml 完善 | P1 | ⬜ 待开发 |

---

**结论**：本补充将最大程度完善体系化迁移方案的工程与物理细节，确保 PROMOTE.md 设计与实现无缝衔接。

---

### 七、设计复审：效率与合理性优化（2026-01-10）

本节对 PROMOTE.md 现有设计进行批判性审查，识别潜在的效率瓶颈与不合理之处，并给出可执行的优化建议。

1) 重复 DFT 计算导致效率低下
- 问题：S1、S3、S4 多次独立做 DFT（OPT/Freq/SP），碎片 SP 未复用父体系的 checkpoint（`.chk`），SCF 重复开销大。
- 建议：对碎片 SP 使用父体系或 TS 的 checkpoint（`Guess=Read` / `%oldchk`），并在文档中给出示例命令/关键词以说明如何生成 fragment.chk 并复用轨道猜测。

2) xTB → DFT 跳跃风险
- 问题：xTB 受限优化与 B3LYP/def2-SVP 势能面差异可能导致 Berny 收敛困难。
- 建议：在 S2 和 S3 之间增加轻量级 DFT 预优化（例如 B3LYP/6-31G* 受限优化 5–10 步）作为中间层；在 `defaults.yaml` 增加 `pre_dft_refine` 配置项。

3) QST2 串行救援不够高效
- 问题：把 QST2 仅作为 Berny 失败后的救援会在失败情况下浪费大量时间。
- 建议：在资源充足环境下采用并行竞赛（Berny 与 QST2 同时提交，先成功者为准），文档中说明适用条件与示例实现思路。

4) 畸变能公式需修正为多碎片形式
- 问题：当前文档只用单一 `E_dist`，但 [5+2] 反应牵涉至少两个片段。
- 建议：将公式修正为
    E_{dist,total} = E_{dist,A} + E_{dist,B}
    E_{int} = ΔE^‡ - E_{dist,total}
并在 Feature_Miner 示例代码中给出双片段计算伪代码。

5) 并行粒度应自适应以避免资源争抢
- 问题：现有分子級并行在处理少量大分子时会造成内存/IO 争抢。
- 建议：实现两级并行策略（分子级 vs 任务级自适应）并在文档中提供策略选择准则与示例。

6) 增加中間結果缓存與斷點續傳規範
- 問题：缺少标准化的状态文件导致失败重跑昂贵。
- 建议：定义工作目录布局与 `.state` 文件（S1_DONE/S2_DONE/...），在批量运行中写入 checkpoint.json/failed.csv，文档列出恢复逻辑示例。

7) IRC 验证应为条件触发而非全量执行
- 建议：提供判断函数（基于虚频幅值、TS 与初猜 RMSD、形成键距离异常）决定是否运行 IRC，从而节省不必要的昂贵 IRC 计算。

8) CREST 参数针对刚性体系优化
- 问题：`energy_window=6.0 kcal/mol` 对刚性桥环过宽，产生大量冗余构象。
- 建议：对刚性产物默认 `energy_window=3.0 kcal/mol` 并启用 `quick` 模式，或在 `defaults.yaml` 提供针对性 profile（rigid vs flexible）。

---

优先级建议（快速清单）：
- P0: 修正畸变能公式并在 Feature_Miner 中实现双片段计算（物理正确性，低复杂度）
- P1: 实现 DFT checkpoint 复用、CREST 参数优化、中间结果缓存（效率收益显著）
- P2: 引入中间层 DFT 预优化、并行竞赛 QST2/Berny、两级并行策略（需要集群资源与更多验证）

实施这些优化可在不改变核心串行架构前提下，显著提高收敛率、降低总计算时间、并提升物理正确性。建议先完成 P0/P1 项，再按需推进 P2 项。

---

### 八、检查测试流程设计（2026-01-10）

完整测试计划详见 **[TEST_PLAN.md](TEST_PLAN.md)**。

#### 8.1 测试架构概览

```
测试金字塔

          ┌─────────────┐
          │ E2E 验收    │  ← 5% (真实数据验证)
          └──────┬──────┘
          ┌──────┴──────┐
          │  集成测试   │  ← 25% (S1→S2→S3→S4)
          └──────┬──────┘
          ┌──────┴──────┐
          │  单元测试   │  ← 70% (每个函数/类)
          └─────────────┘

目标覆盖率: > 85%
```

#### 8.2 单元测试用例汇总

| 模块 | 测试文件 | P0 用例数 | P1 用例数 |
|:---:|---|:---:|:---:|
| S1 ProductAnchor | `test_s1_anchor.py` | 5 | 3 |
| S2 RetroScanner | `test_s2_retro.py` | 8 | 2 |
| S3 TSOptimizer | `test_s3_optimizer.py` | 5 | 4 |
| S4 FeatureMiner | `test_s4_features.py` | 6 | 2 |

关键验证点：
- S1: `verify_no_imaginary()` 确保产物无虚频
- S2: `ts_distance == 2.2 Å`, `break_distance == 3.5 Å`
- S3: `validate_imaginary_freq()` 恰好 1 个虚频, `Guess=Read` checkpoint 复用
- S4: `E_dist = E_dist(substrate) + E_dist(reagent)` 双片段畸变能

#### 8.3 集成测试用例

| 测试用例 | 验证流程 | 通过标准 |
|---|---|---|
| `test_s1_to_s2_flow` | S1 → S2 | 产物 xyz 正确传递 |
| `test_s2_to_s3_flow` | S2 → S3 | ts_guess + reactant 传递 |
| `test_s3_to_s4_flow` | S3 → S4 | ts_final + forming_bonds |
| `test_full_pipeline_mock` | S1→S2→S3→S4 | Mock 完整流程无异常 |
| `test_checkpoint_resume` | 断点恢复 | 从正确步骤继续 |

#### 8.4 端到端验收标准

| 验收项 | 指标 | 容差 |
|---|---|:---:|
| 活化能准确性 | ΔG‡ 与文献值比较 | ±3 kcal/mol |
| 能量守恒 | ΔE‡ = E_dist + E_int | < 0.1 kcal/mol |
| 性能目标 | 单分子完整流程 | < 2h (def2-SVP) |
| checkpoint 加速 | SCF 迭代减少 | > 30% |

#### 8.5 每日测试检查点

| Day | 测试文件 | 通过标准 |
|:---:|---|---|
| 1-2 | `test_s1_anchor.py` + `test_s2_retro.py` | 全部 P0 通过 |
| 3-5 | `test_s3_optimizer.py` | 虚频验证 + checkpoint |
| 6-7 | `test_s4_features.py` | 双片段公式正确 |
| 8 | `test_integration.py` | 完整流程通过 |
| 9 | `test_e2e_acceptance.py` | 误差 < 3 kcal/mol |
| 10 | 断点续传 + 性能 | 恢复成功, 加速 > 30% |

#### 8.6 测试执行命令

```bash
# 快速单元测试（无 QC 依赖）
pytest tests/ -m "not slow and not requires_gaussian" -v

# 完整测试 + 覆盖率
pytest tests/ -v --cov=rph_core --cov-report=html

# 验收测试（需要 Gaussian）
pytest tests/test_e2e_acceptance.py -v -m requires_gaussian
```

---

## 9. 附录: 外部计算后端模板 (Backend Templates)

### 9.1 Step 3.5: L2 高精度单点 (ORCA-RIJCOSX)
**应用场景**: 用于最终能量校正、HBDE 及产物能量确定。
```orca
! M062X def2-TZVPP def2/J RIJCOSX tightSCF 
! Grid4 FinalGrid5
%maxcore 4000
%pal nprocs 16 end
%cpcm 
   smd true 
   SMDsolvent "acetone" 
end

* xyz [CHARGE] [SPIN]
[Coordinates]
*
```

### 9.2 Step 3.5: ORCA 辅助基组匹配规则
- **普通泛函**: `def2-TZVPP` + `def2/J`
- **双杂化泛函 (如 PWPB95)**: `def2-TZVPP` + `def2-TZVPP/C` + `def2/J`

### 9.3 Step 4: Gaussian NBO 7.0 接口
**应用场景**: 提取 NBO 电荷和键级。
```
%chk=ts_final.chk
%mem=32GB
%nproc=16
# B3LYP/def2-SVP Pop=NBO7

NBO Analysis

0 1
```

### 9.4 Step 4: Gaussian NICS(1) Ghost 原子模板
**应用场景**: 计算环的正上方 1.0Å 处的屏蔽常数。
```
%chk=ts_final.chk
# B3LYP/def2-SVP NMR

NICS(1) Calculation with Ghost Atom

0 1
Bq  0.0  0.0  1.0  (Centroid-based coordinates)
...
```

0 1
Bq  0.0  0.0  1.0  (Centroid-based coordinates)
...
```

