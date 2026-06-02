# ReactionProfileHunter (RPH) 全流程研究报告

> **版本**: v2.1.1 | **代码量**: ~54k Python (212 .py 文件) | **生成日期**: 2026-05-29

---

## 目录

1. [项目概述](#1-项目概述)
2. [整体架构](#2-整体架构)
3. [S0: 机理分类 (Mechanism Classifier)](#3-s0-机理分类)
4. [S1: 产物锚定与构象搜索 (Anchor & Conformer Search)](#4-s1-产物锚定与构象搜索)
5. [S2: 逆扫描与TS初猜 (Retro Scan & TS Guess)](#5-s2-逆扫描与ts初猜)
6. [S3: 过渡态优化与验证 (TS Optimization & Validation)](#6-s3-过渡态优化与验证)
7. [S4: 特征提取 (Feature Extraction)](#7-s4-特征提取)
8. [配置体系 (Configuration System)](#8-配置体系)
9. [QC计算基础设施 (QC Computation Infrastructure)](#9-qc计算基础设施)
10. [断点续传与容错 (Checkpoint & Resilience)](#10-断点续传与容错)
11. [批处理与并行 (Batch & Parallel Execution)](#11-批处理与并行)
12. [数据流与步间合同 (Data Flow & Inter-Step Contracts)](#12-数据流与步间合同)
13. [总结与展望](#13-总结与展望)

---

## 1. 项目概述

### 1.1 项目定位

ReactionProfileHunter (RPH) 是一个**产物驱动 (product-driven)** 的自动化反应机理探索管线，专为环加成反应 ([4+3]、[4+2]、[3+2]、[5+2] 等) 的过渡态搜索、几何优化和物理有机特征提取而设计。

### 1.2 核心能力

| 能力 | 描述 |
|------|------|
| **产物驱动策略** | 从产物出发，逆向搜索反应路径 |
| **四步串行管线** | 构象搜索 → 逆扫描 → TS优化 → 特征提取 |
| **双层次计算** | xTB预优化 → B3LYP/M06-2X几何优化 → wB97X-D3BJ/wB97M-V高精度单点能 |
| **多引擎支持** | Gaussian (opt), ORCA (SP), xTB (pre-opt, scan), CREST (conformer), Multiwfn (CDFT) |
| **特征工程** | 活化能、几何参数、电子描述符、构象柔性、NBO分析等 20+ 类特征 |
| **断点续传** | 基于哈希验证的步级checkpoint/resume |
| **批处理** | ProcessPoolExecutor 多进程并行，支持数据集模式 |

### 1.3 技术栈

- **语言**: Python 3.8+
- **量子化学**: Gaussian 16, ORCA 5+, xTB 6+, CREST, Multiwfn, Shermo, ISOSTAT
- **化学信息学**: RDKit (SMILES解析、原子映射)
- **数学**: NumPy, SciPy (Kabsch对齐、Hungarian匹配)
- **数据处理**: pandas, JSON (配置文件、中间数据交换)

---

## 2. 整体架构

### 2.1 核心类关系

```
main() [CLI入口]
  └── ReactionProfileHunter(config_path)
       ├── Load config → normalize_qc_config()
       ├── SmallMoleculeCatalog (小分子库)
       ├── CheckpointManager (断点续传)
       └── Lazy Engines:
            ├── s0_engine: MechanismClassifier (可选，S0机理分类)
            ├── s1_engine: AnchorPhase v3.0 (S1产物锚定)
            ├── s2_engine: RetroScanner (S2逆扫描)
            ├── s3_engine: TSOptimizer (S3 TS优化)
            └── s4_engine: FeatureMiner (S4特征提取，来自rph_features包)
```

### 2.2 主流程 (`run_pipeline()`, orchestrator.py:1844-2804)

```
入口: product_smiles, work_dir
  │
  ├─ [初始化] CheckpointManager(work_dir) + 状态加载/重水合
  │
  ├─ [S0] _run_s0()
  │     └── MechanismClassifier.classify_from_dict(cleaner_data)
  │         → mechanism_graph.json, mechanism_summary.json
  │         → atom_map_smiles.json, dr_branch_plan.json
  │
  ├─ [S1] AnchorPhase.run(molecules={product, precursor, small_molecules})
  │     └── 产物/前体/小分子构象搜索 + DFT优化 + SP
  │         → product_min.xyz, e_sp, fchk/log/thermo
  │         → checkpoint 写入
  │
  ├─ [S1→S2桥] _build_and_save_xyz_mapping()
  │     └── SMILES原子映射 → product XYZ索引映射
  │         → atom_map_xyz.json
  │
  ├─ [S2] run_step2()
  │     ├── 形成键解析 (_resolve_forming_bonds_for_s2)
  │     │    级联回退: atom_map_xyz → mapping annotations → graph edges
  │     │    → mechanism_summary → SMARTS → order_changed
  │     ├── s2_engine.run_retro_scan()
  │     │    └── xTB relaxed scan (concerted/sequential)
  │     │        → 拓扑守护 (topology guard)
  │     │        → 势能面绘图 (scan profile plot)
  │     ├── [可选] s2_engine.run_path_search() (xTB PATH rescue)
  │     └── → ts_guess.xyz, intermediate.xyz, forming_bonds
  │         → scan_profile.json, step2_provenance.json
  │         → checkpoint 写入 (含签名验证)
  │
  ├─ [S3] run_step3()
  │     ├── s3_engine.run()
  │     │    ├── Intermediate DFT opt (intermediate_driver.py)
  │     │    ├── Berny TS opt (berny_driver.py)
  │     │    │    └── 失败→振荡检测→小步长rescue
  │     │    ├── [可选] IRC验证 (irc_driver.py)
  │     │    ├── Reactant SP计算
  │     │    └── SP矩阵组装 (SPMatrixReport)
  │     └── → ts_final.xyz, sp_matrix_report, fchk/log
  │         → checkpoint 写入 (含输入哈希验证)
  │
  ├─ [S3.5] resolve_forming_bonds()  (forming_bonds_resolver.py)
  │     └── 从TS几何差异解析形成键索引 (S2→S4权威源)
  │         → mechanism_meta.json
  │
  ├─ [S4] run_step4()
  │     ├── _resolve_s1_artifacts() → S1构件路径
  │     ├── s4_engine.run()
  │     │    └── 插件管线: 14+ extractors 并行执行
  │     └── → features_raw.csv, features_mlr.csv, feature_meta.json
  │         → checkpoint 写入
  │
  └─ [完成] result.success = True
```

### 2.3 PipelineResult 数据结构

| 字段 | 类型 | 来源 | 流向 |
|------|------|------|------|
| `success` | `bool` | — | 全局状态 |
| `product_smiles` | `str` | 输入 | 日志 |
| `work_dir` | `Path` | 输入 | 所有步骤 |
| `product_xyz` | `Path` | S1→S2→S3→S4 | S2/S3/S4输入 |
| `e_product_l2` | `float` | S1→S3 | S3能量参考 |
| `product_fchk/log/qm_output` | `Path` | S1→S4 | S4特征提取 |
| `product_thermo` | `Path` | S1→S3→S4 | 热力学参考 |
| `product_checkpoint` | `Path` | S1→S3 | 轨道复用 |
| `ts_guess_xyz` | `Path` | S2→S3 | S3输入 |
| `intermediate_xyz` | `Path` | S2→S3→S4 | S3/S4输入 |
| `ts_final_xyz` | `Path` | S3→S4 | S4输入 |
| `forming_bonds` | `Tuple[Tuple[int,int],...]` | S2 or S3.5→S4 | 片段切割 |
| `sp_matrix_report` | `SPMatrixReport` | S3→S4 | 能量特征 |
| `ts_fchk/log/qm_output` | `Path` | S3→S4 | S4特征提取 |
| `intermediate_fchk/log/qm_output` | `Path` | S3→S4 | S4特征提取 |
| `s3_intermediate_xyz` | `Path` | S3→S4 | 中间体几何 |
| `s3_intermediate_l2_energy` | `float` | S3→S4 | 中间体能量 |
| `features_csv` | `Path` | S4 | 最终输出 |
| `error_step/error_message` | `str` | 异常 | 错误追踪 |

### 2.4 目录布局

```
Output/rx_001/
├── pipeline.state              # Checkpoint状态文件
├── S0_Mechanism/
│   ├── mechanism_graph.json    # 反应机理DAG图
│   ├── mechanism_summary.json  # 机理摘要 (反应类型/形成键)
│   ├── atom_map_smiles.json    # SMILES原子映射
│   ├── dr_branch_plan.json     # DR分支计划
│   └── mechanism_graph.png     # 可视化
├── S1_ConfGeneration/
│   ├── product/
│   │   ├── crest/              # CREST构象搜索
│   │   ├── xtb/                # 两阶段xTB (stage1_gfn0/stage2_gfn2)
│   │   ├── cluster/            # ISOSTAT聚类结果
│   │   ├── finalDFT/           # DFT OPT+SP
│   │   └── product_min.xyz     # 最终产物几何
│   ├── precursor/              # 前体构象 (可选)
│   ├── <small_molecule>/       # 小分子参考态
│   ├── atom_map_xyz.json       # XYZ原子映射
│   ├── shermo_summary.json     # 热力学汇总
│   └── conformer_energies.json # 构象能量
├── S2_Retro/
│   ├── ts_guess.xyz            # TS初猜几何
│   ├── intermediate.xyz        # 中间体 (新命名)
│   ├── reactant_complex.xyz    # 中间体 (旧别名)
│   ├── scan_profile.json       # 扫描势能面
│   ├── scan_profile_plot.png   # 势能面图
│   ├── step2_provenance.json   # S2溯源文件
│   └── pes_adapter_fallback/   # PES适配器回退 (可选)
├── S3_TransitionAnalysis/
│   ├── ts_final.xyz            # 优化后的TS
│   ├── ts_opt/                 # TS优化 (berny/)
│   ├── S3_intermediate_opt/    # 中间体DFT优化
│   ├── L2_SP/                  # 高精度单点能
│   ├── sp_matrix_metadata.json # SP矩阵元数据
│   ├── step3_provenance.json   # S3溯源文件
│   └── mechanism_meta.json     # 形成键元数据
└── S4_Data/
    ├── features_raw.csv        # 原始特征表
    ├── features_mlr.csv        # ML就绪特征表
    ├── feature_meta.json       # 特征元数据
    └── qc_nbo.37               # NBO文件 (可选)
```

---

## 3. S0: 机理分类

### 3.1 模块位置与入口

- **目录**: `rph_core/steps/mechanism_classifier/`
- **主类**: `MechanismClassifier` (`classifier.py:40`)
- **入口方法**: `classify_from_dict(cleaner_row)` → `MechanismGraph`

### 3.2 核心组件

| 组件 | 文件 | 功能 |
|------|------|------|
| `MechanismClassifier` | `classifier.py:40` | 主入口，整合Clean适配器和图构建器 |
| `CleanAdapter` | `clean_adapter.py` | 解析Clean程序输出格式（SMILES + 键变化标注） |
| `GraphBuilder` | `graph_builder.py` | 构建MechanismGraph有向无环图 |
| `MechanismGraph` | `models.py:137` | 核心数据模型 |
| `DRCompletion` | `dr_completion.py` | DR (Diradical) 路径补全 |

### 3.3 数据模型

```python
MechanismGraph
├── reaction_id: str           # 反应标识符
├── reaction_type: str         # [4+3], [5+2], [4+2], [3+2]
├── cyclo_mode: CycloMode      # 环化模式枚举
├── topology: TopologyType      # 分子间/分子内
├── nodes: List[GraphNode]     # 反应物→中间体→产物→复合物
├── edges: List[GraphEdge]     # 边属性 (forming_bonds, breaking_bonds, ts_type)
├── pathways: List[PathwayInfo] # 反应路径
├── smiles_atom_mapping         # SMILES原子映射
└── forming_bonds_annotated     # 标注的形成键
```

### 3.4 执行策略

- **可选步骤**: S0 disabled时直接跳过，不影响后续S2
- **失败策略**: `fail_policy: degrade` — S0失败后S2使用SMARTS自动检测回退
- **低置信度处理**: 标记 `low_confidence=True`，S2将defer使用该结果
- **DR补全**: `dr_completion.enabled=true` 时自动生成diradical分支路径

### 3.5 输出物

| 输出 | 文件 | 用途 |
|------|------|------|
| 机理图 | `S0_Mechanism/mechanism_graph.json` | 完整DAG图，供后续分析 |
| 机理摘要 | `S0_Mechanism/mechanism_summary.json` | 反应类型/形成键/低置信度标记 |
| 原子映射 | `S0_Mechanism/atom_map_smiles.json` | SMILES→原子映射号 |
| DR分支计划 | `S0_Mechanism/dr_branch_plan.json` | Diradical分支枚举 |
| 可视化 | `S0_Mechanism/mechanism_graph.png` | 图形化反应机理 |

---

## 4. S1: 产物锚定与构象搜索

### 4.1 模块位置与入口

- **目录**: `rph_core/steps/anchor/`, `rph_core/steps/conformer_search/`
- **主类**: `AnchorPhase` (`anchor/handler.py`) → `ConformerEngine` (`conformer_search/engine.py:64`)
- **入口方法**: `AnchorPhase.run(molecules: dict)` → `AnchorPhaseResult`

### 4.2 构象搜索策略 — UCE v3.1

#### 协议栈 (Protocol Stack)

RPH v2.1 支持四种构象搜索协议，通过 `step1.protocol` 配置切换：

| 协议 | 搜索模式 | 构象数 | DFT计算 | 适用场景 |
|------|---------|--------|---------|---------|
| **ext** (扩展) | 两阶段 GFN0→GFN2 | 3-6 | optimize_all_candidates + shermo_boltzmann | 高精度需求 |
| **full** (完整) | 单阶段 GFN2 | 12-24 | DFT预筛→SP筛选→Boltzmann排序 | 全面构象搜索 |
| **lite** (轻量) | 单阶段 GFN2 | 4-6 | optimize_rank1 + r2SCAN-3c final | 日常计算 (默认) |
| **zero** (最小) | 单阶段 GFN2/skip | 3-4 | optimize_rank1 + r2SCAN-3c | 快速筛选 |

#### 两阶段流程 (ext协议)

```
Stage 1: GFN0 粗筛
  CREST GFN0 (能量窗口 10.0 kcal)
    → ISOSTAT 聚类 (RMSD 0.5 Å, 能量 1.0 kcal)
    → cluster.xyz (候选构象簇)

Stage 2: GFN2 精细优化
  从 cluster.xyz 出发
    → CREST GFN2 (能量窗口 3.0 kcal)
    → ISOSTAT 聚类 (RMSD 1.0 Å, 能量 1.0 kcal)
    → cluster.xyz (最终构象集合)

Stage 3: DFT OPT-SP 耦合
  对每个构象候选
    → DFT 几何优化 (B3LYP/M06-2X/def2-SVP)
    → wB97X-D3BJ/wB97M-V/def2-TZVPP 高精度单点能
    → Shermo 热力学校正 (298.15K, 1 atm)
    → Boltzmann 加权排序
```

### 4.3 分子锚定 (Molecular Anchoring)

`AnchorPhase.run()` 接收分子字典，键为分子角色，值为SMILES：

```python
molecules = {
    "product": "C=C(C)C(=O)O",       # 产物 (必需)
    "precursor": "C=C(C)C(=O)OC",    # 前体 (可选)
    "AcOH": "CC(=O)O",              # 小分子离去基团
    "DMDO": "CC1(C)OO1",            # 氧化剂
}
```

每个分子独立执行完整构象搜索流程，互不干扰（分子自治架构）。

### 4.4 关键输出

| 输出 | 路径 | 流向 |
|------|------|------|
| 产物最优结构 | `S1_ConfGeneration/product_min.xyz` | S2, S3, S4 |
| SP能量 | `e_product_l2` (float) | S3能量参考 |
| fchk文件 | `product/finalDFT/*.fchk` | S4波函数分析 |
| log文件 | `product/finalDFT/*.log` | S4频率分析 |
| 热力学数据 | `conformer_thermo.csv` | S4活化能特征 |
| 构象能量 | `conformer_energies.json` | S4 Boltzmann特征 |
| Shermo汇总 | `shermo_summary.json` | S4 Gibbs自由能 |

---

## 5. S2: 逆扫描与TS初猜

### 5.1 模块位置与入口

- **目录**: `rph_core/steps/step2_retro/`
- **主类**: `RetroScanner` (`retro_scanner.py:47`)
- **入口方法**: `run_retro_scan(product_xyz, output_dir, forming_bonds, scan_config)`

### 5.2 核心组件

| 组件 | 文件 | 功能 |
|------|------|------|
| `RetroScanner` | `retro_scanner.py:47` | 主引擎，协调整体流程 |
| `BondStretcher` | `bond_stretcher.py` | 键拉伸几何操作 |
| `KinematicStretcher` | `kinematic_stretcher.py` | 运动学拉伸（刚体片段） |
| `GeometryGuard` | `geometry_guard.py` | 拓扑守护，检测错误成键 |
| `ScanPolicies` | `scan_policies.py` | 扫描策略 (policy_a/b/c/concerted) |
| `SMARTSMatcher` | `smarts_matcher.py` | SMARTS模式匹配回退 |

### 5.3 形成键解析链 (Forming Bonds Resolution)

S2的核心挑战是确定哪个原子对将形成新键。RPH实现了7级级联回退链 (`_resolve_forming_bonds_for_s2`)：

```
优先级 (高→低):

1️⃣ S1 atom_map_xyz full notation
   └── 从 atom_map_xyz.json 的 forming_bonds_full_notation 提取

2️⃣ S0/S1 atom-map annotations
   └── 从 SMILES映射 + 标注的形成键 → product XYZ索引

3️⃣ S0 mechanism_graph primary edges + S1 atom_map_xyz
   └── 从机理图primary边 + XYZ映射联合推导

4️⃣ S0 mechanism_summary.json legacy fallback
   └── 从 mechanism_summary 的 forming_bonds 字段

5️⃣ SMARTS auto-detection
   └── SMARTSMatcher 模板注册: [5+2], [4+3], [4+2], [3+2]

6️⃣ cleaner data forming_bonds
   └── 从cleaner数据中读取

7️⃣ order_changed derivation
   └── 从 core_bond_changes 的 AROMATIC→SINGLE 等键级变化推导
```

每一级解析后都经过 `_validate_forming_bonds_against_product()` 验证（距离检查 1.2-3.5 Å），失败则降级到下一级。

### 5.4 xTB Relaxed Scan 流程

```
输入: product.xyz + forming_bonds + scan_config

1. BondStretcher: 将形成键拉伸到 scan_start_distance
   └── 可选: clean break (去除小分子离去基团)

2. xTB GFN2 Relaxed Scan (20 steps, 3.5→1.8 Å)
   ├── 扫描模式: concerted (同时扫描) / sequential (依次扫描)
   ├── 扫描策略: policy_a/b/c (不同约束力常数)
   ├── 约束力常数: scan_force_constant (0.5 默认)
   └── 溶剂: acetone (PCM隐式模型)

3. 拓扑守护 (Topology Guard)
   ├── 检测扫描轨迹中是否形成错误化学键
   ├── 自动添加 keep-away 约束
   └── 失败时重试 (topology_retry_once, 增强力常数)

4. 势能面分析
   ├── find_ts_and_dipole_guess(): 定位能量峰值 (TS)
   └── 边界检测: reject_boundary_maximum → 扩展扫描范围重试

5. 输出
   ├── ts_guess.xyz: 能量最高点几何 (TS初猜)
   ├── intermediate.xyz: 键解离态几何 (中间体)
   └── scan_profile.json: 完整扫描数据 + 势能面图
```

### 5.5 xTB Path Search 救援 (S2.2)

当 retro_scan 结果标记为 DEGRADED/FAILED 时，自动触发 xTB PATH 路径搜索：

```
xTB PATH: start.xyz (intermediate) --path→ end.xyz (product)
  ├── npoint: 25 (初始路径点)
  ├── anopt: 10 (精细优化步)
  ├── kpush: 0.003 (推力)
  ├── kpull: -0.015 (拉力)
  └── ppull: 0.05 (优化器拉力)
```

### 5.6 关键输出

| 输出 | 文件 | 流向 |
|------|------|------|
| TS初猜 | `S2_Retro/ts_guess.xyz` | S3输入 |
| 中间体 | `S2_Retro/intermediate.xyz` | S3/S4输入 |
| 形成键 | `(Tuple[int,int], ...)` | S3→S4形成键解析 |
| 扫描数据 | `S2_Retro/scan_profile.json` | S4特征提取 |
| 溯源文件 | `S2_Retro/step2_provenance.json` | Checkpoint验证 |
| 势能面图 | `S2_Retro/scan_profile_plot.png` | 可视化 |

---

## 6. S3: 过渡态优化与验证

### 6.1 模块位置与入口

- **目录**: `rph_core/steps/step3_opt/`
- **主类**: `TSOptimizer` (`ts_optimizer.py:220`)
- **入口方法**: `TSOptimizer.run(ts_guess, intermediate, product, output_dir, ...)`

### 6.2 核心组件

| 组件 | 文件 | 功能 |
|------|------|------|
| `TSOptimizer` | `ts_optimizer.py:220` | 主协调器，编排所有子驱动 |
| `BernyDriver` | `berny_driver.py` | 核心Berny TS优化 |
| `IntermediateDriver` | `intermediate_driver.py` | 中间体DFT优化 |
| `IRCDriver` | `irc_driver.py` | IRC路径跟踪 |
| `QST2Rescue` | `qst2_rescue.py` | QST2救援策略 (默认禁用) |
| `Validator` | `validator.py` | 收敛+虚频一致性验证 |
| `ArtifactResolver` | `artifact_resolver.py` | 目录布局兼容和产物解析 |

### 6.3 TS优化流程

```
输入: ts_guess.xyz + intermediate.xyz + product.xyz

阶段1: Intermediate DFT优化
  └── IntermediateDriver.optimize()
      ├── xTB预优化 (避免原子重叠)
      ├── DFT OPT (M06-2X/def2-SVP)
      └── → intermediate_opt.xyz

阶段2: TS优化 (Berny primary)
  └── BernyDriver.optimize()
      ├── Gaussian关键词: Opt=(TS, CalcFC, NoEigenTest) Freq
      ├── Hessian初始: calcfc (每10步重算)
      ├── 振荡检测: OscillationDetector (窗口10步，能量容忍度0.0001 Ha)
      └── 失败分型:
          ├── 振荡/步数超限 → 小步长Berny rescue (MaxCycles=80, CalcAll)
          │                   从最后几何继续
          ├── 0虚频 → 回原始TS guess重试
          └── 多虚频 → 从失败几何继续

阶段3: [可选] QST2救援 (默认禁用)
  └── QST2Rescue.run()
      关键词: Opt=(QST2, CalcFC) Freq
      起点: intermediate → TS guess
      终点: product

阶段4: [可选] IRC验证
  └── IRCDriver.run()
      关键词: IRC=(CalcFC, MaxPoints=50, StepSize=10)
      验证TS连接正确的反应物和产物

阶段5: 中间体高精度SP
  └── wB97M-V/def2-TZVPP 单点能 (ORCA)

阶段6: SP矩阵组装 (SPMatrixReport)
  └── 综合TS/Reactant/Product/Fragments能量
      → ΔG‡ = G_TS - G_reactant
      → ΔG_rxn = G_product - G_reactant
```

### 6.4 TS救援策略

```
失败分型 → 救援策略映射:

┌─────────────────────┬─────────────────────────────────────┐
│ 失败类型            │ 救援策略                            │
├─────────────────────┼─────────────────────────────────────┤
│ 振荡/步数超限       │ 从最后中断结构出发                  │
│ (>60步触发)         │ 小步长 Berny rescue (MaxCycles=80) │
│                     │ CalcAll (每步重算Hessian)           │
├─────────────────────┼─────────────────────────────────────┤
│ 0虚频 (非TS)       │ 回原始 S2 TS guess 重新开始         │
├─────────────────────┼─────────────────────────────────────┤
│ 多虚频 (>1)        │ 从失败几何继续优化                   │
├─────────────────────┼─────────────────────────────────────┤
│ QST2失败            │ 标记步骤失败，不静默继续            │
└─────────────────────┴─────────────────────────────────────┘
```

### 6.5 SPMatrixReport 数据结构

```python
@dataclass
class SPMatrixReport:
    e_ts: float                    # TS能量 (L2) [Hartree]
    e_reactant: float              # 反应物/中间体能量 (L2)
    e_product: float               # 产物能量 (L2)
    g_ts: float                    # TS Gibbs自由能 [kcal/mol]
    g_reactant: float              # 反应物 Gibbs自由能
    g_product: float               # 产物 Gibbs自由能
    g_ts_source: str               # Gibbs来源 (Shermo/orca/...)
    g_reactant_source: str

    # 派生属性
    @property g_intermediate       # S0语义: g_intermediate ≡ g_reactant

    # 方法
    get_activation_energy() → ΔG‡  # G_TS - G_reactant
    get_reaction_energy() → ΔG_rxn  # G_product - G_reactant
```

### 6.6 关键输出

| 输出 | 路径 | 流向 |
|------|------|------|
| 优化TS | `S3_TransitionAnalysis/ts_final.xyz` | S4 |
| SP矩阵 | `sp_matrix_report` (SPMatrixReport) | S4能量特征 |
| TS fchk | `ts_opt/berny/*.fchk` | S4波函数分析 |
| TS log | `ts_opt/berny/*.log` | S4频率分析 |
| 中间体fchk | `S3_intermediate_opt/*.fchk` | S4中间体分析 |
| 中间体log | `S3_intermediate_opt/*.log` | S4频率分析 |
| metadata | `sp_matrix_metadata.json` | Checkpoint复用 |
| provenance | `step3_provenance.json` | Checkpoint验证 |
| mechanism_meta | `mechanism_meta.json` | 形成键权威源 |

---

## 7. S4: 特征提取

### 7.1 模块位置与入口

- **包**: `rph_features` (独立包，通过 `pip install -e ./rph_features` 安装)
- **主类**: `FeatureMiner` (`rph_features/rph_features/feature_miner.py:29`)
- **入口方法**: `FeatureMiner.run(ts_final, intermediate, product, output_dir, ...)`

### 7.2 插件架构

S4采用**插件管线**设计，基于 `BaseExtractor` 抽象基类：

```python
# 基类接口 (base.py)
class BaseExtractor(ABC):
    @abstractmethod
    def extract(self, context: FeatureContext) -> FeatureResult:
        """从上下文中提取特征，返回 FeatureResult"""

# 自动发现: 子类通过 @register_extractor 装饰器注册
```

### 7.3 提取器清单 (25个文件)

| 提取器 | 文件 | 提取特征 | 前缀 |
|--------|------|---------|------|
| **热力学** | `thermo.py` | dE_activation, dE_reaction, Gibbs能 | `thermo.` |
| **几何** | `geometry.py` | 形成键距离, 不对称度, 紧密接触 | `geom.` |
| **TS质量** | `ts_quality.py` | 虚频数, 虚频值, 收敛状态 | `ts.` |
| **QC检查** | `qc_checks.py` | 计算质量标记, 样品权重 | `qc.` |
| **S1活化** | `step1_activation.py` | 构象活化能, Boltzmann权重, 构象熵 | `s1_` |
| **S2环化** | `step2_cyclization.py` | FMO能量, CDFT指标, GEDT | `s2_` |
| **Multiwfn** | `multiwfn_features.py` | Fukui函数, Dual Descriptor | `mwfn.` |
| **中间体Multiwfn** | `intermediate_multiwfn.py` | 中间体的Multiwfn指标 | `mwfn_int.` |
| **构象布居** | `conformational_population.py` | Boltzmann布居, 有效构象数 | `conf.` |
| **构象预组织** | `conformational_preorganization.py` | 前体预组织能 | `preorg.` |
| **前体几何** | `precursor_geometry.py` | 离去基团几何参数 | `prec_geom.` |
| **前体预组织** | `precursor_preorganization.py` | 前体构象预组织 | `prec_preorg.` |
| **扭转柔性** | `torsion_flexibility.py` | 扭转角柔性分析 | `tors.` |
| **TS能量** | `ts_energy.py` | TS能量分解 | `ts_e.` |
| **相互作用** | `interaction_analysis.py` | 分子间相互作用 | `int.` |
| **NBO E2** | `nbo_e2.py` | NBO二阶稳定化能 | `nbo.` |
| **FMO/CDFT** | `fmo_cdft_dipolar.py` | 前线轨道, 偶极描述符 | `cdft.` |
| **NICS** | `nics.py` | 核独立化学位移 | `nics.` |
| **中间体特征** | `intermediate_features.py` | 中间体几何/电子特征 | `int_` |
| **TS优化元数据** | `ts_optimization_metadata.py` | TS优化过程元数据 | `ts_meta.` |
| **ASM富集** | `asm_enrichment.py` | 活化应变模型 | `asm.` |

### 7.4 特征提取上下文 (FeatureContext)

S4通过 `FeatureContext` 对象传递所有可用数据给每个提取器：

```
FeatureContext:
├── ts_final.xyz         (TS最终结构)
├── intermediate.xyz     (中间体结构)
├── product.xyz          (产物结构)
├── forming_bonds        (形成键索引)
├── sp_matrix_report     (ΔG‡, ΔG_rxn)
├── ts_fchk              (TS波函数)
├── intermediate_fchk    (中间体波函数)
├── product_fchk         (产物波函数)
├── ts_log/out           (TS计算日志)
├── intermediate_log/out (中间体计算日志)
├── product_log/out      (产物计算日志)
├── S1 artifacts:
│   ├── shermo_summary.json
│   ├── conformer_energies.json
│   ├── precursor_xyz
│   ├── atom_map_xyz.json
│   └── small_molecule_gibbs
└── S3 intermediate opt data (V7.1)
```

### 7.5 输出格式

| 输出 | 内容 | 用途 |
|------|------|------|
| `features_raw.csv` | 所有可用特征 (列数随提取器变化) | 完整特征存档 |
| `features_mlr.csv` | 选定列的ML就绪特征表 | ML模型输入 |
| `feature_meta.json` | 特征schema、版本、溯源、状态 | 数据溯源追踪 |

### 7.6 降级处理

当QC产物缺失时（如fchk未生成），提取器自动降级：
- **COMPLETE**: 所有特征成功提取
- **DEGRADED**: 部分特征缺失，填充NaN + warning日志
- **FAILED**: 提取器完全失败，不影响其他提取器

---

## 8. 配置体系

### 8.1 单源真理 — `config/defaults.yaml`

RPH使用单一YAML配置文件（774行），覆盖所有运行参数。禁止分叉配置文件。

### 8.2 配置结构速览

```yaml
executables:        # 量子化学软件路径 (Gaussian/ORCA/xTB/CREST/ISOSTAT/Shermo/Multiwfn)
resources:          # CPU/内存资源 (48GB, 16核, maxcore 2600MB/核)
theory:             # 理论水平
  preoptimization:  # xTB预优化 (GFN2, overlap_threshold=1.0Å)
  optimization:     # 几何优化 (M06-2X/def2-SVP/Gaussian)
  single_point:     # 高精度SP (wB97M-V/def2-TZVPP/ORCA)
optimization_control: # OPT控制 (超时/Hessian/步长/收敛)
step1:              # S1协议栈配置 (protocol, protocol_stack, conformer_search)
step2:              # S2扫描配置 (scan_policy, scan参数, topology_guard, path_search)
step3:              # S3 TS优化配置 (rescue_policy, fragments, intermediate_opt)
step4:              # S4特征提取配置 (enabled_plugins, mlr.columns, multiwfn)
run:                # 运行模式 (dataset/batch/single, resume, output_root)
s0:                 # S0机理分类配置 (enabled, fail_policy, dr_completion)
reaction_profiles:  # 反应类型预置 (扫描起止距离、步数、模式)
reaction_reference_terms: # 小分子参考态计量关系
scheduler:          # V3.0调度架构配置
intra_reaction_parallel: # 单反应内部并行
checkpoint:         # Checkpoint/OldChk复用控制
```

### 8.3 反应Profile驱动

不同反应类型的S2参数通过 `reaction_profiles` 预置：

```yaml
reaction_profiles:
  "[4+3]_default":         # 4+3环加成
    forming_bond_count: 2
    scan:
      scan_start_distance: 3.5    # Å
      scan_end_distance: 1.8
      scan_steps: 20
      scan_mode: concerted
      scan_force_constant: 0.5

  "[5+2]_default":         # 5+2环加成
    forming_bond_count: 2
    scan:
      scan_start_distance: 3.5
      scan_end_distance: 2.0
      scan_steps: 16
      scan_force_constant: 0.8
```

### 8.4 模板系统

`config/templates/` 目录包含 Gaussian 输入模板（`.gjf`/`.com`），运行时动态读取，禁止在Python中硬编码模板字符串。

---

## 9. QC计算基础设施

### 9.1 统一QC接口 — `qc_interface.py`

**核心原则**: 所有量子化学子进程调用必须通过 `rph_core/utils/qc_interface.py`，不允许在step中直接使用 `subprocess.run()`。

### 9.2 架构组件

| 组件 | 文件 | 行数 | 功能 |
|------|------|------|------|
| `QCInterfaceFactory` | `qc_interface.py` | 2036 | 统一QC引擎工厂 |
| `GaussianInterface` | `qc_interface.py` | — | Gaussian输入生成/解析/运行 |
| `XTBInterface` | `qc_interface.py` | — | xTB优化/扫描 (含scan()方法) |
| `CRESTInterface` | `qc_interface.py` | — | CREST构象搜索 |
| `ORCAInterface` | `orca_interface.py` | 1124 | ORCA输入生成/解析/运行 |
| `XTBRunner` | `xtb_runner.py` | 588 | xTB子进程包装 (含run_scan()) |
| `QCTaskRunner` | `qc_task_runner.py` | 989 | OPT-SP耦合循环 |
| `MultiwfnRunner` | `multiwfn_runner.py` | 661 | Multiwfn非交互批量运行 |
| `ShermoRunner` | `shermo_runner.py` | 312 | Shermo热化学计算 |
| `ISOSTATRunner` | `isostat_runner.py` | — | ISOSTAT聚类 |

### 9.3 TaskKind 枚举

```python
class TaskKind(Enum):
    OPTIMIZATION      # 几何优化
    SINGLE_POINT      # 单点能计算
    FREQUENCY         # 频率分析
    TS_OPTIMIZATION   # 过渡态优化 (Berny/QST2)
    IRC              # 内禀反应坐标
    NBO              # 自然键轨道分析
    SCAN             # xTB relaxed scan (NEW: v2.1)
```

### 9.4 沙箱隔离 — LinuxSandbox

```python
class LinuxSandbox(contextlib.ContextManager):
    """安全QC执行环境"""
    # 功能: 磁盘检查、路径隔离、临时目录创建、自动清理
    # 触发条件: is_toxic_path() 检测到路径含空格或 [](){}
    # 用法:
    with LinuxSandbox(work_dir) as sandbox:
        runner.run(sandbox.path, input_content, timeout)
        harvester.harvest(sandbox.path, output_dir)
```

### 9.5 路径安全

```python
# 毒性路径检测
is_path_toxic(path) → bool  # True if path contains spaces or [](){}

# 路径规范化
normalize_path(pathlib.Path) → pathlib.Path

# 所有QC调用前必须检查:
if is_toxic_path(output_dir):
    # 在sandbox中执行
```

### 9.6 可执行文件查找优先级

```
1. config['executables']['gaussian']['path'] (配置文件)
2. 环境变量 (如 ORCA_PATH)
3. 系统 PATH
4. fallback_paths (config中定义的备选路径)
```

---

## 10. 断点续传与容错

### 10.1 CheckpointManager (`checkpoint_manager.py`, 830行)

```python
class CheckpointManager:
    """步级断点续传管理器"""

    # 核心功能:
    - 状态序列化: PipelineState → pipeline.state (JSON)
    - 哈希验证: SHA256[:16] 文件哈希防止输入变化
    - 签名对比: step2/step3 配置+输入复合签名
    - 状态重水合: 从产物文件重建状态 (resume_rehydrate)
    - 文件锁: FileRunLock (过时检测, 6小时)
```

### 10.2 状态模型

```python
PipelineState:
├── product_smiles: str
├── work_dir: str
├── start_time: str (ISO 8601)
├── last_update: str
├── steps: Dict[str, StepCheckpoint]
│   ├── "step_s0": {completed, timestamp, output_files, metadata}
│   ├── "step_s1": {completed, output_files, metadata}
│   ├── "step_s2": {completed, output_files, metadata}
│   ├── "step_s3": {completed, output_files, metadata}
│   └── "step_s4": {completed, output_files, metadata}
└── config_snapshot: Dict
```

### 10.3 各步Checkpoint验证策略

| Step | 验证方法 | 验证内容 |
|------|---------|---------|
| S0 | `is_step_completed("s0")` | 所有3个输出JSON存在 |
| S1 | `is_step_completed("s1")` | product_xyz + e_product_sp 存在 |
| S2 | `is_step_completed("s2")` + 签名验证 | ts_guess + intermediate + step2_signature 匹配 |
| S3 | `is_step3_complete()` | 签名匹配 + 输入哈希匹配 (ts_guess/intermediate/product) + upstream S2签名匹配 |
| S4 | `is_step4_complete()` | features_csv + 机理打包文件存在 |

### 10.4 状态重水合 (Rehydrate)

当 `pipeline.state` 文件丢失时，RPH通过 `rehydrate_state_from_artifacts()` 自动恢复：

```
恢复优先级:
1. Provenance-first: 读取 step2_provenance.json / step3_provenance.json
   └── 精确恢复，包括签名和元数据
2. Best-effort: 扫描磁盘产物文件
   └── 恢复基本完成状态，但签名不可用
```

### 10.5 容错设计

- **S0降级**: 分类失败不阻塞后续步骤 (fail_policy: degrade)
- **S4降级**: 缺失QC产物填充NaN + warning，不抛错
- **部分步骤完成**: `mark_step_failed_partial()` 记录失败阶段，支持重新进入
- **事务性写入**: `atomic_write_json()` 防止写入中断导致状态损坏
- **并发锁**: `FileRunLock` 防止多个进程同时操作同一工作目录

---

## 11. 批处理与并行

### 11.1 运行模式

RPH支持三种运行模式，通过 `run.source` 配置：

| 模式 | 配置 | 说明 |
|------|------|------|
| `dataset` | `run.source: dataset` | CSV/JSON数据集模式 (默认) |
| `batch` | `run.source: batch` | 批量SMILES列表模式 |
| `single` | `run.source: single` | 单分子CLI模式 |

### 11.2 Dataset模式 (主要生产模式)

```
输入: data/reaxys_cleaned.csv
  ├── rx_id                    # 反应ID
  ├── product_smiles_main      # 产物SMILES
  ├── precursor_smiles         # 前体SMILES
  ├── ylide_leaving_group      # 离去基团标识
  └── ...

流程: _run_tasks() (orchestrator.py:3264)
  ├── 1. build_tasks_from_run_config()
  │      └── 读取CSV, 过滤 (filter_ids/max_tasks)
  ├── 2. 按 reaction_id 分组
  ├── 3. 对每个反应组:
  │      ├── 创建 reaction_root + reaction_manifest.json
  │      ├── run_pipeline() (S0-S4全过程)
  │      ├── ConditionThermoCalculator (条件热力学)
  │      ├── ConditionFeatureMerger (条件特征合并)
  │      └── DRAggregator (DR分支聚合)
  └── 4. 输出: 每个反应独立的工作目录
```

### 11.3 并行策略

#### 反应间并行 (`run_batch()`, ProcessPoolExecutor)

```python
# 多进程: 每个反应独立进程, 共享无状态
with ProcessPoolExecutor(max_workers=4) as executor:
    for smiles in smiles_list:
        executor.submit(run_pipeline, smiles, task_dir)
```

#### 反应内并行 (`intra_reaction_parallel`)

```yaml
intra_reaction_parallel:
  enabled: true       # 启用单反应内部并行
  total_cores: 16
  s3:
    parallel_intermediate_ts: true  # S3中间体和TS并行优化
    ts_opt_cores: 9
    intermediate_opt_cores: 7
  s1:
    molecule_parallel: false        # 不同分子并行 (product/precursor/small)
    crest_cores: 16                 # CREST独占全核
    opt_cores_per_job: 16           # DFT OPT独占全核
```

#### V3调度架构 (V3Scheduler)

```
V3Scheduler:
  1. precompute_small_molecules: true
     └── 跨反应共享小分子计算
  2. precursor_s1_once_per_reaction: true
     └── 同一反应的所有条件共享前体S1
  3. branch_product_only_s1: true
     └── DR分支仅重算产物S1
  4. condition_branch_thermo: true
     └── 条件级热力学计算
```

### 11.4 进度追踪

RPH实现了完整的任务进度追踪系统：

```python
V4_TASK_REGISTRY:  # 预定义的15+任务ID
  ├── mechanism       # S0机理分类
  ├── product_anchor  # S1产物锚定
  ├── precursor_anchor # S1前体锚定
  ├── smallmol_anchor # S1小分子锚定
  ├── retro_scan      # S2逆扫描
  ├── ts_opt          # S3 TS优化
  ├── irc_verify      # S3 IRC验证
  ├── reactant_opt    # S3反应物优化
  ├── sp_matrix       # S3 SP矩阵
  ├── thermochemistry # S3热化学
  ├── geom_features   # S4几何特征
  ├── elec_features   # S4电子特征
  ├── thermo_features # S4热力学特征
  └── nbo_features    # S4 NBO特征
```

每个任务追踪状态: `PENDING → RUNNING → COMPLETED/FAILED/SKIPPED/CACHED`

---

## 12. 数据流与步间合同

### 12.1 步间数据流图

```
                  ┌─────────────────────────────────┐
                  │         cleaner_data             │
                  │  (CSV行: SMILES, 键变化, 反应类型)│
                  └──────────┬──────────────────────┘
                             │
                    ┌────────▼────────┐
                    │       S0        │ MechanismClassifier
                    │  机理分类 (可选) │
                    └────────┬────────┘
                             │ mechanism_summary.json (reaction_type, forming_bonds)
                             │ atom_map_smiles.json (SMILES映射)
                             │
        ┌────────────────────┼────────────────────┐
        │ product_smiles     │ precursor_smiles   │ small_molecular_keys
        ▼                    ▼                    ▼
   ┌─────────┐         ┌─────────┐         ┌──────────┐
   │   S1    │         │   S1    │         │    S1    │
   │ product │         │precursor│         │  small   │
   │ anchor  │         │ anchor  │         │ molecule │
   └────┬────┘         └────┬────┘         └────┬─────┘
        │                   │                   │
        │ product_min.xyz   │                   │
        │ e_sp (float)      │                   │
        │ fchk/log/thermo   │                   │ s1_small_molecule_gibbs
        └────────┬──────────┘                   │
                 │                              │
        ┌────────▼────────┐                     │
        │   S1→S2 Bridge   │                    │
        │ atom_map_xyz.json│                    │
        └────────┬────────┘                     │
                 │                              │
        ┌────────▼────────┐                     │
        │       S2        │ RetroScanner        │
        │  逆扫描 (xTB)    │                     │
        └────────┬────────┘                     │
                 │ ts_guess.xyz                  │
                 │ intermediate.xyz              │
                 │ forming_bonds                 │
                 │ scan_profile.json             │
        ┌────────▼────────┐                     │
        │       S3        │ TSOptimizer         │
        │  TS优化 (DFT)    │                     │
        └────────┬────────┘                     │
                 │ ts_final.xyz                  │
                 │ SPMatrixReport                │
                 │ fchk/log (TS+intermediate)    │
                 │                               │
        ┌────────▼────────┐                     │
        │     S3.5        │ Forming Bonds        │
        │  形成键解析      │ Resolver            │
        └────────┬────────┘                     │
                 │ forming_bonds (权威)          │
                 │                               │
        ┌────────▼───────────────────────────────┴──┐
        │                  S4                        │
        │            FeatureMiner                    │
        │  ┌─────────── 插件管线 ───────────┐        │
        │  │ thermo → geom → ts_q → qc →   │        │
        │  │ s1_act → s2_cyc → mwfn →     │        │
        │  │ conf_pop → preorg → nbo → ... │        │
        │  └───────────────────────────────┘        │
        └────────────────────┬──────────────────────┘
                             │
                    ┌────────▼────────┐
                    │   features_raw   │
                    │   features_mlr   │
                    │   feature_meta   │
                    └─────────────────┘
```

### 12.2 步间合同数据类型

```python
# contracts.py — 形式化步间接口

@dataclass
class Step2Artifacts:
    ts_guess_xyz: Path                    # TS初猜
    intermediate_xyz: Path                # 中间体
    forming_bonds: Tuple[Tuple[int,int],...] # 形成键 (0-based)
    generation_method: str               # "retro_scan" / "xtb_path_search"
    status: str                          # "COMPLETE" / "DEGRADED" / "FAILED"
    ts_guess_confidence: str             # 置信度评估
    degraded_reasons: Tuple[str, ...]    # 降级原因
    step2_signature: Dict                # 配置+输入签名
    scan_profile_json: Path              # 扫描势能面数据

@dataclass
class Step3Artifacts:
    ts_final_xyz: Path                   # 优化TS
    sp_report: SPMatrixReport           # ΔG‡, ΔG_rxn
    ts_fchk: Optional[Path]              # TS波函数
    ts_log: Optional[Path]               # TS计算日志
    ts_qm_output: Optional[Path]         # TS QC输出
    intermediate_fchk: Optional[Path]    # 中间体波函数
    intermediate_log: Optional[Path]     # 中间体计算日志
    intermediate_qm_output: Optional[Path] # 中间体QC输出
    intermediate_xyz: Optional[Path]     # 优化后的中间体结构 (V7.1)
    intermediate_l2_energy: Optional[float] # 中间体L2能量 (V7.1)

@dataclass
class Step4Artifacts:
    features_csv: Path                   # 特征CSV路径
```

### 12.3 关键约定

| 约定 | 说明 |
|------|------|
| **路径格式** | 统一使用 `pathlib.Path`，禁止字符串拼接 |
| **索引基准** | 内部使用0-based索引，写入Gaussian约束时转换为1-based |
| **形成键来源** | S2输出为初步形成键，S3.5 `resolve_forming_bonds()` 为权威源 |
| **命名统一** | v6.3+: `intermediate.xyz` (新) / `reactant_complex.xyz` (旧别名) |
| **降级处理** | 缺失QC产物时填充NaN而非抛错 |
| **QC调用** | 全部经过 `qc_interface.py` 沙箱执行 |
| **日志** | `logging.getLogger(__name__)`，禁止 `print()` |
| **导入** | 绝对导入 `from rph_core.utils...`，禁止多级相对导入 |
| **配置** | `config/defaults.yaml` 为单一配置源，禁止分叉 |
| **Checkpoint** | 所有步骤写入 `pipeline.state`，支持哈希验证恢复 |

---

## 13. 总结与展望

### 13.1 RPH核心优势

1. **产物驱动策略**: 从产物逆向搜索，天然适合环加成反应的过渡态定位
2. **双层次计算**: xTB快速扫描 + DFT精确优化 + wB97M-V高精度SP，平衡速度与精度
3. **健壮的容错机制**: 7级形成键解析回退、拓扑守护、TS救援策略、断点续传
4. **丰富的特征工程**: 25类提取器、20+维特征、ML就绪输出格式
5. **多引擎集成**: Gaussian/ORCA/xTB/CREST/Multiwfn无缝协作
6. **可扩展架构**: 插件化提取器、协议栈配置、反应Profile预置

### 13.2 适用反应类型

| 类型 | 支持状态 | 特点 |
|------|---------|------|
| [4+3] 环加成 | ✅ 成熟 | 主反应类型，profile优化完善 |
| [5+2] 环加成 | ✅ 成熟 | legacy backward scan |
| [4+2] Diels-Alder | ✅ 成熟 | 经典反应 |
| [3+2] 1,3-偶极环加成 | ✅ 成熟 | SMARTS模板支持 |

### 13.3 技术指标

| 指标 | 数值 |
|------|------|
| 代码总行数 | ~54,000 行 Python |
| 核心模块数 | 212 .py 文件 |
| 测试文件数 | 53 pytest 文件 |
| 配置文件 | 774行 YAML |
| 支持QC引擎 | 6个 (Gaussian, ORCA, xTB, CREST, Multiwfn, Shermo) |
| S4提取器 | 25个插件文件 |
| Checkpoint/Resume | 步级哈希验证 |

### 13.4 关键文件索引

| 模块 | 文件 | 核心内容 |
|------|------|---------|
| **总指挥** | `rph_core/orchestrator.py` | ReactionProfileHunter类, run_pipeline(), run_batch(), main() CLI |
| **步间协调** | `rph_core/steps/runners.py` | run_step2(), run_step3(), run_step4() |
| **步间合同** | `rph_core/steps/contracts.py` | Step2Artifacts, Step3Artifacts, Step4Artifacts |
| **S0机理** | `rph_core/steps/mechanism_classifier/` | 12个文件, MechanismGraph DAG模型 |
| **S1构象** | `rph_core/steps/conformer_search/engine.py` | ConformerEngine UCE v3.1 (2484行) |
| **S2逆扫描** | `rph_core/steps/step2_retro/retro_scanner.py` | RetroScanner (907行) |
| **S3 TS优化** | `rph_core/steps/step3_opt/ts_optimizer.py` | TSOptimizer (1226行) |
| **S4特征** | `rph_features/rph_features/feature_miner.py` | FeatureMiner (649行) |
| **QC接口** | `rph_core/utils/qc_interface.py` | QCInterfaceFactory, LinuxSandbox, TaskKind (2036行) |
| **Checkpoint** | `rph_core/utils/checkpoint_manager.py` | CheckpointManager, PipelineState (830行) |
| **配置** | `config/defaults.yaml` | 单一配置源 (774行) |
| **形成键** | `rph_core/utils/forming_bonds_resolver.py` | S3.5形成键解析 (282行) |

---

> **报告完成** • 基于 RPH v2.1.1 代码库深度分析 • 2026-05-29
