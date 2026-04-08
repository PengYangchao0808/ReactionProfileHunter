Reaction Profile Hunter (RPH) v2.0.0 - 程序运行流程谱图报告
1. 总览
属性	值
程序名称	Reaction Profile Hunter (RPH)
版本	v2.0.0
任务 SMILES	CC(C)(C)OC(=O)N1CCC[C@@]23C=C[C@@H](CC(=O)[C@H]12)O3
前体 SMILES	C=C=CN(CCCc1ccco1)C(=O)OC(C)(C)C
反应类型	[4+3] 环加成
总耗时	13,521.1 秒 (3 小时 45 分 21 秒)
最终结果摘要
指标	数值
产物 L2 能量	-940.674884 Hartree
前体 L2 能量	-865.331286 Hartree
活化自由能 ΔG‡	26.61 kcal/mol
反应自由能 ΔG_rxn	-24.89 kcal/mol
TS 虚频	-272.3 cm⁻¹
输出文件列表
rph_output/rx_1/
├── S0_Mechanism/
│   ├── mechanism_graph.json
│   └── mechanism_summary.json
├── S1_ConfGeneration/
│   ├── product/
│   │   ├── product_global_min.xyz
│   │   └── dft/conf_000.fchk
│   └── precursor/
│       ├── precursor_global_min.xyz
│       └── dft/conf_000.fchk
├── S2_Retro/
│   ├── ts_guess.xyz
│   ├── intermediate.xyz
│   └── scan_profile.png
├── S3_TS/
│   ├── ts_final.xyz
│   ├── S3_intermediate_opt/
│   └── ASM_SP_Mat/
└── S4_Data/
    ├── features_raw.csv
    ├── features_mlr.csv
    └── feature_meta.json
---
2. 整体架构设计
2.1 模块化管线图
flowchart TD
    subgraph Input["📥 输入层"]
        SMILES[("产物 SMILES\n(S1 锚定目标)")]
        Cleaner[("Cleaner 数据\n(反应元信息)")]
        Config[("defaults.yaml\n(配置参数)")]
    end
    
    subgraph Pipeline["🔄 RPH 四步管线"]
        direction TB
        
        S0["S0: MechanismClassifier<br/>机理分类器"]
        S1["S1: AnchorPhase v3.0<br/>全局最低构象搜索"]
        S2["S2: RetroScanner<br/>TS 初猜构建"]
        S3["S3: TSOptimizer<br/>过渡态优化验证"]
        S4["S4: FeatureMiner<br/>特征提取"]
        
        S0 -->|"forming_bonds<br/>反应类型"| S1
        S1 -->|"product_min.xyz<br/>precursor_min.xyz"| S2
        S2 -->|"ts_guess.xyz<br/>intermediate.xyz"| S3
        S3 -->|"ts_final.xyz<br/>SP 矩阵"| S4
    end
    
    subgraph Engines["⚡ 外部计算引擎"]
        G16["Gaussian 16<br/>几何优化 + 频率"]
        ORCA["ORCA<br/>L2 高精度单点能"]
        XTB["xTB/GFN2<br/>预优化 + 扫描"]
        CREST["CREST<br/>构象搜索"]
        Multiwfn["Multiwfn<br/>波函数分析"]
    end
    
    subgraph Cache["💾 缓存系统"]
        SMC["SmallMoleculeCache<br/>小分子复用"]
        Checkpoint["CheckpointManager<br/>断点续算"]
    end
    
    subgraph Output["📤 输出层"]
        Features[("features_raw.csv<br/>Layer 1+2 特征")]
        QA[("qa_metadata.csv<br/>Layer 3 质量数据")]
        Meta[("feature_meta.json<br/>元数据)")]
    end
    
    SMILES --> S1
    Cleaner --> S0
    Config --> Pipeline
    
    S1 -.->|"查询/存储"| SMC
    S1 -.->|"CREST调用"| CREST
    S1 -.->|"GFN2预优化"| XTB
    S1 -.->|"OPT+Freq"| G16
    S1 -.->|"L2 SP"| ORCA
    
    S2 -.->|"retro_scan<br/>path_search"| XTB
    
    S3 -.->|"Berny TS优化"| G16
    S3 -.->|"L2 SP"| ORCA
    S3 -.->|"Shermo热化学"| S3
    
    S4 -.->|"NBO/Fukui"| Multiwfn
    
    Pipeline -.->|"状态保存"| Checkpoint
    S4 --> Features
    S4 --> QA
    S4 --> Meta
2.2 关键组件职责
组件	类/模块	职责
流程编排	ReactionProfileHunter	总管类，协调 S0-S4 执行，断点续算
机理分类	MechanismClassifier	S0: 反应图构建、键变化检测、环模式识别
构象锚定	AnchorPhase	S1: 两阶段 CREST + DFT OPT-SP 耦合循环
TS 初猜	RetroScanner	S2: 逆扫描 + Path Search，膝点检测
TS 优化	TSOptimizer	S3: Berny/QST2 + 虚频验证 + L2 SP
特征提取	FeatureMiner	S4: 插件化提取器流水线
统一计算中枢	QCTaskRunner	OPT-SP 耦合循环执行器
缓存管理	SmallMoleculeCache	小分子全局缓存，避免重复计算
状态管理	CheckpointManager	步骤级哈希验证断点续算
---
3. 数据流图（DFD）
flowchart LR
    subgraph Input["输入"]
        SMILES["SMILES String"]
    end
    
    subgraph S0_Layer["S0: 机理层"]
        Graph["Reaction Graph<br/>.graphml"]
        Meta["mechanism_summary.json<br/>Forming Bonds: ((12,13),(15,16))"]
    end
    
    subgraph S1_Layer["S1: 构象层"]
        direction TB
        Ensemble["ensemble.xyz<br/>CREST 系综"]
        PM["product_global_min.xyz<br/>E=-940.674884 Ha"]
        PrM["precursor_global_min.xyz<br/>E=-865.331286 Ha"]
        FCHK["*.fchk 检查点文件"]
    end
    
    subgraph S2_Layer["S2: 初猜层"]
        Retro["retro_scan/<br/>xtbscan.log (20 pts)"]
        Path["path_search/<br/>xtbpath.log"]
        TS_Guess["ts_guess.xyz<br/>膝点: 2.337Å"]
        Intermediate["intermediate.xyz<br/>偶极中间体"]
    end
    
    subgraph S3_Layer["S3: 优化层"]
        direction TB
        IntOpt["S3_intermediate_opt/<br/>standard/intermediate_opt.xyz"]
        TS_Opt["ts_opt/berny/<br/>ts_final.xyz"]
        TS_Freq["虚频验证<br/>-272.3 cm⁻¹"]
        L2_SP["L2_SP/<br/>wB97X-D4/def2-TZVPP"]
        SP_Mat["ASM_SP_Mat/<br/>SP 能量矩阵"]
        Thermo["Shermo 热化学<br/>ΔG‡ = 26.61 kcal/mol"]
    end
    
    subgraph S4_Layer["S4: 特征层"]
        Raw["features_raw.csv<br/>(19 columns)"]
        MLR["features_mlr.csv<br/>(ML-ready)"]
        MetaS4["feature_meta.json<br/>(Layer 3 QA)"]
    end
    
    SMILES -->|"解析"| Graph
    Graph -->|"forming_bonds"| S1_Layer
    Graph -->|"reaction_type=[4+3]"| S2_Layer
    
    SMILES -->|"RDKit Embedding<br/>CREST GFN0→GFN2"| Ensemble
    Ensemble -->|"ISOSTAT 聚类<br/>DFT OPT-SP"| PM
    PM --> FCHK
    
    PM -->|"逆扫描起点"| Retro
    Retro -->|"膝点检测"| TS_Guess
    Retro -->|"偶极捕获"| Intermediate
    Intermediate -->|"--path"| Path
    Path -->|"交叉验证"| TS_Guess
    
    TS_Guess -->|"Berny TS"| TS_Opt
    Intermediate -->|"DFT OPT"| IntOpt
    IntOpt -->|"复用为 Reactant"| S3_Layer
    TS_Opt --> TS_Freq
    TS_Opt --> L2_SP
    IntOpt --> L2_SP
    L2_SP --> SP_Mat
    TS_Freq --> Thermo
    L2_SP --> Thermo
    
    SP_Mat --> S4_Layer
    Thermo --> S4_Layer
    PM --> S4_Layer
    PrM --> S4_Layer
    
    S4_Layer --> Raw
    S4_Layer --> MLR
    S4_Layer --> MetaS4
缓存命中点
缓存类型	位置	命中率
SmallMoleculeCache	rph_output/small_molecules/	N/A
SP Cache	内存内	0% (0/1)
Checkpoint Resume	pipeline.state	部分
---
4. 核心算法与工作流详解
4.1 S0 – 反应机制分类
flowchart TD
    A[Cleaner Data<br/>含反应物/产物信息] --> B[构建反应图]
    B --> C[键变化检测]
    C --> D{环模式识别}
    D -->|"检测到 7元环<br/>+ 5元环"| E[[4+3] Cycloaddition]
    D -->|其他模式| F[其他类型]
    E --> G[提取形成键<br/>((12,13),(15,16))]
    G --> H[保存 mechanism_graph.json]
    G --> I[保存 mechanism_summary.json]
算法细节：
- 反应图构建：从 cleaner 数据构建分子图，对比反应物与产物
- 拓扑分类：INTRA_TYPE_I（分子内反应，类型 I）
- 形成键输出：0-based 索引 ((12, 13), (15, 16)) 传递给 S1/S2
---
4.2 S1 – 产物/前体锚定（全局最低构象搜索）
flowchart TD
    subgraph Phase1["阶段 1: CREST 系综搜索"]
        A[RDKit 3D Embedding] --> B[CREST GFN0<br/>快速采样]
        B --> C[ISOSTAT 聚类]
        C --> D[CREST GFN2<br/>精优化]
        D --> E[二次聚类]
    end
    
    subgraph Phase2["阶段 2: DFT OPT-SP 耦合循环"]
        F[选取 Top 3-6 构象] --> G{逐构象优化}
        G -->|conf_000| H[Gaussian OPT<br/>B3LYP/def2-SVP]
        G -->|conf_001| I[Gaussian OPT]
        G -->|conf_002| J[Gaussian OPT]
        G -->|...| K[...]
        
        H --> L[ORCA L2 SP<br/>wB97X-D4/def2-TZVPP]
        I --> L
        J --> L
        K --> L
        
        L --> M[Boltzmann 加权<br/>选择全局最低]
    end
    
    Phase1 --> Phase2
    M --> N[product_global_min.xyz]
性能数据（从日志提取）：
分子	构象数	总耗时
Product	3	~1h
Precursor	6	~2h
优化收敛状态：
- 所有 9 个构象均收敛（无虚频）
- 产物最佳：conf_000 (权重 0.3344, E=-940.674884 Ha)
- 前体最佳：conf_001 (权重 0.1672, E=-865.331286 Ha)
---
4.3 S2 – TS 初猜构建
flowchart TD
    subgraph S2_1["S2.1: Retro Scan (逆扫描)"]
        A[product_min.xyz] --> B[沿形成键拉伸<br/>目标 2.2Å]
        B --> C[xTB GFN2 扫描<br/>20 个扫描点]
        C --> D[能量曲线分析]
        D --> E[膝点检测算法<br/>Knee Point Detection]
        E -->|"TS 初猜"| F[frame_006.xyz<br/>2.337Å]
        E -->|"偶极中间体"| G[frame_010.xyz<br/>2.695Å]
    end
    
    subgraph S2_2["S2.2: Path Search"]
        G --> H[intermediate.xyz] --> I[xTB --path<br/>intermediate → product]
        I --> J[路径优化]
        J --> K[交叉验证 TS 键长<br/>1.586Å]
    end
    
    subgraph Drift["拓扑漂移监控"]
        L[监测新键形成] -->|"(12,17) 键形成<br/>1.397Å → 1.358Å"| M[记录漂移警告]
    end
    
    S2_1 --> S2_2
    S2_1 --> Drift
    S2_2 --> N[ts_guess.xyz]
膝点检测算法 inferred：
- 扫描从产物向反应物方向进行（逆扫描）
- 检测能量曲线的"膝点"（拐点）作为 TS 初猜
- 同时捕获能量峰值和偶极中间体位置
---
4.4 S3 – TS 优化与验证
flowchart TD
    subgraph S3_0["中间体优化"]
        A[intermediate.xyz] --> B[QCTaskRunner<br/>Normal 模式]
        B --> C[Gaussian OPT<br/>B3LYP/def2-SVP]
        C --> D[频率验证]
        D --> E[无虚频 ✓]
        E --> F[ORCA L2 SP]
        F --> G[intermediate_opt.xyz<br/>E=-940.633916 Ha]
    end
    
    subgraph S3_1["TS 优化"]
        H[ts_guess.xyz] --> I[BernyTSDriver]
        I --> J[Gaussian Opt=TS<br/>CalcFC, NoEigenTest]
        J --> K[频率验证]
        K -->|"1 个虚频 ✓<br/>-272.3 cm⁻¹"| L[TS 收敛 ✓]
        K -->|"虚频不符"| M[QST2Rescue<br/>禁用]
        L --> N[ORCA L2 SP<br/>E=-940.589380 Ha]
    end
    
    subgraph S3_2["热化学计算"]
        O[Shermo 分析] --> P[ΔG‡ = 26.61 kcal/mol]
        O --> Q[ΔG_rxn = -25.707 kcal/mol]
    end
    
    G -->|"复用为 Reactant"| S3_2
    N --> S3_2
    S3_0 --> S3_1
关键计算结果：
物种	L2 能量 (Ha)
Reactant (中间体)	-940.633916
TS	-940.589380
ΔG‡	—
---
4.5 S4 – 特征挖掘
flowchart TD
    subgraph Plugins["插件化提取器流水线"]
        A[输入上下文<br/>fchk/xyz/SP矩阵] --> B[thermo<br/>热力学特征]
        A --> C[geometry<br/>几何特征]
        A --> D[qc_checks<br/>QC 验证]
        A --> E[ts_quality<br/>TS 质量指标]
        A --> F[step1_activation<br/>S1 活化特征]
        A --> G[step2_cyclization<br/>S2 环化特征]
        A --> H[multiwfn_features<br/>波函数特征]
    end
    
    subgraph Output["三层输出结构"]
        B --> I[Layer 1<br/>原始特征]
        C --> I
        D --> I
        E --> J[Layer 2<br/>派生特征]
        F --> J
        G --> J
        H --> J
        
        I --> K[deployable_features.csv<br/>19 columns]
        J --> K
        
        B --> L[Layer 3<br/>元数据/质量]
        D --> L
        E --> L
        L --> M[qa_metadata.csv<br/>9 columns]
        
        I --> N[feature_meta.json]
        J --> N
        L --> N
    end
    
    Plugins --> Output
提取器列表：
名称	功能
thermo	热力学量（ΔG, ΔH, ΔS）
geometry	键长、键角、二面角
qc_checks	SCF 收敛、虚频检查
ts_quality	TS 合理性指标
step1_activation	S1 活化能特征
step2_cyclization	S2 环化特征
multiwfn_features	波函数分析
警告记录：
W_S1_MISSING_HOAC_THERMO
W_S1_MISSING_S3_INTERMEDIATE_GIBBS
W_S1_MISSING_THERMO_COMPONENT
---
5. 性能分析与优化建议
5.1 各步骤耗时分析
gantt
    title RPH 执行时间线 (总计 3h 45m)
    dateFormat HH:mm:ss
    axisFormat %H:%M
    
    section S0
    机理分类           :a1, 15:37:43, 1s
    
    section S1
    Product CREST      :b1, 15:37:44, 3m26s
    Product DFT (3×)   :b2, after b1, 56m14s
    Precursor CREST    :b3, after b2, 10m38s
    Precursor DFT (6×) :b4, after b3, 1h49m28s
    
    section S2
    Retro Scan         :c1, 18:37:40, 7s
    Path Search        :c2, after c1, 7s
    
    section S3
    Intermediate OPT   :d1, 18:37:56, 20m39s
    TS Berny OPT       :d2, after d1, 19m27s
    L2 SP & Shermo     :d3, after d2, 3m04s
    
    section S4
    特征提取           :e1, 19:23:03, 1s
5.2 热点分析
排名	热点	耗时
1	Precursor DFT OPT	~109 min
2	Product DFT OPT	~56 min
3	Intermediate OPT	~21 min
4	TS Berny OPT	~19 min
5	CREST 搜索	~14 min
5.3 瓶颈识别
flowchart LR
    A[串行构象优化<br/>瓶颈: 9个构象串行<br/>影响: ~60% 总时间] --> B[建议: 并行化]
    C[SP Cache 命中 0%] --> D[建议: 启用 L2 SP 缓存]
    E[路径含空格<br/>临时目录拷贝开销] --> F[建议: 避免工作目录含特殊字符]
    G[未启用 QST2 Rescue<br/>TS 失败时无备选] --> H[建议: 按需启用]
5.4 优化建议
优先级	建议	预期收益
🔴 高	并行构象优化	减少 S1 时间 60-70%
🔴 高	增强 L2 SP 缓存	避免重复单点能计算
🟡 中	路径规范化	消除临时目录拷贝开销
🟡 中	智能构象筛选	DFT 前用 xTB 能量预筛选
🟢 低	启用 QST2 Rescue	TS 失败时提供备选方案
---
6. 完整谱图（时序图）
sequenceDiagram
    autonumber
    participant Main as Main Process
    participant S0 as S0: Mechanism
    participant S1 as S1: AnchorPhase
    participant S2 as S2: RetroScanner
    participant S3 as S3: TSOptimizer
    participant S4 as S4: FeatureMiner
    participant Cache as SmallMoleculeCache
    participant XTB as xTB/GFN2
    participant G16 as Gaussian 16
    participant ORCA as ORCA
    participant Shermo as Shermo
    
    %% S0
    rect rgb(230, 245, 255)
        Main->>S0: classify_from_dict(cleaner_data)
        S0-->>Main: mechanism_summary.json<br/>forming_bonds: ((12,13),(15,16))
    end
    
    %% S1 - Product
    rect rgb(255, 245, 230)
        Main->>S1: anchor(product)
        S1->>Cache: get_or_create(smiles)
        Cache-->>S1: cache_path
        S1->>XTB: CREST GFN0 Search
        XTB-->>S1: crest_best.xyz
        S1->>XTB: CREST GFN2 Refinement
        XTB-->>S1: ensemble.xyz (3 confs)
        
        loop 3 Conformers
            S1->>G16: Gaussian OPT B3LYP/def2-SVP
            G16-->>S1: conf_xxx.log (converged)
            S1->>ORCA: L2 SP wB97X-D4/def2-TZVPP
            ORCA-->>S1: SP energy
        end
        
        S1-->>Main: product_global_min.xyz<br/>E=-940.674884 Ha
    end
    
    %% S1 - Precursor
    rect rgb(255, 245, 230)
        Main->>S1: anchor(precursor)
        S1->>XTB: CREST Two-Stage
        XTB-->>S1: ensemble.xyz (6 confs)
        
        loop 6 Conformers
            S1->>G16: Gaussian OPT
            G16-->>S1: converged
            S1->>ORCA: L2 SP
            ORCA-->>S1: energy
        end
        
        S1-->>Main: precursor_global_min.xyz<br/>E=-865.331286 Ha
    end
    
    %% S2
    rect rgb(230, 255, 230)
        Main->>S2: run_retro_scan(product_min.xyz)
        S2->>XTB: xtb --opt --input scan.inp
        XTB-->>S2: xtbscan.log (20 points)
        S2->>S2: knee_point_detection()
        S2-->>Main: ts_guess.xyz (2.337Å)<br/>intermediate.xyz
        
        Main->>S2: run_path_search()
        S2->>XTB: xtb --path intermediate product
        XTB-->>S2: path optimized
        S2-->>Main: path validation complete
    end
    
    %% S3
    rect rgb(255, 230, 245)
        Main->>S3: optimize_intermediate()
        S3->>G16: Normal OPT + Freq
        G16-->>S3: intermediate_opt.xyz (no imag)
        S3->>ORCA: L2 SP
        ORCA-->>S3: E=-940.633916 Ha
        
        Main->>S3: optimize_ts(ts_guess.xyz)
        S3->>G16: Berny TS Opt + Freq
        G16-->>S3: ts_final.xyz<br/>imag = -272.3 cm⁻¹
        S3->>ORCA: L2 SP
        ORCA-->>S3: E=-940.589380 Ha
        
        Main->>Shermo: compute_thermochemistry()
        Shermo-->>Main: ΔG‡ = 26.61 kcal/mol
    end
    
    %% S4
    rect rgb(245, 230, 255)
        Main->>S4: extract_features()
        
        par Plugin Execution
            S4->>S4: thermo extractor
            S4->>S4: geometry extractor
            S4->>S4: ts_quality extractor
            S4->>S4: step1_activation extractor
            S4->>S4: step2_cyclization extractor
        end
        
        S4-->>Main: features_raw.csv<br/>feature_meta.json
    end
    
    Main->>Main: Pipeline SUCCESS
---
7. 附录
7.1 关键文件清单
输入文件
文件	路径
CLI 入口	bin/rph_run
配置文件	config/defaults.yaml
运行日志	rph_output/rx_1/rph.log
中间文件
步骤	文件	路径
S0	机理图	S0_Mechanism/mechanism_graph.json
S0	摘要	S0_Mechanism/mechanism_summary.json
S1	产物全局最低	S1_ConfGeneration/product/product_global_min.xyz
S1	前体全局最低	S1_ConfGeneration/precursor/precursor_global_min.xyz
S1	检查点	S1_ConfGeneration/product/dft/conf_000.fchk
S2	TS 初猜	S2_Retro/ts_guess.xyz
S2	中间体	S2_Retro/intermediate.xyz
S2	扫描图	S2_Retro/scan_profile.png
S3	优化后 TS	S3_TS/ts_final.xyz
S3	SP 矩阵	S3_TS/ASM_SP_Mat/
S3	工件索引	S3_TS/artifacts_index.json
输出文件
文件	路径
原始特征	S4_Data/features_raw.csv
ML 特征	S4_Data/features_mlr.csv
特征元数据	S4_Data/feature_meta.json
7.2 环境配置摘要
executables:
  gaussian:
    path: "/opt/software/gaussian/g16/g16"
    root: "/opt/software/gaussian/g16"
  orca:
    path: "/opt/software/orca/orca"
    ld_library_path: "/opt/openmpi418/lib:/opt/software/orca"
  xtb:
    path: "/opt/software/xtb/bin/xtb"
  crest:
    path: "/opt/software/crest/crest"
  multiwfn:
    path: "/opt/software/multiwfn/Multiwfn"
resources:
  mem: "32GB"
  nproc: 16
  orca_maxcore_safety: 0.2
theory:
  optimization:
    method: "B3LYP"
    basis: "def2-SVP"
    dispersion: "GD3BJ"
    engine: "gaussian"
  single_point:
    method: "wB97X-D4"
    basis: "def2-TZVPP"
    engine: "orca"
LD_LIBRARY_PATH:
/opt/openmpi418/lib:/opt/software/orca:/opt/software/gaussian/g16/bsd:/opt/software/gaussian/g16:/opt/software/gaussian/gv/lib
---
报告生成信息
项目	值
分析日期	2026-04-07
RPH 版本	v2.0.0
分析基于	完整运行日志 (13,521.1 秒)
数据来源	rph_output/rx_1/
补充：详细架构验证
背景任务已完成代码库深度探查，以下是与我生成的报告的关键验证点：
模块-文件映射确认
步骤	主文件	关键类
S0	steps/mechanism_classifier/classifier.py	MechanismClassifier
S1	steps/anchor/handler.py	AnchorPhase
S2	steps/step2_retro/retro_scanner.py	RetroScanner
S3	steps/step3_opt/ts_optimizer.py	TSOptimizer
S4	steps/step4_features/feature_miner.py	FeatureMiner
S4 插件清单（完整版）
# 从 extractors/__init__.py 注册
registered_extractors = [
    "thermo",           # 热力学量
    "geometry",         # 几何特征
    "qc_checks",        # QC 验证
    "interaction",      # 相互作用分析
    "nics",             # NICS 芳香性
    "nbo_e2",           # NBO E2 分析
    "ts_quality",       # TS 质量指标
    "asm_enrichment",   # ASM 增强
    "step1_activation", # S1 活化特征
    "step2_cyclization",# S2 环化特征
    "multiwfn_features",# Multiwfn 波函数
    "fmo_cdft_dipolar"  # FMO/CDFT 分析
]
关键调用链验证
Orchestrator.run_pipeline()
├── S0: MechanismClassifier.classify_from_dict()
│   └── 输出: mechanism_summary.json
├── S1: AnchorPhase.run()
│   └── ConformerEngine._step_two_stage_crest()
│       ├── CREST GFN0 快速采样
│       └── CREST GFN2 精优化
│   └── ConformerEngine._step_dft_opt_sp_coupled()
│       ├── Gaussian OPT (B3LYP/def2-SVP)
│       └── ORCA SP (wB97X-D4/def2-TZVPP)
├── S2: RetroScanner.run_retro_scan()
│   └── XTBInterface.scan() [xTB $scan]
├── S3: TSOptimizer.run_with_qctaskrunner()
│   ├── QCTaskRunner.run_opt_sp_cycle() [中间体]
│   └── QCTaskRunner.run_ts_opt_cycle() [TS]
│       ├── BernyTSDriver.optimize()
│       └── Shermo 热化学
└── S4: FeatureMiner.run()
    └── 12 extractors 并行执行