# PlanI — ReactionProfileHunter 通用正向扫描过渡态搜索模块 完整工作计划

**版本**: v3.0-Universal  
**日期**: 2026-03-10  
**状态**: 待实施  
**适用范围**: 任意环加成反应（[4+2]、[3+2]、[4+3]、[5+2]、[2+2]、Claisen/Cope 等）

---

## 一、背景与动机

### 1.1 现有流程的局限

ReactionProfileHunter（RPH）v6.2 的 S1→S4 流水线目前**仅支持 [5+2] 环加成反应**。这一限制体现在以下关键位置：

- **S2 逆向扫描器**（retro_scanner.py，470 行）：采用"产物→过渡态"逆向拉伸策略，将产物中的成键键拉长至 2.2Å 后进行约束优化，获得 TS 猜测结构。此策略的 SMARTS 匹配（smarts_matcher.py，252 行）**完全硬编码为 [3.2.1] 氧杂桥环拓扑**，无法识别其他反应类型的成键位点。
- **分子内片段切割器**（intramolecular_fragmenter.py，565 行）：其 dipole_core_indices 假设偶极核心为 5 个原子（氧化吡喃鎓），这是 [5+2] 特有的。
- **语义切片器**（semantic_slicer.py）、**片段操作模块**（fragment_manipulation.py）、**S4 环化提取器**（step2_cyclization.py）等共计 7 个文件、32 处硬编码引用了 [5+2] 专属逻辑。
- **配置文件**（defaults.yaml，389 行）中无反应类型配置段，所有理论级别和参数均为 [5+2] 默认值。

### 1.2 正向扫描策略的理论依据

对于 [4+3] 等可能为分步机理的反应，传统逆向扫描策略存在根本困难：产物经逆向拉伸可能回到中间体（INT）而非真正的过渡态。

基于以下理论和文献支持，我们采用**正向扫描策略**（从反应物/中间体出发，压缩成键键距寻找能量极大点）：

- **Curtin-Hammett 原理**：对分步反应，中间体到各产物的过渡态之间的能量差决定产物分布，因此从 INT 出发正向搜索 TS₂ 是化学上合理的。
- **GSM 方法**（Zimmerman 2013, 2015）：SE-GSM 使用反应物端的驱动坐标，已在多种周环反应上验证成功。
- **XTBDFT 流程**（Lin 2022）：通用的 xTB→DFT 工作流，不限定反应类型。
- **Houk 组计算研究**：IRC 计算确认协同环加成反应（[4+2]、[3+2] 等）存在单一过渡态，正向/逆向扫描等价。
- **xTB 原生扫描**：xTB 的 $scan 功能块支持多坐标协同扫描（mode=concerted），内部自动处理弛豫优化，无需 Python 循环。

### 1.3 通用化设计理念

本计划的核心思想是：**正向扫描流程本身与反应类型无关**。对于任何环加成反应，模块只需要三个输入：

1. **需要压缩的键列表**（从 cleaner 数据或 SMARTS 匹配获得）
2. **优化后的反应物/中间体结构**（从 S1/S2 获得）
3. **扫描参数**（起始距离、终止距离、步数，从配置文件读取）

反应类型的差异仅体现在：成键数量、理论级别选择、溶剂模型、IRC 验证容忍度等，这些均通过 defaults.yaml 中的 reaction_profiles 配置段参数化，无需修改代码逻辑。

---

## 二、总体架构设计

### 2.1 模块定位

正向扫描模块作为 S2 阶段的**替代执行路径**存在，与现有逆向扫描器（RetroScanner）平行：

- **逆向扫描**（retro_scan）：适用于协同反应（[5+2] 默认），从产物逆向拉伸
- **正向扫描**（forward_scan）：适用于分步或通用场景，从反应物/中间体正向压缩

由 orchestrator.py 根据 reaction_profiles 配置中的 s2_strategy 字段自动路由。

### 2.2 数据流概览

输入阶段（S1 之前）：
- cleaner 项目导出的 TSV/JSON 数据 → cleaner_adapter 转换为 RPH 的 ReactionRecord
- ReactionRecord 中包含 core_bond_changes 字段（如 "1-2:formed;5-6:formed;3-4:broken"）
- 若无 cleaner 数据，回退至 SMARTS 模板匹配

S1 阶段：
- 锚定产物构象搜索 → product_min.xyz、能量、热力学数据

S2 阶段（正向扫描路径）：
- 读取 S1 产物结构 + core_bond_changes 中的 formed 键列表
- 将产物中的 formed 键拉长至扫描起始距离（默认 3.5Å）→ 无约束 xTB 优化 → 反应物复合物
- 以反应物复合物为起点，对所有 formed 键同时压缩至终止距离（默认 1.5Å）
- xTB 原生 $scan（mode=concerted）执行弛豫扫描
- 解析 xtbscan.log，提取能量最高点几何结构作为 ts_guess.xyz
- 同时输出 reactant_complex.xyz

S3 阶段：
- 接收 ts_guess.xyz + reactant_complex.xyz + product_min.xyz
- 根据 reaction_profiles 中的理论级别执行 DFT TS 优化
- Berny → 振荡检测 → QST2 救援（现有流程不变）
- 对分步反应，IRC 验证时需容忍连接到中间体而非最终反应物

S4 阶段：
- 特征提取需适配可变成键数量（不再硬编码 2 键假设）
- 片段切割需参数化核心原子数（不再硬编码 5 原子偶极核心）

### 2.3 配置驱动的反应类型参数化

在 defaults.yaml 中新增 reaction_profiles 配置段。每个配置项定义一种反应类型的完整参数集，包括：

- 反应类型标识符（如 "[4+3]_furan_allenamide"）
- 成键数量（forming_bond_count）—— 可设为自动（从 cleaner 数据推断）
- S2 策略选择（s2_strategy: retro_scan 或 forward_scan）
- 扫描参数：起始距离、终止距离、步数、扫描模式、力常数
- S3 理论级别：优化方法/基组、单点方法/基组
- 溶剂模型及溶剂名称
- IRC 验证配置（是否允许连接到中间体）
- 通用默认配置项（_universal）作为兜底，所有未指定的反应类型均使用此配置

现有 [5+2] 默认行为通过 "[5+2]_default" 配置项保留，确保**向后兼容**。

---

## 三、详细修改方案

### 模块 M1：xTB 扫描基础设施（优先级 P0）

**目标**：为 xTB runner 添加原生 $scan 扫描能力。

**涉及文件**：
- rph_core/utils/xtb_runner.py（328 行）
- rph_core/utils/data_types.py（55 行）
- rph_core/utils/qc_interface.py（1542 行）

**具体改动**：

（一）data_types.py —— 新增 ScanResult 数据类

在现有 QCResult 数据类之后，新增 ScanResult 数据类，包含以下字段：
- success（布尔）：扫描是否成功完成
- energies（浮点数列表）：每个扫描点的能量（Hartree）
- geometries（路径列表或单一 xtbscan.log 路径）：所有扫描点的几何结构
- max_energy_index（整数）：能量最高点的索引
- ts_guess_xyz（路径）：从能量最高点提取的几何结构文件
- scan_log（路径）：xtbscan.log 原始输出文件

（二）xtb_runner.py —— 新增扫描方法

在 XTBRunner 类中新增以下方法：

- run_scan 方法：接收输入结构文件路径、约束-扫描参数字典（包含 force_constant、距离约束列表、扫描区间、步数、扫描模式）、溶剂、电荷、UHF 等参数。内部调用 _write_scan_input 生成输入文件，通过 subprocess 执行 xTB，然后调用 _parse_scan_log 解析结果，返回 ScanResult。

- _write_scan_input 方法：生成包含 $constrain、$scan、$opt 三个功能块的输入文件。约束块写入力常数和所有距离约束（原子编号为 1-based）；扫描块写入 mode=concerted 和每个约束的扫描区间/步数（注意：concerted 模式要求所有约束的步数必须相同）；优化块写入 maxcycle 参数。

- _parse_scan_log 方法：解析 xtbscan.log 文件（XMol 格式）。该文件结构为：每个扫描点由原子数行、注释行（包含 "SCF done" 及能量值）、N 行坐标组成，依次排列。方法需提取所有扫描点的能量，找到能量最大值索引，并将该点的几何结构写出为独立的 xyz 文件。

**关键约束**：
- 原子编号转换：RPH 内部使用 0-based 索引，xTB 约束文件使用 1-based 索引。_write_scan_input 中必须执行 +1 转换，与现有 _write_constraint_input 保持一致。
- 力常数统一：现有 _write_constraint_input 使用 1.0，retro_scanner 使用 0.5。新的扫描方法应从配置读取力常数（默认 0.5），不硬编码。
- xtbscan.log 必须在扫描完成后存在；若文件不存在或为空，run_scan 返回 success=False 并记录错误信息。

（三）qc_interface.py —— 扩展任务类型

在 TaskKind 枚举（第 110-117 行）中新增 SCAN 类型。

在 XTBInterface 类中新增 scan 方法，作为 run_scan 的门面，负责：沙盒路径检查（调用 is_path_toxic）、工作目录创建、调用 XTBRunner.run_scan、结果日志记录。所有 QC 调用必须通过此门面，不得在步骤代码中直接调用 XTBRunner。

### 模块 M2：正向扫描执行器（优先级 P0）

**目标**：在 S2 中实现通用正向扫描 TS 搜索方法。

**涉及文件**：
- rph_core/steps/step2_retro/retro_scanner.py（470 行）
- rph_core/steps/step2_retro/bond_stretcher.py（133 行）

**具体改动**：

（一）retro_scanner.py —— 新增 run_forward_scan 方法

在 RetroScanner 类中新增 run_forward_scan 方法，其签名与现有 run 方法对齐，返回类型相同：(ts_guess_xyz 路径, reactant_complex_xyz 路径, forming_bonds 元组)。

执行流程：

第一步 —— 生成反应物复合物：
- 读取 S1 产物优化结构（product_min.xyz）
- 从 forming_bonds 参数中获取需要形成的键列表（原子索引对）
- 调用 bond_stretcher 将这些键拉长至扫描起始距离（默认 3.5Å，从配置读取）
- 对拉长后的结构进行无约束 xTB 优化 → 得到 reactant_complex.xyz

第二步 —— 执行正向扫描：
- 以 reactant_complex.xyz 为起点
- 构造扫描参数：每个 formed 键从 scan_start_distance 压缩至 scan_end_distance，共 scan_steps 步，mode=concerted
- 通过 qc_interface.scan（而非直接调用 xtb_runner）执行扫描
- 接收 ScanResult

第三步 —— 提取 TS 猜测结构：
- 从 ScanResult.energies 中找到能量极大点
- 若极大点位于扫描端点（首点或末点），记录警告并仍然使用该点，但在日志中标记为"边界极值——建议调整扫描范围"
- 将极大点几何结构写出为 ts_guess.xyz

第四步 —— 输出：
- 返回 (ts_guess.xyz 路径, reactant_complex.xyz 路径, forming_bonds)
- 与现有 run 方法的返回格式完全一致，确保下游 S3 无感切换

（二）bond_stretcher.py —— 泛化键数量

现有 stretch_two_bonds 方法（第 114 行）硬编码了 2 键假设。需改为 stretch_bonds 方法，接受任意数量的（原子索引对, 目标距离）列表。对每对原子，沿键方向线性插值调整距离至目标值。保留原方法名作为兼容别名，内部调用新方法。

### 模块 M3：流水线路由与配置（优先级 P0）

**目标**：使 orchestrator 能根据反应类型配置自动选择 S2 执行路径。

**涉及文件**：
- rph_core/orchestrator.py（1095 行）
- config/defaults.yaml（389 行）

**具体改动**：

（一）defaults.yaml —— 新增 reaction_profiles 配置段

在文件末尾（现有 theory 段之后）新增 reaction_profiles 段，包含以下预定义配置：

"[5+2]_default" 配置：
- forming_bond_count: 2
- s2_strategy: retro_scan
- 其余参数沿用现有默认值（ts_distance: 2.2, break_distance: 3.5）
- 理论级别沿用现有 theory 段
- 此配置确保向后兼容

"[4+3]_default" 配置：
- forming_bond_count: 2
- s2_strategy: forward_scan
- scan_start_distance: 3.5
- scan_end_distance: 1.5
- scan_steps: 20
- scan_mode: concerted
- scan_force_constant: 0.5
- s3_method: B3LYP
- s3_basis: 6-31G(d)
- s3_sp_method: M06-2X
- s3_sp_basis: 6-311+G(d,p)
- solvent: dichloromethane
- solvent_model: cpcm
- irc_allow_intermediate: true

"_universal" 配置（兜底默认）：
- forming_bond_count: auto（从 cleaner 数据推断）
- s2_strategy: forward_scan
- 扫描参数同 [4+3] 默认值
- 理论级别沿用现有 theory 段

（二）orchestrator.py —— S2 路由分发

在 run_pipeline 方法的 S2 调用处（约第 534-541 行），添加路由逻辑：

- 从配置中读取当前反应的 reaction_profile（通过 reaction_type 字段匹配，或使用 _universal 兜底）
- 根据 s2_strategy 字段选择执行路径：
  - "retro_scan" → 调用现有 RetroScanner.run()
  - "forward_scan" → 调用新的 RetroScanner.run_forward_scan()
- 两条路径的返回值格式完全相同：(ts_guess_xyz, reactant_xyz, forming_bonds)
- **关键修复**：现有代码第 541 行使用下划线丢弃了 forming_bonds 返回值。必须修改为显式接收并传递给 S3，而非在 S3 后通过几何比较重新推断。这对正向扫描路径尤为重要，因为正向扫描直接基于已知的 forming_bonds 执行。

（三）orchestrator.py —— 反应类型参数传递

在 CLI 参数解析（main 函数，约第 1036 行）中新增可选参数：
- --reaction-type：反应类型标识符（如 "[4+3]"、"[5+2]"），用于匹配 reaction_profiles 配置
- 默认值为 "[5+2]_default"，确保无参数时行为与现有完全一致

在 ReactionProfileHunter 类初始化时，根据 reaction_type 加载对应的 reaction_profile 配置，合并到运行时配置中。

### 模块 M4：Cleaner 数据适配器（优先级 P0）

**目标**：将 cleaner 项目导出的实验数据转换为 RPH 内部格式，提供通用的成键信息来源。

**涉及文件**：
- rph_core/utils/cleaner_adapter.py（新建文件）
- rph_core/utils/tsv_dataset.py（现有文件，需泛化）

**具体改动**：

（一）cleaner_adapter.py —— 新建 Cleaner→RPH 转换器

该模块负责读取 cleaner 导出的 TSV/JSON 文件，将每条记录转换为 RPH 的 ReactionRecord。

核心逻辑：
- 解析 core_bond_changes 字段（格式："1-2:formed;5-6:formed;3-4:broken"），提取所有 formed 键的原子映射编号对
- 将原子映射编号转换为分子内部索引（通过 SMILES 中的原子映射关系）
- 处理 map_status 字段：仅当 map_status 为 "OK"（置信度≥0.8）或 "LOW_CONFIDENCE" 时使用 cleaner 数据；map_status 为 "FAILED" 时回退到 SMARTS 匹配
- 提取 reaction_type 字段（"[4+3]"、"[5+2]" 等），用于匹配 reaction_profiles 配置
- 提取 product_smiles、reactant_smiles、溶剂信息等

（二）tsv_dataset.py —— 泛化字段名

现有 ReactionRecord 的字段名可能与 cleaner 导出格式不完全对应。需确保以下字段的映射：
- cleaner 的 canonical_product → RPH 的 product_smiles
- cleaner 的 canonical_reactants → RPH 的 reactant_smiles
- cleaner 的 core_bond_changes → RPH 的 forming_bonds + breaking_bonds
- cleaner 的 reaction_type → RPH 的 reaction_type
- cleaner 的 solvent（若有）→ RPH 的 solvent_override

### 模块 M5：SMARTS 模板注册表（优先级 P0）

**目标**：将 SMARTS 匹配从 [5+2] 硬编码改为基于反应类型的模板注册表，作为无 cleaner 数据时的回退方案。

**涉及文件**：
- rph_core/steps/step2_retro/smarts_matcher.py（252 行）

**具体改动**：

将现有的 _topological_core_identification 方法（硬编码 [3.2.1] 氧杂桥环拓扑）重构为模板注册表模式：

- 定义 SMARTSTemplate 数据结构，包含：反应类型标识符、产物 SMARTS 模式、成键原子在 SMARTS 中的位置索引、核心原子数
- 预置模板：
  - "[5+2]"：现有 [3.2.1] 氧杂桥环拓扑匹配逻辑（保持现有行为）
  - "[4+3]"：7 元环产物拓扑，oxyallyl 3 原子核心
  - "[4+2]"：6 元环产物拓扑（Diels-Alder）
  - "[3+2]"：5 元环产物拓扑（1,3-偶极环加成）
- find_reactive_bonds 方法改为：先尝试 cleaner 数据（通过 cleaner_adapter），若不可用则根据 reaction_type 从模板注册表中选择对应模板执行匹配
- 模板注册表存储在 defaults.yaml 的 reaction_profiles 中（SMARTS 模式字符串），而非硬编码在 Python 文件中

**关键原则**：cleaner 数据是成键信息的**主要来源**，SMARTS 匹配仅为回退方案。当 cleaner 提供了高置信度（≥0.8）的 core_bond_changes 时，直接使用，不再执行 SMARTS 匹配。

### 模块 M6：S3 TS 优化适配（优先级 P1）

**目标**：使 S3 阶段能适应不同反应类型的理论级别和验证需求。

**涉及文件**：
- rph_core/steps/step3_opt/ts_optimizer.py（914 行）
- rph_core/steps/step3_opt/berny_driver.py（128 行）
- rph_core/steps/step3_opt/validator.py（217 行）
- rph_core/steps/step3_opt/intramolecular_fragmenter.py（565 行）

**具体改动**：

（一）ts_optimizer.py —— 理论级别参数化

修改 S3 入口方法，使其从 reaction_profile 配置中读取理论级别（方法、基组、溶剂模型），而非使用硬编码的 defaults.yaml 顶层 theory 段。具体而言：

- 优化级别：读取 reaction_profile.s3_method 和 reaction_profile.s3_basis
- 单点级别：读取 reaction_profile.s3_sp_method 和 reaction_profile.s3_sp_basis
- 溶剂：读取 reaction_profile.solvent 和 reaction_profile.solvent_model
- 若 reaction_profile 中未指定，回退到 theory 段默认值

（二）validator.py —— IRC 验证容忍度

对于分步反应（如 [4+3]），IRC 验证可能发现过渡态连接的是中间体而非最终反应物。需新增配置项 irc_allow_intermediate（布尔值，默认 false）：

- 当 irc_allow_intermediate 为 true 时，IRC 只需验证：
  - 正向端连接到产物（或能量更低的结构）
  - 反向端连接到任意局部极小值（中间体或反应物均可）
- 当 irc_allow_intermediate 为 false 时，保持现有严格验证逻辑

（三）intramolecular_fragmenter.py —— 参数化核心原子数

现有 _find_dipole_core_path 方法硬编码搜索 5 原子偶极核心路径。需参数化：

- 从 reaction_profile 配置中读取 core_atom_count（默认 5，[4+3] 为 3，[4+2] 为 4）
- 将路径搜索算法改为根据 core_atom_count 搜索对应长度的最短路径
- 变量名从 dipole_core_indices 泛化为 reactive_core_indices，但保留 dipole_core_indices 作为兼容别名

### 模块 M7：S4 特征提取适配（优先级 P1）

**目标**：使 S4 特征提取能处理可变成键数量和不同片段切割方式。

**涉及文件**：
- rph_core/steps/step4_features/extractors/step2_cyclization.py
- rph_core/steps/step4_features/context.py（294 行）
- rph_core/utils/fragment_cut.py
- rph_core/utils/fragment_manipulation.py
- rph_core/utils/semantic_slicer.py

**具体改动**：

（一）step2_cyclization.py —— 泛化成键数量检查

第 148 行的 2 键硬编码检查改为从上下文中读取预期成键数量（来自 forming_bonds 列表的长度），而非固定为 2。

（二）fragment_cut.py —— 泛化切割逻辑

第 133 行的 2 键假设改为根据实际 forming_bonds 数量执行切割。对于 n 条成键，切割产生 n+1 个或更少的片段（取决于拓扑连通性）。

（三）fragment_manipulation.py —— 泛化电荷分配

第 159 行的氧化吡喃鎓（oxidopyrylium）电荷分配逻辑改为基于反应类型的参数化分配：
- [5+2] oxidopyrylium：+1 电荷在偶极核心
- [4+3] oxyallyl：zwitterion 电荷分布
- 其他类型：中性默认，或从配置读取

（四）semantic_slicer.py —— 泛化核心识别

第 213 行的"偶极核心"识别逻辑改为使用 reactive_core_indices（从 intramolecular_fragmenter 的参数化结果获取），不再假设 5 原子固定长度。

（五）context.py —— 扩展上下文信息

在 FeatureContext 中新增字段：
- reaction_type：当前反应类型标识符
- forming_bond_count：成键数量
- reactive_core_indices：反应核心原子索引列表
- reaction_profile：完整的反应配置字典

### 模块 M8：检查点与批处理适配（优先级 P2）

**目标**：确保检查点系统和批处理模式能正确处理多反应类型场景。

**涉及文件**：
- rph_core/utils/checkpoint_manager.py（603 行）
- rph_core/orchestrator.py（run_batch 方法）

**具体改动**：

（一）checkpoint_manager.py —— 检查点兼容性

现有检查点通过 artifact hash 验证恢复。当 reaction_type 或 s2_strategy 配置变更时，S2 及之后的检查点应失效（因为执行路径不同）。需在检查点元数据中记录 reaction_type 和 s2_strategy，恢复时比对。

（二）orchestrator.py run_batch —— 混合批处理

当批处理中包含不同反应类型的分子时，需在每条记录级别解析 reaction_type 并加载对应 reaction_profile。现有批处理逻辑假设所有记录共享同一配置，需改为逐记录配置解析。

---

## 四、实施优先级与依赖关系

### P0 —— 核心功能（必须首先完成，正向扫描流程才能运行）

实施顺序严格按依赖关系排列：

1. **M1**（xTB 扫描基础设施）：所有上层模块的基础。data_types.py → xtb_runner.py → qc_interface.py，依次完成。
2. **M4**（Cleaner 数据适配器）：M5 依赖其提供成键信息。可与 M1 并行开发。
3. **M5**（SMARTS 模板注册表）：M2 的 run_forward_scan 需要成键信息来源。依赖 M4 定义的接口。
4. **M3**（流水线路由与配置）：defaults.yaml 配置段可先完成；orchestrator 路由逻辑依赖 M2 的 run_forward_scan 方法存在。
5. **M2**（正向扫描执行器）：依赖 M1（xTB 扫描能力）和 M5（成键信息获取）。是正向扫描的核心实现。

**M1 和 M4 可并行开发**，然后 M5，然后 M3 和 M2（M3 的配置部分可提前，路由部分与 M2 同步完成）。

### P1 —— 下游适配（S3/S4 能正确处理新反应类型）

6. **M6**（S3 适配）：依赖 M3 提供的 reaction_profile 配置。可在 M2 完成后立即开始。
7. **M7**（S4 适配）：依赖 M6 完成后的 S3 输出格式确定。但多数改动是独立的泛化重构，可与 M6 并行进行。

### P2 —— 系统健壮性

8. **M8**（检查点与批处理）：最后完成，需要所有其他模块的接口稳定。

### 依赖图

M1（xTB扫描）与 M4（Cleaner适配）可并行 → M5（SMARTS注册表）依赖 M4 → M2（正向扫描执行器）依赖 M1 和 M5 → M3（配置+路由）的配置部分可提前，路由部分与 M2 同步 → M6（S3适配）和 M7（S4适配）依赖 M2/M3 完成 → M8（检查点）最后完成

---

## 五、向后兼容性保证

本计划的所有修改必须确保以下不变量：

1. **无参数运行等价**：当用户不指定 --reaction-type 时，系统默认使用 "[5+2]_default" 配置，执行逆向扫描策略，行为与修改前完全一致。

2. **现有测试通过**：tests/ 目录下的 66 个测试文件（特别是 tests/tmp_v2_2_test/ 中的 Diels-Alder S1 fixture 树）必须全部通过，无需修改测试代码。

3. **RetroScanner.run() 不变**：现有 run 方法的签名、逻辑、返回值均保持不变。run_forward_scan 是新增方法，不影响现有代码路径。

4. **S3/S4 接口兼容**：S3 的 ts_optimizer.optimize() 和 S4 的 feature_miner.extract() 的必选参数不变。新增的 reaction_type、reactive_core_indices 等均为可选参数，缺省时退化为现有行为。

5. **配置文件兼容**：defaults.yaml 中新增的 reaction_profiles 段对现有 step2/step3/step4 段无影响。现有配置文件无需任何修改即可继续使用。

---

## 六、测试策略

### 单元测试

每个模块至少需要以下测试：

**M1 测试**：
- _write_scan_input 生成的输入文件格式是否正确（包含正确的 $constrain/$scan/$opt 块）
- _parse_scan_log 能否正确解析多点 XMol 格式的 xtbscan.log
- 原子编号 0-based → 1-based 转换是否正确
- concerted 模式下步数不一致时是否报错

**M2 测试**：
- 给定已知产物结构和成键信息，run_forward_scan 是否输出合理的 ts_guess.xyz 和 reactant_complex.xyz
- 能量极大点位于边界时是否正确记录警告
- 返回值格式是否与 run 方法一致

**M4 测试**：
- core_bond_changes 字符串解析是否正确（各种格式、边界情况）
- map_status 过滤逻辑（OK → 使用，FAILED → 回退）
- 原子映射编号到内部索引的转换

**M5 测试**：
- 已有 [5+2] 模板的匹配行为与修改前一致
- 新增模板能正确匹配对应反应类型的产物拓扑

### 集成测试

- 端到端正向扫描流程：给定 [4+3] 反应的 SMILES 和 cleaner 数据，能否完成 S1→S2（正向扫描）→S3→S4 全流程
- 端到端逆向扫描流程：给定 [5+2] 反应的 SMILES，确认现有流程不受影响
- 混合批处理：一个批次中同时包含 [5+2] 和 [4+3] 反应，验证分别走入正确的 S2 路径

### 回归测试

- 运行现有全部 66 个测试文件，确认无回归
- 特别关注 tests/test_imports_step4_features.py 等导入烟雾测试

---

## 七、风险与缓解

| 风险 | 影响 | 概率 | 缓解措施 |
|------|------|------|----------|
| xTB concerted scan 对某些反应类型收敛失败 | S2 无法产生 TS 猜测 | 中 | 添加回退逻辑：concerted 失败时改为逐键顺序扫描（sequential scan），取各步最高能量点 |
| 正向扫描能量曲面单调递增/递减，无极大点 | TS 猜测质量差 | 低 | 检测边界极值并扩大扫描范围自动重试（最多一次），仍失败则标记为需人工检查 |
| cleaner 原子映射置信度低导致错误成键信息 | 整个流程基于错误的化学键执行 | 中 | 仅在 map_status 为 OK 时直接使用；LOW_CONFIDENCE 时与 SMARTS 匹配交叉验证 |
| [4+3] 分步机理 IRC 验证失败（连接中间体） | S3 误判 TS 无效 | 高 | irc_allow_intermediate 配置项；分步反应默认启用 |
| 新增 reaction_profiles 配置段与未来 defaults.yaml 更新冲突 | 合并冲突 | 低 | reaction_profiles 独立成段，位于文件末尾，与现有段无交叉引用 |
| IntramolecularFragmenter 泛化后影响 [5+2] 片段切割精度 | S4 特征质量下降 | 中 | core_atom_count 默认值保持 5（[5+2]），确保未显式配置时行为不变 |

---

## 八、交付件清单

完成所有模块后，交付以下内容：

1. 修改后的源文件（按模块列出，共涉及约 18 个文件，新建 1 个文件）
2. 新增/修改的测试文件（每个模块至少一个测试文件）
3. 更新后的 defaults.yaml（含 reaction_profiles 段）
4. 更新后的 AGENTS.md 知识库（反映新模块和新数据流）
5. 全部测试通过的运行日志
6. 至少一个 [4+3] 反应的端到端运行示例输出

---

## 附录 A：已验证的关键技术细节

### xTB 原生扫描功能块语法

输入文件结构（以两键协同扫描为例）：

$constrain 块中写入 force constant 和两条 distance 约束（1-based 原子编号、目标距离）；$scan 块中写入 mode=concerted，以及两条扫描区间（起始距离, 终止距离, 步数）；$opt 块中写入 maxcycle。最后以 $end 结束。

执行命令与普通优化相同（xtb input.xyz --opt --input scan.inp），xTB 通过检测 $scan 块自动切换到扫描模式。

输出文件 xtbscan.log 为 XMol 多帧格式：每帧由原子数行、包含 "SCF done" 和能量值的注释行、以及所有原子的坐标行组成，逐帧排列。

### 已知限制

- concerted 模式要求所有扫描坐标的步数必须相同，否则 xTB 报错退出
- xTB 仅支持弛豫扫描（每个扫描点执行几何优化），不支持刚性扫描
- GitHub issue 1392 报告了二面角扫描在非相邻原子时可能回退到初始坐标的问题，但距离扫描不受影响
- xTB 扫描始终使用 GFN2-xTB 级别，无法指定为 GFN0 或 GFN-FF

### Cleaner 数据格式

core_bond_changes 字段格式为分号分隔的键变化描述，每项由 "原子映射编号1-原子映射编号2:变化类型" 组成，变化类型为 formed 或 broken。

例如 "1-2:formed;5-6:formed;3-4:broken" 表示原子映射编号 1-2 和 5-6 之间形成新键，3-4 之间断键。

仅当 map_status 为 OK（原子映射置信度≥0.8）时完全可信。LOW_CONFIDENCE 时可用但需交叉验证。FAILED 时不可用，需回退到 SMARTS 匹配。

### [4+3] 环加成计算化学参数参考（Houk 组）

- 几何优化：B3LYP/6-31G(d)
- 单点能：M06-2X/6-311+G(d,p)
- 溶剂模型：CPCM（CH₂Cl₂）
- TS 成键距离范围：1.8-2.2 Å，存在显著不对称性
- 电子结构：极性溶剂中为 zwitterion（非双自由基），标准 RB3LYP/M06-2X 即可处理
- 机理：分步，经由 zwitterion 中间体（Reactants → TS₁ → INT → TS₂ → Product）

---

## 附录 B：涉及 [5+2] 硬编码的完整文件清单

以下 7 个文件共 32 处引用了 [5+2] 专属逻辑，均需在对应模块中参数化：

1. rph_core/steps/step2_retro/smarts_matcher.py —— 第 5、129 行（M5 处理）
2. rph_core/steps/step3_opt/intramolecular_fragmenter.py —— 第 5、38、63、69、107、127、175、188-189、195、213、221、244、277、337、339、342、353、363-364 行（M6 处理）
3. rph_core/steps/step4_features/extractors/step2_cyclization.py —— 第 5、148 行（M7 处理）
4. rph_core/utils/semantic_slicer.py —— 第 6、66、213 行（M7 处理）
5. rph_core/utils/fragment_manipulation.py —— 第 159 行（M7 处理）
6. rph_core/steps/step4_features/fragment_extractor.py —— 第 239、342、359 行（M7 处理）
7. rph_core/steps/step3_opt/post_qc_enrichment.py —— 第 242 行（M6 处理）

---

## 附录 C：术语表

| 术语 | 含义 |
|------|------|
| 正向扫描（forward scan） | 从反应物/中间体出发，压缩成键键距搜索能量极大点作为 TS 猜测 |
| 逆向扫描（retro scan） | 从产物出发，拉长成键键距搜索 TS 猜测（现有 S2 策略） |
| 协同扫描（concerted scan） | xTB $scan 的 mode=concerted，所有约束坐标同步扫描 |
| 成键键（forming bond） | 反应中新形成的化学键 |
| 断键键（breaking bond） | 反应中断裂的化学键 |
| Curtin-Hammett 原理 | 快速平衡的中间体之间，产物分布由各自过渡态的能量差决定 |
| Zwitterion | 同时带有正负电荷的中性分子/中间体 |
| Oxyallyl | 氧杂烯丙基，[4+3] 环加成中 3 原子活性核心 |
| IRC（Intrinsic Reaction Coordinate） | 内禀反应坐标，验证 TS 连接的反应物和产物 |
| reaction_profile | defaults.yaml 中定义的反应类型配置参数集 |
| cleaner_adapter | 将 cleaner 导出数据转换为 RPH 内部格式的适配器模块 |
